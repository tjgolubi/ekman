// inset_smooth_corners.cpp
//
// Read XY perimeter (meters), build an inward inset path, detect "real" corners
// by simplifying the ORIGINAL perimeter (0.5 m DP), map those corners onto the
// INSET ring, smooth BETWEEN those corners (endpoints anchored) with a heading
// moving-average, and write:
//   - <output.xy> : smoothed inset ring (CLOSED)
//   - corners.xy  : corner points on the inset ring (for plotting)
//   - deflect.txt : deflection angle (degrees) for each vertex of the
//                   simplified-for-corners perimeter (one value per line)
//
// Usage:
//   inset_smooth_corners <input.xy> <offset_m> <output.xy> [--join miter|round] [--miter 4.0] [--circle 16]
//
// Input:  one "x y" per line; not necessarily closed.
// Output: one "x y" per line; CLOSED ring (last == first).

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/simplify.hpp>
#include <boost/geometry/algorithms/unique.hpp>
#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/algorithms/closest_points.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace bg = boost::geometry;

using Pt       = bg::model::d2::point_xy<double>;
using Polygon  = bg::model::polygon<Pt, /*CW=*/true, /*Closed=*/true>;
using MP       = bg::model::multi_polygon<Polygon>;

// -------------------- TUNABLES --------------------
namespace Tunables {
// Geometry cleanup (applied to geometry we will drive)
constexpr double kSimplifyInputTol    = 0.01; // m; <=0 → skip
constexpr double kSimplifyInsetTol    = 0.01; // m; <=0 → skip
constexpr double kSimplifySmoothedTol = 0.04; // try 0.10–0.30

// Corner detection (on ORIGINAL perimeter via DP 0.5 m)
constexpr double kSimplifyForCorners = 0.50; // m; Douglas–Peucker for corner detection
constexpr double kCornerAngleDeg     = 45.0; // deflection threshold
constexpr double kMinSegLen          = 0.50; // m; ignore tiny edges near vertex

// Smoothing between mapped corner indices on the inset ring
constexpr double kResampleStep       = 1.00; // m
constexpr int    kHeadingHalfWin     = 3;    // moving-average half-window (2*H+1)

// Safety corridor (optional clamp). 0.0 disables.
constexpr double kCorridorShrink     = 0.0;  // m
} // namespace Tunables
// -------------------------------------------------

static double SegLen(const Pt& a, const Pt& b) {
  return std::hypot(b.x() - a.x(), b.y() - a.y());
}

static std::vector<Pt> EnsureClosed(const std::vector<Pt>& ring) {
  if (ring.empty()) return ring;
  if (ring.front().x() == ring.back().x() && ring.front().y() == ring.back().y())
    return ring;
  auto out = ring;
  out.push_back(ring.front());
  return out;
}

static std::vector<Pt> ReadPoints(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open input: " + path);
  std::vector<Pt> pts;
  double x{}, y{};
  while (in >> x >> y) pts.emplace_back(x, y);
  if (pts.size() < 3) throw std::runtime_error("not enough points in: " + path);
  return pts;
}

static void WriteClosed(const std::vector<Pt>& ring, const std::string& path) {
  const auto R = EnsureClosed(ring);
  if (R.size() < 4) throw std::runtime_error("ring degenerate (fewer than 3 unique points)");
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open output: " + path);
  out.setf(std::ios::fixed);
  out << std::setprecision(6);
  for (auto const& p : R) out << p.x() << ' ' << p.y() << '\n';
}

static void WriteCornersXY(const std::vector<Pt>& ring_unique,
                           const std::vector<std::size_t>& idx,
                           const std::string& path = "corners.xy") {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open corners file: " + path);
  out.setf(std::ios::fixed);
  out << std::setprecision(6);
  for (auto i : idx) {
    const auto& c = ring_unique.at(i % ring_unique.size());
    out << c.x() << ' ' << c.y() << '\n';
  }
}

// Build polygon from (possibly open) points
static Polygon MakePolygonFromPoints(const std::vector<Pt>& pts) {
  Polygon poly;
  auto& ring = poly.outer();
  ring.clear();
  ring.reserve(pts.size() + 1);
  for (const auto& p : pts) ring.push_back(p);
  if (ring.size() >= 2 && !bg::equals(ring.front(), ring.back()))
    ring.push_back(ring.front());
  bg::correct(poly);
  return poly;
}

static const Polygon& BiggestByArea(const MP& mp) {
  if (mp.empty()) throw std::runtime_error("buffer produced no polygons (collapsed?)");
  const Polygon* best = &mp.front();
  auto best_area = std::abs(bg::area(*best));
  for (const auto& poly : mp) {
    auto area = std::abs(bg::area(poly));
    if (area > best_area) { best = &poly; best_area = area; }
  }
  return *best;
}

static double Deflection(const Pt& A, const Pt& B, const Pt& C) {
  // Use AB and BC (forward directions), not BA.
  const double ux = B.x() - A.x(), uy = B.y() - A.y(); // AB
  const double vx = C.x() - B.x(), vy = C.y() - B.y(); // BC
  const double nu = std::hypot(ux, uy), nv = std::hypot(vx, vy);
  if (nu == 0.0 || nv == 0.0) return 0.0;
  const double dot = (ux*vx + uy*vy) / (nu*nv);
  return std::acos(std::clamp(dot, -1.0, 1.0)); // 0..π
}

static std::vector<std::size_t>
CornerIdxOnSimplified(const std::vector<Pt>& S,
                      double ang_deg,
                      double min_seg_len)
{
  std::vector<std::size_t> idx;
  const std::size_t n = S.size();
  if (n < 3) return idx;
  const double th = ang_deg * M_PI / 180.0;
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t ip = (i + n - 1) % n, in = (i + 1) % n;
    if (SegLen(S[ip], S[i]) < min_seg_len || SegLen(S[i], S[in]) < min_seg_len) continue;
    if (Deflection(S[ip], S[i], S[in]) >= th) idx.push_back(i);
  }
  return idx;
}

static inline double D2(const Pt& a, const Pt& b) {
  const double dx = a.x() - b.x(), dy = a.y() - b.y();
  return dx*dx + dy*dy;
}

// Map simplified-corner points to indices in the ORIGINAL ring (ordered match)
static std::vector<std::size_t>
MapCornersToOriginal(const std::vector<Pt>& Rorig,
                     const std::vector<Pt>& Rsimp,
                     const std::vector<std::size_t>& simp_corners,
                     double tol /* meters */)
{
  std::vector<std::size_t> out; out.reserve(simp_corners.size());
  if (Rorig.empty() || Rsimp.empty() || simp_corners.empty()) return out;

  const double tol2 = tol * tol;
  std::size_t i0 = 0; // rolling pointer in original ring (unique vertices)

  for (std::size_t k = 0; k < simp_corners.size(); ++k) {
    const Pt& q = Rsimp[simp_corners[k]];
    std::size_t best_i = i0;
    double best_d2 = std::numeric_limits<double>::infinity();

    // scan forward around the ring once
    const std::size_t n = Rorig.size();
    std::size_t i = i0;
    for (std::size_t t = 0; t < n; ++t) {
      const double d2 = D2(Rorig[i], q);
      if (d2 < best_d2) { best_d2 = d2; best_i = i; }
      if (best_d2 <= tol2 && t > 5) break;
      i = (i + 1) % n;
    }
    out.push_back(best_i);
    i0 = (best_i + 1) % Rorig.size();
  }

  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// Map ORIGINAL-corner coordinates to nearest indices on the INSET ring (ordered)
static std::vector<std::size_t>
MapOriginalCornersToInset(const std::vector<Pt>& Rinset,
                          const std::vector<Pt>& Rorig,
                          const std::vector<std::size_t>& orig_corner_idx)
{
  std::vector<std::size_t> out; out.reserve(orig_corner_idx.size());
  if (Rinset.empty() || Rorig.empty() || orig_corner_idx.empty()) return out;

  std::size_t j0 = 0; // rolling pointer along inset ring
  const std::size_t m = Rinset.size();
  for (auto oi : orig_corner_idx) {
    const Pt& q = Rorig[oi];
    std::size_t best_j = j0;
    double best_d2 = std::numeric_limits<double>::infinity();

    // scan forward once around inset ring
    std::size_t j = j0;
    for (std::size_t t = 0; t < m; ++t) {
      const double d2 = D2(Rinset[j], q);
      if (d2 < best_d2) { best_d2 = d2; best_j = j; }
      j = (j + 1) % m;
    }
    out.push_back(best_j);
    j0 = (best_j + 1) % m;
  }

  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// ---------- Heading smoothing (between corners) ----------
static void UnwrapAngles(std::vector<double>& th) {
  for (std::size_t i = 1; i < th.size(); ++i) {
    double d = th[i] - th[i-1];
    while (d >  M_PI) { th[i] -= 2.0*M_PI; d -= 2.0*M_PI; }
    while (d < -M_PI) { th[i] += 2.0*M_PI; d += 2.0*M_PI; }
  }
}

static std::vector<double> MovingAverageCircular(const std::vector<double>& v, int half_win) {
  if (half_win <= 0) return v;
  const int W = 2*half_win + 1;
  const int n = static_cast<int>(v.size());
  std::vector<double> out(n);
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int k = -half_win; k <= half_win; ++k) {
      int j = i + k;
      if (j < 0) j += n;
      if (j >= n) j -= n;
      s += v[j];
    }
    out[i] = s / W;
  }
  return out;
}

static std::vector<Pt> ResamplePolyline(const std::vector<Pt>& seg_pts, double ds) {
  // seg_pts: OPEN polyline (A..B). Returns resampled OPEN list including endpoints.
  std::vector<double> S(seg_pts.size(), 0.0);
  for (std::size_t i = 1; i < seg_pts.size(); ++i)
    S[i] = S[i-1] + SegLen(seg_pts[i-1], seg_pts[i]);
  const double L = S.back();
  if (L == 0.0) return {seg_pts.front(), seg_pts.back()};
  const int N = std::max(2, static_cast<int>(std::round(L / ds)) + 1);
  const double step = L / (N - 1);

  std::vector<Pt> R; R.reserve(N);
  R.push_back(seg_pts.front());
  std::size_t idx = 0;
  for (int k = 1; k < N-1; ++k) {
    const double s = k * step;
    while (idx + 1 < S.size() && S[idx+1] < s) ++idx;
    const double u = (s - S[idx]) / (S[idx+1] - S[idx] + 1e-12);
    const auto& A = seg_pts[idx];
    const auto& B = seg_pts[idx+1];
    R.emplace_back((1-u)*A.x() + u*B.x(), (1-u)*A.y() + u*B.y());
  }
  R.push_back(seg_pts.back());
  return R;
}

static void ClampToCorridor(std::vector<Pt>& P, const Polygon& corridor) {
  for (auto& q : P) {
    if (!bg::within(q, corridor)) {
      bg::model::segment<Pt> seg;
      bg::closest_points(corridor, q, seg);
      q = seg.first;
    }
  }
}

static std::vector<Pt>
SmoothSectionAnchored(const std::vector<Pt>& section_open,
                      double ds, int half_win,
                      const Polygon* corridor /*nullable*/)
{
  if (section_open.size() <= 2) return section_open; // nothing to smooth
  auto R = ResamplePolyline(section_open, ds);       // includes endpoints

  // Headings on internal steps
  std::vector<double> theta; theta.reserve(R.size()-1);
  for (std::size_t i = 0; i + 1 < R.size(); ++i)
    theta.push_back(std::atan2(R[i+1].y() - R[i].y(), R[i+1].x() - R[i].x()));
  UnwrapAngles(theta);

  if (half_win > 0 && theta.size() > static_cast<size_t>(2*half_win+1))
    theta = MovingAverageCircular(theta, half_win);

  // Reconstruct, keeping endpoints fixed
  std::vector<Pt> out(R.size());
  out.front() = R.front();
  for (std::size_t i = 1; i < R.size(); ++i) {
    const double step = SegLen(R[i-1], R[i]);
    const double t = theta[i-1];
    out[i] = Pt(out[i-1].x() + step * std::cos(t),
                out[i-1].y() + step * std::sin(t));
  }
  out.back() = R.back();

  // Optional corridor clamp
  if (corridor) ClampToCorridor(out, *corridor);

  // Drop exact duplicates from numerical noise
  std::vector<Pt> clean; clean.reserve(out.size());
  for (auto const& p : out) {
    if (clean.empty() || SegLen(clean.back(), p) > 1e-6) clean.push_back(p);
  }

  // ---- NEW: simplify this smoothed section as a linestring (endpoints preserved)
  if (Tunables::kSimplifySmoothedTol > 0.0 && clean.size() > 2) {
    using LineString = bg::model::linestring<Pt>;
    LineString ls;
    ls.assign(clean.begin(), clean.end());
    LineString ls_simpl;
    bg::simplify(ls, ls_simpl, Tunables::kSimplifySmoothedTol);
    // ensure at least endpoints remain (DP should, but guard tiny/degenerate cases)
    if (ls_simpl.size() >= 2) {
      std::vector<Pt> simplified;
      simplified.assign(ls_simpl.begin(), ls_simpl.end());
      return simplified;
    }
  }

  return clean;
}

// Assemble full smoothed ring by smoothing each corner-to-corner section (cyclic)
// ring_unique_inset: unique vertices (no duplicate last)
static std::vector<Pt>
SmoothBetweenCorners(const std::vector<Pt>& ring_unique_inset,
                     const std::vector<std::size_t>& corners_on_inset,
                     double ds, int half_win,
                     const Polygon* corridor /*nullable*/)
{
  if (ring_unique_inset.size() < 3 || corners_on_inset.size() < 2)
    return EnsureClosed(ring_unique_inset);

  std::vector<std::size_t> C = corners_on_inset;
  std::sort(C.begin(), C.end());
  C.erase(std::unique(C.begin(), C.end()), C.end());

  std::vector<Pt> out;
  for (std::size_t k = 0; k < C.size(); ++k) {
    const std::size_t i0 = C[k];
    const std::size_t i1 = C[(k + 1) % C.size()];

    std::vector<Pt> sec;
    sec.push_back(ring_unique_inset[i0]);
    std::size_t i = i0;
    while (i != i1) {
      i = (i + 1) % ring_unique_inset.size();
      sec.push_back(ring_unique_inset[i]);
    }

    auto smooth = (sec.size() <= 2)
      ? sec
      : SmoothSectionAnchored(sec, ds, half_win, corridor);

    if (k == 0) {
      out.insert(out.end(), smooth.begin(), smooth.end());
    } else {
      out.insert(out.end(), smooth.begin() + 1, smooth.end());
    }
  }

  out = EnsureClosed(out);
  return out;
}

// ---- NEW: Write per-vertex deflection (degrees) for the simplified perimeter
static void WriteDeflectionsDeg(const std::vector<Pt>& S, const std::string& path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open deflection file: " + path);
  out.setf(std::ios::fixed);
  out << std::setprecision(6);
  const std::size_t n = S.size();
  if (n < 3) {
    for (std::size_t i = 0; i < n; ++i) out << 0.0 << '\n';
    return;
  }
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t ip = (i + n - 1) % n, in = (i + 1) % n;
    const double ang = Deflection(S[ip], S[i], S[in]); // radians
    const double deg = ang * 180.0 / M_PI;
    out << S[i].x() << ' ' << S[i].y() << ' ' << deg << '\n';
  }
}

int main(int argc, char** argv) {
  try {
    if (argc < 4) {
      throw std::runtime_error(
        "usage:\n"
        "  inset_smooth_corners <input.xy> <offset_m> <output.xy>\n"
        "                       [--join miter|round] [--miter 4.0] [--circle 16]\n"
      );
    }
    const std::string in_path  = argv[1];
    const double offset        = std::stod(argv[2]); // inward distance (m)
    const std::string out_path = argv[3];

    std::string join = "miter";
    double miter_limit = 4.0;
    int circle_pts = 16;

    for (int i = 4; i < argc; ++i) {
      std::string key = argv[i];
      if (key == "--join" && i + 1 < argc) {
        join = std::string(argv[++i]);
        if (join != "miter" && join != "round") {
          throw std::runtime_error("invalid --join (miter|round)");
        }
      } else if (key == "--miter" && i + 1 < argc) {
        miter_limit = std::stod(argv[++i]);
        if (miter_limit <= 0.0) throw std::runtime_error("--miter must be > 0");
      } else if (key == "--circle" && i + 1 < argc) {
        circle_pts = std::stoi(argv[++i]);
        if (circle_pts < 4) throw std::runtime_error("--circle must be >= 4");
      } else {
        throw std::runtime_error("unknown or incomplete arg: " + key);
      }
    }
    if (offset <= 0.0) throw std::runtime_error("<offset_m> must be > 0");

    using namespace Tunables;

    // ---- Read ORIGINAL perimeter and build polygon
    const auto raw_pts = ReadPoints(in_path);
    Polygon poly_in = MakePolygonFromPoints(raw_pts);

    // Optional light cleanup for the path we’ll buffer/drive
    if (kSimplifyInputTol > 0.0) {
      Polygon simp; bg::simplify(poly_in, simp, kSimplifyInputTol); bg::correct(simp); poly_in = std::move(simp);
    }

    // ---- Corner detection on ORIGINAL perimeter via DP(0.5 m)
    Polygon poly_dp;
    bg::simplify(poly_in, poly_dp, kSimplifyForCorners);
    bg::correct(poly_dp);

    // Extract unique vertices from original & simplified
    std::vector<Pt> Rorig, Rsimp;
    Rorig.reserve(poly_in.outer().size());
    for (std::size_t i = 0; i + 1 < poly_in.outer().size(); ++i)
      Rorig.push_back(poly_in.outer()[i]);

    Rsimp.reserve(poly_dp.outer().size());
    for (std::size_t i = 0; i + 1 < poly_dp.outer().size(); ++i)
      Rsimp.push_back(poly_dp.outer()[i]);

    // ---- NEW: Write deflection angles (degrees) for simplified perimeter
    WriteDeflectionsDeg(Rsimp, "deflect.txt");

    // Indices on simplified ring (few points)
    auto simp_corners = CornerIdxOnSimplified(Rsimp, kCornerAngleDeg, kMinSegLen);
    // Map those to indices on the ORIGINAL ring
    auto orig_corner_idx = MapCornersToOriginal(Rorig, Rsimp, simp_corners, /*tol=*/kSimplifyForCorners);

    // ---- Inset buffer (negative distance) ----
    bg::strategy::buffer::distance_symmetric<double> distance(-offset);
    bg::strategy::buffer::side_straight              side;
    bg::strategy::buffer::end_flat                   end;
    bg::strategy::buffer::point_circle               circle(circle_pts);

    MP inset_mp;
    if (join == "miter") {
      bg::strategy::buffer::join_miter join_m(miter_limit);
      bg::buffer(poly_in, inset_mp, distance, side, join_m, end, circle);
    } else {
      bg::strategy::buffer::join_round join_r;
      bg::buffer(poly_in, inset_mp, distance, side, join_r, end, circle);
    }
    const Polygon& inset_poly0 = BiggestByArea(inset_mp);
    Polygon inset_poly = inset_poly0;

    if (kSimplifyInsetTol > 0.0) {
      Polygon simp; bg::simplify(inset_poly, simp, kSimplifyInsetTol); bg::correct(simp); inset_poly = std::move(simp);
      bg::unique(inset_poly);
    }

    // Extract INSET unique ring
    std::vector<Pt> Rinset;
    Rinset.reserve(inset_poly.outer().size());
    for (std::size_t i = 0; i + 1 < inset_poly.outer().size(); ++i)
      Rinset.push_back(inset_poly.outer()[i]);

    // Map ORIGINAL-corner positions to INSET ring indices (ordered nearest)
    auto inset_corner_idx = MapOriginalCornersToInset(Rinset, Rorig, orig_corner_idx);

    // ---- Optional corridor (disabled here if kCorridorShrink == 0)
    std::optional<Polygon> corridor;
    if (kCorridorShrink > 0.0) {
      MP corr_mp;
      bg::strategy::buffer::distance_symmetric<double> dist2(-kCorridorShrink);
      if (join == "miter") {
        bg::strategy::buffer::join_miter j2(miter_limit);
        bg::buffer(inset_poly, corr_mp, dist2, side, j2, end, circle);
      } else {
        bg::strategy::buffer::join_round j2;
        bg::buffer(inset_poly, corr_mp, dist2, side, j2, end, circle);
      }
      corridor = BiggestByArea(corr_mp);
    }

    // ---- Smooth BETWEEN those inset corners
    auto smoothed_closed = SmoothBetweenCorners(
      Rinset, inset_corner_idx,
      kResampleStep, kHeadingHalfWin,
      corridor ? &*corridor : nullptr);

    // ---- Emit files
    WriteClosed(smoothed_closed, out_path);
    WriteCornersXY(Rinset, inset_corner_idx, "corners.xy");

    std::cout << "Simplified-for-corners vertices: " << Rsimp.size() << "\n";
    std::cout << "Deflections wrote to deflect.txt (degrees, one per vertex).\n";
    std::cout << "Corners (>= " << kCornerAngleDeg << " deg): " << inset_corner_idx.size() << "\n";
    std::cout << "Wrote smoothed inset: " << out_path
              << "  and corners.xy\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
