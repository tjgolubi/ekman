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
//   inset_bg <input.xy> <offset_m> <output.xy>
//
// Input:  one "x y" per line; not necessarily closed.
// Output: one "x y" per line; CLOSED ring (last == first).

#include <boost/geometry.hpp>

#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <exception>
#include <limits>
#include <cmath>
#include <cstdlib>

namespace bg = boost::geometry;
namespace fs = std::filesystem;

using Pt       = bg::model::d2::point_xy<double>;
using Polygon  = bg::model::polygon<Pt>;
using MP       = bg::model::multi_polygon<Polygon>;

// -------------------- TUNABLES --------------------
namespace Tunables {

constexpr double kMiterLimit   = 4.0;
constexpr int    kCirclePoints = 14;

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
} // Tunables
// -------------------------------------------------

namespace {

constexpr auto Pi = M_PI;
constexpr auto TwoPi = 2 * Pi;
constexpr auto DegPerRad = 180.0 / Pi;
constexpr auto RadPerDeg = Pi / 180.0;

std::vector<Pt> EnsureClosed(const std::vector<Pt>& ring) {
  if (ring.empty()) return ring;
  if (ring.front().x() == ring.back().x()
   && ring.front().y() == ring.back().y())
    return ring;
  auto out = ring;
  out.push_back(ring.front());
  return out;
}

std::vector<Pt> ReadPoints(const fs::path& path) {
  auto in = std::ifstream{path};
  if (!in) throw std::runtime_error{"cannot open input: " + path.string()};
  std::vector<Pt> pts;
  double x, y;
  while (in >> x >> y) pts.emplace_back(x, y);
  if (pts.size() < 3)
    throw std::runtime_error{"not enough points in: " + path.string()};
  return pts;
} // ReadPoints

void WritePoints(const std::vector<Pt>& ring, const fs::path& path) {
  const auto R = EnsureClosed(ring);
  if (R.size() < 4)
    throw std::runtime_error{"ring degenerate (fewer than 3 unique points)"};
  auto out = std::ofstream{path};
  if (!out) throw std::runtime_error{"cannot open output: " + path.string()};
  out.setf(std::ios::fixed);
  out << std::setprecision(6);
  for (auto const& p : R) out << p.x() << ' ' << p.y() << '\n';
} // WritePoints

void WriteCorners(const std::vector<Pt>& ring_unique,
                  const std::vector<std::size_t>& idx,
                  const fs::path& path = "corners.xy")
{
  auto out = std::ofstream{path};
  if (!out)
    throw std::runtime_error{"cannot open corners file: " + path.string()};
  out.setf(std::ios::fixed);
  out << std::setprecision(6);
  for (auto i : idx) {
    const auto& c = ring_unique.at(i % ring_unique.size());
    out << c.x() << ' ' << c.y() << '\n';
  }
} // WriteCorners

// Build polygon from (possibly open) points
Polygon MakePolygon(std::vector<Pt> pts) {
  auto poly = Polygon{};
  poly.outer() = std::move(pts);
  bg::correct(poly);
  return poly;
}

const Polygon& BiggestByArea(const MP& mp) {
  if (mp.empty()) throw std::runtime_error{"buffer produced no polygons (collapsed?)"};
  const Polygon* best = &mp.front();
  auto best_area = std::abs(bg::area(*best));
  for (const auto& poly : mp) {
    auto area = std::abs(bg::area(poly));
    if (area > best_area) { best = &poly; best_area = area; }
  }
  return *best;
}

struct Dir {
  double dx = 0.0;
  double dy = 0.0;
}; // Dir

struct Vec {
  double dx=0.0;
  double dy=0.0;
  constexpr Vec() noexcept = default;
  constexpr Vec(const Vec&) noexcept = default;
  constexpr Vec& operator=(const Vec&) noexcept = default;
  constexpr const Vec& operator+() const noexcept { return *this; }
  constexpr Vec operator-() const noexcept { return Vec{-dx, -dy}; }
  constexpr Vec& operator+=(const Vec& rhs) noexcept
    { dx += rhs.dx; dy += rhs.dy; return *this; }
  constexpr Vec& operator-=(const Vec& rhs) noexcept
    { dx -= rhs.dx; dy -= rhs.dy; return *this; }
  constexpr Vec& operator*=(double s) noexcept
    { dx *= s; dy *= s; return *this; }
  constexpr Vec& operator/=(double s) noexcept
    { dx /= s; dy /= s; return *this; }
  constexpr Vec operator+(const Vec& rhs) const noexcept
    { return Vec{dx+rhs.dx, dy+rhs.dy}; }
  constexpr Vec operator-(const Vec& rhs) const noexcept
    { return Vec{dx-rhs.dx, dy-rhs.dy}; }
  constexpr Vec operator*(double s) const noexcept
    { return Vec{dx*s, dy*s}; }
  constexpr Vec operator/(double s) const noexcept
    { return Vec{dx/s, dy/s}; }
  constexpr friend Vec operator*(double s, const Vec& rhs)
    { return rhs * s; }
  constexpr friend double dot(const Vec& u, const Vec& v)
    { return (u.dx * v.dx + u.dy * v.dy); }
  constexpr friend double operator*(const Vec& u, const Vec& v)
    { return dot(u, v); }
  constexpr friend double cross(const Vec& u, const Vec& v)
    { return (u.dx * v.dy - u.dy * v.dx); }
  constexpr double norm2() const { return dx*dx + dy*dy; }
  constexpr double norm() const { return std::hypot(dx, dy); }
  constexpr Vec unit() const { return *this / norm(); }
  constexpr double angle_wrt(const Vec& ref) const
    { return std::atan2(cross(ref, *this), dot(ref, *this)); }
}; // Vec

constexpr Vec operator-(const Pt& lhs, const Pt& rhs) noexcept
  { return Vec{lhs.x-rhs.x, lhs.y-rhs.y}; }

constexpr Pt operator+(const Pt& p, const Vec& v) noexcept
  { return Pt{p.x() + v.dx, p.y() + v.dy}; }

constexpr Pt operator-(const Pt& p, const Vec& v) noexcept
  { return Pt{p.x() - v.dx, p.y() - v.dy}; }

constexpr double Dist(const Pt& a, const Pt& b) noexcept
  { return (a - b).norm(); }

constexpr double Dist2(const Pt& a, const Pt& b) noexcept
  { return (a - b).norm2(); }

std::vector<gsl::index>
FindCorners(const Ring& R, double ang_deg, double min_seg_len) {
  std::vector<gsl::index> idx;
  const auto n = std::ssize(R);
  if (n < 3) return idx;
  const double th = ang_deg * RadPerDeg;
  auto next = R[0] - R[n-2]; // because R[n-1] == R[0] for closed ring
  for (gsl::index i = 0; i != n-1; ++i) {
    auto curr = next;
    next = R[i+1] - R[i];
    if (std::abs(next.angle_wrt(curr)) >= th) idx.push_back(i);
  }
  return idx;
} // FindCorners

// Map simplified-corner points to indices in the ORIGINAL ring (ordered match)
std::vector<gsl::index>
MapCornersToOriginal(const Ring& Rorig,
                     const Ring& Rsimp,
                     const std::vector<gsl::index>& simp_corners,
                     double tol /* meters */)
{
  std::vector<gsl::index> out;
  out.reserve(simp_corners.size());
  if (Rorig.empty() || Rsimp.empty() || simp_corners.empty()) return out;

  const double tol2 = tol * tol;
  auto i0 = gsl::index{0}; // rolling pointer in original ring (unique vertices)

  for (gsl::index k = 0; k != std::ssize(simp_corners); ++k) {
    const Pt& q = Rsimp[simp_corners[k]];
    auto best_i = i0;
    double best_d2 = std::numeric_limits<double>::infinity();

    // scan forward around the ring once
    const auto n = std::ssize(Rorig);
    auto i = i0;
    for (gsl::index t = 0; t != n; ++t) {
      const double d2 = Dist2(Rorig[i], q);
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
std::vector<std::size_t>
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
      const double d2 = Dist2(Rinset[j], q);
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
void UnwrapAngles(std::vector<double>& th) {
  for (gsl::index i = 1; i != std::ssize(th.size(); ++i) {
    double d = th[i] - th[i-1];
    while (d >  Pi) { th[i] -= TwoPi; d -= TwoPi; }
    while (d < -Pi) { th[i] += TwoPi; d += TwoPi; }
  }
}

std::vector<double>
MovingAverageCircular(const Ring& v, int half_win) {
  if (half_win <= 0) return v;
  const int W = 2*half_win + 1;
  const int n = static_cast<int>(v.size());
  std::vector<double> out(n);
  for (int i = 0; i != n; ++i) {
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
} // MovingAverageCircular

std::vector<Pt> ResamplePolyline(const std::vector<Pt>& seg_pts, double ds) {
  // seg_pts: OPEN polyline (A..B). Returns resampled OPEN list including endpoints.
  std::vector<double> S(seg_pts.size(), 0.0);
  for (std::size_t i = 1; i < seg_pts.size(); ++i)
    S[i] = S[i-1] + Dist(seg_pts[i-1], seg_pts[i]);
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

void ClampToCorridor(std::vector<Pt>& P, const Polygon& corridor) {
  for (auto& q : P) {
    if (!bg::within(q, corridor)) {
      bg::model::segment<Pt> seg;
      bg::closest_points(corridor, q, seg);
      q = seg.first;
    }
  }
}

std::vector<Pt>
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
    const double step = Dist(R[i-1], R[i]);
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
    if (clean.empty() || Dist(clean.back(), p) > 1e-6) clean.push_back(p);
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
std::vector<Pt>
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
void WriteDeflectionsDeg(const Ring& R, const fs::path& path) {
  auto out = std::ofstream{path};
  if (!out)
    throw std::runtime_error{"cannot open deflection file: " + path.string()};
  out.setf(std::ios::fixed);
  out << std::setprecision(6);
  const auto n = std::ssize(R);
  if (n < 3) {
    for (std::size_t i = 0; i < n; ++i) out << 0.0 << '\n';
    return;
  }
  auto curr = R
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t ip = (i + n - 1) % n, in = (i + 1) % n;
    const double ang = Deflection(R[ip], R[i], R[in]); // radians
    const double deg = ang * DegPerRad;
    out << R[i].x() << ' ' << R[i].y() << ' ' << deg << '\n';
  }
}

} // anonymous

int main(int argc, const char* argv[]) {
  try {
    if (argc != 4) {
      throw std::runtime_error{
        "usage:\n"
        "  inset_bg <input.xy> <offset_m> <output.xy>\n"
      };
    }

    const auto in_path  = fs::path{argv[1]};
    const auto offset   = std::stod(argv[2]); // inward distance (m)
    const auto out_path = fs::path{argv[3]};

    if (offset <= 0.0) throw std::runtime_error{"<offset_m> must be > 0"};

    using namespace Tunables;

    // ---- Read ORIGINAL perimeter and build polygon
    auto poly_in = MakePolygon(ReadPoints(in_path));

    // Optional light cleanup for the path we’ll buffer/drive
    if (kSimplifyInputTol > 0.0) {
      auto simp = Polygon{};
      bg::simplify(poly_in, simp, kSimplifyInputTol);
      poly_in = std::move(simp);
    }

    // ---- Corner detection on ORIGINAL perimeter via DP(0.5 m)
    auto poly_dp = Polygon{};
    bg::simplify(poly_in, poly_dp, kSimplifyForCorners);

    // Extract unique vertices from original & simplified
    auto Rorig = std::vector<Pt>{};
    auto Rsimp = std::vector<Pt>{};
    Rorig.reserve(poly_in.outer().size());
    for (std::size_t i = 0; i + 1 < poly_in.outer().size(); ++i)
      Rorig.push_back(poly_in.outer()[i]);

    Rsimp.reserve(poly_dp.outer().size());
    for (std::size_t i = 0; i + 1 < poly_dp.outer().size(); ++i)
      Rsimp.push_back(poly_dp.outer()[i]);

    // ---- NEW: Write deflection angles (degrees) for simplified perimeter
    WriteDeflectionsDeg(Rsimp, "deflect.txt");

    // Indices on simplified ring (few points)
    auto simp_corners = FindCorners(Rsimp, kCornerAngleDeg, kMinSegLen);
    // Map those to indices on the ORIGINAL ring
    auto orig_corner_idx =
          MapCornersToOriginal(Rorig, Rsimp, simp_corners, kSimplifyForCorners);

    // ---- Inset buffer (negative distance) ----
    auto distance = bg::strategy::buffer::distance_symmetric<double>{-offset};
    auto side     = bg::strategy::buffer::side_straight;
    auto join_m   = bg::strategy::buffer::join_miter{kMiterLimit};
    auto end      = bg::strategy::buffer::end_flat;
    auto circle   = bg::strategy::buffer::point_circle{kCirclePoints};

    auto inset_mp = MP{};
    bg::buffer(poly_in, inset_mp, distance, side, join_m, end, circle);
    const Polygon& inset_poly0 = BiggestByArea(inset_mp);
    Polygon inset_poly = inset_poly0;

    if (kSimplifyInsetTol > 0.0) {
      auto simp = Polygon{};
      bg::simplify(inset_poly, simp, kSimplifyInsetTol);
      bg::correct(simp);
      inset_poly = std::move(simp);
      bg::unique(inset_poly);
    }

    // Extract INSET unique ring
    auto Rinset = std::vector<Pt>{};
    Rinset.reserve(inset_poly.outer().size());
    for (std::size_t i = 0; i + 1 < inset_poly.outer().size(); ++i)
      Rinset.push_back(inset_poly.outer()[i]);

    // Map ORIGINAL-corner positions to INSET ring indices (ordered nearest)
    auto inset_corner_idx =
                      MapOriginalCornersToInset(Rinset, Rorig, orig_corner_idx);

    // ---- Optional corridor (disabled here if kCorridorShrink == 0)
    auto corridor = std::optional<Polygon>{};
    if (kCorridorShrink > 0.0) {
      MP corr_mp;
      auto dist2 =
            bg::strategy::buffer::distance_symmetric<double>{-kCorridorShrink};
      auto j2 = bg::strategy::buffer::join_miter{kMiterLimit};
      bg::buffer(inset_poly, corr_mp, dist2, side, j2, end, circle);
      corridor = BiggestByArea(corr_mp);
    }

    // ---- Smooth BETWEEN those inset corners
    auto smoothed_closed = SmoothBetweenCorners(Rinset, inset_corner_idx,
              kResampleStep, kHeadingHalfWin, corridor ? &*corridor : nullptr);

    // ---- Emit files
    WritePoints(smoothed_closed, out_path);
    WriteCorners(Rinset, inset_corner_idx, "corners.xy");

    {
      using namespace std;
      cout << "Simplified-for-corners vertices: " << Rsimp.size()
        << "\nDeflections wrote to deflect.txt (degrees, one per vertex)."
         "\nCorners (>= " << kCornerAngleDeg << " deg): "
        << inset_corner_idx.size()
        << "\nWrote smoothed inset: " << out_path << "  and corners.xy\n";
    }
  }
  catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} // main
