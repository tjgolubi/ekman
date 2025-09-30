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

#include "geom.hpp"

#include <boost/geometry.hpp>
#include <gsl/gsl>

#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <exception>
#include <limits>
#include <cmath>
#include <cstdlib>

namespace bg = boost::geometry;
namespace fs = std::filesystem;

using geom::Vec;
using geom::Pi;
using geom::TwoPi;
using geom::DegPerRad;
using geom::RadPerDeg;

using Pt       = geom::Pt;
using Polygon  = bg::model::polygon<Pt>;
using MP       = bg::model::multi_polygon<Polygon>;
using Ring     = bg::model::ring<Pt>;

#if 0
namespace Tunables {

constexpr double MiterLimit   = 4.0;
constexpr int    CirclePoints = 14;

// Geometry cleanup (applied to geometry we will drive)
constexpr double SimplifyInputTol    = 0.01; // m; <=0 → skip
constexpr double SimplifyInsetTol    = 0.01; // m; <=0 → skip
constexpr double SimplifySmoothedTol = 0.04; // try 0.10–0.30

// Corner detection (on ORIGINAL perimeter via DP 0.5 m)
constexpr double SimplifyForCorners = 0.50; // m; Douglas–Peucker for corner detection
constexpr double CornerAngleDeg     = 45.0; // deflection threshold
constexpr double MinSegLen          = 0.50; // m; ignore tiny edges near vertex

// Smoothing between mapped corner indices on the inset ring
constexpr double ResampleStep       = 1.00; // m
constexpr int    HeadingHalfWin     = 3;    // moving-average half-window (2*H+1)

// Safety corridor (optional clamp). 0.0 disables.
constexpr double CorridorShrink     = 0.0;  // m
} // Tunables
#endif

namespace {

std::ostream& operator<<(std::ostream& os, const Pt& p)
  { return os << p.x << ' ' << p.y; }

std::istream& operator>>(std::istream& os, Pt& p)
  { return os >> p.x >> p.y; }

Ring EnsureClosed(const Ring& ring) {
  if (ring.empty() || ring.front() == ring.back()) return ring;
  auto out = ring;
  out.push_back(ring.front());
  return out;
}

std::vector<Pt> ReadPoints(const fs::path& path) {
  auto in = std::ifstream{path};
  if (!in) throw std::runtime_error{"cannot open input: " + path.string()};
  std::vector<Pt> pts;
  Pt pt;
  while (in >> pt) pts.push_back(pt);
  if (pts.size() < 3)
    throw std::runtime_error{"not enough points in: " + path.string()};
  in.close();
  return pts;
} // ReadPoints

void WritePoints(const std::vector<Pt>& pts, const fs::path& path) {
  auto out = std::ofstream{path};
  if (!out) throw std::runtime_error{"cannot open output: " + path.string()};
  out << std::fixed << std::setprecision(6);
  for (const auto& p : pts) out << p << '\n';
  out.close();
} // WritePoints

Polygon ReadPolygon(const fs::path& path) {
  auto in = std::ifstream{path};
  if (!in) throw std::runtime_error{"cannot open input: " + path.string()};
  auto poly = Polygon{};
  auto* ring = &poly.outer();
  int count = 0;
  int inner = 0;
  std::vector<Pt> pts;
  auto line = std::string{};
  Pt pt;
  while (getline(in, line)) {
    if (line.empty())
      continue;
    if (line != "# inner") {
      std::istringstream iss{line};
      iss >> pt;
      ring->push_back(pt);
      ++count;
      continue;
    }
    // Start a new inner ring.
    poly.inners().emplace_back();
    ring = &poly.inners().back();
    ++inner;
  }
  bg::correct(poly);
  in.close();
  std::cerr << "Read " << count << " points, " << inner << "inners\n";
  return poly;
} // ReadPolygon

void WritePolygon(const Polygon& poly, const fs::path& path) {
    auto out = std::ofstream{path};
    if (!out) throw std::runtime_error{"cannot open output: " + path.string()};
    out << std::fixed << std::setprecision(6);
    for (const auto& p : poly.outer()) out << p << '\n';
    for (const auto& v : poly.inners()) {
      out << "# inner\n";
      for (const auto& p : v) out << p << '\n';
    }
    out.close();
} // WritePolygon

// Build polygon from (possibly open) points
Polygon MakePolygon(std::vector<Pt> pts) {
  auto ring = Ring{pts.begin(), pts.end()};
  auto poly = Polygon{ring};
  bg::correct(poly);
  Expects(poly.outer().front() == poly.outer().back());
  return poly;
}

void WriteDeflections(const Ring& ring, const fs::path& path) {
  Expects(ring.size() >= 3 && ring.front() == ring.back());
  auto out = std::ofstream{path};
  if (!out)
    throw std::runtime_error{"cannot open deflection file: " + path.string()};
  out << std::fixed << std::setprecision(6);
  auto n = std::ssize(ring) - 1;
  auto curr = ring[0] - ring[n-1];
  for (auto i = gsl::index{0}; i != n; ++i) {
    auto prev = curr;
    curr = ring[i+1] - ring[i];
    auto deg  = curr.angle_wrt(prev) * DegPerRad;
    out << ring[i] << ' ' << deg << '\n';
  }
} // WriteDeflections

#if 0
void WriteCorners(const Ring& ring_unique,
                  const std::vector<gsl::index>& idx,
                  const fs::path& path = "corners.xy")
{
  auto out = std::ofstream{path};
  if (!out)
    throw std::runtime_error{"cannot open corners file: " + path.string()};
  out.setf(std::ios::fixed);
  out << std::setprecision(6);
  for (auto i : idx) {
    const auto& c = ring_unique.at(i % ring_unique.size());
    out << c << '\n';
  }
} // WriteCorners

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
      const auto d2 = Dist2(Rorig[i], q);
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
#endif

} // anonymous

int main(int argc, const char* argv[]) {
  try {
    auto arg0 = fs::path{argv[0]};

    if (argc != 4) {
      std::cerr << "usage:\n" << arg0.string()
        << "  <input.xy> <offset_m> <output.xy>\n";
      return EXIT_FAILURE;
    }

    const auto in_path  = fs::path{argv[1]};
    const auto offset   = std::stod(argv[2]); // inward distance (m)
    const auto out_path = fs::path{argv[3]};

    if (offset <= 0.0) throw std::runtime_error{"<offset_m> must be > 0"};

    // ---- Read ORIGINAL perimeter and build polygon
    auto poly_in = ReadPolygon(in_path);

    // ---- NEW: Write deflection angles (degrees) for simplified perimeter
    WriteDeflections(poly_in.outer(), "deflect.txt");

    WritePolygon(poly_in, out_path);

#if 0
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
    WriteDeflections(Rsimp, "deflect.txt");

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
#endif
  }
  catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
} // main
