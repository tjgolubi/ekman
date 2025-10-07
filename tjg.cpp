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

#include "Smooth.hpp"

#include "geom.hpp"

#include <boost/geometry.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <gsl/gsl>

#include <algorithm>
#include <string>
#include <vector>
#include <iterator>
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

template<typename T>
void DbgOut(const char* str, const T& x) {
  std::cerr << str << '=' << x << std::endl;
}

#define DBG(x) DbgOut(BOOST_STRINGIZE(x), (x))

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

using CornerVec = std::vector<gsl::index>;

namespace Tune {

constexpr double MiterLimit   = 4.0;
constexpr int    CirclePoints = 14;

// Geometry cleanup (applied to geometry we will drive)
constexpr double SimplifyInputTol    = 0.01; // m; <=0 → skip
constexpr double SimplifyInsetTol    = 0.01; // m; <=0 → skip
constexpr double SimplifySmoothedTol = 0.10; // try 0.10–0.30
constexpr double SimplifyOutputTol   = 0.10; // try 0.10–0.30

// Corner detection (on ORIGINAL perimeter via DP 0.5 m)
constexpr double SimplifyForCorners = 10.0; // m; Douglas–Peucker for corner detection
constexpr double CornerAngleDeg     = 45.0; // deflection threshold
constexpr double MinSegLen          = 0.50; // m; ignore tiny edges near vertex

// Smoothing between mapped corner indices on the inset ring
constexpr double ResampleStep       = 1.00; // m
constexpr int    HeadingHalfWin     = 3;    // moving-average half-window (2*H+1)

// Safety corridor (optional clamp). 0.0 disables.
constexpr double CorridorShrink     = 0.0;  // m

} // Tune

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
  out << std::fixed << std::setprecision(2);
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
  return poly;
} // ReadPolygon

void WritePolygon(const Polygon& poly, const fs::path& path) {
  auto out = std::ofstream{path};
  if (!out) throw std::runtime_error{"cannot open output: " + path.string()};
  out << std::fixed << std::setprecision(2);
  for (const auto& p : poly.outer()) out << p << '\n';
  out.close();
  if (poly.inners().empty())
    return;
  auto ext  = path.extension();
  int inner = 0;
  for (const auto& v : poly.inners()) {
    auto path2 = path;
    auto ext2 = fs::path{std::to_string(++inner)};
    ext2 += ext;
    path2.replace_extension(ext2);
    out = std::ofstream{path2};
    if (!out)
      throw std::runtime_error{"cannot open output: " + path2.string()};
    out << std::fixed << std::setprecision(2);
    for (const auto& p : v) out << p << '\n';
    out.close();
  }
} // WritePolygon

void WriteMP(const MP& mp, const fs::path& path) {
  if (mp.size() == 1) {
    WritePolygon(mp.front(), path);
    return;
  }
  auto ext = path.extension();
  int num = 0;
  for (const auto& p : mp) {
    auto path2 = path.stem();
    path2 += std::to_string(++num);
    path2 += ext;
    WritePolygon(p, path2);
  }
} // WriteMP

// Build polygon from (possibly open) points
Polygon MakePolygon(std::vector<Pt> pts) {
  auto ring = Ring{pts.begin(), pts.end()};
  auto poly = Polygon{ring};
  bg::correct(poly);
  Expects(poly.outer().front() == poly.outer().back());
  return poly;
}

void WriteDeflections(std::ostream& out, const Ring& ring) {
  Expects(ring.size() >= 3 && ring.front() == ring.back());
  out << std::fixed << std::setprecision(2);
  auto n = std::ssize(ring) - 1;
  auto curr = ring[0] - ring[n-1];
  for (auto i = gsl::index{0}; i != n; ++i) {
    auto prev = curr;
    curr = ring[i+1] - ring[i];
    auto deg  = ToDegrees(curr.angle_wrt(prev));
    out << ring[i] << ' ' << deg << '\n';
  }
} // WriteDeflections

void WriteDeflections(std::ostream& out, const Polygon& poly) {
  WriteDeflections(out, poly.outer());
  for (const auto& r: poly.inners())
    WriteDeflections(out, r);
} // WriteDeflections

void WriteDeflections(std::ostream& out, const MP& mp) {
  for (const auto& poly: mp)
    WriteDeflections(out, poly);
}

template<class Geo>
void WriteDeflections(const Geo& geo, const fs::path& path = "deflect.txt") {
  auto out = std::ofstream{path};
  if (!out)
    throw std::runtime_error{"cannot open deflection file: " + path.string()};
  WriteDeflections(out, geo);
  out.close();
}

void WriteCorners(std::ostream& out, const Ring& ring, const CornerVec& corners)
{
  for (auto i : corners) {
    const auto& c = ring.at(i);
    out << c << '\n';
  }
} // WriteCorners

void WriteCorners(const Ring& ring, const CornerVec& corners,
                  const fs::path& path = "corners.xy")
{
  auto out = std::ofstream{path};
  if (!out)
    throw std::runtime_error{"cannot open corners file: " + path.string()};
  out << std::fixed << std::setprecision(2);
  WriteCorners(out, ring, corners);
  out.close();
} // WriteCorners

void WriteCorners(const Polygon& poly, const std::vector<CornerVec>& allCorners,
                  const fs::path& path = "corners.xy")
{
  auto out = std::ofstream{path};
  if (!out)
    throw std::runtime_error{"cannot open corners file: " + path.string()};
  out << std::fixed << std::setprecision(2);
  auto cp = allCorners.begin();
  WriteCorners(out, poly.outer(), *cp);
  for (const auto& r: poly.inners())
    WriteCorners(out, r, *++cp);
  out.close();
} // WriteCorners

CornerVec FindCornersSimp(const Ring& ring) {
  Expects(ring.front() == ring.back());
  auto corners = CornerVec{};
  constexpr auto Theta = geom::ToRadians(Tune::CornerAngleDeg).value();
  auto n = std::ssize(ring) - 1;
  auto curr = ring[0] - ring[n-1];
  for (auto i = gsl::index{0}; i != n; ++i) {
    auto prev = curr;
    curr = ring[i+1] - ring[i];
    auto th  = curr.angle_wrt(prev);
    if (std::abs(th.value()) >= Theta) corners.push_back(i);
  }
  using namespace std;
  cout << corners.size() << " corners found at:\n";
  for (auto i: corners)
    cout << setw(3) << i << '\t' << ring[i] << '\n';
  return corners;
} // FindCornersSimp

// Map simplified-corner points to indices in the ORIGINAL ring (ordered match)
CornerVec MapCornersToOriginal(const Ring& orig, const Ring& simp,
                               const CornerVec& simp_corners)
{
  auto out = CornerVec{};
  out.reserve(simp_corners.size());
  if (orig.empty() || simp.empty() || simp_corners.empty()) return out;

  auto i0 = gsl::index{0};

  for (auto simp_idx: simp_corners) {
    const Pt& corner = simp[simp_idx];
    auto best_i  = i0;
    auto best_d2 = Dist2(orig[i0], corner);

    // scan forward around the ring once
    const auto n = std::ssize(orig) - 1;
    for (auto i = i0+1; i < n; ++i) {
      const auto d2 = Dist2(orig[i], corner);
      if (d2 < best_d2) { best_d2 = d2; best_i = i; }
    }
    out.push_back(best_i);
    i0 = std::min(gsl::index{0}, best_i + 1);
  }

  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
} // MapCornersToOriginal

CornerVec FindCorners(const Ring& ring) {
  auto simp = Ring{};
  bg::simplify(ring, simp, Tune::SimplifyForCorners);

  auto simp_corners = FindCornersSimp(simp);
  auto corners =  MapCornersToOriginal(ring, simp, simp_corners);
  using namespace std;
  cout << corners.size() << " corners found at:\n";
  for (auto i: corners)
    cout << setw(3) << i << '\t' << ring[i] << '\n';
  return corners;
} // FindCorners

void AdjustCorners(Ring& ring, CornerVec& corners) {
  Expects(ring.size() >= 3);
  Expects(ring.front() == ring.back());
  ring.pop_back();
  if (corners.empty())
    corners.push_back(0);
  if (corners.front() != 0) {
    // Rotate corners so the ring's first coordinate is a corner.
    auto shift1 = corners.front();
    auto shift2 = corners.back() - std::ssize(ring);
    auto mid = (shift1 < -shift2) ? shift1 : shift2;
    std::cout << "Rotating corners: " << mid << '\n';
    if (mid >= 0) {
      for(auto& c: corners)
        c -= mid;
      std::rotate(ring.begin(), ring.begin() + mid, ring.end());
    } else {
      mid = -mid;
      corners.pop_back();
      for (auto& c: corners)
        c += mid;
      corners.insert(corners.begin(), 0);
      std::rotate(ring.begin(), ring.begin() + (ring.size() - mid), ring.end());
    }
  }
  if (corners.size() < 2) {
    std::cout << "Adding distant corner: ";
    // If there's only one corner, add another at the most distant point.
    const auto begin = ring.begin();
    const auto end   = ring.end();
    auto p = begin;
    auto farthest_p = ++p;
    auto farthest_d = Dist(*p, *begin);
    while (++p != end) {
      auto d = Dist(*p, *begin);
      if (d <= farthest_d)
        continue;
      farthest_p = p;
      farthest_d = d;
    }
    auto idx = std::distance(begin, farthest_p);
    std::cout << idx << ' ' << *farthest_p << '\n';
    corners.push_back(idx);
  }
  ring.push_back(ring.front()); // close ring
} // AdjustCorners

MP ComputeInset(const Polygon& in, double offset) {
    // ---- Inset buffer (negative distance) ----
    auto distance = bg::strategy::buffer::distance_symmetric<double>{-offset};
    auto side     = bg::strategy::buffer::side_straight{};
    auto join_m   = bg::strategy::buffer::join_miter{Tune::MiterLimit};
    auto end      = bg::strategy::buffer::end_flat{};
    auto circle   = bg::strategy::buffer::point_circle{Tune::CirclePoints};

    auto inset = MP{};
    bg::buffer(in, inset, distance, side, join_m, end, circle);
    std::cout << "Inset generated " << inset.size() << " polygons.\n";
    return inset;
#if 0
    const Polygon& inset_poly0 = BiggestByArea(inset_mp);
    Polygon inset_poly = inset_poly0;
#endif
} // ComputeInset

} // anonymous

Ring Smooth(Ring ring) {
  auto corners = FindCorners(ring);
  AdjustCorners(ring, corners);
  auto dst = Ring{};
  auto n = corners.size() - 1;
  for (auto i = gsl::index{0}; i != n; ++i) {
    auto line = Smooth(ring, corners[i], corners[i+1]);
    dst.append_range(line);
  }
  bg::correct(dst);
  return dst;
} // Smooth ring

void Smooth(Polygon& poly) {
  poly.outer() = Smooth(poly.outer());
  for (auto& r: poly.inners())
    r = Smooth(r);
} // Smooth polygon

void Smooth(MP& mp) {
  for (auto& p: mp)
    Smooth(p);
} // Smooth multi-polygon

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

    auto reason = std::string{};
    if (!bg::is_valid(poly_in, reason)) {
      std::cerr << "Invalid polygon: " << reason << '\n';
      return EXIT_FAILURE;
    }
    std::cerr << "Polygon is valid\n";

    // ---- NEW: Write deflection angles (degrees) for simplified perimeter
    WriteDeflections(poly_in, "deflect0.txt");

#if 0
    auto allCorners = std::vector<CornerVec>{};
    {
      auto corners = FindCorners(poly_in.outer());
      AdjustCorners(poly_in.outer(), corners);
      allCorners.emplace_back(std::move(corners));
    }
    for (auto& r: poly_in.inners()) {
      auto corners = FindCorners(r);
      AdjustCorners(r, corners);
      allCorners.emplace_back(std::move(corners));
    }
    WriteCorners(poly_in, allCorners);

    auto inset = ComputeInset(poly_in, offset);
#endif

    auto mp_out = MP{};
    mp_out.emplace_back(poly_in);
    // bg::simplify(inset, mp_out, Tune::SimplifyOutputTol);

    Smooth(mp_out);
    DBG((mp_out.size()));
    WriteDeflections(mp_out, "deflect.txt");
    WriteMP(mp_out, out_path);

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
