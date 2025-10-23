/// @file
/// @copyright 2025 Terry Golubiewski, all rights reserved.
/// @author Terry Golubiewski
///
/// Read XY perimeter (meters), build an inward inset path, detect "real" corners
/// by simplifying the ORIGINAL perimeter (0.5 m DP), map those corners onto the
/// INSET ring, smooth BETWEEN those corners (endpoints anchored) with a heading
/// moving-average, and write:
///   - <output.xy> : smoothed inset ring (CLOSED)
///   - corners.xy  : corner points on the inset ring (for plotting)
///   - deflect.txt : deflection angle (degrees) for each vertex of the
///                   simplified-for-corners perimeter (one value per line)
///
/// Usage:
///   inset_bg <input.xy> <offset_m> <output.xy>
///
/// Input:  one "x y" per line; not necessarily closed.
/// Output: one "x y" per line; CLOSED ring (last == first).

#include "Smooth.hpp"
#include "Resample.hpp"

#include "geom_ggl.hpp"

#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/simplify.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/strategies/buffer/cartesian.hpp>
#include <boost/geometry/strategies/agnostic/buffer_distance_symmetric.hpp>
#include <boost/geometry/strategies/cartesian/buffer_side_straight.hpp>
#include <boost/geometry/strategies/cartesian/buffer_join_round.hpp>
#include <boost/geometry/strategies/cartesian/buffer_end_round.hpp>
#include <boost/geometry/strategies/cartesian/buffer_point_circle.hpp>

#include <boost/preprocessor/stringize.hpp>
#include <gsl-lite/gsl-lite.hpp>
namespace gsl = gsl_lite;

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
#include <type_traits>
#include <limits>
#include <cmath>
#include <cstdlib>

template<typename T>
void DbgOut(const char* str, const T& x)
  { std::cerr << str << '=' << x << std::endl; }

#define DBG(x) DbgOut(BOOST_STRINGIZE(x), (x))

namespace ggl = boost::geometry;
namespace fs = std::filesystem;

using Meters = double;
using Degrees = double;

using Pt       = geom::Pt<Meters>;
using Disp     = geom::Vec<Meters>;
using Polygon  = ggl::model::polygon<Pt>;
using MultiPolygon       = ggl::model::multi_polygon<Polygon>;
using Ring     = ggl::model::ring<Pt>;

using CornerVec   = std::vector<gsl::index>;
using PolyCorners = std::vector<CornerVec>;

auto SimpOut = std::ofstream{};

namespace Tune {

#if 0
constexpr Meters MiterLimit   = 4.0;
#endif

constexpr int CirclePoints = 32;

#if 0
// Geometry cleanup (applied to geometry we will drive)
constexpr Meters SimplifyInputTol    = 0.10; // m; <=0 → skip
constexpr Meters SimplifyInsetTol    = 0.10; // m; <=0 → skip
constexpr Meters SimplifySmoothedTol = 0.10; // try 0.10–0.30
#endif

constexpr Meters SimplifyOutputTol   = 0.10; // try 0.10–0.30

// Corner detection (on ORIGINAL perimeter via DP 0.5 m)
constexpr Meters  SimplifyForCorners = 10.0; // m; Douglas–Peucker for corner detection
constexpr Degrees CornerAngleDeg     = 45.0; // deflection threshold

#if 0
constexpr Meters  MinSegLen          = 0.50; // m; ignore tiny edges near vertex

// Smoothing between mapped corner indices on the inset ring
constexpr Meters ResampleStep       = 1.00; // m

constexpr int    HeadingHalfWin     = 3;    // moving-average half-window (2*H+1)

// Safety corridor (optional clamp). 0.0 disables.
constexpr Meters CorridorShrink     = 0.0;  // m
#endif

} // Tune

namespace {

std::ostream& operator<<(std::ostream& os, const Pt& p)
  { return os << p.x() << ' ' << p.y(); }

std::istream& operator>>(std::istream& os, Pt& p)
  { return os >> p.x() >> p.y(); }

std::ostream& operator<<(std::ostream& os, const Ring& r) {
  for (const auto& p: r)
    os << p << '\n';
  return os;
}

#if 0
void EnsureValid(const Ring& ring, bool isInner) {
  auto failure = ggl::validity_failure_type{};
  if (ggl::is_valid(ring, failure)) [[likely]] {
    if (!isInner)
      return;
    failure = ggl::failure_wrong_orientation;
  }
  else {
    if (isInner && failure == ggl::failure_wrong_orientation)
      return;
  }
  auto msg = std::string{"Invalid ring: "};
  msg += ggl::validity_failure_type_message(failure);
  throw std::runtime_error{msg};
} // EnsureValid
#endif

template<class Geo>
requires (!std::is_same_v<Geo, Ring>)
void EnsureValid(const Geo& geo) {
  auto failure = ggl::validity_failure_type{};
  if (ggl::is_valid(geo, failure)) [[likely]]
    return;
  auto msg = std::string{"Invalid geometry: "};
  msg += ggl::validity_failure_type_message(failure);
  throw std::runtime_error{msg};
} // EnsureValid

#if 0
Ring EnsureClosed(const Ring& ring) {
  if (ring.empty() || ring.front() == ring.back()) return ring;
  auto out = ring;
  out.push_back(ring.front());
  return out;
}
std::vector<Pt> ReadPoints(const fs::path& path) {
  auto in = std::ifstream{path, std::ios::binary};
  if (!in) throw std::runtime_error{"cannot open input: " + path.string()};
  auto pts = std::vector<Pt>{};
  auto pt = Pt{};
  while (in >> pt) pts.push_back(pt);
  if (pts.size() < 3)
    throw std::runtime_error{"not enough points in: " + path.string()};
  in.close();
  return pts;
} // ReadPoints

void WritePoints(const std::vector<Pt>& pts, const fs::path& path) {
  auto out = std::ofstream{path, std::ios::binary};
  if (!out) throw std::runtime_error{"cannot open output: " + path.string()};
  out << std::fixed << std::setprecision(2);
  for (const auto& p : pts) out << p << '\n';
  out.close();
} // WritePoints
#endif

Polygon ReadPolygon(const fs::path& path) {
  auto in = std::ifstream{path, std::ios::binary};
  if (!in) throw std::runtime_error{"cannot open input: " + path.string()};
  auto poly = Polygon{};
  auto* ring = &poly.outer();
  auto pts = std::vector<Pt>{};
  auto line = std::string{};
  auto pt = Pt{};
  while (getline(in, line)) {
    if (line.empty())
      continue;
    if (line != "# inner") {
      auto iss = std::istringstream{line};
      iss >> pt;
      ring->push_back(pt);
      continue;
    }
    // Start a new inner ring.
    poly.inners().emplace_back();
    ring = &poly.inners().back();
  }
  in.close();
  ggl::correct(poly);
  EnsureValid(poly);
  return poly;
} // ReadPolygon

void WritePolygon(const Polygon& poly, const fs::path& path) {
  auto out = std::ofstream{path, std::ios::binary};
  if (!out) throw std::runtime_error{"cannot open output: " + path.string()};
  out << std::fixed << std::setprecision(2);
  for (const auto& p : poly.outer()) out << p << '\n';
  out.close();
  if (poly.inners().empty())
    return;
  auto ext   = path.extension();
  auto inner = 0;
  for (const auto& v : poly.inners()) {
    auto path2 = path;
    auto ext2  = fs::path{std::to_string(++inner)};
    ext2 += ext;
    path2.replace_extension(ext2);
    out = std::ofstream{path2, std::ios::binary};
    if (!out)
      throw std::runtime_error{"cannot open output: " + path2.string()};
    out << std::fixed << std::setprecision(2);
    for (const auto& p : v) out << p << '\n';
    out.close();
  }
} // WritePolygon

void WriteMP(const MultiPolygon& mp, const fs::path& path) {
  if (mp.size() == 1) {
    WritePolygon(mp.front(), path);
    return;
  }
  auto ext = path.extension();
  auto num = 0;
  for (const auto& p : mp) {
    auto path2 = path.stem();
    path2 += std::to_string(++num);
    path2 += ext;
    WritePolygon(p, path2);
  }
} // WriteMP

#if 0
// Build polygon from (possibly open) points
Polygon MakePolygon(std::vector<Pt> pts) {
  auto ring = Ring{pts.begin(), pts.end()};
  auto poly = Polygon{ring};
  ggl::correct(poly);
  EnsureValid(poly);
  return poly;
}
#endif

void WriteDeflections(std::ostream& out, const Ring& ring) {
  const auto deflections = MakeDeflections(ring, Meters{1});
  out << std::fixed << std::setprecision(2);
  for (auto d: deflections)
    out << d.degrees() << '\n';
} // WriteDeflections

void WriteDeflections(std::ostream& out, const Polygon& poly) {
  WriteDeflections(out, poly.outer());
  for (const auto& r: poly.inners())
    WriteDeflections(out, r);
} // WriteDeflections

void WriteDeflections(std::ostream& out, const MultiPolygon& mp) {
  for (const auto& poly: mp)
    WriteDeflections(out, poly);
}

template<class Geo>
void WriteDeflections(const Geo& geo, const fs::path& path = "deflect.txt") {
  EnsureValid(geo);
  auto out = std::ofstream{path, std::ios::binary};
  if (!out)
    throw std::runtime_error{"cannot open deflection file: " + path.string()};
  out.exceptions(std::ios::failbit|std::ios::badbit);
  WriteDeflections(out, geo);
  out.close();
} // WriteDeflections (Geo)

void WriteCorners(std::ostream& out, const Ring& ring, const CornerVec& corners)
{
  for (auto i : corners) {
    const auto& c = ring.at(i);
    out << c << '\n';
  }
} // WriteCorners

#if 0
void WriteCorners(const Ring& ring, const CornerVec& corners,
                  const fs::path& path = "corners.xy")
{
  auto out = std::ofstream{path, std::ios::binary};
  if (!out)
    throw std::runtime_error{"cannot open corners file: " + path.string()};
  out << std::fixed << std::setprecision(2);
  WriteCorners(out, ring, corners);
  out.close();
} // WriteCorners
#endif

void WriteCorners(std::ostream& out,
                  const Polygon& poly, const PolyCorners& corners)
{
  gsl_Expects(corners.size() == 1 + poly.inners().size());
  auto iter = corners.begin();
  WriteCorners(out, poly.outer(), *iter);
  for (const auto& r: poly.inners())
    WriteCorners(out, r, *++iter);
} // WriteCorners

#if 0
void WriteCorners(const fs::path& path,
                  const Polygon& poly, const PolyCorners& allCorners)
{
  auto out = std::ofstream{path, std::ios::binary};
  if (!out)
    throw std::runtime_error{"cannot open corners file: " + path.string()};
  out << std::fixed << std::setprecision(2);
  WriteCorners(out, poly, allCorners);
  out.close();
} // WriteCorners
#endif

void WriteCorners(std::ostream& out,
                  const MultiPolygon& polys, const std::vector<PolyCorners>& mpCorners)
{
  gsl_Expects(mpCorners.size() == polys.size());
  auto mp_iter = mpCorners.begin();
  for (const auto& p: polys)
    WriteCorners(out, p, *mp_iter++);
} // WriteCorners

void WriteCorners(const fs::path& path,
                  const MultiPolygon& polys, const std::vector<PolyCorners>& allCorners)
{
  auto out = std::ofstream{path, std::ios::binary};
  if (!out)
    throw std::runtime_error{"cannot open corners file: " + path.string()};
  out << std::fixed << std::setprecision(2);
  WriteCorners(out, polys, allCorners);
  out.close();
} // WriteCorners

CornerVec FindCornersSimp(const Ring& ring) {
  gsl_Expects(ring.size() >= 3);
  gsl_Expects(ring.front() == ring.back());
  auto corners = CornerVec{};
  constexpr auto Theta = -geom::Radians::FromDegrees(Tune::CornerAngleDeg);
  auto n = std::ssize(ring) - 1;
  auto curr = ring[0] - ring[n-1];
  for (auto i = gsl::index{0}; i != n; ++i) {
    auto prev = curr;
    curr = ring[i+1] - ring[i];
    auto th  = curr.angle_wrt(prev);
    if (th <= Theta) corners.push_back(i);
  }
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
    for (auto i = i0+1; i != n; ++i) {
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

template<class Geo>
auto Simplify(const Geo& geo, Meters tolerance = Meters{1}) -> Geo {
  gsl_Expects(tolerance >= Meters{0.01});
  auto simp = Geo{};
  auto failure = ggl::validity_failure_type{};
  while (tolerance >= Meters{0.01}) {
    ggl::simplify(geo, simp, tolerance);
    if (ggl::is_valid(simp, failure) || failure == ggl::failure_wrong_orientation)
      return simp;
    switch (failure) {
      case ggl::failure_self_intersections:
      case ggl::failure_few_points:
        break;
      default: {
        auto msg = std::string{"Simplify: invalid result: "};
        msg += ggl::validity_failure_type_message(failure);
        throw std::runtime_error{msg};
      }
    }
    tolerance /= 2;
    std::cerr
          << "Simplify failed, retrying with tolerance = " << tolerance << '\n';
    simp.clear();
  }
  std::cerr << "Simplify failed: using original geo.\n";
  return geo;
} // Simplify

auto FindCorners(const Ring& ring) -> CornerVec {
  auto simp = Simplify(ring, Tune::SimplifyForCorners);
  if (SimpOut) SimpOut << simp << '\n';
  auto simp_corners = FindCornersSimp(simp);
  auto corners = MapCornersToOriginal(ring, simp, simp_corners);
  return corners;
} // FindCorners

void AdjustCorners(Ring& ring, CornerVec& corners) {
  ring.pop_back();
  if (corners.empty())
    corners.push_back(0);
  if (corners.front() != 0) {
    // Rotate corners so the ring's first coordinate is a corner.
    auto shift1 = corners.front();
    auto shift2 = corners.back() - std::ssize(ring);
    auto mid = (shift1 < -shift2) ? shift1 : shift2;
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
    corners.push_back(idx);
  }
  ring.push_back(ring.front()); // close ring
} // AdjustCorners

MultiPolygon ComputeInset(const Polygon& in, Meters offset) {
  EnsureValid(in);
  gsl_Expects(offset >= Meters{1});

  // ---- Inset buffer (negative distance) ----
  auto distance = ggl::strategy::buffer::distance_symmetric<Meters>{-offset};
  auto side     = ggl::strategy::buffer::side_straight{};
  //  auto join     = ggl::strategy::buffer::join_miter{}; // {Tune::MiterLimit};
  auto join     = ggl::strategy::buffer::join_round{Tune::CirclePoints};
  //  auto end      = ggl::strategy::buffer::end_flat{};
  auto end      = ggl::strategy::buffer::end_round{Tune::CirclePoints};
  auto point    = ggl::strategy::buffer::point_circle{Tune::CirclePoints};
  //  auto point    = ggl::strategy::buffer::point_square{};

  auto inset = MultiPolygon{};
  ggl::buffer(in, inset, distance, side, join, end, point);
  EnsureValid(inset);
  std::cout << "Inset generated " << inset.size() << " polygons.\n";
  return inset;
#if 0
  const Polygon& inset_poly0 = BiggestByArea(inset_mp);
  Polygon inset_poly = inset_poly0;
#endif
} // ComputeInset

} // anonymous

#if 0
Ring Smooth(Ring ring) {
  auto corners = FindCorners(ring);
  AdjustCorners(ring, corners);
  corners.push_back(ring.size()-1);
  auto dst = Ring{};
  const auto n = corners.size() - 1;
  for (auto i = gsl::index{0}; i != n; ++i) {
    auto line = Smooth(ring, corners[i], corners[i+1]);
    dst.append_range(line);
  }
  ggl::correct(dst);
  return dst;
} // Smooth ring

void Smooth(Polygon& poly) {
  poly.outer() = Smooth(poly.outer());
  for (auto& r: poly.inners())
    r = Smooth(r);
} // Smooth polygon

void Smooth(MultiPolygon& mp) {
  for (auto& p: mp)
    Smooth(p);
} // Smooth multi-polygon
#endif

std::vector<CornerVec> FindCorners(Polygon& poly) {
  auto allCorners = std::vector<CornerVec>{};
  {
    auto corners = FindCorners(poly.outer());
    AdjustCorners(poly.outer(), corners);
    allCorners.emplace_back(std::move(corners));
  }
  for (auto& r: poly.inners()) {
    auto corners = FindCorners(r);
    AdjustCorners(r, corners);
    allCorners.emplace_back(std::move(corners));
  }
  return allCorners;
} // FindCorners

int main(int argc, const char* argv[]) {
  try {
    auto arg0 = fs::path{argv[0]};

    if (argc != 4) {
      std::cerr << "usage:\n" << arg0.string()
        << "  <input.xy> <offset_m> <output.xy>\n";
      return EXIT_FAILURE;
    }

    SimpOut = std::ofstream{"simp.xy", std::ios::binary};

    const auto in_path  = fs::path{argv[1]};
    const auto offset   = Meters{std::stod(argv[2])}; // inward distance (m)
    const auto out_path = fs::path{argv[3]};

    if (offset <= Meters{0}) throw std::runtime_error{"<offset_m> must be > 0"};

    // ---- Read ORIGINAL perimeter and build polygon
    auto poly_in = ReadPolygon(in_path);

    //auto allCorners = FindCorners(poly_in); // Modifys poly_in.
    //WriteCorners("corners0.xy", poly_in, allCorners);

    auto inset = ComputeInset(poly_in, offset);
    auto mp_out = Simplify(inset, Tune::SimplifyOutputTol);

    {
      auto allPolyCorners = std::vector<PolyCorners>{};
      for (auto& poly: mp_out) {
        auto corners = FindCorners(poly); // Modifies poly
        allPolyCorners.emplace_back(std::move(corners));
      }
      WriteCorners("corners.xy", mp_out, allPolyCorners);
    }

    // Smooth(mp_out);

    WriteDeflections(mp_out, "deflect.txt");
    WriteMP(mp_out, out_path);

#if 0
    // Optional light cleanup for the path we’ll buffer/drive
    if (kSimplifyInputTol > 0.0) {
      auto simp = Polygon{};
      ggl::simplify(poly_in, simp, kSimplifyInputTol);
      poly_in = std::move(simp);
    }

    // ---- Corner detection on ORIGINAL perimeter via DP(0.5 m)
    auto poly_dp = Polygon{};
    ggl::simplify(poly_in, poly_dp, kSimplifyForCorners);

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
    auto distance = ggl::strategy::buffer::distance_symmetric<double>{-offset};
    auto side     = ggl::strategy::buffer::side_straight;
    auto join_m   = ggl::strategy::buffer::join_miter{kMiterLimit};
    auto end      = ggl::strategy::buffer::end_flat;
    auto circle   = ggl::strategy::buffer::point_circle{kCirclePoints};

    auto inset_mp = MultiPolygon{};
    ggl::buffer(poly_in, inset_mp, distance, side, join_m, end, circle);
    const Polygon& inset_poly0 = BiggestByArea(inset_mp);
    Polygon inset_poly = inset_poly0;

    if (kSimplifyInsetTol > 0.0) {
      auto simp = Polygon{};
      ggl::simplify(inset_poly, simp, kSimplifyInsetTol);
      ggl::correct(simp);
      inset_poly = std::move(simp);
      ggl::unique(inset_poly);
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
      MultiPolygon corr_mp;
      auto dist2 =
            ggl::strategy::buffer::distance_symmetric<double>{-kCorridorShrink};
      auto j2 = ggl::strategy::buffer::join_miter{kMiterLimit};
      ggl::buffer(inset_poly, corr_mp, dist2, side, j2, end, circle);
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
