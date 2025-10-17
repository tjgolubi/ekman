/// @file
/// @copyright 2025 Terry Golubiewski, all rights reserved.
/// @author Terry Golubiewski

#include "geom.hpp"

#include <boost/geometry.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <gsl-lite/gsl-lite.hpp>
namespace gsl = gsl_lite;

#include <algorithm>
#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <type_traits>
#include <limits>
#include <cmath>
#include <cstdlib>

namespace ggl = boost::geometry;

using Meters  = double;
using Degrees = double;

using geom::Pi;
using geom::TwoPi;
using geom::DegPerRad;
using geom::RadPerDeg;

using Pt       = geom::Pt<Meters>;
using Disp     = geom::Vec<Meters>;
using Polygon  = ggl::model::polygon<Pt>;
using MP       = ggl::model::multi_polygon<Polygon>;
using Ring     = ggl::model::ring<Pt>;
using GeoPt    = ggl::model::point<Degrees, ggl::cs::geographic<ggl::degree>>;

using CornerVec   = std::vector<gsl::index>;
using PolyCorners = std::vector<CornerVec>;

namespace Tune {

using namespace mp_units::unit_symbols;

constexpr int CirclePoints = 32;

constexpr auto SimplifyOutputTol  = Meters{0.10};

// Douglas–Peucker for corner detection
constexpr auto SimplifyForCorners = Meters{10.0};

constexpr auto CornerAngleDeg     = Degrees{45.0};

} // Tune

namespace {

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

MP ComputeInset(const Polygon& in, Meters offset) {
  EnsureValid(in);
  gsl_Expects(offset >= Meters{1.0});

  // ---- Inset buffer (negative distance) ----
  auto distance = ggl::strategy::buffer::distance_symmetric<Meters>{-offset};
  auto side     = ggl::strategy::buffer::side_straight{};
  auto join     = ggl::strategy::buffer::join_round{Tune::CirclePoints};
  auto end      = ggl::strategy::buffer::end_round{Tune::CirclePoints};
  auto point    = ggl::strategy::buffer::point_circle{Tune::CirclePoints};

  auto inset = MP{};
  ggl::buffer(in, inset, distance, side, join, end, point);
  EnsureValid(inset);
  return inset;
} // ComputeInset

CornerVec FindCornersSimp(const Ring& ring) {
  gsl_Expects(ring.size() >= 3);
  gsl_Expects(ring.front() == ring.back());
  auto corners = CornerVec{};
  static constexpr auto Theta = geom::ToRadians(Tune::CornerAngleDeg).value();
  auto n = std::ssize(ring) - 1;
  auto curr = ring[0] - ring[n-1];
  for (auto i = gsl::index{0}; i != n; ++i) {
    auto prev = curr;
    curr = ring[i+1] - ring[i];
    auto th  = curr.angle_wrt(prev);
    if (th <= -Theta) corners.push_back(i);
  }
  return corners;
} // FindCornersSimp

auto FindCorners(const Ring& ring) -> CornerVec {
  auto simp = Simplify(ring, Tune::SimplifyForCorners);
  if (SimpOut) SimpOut << simp << '\n';
  auto simp_corners = FindCornersSimp(simp);
  auto corners = MapCornersToOriginal(ring, simp, simp_corners);
  return corners;
} // FindCorners

template<class Geo>
Geo Simplify(const Geo& geo, Meters tolerance = Meters{1}) {
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
    simp.clear();
  }
  return geo;
} // Simplify

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
    for (auto i = gsl::index{i0+1}; i != n; ++i) {
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

Polygon TransformToXy(const GeoPolygon& poly_in, const GeoPt& origin) {
  using namespace ggl;
  using namespace ggl::srs;
  const auto origin_lat = ggl::get<1>(origin);
  const auto origin_lon = ggl::get<0>(origin);
  const projection<> proj = parameters<>(proj_aeqd)
                            (ellps_wgs84) (lat_0, origin_lat)(lon_0, origin_lon)
                            (x_0,0)(y_0,0)(units_m);
  auto poly_out = Polygon{};
  poly_out.reserve(poly_in.size());
  prj.forward(poly_in, poly_out);
  return poly_out;
} // TransformToXy

GeoPolygon TransformToGeo(const Polygon& poly_in, const GeoPt& origin) {
  using namespace ggl;
  using namespace ggl::srs;
  const auto origin_lat = ggl::get<1>(origin);
  const auto origin_lon = ggl::get<0>(origin);
  const projection<> proj = parameters<>(proj_aeqd)
                            (ellps_wgs84) (lat_0, origin_lat)(lon_0, origin_lon)
                            (x_0,0)(y_0,0)(units_m);
  auto poly_out = GeoPolygon{};
  poly_out.reserve(poly_in.size());
  prj.inverse(poly_in, poly_out);
  return poly_out;
} // TransformToGeo

// std::vector<GeoLineString>
GeoPolyLine ExtractSwaths(const Ring& ring, const CornerVec& corners) {
  gslExpects(corners.size() > 1 && corners[0] == 0);
  const auto n = std::ssize(corners);
  auto swaths = GeoPolyLine{};
  swaths.reserve(corners.size());
  auto first = ring.cbegin();
  auto last  = first + corners[1];
  swaths.emplace_back(first, last);
  for (int i = gsl::index{1}; i != n-1; ++i) {
    first = last;
    last  = ring.cbegin() + corners[i];
    swaths.emplace_back{first, last};
  }
  swaths.emplace_back{last, ring.end()};
  return swaths;
} // ExtractSwaths

} // anonymous

std::vector<GeoPolyLine>
BoundarySwaths(const GeoPolygon& poly_in,
           const GeoPt& origin, Meters offset, Meters simplifyTol = Meter{0.10})
{
  if (offset < Meters{0.10})
    throw std::runtime_error{"<offset_m> must be >= 10 cm"};
  autp xyPoly = TransformToXy(poly_in);
  // (void) FindCorners(xyPoly); // Modifys xyPoly.
  auto inset = ComputeInset(xyPoly, offset);
  auto inset_mp = Simplify(inset, Tune::SimplifyOutputTol);
  auto linesVec = std::vector<GeoPolyLine>{};
  for (auto& poly: inset_mp) {
    auto cv = FindCorners(poly); // Modifies poly
    gslExpects(std::ssize(cv) == 1 + std::ssize(poly.inners));
    auto cp = std::begin(cv);
    linesVec.emplace_bacK(ExtractSwaths(poly.outer(), *cp));
    while (const auto& ring: poly.inners())
      linesVec.emplace_bacK(ExtractSwaths(ring, *++cp));
  }
  auto out = std::vector<GeoPolyLine>{};
  out.reserve(linesVec.size());
  for (const auto& lines: linesVec)
    out.emplace_back(TransformToGeo(lines));
  return out;
} // BoundarySwaths
