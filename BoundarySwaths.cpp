/// @file
/// @copyright 2025 Terry Golubiewski, all rights reserved.
/// @author Terry Golubiewski
#include "BoundarySwaths.hpp"

#include <mp-units/systems/si/unit_symbols.h>

#include <boost/geometry/srs/projection.hpp>

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

namespace tjg {

using Disp     = geom::Vec<Meters>;
using MP       = ggl::model::multi_polygon<Polygon>;
using Ring     = ggl::model::ring<Pt>;

using CornerVec   = std::vector<gsl::index>;
using PolyCorners = std::vector<CornerVec>;

namespace Tune {

constexpr int CirclePoints = 32;

// Douglas–Peucker for corner detection
constexpr auto SimplifyForCorners = 10.0 * metre;
constexpr Degrees CornerAngleDeg = 45.0;

} // Tune

namespace detail {

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

MP ComputeInset(const Polygon& in, Distance offset) {
  EnsureValid(in);
  gsl_Expects(offset >= 1.0 * metre);

  // ---- Inset buffer (negative distance) ----
  auto distance = ggl::strategy::buffer::distance_symmetric<Meters>{-offset.numerical_value_in(metre)};
  auto side     = ggl::strategy::buffer::side_straight{};
  auto join     = ggl::strategy::buffer::join_round{Tune::CirclePoints};
  auto end      = ggl::strategy::buffer::end_round{Tune::CirclePoints};
  auto point    = ggl::strategy::buffer::point_circle{Tune::CirclePoints};

  auto inset = MP{};
  ggl::buffer(in, inset, distance, side, join, end, point);
  EnsureValid(inset);
  return inset;
} // ComputeInset

template<class Geo>
Geo Simplify(const Geo& geo, Distance tolerance) {
  gsl_Expects(tolerance >= 0.01 * metre);
  auto simp = Geo{};
  auto failure = ggl::validity_failure_type{};
  while (tolerance >= 0.01 * metre) {
    ggl::simplify(geo, simp, tolerance.numerical_value_in(metre));
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

CornerVec FindCornersSimp(const Ring& ring) {
  gsl_Expects(ring.size() >= 3);
  gsl_Expects(ring.front() == ring.back());
  auto corners = CornerVec{};
  static constexpr auto Theta = geom::Radians::FromDegrees(Tune::CornerAngleDeg);
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
  auto simp_corners = FindCornersSimp(simp);
  auto corners = MapCornersToOriginal(ring, simp, simp_corners);
  return corners;
} // FindCorners

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

PathVec ExtractSwaths(const Ring& ring, const CornerVec& corners) {
  gsl_Expects(corners.size() > 1 && corners[0] == 0);
  const auto n = std::ssize(corners);
  auto swaths = PathVec{};
  swaths.reserve(corners.size());
  auto first = std::cbegin(ring);
  auto last  = first + corners[1];
  gsl_Expects(last != ring.end());
  auto ptr = last;
  if (ptr != ring.end())
    ++ptr;
  swaths.emplace_back(first, ptr);
  for (int i = gsl::index{1}; i != n-1; ++i) {
    first = last;
    auto idx = corners.at(i+1);
    gsl_Expects(idx < ring.size());
    last  = ring.cbegin() + idx;
    gsl_Expects(last != ring.end());
    ptr = last;
    if (ptr != ring.end())
      ++ptr;
    swaths.emplace_back(first, ptr);
  }
  gsl_Expects(last != ring.end());
  swaths.emplace_back(last, ring.end());
  return swaths;
} // ExtractSwaths

template<class XY> struct Geo { using type = void; };
template<class G>  struct Xy  { using type = void; };
template<class XY> using GeoT = Geo<std::remove_cv_t<XY>>::type;
template<class G>  using XyT  =  Xy<std::remove_cv_t<G >>::type;

template<> struct Geo<Path>    { using type = GeoPath; };
template<> struct Geo<PathVec> { using type = GeoPathVec; };
template<> struct Geo<Polygon> { using type = GeoPolygon; };

template<> struct Xy<GeoPath>    { using type = Path; };
template<> struct Xy<GeoPathVec> { using type = PathVec; };
template<> struct Xy<GeoPolygon> { using type = Polygon; };

template<class G>
XyT<G> TransformToXy(const G& geo_in, const GeoPt& origin) {
  using namespace ggl;
  using namespace ggl::srs;
  using namespace ggl::srs::dpar;

  const auto origin_lat = ggl::get<1>(origin);
  const auto origin_lon = ggl::get<0>(origin);
  const projection<> proj = parameters<>(proj_aeqd)
                            (ellps_wgs84) (lat_0, origin_lat)(lon_0, origin_lon)
                            (x_0,0)(y_0,0)(units_m);
  auto geo_out = XyT<G>{};
  proj.forward(geo_in, geo_out);
  return geo_out;
} // TransformToXy

template<class X>
GeoT<X> TransformToGeo(const X& geo_in, const GeoPt& origin) {
  using namespace ggl;
  using namespace ggl::srs;
  using namespace ggl::srs::dpar;

  const auto origin_lat = ggl::get<1>(origin);
  const auto origin_lon = ggl::get<0>(origin);
  const projection<> proj = parameters<>(proj_aeqd)
                            (ellps_wgs84) (lat_0, origin_lat)(lon_0, origin_lon)
                            (x_0,0)(y_0,0)(units_m);
  auto geo_out = GeoT<X>{};
  proj.inverse(geo_in, geo_out);
  return geo_out;
} // TransformToGeo

std::vector<GeoPathVec>
TransformToGeo(const std::vector<PathVec>& in, const GeoPt& origin) {
  auto out = std::vector<GeoPathVec>{};
  out.reserve(in.size());
  for (const auto& pv: in)
    out.emplace_back(TransformToGeo(pv, origin));
  return out;
} // TransformToGeo(vector<PathVec>)

} // detail

std::vector<PathVec>
BoundarySwaths(const Polygon& poly_in, Distance offset, Distance simplifyTol) {
  if (offset < 0.10 * metre)
    throw std::runtime_error{"<offset_m> must be >= 10 cm"};
  // (void) FindCorners(poly_in); // Modifys poly_in.
  auto inset_mp = detail::ComputeInset(poly_in, offset);
  auto simp_mp  = detail::Simplify(inset_mp, simplifyTol);
  auto linesVec = std::vector<PathVec>{};
  for (auto& poly: simp_mp) {
    auto cv = detail::FindCorners(poly); // Modifies poly
    gsl_Expects(std::ssize(cv) == 1 + std::ssize(poly.inners()));
    auto cp = std::begin(cv);
    linesVec.emplace_back(detail::ExtractSwaths(poly.outer(), *cp));
    for (const auto& ring: poly.inners())
      linesVec.emplace_back(detail::ExtractSwaths(ring, *++cp));
  }
  return linesVec;
} // BoundarySwaths

std::vector<GeoPathVec>
BoundarySwaths(const GeoPolygon& poly_in, Distance offset, Distance simplifyTol)
{
  auto box = ggl::model::box<GeoPt>{};
  ggl::envelope(poly_in.outer(), box);
  auto lon = ggl::get<0>(box.min_corner());
  auto lat = ggl::get<1>(box.min_corner());
  auto dLon = (ggl::get<0>(box.max_corner()) - lon) / 2;
  auto dLat = (ggl::get<1>(box.max_corner()) - lat) / 2;
  lon += dLon;
  lat += dLat;
  auto origin = GeoPt{lon, lat};
  auto xyPoly = detail::TransformToXy(poly_in, origin);
  auto xyOut  = BoundarySwaths(xyPoly, offset, simplifyTol);
  return detail::TransformToGeo(xyOut, origin);
} // BoundarySwaths

} // tjg

// ============================================================================

