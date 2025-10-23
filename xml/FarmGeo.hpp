#pragma once
#include "FarmDb.hpp"

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/core/cs.hpp>

namespace boost::geometry::traits {

template<>
struct tag<farm_db::LatLon> { using type = point_tag; };

template<>
struct coordinate_type<farm_db::LatLon> { using type = double; };

template<>
struct coordinate_system<farm_db::LatLon>
  { using type = cs::geographic<degree>; };

template<>
struct dimension<farm_db::LatLon>
  : std::integral_constant<std::size_t, 2> { };

template<std::size_t Dim>
requires (Dim == 0 || Dim == 1)
struct access<farm_db::LatLon, Dim> {
  static constexpr auto deg = mp_units::angular::degree;
  static constexpr double get(const farm_db::LatLon& p) {
    if constexpr (Dim == 0)
      return p.longitude.numerical_value_in(deg);
    else
      return p.latitude.numerical_value_in(deg);
  }
  static constexpr void set(farm_db::LatLon& p, double v) {
    if constexpr (Dim == 0)
      p.longitude = v * deg;
    else
      p.latitude  = v * deg;
  }
}; // access

template<>
struct make<farm_db::LatLon> {
  static constexpr auto deg = mp_units::angular::degree;
  using point_type = farm_db::LatLon;
  static constexpr auto is_specialized = true;
  static constexpr point_type apply(double x, double y)
    { return point_type{y * deg, x * deg}; }
}; // make

} // boost::geometry::traits

namespace farm_db {

namespace ggl = boost::geometry;

namespace geo {
using Point        = LatLon;
using LineString   = ggl::model::linestring<Point>;
using PolyLine     = ggl::model::multi_linestring<LineString>;
using Ring         = ggl::model::ring<Point, true >;
using Hole         = ggl::model::ring<Point, false>;
using Polygon      = ggl::model::polygon<Point>;
using MultiPolygon = ggl::model::multi_polygon<Polygon>;
using Path         = LineString;
using MultiPath    = PolyLine;
} // geo

constexpr geo::Point Geo(const Point& pt) noexcept { return pt.point; }
constexpr Point MakePoint(const geo::Point& pt, Point::Type type) noexcept
  { return Point{pt, type}; }

geo::LineString Geo(const LineString& lstr);
geo::Polygon Geo(const Polygon& poly);

geo::Ring MakeGeoRing(const LineString& lstr);
geo::Hole MakeGeoHole(const LineString& lstr);

Path MakePath(const geo::LineString& lstr);
Path MakePath(const geo::Ring& ring);

} // farm_db
