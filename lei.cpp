#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/srs/projection.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/algorithms/distance.hpp>

#include <mp-units/systems/isq_angle.h>
#include <mp-units/systems/isq/space_and_time.h>
#include <mp-units/systems/angular/units.h>
#include <mp-units/systems/si.h>
#include <mp-units/framework/quantity.h>

#include <iostream>
#include <iomanip>

namespace ggl = boost::geometry;

namespace geo {

constexpr struct latitude  final
  : mp_units::quantity_spec<mp_units::angular::angle> {} latitude;
constexpr struct longitude final
  : mp_units::quantity_spec<mp_units::angular::angle> {} longitude;

} // geo

using GeoPt = ggl::model::point<double, 2, ggl::cs::geographic<ggl::degree>>;
using XyPt  = ggl::model::point<double, 2, ggl::cs::cartesian>;

using Latitude =mp_units::quantity<geo::latitude [mp_units::angular::degree]>;
using Longitude=mp_units::quantity<geo::longitude[mp_units::angular::degree]>;

//using Latitude =mp_units::quantity<geo::latitude,  mp_units::angular::degree>;
//using Longitude=mp_units::quantity<geo::longitude, mp_units::angular::degree>;

struct GeographicCoordinate {
  Latitude  lat;
  Longitude lon;
  GeographicCoordinate() = default;
  GeographicCoordinate(Latitude lat_, Longitude lon_)
    : lat{lat_},  lon{lon_} { }
  explicit GeographicCoordinate(const GeoPt& pt)
    : lat{ggl::get<1>(pt) * mp_units::angular::degree}
    , lon{ggl::get<0>(pt) * mp_units::angular::degree}
    { }
}; // GeographicCoordinate

namespace boost::geometry::traits {

template<>
struct tag<GeographicCoordinate> { using type = point_tag; };

template<>
struct coordinate_type<GeographicCoordinate> { using type = double; };

template<>
struct coordinate_system<GeographicCoordinate>
  { using type = cs::geographic<degree>; };

template<>
struct dimension<GeographicCoordinate>
  : std::integral_constant<std::size_t, 2> { };

template<std::size_t Dim>
requires (Dim == 0 || Dim == 1)
struct access<GeographicCoordinate, Dim> {
  static constexpr auto deg = mp_units::angular::degree;
  static constexpr double get(const GeographicCoordinate& p) {
    if constexpr (Dim == 0)
      return p.lon.numerical_value_in(deg);
    else
      return p.lat.numerical_value_in(deg);
  }
  static constexpr void set(GeographicCoordinate& p, double v) {
    if constexpr (Dim == 0)
      p.lon = v * deg;
    else
      p.lat = v * deg;
  }
}; // access

template<>
struct make<GeographicCoordinate> {
  static constexpr auto deg = mp_units::angular::degree;
  using point_type = GeographicCoordinate;
  static constexpr auto is_specialized = true;
  static constexpr point_type apply(double x, double y)
    { return point_type{y * deg, x * deg}; }
}; // make

} // boost::geometry::traits

std::ostream& operator<<(std::ostream& os, const GeoPt& pt) {
  return os << '[' << ggl::get<0>(pt) << ' ' << ggl::get<1>(pt) << ']';
}

std::ostream& operator<<(std::ostream& os, const XyPt& pt) {
  return os << '(' << ggl::get<0>(pt) << ' ' << ggl::get<1>(pt) << ')';
}

std::ostream& operator<<(std::ostream& os, const GeographicCoordinate& gc) {
  return os << '<' << gc.lat << ", " << gc.lon << '>';
}

constexpr GeoPt Geo(const GeographicCoordinate& gc) noexcept {
  using mp_units::angular::unit_symbols::deg;
  return GeoPt{gc.lon.numerical_value_in(deg), gc.lat.numerical_value_in(deg)};
}

int main(int, char*[]) {
    using namespace boost::geometry;
    using namespace boost::geometry::model;
    using namespace boost::geometry::srs;
    using namespace boost::geometry::srs::dpar;

    using LatLon = GeographicCoordinate;

    using mp_units::angular::unit_symbols::deg;
    auto home  = LatLon{41.907199 * deg, -91.661078 * deg};
    auto hotel = LatLon{41.903958 * deg, -91.655476 * deg};
    std::cout << "home=" << home << ' ' << "hotel=" << hotel << '\n';

    //const auto pt_ll0 = Geo(home);
    //const auto pt_ll  = Geo(hotel);

    const auto pt_xy0 = XyPt{0.0, 0.0};
    std::cout << "pt_xy0=" << pt_xy0 << '\n';

    std::cout << " dist=" << distance(home, hotel) << std::endl;

    using mp_units::angular::unit_symbols::deg;
    auto home_lat = home.lat.numerical_value_in(deg);
    auto home_lon = home.lon.numerical_value_in(deg);
    projection<> prj = parameters<>(proj_aeqd)
              (ellps_wgs84) (lat_0, home_lat)(lon_0, home_lon)
              (x_0,0.0)(y_0,0.0)(units_m);

    auto pt_xy  = XyPt{};
    prj.forward(hotel, pt_xy);

    std::cout << pt_xy << " dist=" << distance(pt_xy, pt_xy0) << std::endl;

    //auto pt_ll2 = GeoPt{0.0, 0.0};
    auto hotel2 = GeoPt{};
    prj.inverse(pt_xy, hotel2);

    std::cout << hotel2 << " dist=" << distance(hotel, hotel2) << std::endl;

    return 0;
} // main
