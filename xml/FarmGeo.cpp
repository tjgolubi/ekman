#include "FarmGeo.hpp"

#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/core/cs.hpp>

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <exception>
#include <utility>

namespace farm_db {

GeoLineString Geo(const LineString& lstr) {
  auto out = GeoLineString{};
  out.reserve(lstr.points.size());
  for (const auto& p: lstr.points)
    out.push_back(Geo(p));
  return out;
} // Geo(LineString)

GeoRing MakeGeoRing(const LineString& lstr) {
  auto ls = Geo(lstr);
  auto out = GeoRing{ls.begin(), ls.end()};
  ggl::correct(out);
  auto msg = std::string{};
  if (!ggl::is_valid(out, msg))
    throw std::runtime_error{"LineString::ring: not a ring: " + msg};
  return out;
} // MakeGeoRing

GeoHole MakeGeoHole(const LineString& lstr) {
  auto ls = Geo(lstr);
  auto out = GeoHole{ls.begin(), ls.end()};
  ggl::correct(out);
  auto msg = std::string{};
  if (!ggl::is_valid(out, msg))
    throw std::runtime_error{"LineString::hole: not a hole: " + msg};
  return out;
} // MakeGeoHole

GeoPolygon Geo(const Polygon& poly) {
  auto out = GeoPolygon{};
  out.outer() = MakeGeoRing(poly.outer);
  for (const auto& p: poly.inners)
    out.inners().push_back(MakeGeoRing(p));
  ggl::correct(out);
  auto msg = std::string{};
  if (!ggl::is_valid(out, msg))
    throw std::runtime_error{"Geo(Polygon): invalid polygon: " + msg};
  return out;
} // Geo(Polygon)

Path MakePath(const GeoLineString& lstr) {
  auto rval = Path{};
  rval.reserve(lstr.size());
  for (const auto& p: lstr)
    rval.push_back(p);
  return rval;
}

Path MakePath(const GeoRing& ring) {
  auto rval = Path{};
  rval.reserve(ring.size());
  for (const auto& p: ring)
    rval.push_back(p);
  return rval;
}

} // farm_db
