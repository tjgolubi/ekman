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

geo::LineString Geo(const LineString& lstr) {
  auto out = geo::LineString{};
  out.reserve(lstr.points.size());
  for (const auto& p: lstr.points)
    out.push_back(Geo(p));
  return out;
} // Geo(LineString)

geo::Ring MakeGeoRing(const LineString& lstr) {
  auto ls = Geo(lstr);
  auto out = geo::Ring{ls.begin(), ls.end()};
  ggl::correct(out);
  auto msg = std::string{};
  if (!ggl::is_valid(out, msg))
    throw std::runtime_error{"LineString::ring: not a ring: " + msg};
  return out;
} // MakeGeoRing

geo::Hole MakeGeoHole(const LineString& lstr) {
  auto ls = Geo(lstr);
  auto out = geo::Hole{ls.begin(), ls.end()};
  ggl::correct(out);
  auto msg = std::string{};
  if (!ggl::is_valid(out, msg))
    throw std::runtime_error{"LineString::hole: not a hole: " + msg};
  return out;
} // MakeGeoHole

geo::Polygon Geo(const Polygon& poly) {
  auto out = geo::Polygon{};
  out.outer() = MakeGeoRing(poly.outer);
  for (const auto& p: poly.inners)
    out.inners().push_back(MakeGeoRing(p));
  ggl::correct(out);
  auto msg = std::string{};
  if (!ggl::is_valid(out, msg))
    throw std::runtime_error{"Geo(Polygon): invalid polygon: " + msg};
  return out;
} // Geo(Polygon)

Path MakePath(const geo::LineString& lstr) {
  auto rval = Path{};
  rval.reserve(lstr.size());
  for (const auto& p: lstr)
    rval.push_back(p);
  return rval;
}

Path MakePath(const geo::Ring& ring) {
  auto rval = Path{};
  rval.reserve(ring.size());
  for (const auto& p: ring)
    rval.push_back(p);
  return rval;
}

} // farm_db
