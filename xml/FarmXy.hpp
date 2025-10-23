#pragma once
#include "FarmDb.hpp"
#include "../geom_ggl.hpp"

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/core/cs.hpp>

namespace farm_db {

namespace ggl = boost::geometry;

namespace xy {

using Meters       = double;
using Point        = geom::Pt<Meters>;
using LineString   = ggl::model::linestring<Point>;
using PolyLine     = ggl::model::multi_linestring<LineString>;
using Ring         = ggl::model::ring<Point, true >;
using Hole         = ggl::model::ring<Point, false>;
using Polygon      = ggl::model::polygon<Point>;
using MultiPolygon = ggl::model::multi_polygon<Polygon>;
using Path         = LineString;
using MultiPath    = PolyLine;

} // xy

} // farm_db
