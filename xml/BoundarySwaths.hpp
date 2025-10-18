/// @file
/// @copyright 2025 Terry Golubiewski, all rights reserved.
/// @author Terry Golubiewski
#pragma once

#include "FarmGeo.hpp"

#include "../geom.hpp"

#include <mp-units/systems/si/units.h>
#include <mp-units/systems/isq/space_and_time.h>

#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/core/cs.hpp>

#include <vector>

namespace farm_db {

using mp_units::quantity;
using mp_units::si::metre;

using Meters   = double;
using XyPt       = geom::Pt<Meters>;
using XyPolygon  = ggl::model::polygon<XyPt>;
using XyPath     = ggl::model::linestring<XyPt>;
using XyMultiPath= ggl::model::multi_linestring<XyPath>;

using GeoPt        = GeoPoint;
using GeoPath      = ggl::model::linestring<GeoPt>;
using GeoMultiPath = ggl::model::multi_linestring<GeoPath>;

using Distance = quantity<mp_units::isq::distance[metre]>;

constexpr Distance DefaultSimplifyTol = 0.10 * metre;

std::vector<XyMultiPath>
BoundarySwaths(const XyPolygon& poly_in, Distance offset,
               Distance simplifyTol = DefaultSimplifyTol);

std::vector<GeoMultiPath>
BoundarySwaths(const GeoPolygon& poly_in, Distance offset,
               Distance simplifyTol = DefaultSimplifyTol);

} // farm_db
