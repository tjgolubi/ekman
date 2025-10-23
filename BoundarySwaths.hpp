/// @file
/// @copyright 2025 Terry Golubiewski, all rights reserved.
/// @author Terry Golubiewski
#pragma once

#include "geom_ggl.hpp"

#include <mp-units/systems/si/units.h>
#include <mp-units/systems/isq/space_and_time.h>

#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/core/cs.hpp>

#include <vector>

namespace tjg {

using mp_units::quantity;
using mp_units::si::metre;

namespace ggl = boost::geometry;
using Meters   = double;
using Pt       = geom::Pt<Meters>;
using Polygon  = ggl::model::polygon<Pt>;
using Path     = ggl::model::linestring<Pt>;
using MultiPath= ggl::model::multi_linestring<Path>;

using Degrees = double;
using GeoPt      = ggl::model::point<Degrees, 2, ggl::cs::geographic<ggl::degree>>;
using GeoPolygon   = ggl::model::polygon<GeoPt>;
using GeoPath      = ggl::model::linestring<GeoPt>;
using GeoMultiPath = ggl::model::multi_linestring<GeoPath>;

using Distance = quantity<mp_units::isq::distance[metre]>;

constexpr Distance DefaultSimplifyTol = 0.10 * metre;

std::vector<MultiPath>
BoundarySwaths(const Polygon& poly_in, Distance offset,
               Distance simplifyTol = DefaultSimplifyTol);

std::vector<GeoMultiPath>
BoundarySwaths(const GeoPolygon& poly_in, Distance offset,
               Distance simplifyTol = DefaultSimplifyTol);

} // tjg
