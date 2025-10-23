/// @file
/// @copyright 2025 Terry Golubiewski, all rights reserved.
/// @author Terry Golubiewski
#pragma once

#include "FarmGeo.hpp"
#include "FarmXy.hpp"

#include <mp-units/systems/si/units.h>
#include <mp-units/systems/isq/space_and_time.h>

#include <vector>

namespace farm_db {

using mp_units::quantity;
using mp_units::si::metre;

using Distance = quantity<mp_units::isq::distance[metre]>;

constexpr Distance DefaultSimplifyTol = 0.10 * metre;

std::vector<xy::MultiPath>
BoundarySwaths(const xy::Polygon& poly_in, Distance offset,
               Distance simplifyTol = DefaultSimplifyTol);

std::vector<geo::MultiPath>
BoundarySwaths(const geo::Polygon& poly_in, Distance offset,
               Distance simplifyTol = DefaultSimplifyTol);

} // farm_db
