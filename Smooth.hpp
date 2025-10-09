/// @file
///   Smooth dense arcs (lots of close points) with uniform-arclength
///   resample → Boost.Math cubic B-splines → simplify. Skip smoothing
///   for sparse arcs (e.g., long straight QGIS edges).
///
/// Contracts:
///   - Ring is closed (r.front()==r.back()).
///   - Caller passes a non-wrapping [start..stop] (contiguous slice).
///   - Pt/Vec algebra as discussed (Pt - Pt -> Vec, .norm(), etc.).

#pragma once
#include "geom.hpp"
#include <boost/geometry.hpp>
#include <gsl-lite/gsl-lite.hpp>

boost::geometry::model::linestring<geom::Pt<double>>
Smooth(const boost::geometry::model::ring<geom::Pt<double>>& r,
       gsl_lite::index start, gsl_lite::index stop);
