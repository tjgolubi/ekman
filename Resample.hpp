#pragma once
#include "geom.hpp"
#include <boost/geometry/geometries/ring.hpp>
#include <vector>

/// Compute uniform-step deflection impulses along a closed ring.
/// The first element is 0 rad (tractor already aligned at ring start).
///
/// The ring is a cyclic polyline P[i]->P[i+1] for i=0..N-2 with
/// P[N-1]==P[0]. The output has one entry per step of length
/// sampleDistance, giving the instantaneous turn (radians) at the start
/// of that step.
///
/// Notes:
/// - The last point is kept to form the final edge P[N-2]->P[N-1]==P[0].
/// - Consecutive duplicate points are removed (Pt::operator==).
/// - Segments with length < sampleDistance are collapsed to reduce
///   jitter during resampling (closure is preserved).
std::vector<geom::Radians>
MakeDeflections(const boost::geometry::model::ring<geom::Pt<double>>& ring,
                double sampleDistance = 1.0 /* meters */);
