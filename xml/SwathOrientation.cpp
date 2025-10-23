#pragma once
#include "SwathOrientation.hpp"

#include <boost/geometry/algorithms/centroid.hpp>
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/length.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/strategies/transform/rotate_transformer.hpp>
#include <boost/geometry/strategies/transform/translate_transformer.hpp>
#include <boost/geometry.hpp>

#include <algorithm>
#include <compare>
#include <vector>
#include <numbers>
#include <cmath>

namespace farm_db {

namespace units {

using HdgDegInt = mp_units::quantity<heading[deg], int>;
using Length    = mp_units::quantity<isq::length[meter]>;
using Width     = mp_units::quantity<isq::width[meter]>;

} // units

namespace xy {

using Point     = geom::Pt<double>;
using Polygon   = ggl::polygon<Point>;
using Line      = ggl::linestring<Point>;
using MultiLine = ggl::multi_linestring<Line>;

}; // xy

namespace detail {

// ---------- Parameters & Result (UpperCamelCase) ----------
namespace Params {
using namespace mp_units;
using namespace units;
constexpr quantity RowWidth = 22 * isq::width[usc::inch];
constexpr int RowsPerSwath  = 36;
constexpr Width SwathWidth = RowsPerSwath * RowWidth;
constexpr Length PruneLen  = 1.2 * SwathWidth;
constexpr double PhaseStepFrac = 1.0 / 16.0; // phase scan step as fraction of SwathWidth (e.g., 1/16)
constexpr HdgDegInt CoarseStep   =  5 * deg; // coarse angle step
constexpr HdgDegInt RefineSpan  =   3 * deg; // refine ±span around best
constexpr HdgDegInt RefineStep  =   1 * deg; // refine step

}; // Params

struct OrientationChoice {
  using namespace mp_units;
  auto swaths   = 0;      // non-empty lattice lines
  auto gaps     = 0;      // Σ(line_segments_on_that_line - 1)
  auto totalLen = Length{0.0};
  auto angle    = HdgDegInt{0}; // (-90 deg, 90 deg]
  auto phase    = Width{0.0};

  int turns() const
    { return (swaths == 0) ? 0 : (swaths - 1) + 2 * gaps; }

  auto operator<=>(const OrientationChoice& rhs) const {
    auto rval = (swaths > 0) <=> (rhs.swaths > 0);
    if (rval != std::strong_ordering_equal)
      return rval;
    rval =  -turns() <=> -rhs.turns();
    if (rval != std::strong_ordering_equal)
      return rval;
    return totalLen <=> rhs.totalLen;
  } // <=>
}; // OrientationChoice

// ---------- Internal: rotate polygon about its centroid by angle (radians) ----------
xy::Polygon RotateAboutCentroid(const xy::Polygon& in, HdgDegInt angle) {
  auto center = ggl::centroid_return<Point>(ggl::return_envelope<Box>(in));

  constexpr auto PtDim = std::size_t{2};
  using PtRep = Point::value_type;
  using Translator =
          ggl::strategy::transform::translate_transformer<PtRep, PtRep, PtDim>;
  using Rotator =
ggl::strategy::transform::rotate_transformer<ggl::degree, PtRep, PtDim, PtDim>;

  auto translator = Translator{-center.x, -center.y};
  auto rotator    = Rotator{angle.numerical_value_in(si::degree)};

  xy::Polygon temp;
  ggl::transform(in, temp, translator);
  xy::Polygon out;
  ggl::transform(temp,  out, rotator);
  return out;
} // RotateAboutCentroid

// ---------- Evaluate a single angle (rotate-by -θ, use horizontal lattice) ----------
OrientationChoice EvalAngle(const xy::Polygon& poly, HdgDegInt angle) {
  // Rotate polygon by -θ so swaths at θ become horizontal (y = const).
  // (Rationale: simpler numerics; reuse same rotated geometry across phase scan.)
  auto rot = RotateAboutCentroid(poly, -angle);

  const auto extents = ggl::return_envelope<Box>(rot);
  const auto minPt   = extents.min_corner();
  const auto maxPt   = extents.max_corner();

  const auto swathWidth = Params::SwathWidth.numerical_value_in(si::metre);

  auto best = OrientationChoice{};
  const int nPhase =
        std::max(1, static_cast<int>(std::round(1.0 / Params::PhaseStepFrac)));

  for (int ip = 0; ip != nPhase; ++ip) {
    const auto baseY = minPt.y + std::fmod(phase + swathWidth, swathWidth);
    const auto lines = static_cast<int>(std::floor((maxPt.y - minPt.y) / swathWidth)) + 4;

    auto cand = OrientationChoice{
      .angle = angle,
      .phase = (ip * Params::PhaseStepFrac) * Params::SwathWidth;
    };

    for (int k = -2; k <= lines + 2; ++k) {
      const auto y = baseY + k * swathWidth;

      auto scan = Line{Point{minPt.x - 10 * swathWidth, y}, Point{maxPt.x + 10 * swathWidth, y}};

      MultiLine fragments;
      (void) ggl::intersection(rot, scan, fragments);

      auto keptOnThisLine = 0;
      for (auto& s : fragments) {
        const auto len = Length{ggl::length(s) * si::metre};
        if (len >= Params::PruneLen) {
          ++keptOnThisLine;
          cand.totalLen += len;
        }
      }

      if (keptOnThisLine > 0) {
        ++cand.swaths;
        cand.gaps += (keptOnThisLine - 1);
      }
    }

    if (cand > best)
      best = cand;
  }

  return best;
} // EvalAngle

// ---------- Public: choose best angle + phase (coarse sweep + refine) ----------
std::pair<HdgDeg, Length> ChooseSwathOrientation(const xy::Polygon& interior) {
  using units::deg;

  auto global = OrientationChoice{};

  // Cardinal sweep
  auto cardinals = std::array{ 0 * deg, 90 * deg, 45 * deg, -45 * deg};
  for (auto angle: cardinals) {
    auto r = EvalAngle(interior, angle);
    if (r > global)
      global = r;
  }
  // Coarse sweep
  for (auto angle = -90 * deg; angle <= 90 * deg; angle += Params::CoarseStep) {
    if (angle % 45 * deg == 0 * deg)
      continue;
    auto r = EvalAngle(interior, angle);
    if (r > global)
      global = r;
  }

  // Refine around best
  const auto lo = global.angle - Params::RefineSpan;
  const auto hi = global.angle + Params::RefineSpan;

  for (auto angle = lo; angle <= hi; angle += Params::RefineStep) {
    auto r = EvalAngle(interior, angle);
    if (r > global)
      global = r;
  }

  return global;
} // ChooseSwathOrientation

} // farm_db
