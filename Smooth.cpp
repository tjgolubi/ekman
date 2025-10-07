/// @file
///   Smooth dense arcs (lots of close points) with uniform-arclength
///   resample → Boost.Math cubic B-splines → simplify. Skip smoothing
///   for sparse arcs (e.g., long straight QGIS edges).
///
/// Contracts:
///   - Ring is closed (r.front()==r.back()).
///   - Caller passes a non-wrapping [start..stop] (contiguous slice).
///   - Pt/Vec algebra as discussed (Pt - Pt -> Vec, .norm(), etc.).
///
/// Tuning guide (dense-only smoothing):
///   Tune::ResampleDs
///     Knot spacing for resampling (m). Typical 0.2–0.5.
///   Tune::EvalDs
///     Output spacing (m). Typical same as ResampleDs.
///   Tune::SimplifyTol
///     Post-spline simplify tol (m). Typical 0.05–0.15.
///   Tune::MinPtsForSmooth
///     Require at least this many input vertices in the arc to smooth.
///     Too small → try smoothing noisy shorts; too big → skip legit arcs.
///     Typical: 6–12.
///   Tune::MaxMeanSpacingForSmooth
///     Mean spacing threshold (m): total_length/(N-1) must be ≤ this
///     to be considered “dense”. Typical: ~0.5–1.0 m.
///   Tune::MaxEdgeForSmooth
///     Largest single segment (m) allowed for smoothing. If any edge
///     exceeds this, the arc is “sparse/straight-ish” → skip smoothing.
///     Typical: 2–5 m.

#include "Smooth.hpp"

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

#include <vector>
#include <cmath>
#include <cstddef>

namespace bg = boost::geometry;

using Pt   = geom::Pt;
using Vec  = geom::Vec;
using Ring = bg::model::ring<Pt>;
using Linestring = bg::model::linestring<Pt>;

// ----------------------- Tunables -----------------------------------
namespace Tune {
  constexpr double ResampleDs            = 0.25;  // m
  constexpr double EvalDs                = 0.25;  // m
  constexpr double SimplifyTol           = 0.10;  // m
  constexpr std::size_t MinPtsForSmooth  = 8;     // input vertices in arc
  constexpr double MaxMeanSpacingForSmooth = 0.75; // m
  constexpr double MaxEdgeForSmooth        = 3.00; // m
} // Tune

namespace {

inline Pt Lerp(const Pt& a, const Pt& b, double t)  { return a + (b - a) * t; }

// ----------------------- Preconditions -------------------------------
inline
void AssertCornerInvariants(const Ring& r, gsl::index start, gsl::index stop) {
  const auto n     = std::ssize(r);
  const auto nuniq = n - 1;
  Expects(n >= 2);
  Expects(r.front() == r.back());
  Expects(nuniq >= 2);
  Expects(start >= 0 && start < nuniq);
  Expects(stop  >= 0 && stop  < nuniq);
  Expects(start <= stop);
} // AssertCornerInvariants

// ----------------------- Utilities ----------------------------------
inline
Linestring ExtractArc(const Ring& r, gsl::index start, gsl::index stop) {
  auto line = Linestring{};
  line.insert(line.end(), r.begin() + start, r.begin() + stop + 1);
  return line;
}

inline double MaxSegmentLen(const Linestring& ls) {
  if (std::ssize(ls) < 2) return 0.0;
  auto m = 0.0;
  for (auto i = 1; i != std::ssize(ls); ++i) {
    const auto d = Dist(ls[i-1], ls[i]);
    if (d > m) m = d;
  }
  return m;
} // MaxSegmentLen

inline double TotalLen(const Linestring& ls) {
  if (std::ssize(ls) < 2) return 0.0;
  auto len = 0.0;
  for (auto i = 1; i != std::ssize(ls); ++i)
    len += Dist(ls[i-1], ls[i]);
  return len;
} // TotalLen

inline bool ShouldSmooth(const Linestring& arc) {
  // Dense-only smoothing predicate.
  if (std::ssize(arc) < Tune::MinPtsForSmooth) return false;

  const auto len = TotalLen(arc);
  if (len <= 0.0) return false;

  const auto mean = len / static_cast<double>(std::ssize(arc) - 1);
  if (mean > Tune::MaxMeanSpacingForSmooth) return false;

  if (MaxSegmentLen(arc) > Tune::MaxEdgeForSmooth) return false;

  return true;
} // ShouldSmooth

} // anonymous

// ----------------------- Smoothing ----------------------------------
Linestring Smooth(const Ring& r, gsl::index start, gsl::index stop) {
  using boost::math::interpolators::cardinal_cubic_b_spline;

  AssertCornerInvariants(r, start, stop);

  // 1) Extract arc
  auto arc = ExtractArc(r, start, stop);
  if (std::ssize(arc) < 2) return arc;

  // If arc is not dense, either return it as-is or optionally simplify.
  if (!ShouldSmooth(arc)) {
    auto out_sparse = Linestring{};
    // Optional: light simplify on sparse arcs (often long straight edges).
    bg::simplify(arc, out_sparse, Tune::SimplifyTol);
    return (std::ssize(out_sparse) >= 2) ? out_sparse : arc;
  }

  // 2) Cumulative arclength on the (dense) arc
  auto s = std::vector<double>{};
  s.reserve(arc.size());
  s.push_back(0.0);
  for (auto i = 1; i != std::ssize(arc); ++i)
    s.push_back(s.back() + Dist(arc[i-1], arc[i]));
  const auto total = s.back();
  if (total <= 0.0) return arc;

  // 3) Resample to uniform knots
  const auto kcount =
    std::max<int>(2, static_cast<int>(std::floor(total/Tune::ResampleDs)) + 1);

  std::vector<double> xs, ys;
  xs.reserve(kcount); ys.reserve(kcount);

  auto seg = 0;
  for (auto k = 0; k != kcount; ++k) {
    const auto sk = (k == kcount - 1) ? total : k * Tune::ResampleDs;

    while ((seg + 1) < std::ssize(s) && s[seg + 1] < sk) ++seg;

    const auto s0 = s[seg];
    const auto s1 = s[seg+1];
    const auto t  = (s1 > s0) ? (sk - s0) / (s1 - s0) : 0.0;

    const auto p = Lerp(arc[seg], arc[seg + 1], t);
    xs.push_back(p.x);
    ys.push_back(p.y);
  }

  // 4) Fit cubic B-splines
  const auto dx =
              (kcount > 1) ? (total / static_cast<double>(kcount - 1)) : total;

  auto sx = cardinal_cubic_b_spline<double>{xs.begin(), xs.end(), 0.0, dx};
  auto sy = cardinal_cubic_b_spline<double>{ys.begin(), ys.end(), 0.0, dx};

  // 5) Evaluate spline
  auto smooth = Linestring{};
  const auto scount =
        std::max<int>(2, static_cast<int>(std::floor(total/Tune::EvalDs)) + 1);

  smooth.reserve(scount);
  for (auto i = 0; i != scount-1; ++i) {
    const auto si = i * Tune::EvalDs;
    smooth.emplace_back(sx(si), sy(si));
  }
  smooth.emplace_back(sx(total), sy(total));

  // 6) Simplify
  auto out = Linestring{};
  bg::simplify(smooth, out, Tune::SimplifyTol);
  return (std::ssize(out) < 2) ? smooth : out;
} // Smooth
