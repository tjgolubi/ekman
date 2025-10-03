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
  inline double ResampleDs            = 0.25;  // m
  inline double EvalDs                = 0.25;  // m
  inline double SimplifyTol           = 0.10;  // m
  inline std::size_t MinPtsForSmooth  = 8;     // input vertices in arc
  inline double MaxMeanSpacingForSmooth = 0.75; // m
  inline double MaxEdgeForSmooth        = 3.00; // m
} // Tune

namespace {

inline Pt Lerp(const Pt& a, const Pt& b, double t)  { return a + (b - a) * t; }

// ----------------------- Preconditions -------------------------------
inline
void AssertCornerInvariants(const Ring& r, gsl::index start, gsl::index stop) {
  const auto n     = static_cast<gsl::index>(r.size());
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
void ExtractArc(const Ring& r, gsl::index start, gsl::index stop,
                Linestring& out)
{
  out.clear();
  const std::size_t a = static_cast<std::size_t>(start);
  const std::size_t b = static_cast<std::size_t>(stop);
  out.insert(out.end(), r.begin() + a, r.begin() + b + 1);
}

inline double MaxSegmentLen(const Linestring& ls) {
  if (ls.size() < 2) return 0.0;
  double m = 0.0;
  for (std::size_t i = 1; i < ls.size(); ++i) {
    const double d = Dist(ls[i - 1], ls[i]);
    if (d > m) m = d;
  }
  return m;
} // MaxSegmentLen

inline double TotalLen(const Linestring& ls) {
  if (ls.size() < 2) return 0.0;
  double L = 0.0;
  for (std::size_t i = 1; i < ls.size(); ++i) L += Dist(ls[i - 1], ls[i]);
  return L;
} // TotalLen

inline bool ShouldSmooth(const Linestring& arc) {
  // Dense-only smoothing predicate.
  if (arc.size() < Tune::MinPtsForSmooth) return false;

  const double L = TotalLen(arc);
  if (L <= 0.0) return false;

  const double mean = L / static_cast<double>(arc.size() - 1);
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
  Linestring arc;
  ExtractArc(r, start, stop, arc);
  if (arc.size() < 2) return arc;

  // If arc is not dense, either return it as-is or optionally simplify.
  if (!ShouldSmooth(arc)) {
    Linestring out_sparse;
    // Optional: light simplify on sparse arcs (often long straight edges).
    bg::simplify(arc, out_sparse, Tune::SimplifyTol);
    if (out_sparse.size() >= 2) return out_sparse;
    return arc;
  }

  // 2) Cumulative arclength on the (dense) arc
  std::vector<double> s;
  s.reserve(arc.size());
  s.push_back(0.0);
  for (std::size_t i = 1; i < arc.size(); ++i)
    s.push_back(s.back() + Dist(arc[i - 1], arc[i]));
  const double total = s.back();
  if (total <= 0.0) return arc;

  // 3) Resample to uniform knots
  const std::size_t kcount =
    std::max<std::size_t>(
      2, static_cast<std::size_t>(std::floor(total / Tune::ResampleDs)) + 1);

  std::vector<double> xs, ys;
  xs.reserve(kcount); ys.reserve(kcount);

  std::size_t seg = 0;
  for (std::size_t k = 0; k < kcount; ++k) {
    const double sk = (k == kcount - 1)
      ? total
      : static_cast<double>(k) * Tune::ResampleDs;

    while ((seg + 1) < s.size() && s[seg + 1] < sk) ++seg;

    const double s0 = s[seg], s1 = s[seg + 1];
    const double t  = (s1 > s0) ? (sk - s0) / (s1 - s0) : 0.0;

    const Pt p = Lerp(arc[seg], arc[seg + 1], t);
    xs.push_back(p.x);
    ys.push_back(p.y);
  }

  // 4) Fit cubic B-splines
  const double dx =
    (kcount > 1) ? (total / static_cast<double>(kcount - 1)) : total;

  cardinal_cubic_b_spline<double> sx(xs.begin(), xs.end(), 0.0, dx);
  cardinal_cubic_b_spline<double> sy(ys.begin(), ys.end(), 0.0, dx);

  // 5) Evaluate spline
  Linestring smooth;
  const std::size_t scount =
    std::max<std::size_t>(
      2, static_cast<std::size_t>(std::floor(total / Tune::EvalDs)) + 1);

  smooth.reserve(scount);
  for (std::size_t i = 0; i < scount; ++i) {
    const double si = (i == scount - 1) ? total
                                        : static_cast<double>(i) * Tune::EvalDs;
    Pt p;
    p.x = sx(si);
    p.y = sy(si);
    smooth.push_back(p);
  }

  // 6) Simplify
  Linestring out;
  bg::simplify(smooth, out, Tune::SimplifyTol);
  if (out.size() < 2) return smooth;
  return out;
} // Smooth
