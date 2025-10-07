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
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace bg = boost::geometry;

using Pt  = geom::Pt;
using Vec = geom::Vec;
using Radians = geom::Radians;
using Ring = bg::model::ring<Pt>;
using Linestring = bg::model::linestring<Pt>;

// ----------------------- Tunables -----------------------------------
namespace Tune {
  constexpr double ResampleDs              = 0.25;  // m
  constexpr double EvalDs                  = 0.25;  // m (kept for uniform step)
  constexpr double SimplifyTol             = 0.10;  // m
  constexpr std::size_t MinPtsForSmooth    = 8;
  constexpr double MaxMeanSpacingForSmooth = 0.75;  // m
  constexpr double MaxEdgeForSmooth        = 3.00;  // m
  constexpr int HeadingHalfWin             = 3;     // try 3–5
} // Tune

namespace {

constexpr Pt Lerp(const Pt& a, const Pt& b, double t)
  { return a + (b - a) * t; }

// ----------------------- Preconditions -------------------------------
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

double MaxSegmentLen(const Linestring& ls) {
  if (std::ssize(ls) < 2) return 0.0;
  auto m = 0.0;
  for (auto i = 1; i != std::ssize(ls); ++i) {
    const auto d = Dist(ls[i-1], ls[i]);
    if (d > m) m = d;
  }
  return m;
} // MaxSegmentLen

double TotalLen(const Linestring& ls) {
  if (std::ssize(ls) < 2) return 0.0;
  auto len = 0.0;
  for (auto i = 1; i != std::ssize(ls); ++i)
    len += Dist(ls[i-1], ls[i]);
  return len;
} // TotalLen

bool ShouldSmooth(const Linestring& arc) {
  // Dense-only smoothing predicate.
  if (std::ssize(arc) < Tune::MinPtsForSmooth) return false;

  const auto len = TotalLen(arc);
  if (len <= 0.0) return false;

  const auto mean = len / static_cast<double>(std::ssize(arc) - 1);
  if (mean > Tune::MaxMeanSpacingForSmooth) return false;

  if (MaxSegmentLen(arc) > Tune::MaxEdgeForSmooth) return false;

  return true;
} // ShouldSmooth

std::vector<Radians> MovingAverage(const std::vector<Radians>& v, int half_win) {
  if (half_win <= 0) return v;
  const auto w = 2*half_win + 1;
  const auto n = static_cast<int>(std::ssize(v));
  auto out = std::vector<Radians>{};
  out.reserve(n);
  for (int i = 0; i != n; ++i) {
    auto s = Radians{0.0};
    for (int k = -half_win; k <= half_win; ++k) {
      int j = std::clamp(i + k, 0, n-1);
      s += v[j];
    }
    out.push_back(s/w);
  }
  return out;
} // MovingAverage

} // anonymous

Linestring Smooth(const Ring& ring, gsl::index start, gsl::index stop) {
  AssertCornerInvariants(ring, start, stop);

  // 1) Extract arc
  auto arc = ExtractArc(ring, start, stop);
  if (std::ssize(arc) < 2) return arc;

  // Sparse → keep (optionally lightly simplify)
  if (!ShouldSmooth(arc)) {
    auto out_sparse = Linestring{};
    bg::simplify(arc, out_sparse, Tune::SimplifyTol);
    return (std::ssize(out_sparse) >= 2) ? out_sparse : arc;
  }

  // 2) Cumulative arclength
  auto s = std::vector<double>{};
  s.reserve(arc.size());
  s.push_back(0.0);
  for (auto i = 1; i != std::ssize(arc); ++i)
    s.push_back(s.back() + Dist(arc[i-1], arc[i]));
  const auto total = s.back();
  if (total <= 0.0) return arc;

  // 3) Resample to uniform spacing
  auto r = std::vector<Pt>{};
  {
    const auto n = std::max<int>(2,
                   static_cast<int>(std::floor(total / Tune::ResampleDs)) + 1);
    r.reserve(n);
    r.emplace_back(arc.front());
    auto seg = 0;
    for (auto k = 1; k != n-1; ++k) {
      const auto sk = k * (total / (n-1));
      while (seg + 1 < std::ssize(s) && s[seg+1] < sk) ++seg;
      const auto s0 = s[seg];
      const auto s1 = s[seg+1];
      const auto t  = (s1 > s0) ? (sk - s0) / (s1 - s0) : 0.0;
      r.emplace_back(Lerp(arc[seg], arc[seg+1], t));
    }
    r.emplace_back(arc.back());
  }

  // 4) Compute per-step headings on the resampled chain
  auto theta = std::vector<Radians>{};
  {
    const auto n = std::ssize(r) - 1;
    theta.reserve(n);
    for (auto i = 0; i != n; ++i) {
      const auto v = r[i+1] - r[i];
      theta.push_back(v.angle());
    }
  }

  // 5) Smooth the headings directly (reduces per-vertex deflection)
  auto theta_s = MovingAverage(theta, Tune::HeadingHalfWin);

  // 6) Reconstruct positions by integrating the smoothed headings
  auto smooth = Linestring{};
  {
    const auto n = std::ssize(r);
    smooth.resize(n);
    smooth[0] = r[0];
    const double step = total / (n - 1);
    for (auto i = 1; i != n; ++i)
      smooth[i] = smooth[i-1] + Vec{step, theta_s[i-1]};
    // force exact endpoint
    smooth.back() = r.back();
  }

  // 7) Optional post-simplify (keeps endpoints)
  if (Tune::SimplifyTol > 0.0 && std::ssize(smooth) > 2) {
    auto out = Linestring{};
    bg::simplify(smooth, out, Tune::SimplifyTol);
    if (std::ssize(out) >= 2) return out;
  }
  return smooth;
}
