// Smooth.cpp
//
// Smooth dense arcs between corners by:
//   1) uniform resample,
//   2) compute per-step deviation angles via Vec::angle_wrt(local tangent),
//   3) moving-average those small signed deviations,
//   4) clamp each step’s deviation so lateral move ≤ LateralTol,
//   5) integrate headings at constant step,
//   6) optional DP simplify.
//
// Contracts:
//   - The input ring is CLOSED: r.front() == r.back().
//   - Caller provides a non-wrapping contiguous slice [start..stop] (indices in the unique range).
//
// Requires: geom.hpp (Pt, Vec, Radians, Dist, operators), GSL for Expects.

#include "Smooth.hpp"

#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/algorithms/simplify.hpp>
#include <boost/geometry.hpp>
#include <gsl-lite/gsl-lite.hpp>
namespace gsl = gsl_lite;
#include <algorithm>
#include <cstddef>
#include <vector>
#include <cmath>

namespace bg = boost::geometry;

using geom::Pt;
using geom::Vec;
using geom::Radians;

using Ring       = bg::model::ring<Pt>;
using Linestring = bg::model::linestring<Pt>;

// ----------------------- Tunables -----------------------------------
namespace Tune {
  // Resampling step (m). Typical: 0.2–0.5
  constexpr double ResampleDs              = 0.50;

  // Post-smoothing DP tolerance (m). Typical: 0.05–0.15
  constexpr double SimplifyTol             = 0.10;

  // Predicate to decide "dense" arcs to smooth
  constexpr int         MinPtsForSmooth    = 8;
  constexpr double      MaxMeanSpacingForSmooth = 0.75;  // m
  constexpr double      MaxEdgeForSmooth        = 3.00;  // m

  // Deviation smoothing: moving-average window half-size
  constexpr int HeadingHalfWin             = 3;          // total window = 2*H+1

  // Safety: per-step lateral deviation limit (meters)
  constexpr double LateralTol              = 0.10;
} // Tune

// ----------------------- Helpers ------------------------------------
namespace {

constexpr Pt Lerp(const Pt& a, const Pt& b, double t) noexcept {
  // a + t*(b-a)
  return a + (b - a) * t;
}

void AssertCornerInvariants(const Ring& r, gsl::index start, gsl::index stop) {
  const auto n     = std::ssize(r);
  const auto nuniq = n - 1;
  gsl_Expects(n >= 2);
  gsl_Expects(r.front() == r.back());
  gsl_Expects(nuniq >= 2);
  gsl_Expects(start >= 0 && start < nuniq);
  gsl_Expects(stop >= 0);
  gsl_Expects(stop <= nuniq);
  gsl_Expects(start <= stop);
}

Linestring ExtractArc(const Ring& r, gsl::index start, gsl::index stop) {
  Linestring line;
  line.insert(line.end(), r.begin() + start, r.begin() + stop + 1);
  return line; // OPEN arc [start..stop]
}

double MaxSegmentLen(const Linestring& ls) {
  if (std::ssize(ls) < 2) return 0.0;
  double m = 0.0;
  for (gsl::index i = 1; i < std::ssize(ls); ++i)
    m = std::max(m, geom::Dist(ls[i-1], ls[i]));
  return m;
}

double TotalLen(const Linestring& ls) {
  if (std::ssize(ls) < 2) return 0.0;
  double L = 0.0;
  for (gsl::index i = 1; i < std::ssize(ls); ++i)
    L += geom::Dist(ls[i-1], ls[i]);
  return L;
}

bool ShouldSmooth(const Linestring& arc) {
  // Dense-only smoothing predicate
  if (std::ssize(arc) < Tune::MinPtsForSmooth) return false;

  const double len  = TotalLen(arc);
  if (len <= 0.0) return false;

  const double mean = len / static_cast<double>(std::ssize(arc) - 1);
  if (mean > Tune::MaxMeanSpacingForSmooth) return false;

  if (MaxSegmentLen(arc) > Tune::MaxEdgeForSmooth) return false;

  return true;
}

std::vector<Radians> MovingAverage(const std::vector<Radians>& v, int half_win) {
  if (half_win <= 0) return v;
  const int W = 2 * half_win + 1;
  const int n = static_cast<int>(v.size());
  std::vector<Radians> out; out.reserve(n);
  for (int i = 0; i < n; ++i) {
    Radians s{0.0};
    for (int k = -half_win; k <= half_win; ++k) {
      int j = std::clamp(i + k, 0, n - 1);
      s += v[j];
    }
    out.push_back(s / static_cast<double>(W));
  }
  return out;
}

} // local

// ----------------------- Smoothing ----------------------------------
Linestring Smooth(const Ring& ring, gsl::index start, gsl::index stop) {
  AssertCornerInvariants(ring, start, stop);

  // 1) Extract arc (OPEN: [start..stop])
  auto arc = ExtractArc(ring, start, stop);
  if (std::ssize(arc) < 2) return arc;

  // 2) If not dense, just DP-simplify lightly and return
  if (!ShouldSmooth(arc)) {
    Linestring sparse;
    bg::simplify(arc, sparse, Tune::SimplifyTol);
    return (std::ssize(sparse) >= 2) ? sparse : arc;
  }

  // 3) Cumulative arclength on arc
  std::vector<double> S; S.reserve(arc.size());
  S.push_back(0.0);
  for (gsl::index i = 1; i < std::ssize(arc); ++i)
    S.push_back(S.back() + geom::Dist(arc[i-1], arc[i]));
  const double total = S.back();
  if (total <= 0.0) return arc;

  // 4) Uniform resample to step ≈ ResampleDs → r (OPEN)
  std::vector<Pt> r;
  const int N = std::max<int>(2, static_cast<int>(std::floor(total / Tune::ResampleDs)) + 1);
  r.reserve(N);
  r.emplace_back(arc.front());
  gsl::index seg = 0;
  for (int k = 1; k < N - 1; ++k) {
    const double sk = k * (total / (N - 1));
    while (seg + 1 < std::ssize(S) && S[seg + 1] < sk) ++seg;
    const double s0 = S[seg], s1 = S[seg + 1];
    const double t  = (s1 > s0) ? (sk - s0) / (s1 - s0) : 0.0;
    r.emplace_back(Lerp(arc[seg], arc[seg + 1], t));
  }
  r.emplace_back(arc.back());

  // Constants
  const double ds = total / (N - 1);         // uniform step
  const auto theta_max = geom::asin(std::clamp(Tune::LateralTol / ds, 0.0, 1.0));

  // 5) Base step headings h[i] for each segment (N-1)
  std::vector<Radians> h; h.reserve(N - 1);
  for (int i = 0; i < N - 1; ++i) {
    const Vec v = r[i + 1] - r[i];
    h.push_back(v.angle());
  }

  // 6) Local reference tangents (centered) and small deviation angles alpha[i]
  auto ref_tangent = [&](int i) -> Vec {
    if (i == 0)     return r[1]     - r[0];
    if (i == N - 2) return r[N - 1] - r[N - 2];
    return (r[i + 1] - r[i - 1]);   // centered difference, not normalized
  };

  std::vector<Radians> alpha; alpha.reserve(N - 1);
  for (int i = 0; i < N - 1; ++i) {
    const Vec v = r[i + 1] - r[i];
    const Vec u = ref_tangent(i);
    alpha.push_back(v.angle_wrt(u));        // small signed deviation
  }

  // 7) Smooth the deviations (reduces per-step deflection)
  auto alpha_s = MovingAverage(alpha, Tune::HeadingHalfWin);

  // 8) Clamp per-step deviation so lateral offset ≤ LateralTol
  for (auto& a : alpha_s)
    a = std::clamp(a, -theta_max, +theta_max);

  // 9) New headings h' = h + alpha_s
  std::vector<Radians> hprime; hprime.reserve(N - 1);
  for (int i = 0; i < N - 1; ++i) hprime.push_back(h[i] + alpha_s[i]);

  // 10) Reconstruct positions by integrating h' at constant step
  Linestring smooth; smooth.resize(N);
  smooth[0] = r[0];
  for (int i = 1; i < N; ++i)
    smooth[i] = smooth[i - 1] + Vec{ds, hprime[i - 1]};  // polar Vec(len, angle)
  // Pin exact endpoint
  smooth.back() = r.back();

  // 11) Optional DP simplify (keeps endpoints)
  if (Tune::SimplifyTol > 0.0 && std::ssize(smooth) > 2) {
    Linestring out;
    bg::simplify(smooth, out, Tune::SimplifyTol);
    if (std::ssize(out) >= 2) return out;
  }
  return smooth;
} // Smooth
