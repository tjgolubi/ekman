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
#include <gsl-lite/gsl-lite.hpp>
namespace gsl = gsl_lite;
#include <algorithm>
#include <cstddef>
#include <vector>
#include <cmath>

namespace bg = boost::geometry;

using geom::Radians;

using Meters     = double;
using Pt         = geom::Pt<Meters>;
using Disp       = geom::Vec<Meters>;
using Ring       = bg::model::ring<Pt>;
using Linestring = bg::model::linestring<Pt>;

// ----------------------- Tunables -----------------------------------
namespace Tune {
  // Resampling step (m). Typical: 0.2–0.5
  constexpr Meters ResampleDs  = 0.50;

  // Post-smoothing DP tolerance (m). Typical: 0.05–0.15
  constexpr Meters SimplifyTol = 0.10;

  // Predicate to decide "dense" arcs to smooth
  constexpr int    MinPtsForSmooth = 8;
  constexpr Meters MaxMeanSpacingForSmooth = 0.75;  // m
  constexpr Meters MaxEdgeForSmooth        = 3.00;  // m

  // Deviation smoothing: moving-average window half-size
  constexpr int HeadingHalfWin             = 3;          // total window = 2*H+1

  // Safety: per-step lateral deviation limit (meters)
  constexpr Meters LateralTol              = 0.10;
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

Meters MaxSegmentLen(const Linestring& ls) {
  if (std::ssize(ls) < 2) return 0.0;
  auto m = Meters{0.0};
  for (gsl::index i = 1; i != std::ssize(ls); ++i)
    m = std::max(m, geom::Dist(ls[i-1], ls[i]));
  return m;
}

Meters TotalLen(const Linestring& ls) {
  if (std::ssize(ls) < 2) return Meters{0};
  auto len = Meters{0};
  for (gsl::index i = 1; i != std::ssize(ls); ++i)
    len += geom::Dist(ls[i-1], ls[i]);
  return len;
}

bool ShouldSmooth(const Linestring& arc) {
  // Dense-only smoothing predicate
  if (std::ssize(arc) < Tune::MinPtsForSmooth) return false;

  const auto len  = TotalLen(arc);
  if (len <= Meters{0}) return false;

  const auto mean = len / (std::ssize(arc) - 1);
  if (mean > Tune::MaxMeanSpacingForSmooth) return false;

  if (MaxSegmentLen(arc) > Tune::MaxEdgeForSmooth) return false;

  return true;
} // ShouldSmooth

std::vector<Radians> MovingAverage(const std::vector<Radians>& v, int half_win) {
  if (half_win <= 0) return v;
  const auto w = 2 * half_win + 1;
  const auto n = static_cast<int>(std::ssize(v));
  std::vector<Radians> out; out.reserve(n);
  for (int i = 0; i != n; ++i) {
    auto s = Radians{0};
    for (int k = -half_win; k <= half_win; ++k) {
      auto j = std::clamp(i + k, 0, n - 1);
      s += v[j];
    }
    out.push_back(s/w);
  }
  return out;
} // MovingAverage

} // local

// ----------------------- Smoothing ----------------------------------
Linestring Smooth(const Ring& ring, gsl::index start, gsl::index stop) {
  AssertCornerInvariants(ring, start, stop);

  // 1) Extract arc (OPEN: [start..stop])
  auto arc = ExtractArc(ring, start, stop);
  if (std::ssize(arc) < 2) return arc;

  // 2) If not dense, just DP-simplify lightly and return
  if (!ShouldSmooth(arc)) {
    auto sparse = Linestring{};
    bg::simplify(arc, sparse, Tune::SimplifyTol);
    return (std::ssize(sparse) >= 2) ? sparse : arc;
  }

  // 3) Cumulative arclength on arc
  std::vector<Meters> S; S.reserve(arc.size());
  S.push_back(Meters{0});
  for (gsl::index i = 1; i != std::ssize(arc); ++i)
    S.push_back(S.back() + geom::Dist(arc[i-1], arc[i]));
  const auto total = S.back();
  if (total <= Meters{0}) return arc;

  // 4) Uniform resample to step ≈ ResampleDs → r (OPEN)
  const auto N =
    std::max<int>(2, static_cast<int>(std::floor(total/Tune::ResampleDs)) + 1);
  auto r = std::vector<Pt>{}; r.reserve(N);
  r.emplace_back(arc.front());
  auto seg = gsl::index{0};
  for (int k = 1; k != N - 1; ++k) {
    const auto sk = k * (total / (N-1));
    while (seg + 1 < std::ssize(S) && S[seg+1] < sk) ++seg;
    const auto s0 = S[seg];
    const auto s1 = S[seg+1];
    const auto t  = (s1 > s0) ? (sk - s0) / (s1 - s0) : 0.0;
    r.emplace_back(Lerp(arc[seg], arc[seg+1], t));
  }
  r.emplace_back(arc.back());

  // Constants
  const auto ds = total / (N-1);         // uniform step
  const auto theta_max = geom::asin(std::clamp(Tune::LateralTol/ds, 0.0, 1.0));

  // 5) Base step headings h[i] for each segment (N-1)
  auto h = std::vector<Radians>{}; h.reserve(N-1);
  for (int i = 0; i != N-1; ++i) {
    const auto v = r[i+1] - r[i];
    h.push_back(v.angle());
  }

  // 6) Local reference tangents (centered) and small deviation angles alpha[i]
  auto ref_tangent = [&](int i) -> Disp {
    if (i == 0)     return r[1] - r[0];
    if (i == N - 2) return r[N-1] - r[N-2];
    return (r[i+1] - r[i-1]);   // centered difference, not normalized
  };

  std::vector<Radians> alpha; alpha.reserve(N - 1);
  for (int i = 0; i < N - 1; ++i) {
    const auto v = r[i+1] - r[i];
    const auto u = ref_tangent(i);
    alpha.push_back(v.angle_wrt(u));        // small signed deviation
  }

  // 7) Smooth the deviations (reduces per-step deflection)
  auto alpha_s = MovingAverage(alpha, Tune::HeadingHalfWin);

  // 8) Clamp per-step deviation so lateral offset ≤ LateralTol
  for (auto& a : alpha_s)
    a = std::clamp(a, -theta_max, +theta_max);

  // 9) New headings h' = h + alpha_s
  auto hprime = std::vector<Radians>{}; hprime.reserve(N-1);
  for (int i = 0; i != N-1; ++i) hprime.push_back(h[i] + alpha_s[i]);

  // 10) Reconstruct positions by integrating h' at constant step
  auto smooth = Linestring{}; smooth.resize(N);
  smooth[0] = r[0];
  for (int i = 1; i != N; ++i)
    smooth[i] = smooth[i-1] + Disp{ds, hprime[i-1]};  // polar Disp(len, angle)
  // Pin exact endpoint
  smooth.back() = r.back();

  // 11) Optional DP simplify (keeps endpoints)
  if (Tune::SimplifyTol > Meters{0} && std::ssize(smooth) > 2) {
    auto out = Linestring{};
    bg::simplify(smooth, out, Tune::SimplifyTol);
    if (std::ssize(out) >= 2) return out;
  }
  return smooth;
} // Smooth
