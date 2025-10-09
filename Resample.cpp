#include "Resample.hpp"

#include <gsl-lite/gsl-lite.hpp>
#include <iterator>  // std::ssize
#include <cmath>     // std::floor

namespace gsl = gsl_lite;

using Meters = double;
using Pt   = geom::Pt<Meters>;
using Disp = geom::Vec<Meters>;
using geom::Radians;
using Ring = boost::geometry::model::ring<Pt>;

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
std::vector<Radians>
MakeDeflections(const Ring& ring, Meters sampleDistance) {
  using gsl::index;

  gsl_Expects(sampleDistance > Meters{0});
  gsl_Expects(std::ssize(ring) >= 2);
  gsl_Expects(ring.front() == ring.back());

  // 1) Copy while collapsing consecutive duplicates. KEEP last point.
  auto P = std::vector<Pt>{};
  for (const Pt& p : ring) {
    if (!P.empty() && p == P.back())
      continue;
    P.push_back(p);
  }
  if (!(P.front() == P.back()))
    P.push_back(P.front());
  gsl_Expects(std::ssize(P) >= 3);

  // 2) Collapse short segments (< sampleDistance), preserving closure.
  auto C = std::vector<Pt>{};
  C.push_back(P.front());
  for (index i = 1; i != std::ssize(P)-1; ++i) {
    const auto& cand = P[i];
    const auto edgeLen = geom::Dist(C.back(), cand);
    if (edgeLen < sampleDistance)
      continue;
    C.push_back(cand);
  }
  if (C.back() != P.back())
    C.push_back(P.back());
  if (C.front() != C.back())
    C.push_back(C.front());
  gsl_Expects(std::ssize(C) >= 3);

  // 3) Build segments C[i]->C[i+1] and their lengths.
  const index Npt = std::ssize(C);
  auto V = std::vector<Disp>{};
  auto L = std::vector<Meters>{};
  for (index i = 0; i != Npt-1; ++i) {
    auto v = C[i+1] - C[i];
    const auto len = v.norm();
    if (len == Meters{0})
      continue;
    V.push_back(v);
    L.push_back(len);
  }
  gsl_Expects(!V.empty());

  // 4) Cumulative arc-length S and perimeter Ltot.
  auto S = std::vector<Meters>{};
  S.push_back(Meters{0});
  auto Ltot = Meters{0};
  for (auto len : L) {
    Ltot += len;
    S.push_back(Ltot);
  }
  gsl_Expects(Ltot > Meters{0});

  // 5) Number of uniform samples; do not place one at the perimeter.
  const auto d = sampleDistance;
  const auto N = static_cast<index>(std::floor(Ltot / d));
  if (N <= 0)
    return {};

  // Helper: locate point at arc-length s. Maintain segment cursor j with
  // S[j] <= s < S[j+1].
  auto locatePoint = [&](Meters s, index& j, Pt& out) {
    constexpr auto eps = 1.0e-12;
    for ( ; j < std::ssize(S)-1; ++j) {
      if (s < S[j+1] - eps * L[j])
        break;
    }
    if (j >= std::ssize(S)-1)
      j = std::ssize(S) - 2;

    const auto segStart = S[j];
    const auto Lj  = L[j];
    const auto& Aj = C[j];
    const auto& Vj = V[j];

    double t = (Lj > Meters{0}) ? (s - segStart) / Lj : 0.0;
    if (t <= eps)
      t = 0.0;
    if (t >= 1.0 - eps) {
      // Snap to vertex; advance for subsequent queries.
      out = C[++j];
    } else {
      out = Aj + (t * Vj);
    }
  }; // locatePoint

  // 6) Resample points and compute deflection impulses.
  auto Qprev = Pt{};
  auto Qcur  = Pt{};
  index j = 0;

  locatePoint(Meters{0}, j, Qprev);
  locatePoint(d, j, Qcur);
  auto Wprev = Qcur - Qprev;

  auto defl = std::vector<Radians>{static_cast<std::size_t>(N)};
  defl[0] = Radians{0.0};
  for (index i = 1; i != N; ++i) {
    Qprev = Qcur;
    locatePoint(i * d, j, Qcur);
    auto Wcur = Qcur - Qprev;
    defl[i] = Wcur.angle_wrt(Wprev);
    Wprev = Wcur;
  }

  return defl;
} // MakeDeflections
