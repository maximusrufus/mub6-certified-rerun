// arith.hpp -- arithmetic adapter for the seed stage, built on the shared
// outward-rounded interval core (../core/interval.hpp).
//
// CERTIFICATE ARITHMETIC: every certified inequality in the seed stage is a
// strict interval comparison on enclosures produced by the core.  There are NO
// additive paper margins (1e-8, 1e-7, 1e-5, .01) and no 1.251 inflation.
// Contract: every thread that evaluates a certificate must hold an
// iv::RoundUp guard (FE_UPWARD); plain double operations performed under that
// guard on NONNEGATIVE operands are valid upper bounds, which the hot paths
// use (marked "up:" in comments).  Build flags: -frounding-math
// -ffp-contract=off, never -ffast-math.
//
// The double-precision helpers at the bottom (PI, expair) are HEURISTIC ONLY:
// they propose normal forms / dual vectors / test Hadamard matrices.  No
// certificate depends on them.
#pragma once
#include <cmath>
#include <cstdint>
#include "../core/interval.hpp"

namespace arith {

using iv::CInterval;
using iv::Interval;
using iv::RoundUp;

// ---- certified building blocks (call under a RoundUp guard) -----------------

// exact point (double) times interval; FE_UPWARD must be active.
inline Interval mulp(const Interval &a, double v) {
  if (v >= 0) return Interval(-((-a.lo) * v), a.hi * v);
  return Interval(-((-a.hi) * v), a.lo * v);
}

// upper bound of |a|
inline double abs_hi(const Interval &a) { return std::max(-a.lo, a.hi); }
// lower bound of |a|
inline double abs_lo(const Interval &a) { return a.lo > 0 ? a.lo : (a.hi < 0 ? -a.hi : 0.0); }

// e((b+1/2)/N) = exp(2 pi i (2b+1)/(2N)) : enclosure of the exact rational midpoint phase.
inline CInterval cis_mid(long b, long N) { return iv::unit_root(2 * b + 1, 2 * N); }

// ---- HEURISTIC double helpers (never used in a certified inequality) --------
static const double PI = 0x1.921fb54442d18p+1;

inline void expair(double x, double &c, double &s) {
  double k = std::nearbyint(4.0 * x);
  double y = 2.0 * PI * (x - k * 0.25);
  double y2 = y * y;
  double cs = 1.0, term = 1.0;
  for (int m = 1; m <= 9; ++m) { term = -term * y2 / (double)((2 * m - 1) * (2 * m)); cs += term; }
  double sn = y, t2 = y;
  for (int m = 1; m <= 9; ++m) { t2 = -t2 * y2 / (double)((2 * m) * (2 * m + 1)); sn += t2; }
  long q = (long)k;
  int r = (int)(((q % 4) + 4) % 4);
  switch (r) {
    case 0: c = cs; s = sn; break;
    case 1: c = -sn; s = cs; break;
    case 2: c = -cs; s = -sn; break;
    default: c = sn; s = -cs; break;
  }
}

}  // namespace arith
