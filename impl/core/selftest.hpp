// Start-up self-test for production binaries (SOUNDNESS_REVIEW S2).  Call mub_startup_selftest("name") first in main().
// Aborts (exit 3) unless: (a) FE_UPWARD really changes results (not constant-folded / fesetround ignored),
// (b) the interval constants are non-degenerate (lo < hi) and enclose the true values.
// Pure core/ dependency; adds no state.  Prints one "selftest ok" line on stderr.
#pragma once
#include <cstdio>
#include <cstdlib>
#include <cfenv>
#include "interval.hpp"
#include "mub_ops.hpp"

namespace iv {
#if defined(__GNUC__)
#define MUB_NOINLINE __attribute__((noinline))
#else
#define MUB_NOINLINE
#endif
MUB_NOINLINE inline double selftest_div(double a, double b) { volatile double x = a, y = b; return x / y; }

inline void mub_startup_selftest(const char* who) {
  const char* fail = nullptr;
  int old = std::fegetround();
  std::fesetround(FE_UPWARD);
  if (std::fegetround() != FE_UPWARD) fail = "fesetround(FE_UPWARD) not honoured";
  double up = selftest_div(1.0, 3.0);
  std::fesetround(FE_DOWNWARD);
  double dn = selftest_div(1.0, 3.0);
  std::fesetround(old);
  if (!fail && !(dn < up)) fail = "1/3 under FE_UPWARD equals FE_DOWNWARD (rounding modes ignored or constant-folded)";
  if (!fail) {
    Interval s = one_sixth(), p = pi_iv(), r = sqrt_iv(6);
    Interval t;
    { RoundUp g; volatile double x = 0.1, y = 0.2; t = Interval(x) + Interval(y); }
    if (!(s.lo < s.hi)) fail = "one_sixth degenerate";
    else if (!(p.lo < p.hi)) fail = "pi_iv degenerate";
    else if (!(r.lo < r.hi)) fail = "sqrt_iv(6) degenerate";
    else if (!(t.lo < t.hi)) fail = "Interval(0.1)+Interval(0.2) degenerate";
    else if (!((long double)p.lo < 3.14159265358979323846264338L && 3.14159265358979323846264338L < (long double)p.hi)) fail = "pi_iv does not enclose pi";
    else if (!((long double)r.lo < 2.44948974278317809819728407L && 2.44948974278317809819728407L < (long double)r.hi)) fail = "sqrt_iv(6) does not enclose sqrt 6";
    else if (!((long double)s.lo < 1.0L/6.0L && 1.0L/6.0L < (long double)s.hi)) fail = "one_sixth does not enclose 1/6";
  }
  if (fail) { std::fprintf(stderr, "%s: STARTUP SELFTEST FAILED: %s\n", who, fail); std::exit(3); }
  std::fprintf(stderr, "%s: startup selftest ok (rounding + interval constants)\n", who);
}
}  // namespace iv
