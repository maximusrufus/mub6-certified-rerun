// Header-only outward-rounded real and complex intervals (binary64).
// Requires: -frounding-math -ffp-contract=off, NEVER -ffast-math.
// All arithmetic assumes the thread's FP rounding mode is FE_UPWARD; create a
// RoundUp guard in EVERY thread (inside each OpenMP parallel region) before use.
// Lower bounds use the negation trick: down(x op y) = -up((-x) op' y).
#pragma once
#include <algorithm>
#include <cassert>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#ifdef __FAST_MATH__
#error "interval.hpp must not be compiled with -ffast-math"
#endif

namespace iv {

// Per-thread rounding guard: sets FE_UPWARD, restores the previous mode.
// fenv state is per-thread, so each OpenMP thread must construct its own.
struct RoundUp {
  int old_;
  RoundUp() : old_(std::fegetround()) { std::fesetround(FE_UPWARD); }
  ~RoundUp() { std::fesetround(old_); }
  RoundUp(const RoundUp&) = delete;
  RoundUp& operator=(const RoundUp&) = delete;
};

inline void require_upward() {
#ifdef IV_CHECK_ROUNDING
  if (std::fegetround() != FE_UPWARD) throw std::logic_error("interval op outside RoundUp guard");
#endif
}

struct Interval {
  double lo, hi;
  Interval() : lo(0), hi(0) {}
  explicit Interval(double x) : lo(x), hi(x) {}  // exact point
  Interval(double l, double h) : lo(l), hi(h) { assert(l <= h); }
  double mid() const { return 0.5 * lo + 0.5 * hi; }  // heuristic only
  double width() const { return hi - lo; }            // heuristic only
  bool contains(double x) const { return lo <= x && x <= hi; }
  bool contains_zero() const { return lo <= 0.0 && 0.0 <= hi; }
  bool subset_of(const Interval& o) const { return o.lo <= lo && hi <= o.hi; }
};

inline Interval operator-(const Interval& a) { return Interval(-a.hi, -a.lo); }
inline Interval operator+(const Interval& a, const Interval& b) {
  require_upward();
  return Interval(-((-a.lo) - b.lo), a.hi + b.hi);
}
inline Interval operator-(const Interval& a, const Interval& b) {
  require_upward();
  return Interval(-(b.hi - a.lo), a.hi - b.lo);
}
inline Interval operator*(const Interval& a, const Interval& b) {
  require_upward();
  double h = std::max(std::max(a.lo * b.lo, a.lo * b.hi), std::max(a.hi * b.lo, a.hi * b.hi));
  double nl = std::max(std::max((-a.lo) * b.lo, (-a.lo) * b.hi), std::max((-a.hi) * b.lo, (-a.hi) * b.hi));
  return Interval(-nl, h);
}
inline Interval sqr(const Interval& a) {
  require_upward();
  if (a.lo >= 0) return Interval(-((-a.lo) * a.lo), a.hi * a.hi);
  if (a.hi <= 0) return Interval(-((-a.hi) * a.hi), a.lo * a.lo);
  return Interval(0.0, std::max(a.lo * a.lo, a.hi * a.hi));
}
inline Interval recip(const Interval& b) {
  require_upward();
  if (b.lo <= 0 && b.hi >= 0) throw std::domain_error("interval division by interval containing 0");
  // down(1/hi) = -up((-1)/hi); up(1/lo)
  return Interval(-((-1.0) / b.hi), 1.0 / b.lo);
}
inline Interval operator/(const Interval& a, const Interval& b) { return a * recip(b); }
inline Interval sqrt(const Interval& a) {
  require_upward();
  if (a.lo < 0) throw std::domain_error("interval sqrt of negative");
  double h = std::sqrt(a.hi);  // upward-rounded
  double l = std::sqrt(a.lo);  // upward-rounded; true floor is at most 1 ulp lower
  if (l > 0) l = std::nextafter(l, 0.0);
  return Interval(l, h);
}
inline Interval abs(const Interval& a) {
  if (a.lo >= 0) return a;
  if (a.hi <= 0) return -a;
  return Interval(0.0, std::max(-a.lo, a.hi));
}

// Certain comparisons: true only if the relation holds for every pair of points.
inline bool certainly_lt(const Interval& a, const Interval& b) { return a.hi < b.lo; }
inline bool certainly_le(const Interval& a, const Interval& b) { return a.hi <= b.lo; }
inline bool certainly_gt(const Interval& a, const Interval& b) { return a.lo > b.hi; }
inline bool certainly_ge(const Interval& a, const Interval& b) { return a.lo >= b.hi; }
inline bool disjoint(const Interval& a, const Interval& b) { return a.hi < b.lo || b.hi < a.lo; }
inline bool overlaps(const Interval& a, const Interval& b) { return !disjoint(a, b); }

// ---- complex ----
struct CInterval {
  Interval re, im;
  CInterval() {}
  CInterval(const Interval& r, const Interval& i) : re(r), im(i) {}
  explicit CInterval(const Interval& r) : re(r), im(0.0) {}
  CInterval(double r, double i) : re(r), im(i) {}
};
inline CInterval operator-(const CInterval& a) { return CInterval(-a.re, -a.im); }
inline CInterval conj(const CInterval& a) { return CInterval(a.re, -a.im); }
inline CInterval operator+(const CInterval& a, const CInterval& b) { return CInterval(a.re + b.re, a.im + b.im); }
inline CInterval operator-(const CInterval& a, const CInterval& b) { return CInterval(a.re - b.re, a.im - b.im); }
inline CInterval operator*(const CInterval& a, const CInterval& b) {
  return CInterval(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
inline CInterval operator*(const CInterval& a, const Interval& s) { return CInterval(a.re * s, a.im * s); }
inline Interval abs2(const CInterval& a) { return sqr(a.re) + sqr(a.im); }
inline CInterval operator/(const CInterval& a, const CInterval& b) {
  Interval d = abs2(b);
  if (!(d.lo > 0)) throw std::domain_error("complex interval division by interval containing 0");
  CInterval n = a * conj(b);
  Interval r = recip(d);
  return CInterval(n.re * r, n.im * r);
}
inline bool disjoint(const CInterval& a, const CInterval& b) { return disjoint(a.re, b.re) || disjoint(a.im, b.im); }

// ---- rigorous constants (computed, never rounded doubles) ----
inline Interval pi_iv() {
  // pi_d = 0x1.921fb54442d18p+1 < pi with error < 1.23e-16 < ulp, so [pi_d, nextafter(pi_d,+inf)] encloses pi.
  const double pd = 0x1.921fb54442d18p+1;
  return Interval(pd, std::nextafter(pd, 1e300));
}
inline Interval from_int(long n) { return Interval((double)n); }  // |n| < 2^53
inline Interval ratio(long p, long q) { RoundUp g; return from_int(p) / from_int(q); }
inline Interval one_sixth() { return ratio(1, 6); }
inline Interval sqrt_iv(long n) { RoundUp g; return iv::sqrt(from_int(n)); }
inline Interval inv_sqrt6() { RoundUp g; return Interval(1.0) / sqrt_iv(6); }

namespace detail {
constexpr int TAYLOR_M = 12;
// 1/k! as intervals, k = 0..2M+3 (recursive interval division).
struct Fact {
  Interval inv[2 * TAYLOR_M + 4];
  Fact() {
    RoundUp g;
    inv[0] = Interval(1.0);
    for (int k = 1; k < 2 * TAYLOR_M + 4; ++k) inv[k] = inv[k - 1] / Interval((double)k);
  }
};
inline const Fact& fact() { static const Fact f; return f; }
}  // namespace detail

// cos(y), sin(y) enclosures for an interval y with |y| <= 2 (Taylor + alternating-series remainder).
inline void cos_sin_small(const Interval& y, Interval& c, Interval& s) {
  RoundUp g;
  const auto& F = detail::fact();
  const int M = detail::TAYLOR_M;
  Interval y2 = sqr(y);
  Interval cc(0.0), ss(0.0);
  for (int n = M; n >= 0; --n) {
    Interval ck = F.inv[2 * n], sk = F.inv[2 * n + 1];
    if (n & 1) { ck = -ck; sk = -sk; }
    cc = ck + y2 * cc;
    ss = sk + y2 * ss;
  }
  ss = y * ss;
  double am = abs(y).hi;
  double p = 1.0;  // upper bound of am^(2M+2)
  for (int i = 0; i < 2 * M + 2; ++i) p = p * am;
  double rc = p * F.inv[2 * M + 2].hi;
  double rs = p * am * F.inv[2 * M + 3].hi;
  c = cc + Interval(-rc, rc);
  s = ss + Interval(-rs, rs);
}

// exp(2 pi i t) for interval t in turns; reduce by nearest quarter turn.
inline CInterval cis_turns(const Interval& t) {
  RoundUp g;
  long q = (long)std::floor(4.0 * t.mid() + 0.5);  // any integer is valid
  Interval y = (t - Interval(0.25 * (double)q)) * (pi_iv() * Interval(2.0));
  Interval c, s;
  cos_sin_small(y, c, s);
  switch (((q % 4) + 4) % 4) {
    case 0: return CInterval(c, s);
    case 1: return CInterval(-s, c);
    case 2: return CInterval(-c, -s);
    default: return CInterval(s, -c);
  }
}
// exp(2 pi i k / N)
inline CInterval unit_root(long k, long N) {
  RoundUp g;
  return cis_turns(Interval((double)k) / Interval((double)N));
}

}  // namespace iv
