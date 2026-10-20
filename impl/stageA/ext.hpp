// ext.hpp -- small helpers on top of ../core/interval.hpp for Stage A.
// Everything here must run inside an iv::RoundUp guard (FE_UPWARD) per thread.
// Convention: a "double" that is used as an upper bound carries the suffix U
// or is documented as such; plain doubles fed into interval ops are EXACT.
#pragma once
#include "../core/interval.hpp"
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace sa {
using iv::CInterval;
using iv::Interval;
typedef Interval I;
typedef CInterval Z;

struct Fail : std::runtime_error {
  explicit Fail(const char* m) : std::runtime_error(m) {}
};

// exact-double * interval
inline I mulp(double x, const I& a) {
  if (x >= 0) return I(-((-x) * a.lo), x * a.hi);
  return I(-((-x) * a.hi), x * a.lo);
}
inline Z mulp(double x, const Z& a) { return Z(mulp(x, a.re), mulp(x, a.im)); }
inline I addp(const I& a, double x) { return a + I(x); }
// sqrt of an interval that is known to be mathematically >= 0 (clamps rounding artefacts)
inline I sqrt0(const I& a) {
  if (!(a.hi >= 0)) throw Fail("sqrt0: negative/NaN");
  return iv::sqrt(I(a.lo < 0 ? 0.0 : a.lo, a.hi));
}
inline I zabs(const Z& a) { return sqrt0(iv::abs2(a)); }
inline I imax(const I& a, const I& b) { return I(std::max(a.lo, b.lo), std::max(a.hi, b.hi)); }
inline I imin(const I& a, const I& b) { return I(std::min(a.lo, b.lo), std::min(a.hi, b.hi)); }
inline Z zmulI(const Z& a, const I& s) { return Z(a.re * s, a.im * s); }
inline Z zi_mul(const Z& a) { return Z(-a.im, a.re); }  // i * a
inline bool finite(const I& a) { return std::isfinite(a.lo) && std::isfinite(a.hi); }

// ---- constants (valid after threads hold RoundUp; function statics) ----
inline const I& PI() { static const I v = iv::pi_iv(); return v; }
inline const I& TWOPI() { static const I v = [] { iv::RoundUp g; return iv::pi_iv() * I(2.0); }(); return v; }
inline const I& SQRT6() { static const I v = iv::sqrt_iv(6); return v; }
inline const I& SQRT8() { static const I v = iv::sqrt_iv(8); return v; }
inline const I& SQRT252() { static const I v = iv::sqrt_iv(252); return v; }
inline const I& SIXTH() { static const I v = iv::one_sixth(); return v; }

// fixed-size complex interval matrices
template <int R, int C>
struct CM {
  Z a[R][C];
};
template <int R, int K, int C>
inline CM<R, C> mm(const CM<R, K>& A, const CM<K, C>& B) {
  CM<R, C> O;
  for (int i = 0; i < R; ++i)
    for (int j = 0; j < C; ++j) {
      Z s = A.a[i][0] * B.a[0][j];
      for (int k = 1; k < K; ++k) s = s + A.a[i][k] * B.a[k][j];
      O.a[i][j] = s;
    }
  return O;
}
template <int R, int C>
inline CM<C, R> adj(const CM<R, C>& A) {
  CM<C, R> O;
  for (int i = 0; i < R; ++i)
    for (int j = 0; j < C; ++j) O.a[j][i] = iv::conj(A.a[i][j]);
  return O;
}
template <int R, int C>
inline CM<R, C> madd(const CM<R, C>& A, const CM<R, C>& B) {
  CM<R, C> O;
  for (int i = 0; i < R; ++i)
    for (int j = 0; j < C; ++j) O.a[i][j] = A.a[i][j] + B.a[i][j];
  return O;
}
template <int R, int C>
inline CM<R, C> msub(const CM<R, C>& A, const CM<R, C>& B) {
  CM<R, C> O;
  for (int i = 0; i < R; ++i)
    for (int j = 0; j < C; ++j) O.a[i][j] = A.a[i][j] - B.a[i][j];
  return O;
}
template <int R, int C>
inline CM<R, C> mneg(const CM<R, C>& A) {
  CM<R, C> O;
  for (int i = 0; i < R; ++i)
    for (int j = 0; j < C; ++j) O.a[i][j] = -A.a[i][j];
  return O;
}
template <int R, int C>
inline CM<R, C> mzero() {
  CM<R, C> O;
  for (int i = 0; i < R; ++i)
    for (int j = 0; j < C; ++j) O.a[i][j] = Z(0.0, 0.0);
  return O;
}
// 2x2 inverse (throws if the determinant interval touches 0)
inline CM<2, 2> inv2(const CM<2, 2>& A) {
  Z det = A.a[0][0] * A.a[1][1] - A.a[0][1] * A.a[1][0];
  I d2 = abs2(det);
  if (!(d2.lo > 0)) throw Fail("inv2: singular");
  I r = recip(d2);
  Z id = zmulI(iv::conj(det), r);  // 1/det
  CM<2, 2> O;
  O.a[0][0] = A.a[1][1] * id;
  O.a[0][1] = (-A.a[0][1]) * id;
  O.a[1][0] = (-A.a[1][0]) * id;
  O.a[1][1] = A.a[0][0] * id;
  return O;
}
inline double midp(const I& a) { return 0.5 * a.lo + 0.5 * a.hi; }

}  // namespace sa
