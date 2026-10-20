// Tests of the certificate arithmetic (arith.hpp over ../core/interval.hpp): every enclosure
// must CONTAIN the libquadmath value; widths and lower/upper bound helpers are checked.
#include <quadmath.h>
#include <random>
#include "check.hpp"
#include "../cover.hpp"

using namespace seedcov;

int main() {
  iv::RoundUp guard;
  std::mt19937_64 rng(7);
  const __float128 PI128 = M_PIq;
  for (int N : {112, 64, 48}) {
    Tables T(N);
    double widest = 0, wPP = 0, wc = 0;
    for (int b = 0; b < N; ++b) {
      __float128 a = 2 * PI128 * ((__float128)(2 * b + 1) / (__float128)(2 * N));
      __float128 cc = cosq(a), ss = sinq(a);
      CHECK((__float128)T.ec[b].lo <= cc && cc <= (__float128)T.ec[b].hi);
      CHECK((__float128)T.es[b].lo <= ss && ss <= (__float128)T.es[b].hi);
      widest = std::max(widest, std::max(T.ec[b].width(), T.es[b].width()));
      // PP / SP enclosures of exact-double directions against all J_b'
      for (int c = 0; c < N; ++c) {
        __float128 a2 = 2 * PI128 * ((__float128)(2 * c + 1) / (__float128)(2 * N));
        __float128 pp = (__float128)T.pv0[b] * cosq(a2) + (__float128)T.pv1[b] * sinq(a2);
        __float128 sp = (__float128)T.sv0[b] * cosq(a2) + (__float128)T.sv1[b] * sinq(a2);
        const iv::Interval &P = T.PP[(size_t)b * N + c], &S = T.SP[(size_t)b * N + c];
        CHECK((__float128)P.lo <= pp && pp <= (__float128)P.hi);
        CHECK((__float128)S.lo <= sp && sp <= (__float128)S.hi);
        CHECK((__float128)T.PPa[(size_t)b * N + c] >= fabsq(pp));
        CHECK((__float128)T.SPa[(size_t)b * N + c] >= fabsq(sp));
        wPP = std::max(wPP, P.width());
      }
      // direction norms are upper bounds
      __float128 n = sqrtq((__float128)T.pv0[b] * T.pv0[b] + (__float128)T.pv1[b] * T.pv1[b]);
      CHECK((__float128)T.pnv[b] >= n);
    }
    // list data c encloses the quad value
    std::uniform_int_distribution<int> B(0, N - 1);
    for (int it = 0; it < 100000; ++it) {
      int b[5]; for (int j = 0; j < 5; ++j) b[j] = B(rng);
      ListData d; make_ld(T, b, d);
      __float128 re = 1, im = 0;
      for (int j = 0; j < 5; ++j) {
        __float128 a = 2 * PI128 * ((__float128)(2 * b[j] + 1) / (__float128)(2 * N));
        re += cosq(a); im += sinq(a);
      }
      __float128 c0 = im / (2 * PI128), c1 = -re / (2 * PI128);
      CHECK((__float128)d.c0.lo <= c0 && c0 <= (__float128)d.c0.hi);
      CHECK((__float128)d.c1.lo <= c1 && c1 <= (__float128)d.c1.hi);
      wc = std::max(wc, std::max(d.c0.width(), d.c1.width()));
    }
    // constants: RN2lo <= (0.032 N)^2 and tight; Eg encloses 5pi/(4N); pi enclosure
    __float128 RN = (__float128)0.032Q * N, ideal = RN * RN;
    CHECK((__float128)T.RN2lo <= ideal);
    CHECK(ideal - (__float128)T.RN2lo < 1e-12Q);
    CHECK((__float128)T.Eg.lo <= 5 * PI128 / (4 * (__float128)N) && 5 * PI128 / (4 * (__float128)N) <= (__float128)T.Eg.hi);
    CHECK((__float128)iv::pi_iv().lo <= PI128 && PI128 <= (__float128)iv::pi_iv().hi);
    printf("  N=%d  widest e(u_b) enclosure %.3g, widest PP entry %.3g, widest c %.3g, RN2lo %.12f\n", N, widest, wPP, wc, T.RN2lo);
    CHECK(widest < 1e-14);
    CHECK(wc < 1e-13);
  }
  // helper semantics
  {
    iv::Interval a(-2.0, 3.0), b(1.0, 2.0);
    CHECK(abs_lo(a) == 0.0 && abs_hi(a) == 3.0);
    CHECK(abs_lo(b) == 1.0 && abs_hi(b) == 2.0);
    CHECK(abs_lo(iv::Interval(-5.0, -4.0)) == 4.0);
    iv::Interval m = mulp(iv::Interval(1.0, 2.0), -3.0);
    CHECK(m.lo <= -6.0 && m.hi >= -3.0);
  }
  return finish("test_arith");
}
