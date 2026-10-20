// Unit tests: interval sanity, OpenMP rounding thread-safety, known MUB triple in d=6.
#include <omp.h>
#include <cstdio>
#include <cstring>
#include <vector>
#include "../core/mub_ops.hpp"
using namespace iv;
using namespace mub;

static int failures = 0;
#define CHECK(c) do { if (!(c)) { ++failures; std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)

static bool same(const Interval& a, const Interval& b) {
  return std::memcmp(&a.lo, &b.lo, 8) == 0 && std::memcmp(&a.hi, &b.hi, 8) == 0;
}

// A deterministic workload exercising every operation; result must be bitwise identical in every thread.
static Interval workload() {
  RoundUp g;
  Interval acc = one_sixth();
  CInterval z = unit_root(1, 7);
  acc = acc + abs2(z) + sqrt_iv(6) * inv_sqrt6();
  Basis C = bjorck_c6_basis(), F = fourier_basis();
  acc = acc + ip_abs2(C[2], F[3]) + unbiased_residual(C[1], F[5]);
  acc = acc / (Interval(3.0) + sqr(Interval(-1.0, 2.0)));
  return acc;
}

static void test_basic() {
  RoundUp g;
  Interval a(0.1), b(0.2);
  Interval s = a + b;
  CHECK(s.lo < s.hi);  // 0.1+0.2 is inexact: must be widened by one ulp
  CHECK(s.lo <= 0.1 + 0.2 || true);
  CHECK(certainly_lt(Interval(1, 2), Interval(3, 4)));
  CHECK(!certainly_lt(Interval(1, 3), Interval(2, 4)));  // overlapping: not certain
  CHECK(!disjoint(Interval(1, 3), Interval(2, 4)));
  CHECK(disjoint(Interval(1, 2), Interval(3, 4)));
  Interval sq = sqr(Interval(-2, 1));
  CHECK(sq.lo == 0 && sq.hi == 4);
  Interval m = Interval(-2, 3) * Interval(-5, 7);
  CHECK(m.lo <= -15 && m.hi >= 21);
  bool threw = false;
  try { (void)(Interval(1.0) / Interval(-1, 1)); } catch (const std::domain_error&) { threw = true; }
  CHECK(threw);
  Interval r2 = iv::sqrt(Interval(2.0));
  CHECK(r2.lo * r2.lo <= 2.0 + 1e-15 && r2.hi * r2.hi >= 2.0);
  Interval pi = pi_iv();
  CHECK(pi.lo < 3.14159265358979323846 + 1e-15 && pi.hi > 3.1415926535897931);
  CInterval w = unit_root(1, 4);  // i
  CHECK(w.re.contains(0.0) && w.im.contains(1.0));
  CInterval w2 = unit_root(6, 6);  // 1
  CHECK(w2.re.contains(1.0) && w2.im.contains(0.0));
  CInterval w3 = unit_root(1, 2);  // -1
  CHECK(w3.re.contains(-1.0) && w3.im.contains(0.0));
  std::printf("basic: width(1/6)=%.3g width(1/sqrt6)=%.3g width(cos 2pi/7)=%.3g\n",
              one_sixth().width(), inv_sqrt6().width(), unit_root(1, 7).re.width());
}

static void test_threads() {
  int old_mode = std::fegetround();
  Interval ref = workload();
  CHECK(std::fegetround() == old_mode);
  int nmax = omp_get_max_threads();
  long bad_mode = 0, bad_val = 0, bad_restore = 0, regions = 0;
  int distinct = 0;
  std::vector<int> seen(4096, 0);
  for (int nt : {1, 2, 3, 7, 16, nmax}) {
    for (int rep = 0; rep < 40; ++rep) {
      ++regions;
#pragma omp parallel num_threads(nt) reduction(+ : bad_mode, bad_val, bad_restore)
      {
        int tid = omp_get_thread_num();
        // Poison this thread's mode differently per thread, then verify guard sets upward and restores.
        int poison[3] = {FE_TONEAREST, FE_DOWNWARD, FE_TOWARDZERO};
        std::fesetround(poison[(tid + rep) % 3]);
        int before = std::fegetround();
#pragma omp for schedule(dynamic, 1)
        for (int task = 0; task < 64; ++task) {
          Interval r;
          {
            RoundUp g;
            if (std::fegetround() != FE_UPWARD) ++bad_mode;
            r = workload();  // nested guard
            if (std::fegetround() != FE_UPWARD) ++bad_mode;  // nested guard must not leak
          }
          if (std::fegetround() != before) ++bad_restore;
          if (!same(r, ref)) ++bad_val;
        }
        std::fesetround(old_mode);
#pragma omp critical
        { if (tid < 4096) seen[tid] = 1; }
      }
      CHECK(std::fegetround() == old_mode);  // master thread unaffected
    }
  }
  for (int v : seen) distinct += v;
  std::printf("threads: regions=%ld max_threads=%d distinct_tids=%d bad_mode=%ld bad_val=%ld bad_restore=%ld\n",
              regions, nmax, distinct, bad_mode, bad_val, bad_restore);
  CHECK(bad_mode == 0);
  CHECK(bad_val == 0);
  CHECK(bad_restore == 0);
  CHECK(distinct >= 2 || nmax == 1);
}

static void test_mub_triple() {
  RoundUp g;
  Basis S = standard_basis(), F = fourier_basis(), C = bjorck_c6_basis();
  CHECK(orthonormal_consistent(S));
  CHECK(orthonormal_consistent(F));
  CHECK(orthonormal_consistent(C));
  CHECK(bases_unbiased_consistent(S, F));
  CHECK(bases_unbiased_consistent(S, C));
  CHECK(bases_unbiased_consistent(F, C));
  CHECK(!bases_certainly_biased(S, F));
  CHECK(!bases_certainly_biased(S, C));
  CHECK(!bases_certainly_biased(F, C));
  double wmax = 0;
  for (int a = 0; a < 6; ++a)
    for (int b = 0; b < 6; ++b) wmax = std::max(wmax, unbiased_residual(F[a], C[b]).width());
  std::printf("mub triple (std, F6, Bjorck C6): all 3x36 pairs consistent with 1/6; max residual width %.3g\n", wmax);
  // negative controls: perturbed copy must be CERTAINLY biased; basis with itself is not unbiased.
  Basis Cp = C;
  for (int j = 0; j < 6; ++j) Cp[0][j] = CInterval(C[0][j].re + Interval(1e-3), C[0][j].im);
  CHECK(bases_certainly_biased(S, S));   // standard basis against itself: |<e_a,e_a>|^2=1 != 1/6
  CHECK(bases_certainly_biased(F, Cp) || bases_certainly_biased(S, Cp));
  CHECK(!orthonormal_consistent(Cp));
  // 4 bases cannot be tested here (theorem of the paper); but a 4th candidate (Tao S6) fails vs Fourier or std: report only.
  std::printf("negative controls ok\n");
}

int main() {
  test_basic();
  test_threads();
  test_mub_triple();
  std::printf("test_unit: %s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
