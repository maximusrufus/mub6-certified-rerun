// Tests of the seed cover:
//  A. filters never reject an exact solution (sumok, extend, shared interval);
//  B. exact Hadamard matrices -> seed:arrangement -> every filter passes, the
//     pair is enumerated, and a center of the cover lies within .0320001
//     (the coverage self-check; see COVERAGE.md);
//  C. the center-reuse bounds (box and dual) dominate the true distance;
//  D. optimised enumeration == naive long-double re-implementation at small N.
//
//   test_cover               N=64 in-process cover + all of the above
//   test_cover 112 seeds.txt also checks the coverage of a written N=112 seed file
#include <quadmath.h>
#include <cstring>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include "check.hpp"
#include "hadamard.hpp"
#include "../cover.hpp"
#include "../normal.hpp"
#include "../reuse.hpp"
#include "../seedgen.hpp"

using namespace seedcov;

static const double TWOPI = 2.0 * M_PI;

// ---- exact random lists -----------------------------------------------------
// phases (cycles) p[0..5] with p[0]=0 summing to zero: pick `free` phases,
// solve the remaining two.  Returns false if infeasible.
static bool solve_two(double a_re, double a_im, double &p, double &q) {
  // need e(p)+e(q) = (a_re, a_im) = -z; |a| <= 2
  double m = std::hypot(a_re, a_im);
  if (m > 2.0) return false;
  double th = std::atan2(a_im, a_re) / TWOPI;
  double dlt = std::acos(m / 2.0) / TWOPI;
  p = th + dlt; q = th - dlt;
  return true;
}

// list (0, v0..v4) with v[fixedIdx] given: choose three of the other four at random, solve the remaining two
static bool random_exact_list(std::mt19937_64 &rng, double *v /*5*/, int fixedIdx, double fixedVal) {
  std::uniform_real_distribution<double> U(0.0, 1.0);
  int idx[4], n = 0;
  for (int j = 0; j < 5; ++j) if (j != fixedIdx) idx[n++] = j;
  // random order of the four free indices
  for (int i = 3; i > 0; --i) std::swap(idx[i], idx[rng() % (i + 1)]);
  v[fixedIdx] = fixedVal;
  double re = 1.0, im = 0.0;
  re += std::cos(TWOPI * fixedVal); im += std::sin(TWOPI * fixedVal);
  for (int t = 0; t < 2; ++t) { v[idx[t]] = U(rng); re += std::cos(TWOPI * v[idx[t]]); im += std::sin(TWOPI * v[idx[t]]); }
  double p, q;
  if (!solve_two(-re, -im, p, q)) return false;
  v[idx[2]] = p - std::floor(p); v[idx[3]] = q - std::floor(q);
  return true;
}

// ---- A: filters keep exact solutions ----------------------------------------
static void test_filters_exact(std::mt19937_64 &rng, int N) {
  iv::RoundUp guard;
  Tables T(N);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  long n = 0, nsum = 0;
  for (long it = 0; it < 400000 && n < 100000; ++it) {
    double x[5], y[5];
    // x list: (0,x1..x5); y list: (0, xi, y2..y5) with xi = x[s]
    int s = (int)(rng() % 5);
    if (!random_exact_list(rng, x, s, U(rng))) continue;
    if (!random_exact_list(rng, y, 0, x[s])) continue;
    // y list index 0 is the shared entry; y[1..4] tail
    int bx[5], by[5];
    for (int j = 0; j < 5; ++j) { bx[j] = std::min(N - 1, (int)std::floor(x[j] * N)); by[j] = std::min(N - 1, (int)std::floor(y[j] * N)); }
    // the pipeline shares the bin of xi exactly
    by[0] = bx[s];
    ++n;
    // sumok keeps exact solutions (any ordering of the list)
    CHECKM(sumok_bins(T, bx), "sumok rejected exact X");
    CHECKM(sumok_bins(T, by), "sumok rejected exact Y");
    ++nsum;
    ListData dx, dy; make_ld(T, bx, dx); make_ld(T, by, dy);
    Ivl ix = extend(T, dx, s), iy = extend(T, dy, 0);
    double h = (double)N * (x[s] - T.mid(bx[s]));      // sample is a double: allow 1e-12 for ITS rounding only
    CHECKM(h >= ix.lo - 1e-12 && h <= ix.hi + 1e-12, "true h=%g outside X interval [%g,%g]", h, ix.lo, ix.hi);
    CHECKM(h >= iy.lo - 1e-12 && h <= iy.hi + 1e-12, "true h=%g outside Y interval [%g,%g]", h, iy.lo, iy.hi);
    CHECK(!ivl_disjoint(ix, iy));
  }
  printf("  A: %ld exact lists passed sumok/extend\n", nsum);
  CHECK(nsum > 20000);
}

// ---- C: reuse bounds dominate the true distance --------------------------------
static void test_reuse_bounds(std::mt19937_64 &rng, int N) {
  iv::RoundUp guard;
  Tables T(N);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  long n = 0, nbox = 0, ndual = 0, ndualWins = 0, nsucc = 0;
  for (long it = 0; it < 2000000 && n < 150000; ++it) {
    double xs[5], y[5];
    int s = (int)(rng() % 5);
    if (!random_exact_list(rng, xs, s, U(rng))) continue;
    if (!random_exact_list(rng, y, 0, xs[s])) continue;
    // sort x1..x5 (the pipeline lists x nondecreasing); selected position follows the entry
    double xi = xs[s];
    int order[5] = {0, 1, 2, 3, 4};
    std::sort(order, order + 5, [&](int a, int b) { return xs[a] < xs[b]; });
    double xsort[5]; int sPos = 0;
    for (int m = 0; m < 5; ++m) { xsort[m] = xs[order[m]]; if (order[m] == s) sPos = m; }
    double ytail[4] = {y[1], y[2], y[3], y[4]}; std::sort(ytail, ytail + 4);
    int bx[5], bt[4];
    for (int j = 0; j < 5; ++j) bx[j] = std::min(N - 1, (int)std::floor(xsort[j] * N));
    for (int j = 0; j < 4; ++j) bt[j] = std::min(N - 1, (int)std::floor(ytail[j] * N));
    int sb = bx[sPos];
    int yb[5] = {sb, bt[0], bt[1], bt[2], bt[3]};
    ListData dx, dy; make_ld(T, bx, dx); make_ld(T, yb, dy);
    PairCtx ctx; build_ctx(T, dx, dy, sPos, ctx);
    // nine-order exact point and midpoint
    double q9[9], u9[9]; int b9[9];
    q9[0] = xi; b9[0] = sb; { int p = 1; for (int m = 0; m < 5; ++m) if (m != sPos) { q9[p] = xsort[m]; b9[p] = bx[m]; ++p; } }
    for (int m = 0; m < 4; ++m) { q9[5 + m] = ytail[m]; b9[5 + m] = bt[m]; }
    for (int m = 0; m < 9; ++m) u9[m] = T.mid(b9[m]);
    // random integer displacement k (centre u' = u - k/N)
    int k[9]; int kn2 = 0, kl1 = 0;
    for (int m = 0; m < 9; ++m) { int r = (int)(rng() % 7); k[m] = (r < 3) ? 0 : (r == 3 ? 1 : (r == 4 ? -1 : (r == 5 ? 2 : -2))); }
    if (n % 3 == 0) for (int m = 0; m < 9; ++m) k[m] = (int)(rng() % 3) - 1;   // many small k
    for (int m = 0; m < 9; ++m) { kn2 += k[m] * k[m]; kl1 += std::abs(k[m]); }
    // true N^2 ||q - u'||^2 in long double
    long double d2 = 0;
    for (int m = 0; m < 9; ++m) { long double d = (long double)N * ((long double)q9[m] - (long double)u9[m]) + k[m]; d2 += d * d; }
    // box bound
    long double boxb = (long double)kn2 + 2.25L + kl1;
    CHECKM(d2 <= boxb + 1e-9L, "box bound violated: true %Lg bound %Lg", d2, boxb);
    // dual bound at the heuristic v: CERTIFIED interval bound
    double v[4];
    propose_v(T, ctx, k, v);
    double dsup = dual_cert_hi(T, ctx, k, v);
    long double dualb = (long double)kn2 + 2.25L + dsup;
    CHECKM(d2 <= dualb + 1e-9L, "dual bound violated: true %Lg bound %Lg (kn2=%d)", d2, dualb, kn2);
    // dual bound at random v: 2k.h <= formula for any v
    double vr[4]; for (double &tt : vr) tt = (U(rng) - 0.5) * 4.0;
    double sr = dual_cert_hi(T, ctx, k, vr);
    long double dualr = (long double)kn2 + 2.25L + sr;
    CHECKM(d2 <= dualr + 1e-9L, "dual(random v) violated: true %Lg bound %Lg", d2, dualr);
    ++n;
    if (d2 <= boxb) ++nbox;
    if (dsup < 1e9) ++ndual;
    if (dualb < boxb) ++ndualWins;
    bool succ = box_success(T, kn2, kl1) || dual_success(T, ctx, k, kn2);
    if (succ) {
      ++nsucc;
      // success must imply the true distance is under the threshold
      long double RN = (long double)N * 0.032L;
      CHECKM(d2 < RN * RN, "success but true %Lg >= (RN)^2 %Lg", d2, RN * RN);
    }
  }
  printf("  C: %ld exact points; dual proposals used %ld; dual beats box %ld; reuse success %ld\n", n, ndual, ndualWins, nsucc);
  CHECK(n > 50000);
  CHECK(ndualWins > 100);   // the dual path is exercised
}

// ---- D: interval filters vs a naive __float128 re-implementation, no margins ----
// The naive code uses the SAME exact-double directions (perpendiculars from the tables) so the
// comparison is item-by-item: an interval rejection must be a true violation (up to 1e-12, the
// direction difference for the residual direction), a true violation larger than 1e-12 must be rejected.
typedef __float128 Q;
struct QList { Q cx, cy; Q Jc[5], Js[5]; };   // c' = 2 pi c = (sum sin, -(1+sum cos))
static QList qlist(int N, const int *b) {
  QList L; Q re = 1, im = 0;
  for (int j = 0; j < 5; ++j) {
    Q a = 2 * M_PIq * ((Q)(2 * b[j] + 1) / (Q)(2 * N));
    L.Jc[j] = cosq(a); L.Js[j] = sinq(a); re += L.Jc[j]; im += L.Js[j];
  }
  L.cx = im; L.cy = -re;
  return L;
}
static double naive_violation(const Tables &T, int N, const int *b) {
  QList L = qlist(N, b);
  Q tp = 2 * M_PIq, w = (Q)1 / N, coef = 5 * M_PIq * w * w / 4;
  Q worst = -1e30Q;
  for (int dir = 0; dir < 6; ++dir) {
    Q v0, v1;
    if (dir < 5) { v0 = T.pv0[b[dir]]; v1 = T.pv1[b[dir]]; }
    else { Q n = sqrtq(L.cx * L.cx + L.cy * L.cy); if (n == 0) continue; v0 = L.cx / n; v1 = L.cy / n; }
    Q nv = sqrtq(v0 * v0 + v1 * v1), sup = 0;
    for (int j = 0; j < 5; ++j) sup += fabsq(v0 * L.Jc[j] + v1 * L.Js[j]);
    Q lhs = fabsq(v0 * L.cx + v1 * L.cy), rhs = tp * (coef * nv + w / 2 * sup);
    if (lhs - rhs > worst) worst = lhs - rhs;
  }
  return (double)worst;
}
static void naive_extend(const Tables &T, int N, const int *b, int t, Q &lo, Q &hi) {
  QList L = qlist(N, b);
  Q tp = 2 * M_PIq, E = 5 * M_PIq / (4 * (Q)N);
  lo = -0.5Q; hi = 0.5Q;
  for (int dir = 0; dir < 5; ++dir) {
    Q v0, v1;
    if (dir == 0) { v0 = T.sv0[b[t]]; v1 = T.sv1[b[t]]; }
    else { int j = dir - 1; if (j >= t) ++j; v0 = T.pv0[b[j]]; v1 = T.pv1[b[j]]; }
    Q a = v0 * L.Jc[t] + v1 * L.Js[t];
    if (fabsq(a) < 1e-6Q) continue;
    Q r = E * sqrtq(v0 * v0 + v1 * v1);
    for (int j = 0; j < 5; ++j) if (j != t) r += fabsq(v0 * L.Jc[j] + v1 * L.Js[j]) / 2;
    Q vc = (v0 * L.cx + v1 * L.cy) / tp;
    Q ctr = -(Q)N * vc / a, rad = r / fabsq(a);
    if (ctr - rad > lo) lo = ctr - rad;
    if (ctr + rad < hi) hi = ctr + rad;
  }
}

static void test_naive(int N) {
  iv::RoundUp guard;
  Tables T(N);
  std::vector<XList> X = enumerate_x(T);
  TailTable tt = build_tails(T, X);
  PairCounts pc = count_pairs(T, X, tt);
  long nxN = 0, falseRej = 0, missed = 0, listsChecked = 0, extChecked = 0;
  double widest = 0;
  auto compare_list = [&](const int *b, int t, bool ivKeeps) {
    double viol = naive_violation(T, N, b);
    ++listsChecked;
    if (!ivKeeps && viol < -1e-12) ++falseRej;     // interval rejected something that is not violated
    if (ivKeeps && viol > 1e-12) ++missed;         // interval kept a clearly violated list
    if (ivKeeps) {
      ListData d; make_ld(T, b, d);
      Ivl I = extend(T, d, t);
      Q qlo, qhi; naive_extend(T, N, b, t, qlo, qhi);
      ++extChecked;
      // outward containment: interval bounds enclose the quad value (no margin)
      CHECK((Q)I.lo <= qlo + 1e-26Q || qlo > qhi);
      CHECK((Q)I.hi >= qhi - 1e-26Q || qlo > qhi);
      if (qlo <= qhi) widest = std::max(widest, std::max((double)((Q)qlo - (Q)I.lo), (double)((Q)I.hi - qhi)));
    }
  };
  std::vector<std::array<int, 5>> naiveX;
  for (int a = 0; a < N; ++a) for (int b = a; b < N; ++b) for (int c = b; c < N; ++c) for (int d = c; d < N; ++d) for (int l = d; l < N; ++l) {
    int bx[5] = {a, b, c, d, l};
    if (!x_orient_ok(bx) || !x_gap_ok(T, bx)) continue;
    bool k = sumok_bins(T, bx);
    double viol = naive_violation(T, N, bx);
    compare_list(bx, select_s(N, bx), k);
    if (viol <= 0) ++nxN;
  }
  for (int sb = 0; sb < N; sb += 3)
    for (int t1 = 0; t1 < N; ++t1) for (int t2 = t1; t2 < N; ++t2) for (int t3 = t2; t3 < N; ++t3) for (int t4 = t3; t4 < N; ++t4) {
      int b[5] = {sb, t1, t2, t3, t4};
      compare_list(b, 0, sumok_bins(T, b));
    }
  printf("  D: N=%d  interval counts (x %zu, sumok-only %ld, gap2 %ld, extend %ld); %ld lists vs quad: false rejections %ld, clearly-violated kept %ld, "
         "extend containment checked on %ld lists (widest outward slack %.2e); naive-quad x-lists %ld\n",
         N, X.size(), pc.sumok_only, pc.gap2, pc.extend, listsChecked, falseRej, missed, extChecked, widest, nxN);
  CHECK(falseRej == 0);
  CHECK(missed == 0);
  CHECK((long)X.size() == nxN);
}

// ---- B: Hadamard -> arrangement -> filters + coverage -------------------------
struct CoverIndex {
  int N;
  std::vector<std::array<uint8_t, 9>> c[5];
};

// distance test: double prefilter, then a strict interval confirmation
// (centre = exact rational (2b+1)/(2N); sample point q9 is an exact double; radius^2 lower bound).
static bool covered(const CoverIndex &ci, int s, const double *q9, double radius, double *bestd = nullptr) {
  double r2 = radius * radius * (1 + 1e-12), best = 1e9;
  iv::Interval rr = iv::Interval(320001.0) / iv::Interval(1e7);   // radius .0320001 (the tested radius is this value)
  double r2lo = iv::sqr(rr).lo;
  (void)radius;
  for (auto &b : ci.c[s]) {
    double acc = 0;
    for (int m = 0; m < 9; ++m) { double d = q9[m] - ((double)b[m] + 0.5) / ci.N; acc += d * d; if (acc > r2) break; }
    best = std::min(best, acc);
    if (acc <= r2) {
      iv::Interval S(0.0);
      for (int m = 0; m < 9; ++m) {
        iv::Interval c = iv::Interval((double)(2 * b[m] + 1)) / iv::Interval((double)(2 * ci.N));
        S = S + iv::sqr(iv::Interval(q9[m]) - c);
      }
      if (S.hi <= r2lo) { if (bestd) *bestd = std::sqrt(acc); return true; }
    }
  }
  if (bestd) *bestd = std::sqrt(best);
  return false;
}

// tails are stored in lexicographic (t1,t2,t3,t4) order, t1 in the lowest byte
static bool tail_less(uint32_t a, uint32_t b) {
  for (int sh = 0; sh < 32; sh += 8) { uint32_t x = (a >> sh) & 255, y = (b >> sh) & 255; if (x != y) return x < y; }
  return false;
}
static bool in_sorted_tails(const std::vector<uint32_t> &v, uint32_t key) { return std::binary_search(v.begin(), v.end(), key, tail_less); }

static void check_sample(const double *P, int N, const Tables &T, const std::vector<XList> &X, const TailTable &tails,
                         const CoverIndex &ci, int pc, int gt, bool bc, long &nCovered, double &worst, const char *tag) {
  nf::NormalForm nfm;
  if (!nf::arrange(P, N, pc, gt, bc, nfm)) { CHECKM(false, "%s arrange failed", tag); return; }
  int bx[5]; for (int m = 0; m < 5; ++m) bx[m] = nfm.bx[m];
  int s = nfm.s, sb = bx[s];
  // F1, F2, F3
  CHECKM(x_orient_ok(bx), "%s F1 rejected an exact Hadamard", tag);
  CHECKM(x_gap_ok(T, bx), "%s F2 rejected an exact Hadamard", tag);
  CHECKM(sumok_bins(T, bx), "%s F3 rejected an exact Hadamard", tag);
  int tail[4] = {nfm.by[2], nfm.by[3], nfm.by[4], nfm.by[5]};
  int yb[5] = {sb, tail[0], tail[1], tail[2], tail[3]};
  CHECKM(sumok_bins(T, yb), "%s Y sumok rejected an exact Hadamard", tag);
  CHECKM(y_gap_ok(T, bx[4], sb, tail), "%s gap2 rejected an exact Hadamard", tag);
  ListData dx, dy; make_ld(T, bx, dx); make_ld(T, yb, dy);
  Ivl ix = extend(T, dx, s), iy = extend(T, dy, 0);
  double h = (double)N * (nfm.x[s + 1] - T.mid(sb));   // sample is a double: 1e-12 for its own rounding only
  CHECKM(h >= ix.lo - 1e-12 && h <= ix.hi + 1e-12 && h >= iy.lo - 1e-12 && h <= iy.hi + 1e-12 && !ivl_disjoint(ix, iy), "%s shared interval rejected an exact Hadamard", tag);
  // the pair is in the enumeration
  XList key; for (int m = 0; m < 5; ++m) key.bx[m] = (uint8_t)bx[m];
  bool foundX = std::binary_search(X.begin(), X.end(), key, [](const XList &a, const XList &b) { return std::lexicographical_compare(a.bx, a.bx + 5, b.bx, b.bx + 5); });
  CHECKM(foundX, "%s X-list not enumerated", tag);
  uint32_t tk = (uint32_t)tail[0] | ((uint32_t)tail[1] << 8) | ((uint32_t)tail[2] << 16) | ((uint32_t)tail[3] << 24);
  CHECKM(in_sorted_tails(tails.tail[sb], tk), "%s Y tail not enumerated", tag);
  // coverage: centre within .0320001 (cycles) with the same s
  double dist = 0;
  bool cv = covered(ci, s, nfm.q9, 0.0320001, &dist);
  CHECKM(cv, "%s NOT COVERED (s=%d, nearest centre %g)", tag, s, dist);
  if (cv) { ++nCovered; worst = std::max(worst, dist); }
}

static CoverIndex index_from_result(int N, const std::vector<XList> &X, const CoverResult &R) {
  CoverIndex ci; ci.N = N;
  for (const Center &c : R.centers) {
    int b[9]; center_bins(X[c.xi], c.tail, b);
    std::array<uint8_t, 9> a; for (int m = 0; m < 9; ++m) a[m] = (uint8_t)b[m];
    ci.c[X[c.xi].s].push_back(a);
  }
  return ci;
}
static CoverIndex index_from_file(int N, const char *fn) {
  CoverIndex ci; ci.N = N;
  std::ifstream f(fn);
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream is(line);
    int s; is >> s; std::array<uint8_t, 9> a; int v;
    for (int m = 0; m < 9; ++m) { is >> v; a[m] = (uint8_t)v; }
    ci.c[s].push_back(a);
  }
  return ci;
}

static void test_coverage(std::mt19937_64 &rng, int N, const CoverIndex &ci, const std::vector<XList> &X, const TailTable &tails, int nLM) {
  iv::RoundUp guard;
  Tables T(N);
  long nCov = 0; double worst = 0;
  // F6 family on a grid, including ties and boundaries
  for (int ia = 0; ia < 16; ++ia) for (int ib = 0; ib < 16; ++ib) {
    double P[36]; hd::fourier_family(ia / 16.0 + 0.0031 * ib, ib / 16.0 + 0.0017 * ia, P);
    for (int pc = 0; pc < 3; ++pc) for (int gt = 0; gt < 2; ++gt) check_sample(P, N, T, X, tails, ci, pc, gt, (pc + gt) & 1, nCov, worst, "F6");
  }
  { double P[36]; hd::fourier_family(0, 0, P);
    for (int pc = 0; pc < 30; ++pc) for (int gt = 0; gt < 3; ++gt) for (int bc = 0; bc < 2; ++bc) check_sample(P, N, T, X, tails, ci, pc, gt, bc, nCov, worst, "F6exact"); }
  int got = 0, tries = 0;
  while (got < nLM && tries < 20 * nLM) {
    ++tries;
    double P[36];
    if (!hd::random_hadamard(rng, P)) continue;
    ++got;
    for (int pc = 0; pc < 2; ++pc) for (int bc = 0; bc < 2; ++bc) check_sample(P, N, T, X, tails, ci, pc, pc, bc, nCov, worst, "LM");
  }
  printf("  B: N=%d  %ld normal-form samples covered (largest distance to the first covering centre: %.5f <= .0320001); %d LM Hadamards\n", N, nCov, worst, got);
  CHECK(nCov > 1000);
}

int main(int argc, char **argv) {
  std::mt19937_64 rng(20261019);
  iv::RoundUp mainGuard;
  printf("test_cover\n");
  if (argc >= 3) {   // file mode: coverage check of a written seed file only
    int N = atoi(argv[1]);
    Tables T(N);
    std::vector<XList> X = enumerate_x(T);
    TailTable tails = build_tails(T, X);
    CoverIndex ci = index_from_file(N, argv[2]);
    size_t tot = 0;
    for (int s = 0; s < 5; ++s) tot += ci.c[s].size();
    printf("  file %s: %zu centres (N=%d)\n", argv[2], tot, N);
    test_coverage(rng, N, ci, X, tails, 3000);
    return finish("test_cover(file)");
  }
  test_filters_exact(rng, 112);
  test_filters_exact(rng, 64);
  test_reuse_bounds(rng, 112);
  test_reuse_bounds(rng, 64);
  test_naive(20);
  test_naive(24);

  // in-process N=48 cover without reuse, then coverage of Hadamard normal forms
  {
    int N = 48;
    Tables T(N);
    std::vector<XList> X = enumerate_x(T);
    TailTable tails = build_tails(T, X);
    Options o; o.N = N;
    o.cut = 0;   // no reuse: every pair becomes a centre (radius 3/(2N) <= .032 for N >= 47); reuse is covered by C and by the N=112 file check
    CoverResult R = run_cover(T, X, tails, o);
    printf("  in-process N=%d cover: %zu x-lists, %ld pairs, %zu centres\n", N, X.size(), R.pairs, R.centers.size());
    CoverIndex ci = index_from_result(N, X, R);
    test_coverage(rng, N, ci, X, tails, 400);
  }
  return finish("test_cover");
}
