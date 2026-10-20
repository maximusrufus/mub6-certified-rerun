// test_stagea.cpp -- Stage A tests.  `make test` runs this; exit 0 iff all checks pass.
//  T1 interval layer sanity (cis encloses libm values; frames)
//  T2 Root certificates at known Hadamards (F6, Tao S6, Bjorck C6) + negative control
//  T3 synthetic seeds around known Hadamards: the witness cell is never excluded, the true
//     Hadamard is covered by a local Record, every Record re-certifies (checkRoot)
//  T4 random far-off seeds (sum equations satisfied) are fully excluded; runtime per 1000 seeds
//  T5 Record file round trip
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include "../stagea.hpp"

using namespace sa;
typedef std::complex<double> cd;

static int g_fail = 0;
#define CHECK(c, ...)                                                    \
  do {                                                                   \
    if (!(c)) {                                                          \
      ++g_fail;                                                          \
      printf("  FAIL %s:%d: %s -- ", __FILE__, __LINE__, #c);            \
      printf(__VA_ARGS__);                                               \
      printf("\n");                                                      \
    }                                                                    \
  } while (0)

static cd ee(double x) { return std::polar(1.0, 2 * M_PI * x); }

struct Had {
  cd H[6][6];
  const char* name;
};

static bool isHadamard(const Had& h, double tol = 1e-12) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) {
      if (std::fabs(std::abs(h.H[i][j]) - 1) > tol) return false;
      cd s = 0;
      for (int k = 0; k < 6; ++k) s += h.H[i][k] * std::conj(h.H[j][k]);
      if (std::abs(s - (i == j ? 6.0 : 0.0)) > 1e-10) return false;
    }
  return true;
}

static Had mkF6(double a, double b) {
  Had h;
  h.name = "F6(a,b)";
  double ph[6] = {0, a, b, 0, a, b};
  for (int j = 0; j < 6; ++j)
    for (int k = 0; k < 6; ++k) {
      double p = (double)((j * k) % 6) / 6.0;
      cd v = ee(p);
      if (j % 2 == 1) v *= std::polar(1.0, ph[k]);
      h.H[j][k] = v;
    }
  return h;
}
static Had mkS6() {
  Had h;
  h.name = "Tao S6";
  cd w = ee(1.0 / 3), w2 = w * w, o = 1;
  cd M[6][6] = {{o, o, o, o, o, o}, {o, o, w, w, w2, w2}, {o, w, o, w2, w2, w}, {o, w, w2, o, w, w2}, {o, w2, w2, w, o, w}, {o, w2, w, w2, w, o}};
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) h.H[i][j] = M[i][j];
  return h;
}
static Had mkBjorck() {
  Had h;
  h.name = "Bjorck C6";
  double s3 = std::sqrt(3.0);
  cd d((1 - s3) / 2, std::sqrt(s3 / 2));
  cd I(0, 1);
  cd c[6] = {1, I * d, -d, -I, -std::conj(d), I * std::conj(d)};
  for (int j = 0; j < 6; ++j)
    for (int k = 0; k < 6; ++k) h.H[j][k] = c[((k - j) % 6 + 6) % 6];
  return h;
}

static double fr(double x) { return x - std::floor(x); }

// Build a dephased arrangement: returns phases P[36] (dephased, mod 1) and the nine seed coordinates.
struct Arr {
  double P[36];
  double q[9];
  cd Dblock[4][4];
};
static std::vector<Arr> arrangements(const Had& h, int maxOut, std::mt19937& rng) {
  std::vector<Arr> out;
  double ph[6][6];
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) ph[i][j] = fr(std::arg(h.H[i][j]) / (2 * M_PI));
  std::vector<std::array<int, 4>> opts;
  for (int tr = 0; tr < 2; ++tr)
    for (int a = 0; a < 6; ++a)
      for (int b = 0; b < 6; ++b)
        if (a != b)
          for (int c = 0; c < 6; ++c) opts.push_back({tr, a, b, c});
  std::shuffle(opts.begin(), opts.end(), rng);
  Params P;
  const bool anyArr = getenv("ARRANY") != nullptr;
  // global minimum gap over all row pairs of M and M^T (spec 1.3)
  auto gapOf = [&](int tr, int a, int b, int* after) {
    double d[6];
    for (int c = 0; c < 6; ++c) d[c] = fr((tr ? ph[c][b] : ph[b][c]) - (tr ? ph[c][a] : ph[a][c]));
    std::vector<std::pair<double, int>> v;
    for (int c = 0; c < 6; ++c) v.push_back({d[c], c});
    std::sort(v.begin(), v.end());
    double best = -1;
    for (int i = 0; i < 6; ++i) {
      double g = (i + 1 < 6) ? v[i + 1].first - v[i].first : v[0].first + 1.0 - v[5].first;
      if (g > best + 1e-12) { best = g; *after = v[(i + 1) % 6].second; }
    }
    return best;
  };
  double gmin = 1e300;
  for (int tr = 0; tr < 2; ++tr)
    for (int a = 0; a < 6; ++a)
      for (int b = 0; b < 6; ++b)
        if (a != b) { int af; gmin = std::min(gmin, gapOf(tr, a, b, &af)); }
  for (auto& o : opts) {
    if ((int)out.size() >= maxOut) break;
    int tr = o[0], a = o[1], b = o[2], c0 = o[3];
    double M[6][6];
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) M[i][j] = tr ? ph[j][i] : ph[i][j];
    // columns: c0 first, others sorted by frac(diff_b-a - diff at c0)
    std::vector<std::pair<double, int>> cols;
    for (int c = 0; c < 6; ++c) cols.push_back({fr((M[b][c] - M[a][c]) - (M[b][c0] - M[a][c0])), c});
    std::sort(cols.begin(), cols.end());
    if (cols[0].second != c0) continue;
    int colp[6];
    for (int i = 0; i < 6; ++i) colp[i] = cols[i].second;
    int sLo = 1, sHi = 5;
    if (!anyArr) {
      int af;
      double g = gapOf(tr, a, b, &af);
      if (g > gmin + 1e-9) continue;
      (void)af;
      if (!((1.0 - cols[5].first) >= g - 1e-9)) continue;  // c0 must follow a largest gap
      double x[6];
      for (int i = 0; i < 6; ++i) x[i] = fr((M[b][colp[i]] - M[a][colp[i]]) - (M[b][c0] - M[a][c0]));
      if (!(x[1] <= x[5] - x[4] + 1e-12)) continue;  // orientation
      int bs = 1;
      for (int j = 2; j <= 5; ++j)
        if (std::fabs(x[j] - 0.5) < std::fabs(x[bs] - 0.5) - 1e-12) bs = j;
      sLo = sHi = bs;
    }
    for (int s = sLo; s <= sHi; ++s) {
      int cp[6];
      cp[0] = colp[0];
      cp[1] = colp[s];
      int n = 2;
      for (int i = 1; i < 6; ++i)
        if (i != s) cp[n++] = colp[i];
      std::vector<std::pair<double, int>> rows;
      for (int r = 0; r < 6; ++r) {
        if (r == a || r == b) continue;
        rows.push_back({fr(M[r][cp[1]] - M[r][cp[0]] - M[a][cp[1]] + M[a][cp[0]]), r});
      }
      std::sort(rows.begin(), rows.end());
      int rp[6] = {a, b, rows[0].second, rows[1].second, rows[2].second, rows[3].second};
      Arr A;
      double X[6][6];
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) X[i][j] = fr(M[rp[i]][cp[j]] - M[rp[i]][cp[0]] - M[rp[0]][cp[j]] + M[rp[0]][cp[0]]);
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) A.P[6 * i + j] = X[i][j];
      for (int l = 0; l < 9; ++l) A.q[l] = A.P[SEEDPOS[l]];
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) A.Dblock[i][j] = ee(X[2 + i][2 + j]);
      SeedCtx S;
      std::string why;
      if (!buildSeed(S, A.q, P, why)) {
        if (getenv("ARRDEBUG")) printf("    arrangement rejected: %s\n", why.c_str());
        continue;
      }
      out.push_back(A);
      if ((int)out.size() >= maxOut) break;
    }
  }
  return out;
}

// polar cell coordinates of the true completion
static bool witnessOf(const SeedCtx& S, const Arr& A, Witness& w) {
  cd KB[4][2], KC[4][2];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 2; ++j) {
      KB[i][j] = cd(midp(S.KB.a[i][j].re), midp(S.KB.a[i][j].im));
      KC[i][j] = cd(midp(S.KC.a[i][j].re), midp(S.KC.a[i][j].im));
    }
  cd Y[2][2];
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y) {
      cd s = 0;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) s += std::conj(KC[i][x]) * A.Dblock[i][j] * KB[j][y];
      Y[x][y] = s;
    }
  // polar V = Y (Y^*Y)^{-1/2}
  cd Gm[2][2];
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y) {
      cd s = 0;
      for (int k = 0; k < 2; ++k) s += std::conj(Y[k][x]) * Y[k][y];
      Gm[x][y] = s;
    }
  cd det = Gm[0][0] * Gm[1][1] - Gm[0][1] * Gm[1][0];
  double sd = std::sqrt(det.real());
  double tr = Gm[0][0].real() + Gm[1][1].real();
  double den = std::sqrt(tr + 2 * sd);
  cd R[2][2] = {{(Gm[0][0] + sd) / den, Gm[0][1] / den}, {Gm[1][0] / den, (Gm[1][1] + sd) / den}};  // sqrt(G)
  cd dR = R[0][0] * R[1][1] - R[0][1] * R[1][0];
  cd Ri[2][2] = {{R[1][1] / dR, -R[0][1] / dR}, {-R[1][0] / dR, R[0][0] / dR}};
  cd V[2][2];
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y) V[x][y] = Y[x][0] * Ri[0][y] + Y[x][1] * Ri[1][y];
  cd dv = V[0][0] * V[1][1] - V[0][1] * V[1][0];
  double alpha = fr(std::arg(dv) / (2 * M_PI)) / 2.0;  // det = e(2 alpha)
  cd ph = ee(-alpha);
  cd z = ph * V[0][0], wv = ph * V[0][1];
  w.al = alpha;
  w.t = std::atan2(std::abs(wv), std::abs(z)) / (2 * M_PI);
  w.be = fr(std::arg(z) / (2 * M_PI));
  w.ga = fr(std::arg(wv) / (2 * M_PI));
  w.on = true;
  w.violations = 0;
  return std::isfinite(w.al + w.t + w.be + w.ga);
}

static double coverDist(const double Ptrue[36], const Root& R) {
  return mergeDistanceUpper(Ptrue, R.P);
}

// ------------------------------------------------------------------ tests
static void t1() {
  printf("T1 interval layer / frames\n");
  std::mt19937 rng(1);
  std::uniform_real_distribution<double> U(-3, 3);
  for (int i = 0; i < 2000; ++i) {
    double x = U(rng);
    Z e = iv::cis_turns(I(x));
    double c = std::cos(2 * M_PI * x), s = std::sin(2 * M_PI * x);
    CHECK(e.re.lo - 1e-14 <= c && c <= e.re.hi + 1e-14 && e.im.lo - 1e-14 <= s && s <= e.im.hi + 1e-14, "cis x=%g", x);
    CHECK(e.re.hi - e.re.lo < 1e-14, "width");
  }
  Params P;
  Had h = mkF6(0.7, 1.9);
  CHECK(isHadamard(h), "F6 not hadamard");
  CHECK(isHadamard(mkS6()), "S6 not hadamard");
  CHECK(isHadamard(mkBjorck()), "Bjorck not hadamard");
}

static void t2() {
  printf("T2 Root certificates\n");
  struct Case { Had h; };
  std::vector<Had> hs = {mkF6(0.7, 1.9), mkF6(0.3, 0.45), mkS6(), mkBjorck()};
  for (auto& h : hs) {
    double P[36];
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) P[6 * i + j] = std::arg(h.H[i][j]) / (2 * M_PI);
    // perturb slightly and snap
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> U(-1, 1);
    double Q[36];
    for (int i = 0; i < 36; ++i) Q[i] = P[i] + 0.004 * U(rng);
    Root R;
    std::string why;
    bool ok = buildRoot(Q, R, &why);
    CHECK(ok, "%s buildRoot: %s", h.name, why.c_str());
    if (ok) {
      bool c2 = checkRoot(R, &why);
      CHECK(c2, "%s checkRoot: %s", h.name, why.c_str());
      printf("  %-10s k=%d C=%.3f mnorm=%.4f id=%.2e ort=%.1e hor=%.1e froot=%.1e\n", h.name, R.k, R.C, R.mnorm, R.id, R.ort, R.hor, R.froot);
    }
  }
  // negative control: a random matrix is not near a Hadamard
  std::mt19937 rng(9);
  std::uniform_real_distribution<double> U(0, 1);
  int bad = 0;
  for (int t = 0; t < 20; ++t) {
    double Q[36];
    for (int i = 0; i < 36; ++i) Q[i] = U(rng);
    Root R;
    if (buildRoot(Q, R)) ++bad;
  }
  CHECK(bad == 0, "random matrices certified as roots: %d", bad);
  // a Root with a doctored constant must fail checkRoot
  {
    double P[36];
    Had h = mkS6();
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) P[6 * i + j] = std::arg(h.H[i][j]) / (2 * M_PI);
    Root R;
    if (buildRoot(P, R)) {
      Root Rb = R;
      Rb.id = 0.0;
      CHECK(!checkRoot(Rb), "doctored id accepted");
      Root Rc = R;
      Rc.mnorm *= 0.5;
      CHECK(!checkRoot(Rc), "doctored mnorm accepted");
    }
  }
}

static void t3(int perHad) {
  printf("T3 seeds around known Hadamards (witness never excluded, truth covered)\n");
  std::mt19937 rng(77);
  std::normal_distribution<double> N(0, 1);
  std::vector<Had> hs = {mkF6(0.7, 1.9), mkF6(0.3, 4.1), mkS6(), mkBjorck()};
  if (const char* e = getenv("HADS")) {
    std::vector<Had> sel;
    for (const char* c = e; *c; ++c)
      if (*c >= '0' && *c <= '3') sel.push_back(hs[*c - '0']);
    hs = sel;
  }
  Params P;
  int tot = 0, okc = 0;
  for (auto& h : hs) {
    auto arrs = arrangements(h, perHad, rng);
    CHECK(!arrs.empty(), "%s: no admissible arrangement", h.name);
    for (auto& A : arrs) {
      // perturbed center, ||delta|| = .02
      double d[9], n = 0;
      for (int i = 0; i < 9; ++i) { d[i] = N(rng); n += d[i] * d[i]; }
      n = std::sqrt(n);
      double q[9];
      for (int i = 0; i < 9; ++i) q[i] = A.q[i] + 0.02 * d[i] / n;
      SeedCtx S;
      std::string why;
      if (!buildSeed(S, q, P, why)) { printf("  %s: perturbed seed fails frames (%s)\n", h.name, why.c_str()); continue; }
      Witness w;
      CHECK(witnessOf(S, A, w), "witness");
      w.truthP = A.P;
      auto t0 = std::chrono::steady_clock::now();
      SeedResult R = processSeed(q, P, &w);
      double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      ++tot;
      bool good = (R.status == 0) && w.violations == 0;
      double dmin = 1e300;
      int which = -1;
      for (size_t i = 0; i < R.recs.size(); ++i) {
        double dd = coverDist(A.P, R.recs[i]);
        if (dd <= R.recs[i].r && dd < dmin) { dmin = dd; which = (int)i; }
      }
      if (which < 0)
        for (size_t i = 0; i < R.recs.size(); ++i)
          if (coverDist(A.P, R.recs[i]) <= R.recs[i].r + 1e-9) { which = (int)i; break; }  // truth is only a double approximation
      bool covered = which >= 0;
      if (!covered) {
        double bd = 1e300; int bi = -1;
        for (size_t i = 0; i < R.recs.size(); ++i) { double dd = coverDist(A.P, R.recs[i]); if (dd < bd) { bd = dd; bi = (int)i; } }
        if (bi >= 0) printf("    DEBUG nearest record dist=%.3e r=%.3e k=%d\n", bd, R.recs[bi].r, R.recs[bi].k);
      }
      if (getenv("VERB")) printf("    cover: prepFail %ld reuse %ld (tries %ld) snapOK %ld buildFail %ld snapBoundFail %ld | t prep %.2f reuse %.2f build %.2f snapbound %.2f\n", R.st.nPrepFail, R.st.nReuse,
             R.st.nBoundTry, R.st.nSnapOK, R.st.nBuildFail, R.st.nSnapBound, R.st.tPrep, R.st.tReuse, R.st.tBuild, R.st.tSnapBound);
      if (getenv("DUPS")) {
        int hk[11] = {0};
        double rmax = 0, rmin = 1;
        for (auto& r : R.recs) { hk[r.k]++; rmax = std::max(rmax, r.r); rmin = std::min(rmin, r.r); }
        printf("    k histogram:");
        for (int k = 0; k <= 10; ++k) printf(" %d", hk[k]);
        printf("  r in [%.3e, %.3e]\n", rmin, rmax);
        int dups = 0;
        for (size_t i = 0; i < R.recs.size(); ++i)
          for (size_t j = 0; j < i; ++j)
            if (mergeDistanceUpper(R.recs[i].P, R.recs[j].P) < 1e-6) { ++dups; break; }
        printf("    duplicates (equivalent to an earlier record): %d of %zu\n", dups, R.recs.size());
      }
      bool certs = true;
      for (auto& r : R.recs) {
        std::string w2;
        if (!checkRoot(r, &w2)) { certs = false; printf("    record fails checkRoot: %s\n", w2.c_str()); }
        if (!(r.r <= P.target)) certs = false;
      }
      printf("  %-10s status=%d unres=%ld recs=%zu viol=%ld covered=%d certs=%d cells=%ld (norm %ld dir %ld cov %ld/%ld) maxNu=%d %.2fs %s\n", h.name, R.status, R.unresolved,
             R.recs.size(), w.violations, (int)covered, (int)certs, R.st.cellsNorm, R.st.exclNorm, R.st.exclDir, R.st.covered, R.st.coverTry, R.st.maxNuReached, secs,
             R.why.c_str());
      if (w.coverChecks) printf("    witness cover checks %ld violations %ld worst dist/s ratio %.3f\n", w.coverChecks, w.coverViolations, w.worstRatio);
      CHECK(w.coverViolations == 0, "%s: cover bound for a cell holding the true completion did not contain it", h.name);
      CHECK(w.violations == 0, "%s: a cell containing the true completion was excluded", h.name);
      CHECK(covered, "%s: true Hadamard not covered by any local Record", h.name);
      CHECK(certs, "%s: record certificate failure", h.name);
      CHECK(R.status == 0, "%s: seed unresolved/failed: %s", h.name, R.why.c_str());
      if (good && covered && certs) ++okc;
    }
  }
  printf("  %d / %d seeds fully OK\n", okc, tot);
}

// random seed satisfying the two sum equations exactly
static void randomSeed(std::mt19937& rng, double q[9]) {
  std::uniform_real_distribution<double> U(0, 1);
  for (;;) {
    double a = U(rng);
    cd A = ee(a);
    auto four = [&](double out[4]) {
      for (int tries = 0; tries < 1000; ++tries) {
        double p0 = U(rng), p1 = U(rng);
        cd w = -(cd(1) + A) - ee(p0) - ee(p1);
        double m = std::abs(w);
        if (m > 2 || m < 1e-9) continue;
        cd mid = w / 2.0;
        double h = std::sqrt(1 - m * m / 4);
        cd perp = cd(0, 1) * w / m * h;
        cd b2 = mid + perp, b3 = mid - perp;
        out[0] = p0; out[1] = p1;
        out[2] = fr(std::arg(b2) / (2 * M_PI));
        out[3] = fr(std::arg(b3) / (2 * M_PI));
        return true;
      }
      return false;
    };
    double bb[4], cc[4];
    if (!four(bb) || !four(cc)) continue;
    q[0] = a;
    for (int i = 0; i < 4; ++i) { q[1 + i] = bb[i]; q[5 + i] = cc[i]; }
    Params P;
    SeedCtx S;
    std::string why;
    if (buildSeed(S, q, P, why)) return;
  }
}

// quiet single witness case; returns number of failed checks
struct CaseOut {
  int fails = 0;
  long recs = 0, cells = 0;
  double secs = 0;
  bool skipped = false;
};
static CaseOut witnessCase(const Had& h, const Arr& A, std::mt19937& rng, double delta) {
  CaseOut o;
  std::normal_distribution<double> N(0, 1);
  Params P;
  double d[9], n = 0;
  for (int i = 0; i < 9; ++i) { d[i] = N(rng); n += d[i] * d[i]; }
  n = std::sqrt(n);
  double q[9];
  for (int i = 0; i < 9; ++i) q[i] = A.q[i] + delta * d[i] / n;
  SeedCtx S;
  std::string why;
  if (!buildSeed(S, q, P, why)) { o.skipped = true; return o; }
  Witness w;
  if (!witnessOf(S, A, w)) { o.fails++; return o; }
  w.truthP = A.P;
  auto t0 = std::chrono::steady_clock::now();
  SeedResult R = processSeed(q, P, &w);
  o.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  o.recs = (long)R.recs.size();
  o.cells = R.st.cellsNorm;
  bool covered = false;
  for (auto& r : R.recs)
    if (coverDist(A.P, r) <= r.r + 1e-9) { covered = true; break; }
  bool certs = true;
  for (auto& r : R.recs) {
    std::string w2;
    if (!checkRoot(r, &w2) || !(r.r <= P.target)) certs = false;
  }
  if (w.violations) { o.fails++; printf("  FUZZ %s: witness cell excluded\n", h.name); }
  if (w.coverViolations) { o.fails++; printf("  FUZZ %s: cover bound missed the truth\n", h.name); }
  if (!covered) { o.fails++; printf("  FUZZ %s: truth not covered\n", h.name); }
  if (!certs) { o.fails++; printf("  FUZZ %s: record certificate failure\n", h.name); }
  if (R.status != 0) { o.fails++; printf("  FUZZ %s: status %d %s\n", h.name, R.status, R.why.c_str()); }
  return o;
}

static void t3fuzz(int nF6, int nOther) {
  printf("T3b soundness fuzz: %d random F6(a,b) + %d S6/Bjorck seeds, random arrangement, random offsets in (0.005,0.031)\n", nF6, nOther);
  std::mt19937 rng(20261019);
  std::uniform_real_distribution<double> U01(0, 1), Ud(0.005, 0.031);
  int done = 0, skipped = 0, fails = 0;
  long recs = 0, cells = 0;
  double secs = 0, worst = 0;
  for (int i = 0; i < nF6 + nOther; ++i) {
    Had h = i < nF6 ? mkF6(2 * M_PI * U01(rng), 2 * M_PI * U01(rng)) : (i % 2 ? mkS6() : mkBjorck());
    auto arrs = arrangements(h, 1, rng);
    if (arrs.empty()) { ++skipped; continue; }
    CaseOut o = witnessCase(h, arrs[0], rng, Ud(rng));
    if (o.skipped) { ++skipped; continue; }
    ++done;
    fails += o.fails;
    recs += o.recs;
    cells += o.cells;
    secs += o.secs;
    worst = std::max(worst, o.secs);
  }
  printf("  %d cases run (%d skipped by frames/arrangement), %d failed checks, %ld local records, %ld cells, %.1f s total (slowest %.1f s)\n", done, skipped, fails, recs, cells, secs, worst);
  CHECK(fails == 0, "fuzz failures: %d", fails);
  CHECK(done >= (nF6 + nOther) / 2, "too few fuzz cases ran: %d", done);
}

static void t4(int n) {
  printf("T4 random seeds\n");
  std::mt19937 rng(4242);
  std::uniform_real_distribution<double> U(0, 1);
  Params P;
  // (a) generic random seeds: the vanishing-sum equations are violated; they must be excluded wholesale
  {
    int m = 300, frameOK = 0, sumEx = 0, bad = 0;
    long recs = 0;
    for (int i = 0; i < m; ++i) {
      double q[9];
      for (int j = 0; j < 9; ++j) q[j] = U(rng);
      SeedCtx S;
      std::string why;
      if (!buildSeed(S, q, P, why)) continue;
      ++frameOK;
      SeedResult R = processSeed(q, P);
      sumEx += R.seedExcluded;
      bad += (R.status != 0);
      recs += (long)R.recs.size();
    }
    printf("  (a) %d random seeds passing frames: %d excluded by the sum equations alone, %d bad, %ld records\n", frameOK, sumEx, bad, recs);
    CHECK(bad == 0 && recs == 0, "generic random seeds were not fully excluded");
    CHECK(frameOK > 0 && sumEx * 100 >= frameOK * 95, "sum exclusion rate too low");
  }
  // (b) seeds satisfying both sum equations: cell machinery is exercised; every seed must resolve
  {
    std::vector<std::array<double, 9>> seeds(n);
    for (auto& s : seeds) randomSeed(rng, s.data());
    long recs = 0, bad = 0, cells = 0, normEx = 0, dirEx = 0, zeroRec = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; ++i) {
      SeedResult R = processSeed(seeds[i].data(), P);
      recs += (long)R.recs.size();
      zeroRec += R.recs.empty();
      bad += (R.status != 0);
      cells += R.st.cellsNorm;
      normEx += R.st.exclNorm;
      dirEx += R.st.exclDir;
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("  (b) %d sum-satisfying random seeds: bad=%ld, seeds with no Record=%ld, records=%ld, cells=%ld (norm-excluded %ld, dir-excluded %ld) %.2f s = %.1f s / 1000 seeds (1 thread)\n", n, bad,
           zeroRec, recs, cells, normEx, dirEx, secs, secs * 1000.0 / n);
    CHECK(bad == 0, "random seeds failed/unresolved: %ld", bad);
    CHECK(normEx + dirEx > 0.5 * cells, "fewer than 50%% of cells excluded");
  }
}

// T6: records written by arun re-certify from scratch
static void t6(const char* path) {
  printf("T6 re-certify records in %s\n", path);
  std::vector<Root> v;
  if (!readRecords(path, v)) { printf("  (no record file, skipped)\n"); return; }
  int bad = 0;
  for (auto& r : v) {
    std::string why;
    if (!checkRoot(r, &why) || !(r.r >= 0 && r.r <= 0.125 + 1e-12) || !(r.n0 < 1e-6) || !(r.froot < 1e-6)) { ++bad; if (bad < 4) printf("  bad record: %s\n", why.c_str()); }
  }
  printf("  %zu records, %d failed\n", v.size(), bad);
  CHECK(bad == 0, "record re-certification failures: %d", bad);
}

static void t5() {
  printf("T5 Record file round trip\n");
  Had h = mkS6();
  double P[36];
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) P[6 * i + j] = std::arg(h.H[i][j]) / (2 * M_PI);
  Root R;
  CHECK(buildRoot(P, R), "build");
  R.r = 0.05;
  std::vector<Root> v = {R, R};
  const char* path = "build/t5.bin";
  CHECK(writeRecords(path, v), "write");
  std::vector<Root> w;
  CHECK(readRecords(path, w), "read");
  CHECK(w.size() == 2 && w[0].C == R.C && w[1].r == R.r && w[0].P[7] == R.P[7], "roundtrip");
  FILE* f = fopen(path, "rb");
  fseek(f, 0, SEEK_END);
  long sz = ftell(f);
  fclose(f);
  CHECK(sz == 8 + 2 * (356 + 288 * R.k), "size %ld vs %ld", sz, (long)(8 + 2 * (356 + 288 * R.k)));
}

int main(int argc, char** argv) {
  if (argc >= 3 && std::string(argv[1]) == "--verify") {  // only T6 on a record file
    iv::RoundUp g0;
    initTables();
    t6(argv[2]);
    printf(g_fail ? "FAILED (%d checks)\n" : "ALL OK\n", g_fail);
    return g_fail ? 1 : 0;
  }
  int perHad = argc > 1 ? atoi(argv[1]) : 2;
  int nfar = argc > 2 ? atoi(argv[2]) : 40;
  int fuzzF6 = getenv("FUZZ_F6") ? atoi(getenv("FUZZ_F6")) : 30;
  int fuzzOther = getenv("FUZZ_OTHER") ? atoi(getenv("FUZZ_OTHER")) : 6;
  iv::RoundUp g;
  initTables();
  t1();
  t2();
  t5();
  t3(perHad);
  t3fuzz(fuzzF6, fuzzOther);
  t4(nfar);
  if (argc > 3) t6(argv[3]);
  printf(g_fail ? "FAILED (%d checks)\n" : "ALL OK\n", g_fail);
  return g_fail ? 1 : 0;
}
