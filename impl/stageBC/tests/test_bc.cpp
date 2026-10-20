// Stage B / Stage C unit and integration tests.  `test_bc` runs everything (exit 0 iff all pass).
//   test_bc [name...]    names: io root rel stageB graph cluster single edge triple soundness
#include "../stageB.hpp"
#include "../stageC.hpp"
#include "fixtures.hpp"
#include <chrono>
#include <cstdio>
#include <set>

using namespace mb;
static int failures = 0, checks = 0;
#define CHECK(c) do { ++checks; if (!(c)) { ++failures; std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)
#define CHECKM(c, ...) do { ++checks; if (!(c)) { ++failures; std::printf("FAIL %s:%d: %s :: ", __FILE__, __LINE__, #c); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static std::mt19937_64 rng(20261019);
static double urand(double a, double b) { return std::uniform_real_distribution<double>(a, b)(rng); }

static Record rootAt(const double* P0) {
  double P[36];
  std::memcpy(P, P0, sizeof P);
  Record R;
  RootInfo info;
  bool s = snapSolution(P);
  bool c = s && certifyRoot(P, R, &info);
  CHECKM(c, "root at fixture failed: snap=%d why=%s", (int)s, info.why.c_str());
  return R;
}
static Record fourierRoot(double a, double b, double r) {
  double P[36];
  fx::fourierPhases(a, b, P);
  Record R = rootAt(P);
  R.r = r;
  return R;
}

// ------------------------------------------------------------------ io
static void test_io() {
  Record R = fourierRoot(0.3, 0.5, 0.05);
  Record Z = fourierRoot(0.0, 0.0, 0.0);
  std::vector<Record> v{R, Z};
  CHECK(writeRecords("build/t_io.bin", v));
  std::vector<Record> back;
  std::string err;
  CHECK(readRecords("build/t_io.bin", back, err));
  CHECK(back.size() == 2 && back[0].k == R.k && back[1].k == Z.k);
  CHECK(std::memcmp(back[0].p, R.p, sizeof R.p) == 0 && back[0].T == R.T && back[0].r == R.r && back[0].C == R.C);
  // size = 8 + sum(356 + 288k)
  FILE* f = std::fopen("build/t_io.bin", "rb");
  std::fseek(f, 0, SEEK_END);
  long sz = std::ftell(f);
  std::fclose(f);
  CHECK(sz == 8 + (356 + 288L * R.k) + (356 + 288L * Z.k));
  // zero-count file = 8 bytes
  CHECK(writeRecords("build/t_io0.bin", std::vector<Record>{}));
  f = std::fopen("build/t_io0.bin", "rb");
  std::fseek(f, 0, SEEK_END);
  CHECK(std::ftell(f) == 8);
  std::fclose(f);
  // trailing byte => rejected
  f = std::fopen("build/t_io.bin", "ab");
  std::fputc(1, f);
  std::fclose(f);
  std::vector<Record> bad;
  CHECK(!readRecords("build/t_io.bin", bad, err));
  std::string why;
  CHECK(validateRecord(R, 0.126, why));
  Record B2 = R;
  B2.r = 0.2;
  CHECK(!validateRecord(B2, 0.126, why));
  B2 = R;
  B2.froot = 1e-3;
  CHECK(!validateRecord(B2, 0.126, why));
  B2 = R;
  B2.p[3] = 5;
  CHECK(!validateRecord(B2, 0.126, why));
  std::printf("io: ok\n");
}

// ------------------------------------------------------------------ root certificates
static void test_root() {
  const int N = 6;
  for (int t = 0; t < N; ++t) {
    double a = urand(-1.2, 1.2), b = urand(-1.2, 1.2);
    Record R = fourierRoot(a, b, 0.0);
    CHECK(R.k >= 1 && R.k <= 10);
    CHECK(R.id < 1e-4 && R.ort < 1e-5 && R.hor < 1e-5 && R.C < 20 && R.mnorm < 2 && R.froot < 1e-9);
    // independent check of the Root inequalities on exact solutions of the family near the centre
    for (int s = 0; s < 200; ++s) {
      double H[36];
      fx::fourierPhases(a + urand(-0.12, 0.12), b + urand(-0.12, 0.12), H);
      double d[36];
      fx::aligned(R.p, H, d);
      double sn = fx::norm36(d);
      if (sn > 0.125) continue;
      double eta[10] = {0};
      for (int c = 0; c < R.k; ++c)
        for (int m = 0; m < 36; ++m) eta[c] += R.col(c)[m] * d[m];
      double nu[36];
      for (int m = 0; m < 36; ++m) { nu[m] = d[m]; for (int c = 0; c < R.k; ++c) nu[m] -= R.col(c)[m] * eta[c]; }
      double nn = fx::norm36(nu);
      double nP = R.n0 + R.id * sn + R.C * sn * sn;
      CHECKM(nn <= nP + 1e-12, "root-normal violated: %g > %g (s=%g)", nn, nP, sn);
      double e2 = 0;
      for (int c = 0; c < R.k; ++c) e2 += eta[c] * eta[c];
      CHECK(sn * sn <= (1 + R.ort) * e2 + nn * nn + 1e-12);
    }
  }
  // the Bjorck circulant: a defect-0 point should give k=0 or small k; just certify
  double B[36];
  fx::bjorckPhases(B);
  double Bs[36];
  std::memcpy(Bs, B, sizeof B);
  Interval Fn = residualNorm(B);
  CHECKM(Fn.hi < 1e-12, "bjorck not Hadamard: %g", Fn.hi);
  Record R = rootAt(Bs);
  std::printf("root: ok (Bjorck C6 k=%d mnorm=%.3f C=%.3f)\n", R.k, R.mnorm, R.C);
  // a random phase matrix is not a solution: snapping must fail or certification must refuse
  double Rnd[36];
  for (int m = 0; m < 36; ++m) Rnd[m] = urand(-0.5, 0.5);
  Record Rr;
  RootInfo info;
  bool s = snapSolution(Rnd);
  bool c = s && certifyRoot(Rnd, Rr, &info);
  CHECK(!c || residualNorm(Rr.p).hi < 1e-9);  // if it did converge to some Hadamard, the certificate is still honest
  double Rnd2[36];
  for (int m = 0; m < 36; ++m) Rnd2[m] = urand(-0.5, 0.5);
  CHECK(!certifyRoot(Rnd2, Rr, &info));  // unsnapped random matrix has froot >> 1e-9
}

// ------------------------------------------------------------------ relation soundness (the heart of Stage B)
// ------------------------------------------------------------------ symmetry maps Phi (phase:transport)
static void test_phi() {
  for (int t = 0; t < 4; ++t) {
    double a = urand(-1.2, 1.2), b = urand(-1.2, 1.2);
    Record R = fourierRoot(a, b, 0.0);
    std::vector<Phi> vs = reuseVariants(R.p);
    CHECK(!vs.empty() && vs[0].isIdentity());
    int nontriv = 0;
    for (auto& f : vs) {
      Phi g = f.inverse();
      double x[36], y[36], z[36];
      for (int m = 0; m < 36; ++m) x[m] = urand(-1, 1);
      applyPhiVec(f, x, y);
      applyPhiVec(g, y, z);
      double e = 0, nx = 0, ny = 0;
      for (int m = 0; m < 36; ++m) { e += (z[m] - x[m]) * (z[m] - x[m]); nx += x[m] * x[m]; ny += y[m] * y[m]; }
      CHECK(e == 0);                       // inverse is exact (a signed permutation)
      CHECK(std::fabs(nx - ny) < 1e-12);  // isometry
      // Pi commutes with Phi
      double px[36], ppx[36], pxp[36];
      PiD(x, px); applyPhiVec(f, px, ppx); PiD(y, pxp);
      double ec = 0;
      for (int m = 0; m < 36; ++m) ec += (ppx[m] - pxp[m]) * (ppx[m] - pxp[m]);
      CHECK(ec < 1e-24);
      if (f.isIdentity()) continue;
      ++nontriv;
      Record Q = applyPhi(R, f);
      // transported Root: residual invariance and the Root inequalities on transported exact solutions
      CHECKM(residualNorm(Q.p).hi <= R.froot + 1e-12, "residual not invariant: %g > %g", residualNorm(Q.p).hi, R.froot);
      Record F;
      RootInfo info;
      bool c = certifyRoot(Q.p, F, &info);
      CHECKM(c && F.k == R.k, "transported root does not recertify: %s", info.why.c_str());
      for (int s = 0; s < 100; ++s) {
        double H[36], Hf[36];
        fx::fourierPhases(a + urand(-0.12, 0.12), b + urand(-0.12, 0.12), H);
        applyPhiVec(f, H, Hf);
        double d[36];
        fx::aligned(Q.p, Hf, d);
        double sn = fx::norm36(d);
        if (sn > 0.125) continue;
        CHECK(residualNorm(Hf).hi < 1e-10);  // Phi H is an exact solution
        double eta[10] = {0};
        for (int cc = 0; cc < Q.k; ++cc)
          for (int m = 0; m < 36; ++m) eta[cc] += Q.col(cc)[m] * d[m];
        double nu[36];
        for (int m = 0; m < 36; ++m) { nu[m] = d[m]; for (int cc = 0; cc < Q.k; ++cc) nu[m] -= Q.col(cc)[m] * eta[cc]; }
        double nn = fx::norm36(nu);
        CHECKM(nn <= Q.n0 + Q.id * sn + Q.C * sn * sn + 1e-12, "transported root-normal violated: %g (s=%g)", nn, sn);
      }
      // a relation from the transported parent to the original is the relation of the pulled-back root (isometry)
      Record Qr = Q; Qr.r = 0.05;
      Record Rr = R; Rr.r = 0.05;
      Rel r1 = makeRel(Qr, R, 0.249);
      Rel r2 = makeRel(Rr, applyPhi(R, g), 0.249);
      CHECK(r1.ok == r2.ok);
      if (r1.ok && r2.ok) CHECK(std::fabs(r1.eval(0.02) - r2.eval(0.02)) < 1e-9 && std::fabs(r1.dqNorm - r2.dqNorm) < 1e-9);
    }
    std::printf("phi: %zu variants (%d non-identity) at fourier(%.2f,%.2f)\n", vs.size(), nontriv, a, b);
  }
}

// ------------------------------------------------------------------ tensor improvement of C (core/tensor_c.hpp)
static void test_tensor() {
  int better = 0, total = 0;
  for (int t = 0; t < 12; ++t) {
    int g = 15 + (int)(urand(0, 11));
    std::vector<double> HV((size_t)36 * g), LV((size_t)30 * g);
    double sc = urand(0.2, 0.6);
    for (auto& x : HV) x = urand(-sc, sc);
    for (auto& x : LV) x = urand(-sc, sc);
    // explicit S (30 x 216): row a = M_{:,a} (x) (e_i - e_j),  M = HV LV^T (36 x 30)
    std::vector<double> M((size_t)36 * 30, 0.0);
    for (int m = 0; m < 36; ++m)
      for (int a = 0; a < 30; ++a) { double s = 0; for (int c = 0; c < g; ++c) s += HV[(size_t)m * g + c] * LV[(size_t)a * g + c]; M[(size_t)m * 30 + a] = s; }
    std::vector<double> S((size_t)30 * 216, 0.0);
    for (int a = 0; a < 30; ++a) {
      int i, j; tc::edgeOf(a, i, j);
      for (int m = 0; m < 36; ++m)
        for (int q = 0; q < 6; ++q) S[(size_t)a * 216 + 6 * m + q] = M[(size_t)m * 30 + a] * ((q == i) - (q == j));
    }
    std::vector<double> G((size_t)30 * 30), V, ev;
    for (int a = 0; a < 30; ++a)
      for (int b = 0; b < 30; ++b) { double s = 0; for (int c = 0; c < 216; ++c) s += S[(size_t)a * 216 + c] * S[(size_t)b * 216 + c]; G[(size_t)a * 30 + b] = s; }
    jacobiEig(G, 30, V, ev);
    double Snorm = std::sqrt(std::max(0.0, ev[29]));
    double truth = cn().pi.hi * std::sqrt(6.0) * Snorm;
    tc::TensorInfo info;
    double C = tc::tensorC(HV, LV, g, cn().pi, &info);
    CHECKM(std::isfinite(C) && C >= truth * (1 - 1e-9), "tensor C %g below pi sqrt6 ||S|| = %g", C, truth);
    CHECKM(C <= truth * 1.5 + 1e-9, "tensor C %g is far above the truth %g (beta=%g)", C, truth, info.beta);
    // plain bound for comparison
    std::vector<double> MM((size_t)30 * 30), V2, ev2;
    for (int a = 0; a < 30; ++a)
      for (int b = 0; b < 30; ++b) { double s = 0; for (int m = 0; m < 36; ++m) s += M[(size_t)m * 30 + a] * M[(size_t)m * 30 + b]; MM[(size_t)a * 30 + b] = s; }
    jacobiEig(MM, 30, V2, ev2);
    double plain = cn().pi.hi * std::sqrt(8.0 * std::max(0.0, ev2[29]));
    ++total; if (C < plain) ++better;
    // the bound is monotone: scaling the factors by 1.1 scales both sides by 1.21
    std::vector<double> HV2 = HV; for (auto& x : HV2) x *= 1.1;
    double C2 = tc::tensorC(HV2, LV, g, cn().pi);
    CHECK(C2 >= C * 1.1 * (1 - 1e-9));
  }
  // degenerate proposal: a factor pair whose Gram is far from diagonalisable within the entry bound -> +inf, never small
  std::vector<double> HVb((size_t)36 * 3, 3.0), LVb((size_t)30 * 3, 3.0);
  double Cb = tc::tensorC(HVb, LVb, 3, cn().pi);
  CHECK(std::isfinite(Cb) ? Cb >= cn().pi.hi * std::sqrt(6.0) * 1.0 : true);
  std::printf("tensor: bound >= pi sqrt6 ||S|| on %d random factor pairs; tensor < plain in %d/%d\n", total, better, total);
}

static void test_rel() {
  long tested = 0;
  for (int trial = 0; trial < 4; ++trial) {
    double a0 = urand(-1, 1), b0 = urand(-1, 1);
    Record P = fourierRoot(a0, b0, 0.0);
    P.r = (trial == 0) ? 0.05 : (trial == 1 ? 0.09 : 0.125);
    // second root Q at a nearby family point
    double da = urand(-0.06, 0.06), db = urand(-0.06, 0.06);
    Record Q = fourierRoot(a0 + da, b0 + db, 0.0);
    for (double cut : {0.249}) {
      Rel rel = makeRel(P, Q, cut);
      if (!rel.ok) { std::printf("  (relation declined, dq=%g)\n", rel.dqNorm); continue; }
      // exact solutions near P
      for (int s = 0; s < 400; ++s) {
        double H[36];
        fx::fourierPhases(a0 + urand(-0.15, 0.15), b0 + urand(-0.15, 0.15), H);
        double d[36];
        fx::aligned(P.p, H, d);
        if (fx::norm36(d) > P.r) continue;
        double dq[36];
        {
          Interval di[36];
          alignedDisp(P.p, Q.p, di);
          for (int m = 0; m < 36; ++m) dq[m] = di[m].mid();
        }
        double xi2 = 0;
        for (int c = 0; c < P.k; ++c) {
          double e = 0, eq = 0;
          for (int m = 0; m < 36; ++m) { e += P.col(c)[m] * d[m]; eq += P.col(c)[m] * dq[m]; }
          xi2 += (e - eq) * (e - eq);
        }
        double t = std::sqrt(xi2) * (1 + 1e-9) + 1e-13;
        double hp[36];
        for (int m = 0; m < 36; ++m) hp[m] = d[m] - dq[m];
        double dtrue = fx::norm36(hp);
        double bound = rel.eval(t);
        ++tested;
        CHECKM(dtrue <= bound + 1e-12, "Dhat violated: true %g > bound %g (t=%g)", dtrue, bound, t);
        // rho certificate: if t <= rho(a) then dtrue <= a
        double a = 0.07;
        double rho = rel.rho(a);
        if (rho > 0 && t <= rho) CHECK(dtrue <= a);
      }
    }
  }
  std::printf("rel: %ld exact solutions checked against Dhat\n", tested);
}

// ------------------------------------------------------------------ Stage B end to end
// H is covered if some atlas ball contains H or Phi H for one of the maps the refiner actually used in a merge
// (refine:atlas-assignment: a merge through Phi certifies Phi Q, i.e. the Phi-image of the solutions of the box).
static bool covered(const std::vector<Record>& atlas, const double* H, const std::vector<Phi>* maps = nullptr) {
  std::vector<Phi> ms;
  ms.push_back(Phi::identity());
  if (maps) ms.insert(ms.end(), maps->begin(), maps->end());
  for (auto& f : ms) {
    double Hf[36];
    applyPhiVec(f, H, Hf);
    for (auto& Z : atlas) {
      Interval d[36];
      alignedDisp(Z.p, Hf, d);
      double s = 0;
      for (int m = 0; m < 36; ++m) s += d[m].mid() * d[m].mid();
      if (std::sqrt(s) <= Z.r + 1e-12) return true;
    }
  }
  return false;
}
static void run_stageB(double a0, double b0, double R, double tau, double w, int samples, const char* label, bool noPhi = false) {
  Record P = fourierRoot(a0, b0, R);
  Refiner rf(tau, w);
  auto t0 = std::chrono::steady_clock::now();
  bool ok = rf.processParent(P);
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("stageB %s: R=%.3f tau=%.3f w=%.4f -> ok=%d roots=%zu (merged %ld, direct %ld, marked %ld) %.1fs\n", label, R, tau, w, (int)ok,
              rf.atlas.recs.size(), rf.st.merged, rf.st.direct, rf.st.marked, secs);
  CHECK(ok);
  if (!ok) return;
  for (auto& Z : rf.atlas.recs) {
    std::string why;
    CHECK(Z.r <= tau && Z.r >= 0);
    Record Zv = Z;
    CHECK(validateRecord(Zv, 0.126, why));
    CHECK(residualNorm(Z.p).hi <= Z.froot);  // stored froot is a true upper bound
  }
  int cov = 0, tested = 0, bad = 0;
  for (int s = 0; s < samples * 20 && tested < samples; ++s) {
    double H[36];
    fx::fourierPhases(a0 + urand(-1.5 * R, 1.5 * R), b0 + urand(-1.5 * R, 1.5 * R), H);
    double d[36];
    fx::aligned(P.p, H, d);
    if (fx::norm36(d) > R) continue;
    ++tested;
    if (covered(rf.atlas.recs, H, &rf.usedMaps)) ++cov; else ++bad;
  }
  CHECKM(bad == 0, "%d of %d exact solutions NOT covered by the child balls", bad, tested);
  std::printf("  cover check: %d/%d exact solutions covered (%zu Phi maps used in merges)\n", cov, tested, rf.usedMaps.size());
  if (noPhi) {  // same parent without symmetry maps: must also succeed (Phi is reuse only), typically with more roots
    Refiner r2(tau, w);
    r2.usePhi = false;
    bool ok2 = r2.processParent(P);
    CHECK(ok2);
    std::printf("  without Phi: roots=%zu (with Phi: %zu)\n", r2.atlas.recs.size(), rf.atlas.recs.size());
  }
}
static void test_stageB() {  // quick: runs in `make test`
  run_stageB(0.3, 0.5, 0.05, 0.07, 0.013, 400, "small", true);
  run_stageB(-0.6, 0.2, 0.03, 0.016, 0.005, 400, "fine-small");
}
static void test_stageBslow() {  // `make test-slow`: paper-scale radii (R=.125 -> .07 and .07 -> .016)
  run_stageB(-0.6, 0.2, 0.125, 0.07, 0.013, 600, "mid");
  run_stageB(0.9, -0.4, 0.07, 0.016, 0.005, 400, "fine");
}

// ------------------------------------------------------------------ markRec: soundness (every marked box is within rho) and completeness
static void test_mark() {
  Record P = fourierRoot(0.3, 0.5, 0.06);
  Record Q = fourierRoot(0.32, 0.49, 0.0);
  double w = 0.008, tau = 0.04;
  Refiner rf(tau, w);
  Refiner::Grid g;
  CHECK(rf.setupGrid(P, g));
  Rel rel = makeRel(P, Q, 0.249);
  CHECK(rel.ok);
  double rho = rel.rho(tau);
  CHECK(rho > 0 && rho < tau);
  double tmax = 0;
  long cnt = rf.markRec(rel, rho, g, true, tmax);
  CHECK(cnt > 0 && tmax <= rho);
  long marked = 0, shouldBe = 0, wrong = 0, missed = 0;
  const int k = g.k;
  for (uint64_t idx = 0; idx < g.total; ++idx) {
    uint64_t u = idx;
    double far2 = 0;
    for (int a = 0; a < k; ++a) {
      long b = (long)(u % (uint64_t)g.n) - g.L;
      u /= (uint64_t)g.n;
      double d = std::fabs(b * w - rel.etaq[a].mid()) + w / 2;
      far2 += d * d;
    }
    double far = std::sqrt(far2);
    bool m = g.isMarked(idx);
    marked += m;
    if (m && far > rho) ++wrong;                  // soundness: a marked box must have farthest corner <= rho
    if (!m && far < rho - 1e-5) ++missed;         // completeness (not required, but expected)
    if (far < rho - 1e-5) ++shouldBe;
  }
  CHECKM(wrong == 0, "%ld marked boxes violate the farthest-corner bound", wrong);
  CHECKM(missed == 0, "%ld boxes within rho not marked", missed);
  std::printf("mark: rho=%.4f, %ld boxes marked (expected >= %ld), tmax=%.4f\n", rho, marked, shouldBe, tmax);
}

// ------------------------------------------------------------------ graph
static bool bruteExist(const Graph& G) {
  int n = G.n;
  auto O = [&](int i, int j) { return (G.o(i)[j >> 6] >> (j & 63)) & 1; };
  auto U = [&](int i, int j) { return (G.u(i)[j >> 6] >> (j & 63)) & 1; };
  std::vector<std::vector<int>> cl;
  std::vector<int> cur;
  std::function<void(int)> rec = [&](int s) {
    if (cur.size() == 6) { cl.push_back(cur); return; }
    for (int v = s; v < n; ++v) {
      bool ok = true;
      for (int u : cur) if (!O(u, v)) { ok = false; break; }
      if (!ok) continue;
      cur.push_back(v);
      rec(v + 1);
      cur.pop_back();
    }
  };
  rec(0);
  for (size_t x = 0; x < cl.size(); ++x)
    for (size_t y = 0; y < cl.size(); ++y) {
      if (x == y) continue;
      bool ok = true;
      for (int a : cl[x]) for (int b : cl[y]) if (a == b || !U(a, b)) ok = false;
      if (ok) return true;
    }
  return false;
}
static void test_graph() {
  // planted configuration is found
  for (int t = 0; t < 20; ++t) {
    int n = 30 + (int)(rng() % 40);
    Graph G;
    G.init(n);
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) { if (urand(0, 1) < 0.25) G.addO(i, j); if (urand(0, 1) < 0.25) G.addU(i, j); }
    std::vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i;
    std::shuffle(perm.begin(), perm.end(), rng);
    for (int a = 0; a < 6; ++a) for (int b = a + 1; b < 6; ++b) { G.addO(perm[a], perm[b]); G.addO(perm[6 + a], perm[6 + b]); }
    for (int a = 0; a < 6; ++a) for (int b = 0; b < 6; ++b) G.addU(perm[a], perm[6 + b]);
    GraphSearch gs(G, 50000000);
    CHECK(gs.exist() && !gs.gaveUp);
    std::vector<int> keep = pruneVertices(G);
    Graph H = induced(G, keep);
    GraphSearch gh(H, 50000000);
    CHECK(gh.exist() && !gh.gaveUp);
  }
  // exactness vs brute force on small random graphs, with and without pruning
  int mism = 0, pos = 0;
  for (int t = 0; t < 400; ++t) {
    int n = 12 + (int)(rng() % 7);
    double pO = urand(0.4, 0.95), pU = urand(0.3, 0.95);
    Graph G;
    G.init(n);
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) { if (urand(0, 1) < pO) G.addO(i, j); if (urand(0, 1) < pU) G.addU(i, j); }
    bool br = bruteExist(G);
    GraphSearch gs(G, 100000000);
    bool r1 = gs.exist();
    std::vector<int> keep = pruneVertices(G);
    Graph H = induced(G, keep);
    GraphSearch gh(H, 100000000);
    bool r2 = gh.exist();
    if (br) ++pos;
    if (r1 != br || r2 != br || gs.gaveUp || gh.gaveUp) ++mism;
  }
  CHECKM(mism == 0, "%d mismatches vs brute force (positives=%d)", mism, pos);
  CHECK(pos > 5);
  // cap exhaustion never yields a negative: dense graph, tiny cap => gaveUp
  {
    Graph G;
    G.init(40);
    for (int i = 0; i < 40; ++i) for (int j = i + 1; j < 40; ++j) { if (urand(0, 1) < 0.9) G.addO(i, j); if (urand(0, 1) < 0.1) G.addU(i, j); }
    GraphSearch gs(G, 3);
    bool r = gs.exist();
    CHECK(gs.gaveUp && r);  // possibility reported
  }
  // a true negative with a generous cap completes without giving up
  {
    Graph G;
    G.init(50);
    for (int i = 0; i < 50; ++i) for (int j = i + 1; j < 50; ++j) { if (urand(0, 1) < 0.5) G.addO(i, j); }
    GraphSearch gs(G, 100000000);
    CHECK(!gs.exist() && !gs.gaveUp);
  }
  std::printf("graph: ok (%d positives among 400 random small graphs)\n", pos);
}

// ------------------------------------------------------------------ clustering covers every bin point
static void test_cluster() {
  for (int trial = 0; trial < 3; ++trial) {
    int M = 96;
    std::vector<std::array<uint8_t, 5>> bins;
    int clumps = 4 + trial * 3;
    for (int c = 0; c < clumps; ++c) {
      std::array<int, 5> base;
      for (int j = 0; j < 5; ++j) base[j] = (int)(rng() % 90);
      for (int t = 0; t < 120; ++t) {
        std::array<uint8_t, 5> b;
        for (int j = 0; j < 5; ++j) b[j] = (uint8_t)std::min(95, std::max(0, base[j] + (int)(rng() % 7) - 3));
        bins.push_back(b);
      }
    }
    ClusterOut cl = clusterBins(bins, M);
    CHECKM(cl.ok, "cluster failed: %s", cl.why.c_str());
    if (!cl.ok) continue;
    int bad = 0, total = 0;
    for (auto& b : bins)
      for (int s = 0; s < 20; ++s) {
        double beta[6] = {0};
        for (int j = 1; j < 6; ++j) beta[j] = (b[j - 1] + urand(0, 1)) / M;
        bool cov = false;
        for (size_t c = 0; c < cl.v.size() && !cov; ++c) {
          double x[6], mean = 0;
          for (int j = 0; j < 6; ++j) { double d = beta[j] - cl.v[c][j]; x[j] = d - rnd(d); mean += x[j]; }
          mean /= 6;
          double d2 = 0;
          for (int j = 0; j < 6; ++j) d2 += (x[j] - mean) * (x[j] - mean);
          if (std::sqrt(d2) <= cl.s[c] + 1e-12) cov = true;
        }
        ++total;
        if (!cov) ++bad;
      }
    CHECKM(bad == 0, "%d/%d points of retained bins not covered by cluster balls", bad, total);
    std::printf("cluster: %zu bins -> %zu balls, max radius %.4f\n", bins.size(), cl.v.size(), *std::max_element(cl.s.begin(), cl.s.end()));
  }
}

// ------------------------------------------------------------------ single-vector tests keep every exact unbiased vector
static void test_single() {
  Record P0 = fourierRoot(0.3, 0.5, 0.016);
  for (int trial = 0; trial < 2; ++trial) {
    // exact Hadamard inside the ball: P0 + delta, delta from a family point at distance <= .8 r
    double Hh[36];
    if (trial == 0) std::memcpy(Hh, P0.p, sizeof Hh);
    else {
      for (;;) {
        double H[36], d[36];
        fx::fourierPhases(0.3 + urand(-0.03, 0.03), 0.5 + urand(-0.03, 0.03), H);
        fx::aligned(P0.p, H, d);
        double s = fx::norm36(d);
        if (s > 0.008 && s < 0.8 * P0.r) { for (int m = 0; m < 36; ++m) Hh[m] = P0.p[m] + d[m]; break; }
      }
    }
    CHECK(residualNorm(Hh).hi < 1e-12);
    std::vector<fx::UVec> V = fx::findUnbiasedVectors(Hh, 400, rng);
    CHECK(V.size() > 0);
    CRec C;
    initCRec(C, P0);
    Single s;
    int retainedAll = 0;
    for (auto& v : V) {
      bool ok = true;
      for (int M : {12, 24, 48, 96}) {
        Layer L;
        L.init(C, M);
        uint8_t b[5];
        for (int j = 1; j < 6; ++j) b[j - 1] = (uint8_t)std::min(M - 1, (int)std::floor(v.beta[j] * M));
        if (binExcluded(C, L, b, s)) { ok = false; CHECKM(false, "exact unbiased vector excluded at M=%d (trial %d)", M, trial); break; }
      }
      retainedAll += ok;
    }
    // power: random vectors are almost always excluded at M=12 (a bin is only excluded when it contains NO unbiased vector;
    // random bins contain none with overwhelming probability)
    Layer L12;
    L12.init(C, 12);
    int excl = 0, tot = 600;
    for (int t = 0; t < tot; ++t) {
      uint8_t b[5];
      for (int j = 0; j < 5; ++j) b[j] = (uint8_t)(rng() % 12);
      if (binExcluded(C, L12, b, s)) ++excl;
    }
    std::printf("single: trial %d: %zu exact vectors, all retained: %d; random bins excluded %d/%d at M=12\n", trial, V.size(), retainedAll, excl, tot);
    CHECK(retainedAll == (int)V.size());
    CHECK(excl > tot / 2);
  }
}

// ------------------------------------------------------------------ edge and joint tests keep every exact relation
static void test_edge() {
  Record P0 = fourierRoot(0.3, 0.5, 0.016);
  CRec C;
  initCRec(C, P0);
  long nO = 0, nU = 0, nOkept = 0, nUkept = 0, nOexcl = 0, nTot = 0;
  for (int trial = 0; trial < 3; ++trial) {
    double Hh[36];
    if (trial == 0) std::memcpy(Hh, P0.p, sizeof Hh);
    else {
      for (;;) {
        double H[36], d[36];
        fx::fourierPhases(0.3 + urand(-0.03, 0.03), 0.5 + urand(-0.03, 0.03), H);
        fx::aligned(P0.p, H, d);
        double s = fx::norm36(d);
        if (s > 0.005 && s < 0.9 * P0.r) { for (int m = 0; m < 36; ++m) Hh[m] = P0.p[m] + d[m]; break; }
      }
    }
    std::vector<fx::UVec> V = fx::findUnbiasedVectors(Hh, 300, rng);
    const double s = 0.03;
    std::vector<std::unique_ptr<Ball>> balls;
    for (auto& v : V) {
      // centre = beta - x0 with x0 mean-zero, ||x0|| <= 0.9 s (the true vector lies inside the ball)
      double x0[6], mean = 0, nn = 0;
      for (int j = 0; j < 6; ++j) { x0[j] = urand(-1, 1); mean += x0[j]; }
      mean /= 6;
      for (int j = 0; j < 6; ++j) { x0[j] -= mean; nn += x0[j] * x0[j]; }
      double sc = urand(0.3, 0.9) * s / std::sqrt(nn);
      double c[6];
      for (int j = 0; j < 6; ++j) c[j] = v.beta[j] - x0[j] * sc;
      balls.emplace_back(new Ball());
      makeBall(C, c, s, *balls.back());
    }
    for (size_t i = 0; i < V.size(); ++i)
      for (size_t j = i + 1; j < V.size(); ++j) {
        double ip = fx::ipAbs2(V[i].beta, V[j].beta);
        bool isO = ip < 1e-9, isU = std::fabs(ip - 6.0) < 1e-6;
        EdgeRes er = testEdge(C, *balls[i], *balls[j]);
        ++nTot;
        if (isO) { ++nO; nOkept += er.O; CHECKM(er.O, "orthogonal pair lost its O edge (trial %d, pair %zu,%zu)", trial, i, j); }
        else if (er.O) {}
        else ++nOexcl;
        if (isU) { ++nU; nUkept += er.U; CHECKM(er.U, "unbiased pair lost its U edge (trial %d, pair %zu,%zu)", trial, i, j); }
      }
  }
  std::printf("edge: %ld pairs; orthogonal %ld (kept %ld), unbiased %ld (kept %ld), non-O pairs with O edge removed: %ld\n", nTot, nO, nOkept, nU, nUkept, nOexcl);
  CHECK(nO > 0);
}

// ------------------------------------------------------------------ Taylor-remainder checks on arbitrary (non-solution) points
// For ANY centres, displacements within the radii and any directions y, the exact quantities must satisfy the inequalities
// that define the support functions (these are the pure analytic lemmas of exclusion.tex, tested numerically in doubles).
static void test_taylor() {
  const double TWO_PI = 6.283185307179586476925286766559;
  Record P0 = fourierRoot(0.3, 0.5, 0.016);
  CRec C;
  initCRec(C, P0);
  CRec C0 = C;
  C0.nhat = Interval(0.0);
  long checked = 0, tight = 0;
  for (int sample = 0; sample < 400; ++sample) {
    double Hh[36], delta[36] = {0};
    if (sample % 2) {
      for (;;) {
        double H[36];
        fx::fourierPhases(0.3 + urand(-0.03, 0.03), 0.5 + urand(-0.03, 0.03), H);
        fx::aligned(P0.p, H, delta);
        if (fx::norm36(delta) < 0.95 * P0.r) break;
      }
    }
    for (int m = 0; m < 36; ++m) Hh[m] = P0.p[m] + delta[m];
    double eta[10] = {0}, nu[36];
    for (int c = 0; c < P0.k; ++c) for (int m = 0; m < 36; ++m) eta[c] += P0.col(c)[m] * delta[m];
    for (int m = 0; m < 36; ++m) { nu[m] = delta[m]; for (int c = 0; c < P0.k; ++c) nu[m] -= P0.col(c)[m] * eta[c]; }
    std::unique_ptr<Ball> bl[2];
    double v[2][6], x[2][6], s[2], y[2][6], ye[2];
    for (int t = 0; t < 2; ++t) {
      s[t] = urand(0.005, 0.05);
      for (int j = 0; j < 6; ++j) { v[t][j] = urand(0, 1); y[t][j] = urand(-1, 1); }
      bl[t].reset(new Ball());
      makeBall(C, v[t], s[t], *bl[t]);
    }
    ye[0] = urand(-1, 1);
    ye[1] = urand(-1, 1);
    for (int rel = 0; rel < 2; ++rel) {
      const bool isO = rel == 0;
      EdgeData e;
      edgeData(*bl[0], *bl[1], e);
      double DE[2][6];
      for (int j = 0; j < 6; ++j) {
        if (isO) { DE[0][j] = e.DO[0][j].mid(); DE[1][j] = e.DO[1][j].mid(); }
        else { DE[0][j] = e.DU[j].mid(); DE[1][j] = 0; }
      }
      const int l = isO ? 2 : 1;
      double t0[6], t1[6];
      for (int j = 0; j < 6; ++j) {
        t0[j] = 0; t1[j] = 0;
        for (int i = 0; i < 6; ++i) { t0[j] += y[0][i] * bl[0]->single.J[i][j].mid(); t1[j] += y[1][i] * bl[1]->single.J[i][j].mid(); }
        for (int q = 0; q < l; ++q) { t0[j] += ye[q] * DE[q][j]; t1[j] -= ye[q] * DE[q][j]; }
      }
      bool worst = sample % 3 == 0;  // worst-case direction for the linear term (tests the (-1)^t sign convention)
      for (int t = 0; t < 2; ++t) {
        double nn = 0;
        const double* tv = t == 0 ? t0 : t1;
        for (int j = 0; j < 6; ++j) { x[t][j] = worst ? -tv[j] : urand(-1, 1); nn += x[t][j] * x[t][j]; }
        double sc = (worst ? 1.0 : urand(0.2, 1.0)) * s[t] / std::sqrt(nn);
        for (int j = 0; j < 6; ++j) x[t][j] *= sc;
      }
      double phi = 0, cst = 0;
      for (int t = 0; t < 2; ++t) {
        double beta[6];
        for (int j = 0; j < 6; ++j) beta[j] = v[t][j] + x[t][j];
        double lin = 0, nuTerm = 0, Ny2 = 0, c0 = 0, eb = 0;
        for (int i = 0; i < 6; ++i) {
          double sr = 0, si = 0;
          for (int j = 0; j < 6; ++j) {
            double ang = TWO_PI * (beta[j] - Hh[6 * j + i]);
            sr += std::cos(ang) / std::sqrt(6.0);
            si += std::sin(ang) / std::sqrt(6.0);
          }
          eb += y[t][i] * (sr * sr + si * si - 1) / (4 * M_PI);
          c0 += y[t][i] * bl[t]->single.c[i].mid();
          for (int j = 0; j < 6; ++j) {
            double Jij = bl[t]->single.J[i][j].mid();
            lin += y[t][i] * Jij * (x[t][j] - delta[6 * j + i]);
            nuTerm += y[t][i] * Jij * nu[6 * j + i];
            Ny2 += y[t][i] * y[t][i] * Jij * Jij;
          }
        }
        phi += eb;
        cst += c0;
        double rem = std::fabs(eb - c0 - lin);
        double bound = weightedError(C0, bl[t]->single, y[t], Interval(s[t]), Interval(0.0), false).hi;
        ++checked;
        CHECKM(rem <= bound * (1 + 1e-9) + 1e-12, "ball Taylor remainder %g > bound %g (sample %d t %d)", rem, bound, sample, t);
        if (rem > 0.2 * bound) ++tight;
        double nb = (C.nhat * Interval(std::sqrt(Ny2))).hi;  // |sum y_i J_ij nu_ji| <= nhat * N_y
        CHECKM(std::fabs(nuTerm) <= nb * (1 + 1e-9) + 1e-14, "normal part %g > %g", std::fabs(nuTerm), nb);
      }
      double sr = 0, si = 0;
      for (int j = 0; j < 6; ++j) {
        double ang = TWO_PI * ((v[0][j] + x[0][j]) - (v[1][j] + x[1][j]));
        sr += std::cos(ang);
        si += std::sin(ang);
      }
      double Sr0 = e.S.re.mid(), Si0 = e.S.im.mid();
      double z[6];
      for (int j = 0; j < 6; ++j) z[j] = x[0][j] - x[1][j];
      double ePhi, eConst;
      if (isO) {
        double dr = 0, di = 0;
        for (int j = 0; j < 6; ++j) { dr += DE[0][j] * z[j]; di += DE[1][j] * z[j]; }
        double r0 = sr / TWO_PI - (Sr0 / TWO_PI + dr), r1 = si / TWO_PI - (Si0 / TWO_PI + di);
        double rem = std::sqrt(r0 * r0 + r1 * r1);
        double EO = (cn().pi * iv::sqr(e.re)).hi;
        ++checked;
        CHECKM(rem <= EO * (1 + 1e-9) + 1e-14, "orth edge remainder %g > E_O %g", rem, EO);
        ePhi = ye[0] * sr / TWO_PI + ye[1] * si / TWO_PI;
        eConst = ye[0] * Sr0 / TWO_PI + ye[1] * Si0 / TWO_PI;
      } else {
        double dd = 0;
        for (int j = 0; j < 6; ++j) dd += DE[0][j] * z[j];
        double qU = ((Sr0 * Sr0 + Si0 * Si0) / 6 - 1) / (4 * M_PI);
        double aU = ((sr * sr + si * si) / 6 - 1) / (4 * M_PI);
        double rem = std::fabs(aU - (qU + dd));
        double EU = (cn().pi * iv::sqr(e.re) * (Interval(1.0) + e.absS * cn().sixth)).hi;
        ++checked;
        CHECKM(rem <= EU * (1 + 1e-9) + 1e-14, "unbiased edge remainder %g > E_U %g", rem, EU);
        ePhi = ye[0] * aU;
        eConst = ye[0] * qU;
      }
      // joint inequality: |const| <= |Phi| + RHS for every sample (Phi = 0 at an exact configuration)
      double Phi = phi + ePhi;
      double cstT = cst + eConst;
      double ycat[14];
      for (int i = 0; i < 6; ++i) { ycat[i] = y[0][i]; ycat[6 + i] = y[1][i]; }
      ycat[12] = ye[0];
      ycat[13] = ye[1];
      SupportOut so = jointSupport(C, *bl[0], *bl[1], e, isO, ycat);
      ++checked;
      CHECKM(std::fabs(cstT) <= std::fabs(Phi) + so.rhs.hi + 1e-12, "joint inequality violated: |const|=%g > |Phi|=%g + rhs=%g (rel %d)", std::fabs(cstT), std::fabs(Phi), so.rhs.hi, rel);
    }
  }
  std::printf("taylor: %ld remainder/joint checks, %ld ball remainders within 5x of their bound\n", checked, tight);
}

// ------------------------------------------------------------------ near-MUB: {I, F6, X} must not be excluded as a triple
static void test_triple() {
  Record P = fourierRoot(0.0, 0.0, 0.01);
  CParams cp;
  cp.tripleMode = true;
  CStats st;
  COutcome o = processRecordC(P, cp, st);
  std::printf("triple (F6, r=.01): %s, %zu balls, O-edges %zu\n", outcomeName(o), st.balls, st.edgesO);
  CHECK(o == COutcome::TripleFound);  // Grassl: {I,F6} extends to a third basis
  CParams cp2;
  CStats st2;
  COutcome o2 = processRecordC(P, cp2, st2);
  std::printf("two further bases (F6, r=.01): %s (theorem: no 4 MUBs)\n", outcomeName(o2));
  CHECK(isNegative(o2));
  // large radius: capped records must come back UNRESOLVED, never negative-by-accident
  Record Pc = fourierRoot(0.3, 0.5, 0.03);
  CParams cp3;
  cp3.graphCap = 2000;
  CStats st3;
  COutcome o3 = processRecordC(Pc, cp3, st3);
  std::printf("capped run (r=.03, cap 2000): %s\n", outcomeName(o3));
  CHECK(!isNegative(o3));
}

// ------------------------------------------------------------------ closest-to-4-MUB configuration built from exact vectors
// Take the exact vectors unbiased to I and F6 (Grassl's 48), group them into orthonormal bases (O-cliques), and pick the two
// disjoint bases whose cross overlaps are closest to unbiased.  NOTE: an exact (or genuinely near) 4-MUB cannot be built - the best
// pair here has defect ~1 (printed), so the cross U edges are legitimately rejected for the bad pairs.  What is tested is the
// CONSISTENCY chain on real vectors: orthogonal pairs always keep their O edge; a configuration all of whose necessary tests pass
// is found by the graph search and never excluded; the search answers exactly the planted-pattern question.  The genuinely
// "near-MUB" retention test is `triple` (F6 keeps its third basis) together with `single`/`edge` (exact vectors/pairs survive).
static void test_nearmub() {
  Record P0 = fourierRoot(0.0, 0.0, 0.01);
  CRec C;
  initCRec(C, P0);
  std::vector<fx::UVec> V = fx::findUnbiasedVectors(P0.p, 900, rng);
  const int n = (int)V.size();
  std::printf("nearmub: %d exact vectors unbiased to I and F6\n", n);
  CHECK(n >= 12);
  std::vector<std::vector<char>> orth(n, std::vector<char>(n, 0));
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (i != j) orth[i][j] = fx::ipAbs2(V[i].beta, V[j].beta) < 1e-9;
  std::vector<std::vector<int>> bases;
  std::vector<int> cur;
  std::function<void(int)> rec = [&](int s) {
    if (cur.size() == 6) { bases.push_back(cur); return; }
    for (int v = s; v < n; ++v) {
      bool ok = true;
      for (int u : cur) if (!orth[u][v]) { ok = false; break; }
      if (ok) { cur.push_back(v); rec(v + 1); cur.pop_back(); }
    }
  };
  rec(0);
  std::printf("nearmub: %zu orthonormal bases among them\n", bases.size());
  CHECK(bases.size() >= 2);
  double best = 1e9;
  int bx = -1, by = -1;
  for (size_t x = 0; x < bases.size(); ++x)
    for (size_t y = x + 1; y < bases.size(); ++y) {
      bool disjoint = true;
      for (int a : bases[x]) for (int b : bases[y]) if (a == b) disjoint = false;
      if (!disjoint) continue;
      double dmax = 0;
      for (int a : bases[x]) for (int b : bases[y]) dmax = std::max(dmax, std::fabs(fx::ipAbs2(V[a].beta, V[b].beta) - 6.0) / 6.0);
      if (dmax < best) { best = dmax; bx = (int)x; by = (int)y; }
    }
  CHECK(bx >= 0);
  if (bx < 0) return;
  std::printf("nearmub: best pair of disjoint bases has max relative unbiasedness defect %.4f\n", best);
  std::vector<int> ids = bases[bx];
  ids.insert(ids.end(), bases[by].begin(), bases[by].end());
  const double s = 0.03;
  std::vector<std::unique_ptr<Ball>> balls;
  for (int id : ids) {
    double x0[6], mean = 0, nn = 0;
    for (int j = 0; j < 6; ++j) { x0[j] = urand(-1, 1); mean += x0[j]; }
    mean /= 6;
    for (int j = 0; j < 6; ++j) { x0[j] -= mean; nn += x0[j] * x0[j]; }
    double sc = urand(0.3, 0.9) * s / std::sqrt(nn);
    double c[6];
    for (int j = 0; j < 6; ++j) c[j] = V[id].beta[j] - x0[j] * sc;
    balls.emplace_back(new Ball());
    makeBall(C, c, s, *balls.back());
  }
  Graph G;
  G.init(12);
  int needU = 0, haveU = 0, haveO = 0;
  for (int a = 0; a < 12; ++a)
    for (int b = a + 1; b < 12; ++b) {
      EdgeRes er = testEdge(C, *balls[a], *balls[b]);
      bool sameBasis = (a < 6) == (b < 6);
      if (er.O) G.addO(a, b);
      if (er.U) G.addU(a, b);
      if (sameBasis) { CHECK(er.O); haveO += er.O; }  // true orthogonal pairs always keep their O edge
      else { ++needU; haveU += er.U; }
    }
  GraphSearch gs(G, 5000000);
  bool found = gs.exist();
  std::printf("nearmub: O edges kept %d/30, cross U edges kept %d/%d -> graph search %s\n", haveO, haveU, needU, found ? "FINDS the configuration" : "negative");
  CHECK(!gs.gaveUp);
  CHECK(found == (haveU == needU && haveO == 30));  // the search answers exactly the planted-pattern question
  if (haveU == needU) CHECK(found);                 // a configuration that survives all tests is never excluded
}

// ------------------------------------------------------------------ never silently drop: failure categories are UNRESOLVED
static void test_soundness() {
  Record P = fourierRoot(0.3, 0.5, 0.016);
  CParams cp;
  cp.binCap = 50;  // force the overflow path
  CStats st;
  COutcome o = processRecordC(P, cp, st);
  CHECK(o == COutcome::Overflow && !isNegative(o));
  CParams cp2;
  cp2.radiusParam = 0.0005;  // absurdly small cluster radius => > 2048 centres or radius check, never negative-by-accident
  CStats st2;
  COutcome o2 = processRecordC(P, cp2, st2);
  std::printf("soundness: overflow path %s; tiny cluster radius -> %s\n", outcomeName(o), outcomeName(o2));
  // Stage B: a parent that cannot be refined is reported failed (never dropped): huge k grid
  Record Q = fourierRoot(0.3, 0.5, 0.125);
  Refiner rf(0.0001, 0.0001);  // grid limit 50M boxes exceeded (k=4, L ~ 1250 -> 2501^4)
  bool ok = rf.processParent(Q);
  CHECK(!ok);
  CHECK(rf.lastFail.find("50,000,000") != std::string::npos);
  std::printf("soundness: oversize grid reported as failed parent: %s\n", rf.lastFail.c_str());
}

int main(int argc, char** argv) {
  cn();
  iv::RoundUp guard;
  std::set<std::string> want;
  for (int i = 1; i < argc; ++i) want.insert(argv[i]);
  auto run = [&](const char* name, void (*f)()) {
    if (!want.empty() && !want.count(name)) return;
    auto t0 = std::chrono::steady_clock::now();
    int before = failures;
    f();
    std::printf("[%s] %s (%.1fs)\n", name, failures == before ? "PASS" : "FAIL", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  };
  run("io", test_io);
  run("root", test_root);
  run("rel", test_rel);
  run("phi", test_phi);
  run("tensor", test_tensor);
  run("mark", test_mark);
  run("graph", test_graph);
  run("cluster", test_cluster);
  run("single", test_single);
  run("edge", test_edge);
  run("taylor", test_taylor);
  run("nearmub", test_nearmub);
  run("triple", test_triple);
  run("soundness", test_soundness);
  run("stageB", test_stageB);
  if (!want.empty() && want.count("stageBslow")) run("stageBslow", test_stageBslow);
  std::printf("checks: %d, failures: %d\n", checks, failures);
  return failures ? 1 : 0;
}
