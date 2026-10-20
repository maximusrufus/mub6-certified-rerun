// Tests of the normal forms: seed:equivalences / seed:arrangement /
// phase:alignment / phase:transport.
#include <random>
#include "check.hpp"
#include "hadamard.hpp"

using namespace nf;

static double circ_gap_of(const double *v6) {
  double s[6]; for (int i = 0; i < 6; ++i) s[i] = frac01(v6[i]);
  std::sort(s, s + 6);
  double g = 0;
  for (int i = 0; i < 6; ++i) { double d = (i == 5 ? s[0] + 1.0 : s[i + 1]) - s[i]; g = std::max(g, d); }
  return g;
}

static void test_projector(std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> U(-3, 3);
  for (int t = 0; t < 50; ++t) {
    double z[36], p[36], pp[36];
    for (double &v : z) v = U(rng);
    Pi(z, p); Pi(p, pp);
    double rs, cs, nz = 0, np = 0, idem = 0, orth = 0;
    for (int i = 0; i < 6; ++i) { rs = cs = 0; for (int j = 0; j < 6; ++j) { rs += p[6 * i + j]; cs += p[6 * j + i]; } CHECK(std::fabs(rs) < 1e-12); CHECK(std::fabs(cs) < 1e-12); }
    for (int m = 0; m < 36; ++m) { nz += z[m] * z[m]; np += p[m] * p[m]; idem = std::max(idem, std::fabs(p[m] - pp[m])); orth += (z[m] - p[m]) * p[m]; }
    CHECK(idem < 1e-12);
    CHECK(np <= nz + 1e-12);        // seed:horizontal norm bound
    CHECK(std::fabs(orth) < 1e-10); // z - Pi z orthogonal to range
  }
}

static void test_alignment(std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> U(-2, 2);
  for (int t = 0; t < 200; ++t) {
    double P[36], Q[36], Nint[36], d[36], dl[36], Qd[36], Pl[36], Ql[36], hz[36];
    for (int m = 0; m < 36; ++m) { P[m] = U(rng); Q[m] = U(rng); }
    // arbitrary integer choice (not only the rounded one)
    for (int m = 0; m < 36; ++m) Nint[m] = (t % 2) ? std::nearbyint(U(rng) * 3) : 0.0;
    if (t % 3 == 0) alignment_integers(P, Q, Nint);
    alignment_delta(P, Q, Nint, d);
    // d is horizontal
    double hs[36]; Pi(d, hs); for (int m = 0; m < 36; ++m) CHECK(std::fabs(hs[m] - d[m]) < 1e-12);
    for (int m = 0; m < 36; ++m) Qd[m] = Q[m] + d[m];
    CHECK(rephasing_equivalent(Qd, P, 1e-9));            // phase:alignment
    // and for every horizontal delta
    for (int m = 0; m < 36; ++m) dl[m] = U(rng);
    Pi(dl, hz);
    for (int m = 0; m < 36; ++m) { Pl[m] = P[m] + hz[m]; Ql[m] = Q[m] + d[m] + hz[m]; }
    CHECK(rephasing_equivalent(Ql, Pl, 1e-9));
  }
}

static Sym random_sym(std::mt19937_64 &rng) {
  Sym S;
  S.transpose = rng() & 1; S.flip = rng() & 1;
  for (int i = 5; i > 0; --i) { std::swap(S.rp[i], S.rp[rng() % (i + 1)]); std::swap(S.cp[i], S.cp[rng() % (i + 1)]); }
  return S;
}

static void test_transport(std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> U(-1, 1);
  for (int t = 0; t < 200; ++t) {
    double P[36], Pp[36], z[36], hz[36], phz[36], Phz[36], hphz[36];
    for (double &v : P) v = U(rng);
    Sym S = random_sym(rng);
    apply_sym(S, P, Pp);
    // residual invariance: ||F(Phi P)|| = ||F(P)|| and = Gram identity
    double f0 = normF(P), f1 = normF(Pp);
    CHECKM(std::fabs(f0 - f1) < 1e-10 * (1 + f0), "%g %g", f0, f1);
    double g = normF2_via_gram(P);
    CHECKM(std::fabs(g - f0 * f0) < 1e-9 * (1 + f0 * f0), "%g %g", g, f0 * f0);
    // Phi is a linear isometry commuting with Pi
    for (double &v : z) v = U(rng);
    Pi(z, hz);
    apply_sym(S, z, phz);
    Pi(phz, hphz);
    apply_sym(S, hz, Phz);
    double nz = 0, nphz = 0;
    for (int m = 0; m < 36; ++m) { CHECK(std::fabs(hphz[m] - Phz[m]) < 1e-12); nz += z[m] * z[m]; nphz += phz[m] * phz[m]; }
    CHECK(std::fabs(nz - nphz) < 1e-10);
    // tangent frame transport keeps column norms / inner products
    double T[3 * 36], T2[3 * 36];
    for (double &v : T) v = U(rng);
    transport_T(S, T, 3, T2);
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
      double i1 = 0, i2 = 0; for (int m = 0; m < 36; ++m) { i1 += T[36 * a + m] * T[36 * b + m]; i2 += T2[36 * a + m] * T2[36 * b + m]; }
      CHECK(std::fabs(i1 - i2) < 1e-10);
    }
  }
}

// one arrangement result, all lemma conditions
static void check_arrangement(const double *P0, int N, int pairIdx, int gapTie, bool binCeil, const char *tag) {
  NormalForm nfm;
  bool ok = arrange(P0, N, pairIdx, gapTie, binCeil, nfm);
  CHECKM(ok, "%s arrange failed", tag);
  if (!ok) return;
  const double tol = 1e-8;
  // E is the dephased permuted (maybe transposed) original
  double Pt[36], Pp[36], D[36];
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Pt[6 * i + j] = nfm.transposed ? P0[6 * j + i] : P0[6 * i + j];
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Pp[6 * i + j] = Pt[6 * nfm.rowperm[i] + nfm.colperm[j]];
  D0(Pp, D);
  for (int m = 0; m < 36; ++m) { double d = D[m] - nfm.E[m]; CHECKM(std::fabs(d - std::nearbyint(d)) < 1e-9, "%s equivalence", tag); }
  CHECKM(hd::unitarity_error(nfm.E) < 1e-8, "%s unitarity %g", tag, hd::unitarity_error(nfm.E));
  // first row/col zero
  for (int j = 0; j < 6; ++j) { CHECK(nfm.E[j] == 0.0); CHECK(nfm.E[6 * j] == 0.0); }
  // x arrangement
  CHECK(nfm.x[0] == 0.0);
  for (int m = 1; m < 6; ++m) CHECKM(nfm.x[m] >= nfm.x[m - 1], "%s x sorted", tag);
  CHECKM(nfm.x[1] <= nfm.x[5] - nfm.x[4] + tol, "%s orientation x1=%g x5-x4=%g", tag, nfm.x[1], nfm.x[5] - nfm.x[4]);
  double mg = min_gap_pairs(P0, nullptr);
  CHECKM(std::fabs(nfm.gx - mg) < tol, "%s gx %g vs min gap %g", tag, nfm.gx, mg);
  CHECKM(std::fabs(circ_gap_of(nfm.x) - nfm.gx) < tol, "%s gx is the gap of x", tag);
  // seed:local-conditions
  double sr = 1, si = 0, tr = 1, ti = 0;
  for (int j = 1; j < 6; ++j) { sr += std::cos(2 * M_PI * nfm.x[j]); si += std::sin(2 * M_PI * nfm.x[j]); tr += std::cos(2 * M_PI * nfm.y[j]); ti += std::sin(2 * M_PI * nfm.y[j]); }
  CHECKM(std::hypot(sr, si) < tol, "%s sum x %g", tag, std::hypot(sr, si));
  CHECKM(std::hypot(tr, ti) < tol, "%s sum y %g", tag, std::hypot(tr, ti));
  CHECKM(circ_gap_of(nfm.y) >= nfm.gx - tol, "%s g_y >= g_x", tag);
  // y arrangement
  CHECK(nfm.y[0] == 0.0);
  CHECK(std::fabs(nfm.y[1] - nfm.x[nfm.s + 1]) < 1e-9);
  for (int m = 3; m < 6; ++m) CHECKM(nfm.y[m] >= nfm.y[m - 1], "%s y sorted", tag);
  // selection rule on the recorded bins
  int s = 0, bv = std::abs(2 * nfm.bx[0] + 1 - N);
  for (int m = 1; m < 5; ++m) { int v = std::abs(2 * nfm.bx[m] + 1 - N); if (v < bv) { bv = v; s = m; } }
  CHECK(s == nfm.s);
  CHECK(nfm.by[1] == nfm.bx[nfm.s]);
  // bins really contain the phases (closed bins)
  for (int m = 0; m < 5; ++m) { double lo = (double)nfm.bx[m] / N - 1e-9, hi = (double)(nfm.bx[m] + 1) / N + 1e-9; CHECK(nfm.x[m + 1] >= lo && nfm.x[m + 1] <= hi); }
  for (int m = 2; m < 6; ++m) { double lo = (double)nfm.by[m] / N - 1e-9, hi = (double)(nfm.by[m] + 1) / N + 1e-9; CHECK(nfm.y[m] >= lo && nfm.y[m] <= hi); }
}

int main() {
  std::mt19937_64 rng(20261019);
  test_projector(rng);
  test_alignment(rng);
  test_transport(rng);

  // Fourier family: unitarity (validates the formula) + arrangement on a grid
  int nF = 0;
  for (int ia = 0; ia < 12; ++ia) for (int ib = 0; ib < 12; ++ib) {
    double P[36]; hd::fourier_family(ia / 12.0 + 0.013, ib / 12.0 + 0.007 * ia, P);
    CHECKM(hd::unitarity_error(P) < 1e-12, "F6(a,b) unitarity %g", hd::unitarity_error(P));
    for (int pc = 0; pc < 4; ++pc) for (int gt = 0; gt < 2; ++gt) check_arrangement(P, 112, pc, gt, (pc + gt) & 1, "F6");
    ++nF;
  }
  // F6 itself (a=b=0) has many tied pairs and boundary values
  {
    double P[36]; hd::fourier_family(0, 0, P);
    for (int pc = 0; pc < 30; ++pc) for (int gt = 0; gt < 3; ++gt) for (int bc = 0; bc < 2; ++bc) check_arrangement(P, 112, pc, gt, bc, "F6exact");
    for (int N : {56, 64, 100}) for (int pc = 0; pc < 30; ++pc) check_arrangement(P, N, pc, 0, pc & 1, "F6exactN");
  }
  // random exact Hadamards from LM
  int got = 0, tries = 0;
  while (got < 150 && tries < 6000) {
    ++tries;
    double P[36];
    if (!hd::random_hadamard(rng, P)) continue;
    CHECKM(hd::unitarity_error(P) < 1e-10, "LM unitarity");
    ++got;
    for (int pc = 0; pc < 3; ++pc) check_arrangement(P, 112, pc, pc, pc & 1, "LM");
  }
  printf("random exact Hadamards: %d of %d LM starts converged\n", got, tries);
  CHECK(got >= 100);
  return finish("test_normal");
}
