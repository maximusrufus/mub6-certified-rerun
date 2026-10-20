// Root certificates (SPEC 3.6; phase-cover.tex phase:root-bounds, phase:stored-root-certificate).
// The eigen-decomposition / Newton snapping are PROPOSALS in plain doubles.  Every constant stored in the
// resulting Record is an outward (interval) upper bound computed from the exact stored doubles P, T, H_V, L_V.
// Differences from the paper's buildH::patch (documented in DECISIONS.md): interval arithmetic replaces all
// additive margins.  C = min(pi*sqrt(8)*mnorm, tensor bound) -- both certified upper bounds (core/tensor_c.hpp).
#pragma once
#include "common.hpp"
#include "../core/tensor_c.hpp"
#include <memory>

namespace mb {

// L_P as a 30x36 interval matrix, ordered pairs (i>j) by increasing i then j, rows (Re, Im).
struct PhaseJac {
  Interval L[30][36];
  CInterval Z[36];  // e(P_ij)
  double Ld[30][36];
};

inline void buildPhaseJac(const double* P, PhaseJac& J) {
  for (int m = 0; m < 36; ++m) J.Z[m] = iv::cis_turns(Interval(P[m]));
  for (int a = 0; a < 30; ++a)
    for (int b = 0; b < 36; ++b) J.L[a][b] = Interval(0.0);
  int m = 0;
  for (int i = 1; i < 6; ++i)
    for (int j = 0; j < i; ++j, ++m)
      for (int g = 0; g < 6; ++g) {
        CInterval c = J.Z[ix(i, g)] * iv::conj(J.Z[ix(j, g)]);  // e(P_ig - P_jg)
        J.L[2 * m][ix(i, g)] = J.L[2 * m][ix(i, g)] + c.re;
        J.L[2 * m][ix(j, g)] = J.L[2 * m][ix(j, g)] - c.re;
        J.L[2 * m + 1][ix(i, g)] = J.L[2 * m + 1][ix(i, g)] + c.im;
        J.L[2 * m + 1][ix(j, g)] = J.L[2 * m + 1][ix(j, g)] - c.im;
      }
  for (int a = 0; a < 30; ++a)
    for (int b = 0; b < 36; ++b) J.Ld[a][b] = J.L[a][b].mid();
}

// F(P): 30 real coordinates (SPEC 3.6), interval enclosure.
inline void phaseResidual(const double* P, Interval* F) {
  const Consts& c = cn();
  CInterval Z[36];
  for (int m = 0; m < 36; ++m) Z[m] = iv::cis_turns(Interval(P[m]));
  Interval inv2pi = Interval(1.0) / c.twopi;
  int m = 0;
  for (int i = 1; i < 6; ++i)
    for (int j = 0; j < i; ++j, ++m) {
      CInterval s = Z[ix(i, 0)] * iv::conj(Z[ix(j, 0)]);
      for (int g = 1; g < 6; ++g) s = s + Z[ix(i, g)] * iv::conj(Z[ix(j, g)]);
      // s/(2 pi i) = (Im s - i Re s)/(2 pi)
      F[2 * m] = s.im * inv2pi;
      F[2 * m + 1] = (-s.re) * inv2pi;
    }
}
inline Interval residualNorm(const double* P) {
  Interval F[30];
  phaseResidual(P, F);
  return inorm(F, 30);
}

// double helper: horizontal orthonormal basis (36x25), Helmert (x) Helmert
inline void horizBasis(double B[36][25]) {
  double H[6][5] = {{0}};
  for (int a = 1; a <= 5; ++a) {
    double nrm = std::sqrt((double)a * (a + 1));
    for (int i = 0; i < a; ++i) H[i][a - 1] = 1.0 / nrm;
    H[a][a - 1] = -(double)a / nrm;
  }
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      for (int a = 0; a < 5; ++a)
        for (int b = 0; b < 5; ++b) B[ix(i, j)][5 * a + b] = H[i][a] * H[j][b];
}

struct HB { double B[36][25]; HB() { horizBasis(B); } };
inline const HB& hb() { static const HB h; return h; }

// Orthonormal basis (m x 36, row = vector) of the horizontal vectors orthogonal to the columns of T.  Heuristic.
inline int complementBasis(const double* T, int k, std::vector<double>& out) {
  const double(*B)[25] = hb().B;
  std::vector<std::vector<double>> acc;
  auto addvec = [&](std::vector<double> v, bool keep) {
    for (int pass = 0; pass < 2; ++pass)
      for (auto& q : acc) {
        double d = 0;
        for (int m = 0; m < 36; ++m) d += q[m] * v[m];
        for (int m = 0; m < 36; ++m) v[m] -= d * q[m];
      }
    double n = 0;
    for (double x : v) n += x * x;
    n = std::sqrt(n);
    if (n < 1e-3) return;
    for (double& x : v) x /= n;
    acc.push_back(v);
    (void)keep;
  };
  for (int a = 0; a < k; ++a) addvec(std::vector<double>(T + 36 * a, T + 36 * a + 36), true);
  size_t nT = acc.size();
  for (int c = 0; c < 25; ++c) {
    std::vector<double> v(36);
    for (int m = 0; m < 36; ++m) v[m] = B[m][c];
    addvec(v, false);
  }
  out.clear();
  for (size_t i = nT; i < acc.size(); ++i) out.insert(out.end(), acc[i].begin(), acc[i].end());
  return (int)(acc.size() - nT);
}

// double-precision residual / Jacobian used ONLY by the snapping heuristic (no rigor needed: certifyRoot re-verifies).
inline double residD(const double* X, double* F, double* cr = nullptr, double* ci = nullptr) {
  const double TWO_PI = 6.283185307179586476925286766559;
  double zr[36], zi[36];
  for (int m = 0; m < 36; ++m) { zr[m] = std::cos(TWO_PI * X[m]); zi[m] = std::sin(TWO_PI * X[m]); }
  if (cr) { std::memcpy(cr, zr, sizeof zr); std::memcpy(ci, zi, sizeof zi); }
  double n2 = 0;
  int m = 0;
  for (int i = 1; i < 6; ++i)
    for (int j = 0; j < i; ++j, ++m) {
      double sr = 0, si = 0;
      for (int g = 0; g < 6; ++g) {
        sr += zr[6 * i + g] * zr[6 * j + g] + zi[6 * i + g] * zi[6 * j + g];
        si += zi[6 * i + g] * zr[6 * j + g] - zr[6 * i + g] * zi[6 * j + g];
      }
      F[2 * m] = si / TWO_PI;
      F[2 * m + 1] = -sr / TWO_PI;
      n2 += F[2 * m] * F[2 * m] + F[2 * m + 1] * F[2 * m + 1];
    }
  return std::sqrt(n2);
}

// Newton snapping onto F(P)=0 with minimum-norm Gauss-Newton steps inside the given basis (rows of `basis`, nb vectors
// of length 36); basis==nullptr uses all 25 horizontal directions.  Heuristic only.
inline bool snapSolution(double* P, const std::vector<double>* basis = nullptr, int maxit = 60) {
  std::vector<double> Bm;  // nb x 36
  int nb;
  if (basis) { Bm = *basis; nb = (int)(Bm.size() / 36); }
  else {
    const double(*B)[25] = hb().B;
    nb = 25;
    Bm.assign(25 * 36, 0.0);
    for (int c = 0; c < 25; ++c) for (int m = 0; m < 36; ++m) Bm[c * 36 + m] = B[m][c];
  }
  if (nb == 0) { double F[30]; return residD(P, F) < 1e-13; }
  for (int it = 0; it < maxit; ++it) {
    double F[30], zr[36], zi[36];
    double nf = residD(P, F, zr, zi);
    if (!std::isfinite(nf)) return false;
    if (nf < 3e-16) return true;
    // Jacobian rows (Re, Im) of L_P in the 36 phase variables
    static thread_local double L[30][36];
    std::memset(L, 0, sizeof L);
    int m = 0;
    for (int i = 1; i < 6; ++i)
      for (int j = 0; j < i; ++j, ++m)
        for (int g = 0; g < 6; ++g) {
          double cre = zr[6 * i + g] * zr[6 * j + g] + zi[6 * i + g] * zi[6 * j + g];
          double cim = zi[6 * i + g] * zr[6 * j + g] - zr[6 * i + g] * zi[6 * j + g];
          L[2 * m][6 * i + g] += cre;  L[2 * m][6 * j + g] -= cre;
          L[2 * m + 1][6 * i + g] += cim;  L[2 * m + 1][6 * j + g] -= cim;
        }
    std::vector<double> Jb((size_t)30 * nb, 0.0), G((size_t)nb * nb, 0.0), rhs(nb, 0.0);
    for (int a = 0; a < 30; ++a)
      for (int c = 0; c < nb; ++c) {
        double s = 0;
        for (int q = 0; q < 36; ++q) s += L[a][q] * Bm[c * 36 + q];
        Jb[(size_t)a * nb + c] = s;
      }
    for (int c = 0; c < nb; ++c)
      for (int d = c; d < nb; ++d) {
        double s = 0;
        for (int a = 0; a < 30; ++a) s += Jb[(size_t)a * nb + c] * Jb[(size_t)a * nb + d];
        G[(size_t)c * nb + d] = s;
        G[(size_t)d * nb + c] = s;
      }
    for (int c = 0; c < nb; ++c) {
      double s = 0;
      for (int a = 0; a < 30; ++a) s += Jb[(size_t)a * nb + c] * F[a];
      rhs[c] = -s;
    }
    std::vector<double> step;
    {
      std::vector<double> Gr = G;
      double tr = 0;
      for (int c = 0; c < nb; ++c) tr += G[(size_t)c * nb + c];
      for (int c = 0; c < nb; ++c) Gr[(size_t)c * nb + c] += 1e-12 * tr;  // tiny ridge: solves stay finite if the basis is not a clean complement
      if (!solveLin(Gr, rhs, nb, step)) return false;
    }
    double Pn[36];
    double damp = 1.0;
    bool improved = false;
    for (int tr = 0; tr < 12; ++tr) {
      for (int q = 0; q < 36; ++q) {
        double s = 0;
        for (int c = 0; c < nb; ++c) s += Bm[c * 36 + q] * step[c];
        Pn[q] = P[q] + damp * s;
      }
      double Fn[30];
      if (residD(Pn, Fn) < nf) { improved = true; break; }
      damp *= 0.5;
    }
    if (!improved) return nf < 1e-13;
    std::memcpy(P, Pn, sizeof Pn);
  }
  double F[30];
  return residD(P, F) < 1e-13;
}

struct RootInfo {
  std::string why;
  double minEig = 0;
  double Cplain = 0, Ctensor = 0;  // the two candidate remainder constants (diagnostics)
};
// Tensor improvement of C (core/tensor_c.hpp): on by default; MUB_NO_TENSOR_C=1 reverts to C = pi sqrt8 mnorm
// (measurement switch only: both are valid upper bounds).
inline bool useTensorC() { static const bool v = std::getenv("MUB_NO_TENSOR_C") == nullptr; return v; }

// Build and certify a Root centred at Pin (radius left 0).  Returns false (with reason) on any failed check.
inline bool certifyRoot(const double* Pin, Record& out, RootInfo* info = nullptr) {
  auto fail = [&](const char* s) { if (info) info->why = s; return false; };
  const Consts& cc = cn();
  Record R;
  for (int m = 0; m < 36; ++m) {
    if (!(std::fabs(Pin[m]) <= 4)) return fail("|P| > 4");
    R.p[m] = Pin[m];
  }
  const double(*B)[25] = hb().B;
  PhaseJac* J = new PhaseJac();
  struct Del { PhaseJac* p; ~Del() { delete p; } } del{J};
  buildPhaseJac(R.p, *J);

  // J_b = L B (30x25), Gram, eigen
  std::vector<double> Jb(30 * 25, 0.0), G(25 * 25, 0.0);
  for (int a = 0; a < 30; ++a)
    for (int c = 0; c < 25; ++c) {
      double s = 0;
      for (int m = 0; m < 36; ++m) s += J->Ld[a][m] * B[m][c];
      Jb[a * 25 + c] = s;
    }
  for (int c = 0; c < 25; ++c)
    for (int d = 0; d < 25; ++d) {
      double s = 0;
      for (int a = 0; a < 30; ++a) s += Jb[a * 25 + c] * Jb[a * 25 + d];
      G[c * 25 + d] = s;
    }
  std::vector<double> V, ev;
  jacobiEig(G, 25, V, ev);
  int k = 0;
  while (k < 25 && ev[k] < 0.64) ++k;
  if (k > 10) return fail("more than 10 tangent columns");
  int g = 25 - k;
  if (info) info->minEig = ev[0];
  // T columns = B V_c (c < k); H_V, L_V columns for c >= k
  R.k = k;
  R.T.assign((size_t)k * 36, 0.0);
  for (int a = 0; a < k; ++a)
    for (int m = 0; m < 36; ++m) {
      double s = 0;
      for (int c = 0; c < 25; ++c) s += B[m][c] * V[c * 25 + a];
      R.T[(size_t)a * 36 + m] = s;
    }
  std::vector<double> HV((size_t)36 * g), LV((size_t)30 * g);  // row-major [row][col]
  for (int a = 0; a < g; ++a) {
    int e = k + a;
    double sc = 1.0 / std::sqrt(ev[e]);
    for (int m = 0; m < 36; ++m) {
      double s = 0;
      for (int c = 0; c < 25; ++c) s += B[m][c] * V[c * 25 + e];
      HV[(size_t)m * g + a] = s * sc;
    }
    for (int r = 0; r < 30; ++r) {
      double s = 0;
      for (int c = 0; c < 25; ++c) s += Jb[r * 25 + c] * V[c * 25 + e];
      LV[(size_t)r * g + a] = s * sc;
    }
  }
  for (double x : HV) if (!(std::fabs(x) <= 2)) return fail("H_V entry > 2");
  for (double x : LV) if (!(std::fabs(x) <= 2)) return fail("L_V entry > 2");
  for (double x : R.T) if (!(std::fabs(x) <= 1.01)) return fail("T entry > 1.01");

  // ---- rigorous constants
  // ort = ||T^T T - I||_F
  Interval ortSq(0.0);
  for (int a = 0; a < k; ++a)
    for (int b = 0; b < k; ++b) {
      Interval s(0.0);
      for (int m = 0; m < 36; ++m) s = s + Interval(R.T[(size_t)a * 36 + m]) * Interval(R.T[(size_t)b * 36 + m]);
      if (a == b) s = s - Interval(1.0);
      ortSq = ortSq + iv::sqr(s);
    }
  double ort = iv::sqrt(ortSq).hi;
  // hor = ||(I-Pi)T||_F
  Interval horSq(0.0);
  for (int a = 0; a < k; ++a) {
    Interval t[36], pt[36];
    for (int m = 0; m < 36; ++m) t[m] = Interval(R.T[(size_t)a * 36 + m]);
    PiI(t, pt);
    for (int m = 0; m < 36; ++m) horSq = horSq + iv::sqr(t[m] - pt[m]);
  }
  double hor = iv::sqrt(horSq).hi;
  // froot
  double froot = residualNorm(R.p).hi;
  // mnorm from absolute row sums of the Gram matrices of H_V, L_V
  auto gramRowSum = [&](const std::vector<double>& A, int rows) {
    double best = 0;
    for (int a = 0; a < g; ++a) {
      Interval rs(0.0);
      for (int b = 0; b < g; ++b) {
        Interval s(0.0);
        for (int r = 0; r < rows; ++r) s = s + Interval(A[(size_t)r * g + a]) * Interval(A[(size_t)r * g + b]);
        rs = rs + iv::abs(s);
      }
      best = std::max(best, rs.hi);
    }
    return best;
  };
  double ah = gramRowSum(HV, 36), al = gramRowSum(LV, 30);
  double mnorm = iv::sqrt(Interval(ah) * Interval(al)).hi;
  double n0 = (Interval(mnorm) * Interval(froot)).hi;
  double Cplain = (cc.pi * cc.sqrt8 * Interval(mnorm)).hi;
  double Ctensor = tc::tensorC(HV, LV, g, cc.pi);
  double C = useTensorC() ? std::min(Cplain, Ctensor) : Cplain;
  if (info) { info->Cplain = Cplain; info->Ctensor = Ctensor; }

  // id = || Pi - T T^T - H_V (L_V^T L_P) Pi ||_F, full 36x36
  // W = L_V^T L_P (g x 36), then Pi applied to each row, then H_V * (W Pi)
  std::vector<Interval> WP((size_t)g * 36);
  for (int a = 0; a < g; ++a) {
    Interval w[36], wp[36];
    for (int m = 0; m < 36; ++m) {
      Interval s(0.0);
      for (int r = 0; r < 30; ++r) s = s + Interval(LV[(size_t)r * g + a]) * J->L[r][m];
      w[m] = s;
    }
    PiI(w, wp);  // (W Pi) row = Pi applied to the row (Pi symmetric)
    for (int m = 0; m < 36; ++m) WP[(size_t)a * 36 + m] = wp[m];
  }
  Interval idSq(0.0);
  for (int n = 0; n < 36; ++n) {
    Interval e[36], pe[36];
    for (int m = 0; m < 36; ++m) e[m] = Interval(m == n ? 1.0 : 0.0);
    PiI(e, pe);  // column n of Pi (Pi symmetric)
    for (int m = 0; m < 36; ++m) {
      Interval tt(0.0);
      for (int a = 0; a < k; ++a) tt = tt + Interval(R.T[(size_t)a * 36 + m]) * Interval(R.T[(size_t)a * 36 + n]);
      Interval ml(0.0);
      for (int a = 0; a < g; ++a) ml = ml + Interval(HV[(size_t)m * g + a]) * WP[(size_t)a * 36 + n];
      idSq = idSq + iv::sqr(pe[m] - tt - ml);
    }
  }
  double id = iv::sqrt(idSq).hi;

  if (!(id < 1e-4)) return fail("id >= 1e-4");
  if (!(ort < 1e-5)) return fail("ort >= 1e-5");
  if (!(hor < 1e-5)) return fail("hor >= 1e-5");
  if (!(C < 20)) return fail("C >= 20");
  if (!(mnorm < 2)) return fail("mnorm >= 2");
  if (!(froot < 1e-9)) return fail("froot >= 1e-9 (not a snapped solution)");
  if (!(n0 < 1e-6)) return fail("n0 >= 1e-6");
  R.C = C; R.n0 = n0; R.id = id; R.ort = ort; R.hor = hor; R.mnorm = mnorm; R.froot = froot; R.r = 0;
  out = std::move(R);
  return true;
}

}  // namespace mb
