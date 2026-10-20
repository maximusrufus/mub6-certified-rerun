// reuse.hpp -- certified center reuse (seed.tex sec 2.4: testBucket / supdual)
// with interval certificates.  Heuristics (the proposal of v) run in plain
// doubles; the certificate `dual_cert_hi` is evaluated with enclosures and an
// outward-rounded upper bound, and success is a strict comparison against the
// lower bound of (R N)^2.  No additive margins, no 1.251 inflation.
// Call under an iv::RoundUp guard.
//
// Nine-coordinate order: position 0 = shared entry xi = x_{s+1}; 1..4 = the
// other four x entries in their former order; 5..8 = y_2..y_5.
#pragma once
#include <cmath>
#include <cstdint>
#include "cover.hpp"

namespace seedcov {

struct PairCtx {
  int s;
  int xb[9], yb[9];          // bin of the X / Y list entry at each nine position (-1 if none)
  Interval p[4];             // p_* = (-N c_X ; -N c_Y), enclosures
  // heuristic (double) copies
  double Jt[9][4];
  double pd[4];
  double M[4][4];            // Cholesky factor of Jt^T Jt + ridge
  bool hasM;
};

inline void build_ctx(const Tables &T, const ListData &dx, const ListData &dy, int s, PairCtx &c) {
  c.s = s;
  for (int j = 0; j < 9; ++j) { c.xb[j] = -1; c.yb[j] = -1; for (int a = 0; a < 4; ++a) c.Jt[j][a] = 0.0; }
  int pos = 1;
  for (int m = 0; m < 5; ++m) {
    int q = (m == s) ? 0 : pos++;
    c.xb[q] = dx.b[m];
  }
  c.yb[0] = dy.b[0];
  for (int m = 1; m < 5; ++m) c.yb[4 + m] = dy.b[m];
  for (int j = 0; j < 9; ++j) {
    if (c.xb[j] >= 0) { c.Jt[j][0] = T.ec[c.xb[j]].mid(); c.Jt[j][1] = T.es[c.xb[j]].mid(); }
    if (c.yb[j] >= 0) { c.Jt[j][2] = T.ec[c.yb[j]].mid(); c.Jt[j][3] = T.es[c.yb[j]].mid(); }
  }
  Interval Nn((double)T.N);
  c.p[0] = -(Nn * dx.c0); c.p[1] = -(Nn * dx.c1);
  c.p[2] = -(Nn * dy.c0); c.p[3] = -(Nn * dy.c1);
  for (int a = 0; a < 4; ++a) c.pd[a] = c.p[a].mid();
  double A[4][4];
  for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) {
    double sacc = 0.0;
    for (int j = 0; j < 9; ++j) sacc += c.Jt[j][a] * c.Jt[j][b];
    A[a][b] = sacc + (a == b ? 1e-9 : 0.0);
  }
  c.hasM = true;
  for (int i = 0; i < 4 && c.hasM; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sacc = A[i][j];
      for (int m = 0; m < j; ++m) sacc -= c.M[i][m] * c.M[j][m];
      if (i == j) {
        if (!(sacc > 1e-12)) { c.hasM = false; break; }
        c.M[i][i] = std::sqrt(sacc);
      } else {
        c.M[i][j] = sacc / c.M[j][j];
      }
    }
  }
}

// CERTIFIED upper bound (seed:dual-formula) of  2 k.h  over the feasible set:
//   sum_j |k_j - (J^T v)_j| + 2 v.p + 2 E (|v_X| + |v_Y|)       with E = 5 pi/(4N).
// v is any finite exact double vector.  Returns +inf for unusable v.
inline double dual_cert_hi(const Tables &T, const PairCtx &c, const int *k, const double *v) {
  for (int a = 0; a < 4; ++a) if (!std::isfinite(v[a])) return INFINITY;
  double sacc = 0.0;
  for (int j = 0; j < 9; ++j) {
    Interval g(0.0);
    if (c.xb[j] >= 0) g = g + mulp(T.ec[c.xb[j]], v[0]) + mulp(T.es[c.xb[j]], v[1]);
    if (c.yb[j] >= 0) g = g + mulp(T.ec[c.yb[j]], v[2]) + mulp(T.es[c.yb[j]], v[3]);
    Interval d = Interval((double)k[j]) - g;
    sacc += abs_hi(d);       // up
  }
  Interval lin = mulp(c.p[0], v[0]) + mulp(c.p[1], v[1]) + mulp(c.p[2], v[2]) + mulp(c.p[3], v[3]);
  double nx = std::sqrt(v[0] * v[0] + v[1] * v[1]);   // up
  double ny = std::sqrt(v[2] * v[2] + v[3] * v[3]);   // up
  double tot = sacc + 2.0 * lin.hi + 2.0 * (T.Eg_hi * (nx + ny));
  // lin.hi may be negative: the sum is then an upper bound of a sum with a negative term, still valid
  // under upward rounding of every partial result.
  return tot;
}

// ---- heuristic (double) objective used only to PROPOSE v -----------------------
inline double dual_est(const Tables &T, const PairCtx &c, const int *k, const double *v) {
  double sacc = 0.0;
  for (int j = 0; j < 9; ++j) {
    double g = c.Jt[j][0] * v[0] + c.Jt[j][1] * v[1] + c.Jt[j][2] * v[2] + c.Jt[j][3] * v[3];
    sacc += std::fabs((double)k[j] - g);
  }
  double lin = 0.0;
  for (int a = 0; a < 4; ++a) lin += c.pd[a] * v[a];
  double Ed = 1.25 * arith::PI / (double)T.N;
  return sacc + 2.0 * lin + 2.0 * Ed * (std::sqrt(v[0] * v[0] + v[1] * v[1]) + std::sqrt(v[2] * v[2] + v[3] * v[3]));
}

static inline void chol_solve4(const double L[4][4], const double *b, double *x) {
  double y[4];
  for (int i = 0; i < 4; ++i) {
    double sacc = b[i];
    for (int m = 0; m < i; ++m) sacc -= L[i][m] * y[m];
    y[i] = sacc / L[i][i];
  }
  for (int i = 3; i >= 0; --i) {
    double sacc = y[i];
    for (int m = i + 1; m < 4; ++m) sacc -= L[m][i] * x[m];
    x[i] = sacc / L[i][i];
  }
}

// Heuristic proposal of v (SPEC A4): least squares start, then a few damped Newton
// steps on a smoothed dual objective.  Returns the best estimated bound; vbest = proposal.
inline double propose_v(const Tables &T, const PairCtx &c, const int *k, double *vbest) {
  double best = 1e300;
  for (int a = 0; a < 4; ++a) vbest[a] = 0.0;
  if (!c.hasM) return best;
  double v[4];
  {
    double rhs[4];
    for (int a = 0; a < 4; ++a) {
      double sacc = 0.0;
      for (int j = 0; j < 9; ++j) sacc += c.Jt[j][a] * (double)k[j];
      rhs[a] = sacc;
    }
    chol_solve4(c.M, rhs, v);
  }
  double cur = dual_est(T, c, k, v);
  if (cur < best) { best = cur; for (int a = 0; a < 4; ++a) vbest[a] = v[a]; }
  const double Ed = 1.25 * arith::PI / (double)T.N;
  static const double EPS[3] = {0.3, 0.05, 0.01};
  for (int e = 0; e < 3; ++e) {
    double eps = EPS[e], eps2 = eps * eps;
    for (int it = 0; it < 3; ++it) {
      double g[4] = {2.0 * c.pd[0], 2.0 * c.pd[1], 2.0 * c.pd[2], 2.0 * c.pd[3]};
      double H[4][4] = {{0}};
      for (int j = 0; j < 9; ++j) {
        double gj = c.Jt[j][0] * v[0] + c.Jt[j][1] * v[1] + c.Jt[j][2] * v[2] + c.Jt[j][3] * v[3];
        double r = (double)k[j] - gj;
        double rho = std::sqrt(r * r + eps2);
        double wg = r / rho;
        double wh = eps2 / (rho * rho * rho);
        for (int a = 0; a < 4; ++a) {
          g[a] -= wg * c.Jt[j][a];
          for (int b = 0; b < 4; ++b) H[a][b] += wh * c.Jt[j][a] * c.Jt[j][b];
        }
      }
      double sx = std::sqrt(v[0] * v[0] + v[1] * v[1] + eps2);
      double sy = std::sqrt(v[2] * v[2] + v[3] * v[3] + eps2);
      g[0] += 2.0 * Ed * v[0] / sx; g[1] += 2.0 * Ed * v[1] / sx;
      g[2] += 2.0 * Ed * v[2] / sy; g[3] += 2.0 * Ed * v[3] / sy;
      double hx = 2.0 * Ed * eps2 / (sx * sx * sx), hy = 2.0 * Ed * eps2 / (sy * sy * sy);
      H[0][0] += hx; H[1][1] += hx; H[2][2] += hy; H[3][3] += hy;
      double L[4][4]; bool ok = true;
      for (int i = 0; i < 4 && ok; ++i) for (int j = 0; j <= i; ++j) {
        double sacc = H[i][j] + (i == j ? 1e-9 : 0.0);
        for (int m = 0; m < j; ++m) sacc -= L[i][m] * L[j][m];
        if (i == j) { if (!(sacc > 1e-14)) { ok = false; break; } L[i][i] = std::sqrt(sacc); }
        else L[i][j] = sacc / L[j][j];
      }
      if (!ok) break;
      double rhs[4] = {-g[0], -g[1], -g[2], -g[3]}, d[4];
      chol_solve4(L, rhs, d);
      double step = 1.0, fnew = cur, vn[4];
      bool improved = false;
      for (int ls = 0; ls < 4; ++ls) {
        for (int a = 0; a < 4; ++a) vn[a] = v[a] + step * d[a];
        fnew = dual_est(T, c, k, vn);
        if (fnew < cur) { improved = true; break; }
        step *= 0.5;
      }
      if (!improved) break;
      for (int a = 0; a < 4; ++a) v[a] = vn[a];
      cur = fnew;
      if (cur < best) { best = cur; for (int a = 0; a < 4; ++a) vbest[a] = v[a]; }
    }
  }
  return best;
}

// testBucket for one previous center.  kn2 = ||k||^2 (exact integer), kl1 = sum|k_j|.
// Success iff the CERTIFIED upper bound of ||k+h||^2 is strictly below the lower bound of
// (R N)^2.  Box bound first (exact: kn2 + 9/4 + kl1), then the dual bound.
// `ctx` is built lazily through build() (called at most once per pair) by the caller.
inline bool box_success(const Tables &T, int kn2, int kl1) {
  return (double)kn2 + 2.25 + (double)kl1 < T.RN2lo;     // all terms exactly representable
}
inline bool dual_success(const Tables &T, const PairCtx &c, const int *k, int kn2) {
  double v[4];
  double est = propose_v(T, c, k, v);
  // pre-filter on the heuristic estimate only (saves interval work); NOT a certificate
  if (!(est + (double)kn2 + 2.25 < T.RN2lo + 0.5)) return false;
  double b = dual_cert_hi(T, c, k, v);
  double tot = ((double)kn2 + 2.25) + b;      // up
  return tot < T.RN2lo;
}

}  // namespace seedcov
