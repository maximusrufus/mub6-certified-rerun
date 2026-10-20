// tests/hadamard.hpp -- sources of exact complex Hadamard matrices (phase form,
// cycles) for tests: the Fourier family F6(a,b) and random Levenberg-Marquardt
// solutions of F(P)=0.  Test-only: nothing certified depends on this.
#pragma once
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>
#include "../arith.hpp"
#include "../normal.hpp"

namespace hd {

// Fourier family F6^(2)(a,b) = F6 o exp(i R(a,b)), R = [0;0 a b 0 a b;0;0 a b 0 a b;0;0 a b 0 a b]
// (a,b in cycles).  Verified unitary numerically in tests (unitarity_error).
inline void fourier_family(double a, double b, double *P) {
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
    double ph = (double)((i * j) % 6) / 6.0;
    if (i % 2 == 1) {
      if (j == 1 || j == 4) ph += a;
      if (j == 2 || j == 5) ph += b;
    }
    P[6 * i + j] = ph;
  }
}

// max |(H H*)_{ij} - 6 delta_ij| for H = e(P)
inline double unitarity_error(const double *P) {
  double Ac[36], As[36];
  for (int m = 0; m < 36; ++m) { Ac[m] = std::cos(2 * M_PI * P[m]); As[m] = std::sin(2 * M_PI * P[m]); }
  double err = 0;
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
    double re = 0, im = 0;
    for (int g = 0; g < 6; ++g) {
      re += Ac[6 * i + g] * Ac[6 * j + g] + As[6 * i + g] * As[6 * j + g];
      im += As[6 * i + g] * Ac[6 * j + g] - Ac[6 * i + g] * As[6 * j + g];
    }
    if (i == j) re -= 6.0;
    err = std::max(err, std::max(std::fabs(re), std::fabs(im)));
  }
  return err;
}

// Levenberg-Marquardt on the 30 real residuals F(P); returns true on convergence
// (||F|| < tol).  P is modified in place.
inline bool lm_solve(double *P, int maxit = 200, double tol = 1e-14) {
  double mu = 1e-2;
  double F[30];
  nf::residualF(P, F);
  double f2 = 0; for (double v : F) f2 += v * v;
  for (int it = 0; it < maxit; ++it) {
    if (std::sqrt(f2) < tol) return true;
    // Jacobian 30 x 36
    static double J[30][36];
    for (auto &r : J) for (double &v : r) v = 0;
    int n = 0;
    for (int i = 1; i < 6; ++i) for (int j = 0; j < i; ++j) {
      for (int g = 0; g < 6; ++g) {
        double c = std::cos(2 * M_PI * (P[6 * i + g] - P[6 * j + g]));
        double s = std::sin(2 * M_PI * (P[6 * i + g] - P[6 * j + g]));
        J[n][6 * i + g] += c; J[n][6 * j + g] -= c;
        J[n + 1][6 * i + g] += s; J[n + 1][6 * j + g] -= s;
      }
      n += 2;
    }
    double A[36][36], g[36];
    for (int a = 0; a < 36; ++a) {
      g[a] = 0;
      for (int r = 0; r < 30; ++r) g[a] += J[r][a] * F[r];
      for (int b = 0; b < 36; ++b) { double sacc = 0; for (int r = 0; r < 30; ++r) sacc += J[r][a] * J[r][b]; A[a][b] = sacc; }
    }
    bool stepped = false;
    for (int tries = 0; tries < 12 && !stepped; ++tries) {
      double L[36][36];
      bool ok = true;
      for (int i = 0; i < 36 && ok; ++i) for (int j = 0; j <= i; ++j) {
        double sacc = A[i][j] + (i == j ? mu : 0.0);
        for (int m = 0; m < j; ++m) sacc -= L[i][m] * L[j][m];
        if (i == j) { if (!(sacc > 1e-300)) { ok = false; break; } L[i][i] = std::sqrt(sacc); }
        else L[i][j] = sacc / L[j][j];
      }
      if (!ok) { mu *= 10; continue; }
      double y[36], d[36];
      for (int i = 0; i < 36; ++i) { double sacc = -g[i]; for (int m = 0; m < i; ++m) sacc -= L[i][m] * y[m]; y[i] = sacc / L[i][i]; }
      for (int i = 35; i >= 0; --i) { double sacc = y[i]; for (int m = i + 1; m < 36; ++m) sacc -= L[m][i] * d[m]; d[i] = sacc / L[i][i]; }
      double Pn[36], Fn[30];
      for (int m = 0; m < 36; ++m) Pn[m] = P[m] + d[m];
      nf::residualF(Pn, Fn);
      double fn2 = 0; for (double v : Fn) fn2 += v * v;
      if (fn2 < f2) {
        for (int m = 0; m < 36; ++m) P[m] = Pn[m];
        for (int m = 0; m < 30; ++m) F[m] = Fn[m];
        f2 = fn2; mu = std::max(mu / 4, 1e-15); stepped = true;
      } else mu *= 4;
    }
    if (!stepped) break;
  }
  return std::sqrt(f2) < tol;
}

// A random exact Hadamard matrix (phases in cycles, with random row/column
// gauge phases).  Returns false if the LM start did not converge.
inline bool random_hadamard(std::mt19937_64 &rng, double *P) {
  std::uniform_real_distribution<double> U(0.0, 1.0);
  for (int m = 0; m < 36; ++m) P[m] = U(rng);
  // Seed near a dephased structure: first row/col zero is a valid gauge
  for (int j = 0; j < 6; ++j) { P[j] = 0; P[6 * j] = 0; }
  if (!lm_solve(P)) return false;
  // re-gauge: random row and column phases
  double r[6], c[6];
  for (int i = 0; i < 6; ++i) { r[i] = U(rng); c[i] = U(rng); }
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) P[6 * i + j] += r[i] + c[j];
  return true;
}

}  // namespace hd
