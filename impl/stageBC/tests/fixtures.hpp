// Test fixtures: exact complex Hadamard matrices (as phase matrices in cycles), exact unbiased vectors by Newton.
#pragma once
#include "../root.hpp"
#include <random>

namespace fx {
using namespace mb;

// Fourier family F6(a,b), phases in cycles (a, b in radians), flattened 6i+j.
inline void fourierPhases(double a, double b, double* P) {
  const double TWO_PI = 6.283185307179586476925286766559;
  for (int j = 0; j < 6; ++j)
    for (int k = 0; k < 6; ++k) {
      double R = 0;
      if (j % 2 == 1 && (k % 3) == 1) R = a;
      if (j % 2 == 1 && (k % 3) == 2) R = b;
      P[6 * j + k] = (double)((j * k) % 6) / 6.0 + R / TWO_PI;
    }
}

// Bjorck circulant C6 phases: first row (1, i d, -d, -i, -conj d, i conj d), d = (1-sqrt3)/2 + i sqrt(sqrt3/2).
inline void bjorckPhases(double* P) {
  const double s3 = std::sqrt(3.0);
  double dre = (1 - s3) / 2, dim = std::sqrt(s3 / 2);
  double cre[6] = {1, -dim, -dre, 0, dre, -dim};
  double cim[6] = {0, dre, -dim, -1, dim, dre};
  // i d = (-dim, dre); -d = (-dre,-dim); -i = (0,-1); -conj d = (-dre, dim) ; i conj d = (dim, dre)
  cre[4] = -dre; cim[4] = dim;
  cre[5] = dim;  cim[5] = dre;
  const double TWO_PI = 6.283185307179586476925286766559;
  for (int j = 0; j < 6; ++j)
    for (int k = 0; k < 6; ++k) {
      int t = ((k - j) % 6 + 6) % 6;
      P[6 * j + k] = std::atan2(cim[t], cre[t]) / TWO_PI;
    }
}

// Horizontal displacement delta = Pi(D0 H - D0 P - N) (double) and the aligned Hadamard  P + delta.
inline void aligned(const double* P, const double* H, double* delta) {
  double dp[36], dh[36], v[36];
  D0D(P, dp);
  D0D(H, dh);
  for (int m = 0; m < 36; ++m) v[m] = dh[m] - dp[m] - rnd(dh[m] - dp[m]);
  PiD(v, delta);
}
inline double norm36(const double* x) {
  double s = 0;
  for (int m = 0; m < 36; ++m) s += x[m] * x[m];
  return std::sqrt(s);
}

// exact (to ~1e-15) unit vectors e(beta)/sqrt6, beta_0 = 0, unbiased to the standard basis and to the columns of e(Hphase)/sqrt6.
struct UVec { double beta[6]; };

inline bool newtonUnbiased(const double* H, double* beta) {
  const double TWO_PI = 6.283185307179586476925286766559;
  double lambda = 1e-3;
  for (int it = 0; it < 200; ++it) {
    double f[6], Jm[6][5];
    double res = 0;
    for (int i = 0; i < 6; ++i) {
      double sr = 0, si = 0;
      double qr[6], qi[6];
      for (int j = 0; j < 6; ++j) {
        double ang = TWO_PI * (beta[j] - H[6 * j + i]);
        qr[j] = std::cos(ang) / std::sqrt(6.0);
        qi[j] = std::sin(ang) / std::sqrt(6.0);
        sr += qr[j];
        si += qi[j];
      }
      f[i] = sr * sr + si * si - 1;
      res += f[i] * f[i];
      for (int j = 1; j < 6; ++j) {
        // d|S|^2/d beta_j = -4 pi Im(conj(S) q) = -4 pi (sr*qi - si*qr)
        Jm[i][j - 1] = -2 * TWO_PI * (sr * qi[j] - si * qr[j]);
      }
    }
    if (std::sqrt(res) < 1e-15) return true;
    // LM step on 5 unknowns
    for (int tr = 0; tr < 20; ++tr) {
      std::vector<double> A(25, 0.0), b(5, 0.0), x;
      for (int p = 0; p < 5; ++p) {
        for (int q = 0; q < 5; ++q) {
          double s = 0;
          for (int i = 0; i < 6; ++i) s += Jm[i][p] * Jm[i][q];
          A[p * 5 + q] = s + (p == q ? lambda : 0);
        }
        double s = 0;
        for (int i = 0; i < 6; ++i) s += Jm[i][p] * f[i];
        b[p] = -s;
      }
      if (!solveLin(A, b, 5, x)) { lambda *= 10; continue; }
      double nb[6] = {0};
      for (int j = 1; j < 6; ++j) nb[j] = beta[j] + x[j - 1];
      double res2 = 0;
      for (int i = 0; i < 6; ++i) {
        double sr = 0, si = 0;
        for (int j = 0; j < 6; ++j) {
          double ang = TWO_PI * (nb[j] - H[6 * j + i]);
          sr += std::cos(ang) / std::sqrt(6.0);
          si += std::sin(ang) / std::sqrt(6.0);
        }
        double fi = sr * sr + si * si - 1;
        res2 += fi * fi;
      }
      if (res2 < res) { for (int j = 1; j < 6; ++j) beta[j] = nb[j]; lambda = std::max(lambda / 3, 1e-12); break; }
      lambda *= 4;
    }
    if (lambda > 1e12) return false;
  }
  return false;
}

inline double ipAbs2(const double* b1, const double* b2) {
  const double TWO_PI = 6.283185307179586476925286766559;
  double sr = 0, si = 0;
  for (int j = 0; j < 6; ++j) { sr += std::cos(TWO_PI * (b2[j] - b1[j])); si += std::sin(TWO_PI * (b2[j] - b1[j])); }
  return (sr * sr + si * si);  // unscaled |sum|^2 ; orthogonal: 0, unbiased: 6
}

inline std::vector<UVec> findUnbiasedVectors(const double* H, int starts, std::mt19937_64& rng) {
  std::vector<UVec> out;
  std::uniform_real_distribution<double> U(0, 1);
  for (int s = 0; s < starts; ++s) {
    UVec v;
    v.beta[0] = 0;
    for (int j = 1; j < 6; ++j) v.beta[j] = U(rng);
    if (!newtonUnbiased(H, v.beta)) continue;
    for (int j = 1; j < 6; ++j) v.beta[j] -= std::floor(v.beta[j]);  // lift into [0,1)
    bool dup = false;
    for (auto& w : out)
      if (ipAbs2(w.beta, v.beta) > 36 - 1e-8) { dup = true; break; }
    if (!dup) out.push_back(v);
  }
  return out;
}

}  // namespace fx
