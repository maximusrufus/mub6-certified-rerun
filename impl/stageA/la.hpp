// la.hpp -- plain binary64 dense helpers used ONLY to PROPOSE things
// (eigenvectors, solves, directions).  Nothing here is trusted: every proposal
// is re-checked in interval arithmetic by the caller.
#pragma once
#include <cmath>
#include <vector>

namespace la {

// Cyclic Jacobi for a symmetric n x n matrix A (row-major, destroyed).
// On return ev[i] are eigenvalues, V (row-major n x n) has eigenvectors in columns.
inline void jacobi(std::vector<double>& A, int n, std::vector<double>& ev, std::vector<double>& V) {
  V.assign((size_t)n * n, 0.0);
  for (int i = 0; i < n; ++i) V[(size_t)i * n + i] = 1.0;
  for (int sweep = 0; sweep < 60; ++sweep) {
    double off = 0;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) off += A[(size_t)i * n + j] * A[(size_t)i * n + j];
    if (off < 1e-30) break;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) {
        double apq = A[(size_t)p * n + q];
        if (std::fabs(apq) < 1e-300) continue;
        double app = A[(size_t)p * n + p], aqq = A[(size_t)q * n + q];
        double theta = (aqq - app) / (2.0 * apq);
        double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (int k = 0; k < n; ++k) {
          double akp = A[(size_t)k * n + p], akq = A[(size_t)k * n + q];
          A[(size_t)k * n + p] = c * akp - s * akq;
          A[(size_t)k * n + q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          double apk = A[(size_t)p * n + k], aqk = A[(size_t)q * n + k];
          A[(size_t)p * n + k] = c * apk - s * aqk;
          A[(size_t)q * n + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          double vkp = V[(size_t)k * n + p], vkq = V[(size_t)k * n + q];
          V[(size_t)k * n + p] = c * vkp - s * vkq;
          V[(size_t)k * n + q] = s * vkp + c * vkq;
        }
      }
  }
  ev.resize(n);
  for (int i = 0; i < n; ++i) ev[i] = A[(size_t)i * n + i];
}

// Solve SPD system (A + lam I) x = b by Cholesky; returns false if not PD.
inline bool chol_solve(const double* A, int n, double lam, const double* b, double* x) {
  std::vector<double> L((size_t)n * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = A[(size_t)i * n + j] + (i == j ? lam : 0.0);
      for (int k = 0; k < j; ++k) s -= L[(size_t)i * n + k] * L[(size_t)j * n + k];
      if (i == j) {
        if (!(s > 1e-300)) return false;
        L[(size_t)i * n + i] = std::sqrt(s);
      } else
        L[(size_t)i * n + j] = s / L[(size_t)j * n + j];
    }
  std::vector<double> y(n);
  for (int i = 0; i < n; ++i) {
    double s = b[i];
    for (int k = 0; k < i; ++k) s -= L[(size_t)i * n + k] * y[k];
    y[i] = s / L[(size_t)i * n + i];
  }
  for (int i = n - 1; i >= 0; --i) {
    double s = y[i];
    for (int k = i + 1; k < n; ++k) s -= L[(size_t)k * n + i] * x[k];
    x[i] = s / L[(size_t)i * n + i];
  }
  return true;
}

}  // namespace la
