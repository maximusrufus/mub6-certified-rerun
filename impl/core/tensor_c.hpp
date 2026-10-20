// tensor_c.hpp -- the optional "tensor improvement" of the Root remainder constant C
// (phase-cover.tex phase:remainder-constants second choice, phase:tensor-gram; rounding.tex round:tensor-c).
//
//   C >= pi * b,   b = sqrt(6) * ||S||,   S (30 x 216) has rows  M_{:,a}^T (x) (e_i - e_j)^T
//   (S S^T)_{ab} = ((e_i - e_j) . (e_i' - e_j')) * (M_{:,a} . M_{:,b}),      M = H_V L_V^T (exact stored doubles)
//
// Certificate (paper, Section phase:stored-root-certificate / round:tensor-c, margins replaced by intervals):
//   a real change of basis U (proposal: Jacobi eigenvectors of mid(S S^T)) with  beta >= ||U^T U - I||_F < .1
//   gives ||U^{-1}||^2 <= 1/(1-beta)  (U^T U >= (1-beta) I), and for the exact symmetric A = S S^T,
//   ||A|| = ||U^{-T} (U^T A U) U^{-1}|| <= row / (1-beta)   with  row >= max_a sum_b |(U^T A U)_{ab}|
//   (absolute row sum bounds the spectral norm of a symmetric matrix; the true U^T A U lies in the enclosure).
// Everything after the proposal is interval arithmetic on the exact stored doubles (H_V, L_V, U); the result is
// an upper bound, so min(pi sqrt8 mnorm, tensor) is a valid C (both are upper bounds for a remainder constant).
// Must run inside an iv::RoundUp guard.  Returns +inf when the proposal is unusable (never a smaller C).
#pragma once
#include "interval.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace tc {
using iv::Interval;

// row index a (0..29) -> edge (i, j), i > j, ordered by increasing i then j; a = 2m (Re) or 2m+1 (Im)
inline void edgeOf(int a, int& i, int& j) {
  int m = a / 2, e = 0;
  for (i = 1; i < 6; ++i)
    for (j = 0; j < i; ++j, ++e)
      if (e == m) return;
}
inline int edgeSign(int a, int b) {  // (e_i - e_j) . (e_i' - e_j')
  int i, j, ii, jj;
  edgeOf(a, i, j);
  edgeOf(b, ii, jj);
  return (i == ii) - (i == jj) - (j == ii) + (j == jj);
}

// Cyclic Jacobi (proposal only).  A n x n row-major (destroyed); V columns = eigenvectors.
inline void jacobiProp(std::vector<double>& A, int n, std::vector<double>& V) {
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
        double th = (A[(size_t)q * n + q] - A[(size_t)p * n + p]) / (2 * apq);
        double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
        double c = 1 / std::sqrt(t * t + 1), s = t * c;
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
}

struct TensorInfo {
  double beta = 0, row = 0;
  double plainRowBound = 0;  // abs-row-sum bound of S S^T itself (no change of basis), for diagnostics
};

// HV: 36 x g row-major, LV: 30 x g row-major (exact stored doubles).  pi: enclosure of pi.
// Returns hi(pi * sqrt(6 * row / (1 - beta))) or +inf.
inline double tensorC(const std::vector<double>& HV, const std::vector<double>& LV, int g, const Interval& pi, TensorInfo* info = nullptr) {
  iv::require_upward();
  const double INF = std::numeric_limits<double>::infinity();
  if (g <= 0) return INF;
  // Gh = H_V^T H_V (g x g)
  std::vector<Interval> Gh((size_t)g * g);
  for (int a = 0; a < g; ++a)
    for (int b = 0; b < g; ++b) {
      Interval s(0.0);
      for (int m = 0; m < 36; ++m) s = s + Interval(HV[(size_t)m * g + a]) * Interval(HV[(size_t)m * g + b]);
      Gh[(size_t)a * g + b] = s;
    }
  // Y = L_V Gh (30 x g), GM = Y L_V^T (30 x 30) = M^T M
  std::vector<Interval> Y((size_t)30 * g), A((size_t)30 * 30);
  for (int r = 0; r < 30; ++r)
    for (int b = 0; b < g; ++b) {
      Interval s(0.0);
      for (int a = 0; a < g; ++a) s = s + Interval(LV[(size_t)r * g + a]) * Gh[(size_t)a * g + b];
      Y[(size_t)r * g + b] = s;
    }
  for (int r = 0; r < 30; ++r)
    for (int c = 0; c < 30; ++c) {
      int sg = edgeSign(r, c);
      if (sg == 0) { A[(size_t)r * 30 + c] = Interval(0.0); continue; }
      Interval s(0.0);
      for (int b = 0; b < g; ++b) s = s + Y[(size_t)r * g + b] * Interval(LV[(size_t)c * g + b]);
      A[(size_t)r * 30 + c] = s * Interval((double)sg);
    }
  if (info) {
    double best = 0;
    for (int r = 0; r < 30; ++r) {
      Interval rs(0.0);
      for (int c = 0; c < 30; ++c) rs = rs + iv::abs(A[(size_t)r * 30 + c]);
      best = std::max(best, rs.hi);
    }
    info->plainRowBound = best;
  }
  // proposal U
  std::vector<double> Am((size_t)30 * 30), U;
  for (size_t t = 0; t < Am.size(); ++t) Am[t] = A[t].mid();
  for (int r = 0; r < 30; ++r)
    for (int c = r + 1; c < 30; ++c) { double v = 0.5 * (Am[(size_t)r * 30 + c] + Am[(size_t)c * 30 + r]); Am[(size_t)r * 30 + c] = Am[(size_t)c * 30 + r] = v; }
  jacobiProp(Am, 30, U);
  for (double x : U) if (!(std::fabs(x) <= 2)) return INF;
  // beta = ||U^T U - I||_F (interval)
  Interval bsq(0.0);
  for (int a = 0; a < 30; ++a)
    for (int b = 0; b < 30; ++b) {
      Interval s(0.0);
      for (int m = 0; m < 30; ++m) s = s + Interval(U[(size_t)m * 30 + a]) * Interval(U[(size_t)m * 30 + b]);
      if (a == b) s = s - Interval(1.0);
      bsq = bsq + iv::sqr(s);
    }
  double beta = iv::sqrt(bsq).hi;
  if (info) info->beta = beta;
  if (!(beta < 0.1)) return INF;
  // X = A U (30 x 30), B = U^T X; row = max abs row sum of B
  std::vector<Interval> X((size_t)30 * 30);
  for (int r = 0; r < 30; ++r)
    for (int c = 0; c < 30; ++c) {
      Interval s(0.0);
      for (int m = 0; m < 30; ++m) s = s + A[(size_t)r * 30 + m] * Interval(U[(size_t)m * 30 + c]);
      X[(size_t)r * 30 + c] = s;
    }
  double row = 0;
  for (int r = 0; r < 30; ++r) {
    Interval rs(0.0);
    for (int c = 0; c < 30; ++c) {
      Interval s(0.0);
      for (int m = 0; m < 30; ++m) s = s + Interval(U[(size_t)m * 30 + r]) * X[(size_t)m * 30 + c];
      rs = rs + iv::abs(s);
    }
    row = std::max(row, rs.hi);
  }
  if (info) info->row = row;
  if (!(row >= 0) || !std::isfinite(row)) return INF;
  Interval val = Interval(6.0) * Interval(row) / (Interval(1.0) - Interval(beta));
  if (!(val.lo >= 0)) val = Interval(0.0, val.hi);
  double C = (pi * iv::sqrt(val)).hi;
  if (!std::isfinite(C)) return INF;
  return C;
}

}  // namespace tc
