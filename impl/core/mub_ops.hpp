// 6-dimensional complex interval vectors, inner products, unbiasedness / orthonormality checks.
// Every function here must be called inside an iv::RoundUp guard (per thread).
#pragma once
#include <array>
#include "interval.hpp"

namespace mub {
using iv::CInterval;
using iv::Interval;

constexpr int D = 6;
using CVec6 = std::array<CInterval, D>;
using Basis = std::array<CVec6, D>;  // Basis[k] = k-th basis vector

// <u,v> = sum_j conj(u_j) v_j
inline CInterval inner(const CVec6& u, const CVec6& v) {
  CInterval s = iv::conj(u[0]) * v[0];
  for (int j = 1; j < D; ++j) s = s + iv::conj(u[j]) * v[j];
  return s;
}
inline Interval norm2(const CVec6& u) {
  Interval s = iv::abs2(u[0]);
  for (int j = 1; j < D; ++j) s = s + iv::abs2(u[j]);
  return s;
}
// |<u,v>|^2
inline Interval ip_abs2(const CVec6& u, const CVec6& v) { return iv::abs2(inner(u, v)); }
// |<u,v>|^2 - 1/6  (the unbiasedness residual; an exact MUB pair gives an interval containing 0)
inline Interval unbiased_residual(const CVec6& u, const CVec6& v) { return ip_abs2(u, v) - iv::one_sixth(); }

// Consistent with unbiasedness: the interval |<u,v>|^2 overlaps 1/6.
inline bool unbiased_consistent(const CVec6& u, const CVec6& v) {
  return unbiased_residual(u, v).contains_zero();
}
// Certified biased: |<u,v>|^2 is disjoint from 1/6, so no exact vectors inside can be unbiased.
inline bool certainly_biased(const CVec6& u, const CVec6& v) {
  return iv::disjoint(ip_abs2(u, v), iv::one_sixth());
}
// Consistent with unit norm.
inline bool unit_consistent(const CVec6& u) { return iv::overlaps(norm2(u), Interval(1.0)); }
// Consistent with orthogonality: both Re and Im of <u,v> contain 0.
inline bool orthogonal_consistent(const CVec6& u, const CVec6& v) {
  CInterval p = inner(u, v);
  return p.re.contains_zero() && p.im.contains_zero();
}
// Consistent with orthonormal basis.
inline bool orthonormal_consistent(const Basis& B) {
  for (int a = 0; a < D; ++a) {
    if (!unit_consistent(B[a])) return false;
    for (int b = a + 1; b < D; ++b)
      if (!orthogonal_consistent(B[a], B[b])) return false;
  }
  return true;
}
// Consistent with the two bases being unbiased: all 36 pairs.
inline bool bases_unbiased_consistent(const Basis& A, const Basis& B) {
  for (int a = 0; a < D; ++a)
    for (int b = 0; b < D; ++b)
      if (!unbiased_consistent(A[a], B[b])) return false;
  return true;
}
// Certified: some pair is certainly biased.
inline bool bases_certainly_biased(const Basis& A, const Basis& B) {
  for (int a = 0; a < D; ++a)
    for (int b = 0; b < D; ++b)
      if (certainly_biased(A[a], B[b])) return true;
  return false;
}

// ---- known bases (entries are rigorous intervals) ----
inline Basis standard_basis() {
  Basis B;
  for (int k = 0; k < D; ++k)
    for (int j = 0; j < D; ++j) B[k][j] = CInterval(j == k ? 1.0 : 0.0, 0.0);
  return B;
}
// Fourier basis: component j of vector k is exp(2 pi i jk/6)/sqrt(6)
inline Basis fourier_basis() {
  iv::RoundUp g;
  Interval s = iv::inv_sqrt6();
  Basis B;
  for (int k = 0; k < D; ++k)
    for (int j = 0; j < D; ++j) B[k][j] = iv::unit_root((j * k) % 6, 6) * s;
  return B;
}
// Bjorck circulant C6 (Tadej-Zyczkowski catalogue): first row (1, i d, -d, -i, -conj d, i conj d),
// d = (1 - sqrt3)/2 + i sqrt(sqrt3/2). Circulant H_{jk} = c[(k-j) mod 6]. A circulant complex Hadamard
// is diagonalised by the Fourier matrix, hence unbiased to both the standard and the Fourier basis.
inline Basis bjorck_c6_basis() {
  iv::RoundUp g;
  Interval s3 = iv::sqrt_iv(3);
  Interval dre = (Interval(1.0) - s3) / Interval(2.0);
  Interval dim = iv::sqrt(s3 / Interval(2.0));
  CInterval d(dre, dim), dc = iv::conj(d);
  CInterval I(0.0, 1.0);
  CInterval c[6] = {CInterval(1.0, 0.0), I * d, -d, CInterval(0.0, -1.0), -dc, I * dc};
  Interval s = iv::inv_sqrt6();
  Basis B;
  for (int k = 0; k < D; ++k)
    for (int j = 0; j < D; ++j) B[k][j] = c[((k - j) % 6 + 6) % 6] * s;
  return B;
}

}  // namespace mub
