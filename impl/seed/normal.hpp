// normal.hpp -- normal forms and symmetry reductions (SPEC sec 1; seed.tex
// seed:equivalences, seed:arrangement; phase-cover.tex phase:alignment,
// phase:transport).  Phase matrices are 6x6 row-major, in CYCLES, entry (i,j)
// at position 6i+j.  Everything here is plain binary64: these routines
// PROPOSE normal forms / integer lifts / symmetries; no certificate depends on
// their rounding (every identity used downstream is exact for any integers N,
// see phase:alignment), and the tests check each identity numerically.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include "arith.hpp"

namespace nf {

typedef double Mat[36];

inline double frac01(double x) { double f = x - std::floor(x); return f >= 1.0 ? 0.0 : f; }

// ---- seed:projection-formula : double centering ----------------------------
inline void Pi(const double *z, double *out) {
  double rs[6] = {0}, cs[6] = {0}, tot = 0;
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) { rs[i] += z[6 * i + j]; cs[j] += z[6 * i + j]; tot += z[6 * i + j]; }
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j)
    out[6 * i + j] = z[6 * i + j] - rs[i] / 6.0 - cs[j] / 6.0 + tot / 36.0;
}

// ---- seed:dephasing-formula : unwrapped dephasing --------------------------
inline void D0(const double *P, double *out) {
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j)
    out[6 * i + j] = P[6 * i + j] - P[6 * i] - P[j] + P[0];
}

// ---- phase:alignment --------------------------------------------------------
// Integer matrix N chosen by rounding (any choice is valid): N = nint(D0(P)-D0(Q)).
inline void alignment_integers(const double *P, const double *Q, double *Nint) {
  double a[36], b[36];
  D0(P, a); D0(Q, b);
  for (int m = 0; m < 36; ++m) Nint[m] = std::nearbyint(a[m] - b[m]);
}
// d = Pi(D0(P) - D0(Q) - N).  Then e(Q+d) is rephasing-equivalent to e(P), and
// d works for e(Q+d+delta) vs e(P+delta) for every horizontal delta.
inline void alignment_delta(const double *P, const double *Q, const double *Nint, double *d) {
  double a[36], b[36], w[36];
  D0(P, a); D0(Q, b);
  for (int m = 0; m < 36; ++m) w[m] = a[m] - b[m] - Nint[m];
  Pi(w, d);
}
// Convenience: d with the rounded N; returns ||d||_F.
inline double alignment_d(const double *P, const double *Q, double *d) {
  double Nint[36]; alignment_integers(P, Q, Nint); alignment_delta(P, Q, Nint, d);
  double s = 0; for (int m = 0; m < 36; ++m) s += d[m] * d[m];
  return std::sqrt(s);
}
// Rephasing equivalence of e(A) and e(B): D0(A)-D0(B) is an integer matrix.
inline bool rephasing_equivalent(const double *A, const double *B, double tol) {
  double a[36], b[36]; D0(A, a); D0(B, b);
  for (int m = 0; m < 36; ++m) { double d = a[m] - b[m]; if (std::fabs(d - std::nearbyint(d)) > tol) return false; }
  return true;
}

// ---- phase:residual / phase:residual-invariance ----------------------------
// F_ij(P) = (1/2 pi i) sum_g e(P_ig - P_jg), 0 <= j < i <= 5; (Re, Im) order,
// i increasing then j: 30 reals.
inline void residualF(const double *P, double *F) {
  int n = 0;
  for (int i = 1; i < 6; ++i) for (int j = 0; j < i; ++j) {
    double re = 0, im = 0;
    for (int g = 0; g < 6; ++g) { double c, s; arith::expair(P[6 * i + g] - P[6 * j + g], c, s); re += c; im += s; }
    // (re + i im)/(2 pi i) = (im - i re)/(2 pi)
    F[n++] = im / (2.0 * arith::PI);
    F[n++] = -re / (2.0 * arith::PI);
  }
}
inline double normF(const double *P) {
  double F[30], s = 0; residualF(P, F); for (double v : F) s += v * v; return std::sqrt(s);
}
// ||F(P)||^2 = (||A A*||_F^2 - 216)/(8 pi^2), A = e(P)
inline double normF2_via_gram(const double *P) {
  double Ac[36], As[36];
  for (int m = 0; m < 36; ++m) arith::expair(P[m], Ac[m], As[m]);
  double fro2 = 0;
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
    double re = 0, im = 0;
    for (int g = 0; g < 6; ++g) {   // (A A*)_ij = sum_g A_ig conj(A_jg)
      re += Ac[6 * i + g] * Ac[6 * j + g] + As[6 * i + g] * As[6 * j + g];
      im += As[6 * i + g] * Ac[6 * j + g] - Ac[6 * i + g] * As[6 * j + g];
    }
    fro2 += re * re + im * im;
  }
  return (fro2 - 216.0) / (8.0 * arith::PI * arith::PI);
}

// ---- phase:transport --------------------------------------------------------
// Phi: optional transpose, row perm, col perm, optional global sign flip.
//   (Phi P)_{ij} = sgn * Pt_{rp[i], cp[j]},  Pt = transpose? P^T : P.
struct Sym {
  bool transpose = false;
  bool flip = false;
  int rp[6] = {0, 1, 2, 3, 4, 5};
  int cp[6] = {0, 1, 2, 3, 4, 5};
};
inline void apply_sym(const Sym &S, const double *P, double *out) {
  double Pt[36];
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) Pt[6 * i + j] = S.transpose ? P[6 * j + i] : P[6 * i + j];
  double sg = S.flip ? -1.0 : 1.0;
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) out[6 * i + j] = sg * Pt[6 * S.rp[i] + S.cp[j]];
}
// Transport of a Root's tangent frame T (k columns of R^36, stored as k rows of 36):
// the same linear isometry acts on every column.
inline void transport_T(const Sym &S, const double *T, int k, double *out) {
  for (int a = 0; a < k; ++a) apply_sym(S, T + 36 * a, out + 36 * a);
}

// ---- seed:arrangement ------------------------------------------------------
// circular gap of six values mod 1 (largest open empty arc); also returns the
// sorted order of the six points and the index (in sorted order) of the point
// that starts right after a largest gap (first largest gap, ties by index).
struct GapInfo { double gap; int startIdx; };
inline GapInfo circ_gap(const double v[6], const int *order /*sorted indices*/) {
  GapInfo g{-1.0, 0};
  for (int a = 0; a < 6; ++a) {
    int i0 = order[a], i1 = order[(a + 1) % 6];
    double d = v[i1] - v[i0];
    if (a == 5) d += 1.0;
    if (d > g.gap) { g.gap = d; g.startIdx = (a + 1) % 6; }
  }
  return g;
}

struct NormalForm {
  bool transposed = false;
  int rowperm[6], colperm[6];  // E_{ij} = D0(P')_{rowperm[i],colperm[j]} mod 1, P' = transposed? P^T : P
  double E[36];                // arranged, dephased, in [0,1)
  double x[6];                 // (0,x1..x5) sorted, BEFORE moving column s+1 (E row 1 is permuted)
  double y[6];                 // E column 1: (0,xi,y2..y5)
  int s;                       // selected position 0..4
  double gx;                   // 1 - x5
  double q9[9];                // nine seed coordinates in cycles, nine-order
  int bx[5], by[6];            // bins of x1..x5 (original order) / y list (by[1] = shared)
};

// Minimum gap over all 15 row pairs and 15 column pairs of P (phases in cycles).
// pairs attaining the minimum are returned (type 0 = rows, 1 = cols).
struct PairRef { int type, a, b; };
inline double min_gap_pairs(const double *P, std::vector<PairRef> *argmins, double tol = 1e-12) {
  struct Cand { double g; PairRef r; };
  std::vector<Cand> all;
  for (int type = 0; type < 2; ++type)
    for (int a = 0; a < 6; ++a) for (int b = a + 1; b < 6; ++b) {
      double d[6];
      for (int j = 0; j < 6; ++j) d[j] = frac01(type == 0 ? P[6 * b + j] - P[6 * a + j] : P[6 * j + b] - P[6 * j + a]);
      int ord[6] = {0, 1, 2, 3, 4, 5};
      std::sort(ord, ord + 6, [&](int u, int v) { return d[u] < d[v]; });
      GapInfo g = circ_gap(d, ord);
      all.push_back({g.gap, {type, a, b}});
    }
  double m = 1e9;
  for (auto &c : all) m = std::min(m, c.g);
  if (argmins) { argmins->clear(); for (auto &c : all) if (c.g <= m + tol) argmins->push_back(c.r); }
  return m;
}

// value mod 1 in [0,1), with values within 1e-9 below 1 snapped to 0: for a
// genuine Hadamard matrix every circular gap is >= 1/6 so such a value is a
// rounding artefact of a repeated point (used only to PROPOSE the normal form).
inline double fsnap(double x) { double f = frac01(x); return f > 1.0 - 1e-9 ? 0.0 : f; }

// bin of t in [0,1] for N bins with the closed-bin boundary convention:
// within 1e-9 of a boundary both neighbours are valid; binCeil picks the lower one.
inline int bin_of(double t, int N, bool binCeil) {
  double u = t * (double)N, r = std::nearbyint(u);
  int b;
  if (std::fabs(u - r) < 1e-9) b = binCeil ? (int)r - 1 : (int)r;
  else b = (int)std::floor(u);
  if (b < 0) b = 0;
  if (b > N - 1) b = N - 1;
  return b;
}

// Construct the seed arrangement of an (approximately exact) Hadamard phase
// matrix P (cycles).  Choices: pairIdx selects among minimum-gap pairs,
// gapTie among the (near-)largest gaps of that pair, binCeil the boundary
// convention of the closed bins.  Returns false if some internal consistency
// check fails (bins not monotone).
inline bool arrange(const double *P0, int N, int pairIdx, int gapTie, bool binCeil, NormalForm &nfm) {
  std::vector<PairRef> am;
  min_gap_pairs(P0, &am);
  PairRef pr = am[pairIdx % am.size()];
  double P[36];
  nfm.transposed = (pr.type == 1);
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) P[6 * i + j] = nfm.transposed ? P0[6 * j + i] : P0[6 * i + j];
  int r0 = pr.a, r1 = pr.b;
  double d[6];
  for (int j = 0; j < 6; ++j) d[j] = fsnap(P[6 * r1 + j] - P[6 * r0 + j]);
  int ord[6] = {0, 1, 2, 3, 4, 5};
  std::sort(ord, ord + 6, [&](int u, int v) { return d[u] < d[v]; });
  std::vector<int> starts;
  double best = -1;
  for (int a = 0; a < 6; ++a) {
    double dd = d[ord[(a + 1) % 6]] - d[ord[a]]; if (a == 5) dd += 1.0;
    if (dd > best + 1e-12) { best = dd; starts.assign(1, (a + 1) % 6); }
    else if (std::fabs(dd - best) <= 1e-12) starts.push_back((a + 1) % 6);
  }
  int st = starts[gapTie % starts.size()];
  int start = ord[st];
  // order the columns along the arc beginning at `start` (start first)
  auto arcorder = [&](int rA, int rB, int c0, int *cols, double *xv) {
    double xx[6];
    for (int j = 0; j < 6; ++j) xx[j] = fsnap((P[6 * rB + j] - P[6 * rA + j]) - (P[6 * rB + c0] - P[6 * rA + c0]));
    int o[6] = {0, 1, 2, 3, 4, 5};
    std::sort(o, o + 6, [&](int u, int v) {
      if (u == c0 || v == c0) return u == c0 && v != c0;   // start column first
      return xx[u] < xx[v];
    });
    for (int m = 0; m < 6; ++m) { cols[m] = o[m]; xv[m] = xx[o[m]]; }
  };
  int cols[6]; double xv[6];
  arcorder(r0, r1, start, cols, xv);
  if (!(xv[1] <= xv[5] - xv[4] + 1e-12)) {
    // orientation fails: swap rows, start from the old x5 column, reverse
    int newstart = cols[5];
    std::swap(r0, r1);
    arcorder(r0, r1, newstart, cols, xv);
  }
  int bxs[5];
  for (int m = 0; m < 5; ++m) bxs[m] = bin_of(xv[m + 1], N, binCeil);
  for (int m = 1; m < 5; ++m) if (bxs[m] < bxs[m - 1]) return false;
  int s = 0, bv = std::abs(2 * bxs[0] + 1 - N);
  for (int m = 1; m < 5; ++m) { int v = std::abs(2 * bxs[m] + 1 - N); if (v < bv) { bv = v; s = m; } }
  // new column order: col0, col s+1, the other four in former order
  int cn[6]; cn[0] = cols[0]; cn[1] = cols[s + 1];
  { int p = 2; for (int m = 1; m < 6; ++m) if (m != s + 1) cn[p++] = cols[m]; }
  int others[4]; { int p = 0; for (int i = 0; i < 6; ++i) if (i != r0 && i != r1) others[p++] = i; }
  auto ent = [&](int i, int j) {   // dephased phase w.r.t. row r0 and column cn[0]
    return fsnap(P[6 * i + j] - P[6 * i + cn[0]] - P[6 * r0 + j] + P[6 * r0 + cn[0]]);
  };
  std::sort(others, others + 4, [&](int a, int b) { return ent(a, cn[1]) < ent(b, cn[1]); });
  int rn[6] = {r0, r1, others[0], others[1], others[2], others[3]};
  for (int i = 0; i < 6; ++i) { nfm.rowperm[i] = rn[i]; nfm.colperm[i] = cn[i]; }
  for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) nfm.E[6 * i + j] = ent(rn[i], cn[j]);
  for (int j = 0; j < 6; ++j) { nfm.E[j] = 0.0; nfm.E[6 * j] = 0.0; }   // exact dephasing
  for (int j = 0; j < 6; ++j) nfm.x[j] = xv[j];   // original (sorted) order, before the column move
  for (int i = 0; i < 6; ++i) nfm.y[i] = nfm.E[6 * i + 1];
  nfm.s = s;
  nfm.gx = 1.0 - nfm.x[5];
  static const int pos9[9] = {7, 8, 9, 10, 11, 13, 19, 25, 31};
  for (int m = 0; m < 9; ++m) nfm.q9[m] = nfm.E[pos9[m]];
  for (int m = 0; m < 5; ++m) nfm.bx[m] = bxs[m];
  nfm.by[0] = 0; nfm.by[1] = bxs[s];   // shared entry: exactly the x bin
  for (int m = 2; m < 6; ++m) nfm.by[m] = bin_of(nfm.y[m], N, binCeil);
  return true;
}

}  // namespace nf
