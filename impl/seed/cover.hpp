// cover.hpp -- seed cover primitives (seed.tex sec 2.2-2.4) with interval
// certificates: every rejection below is "the inequality is certainly violated
// for every point of the enclosures".  No additive paper margins.
// All functions must be called by a thread holding an iv::RoundUp guard.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>
#include "arith.hpp"

namespace seedcov {

using arith::abs_hi;
using arith::abs_lo;
using arith::mulp;
using iv::Interval;

static const double REUSE_NUM = 32.0, REUSE_DEN = 1000.0;   // R = .032 exactly
static const double EXT_MIN_DIV = 1e-6;                      // extend: skip |a| below this (drops a constraint: sound)

struct Tables {
  int N;
  std::vector<Interval> ec, es;           // e(u_b) = (cos, sin), enclosures of the exact rational phases
  std::vector<double> pv0, pv1;           // exact doubles: direction perpendicular to mid(J_b)
  std::vector<double> sv0, sv1;           // exact doubles: direction mid(J_b)
  std::vector<Interval> PP, SP;           // PP[b*N+b'] = pv_b . J_b' ; SP = sv_b . J_b'  (signed enclosures)
  std::vector<double> PPa, SPa;           // upper bounds of |.|
  std::vector<double> pnv, snv;           // upper bounds of the norms of pv_b / sv_b
  Interval tp, inv2pi, Eg, whalf, coefT;  // 2pi, 1/(2pi), 5pi/(4N), 1/(2N), 5 pi/(4N^2)
  double tp_hi, Eg_hi, whalf_hi, coef_hi; // upper bounds
  double RN2lo;                           // lower bound of (R N)^2, R = .032
  double mid(int b) const { return ((double)b + 0.5) / (double)N; }   // heuristic/test only

  explicit Tables(int N_) : N(N_) {
    iv::RoundUp g;
    ec.resize(N); es.resize(N); pv0.resize(N); pv1.resize(N); sv0.resize(N); sv1.resize(N);
    pnv.resize(N); snv.resize(N);
    for (int b = 0; b < N; ++b) {
      iv::CInterval z = arith::cis_mid(b, N);
      ec[b] = z.re; es[b] = z.im;
      sv0[b] = ec[b].mid(); sv1[b] = es[b].mid();
      pv0[b] = -sv1[b]; pv1[b] = sv0[b];
      pnv[b] = iv::sqrt(Interval(pv0[b] * pv0[b]) + Interval(pv1[b] * pv1[b])).hi;
      snv[b] = pnv[b];
    }
    PP.resize((size_t)N * N); SP.resize((size_t)N * N); PPa.resize((size_t)N * N); SPa.resize((size_t)N * N);
    for (int b = 0; b < N; ++b) for (int c = 0; c < N; ++c) {
      Interval p = mulp(ec[c], pv0[b]) + mulp(es[c], pv1[b]);
      Interval s = mulp(ec[c], sv0[b]) + mulp(es[c], sv1[b]);
      PP[(size_t)b * N + c] = p; SP[(size_t)b * N + c] = s;
      PPa[(size_t)b * N + c] = abs_hi(p); SPa[(size_t)b * N + c] = abs_hi(s);
    }
    Interval pi = iv::pi_iv(), Nn((double)N);
    tp = pi * Interval(2.0);
    inv2pi = Interval(1.0) / tp;
    Eg = Interval(5.0) * pi / (Interval(4.0) * Nn);
    whalf = Interval(1.0) / (Interval(2.0) * Nn);
    coefT = Interval(5.0) * pi / (Interval(4.0) * Nn * Nn);
    tp_hi = tp.hi; Eg_hi = Eg.hi; whalf_hi = whalf.hi; coef_hi = coefT.hi;
    Interval R = Interval(REUSE_NUM) / Interval(REUSE_DEN);
    Interval RN = R * Nn;
    RN2lo = iv::sqr(RN).lo;
  }
};

// Data of a five-entry list (0,u_0..u_4): bins and the enclosure of
// c = (1+sum e(u_j))/(2 pi i)  = (sum sin, -(1+sum cos))/(2 pi).
struct ListData {
  int b[5];
  Interval c0, c1;
};

inline void make_ld(const Tables &T, const int *b, ListData &d) {
  Interval re(1.0), im(0.0);
  for (int j = 0; j < 5; ++j) { d.b[j] = b[j]; re = re + T.ec[b[j]]; im = im + T.es[b[j]]; }
  d.c0 = im * T.inv2pi;
  d.c1 = -(re * T.inv2pi);
}

// sumok (seed:sum-support): in direction v (exact doubles), certainly violated iff
//   lb|v.c'| > (2pi)^up * ( coef^up * ||v||^up + (w/2)^up * S^up ),  c' = 2 pi c.
// S^up = upper bound of sum_k |v.J_k|.  Any exact direction is valid.
inline bool sum_certainly_violated(const Tables &T, const Interval &cx, const Interval &cy, double v0, double v1,
                                   double nv_hi, double S_hi) {
  Interval X = mulp(cx, v0) + mulp(cy, v1);
  double lhs = abs_lo(X);
  double rhs = T.tp_hi * (T.coef_hi * nv_hi + T.whalf_hi * S_hi);   // up: nonnegative operands
  return lhs > rhs;
}

// returns true iff the list is NOT certainly infeasible (kept).
inline bool sumok_bins(const Tables &T, const int *b) {
  Interval sx(0.0), sy(1.0);
  for (int j = 0; j < 5; ++j) { sx = sx + T.es[b[j]]; sy = sy + T.ec[b[j]]; }
  Interval cx = sx, cy = -sy;    // c' = (sum sin, -(1+sum cos))
  for (int j = 0; j < 5; ++j) {
    const double *row = &T.PPa[(size_t)b[j] * T.N];
    double S = 0.0;
    for (int k = 0; k < 5; ++k) S += row[b[k]];     // up
    if (sum_certainly_violated(T, cx, cy, T.pv0[b[j]], T.pv1[b[j]], T.pnv[b[j]], S)) return false;
  }
  // residual direction: midpoint of c', normalised (exact doubles); skipped if zero
  double m0 = cx.mid(), m1 = cy.mid();
  double nm = std::sqrt(m0 * m0 + m1 * m1);
  if (nm > 0.0 && std::isfinite(nm)) {
    double v0 = m0 / nm, v1 = m1 / nm;
    double nv = std::sqrt(v0 * v0 + v1 * v1);      // up (nonnegative operands)
    double S = 0.0;
    for (int k = 0; k < 5; ++k) {
      Interval t = mulp(T.ec[b[k]], v0) + mulp(T.es[b[k]], v1);
      S += abs_hi(t);
    }
    if (sum_certainly_violated(T, cx, cy, v0, v1, nv, S)) return false;
  }
  return true;
}

// ---- gap inequalities: exact integer arithmetic ------------------------------
// With G = 2N * gap the midpoint gap of sorted bins q0..q4 (plus fixed zero) is the
// integer max{2q0+1, 2(q_{i+1}-q_i), 2N-(2q4+1)}; all comparisons are exact.
inline int mid_gap2N(int N, const int q[5]) {
  int g = 2 * q[0] + 1;
  for (int i = 1; i < 5; ++i) g = std::max(g, 2 * (q[i] - q[i - 1]));
  g = std::max(g, 2 * N - (2 * q[4] + 1));
  return g;
}

inline int select_s(int N, const int bx[5]) {
  int best = 0, bv = std::abs(2 * bx[0] + 1 - N);
  for (int s = 1; s < 5; ++s) {
    int v = std::abs(2 * bx[s] + 1 - N);
    if (v < bv) { bv = v; best = s; }
  }
  return best;
}

// (F1) a <= l-d+1; (F2) gap_mid <= 1-(l-1)/N  <=>  G <= 2N - 2(l-1)
inline bool x_orient_ok(const int bx[5]) { return bx[0] <= bx[4] - bx[3] + 1; }
inline bool x_gap_ok(const Tables &T, const int bx[5]) {
  return mid_gap2N(T.N, bx) <= 2 * T.N - 2 * (bx[4] - 1);
}
// second gap inequality: gap_mid(Y sorted copy) >= 1-(l+2)/N  <=>  G >= 2N - 2(l+2)
inline bool y_gap_ok(const Tables &T, int l, int shared, const int tail[4]) {
  int q[5] = {shared, tail[0], tail[1], tail[2], tail[3]};
  std::sort(q, q + 5);
  return mid_gap2N(T.N, q) >= 2 * T.N - 2 * (l + 2);
}

// ---- extend: shared-coordinate interval (seed:clip-interval), outward ---------
struct Ivl { double lo, hi; };   // rigorous bounds of an interval of h_t values

// t = position of the shared entry.  Directions: J_t itself (exact double direction mid(J_t))
// and the perpendiculars to the other four columns.  Starts from [-1/2,1/2].
inline Ivl extend(const Tables &T, const ListData &d, int t) {
  Ivl I{-0.5, 0.5};
  const int N = T.N;
  const Interval Nn((double)N);
  const int bt = d.b[t];
  for (int dir = 0; dir < 5; ++dir) {
    double v0, v1, nv;
    Interval a;
    const double *absrow;
    int jd = -1;
    if (dir == 0) {
      v0 = T.sv0[bt]; v1 = T.sv1[bt]; nv = T.snv[bt];
      a = T.SP[(size_t)bt * N + bt]; absrow = &T.SPa[(size_t)bt * N];
    } else {
      jd = dir - 1; if (jd >= t) ++jd;
      int bj = d.b[jd];
      v0 = T.pv0[bj]; v1 = T.pv1[bj]; nv = T.pnv[bj];
      a = T.PP[(size_t)bj * N + bt]; absrow = &T.PPa[(size_t)bj * N];
    }
    double amin = abs_lo(a);
    if (!(amin >= EXT_MIN_DIV)) continue;           // skipping a constraint is always sound
    // r^up = E ||v|| + (1/2) sum_{j != t} |v.J_j|
    double S = 0.0;
    for (int j = 0; j < 5; ++j) if (j != t) S += absrow[d.b[j]];        // up
    double r = T.Eg_hi * nv + 0.5 * S;                                 // up
#ifdef DIAG_PAPER_MARGINS   // diagnostic only (explains the count delta); never in a certificate build
    r += 1e-7;
#endif
    Interval vc = mulp(d.c0, v0) + mulp(d.c1, v1);
    Interval ctr = -((Nn * vc) / a);                                   // encloses -N v.c/a
    double rad = r / amin;                                             // up
    double lo = -(rad - ctr.lo);                                       // down(ctr.lo - rad)
    double hi = ctr.hi + rad;                                          // up
#ifdef DIAG_PAPER_MARGINS
    lo -= 1e-5; hi += 1e-5;
#endif
    if (lo > I.lo) I.lo = lo;
    if (hi < I.hi) I.hi = hi;
  }
  return I;
}

// reject only if certainly (strictly) disjoint
inline bool ivl_disjoint(const Ivl &a, const Ivl &b) { return a.lo > b.hi || b.lo > a.hi; }

}  // namespace seedcov
