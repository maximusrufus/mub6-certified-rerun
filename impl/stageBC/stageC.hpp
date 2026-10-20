// Stage C: vector-bin cover and graph exclusion (SPEC section 5; exclusion.tex).
// EVERY exclusion inequality is evaluated in interval arithmetic (core/interval.hpp, outward rounding) from the exact
// stored doubles; a test excludes only when the certified lower bound of the left side exceeds the certified upper
// bound of the right side.  The paper's additive margins (1e-12, 2e-7, 3e-7, ...) are therefore NOT used (they exist
// only to cover floating-point error, which the interval layer accounts for - see DECISIONS.md, A23).
// All functions run inside an iv::RoundUp guard.  Directions y, clustering, vertex order are heuristic proposals.
#pragma once
#include "common.hpp"
#include "graph.hpp"
#include <memory>

namespace mb {

// ------------------------------------------------------------------ per-record context
struct CRec {
  const Record* R = nullptr;
  int k = 0;
  Interval r, rH, nhat;  // r (exact), sqrt(1+ort) r, n0 + id r + C r^2
  CInterval Fneg[6][6];  // q_ji = e(u_j) * Fneg[j][i]; Fneg[j][i] = e(-P_ji)/sqrt6
  CInterval PhaseCis[6][6];
  CRec() {}
};

inline void initCRec(CRec& C, const Record& R) {
  const Consts& c = cn();
  C.R = &R;
  C.k = R.k;
  C.r = Interval(R.r);
  C.rH = iv::sqrt(Interval(1.0) + Interval(R.ort)) * C.r;
  C.nhat = Interval(R.n0) + Interval(R.id) * C.r + Interval(R.C) * iv::sqr(C.r);
  for (int j = 0; j < 6; ++j)
    for (int i = 0; i < 6; ++i) C.Fneg[j][i] = iv::cis_turns(-Interval(R.p[6 * j + i])) * c.isqrt6;
}

// ------------------------------------------------------------------ single-vector data
struct Single {
  CInterval q[6][6];  // q[j][i]
  CInterval S[6];
  Interval a2[6], absS[6], c[6], J[6][6];  // J[i][j]
  Interval Smax;
  Interval B[6][10];  // B[i][k]
  bool haveCJ = false, haveB = false;
};

inline void fillS(Single& s) {
  for (int i = 0; i < 6; ++i) {
    CInterval t = s.q[0][i];
    for (int j = 1; j < 6; ++j) t = t + s.q[j][i];
    s.S[i] = t;
    s.a2[i] = iv::abs2(t);
  }
  s.haveCJ = false;
  s.haveB = false;
}
inline void fillCJ(Single& s) {
  const Consts& cc = cn();
  s.Smax = Interval(0.0);
  for (int i = 0; i < 6; ++i) {
    s.absS[i] = iv::sqrt(s.a2[i]);
    s.Smax = imax(s.Smax, s.absS[i]);
    s.c[i] = (s.a2[i] - Interval(1.0)) / cc.c4pi;
    for (int j = 0; j < 6; ++j) s.J[i][j] = s.S[i].im * s.q[j][i].re - s.S[i].re * s.q[j][i].im;  // -Im(conj(S_i) q_ji)
  }
  s.haveCJ = true;
}
inline void fillB(Single& s, const CRec& C) {
  for (int i = 0; i < 6; ++i)
    for (int k = 0; k < C.k; ++k) {
      Interval t(0.0);
      for (int j = 0; j < 6; ++j) t = t + s.J[i][j] * Interval(C.R->col(k)[6 * j + i]);
      s.B[i][k] = t;
    }
  s.haveB = true;
}

// Direct tests of exclude:direct-tests (strictly violated => excluded).  a,X,Lm are the bin (or ball) constants.
inline bool directExcluded(const CRec& C, const Single& s, const Interval& Lm, const Interval& X) {
  const Consts& cc = cn();
  Interval tau = cc.twopi * (Lm * cc.isqrt6 + C.r);
  Interval up = iv::sqr(Interval(1.0) + tau);
  Interval base = Interval(1.0) - tau;
  double lowThr = 0;
  if (base.lo > 0) lowThr = iv::sqr(Interval(base.lo)).lo;
  for (int i = 0; i < 6; ++i) {
    if (s.a2[i].lo > up.hi) return true;
    if (lowThr > 0 && s.a2[i].hi < lowThr) return true;
  }
  Interval agg(0.0);
  for (int i = 0; i < 6; ++i) agg = agg + iv::sqr(iv::sqrt(s.a2[i]) - Interval(1.0));
  Interval thr = iv::sqr(cc.twopi * (cc.K * X + C.r));
  return agg.lo > thr.hi;
}

// scalar linear test (exclude:scalar-linear), needs fillCJ.
inline bool scalarExcluded(const CRec& C, const Single& s, const Interval& a, const Interval& Lm, const Interval& X, const Interval& I) {
  const Consts& cc = cn();
  Interval KX = cc.K * X;
  Interval rho1 = cc.pi * (iv::sqr(KX + C.r) + s.Smax * (cc.K * I * X + (Interval(2.0) * X * C.r + iv::sqr(C.r)) * cc.isqrt6));
  Interval rho2 = cc.pi * (iv::sqr(Lm * cc.isqrt6 + C.r) + s.Smax * iv::sqr(X + C.r) * cc.isqrt6);
  Interval rho = imin(rho1, rho2);
  for (int i = 0; i < 6; ++i) {
    Interval sum1(0.0);
    for (int j = 1; j < 6; ++j) sum1 = sum1 + iv::abs(s.J[i][j]);
    Interval nrm = inorm(s.J[i], 6);
    Interval rhs = a * sum1 + C.r * nrm + rho;
    if (iv::abs(s.c[i]).lo > rhs.hi) return true;
  }
  return false;
}

struct SupportOut {
  Interval lhs, rhs;
  bool excluded() const { return lhs.lo > rhs.hi; }
};

// pure-single remainder  E(y; X, I)  (exclude:weighted-error).  useMin: bin version (min of two bounds); else the ball version
// (first bound only).  Returns the interval of  n N_y + pi[...].
inline Interval weightedError(const CRec& C, const Single& s, const double* y, const Interval& X, const Interval& I, bool useMin) {
  const Consts& cc = cn();
  Interval yinf(0.0), Sy2(0.0), Syinf(0.0), Ny2(0.0);
  for (int i = 0; i < 6; ++i) {
    Interval yi(y[i]);
    yinf = imax(yinf, iv::abs(yi));
    Sy2 = Sy2 + iv::sqr(yi) * s.a2[i];
    Syinf = imax(Syinf, iv::abs(yi) * s.absS[i]);
    for (int j = 0; j < 6; ++j) Ny2 = Ny2 + iv::sqr(yi) * iv::sqr(s.J[i][j]);
  }
  Interval Sy = iv::sqrt(Sy2), Ny = iv::sqrt(Ny2);
  Interval Gy(0.0);
  for (int j = 0; j < 6; ++j) {
    CInterval t(Interval(0.0), Interval(0.0));
    for (int i = 0; i < 6; ++i) t = t + iv::conj(s.S[i]) * s.q[j][i] * Interval(y[i]);
    Gy = imax(Gy, iv::sqrt(iv::abs2(t)));
  }
  Interval gx2 = Gy * iv::sqr(X);
  Interval pureX = gx2;
  if (useMin) pureX = imin(gx2, cc.K * Sy * I * X);
  Interval br = yinf * iv::sqr(cc.K * X + C.r) + pureX + Interval(2.0) * cc.isqrt6 * Sy * X * C.r + cc.isqrt6 * Syinf * iv::sqr(C.r);
  return C.nhat * Ny + cc.pi * br;
}

// bin support inequality (exclude:bin-support)
inline SupportOut binSupport(const CRec& C, const Single& s, const double* y, const Interval& a, const Interval& X, const Interval& I) {
  Interval yc(0.0);
  for (int i = 0; i < 6; ++i) yc = yc + Interval(y[i]) * s.c[i];
  Interval sumg(0.0);
  for (int j = 1; j < 6; ++j) {
    Interval g(0.0);
    for (int i = 0; i < 6; ++i) g = g + Interval(y[i]) * s.J[i][j];
    sumg = sumg + iv::abs(g);
  }
  Interval bv[10];
  for (int k = 0; k < C.k; ++k) {
    Interval t(0.0);
    for (int i = 0; i < 6; ++i) t = t + Interval(y[i]) * s.B[i][k];
    bv[k] = t;
  }
  Interval nB = inorm(bv, C.k);
  SupportOut o;
  o.lhs = iv::abs(yc);
  o.rhs = a * sumg + C.rH * nB + weightedError(C, s, y, X, I, true);
  return o;
}

// -------- heuristic directions
// Solve  H w = g  for the quadratic surrogate of the support function, normalising to max |entry| = 1.
inline bool normalizeDir(std::vector<double>& y) {
  double mx = 0;
  for (double v : y) {
    if (!std::isfinite(v)) return false;
    mx = std::max(mx, std::fabs(v));
  }
  if (!(mx > 1e-100 && mx < 1e100)) return false;
  for (double& v : y) v /= mx;
  return true;
}

// proposals for the single-vector directional test: L2 surrogate then a few IRLS passes
inline void proposeSingle(const CRec& C, const Single& s, double a, std::vector<std::vector<double>>& cands) {
  double Jm[6][6], Bm[6][10], cm[6];
  for (int i = 0; i < 6; ++i) {
    cm[i] = s.c[i].mid();
    for (int j = 0; j < 6; ++j) Jm[i][j] = s.J[i][j].mid();
    for (int k = 0; k < C.k; ++k) Bm[i][k] = s.B[i][k].mid();
  }
  double rH = C.rH.hi;
  std::vector<double> wj(6, a * a), wk(std::max(C.k, 1), rH * rH);
  std::vector<double> y(6);
  for (int pass = 0; pass < 5; ++pass) {
    std::vector<double> H(36, 0.0);
    for (int j = 1; j < 6; ++j)
      for (int p = 0; p < 6; ++p)
        for (int q = 0; q < 6; ++q) H[p * 6 + q] += wj[j] * Jm[p][j] * Jm[q][j];
    for (int k = 0; k < C.k; ++k)
      for (int p = 0; p < 6; ++p)
        for (int q = 0; q < 6; ++q) H[p * 6 + q] += wk[k] * Bm[p][k] * Bm[q][k];
    double tr = 0;
    for (int p = 0; p < 6; ++p) tr += H[p * 6 + p];
    for (int p = 0; p < 6; ++p) H[p * 6 + p] += 1e-9 * tr + 1e-18;
    std::vector<double> b(cm, cm + 6);
    if (!solveLin(H, b, 6, y)) break;
    std::vector<double> yn = y;
    if (!normalizeDir(yn)) break;
    cands.push_back(yn);
    // IRLS weights: a|g_j| ~ a g_j^2 / |g_j|
    for (int j = 1; j < 6; ++j) {
      double g = 0;
      for (int i = 0; i < 6; ++i) g += yn[i] * Jm[i][j];
      wj[j] = a / std::max(std::fabs(g), 1e-6);
    }
    if (C.k) {
      for (int k = 0; k < C.k; ++k) {
        double t = 0;
        for (int i = 0; i < 6; ++i) t += yn[i] * Bm[i][k];
        (void)t;
      }
      double nb = 0;
      for (int k = 0; k < C.k; ++k) {
        double t = 0;
        for (int i = 0; i < 6; ++i) t += yn[i] * Bm[i][k];
        nb += t * t;
      }
      nb = std::sqrt(nb);
      for (int k = 0; k < C.k; ++k) wk[k] = rH / std::max(nb, 1e-6);
    }
  }
  {  // plain residual direction
    std::vector<double> yc(cm, cm + 6);
    if (normalizeDir(yc)) cands.push_back(yc);
  }
}

// ------------------------------------------------------------------ bin layer
struct Layer {
  int M = 0;
  Interval a, X, Lm, I;
  std::vector<CInterval> tq;  // [j=1..5][b][i] flattened: ((j-1)*M + b)*6 + i
  CInterval q0[6];
  void init(const CRec& C, int M_) {
    const Consts& cc = cn();
    M = M_;
    a = Interval(0.5) / Interval((double)M);
    X = cc.sqrt29_6 * a;
    Lm = Interval(5.0) * a;
    I = Interval(1.5) * a;
    for (int i = 0; i < 6; ++i) q0[i] = C.Fneg[0][i];
    tq.assign((size_t)5 * M * 6, CInterval());
    for (int j = 1; j < 6; ++j)
      for (int b = 0; b < M; ++b) {
        Interval u = Interval((double)(2 * b + 1)) / Interval((double)(2 * M));
        for (int i = 0; i < 6; ++i) tq[((size_t)(j - 1) * M + b) * 6 + i] = iv::cis_turns(u - Interval(C.R->p[6 * j + i])) * cc.isqrt6;
      }
  }
  const CInterval* row(int j, int b) const { return j == 0 ? q0 : &tq[((size_t)(j - 1) * M + b) * 6]; }
};

// load q rows for bin b[0..4] (coordinates j=1..5)
inline void loadBin(const Layer& L, const uint8_t* b, Single& s) {
  for (int i = 0; i < 6; ++i) s.q[0][i] = L.q0[i];
  for (int j = 1; j < 6; ++j) {
    const CInterval* r = L.row(j, b[j - 1]);
    for (int i = 0; i < 6; ++i) s.q[j][i] = r[i];
  }
}

// true iff the bin is certified to contain no unbiased vector for any admissible Hadamard displacement
inline bool binExcluded(const CRec& C, const Layer& L, const uint8_t* b, Single& s, bool dbgNoDirectional = false) {
  loadBin(L, b, s);
  fillS(s);
  if (directExcluded(C, s, L.Lm, L.X)) return true;
  fillCJ(s);
  if (scalarExcluded(C, s, L.a, L.Lm, L.X, L.I)) return true;
  if (dbgNoDirectional) return false;
  fillB(s, C);
  std::vector<std::vector<double>> cands;
  proposeSingle(C, s, L.a.hi, cands);
  for (auto& y : cands)
    if (binSupport(C, s, y.data(), L.a, L.X, L.I).excluded()) return true;
  return false;
}

// ------------------------------------------------------------------ clustering
struct Ball {
  double v[6];
  double s = 0;  // certified upper bound for the radius
  CInterval E[6];  // e(v_j)
  Single single;
  bool ready = false;
};

struct ClusterOut {
  bool ok = false;
  std::string why;
  std::vector<std::array<double, 6>> v;
  std::vector<double> s;
};

inline ClusterOut clusterBins(const std::vector<std::array<uint8_t, 5>>& bins, int M, double radiusParam = 0.042, int maxCenters = 2048) {
  const Consts& cc = cn();
  ClusterOut out;
  const size_t n = bins.size();
  std::vector<std::array<double, 6>> U(n);
  for (size_t t = 0; t < n; ++t) {
    U[t][0] = 0;
    for (int j = 1; j < 6; ++j) U[t][j] = (2.0 * bins[t][j - 1] + 1.0) / (2.0 * M);
  }
  double thr = (radiusParam - 1.11 / M);
  thr *= thr;
  std::vector<std::array<double, 6>> V;
  std::vector<int> cid(n, -1);
  std::vector<std::array<int, 6>> mv(n);
  for (int pass = 0; pass < 2; ++pass) {
    std::vector<std::array<double, 6>> sum(V.size());
    std::vector<int> cnt(V.size(), 0);
    for (auto& s : sum) s.fill(0.0);
    for (size_t t = 0; t < n; ++t) {
      int bestc = -1;
      double bestd = 1e300;
      std::array<int, 6> bm{};
      for (size_t c = 0; c < V.size(); ++c) {
        double x[6];
        std::array<int, 6> m{};
        for (int j = 0; j < 6; ++j) {
          double d = U[t][j] - V[c][j];
          double mm = rnd(d);
          m[j] = (int)mm;
          x[j] = d - mm;
        }
        double mean = 0;
        for (int j = 0; j < 6; ++j) mean += x[j];
        mean /= 6;
        double d2 = 0;
        for (int j = 0; j < 6; ++j) d2 += (x[j] - mean) * (x[j] - mean);
        if (d2 <= thr && d2 < bestd) { bestd = d2; bestc = (int)c; bm = m; }
      }
      if (bestc < 0) {
        if ((int)V.size() >= maxCenters) { out.why = "more than 2048 centers"; return out; }
        V.push_back(U[t]);
        sum.push_back(std::array<double, 6>{0, 0, 0, 0, 0, 0});
        cnt.push_back(0);
        bestc = (int)V.size() - 1;
        bm = {0, 0, 0, 0, 0, 0};
      }
      cid[t] = bestc;
      mv[t] = bm;
      ++cnt[bestc];
      for (int j = 0; j < 6; ++j) sum[bestc][j] += U[t][j] - bm[j];
    }
    if (pass == 0) {  // move centres to centroids; the last pass does not move them
      for (size_t c = 0; c < V.size(); ++c)
        if (cnt[c]) for (int j = 0; j < 6; ++j) V[c][j] = sum[c][j] / cnt[c];
    } else {
      // keep only non-empty centres
      std::vector<int> remap(V.size(), -1);
      std::vector<std::array<double, 6>> V2;
      for (size_t c = 0; c < V.size(); ++c)
        if (cnt[c]) { remap[c] = (int)V2.size(); V2.push_back(V[c]); }
      for (size_t t = 0; t < n; ++t) cid[t] = remap[cid[t]];
      V = V2;
    }
  }
  // certified radii: s_c >= max ||h|| + X_M  (exclude:cluster-cover)
  Interval XM = cc.sqrt29_6 * (Interval(0.5) / Interval((double)M));
  std::vector<double> rad(V.size(), 0.0);
  for (size_t t = 0; t < n; ++t) {
    int c = cid[t];
    Interval h[6];
    Interval mean(0.0);
    for (int j = 0; j < 6; ++j) {
      Interval uj = (j == 0) ? Interval(0.0) : Interval((double)(2 * bins[t][j - 1] + 1)) / Interval((double)(2 * M));
      h[j] = uj - Interval(V[c][j]) - Interval((double)mv[t][j]);
      mean = mean + h[j];
    }
    mean = mean * cc.sixth;
    for (int j = 0; j < 6; ++j) h[j] = h[j] - mean;
    double sr = (inorm(h, 6) + XM).hi;
    rad[c] = std::max(rad[c], sr);
  }
  for (size_t c = 0; c < V.size(); ++c) {
    for (int j = 0; j < 6; ++j)
      if (!(std::fabs(V[c][j]) <= 2)) { out.why = "|v_j| > 2"; return out; }
    if (!(rad[c] < 0.1)) { out.why = "cluster radius >= .1"; return out; }
  }
  out.v = V;
  out.s = rad;
  out.ok = true;
  return out;
}

// ------------------------------------------------------------------ balls and edges
inline void makeBall(const CRec& C, const double* v, double s, Ball& b) {
  const Consts& cc = cn();
  for (int j = 0; j < 6; ++j) {
    b.v[j] = v[j];
    b.E[j] = iv::cis_turns(Interval(v[j]));
    for (int i = 0; i < 6; ++i) b.single.q[j][i] = iv::cis_turns(Interval(v[j]) - Interval(C.R->p[6 * j + i])) * cc.isqrt6;
  }
  b.s = s;
  fillS(b.single);
  fillCJ(b.single);
  fillB(b.single, C);
  b.ready = true;
}

struct EdgeData {
  CInterval Q[6], S;
  Interval a2, absS, re;
  Interval DO[2][6], DU[6], nDU;
};

inline void edgeData(const Ball& A, const Ball& B, EdgeData& e) {
  const Consts& cc = cn();
  e.S = CInterval(Interval(0.0), Interval(0.0));
  for (int j = 0; j < 6; ++j) {
    e.Q[j] = A.E[j] * iv::conj(B.E[j]);
    e.S = e.S + e.Q[j];
  }
  e.a2 = iv::abs2(e.S);
  e.absS = iv::sqrt(e.a2);
  e.re = Interval(A.s) + Interval(B.s);
  for (int j = 0; j < 6; ++j) {
    e.DO[0][j] = -e.Q[j].im;
    e.DO[1][j] = e.Q[j].re;
    e.DU[j] = (e.S.re * e.DO[0][j] + e.S.im * e.DO[1][j]) * cc.sixth;
  }
  e.nDU = inorm(e.DU, 6);
}

// ball-level pure single remainder for the joint test (exclude:ball-error); the maximum coefficient is formed from the
// pure single data (before edge derivatives are added), as the paper requires.
inline Interval ballError(const CRec& C, const Ball& b, const double* y) {
  Interval X(b.s);
  Interval I(0.0);  // unused for the ball version
  return weightedError(C, b.single, y, X, I, false);
}

// joint inequality (exclude:joint-inequality) for relation `isO`.  y has 12 + l entries.
inline SupportOut jointSupport(const CRec& C, const Ball& A, const Ball& B, const EdgeData& e, bool isO, const double* y) {
  const Consts& cc = cn();
  const int l = isO ? 2 : 1;
  const double* y0 = y;
  const double* y1 = y + 6;
  const double* ye = y + 12;
  Interval lhs(0.0);
  for (int i = 0; i < 6; ++i) lhs = lhs + Interval(y0[i]) * A.single.c[i] + Interval(y1[i]) * B.single.c[i];
  Interval EE, nye;
  if (isO) {
    Interval q0 = e.S.re / cc.twopi, q1 = e.S.im / cc.twopi;
    lhs = lhs + Interval(ye[0]) * q0 + Interval(ye[1]) * q1;
    EE = cc.pi * iv::sqr(e.re);
    Interval t[2] = {Interval(ye[0]), Interval(ye[1])};
    nye = inorm(t, 2);
  } else {
    Interval qU = (e.a2 * cc.sixth - Interval(1.0)) / cc.c4pi;
    lhs = lhs + Interval(ye[0]) * qU;
    EE = cc.pi * iv::sqr(e.re) * (Interval(1.0) + e.absS * cc.sixth);
    nye = iv::abs(Interval(ye[0]));
  }
  // edge terms
  Interval t0[6], t1[6];
  for (int j = 0; j < 6; ++j) {
    Interval g0(0.0), g1(0.0), de(0.0);
    for (int i = 0; i < 6; ++i) { g0 = g0 + Interval(y0[i]) * A.single.J[i][j]; g1 = g1 + Interval(y1[i]) * B.single.J[i][j]; }
    if (isO) de = Interval(ye[0]) * e.DO[0][j] + Interval(ye[1]) * e.DO[1][j];
    else de = Interval(ye[0]) * e.DU[j];
    t0[j] = g0 + de;
    t1[j] = g1 - de;
  }
  Interval rhs = Interval(A.s) * inorm(t0, 6) + Interval(B.s) * inorm(t1, 6);
  rhs = rhs + ballError(C, A, y0) + ballError(C, B, y1);
  Interval bv[10];
  for (int k = 0; k < C.k; ++k) {
    Interval t(0.0);
    for (int i = 0; i < 6; ++i) t = t + Interval(y0[i]) * A.single.B[i][k] + Interval(y1[i]) * B.single.B[i][k];
    bv[k] = t;
  }
  rhs = rhs + C.rH * inorm(bv, C.k) + EE * nye;
  (void)l;
  SupportOut o;
  o.lhs = iv::abs(lhs);
  o.rhs = rhs;
  return o;
}

// heuristic joint direction (one proposal; second attempt with IRLS-like reweighting)
inline void proposeJoint(const CRec& C, const Ball& A, const Ball& B, const EdgeData& e, bool isO, std::vector<std::vector<double>>& cands) {
  const int l = isO ? 2 : 1;
  const int n = 12 + l;
  double J0[6][6], J1[6][6], B0[6][10], B1[6][10], c0[6], c1[6], DE[2][6];
  for (int i = 0; i < 6; ++i) {
    c0[i] = A.single.c[i].mid();
    c1[i] = B.single.c[i].mid();
    for (int j = 0; j < 6; ++j) { J0[i][j] = A.single.J[i][j].mid(); J1[i][j] = B.single.J[i][j].mid(); }
    for (int k = 0; k < C.k; ++k) { B0[i][k] = A.single.B[i][k].mid(); B1[i][k] = B.single.B[i][k].mid(); }
  }
  for (int j = 0; j < 6; ++j) {
    if (isO) { DE[0][j] = e.DO[0][j].mid(); DE[1][j] = e.DO[1][j].mid(); }
    else { DE[0][j] = e.DU[j].mid(); DE[1][j] = 0; }
  }
  double s0 = A.s, s1 = B.s, rH = C.rH.hi;
  const Consts& cc = cn();
  double EE = (isO ? (cc.pi * iv::sqr(e.re)) : (cc.pi * iv::sqr(e.re) * (Interval(1.0) + e.absS * cc.sixth))).hi;
  double qE[2];
  if (isO) { qE[0] = (e.S.re / cc.twopi).mid(); qE[1] = (e.S.im / cc.twopi).mid(); }
  else { qE[0] = ((e.a2 * cc.sixth - Interval(1.0)) / cc.c4pi).mid(); qE[1] = 0; }
  std::vector<double> g(n);
  for (int i = 0; i < 6; ++i) { g[i] = c0[i]; g[6 + i] = c1[i]; }
  for (int t = 0; t < l; ++t) g[12 + t] = qE[t];
  double wnorm0 = 1, wnorm1 = 1, wB = 1;
  for (int pass = 0; pass < 4; ++pass) {
    std::vector<double> H((size_t)n * n, 0.0);
    auto addRow = [&](const std::vector<double>& row, double wt) {
      for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) H[(size_t)p * n + q] += wt * row[p] * row[q];
    };
    for (int j = 0; j < 6; ++j) {
      std::vector<double> r0(n, 0.0), r1(n, 0.0);
      for (int i = 0; i < 6; ++i) { r0[i] = J0[i][j]; r1[6 + i] = J1[i][j]; }
      for (int t = 0; t < l; ++t) { r0[12 + t] = DE[t][j]; r1[12 + t] = -DE[t][j]; }
      addRow(r0, s0 * s0 * wnorm0);
      addRow(r1, s1 * s1 * wnorm1);
    }
    for (int k = 0; k < C.k; ++k) {
      std::vector<double> r(n, 0.0);
      for (int i = 0; i < 6; ++i) { r[i] = B0[i][k]; r[6 + i] = B1[i][k]; }
      addRow(r, rH * rH * wB);
    }
    for (int t = 0; t < l; ++t) {
      std::vector<double> r(n, 0.0);
      r[12 + t] = 1;
      addRow(r, EE * EE);
    }
    double tr = 0;
    for (int p = 0; p < n; ++p) tr += H[(size_t)p * n + p];
    for (int p = 0; p < n; ++p) H[(size_t)p * n + p] += 1e-9 * tr + 1e-18;
    std::vector<double> y;
    if (!solveLin(H, g, n, y)) return;
    if (!normalizeDir(y)) return;
    cands.push_back(y);
    // reweight: L2 surrogate of a norm  s*||v|| ~ s*(||v||^2/||v_old|| + ||v_old||)/2
    double n0 = 0, n1 = 0, nb = 0;
    for (int j = 0; j < 6; ++j) {
      double a0 = 0, a1 = 0;
      for (int i = 0; i < 6; ++i) { a0 += y[i] * J0[i][j]; a1 += y[6 + i] * J1[i][j]; }
      for (int t = 0; t < l; ++t) { a0 += y[12 + t] * DE[t][j]; a1 -= y[12 + t] * DE[t][j]; }
      n0 += a0 * a0; n1 += a1 * a1;
    }
    for (int k = 0; k < C.k; ++k) {
      double a = 0;
      for (int i = 0; i < 6; ++i) a += y[i] * B0[i][k] + y[6 + i] * B1[i][k];
      nb += a * a;
    }
    wnorm0 = 1.0 / std::max(std::sqrt(n0) * s0, 1e-9) * (s0 > 0 ? 1.0 : 0.0) + 1e-12;
    wnorm1 = 1.0 / std::max(std::sqrt(n1) * s1, 1e-9) * (s1 > 0 ? 1.0 : 0.0) + 1e-12;
    wB = 1.0 / std::max(std::sqrt(nb) * rH, 1e-9);
    wnorm0 *= s0 > 0 ? 1.0 : 1.0;
  }
}

struct EdgeRes { bool O = false, U = false; };

// An edge is kept unless some test is strictly violated.  Returns which relations remain possible.
inline EdgeRes testEdge(const CRec& C, const Ball& A, const Ball& B, bool useJoint = true) {
  const Consts& cc = cn();
  EdgeData e;
  edgeData(A, B, e);
  EdgeRes res;
  res.O = true;
  res.U = true;
  // distance tests
  if (((e.absS / (cc.twopi * cc.sqrt6)).lo) > e.re.hi) res.O = false;
  if ((iv::abs(e.absS * cc.isqrt6 - Interval(1.0)) / cc.twopi).lo > e.re.hi) res.U = false;
  if (!res.O && !res.U) return res;
  if (res.O) {
    Interval EO = cc.pi * iv::sqr(e.re);
    Interval lhs = e.a2 / cc.twopi;
    Interval rhs = Interval(6.0) * e.re * e.nDU + EO * e.absS;
    if (lhs.lo > rhs.hi) res.O = false;
  }
  if (res.U) {
    Interval qU = (e.a2 * cc.sixth - Interval(1.0)) / cc.c4pi;
    Interval EU = cc.pi * iv::sqr(e.re) * (Interval(1.0) + e.absS * cc.sixth);
    Interval rhs = e.re * e.nDU + EU;
    if (iv::abs(qU).lo > rhs.hi) res.U = false;
  }
  if (!useJoint) return res;
  for (int rel = 0; rel < 2; ++rel) {
    bool isO = rel == 0;
    if ((isO && !res.O) || (!isO && !res.U)) continue;
    std::vector<std::vector<double>> cands;
    proposeJoint(C, A, B, e, isO, cands);
    for (auto& y : cands) {
      if (jointSupport(C, A, B, e, isO, y.data()).excluded()) {
        if (isO) res.O = false; else res.U = false;
        break;
      }
    }
  }
  return res;
}

// ------------------------------------------------------------------ record driver
enum class COutcome { NegEmptyLayer, NegGraph, Overflow, ClusterFail, GraphCap, GraphCandidate, TripleFound, TripleAbsent };
inline const char* outcomeName(COutcome o) {
  switch (o) {
    case COutcome::NegEmptyLayer: return "NEGATIVE(empty-layer)";
    case COutcome::NegGraph: return "NEGATIVE(graph)";
    case COutcome::Overflow: return "UNRESOLVED(bin-overflow)";
    case COutcome::ClusterFail: return "UNRESOLVED(cluster-failure)";
    case COutcome::GraphCap: return "UNRESOLVED(graph-cap)";
    case COutcome::GraphCandidate: return "UNRESOLVED(graph-candidate)";
    case COutcome::TripleFound: return "TRIPLE-FOUND";
    case COutcome::TripleAbsent: return "TRIPLE-ABSENT";
  }
  return "?";
}
inline bool isNegative(COutcome o) { return o == COutcome::NegEmptyLayer || o == COutcome::NegGraph; }

struct CParams {
  int M0 = 12, layers = 4;
  double radiusParam = 0.042;
  uint64_t graphCap = 5000000;
  size_t binCap = 2000000;
  bool tripleMode = false;  // test hook: look for ONE O-6-clique (single further basis) instead of the 2-clique pattern
};

struct CStats {
  std::vector<size_t> layerBins;
  size_t balls = 0, edgesO = 0, edgesU = 0, pruned = 0;
  uint64_t steps = 0;
  std::string note;
};

inline COutcome processRecordC(const Record& R, const CParams& P, CStats& st) {
  CRec C;
  initCRec(C, R);
  std::vector<std::array<uint8_t, 5>> cur, nxt;
  Single s;
  int M = P.M0;
  for (int layer = 0; layer < P.layers; ++layer, M *= 2) {
    Layer L;
    L.init(C, M);
    nxt.clear();
    if (layer == 0) {
      uint8_t b[5];
      for (b[0] = 0; b[0] < M; ++b[0])
        for (b[1] = 0; b[1] < M; ++b[1])
          for (b[2] = 0; b[2] < M; ++b[2])
            for (b[3] = 0; b[3] < M; ++b[3])
              for (b[4] = 0; b[4] < M; ++b[4])
                if (!binExcluded(C, L, b, s)) nxt.push_back({b[0], b[1], b[2], b[3], b[4]});
    } else {
      for (auto& par : cur)
        for (int ch = 0; ch < 32; ++ch) {
          uint8_t b[5];
          for (int j = 0; j < 5; ++j) b[j] = (uint8_t)(2 * par[j] + ((ch >> j) & 1));
          if (!binExcluded(C, L, b, s)) {
            nxt.push_back({b[0], b[1], b[2], b[3], b[4]});
            if (nxt.size() > P.binCap) { st.layerBins.push_back(nxt.size()); return COutcome::Overflow; }
          }
        }
    }
    st.layerBins.push_back(nxt.size());
    if (nxt.empty()) return COutcome::NegEmptyLayer;
    cur.swap(nxt);
  }
  M /= 2;  // final layer resolution
  ClusterOut cl = clusterBins(cur, M, P.radiusParam);
  if (!cl.ok) { st.note = cl.why; return COutcome::ClusterFail; }
  const int nb = (int)cl.v.size();
  st.balls = nb;
  std::vector<std::unique_ptr<Ball>> balls(nb);
  for (int t = 0; t < nb; ++t) {
    balls[t].reset(new Ball());
    makeBall(C, cl.v[t].data(), cl.s[t], *balls[t]);
  }
  Graph G;
  G.init(nb);
  for (int a = 0; a < nb; ++a)
    for (int b = a + 1; b < nb; ++b) {
      EdgeRes er = testEdge(C, *balls[a], *balls[b]);
      if (er.O) { G.addO(a, b); ++st.edgesO; }
      if (er.U) { G.addU(a, b); ++st.edgesU; }
    }
  if (P.tripleMode) {
    GraphSearch gs(G, P.graphCap);
    bool f = gs.hasO6();
    st.steps = gs.steps;
    if (gs.gaveUp) return COutcome::TripleFound;  // possibility
    return f ? COutcome::TripleFound : COutcome::TripleAbsent;
  }
  std::vector<int> keep = pruneVertices(G);
  st.pruned = nb - keep.size();
  Graph H = induced(G, keep);
  GraphSearch gs(H, P.graphCap);
  bool found = gs.exist();
  st.steps = gs.steps;
  if (gs.gaveUp) return COutcome::GraphCap;
  if (found) return COutcome::GraphCandidate;
  return COutcome::NegGraph;
}

}  // namespace mb
