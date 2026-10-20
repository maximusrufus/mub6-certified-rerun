// phi.hpp -- symmetry maps Phi for Stage B (SPEC 1.4 phase:transport, 3.8, 4.5 refine:atlas-assignment,
// refine:mapped-coordinates).  A Phi is a signed row/column permutation with optional transpose and optional
// global sign flip: a linear isometry of R^36 commuting with Pi.  By phase:transport it maps a certified Root
// (P, T, consts) to the certified Root (Phi P, Phi T, same consts), and the solution spheres are isometric.
//
// Soundness: ANY Phi is valid; the choice (seed-arrangement signature, as in stageA/atlas.cpp) is a heuristic
// that only decides which merge/reuse attempts are tried.  A box is never declared covered by a Phi alone:
// every use feeds a pulled-back Root  Phi^{-1} Z  into the ordinary interval relation (makeRel) or the ordinary
// aligned-distance merge (assign), which are the certificates.
#pragma once
#include "common.hpp"
#include <algorithm>
#include <vector>

namespace mb {

struct Phi {
  int map[36];  // (Phi x)[i] = sign * x[map[i]]
  int sign = 1;
  static Phi identity() { Phi f; for (int i = 0; i < 36; ++i) f.map[i] = i; f.sign = 1; return f; }
  bool isIdentity() const { if (sign != 1) return false; for (int i = 0; i < 36; ++i) if (map[i] != i) return false; return true; }
  Phi inverse() const { Phi g; g.sign = sign; for (int i = 0; i < 36; ++i) g.map[map[i]] = i; return g; }
};

inline void applyPhiVec(const Phi& f, const double* x, double* y) {
  for (int i = 0; i < 36; ++i) y[i] = f.sign * x[f.map[i]];
}
// Transported Root (phase:transport): P and every T column mapped, constants and radius unchanged.
inline Record applyPhi(const Record& R, const Phi& f) {
  Record O = R;
  applyPhiVec(f, R.p, O.p);
  for (int a = 0; a < R.k; ++a) applyPhiVec(f, R.col(a), O.T.data() + 36 * (size_t)a);
  return O;
}

namespace phidet {
inline double frac(double x) { return x - std::floor(x); }
inline void dephase(const double* P, double* X) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) X[6 * i + j] = frac(P[6 * i + j] - P[6 * i] - P[j] + P[0]);
}
inline double pairGap(const double* d, int n) {
  std::vector<double> v(d, d + n);
  std::sort(v.begin(), v.end());
  double best = -1;
  for (int i = 0; i < n; ++i) {
    double g = (i + 1 < n) ? v[i + 1] - v[i] : v[0] + 1.0 - v[n - 1];
    best = std::max(best, g);
  }
  return best;
}
}  // namespace phidet

// Candidate maps by the seed-arrangement signature (same rule as stageA/atlas.cpp `variants`): over
// transpose x sign x ordered row pairs (a,b), keep pairs whose largest circular gap of the six phase differences
// is within `tol` of the global minimum; columns start after any gap within 1e-3 of the pair's largest gap;
// the remaining four rows sorted by their dephased entry in column position 1.  gmin = the invariant minimum gap.
inline std::vector<Phi> phiVariants(const double* P, double tol, double& gmin) {
  using phidet::frac;
  double X[36];
  phidet::dephase(P, X);
  struct Cand { int tr, sg, a, b; double gap; };
  std::vector<Cand> cs;
  gmin = 1e300;
  for (int tr = 0; tr < 2; ++tr)
    for (int sg = 0; sg < 2; ++sg) {
      double M[36];
      for (int a = 0; a < 6; ++a)
        for (int b = 0; b < 6; ++b) M[6 * a + b] = frac((sg ? -1.0 : 1.0) * (tr ? X[6 * b + a] : X[6 * a + b]));
      for (int a = 0; a < 6; ++a)
        for (int b = 0; b < 6; ++b) {
          if (a == b) continue;
          double d[6];
          for (int c = 0; c < 6; ++c) d[c] = frac(M[6 * b + c] - M[6 * a + c]);
          double g = phidet::pairGap(d, 6);
          cs.push_back({tr, sg, a, b, g});
          gmin = std::min(gmin, g);
        }
    }
  std::vector<Phi> out;
  for (auto& c : cs) {
    if (c.gap > gmin + tol) continue;
    double M[36];
    for (int a = 0; a < 6; ++a)
      for (int b = 0; b < 6; ++b) M[6 * a + b] = frac((c.sg ? -1.0 : 1.0) * (c.tr ? X[6 * b + a] : X[6 * a + b]));
    double d[6];
    for (int cc = 0; cc < 6; ++cc) d[cc] = frac(M[6 * c.b + cc] - M[6 * c.a + cc]);
    std::vector<std::pair<double, int>> v;
    for (int i = 0; i < 6; ++i) v.push_back({d[i], i});
    std::sort(v.begin(), v.end());
    std::vector<int> starts;
    for (int i = 0; i < 6; ++i) {
      double g = (i + 1 < 6) ? v[i + 1].first - v[i].first : v[0].first + 1.0 - v[5].first;
      if (g >= c.gap - 1e-3) starts.push_back(v[(i + 1) % 6].second);
    }
    for (int s0 : starts) {
      int colp[6], rowp[6];
      std::vector<std::pair<double, int>> cols;
      for (int cc = 0; cc < 6; ++cc) cols.push_back({frac(d[cc] - d[s0]), cc});
      std::sort(cols.begin(), cols.end());
      for (int i = 0; i < 6; ++i) colp[i] = cols[i].second;
      rowp[0] = c.a;
      rowp[1] = c.b;
      std::vector<std::pair<double, int>> rows;
      for (int r = 0; r < 6; ++r) {
        if (r == c.a || r == c.b) continue;
        double val = frac(M[6 * r + colp[1]] - M[6 * r + colp[0]] - M[6 * c.a + colp[1]] + M[6 * c.a + colp[0]]);
        rows.push_back({val, r});
      }
      std::sort(rows.begin(), rows.end());
      for (int i = 0; i < 4; ++i) rowp[2 + i] = rows[i].second;
      Phi f;
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) f.map[6 * i + j] = c.tr ? 6 * colp[j] + rowp[i] : 6 * rowp[i] + colp[j];
      f.sign = c.sg ? -1 : 1;
      out.push_back(f);
    }
  }
  return out;
}

// Variants used for reuse / merging: identity first, then the signature maps (duplicates removed).
inline std::vector<Phi> reuseVariants(const double* P, double tol = 0.02) {
  double gmin;
  std::vector<Phi> vs = phiVariants(P, tol, gmin);
  std::vector<Phi> out;
  out.push_back(Phi::identity());
  for (auto& f : vs) {
    bool dup = false;
    for (auto& g : out) {
      if (g.sign != f.sign) continue;
      bool same = true;
      for (int i = 0; i < 36 && same; ++i) same = g.map[i] == f.map[i];
      if (same) { dup = true; break; }
    }
    if (!dup) out.push_back(f);
  }
  return out;
}

// heuristic aligned distance (double) for ranking candidates
inline double approxAligned(const double* P1, const double* P2) {
  double W[36], mr[6] = {0}, mc[6] = {0}, tm = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) {
      double a = P1[6 * i + j] - P1[6 * i] - P1[j] + P1[0];
      double b = P2[6 * i + j] - P2[6 * i] - P2[j] + P2[0];
      double N = rnd(a - b);
      W[6 * i + j] = a - b - N;
      mr[i] += W[6 * i + j] / 6;
      mc[j] += W[6 * i + j] / 6;
      tm += W[6 * i + j] / 36;
    }
  double s = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) {
      double v = W[6 * i + j] - mr[i] - mc[j] + tm;
      s += v * v;
    }
  return std::sqrt(s);
}

}  // namespace mb
