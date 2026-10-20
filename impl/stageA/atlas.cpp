// atlas.cpp -- atlas merge (SPEC 3.8) and the binary Record format (SPEC 7).
// Variants Phi (A13): signed row/column permutations with optional transpose, chosen by the
// seed-arrangement rule (smallest-gap row pair, start after the largest gap, remaining rows
// sorted by their entry in column 1).  Any Phi is sound; the signature is only a heuristic.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include "stagea.hpp"

namespace sa {

Root applyPhi(const Root& R, const Phi& f) {
  Root O = R;
  for (int i = 0; i < 36; ++i) O.P[i] = f.sign * R.P[f.map[i]];
  for (int a = 0; a < R.k; ++a)
    for (int i = 0; i < 36; ++i) O.T[a][i] = f.sign * R.T[a][f.map[i]];
  return O;
}

namespace {

double frac(double x) {
  double y = x - std::floor(x);
  return y;
}

struct Var {
  Phi phi;
  double gap;  // largest gap of the pivot pair
};

// dephased phases (double) of P
void dephase(const double P[36], double X[36]) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) X[6 * i + j] = frac(P[6 * i + j] - P[6 * i] - P[j] + P[0]);
}

double pairGap(const double* d, int n, int* startAfterLargest, double* gapsOut = nullptr) {
  // d: n phases mod 1; return largest circular gap; startAfterLargest = index (in d) of the point just after it
  std::vector<std::pair<double, int>> v;
  for (int i = 0; i < n; ++i) v.push_back({d[i], i});
  std::sort(v.begin(), v.end());
  double best = -1;
  int bi = 0;
  for (int i = 0; i < n; ++i) {
    double g = (i + 1 < n) ? v[i + 1].first - v[i].first : v[0].first + 1.0 - v[n - 1].first;
    if (g > best) {
      best = g;
      bi = (i + 1) % n;
    }
  }
  *startAfterLargest = v[bi].second;
  (void)gapsOut;
  return best;
}

// all variants whose pivot gap is within `tol` of the minimum; gmin returned
std::vector<Var> variants(const double P[36], double tol, double& gmin) {
  double X[36];
  dephase(P, X);
  struct Cand {
    int tr, sg, a, b;
    double gap;
  };
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
          int st;
          double g = pairGap(d, 6, &st);
          cs.push_back({tr, sg, a, b, g});
          gmin = std::min(gmin, g);
        }
    }
  std::vector<Var> out;
  for (auto& c : cs) {
    if (c.gap > gmin + tol) continue;
    double M[36];
    for (int a = 0; a < 6; ++a)
      for (int b = 0; b < 6; ++b) M[6 * a + b] = frac((c.sg ? -1.0 : 1.0) * (c.tr ? X[6 * b + a] : X[6 * a + b]));
    double d[6];
    for (int cc = 0; cc < 6; ++cc) d[cc] = frac(M[6 * c.b + cc] - M[6 * c.a + cc]);
    // start columns: after any gap within 1e-3 of the largest
    std::vector<std::pair<double, int>> v;
    for (int i = 0; i < 6; ++i) v.push_back({d[i], i});
    std::sort(v.begin(), v.end());
    std::vector<int> starts;
    for (int i = 0; i < 6; ++i) {
      double g = (i + 1 < 6) ? v[i + 1].first - v[i].first : v[0].first + 1.0 - v[5].first;
      if (g >= c.gap - 1e-3) starts.push_back(v[(i + 1) % 6].second);
    }
    for (int s0 : starts) {
      Var V;
      V.gap = c.gap;
      int colp[6];
      std::vector<std::pair<double, int>> cols;
      for (int cc = 0; cc < 6; ++cc) cols.push_back({frac(d[cc] - d[s0]), cc});
      std::sort(cols.begin(), cols.end());
      // the start column has key 0 -> first (ties with repeated points are harmless)
      for (int i = 0; i < 6; ++i) colp[i] = cols[i].second;
      int rowp[6];
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
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) V.phi.map[6 * i + j] = c.tr ? 6 * colp[j] + rowp[i] : 6 * rowp[i] + colp[j];
      V.phi.sign = c.sg ? -1 : 1;
      out.push_back(V);
    }
  }
  return out;
}

}  // namespace

// Pi(D0(P1) - D0(P2) - N) norm, interval-evaluated; N = nearest integers.
double mergeDistanceUpper(const double P1[36], const double P2[36]) {
  try {
    I W[36], Dl[36];
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) {
        I a = I(P1[6 * i + j]) - I(P1[6 * i]) - I(P1[j]) + I(P1[0]);
        I b = I(P2[6 * i + j]) - I(P2[6 * i]) - I(P2[j]) + I(P2[0]);
        double N = std::floor(midp(a) - midp(b) + 0.5);
        W[6 * i + j] = a - b - I(N);
      }
    piApply(W, Dl);
    I s(0.0);
    for (int p = 0; p < 36; ++p) s = s + sqr(Dl[p]);
    return sqrt0(s).hi;
  } catch (...) {
    return 1e300;
  }
}

namespace {
double approxMerge(const double P1[36], const double P2[36]) {
  double W[36], mr[6] = {0}, mc[6] = {0}, tm = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) {
      double a = P1[6 * i + j] - P1[6 * i] - P1[j] + P1[0];
      double b = P2[6 * i + j] - P2[6 * i] - P2[j] + P2[0];
      double N = std::floor(a - b + 0.5);
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
}  // namespace

bool Atlas::add(const Root& R) {
  double gmin;
  std::vector<Var> vs = variants(R.P, 0.02, gmin);
  if (vs.empty()) {  // degenerate: insert as is
    recs.push_back(R);
    return false;
  }
  int b0 = (int)std::floor(gmin * 100.0);
  if ((int)buckets.size() < 120) buckets.resize(120);
  struct Best {
    double d = 1e300;
    int idx = -1;
    int v = -1;
  } best;
  for (int b = std::max(0, b0 - 1); b <= std::min(b0 + 1, (int)buckets.size() - 1); ++b)
    for (int idx : buckets[b]) {
      for (size_t v = 0; v < vs.size(); ++v) {
        double Pv[36];
        for (int i = 0; i < 36; ++i) Pv[i] = vs[v].phi.sign * R.P[vs[v].phi.map[i]];
        double d = approxMerge(Pv, recs[idx].P);
        if (d < best.d) best = {d, idx, (int)v};
      }
    }
  if (best.idx >= 0 && best.d < cap) {
    double Pv[36];
    const Phi& f = vs[best.v].phi;
    for (int i = 0; i < 36; ++i) Pv[i] = f.sign * R.P[f.map[i]];
    double d = mergeDistanceUpper(Pv, recs[best.idx].P);
    I need = I(R.r) + I(d) + I(1e-8);
    if (need.hi + 2e-8 <= cap) {
      if (need.hi > recs[best.idx].r) recs[best.idx].r = need.hi;
      return true;
    }
  }
  Root T = applyPhi(R, vs[0].phi);
  {
    bool ident = vs[0].phi.sign == 1;
    for (int i = 0; i < 36 && ident; ++i) ident = vs[0].phi.map[i] == i;
    if (!ident) recertifyTransported(T);  // fresh constants in the stored orientation (transported ones kept on failure)
  }
  int b = std::min(std::max(0, b0), (int)buckets.size() - 1);
  buckets[b].push_back((int)recs.size());
  recs.push_back(T);
  return false;
}

bool writeRecords(const char* path, const std::vector<Root>& v) {
  FILE* f = fopen(path, "wb");
  if (!f) return false;
  uint64_t n = v.size();
  fwrite(&n, 8, 1, f);
  for (const Root& R : v) {
    int32_t k = R.k;
    fwrite(&k, 4, 1, f);
    fwrite(R.P, 8, 36, f);
    for (int a = 0; a < R.k; ++a) fwrite(R.T[a], 8, 36, f);
    double c[8] = {R.C, R.n0, R.id, R.ort, R.hor, R.mnorm, R.froot, R.r};
    fwrite(c, 8, 8, f);
  }
  fclose(f);
  return true;
}

bool readRecords(const char* path, std::vector<Root>& v) {
  FILE* f = fopen(path, "rb");
  if (!f) return false;
  uint64_t n;
  if (fread(&n, 8, 1, f) != 1) { fclose(f); return false; }
  v.clear();
  for (uint64_t i = 0; i < n; ++i) {
    Root R;
    int32_t k;
    if (fread(&k, 4, 1, f) != 1 || k < 0 || k > 10) { fclose(f); return false; }
    R.k = k;
    if (fread(R.P, 8, 36, f) != 36) { fclose(f); return false; }
    for (int a = 0; a < k; ++a)
      if (fread(R.T[a], 8, 36, f) != 36) { fclose(f); return false; }
    double c[8];
    if (fread(c, 8, 8, f) != 8) { fclose(f); return false; }
    R.C = c[0]; R.n0 = c[1]; R.id = c[2]; R.ort = c[3]; R.hor = c[4]; R.mnorm = c[5]; R.froot = c[6]; R.r = c[7];
    v.push_back(R);
  }
  uint8_t extra;
  bool eof = fread(&extra, 1, 1, f) == 0;
  fclose(f);
  return eof;
}

}  // namespace sa
