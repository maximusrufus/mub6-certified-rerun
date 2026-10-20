// covcheck.cpp -- independent coverage check of the seed cover (pre-registration item 5).
//
// Claim checked (SPEC.md sec 2, seed:cover): every bin pair (bX, bY) of the paper's seed
// enumeration that is not excluded by one of the paper's necessary conditions (F1, F2, sumok,
// the second gap inequality, the shared-interval test) has its FEASIBLE part (the points of the
// closed 9-dimensional bin that satisfy the linearised vanishing-sum constraints of
// seed:grid-error) within Euclidean distance R = 0.032 of a written centre with the same s.
//
// Independence: this program shares no code with the seed generator (impl/seed).  It has its own
// interval type (round-to-nearest + nextafter outward widening; no dependence on the FE rounding
// mode, needs -ffp-contract=off), its own unit-root table (roots_112.hpp, python-flint arb at 200
// bits, outward rounded), its own enumeration, filters and certificate.  The only shared inputs are
// the paper's formulas and the file seeds.txt that is being checked.
//
// Soundness of the certificate (seed:distance-bound, seed:dual-formula, any v):
//   N^2 ||q - u'||^2 = ||k + h||^2 <= ||k||^2 + 9/4 + 2 k.h,
//   2 k.h <= sum_j |k_j - (J^T v)_j| + 2 v.p + 2E(||v_X|| + ||v_Y||)  (also <= sum_j |k_j|, v = 0),
// for every feasible h (|h_j| <= 1/2, J h = p + r, ||r_X||,||r_Y|| <= E = 5 pi/(4N)).
// Success iff the interval UPPER bound of the right-hand side is < RN2_LO <= (R N)^2.
// A pair that is itself a centre is covered by its midpoint (half-diagonal 3/(2N) < R).
//
// Output: counts and, for every uncovered pair, its bins (the exact uncovered region).
// Exit 0 iff no uncovered pair; exit 1 otherwise; exit 2 on input/format problems.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>
#include "roots_112.hpp"

// ---------------------------------------------------------------- interval arithmetic
struct Iv {
  double lo, hi;
  Iv() : lo(0), hi(0) {}
  Iv(double a) : lo(a), hi(a) {}
  Iv(double a, double b) : lo(a), hi(b) {}
};
static inline double dn(double x) { return std::nextafter(x, -INFINITY); }
static inline double up(double x) { return std::nextafter(x, INFINITY); }
static inline Iv operator+(Iv a, Iv b) { return Iv(dn(a.lo + b.lo), up(a.hi + b.hi)); }
static inline Iv operator-(Iv a, Iv b) { return Iv(dn(a.lo - b.hi), up(a.hi - b.lo)); }
static inline Iv operator-(Iv a) { return Iv(-a.hi, -a.lo); }
static inline Iv operator*(Iv a, Iv b) {
  double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
  return Iv(dn(std::min(std::min(p1, p2), std::min(p3, p4))), up(std::max(std::max(p1, p2), std::max(p3, p4))));
}
static inline Iv operator/(Iv a, Iv b) {  // requires 0 not in b
  if (b.lo <= 0 && b.hi >= 0) { std::fprintf(stderr, "division by interval containing 0\n"); std::exit(2); }
  double p1 = a.lo / b.lo, p2 = a.lo / b.hi, p3 = a.hi / b.lo, p4 = a.hi / b.hi;
  return Iv(dn(std::min(std::min(p1, p2), std::min(p3, p4))), up(std::max(std::max(p1, p2), std::max(p3, p4))));
}
static inline Iv iabs(Iv a) {
  if (a.lo >= 0) return a;
  if (a.hi <= 0) return Iv(-a.hi, -a.lo);
  return Iv(0, std::max(-a.lo, a.hi));
}
static inline Iv isqrt(Iv a) {  // a.lo >= 0
  return Iv(dn(std::sqrt(std::max(0.0, a.lo))), up(std::sqrt(a.hi)));
}
static const Iv PI(PI_LO, PI_HI);

// ---------------------------------------------------------------- problem constants
static const int N = 112;
static const double HALF = 0.5;
// RN2_LO: a rigorous lower bound of (0.032*112)^2 = 12.845056 (exact decimal); 12.845 < 12.845056.
static const double RN2_LO = 12.845;
// E = 5 pi / (4 N) as an interval; w = 1/N.
static Iv E_IV() { return Iv(5.0) * PI / Iv(4.0 * N); }
static Iv COEF_IV() { return Iv(5.0) * PI / (Iv(4.0) * Iv((double)N) * Iv((double)N)); }  // 5 pi w^2 / 4
static Iv W2_IV() { return Iv(1.0) / Iv(2.0 * N); }                                       // w/2

// unit root of bin b as (cos, sin) enclosure
static inline Iv rc(int b) { return Iv(ROOT_COS[b][0], ROOT_COS[b][1]); }
static inline Iv rs(int b) { return Iv(ROOT_SIN[b][0], ROOT_SIN[b][1]); }

// Residual data for a list (0, u_0..u_4): c = (1 + sum e(u_j)) / (2 pi i) realified = (Im S, -Re S)/(2 pi);
// J_j = (cos, sin)(2 pi u_j).
struct ListData {
  int b[5];
  Iv c0, c1;        // c
  Iv J0[5], J1[5];  // columns
  double Jm0[5], Jm1[5], cm0, cm1;  // midpoints (heuristic directions only)
};
static void list_data(const int b[5], ListData &d) {
  Iv re(1.0), im(0.0);
  for (int j = 0; j < 5; ++j) {
    d.b[j] = b[j];
    d.J0[j] = rc(b[j]); d.J1[j] = rs(b[j]);
    re = re + d.J0[j]; im = im + d.J1[j];
    d.Jm0[j] = 0.5 * (d.J0[j].lo + d.J0[j].hi); d.Jm1[j] = 0.5 * (d.J1[j].lo + d.J1[j].hi);
  }
  Iv twopi = Iv(2.0) * PI;
  d.c0 = im / twopi; d.c1 = -(re / twopi);
  d.cm0 = 0.5 * (d.c0.lo + d.c0.hi); d.cm1 = 0.5 * (d.c1.lo + d.c1.hi);
}

// seed:sum-support in direction v (exact doubles): |v.c| <= coef ||v|| + (w/2) sum |v.J_j|.
// Reject iff the inequality is certainly violated.
static bool violates(const ListData &d, double v0, double v1) {
  if (!std::isfinite(v0) || !std::isfinite(v1)) return false;
  Iv V0(v0), V1(v1);
  Iv lhs = iabs(V0 * d.c0 + V1 * d.c1);
  Iv nv = isqrt(V0 * V0 + V1 * V1);
  Iv rhs = COEF_IV() * nv;
  Iv w2 = W2_IV();
  for (int j = 0; j < 5; ++j) rhs = rhs + w2 * iabs(V0 * d.J0[j] + V1 * d.J1[j]);
  return lhs.lo > rhs.hi;
}
static bool sumok(const ListData &d) {
  // directions: residual c (normalised midpoint), perpendiculars to each J_j midpoint
  double nc = std::hypot(d.cm0, d.cm1);
  if (nc > 0) { if (violates(d, d.cm0 / nc, d.cm1 / nc)) return false; }
  for (int j = 0; j < 5; ++j) if (violates(d, -d.Jm1[j], d.Jm0[j])) return false;
  return true;
}

// seed:clip-interval: enclosure of the admissible shared displacement h_t for list d at position t.
// Directions: J_t midpoint and perpendiculars to the other four columns; skip lower|a| < 1e-6.
static void extend(const ListData &d, int t, double &lo, double &hi) {
  lo = -HALF; hi = HALF;
  Iv E = E_IV(), Nn((double)N);
  for (int dir = 0; dir < 5; ++dir) {
    double v0, v1;
    if (dir == 0) { v0 = d.Jm0[t]; v1 = d.Jm1[t]; }
    else {
      int j = dir - 1; if (j >= t) ++j;  // the other four columns
      v0 = -d.Jm1[j]; v1 = d.Jm0[j];
    }
    Iv V0(v0), V1(v1);
    Iv a = V0 * d.J0[t] + V1 * d.J1[t];
    Iv aa = iabs(a);
    if (aa.lo < 1e-6) continue;
    Iv r = E * isqrt(V0 * V0 + V1 * V1);
    for (int j = 0; j < 5; ++j) if (j != t) r = r + Iv(HALF) * iabs(V0 * d.J0[j] + V1 * d.J1[j]);
    Iv vc = Nn * (V0 * d.c0 + V1 * d.c1);
    Iv centre = -(vc / a);
    Iv rad = r / aa;
    Iv L = centre - rad, U = centre + rad;
    lo = std::max(lo, L.lo); hi = std::min(hi, U.hi);  // outward: keep lower of L, upper of U
  }
}

// ---------------------------------------------------------------- enumeration
static int select_s(const int bx[5]) {
  int best = -1, bs = 0;
  for (int s = 0; s < 5; ++s) { int v = std::abs(2 * bx[s] + 1 - N); if (best < 0 || v < best) { best = v; bs = s; } }
  return bs;
}
// 2N * midpoint gap of the sorted bins q[0..4] together with the fixed 0 (exact integer)
static int gap2N(int q[5]) {
  std::sort(q, q + 5);
  int g = 2 * q[0] + 1;
  for (int i = 1; i < 5; ++i) g = std::max(g, 2 * (q[i] - q[i - 1]));
  g = std::max(g, 2 * N - 2 * q[4] - 1);
  return g;
}

struct XL { int bx[5]; int s; };

// ---------------------------------------------------------------- centres
struct Key9 {
  uint64_t a;
  bool operator==(const Key9 &o) const { return a == o.a; }
};
struct Key9H { size_t operator()(const Key9 &k) const { return std::hash<uint64_t>()(k.a * 0x9E3779B97F4A7C15ULL); } };
static inline uint64_t pack5(const int b[5]) { uint64_t a = 0; for (int i = 0; i < 5; ++i) a = a * 128 + b[i]; return a; }
static inline uint32_t pack4(const int b[4]) { uint32_t a = 0; for (int i = 0; i < 4; ++i) a = a * 128 + b[i]; return a; }

struct CentreX {              // one X-tuple of the centre file, with its tails
  int bx[5]; int s; uint64_t key;
  std::vector<uint32_t> tails;  // packed (y2..y5) sorted
};

// ---------------------------------------------------------------- certificate
struct PairCtx {
  // h order: 0..4 = x1..x5 displacements, 5..8 = y2..y5; the shared y entry is h_s.
  int s;
  Iv JT[9][4];       // J_*^T: row j = coefficients of h_j in (X row0, X row1, Y row0, Y row1)
  double JTm[9][4];
  Iv p[4]; double pm[4];
};
static void build_ctx(const ListData &dx, const ListData &dy, int s, PairCtx &c) {
  c.s = s;
  for (int j = 0; j < 9; ++j) for (int a = 0; a < 4; ++a) c.JT[j][a] = Iv(0.0);
  for (int j = 0; j < 5; ++j) { c.JT[j][0] = dx.J0[j]; c.JT[j][1] = dx.J1[j]; }
  // Y list (0, u_shared, y2..y5): varying index 0 = shared -> h_s; indices 1..4 -> h_5..h_8
  c.JT[s][2] = c.JT[s][2] + dy.J0[0]; c.JT[s][3] = c.JT[s][3] + dy.J1[0];
  for (int m = 1; m < 5; ++m) { c.JT[4 + m][2] = dy.J0[m]; c.JT[4 + m][3] = dy.J1[m]; }
  Iv Nn((double)N);
  c.p[0] = -(Nn * dx.c0); c.p[1] = -(Nn * dx.c1); c.p[2] = -(Nn * dy.c0); c.p[3] = -(Nn * dy.c1);
  for (int j = 0; j < 9; ++j) for (int a = 0; a < 4; ++a) c.JTm[j][a] = 0.5 * (c.JT[j][a].lo + c.JT[j][a].hi);
  for (int a = 0; a < 4; ++a) c.pm[a] = 0.5 * (c.p[a].lo + c.p[a].hi);
}

// rigorous upper bound of  sum_j |k_j - (J^T v)_j| + 2 v.p + 2E(||v_X|| + ||v_Y||)
static double dual_hi(const PairCtx &c, const int k[9], const double v[4]) {
  Iv V[4]; for (int a = 0; a < 4; ++a) V[a] = Iv(v[a]);
  Iv acc(0.0);
  for (int j = 0; j < 9; ++j) {
    Iv t((double)k[j]);
    for (int a = 0; a < 4; ++a) t = t - c.JT[j][a] * V[a];
    acc = acc + iabs(t);
  }
  for (int a = 0; a < 4; ++a) acc = acc + Iv(2.0) * V[a] * c.p[a];
  Iv nx = isqrt(V[0] * V[0] + V[1] * V[1]), ny = isqrt(V[2] * V[2] + V[3] * V[3]);
  acc = acc + Iv(2.0) * E_IV() * (nx + ny);
  return acc.hi;
}

// heuristic proposal of v: least squares then IRLS on the smoothed L1 objective (doubles only)
static bool solve4(double A[4][4], double b[4], double x[4]) {
  // Gaussian elimination with partial pivoting
  double M[4][5];
  for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) M[i][j] = A[i][j]; M[i][4] = b[i]; }
  for (int col = 0; col < 4; ++col) {
    int piv = col;
    for (int r = col + 1; r < 4; ++r) if (std::fabs(M[r][col]) > std::fabs(M[piv][col])) piv = r;
    if (std::fabs(M[piv][col]) < 1e-14) return false;
    if (piv != col) for (int j = 0; j < 5; ++j) std::swap(M[piv][j], M[col][j]);
    for (int r = 0; r < 4; ++r) if (r != col) {
      double f = M[r][col] / M[col][col];
      for (int j = col; j < 5; ++j) M[r][j] -= f * M[col][j];
    }
  }
  for (int i = 0; i < 4; ++i) { x[i] = M[i][4] / M[i][i]; if (!std::isfinite(x[i])) return false; }
  return true;
}
static double obj_mid(const PairCtx &c, const int k[9], const double v[4]) {
  double acc = 0;
  for (int j = 0; j < 9; ++j) { double t = k[j]; for (int a = 0; a < 4; ++a) t -= c.JTm[j][a] * v[a]; acc += std::fabs(t); }
  for (int a = 0; a < 4; ++a) acc += 2 * v[a] * c.pm[a];
  double E = 0.5 * (E_IV().lo + E_IV().hi);
  acc += 2 * E * (std::hypot(v[0], v[1]) + std::hypot(v[2], v[3]));
  return acc;
}
static bool propose_v(const PairCtx &c, const int k[9], double v[4]) {
  // Damped Newton on the smoothed objective  sum_j sqrt(r_j^2 + eps^2) + 2 v.p + 2E(sqrt(|vX|^2+eps^2) + sqrt(|vY|^2+eps^2)),
  // r = k - J^T v, with backtracking on the true (double) objective; eps decreases 1e-1 -> 1e-4.
  // Start: least squares  J^T v ~ k.  Heuristic only: any v is sound.
  double E = 0.5 * (E_IV().lo + E_IV().hi);
  double cur[4] = {0, 0, 0, 0};
  {
    double A[4][4], b[4];
    for (int a = 0; a < 4; ++a) {
      for (int bb = 0; bb < 4; ++bb) { double s = 0; for (int j = 0; j < 9; ++j) s += c.JTm[j][a] * c.JTm[j][bb]; A[a][bb] = s + (a == bb ? 1e-9 : 0); }
      double s = 0; for (int j = 0; j < 9; ++j) s += c.JTm[j][a] * k[j]; b[a] = s;
    }
    if (!solve4(A, b, cur)) for (int a = 0; a < 4; ++a) cur[a] = 0;
  }
  double best[4] = {0, 0, 0, 0}; double bestf = obj_mid(c, k, best);
  double fcur = obj_mid(c, k, cur);
  if (fcur < bestf) { bestf = fcur; for (int a = 0; a < 4; ++a) best[a] = cur[a]; }
  for (double eps = 1.0; eps >= 1e-6 * 0.5; eps *= 0.1) {
    for (int it = 0; it < 40; ++it) {
      double g[4] = {0, 0, 0, 0}, H[4][4] = {{0}};
      for (int j = 0; j < 9; ++j) {
        double r = k[j]; for (int a = 0; a < 4; ++a) r -= c.JTm[j][a] * cur[a];
        double q = std::sqrt(r * r + eps * eps);
        double d1 = r / q, d2 = eps * eps / (q * q * q);
        for (int a = 0; a < 4; ++a) {
          g[a] += -c.JTm[j][a] * d1;
          for (int bb = 0; bb < 4; ++bb) H[a][bb] += c.JTm[j][a] * c.JTm[j][bb] * d2;
        }
      }
      for (int a = 0; a < 4; ++a) g[a] += 2 * c.pm[a];
      for (int blk = 0; blk < 2; ++blk) {
        double x = cur[2 * blk], y = cur[2 * blk + 1];
        double q = std::sqrt(x * x + y * y + eps * eps);
        g[2 * blk] += 2 * E * x / q; g[2 * blk + 1] += 2 * E * y / q;
        double q3 = q * q * q;
        H[2 * blk][2 * blk] += 2 * E * (y * y + eps * eps) / q3;
        H[2 * blk + 1][2 * blk + 1] += 2 * E * (x * x + eps * eps) / q3;
        H[2 * blk][2 * blk + 1] += -2 * E * x * y / q3; H[2 * blk + 1][2 * blk] += -2 * E * x * y / q3;
      }
      for (int a = 0; a < 4; ++a) H[a][a] += 1e-12;
      double step[4], rhs[4]; for (int a = 0; a < 4; ++a) rhs[a] = -g[a];
      if (!solve4(H, rhs, step)) break;
      double t = 1.0; bool moved = false;
      for (int ls = 0; ls < 30; ++ls, t *= 0.5) {
        double nv[4]; for (int a = 0; a < 4; ++a) nv[a] = cur[a] + t * step[a];
        double f = obj_mid(c, k, nv);
        if (f < fcur) { fcur = f; for (int a = 0; a < 4; ++a) cur[a] = nv[a]; moved = true; break; }
      }
      if (fcur < bestf) { bestf = fcur; for (int a = 0; a < 4; ++a) best[a] = cur[a]; }
      if (!moved) break;
    }
  }
  for (int a = 0; a < 4; ++a) { v[a] = best[a]; if (!std::isfinite(v[a]) || std::fabs(v[a]) > 1e15) return false; }
  return true;
}

// Certified: the feasible part of the pair with k = N(u - u') (pair midpoint minus centre) is within R of u'.
// Sign matters for the dual bound: ||k + h||^2 = N^2 ||q - u'||^2 with h = N(q - u).
static bool certified(const PairCtx &c, const int k[9], long &nbox, long &ndual) {
  int kn2 = 0, k1 = 0;
  for (int j = 0; j < 9; ++j) { kn2 += k[j] * k[j]; k1 += std::abs(k[j]); }
  if (kn2 + 2.25 + k1 < RN2_LO) { ++nbox; return true; }
  double v[4];
  if (!propose_v(c, k, v)) return false;
  double B = dual_hi(c, k, v);
  Iv tot = Iv((double)kn2) + Iv(2.25) + Iv(B);
  if (tot.hi < RN2_LO) { ++ndual; return true; }
  return false;
}

// ---------------------------------------------------------------- main
int main(int argc, char **argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: covcheck seeds.txt [--limit-x M]\n"); return 2; }
  long limitx = 0; bool debug = false, dump = false;
  for (int i = 2; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--limit-x") && i + 1 < argc) limitx = std::atol(argv[++i]);
    else if (!std::strcmp(argv[i], "--debug")) debug = true;
    else if (!std::strcmp(argv[i], "--dump-pairs")) dump = true;
  }
  if (ROOTS_N != N) { std::fprintf(stderr, "root table N mismatch\n"); return 2; }

  // ---- read centres: "s b0 b1 b2 b3 b4 b5 b6 b7 b8", nine order (xi, other four x, y2..y5)
  FILE *f = std::fopen(argv[1], "r");
  if (!f) { std::fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
  std::unordered_map<Key9, int, Key9H> xidx;
  std::vector<CentreX> cx;
  long ncent = 0;
  {
    int s, b[9];
    while (std::fscanf(f, "%d %d %d %d %d %d %d %d %d %d", &s, &b[0], &b[1], &b[2], &b[3], &b[4], &b[5], &b[6], &b[7], &b[8]) == 10) {
      ++ncent;
      if (s < 0 || s > 4) { std::fprintf(stderr, "bad s on line %ld\n", ncent); return 2; }
      for (int j = 0; j < 9; ++j) if (b[j] < 0 || b[j] >= N) { std::fprintf(stderr, "bad bin on line %ld\n", ncent); return 2; }
      int bx[5]; int p = 1;
      for (int m = 0; m < 5; ++m) bx[m] = (m == s) ? b[0] : b[p++];
      for (int m = 1; m < 5; ++m) if (bx[m] < bx[m - 1]) { std::fprintf(stderr, "x-tuple not sorted, line %ld\n", ncent); return 2; }
      if (select_s(bx) != s) { std::fprintf(stderr, "s inconsistent with selection rule, line %ld\n", ncent); return 2; }
      for (int m = 6; m < 9; ++m) if (b[m] < b[m - 1]) { std::fprintf(stderr, "y-tail not sorted, line %ld\n", ncent); return 2; }
      Key9 key{pack5(bx)};
      auto it = xidx.find(key);
      int id;
      if (it == xidx.end()) { id = (int)cx.size(); xidx.emplace(key, id); CentreX c; std::memcpy(c.bx, bx, sizeof bx); c.s = s; c.key = key.a; cx.push_back(c); }
      else id = it->second;
      cx[id].tails.push_back(pack4(b + 5));
    }
    std::fclose(f);
  }
  long dup = 0;
  for (auto &c : cx) { std::sort(c.tails.begin(), c.tails.end()); size_t n = c.tails.size(); c.tails.erase(std::unique(c.tails.begin(), c.tails.end()), c.tails.end()); dup += n - c.tails.size(); }
  std::fprintf(stderr, "centres %ld (x-tuples %zu, duplicates %ld)\n", ncent, cx.size(), dup);

  // ---- X-lists: nondecreasing 5-tuples passing F1, F2, sumok
  std::vector<XL> X;
  for (int a = 0; a < N; ++a) for (int b = a; b < N; ++b) for (int c = b; c < N; ++c) for (int d = c; d < N; ++d) for (int l = d; l < N; ++l) {
    if (!(a <= l - d + 1)) continue;                     // F1
    int q[5] = {a, b, c, d, l};
    if (!(gap2N(q) <= 2 * N - 2 * (l - 1))) continue;    // F2  (g_mid <= 1 - (l-1)/N)
    int bx[5] = {a, b, c, d, l};
    ListData dx; list_data(bx, dx);
    if (!sumok(dx)) continue;                            // F3
    XL x; std::memcpy(x.bx, bx, sizeof bx); x.s = select_s(bx); X.push_back(x);
  }
  std::fprintf(stderr, "x-lists %zu\n", X.size());
  if (limitx > 0 && (long)X.size() > limitx) X.resize(limitx);

  // ---- Y tails per shared bin: nondecreasing 4-tuples passing sumok on (0, u_shared, u_1..u_4)
  std::vector<std::vector<std::array<int, 4>>> tails(N);
  #pragma omp parallel for schedule(dynamic)
  for (int sb = 0; sb < N; ++sb) {
    std::vector<std::array<int, 4>> &T = tails[sb];
    for (int a = 0; a < N; ++a) for (int b = a; b < N; ++b) for (int c = b; c < N; ++c) for (int d = c; d < N; ++d) {
      int by[5] = {sb, a, b, c, d};
      ListData dy; list_data(by, dy);
      if (!sumok(dy)) continue;
      T.push_back({a, b, c, d});
    }
  }
  long ntails = 0; for (auto &T : tails) ntails += T.size();
  std::fprintf(stderr, "y-tails total %ld\n", ntails);

  // ---- neighbour X-tuples among centres (same s, ||kx||^2 <= CUT)
  const int CUT = 16;
  std::vector<std::vector<std::pair<int, int>>> nbr(X.size());  // (centre x id, ||kx||^2)
  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < X.size(); ++i) {
    for (size_t j = 0; j < cx.size(); ++j) {
      if (cx[j].s != X[i].s) continue;
      int d2 = 0; for (int m = 0; m < 5 && d2 <= CUT; ++m) { int d = cx[j].bx[m] - X[i].bx[m]; d2 += d * d; }
      if (d2 <= CUT) nbr[i].push_back({(int)j, d2});
    }
    std::sort(nbr[i].begin(), nbr[i].end(), [](auto &a, auto &b) { return a.second < b.second; });
  }

  // ---- pairs
  long n_sum = 0, n_gap = 0, n_ext = 0, n_self = 0, n_box = 0, n_dual = 0, n_unc = 0;
  std::vector<std::string> uncovered;
  #pragma omp parallel for schedule(dynamic) reduction(+:n_sum,n_gap,n_ext,n_self,n_box,n_dual,n_unc)
  for (size_t i = 0; i < X.size(); ++i) {
    const XL &x = X[i];
    int sb = x.bx[x.s], l = x.bx[4];
    ListData dx; list_data(x.bx, dx);
    double xlo, xhi; extend(dx, x.s, xlo, xhi);
    Key9 selfkey{pack5(x.bx)};
    auto selfit = xidx.find(selfkey);
    const CentreX *self = (selfit == xidx.end()) ? nullptr : &cx[selfit->second];
    std::vector<std::string> local_unc;
    for (const auto &t : tails[sb]) {
      ++n_sum;
      int q[5] = {sb, t[0], t[1], t[2], t[3]};
      if (!(gap2N(q) >= 2 * N - 2 * (l + 2))) continue;   // second gap inequality  g_y,mid >= 1 - (l+2)/N
      ++n_gap;
      int by[5] = {sb, t[0], t[1], t[2], t[3]};
      ListData dy; list_data(by, dy);
      double ylo, yhi; extend(dy, 0, ylo, yhi);
      if (xhi < ylo || yhi < xlo) continue;                 // strictly disjoint -> excluded
      ++n_ext;
      if (dump) {
        #pragma omp critical
        std::printf("PAIR s=%d x=%d,%d,%d,%d,%d y=%d,%d,%d,%d\n", x.s, x.bx[0], x.bx[1], x.bx[2], x.bx[3], x.bx[4], t[0], t[1], t[2], t[3]);
        continue;
      }
      // coverage
      uint32_t tk = pack4(t.data());
      if (self && std::binary_search(self->tails.begin(), self->tails.end(), tk)) { ++n_self; continue; }
      PairCtx ctx; build_ctx(dx, dy, x.s, ctx);
      bool ok = false;
      for (const auto &nb : nbr[i]) {
        const CentreX &c = cx[nb.first];
        int rem = CUT - nb.second;
        int k[9];
        for (int m = 0; m < 5; ++m) k[m] = x.bx[m] - c.bx[m];
        for (uint32_t ptk : c.tails) {
          int cy[4] = {(int)(ptk >> 21) & 127, (int)(ptk >> 14) & 127, (int)(ptk >> 7) & 127, (int)ptk & 127};
          int d2 = 0; for (int m = 0; m < 4; ++m) { k[5 + m] = t[m] - cy[m]; d2 += k[5 + m] * k[5 + m]; }
          if (d2 > rem) continue;
          if (certified(ctx, k, n_box, n_dual)) { ok = true; break; }
        }
        if (ok) break;
      }
      if (!ok) {
        ++n_unc;
        if (debug && n_unc <= 20) {
          #pragma omp critical
          {
            std::printf("DEBUG pair s=%d x=%d,%d,%d,%d,%d y=%d,%d,%d,%d  xint=[%g,%g] yint=[%g,%g]\n", x.s, x.bx[0], x.bx[1], x.bx[2], x.bx[3], x.bx[4], t[0], t[1], t[2], t[3], xlo, xhi, ylo, yhi);
            for (const auto &nb : nbr[i]) {
              const CentreX &c = cx[nb.first];
              int k[9]; for (int m = 0; m < 5; ++m) k[m] = x.bx[m] - c.bx[m];
              for (uint32_t ptk : c.tails) {
                int cy[4] = {(int)(ptk >> 21) & 127, (int)(ptk >> 14) & 127, (int)(ptk >> 7) & 127, (int)ptk & 127};
                int d2 = nb.second; for (int m = 0; m < 4; ++m) { k[5 + m] = t[m] - cy[m]; d2 += k[5 + m] * k[5 + m]; }
                if (d2 > CUT) continue;
                double v[4]; bool pv = propose_v(ctx, k, v);
                double B = pv ? dual_hi(ctx, k, v) : 1e300;
                int k1 = 0; for (int m = 0; m < 9; ++m) k1 += std::abs(k[m]);
                std::printf("  cand k=%d,%d,%d,%d,%d,%d,%d,%d,%d kn2=%d box=%g dual=%g v=%g,%g,%g,%g\n", k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],k[8], d2, d2+2.25+k1, d2+2.25+B, v[0],v[1],v[2],v[3]);
              }
            }
          }
        }
        char buf[128];
        std::snprintf(buf, sizeof buf, "s=%d x=%d,%d,%d,%d,%d y=%d,%d,%d,%d", x.s, x.bx[0], x.bx[1], x.bx[2], x.bx[3], x.bx[4], t[0], t[1], t[2], t[3]);
        local_unc.push_back(buf);
      }
    }
    if (!local_unc.empty()) {
      #pragma omp critical
      uncovered.insert(uncovered.end(), local_unc.begin(), local_unc.end());
    }
  }
  std::printf("centres %ld\nx-lists %zu\npairs sumok-only %ld\npairs +gap2 %ld\npairs +extend %ld\n", ncent, X.size(), n_sum, n_gap, n_ext);
  std::printf("covered by own centre %ld\ncovered by box bound %ld\ncovered by dual bound %ld\nuncovered %ld\n", n_self, n_box, n_dual, n_unc);
  std::sort(uncovered.begin(), uncovered.end());
  for (auto &u : uncovered) std::printf("UNCOVERED %s\n", u.c_str());
  std::printf(n_unc == 0 ? "COVERAGE-PROVEN\n" : "COVERAGE-FAILED\n");
  return n_unc == 0 ? 0 : 1;
}
