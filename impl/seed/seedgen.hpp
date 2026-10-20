// seedgen.hpp -- enumeration of X-lists, Y-tails, pairs and the reuse traversal.
// Every OpenMP thread that evaluates a certificate holds its own iv::RoundUp guard.
#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include "cover.hpp"
#include "reuse.hpp"

namespace seedcov {

struct Options {
  int N = 112;
  int cut = -1;          // reach: only try previous centers with ||k||^2 <= cut (-1: floor(RN2)+2 = 14 at N=112)
  long maxx = 0;
  bool count_only = false;
  bool verbose = true;
};

struct XList {
  uint8_t bx[5];   // bins of x1..x5 in original (nondecreasing) order
  uint8_t s;       // selected position
  uint8_t sb;      // shared bin = bx[s]
  uint8_t nx[5];   // nine-order X part: (xi, other four in former order)
};

inline void make_nine_x(XList &x) {
  x.nx[0] = x.bx[x.s];
  int p = 1;
  for (int m = 0; m < 5; ++m) if (m != x.s) x.nx[p++] = x.bx[m];
}

// All weakly increasing 5-bin tuples passing F1 (orientation), F2 (midpoint
// gap) and F3 (sumok); lexicographic order.
inline std::vector<XList> enumerate_x(const Tables &T) {
  int N = T.N;
  std::vector<std::vector<XList>> per(N);
#pragma omp parallel
  {
    iv::RoundUp guard;
#pragma omp for schedule(dynamic, 1)
    for (int a = 0; a < N; ++a) {
      std::vector<XList> loc;
      for (int b = a; b < N; ++b) for (int c = b; c < N; ++c) for (int d = c; d < N; ++d) for (int l = d; l < N; ++l) {
        int bx[5] = {a, b, c, d, l};
        if (!x_orient_ok(bx)) continue;
        if (!x_gap_ok(T, bx)) continue;
        if (!sumok_bins(T, bx)) continue;
        XList x;
        for (int m = 0; m < 5; ++m) x.bx[m] = (uint8_t)bx[m];
        x.s = (uint8_t)select_s(N, bx);
        x.sb = x.bx[x.s];
        make_nine_x(x);
        loc.push_back(x);
      }
      per[a] = std::move(loc);
    }
  }
  std::vector<XList> X;
  for (auto &v : per) X.insert(X.end(), v.begin(), v.end());
  return X;
}

// Y tails per shared bin: weakly increasing (t1..t4) with sumok on (shared,t1..t4).
// Packed little-endian into uint32 (t1 lowest byte); lexicographic order.  For each
// tail the outward interval of the shared grid displacement (extend with t=0) is
// precomputed: it depends only on (shared, tail), not on the X-list.
struct TailTable {
  std::vector<std::vector<uint32_t>> tail;
  std::vector<std::vector<double>> ylo, yhi;
};

inline void unpack_tail(uint32_t tl, int *t) {
  t[0] = (int)(tl & 255); t[1] = (int)((tl >> 8) & 255); t[2] = (int)((tl >> 16) & 255); t[3] = (int)(tl >> 24);
}

inline TailTable build_tails(const Tables &T, const std::vector<XList> &X) {
  int N = T.N;
  std::vector<char> need(N, 0);
  for (const XList &x : X) need[x.sb] = 1;
  TailTable tt;
  tt.tail.resize(N); tt.ylo.resize(N); tt.yhi.resize(N);
#pragma omp parallel
  {
    iv::RoundUp guard;
#pragma omp for schedule(dynamic, 1)
    for (int sb = 0; sb < N; ++sb) {
      if (!need[sb]) continue;
      auto &v = tt.tail[sb]; auto &lo = tt.ylo[sb]; auto &hi = tt.yhi[sb];
      for (int t1 = 0; t1 < N; ++t1) for (int t2 = t1; t2 < N; ++t2) for (int t3 = t2; t3 < N; ++t3) for (int t4 = t3; t4 < N; ++t4) {
        int b[5] = {sb, t1, t2, t3, t4};
        if (!sumok_bins(T, b)) continue;
        ListData dy; make_ld(T, b, dy);
        Ivl iy = extend(T, dy, 0);
        v.push_back((uint32_t)t1 | ((uint32_t)t2 << 8) | ((uint32_t)t3 << 16) | ((uint32_t)t4 << 24));
        lo.push_back(iy.lo); hi.push_back(iy.hi);
      }
    }
  }
  return tt;
}

struct PairCounts { long sumok_only = 0, gap2 = 0, extend = 0; };

inline PairCounts count_pairs(const Tables &T, const std::vector<XList> &X, const TailTable &tt) {
  long p0 = 0, p1 = 0, p2 = 0;
#pragma omp parallel reduction(+ : p0, p1, p2)
  {
    iv::RoundUp guard;
#pragma omp for schedule(dynamic, 16)
    for (size_t xi = 0; xi < X.size(); ++xi) {
      const XList &x = X[xi];
      int bx[5]; for (int m = 0; m < 5; ++m) bx[m] = x.bx[m];
      ListData dx; make_ld(T, bx, dx);
      Ivl ix = extend(T, dx, x.s);
      const auto &tv = tt.tail[x.sb]; const auto &lo = tt.ylo[x.sb]; const auto &hi = tt.yhi[x.sb];
      for (size_t i = 0; i < tv.size(); ++i) {
        int tail[4]; unpack_tail(tv[i], tail);
        ++p0;
        if (!y_gap_ok(T, x.bx[4], x.sb, tail)) continue;
        ++p1;
        if (ivl_disjoint(ix, Ivl{lo[i], hi[i]})) continue;
        ++p2;
      }
    }
  }
  PairCounts pc; pc.sumok_only = p0; pc.gap2 = p1; pc.extend = p2;
  return pc;
}

// ---- reuse traversal --------------------------------------------------------
struct Center { uint32_t xi; uint32_t tail; };

struct CoverResult {
  std::vector<Center> centers;   // traversal (creation) order
  long pairs = 0, n_box = 0, n_dual = 0, n_cache = 0;
};

inline void center_bins(const XList &x, uint32_t tail, int *b) {
  for (int m = 0; m < 5; ++m) b[m] = x.nx[m];
  int t[4]; unpack_tail(tail, t);
  for (int m = 0; m < 4; ++m) b[5 + m] = t[m];
}

// Traversal order (SPEC A2): X-lists lexicographic, Y-tails lexicographic.  The
// selected position s partitions the problem (reuse needs equal s), so the five
// streams are independent and run in parallel; the result is deterministic.
inline CoverResult run_cover(const Tables &T, const std::vector<XList> &X, const TailTable &tt, const Options &o) {
  const int cut = o.cut >= 0 ? o.cut : (int)T.RN2lo + 2;
  std::vector<int> grp[5];
  for (size_t i = 0; i < X.size(); ++i) grp[X[i].s].push_back((int)i);
  std::vector<Center> res[5];
  long pairs[5] = {0}, nbox[5] = {0}, ndual[5] = {0}, ncache[5] = {0};

#pragma omp parallel for schedule(dynamic, 1)
  for (int s = 0; s < 5; ++s) {
    iv::RoundUp guard;
    const std::vector<int> &xs = grp[s];
    int n = (int)xs.size();
    std::vector<std::vector<uint32_t>> cen(n);
    int cl = -1; uint32_t ct = 0;
    std::vector<std::pair<int, int>> nb;
    for (int li = 0; li < n; ++li) {
      const XList &x = X[xs[li]];
      if (getenv("SEED_PROGRESS") && li % 20 == 0) fprintf(stderr, "s=%d li=%d/%d centers=%zu pairs=%ld\n", s, li, n, res[s].size(), pairs[s]);
      nb.clear();
      for (int lj = 0; lj <= li; ++lj) {
        const XList &y = X[xs[lj]];
        int kx2 = 0;
        for (int m = 0; m < 5; ++m) { int d = (int)x.nx[m] - (int)y.nx[m]; kx2 += d * d; }
        if (kx2 <= cut) nb.push_back({kx2, -lj});
      }
      std::sort(nb.begin(), nb.end());
      int bxi[5]; for (int m = 0; m < 5; ++m) bxi[m] = x.bx[m];
      ListData dx; make_ld(T, bxi, dx);
      Ivl ix = extend(T, dx, s);
      const auto &tv = tt.tail[x.sb]; const auto &lov = tt.ylo[x.sb]; const auto &hiv = tt.yhi[x.sb];
      for (size_t ti = 0; ti < tv.size(); ++ti) {
        uint32_t tl = tv[ti];
        int tail[4]; unpack_tail(tl, tail);
        if (!y_gap_ok(T, x.bx[4], x.sb, tail)) continue;
        if (ivl_disjoint(ix, Ivl{lov[ti], hiv[ti]})) continue;
        ++pairs[s];
        PairCtx ctx; bool ctxBuilt = false;
        int k[9];
        auto attempt = [&](int lj, uint32_t tl2) -> bool {
          const XList &y = X[xs[lj]];
          int kn2 = 0, kl1 = 0;
          for (int m = 0; m < 5; ++m) { k[m] = (int)x.nx[m] - (int)y.nx[m]; kn2 += k[m] * k[m]; kl1 += std::abs(k[m]); }
          int t2[4]; unpack_tail(tl2, t2);
          for (int m = 0; m < 4; ++m) { k[5 + m] = tail[m] - t2[m]; kn2 += k[5 + m] * k[5 + m]; kl1 += std::abs(k[5 + m]); }
          if (kn2 > cut) return false;
          if (box_success(T, kn2, kl1)) { ++nbox[s]; return true; }
          if (!ctxBuilt) {
            int yb[5] = {x.sb, tail[0], tail[1], tail[2], tail[3]};
            ListData dy; make_ld(T, yb, dy);
            build_ctx(T, dx, dy, s, ctx); ctxBuilt = true;
          }
          if (dual_success(T, ctx, k, kn2)) { ++ndual[s]; return true; }
          return false;
        };
        bool found = false;
        if (cl >= 0 && attempt(cl, ct)) { found = true; ++ncache[s]; }
        if (!found) {
          for (auto &pr : nb) {
            int lj = -pr.second;
            // centres of an X-list are stored in traversal (= lexicographic tail)
            // order, so the t1 band |t1-t1'| <= floor(sqrt(cut-kx2)) is contiguous
            const std::vector<uint32_t> &cv = cen[lj];
            int rr = 0; while ((rr + 1) * (rr + 1) <= cut - pr.first) ++rr;
            int lo1 = tail[0] - rr, hi1 = tail[0] + rr;
            size_t a0 = std::lower_bound(cv.begin(), cv.end(), lo1, [](uint32_t v, int key) { return (int)(v & 255) < key; }) - cv.begin();
            for (size_t ci = a0; ci < cv.size(); ++ci) {
              uint32_t tl2 = cv[ci];
              if ((int)(tl2 & 255) > hi1) break;
              if (lj == cl && tl2 == ct) continue;
              if (attempt(lj, tl2)) { found = true; cl = lj; ct = tl2; break; }
            }
            if (found) break;
          }
        }
        if (!found) {
          cen[li].push_back(tl);
          res[s].push_back(Center{(uint32_t)xs[li], tl});
          cl = li; ct = tl;
        }
      }
    }
  }
  CoverResult R;
  for (int s = 0; s < 5; ++s) {
    R.centers.insert(R.centers.end(), res[s].begin(), res[s].end());
    R.pairs += pairs[s]; R.n_box += nbox[s]; R.n_dual += ndual[s]; R.n_cache += ncache[s];
  }
  std::stable_sort(R.centers.begin(), R.centers.end(), [](const Center &a, const Center &b) { return a.xi < b.xi; });
  return R;
}

}  // namespace seedcov
