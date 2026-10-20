// Graph search of SPEC 5.5 (exclusion.tex exclude:graph-search): two disjoint 6-cliques A,B in the O-graph with
// A x B inside the U-graph.  Exact bitset implementation; a step-cap exhaustion is reported through `gaveUp`
// and callers must treat it as POSSIBLE (never as a negative certificate).
#pragma once
#include <cstdint>
#include <vector>

namespace mb {

struct Graph {
  int n = 0, W = 0;
  std::vector<uint64_t> O, U;  // n*W words each
  void init(int n_) {
    n = n_;
    W = (n + 63) / 64;
    O.assign((size_t)n * W, 0);
    U.assign((size_t)n * W, 0);
  }
  const uint64_t* o(int i) const { return &O[(size_t)i * W]; }
  const uint64_t* u(int i) const { return &U[(size_t)i * W]; }
  void addO(int i, int j) { if (i == j) return; O[(size_t)i * W + (j >> 6)] |= 1ull << (j & 63); O[(size_t)j * W + (i >> 6)] |= 1ull << (i & 63); }
  void addU(int i, int j) { if (i == j) return; U[(size_t)i * W + (j >> 6)] |= 1ull << (j & 63); U[(size_t)j * W + (i >> 6)] |= 1ull << (i & 63); }
};

struct GraphSearch {
  const Graph& G;
  uint64_t cap, steps = 0;
  bool gaveUp = false;
  GraphSearch(const Graph& g, uint64_t cap_) : G(g), cap(cap_) {}

  static int pc(const uint64_t* a, int W) { int c = 0; for (int i = 0; i < W; ++i) c += __builtin_popcountll(a[i]); return c; }
  static int lowest(const uint64_t* a, int W) {
    for (int i = 0; i < W; ++i) if (a[i]) return i * 64 + __builtin_ctzll(a[i]);
    return -1;  // callers only call on non-empty sets
  }
  static void clear1(uint64_t* a, int i) { a[i >> 6] &= ~(1ull << (i & 63)); }

  bool tick() {
    if (++steps > cap) { gaveUp = true; return false; }
    return true;
  }

  // clique(P, m): is there an m-clique inside P in the O-graph?
  bool clique(uint64_t* P, int m) {
    if (m == 0) return true;
    const int W = G.W;
    if (m == 1) return pc(P, W) > 0;
    while (pc(P, W) >= m) {
      if (!tick()) return true;  // exhaustion = possibility
      int i = lowest(P, W);
      clear1(P, i);
      uint64_t Q[32];
      const uint64_t* oi = G.o(i);
      for (int w = 0; w < W; ++w) Q[w] = P[w] & oi[w];
      if (clique(Q, m - 1)) return true;
      if (gaveUp) return true;
    }
    return false;
  }

  // recursive part of exist() (mode 1).  first = candidates for the next A vertex (> max S, O-adjacent to S),
  // second = candidates for B (> a, U-adjacent to all of S).
  bool rec(int chosen, const uint64_t* first, const uint64_t* second, bool root) {
    const int W = G.W;
    if (!tick()) return true;
    int need = 6 - chosen;
    uint64_t F[32], S2[32];
    for (int w = 0; w < W; ++w) { F[w] = first[w]; S2[w] = second[w]; }
    if (pc(F, W) < need) return false;
    if (pc(S2, W) < 6) return false;
    {  // mode-1 prune: exact 6-clique test on the second pool
      uint64_t T[32];
      for (int w = 0; w < W; ++w) T[w] = S2[w];
      if (!clique(T, 6)) return false;
      if (gaveUp) return true;
    }
    if (need == 0) return true;
    while (pc(F, W) >= need) {
      int i = lowest(F, W);
      clear1(F, i);
      uint64_t F2[32], S3[32];
      const uint64_t* oi = G.o(i);
      const uint64_t* ui = G.u(i);
      for (int w = 0; w < W; ++w) { F2[w] = F[w] & oi[w]; S3[w] = S2[w] & ui[w]; }
      if (root) clear1(S2, i);  // root: delete each tried vertex from both pools
      if (rec(chosen + 1, F2, S3, false)) return true;
      if (gaveUp) return true;
    }
    return false;
  }

  // true iff a configuration exists OR the step cap was hit.
  bool exist() {
    uint64_t all[32] = {0};
    for (int i = 0; i < G.n; ++i) all[i >> 6] |= 1ull << (i & 63);
    return rec(0, all, all, true);
  }
  // pure O 6-clique existence (the "one further basis" sanity check)
  bool hasO6() {
    uint64_t all[32] = {0};
    for (int i = 0; i < G.n; ++i) all[i >> 6] |= 1ull << (i & 63);
    return clique(all, 6);
  }
};

// Remove vertices that cannot belong to any configuration: every vertex of A needs >= 5 O-neighbours and >= 6
// U-neighbours; every vertex of B needs >= 5 O-neighbours (inside B) and >= 6 U-neighbours (all of A).  Iterated.
// Exact (deleting a vertex that lies in no configuration cannot change the answer).  Returns surviving vertex ids.
inline std::vector<int> pruneVertices(const Graph& G) {
  std::vector<char> alive(G.n, 1);
  bool changed = true;
  while (changed) {
    changed = false;
    uint64_t mask[32] = {0};
    for (int i = 0; i < G.n; ++i) if (alive[i]) mask[i >> 6] |= 1ull << (i & 63);
    for (int i = 0; i < G.n; ++i) {
      if (!alive[i]) continue;
      int dO = 0, dU = 0;
      for (int w = 0; w < G.W; ++w) {
        dO += __builtin_popcountll(G.o(i)[w] & mask[w]);
        dU += __builtin_popcountll(G.u(i)[w] & mask[w]);
      }
      if (dO < 5 || dU < 6) { alive[i] = 0; changed = true; mask[i >> 6] &= ~(1ull << (i & 63)); }
    }
  }
  std::vector<int> out;
  for (int i = 0; i < G.n; ++i) if (alive[i]) out.push_back(i);
  return out;
}

inline Graph induced(const Graph& G, const std::vector<int>& keep) {
  Graph H;
  H.init((int)keep.size());
  for (size_t a = 0; a < keep.size(); ++a)
    for (size_t b = a + 1; b < keep.size(); ++b) {
      int i = keep[a], j = keep[b];
      if ((G.o(i)[j >> 6] >> (j & 63)) & 1) H.addO((int)a, (int)b);
      if ((G.u(i)[j >> 6] >> (j & 63)) & 1) H.addU((int)a, (int)b);
    }
  return H;
}

}  // namespace mb
