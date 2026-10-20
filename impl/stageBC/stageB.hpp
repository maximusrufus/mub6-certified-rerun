// Stage B: tangent-grid refinement (SPEC section 4; refinement.tex).
// Soundness model: every bound that gates a cover decision is an UPPER bound computed either with Interval
// arithmetic or with plain doubles under FE_UPWARD on non-negative, monotone expressions (all functions here
// must run inside iv::RoundUp).  Proposals (snapping, eigen-solves, hash lookups) are heuristic and are never
// trusted: a proposed Root must pass certifyRoot() and a relation test.
#pragma once
#include "common.hpp"
#include "root.hpp"
#include "phi.hpp"
#include <cstdint>
#include <unordered_map>
#include <functional>

namespace mb {

// ------------------------------------------------------------------ relation between a parent and a second Root
struct Rel {
  bool ok = false;
  int kP = 0;
  double oP = 0, aP = 0, gP = 0, muP = 0;
  double oQ = 0, aQ = 0, gQ = 0, n0Q = 0;
  double D = 0, v0 = 0, lam = 0, e0 = 0, sigma = 0, beta = 0;
  double dqNorm = 0;
  Interval etaq[10];
  mutable int lastSteps = 0;

  // Dhat(t): upper bound on ||delta - d_q|| for every parent solution with ||eta - eta_q|| <= t (refine:iteration-*).
  // All operands are non-negative and the operations monotone, so FE_UPWARD evaluation yields an upper bound.
  double eval(double t) const {
    iv::require_upward();
    if (!(t >= 0 && t <= 2)) return 10;
    double z = std::sqrt(1 + oP) * t;
    double b = std::min(D, std::sqrt(z * z + v0 * v0));
    int steps = 0;
    for (int it = 0; it < 30; ++it) {
      double v = std::min(v0, std::min((1 + oP) * b, lam * b + e0));
      double psiP = std::sqrt(z * z + v * v);
      double ep = beta * t + sigma * v;
      double nq = n0Q + aQ * b + gQ * b * b;
      double psiQ = std::sqrt((1 + oQ) * ep * ep + nq * nq);
      double nb = std::min(b, std::min(psiP, psiQ));
      if (!(nb < b)) break;
      b = nb;
      ++steps;
    }
    lastSteps = steps;
    return b;
  }
  // refine:rho-trial / rho-return.  Returns 0 if unusable.
  double rho(double a) const {
    iv::require_upward();
    if (eval(0.0) > a - 1e-7) return 0;
    double l = 0, u = 1.01 * a;
    for (int i = 0; i < 25; ++i) {
      double m = 0.5 * (l + u);
      if (eval(m) < a - 2e-7) l = m; else u = m;
    }
    double r = std::max(0.0, l - 1e-8);
    return r > 1e-7 ? r : 0;
  }
};

// spectralUpper (refine:sigma-gram): min(max abs row sum, sqrt ||G^2||_F), for the true symmetric G inside the intervals.
inline double spectralUpper(const std::vector<Interval>& G, int n) {
  double rowmax = 0;
  for (int i = 0; i < n; ++i) {
    Interval rs(0.0);
    for (int j = 0; j < n; ++j) rs = rs + iv::abs(G[(size_t)i * n + j]);
    rowmax = std::max(rowmax, rs.hi);
  }
  Interval fro(0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      Interval s(0.0);
      for (int m = 0; m < n; ++m) s = s + G[(size_t)i * n + m] * G[(size_t)m * n + j];
      fro = fro + iv::sqr(s);
    }
  double alt = iv::sqrt(iv::sqrt(fro)).hi;  // sqrt(||G^2||_F)
  alt = std::max(alt, 0.0);
  return std::min(rowmax, alt);
}

// Build the relation from parent P (radius P.r) to Root Q.  cut: .249 (direct proposal) / .23 (reuse).
inline Rel makeRel(const Record& P, const Record& Q, double cut) {
  Rel R;
  R.kP = P.k;
  Interval d[36];
  alignedDisp(P.p, Q.p, d);
  Interval dn = inorm(d, 36);
  R.dqNorm = dn.hi;
  if (!(R.dqNorm <= cut + 1e-6)) return R;  // declined
  const int k = P.k;
  Interval eta[10];
  for (int a = 0; a < k; ++a) {
    Interval s(0.0);
    for (int m = 0; m < 36; ++m) s = s + Interval(P.col(a)[m]) * d[m];
    eta[a] = s;
    R.etaq[a] = s;
  }
  // A d = d - T T^T d
  Interval Ad[36];
  for (int m = 0; m < 36; ++m) {
    Interval s = d[m];
    for (int a = 0; a < k; ++a) s = s - Interval(P.col(a)[m]) * eta[a];
    Ad[m] = s;
  }
  double nq = inorm(Ad, 36).hi;
  // sigma >= ||S_h^T A||,  G = (A S_h)^T (A S_h)
  double sigma = 0;
  if (Q.k > 0) {
    int kq = Q.k;
    std::vector<Interval> W((size_t)kq * 36);
    for (int c = 0; c < kq; ++c) {
      Interval col[36], pc[36];
      for (int m = 0; m < 36; ++m) col[m] = Interval(Q.col(c)[m]);
      PiI(col, pc);
      Interval coef[10];
      for (int a = 0; a < k; ++a) {
        Interval s(0.0);
        for (int m = 0; m < 36; ++m) s = s + Interval(P.col(a)[m]) * pc[m];
        coef[a] = s;
      }
      for (int m = 0; m < 36; ++m) {
        Interval s = pc[m];
        for (int a = 0; a < k; ++a) s = s - Interval(P.col(a)[m]) * coef[a];
        W[(size_t)c * 36 + m] = s;
      }
    }
    std::vector<Interval> G((size_t)kq * kq);
    for (int c = 0; c < kq; ++c)
      for (int e = 0; e < kq; ++e) {
        Interval s(0.0);
        for (int m = 0; m < 36; ++m) s = s + W[(size_t)c * 36 + m] * W[(size_t)e * 36 + m];
        G[(size_t)c * kq + e] = s;
      }
    double su = spectralUpper(G, kq);
    sigma = iv::sqrt(Interval(su)).hi;
  }
  R.oP = P.ort; R.aP = P.id; R.gP = P.C; R.muP = P.mnorm;
  R.oQ = Q.ort; R.aQ = Q.id; R.gQ = Q.C; R.n0Q = Q.n0;
  const double Rr = P.r;
  R.D = (Interval(Rr) + Interval(R.dqNorm)).hi;
  double nPR = (Interval(P.n0) + Interval(P.id) * Interval(Rr) + Interval(P.C) * iv::sqr(Interval(Rr))).hi;
  R.v0 = (Interval(nPR) + Interval(nq)).hi;
  R.lam = (Interval(P.id) + Interval(P.C) * Interval(R.D)).hi;
  R.e0 = (Interval(P.mnorm) * Interval(Q.froot)).hi;
  R.sigma = sigma;
  Interval b1 = Interval(1.0) + Interval(P.ort), b2 = Interval(1.0) + Interval(Q.ort);
  R.beta = (b1 * iv::sqrt(b1 * b2)).hi;
  R.ok = true;
  return R;
}

// ------------------------------------------------------------------ atlas of certified child Roots
struct Atlas {
  std::vector<Record> recs;
  std::unordered_map<uint64_t, std::vector<int>> cells;
  static constexpr int NK = 5;
  static void keyEntries(const double* P, double* e) {
    static const int ent[NK] = {7, 8, 13, 15, 22};  // (1,1),(1,2),(2,1),(2,3),(3,4)
    double d[36];
    D0D(P, d);
    for (int i = 0; i < NK; ++i) e[i] = d[ent[i]];
  }
  static int cellOf(double x) {
    double f = x - std::floor(x);
    int c = (int)(f * 8.0);
    return c < 0 ? 0 : (c > 7 ? 7 : c);
  }
  static uint64_t keyOf(const int* c) {
    uint64_t k = 0;
    for (int i = 0; i < NK; ++i) k = k * 8 + (uint64_t)c[i];
    return k;
  }
  int insert(const Record& R) {
    int idx = (int)recs.size();
    recs.push_back(R);
    double e[NK];
    keyEntries(R.p, e);
    int c[NK];
    for (int i = 0; i < NK; ++i) c[i] = cellOf(e[i]);
    cells[keyOf(c)].push_back(idx);
    return idx;
  }
  // candidates within +-1 cell (heuristic neighbour finder)
  void near(const double* P, std::vector<int>& out) const {
    out.clear();
    double e[NK];
    keyEntries(P, e);
    int c[NK];
    for (int i = 0; i < NK; ++i) c[i] = cellOf(e[i]);
    int off[NK];
    long total = 1;
    for (int i = 0; i < NK; ++i) total *= 3;
    for (long t = 0; t < total; ++t) {
      long u = t;
      int cc[NK];
      for (int i = 0; i < NK; ++i) { off[i] = (int)(u % 3) - 1; u /= 3; cc[i] = (c[i] + off[i] + 8) & 7; }
      auto it = cells.find(keyOf(cc));
      if (it != cells.end()) out.insert(out.end(), it->second.begin(), it->second.end());
    }
  }
};

// ------------------------------------------------------------------ refinement driver
struct BStats {
  long reused = 0, parents = 0, ok = 0, failed = 0, boxes = 0, outside = 0, marked = 0, direct = 0, directFail = 0, merged = 0, inserted = 0, subdiv = 0;
};

struct Refiner {
  double tau, w;
  Atlas atlas;
  std::vector<Record> blockers;
  BStats st;
  std::vector<Phi> usedMaps;  // distinct non-identity maps used by merges (diagnostics / tests)
  std::string lastFail;
  bool verbose = false;
  bool usePhi = std::getenv("MUB_NO_PHI") == nullptr;  // measurement switch only (Phi reuse is optional, never load-bearing)

  Refiner(double tau_, double w_) : tau(tau_), w(w_) {}

  // per-parent grid context
  struct Grid {
    const Record* P = nullptr;
    int k = 0;
    long L = 0, n = 0;
    uint64_t total = 1;
    double Ed = 0;  // certified upper bound for sqrt(1+ort)*R (inflated by 1e-7 like the paper)
    std::vector<uint64_t> marks;
    // Parent frames (refine:mapped-coordinates): PV[0] = P itself, PV[v] = Phi_v P (same grid coordinates, since
    // (Phi T)^T (Phi delta) = T^T delta).  A relation from Phi_v P to a stored child/blocker Z is the relation from
    // P to the pulled-back Root Phi_v^{-1} Z (isometry), so every mark is still certified by makeRel + markRec.
    std::vector<Record> PV;
    std::unordered_map<uint64_t, Rel> relCache;  // key (variant v, child zi): relation PV[v] -> child (box independent)
    std::unordered_map<uint64_t, char> childDone; // markByChild already run for (v, zi)
    static uint64_t key(int v, int zi) { return ((uint64_t)(uint32_t)zi << 8) | (uint64_t)(v & 255); }
    std::vector<double> comp;  // complement-of-T basis used by the snapper (normal directions only)
    bool isMarked(uint64_t idx) const { return (marks[idx >> 6] >> (idx & 63)) & 1; }
    void mark(uint64_t idx) { marks[idx >> 6] |= (1ull << (idx & 63)); }
  };

  bool setupGrid(const Record& P, Grid& g) {
    g.P = &P;
    g.k = P.k;
    complementBasis(P.T.data(), P.k, g.comp);
    g.PV.clear();
    if (usePhi) {
      std::vector<Phi> vs = reuseVariants(P.p);
      if (vs.size() > 200) vs.resize(200);
      for (auto& f : vs) g.PV.push_back(applyPhi(P, f));
    } else {
      g.PV.push_back(P);
    }
    Interval E = iv::sqrt(Interval(1.0) + Interval(P.ort)) * Interval(P.r);
    g.Ed = (E + Interval(1e-7)).hi;
    double q = (Interval(g.Ed) / Interval(w)).hi + 0.500001;
    long L = (long)std::ceil(q);
    for (;;) {  // verify (L + 1/2) w >= Ed rigorously
      Interval chk = (Interval((double)L) + Interval(0.5)) * Interval(w);
      if (chk.lo >= g.Ed) break;
      ++L;
    }
    g.L = L;
    g.n = 2 * L + 1;
    uint64_t total = 1;
    for (int a = 0; a < g.k; ++a) {
      if (total > 50000000ull / (uint64_t)g.n) { lastFail = "base grid exceeds 50,000,000 boxes"; return false; }
      total *= (uint64_t)g.n;
    }
    g.total = total;
    g.marks.assign((total + 63) / 64, 0);
    return true;
  }

  // markRec: mark every base box whose farthest-corner tangent distance to eta_q is certified <= rho.
  // Returns number of boxes (marked or, if !commit, that would be marked) and tmax (certified upper bound).
  long markRec(const Rel& rel, double rho, Grid& g, bool commit, double& tmax) {
    iv::require_upward();
    const int k = g.k;
    const long n = g.n, L = g.L;
    tmax = 0;
    std::vector<std::vector<double>> tab(k, std::vector<double>(n));
    Interval hw = Interval(w) * Interval(0.5);
    for (int a = 0; a < k; ++a)
      for (long b = -L; b <= L; ++b) {
        Interval u = Interval(w) * Interval((double)b);
        Interval dev = iv::abs(u - rel.etaq[a]) + hw + Interval(1e-8);
        tab[a][b + L] = iv::sqr(dev).hi;
      }
    double thr = (rho - 5e-8) * (rho - 5e-8);
    if (!(rho - 5e-8 > 0)) return 0;
    long count = 0;
    std::vector<long> bsel(k, 0);
    std::vector<uint64_t> pend;
    // iterative DFS
    struct Rec {
      static void go(int a, int k, long n, double sum, double thr, const std::vector<std::vector<double>>& tab, std::vector<long>& bsel,
                     const std::function<void(double)>& leaf) {
        if (a == k) { leaf(sum); return; }
        for (long bi = 0; bi < n; ++bi) {
          double s2 = sum + tab[a][bi];
          if (s2 <= thr) { bsel[a] = bi; go(a + 1, k, n, s2, thr, tab, bsel, leaf); }
        }
      }
    };
    std::vector<std::pair<uint64_t, double>> hits;
    auto leaf = [&](double s) {
      double t = std::sqrt(s) + 2e-8;
      if (!(t <= rho)) return;
      uint64_t idx = 0, mul = 1;
      for (int a = 0; a < k; ++a) { idx += (uint64_t)bsel[a] * mul; mul *= (uint64_t)n; }
      hits.emplace_back(idx, t);
    };
    Rec::go(0, k, n, 0.0, thr, tab, bsel, leaf);
    for (auto& h : hits) {
      if (g.isMarked(h.first)) continue;
      ++count;
      tmax = std::max(tmax, h.second);
      if (commit) g.mark(h.first);
    }
    return count;
  }

  // Mark by a blocker relation (no radius growth).
  void markByBlocker(const Record& Pr, const Record& B, Grid& g) {
    Rel rel = makeRel(Pr, B, 0.23);
    if (!rel.ok) return;
    double rho = rel.rho(B.r);
    if (rho <= 0) return;
    double tmax;
    long c = markRec(rel, rho, g, true, tmax);
    st.marked += c;
  }
  const Rel& getRel(const Record& Pf, int v, int zi, Grid& g, double cut) {
    uint64_t key = Grid::key(v, zi);
    auto it = g.relCache.find(key);
    if (it == g.relCache.end()) it = g.relCache.emplace(key, makeRel(Pf, atlas.recs[zi], cut)).first;
    return it->second;
  }
  // Mark by an atlas child Z through the relation `rel` (parent frame -> Z); Z's radius grows to the certified
  // Dhat(tmax)+1e-8, which must be <= tau.  Returns boxes marked.
  long markByRel(const Rel& rel, int zi, Grid& g) {
    if (!rel.ok) return 0;
    double rho = rel.rho(tau);
    if (verbose) std::fprintf(stderr, "  child %d: dq=%.4f sigma=%.3f lam=%.3f v0=%.4f D=%.4f Dhat(.01)=%.4f Dhat(.02)=%.4f Dhat(.04)=%.4f rho=%.4f\n", zi, rel.dqNorm, rel.sigma, rel.lam, rel.v0, rel.D, rel.eval(.01), rel.eval(.02), rel.eval(.04), rho);
    if (rho <= 0) return 0;
    double tmax;
    long c = markRec(rel, rho, g, false, tmax);
    if (verbose) std::fprintf(stderr, "    boxes %ld tmax=%.4f newR=%.4f\n", c, tmax, c ? rel.eval(tmax) + 1e-8 : 0.0);
    if (c == 0) return 0;
    double newR = rel.eval(tmax) + 1e-8;
    if (!(tmax <= rho && newR <= tau)) return 0;  // refuse (the paper aborts; we simply do not mark)
    markRec(rel, rho, g, true, tmax);
    atlas.recs[zi].r = std::max(atlas.recs[zi].r, newR);
    st.marked += c;
    return c;
  }
  // Mark by atlas child zi seen from parent frame v (cached relation; run once per (v, zi)).
  void markByChild(int v, int zi, Grid& g, double cut) {
    uint64_t key = Grid::key(v, zi);
    if (g.childDone.count(key)) return;
    g.childDone[key] = 1;
    markByRel(getRel(g.PV[v], v, zi, g, cut), zi, g);
  }
  // Mark by atlas child zi pulled back by Phi^{-1} into the parent's own frame (the paper's inverseTransform /
  // addNear): a fresh relation P -> Phi^{-1} Z.  Not cached (used once, after a merge through Phi).
  void markByPulledBack(const Record& P, int zi, const Phi& f, Grid& g, double cut) {
    if (f.isIdentity()) { markByChild(0, zi, g, cut); return; }
    Record Zb = applyPhi(atlas.recs[zi], f.inverse());
    Rel rel = makeRel(P, Zb, cut);
    markByRel(rel, zi, g);
  }

  // assign (refine:atlas-assignment): merge into an existing centre or insert.  Returns atlas index.
  // record a non-identity map used by assign (tests re-check coverage through exactly these maps)
  void noteMap(const Phi& f) {
    if (f.isIdentity()) return;
    for (auto& m : usedMaps) { bool same = m.sign == f.sign; for (int i = 0; i < 36 && same; ++i) same = m.map[i] == f.map[i]; if (same) return; }
    usedMaps.push_back(f);
  }
  // The map used for the merge / insertion (identity if none) is returned through `used`.
  int assign(const Record& Q, double r, Phi* used = nullptr) {
    std::vector<Phi> vs;
    if (usePhi) vs = reuseVariants(Q.p); else vs.push_back(Phi::identity());
    int best = -1, bestV = 0;
    double bestL = 1e300;
    std::vector<int> cand;
    for (size_t v = 0; v < vs.size(); ++v) {
      double Qp[36];
      applyPhiVec(vs[v], Q.p, Qp);
      atlas.near(Qp, cand);
      for (int zi : cand) {
        double a = approxAligned(atlas.recs[zi].p, Qp);
        if (a < bestL) { bestL = a; best = zi; bestV = (int)v; }
      }
    }
    if (best >= 0) {
      double Qp[36];
      applyPhiVec(vs[bestV], Q.p, Qp);
      Interval d[36];
      alignedDisp(atlas.recs[best].p, Qp, d);
      double l = inorm(d, 36).hi;  // certified aligned distance of Phi Q from Z
      double tot = (Interval(r) + Interval(l) + Interval(2e-8)).hi;
      if (tot < tau) {
        double nr = (Interval(r) + Interval(l) + Interval(1e-8)).hi;
        atlas.recs[best].r = std::max(atlas.recs[best].r, nr);
        ++st.merged;
        if (used) *used = vs[bestV];
        noteMap(vs[bestV]);
        return best;
      }
    }
    // Insert the CANONICAL transformed Root (first signature map; identity if none), as the paper does: stored
    // centres then share the signature frame, so a later query's variants can meet them.  Any Phi is allowed
    // (phase:transport); the caller pulls the stored root back by Phi^{-1} (exact) for its own relation.
    Phi canon = (vs.size() > 1) ? vs[1] : Phi::identity();
    Record Z = applyPhi(Q, canon);
    if (!canon.isIdentity()) {
      // Re-certify in the stored orientation so that every stored constant is a fresh interval bound evaluated on
      // the stored doubles (the transported constants are valid by phase:transport; fresh ones are kept only
      // when the certificate succeeds with the same k, otherwise the transported ones stay).
      Record F;
      if (certifyRoot(Z.p, F) && F.k == Z.k) Z = F;
    }
    Z.r = r;
    ++st.inserted;
    if (used) *used = canon;
    noteMap(canon);
    return atlas.insert(Z);
  }

  // Box-level reuse of an existing child: Dhat at the farthest corner of the box, radius grown (<= tau).
  bool tryReuse(const Record& Pr, const std::vector<double>& s, double hs, Grid& g) {
    const int k = Pr.k;
    double Q0[36];
    for (int m = 0; m < 36; ++m) Q0[m] = Pr.p[m];
    for (int a = 0; a < k; ++a) {
      double u = w * s[a];
      for (int m = 0; m < 36; ++m) Q0[m] += Pr.col(a)[m] * u;
    }
    std::vector<int> cand;
    atlas.near(Q0, cand);
    Interval hbox = Interval(w) * Interval(hs);
    int best = -1;
    double bestR = 1e300;
    for (int zi : cand) {
      const Rel& rel = getRel(g.PV[0], 0, zi, g, 0.249);
      if (!rel.ok) continue;
      Interval sq(0.0);
      for (int a = 0; a < k; ++a) {
        Interval u = Interval(w) * Interval(s[a]);
        sq = sq + iv::sqr(iv::abs(u - rel.etaq[a]) + hbox);
      }
      double tB = iv::sqrt(sq).hi;
      double nr = rel.eval(tB) + 1e-8;
      if (nr <= tau && nr < bestR) { bestR = nr; best = zi; }
    }
    if (best < 0) return false;
    atlas.recs[best].r = std::max(atlas.recs[best].r, bestR);
    ++st.reused;
    return true;
  }

  // direct proposal for a box (centre s_a * w, half-width hs * w).  Returns true if the box is covered.
  bool direct(const Record& Pr, const std::vector<double>& s, double hs, Grid& g) {
    ++st.direct;
    const int k = Pr.k;
    double Q0[36];
    for (int m = 0; m < 36; ++m) Q0[m] = Pr.p[m];
    for (int a = 0; a < k; ++a) {
      double u = w * s[a];
      for (int m = 0; m < 36; ++m) Q0[m] += Pr.col(a)[m] * u;
    }
    if (!snapSolution(Q0, &g.comp)) { ++st.directFail; if (verbose) std::fprintf(stderr, "  direct: snap failed\n"); return false; }
    for (int m = 0; m < 36; ++m) Q0[m] -= rnd(Q0[m]);  // equivalent phases in [-1/2,1/2]
    Record Q;
    RootInfo info;
    if (!certifyRoot(Q0, Q, &info)) { ++st.directFail; if (verbose) std::fprintf(stderr, "  direct: certify failed: %s\n", info.why.c_str()); return false; }
    Rel rel = makeRel(Pr, Q, 0.249);
    if (!rel.ok) { ++st.directFail; if (verbose) std::fprintf(stderr, "  direct: relation declined dq=%g\n", rel.dqNorm); return false; }
    Interval hbox = Interval(w) * Interval(hs);
    Interval sq(0.0);
    for (int a = 0; a < k; ++a) {
      Interval u = Interval(w) * Interval(s[a]);
      Interval dv = iv::abs(u - rel.etaq[a]) + hbox;
      sq = sq + iv::sqr(dv);
    }
    double tB = iv::sqrt(sq).hi;
    double r = rel.eval(tB) + 2e-8;
    if (!(r < tau - 1e-7)) { ++st.directFail; if (verbose) std::fprintf(stderr, "  direct: r=%g tB=%g tau=%g D=%g v0=%g lam=%g sigma=%g beta=%g\n", r, tB, tau, rel.D, rel.v0, rel.lam, rel.sigma, rel.beta); return false; }
    Phi used;
    int zi = assign(Q, r, &used);
    markByPulledBack(Pr, zi, used, g, 0.249);
    return true;
  }

  bool outside(const Grid& g, const std::vector<double>& s, double hs) {
    Interval hbox = Interval(w) * Interval(hs);
    Interval sum(0.0);
    for (int a = 0; a < g.k; ++a) {
      Interval u = Interval(w) * Interval(s[a]);
      Interval x = iv::abs(u) - hbox;
      double lo = std::max(0.0, x.lo);
      sum = sum + iv::sqr(Interval(lo));
    }
    return sum.lo > (Interval(g.Ed) * Interval(g.Ed)).hi;
  }

  bool solveBox(const Record& Pr, std::vector<double> s, double hs, int depth, Grid& g) {
    ++st.boxes;
    if (g.k > 0 && outside(g, s, hs)) { ++st.outside; return true; }
    if (tryReuse(Pr, s, hs, g)) return true;
    if (direct(Pr, s, hs, g)) return true;
    if (g.k == 0) { lastFail = "unresolved cell with k=0"; return false; }
    if (depth >= 6) { lastFail = "unresolved cell at depth 6"; return false; }
    ++st.subdiv;
    int k = g.k;
    for (int c = 0; c < (1 << k); ++c) {
      std::vector<double> s2(k);
      for (int a = 0; a < k; ++a) s2[a] = s[a] + ((c >> a) & 1 ? 0.5 : -0.5) * hs;
      if (!solveBox(Pr, s2, hs * 0.5, depth + 1, g)) return false;
    }
    return true;
  }

  // Process one parent.  Returns true iff its whole sphere is now covered by the atlas (or a blocker).
  bool processParent(const Record& P) {
    ++st.parents;
    lastFail.clear();
    Grid g;
    if (!setupGrid(P, g)) { ++st.failed; return false; }
    for (size_t v = 0; v < g.PV.size(); ++v)
      for (const Record& B : blockers) markByBlocker(g.PV[v], B, g);
    // existing children: scan every parent frame with a cheap entry prefilter
    for (size_t v = 0; v < g.PV.size(); ++v) {
      double ep[Atlas::NK];
      Atlas::keyEntries(g.PV[v].p, ep);
      for (size_t zi = 0; zi < atlas.recs.size(); ++zi) {
        double ez[Atlas::NK];
        Atlas::keyEntries(atlas.recs[zi].p, ez);
        bool far = false;
        for (int i = 0; i < Atlas::NK; ++i) {
          double x = ez[i] - ep[i];
          x -= rnd(x);
          if (std::fabs(x) > 0.4) { far = true; break; }
        }
        if (!far) markByChild((int)v, (int)zi, g, 0.249);
      }
    }
    std::vector<double> s(g.k);
    for (uint64_t idx = 0; idx < g.total; ++idx) {
      if (g.isMarked(idx)) continue;
      uint64_t u = idx;
      for (int a = 0; a < g.k; ++a) { s[a] = (double)((long)(u % (uint64_t)g.n) - g.L); u /= (uint64_t)g.n; }
      if (!solveBox(P, s, 0.5, 0, g)) { ++st.failed; return false; }
    }
    ++st.ok;
    blockers.push_back(P);
    return true;
  }
};

}  // namespace mb
