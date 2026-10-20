// bdiag: measurement tool for PERF_FIX.md.  ./bdiag records.bin [tau]
#include "../phi.hpp"
// For every record: fresh certifyRoot at P -> plain C vs tensor C; and the relation-iteration shape (rho for a
// self-relation at radius tau) under both constants.  Diagnostics only (no certificate is emitted).
#include "../stageB.hpp"
#include "../phi.hpp"
#include <cstdio>

using namespace mb;

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: bdiag records.bin [tau]\n"); return 2; }
  double tau = argc > 2 ? std::atof(argv[2]) : 0.07;
  iv::RoundUp g;
  std::vector<Record> in;
  std::string err;
  if (!readRecords(argv[1], in, err)) { std::fprintf(stderr, "%s\n", err.c_str()); return 2; }
  if (argc > 3 && std::string(argv[3]) == "pairs") {
    // For each record j: nearest earlier record i under identity vs under the Phi variants (both sides),
    // and how many j would be within tau (merge/reuse range) of an earlier one.
    std::vector<std::vector<Phi>> V(in.size());
    for (size_t j = 0; j < in.size(); ++j) V[j] = reuseVariants(in[j].p);
    int nId = 0, nPhi = 0, nIdR = 0, nPhiR = 0; double sId = 0, sPhi = 0; size_t nv = 0;
    for (size_t j = 1; j < in.size(); ++j) {
      double bId = 1e300, bPhi = 1e300;
      nv += V[j].size();
      for (size_t i = 0; i < j; ++i) {
        bId = std::min(bId, approxAligned(in[j].p, in[i].p));
        for (auto& f : V[j]) { double Pj[36]; applyPhiVec(f, in[j].p, Pj);
          for (auto& h : V[i]) { double Pi[36]; applyPhiVec(h, in[i].p, Pi); bPhi = std::min(bPhi, approxAligned(Pj, Pi)); } }
      }
      sId += bId; sPhi += bPhi;
      if (bId < tau) ++nId;
      if (bPhi < tau) ++nPhi;
      if (bId < in[j].r + 0.125) ++nIdR;   // balls overlap (R_j + R_i <= .25)
      if (bPhi < in[j].r + 0.125) ++nPhiR;
    }
    std::printf("pairs: n=%zu avg variants/record %.2f; nearest-earlier mean dist identity %.4f / with Phi %.4f; within tau=%.3f: identity %d, Phi %d; balls overlapping (d<r_j+.125): identity %d, Phi %d\n",
                in.size(), (double)nv / (in.size() - 1), sId / (in.size() - 1), sPhi / (in.size() - 1), tau, nId, nPhi, nIdR, nPhiR);
    return 0;
  }
  if (argc > 3 && std::string(argv[3]) == "rel") {
    // pairwise relations among the first records (identity Phi): dq, sigma (certified) vs plain-double principal-angle
    // sine of the tangent spaces, lam, rho at the blocker radius, and rho under the Phi variants.
    size_t nn = std::min<size_t>(in.size(), argc > 4 ? std::atoi(argv[4]) : 12);
    std::printf("# j i : dq sigma sigma_true lam v0 rho(a=R_i) | best over Phi: dq rho\n");
    for (size_t j = 1; j < nn; ++j)
      for (size_t i = 0; i < j; ++i) {
        Rel r = makeRel(in[j], in[i], 0.23);
        double st = 0;
        if (r.ok) {
          // power iteration on S^T A A S with A = I - T T^T (heuristic check of sigma)
          int kq = in[i].k, k = in[j].k;
          std::vector<double> G((size_t)kq * kq, 0.0);
          for (int c = 0; c < kq; ++c)
            for (int e = 0; e < kq; ++e) {
              double sum = 0;
              for (int m = 0; m < 36; ++m) {
                double wc = in[i].col(c)[m], we = in[i].col(e)[m];
                for (int a = 0; a < k; ++a) {
                  double tc = 0, te = 0;
                  for (int q = 0; q < 36; ++q) { tc += in[j].col(a)[q] * in[i].col(c)[q]; te += in[j].col(a)[q] * in[i].col(e)[q]; }
                  wc -= in[j].col(a)[m] * tc; we -= in[j].col(a)[m] * te;
                }
                sum += wc * we;
              }
              G[(size_t)c * kq + e] = sum;
            }
          std::vector<double> V, ev;
          jacobiEig(G, kq, V, ev);
          st = std::sqrt(std::max(0.0, ev[kq - 1]));
        }
        double bestRho = -1, bestDq = 9;
        for (auto& f : reuseVariants(in[j].p)) {
          Record Pj = applyPhi(in[j], f);
          Rel rr = makeRel(Pj, in[i], 0.23);
          if (!rr.ok) continue;
          double rh = rr.rho(in[i].r);
          if (rh > bestRho) { bestRho = rh; bestDq = rr.dqNorm; }
        }
        std::printf("%zu %zu : %.4f %.3f %.3f %.3f %.4f %.4f | %.4f %.4f\n", j, i, r.dqNorm, r.sigma, st, r.lam, r.v0, r.ok ? r.rho(in[i].r) : -1, bestDq, bestRho);
      }
    return 0;
  }
  double sumRatio = 0; int n = 0, nT = 0;
  std::vector<double> ratios;
  std::printf("# idx k r C_stored C_plain C_tensor beta row plainRow mnorm | rho(self,tau) with C_stored and with C_tensor\n");
  for (size_t i = 0; i < in.size(); ++i) {
    Record R;
    RootInfo info;
    bool ok = certifyRoot(in[i].p, R, &info);
    if (!ok) { std::printf("%zu certify failed: %s\n", i, info.why.c_str()); continue; }
    // relation of the parent to a copy of itself: Dhat(t) shape and rho at radius tau
    Record P = in[i];
    Record Pt = in[i]; Pt.C = std::min(info.Cplain, info.Ctensor);
    Rel r1 = makeRel(P, P, 0.249), r2 = makeRel(Pt, Pt, 0.249);
    double rho1 = r1.ok ? r1.rho(tau) : -1, rho2 = r2.ok ? r2.rho(tau) : -1;
    std::printf("%zu %d %.4f %.4f %.4f %.4f %.3g %.4g %.4g %.4f | %.4f %.4f\n", i, in[i].k, in[i].r, in[i].C, info.Cplain, info.Ctensor, 0.0, 0.0, 0.0, in[i].mnorm, rho1, rho2);
    if (std::isfinite(info.Ctensor)) { sumRatio += info.Ctensor / info.Cplain; ++nT; ratios.push_back(info.Ctensor / info.Cplain); }
    ++n;
  }
  std::sort(ratios.begin(), ratios.end());
  if (nT) std::printf("# tensor usable %d/%d, C_tensor/C_plain mean %.3f median %.3f min %.3f max %.3f\n", nT, n, sumRatio / nT, ratios[ratios.size() / 2], ratios.front(), ratios.back());
  return 0;
}
