// cover.cpp -- from a retained U(2) cell to a certified phase ball (SPEC 3.7).
// Candidate::prep (arc bounds, center P0, eps), covercompute (direct bound s0,
// tangent bounds via ellipsoids + ell), bootstrap, local::test.  All interval.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include "cover.hpp"
#include "la.hpp"
#include "stagea.hpp"

namespace sa {

namespace {

struct CoverData {
  I mj[16];
  Z Dc[16];
  I Xd[16][13];
  double P0[36];
  double eps = 0, arc = 0, weighted = 0;
  double invmag[16];
};

bool prepCover(const SeedCtx& S, const Params& P, CellWork& w, CoverData& cd) {
  (void)P;
  ensureJ(S, w);
  const LevelK& lk = S.L(w.Nu);
  I sumSq(0.0);
  double mi = 1e300;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      int f = 4 * i + j;
      cd.mj[f] = zabs(w.D.a[i][j]);
      if (!(cd.mj[f].lo > 0.15)) return false;  // magnitude guard (A24)
      mi = std::min(mi, cd.mj[f].lo);
      sumSq = sumSq + sqr(cd.mj[f] - I(1.0));
    }
  I slackI = I(lk.z2) - I(sumSq.lo);
  double slack = std::max(0.0, slackI.hi);
  double qd = (I(slack) / I(mi)).hi;
  if (!(qd < 3.8)) return false;
  I den = sqrt0(I(1.0) - mulp(0.25, I(qd)));
  cd.arc = (sqrt0(I(qd)) / den).hi;
  cd.weighted = (sqrt0(I(slack)) / den).hi;
  if (!(cd.arc < 8.8) || !(cd.weighted < 10.0)) return false;
  for (int p = 0; p < 36; ++p) cd.P0[p] = 0.0;
  for (int l = 0; l < 9; ++l) cd.P0[SEEDPOS[l]] = S.q[l];
  I e2(0.0);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      int f = 4 * i + j;
      cd.invmag[f] = recip(I(cd.mj[f].lo)).hi;
      cd.Dc[f] = zmulI(w.D.a[i][j], recip(cd.mj[f]));
      double ang = std::atan2(midp(w.D.a[i][j].im), midp(w.D.a[i][j].re)) / (2.0 * M_PI);
      cd.P0[dpos(f)] = ang;
      Z e = iv::cis_turns(I(ang));
      I kap = sqrt0(abs2(e - cd.Dc[f]));
      e2 = e2 + sqr(mulp(0.25, kap));
    }
  cd.eps = sqrt0(e2).hi;
  if (!(cd.eps < 1e-3)) return false;
  for (int f = 0; f < 16; ++f) {
    int i = f / 4, j = f % 4;
    for (int l = 0; l < 13; ++l) {
      const Z& J = w.Jc[l].a[i][j];
      cd.Xd[f][l] = cd.Dc[f].re * J.im - cd.Dc[f].im * J.re;
    }
  }
  return true;
}

// ---- spectral upper bound of a symmetric exact-double k x k PSD matrix
double spectralUpper(const double* G, int k) {
  double best = 1e300;
  {
    double mx = 0;
    for (int a = 0; a < k; ++a) {
      I s(0.0);
      for (int b = 0; b < k; ++b) s = s + iv::abs(I(G[a * k + b]));
      mx = std::max(mx, s.hi);
    }
    best = mx;
  }
  // tight bound through an approximate eigenbasis: lam <= rowsum(V^T G V)/(1-beta)
  try {
    std::vector<double> A(G, G + k * k), ev, V;
    la::jacobi(A, k, ev, V);
    std::vector<I> GV((size_t)k * k), VGV((size_t)k * k);
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) {
        I s(0.0);
        for (int c = 0; c < k; ++c) s = s + mulp(G[a * k + c], I(V[c * k + b]));
        GV[(size_t)a * k + b] = s;
      }
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) {
        I s(0.0);
        for (int c = 0; c < k; ++c) s = s + mulp(V[c * k + a], GV[(size_t)c * k + b]);
        VGV[(size_t)a * k + b] = s;
      }
    double r = 0;
    for (int a = 0; a < k; ++a) {
      I s(0.0);
      for (int b = 0; b < k; ++b) s = s + iv::abs(VGV[(size_t)a * k + b]);
      r = std::max(r, s.hi);
    }
    I bsum(0.0);
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) {
        I s(0.0);
        for (int c = 0; c < k; ++c) s = s + mulp(V[c * k + a], I(V[c * k + b]));
        if (a == b) s = s - I(1.0);
        bsum = bsum + sqr(s);
      }
    double beta = sqrt0(bsum).hi;
    if (beta < 0.5) {
      double t = (I(r) / (I(1.0) - I(beta))).hi;
      best = std::min(best, t);
    }
  } catch (...) {
  }
  return best;
}

// certified upper bound on max over the offsets of ||U + G^{1/2} w||, |w|<=1.
double ell(const double* G, int k, const std::vector<std::array<I, 10>>& corners) {
  double gu = spectralUpper(G, k);
  double mx = 0;
  for (auto& u : corners) {
    I s(0.0);
    for (int a = 0; a < k; ++a) s = s + sqr(u[a]);
    mx = std::max(mx, sqrt0(s).hi);
  }
  double best = (I(mx) + sqrt0(I(gu))).hi;
  if (!(gu < 400.0)) return best;
  const double xs[] = {1e-4, 1e-3, 1e-2, 0.05, 0.2, 1.0, 4.0};
  for (double x : xs) {
    double gam = (I(gu) * I(1.0 + x)).hi + 2e-6;
    if (!(gam < 418.0)) continue;
    I gapI = I(gam) - I(gu);
    if (!(gapI.lo >= 1e-6)) continue;
    // matrix A = gam I - G (doubles) for the heuristic solve
    std::vector<double> Am(k * k);
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) Am[a * k + b] = (a == b ? gam : 0.0) - G[a * k + b];
    double high = 0;
    bool ok = true;
    for (auto& u : corners) {
      double um[10], z[10];
      for (int a = 0; a < k; ++a) um[a] = midp(u[a]);
      if (!la::chol_solve(Am.data(), k, 0.0, um, z)) { ok = false; break; }
      double sz = 1.0;
      for (int a = 0; a < k; ++a) {
        if (!(std::fabs(z[a]) < 1e8)) ok = false;
        sz = std::max(sz, std::fabs(z[a]));
      }
      if (!ok) break;
      I dot(0.0), un(0.0), res(0.0);
      for (int a = 0; a < k; ++a) {
        dot = dot + mulp(z[a], u[a]);
        un = un + sqr(u[a]);
        I Az = mulp(gam, I(z[a]));
        for (int b = 0; b < k; ++b) Az = Az - mulp(G[a * k + b], I(z[b]));
        res = res + sqr(u[a] - Az);
      }
      I bound = dot + sqrt0(un) * sqrt0(res) / gapI;
      high = std::max(high, bound.hi);
    }
    if (!ok) continue;
    double cand = sqrt0(I(gam) * (I(1.0) + I(high))).hi;
    best = std::min(best, cand);
  }
  return best;
}

// symmetric exact-double matrix dominating an interval matrix (entries a<=b given)
void dominate(const std::vector<I>& Gi, int k, double* Gd) {
  double rho = 0;
  std::vector<double> dev(k * k);
  for (int a = 0; a < k; ++a)
    for (int b = a; b < k; ++b) {
      const I& e = Gi[(size_t)a * k + b];
      double m = midp(e);
      Gd[a * k + b] = Gd[b * k + a] = m;
      double d = std::max((I(e.hi) - I(m)).hi, (I(m) - I(e.lo)).hi);
      dev[a * k + b] = dev[b * k + a] = d;
    }
  for (int a = 0; a < k; ++a) {
    I s(0.0);
    for (int b = 0; b < k; ++b) s = s + I(dev[a * k + b]);
    rho = std::max(rho, s.hi);
  }
  for (int a = 0; a < k; ++a) Gd[a * k + a] = (I(Gd[a * k + a]) + I(rho)).hi;
}

// Certified radius bound for the cell w.r.t. Root R.  Returns false on any guard failure.
bool coverBound(const Root& R, const CoverData& cd, const SeedCtx& S, const LevelK& lk, double target, double& sOut) {
  try {
    if (!(R.C < 20.0)) return false;
    const I RsI(S.Rs), hI(lk.h), ErI(lk.Er), arcI(cd.arc), epsI(cd.eps);
    I Delta[36];
    alignDelta(cd.P0, R.P, Delta);
    I nd(0.0);
    for (int p = 0; p < 36; ++p) nd = nd + sqr(Delta[p]);
    I Ds(0.0), Dd(0.0);
    for (int l = 0; l < 9; ++l) Ds = Ds + sqr(Delta[SEEDPOS[l]]);
    double dmaxabs = 0;
    for (int f = 0; f < 16; ++f) {
      Dd = Dd + sqr(Delta[dpos(f)]);
      dmaxabs = std::max(dmaxabs, iv::abs(Delta[dpos(f)]).hi);
    }
    I Rs2arc2 = sqr(RsI) + sqr(arcI);
    I arc3_6 = (sqr(arcI) * arcI) * recip(I(6.0));
    // XtDelta
    I g[13];
    for (int l = 0; l < 13; ++l) {
      I s = (l < 9) ? Delta[SEEDPOS[l]] : I(0.0);
      for (int f = 0; f < 16; ++f) s = s + cd.Xd[f][l] * Delta[dpos(f)];
      g[l] = s;
    }
    I gs(0.0), gt(0.0);
    for (int l = 0; l < 9; ++l) gs = gs + sqr(g[l]);
    for (int l = 9; l < 13; ++l) gt = gt + iv::abs(g[l]);
    I S1 = RsI * sqrt0(gs) + hI * gt + ErI * sqrt0(Dd) + I(dmaxabs) * arc3_6;
    I w2(0.0);
    for (int f = 0; f < 16; ++f) w2 = w2 + sqr(Delta[dpos(f)]) * I(cd.invmag[f]);
    I S2 = RsI * sqrt0(Ds) + I(cd.weighted) * sqrt0(w2);
    double supp = std::min(S1.hi, S2.hi);
    I d2 = Rs2arc2 / sqr(TWOPI());
    I inner = nd + d2 + I(supp) / PI();
    I s0I = sqrt0(I(inner.lo < 0 ? 0.0 : inner.lo, inner.hi)) + epsI;
    I simple = sqrt0(nd) + epsI + sqrt0(Rs2arc2) / TWOPI();
    double s = std::min(s0I.hi, simple.hi);
    if (!(s < 3.0)) return false;
    // Hopeless-candidate gate (speed only; returning false is never an exclusion): the bootstrap can only
    // shrink s when n_P(s) < s; if n_P(s)/s >= 1 already, nothing improves, so s > target fails.
    if (s > target) {
      I ratio = I(R.C) * I(s) + I(R.id) + I(R.n0) / I(s);
      if (ratio.lo >= 1.0) { sOut = s; return true; }
    }
    // tangent
    double t0 = 0;
    const int k = R.k;
    if (k > 0) {
      std::vector<I> AB((size_t)k * 13);
      for (int a = 0; a < k; ++a)
        for (int l = 0; l < 13; ++l) {
          I v = (l < 9) ? I(R.T[a][SEEDPOS[l]]) : I(0.0);
          for (int f = 0; f < 16; ++f) v = v + mulp(R.T[a][dpos(f)], cd.Xd[f][l]);
          AB[(size_t)a * 13 + l] = v;
        }
      I U0[10];
      for (int a = 0; a < k; ++a) {
        I v(0.0);
        for (int p = 0; p < 36; ++p) v = v + mulp(R.T[a][p], Delta[p]);
        U0[a] = v * TWOPI();
      }
      double maxrow = 0;
      for (int f = 0; f < 16; ++f) {
        I r2(0.0);
        for (int a = 0; a < k; ++a) r2 = r2 + sqr(I(R.T[a][dpos(f)]));
        maxrow = std::max(maxrow, sqrt0(r2).hi);
      }
      I sine = I(maxrow) * arc3_6;
      I holterm = I(R.hor) * sqrt0(Rs2arc2);
      I sq1ort = sqrt0(I(1.0) + I(R.ort));
      double bestT = 1e300;
      // ---- linear variant (beta = .45)
      {
        const double beta = 0.45;
        double c1 = ((I(1.0) + I(beta)) * sqr(RsI)).hi;
        double c2 = ((I(1.0) + recip(I(beta))) * sqr(ErI)).hi;
        std::vector<I> Gi((size_t)k * k);
        for (int a = 0; a < k; ++a)
          for (int b = a; b < k; ++b) {
            I s1(0.0), s2(0.0);
            for (int l = 0; l < 9; ++l) s1 = s1 + AB[(size_t)a * 13 + l] * AB[(size_t)b * 13 + l];
            for (int f = 0; f < 16; ++f) s2 = s2 + mulp(R.T[a][dpos(f)], I(R.T[b][dpos(f)]));
            Gi[(size_t)a * k + b] = mulp(c1, s1) + mulp(c2, s2);
          }
        std::vector<double> Gd((size_t)k * k);
        dominate(Gi, k, Gd.data());
        std::vector<std::array<I, 10>> corners;
        for (int m = 0; m < 16; ++m) {
          std::array<I, 10> u;
          for (int a = 0; a < k; ++a) {
            I v = U0[a];
            for (int b = 0; b < 4; ++b) {
              double sg = (m >> b & 1) ? 1.0 : -1.0;
              v = v + mulp(sg, hI * AB[(size_t)a * 13 + 9 + b]);
            }
            u[a] = v;
          }
          corners.push_back(u);
        }
        double E = ell(Gd.data(), k, corners);
        I t = (I(E) + sine + holterm) / TWOPI() + sq1ort * epsI;
        bestT = std::min(bestT, t.hi);
      }
      // ---- direct variant (beta = .8)
      {
        const double beta = 0.8;
        double c1 = ((I(1.0) + I(beta)) * sqr(RsI)).hi;
        double c2 = ((I(1.0) + recip(I(beta))) * sqr(I(cd.weighted))).hi;
        std::vector<I> Gi((size_t)k * k);
        for (int a = 0; a < k; ++a)
          for (int b = a; b < k; ++b) {
            I s1(0.0), s2(0.0);
            for (int l = 0; l < 9; ++l) s1 = s1 + mulp(R.T[a][SEEDPOS[l]], I(R.T[b][SEEDPOS[l]]));
            for (int f = 0; f < 16; ++f) s2 = s2 + mulp(R.T[a][dpos(f)], I(R.T[b][dpos(f)])) * I(cd.invmag[f]);
            Gi[(size_t)a * k + b] = mulp(c1, s1) + mulp(c2, s2);
          }
        std::vector<double> Gd((size_t)k * k);
        dominate(Gi, k, Gd.data());
        std::vector<std::array<I, 10>> corners(1);
        for (int a = 0; a < k; ++a) corners[0][a] = U0[a];
        double E = ell(Gd.data(), k, corners);
        I t = (I(E) + holterm) / TWOPI() + sq1ort * epsI;
        bestT = std::min(bestT, t.hi);
      }
      t0 = bestT;
      if (!(t0 < 10.0)) return false;
    }
    // bootstrap
    I ort1 = I(1.0) + I(R.ort);
    double cur = s;
    for (int it = 0; it < 30; ++it) {
      I n = I(R.n0) + I(R.id) * I(cur) + I(R.C) * sqr(I(cur));
      I nxt = sqrt0(sqr(n) + ort1 * sqr(I(t0)));
      if (nxt.hi < cur) cur = nxt.hi;
      else break;
    }
    sOut = cur;
    return std::isfinite(cur);
  } catch (const std::exception&) {
    return false;
  }
}

double approxDist(const double P0[36], const double P[36]) {
  double W[36], mr[6] = {0}, mc[6] = {0}, tm = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) {
      double d0 = P[6 * i + j] - P[6 * i] - P[j] + P[0];
      double N = std::floor(P0[6 * i + j] - d0 + 0.5);
      W[6 * i + j] = P0[6 * i + j] - d0 - N;
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

bool tryCover(const SeedCtx& S, const Params& P, CellWork& w, CoverCtx& cc) {
  using clk = std::chrono::steady_clock;
  auto el = [](clk::time_point a) { return std::chrono::duration<double>(clk::now() - a).count(); };
  auto t0 = clk::now();
  CoverData cd;
  if (!prepCover(S, P, w, cd)) { cc.nPrepFail++; cc.tPrep += el(t0); return false; }
  cc.tPrep += el(t0);
  const LevelK& lk = S.L(w.Nu);
  t0 = clk::now();
  std::vector<std::pair<double, int>> cand;
  for (size_t i = 0; i < cc.roots.size(); ++i) cand.push_back({approxDist(cd.P0, cc.roots[i].P), (int)i});
  std::sort(cand.begin(), cand.end());
  if ((int)cand.size() > P.topK) cand.resize(P.topK);
  auto wcheck = [&](const Root& Rt, double sb) {
    if (!cc.truth || !cc.cellHasWitness) return;
    double d = mergeDistanceUpper(cc.truth, Rt.P);
    cc.witChecks++;
    if (sb > 0) cc.witWorstRatio = std::max(cc.witWorstRatio, d / sb);
    if (d > sb + 1e-9) cc.witViol++;
  };
  for (auto& c : cand) {
    double s;
    cc.nBoundTry++;
    bool okb = coverBound(cc.roots[c.second], cd, S, lk, P.target, s);
    if (okb) wcheck(cc.roots[c.second], s);
    if (okb && s <= P.target) {
      if (s > cc.roots[c.second].r) cc.roots[c.second].r = s;
      cc.nReuse++;
      cc.tReuse += el(t0);
      return true;
    }
  }
  cc.tReuse += el(t0);
  t0 = clk::now();
  Root R;
  bool okb = buildRoot(cd.P0, R);
  cc.tBuild += el(t0);
  if (!okb) { cc.nBuildFail++; return false; }
  t0 = clk::now();
  double s;
  bool okc0 = coverBound(R, cd, S, lk, P.target, s);
  if (okc0) wcheck(R, s);
  bool okc = okc0 && (s <= P.target);
  cc.tSnapBound += el(t0);
  if (!okc) { cc.nSnapBound++; return false; }
  R.r = s;
  cc.roots.push_back(R);
  cc.nSnapOK++;
  return true;
}

}  // namespace sa
