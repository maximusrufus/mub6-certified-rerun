// seedproc.cpp -- Stage A seed blocks, frames, U(2) cells, norm/directional tests.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include "cover.hpp"
#include "la.hpp"
#include "stagea.hpp"

namespace sa {

thread_local double g_dirCertY[20];

// ------------------------------------------------------------------ tables
struct LevelTab {
  int Nu = 0;
  std::vector<Z> eA, eB, eG, eT;  // e(alpha), e(beta), e(gamma), cos+i sin of 2 pi t
};
static LevelTab g_tabs[6];
static bool g_tabsReady = false;

static int levelIndex(int Nu) {
  int i = 0;
  while ((16 << i) < Nu) ++i;
  return i;
}

void initTables() {
  if (g_tabsReady) return;
  for (int li = 0; li < 6; ++li) {
    LevelTab& T = g_tabs[li];
    int Nu = 16 << li;
    T.Nu = Nu;
    auto mid = [&](int i) { return (2.0 * i + 1.0) / (2.0 * Nu); };  // exact dyadic
    for (int i = 0; i < Nu / 2; ++i) T.eA.push_back(iv::cis_turns(I(mid(i))));
    for (int i = 0; i < Nu; ++i) T.eB.push_back(iv::cis_turns(I(mid(i))));
    for (int i = 0; i < Nu; ++i) T.eG.push_back(iv::cis_turns(I(mid(i))));
    for (int i = 0; i < Nu / 4; ++i) T.eT.push_back(iv::cis_turns(I(mid(i))));
  }
  g_tabsReady = true;
}

const LevelK& SeedCtx::L(int Nu) const { return lev[levelIndex(Nu)]; }

// ------------------------------------------------------------------ frames
// Orthonormal 4x2 K with range ker M^*, M = (1 p), by Gram-Schmidt in intervals.
// Selection of coordinate vectors by the largest (midpoint) residual norm (A7).
static void frameK(const Z p[4], CM<4, 2>& K, I& Sabs, int fc[2]) {
  Z S = p[0] + p[1] + p[2] + p[3];
  Sabs = zabs(S);
  Z q1[4], q2[4];
  for (int k = 0; k < 4; ++k) q1[k] = Z(I(0.5), I(0.0));
  Z c1 = zmulI(S, I(0.5));  // <q1,p> = S/2
  Z r2[4];
  I n2s(0.0);
  for (int k = 0; k < 4; ++k) {
    r2[k] = p[k] - q1[k] * c1;
    n2s = n2s + abs2(r2[k]);
  }
  I n2 = sqrt0(n2s);
  if (!(n2.lo > 0.449999)) throw Fail("frame: GS residual 2 small");
  I in2 = recip(n2);
  for (int k = 0; k < 4; ++k) q2[k] = zmulI(r2[k], in2);
  // third vector
  int bk = -1;
  double bn = -1;
  Z res[4][4];
  I nrm[4];
  for (int k = 0; k < 4; ++k) {
    // e_k - q1 <q1,e_k> - q2 <q2,e_k>; <q1,e_k>=1/2 ; <q2,e_k>=conj(q2_k)
    I s(0.0);
    for (int j = 0; j < 4; ++j) {
      Z v(j == k ? 1.0 : 0.0, 0.0);
      v = v - q1[j] * Z(I(0.5), I(0.0)) - q2[j] * iv::conj(q2[k]);
      res[k][j] = v;
      s = s + abs2(v);
    }
    nrm[k] = sqrt0(s);
    if (midp(nrm[k]) > bn) {
      bn = midp(nrm[k]);
      bk = k;
    }
  }
  if (!(nrm[bk].lo > 0.449999)) throw Fail("frame: GS residual 3 small");
  fc[0] = bk;
  I in3 = recip(nrm[bk]);
  Z q3[4];
  for (int j = 0; j < 4; ++j) q3[j] = zmulI(res[bk][j], in3);
  // fourth vector
  int bl = -1;
  double bn2 = -1;
  Z res2[4][4];
  I nrm2[4];
  for (int l = 0; l < 4; ++l) {
    if (l == bk) continue;
    I s(0.0);
    for (int j = 0; j < 4; ++j) {
      Z v(j == l ? 1.0 : 0.0, 0.0);
      v = v - q1[j] * Z(I(0.5), I(0.0)) - q2[j] * iv::conj(q2[l]) - q3[j] * iv::conj(q3[l]);
      res2[l][j] = v;
      s = s + abs2(v);
    }
    nrm2[l] = sqrt0(s);
    if (midp(nrm2[l]) > bn2) {
      bn2 = midp(nrm2[l]);
      bl = l;
    }
  }
  if (!(nrm2[bl].lo > 0.449999)) throw Fail("frame: GS residual 4 small");
  fc[1] = bl;
  I in4 = recip(nrm2[bl]);
  for (int j = 0; j < 4; ++j) {
    K.a[j][0] = q3[j];
    K.a[j][1] = zmulI(res2[bl][j], in4);
  }
}

bool buildSeed(SeedCtx& S, const double q[9], const Params& P, std::string& why) {
  try {
    for (int i = 0; i < 9; ++i) S.q[i] = q[i];
    S.a = iv::cis_turns(I(q[0]));
    for (int j = 0; j < 4; ++j) S.b[j] = iv::cis_turns(I(q[1 + j]));
    for (int i = 0; i < 4; ++i) S.c[i] = iv::cis_turns(I(q[5 + i]));
    const Z one(1.0, 0.0), zero(0.0, 0.0);
    S.E0.a[0][0] = one; S.E0.a[0][1] = one; S.E0.a[1][0] = one; S.E0.a[1][1] = S.a;
    for (int j = 0; j < 4; ++j) { S.B0.a[0][j] = one; S.B0.a[1][j] = S.b[j]; }
    for (int i = 0; i < 4; ++i) { S.C0.a[i][0] = one; S.C0.a[i][1] = S.c[i]; }
    // Rs: upper bound of 2 pi rs, rs >= .0320001
    double rsU = std::nextafter(P.rsCycles, 1.0);
    S.Rs = (TWOPI() * I(rsU)).hi;
    I RsI(S.Rs);
    // frames
    Z pb[4], pc[4];
    for (int j = 0; j < 4; ++j) { pb[j] = iv::conj(S.b[j]); pc[j] = S.c[j]; }
    I SbAbs, ScAbs;
    frameK(pb, S.KB, SbAbs, S.fc);
    frameK(pc, S.KC, ScAbs, S.fc + 2);
    S.mB = sqrt0(I(4.0) - SbAbs);
    S.mC = sqrt0(I(4.0) - ScAbs);
    S.m = imin(S.mB, S.mC);
    if (!(S.m.lo > P.mFrameMin)) throw Fail("frame: m too small");
    S.Cop = sqrt0(I(4.0) + ScAbs);
    // G_B = B0^* (B0 B0^*)^{-1}, G_C = C0 (C0^* C0)^{-1}
    CM<2, 4> B0 = S.B0;
    CM<4, 2> B0s = adj(B0);
    CM<2, 2> BBs = mm(B0, B0s);
    S.GB = mm(B0s, inv2(BBs));
    CM<2, 4> C0s = adj(S.C0);
    CM<2, 2> CsC = mm(C0s, S.C0);
    S.GC = mm(S.C0, inv2(CsC));
    CM<2, 4> KBs = adj(S.KB);
    S.KBs = KBs;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 2; ++j) S.KC6.a[i][j] = Z(S.KC.a[i][j].re * SQRT6(), S.KC.a[i][j].im * SQRT6());
    S.PK = mm(S.KB, KBs);
    CM<2, 4> GBs = adj(S.GB);
    CM<2, 2> E0s = adj(S.E0);
    S.Dspan = mneg(mm(mm(S.C0, E0s), GBs));
    // constants A0, A*
    S.A0 = sqrt0(I(2.0) + zabs(S.E0.a[0][1] + S.a));  // |1+a|
    {
      I t = S.A0 + RsI;
      S.Astar = I(std::min(2.0, t.lo), std::min(2.0, t.hi));
    }
    // seed variation matrices
    for (int l = 0; l < 9; ++l) {
      CM<2, 2> E1 = mzero<2, 2>();
      CM<2, 4> B1 = mzero<2, 4>();
      CM<4, 2> C1 = mzero<4, 2>();
      if (l == 0) E1.a[1][1] = zi_mul(S.a);
      else if (l <= 4) B1.a[1][l - 1] = zi_mul(S.b[l - 1]);
      else C1.a[l - 5][1] = zi_mul(S.c[l - 5]);
      CM<2, 2> E1s = adj(E1);
      CM<4, 2> B1s = adj(B1);
      CM<2, 4> C1s = adj(C1);
      // Cst = -C1 E0* GB* - C0 E1* GB* - GC E0* B1 PK
      CM<4, 4> t1 = mm(mm(C1, E0s), GBs);
      CM<4, 4> t2 = mm(mm(S.C0, E1s), GBs);
      CM<4, 4> t3 = mm(mm(mm(S.GC, E0s), B1), S.PK);
      S.Cst[l] = mneg(madd(madd(t1, t2), t3));
      S.Lft[l] = mm(B1s, GBs);
      S.Lc[l] = mm(S.GC, C1s);
    }
    // sums: c_X = (1 + a + sum b)/i, c_Y
    {
      Z sx = one + S.a, sy = one + S.a;
      for (int j = 0; j < 4; ++j) { sx = sx + S.b[j]; sy = sy + S.c[j]; }
      S.cX[0] = sx.im; S.cX[1] = -sx.re;
      S.cY[0] = sy.im; S.cY[1] = -sy.re;
    }
    // level constants
    S.lev.clear();
    for (int li = 0; li < 6; ++li) {
      LevelK lk;
      lk.Nu = 16 << li;
      I h = PI() * recip(I((double)lk.Nu));
      lk.h = h.hi;
      I Rs2 = sqr(RsI);
      I mC2 = sqr(S.mC);
      I qq = Rs2 * (I(1.0) + (sqr(S.Astar) + I(6.0)) / mC2);
      if (!(qq.hi < 0.32)) throw Fail("scale: q >= .32");
      lk.q = qq.hi;
      I pp = I(qq.hi) / (SQRT6() + sqrt0(I(6.0) - I(qq.hi)));
      if (!(pp.hi < 0.068)) throw Fail("scale: pol >= .068");
      lk.p = pp.hi;
      I off = RsI * (SQRT6() + S.Astar) / S.m;
      if (li == 0 && !(off.hi < 0.74)) throw Fail("scale: off >= .74");
      I z2 = sqr(off) + sqr(mulp(6.0, h) + I(lk.p));
      lk.z2 = z2.hi;
      if (li == 0 && !(sqrt0(z2).hi < 1.5)) throw Fail("scale: z >= 1.5");
      I zz = sqrt0(I(lk.z2));
      I L = imax(SQRT6(), S.Cop * sqrt0(I(1.0) + sqr(S.A0 / S.mB)));
      I eoff = Rs2 * (L + S.A0 + I(1.0)) / (S.m * I(2.0)) + zz * RsI / S.m;
      I e3 = SQRT252() * I(0.5) * sqr(h) + I(lk.p);
      I Er = sqrt0(sqr(eoff) + sqr(e3));
      if (li == 0 && !(Er.hi < 1.0)) throw Fail("scale: E >= 1");
      lk.Er = Er.hi;
      S.lev.push_back(lk);
    }
    return true;
  } catch (const std::exception& e) {
    why = e.what();
    return false;
  }
}

// ------------------------------------------------------------------- cells
static void buildU(const LevelTab& T, const Cell& c, Z U[2][2]) {
  const Z& eA = T.eA[c.ia];
  const Z& eB = T.eB[c.ib];
  const Z& eG = T.eG[c.ig];
  const I& cT = T.eT[c.it].re;
  const I& sT = T.eT[c.it].im;
  Z ab = eA * eB, ag = eA * eG;
  Z abc = eA * iv::conj(eB), agc = eA * iv::conj(eG);
  U[0][0] = zmulI(ab, cT);
  U[0][1] = zmulI(ag, sT);
  U[1][0] = -zmulI(agc, sT);
  U[1][1] = zmulI(abc, cT);
}

void cellD(const SeedCtx& S, int Nu, const Cell& c, CellWork& w) {
  const LevelTab& T = g_tabs[levelIndex(Nu)];
  w.c = c;
  w.Nu = Nu;
  buildU(T, c, w.U);
  CM<2, 2> U;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) U.a[i][j] = w.U[i][j];
  CM<4, 2> W = mm(S.KC6, U);
  w.D = madd(S.Dspan, mm(W, S.KBs));
  w.haveJ = false;
}

bool normExcluded(const SeedCtx&, const LevelK& lk, const CellWork& w) {
  I s(0.0);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      I m = zabs(w.D.a[i][j]);
      s = s + sqr(m - I(1.0));
    }
  return s.lo > lk.z2;
}

void ensureJ(const SeedCtx& S, CellWork& w) {
  if (w.haveJ) return;
  const LevelTab& T = g_tabs[levelIndex(w.Nu)];
  const Cell& c = w.c;
  CM<4, 4> DPK = mm(w.D, S.PK);
  for (int l = 0; l < 9; ++l) {
    CM<4, 4> t = msub(S.Cst[l], mm(w.D, S.Lft[l]));
    w.Jc[l] = msub(t, mm(S.Lc[l], DPK));
  }
  // angular derivatives of U
  const Z& eA = T.eA[c.ia];
  const Z& eB = T.eB[c.ib];
  const Z& eG = T.eG[c.ig];
  const I& cT = T.eT[c.it].re;
  const I& sT = T.eT[c.it].im;
  CM<2, 2> dU[4];
  for (int a = 0; a < 4; ++a) dU[a] = mzero<2, 2>();
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) dU[0].a[i][j] = zi_mul(w.U[i][j]);  // alpha
  dU[1].a[0][0] = -zmulI(eA * eB, sT);                                // t
  dU[1].a[0][1] = zmulI(eA * eG, cT);
  dU[1].a[1][0] = -zmulI(eA * iv::conj(eG), cT);
  dU[1].a[1][1] = -zmulI(eA * iv::conj(eB), sT);
  dU[2].a[0][0] = zi_mul(w.U[0][0]);                                  // beta
  dU[2].a[1][1] = -zi_mul(w.U[1][1]);
  dU[3].a[0][1] = zi_mul(w.U[0][1]);                                  // gamma
  dU[3].a[1][0] = -zi_mul(w.U[1][0]);
  // theta ordering of J columns: (alpha, t, beta, gamma) = columns 9..12
  for (int a = 0; a < 4; ++a) w.Jc[9 + a] = mm(mm(S.KC6, dU[a]), S.KBs);
  w.haveJ = true;
}

// ----------------------------------------------------- directional machinery
namespace {

struct DirData {
  I J[20][13];
  I c[20];
  I D2[16];
  double Jm[20][13], cm[20], D2m[16];
};

// V(y) in doubles (heuristic); returns V and fills grad (subgradient of V).
double approxV(const double y[20], const DirData& d, double Rs, double h, double z2, double Er, double grad[20]) {
  double g[13];
  for (int l = 0; l < 13; ++l) {
    double s = 0;
    for (int i = 0; i < 20; ++i) s += y[i] * d.Jm[i][l];
    g[l] = s;
  }
  double gs2 = 0;
  for (int l = 0; l < 9; ++l) gs2 += g[l] * g[l];
  double gsn = std::sqrt(gs2);
  double wl[13];
  for (int l = 0; l < 9; ++l) wl[l] = gsn > 0 ? Rs * g[l] / gsn : 0.0;
  for (int l = 9; l < 13; ++l) wl[l] = h * (g[l] > 0 ? 1.0 : (g[l] < 0 ? -1.0 : 0.0));
  double rhs = Rs * gsn;
  for (int l = 9; l < 13; ++l) rhs += h * std::fabs(g[l]);
  // E_r ||D o y_D||
  double dd = 0;
  for (int f = 0; f < 16; ++f) dd += d.D2m[f] * y[4 + f] * y[4 + f];
  double ddn = std::sqrt(dd);
  rhs += Er * ddn;
  double mx = 0;
  int imx = -1;
  for (int f = 0; f < 16; ++f)
    if (-y[4 + f] > mx) { mx = -y[4 + f]; imx = f; }
  rhs += 0.5 * z2 * mx;
  double yx = std::hypot(y[0], y[1]), yy = std::hypot(y[2], y[3]);
  rhs += 0.5 * Rs * Rs * (yx + yy);
  double yc = 0;
  for (int i = 0; i < 20; ++i) yc += y[i] * d.cm[i];
  // gradient
  for (int i = 0; i < 20; ++i) {
    double s = 0;
    for (int l = 0; l < 9; ++l) s += wl[l] * d.Jm[i][l];
    for (int l = 9; l < 13; ++l) s += wl[l] * d.Jm[i][l];
    grad[i] = d.cm[i] - s;
  }
  if (ddn > 0)
    for (int f = 0; f < 16; ++f) grad[4 + f] -= Er * d.D2m[f] * y[4 + f] / ddn;
  if (imx >= 0) grad[4 + imx] += 0.5 * z2;
  if (yx > 0) { grad[0] -= 0.5 * Rs * Rs * y[0] / yx; grad[1] -= 0.5 * Rs * Rs * y[1] / yx; }
  if (yy > 0) { grad[2] -= 0.5 * Rs * Rs * y[2] / yy; grad[3] -= 0.5 * Rs * Rs * y[3] / yy; }
  return yc - rhs;
}

// rigorous check of the directional inequality for exact doubles y.
bool certifyDir(const double y[20], const DirData& d, double Rs, double h, double z2, double Er) {
  I lhs(0.0);
  for (int i = 0; i < 20; ++i) lhs = lhs + mulp(y[i], d.c[i]);
  I g[13];
  for (int l = 0; l < 13; ++l) {
    I s(0.0);
    for (int i = 0; i < 20; ++i) s = s + mulp(y[i], d.J[i][l]);
    g[l] = s;
  }
  I gs(0.0);
  for (int l = 0; l < 9; ++l) gs = gs + sqr(g[l]);
  I rhs = I(Rs) * sqrt0(gs);
  I gt(0.0);
  for (int l = 9; l < 13; ++l) gt = gt + iv::abs(g[l]);
  rhs = rhs + I(h) * gt;
  I dd(0.0);
  for (int f = 0; f < 16; ++f) dd = dd + d.D2[f] * sqr(I(y[4 + f]));
  rhs = rhs + I(Er) * sqrt0(dd);
  double mx = 0;
  for (int f = 0; f < 16; ++f) mx = std::max(mx, -y[4 + f]);
  rhs = rhs + I(0.5) * I(z2) * I(mx);
  I yx = sqrt0(sqr(I(y[0])) + sqr(I(y[1]))), yy = sqrt0(sqr(I(y[2])) + sqr(I(y[3])));
  rhs = rhs + I(0.5) * sqr(I(Rs)) * (yx + yy);
  return lhs.lo > rhs.hi;
}

}  // namespace

static void buildDirData(const SeedCtx& S, const CellWork& w, DirData& d) {
  for (int i = 0; i < 20; ++i)
    for (int l = 0; l < 13; ++l) d.J[i][l] = I(0.0);
  d.J[0][0] = S.a.re; d.J[1][0] = S.a.im;
  d.J[2][0] = S.a.re; d.J[3][0] = S.a.im;
  for (int j = 0; j < 4; ++j) {
    d.J[0][1 + j] = S.b[j].re; d.J[1][1 + j] = S.b[j].im;
    d.J[2][5 + j] = S.c[j].re; d.J[3][5 + j] = S.c[j].im;
  }
  d.c[0] = S.cX[0]; d.c[1] = S.cX[1]; d.c[2] = S.cY[0]; d.c[3] = S.cY[1];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      int f = 4 * i + j;
      const Z& Dij = w.D.a[i][j];
      d.D2[f] = abs2(Dij);
      d.c[4 + f] = mulp(0.5, d.D2[f] - I(1.0));
      for (int l = 0; l < 13; ++l) {
        const Z& Jl = w.Jc[l].a[i][j];
        d.J[4 + f][l] = Dij.re * Jl.re + Dij.im * Jl.im;
      }
    }
  for (int i = 0; i < 20; ++i) {
    d.cm[i] = midp(d.c[i]);
    for (int l = 0; l < 13; ++l) d.Jm[i][l] = midp(d.J[i][l]);
  }
  for (int f = 0; f < 16; ++f) d.D2m[f] = midp(d.D2[f]);
}

bool dirExcluded(const SeedCtx& S, const LevelK& lk, CellWork& w) {
  ensureJ(S, w);
  DirData d;
  buildDirData(S, w, d);
  const double Rs = S.Rs, h = lk.h, z2 = lk.z2, Er = lk.Er;
  // candidate directions: c - J (J^T J + lam I)^{-1} J^T c
  double N[13 * 13], rhs[13];
  for (int a = 0; a < 13; ++a) {
    for (int b = 0; b < 13; ++b) {
      double s = 0;
      for (int i = 0; i < 20; ++i) s += d.Jm[i][a] * d.Jm[i][b];
      N[a * 13 + b] = s;
    }
    double s = 0;
    for (int i = 0; i < 20; ++i) s += d.Jm[i][a] * d.cm[i];
    rhs[a] = s;
  }
  double best[20];
  double bestV = -1e300;
  const double lams[] = {1e-9, 1e-6, 1e-4, 1e-2, 1e-1};
  double grad[20];
  auto consider = [&](double y[20]) {
    double nn = 0;
    for (int i = 0; i < 20; ++i) nn += y[i] * y[i];
    nn = std::sqrt(nn);
    if (!(nn > 1e-300)) return;
    double yn[20];
    for (int i = 0; i < 20; ++i) yn[i] = y[i] / nn;
    double V = approxV(yn, d, Rs, h, z2, Er, grad);
    if (V > bestV) { bestV = V; for (int i = 0; i < 20; ++i) best[i] = yn[i]; }
  };
  for (double lam : lams) {
    double wv[13];
    if (!la::chol_solve(N, 13, lam, rhs, wv)) continue;
    double y[20];
    for (int i = 0; i < 20; ++i) {
      double s = d.cm[i];
      for (int l = 0; l < 13; ++l) s -= d.Jm[i][l] * wv[l];
      y[i] = s;
    }
    consider(y);
  }
  {
    double y[20];
    for (int i = 0; i < 20; ++i) y[i] = d.cm[i];
    consider(y);
  }
  if (!(bestV > -1e299)) return false;
  // projected subgradient ascent on the sphere
  double cur[20];
  for (int i = 0; i < 20; ++i) cur[i] = best[i];
  static const int DIRIT = getenv("DIRITERS") ? atoi(getenv("DIRITERS")) : 30;
  for (int it = 0; it < DIRIT; ++it) {
    double V = approxV(cur, d, Rs, h, z2, Er, grad);
    if (V > bestV) { bestV = V; for (int i = 0; i < 20; ++i) best[i] = cur[i]; }
    double gy = 0;
    for (int i = 0; i < 20; ++i) gy += grad[i] * cur[i];
    double gn = 0;
    for (int i = 0; i < 20; ++i) { grad[i] -= gy * cur[i]; gn += grad[i] * grad[i]; }
    gn = std::sqrt(gn);
    if (!(gn > 1e-14)) break;
    double eta = 0.5 / (1.0 + it * 0.5);
    double nn = 0;
    for (int i = 0; i < 20; ++i) { cur[i] += eta * grad[i] / gn * 0.5; nn += cur[i] * cur[i]; }
    nn = std::sqrt(nn);
    for (int i = 0; i < 20; ++i) cur[i] /= nn;
  }
  if (!(bestV > 0)) return false;
  // finalize: exact doubles, divide by max abs
  double mxa = 0;
  for (int i = 0; i < 20; ++i) mxa = std::max(mxa, std::fabs(best[i]));
  if (!(mxa > 1e-100 && mxa < 1e100)) return false;
  double yf[20];
  for (int i = 0; i < 20; ++i) yf[i] = best[i] / mxa;
  bool ok = certifyDir(yf, d, Rs, h, z2, Er);
  if (ok) for (int i = 0; i < 20; ++i) g_dirCertY[i] = yf[i];
  return ok;
}

bool seedSumExcluded(const SeedCtx& S) {
  // sum equations only: y = (yX, yY, 0)
  I RsI(S.Rs);
  for (int mode = 0; mode < 3; ++mode) {
    double yx[2] = {0, 0}, yy[2] = {0, 0};
    if (mode != 1) { yx[0] = midp(S.cX[0]); yx[1] = midp(S.cX[1]); }
    if (mode != 0) { yy[0] = midp(S.cY[0]); yy[1] = midp(S.cY[1]); }
    double nx = std::hypot(yx[0], yx[1]), ny = std::hypot(yy[0], yy[1]);
    if (nx > 0) { yx[0] /= nx; yx[1] /= nx; }
    if (ny > 0) { yy[0] /= ny; yy[1] /= ny; }
    I lhs = mulp(yx[0], S.cX[0]) + mulp(yx[1], S.cX[1]) + mulp(yy[0], S.cY[0]) + mulp(yy[1], S.cY[1]);
    I gs(0.0);
    {
      I g0 = mulp(yx[0], S.a.re) + mulp(yx[1], S.a.im) + mulp(yy[0], S.a.re) + mulp(yy[1], S.a.im);
      gs = sqr(g0);
      for (int j = 0; j < 4; ++j) {
        I gb = mulp(yx[0], S.b[j].re) + mulp(yx[1], S.b[j].im);
        I gc = mulp(yy[0], S.c[j].re) + mulp(yy[1], S.c[j].im);
        gs = gs + sqr(gb) + sqr(gc);
      }
    }
    I yxn = sqrt0(sqr(I(yx[0])) + sqr(I(yx[1]))), yyn = sqrt0(sqr(I(yy[0])) + sqr(I(yy[1])));
    I rhs = RsI * sqrt0(gs) + I(0.5) * sqr(RsI) * (yxn + yyn);
    if (lhs.lo > rhs.hi) return true;
  }
  return false;
}

// ------------------------------------------------------------ seed driver
namespace {
bool witnessContains(const Witness& w, int Nu, const Cell& c) {
  const double tol = 1e-9;
  auto in01 = [&](double x, int idx, int n) {  // mod 1 closed bin [idx/n,(idx+1)/n]
    double lo = (double)idx / n, hi = (double)(idx + 1) / n;
    for (int s = -1; s <= 1; ++s) {
      double xx = x + s;
      if (xx >= lo - tol && xx <= hi + tol) return true;
    }
    return false;
  };
  auto inr = [&](double x, int idx, int n) {
    double lo = (double)idx / n, hi = (double)(idx + 1) / n;
    return x >= lo - tol && x <= hi + tol;
  };
  return inr(w.al, c.ia, Nu) && inr(w.t, c.it, Nu) && in01(w.be, c.ib, Nu) && in01(w.ga, c.ig, Nu);
}
}  // namespace

static uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// S6: deterministic seeded sample (hash of the seed bits, level, cell, kind) of excluded cells
static void pushSample(SeedResult& R, const Params& P, const SeedCtx& S, const LevelK& lk, int Nu, const Cell& c,
                       const double q[9], int kind, const double* y) {
  uint64_t key = P.sampleSeed;
  for (int i = 0; i < 9; ++i) { uint64_t b; std::memcpy(&b, &q[i], 8); key = splitmix64(key ^ b); }
  key = splitmix64(key ^ ((uint64_t)Nu << 40) ^ ((uint64_t)c.ia << 30) ^ ((uint64_t)c.it << 20) ^ ((uint64_t)c.ib << 10) ^ (uint64_t)c.ig);
  key = splitmix64(key ^ (uint64_t)kind);
  if (key % P.sampleMod != 0) return;
  ExclSample e;
  e.kind = kind; e.Nu = Nu; e.ia = c.ia; e.it = c.it; e.ib = c.ib; e.ig = c.ig;
  for (int i = 0; i < 4; ++i) e.fc[i] = S.fc[i];
  e.h = lk.h; e.z2 = lk.z2; e.Er = lk.Er; e.Rs = S.Rs;
  if (y) for (int i = 0; i < 20; ++i) e.y[i] = y[i];
  R.samples.push_back(e);
}

SeedResult processSeed(const double q[9], const Params& P, Witness* wit) {
  SeedResult R;
  SeedCtx S;
  if (!buildSeed(S, q, P, R.why)) {
    R.status = 1;
    return R;
  }
  try {
    if (seedSumExcluded(S)) {
      R.seedExcluded = true;
      if (wit && wit->on) wit->violations++;  // a witness inside an excluded seed ball is a violation
      return R;
    }
  } catch (const std::exception& e) {
    R.status = 1;
    R.why = e.what();
    return R;
  }
  struct Node {
    Cell c;
    int Nu;
  };
  std::vector<Node> stack;
  for (int ia = 0; ia < 8; ++ia)
    for (int it = 0; it < 4; ++it)
      for (int ib = 0; ib < 16; ++ib)
        for (int ig = 0; ig < 16; ++ig) stack.push_back({{(uint16_t)ia, (uint16_t)it, (uint16_t)ib, (uint16_t)ig}, 16});
  long extra = 0;
  bool giveUp = false;
  CoverCtx cc;
  if (wit && wit->on && wit->truthP) cc.truth = wit->truthP;
  while (!stack.empty()) {
    Node nd = stack.back();
    stack.pop_back();
    R.st.maxNuReached = std::max(R.st.maxNuReached, nd.Nu);
    if (nd.Nu >= 128) {
      if (++extra > P.extraCap) {
        giveUp = true;
        R.unresolved++;
        stack.clear();
        R.why = "extra node cap";
        break;
      }
    }
    try {
      CellWork w;
      auto tA = std::chrono::steady_clock::now();
      cellD(S, nd.Nu, nd.c, w);
      const LevelK& lk = S.L(nd.Nu);
      bool wc = wit && wit->on && witnessContains(*wit, nd.Nu, nd.c);
      R.st.cellsNorm++;
      bool nx = normExcluded(S, lk, w);
      auto tB = std::chrono::steady_clock::now();
      R.st.tNorm += std::chrono::duration<double>(tB - tA).count();
      if (nx) {
        R.st.exclNorm++;
        if (wc) wit->violations++;
        if (P.sampleMod) pushSample(R, P, S, lk, nd.Nu, nd.c, q, 0, nullptr);
        continue;
      }
      if (nd.Nu >= P.firstDirNu) {
        R.st.cellsDir++;
        bool dx = dirExcluded(S, lk, w);
        auto tC = std::chrono::steady_clock::now();
        R.st.tDir += std::chrono::duration<double>(tC - tB).count();
        tB = tC;
        if (dx) {
          R.st.exclDir++;
          if (wc) wit->violations++;
          if (P.sampleMod) pushSample(R, P, S, lk, nd.Nu, nd.c, q, 1, g_dirCertY);
          continue;
        }
      }
      if (nd.Nu >= P.firstCoverNu) {
        R.st.coverTry++;
        cc.cellHasWitness = wc;
        bool cv = tryCover(S, P, w, cc);
        R.st.tCover += std::chrono::duration<double>(std::chrono::steady_clock::now() - tB).count();
        if (cv) {
          R.st.covered++;
          continue;
        }
      }
      if (nd.Nu >= P.maxNu) {
        R.unresolved++;
        continue;
      }
      R.st.split++;
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          for (int c = 0; c < 2; ++c)
            for (int d = 0; d < 2; ++d)
              stack.push_back({{(uint16_t)(2 * nd.c.ia + a), (uint16_t)(2 * nd.c.it + b), (uint16_t)(2 * nd.c.ib + c),
                                (uint16_t)(2 * nd.c.ig + d)},
                               nd.Nu * 2});
    } catch (const std::exception& e) {
      // numerical breakdown: never an exclusion; refine or leave unresolved
      if (nd.Nu >= P.maxNu) {
        R.unresolved++;
      } else {
        for (int a = 0; a < 2; ++a)
          for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c)
              for (int d = 0; d < 2; ++d)
                stack.push_back({{(uint16_t)(2 * nd.c.ia + a), (uint16_t)(2 * nd.c.it + b), (uint16_t)(2 * nd.c.ib + c),
                                  (uint16_t)(2 * nd.c.ig + d)},
                                 nd.Nu * 2});
      }
    }
  }
  (void)giveUp;
  R.recs = cc.roots;
  if (wit) { wit->coverChecks += cc.witChecks; wit->coverViolations += cc.witViol; wit->worstRatio = std::max(wit->worstRatio, cc.witWorstRatio); }
  R.st.nPrepFail = cc.nPrepFail; R.st.nReuse = cc.nReuse; R.st.nSnapOK = cc.nSnapOK; R.st.nBuildFail = cc.nBuildFail;
  R.st.nSnapBound = cc.nSnapBound; R.st.nBoundTry = cc.nBoundTry;
  R.st.tPrep = cc.tPrep; R.st.tReuse = cc.tReuse; R.st.tBuild = cc.tBuild; R.st.tSnapBound = cc.tSnapBound;
  R.status = R.unresolved ? 2 : 0;
  return R;
}

}  // namespace sa
