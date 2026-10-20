// Emits JSON cases (hex floats = exact doubles) for tests/xcheck_bc.py, which re-derives every quantity with python-flint
// (arb/acb, >=200 bits) and checks CONTAINMENT of the flint ball inside the C++ interval (one-sided for the double upper
// bounds of Stage B).  usage: xgen N seed out.json
#include "../stageB.hpp"
#include "../stageC.hpp"
#include "fixtures.hpp"
#include <cstdio>
#include <string>

using namespace mb;
static std::mt19937_64 rng;
static double urand(double a, double b) { return std::uniform_real_distribution<double>(a, b)(rng); }

static std::string hx(double x) { char b[64]; std::snprintf(b, sizeof b, "\"%a\"", x); return b; }
static std::string hxi(const Interval& i) { return "[" + hx(i.lo) + "," + hx(i.hi) + "]"; }
static std::string hvec(const double* v, int n) {
  std::string s = "[";
  for (int i = 0; i < n; ++i) s += (i ? "," : "") + hx(v[i]);
  return s + "]";
}
static std::string recJson(const Record& R) {
  std::string s = "{\"k\":" + std::to_string(R.k) + ",\"p\":" + hvec(R.p, 36) + ",\"T\":" + hvec(R.T.data(), R.T.size());
  s += ",\"C\":" + hx(R.C) + ",\"n0\":" + hx(R.n0) + ",\"id\":" + hx(R.id) + ",\"ort\":" + hx(R.ort) + ",\"hor\":" + hx(R.hor);
  s += ",\"mnorm\":" + hx(R.mnorm) + ",\"froot\":" + hx(R.froot) + ",\"r\":" + hx(R.r) + "}";
  return s;
}
static Record fourierRoot(double a, double b, double r) {
  double P[36];
  fx::fourierPhases(a, b, P);
  snapSolution(P);
  Record R;
  certifyRoot(P, R);
  R.r = r;
  return R;
}

int main(int argc, char** argv) {
  if (argc < 4) return 2;
  int N = std::atoi(argv[1]);
  rng.seed(std::atoll(argv[2]));
  cn();
  iv::RoundUp g;
  FILE* f = std::fopen(argv[3], "w");
  std::fprintf(f, "[\n");
  bool first = true;
  auto emit = [&](const std::string& s) { std::fprintf(f, "%s%s\n", first ? "" : ",", s.c_str()); first = false; };
  for (int t = 0; t < N; ++t) {
    Record P = fourierRoot(urand(-1.5, 1.5), urand(-1.5, 1.5), urand(0.0, 0.03));
    CRec C;
    initCRec(C, P);
    // ---- single-vector bin case
    {
      int Ms[4] = {12, 24, 48, 96};
      int M = Ms[rng() % 4];
      Layer L;
      L.init(C, M);
      uint8_t b[5];
      for (int j = 0; j < 5; ++j) b[j] = (uint8_t)(rng() % M);
      Single s;
      loadBin(L, b, s);
      fillS(s);
      fillCJ(s);
      fillB(s, C);
      double y[6];
      for (int i = 0; i < 6; ++i) y[i] = urand(-1, 1);
      SupportOut so = binSupport(C, s, y, L.a, L.X, L.I);
      // direct-test aggregate and scalar-test values
      const Consts& cc = cn();
      Interval agg(0.0);
      for (int i = 0; i < 6; ++i) agg = agg + iv::sqr(iv::sqrt(s.a2[i]) - Interval(1.0));
      Interval thr = iv::sqr(cc.twopi * (cc.K * L.X + C.r));
      Interval tau = cc.twopi * (L.Lm * cc.isqrt6 + C.r);
      std::string line = "{\"type\":\"single\",\"M\":" + std::to_string(M) + ",\"bin\":[";
      for (int j = 0; j < 5; ++j) line += (j ? "," : "") + std::to_string((int)b[j]);
      line += "],\"rec\":" + recJson(P) + ",\"y\":" + hvec(y, 6);
      line += ",\"lhs\":" + hxi(so.lhs) + ",\"rhs\":" + hxi(so.rhs) + ",\"agg\":" + hxi(agg) + ",\"thr\":" + hxi(thr) + ",\"tau\":" + hxi(tau);
      line += ",\"absS2\":[";
      for (int i = 0; i < 6; ++i) line += (i ? "," : "") + hxi(s.a2[i]);
      line += "],\"c\":[";
      for (int i = 0; i < 6; ++i) line += (i ? "," : "") + hxi(s.c[i]);
      line += "]}";
      emit(line);
    }
    // ---- edge + joint case
    {
      double v[2][6], sr[2], y[14];
      std::unique_ptr<Ball> bl[2];
      for (int q = 0; q < 2; ++q) {
        for (int j = 0; j < 6; ++j) v[q][j] = urand(-1, 2);
        sr[q] = urand(0.005, 0.06);
        bl[q].reset(new Ball());
        makeBall(C, v[q], sr[q], *bl[q]);
      }
      for (int i = 0; i < 14; ++i) y[i] = urand(-1, 1);
      EdgeData e;
      edgeData(*bl[0], *bl[1], e);
      const Consts& cc = cn();
      SupportOut so[2];
      for (int rel = 0; rel < 2; ++rel) so[rel] = jointSupport(C, *bl[0], *bl[1], e, rel == 0, y);
      Interval EO = cc.pi * iv::sqr(e.re);
      Interval oLhs = e.a2 / cc.twopi, oRhs = Interval(6.0) * e.re * e.nDU + EO * e.absS;
      Interval qU = (e.a2 * cc.sixth - Interval(1.0)) / cc.c4pi;
      Interval EU = cc.pi * iv::sqr(e.re) * (Interval(1.0) + e.absS * cc.sixth);
      Interval uLhs = iv::abs(qU), uRhs = e.re * e.nDU + EU;
      std::string line = "{\"type\":\"joint\",\"rec\":" + recJson(P) + ",\"v0\":" + hvec(v[0], 6) + ",\"v1\":" + hvec(v[1], 6) + ",\"s\":" + hvec(sr, 2) + ",\"y\":" + hvec(y, 14);
      line += ",\"absS\":" + hxi(e.absS) + ",\"jointO\":{\"lhs\":" + hxi(so[0].lhs) + ",\"rhs\":" + hxi(so[0].rhs) + "},\"jointU\":{\"lhs\":" + hxi(so[1].lhs) + ",\"rhs\":" + hxi(so[1].rhs) + "}";
      line += ",\"orthSupport\":{\"lhs\":" + hxi(oLhs) + ",\"rhs\":" + hxi(oRhs) + "},\"unbSupport\":{\"lhs\":" + hxi(uLhs) + ",\"rhs\":" + hxi(uRhs) + "}}";
      emit(line);
    }
    // ---- Stage B relation (one-sided upper bounds)
    if (t % 4 == 0) {
      double a0 = urand(-1, 1), b0 = urand(-1, 1);
      Record Pp = fourierRoot(a0, b0, urand(0.03, 0.125));
      Record Q = fourierRoot(a0 + urand(-0.06, 0.06), b0 + urand(-0.06, 0.06), 0.0);
      Rel rel = makeRel(Pp, Q, 0.249);
      if (rel.ok) {
        std::string line = "{\"type\":\"rel\",\"P\":" + recJson(Pp) + ",\"Q\":" + recJson(Q) + ",\"dq\":" + hx(rel.dqNorm) + ",\"D\":" + hx(rel.D) + ",\"v0\":" + hx(rel.v0) + ",\"lam\":" + hx(rel.lam);
        line += ",\"e0\":" + hx(rel.e0) + ",\"sigma\":" + hx(rel.sigma) + ",\"beta\":" + hx(rel.beta) + ",\"etaq\":[";
        for (int a = 0; a < Pp.k; ++a) line += (a ? "," : "") + hxi(rel.etaq[a]);
        line += "],\"evals\":[";
        double ts[4] = {0.0, 0.01, 0.03, 0.06};
        for (int i = 0; i < 4; ++i) {
          double val = rel.eval(ts[i]);
          line += (i ? "," : "") + std::string("[") + hx(ts[i]) + "," + hx(val) + "," + hx((double)rel.lastSteps) + "]";
        }
        line += "]}";
        emit(line);
      }
    }
  }
  std::fprintf(f, "]\n");
  std::fclose(f);
  return 0;
}
