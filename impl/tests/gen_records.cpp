// Generates JSON records (one per line, hex-float strings = exact) of random inputs and the C++ interval
// outputs, for harness/flint_check.py to re-evaluate with python-flint acb at 128+ bits.
// usage: gen_records N seed out.json
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include "../core/mub_ops.hpp"
using namespace iv;
using namespace mub;

static std::mt19937_64 rng;
static double U01() { return std::uniform_real_distribution<double>(0.0, 1.0)(rng); }

struct Pt { Interval I; double x; };  // interval and an exact double sample inside it

static Pt rand_pt(bool positive = false, bool nonzero = false) {
  // magnitude spread over a few decades; 40% exact points, else a random-width interval
  double scale = std::pow(10.0, std::uniform_real_distribution<double>(-3.0, 2.0)(rng));
  for (;;) {
    double c = (U01() * 2 - 1) * scale;
    if (positive) c = std::fabs(c) + 1e-3 * scale;
    double lo, hi;
    if (U01() < 0.4) { lo = hi = c; }
    else {
      double r = std::pow(10.0, std::uniform_real_distribution<double>(-14.0, -1.0)(rng)) * scale * U01();
      lo = c - r; hi = c + r;
      if (positive && lo < 0) lo = 0.0 + 1e-300;
    }
    if (nonzero && lo <= 0 && hi >= 0) continue;
    double x = lo + U01() * (hi - lo);
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return Pt{Interval(lo, hi), x};
  }
}
struct CPt { CInterval I; double xr, xi; };
static CPt rand_cpt(bool nonzero = false) {
  for (;;) {
    Pt r = rand_pt(), i = rand_pt();
    CInterval z(r.I, i.I);
    if (nonzero) { RoundUp g; if (!(abs2(z).lo > 0)) continue; }
    return CPt{z, r.x, i.x};
  }
}

static FILE* out;
static bool first = true;
static void h(double x) { std::fprintf(out, "\"%a\"", x); }
static void pt(const Pt& p) { std::fprintf(out, "["); h(p.I.lo); std::fprintf(out, ","); h(p.I.hi); std::fprintf(out, ","); h(p.x); std::fprintf(out, "]"); }
static void iv2(const Interval& a) { std::fprintf(out, "["); h(a.lo); std::fprintf(out, ","); h(a.hi); std::fprintf(out, "]"); }
static void cpt(const CPt& c) {
  std::fprintf(out, "{\"re\":"); pt(Pt{c.I.re, c.xr}); std::fprintf(out, ",\"im\":"); pt(Pt{c.I.im, c.xi}); std::fprintf(out, "}");
}
static void civ(const CInterval& a) { std::fprintf(out, "{\"re\":"); iv2(a.re); std::fprintf(out, ",\"im\":"); iv2(a.im); std::fprintf(out, "}"); }
static void begin(const char* op) { std::fprintf(out, "%s{\"op\":\"%s\"", first ? "" : ",\n", op); first = false; }
static void end() { std::fprintf(out, "}"); }

static void rec_real(int kind) {
  static const char* names[] = {"add", "sub", "mul", "div", "sqr", "sqrt"};
  RoundUp g;
  Pt a = rand_pt(kind == 5), b = rand_pt(false, kind == 3);
  Interval r;
  switch (kind) {
    case 0: r = a.I + b.I; break;
    case 1: r = a.I - b.I; break;
    case 2: r = a.I * b.I; break;
    case 3: r = a.I / b.I; break;
    case 4: r = sqr(a.I); break;
    default: r = iv::sqrt(a.I); break;
  }
  begin(names[kind]);
  std::fprintf(out, ",\"a\":"); pt(a);
  if (kind <= 3) { std::fprintf(out, ",\"b\":"); pt(b); }
  std::fprintf(out, ",\"out\":"); iv2(r);
  end();
}
static void rec_complex(int kind) {
  static const char* names[] = {"cmul", "cdiv", "cabs2"};
  RoundUp g;
  CPt a = rand_cpt(), b = rand_cpt(kind == 1);
  begin(names[kind]);
  std::fprintf(out, ",\"a\":"); cpt(a);
  if (kind < 2) { std::fprintf(out, ",\"b\":"); cpt(b); }
  std::fprintf(out, ",\"out\":");
  if (kind == 0) civ(a.I * b.I);
  else if (kind == 1) civ(a.I / b.I);
  else iv2(abs2(a.I));
  end();
}
static void vec_out(const CPt* v) {
  std::fprintf(out, "[");
  for (int j = 0; j < 6; ++j) { if (j) std::fprintf(out, ","); cpt(v[j]); }
  std::fprintf(out, "]");
}
static void rec_ip() {
  RoundUp g;
  CPt u[6], v[6];
  CVec6 U, V;
  for (int j = 0; j < 6; ++j) {
    u[j] = rand_cpt(); v[j] = rand_cpt();
    U[j] = u[j].I; V[j] = v[j].I;
  }
  CInterval p = inner(U, V);
  begin("ip");
  std::fprintf(out, ",\"u\":"); vec_out(u);
  std::fprintf(out, ",\"v\":"); vec_out(v);
  std::fprintf(out, ",\"ip\":"); civ(p);
  std::fprintf(out, ",\"abs2\":"); iv2(iv::abs2(p));
  std::fprintf(out, ",\"resid\":"); iv2(unbiased_residual(U, V));
  std::fprintf(out, ",\"norm2u\":"); iv2(norm2(U));
  end();
}
static void rec_const() {
  RoundUp g;
  static const long Ns[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 18, 24, 36, 60, 100, 360, 1000};
  long N = Ns[rng() % (sizeof(Ns) / sizeof(Ns[0]))];
  long k = (long)(rng() % (unsigned long)(3 * N + 1)) - N;  // includes negatives and k > N
  CInterval z = unit_root(k, N);
  begin("unit_root");
  std::fprintf(out, ",\"k\":%ld,\"N\":%ld,\"out\":", k, N);
  civ(z);
  end();
}
static void rec_named_consts() {
  RoundUp g;
  begin("consts");
  std::fprintf(out, ",\"one_sixth\":"); iv2(one_sixth());
  std::fprintf(out, ",\"inv_sqrt6\":"); iv2(inv_sqrt6());
  std::fprintf(out, ",\"sqrt6\":"); iv2(sqrt_iv(6));
  std::fprintf(out, ",\"sqrt3\":"); iv2(sqrt_iv(3));
  std::fprintf(out, ",\"pi\":"); iv2(pi_iv());
  end();
}
static void rec_mub() {
  // C++ enclosures of the three known bases and of all cross |<u,v>|^2.
  RoundUp g;
  Basis S = standard_basis(), F = fourier_basis(), C = bjorck_c6_basis();
  const Basis* Bs[3] = {&S, &F, &C};
  begin("mub");
  std::fprintf(out, ",\"bases\":[");
  for (int b = 0; b < 3; ++b) {
    if (b) std::fprintf(out, ",");
    std::fprintf(out, "[");
    for (int k = 0; k < 6; ++k) {
      if (k) std::fprintf(out, ",");
      std::fprintf(out, "[");
      for (int j = 0; j < 6; ++j) { if (j) std::fprintf(out, ","); civ((*Bs[b])[k][j]); }
      std::fprintf(out, "]");
    }
    std::fprintf(out, "]");
  }
  std::fprintf(out, "],\"abs2\":[");  // pairs (0,1),(0,2),(1,2) x 6 x 6
  int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  for (int p = 0; p < 3; ++p) {
    if (p) std::fprintf(out, ",");
    std::fprintf(out, "[");
    for (int a = 0; a < 6; ++a)
      for (int b = 0; b < 6; ++b) { if (a || b) std::fprintf(out, ","); iv2(ip_abs2((*Bs[pairs[p][0]])[a], (*Bs[pairs[p][1]])[b])); }
    std::fprintf(out, "]");
  }
  std::fprintf(out, "]");
  end();
}

int main(int argc, char** argv) {
  if (argc < 4) { std::fprintf(stderr, "usage: gen_records N seed out.json\n"); return 2; }
  long N = std::atol(argv[1]);
  rng.seed((unsigned long)std::atol(argv[2]));
  out = std::fopen(argv[3], "w");
  if (!out) return 2;
  std::fprintf(out, "[\n");
  rec_named_consts();
  rec_mub();
  for (long i = 0; i < N; ++i) {
    unsigned r = rng() % 100;
    if (r < 40) rec_real((int)(rng() % 6));
    else if (r < 65) rec_complex((int)(rng() % 3));
    else if (r < 85) rec_ip();
    else rec_const();
  }
  std::fprintf(out, "\n]\n");
  std::fclose(out);
  return 0;
}
