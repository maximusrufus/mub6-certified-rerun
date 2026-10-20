// root.cpp -- Root certificates (SPEC 3.6, phase-cover.tex phase:root-bounds).
// Construction (A9, decided in DECISIONS.md): eigen-decomposition of the 25x25 Gram
// A^T A of the horizontal Jacobian A = L_P Bh, with Bh an orthonormal basis of the
// horizontal space.  Eigenvalues < .64 -> columns of T; the others give
//   H_V = Bh v / sqrt(lam),  L_V = A v / sqrt(lam),   M := H_V L_V^T  (exact product
// of the STORED doubles).  All constants are interval-evaluated; every proposal is
// re-checked by the Proposition phase:root-bounds inequalities.
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include "la.hpp"
#include "../core/tensor_c.hpp"
#include "stagea.hpp"

namespace sa {

// ---------------------------------------------------------------- basics
void piApply(const I in[36], I out[36]) {
  I rm[6], cm[6];
  for (int i = 0; i < 6; ++i) {
    I s(0.0);
    for (int j = 0; j < 6; ++j) s = s + in[6 * i + j];
    rm[i] = s * SIXTH();
  }
  for (int j = 0; j < 6; ++j) {
    I s(0.0);
    for (int i = 0; i < 6; ++i) s = s + in[6 * i + j];
    cm[j] = s * SIXTH();
  }
  I tm(0.0);
  for (int i = 0; i < 6; ++i) tm = tm + rm[i];
  tm = tm * SIXTH();
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) out[6 * i + j] = in[6 * i + j] - rm[i] - cm[j] + tm;
}

void alignDelta(const double P0[36], const double P[36], I Delta[36]) {
  I W[36];
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) {
      I d0 = I(P[6 * i + j]) - I(P[6 * i]) - I(P[j]) + I(P[0]);
      double mid = midp(d0);
      double N = std::floor(P0[6 * i + j] - mid + 0.5);
      W[6 * i + j] = I(P0[6 * i + j]) - d0 - I(N);
    }
  piApply(W, Delta);
}

namespace {

const double* Bh() {  // 36 x 25 row-major; thread-safe init (C++11 magic static), no data race under OpenMP
  static const std::array<double, 36 * 25> B = [] {
    std::array<double, 36 * 25> Bv{};
    double H[6][5];
    for (int u = 0; u < 5; ++u) {
      double nrm = std::sqrt((u + 1.0) * (u + 2.0));
      for (int i = 0; i < 6; ++i) H[i][u] = (i <= u) ? 1.0 / nrm : (i == u + 1 ? -(u + 1.0) / nrm : 0.0);
    }
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
        for (int u = 0; u < 5; ++u)
          for (int v = 0; v < 5; ++v) Bv[(6 * i + j) * 25 + 5 * u + v] = H[i][u] * H[j][v];
    return Bv;
  }();
  return B.data();
}

// doubles: F (30) and L (30x36)
void FL_double(const double P[36], double F[30], std::vector<double>* L) {
  if (L) L->assign(30 * 36, 0.0);
  const double TW = 2.0 * M_PI;
  int e = 0;
  for (int i = 1; i < 6; ++i)
    for (int j = 0; j < i; ++j, ++e) {
      double sx = 0, sy = 0;
      for (int g = 0; g < 6; ++g) {
        double ph = TW * (P[6 * i + g] - P[6 * j + g]);
        double c = std::cos(ph), s = std::sin(ph);
        sx += c;
        sy += s;
        if (L) {
          (*L)[(2 * e) * 36 + 6 * i + g] += c;
          (*L)[(2 * e) * 36 + 6 * j + g] -= c;
          (*L)[(2 * e + 1) * 36 + 6 * i + g] += s;
          (*L)[(2 * e + 1) * 36 + 6 * j + g] -= s;
        }
      }
      F[2 * e] = sy / TW;
      F[2 * e + 1] = -sx / TW;
    }
}

struct EdgeTab {
  Z E[15][6];  // e(P_ig - P_jg) enclosures
};
void edgeTab(const double P[36], EdgeTab& T) {
  int e = 0;
  for (int i = 1; i < 6; ++i)
    for (int j = 0; j < i; ++j, ++e)
      for (int g = 0; g < 6; ++g) T.E[e][g] = iv::cis_turns(I(P[6 * i + g]) - I(P[6 * j + g]));
}

// interval residual: returns sum of squares (upper) and F as intervals
I residualIv(const EdgeTab& T, I F[30]) {
  I val(0.0);
  int e = 0;
  I itp = recip(TWOPI());
  for (int i = 1; i < 6; ++i)
    for (int j = 0; j < i; ++j, ++e) {
      I sx(0.0), sy(0.0);
      for (int g = 0; g < 6; ++g) {
        sx = sx + T.E[e][g].re;
        sy = sy + T.E[e][g].im;
      }
      F[2 * e] = sy * itp;
      F[2 * e + 1] = -(sx * itp);
      val = val + sqr(F[2 * e]) + sqr(F[2 * e + 1]);
    }
  return val;
}

}  // namespace

double residualUpper(const double P[36]) {
  try {
    EdgeTab T;
    edgeTab(P, T);
    I F[30];
    return sqrt0(residualIv(T, F)).hi;
  } catch (...) {
    return 1e300;
  }
}

bool useTensorC() { static const bool v = std::getenv("MUB_NO_TENSOR_C") == nullptr; return v; }

namespace {

struct Cert {
  Root R;
  bool ok = false;
};

// Certify center P.  If Tin != nullptr use it (k columns), else propose T by the eigen rule.
bool certify(const double P[36], const double (*Tin)[36], int kin, Root& out, std::string* why, bool inflate) {
  auto fail = [&](const char* m) {
    if (why) *why = m;
    return false;
  };
  for (int i = 0; i < 36; ++i)
    if (!(std::fabs(P[i]) <= 4.0)) return fail("|P_j|>4");
  std::vector<double> Ld;
  double Fd[30];
  FL_double(P, Fd, &Ld);
  const double* B = Bh();
  // A = L Bh (30x25), G = A^T A
  std::vector<double> A(30 * 25, 0.0), G(25 * 25, 0.0);
  for (int r = 0; r < 30; ++r)
    for (int c = 0; c < 25; ++c) {
      double s = 0;
      for (int p = 0; p < 36; ++p) s += Ld[r * 36 + p] * B[p * 25 + c];
      A[r * 25 + c] = s;
    }
  for (int a = 0; a < 25; ++a)
    for (int b = a; b < 25; ++b) {
      double s = 0;
      for (int r = 0; r < 30; ++r) s += A[r * 25 + a] * A[r * 25 + b];
      G[a * 25 + b] = G[b * 25 + a] = s;
    }
  std::vector<double> ev, V;
  la::jacobi(G, 25, ev, V);
  std::vector<int> tIdx, vIdx;
  for (int i = 0; i < 25; ++i) (ev[i] < 0.64 ? tIdx : vIdx).push_back(i);
  int k = (int)tIdx.size();
  if (k > 10) return fail("k>10");
  if (Tin && kin != k) return fail("k mismatch");
  int g = (int)vIdx.size();
  Root R;
  R.k = k;
  for (int i = 0; i < 36; ++i) R.P[i] = P[i];
  // T columns
  for (int a = 0; a < k; ++a)
    for (int p = 0; p < 36; ++p) {
      if (Tin) {
        R.T[a][p] = Tin[a][p];
      } else {
        double s = 0;
        for (int c = 0; c < 25; ++c) s += B[p * 25 + c] * V[c * 25 + tIdx[a]];
        R.T[a][p] = s;
      }
    }
  // H_V (36 x g), L_V (30 x g)
  std::vector<double> HV((size_t)36 * g), LV((size_t)30 * g);
  for (int a = 0; a < g; ++a) {
    int col = vIdx[a];
    double sl = std::sqrt(ev[col]);
    if (!(sl > 0.5)) return fail("tiny eigenvalue in V part");
    for (int p = 0; p < 36; ++p) {
      double s = 0;
      for (int c = 0; c < 25; ++c) s += B[p * 25 + c] * V[c * 25 + col];
      HV[(size_t)p * g + a] = s / sl;
      if (!(std::fabs(HV[(size_t)p * g + a]) <= 2.0)) return fail("H_V entry >2");
    }
    for (int r = 0; r < 30; ++r) {
      double s = 0;
      for (int c = 0; c < 25; ++c) s += A[r * 25 + c] * V[c * 25 + col];
      LV[(size_t)r * g + a] = s / sl;
      if (!(std::fabs(LV[(size_t)r * g + a]) <= 2.0)) return fail("L_V entry >2");
    }
  }
  for (int a = 0; a < k; ++a)
    for (int p = 0; p < 36; ++p)
      if (!(std::fabs(R.T[a][p]) <= 2.0)) return fail("T entry >2");
  // ---- interval certificates
  EdgeTab ET;
  edgeTab(P, ET);
  I F[30];
  I val = residualIv(ET, F);
  if (!(val.hi < 1e-18)) return fail("residual not small");
  R.froot = sqrt0(val).hi;
  // ort, hor
  {
    I s(0.0);
    for (int a = 0; a < k; ++a)
      for (int b = 0; b < k; ++b) {
        I d(0.0);
        for (int p = 0; p < 36; ++p) d = d + mulp(R.T[a][p], I(R.T[b][p]));
        if (a == b) d = d - I(1.0);
        s = s + sqr(d);
      }
    R.ort = sqrt0(s).hi;
    I sh(0.0);
    for (int a = 0; a < k; ++a) {
      I col[36], pc[36];
      for (int p = 0; p < 36; ++p) col[p] = I(R.T[a][p]);
      piApply(col, pc);
      for (int p = 0; p < 36; ++p) sh = sh + sqr(col[p] - pc[p]);
    }
    R.hor = sqrt0(sh).hi;
  }
  // mnorm
  {
    auto rowsum = [&](const std::vector<double>& X, int rows) {
      double best = 0;
      std::vector<I> gram((size_t)g * g);
      for (int a = 0; a < g; ++a)
        for (int b = 0; b < g; ++b) {
          I d(0.0);
          for (int r = 0; r < rows; ++r) d = d + mulp(X[(size_t)r * g + a], I(X[(size_t)r * g + b]));
          gram[(size_t)a * g + b] = d;
        }
      for (int a = 0; a < g; ++a) {
        I s(0.0);
        for (int b = 0; b < g; ++b) s = s + iv::abs(gram[(size_t)a * g + b]);
        best = std::max(best, s.hi);
      }
      return best;
    };
    double ah = rowsum(HV, 36), al = rowsum(LV, 30);
    R.mnorm = sqrt0(I(ah) * I(al)).hi;
  }
  if (!(R.mnorm < 2.0)) return fail("mnorm>=2");
  R.n0 = (I(R.mnorm) * I(R.froot)).hi;
  // C = pi * min(sqrt8 ||M||, sqrt6 ||S||): both are certified upper bounds of a remainder constant
  // (phase:remainder-constants); the tensor bound is core/tensor_c.hpp (interval, no paper margins).
  double Cplain = (PI() * SQRT8() * I(R.mnorm)).hi;
  double Ctensor = useTensorC() ? tc::tensorC(HV, LV, g, PI()) : Cplain;
  R.C = std::min(Cplain, Ctensor);
  if (!(R.C < 20.0)) return fail("C>=20");
  // id = || Pi - T T^* - M L_P Pi ||_F, M = H_V L_V^T
  {
    // L Pi (30 x 36)
    std::vector<I> LP(30 * 36);
    int e = 0;
    for (int i = 1; i < 6; ++i)
      for (int j = 0; j < i; ++j, ++e) {
        I rowc[36], rows[36], oc[36], os[36];
        for (int p = 0; p < 36; ++p) rowc[p] = rows[p] = I(0.0);
        for (int gg = 0; gg < 6; ++gg) {
          rowc[6 * i + gg] = ET.E[e][gg].re;
          rowc[6 * j + gg] = -ET.E[e][gg].re;
          rows[6 * i + gg] = ET.E[e][gg].im;
          rows[6 * j + gg] = -ET.E[e][gg].im;
        }
        piApply(rowc, oc);
        piApply(rows, os);
        for (int p = 0; p < 36; ++p) {
          LP[(2 * e) * 36 + p] = oc[p];
          LP[(2 * e + 1) * 36 + p] = os[p];
        }
      }
    std::vector<I> W((size_t)g * 36);
    for (int a = 0; a < g; ++a)
      for (int c = 0; c < 36; ++c) {
        I s(0.0);
        for (int r = 0; r < 30; ++r) s = s + mulp(LV[(size_t)r * g + a], LP[r * 36 + c]);
        W[(size_t)a * 36 + c] = s;
      }
    I fro(0.0);
    I m6 = I(1.0) - SIXTH();   // 1 - 1/6
    I m6n = -SIXTH();
    for (int r = 0; r < 36; ++r)
      for (int c = 0; c < 36; ++c) {
        int ri = r / 6, rj = r % 6, ci = c / 6, cj = c % 6;
        I fi = (ri == ci) ? m6 : m6n, fj = (rj == cj) ? m6 : m6n;
        I ent = fi * fj;  // Pi entry
        for (int a = 0; a < k; ++a) ent = ent - mulp(R.T[a][r], I(R.T[a][c]));
        I mp(0.0);
        for (int a = 0; a < g; ++a) mp = mp + mulp(HV[(size_t)r * g + a], W[(size_t)a * 36 + c]);
        ent = ent - mp;
        fro = fro + sqr(ent);
      }
    R.id = sqrt0(fro).hi;
  }
  if (!(R.id < 1e-4)) return fail("id>=1e-4");
  if (!(R.ort < 1e-5)) return fail("ort>=1e-5");
  if (!(R.hor < 1e-5)) return fail("hor>=1e-5");
  if (inflate) {
    // Re-verification robustness (NOT needed for soundness: every constant is only an upper bound, and larger
    // constants are still valid).  A Root that is re-certified from scratch after a symmetry Phi (different
    // rounding in the eigensolve / interval sums) must still dominate the fresh constants.
    R.froot = (I(R.froot) + I(1e-13)).hi;
    R.id = (I(R.id) + I(1e-11)).hi;
    R.ort = (I(R.ort) + I(1e-12)).hi;
    R.hor = (I(R.hor) + I(1e-12)).hi;
    R.mnorm = (I(R.mnorm) * I(1.0 + 1e-9)).hi;
    R.n0 = (I(R.mnorm) * I(R.froot)).hi;
    R.C = std::min((PI() * SQRT8() * I(R.mnorm)).hi, (I(Ctensor) * I(1.0 + 1e-9)).hi);
  }
  out = R;
  return true;
}

}  // namespace

bool recertifyTransported(Root& R, std::string* why) {
  Root fresh;
  try {
    if (!certify(R.P, R.k ? R.T : nullptr, R.k, fresh, why, true)) return false;
  } catch (const std::exception& e) {
    if (why) *why = e.what();
    return false;
  }
  if (fresh.k != R.k) { if (why) *why = "k differs"; return false; }
  double r = R.r;
  R = fresh;
  R.r = r;
  return true;
}

bool checkRoot(const Root& R, std::string* why) {
  Root fresh;
  try {
    if (!certify(R.P, R.k ? R.T : nullptr, R.k, fresh, why, false)) return false;
  } catch (const std::exception& e) {
    if (why) *why = e.what();
    return false;
  }
  // when k==0 certify(Tin=nullptr) proposes T which must also be empty (k mismatch handled via tIdx)
  if (fresh.k != R.k) {
    if (why) *why = "k differs";
    return false;
  }
  try {
    if (!(R.froot >= fresh.froot) || !(R.ort >= fresh.ort) || !(R.hor >= fresh.hor) || !(R.id >= fresh.id) ||
        !(R.mnorm >= fresh.mnorm) || !(R.C >= fresh.C)) {
      if (why) *why = "stored constant smaller than recomputed";
      return false;
    }
    I prod = I(R.mnorm) * I(R.froot);
    if (!(R.n0 >= prod.hi)) {
      if (why) *why = "n0 < mnorm*froot";
      return false;
    }
    // (C >= fresh.C above: fresh.C is the certified min of the plain and the tensor remainder constants)
  } catch (...) {
    return false;
  }
  return true;
}

bool buildRoot(const double P0[36], Root& out, std::string* why) {
  try {
    double P[36];
    for (int i = 0; i < 36; ++i) P[i] = P0[i];
    const double* B = Bh();
    double lastF = 1e300;
    for (int it = 0; it < 25; ++it) {
      double F[30];
      std::vector<double> L;
      FL_double(P, F, &L);
      double fn = 0;
      for (int i = 0; i < 30; ++i) fn += F[i] * F[i];
      fn = std::sqrt(fn);
      if (fn < 2e-15) break;
      if (!(fn < 0.5)) {
        if (why) *why = "snap: residual too large";
        return false;
      }
      if (it > 3 && fn > 0.9 * lastF && fn > 1e-9) {
        if (why) *why = "snap: stagnation";
        return false;
      }
      lastF = fn;
      std::vector<double> A(30 * 25), G(25 * 25);
      for (int r = 0; r < 30; ++r)
        for (int c = 0; c < 25; ++c) {
          double s = 0;
          for (int p = 0; p < 36; ++p) s += L[r * 36 + p] * B[p * 25 + c];
          A[r * 25 + c] = s;
        }
      for (int a = 0; a < 25; ++a)
        for (int b = a; b < 25; ++b) {
          double s = 0;
          for (int r = 0; r < 30; ++r) s += A[r * 25 + a] * A[r * 25 + b];
          G[a * 25 + b] = G[b * 25 + a] = s;
        }
      std::vector<double> ev, V;
      la::jacobi(G, 25, ev, V);
      double y[25] = {0};
      for (int c = 0; c < 25; ++c) {
        if (ev[c] < 0.64) continue;
        double s = 0;  // v_c^T A^T F
        for (int r = 0; r < 30; ++r) {
          double av = 0;
          for (int a = 0; a < 25; ++a) av += A[r * 25 + a] * V[a * 25 + c];
          s += av * F[r];
        }
        for (int a = 0; a < 25; ++a) y[a] -= V[a * 25 + c] * s / ev[c];
      }
      for (int p = 0; p < 36; ++p) {
        double s = 0;
        for (int a = 0; a < 25; ++a) s += B[p * 25 + a] * y[a];
        P[p] += s;
      }
    }
    return certify(P, nullptr, 0, out, why, true);
  } catch (const std::exception& e) {
    if (why) *why = e.what();
    return false;
  }
}

}  // namespace sa
