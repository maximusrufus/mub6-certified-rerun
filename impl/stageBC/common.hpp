// Shared types for Stage B / Stage C: Record I/O (SPEC section 7), interval helpers, small dense linear algebra.
// Built on ../core/interval.hpp.  Every function that touches intervals must run inside an iv::RoundUp guard
// (one per thread).  Plain-double code in here is HEURISTIC ONLY (proposals); soundness lives in Interval code.
#pragma once
#include "../core/interval.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace mb {
using iv::CInterval;
using iv::Interval;

// ---------------------------------------------------------------- Record (SPEC section 7)
struct Record {
  int32_t k = 0;
  double p[36];
  std::vector<double> T;  // k*36, row a = column a of T (flattened 6i+j)
  double C = 0, n0 = 0, id = 0, ort = 0, hor = 0, mnorm = 0, froot = 0, r = 0;
  Record() { std::memset(p, 0, sizeof p); }
  const double* col(int a) const { return T.data() + 36 * (size_t)a; }
};

inline bool readRecords(const std::string& path, std::vector<Record>& out, std::string& err) {
  FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) { err = "cannot open " + path; return false; }
  uint64_t n = 0;
  if (std::fread(&n, 8, 1, f) != 1) { std::fclose(f); err = "short header " + path; return false; }
  if (n > 2000000ull) { std::fclose(f); err = "more than 2,000,000 records in " + path; return false; }
  out.reserve(out.size() + n);
  for (uint64_t t = 0; t < n; ++t) {
    Record R;
    int32_t k;
    if (std::fread(&k, 4, 1, f) != 1 || k < 0 || k > 10) { std::fclose(f); err = "bad k in " + path; return false; }
    R.k = k;
    if (std::fread(R.p, 8, 36, f) != 36) { std::fclose(f); err = "short p in " + path; return false; }
    R.T.resize((size_t)k * 36);
    if (k && std::fread(R.T.data(), 8, (size_t)k * 36, f) != (size_t)k * 36) { std::fclose(f); err = "short T in " + path; return false; }
    double c[8];
    if (std::fread(c, 8, 8, f) != 8) { std::fclose(f); err = "short consts in " + path; return false; }
    R.C = c[0]; R.n0 = c[1]; R.id = c[2]; R.ort = c[3]; R.hor = c[4]; R.mnorm = c[5]; R.froot = c[6]; R.r = c[7];
    out.push_back(std::move(R));
  }
  unsigned char extra;
  if (std::fread(&extra, 1, 1, f) != 0) { std::fclose(f); err = "trailing bytes in " + path; return false; }
  std::fclose(f);
  return true;
}

inline bool writeRecords(const std::string& path, const std::vector<const Record*>& recs) {
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) return false;
  uint64_t n = recs.size();
  bool ok = std::fwrite(&n, 8, 1, f) == 1;
  for (const Record* R : recs) {
    int32_t k = R->k;
    ok = ok && std::fwrite(&k, 4, 1, f) == 1 && std::fwrite(R->p, 8, 36, f) == 36;
    if (k) ok = ok && std::fwrite(R->T.data(), 8, (size_t)k * 36, f) == (size_t)k * 36;
    double c[8] = {R->C, R->n0, R->id, R->ort, R->hor, R->mnorm, R->froot, R->r};
    ok = ok && std::fwrite(c, 8, 8, f) == 8;
  }
  ok = (std::fclose(f) == 0) && ok;
  return ok;
}
inline bool writeRecords(const std::string& path, const std::vector<Record>& recs) {
  std::vector<const Record*> v;
  for (auto& r : recs) v.push_back(&r);
  return writeRecords(path, v);
}

// SPEC 4.1 validation.  maxR: .126 for B; .03 for C.
inline bool validateRecord(const Record& R, double maxR, std::string& why) {
  auto bad = [&](const char* s) { why = s; return false; };
  if (!(R.r >= 0 && R.r <= maxR)) return bad("radius out of range");
  if (!(R.C >= 0 && R.C < 20)) return bad("C out of range");
  if (!(R.n0 >= 0 && R.n0 < 1e-6)) return bad("n0 out of range");
  if (!(R.froot >= 0 && R.froot < 1e-6)) return bad("froot out of range");
  if (!(R.id >= 0 && R.id < 1e-4)) return bad("id out of range");
  if (!(R.ort >= 0 && R.ort < 1e-5)) return bad("ort out of range");
  if (!(R.hor >= 0 && R.hor < 1e-5)) return bad("hor out of range");
  if (!(R.mnorm >= 0 && R.mnorm < 2)) return bad("mnorm out of range");
  if (R.k < 0 || R.k > 10) return bad("k out of range");
  for (int j = 0; j < 36; ++j) if (!(std::fabs(R.p[j]) <= 4)) return bad("|p| > 4 or NaN");
  for (double t : R.T) if (!(std::fabs(t) <= 1.01)) return bad("|T| > 1.01 or NaN");
  return true;
}

// ---------------------------------------------------------------- interval helpers
inline Interval imax(const Interval& a, const Interval& b) { return Interval(std::max(a.lo, b.lo), std::max(a.hi, b.hi)); }
inline Interval imin(const Interval& a, const Interval& b) { return Interval(std::min(a.lo, b.lo), std::min(a.hi, b.hi)); }
// sqrt(sum of squares) of n intervals
inline Interval inorm(const Interval* x, int n) {
  Interval s(0.0);
  for (int i = 0; i < n; ++i) s = s + iv::sqr(x[i]);
  return iv::sqrt(s);
}
inline Interval isumsq(const Interval* x, int n) {
  Interval s(0.0);
  for (int i = 0; i < n; ++i) s = s + iv::sqr(x[i]);
  return s;
}

struct Consts {
  Interval pi, twopi, sqrt2, sqrt6, isqrt6, sqrt8, sixth, th36, K, one;
  Interval sqrt29_6, c4pi, sqrt6o2;
};
inline const Consts& cn() {
  static const Consts c = [] {
    iv::RoundUp g;
    Consts c;
    c.pi = iv::pi_iv();
    c.twopi = c.pi * Interval(2.0);
    c.sqrt2 = iv::sqrt_iv(2);
    c.sqrt6 = iv::sqrt_iv(6);
    c.isqrt6 = iv::inv_sqrt6();
    c.sqrt8 = iv::sqrt_iv(8);
    c.sixth = iv::ratio(1, 6);
    c.th36 = iv::ratio(1, 36);
    c.K = iv::ratio(100001, 100000);  // K := 1.00001 (exclusion.tex exclude:operator-bound)
    c.one = Interval(1.0);
    c.sqrt29_6 = iv::sqrt(Interval(29.0) / Interval(6.0));
    c.c4pi = c.pi * Interval(4.0);
    c.sqrt6o2 = c.sqrt6 / Interval(2.0);
    return c;
  }();
  return c;
}

// ---------------------------------------------------------------- phase-matrix helpers (SPEC 1.2, 3.6)
inline int ix(int i, int j) { return 6 * i + j; }

// Pi z = z - rowmean - colmean + totalmean (SPEC 1.2), interval version.
inline void PiI(const Interval* z, Interval* o) {
  const Consts& c = cn();
  Interval rs[6], cs[6], tot(0.0);
  for (int i = 0; i < 6; ++i) {
    rs[i] = z[ix(i, 0)];
    for (int j = 1; j < 6; ++j) rs[i] = rs[i] + z[ix(i, j)];
  }
  for (int j = 0; j < 6; ++j) {
    cs[j] = z[ix(0, j)];
    for (int i = 1; i < 6; ++i) cs[j] = cs[j] + z[ix(i, j)];
  }
  for (int i = 0; i < 6; ++i) tot = tot + rs[i];
  Interval t36 = tot * c.th36;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) o[ix(i, j)] = z[ix(i, j)] - rs[i] * c.sixth - cs[j] * c.sixth + t36;
}
// double version (heuristic)
inline void PiD(const double* z, double* o) {
  double rs[6] = {0}, cs[6] = {0}, tot = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) { rs[i] += z[ix(i, j)]; cs[j] += z[ix(i, j)]; tot += z[ix(i, j)]; }
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) o[ix(i, j)] = z[ix(i, j)] - rs[i] / 6 - cs[j] / 6 + tot / 36;
}
// unwrapped dephasing D0(P)_ij = P_ij - P_i0 - P_0j + P_00, exact doubles -> intervals
inline void D0I(const double* P, Interval* o) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) o[ix(i, j)] = Interval(P[ix(i, j)]) - Interval(P[ix(i, 0)]) - Interval(P[ix(0, j)]) + Interval(P[ix(0, 0)]);
}
inline void D0D(const double* P, double* o) {
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) o[ix(i, j)] = P[ix(i, j)] - P[ix(i, 0)] - P[ix(0, j)] + P[ix(0, 0)];
}
inline double rnd(double x) { return std::floor(x + 0.5); }  // any integer is valid (phase:alignment); avoids mode-dependent nearbyint

// Aligned displacement d = Pi(D0 Q - D0 P - N), N = rounded wraps.  Returns 36 intervals.
inline void alignedDisp(const double* P, const double* Q, Interval* d) {
  double dp[36], dq[36];
  D0D(P, dp);
  D0D(Q, dq);
  Interval ip[36], iq[36], v[36];
  D0I(P, ip);
  D0I(Q, iq);
  for (int m = 0; m < 36; ++m) {
    double N = rnd(dq[m] - dp[m]);
    v[m] = iq[m] - ip[m] - Interval(N);
  }
  PiI(v, d);
}

// ---------------------------------------------------------------- tiny dense linear algebra (double, heuristic)
// Symmetric Jacobi eigen-decomposition.  A is n*n row-major (destroyed); V (n*n) columns are eigenvectors; ev ascending.
inline void jacobiEig(std::vector<double>& A, int n, std::vector<double>& V, std::vector<double>& ev) {
  V.assign((size_t)n * n, 0.0);
  for (int i = 0; i < n; ++i) V[(size_t)i * n + i] = 1.0;
  for (int sweep = 0; sweep < 60; ++sweep) {
    double off = 0;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) off += A[(size_t)i * n + j] * A[(size_t)i * n + j];
    if (off < 1e-32) break;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) {
        double apq = A[(size_t)p * n + q];
        if (std::fabs(apq) < 1e-300) continue;
        double th = (A[(size_t)q * n + q] - A[(size_t)p * n + p]) / (2 * apq);
        double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
        double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < n; ++k) {
          double akp = A[(size_t)k * n + p], akq = A[(size_t)k * n + q];
          A[(size_t)k * n + p] = c * akp - s * akq;
          A[(size_t)k * n + q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          double apk = A[(size_t)p * n + k], aqk = A[(size_t)q * n + k];
          A[(size_t)p * n + k] = c * apk - s * aqk;
          A[(size_t)q * n + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          double vkp = V[(size_t)k * n + p], vkq = V[(size_t)k * n + q];
          V[(size_t)k * n + p] = c * vkp - s * vkq;
          V[(size_t)k * n + q] = s * vkp + c * vkq;
        }
      }
  }
  std::vector<int> ord(n);
  for (int i = 0; i < n; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](int a, int b) { return A[(size_t)a * n + a] < A[(size_t)b * n + b]; });
  std::vector<double> V2((size_t)n * n);
  ev.resize(n);
  for (int c = 0; c < n; ++c) {
    ev[c] = A[(size_t)ord[c] * n + ord[c]];
    for (int r = 0; r < n; ++r) V2[(size_t)r * n + c] = V[(size_t)r * n + ord[c]];
  }
  V = V2;
}

// Solve A x = b by Gaussian elimination with pivoting.  Returns false if singular / non-finite.
inline bool solveLin(std::vector<double> A, std::vector<double> b, int n, std::vector<double>& x) {
  for (int c = 0; c < n; ++c) {
    int piv = c;
    for (int r = c + 1; r < n; ++r)
      if (std::fabs(A[(size_t)r * n + c]) > std::fabs(A[(size_t)piv * n + c])) piv = r;
    if (!(std::fabs(A[(size_t)piv * n + c]) > 1e-300)) return false;
    if (piv != c) {
      for (int k = 0; k < n; ++k) std::swap(A[(size_t)c * n + k], A[(size_t)piv * n + k]);
      std::swap(b[c], b[piv]);
    }
    for (int r = c + 1; r < n; ++r) {
      double f = A[(size_t)r * n + c] / A[(size_t)c * n + c];
      if (f == 0) continue;
      for (int k = c; k < n; ++k) A[(size_t)r * n + k] -= f * A[(size_t)c * n + k];
      b[r] -= f * b[c];
    }
  }
  x.assign(n, 0.0);
  for (int r = n - 1; r >= 0; --r) {
    double s = b[r];
    for (int k = r + 1; k < n; ++k) s -= A[(size_t)r * n + k] * x[k];
    x[r] = s / A[(size_t)r * n + r];
  }
  for (int i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
  return true;
}

}  // namespace mb
