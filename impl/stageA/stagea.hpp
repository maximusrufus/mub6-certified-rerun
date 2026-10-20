// stagea.hpp -- Stage A (arun) declarations. SPEC.md section 3.
// Every certificate inequality is evaluated in outward-rounded interval
// arithmetic (../core/interval.hpp); heuristics (eigensolves, Newton, direction
// search) only PROPOSE.  Call everything inside an iv::RoundUp guard.
#pragma once
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>
#include "ext.hpp"

namespace sa {

constexpr int SEEDPOS[9] = {7, 8, 9, 10, 11, 13, 19, 25, 31};
inline int dpos(int f) { return 6 * (2 + f / 4) + 2 + (f % 4); }  // f = 4*(i-2)+(j-2)

// ---------------------------------------------------------------- Root/Record
struct Root {
  int k = 0;
  double P[36] = {0};
  double T[10][36] = {{0}};  // T[a] = column a of T
  double C = 0, n0 = 0, id = 0, ort = 0, hor = 0, mnorm = 0, froot = 0;
  double r = 0;  // Record radius (cycles)
};

// Build a certified Root near P0 (Gauss-Newton snap + interval certificate).
// Returns false (never throws) if any check fails.
bool buildRoot(const double P0[36], Root& out, std::string* why = nullptr);
// Re-certify the data of an existing Root with fresh interval evaluations
// (independent of how it was built).  Used by tests and by the output check.
bool checkRoot(const Root& R, std::string* why = nullptr);
// Re-certify a Phi-transported Root from scratch (certify with its stored T) and replace its constants by the
// fresh (inflated) ones, so that checkRoot on the stored orientation is deterministic.  The transported constants
// are valid by phase:transport anyway; on failure they are kept.  Radius unchanged.
bool recertifyTransported(Root& R, std::string* why = nullptr);
// MUB_NO_TENSOR_C=1 disables the tensor improvement of C (measurement switch; both choices are valid)
bool useTensorC();
// ||F(P)|| upper bound
double residualUpper(const double P[36]);

// ----------------------------------------------------------------- seed data
struct Params {
  double target = 0.103;     // local target radius
  double cap = 0.125;        // global atlas cap
  double rsCycles = 0.0320001;
  int maxNu = 512;
  int baseNu = 16;
  int firstDirNu = 32;       // directional test starts here (A6)
  int firstCoverNu = 64;     // phase-ball attempts start here
  long extraCap = 1000000;   // nodes beyond the list-held levels
  int topK = 4;
  double mFrameMin = 1.24;   // lower bound on m (A24)
  bool verbose = false;
  // S6: deterministic seeded sample of EXCLUDED cells for independent python-flint re-verification (0 = off)
  unsigned sampleMod = 0;   // keep a cell iff splitmix64(sampleSeed ^ key) % sampleMod == 0
  uint64_t sampleSeed = 0x4d55425f53365f31ULL;
};

struct LevelK {
  int Nu = 0;
  double h = 0;   // upper bound of pi/Nu (radians)
  double p = 0, q = 0, z2 = 0, Er = 0;
};

struct SeedCtx {
  double q[9];
  Z a, b[4], c[4];
  CM<2, 2> E0;
  CM<2, 4> B0;
  CM<4, 2> C0;
  CM<4, 2> KB, KC, GB, GC;
  CM<2, 4> KBs;   // K_B^*
  CM<4, 2> KC6;   // sqrt6 K_C
  CM<4, 4> PK, Dspan;
  CM<4, 4> Cst[9], Lft[9], Lc[9];
  I mB, mC, m, A0, Astar, Cop;
  int fc[4] = {0, 0, 0, 0};  // Gram-Schmidt coordinate choices (bk, bl) of the B and C frames (S6 dump)
  double Rs = 0;  // upper bound 2*pi*rs
  std::vector<LevelK> lev;  // indexed by log2(Nu/16)
  // seed-sum rows: c_X, c_Y (real pairs) as intervals
  I cX[2], cY[2];
  const LevelK& L(int Nu) const;
};

bool buildSeed(SeedCtx& S, const double q[9], const Params& P, std::string& why);

struct Cell {
  uint16_t ia, it, ib, ig;
};

struct CellWork {
  Cell c;
  int Nu = 0;
  Z U[2][2];
  CM<4, 4> D;
  bool haveJ = false;
  CM<4, 4> Jc[13];  // linearisation columns (9 seed + 4 angles)
};

void initTables();  // once, before threads (inside RoundUp)
void cellD(const SeedCtx& S, int Nu, const Cell& c, CellWork& w);
bool normExcluded(const SeedCtx& S, const LevelK& lk, const CellWork& w);
void ensureJ(const SeedCtx& S, CellWork& w);
bool dirExcluded(const SeedCtx& S, const LevelK& lk, CellWork& w);
bool seedSumExcluded(const SeedCtx& S);  // whole-seed-ball test from the X/Y sums only

// ----------------------------------------------------- phase-ball covering
struct CoverResult {
  bool ok = false;
  double s = 0;   // certified radius bound for the (possibly new) Root
  int rootIndex = -1;
};

// ---------------------------------------------------------- seed processing
struct Witness {
  bool on = false;
  const double* truthP = nullptr;        // exact dephased phases (36) of the true completion
  long coverChecks = 0, coverViolations = 0;
  double worstRatio = 0;
  double al = 0, t = 0, be = 0, ga = 0;  // true polar cell coordinates (cycles)
  long violations = 0;                   // times a cell holding the witness was excluded
};

struct SeedStats {
  long cellsNorm = 0, cellsDir = 0, exclNorm = 0, exclDir = 0, covered = 0, coverTry = 0, split = 0;
  long deep = 0;
  int maxNuReached = 0;
  long nPrepFail = 0, nReuse = 0, nSnapOK = 0, nBuildFail = 0, nSnapBound = 0, nBoundTry = 0;
  double tPrep = 0, tReuse = 0, tBuild = 0, tSnapBound = 0, tNorm = 0, tDir = 0, tCover = 0;
};

// one sampled excluded cell (kind 0 = norm test, 1 = directional test with certificate y[20])
struct ExclSample {
  int kind = 0, Nu = 0;
  uint16_t ia = 0, it = 0, ib = 0, ig = 0;
  int fc[4] = {0, 0, 0, 0};
  double h = 0, z2 = 0, Er = 0, Rs = 0;  // constants the C++ test used (claims; checked >= recomputed bounds)
  double y[20] = {0};
};
extern thread_local double g_dirCertY[20];  // certificate direction of the last successful dirExcluded (S6 dump)

struct SeedResult {
  int status = 0;  // 0 ok, 1 setup failure (frames/guards), 2 unresolved cells
  std::string why;
  long unresolved = 0;
  std::vector<Root> recs;  // local Records (identity gauge), radius = certified bound
  SeedStats st;
  bool seedExcluded = false;
  std::vector<ExclSample> samples;  // empty unless Params::sampleMod > 0
};

SeedResult processSeed(const double q[9], const Params& P, Witness* wit = nullptr);

// ------------------------------------------------------------------- atlas
struct Atlas {
  std::vector<Root> recs;
  double cap = 0.125;
  // Try to merge a local Record; returns true if merged into an existing one.
  // Otherwise inserts a transformed (canonical gauge) copy.
  bool add(const Root& R);
  struct Impl;
  std::vector<int> bucketIdx;  // unused placeholder for hashing
  std::vector<std::vector<int>> buckets;
};

// distance bound for merge: min over variants Phi and wrap integers N of ||Pi(D0(Phi P) - D0(Q) - N)||
// (upper bound, interval evaluated).  Exposed for tests.
double mergeDistanceUpper(const double P[36], const double Q[36]);
// apply symmetry Phi (index map + sign) to a Root
struct Phi {
  int map[36];  // new[i] = old[map[i]]
  int sign = 1;
};
Root applyPhi(const Root& R, const Phi& f);

// binary Record I/O (SPEC section 7)
bool writeRecords(const char* path, const std::vector<Root>& v);
bool readRecords(const char* path, std::vector<Root>& v);

// ----------------------------------------------------------------- utilities
// Pi (double centering) applied to 36 intervals
void piApply(const I in[36], I out[36]);
// D0(P) wrapped difference: returns Delta = Pi(P0 - D0(P) - N) as intervals, N chosen nearest.
void alignDelta(const double P0[36], const double P[36], I Delta[36]);

}  // namespace sa
