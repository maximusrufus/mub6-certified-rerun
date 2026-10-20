// cover.hpp -- internal interface between the seed driver and the phase-ball covering.
#pragma once
#include <vector>
#include "stagea.hpp"

namespace sa {
struct CoverCtx {
  std::vector<Root> roots;  // seed-local Records (identity gauge)
  long nPrepFail = 0, nReuse = 0, nSnapOK = 0, nBuildFail = 0, nSnapBound = 0, nBoundTry = 0;
  const double* truth = nullptr;  // test hook: exact phases of a known completion
  bool cellHasWitness = false;
  long witChecks = 0, witViol = 0;
  double witWorstRatio = 0;       // max dist/s over checked bounds
  double tPrep = 0, tReuse = 0, tBuild = 0, tSnapBound = 0;
};
// Try to cover the (non-excluded) cell by a certified phase ball.  On success the
// radius of an existing local Root is grown or a new Root is appended.
bool tryCover(const SeedCtx& S, const Params& P, CellWork& w, CoverCtx& cc);
}  // namespace sa
