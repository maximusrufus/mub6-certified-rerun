// seed.cpp -- seed cover driver ("covrun" re-implementation), SPEC sec 2.
//
//   ./seed --count-only            prints the two deterministic checkpoints
//                                  (x-lists, pairs) on stdout, one per line
//   ./seed --out seeds.txt         full cover; one center per line (SPEC A1 (b))
//   options: --N n (default 112), --cut c (reuse reach, integer ||k||^2),
//            --maxx m (debug: only the first m x-lists)
//
// Output line: "s b0 b1 b2 b3 b4 b5 b6 b7 b8" -- selected position s (0..4) and
// the nine bin indices in nine-coordinate order (xi, other four x, y2..y5).
// Center phases are (b+1/2)/N (division-rounded in binary64).
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <string>
#include <omp.h>
#include "../core/selftest.hpp"
#include "cover.hpp"
#include "reuse.hpp"
#include "seedgen.hpp"

using namespace seedcov;

int main(int argc, char **argv) {
  iv::mub_startup_selftest("seed");
  Options o;
  const char *out = "seeds.txt";
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--count-only") o.count_only = true;
    else if (a == "--N" && i + 1 < argc) o.N = atoi(argv[++i]);
    else if (a == "--cut" && i + 1 < argc) o.cut = atoi(argv[++i]);
    else if (a == "--maxx" && i + 1 < argc) o.maxx = atol(argv[++i]);
    else if (a == "--out" && i + 1 < argc) out = argv[++i];
    else if (a == "--quiet") o.verbose = false;
    else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 2; }
  }
  Tables T(o.N);   // builds the interval tables under its own rounding guard
  auto t0 = std::chrono::steady_clock::now();
  auto el = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };

  iv::RoundUp mainGuard;   // main thread evaluates certificates too (workers hold their own)
  std::vector<XList> X = enumerate_x(T);
  if (o.maxx > 0 && (long)X.size() > o.maxx) X.resize(o.maxx);
  if (o.verbose) fprintf(stderr, "x-lists %zu (%.1fs)\n", X.size(), el());
  TailTable tails = build_tails(T, X);
  if (o.verbose) fprintf(stderr, "tails built (%.1fs)\n", el());

  if (o.count_only) {
    PairCounts pc = count_pairs(T, X, tails);
    if (o.verbose)
      fprintf(stderr, "pairs: sumok-only %ld, +gap2 %ld, +extend %ld (%.1fs)\n", pc.sumok_only, pc.gap2, pc.extend, el());
    printf("%zu\n%ld\n", X.size(), pc.extend);
    return 0;
  }

  CoverResult R = run_cover(T, X, tails, o);
  if (o.verbose)
    fprintf(stderr, "pairs %ld  centers %zu  box %ld dual %ld cache %ld (%.1fs)\n", R.pairs, R.centers.size(), R.n_box,
            R.n_dual, R.n_cache, el());
  FILE *f = fopen(out, "w");
  if (!f) { perror("open out"); return 1; }
  for (const Center &c : R.centers) {
    int b[9]; center_bins(X[c.xi], c.tail, b);
    fprintf(f, "%d", X[c.xi].s);
    for (int j = 0; j < 9; ++j) fprintf(f, " %d", b[j]);
    fputc('\n', f);
  }
  fclose(f);
  printf("%zu\n", R.centers.size());
  return 0;
}
