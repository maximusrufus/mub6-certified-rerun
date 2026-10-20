// crun: Stage C driver.
//   ./crun M layers dirs radius arg5 mode graphcap shard parts badout.bin input.bin
// (SPEC 5; the paper's invocation is `./crun 12 4 1 .042 1 1 5000000 i parts badout.bin input.bin`).
// Processes records with index = shard (mod parts).  Exit 0 iff every processed record is certified NEGATIVE;
// otherwise exit 1 and badout.bin holds the unresolved records in input order (an ordered subsequence of the input).
// Exit 2 on invalid input.  Extra flag --triple as 12th argument: sanity mode (look for one O-6-clique).
#include "stageC.hpp"
#include "../core/selftest.hpp"
#include <omp.h>
#include <chrono>

using namespace mb;

int main(int argc, char** argv) {
  iv::mub_startup_selftest("crun");
  if (argc < 12) {
    std::fprintf(stderr, "usage: crun M layers dirs radius arg5 mode graphcap shard parts badout.bin input.bin\n");
    return 2;
  }
  CParams P;
  P.M0 = std::atoi(argv[1]);
  P.layers = std::atoi(argv[2]);
  P.radiusParam = std::atof(argv[4]);
  int mode = std::atoi(argv[6]);
  P.graphCap = (uint64_t)std::atoll(argv[7]);
  int shard = std::atoi(argv[8]), parts = std::atoi(argv[9]);
  std::string badp = argv[10], inp = argv[11];
  if (argc > 12 && std::string(argv[12]) == "--triple") P.tripleMode = true;
  if (mode != 1 || P.M0 < 1 || P.M0 > 100 || P.layers < 1 || parts < 1 || shard < 0 || shard >= parts) {
    std::fprintf(stderr, "unsupported parameters (only graph mode 1 is implemented)\n");
    return 2;
  }
  if (P.M0 << (P.layers - 1) > 255) { std::fprintf(stderr, "final resolution exceeds 255\n"); return 2; }
  std::vector<Record> in;
  std::string err;
  if (!readRecords(inp, in, err)) { std::fprintf(stderr, "input error: %s\n", err.c_str()); return 2; }
  std::vector<size_t> mine;
  for (size_t i = 0; i < in.size(); ++i) {
    std::string why;
    if (!validateRecord(in[i], 0.03, why)) { std::fprintf(stderr, "record %zu invalid: %s\n", i, why.c_str()); return 2; }
    if ((long)(i % (size_t)parts) == shard) mine.push_back(i);
  }
  cn();  // construct constants before the parallel region
  std::vector<COutcome> out(mine.size());
  std::vector<CStats> stats(mine.size());
  auto t0 = std::chrono::steady_clock::now();
#pragma omp parallel
  {
    iv::RoundUp g;  // every thread needs its own upward-rounding guard
#pragma omp for schedule(dynamic, 1)
    for (long t = 0; t < (long)mine.size(); ++t) {
      try {
        out[t] = processRecordC(in[mine[t]], P, stats[t]);
      } catch (const std::exception& e) {
        out[t] = COutcome::ClusterFail;  // never an exclusion
        stats[t].note = std::string("exception: ") + e.what();
      }
    }
  }
  size_t nneg = 0;
  std::vector<const Record*> bad;
  for (size_t t = 0; t < mine.size(); ++t) {
    bool neg = isNegative(out[t]);
    if (P.tripleMode) neg = (out[t] == COutcome::TripleAbsent);
    if (neg) ++nneg; else bad.push_back(&in[mine[t]]);
    std::fprintf(stderr, "record %zu: %s balls=%zu edgesO=%zu edgesU=%zu pruned=%zu steps=%llu layers=", mine[t], outcomeName(out[t]), stats[t].balls,
                 stats[t].edgesO, stats[t].edgesU, stats[t].pruned, (unsigned long long)stats[t].steps);
    for (size_t b : stats[t].layerBins) std::fprintf(stderr, "%zu,", b);
    std::fprintf(stderr, " %s\n", stats[t].note.c_str());
  }
  if (!writeRecords(badp, bad)) { std::fprintf(stderr, "cannot write %s\n", badp.c_str()); return 2; }
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::fprintf(stderr, "[C] shard %d/%d: %zu records, %zu negative, %zu unresolved, %.2fs\n", shard, parts, mine.size(), nneg, bad.size(), secs);
  std::printf("CRUN processed=%zu negative=%zu bad=%zu\n", mine.size(), nneg, bad.size());
  return bad.empty() ? 0 : 1;
}
