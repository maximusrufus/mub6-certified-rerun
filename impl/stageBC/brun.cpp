// brun: Stage B driver.   ./brun tau w out.bin bad.bin batchprefix|- inputs...
// Exit 0 iff every parent was fully refined (bad.bin then holds the 8-byte zero count); exit 1 if some parents
// failed (they are written to bad.bin: the residual obligation); exit 2 on invalid input / I/O error.
#include "stageB.hpp"
#include "../core/selftest.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>

using namespace mb;

int main(int argc, char** argv) {
  iv::mub_startup_selftest("brun");
  if (argc < 7) {
    std::fprintf(stderr, "usage: brun tau w out.bin bad.bin batchprefix|- inputs...\n");
    return 2;
  }
  double tau = std::atof(argv[1]), w = std::atof(argv[2]);  // parsed in round-to-nearest, before the guard
  iv::RoundUp guard;
  std::string outp = argv[3], badp = argv[4], prefix = argv[5];
  if (!(tau >= 0.0001 && tau <= 0.126) || !(w >= 0.0001 && w <= 0.05)) {
    std::fprintf(stderr, "parameters out of range\n");
    return 2;
  }
  std::vector<Record> in;
  for (int a = 6; a < argc; ++a) {
    std::string err;
    if (!readRecords(argv[a], in, err)) { std::fprintf(stderr, "input error: %s\n", err.c_str()); return 2; }
  }
  if (in.size() > 8000000) { std::fprintf(stderr, "combined input exceeds 8,000,000\n"); return 2; }
  for (size_t i = 0; i < in.size(); ++i) {
    std::string why;
    if (!validateRecord(in[i], 0.126, why)) { std::fprintf(stderr, "record %zu invalid: %s\n", i, why.c_str()); return 2; }
  }
  auto t0 = std::chrono::steady_clock::now();
  Refiner R(tau, w);
  R.verbose = std::getenv("BDEBUG") != nullptr;
  std::vector<const Record*> bad;
  for (size_t i = 0; i < in.size(); ++i) {
    bool ok = R.processParent(in[i]);
    if (!ok) {
      bad.push_back(&in[i]);
      std::fprintf(stderr, "parent %zu FAILED: %s\n", i, R.lastFail.c_str());
    }
    if (R.atlas.recs.size() > 2000000) { std::fprintf(stderr, "atlas exceeds 2,000,000 roots\n"); return 2; }
    if (((i + 1) % 10 == 0 || i + 1 == in.size() || std::getenv("BPROGRESS")))
      std::fprintf(stderr, "[B] %zu/%zu parents, atlas %zu, blockers %zu, bad %zu  (cum: boxes %ld outside %ld marked %ld direct %ld merged %ld inserted %ld, %.0fs)\n", i + 1, in.size(), R.atlas.recs.size(), R.blockers.size(), bad.size(), R.st.boxes, R.st.outside, R.st.marked, R.st.direct, R.st.merged, R.st.inserted, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  }
  // outputs
  std::vector<const Record*> outv;
  for (auto& r : R.atlas.recs) outv.push_back(&r);
  if (!writeRecords(outp, outv)) { std::fprintf(stderr, "cannot write %s\n", outp.c_str()); return 2; }
  if (!writeRecords(badp, bad)) { std::fprintf(stderr, "cannot write %s\n", badp.c_str()); return 2; }
  if (prefix != "-") {
    // batches of <=256 roots in insertion order; written radius = tau; atomic rename; marker fdone.txt
    size_t nb = 0;
    std::vector<Record> tmp;
    for (size_t s = 0; s < R.atlas.recs.size(); s += 256) {
      tmp.clear();
      for (size_t i = s; i < std::min(R.atlas.recs.size(), s + 256); ++i) {
        Record r = R.atlas.recs[i];
        r.r = tau;
        tmp.push_back(std::move(r));
      }
      std::string name = prefix + "f" + std::to_string(nb) + ".bin";
      if (!writeRecords(name + ".tmp", tmp)) { std::fprintf(stderr, "cannot write batch\n"); return 2; }
      std::filesystem::rename(name + ".tmp", name);
      ++nb;
    }
    std::ofstream m(prefix + "fdone.txt.tmp");
    m << nb << " " << R.atlas.recs.size() << "\n";
    m.close();
    std::filesystem::rename(prefix + "fdone.txt.tmp", prefix + "fdone.txt");
  }
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::fprintf(stderr,
               "[B] done: %zu parents -> %zu roots (merged %ld, inserted %ld), failed parents %zu, boxes %ld outside %ld marked %ld reused %ld direct %ld/%ld fail, %.2fs\n",
               in.size(), R.atlas.recs.size(), R.st.merged, R.st.inserted, bad.size(), R.st.boxes, R.st.outside, R.st.marked, R.st.reused, R.st.direct, R.st.directFail, secs);
  std::printf("BRUN in=%zu out=%zu bad=%zu\n", in.size(), R.atlas.recs.size(), bad.size());
  return bad.empty() ? 0 : 1;
}
