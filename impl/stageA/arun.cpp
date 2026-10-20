// arun -- Stage A driver.
//   ./arun centers.txt i n out.bin bad.txt [--threads T] [--limit L] [-v]
// Seed file (SPEC sec 7 / A1): per line either "s b0 ... b8" (selected position + nine bin
// indices at N=112, as written by ../seed/seed) or nine decimal phases (+ optional s).
// Shard i processes lines with index = i (mod n).  Output: Record file (SPEC sec 7) and a text
// "bad" file listing every seed that failed setup or has unresolved cells.  Exit 0 iff the bad
// file is empty.
#include <omp.h>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include "stagea.hpp"
#include "../core/selftest.hpp"

using namespace sa;

static bool parseLine(const std::string& ln, double q[9]) {
  std::istringstream is(ln);
  std::vector<std::string> tok;
  std::string t;
  while (is >> t) tok.push_back(t);
  if (tok.size() < 9) return false;
  bool dec = false;
  for (auto& s : tok)
    if (s.find('.') != std::string::npos || s.find('e') != std::string::npos) dec = true;
  if (dec) {
    for (int i = 0; i < 9; ++i) q[i] = strtod(tok[i].c_str(), nullptr);
    return true;
  }
  if (tok.size() != 10) return false;
  for (int i = 0; i < 9; ++i) {
    long b = atol(tok[1 + i].c_str());
    q[i] = ((double)b + 0.5) / 112.0;  // division-rounded midpoint, N = 112
  }
  return true;
}

int main(int argc, char** argv) {
  iv::mub_startup_selftest("arun");
  if (argc < 6) {
    fprintf(stderr, "usage: arun centers.txt i n out.bin bad.txt [--threads T] [--limit L] [-v]\n");
    return 2;
  }
  const char* cfile = argv[1];
  int shard = atoi(argv[2]), nsh = atoi(argv[3]);
  const char* outf = argv[4];
  const char* badf = argv[5];
  int threads = omp_get_max_threads();
  long limit = -1;
  const char* dumpf = nullptr;
  std::vector<std::string> dumpLines;
  Params P;
  for (int i = 6; i < argc; ++i) {
    if (!strcmp(argv[i], "--threads") && i + 1 < argc) threads = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--limit") && i + 1 < argc) limit = atol(argv[++i]);
    else if (!strcmp(argv[i], "-v")) P.verbose = true;
    else if (!strcmp(argv[i], "--sample-excl") && i + 2 < argc) { P.sampleMod = (unsigned)atoi(argv[i + 1]); dumpf = argv[i + 2]; i += 2; }
    else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 2; }
  }
  std::ifstream in(cfile);
  if (!in) { fprintf(stderr, "cannot open %s\n", cfile); return 2; }
  std::vector<std::array<double, 9>> seeds;
  std::vector<long> ids;
  std::string ln;
  long idx = 0;
  while (std::getline(in, ln)) {
    if (ln.empty() || ln[0] == '#') continue;
    if (idx % nsh == shard) {
      std::array<double, 9> q;
      if (!parseLine(ln, q.data())) { fprintf(stderr, "bad line %ld\n", idx); return 2; }
      seeds.push_back(q);
      ids.push_back(idx);
      if (limit >= 0 && (long)seeds.size() >= limit) break;
    }
    ++idx;
  }
  { iv::RoundUp g; initTables(); }
  Atlas fin;
  fin.cap = P.cap;
  std::vector<std::string> bad;
  long nbad = 0, totUnres = 0;
  double tN = 0, tD = 0, tC = 0, tBu = 0;
  long cells = 0, covered = 0, exclN = 0, exclD = 0, localRecs = 0, mergedIn = 0;
  long hist[6] = {0, 0, 0, 0, 0, 0};  // local record radii: <.02 <.04 <.06 <.08 <.1 >=.1
  auto t0 = std::chrono::steady_clock::now();
  // Seeds are processed in parallel chunks (dynamic schedule); the atlas merge is then done serially in
  // seed-index order, so the output does not depend on the thread count or on scheduling.
  const long chunk = std::max<long>(256, 8L * threads);
  double secsMerge = 0;
  for (long c0 = 0; c0 < (long)seeds.size(); c0 += chunk) {
    long c1 = std::min<long>(c0 + chunk, (long)seeds.size());
    std::vector<SeedResult> res(c1 - c0);
#pragma omp parallel num_threads(threads)
    {
      iv::RoundUp g;
#pragma omp for schedule(dynamic, 1)
      for (long j = c0; j < c1; ++j) res[j - c0] = processSeed(seeds[j].data(), P);
    }
    auto m0 = std::chrono::steady_clock::now();
    {
      iv::RoundUp g;
      for (long j = c0; j < c1; ++j) {
        SeedResult& R = res[j - c0];
        long mg = 0;
        for (const Root& r : R.recs) {
          if (r.r < 0.02) ++hist[0]; else if (r.r < 0.04) ++hist[1]; else if (r.r < 0.06) ++hist[2]; else if (r.r < 0.08) ++hist[3]; else if (r.r < 0.1) ++hist[4]; else ++hist[5];
          mg += fin.add(r) ? 1 : 0;
        }
        localRecs += (long)R.recs.size();
        mergedIn += mg;
        tN += R.st.tNorm; tD += R.st.tDir; tC += R.st.tCover; tBu += R.st.tBuild;
        for (const ExclSample& e : R.samples) {
          char b[1024]; int n = snprintf(b, sizeof b, "%ld %d %d %d %d %d %d %a %a %a %a %d%d%d%d", ids[j], e.kind, e.Nu, (int)e.ia, (int)e.it, (int)e.ib, (int)e.ig, e.h, e.z2, e.Er, e.Rs, e.fc[0], e.fc[1], e.fc[2], e.fc[3]);
          if (e.kind == 1) for (int t = 0; t < 20; ++t) n += snprintf(b + n, sizeof b - n, " %a", e.y[t]);
          dumpLines.push_back(b);
        }
        cells += R.st.cellsNorm;
        covered += R.st.covered;
        exclN += R.st.exclNorm;
        exclD += R.st.exclDir;
        if (R.status != 0) {
          ++nbad;
          totUnres += R.unresolved;
          char buf[256];
          snprintf(buf, sizeof buf, "%ld status=%d unresolved=%ld %s", ids[j], R.status, R.unresolved, R.why.c_str());
          bad.push_back(buf);
        }
      }
    }
    secsMerge += std::chrono::duration<double>(std::chrono::steady_clock::now() - m0).count();
  }
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  fprintf(stderr, "serial atlas merge time %.2fs of %.2fs\n", secsMerge, secs);
  fprintf(stderr, "local record radii: <.02 %ld  <.04 %ld  <.06 %ld  <.08 %ld  <.1 %ld  >=.1 %ld\n", hist[0], hist[1], hist[2], hist[3], hist[4], hist[5]);
  if (!writeRecords(outf, fin.recs)) { fprintf(stderr, "cannot write %s\n", outf); return 2; }
  fprintf(stderr, "cpu-seconds: norm %.1f dir %.1f cover %.1f (of which root builds %.1f)\n", tN, tD, tC, tBu);
  if (dumpf) {
    std::ofstream df(dumpf);
    df << "# idx kind Nu ia it ib ig h z2 Er Rs fc [y20]  sampleMod=" << P.sampleMod << "\n";
    for (auto& l : dumpLines) df << l << "\n";
    fprintf(stderr, "sampled %zu excluded cells (1/%u) -> %s\n", dumpLines.size(), P.sampleMod, dumpf);
  }
  std::sort(bad.begin(), bad.end());
  { std::ofstream bf(badf); for (auto& s : bad) bf << s << "\n"; }
  fprintf(stderr, "seeds %zu  local-records %ld merged-into-existing %ld  final-records %zu  bad %ld  unresolved %ld  cells %ld norm-excl %ld dir-excl %ld covered %ld  wall %.2fs (%.2fs / 1000 seeds, %d threads)\n",
          seeds.size(), localRecs, mergedIn, fin.recs.size(), nbad, totUnres, cells, exclN, exclD, covered, secs,
          seeds.empty() ? 0.0 : secs * 1000.0 / seeds.size(), threads);
  return nbad ? 1 : 0;
}
