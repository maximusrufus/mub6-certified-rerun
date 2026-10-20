// Fixture generator: certified Roots (Records) at known Hadamard matrices.
//   mkrec fourier a b radius out.bin     (Fourier family F6(a,b), a,b in radians, SPEC-format record)
// Also provides fourierPhases() for the unit tests.
#include "../root.hpp"
#include <cstdio>
using namespace mb;

static void fourierPhases(double a, double b, double* P) {
  const double TWO_PI = 6.283185307179586476925286766559;
  for (int j = 0; j < 6; ++j)
    for (int k = 0; k < 6; ++k) {
      double R = 0;
      if (j % 2 == 1 && (k % 3) == 1) R = a;
      if (j % 2 == 1 && (k % 3) == 2) R = b;
      P[6 * j + k] = (double)((j * k) % 6) / 6.0 + R / TWO_PI;
    }
}

int main(int argc, char** argv) {
  if (argc < 6) { std::fprintf(stderr, "usage: mkrec fourier a b radius out.bin\n"); return 2; }
  double a = std::atof(argv[2]), b = std::atof(argv[3]), r = std::atof(argv[4]);
  iv::RoundUp g;
  double P[36];
  fourierPhases(a, b, P);
  if (!snapSolution(P)) { std::fprintf(stderr, "snap failed\n"); return 1; }
  Record R;
  RootInfo info;
  if (!certifyRoot(P, R, &info)) { std::fprintf(stderr, "certify failed: %s\n", info.why.c_str()); return 1; }
  R.r = r;
  std::printf("k=%d C=%.4g n0=%.3g id=%.3g ort=%.3g hor=%.3g mnorm=%.4g froot=%.3g minEig=%.3g\n", R.k, R.C, R.n0, R.id, R.ort, R.hor, R.mnorm, R.froot, info.minEig);
  return writeRecords(argv[5], std::vector<Record>{R}) ? 0 : 1;
}
