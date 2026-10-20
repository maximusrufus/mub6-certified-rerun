#!/bin/bash
# Smoke test from a clean clone (about 10-15 minutes on 8 cores).  Exit 0 iff everything passes.
# Requirements: g++ (GCC >= 12, OpenMP), GNU make, python3 with python-flint and numpy, bc, sha256sum.
# Steps: (1) core interval library unit tests + python-flint cross-check; (2) seed-stage unit tests (arithmetic
# budgets, normal forms, filters on exact Hadamards, reuse-bound soundness, N=48 cover);
# (3) Stage A tests + a 25-seed run (checked build vs plain build vs single thread, byte-identical);
# (4) Stage B/C unit tests + the synthetic end-to-end B/C pipeline test; (5) the independent coverage checker on
# the first 200 X-lists (must print COVERAGE-PROVEN) and the hash of the committed seed file.
# The full 1% end-to-end sample (seed regeneration, Stage A on every 100th seed, S5/S6 verification, B/C with a
# 600 s budget) is `bash impl/pipeline/smoke_test.sh` (about 30 minutes on 16 threads); the full run is
# `THREADS=16 bash impl/pipeline/run_all.sh` (12.4 h on 16 threads).
set -eu
cd "$(dirname "$0")"
T=${THREADS:-8}
export OMP_NUM_THREADS=$T
echo "== [1/5] core interval library"
( cd impl && make -s test N=20000 )
echo "== [2/5] seed stage"
( cd impl/seed && make -s test )
echo "== [3/5] stage A"
( cd impl/stageA && make -s test THREADS=$T )
echo "== [4/5] stage B/C"
( cd impl/stageBC && make -s all && ./build/test_bc && bash tests/pipeline_test.sh | tail -2 )
echo "== [5/5] coverage checker (first 200 X-lists) + seed file hash"
( cd coverage && make -s covcheck && ./covcheck ../impl/seed/seeds.txt --limit-x 200 2>/dev/null | tail -1 | grep -q '^COVERAGE-PROVEN$' && echo COVERAGE-PROVEN )
echo "22cd76ed4e211e842cbf03530ed46fba3416fb9f704f0342d020d4e838d46311  impl/seed/seeds.txt" | sha256sum -c -
echo "SMOKE OK"
