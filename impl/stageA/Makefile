CXX      ?= g++
# NEVER add -ffast-math. -frounding-math and -ffp-contract=off are mandatory for outward rounding.
BASE     ?= -std=c++20 -O2 -frounding-math -ffp-contract=off -fopenmp -Wall -Wextra
SRC      = seedproc.cpp root.cpp cover.cpp atlas.cpp
HDRS     = ext.hpp la.hpp stagea.hpp cover.hpp ../core/interval.hpp
THREADS ?= 4
PY      ?= python3

all: build/arun build/arun_chk build/test_stagea

# production driver: no per-operation rounding-mode assertion (2.6x faster)
build/arun: arun.cpp $(SRC) $(HDRS)
	@mkdir -p build
	$(CXX) $(BASE) -o $@ arun.cpp $(SRC)

# checked driver: every interval operation asserts FE_UPWARD (catches a missing RoundUp guard)
build/arun_chk: arun.cpp $(SRC) $(HDRS)
	@mkdir -p build
	$(CXX) $(BASE) -DIV_CHECK_ROUNDING -o $@ arun.cpp $(SRC)

build/test_stagea: tests/test_stagea.cpp $(SRC) $(HDRS)
	@mkdir -p build
	$(CXX) $(BASE) -DIV_CHECK_ROUNDING -o $@ tests/test_stagea.cpp $(SRC)

# make test:
#   1. unit + soundness tests (witness cells never excluded, covers contain the truth, far-off seeds excluded)
#   2. checked driver on tests/sample_seeds.txt (real seeds, SPEC A1 format): no bad/unresolved seeds
#   3. production driver on the same sample must give a byte-identical record file (checked build, 1 thread vs N threads)
#   4. every record re-certified from scratch (C++ checkRoot) and cross-checked with python-flint
test: build/test_stagea build/arun build/arun_chk
	FUZZ_F6=20 FUZZ_OTHER=4 ./build/test_stagea 2 8
	./build/arun_chk tests/sample_seeds.txt 0 1 build/sample_chk.bin build/sample.bad --threads $(THREADS)
	test ! -s build/sample.bad
	./build/arun tests/sample_seeds.txt 0 1 build/sample.bin build/sample.bad --threads $(THREADS)
	test ! -s build/sample.bad
	cmp build/sample.bin build/sample_chk.bin
	./build/arun tests/sample_seeds.txt 0 1 build/sample_1t.bin build/sample.bad --threads 1
	cmp build/sample.bin build/sample_1t.bin
	./build/test_stagea --verify build/sample.bin
	$(PY) tools/verify_records.py build/sample.bin 20

clean:
	rm -rf build

.PHONY: all test clean
