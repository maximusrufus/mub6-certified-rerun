# Stage B / Stage C (MUB d=6 re-implementation).  Build/test inside WSL:  make test
# NEVER add -ffast-math.  -frounding-math and -ffp-contract=off are mandatory for the outward-rounded interval layer.
CXX      ?= g++
PY       ?= python3
BASE     = -std=c++20 -O2 -frounding-math -ffp-contract=off -fopenmp -Wall -Wextra
# test builds check at runtime that EVERY interval operation executes under FE_UPWARD
TESTFLAGS = $(BASE) -DIV_CHECK_ROUNDING
PRODFLAGS = $(BASE) -DNDEBUG

HDRS = ../core/interval.hpp common.hpp root.hpp stageB.hpp stageC.hpp graph.hpp

all: build/brun build/crun build/mkrec build/test_bc build/xgen

build/brun: brun.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(TESTFLAGS) -o $@ $<
build/crun: crun.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(TESTFLAGS) -o $@ $<
build/mkrec: tests/mkrec.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(TESTFLAGS) -o $@ $<
build/test_bc: tests/test_bc.cpp tests/fixtures.hpp $(HDRS)
	@mkdir -p build
	$(CXX) $(TESTFLAGS) -o $@ $<
build/xgen: tests/xgen.cpp tests/fixtures.hpp $(HDRS)
	@mkdir -p build
	$(CXX) $(TESTFLAGS) -o $@ $<

# production binaries (no per-operation rounding-mode assertion)
prod: build/prod/brun build/prod/crun
build/prod/brun: brun.cpp $(HDRS)
	@mkdir -p build/prod
	$(CXX) $(PRODFLAGS) -o $@ $<
build/prod/crun: crun.cpp $(HDRS)
	@mkdir -p build/prod
	$(CXX) $(PRODFLAGS) -o $@ $<

N    ?= 120
SEED ?= 20261019

test: all
	./build/test_bc
	./build/xgen $(N) $(SEED) build/xcases.json
	$(PY) tests/xcheck_bc.py build/xcases.json
	bash tests/pipeline_test.sh

# the slow Stage B end-to-end cover tests (large parents); not part of `make test`
test-slow: all
	./build/test_bc stageBslow

clean:
	rm -rf build

.PHONY: all prod test test-slow clean
