CXX      ?= g++
# NEVER add -ffast-math. -frounding-math and -ffp-contract=off are mandatory for outward rounding.
CXXFLAGS ?= -std=c++20 -O2 -frounding-math -ffp-contract=off -fopenmp -Wall -Wextra -DIV_CHECK_ROUNDING
PY       ?= python3
N        ?= 100000
SEED     ?= 20261019

HDRS = core/interval.hpp core/mub_ops.hpp

all: build/test_unit build/gen_records

build/test_unit: tests/test_unit.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -o $@ $<

build/gen_records: tests/gen_records.cpp $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -o $@ $<

test: all
	./build/test_unit
	./build/gen_records $(N) $(SEED) build/records.json
	$(PY) harness/flint_check.py build/records.json

clean:
	rm -rf build

.PHONY: all test clean
