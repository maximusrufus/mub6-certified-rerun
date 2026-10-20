CXX      ?= g++
CXXFLAGS ?= -std=gnu++17 -O2 -ffp-contract=off -frounding-math -fopenmp -Wall -Wextra
# NEVER add -ffast-math / -Ofast / -ffp-contract=fast: the margins assume IEEE
# binary64 without contraction (SPEC sec 8).

HDRS  = arith.hpp cover.hpp reuse.hpp seedgen.hpp normal.hpp ../core/interval.hpp
TESTS = tests/test_arith tests/test_normal tests/test_cover

all: seed

seed: seed.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) seed.cpp -o seed

tests/test_arith: tests/test_arith.cpp $(HDRS) tests/check.hpp
	$(CXX) $(CXXFLAGS) tests/test_arith.cpp -o $@ -lquadmath

tests/test_normal: tests/test_normal.cpp $(HDRS) tests/check.hpp tests/hadamard.hpp
	$(CXX) $(CXXFLAGS) tests/test_normal.cpp -o $@

tests/test_cover: tests/test_cover.cpp $(HDRS) tests/check.hpp tests/hadamard.hpp
	$(CXX) $(CXXFLAGS) tests/test_cover.cpp -o $@ -lquadmath

# Fast unit tests (about a minute): arithmetic budgets, normal forms,
# filters on exact solutions, reuse-bound soundness, naive-vs-optimised counts,
# arrangement -> enumeration -> coverage at N=48.
test: $(TESTS)
	./tests/test_arith
	./tests/test_normal
	./tests/test_cover

# The full coverage self-check (about 4 minutes): generates the N=112 seed file
# and checks that exact Hadamard normal forms lie within .0320001 of a centre.
seeds.txt: seed
	./seed --out seeds.txt

coverage: seeds.txt tests/test_cover
	./tests/test_cover 112 seeds.txt

clean:
	rm -f seed $(TESTS)

.PHONY: all test coverage clean
