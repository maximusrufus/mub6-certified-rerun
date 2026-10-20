CXX      ?= g++
# -ffp-contract=off is mandatory (the outward widening assumes separately rounded + - * /); never -ffast-math.
CXXFLAGS ?= -std=c++17 -O2 -ffp-contract=off -fopenmp -Wall -Wextra

covcheck: covcheck.cpp roots_112.hpp
	$(CXX) $(CXXFLAGS) covcheck.cpp -o covcheck

roots_112.hpp: gen_roots.py
	python3 gen_roots.py 112 roots_112.hpp

check: covcheck
	./covcheck ../impl/seed/seeds.txt | tee COVERAGE_RUN.out

clean:
	rm -f covcheck
.PHONY: check clean
