#!/bin/bash
# usage: build.sh target.cpp out
cd "$(dirname "$0")"
g++ -std=c++20 -O2 -frounding-math -ffp-contract=off -fopenmp -Wall -Wextra -DIV_CHECK_ROUNDING -o "$2" "$1"
