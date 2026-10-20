# COVERAGE.md -- the seed cells cover the normal-form space

Claim (seed:cover): for every exact complex Hadamard matrix H, after the
symmetry reduction of seed:arrangement, its nine seed coordinates
q = (xi, other four x, y2..y5) lie within radius .0320001 (cycles, Euclidean)
of a center in `seeds.txt` that has the same selected position s.

## Argument
1. **Normal form exists** (seed:arrangement, implemented in `nf::arrange`): any
   exact H is equivalent (row/column permutation, rephasing, transposition) to
   a dephased matrix with `x=(0,x1<=..<=x5)`, `g_x=1-x5` the minimum gap over all
   30 row/column pairs, `x1 <= x5-x4`, a selected column moved second, rows 2..5
   sorted by it, `sum e(x_j)=0`, `sum e(y_j)=0`, `g_y >= g_x`.
2. **Bins.** Closed bins `[b/N,(b+1)/N]` cover [0,1] (both ends). Choose
   nondecreasing bins `bX=(a,b,c,d,l)` for x1..x5; s = first minimiser of
   `|2 b_s+1-N|` (a function of bX alone); the shared y entry uses bin `bX_s`
   exactly; the four tail bins are nondecreasing. The pair (bX, bY) is in the
   enumeration of `enumerate_x` x `build_tails` unless a filter rejects it.
3. **No filter rejects an exact point** (each is a proved necessary condition,
   evaluated as strict interval comparisons on enclosures of the exact quantities,
   tests/test_arith and tests/test_cover D): F1 `a <= l-d+1` is integral; F2 and the second gap
   inequality follow from seed:gap-stability and `g_y >= g_x`; `sumok` is
   seed:sum-support; `extend`/shared-interval is seed:clip-interval (the same
   grid displacement `h_s` lies in both intervals); a pair is dropped only when
   the intervals are strictly disjoint.
4. **Traversal induction.** Process surviving pairs in order. Each either
   (a) becomes a center (its midpoint; every point of its closed 9-dimensional
   bin is within `3/(2N) = .01339 < .032` of it), or (b) is certified by
   `test_center` against an already retained center u' with the same s: success
   requires `||k||^2 + 9/4 + B < RN2lo` with `B` either the box bound
   `sum|k_j|` (exact integers) or the interval upper bound of the dual bound (seed:dual-support, any v,
   `dual_cert_hi`), and `RN2lo <= (R N)^2` rigorously. By seed:distance-bound every feasible
   point of the pair then satisfies `N||q-u'|| < R N`. Retained centers are never
   removed, so by induction every feasible point of every surviving pair is
   within R = .032 of some retained center.
5. **Permutation.** Nine-coordinate order is the same for the center and the point (same s); distances
   are unchanged. Centers are the exact rationals `(b+1/2)/N`; a consumer that converts them to binary64
   moves each by < 1e-14, so radius `.0320001` covers that case too (the coverage test checks `.0320001`).

The argument depends on the generator only through "success is certified by the
checked inequality"; candidate order, the reach cut and the v-proposal are
unconstrained heuristics.

## Tests of the argument (`make test`, `make coverage`)
- `tests/test_cover.cpp` A: 2x100,000 random EXACT lists (N=112 and 64) pass
  sumok and extend with the true `h_s` inside both intervals (step 3).
- `tests/test_cover.cpp` C: 2x150,000 exact points: the true `N^2||q-u'||^2`
  is <= the box bound, <= the dual bound at the heuristic v, <= the dual
  bound at random v, and `test_center` success implies `< RN2` (step 4).
  A deliberately broken sign in `p_*` made the earlier (double) version of this test fail; the interval version keeps the same assertions.
- `tests/test_cover.cpp` D: at N=20 and N=24 every list is compared with a margin-free `__float128`
  re-implementation (same exact-double directions): 0 false rejections, 0 clearly-violated lists kept,
  `extend` bounds contain the quad interval (outward slack <= 1.5e-13), x-list counts equal.
- `tests/test_normal.cpp`: arrangement lemma conditions (equivalence of the
  arranged matrix, unitarity, sorted x/y, orientation, `g_x` = min gap,
  vanishing sums, `g_y >= g_x`, selection rule, bins contain the phases) on
  the Fourier family F6(a,b) (incl. the exact F6 with all tied pairs and bin
  boundaries) and 150 random exact Hadamards from Levenberg-Marquardt;
  alignment and transport identities.
- `make test` also runs an in-process N=48 cover without reuse
  (every pair a center) and checks 3,316 arranged samples are covered.
- `make coverage` (`tests/test_cover 112 seeds.txt`): generates the full N=112
  file (249,976 centers) and checks **13,716** arrangement samples (F6 grid,
  exact F6 with every minimum-gap pair, gap tie and bin-boundary convention,
  3000 random Hadamards x 4 choices): every filter passes, the pair is
  enumerated, and a same-s center lies within .0320001 (distance confirmed by interval comparison). Result: 0 failures.

Limits (untested, not claimed): the sample is finite (F6 family, LM-found Hadamards) and the sample points are
doubles (rounding ~1e-16, inside the 1e-7 slack of the tested radius .0320001); the proof is the argument above.
