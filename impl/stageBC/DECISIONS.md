# DECISIONS.md — Stage B / Stage C re-implementation (MUB d=6)

Scope: `brun` (Stage B: relation iteration, rho search, tangent grid, markRec, assign, batching, the three (tau,w) passes) and
`crun` (Stage C: bin tests, clustering, O/U edge tests, joint test, mode-1 two-disjoint-6-clique search with the 5M-step cap),
plus `run_pipeline.sh` (mid -> fine+batches -> C per batch -> fallback -> final C). Spec: `../../SPEC.md` §4–§9, LaTeX in `../../src/build/sections/`.
Interval layer: `../core/interval.hpp` (outward-rounded binary64, FE_UPWARD negation trick). Build/test inside WSL: `make test`.

## 0. Invariant
Every exclusion decision is a certified inequality evaluated in interval arithmetic from the exact stored doubles. Anything not
certified is carried forward: Stage B failed parents -> `bad.bin` (exit 1); Stage C unresolved records (bin overflow, cluster failure,
graph step cap, graph candidate, any exception) -> `badout.bin` (exit 1). The pipeline ends with `RESULT: 0 unresolved` or an explicit
residual list (`RESIDUAL_*.bin`). The paper's counts (600, 22,379, 163, 1127, 282/282/282/281) are NOT used anywhere as acceptance.

## 1. A23 (S): expair / operation orderings -> interval evaluation, ordering irrelevant
The paper's additive margins (1e-12, 2e-7, 3e-7, .0001, 1e-8 ... in exclusion.tex §5.x and rounding.tex) exist only to absorb binary64 error
of one specific operation order. This code does not reproduce that order and therefore uses NO such margin:
* every quantity on either side of an inequality is an `Interval` (or, in Stage B hot paths, a plain double computed under FE_UPWARD from
  non-negative operands with monotone operations, which is an upper bound by induction); a test excludes iff `lhs.lo > rhs.hi`
  (strict, so equality and overlap retain the candidate);
* `e(t)` is `iv::cis_turns` (Taylor with remainder, enclosure of the true value), so no `expair` accuracy premise is used (A23);
* constants `pi, sqrt6, 1/sqrt6, 1/6, K=1.00001, sqrt(29/6)` are enclosures (`common.hpp: cn()`), never rounded doubles;
* soundness therefore does not depend on `-ffp-contract`, `pow(x,2)` lowering, FLT_EVAL_METHOD or operation order, only on the interval layer
  (which requires `-frounding-math`, no `-ffast-math`, and a `RoundUp` guard in EVERY thread: `crun` constructs one inside the OpenMP region;
  test binaries are built with `-DIV_CHECK_ROUNDING`, which asserts FE_UPWARD on every interval operation).
* Decimal command-line parameters (`tau`, `w`, `0.03`) are parsed BEFORE the guard (strtod honours the rounding mode: under FE_UPWARD `0.03`
  parses to the next double). `printf("%.3f")` inside a guard also rounds upward: cosmetic only.
* Heuristics (Newton snapping, Jacobi eigen-decompositions, linear solves for directions y, Lloyd clustering, hash neighbours, vertex order,
  integer wraps via `floor(x+.5)`) run in plain doubles and are never trusted: validity rests on the final checked inequality
  (introduction.tex L141–148). `rnd()` replaces `nearbyint` because the latter is mode dependent.

Consequence: the paper's "K" operator bound uses `K = 1.00001` as an interval; Stage C additionally requires the stored `froot < 1e-6`, `ort < 1e-5`
(validation) which is what `exclude:operator-bound` needs.

## 2. Stage B choices
* Input validation exactly SPEC 4.1 (ranges, k<=10, |p|<=4, |T|<=1.01, tau,w ranges, <=2M records/file, <=8M combined); invalid input => exit 2.
* `Rel` (relation parent->Q): `d_q = Pi(D0Q - D0P - N)` with N = `floor(x+.5)` per entry, all in intervals; `eta_q`, `n_q`, `sigma`
  (via `spectralUpper` = min(abs row sum, sqrt||G^2||_F), square-rooted), `D, v0, lambda, e0, beta` as upper bounds; `eval(t)` is the iteration
  of refine:iteration-* (<=30 improvements, returns 10 if t not in [0,2]). The `1e-8` enlargements of the paper are replaced by upward rounding.
* rho search = refine:rho-trial/return verbatim (25 midpoint trials, `D(0) > a-1e-7` reject, `D(m) < a-2e-7`, `rho = max(0,l-1e-8)`, usable iff `> 1e-7`).
  Cuts on `||d_q||`: .249 (children/direct) and .23 (blockers), +1e-6 allowance (A15 reading: norm of d_q).
* Grid: `E = fl(sqrt(1+ort)R)+1e-7` (upper bound), `L = ceil(E/w+.500001)` then verified `(L+1/2)w >= E` rigorously; <=50M boxes, depth <= 6.
  `outside` uses lower bounds of `max(0,|u|-h)` per coordinate vs `E^2` (upper). Corner bound `t_B` = interval norm upper bound.
* `markRec`: DFS over per-axis tables of `(|b w - eta_q| + w/2 + 1e-8)^2` (upper), threshold `(rho-5e-8)^2`, marked bound `sqrt(sum)+2e-8`, and a mark
  is only committed when that bound is `<= rho` (checked explicitly, independent of the pruning threshold). Test `mark` verifies soundness
  (every marked box has farthest corner <= rho) and completeness against brute force.
* Atlas children: radius growth to `Dhat(t_max)+1e-8` must stay `<= tau`, else the relation is simply not used (the paper aborts; we only decline).
  `assign` = merge if `r+l+2e-8 < tau` (radius `max(old, r+l+1e-8)`), else insert; l = interval upper bound of the aligned distance.
* Direct proposals: Newton snap of `P + T u` restricted to the horizontal complement of T (keeps eta_q = u; snapping inside all 25 directions
  drifts along near-singular tangent directions and was measured to give t_B ~ .08 at |u| = .05), wrap entries to [-1/2,1/2], `certifyRoot`
  (interval re-verification of ort, hor, froot, mnorm, n0, C, id), relation (cut .249), `r = Dhat(t_B)+2e-8 < tau-1e-7`.
* A9 (Root construction): eigen-decomposition of the 25x25 horizontal Gram `(L_P B)^T (L_P B)`; eigenvalues `< .64` -> T (more than 10 => fail);
  `H_V = B V Lambda^{-1/2}`, `L_V = L_P B V Lambda^{-1/2}`, `M = H_V L_V^T`; `mnorm = sqrt(ah*al)` from certified abs-row-sums of the column Grams;
  `C = pi*sqrt(8)*mnorm`; `id = ||Pi - TT^T - H_V (L_V^T L_P) Pi||_F` on the full 36x36 matrix; all in intervals. Acceptance: id<1e-4, ort,hor<1e-5,
  C<20, mnorm<2, froot<1e-9, n0<1e-6. `C = min(pi sqrt8 mnorm, tensor)` with the tensor improvement of round:tensor-c in
  `../core/tensor_c.hpp` (interval congruence certificate, see stageA/DECISIONS.md 3); 6% effect on real records; `MUB_NO_TENSOR_C=1` disables.
* Symmetry maps Phi (`phi.hpp`, implemented 2026-10-19; A13 signature rule copied from stageA/atlas.cpp, identity always included):
  - parent frames: `Grid::PV[v] = Phi_v P` for the parent's variants; blockers and atlas children are related from EVERY frame
    (refine:mapped-coordinates: grid indices are unchanged because `(Phi T)^T (Phi delta) = T^T delta`); relation cache keyed by (v, child);
  - `assign`: variants of the proposed child Q are tried against hash neighbours; merge = certified aligned distance of `Phi Q` from Z
    (refine:atlas-assignment); otherwise the CANONICAL `Phi Q` (first signature map) is inserted, re-certified from scratch in the stored
    orientation (fresh constants; transported ones kept if that fails);
  - after a direct proposal the stored child is pulled back by `Phi^{-1}` (exact signed permutation) for the parent's own relation
    (the paper's inverseTransform/addNear).
  Every mark is still a `makeRel` + `markRec` certificate on a transported Root (phase:transport); a Phi never covers a box by itself.
  Measured effect (../pipeline/PERF_FIX.md): mid roots for the first 3 real Stage A records 587 -> 18, for 20 records 3153 -> see file.
  `MUB_NO_PHI=1` disables (measurement switch). Tests: `phi` (inverse, isometry, Pi-commutation, residual invariance, transported
  root-normal inequality on exact solutions, relation isometry), `stageB` cover check through exactly the maps used (`usedMaps`),
  and the same parent without Phi.
* Parent order = input order (A14). Blocker list = earlier successful parents, queried at their fixed input radius. Children created during a failed
  parent remain in the atlas. Box-level reuse (`tryReuse`) and `markByChild` are caches of relations parent->child (relations do not depend on the box).
* Outputs: `out.bin` = atlas in insertion order with achieved radii; with a batch prefix, `f{j}.bin` (<=256 roots, written radius = tau, temp name +
  atomic rename) and `fdone.txt` ("<#batches> <#roots>") are written at the END of the run (the paper may emit them while running; final state identical).
  Exit 0 iff no failed parent; failed parents are written to `bad.bin`.

## 3. Stage C choices
* Bin tests in the paper's order: direct (exclude:direct-tests), scalar-linear, directional support with proposed y. Several candidate y per bin
  (L2 surrogate + 4 IRLS passes + residual direction) instead of "one proposed direction" (A8/A19): any y is sound, more candidates only raise power.
* Joint test: candidate y from a 14x14 (13x13 for U) regularised surrogate with reweighting; ball error formed from the pure single data BEFORE adding edge
  derivatives (as exclusion.tex requires); tangent term shared and accumulated over both balls.
* Clustering: two passes (second does not move centres), threshold `(.042-1.11/M)^2`, radius = interval upper bound of `||h|| + X_M` per assigned bin.
  Failure (more than 2048 centres, |v_j|>2, radius >= .1) => UNRESOLVED.
* Graph search: exact bitset implementation of mode 1. A18: one "step" = one recursive call of `rec` or `clique` (including calls from the 6-clique prune).
  Exhaustion sets a persistent give-up flag and is reported as UNRESOLVED(graph-cap), never negative. ADDITION (exact, not in the paper): before the
  search, vertices with O-degree < 5 or U-degree < 6 are deleted iteratively; such a vertex lies in no configuration, so this cannot change a negative
  answer (tested against brute force on 400 random graphs, with and without pruning).
* Outcomes: NEGATIVE(empty-layer) / NEGATIVE(graph) are the only exclusions. UNRESOLVED(bin-overflow | cluster-failure | graph-cap | graph-candidate).
* `--triple` (test hook): looks for ONE O 6-clique (a single further basis). F6 must keep it ({I,F6} extends to a third basis; SPEC §10.6).

## 4. Tests (`make test`; all pass in WSL g++ 15)
* `test_bc` (names: io root mark rel graph cluster single edge taylor triple soundness stageB; `stageBslow` via `make test-slow`):
  - root: certificates on the Fourier family; root-normal/root-square inequalities hold on exact family solutions;
  - rel: Dhat >= true ||delta - d_q|| on ~1.6k exact solutions; rho certificate checked;
  - stageB: every output root revalidated, stored froot >= interval norm, and all sampled exact family solutions inside the parent ball are covered by a child ball;
  - single/edge: all 36–39 exact unbiased vectors of an exact Hadamard in the ball survive all four layers; all exact orthogonal pairs keep their O edge (incl. joint test);
  - taylor: 3200 checks of the analytic lemmas on arbitrary (non-solution) points, including worst-case displacement directions (sign convention of the joint test;
    a deliberate sign flip makes this test fail - verified by mutation);
  - triple: F6 keeps a third basis; two further bases on F6 (r=.01) are certified NEGATIVE (consistent with Grassl / the theorem); capped runs are UNRESOLVED;
  - soundness: overflow path, tiny cluster radius, oversize grid all reported as UNRESOLVED / failed parent.
* `xgen` + `tests/xcheck_bc.py`: randomised cases re-derived independently with python-flint (arb/acb, 220 bits); C++ enclosures must contain the flint balls
  (exact rational test, `harness/flint_check.py: contains` when importable). Stage B double bounds must be >= flint upper bounds. Mutated expectations are detected.
* `tests/pipeline_test.sh`: synthetic Stage-A-like records -> full pipeline with a small initial graph cap to force the fallback path; `pipeline_check.py`
  verifies file formats, batch size/radius/marker, batches == fine.bin up to radius, bad files are ordered subsequences.

## 5. Not done / untested (keep as a to-do, these are NOT "gated")
* No real Stage A output exists yet: all tests use certified Roots at Fourier-family (and Bjorck) Hadamards. Behaviour on the real 3083 records (their k, C, mnorm,
  radii up to .125) is UNTESTED; the success of the three B passes on them is unknown (C without the tensor improvement and no Phi maps are the likely risks).
* U-edge joint soundness on exact unbiased PAIRS is covered only through the Taylor test (no exact mutually-unbiased pair among vectors unbiased to I and F6 exists,
  as expected: that would be a 4th MUB); no exact configuration of that kind can exist.
* Performance: Stage B direct proposals cost ~0.6 ms each (certifyRoot dominates); Stage B is single-threaded (the atlas is sequential). Stage C is OpenMP over records.
  Full-scale runtime unmeasured.
* No paid API, no network; nothing pushed.
