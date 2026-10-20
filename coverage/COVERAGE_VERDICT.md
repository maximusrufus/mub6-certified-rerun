# Seed-cover coverage check: verdict

**COVERAGE-PROVEN** (run of 2026-10-20, `covcheck ../impl/seed/seeds.txt`, exit 0; raw output in `COVERAGE_RUN.out`).

## Claim proved
For every bin pair (X-tuple, Y-tail) of the paper's seed enumeration (N = 112) that is not excluded by one of the
paper's necessary conditions (orientation filter, the two midpoint-gap inequalities, the vanishing-sum test, the
shared-interval test), the feasible part of the closed 9-dimensional bin (the points that satisfy the linearised
vanishing-sum constraints `||J h - p|| <= E`, `|h_j| <= 1/2`) lies within Euclidean distance R = 0.032 of a centre
in `impl/seed/seeds.txt` (249,976 centres) with the same selected position s. With the paper's arrangement lemma
and the necessary-condition lemmas, this is the paper's claim `seed:cover` for this seed file, including the
symmetry reduction (the normal form fixes the row/column permutation, rephasing and transposition; the selection
rule fixes the column order; the checker recomputes s from the X-tuple and rejects any centre whose stated s
disagrees).

## Method (independent of the seed generator)
- Own interval type: round-to-nearest double arithmetic widened by one `nextafter` step per operation (no
  dependence on the FE rounding mode; built with `-ffp-contract=off`), own enclosure of pi.
- Own table of unit roots e((2b+1)/224), b = 0..111: python-flint `arb`, 200 bits, outward-rounded to binary64
  (`gen_roots.py` -> `roots_112.hpp`, enclosure widths < 1e-15).
- Own enumeration of the X-tuples (C(116,5) candidates) and Y-tails (C(115,4) per shared bin), own interval
  implementation of each necessary condition (same direction sets as the paper; exact integer arithmetic for the
  gap inequalities).
- Certificate per pair: the pair is itself a centre (its closed bin has half-diagonal 3/224 < 0.032), or for some
  same-s centre u' with k = N(u - u'), ||k||^2 <= 16: box bound `||k||^2 + 9/4 + sum|k_j| < 12.845` (exact
  integers), or dual bound `||k||^2 + 9/4 + upper(sum_j |k_j - (J^T v)_j| + 2 v.p + 2E(||v_X|| + ||v_Y||)) < 12.845`
  with E = 5 pi/(4N) and 12.845 < (0.032 * 112)^2 = 12.845056. The multiplier v is a heuristic (least squares
  followed by damped Newton on a smoothed objective); only the interval upper bound decides.

## Counts
| quantity | checker | seed generator |
|---|---|---|
| X-lists after the three X filters | 22,829 | 22,829 |
| pairs after sumok | 298,638,554 | 298,638,554 |
| pairs after the second gap inequality | 190,192,620 | 190,192,620 |
| pairs after the shared-interval test | 151,742,926 | 151,742,926 |
| covered by their own centre | 249,976 | (249,976 centres) |
| covered by the box bound | 86,187,558 | |
| covered by the dual bound | 65,305,392 | |
| uncovered | **0** | |

The exact agreement of the exclusion counts is a consistency check between two implementations of the same
necessary conditions; the proof is the certificate for every kept pair.

## Negative control
`NEGATIVE_CONTROL.out`: with every second centre removed from the seed file, the checker (first 300 X-lists)
reports 3,020 uncovered pairs and `COVERAGE-FAILED` (exit 1). Removing every 250th centre is still covered
(the cover is redundant), which is expected.

## Run
8 OpenMP threads, 25 min 34 s wall, 12,193 CPU-seconds, 133 MB resident. Inputs hashed in
`COVERAGE_INPUTS.sha256` (seed file 22cd76ed..., `covcheck.cpp`, `roots_112.hpp`).

## Not covered by this check
The arrangement lemma and the necessary-condition lemmas are the paper's theorems and are used as stated. The
binary64 rounding of the rational centres (< 1e-14 per coordinate, absorbed by the search radius 0.0320001) is
argued, not machine-checked, here.
