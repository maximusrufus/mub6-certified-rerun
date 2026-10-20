# Independent certified re-run of the computer search for N(6) = 3

This repository is an independent, from-the-text re-implementation and certified re-run of the
computer-assisted part of

> *The maximum number of mutually unbiased bases in dimension six* (openai/math preprint family 266,
> 2026-09-24), Theorem 1.1: there are no four pairwise mutually unbiased orthonormal bases in C^6, so N(6) = 3.

Credit for the result and the method belongs to the original authors. Nothing here is new mathematics. What is
new is a second, independent execution of their search with outward-rounded interval arithmetic, a second
implementation of the exclusion test (python-flint, 128+ bits) applied to every surviving record and to a 1 %
sample of excluded cells, and a separate machine proof that the seed cover is complete.

**Words.** The outcome is "independently reproduced by a certified computation". It is not "verified" and not
"proved" in the sense of a formal proof assistant, and it is not a peer review of the paper's normal-form
reduction (see Caveats).

## Outcome

Under the pre-registered decision rules ([preregistration/PREREGISTRATION.md](preregistration/PREREGISTRATION.md)):
**REPRODUCED**, with the caveats below. [RESULT.md](RESULT.md) is the result note as written when the run
finished.

| stage | result |
|---|---|
| Seeds | 249,976 seed cells regenerated from the paper's text; sha256 `22cd76ed...` (identical to the committed `impl/seed/seeds.txt`) |
| Seed cover | the union of the seed cells covers the normal-form parameter space: **COVERAGE-PROVEN** by an independent checker (see below) |
| Stage A | 7,726,523,584 cells processed; 1,166 records; 0 bad, 0 unresolved |
| S5 | all 1,166 records re-certified independently (0 failures); 40 sampled records pass exact-solution checks |
| S6 | 72,265,439 excluded cells (1/100 sample) re-checked with python-flint: 0 disagreements |
| Stages B, C | 61 batches, 15,595 roots, 5 batches through the fallback; final `0 unresolved`, pipeline_check 0 failures |
| Wall time | 44,704 s (12.4 h) on 16 threads |

The paper reports 249,878 seeds, Stage B records 3083 -> 600 -> 22,379 and 163 records resolved by the fallback.
Our counts differ (249,976 seeds, 1,166 Stage A records, 15,595 Stage B fine roots). The counts depend on
implementation choices (symmetry merging, interval margins), not on correctness; they were pre-registered as
diagnostics, not as pass criteria. Re-adding the paper's margins to the seed stage reproduces its pair count exactly
(`impl/seed/DECISIONS.md`).

## Caveats and what was not checked

(verbatim from RESULT.md, except that the coverage item is now resolved)

- This reproduces the computer search under the paper's normal-form reduction. It is not a peer review of that
  reduction (pre-registration, "what we will claim").
- Coverage check (pre-registration item 5): resolved after RESULT.md was written. The independent checker in
  `coverage/` proves that the union of the seed cells covers the normal-form parameter space; see "Seed cover"
  below for exactly what is and is not covered by that proof.
- The fallback stage still uses the pipeline's own interval code. S6's independent check samples 1 % of excluded
  cells, not all of them.
- No outside party has reproduced our run.
- Words to use: "independently reproduced by a certified computation", never "verified" or "proved".

## Method

The pipeline follows the paper's stages, each re-implemented from the LaTeX alone (no original code was used or
requested):

1. **Seed cover** (`impl/seed`): bins of width 1/112 on the nine seed coordinates, the paper's necessary
   conditions (orientation, midpoint-gap inequalities, vanishing-sum test, shared-interval test) and the
   centre-reuse certificate (box bound and dual bound). Every rejection and every reuse is a strict comparison on
   interval enclosures; the paper's additive margins are not used.
2. **Stage A** (`impl/stageA`): completion of each seed ball to candidate Hadamard matrices, U(2) cell
   subdivision, root certificates, phase-ball coverage and the atlas merge.
3. **Stages B and C** (`impl/stageBC`): refinement of the records to radius .016 (batched), the vector-bin cover
   and graph exclusion, the fallback refinement on unresolved batches, and the final exclusion.
4. **Verification layers** that are not in the paper: `S2` start-up self-tests of the rounding mode and interval
   constants in every binary; `S5` independent re-certification of every Stage A record in C++ and in
   python-flint; `S6` python-flint re-verification of a deterministic 1/100 sample of excluded Stage A cells;
   hand-off format checks between stages; a structural check of the B/C pipeline output.

Arithmetic: outward-rounded intervals in C++ (`impl/core/interval.hpp`; FE_UPWARD per thread, `-frounding-math`,
`-ffp-contract=off`, never `-ffast-math`; every official binary is built with `-DIV_CHECK_ROUNDING`, a
per-operation assertion on the rounding mode). Constants are held as intervals.

### Seed cover (the coverage proof)

`coverage/covcheck.cpp` is written independently of the seed generator: its own interval type (round-to-nearest
plus `nextafter` widening, no dependence on the FE rounding mode), its own table of unit roots
(`coverage/roots_112.hpp`, python-flint `arb` at 200 bits, outward-rounded, generated by `coverage/gen_roots.py`),
its own enumeration, filters and certificate. The only shared inputs are the paper's formulas and the seed file
being checked.

What it proves (exit 0, last line `COVERAGE-PROVEN`): for every bin pair (X-tuple, Y-tail) of the paper's
enumeration that is not excluded by one of the paper's necessary conditions, the feasible part of the closed
9-dimensional bin (the points satisfying the linearised vanishing-sum constraints of the paper's `grid-error`
lemma) lies within Euclidean distance 0.032 of a centre in `impl/seed/seeds.txt` with the same selected position.
Together with the paper's arrangement lemma (every exact Hadamard matrix has a normal form whose nine coordinates
fall in some bin pair, with the stated selection rule) and its necessary-condition lemmas (none of which rejects an
exact point), this is the paper's claim `seed:cover` for this seed file. The conversion of the rational centres to
binary64 moves each by less than 1e-14, which the search radius 0.0320001 absorbs.

What it does not prove: the arrangement lemma and the necessary-condition lemmas are theorems in the paper; they
are used, not re-proved. The checker's exclusions coincide exactly with the generator's (same 22,829 X-lists and
the same kept pairs on every tested prefix), which is a consistency check, not part of the proof.

Result of the full run (`coverage/COVERAGE_VERDICT.md`, raw output `coverage/COVERAGE_RUN.out`): 22,829 X-lists,
151,742,926 kept pairs (both counts identical to the generator's), 249,976 pairs covered by their own centre,
86,187,558 by the box bound, 65,305,392 by the dual bound, **0 uncovered**; 25.5 min on 8 threads. Negative
control: with every second centre removed, the checker reports 3,020 uncovered pairs and fails
(`coverage/NEGATIVE_CONTROL.out`).

## How to rerun

Requirements: Linux (tested on Ubuntu, GCC 15.2, 32 threads), GNU make, `python3` with `python-flint` (>= 0.9) and
`numpy`, `bc`, `sha256sum`.

```
bash smoke.sh                       # unit tests of every stage + coverage checker on 200 X-lists (~10-15 min, 8 cores)
bash impl/pipeline/smoke_test.sh    # 1 % end-to-end sample with all verification layers (~30 min, 16 threads)
THREADS=16 bash impl/pipeline/run_all.sh impl/pipeline/work   # the full run (12.4 h on 16 threads)
(cd coverage && make check)         # the full coverage proof (see coverage/COVERAGE_VERDICT.md for the run time)
```

`run_all.sh` ends with exactly `0 unresolved` (exit 0) or a residual list with reasons (exit 1); exit 2 is a
tool, format or verification failure. It writes a SHA-256 manifest of every input and output.

## Certificates and hashes

- `certificates/`: the full run's outputs. `A.bin` (Stage A records), `mid.bin`, `fine.bin`, `batches/` (Stage B),
  `fb1.bin` (fallback), the `*bad*` files (all empty), `MANIFEST.txt` (SHA-256 of every binary, input and output
  of the run), `fullrun.out` (the driver's output) and `logs/`. The 10 GB file `excl_sample.txt` (the S6 sample)
  is not included; its hash is in `MANIFEST.txt`.
- `certificates/CODE_MANIFEST.sha256`: SHA-256 of every source file of the run (the files under `impl/` match it;
  check with `cd impl && sha256sum -c ../certificates/CODE_MANIFEST.sha256`).
- `certificates/OUTPUT_MANIFEST.sha256`: SHA-256 of `MANIFEST.txt` and of the original `fullrun.out`. The shipped
  `fullrun.out` and the shipped logs differ from the originals in exactly one way: the absolute path of the work
  directory is replaced by `<work>` / `<pipeline>`.
- `coverage/COVERAGE_INPUTS.sha256`: hashes of the seed file and the checker sources at the time of the proof run.
- `preregistration/`: the pre-registration text, its SHA-256 and the RFC 3161 timestamp (FreeTSA, 2026-10-19
  22:29:54 GMT; request `.tsq`, response `.tsr`). The text was frozen before any search code ran. One line of the
  original names the author's internal registration tooling and is redacted in the shipped copy; the hash and the
  token cover the unredacted original. Verify with
  `openssl ts -verify -in PREREGISTRATION.tsr -queryfile PREREGISTRATION.tsq -CAfile <FreeTSA cacert.pem>`.

## Licences

Code: MIT (`LICENSE`). Data (certificates, seed file, manifests, pre-registration files, root table): CC0 1.0
(`LICENSE-DATA`).
