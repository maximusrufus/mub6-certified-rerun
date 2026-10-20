# Result: independent certified re-run of N(6) = 3 (openai/math entry 266)

Run finished 2026-10-20. Pre-registration: preregistration/PREREGISTRATION.md (sha256 ae6448ca…, FreeTSA 2026-10-19 22:29:54 GMT). Code: certificates/CODE_MANIFEST.sha256 (63f67c49…). Outputs: certificates/OUTPUT_MANIFEST.sha256.

## Outcome under the pre-registered rules: REPRODUCED (with the caveats below)
- **Seeds:** 249,976 seeds were regenerated from the paper's text. The regenerated set is identical to the committed seed file (sha256 22cd76ed…).
- **Stage A:** 7,726,523,584 cells processed, ending in 1,166 records, 0 bad and 0 unresolved. Arithmetic was outward-rounded intervals with a rounding check.
- **S5 (independent re-certification):** all 1,166 records were re-certified with 0 failures; 40 sampled records passed the exact-solution normal checks.
- **S6 (independent python-flint check of excluded cells):** 72,265,439 cells were sampled (1/100) with 0 disagreements.
- **Stages B and C:** 61 batches and 15,595 roots; 5 batches went to the fallback stage. Final: RESULT 0 unresolved, and pipeline_check reported 0 failures.
- **Wall time:** 44,704 s (12.4 h), on 16 threads, at $0.

## Checkpoints (diagnostic only, not pass criteria)
The paper reports 249,878 seeds and Stage B records 3083 → 600 → 22,379, with 163 records resolved by the fallback. Our counts differ: 249,976 seeds (the paper's margins were removed; re-adding them reproduces its pair count), 1,166 Stage A records, and 15,595 Stage B fine roots. The counts depend on implementation choices (symmetry merging, interval margins), not on correctness.

## Caveats and what was not checked
- This reproduces the computer search under the paper's normal-form reduction. It is not a peer review of that reduction (pre-registration, "what we will claim").
- Coverage check (pre-registration item 5): the seed cells match the paper's construction exactly, but we have no separate machine proof that the seed union covers the normal-form space. UNTESTED.
- The fallback stage still uses the pipeline's own interval code. S6's independent check samples 1% of excluded cells, not all of them.
- No outside party has reproduced our run. The code and certificates are unreleased (founder approval required).
- Words to use: "independently reproduced by a certified computation", never "verified" or "proved".
