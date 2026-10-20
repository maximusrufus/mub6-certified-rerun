# Pre-registration (DRAFT): independent certified re-implementation of N(6) = 3

Draft written 2026-10-19 by research-dc, before any search code has run. The founder registers it, using CommitmentNotary register plus a freetsa anchor. This text is frozen at registration; any later change goes in PREREG_AMENDMENTS.md with a timestamp.

## Target claim
openai/math family 266, "The maximum number of mutually unbiased bases in dimension six" (2026-09-24), Theorem 1.1: there are no four pairwise mutually unbiased orthonormal bases in C^6, so N(6) = 3. The original certificate is one recorded binary64 execution of an unreleased pipeline (MUB_ASSESS.md).

## Method (fixed in advance)
1. Re-implement the paper's certified search from its LaTeX alone (SPEC.md). We do not use OpenAI code; the founder decided not to request it.
2. Arithmetic: outward-rounded interval arithmetic in C++ (FE_UPWARD, set per OpenMP thread, -frounding-math, no -ffast-math), with constants held as intervals. Every surviving cell, and a random 1% sample of excluded cells, is re-checked in python-flint acb at 128 bits or more.
3. Checkpoints, i.e. the paper's stated counts: 249,878 seeds; Stage B records 3083 -> 600 -> 22,379; 163 capped records resolved by the fallback. Matching a checkpoint is a diagnostic of implementation fidelity, not a pass criterion. Our rigorous search may legitimately produce different counts.
4. Decision rules and outcomes:
   - **REPRODUCED:** our interval-arithmetic search covers the full parameter space defined by the paper's normal form and excludes every cell, and an independent second implementation of the exclusion test (python-flint) agrees on all surviving and sampled cells.
   - **REFUTED:** we find a candidate quadruple of bases whose mutual-unbiasedness defects are certified, by interval Newton / Krawczyk or an exact algebraic check, to vanish within a certified neighbourhood. Any claimed refutation must also survive an independent blind confirmation before it is reported anywhere.
   - **SPEC-GAP:** the paper's text is insufficient to define a step (an AMBIGUOUS item in SPEC.md that cannot be resolved from the paper or standard results). We report which step and why.
   - **INCONCLUSIVE:** the cover cannot be completed within the compute budget, or cells remain that neither exclusion nor a certified solution resolves.
5. Coverage check: the union of seed cells is proven, separately from the search, to cover the normal-form parameter space, including the symmetry reductions as stated and proved in the paper.
6. Compute budget: at most 72 CPU-hours on this workstation under a CPU lease, and no paid compute. If the budget runs out, the outcome is INCONCLUSIVE with the residual cells reported.

## What we will and will not claim
- REPRODUCED means an independent, rigorous computer-assisted reproduction of the exclusion, under the paper's normal-form reduction. It is not a peer review of that reduction, which we audit separately and report.
- We will not claim novelty for N(6) = 3 itself. Credit belongs to the original authors, and the prior work (Jaming–Matolcsi–Móra–Szöllősi–Weiner 2009 and the discretization programme) is cited.

## Analysis plan
Report the outcome class, the counts against the checkpoints, the full list of residual cells (if any), runtimes, and hashes of the code and outputs. Release code and certificates only after founder approval, through the posting gate.
