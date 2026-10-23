# Reproduction report: OpenAI's claimed N(6)=3 (mutually unbiased bases in dimension six)

Date 2026-10-23. Scope: one claim (openai/math entry 266), five fixed questions. Credit: the claim and method are OpenAI's; the published data we cross-checked is Manvir Jaswal's; the reduction lemmas are in Jaswal's paper; carsoncrocker has been invited to read them; Carlos Toledo's cert-machine audits the same release at entry level. This is not peer review and not a formal proof.

## 1. What exactly was claimed?

Preprint "The maximum number of mutually unbiased bases in dimension six", author OpenAI, dated 2026-09-24, family 266 of github.com/openai/math.
- Pin: openai/math commit `1581d4cbdd69a0faf34b2089614fe2a814316a6a` (2026-10-12). PDF sha256 `b4cc74c840c902494d589b51c8e9706490f37cf2735e69a0486122131c6b7446`.
- Statement (Theorem 1.1, verbatim): "Under the arithmetic and compiler conditions in Appendix A, the fresh complete execution documented in Section 7.5 certifies that there are no four pairwise mutually unbiased orthonormal bases in C^6. Consequently, N(6) = 3." The second basis is unrestricted.
- Lean scope: `lean/docs/266.md` in the same repository says the Lean formalisation covers a weaker bound (any family has at most five members) plus a Fourier-vanishing lemma, and "does not establish the paper's upper bound of three or its computer-assisted exclusion of four arbitrary bases." So the headline is not formalised.

## 2. What data and code are required, and what was released?

- Required by the paper (Section 7, Appendix A): the verification/code directory (seven C++ sources including stageC.cpp), extraction.json with SHA-256 digests, proof_run.py, build manifest, REPRODUCE.md, and the execution archive full-run-01/.
- Released in openai/math at the pin: the PDF, LaTeX sources (build/sections/*.tex) and README. None of the files above is present (no matching path among the 133,956 files at the pin).
- the search code and run archive were not in the release. cert-machine's `corpus/openai-math/findings.json` (finding "missing-artifacts", checked 2026-10-23) records that the paper "cites a verification/code directory of seven programs and an extraction.json — absent, so the computer search behind its headline cannot be re-run" (dash normalised) (https://github.com/carlostoledo1891/cert-machine/blob/main/corpus/openai-math/findings.json).
- Consequence: anyone re-running it must rebuild the pipeline from the paper text.

## 3. Can it be reproduced from the paper's text?

Yes, as a computation, from the paper text alone (no OpenAI code was used or requested). Public repo: https://github.com/maximusrufus/mub6-certified-rerun. Evidence: `RESULT.md` and `coverage/COVERAGE_VERDICT.md` in that repository.
- Pre-registration: `preregistration/PREREGISTRATION.md`, sha256 `ae6448ca48ff62a2f03460f8451cbd59d30df31f749266276469a8d07c65b30e`, FreeTSA timestamp 2026-10-19 22:29:54 GMT, before the run (finished 2026-10-20).
- Search: 7,726,523,584 cells processed; 1,166 Stage A records, 0 bad, 0 unresolved; Stages B and C ended with 0 unresolved and pipeline_check 0 failures. 12.4 h on 16 threads, $0.
- Separate seed-cover check (own interval type, own root table, own enumeration): all 151,742,926 kept bin pairs certified, 0 uncovered; negative control (every second centre removed) fails as it must.
- Cross-check: python-flint re-check of a 1/100 sample of excluded cells (72,265,439 cells), 0 disagreements; all 1,166 Stage A records re-certified, 0 failures.
- Differences from the paper's counts (249,976 seeds vs 249,878; 1,166 Stage A records; 15,595 Stage B roots) were pre-registered as diagnostics only and trace to implementation choices.
- Outcome under the pre-registered rules: REPRODUCED, with the caveats in section 5.

## 4. Does the verification establish the full claim?

Only the computation, under the paper's hand-proved normal-form reduction. It does not by itself establish N(6)=3.
- The reduction lemmas are proved in M. Jaswal's paper v2, Appendix C.2, Lemmas C.4 and C.5 (https://zenodo.org/records/23533384). carsoncrocker has been invited to read them (https://github.com/google-deepmind/formal-conjectures/issues/6918). We have not reviewed them.
- Cross-check of Jaswal's published data (Zenodo 10.5281/zenodo.23255817, v1, CC BY 4.0, read as data, none of his code run): all 57,480 stored root records pass our checker (root conditions only): 57,480 in 200-bit python-flint, 0 rejected; 1,718 are valid only in transposed coordinates, covered by his Lemma C.24 (residual equality checked numerically, max relative difference 1.1e-43). 26,677 exact-solution tests, 0 failures. This is not a replay of cell-level exclusions: neither run stores them. The cross-check scripts and outputs are available on request.

## 5. What could not be established, and why?

- The normal-form reduction is not reviewed by us (hand-proved in the cited lemmas; no formalisation).
- Arithmetic trust base: our C++ interval code (outward rounding, rounding-mode assertion in every binary) plus the python-flint cross-checks. The fallback stage (5 of 61 batches) uses the pipeline's own interval code, and the separate excluded-cell check samples 1%, not all cells.
- The binary64 rounding of the rational seed centres is argued, not machine-checked (`coverage/COVERAGE_VERDICT.md`, "Not covered").
- The paper's own Appendix A conditions (binary64 arithmetic, compiler) apply to OpenAI's run; ours is a different implementation, so we say nothing about their execution archive.
- Nobody outside our side has reproduced our run; nothing here is peer reviewed.
- Nothing beyond the Lean scope note in section 1 is formalised.
