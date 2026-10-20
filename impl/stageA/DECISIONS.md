# Stage A decisions (arun re-implementation)

Scope: SPEC section 3 (frames, polar bound, U(2) cells 16 -> 512, z/E_r, norm + directional tests,
Root certificates, Candidate arc bounds, ellipsoid/ell/bootstrap, atlas merge). Interval layer:
`../core/interval.hpp` (per-thread `iv::RoundUp`). Seed input: `../seed/seed` output (A1 option (b)).

## 0. Soundness policy (what is interval-evaluated)

No margin from rounding.tex is copied. A certificate inequality is "A rejects/covers iff lo(LHS) > hi(RHS)"
(or the bound is `hi` of an interval expression); proposals (eigenvectors, Gauss-Newton, direction search,
Cholesky solves, atan2, libm) are plain binary64 and are only ever fed back as EXACT doubles into interval
expressions. NaN/inf/negative-sqrt/zero-division make a test fail (never exclude).

| Inequality | Where | Evaluation |
|---|---|---|
| seed center -> e(q) | `buildSeed` | `iv::cis_turns` enclosures of the stored binary64 phases (treated as exact) |
| frames K_B, K_C, G_B, G_C, P_K, D_span, m, A_0, A_* | `buildSeed` | interval Gram-Schmidt; coordinate-vector choice by midpoint only |
| q, p, z^2, L, e_off, E_r per level | `buildSeed` | interval, stored as upper doubles; `h = hi(pi/N_u)` |
| norm test sum(|D|-1)^2 > z^2 | `normExcluded` | `lo(sum) > z2U`, strict, no 1e-5 margin |
| directional test | `certifyDir` | J_all, c_all, D, `y` all intervals/exact doubles; exclude iff `lo(y.c) > hi(RHS)`, no .0001 margin |
| whole-ball sum test | `seedSumExcluded` | same inequality with y supported on the X/Y rows |
| Root: froot, ort, hor, mnorm, id, C, n0 | `root.cpp certify` | interval (see 3) |
| arc/weighted/eps, direct bound s0, S_1, S_2 | `cover.cpp` | interval, hi of each |
| ellipsoid G, offsets, ell, t0, bootstrap | `cover.cpp` | interval; G dominated by `dominate()` (see 4) |
| merge distance, radius growth, cap | `atlas.cpp` | interval, `hi` |

Radii are `hi` values; isolated Roots (k = 0) therefore get radius ~1e-15 (the certified bound really is
that small). Stage B must accept `R >= 0` (SPEC 4.1 validation does).

## 1. Seeds and parameters

* **A1** input line = `s b0..b8` (selected position + nine bin indices, N = 112) exactly as `../seed/seed`
  writes it; phase = `((double)b + .5)/112`. Nine decimal phases are also accepted (tests). `s` is not used by A.
  The double differs from the rational midpoint by < 1e-16, absorbed by .0320001 vs .032 (gap 1e-7).
* Seed radius `r_s`: `nextafter(.0320001)`; `R_s = hi(2*pi*r_s)` (no extra +1e-9; nothing to absorb).
* **A6** `levels=3` = grids 16/32/64 are the list-held levels (we use an explicit DFS stack, no lists);
  `iters=2` is not needed (direction search below is stronger: 5 regularised starts + <=30 ascent steps);
  `first=1` = first level at which the directional test is applied is the SECOND grid (N_u = 32). Only the norm
  test runs at 16. Running more tests than the paper is sound. Cover attempts start at 64 (`firstCoverNu`).
  `topK = 4` nearest local Roots (ranked by approximate ||Delta||). `maxNu = 512`. Node cap: cells with
  N_u >= 128 counted; > 1,000,000 => seed unresolved (`extra node cap`), no exclusion.
* **A7** coordinate vectors in Gram-Schmidt: largest midpoint residual first; interval residual norm must have
  `lo > .449999` (paper's threshold, kept).
* **A8** directions: y = c - J (J^T J + lam I)^-1 J^T c for lam in {1e-9,1e-6,1e-4,1e-2,1e-1} plus y = c; best
  by approximate V(y); then <= 30 projected subgradient steps on the sphere (V is concave, positively
  homogeneous). Final y: divide by max |entry| (exact doubles), reject if max not in (1e-100, 1e100).
  Experiment: 0 / 30 / 300 ascent steps -> 60.4k / 33.6k / 32.9k covered (cover-attempted) cells on 51 real seeds (saturated by 30).
* Level-16 cells: norm test only (SPEC). Directional at >= 32.

## 2. A24 scale guards (soundness-relevant ambiguity)

With interval arithmetic every error term is tracked, so the paper's scale bounds are NOT needed for validity.
Decision: keep every published guard that is an input/scale condition (extra failures only cost work, never
soundness) and enforce the lemma preconditions rigorously:

* frames: lower bound of `m = min(m_B,m_C)` > **1.24** (paper: stored m > 1.25; its ideal bound is 1.24 -- we use
  the ideal bound because our `m` is an enclosure of the ideal value); residual norms lo > .449999.
* `q < .32`, `p < .068` (hi values, every level), `off = R_s(sqrt6+A_*)/m < .74`, `z < 1.5`, `E_r < 1` (level 16).
  Lemma preconditions `q < 6`, `|S| < 4` follow.
* arc stage: `lo(m_j) > .15`, `q_d < 3.8`, `arc < 8.8`, `weighted < 10`, `eps < 1e-3`, initial radius `< 3`.
* Root: `|P_j| <= 4`; `|T|,|H_V|,|L_V|` entries <= 2; k <= 10; eigenvalue cut .64; `val = ||F||^2 < 1e-18`;
  `id < 1e-4`, `ort, hor < 1e-5`, `C < 20`, `mnorm < 2`.
* ell: `gamma < 418`, `gap >= 1e-6` (gap is a certified lower bound of gamma - lambda_max), `|z_i| < 1e8`.
* NOT enforced: the paper's `nd <= .145^2` direct-path guard (it only protects the paper's margin; our bound is
  rigorous, and a successful cover has `s >= sqrt(nd)` with `s <= .103`, so it holds anyway); the paper's
  per-operation +1e-k margins (replaced by intervals).
* Hopeless-candidate gate in `coverBound` (speed only): if `s > target` and `n_P(s)/s >= 1` the bootstrap cannot
  improve, so the candidate is rejected without the tangent computation. Returning "not covered" is never an exclusion.

## 3. A9 Root construction (decided and why it is sound)

For a center P (after Gauss-Newton snapping `buildRoot`, horizontal steps only, stops at ||F|| < 2e-15):
`Bh` = Helmert (x) Helmert, 36x25 orthonormal horizontal basis; `A = L_P Bh`; Jacobi eigendecomposition of
`A^T A`; eigenvalues `< .64` -> columns of T (`T = Bh v`); the rest give
`H_V = Bh v / sqrt(lam)` and `L_V = A v / sqrt(lam)`; **`M := H_V L_V^T` is defined as the exact real product of the
stored doubles** (so there is no "accumulated M vs exact product" gap to bound). Then, all in intervals with
exact-double inputs:

* `mnorm = sqrt(ah*al)` where `ah, al` = max absolute row sum of the interval Grams `H_V^T H_V`, `L_V^T L_V`
  (each bounds a squared operator norm) -- so `||M|| <= ||H_V|| ||L_V|| <= mnorm`;
* `froot = hi ||F(P)||`, `n0 = hi(mnorm*froot)` (so `mnorm*froot <= n0`);
* `C = hi(pi * sqrt8 * mnorm)` (Lemma phase:remainder-constants, first choice);
* `id = hi ||Pi - T T^T - M L_P Pi||_F` over the full 36x36 identity (Frobenius >= operator);
* `ort = hi ||T^T T - I||_F`, `hor = hi ||(I-Pi) T||_F`.

These are exactly the hypotheses of Proposition phase:root-bounds, so ANY T/H_V/L_V passing them is a valid
Root; the eigen-decomposition quality only affects acceptance. **Tensor improvement of C** (round:tensor-c,
implemented 2026-10-19 in `../core/tensor_c.hpp`, shared with Stage B): `C = min(pi sqrt8 mnorm, pi sqrt(6 row/(1-beta)))`
where `A = (M^T M) o E` (E = edge-sign matrix of phase:tensor-gram, `M^T M = L_V (H_V^T H_V) L_V^T` in intervals),
`U` = Jacobi eigenvectors of mid(A) (proposal, entries <= 2), `beta = hi ||U^T U - I||_F < .1` (else the tensor
candidate is +inf), `row = max_a hi sum_b |(U^T A U)_ab|`. The certificate is the paper's congruence argument
(`||A|| <= ||U^-1||^2 ||U^T A U|| <= row/(1-beta)`), evaluated in intervals with no `.0002`/`.00001` margins.
Measured on 793 real Stage A records: `C_tensor/C_plain` median .954 (min .913, max 1.08 -> plain kept); a
6% effect only (see ../pipeline/PERF_FIX.md). `MUB_NO_TENSOR_C=1` reverts to the plain C (measurement switch).
`checkRoot` re-derives M from (P, T) with a fresh eigendecomposition and accepts a stored Record only if every
stored constant is >= the freshly certified one (fresh `C` is the same min) and `n0 >= mnorm*froot`.

## 4. Cover computation (Candidate::prep, covercompute, ell, bootstrap)

* Candidate center P0: first row/column 0; seed entries = stored q; lower-right 16 = `atan2(mid D)/2pi`;
  `eps = hi sqrt(sum (kappa_j/4)^2)`, `kappa_j = hi |e(P0_j) - D_j/|D_j||` (chord >= 4 * phase distance), replacing
  the paper's `|e - D/m| < 1e-7` test. No other check is needed because kappa is bounded by interval arithmetic.
* Alignment: `N = floor(P0 - mid(D0(P)) + .5)` (A11), `Delta = Pi(P0 - D0(P) - N)` in intervals (any integers valid).
* Direct bound: `s0 = hi( sqrt(|Delta|^2 + (R_s^2+arc^2)/(2pi)^2 + min(S_1,S_2)/pi) + eps )`; also the always-valid
  `|Delta| + eps + sqrt(R_s^2+arc^2)/(2pi)`; the min is the initial radius.
* Tangent (k > 0): `G_lin` (beta = .45, 16 corner offsets, sine term) and `G_dir` (beta = .8, `invmag` = hi(1/lo(m_j)));
  `t0 = min` of both. Interval entries of G are dominated by an exact symmetric double matrix
  `G_hat = mid(G) + rho I`, `rho` = max row sum of entry deviations (so `G_hat - G_true` has spectral norm <= rho
  and `G_hat >= G_true` as PSD matrices) -- this replaces the paper's `+1e-7 I`.
* `ell`: lambda_max bound = min(row sum, `rowsum(V^T G V)/(1-beta)` with `beta >= ||V^T V - I||_F < .5`, V from a
  Jacobi eigendecomposition); gamma tried at `gu(1+x)+2e-6`, x in {1e-4,..,4}; residual bound
  `U^T z + ||U|| ||r|| / gap` with U the interval offset, z an exact double Cholesky proposal, r = U - (gamma I - G)z
  in intervals; also the triangle bound `max||U|| + sqrt(gu)`; the minimum is returned.
* Bootstrap: <= 30 updates, `s <- sqrt(n_P(s)^2 + (1+ort) t0^2)` only when smaller. Target .103.
* `local::test`: top-4 local Roots by approximate ||Delta||; else `buildRoot(P0)` and the same bound. Failure
  => the cell is split (never excluded). A Record's radius grows to `max(r, s)`.

## 5. Atlas (A13)

Gauge Phi = (transpose?, sign, row permutation, column permutation) = the seed-arrangement rule: for all ordered
row pairs of M or M^T (and of -M), keep pairs whose largest circular gap is within .02 of the global minimum;
start columns after any largest gap (within 1e-3); remaining four rows sorted by their dephased entry in column
position 1. Buckets by `floor(100 * min gap)` (the min gap is a Phi-invariant), +-1 neighbouring bucket. The best
variant/atlas center is chosen with a double distance, then `d = hi ||Pi(D0(Phi P) - D0(Q) - N)||` is recomputed in
intervals; merge iff `hi(r + d + 1e-8) + 2e-8 <= .125`, new radius `hi(r+d+1e-8)`. Otherwise the transformed Root
(P and T columns permuted/negated) is inserted with its radius; since 2026-10-19 a non-identity transform is
**re-certified from scratch in the stored orientation** (`recertifyTransported`: `certify` with the stored T,
inflated constants) and those fresh constants are stored (the transported ones, valid by phase:transport, are
kept only if the re-certification fails), so `checkRoot` on the output file is deterministic for every record.
`arun` processes seeds in parallel chunks (dynamic schedule) and merges the resulting local Records SERIALLY in
seed-index order into one atlas, so the output file is byte-identical for any thread count (checked: 1 vs 8
threads, checked vs production build). A per-thread atlas followed by a final merge was tried first and gave
2.4x more Records on the sample (1226 vs 517) because radii grown in one atlas block later merges (r + d <= cap).
Roots whose seeds contain a family of Hadamard matrices produce many local Records; counts differ from the paper's 3083.

Re-verification margins: `buildRoot` inflates the stored constants by tiny absolute amounts (froot +1e-13,
id +1e-11, ort/hor +1e-12, mnorm x(1+1e-9), then n0 and C recomputed from them). This is NOT a soundness margin (all
constants are only upper bounds, larger is still valid); it makes `checkRoot` on a symmetry-transformed Record
(different rounding in the fresh eigensolve) robust. `checkRoot` itself compares against UN-inflated fresh constants.

## 6. Not done / untested (explicit)

* The full 250,490-seed list has not been processed here (sampled only); see the report for the sample numbers.
* No `Usearch`/`scalar` tests (paper does not use them). Tensor C: implemented (3), 6% effect.
* Local record radii (measured, `arun` prints a histogram): all >= .06, most in [.08, .125) -- the within-seed
  family patch is ~.1 wide, so records are fat; the atlas cap .125 then makes the 793 smoke records overlap many
  times. This is the paper's design too (local target .103, cap .125).
* Atlas size cap (2,000,000) is not enforced in `arun`; sizes seen are far below it.
* The seed-level whole-ball sum test is an addition (sound; seeds from the cover satisfy `sumok` so it rarely fires).
* Independent cross-check: `tools/verify_records.py` (python-flint 160-bit) re-evaluates ||F||, ort, hor and checks
  the Root assertions on exact solutions found by high-precision Newton; it is a sanity check, not a proof.
