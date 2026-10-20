# DECISIONS.md -- seed stage (SPEC sec 9 items touching this stage)

Status key: S = soundness-relevant, P = affects only counts/performance.
Evidence for every statement is a test or a measured number in this directory.

## Arithmetic (pre-registration: outward-rounded intervals, no paper margins)
All certificates use `../core/interval.hpp` through `arith.hpp`/`cover.hpp`/`reuse.hpp`.
Every rejection or success is a strict interval comparison:
- unit exponentials `e((b+1/2)/N)` = `iv::unit_root(2b+1, 2N)` (enclosure of the exact rational phase);
- `sumok`: reject iff `lower(|v.c'|) > upper((2pi)(coef*|v| + (w/2) sum|v.J_k|))`, `c' = 2 pi c`, for
  EXACT-double directions v (perpendiculars to mid(J_b), and the normalised midpoint of c');
- gap inequalities: exact integer arithmetic (`2N*gap` is an integer), so no tolerance is needed;
- `extend`: outward interval for `h_s` from enclosures of c, a, r and E = 5pi/(4N); constraints with
  `lower|a| < 1e-6` are skipped (drops a constraint: sound); pair rejected iff the X and Y bounds are
  strictly disjoint; no `1e-7`/`1e-5` additions;
- reuse: box bound exact in integers; dual bound = upper bound of `sum|k-J^T v| + 2 v.p + 2E(|vX|+|vY|)`
  with interval J, p and E = 5pi/(4N) (no 1.251, no +.01); success iff `kn2 + 9/4 + bound < lower((.032 N)^2)`
  (RN2lo = 12.845056, no -.0001). The proposal of v (and its estimated-bound pre-filter) is plain double
  heuristics and never decides success;
- every thread that evaluates a certificate holds an `iv::RoundUp` guard (OpenMP regions, main, tests).
Tests: interval enclosures contain the libquadmath values (tests/test_arith), filters are compared item
by item with a margin-free `__float128` re-implementation (tests/test_cover D).

## Checkpoints (measured, `./seed --count-only`, interval arithmetic)
- X-lists after F1-F3: **22829** (exact match with the paper, unchanged by the move to intervals).
- Pairs entering reuse: **151,742,926** (paper: 151,744,086; **1,160 fewer, -0.0008%**). Intermediate counts are
  unchanged: sumok-only 298,638,554, +gap2 190,192,620.
  EXPLANATION (verified by experiment): the only filter that changed is the shared-interval test `extend`.
  The paper widens its numerator by 1e-7 and each endpoint by 1e-5; those margins keep 1,160 pairs
  whose X and Y intervals are disjoint by less than about 2e-5 but are certainly disjoint under
  rigorous enclosures (widths about 1e-13). Re-adding exactly those two margins to the interval engine
  (diagnostic build `-DDIAG_PAPER_MARGINS`, not a certificate build) reproduces **151,744,086**
  exactly. `sumok` and the gap filters are unchanged because the paper's 1e-8 margins never decided a
  case at N=112 (gaps are exact integers; no sumok list lies within 1e-8 of its boundary).
- Centers: **249,976** vs the paper's 249,878 (+0.04%; the margin-free thresholds give slightly more reuse
  than the margined ones: RN2 12.845056 vs 12.844956, E = 5pi/(4N) vs 1.251 pi/N). Not deterministic in the
  paper (SPEC sec 6); my traversal/reuse heuristics (A2, A4, A26) differ from the source. Runs are
  deterministic (box 92,246,600 / dual 59,246,350 / cache 115,446,376 successes).

## Items
- **A1 (P)** Seed file = SPEC choice (b): one line per center, `s b0 ... b8`,
  `s` in 0..4 (selected position), `b0..b8` the nine bin indices in
  nine-coordinate order (xi, other four x in former order, y2..y5). Stage A
  forms phases as `(b+0.5)/112` in binary64 (division-rounded; error < 1e-14,
  absorbed by radius .0320001).
- **A2 (P)** X-lists lexicographic, Y-tails lexicographic. Reuse needs equal
  `s`, so the five `s` streams are independent and run in parallel; output is
  merged by X index (deterministic). Any order is sound (induction in
  COVERAGE.md uses only "earlier retained centers stay retained").
- **A3 (P)** `sumok` directions: the residual `c` normalised (exact doubles from the midpoint of the
  enclosure; skipped if zero) plus the five perpendiculars to mid(`J_j`). Any exact direction set is sound
  (seed:sum-support holds for every v). Reproduces 22829 x-lists exactly.
- **A4 (P)** `supdual` proposals: least-squares start, then <= 9 damped Newton steps on a smoothed dual
  objective (doubles). The certificate is `dual_cert_hi()` (interval upper bound); non-finite v is rejected.
  A wrong v can only fail to certify.
- **A5 (P)** Y midpoint gap = sorted copy of the five midpoint bins (shared
  bin included) together with the fixed 0, exactly as the paper.
- **A23 (S)** RESOLVED by interval arithmetic: no `expair` error budget is relied on; `cis_turns` of the
  core encloses the exact phase (widest enclosure 3.1e-15 at N=112). The double `expair` in `arith.hpp` is
  heuristic-only (normal-form proposals, test Hadamard generation). `RN2lo` is a rigorous lower bound
  of (0.032 N)^2, computed from the exact ratio 32/1000, so no `pow`/contraction assumption exists;
  the core refuses to compile with `-ffast-math`.
- **A24 (S)** RESOLVED: no scale-dependent margins remain. Guards that remain: finite `v`; `lower|a| >= 1e-6`
  in `extend` (skipping only drops a constraint). Division by an enclosure containing 0 throws in the core
  (never reached: guarded by the `1e-6` test).
- **A26 (P)** Reuse across different X tuples is attempted against earlier
  X-lists with the same `s` and `||k_x||^2 <= cut`; within an X-list the t1
  band `|t1 - t1'| <= floor(sqrt(cut - ||k_x||^2))` is searched (exact filter,
  verified: identical output with and without it).
- **Not in SPEC list: empty one-sided interval.** If the X (or Y) interval is
  empty by itself (lo > hi) the paper's rule "reject only strictly disjoint"
  keeps the pair when it is not disjoint from the other list's bounds. I
  follow the paper literally: conservative (more pairs kept), and it is
  required to reproduce 151744086.
- **Reuse reach cut.** Candidates with `||k||^2 > floor(RN2)+2` (14 at N=112)
  are not tried. This only skips reuse attempts (sound: a skipped attempt
  retains the current midpoint). Raising it can only reduce the center count.
- **Normal forms (nf namespace).** `arrange` implements seed:arrangement
  (min-gap pair among all 30, transpose, largest-gap start, orientation swap,
  selection rule, row sort); repeated points are handled by snapping values
  within 1e-9 below 1 to 0 (valid because every gap of a Hadamard matrix is
  >= 1/6). It is a PROPOSER used by tests; nothing certified depends on it.
  `alignment_*` and `Sym`/`apply_sym`/`transport_T` implement phase:alignment
  and phase:transport; identities are tested numerically (tests/test_normal).
