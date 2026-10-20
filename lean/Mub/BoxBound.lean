import Mathlib
/-!
# Soundness of the seed-cover box-bound centre-reuse test

Source: `../LEAN_PATH_PLAN.md` section 4 and `../SPEC.md` section 2.7
(`seed:distance-bound`, `seed:box-support`).

Setting. `N = 112`; centres `u, u' ∈ (1/N) ℤ^9`, so `k = N (u - u')` is an
integer vector; `q ∈ ℝ^9` lies in the closed bin around `u`:
`|N (q - u)_j| ≤ 1/2` for every `j`.  Claim:

  `N² ‖q - u'‖² ≤ ‖k‖² + Σ_j |k_j| + 9/4`.

Corollary used by `covcheck.cpp::certified` (line 303):
if `‖k‖² + 9/4 + Σ|k_j| < 12.845` then `‖q - u'‖ < 0.032`
(because `12.845 < (0.032 · 112)² = 12.845056`).

Nothing here depends on `N = 112` except the corollary.
-/

open Finset BigOperators

namespace Mub

/-- One coordinate: if `|h| ≤ 1/2` then `(k + h)² ≤ k² + |k| + 1/4`. -/
theorem coord_bound (k h : ℝ) (hh : |h| ≤ 1 / 2) :
    (k + h) ^ 2 ≤ k ^ 2 + |k| + 1 / 4 := by
  have hk0 : 0 ≤ |k| := abs_nonneg k
  have hkh : k * h ≤ |k| * |h| := by rw [← abs_mul]; exact le_abs_self _
  have h1 : |k| * |h| ≤ |k| * (1 / 2) := mul_le_mul_of_nonneg_left hh hk0
  have h2 : h ^ 2 ≤ (1 / 2 : ℝ) ^ 2 := by
    rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg h) hh 2
  have e : (k + h) ^ 2 = k ^ 2 + 2 * (k * h) + h ^ 2 := by ring
  rw [e]; norm_num at h2; linarith

/-- The box bound, in the exact form of the plan: `N` any real, `u u' q : Fin 9 → ℝ`,
`k : Fin 9 → ℤ` with `k = N (u - u')`, and `q` within half a bin of `u`. -/
theorem box_bound_sound (N : ℝ) (u u' q : Fin 9 → ℝ) (k : Fin 9 → ℤ)
    (hk : ∀ j, (k j : ℝ) = N * (u j - u' j))
    (hq : ∀ j, |N * (q j - u j)| ≤ 1 / 2) :
    N ^ 2 * ∑ j, (q j - u' j) ^ 2 ≤ ∑ j, ((k j : ℝ) ^ 2) + ∑ j, |(k j : ℝ)| + 9 / 4 := by
  have hcoord : ∀ j, N ^ 2 * (q j - u' j) ^ 2 ≤ (k j : ℝ) ^ 2 + |(k j : ℝ)| + 1 / 4 := by
    intro j
    have e : N * (q j - u' j) = (k j : ℝ) + N * (q j - u j) := by rw [hk j]; ring
    have : N ^ 2 * (q j - u' j) ^ 2 = (N * (q j - u' j)) ^ 2 := by ring
    rw [this, e]
    exact coord_bound _ _ (hq j)
  calc N ^ 2 * ∑ j, (q j - u' j) ^ 2
      = ∑ j, N ^ 2 * (q j - u' j) ^ 2 := by rw [Finset.mul_sum]
    _ ≤ ∑ j : Fin 9, ((k j : ℝ) ^ 2 + |(k j : ℝ)| + 1 / 4) := Finset.sum_le_sum (fun j _ => hcoord j)
    _ = ∑ j, ((k j : ℝ) ^ 2) + ∑ j, |(k j : ℝ)| + 9 / 4 := by
        rw [Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ,
          Fintype.card_fin, nsmul_eq_mul]
        norm_num

/-- Exact-integer form of the checker's test (`covcheck.cpp` line 303, `RN2_LO = 12.845`). -/
def boxCert (k : Fin 9 → ℤ) : Bool :=
  decide ((∑ j, k j ^ 2) * 1000 + 2250 + (∑ j, |k j|) * 1000 < 12845)

/-- Corollary used by the checker: with `N = 112`, a passing `boxCert` implies
`‖q - u'‖ < 0.032`, i.e. `‖q - u'‖² < 0.032²`. -/
theorem boxCert_sound (u u' q : Fin 9 → ℝ) (k : Fin 9 → ℤ)
    (hk : ∀ j, (k j : ℝ) = 112 * (u j - u' j))
    (hq : ∀ j, |112 * (q j - u j)| ≤ 1 / 2)
    (hc : boxCert k = true) :
    ∑ j, (q j - u' j) ^ 2 < (0.032 : ℝ) ^ 2 := by
  have hb := box_bound_sound 112 u u' q k hk hq
  unfold boxCert at hc
  have hz : (∑ j, k j ^ 2) * 1000 + 2250 + (∑ j, |k j|) * 1000 < 12845 :=
    of_decide_eq_true hc
  have hr : (∑ j, ((k j : ℝ) ^ 2)) + ∑ j, |(k j : ℝ)| + 9 / 4 < 12.845 := by
    have := (Int.cast_lt (R := ℝ)).mpr hz
    push_cast at this
    linarith
  have h112 : (112 : ℝ) ^ 2 * (0.032 : ℝ) ^ 2 = 12.845056 := by norm_num
  have hlt : (112 : ℝ) ^ 2 * ∑ j, (q j - u' j) ^ 2 < (112 : ℝ) ^ 2 * (0.032 : ℝ) ^ 2 := by
    rw [h112]; linarith
  exact lt_of_mul_lt_mul_left hlt (by positivity)

end Mub
