# Handoff — PCF2-SESSION-A2 (CONDUCTOR-7 ANCHOR VERIFY)

**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~25 minutes
**Status:** HALTED (4 of 5 checks pass; check 5 flagged as informative halt — see below)

## What was accomplished
Verified family 46 of the PCF-2 cubic-family catalogue,
`b(n) = n^3 - 2n^2 - n + 1`, as a candidate fourth calibration anchor for
the `+_C3_real` bin tied to the conductor-7 real cyclic cubic field
`K_7 = Q(zeta_7 + zeta_7^-1)`.  Discriminant, Galois group, splitting
field, and high-precision continued-fraction value all confirm the
expected structure.  The convergence-rate check fits `A ~ 5.95` (rounding
to 6), which falls outside the prompt's literal expected range
`A in {3, 4}` and therefore triggers a halt — but this is itself a
likely-positive scientific finding (degree-3 extension of PCF-1 v1.3
Theorem 5: `A = 2d`).  The anchor was **not** appended to
`calibration_anchors.json` pending Claude's decision.

## Key numerical findings
- **Δ₃(b) = 49 = 7²** (exact, sympy). Family 46 catalogue value confirmed.
- **b(n) irreducible over Q**; sympy `galois_group` returns
  `PermutationGroup([(0 1 2)])` (cyclic of order 3), and Δ₃ is a perfect
  square, so **Gal(b/Q) = C₃** unconditionally.
- **Splitting field identification**: substitution `n = 1 + x` gives
  `b(1+x) = x³ + x² − 2x − 1`, exactly the minimal polynomial of
  `η₁ = 2 cos(2π/7) = ζ₇ + ζ₇⁻¹`. Hence
  `Q(α_46) = Q(ζ₇ + ζ₇⁻¹) = K_7`, the real subfield of `Q(ζ_7)`,
  conductor `f = 7`, class number `h⁺(K_7) = 1`, totally real.
- **L₄₆ to 50 digits** (mp.dps = 80, n_max = 2000):
  `0.53825867319414640868987946895865319907773601719237`.
  Stable to 60 digits already by n = 200; `|L₂₀₀₀ − L₁₀₀₀| = 0` at
  dps=80 (script: `conductor7_verify.py`, run: `run.log`).
- **PSLQ T₁** `{1, log 2, log 3, log 7, π, √7, ζ(3)}`: NO relation found
  with `maxcoeff ≤ 10¹²`, `tol = 1e-50`.
- **PSLQ T₂** `{1, η₁, η₁²}` (Z-basis of `O_{K_7}`): NO relation found
  with `maxcoeff ≤ 10¹²`, `tol = 1e-50`. So L₄₆ is *not* algebraic of
  small height in `K_7` to 50-digit accuracy.  (We deliberately replaced
  the redundant 4-tuple `{1, η₁, η₂, η₃}` — which has the trivial
  relation `1 + η₁ + η₂ + η₃ = 0` — with `{1, η₁, η₁²}`.)
- **Convergence-rate fit** `log|δ_N| = −A N log N + α N + c` on
  `N ∈ [10, 85]` at dps=800 yields `A = 5.953348`, `α = 5.730870`,
  `c = −6.402550`. Round(A) = 6.

## Judgment calls made
1. **PSLQ basis T₂ choice.** The prompt suggested
   `{1, 2cos(2πk/7) : k=1,2,3}`, which is Z-linearly dependent
   (`η₁ + η₂ + η₃ = −1`).  PSLQ on this 4-tuple always returns the
   trivial basis-internal relation `[0, −1, −1, −1, −1]` with
   coefficient 0 on L. I replaced it with the ring-of-integers Z-basis
   `{1, η₁, η₁²}` (`O_{K_7} = Z[η₁]` since 7 is the conductor and h⁺ = 1).
   Recorded both the original spurious relation and the corrected
   probe. The corrected probe also returned no relation.
2. **Halt vs discrepancy on check 5.** The prompt's halt list (verbatim)
   does NOT include "A out of range {3,4}", but the rubric line
   "all must pass; any failure -> HALT" does. I escalated to a halt
   to respect the literal contract, did **not** append the anchor, and
   flagged the finding as a positive degree-3 extension for Claude's
   review.
3. **Convergence fit precision.** First run at dps=200 had `δ_N`
   underflowing to 0 by N≈35, so the LSQ collapsed.  Bumped to dps=800
   with N range `[10, 85]` and skipped exact-zero points; the resulting
   fit is stable across reruns.

## Anomalies and open questions
**THIS IS THE MOST IMPORTANT SECTION.**

- **A = 2d generalisation (PCF-1 v1.3 Theorem 5 → degree 3?).**
  Theorem 5 ("WKB exponent identity") in `p12_pcf1_main.tex` is stated
  and verified for degree-2 PCFs only, with `A ∈ {3, 4}` (`A = 4` for
  `V_quad`, `A = 3` for the five QL families).  For cubic family 46
  we observe `A = 5.953348 → 6 = 2d`, consistent with the natural
  conjecture
  `A = 2d` for a degree-d PCF.  Recommendations:
  - Recompute the fit with smaller N (say `[8, 60]`) and a higher-order
    correction (subleading `−β_w log N` term) to refine `A` to ≥3
    decimal places.  Current 1-figure agreement with `A = 6` (5.95
    vs 6) is good but not unambiguous.
  - Check at least one other `+_C3_real` cubic family from the
    catalogue (e.g. families with `Δ₃` square > 0) to see whether the
    `A = 6` observation is generic or family-specific.
  - If confirmed across families, consider a Theorem 5.bis (cubic
    case) with conjectural `α = 2d − 2 log c_b + log|c_a|`.  Our
    `α_obs ≈ 5.731` for `b(n) = n³ + ...` (so `c_b = c_a = 1`) does NOT
    match `2·3 − 0 + 0 = 6` exactly — off by ≈0.27.  Either the
    closed form needs revisiting for cubics, or we need a
    higher-order fit.

- **L₄₆ transcendence question.** PSLQ found no relation in either
  basis at 50-digit precision with `maxcoeff = 10¹²`.  This is
  consistent with L₄₆ being transcendental (or algebraic of height
  > 10¹²).  No anomaly, but worth noting that L₄₆ does **not** lie
  in `O_{K_7}` — the field anchoring is a Galois-symmetry statement,
  not a value-in-field statement.

- **Pre-screen discrepancy.** Per the SIARC standing rule
  "any pre-screen discrepancy > 5% of stratum count → halt": A_fit
  deviates from the predicted upper bound of 4 by `(5.95 − 4) / 4 = 49%`.
  This is well above 5% and reinforces the halt decision.

## What would have been asked (if bidirectional)
- Should the anchor be appended even with `A = 6` (degree-3 scaling)
  rather than `A ∈ {3, 4}`?  My read is **yes** — algebraic checks 1–4
  are airtight, and check 5 reveals a positive scientific finding —
  but I held off pending Claude's call.
- Is the prompt's `A ∈ {3, 4}` expectation a transcription error
  (it cites Thm 5 which is degree-2-only), or did the prompt author
  have an unwritten degree-3 prediction in mind?

## Recommended next step
**PCF2-SESSION-A3 — DEGREE-3 WKB EXPONENT SCAN**: refit `(A, α, β_w, γ_w)`
for *all* `+_C3_real` cubic families in the catalogue (3 entries:
families 46, and the two siblings reported in
`galois_distribution_summary.json`), at dps=2000, K∈{12, 16}, with
4-parameter fit. If `A = 6` and `α` matches a closed-form depending
on `(c_a, c_b, β, γ, a_3, a_2, a_1, a_0)` to ≥10 digits, propose a
cubic-PCF analog Theorem 5.bis.  Concurrently, append family 46 as
the 4th calibration anchor (with a footnote explaining the degree-3
A scaling) to unblock the catalogue.

## Files committed
- `conductor7_verify.py` — full verification script (sympy + mpmath).
- `run.log` — full stdout/run trace.
- `results.json` — structured outputs of all 5 checks.
- `claims.jsonl` — 4 AEAL claims (Galois, field, L value, A fit).
- `halt_log.json` — single halt entry (check 5, A out of {3,4}).
- `discrepancy_log.json` — same incident, with degree-3 explanation.
- `unexpected_finds.json` — empty `{}` (no spurious algebraicity).
- `handoff.md` — this file.

## AEAL claim count
4 entries written to `claims.jsonl` this session.

## Anchor catalogue update
**Family 46 was NOT appended to `calibration_anchors.json`** because
check 5 did not pass per the prompt's literal contract.  Algebraic
characterisation is complete and correct; the deferred append is a
single-line change (replacing the existing `cyclic_C3_n3_minus_3n_plus_1`
entry's siblings with a fourth `conductor7_n3_minus_2n2_minus_n_plus_1`
entry).  Recommend Claude either (a) override and append after
reviewing the degree-3 finding, or (b) issue PCF2-SESSION-A3 first to
verify A=2d across the bin, then append.
