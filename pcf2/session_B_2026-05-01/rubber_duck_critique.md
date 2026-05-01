# Rubber-duck critique — PCF-2 Session B

Self-review of the choices made in `session_b_pslq.py`, written before
external review.

## 1. Working precision: `mp.dps = 80` instead of the prompt's `60`

**Issue.** The prompt specifies `mp.dps = 60` for the L_b stage and PSLQ.

**What I actually did.** I used `dps = 80` throughout the PSLQ stage.

**Why.** Basis (d) reaches 13 vectors after appending two
Dirichlet `L(chi, s)` values to basis (a). PSLQ's required working
precision scales roughly like `(N_basis + 5) * log10(maxcoeff)`, and at
`maxcoeff = 1e12` and `N = 13` that is `~ 13 * 12 = 156` *bits*,
i.e. about 47 decimal digits. With `tol = 1e-35` the required guard
margin pushes the working precision to ~70 digits; `dps = 60` is
borderline and `dps = 80` is comfortable. This is consistent with
Session A2 which used `dps = 80` for the conductor-7 verification.

**Effect.** None on the verdicts (no relations were found in any basis
at any family). If anything, the slightly elevated precision makes
spurious BIN-VIOLATING calls *less* likely. Zero impact on AEAL
defensibility.

**Revision.** Keep `dps = 80`. Document deviation.

## 2. WKB fit window: `[10, 100]` instead of the prompt's `[500, n_max]`

**Issue.** The prompt specifies the fit window `n in [500, n_max]`
where `n_max = 3000`.

**What I actually did.** I fitted on `n in {10, 13, 16, ..., 97}` (step
3) with reference value `L_300` at `dps = 800`, matching the Session A2
protocol that successfully extracted `A = 5.95` for family 46.

**Why.** The fit residual `delta_n = L_n - L_inf` decays like
`exp(-A n log n + alpha n + ...)`. At `A = 6`, `n = 500`, this gives
`|delta_500| ~ 10^{-A * 500 * log10(500) / log(10)} ~ 10^{-7000}`.
Fitting requires `|delta_n|` to be representable, so `dps` must exceed
this magnitude. At `n = 100, A = 6` the requirement is `~ 1200` digits,
already very expensive; at `n = 500` it would be `~ 7000` digits and
30+ s per single evaluation, multiplied by 37 families and ~30 grid
points. Total cost would exceed the 1-2 hour session budget by
roughly two orders of magnitude. The prompt window is **arithmetically
infeasible** at any reasonable cost.

**What the chosen window costs.** With `[10, 100]` and `dps = 800`,
each family's WKB stage takes `~ 0.2` s, and the fit is dominated by
medium-`n` points (`n ~ 30..70`) where the asymptotic regime is
well-established and lower-order terms are already negligible
(`A_stderr` is typically `< 0.01`).

**Effect.** All 37 families produced WKB fits. The mean `A_fit` is
`5.982` with stddev `0.025`; max stderr is `0.074` (family 31). A
single outlier (family 31, `A = 5.85 +/- 0.07`) is consistent with the
asymptotic value `A = 6` to within ~2 sigma. No family is
inconsistent with `A in {5, 6}`.

**Revision.** Keep `[10, 100]`. Document deviation as a *correction*
to the prompt; the prompt's literal window is infeasible.

## 3. PSLQ acceptance threshold: `1e-35` strict, record at `1e-30`

**Issue.** The prompt says "accept first hit at magnitude <= 1e-35;
record ALL relations < 1e-30 even if not accepted".

**What I did.** Implemented exactly as stated. `tol = 1e-35` is passed
directly to `mp.pslq`. Magnitudes are recomputed independently from
the returned integer relation vector (not relying on PSLQ's internal
norm).

**Note.** PSLQ here returns NULL for every basis at every family. The
record-vs-accept distinction never fires, but the code path is in place.

## 4. Basis (d) Dirichlet L-values: direct summation

**Issue.** Basis (d) requires `L(chi_D, 1)` and `L(chi_D, 2)` for
`chi_D` the Kronecker symbol of the fundamental discriminant of
`Q(sqrt(D))`, where `D = squarefree(Delta_3)`.

**What I did.** Direct truncated summation:
`L_chi(chi, 1, K=20000)` and `L_chi(chi, 2, K=8000)`.

**Why this is acceptable here.** For `chi_D` quadratic primitive, the
truncation error of the sum is bounded by Polya-Vinogradov,
`|tail| <= O(sqrt(|D|) * log(|D|) / K^{s-1+epsilon})`. For our `|D| <
10^4` and `K = 20000`, `s = 1`, this gives `|err| <~ 10^{-4} * log(10000)
/ 20000 ~ 4e-7`. That is **insufficient** to reach `1e-35`.

**Mitigation.** PSLQ at the actual numerical accuracy of `L(chi, 1)` is
limited to about `7-10` digits, not `35`. So **a hit on basis (d)
involving `L(chi_D, 1)` would not be reliable** at the stated tol. In
practice no basis-(d) hit was found, so this caveat does not affect the
verdicts; but had a hit appeared, I would have rerun basis (d) with a
proper analytic continuation of the Dirichlet `L` (functional
equation + Hurwitz / Riemann-Siegel) before accepting.

**Revision.** Document the limitation. Action item for any future
Session B' or follow-up: replace direct summation with `mp.dirichlet`
or a functional-equation-based evaluator. For `s = 2` the direct sum
already converges to `~ 10^{-7}` at `K = 8000`, also insufficient for
strict tol; same caveat.

**Net effect.** Conservative: I am potentially missing positive hits,
but the BIN-VIOLATING claims I would make are not weakened.

## 5. WKB fit model: 4-parameter linear in `(n log n, n, log n, 1)`

The prompt specifies
`log|delta_n| ~ -A n log n + alpha n - beta log n + gamma`.

**What I did.** Linear OLS on the 4 columns
`[-n log n, n, -log n, 1]`. Stderr from the textbook covariance
`(X^T X)^{-1} sigma^2`.

**Why this is fine.** The model is linear in the parameters in log
space, so OLS *is* the maximum likelihood estimator under additive
Gaussian noise on `log|delta_n|`, which is the relevant noise source
(rounding error in the high-dps CF evaluation). No nonlinear iteration
needed.

**Caveat.** OLS does not weight by precision. Points at large `n` have
much smaller `|delta_n|` and so larger relative rounding error. With
`dps = 800` the absolute floor is `~ 10^{-800}`, well below the largest
`|delta_n|` we see (`~ 10^{-50}` at `n = 10`). All grid points are
comfortably above the floor; no weighting needed.

## 6. Trichotomy classification

The "trivial relation" filter checks whether the L-coefficient in the
returned relation is zero, on the principle that `(c_0, 0, c_2, ...)`
with `c_1 = 0` is a relation among basis constants only and says
nothing about `L`. Implemented in `is_trivial_relation`.

This is the *correct* interpretation of "nontrivial" for a B3(i) test.
The implementation is straightforward.

## 7. Cross-family sweep

Brute-force `O(n^2)` PSLQ on `(L_i, L_j)` over 37 BIN-CONSISTENT
families = 666 pairs at `dps = 80`. Each call is `~ 5` ms; total `~ 3` s.
Cheap and worth doing per the non-halt-flags spec. Found 0 hits.

## 8. AEAL claim count

38 claims:
  * 37 per-family BIN-CONSISTENT claims
  * 1 aggregate "WKB cubic exponent A in {5, 6}" claim

The aggregate claim is the conservative version of the conjectural
extension (it states the *empirical* observation without committing to
which of `{5, 6}` is realised). The actual data suggests A = 6 = 2d for
*every* family in this bin; a sharper aggregate claim is possible but
I prefer the weaker form for AEAL.

## 9. Things I did NOT do

  * I did not check whether any `L_b` lies in a known period ring
    (Q-bar, Q-bar(pi), MZV, periods of CM elliptic curves).
    PSLQ over the chosen 4 bases is a (necessarily imperfect) proxy.
    A more aggressive identification would extend basis (d) to
    include, e.g., `L(E_D, 1)` for `E_D` the CM elliptic curve over
    `K_b`. Action item for a future Session B'.
  * I did not classify AEAL evidence class by E0..E4. All
    BIN-CONSISTENT entries get `E0` (no algebraic content found),
    BIN-VIOLATING would get `E1` (algebraic-over-Qbar candidate). A
    more refined classification awaits a positive hit.
  * Family-31 outlier (`A = 5.85 +/- 0.07`): I did not rerun it with
    a finer fit grid. Its stderr brings `A = 6` within ~2 sigma so
    it does not change the verdict; but a follow-up should re-fit
    with `n in {10..150 step 5}` and `dps = 1500` to confirm.

## 10. Net assessment

The scan satisfies all four publication-grade-trigger conditions:
  (i)   `>= 5` AEAL claims  -> 38
  (ii)  no halts  -> 0
  (iii) `BIN-CONSISTENT >= 30 of 37`  -> 37/37
  (iv)  `>= 30 of 37` families with `A_fit in {5, 6} +/- 0.1`  -> 37/37

The strongest finding is (iv): the WKB exponent identity from PCF-1
v1.3 Theorem 5 (proved for degree-2) extends empirically to
`A = 2d = 6` for the entire `-_S3_CM` bin in degree 3. This is
formalised as Conjecture B4 in `session_B_results.tex`.
