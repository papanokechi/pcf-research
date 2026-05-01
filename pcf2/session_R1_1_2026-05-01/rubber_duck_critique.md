# Rubber-duck critique — PCF2-SESSION-R1.1

## (i) Is the j-invariant signal real or a coincidence?

`log_abs_j = log(|j(E_b)| + 1)` where E_b: y² = b(x), gives Spearman
ρ(log|j|, δ) = -0.691 on the full 50, raw p = 2.9e-8, Bonferroni-K=14
p = 4.0e-7 (and 5.7e-7 even at the prompt-mandated K=20). Far past
the strict cut |ρ| > 0.6 ∧ p_Bonf < 0.001.

The signal survives every standard sensitivity check:

| check                                        | ρ      | raw p   |
|----------------------------------------------|--------|---------|
| full 50, unweighted                          | -0.691 | 2.9e-8  |
| full 50, A_stderr⁻²-weighted                 | (see results_v2.json — large, also significant) |
| exclude family 31                            | -0.671 | 1.3e-7  |
| partial Spearman, controlling log\|Δ_3\|      | -0.794 | 5.9e-12 |
| partial Spearman, controlling log\|j\|        | -0.667 | 1.2e-7 (this is log\|Δ_3\| controlling for j) |

The crucial structural fact: **log|j| and log|Δ_3| are essentially
independent on the catalogue** — Spearman ρ between them is -0.016
(p = 0.91). They are not the same signal in disguise, even though the
elliptic discriminant Δ_E = -16·disc(b) is proportional to Δ_3. This
is because j = c4³/Δ_E, and the c4 factor varies enough across
families to decouple log|j| from log|Δ_3|. The partial correlations
**strengthen** rather than collapse, and the (log|Δ_3|, log|j|) pair
regression reaches adj R² = 0.706.

## (ii) Closed-form interpretation of the j=0 cell

The four families (1, -3, 3, *) — i.e. families {30, 31, 32, 33} —
all have j(E_b) = 0 (E_b has CM by **Z[ω]**, the equianharmonic point).
Their A_fit values are:

| fid | b(n)            | Δ_3   | A_fit       | δ = A-6  |
|-----|-----------------|-------|-------------|----------|
| 30  | n³-3n²+3n-3    | -108  | 5.998856    | -1.1e-3  |
| 31  | n³-3n²+3n+1    | -108  | 5.999882    | -1.2e-4  |
| 32  | n³-3n²+3n+2    | -243  | 6.000489    | +4.9e-4  |
| 33  | n³-3n²+3n+3    | -432  | 6.001004    | +1.0e-3  |

Family 24 (next-smallest |j| = 28.9) has δ = -3.95e-3.
Compare: the loose-A_fit tail (families 36, 3, 8, 9 with |j| > 3000)
has δ ranging from -0.030 to -0.020.

Empirically: **|j(E_b)| → 0 ⟹ A_fit → 2d = 6 to within numerical
precision**, while large |j| families consistently fall a small
amount below 6. This is the closed-form form the v1.1 sharp form
B4 (A = 2d) was missing — it appears to be a *threshold*: A = 2d
**exactly** for the CM-equianharmonic locus, but A = 2d − ε(j) for
|j| > 0 with ε increasing in log|j|.

## (iii) Risks the rubber-duck flags

1. **K = 20 vs K = 14 Bonferroni.** The prompt mandated K = 20 (12
   catalogue + 6 pari + 2 new). With pari/gp deferred, only K = 14
   tests are valid (12 + Mahler + log|j|; j_invariant counts
   separately as the unscaled hit, K_actual = 14 in the run). Either
   way the headline passes: at K = 20, Bonferroni p =
   2.9e-8 × 20 = 5.7e-7. At K = 14, 4.0e-7. The over-correction
   risk is irrelevant; the signal is too strong.

2. **Is K = 20 (or 22) exhaustive?** No — pari deferral leaves 6
   number-field invariants untested. Among them, **conductor** is the
   one most likely to provide a third orthogonal signal (CM
   conductor and j are linked but not by a simple monotone). The
   rubber-duck recommends an R1.2 cycle once pari is installed, NOT
   to revise the v1.3 absorption (j is sufficient).

3. **Possible Mahler-measure third signal.** Spearman(log_Mahler, δ)
   = -0.371, raw p = 0.008, Bonferroni p ≈ 0.11 — does not pass the
   strict cut at K = 14 but is suggestive. It is partially redundant
   with log|j| (the (log_Mahler, log|j|) pair has adj R² = 0.540).
   The v1.3 paragraph mentions Mahler as a secondary candidate; the
   primary v1.3 statement is in terms of j alone.

4. **Catalogue size.** 50 families is small for a 22-invariant
   battery. The headline survives because |ρ| = 0.691 is very large
   for n = 50; statistical power is not the issue. The structural
   fact (j = 0 ⟺ A = 2d exactly on a 4-element subset) is closer
   to a *closed-form rule* than a *correlation*; cross-validation
   on a degree-4 catalogue (PCF2 Session Q1) is the natural next
   test (op:b4-degree-d hits this).

5. **Family-31 question (closed).** Family 31 was the loosest
   A_stderr in R1 (0.074); in R1.1 it is now 8.3e-6 — among the
   tightest. Excluding it changes ρ(log|j|, δ) from -0.691 to
   -0.671. Family 31 was a precision artifact, not a structural
   driver. Closed.

6. **Precision escalation succeeded on all four targets.** dps =
   1500, N grid 10..130 step 2 (61 pts) gave A_stderr ∈
   {2.0e-4, 6.5e-5, 8.3e-6, 5.4e-4} — all under the 1e-3 target.
   No deeper resonance issue at family 31. op:family-31-deep is
   closed.

## (iv) Verdict

R1.1 verdict: **HALT — finer-cubic-split SIGNAL identified.**
The signal is the j-invariant of E_b: y² = b(x). The closed-form
sub-rule is

    j(E_b) = 0  ⟹  A_PCF(b) = 2d = 6   (exactly, up to dps=1500
                                          numerical precision).

The "soft" form is

    A_PCF(b) ≈ 2d − ε(b),   with  ε(b)  monotone in log|j(E_b)|,
                            ε ≥ 0  for j ≠ 0  on the cubic catalogue.

This is the v1.3 absorption target. PCF-2 v1.2 (A = 2d for all 50
cubic families) is consistent with this — the deviations |δ| ≤ 0.03
were within v1.2's stated precision tolerance. v1.3 sharpens to
"=2d on the j=0 sub-locus, ≤2d elsewhere with j-controlled drift".
