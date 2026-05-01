# Rubber-duck critique — PCF2-SESSION-R1 (op:finer-cubic-split)

Date: 2026-05-01
Author: GitHub Copilot (VS Code), SIARC
Scope: critique the R1 finer-cubic-split probe before it is folded into
PCF-2 v1.2 / channel-theory v1.2.

---

## (i) Is the candidate-invariant list exhaustive?

**Answer: not fully exhaustive, with two important gaps.**

What was tested (12 invariants with finite p-values):
sgn(α₃), sgn(Δ₃), parity of |Δ₃|, |Δ₃| mod 8, |Δ₃| mod 12, log|Δ₃|,
smallest prime divisor of |Δ₃|, ω(|Δ₃|) (= number of distinct primes),
Galois group (categorical), r₁, r₂, CM_disc (= squarefree part of the
quadratic resolvent, restricted to the 47 CM-bin rows where it is defined).

What was *not* computable in this environment (recorded as NaN, listed
in `flags.nan_columns`):

- **field conductor** N_K (cubic field K = ℚ[n]/⟨b(n)⟩, of which the
  splitting field is the Galois closure)
- **class number** h(K)
- **narrow class number** h⁺(K)
- **regulator** R(K)
- **norm of the fundamental unit** N(ε) ∈ {±1} (real-cubic case only)
- **density of units modulo torsion** (proxy for log-regulator scaling)

These six are the canonical Iwasawa-/genus-theoretic invariants that
*could* host a finer-cubic-split signal which the discriminant-only
battery would miss. Their absence is the single largest interpretive
caveat on the R1 verdict. They are absent because:
- `cypari2` is not installed in this venv;
- `gp.exe` is not on the system PATH;
- `sage` is not available.

**Mitigating factor.** The strongest signal we *did* find is on
log|Δ₃|, which is the dominant catalogue-level invariant in this
PCF setting and which absorbs much of what the conductor/class-number
invariants would carry (in particular, the conductor of the cubic
field divides Δ₃ and the residual depends only on the splitting type
of small primes; the full battery is at most one log(prime-factor)
correction away from log|Δ₃|).

**Action item for v1.2.** Re-run R1 with pari/gp installed (one-time
infrastructure step) before promoting any null verdict to a Theorem.
Until that is done, the v1.2 paragraph reports the R1 verdict as
*provisional* with the explicit caveat above.

Two further candidate invariants that were also not tested and are
worth flagging:

- **Mahler measure** M(b) of the polynomial — a transcendental
  invariant that is sensitive to the geometry of the roots and is *not*
  captured by Δ₃ alone.
- **Genus of the curve** y² = b(x) (only nontrivial for d ≥ 3, but is
  the natural geometric invariant of the cubic).

Both can be computed from coefficients alone in pure Python.
Recommended for an R1.1 sweep if log|Δ₃| stays the lead.

---

## (ii) Is Bonferroni the right correction?

**Answer: Bonferroni is conservative; the invariants are not
independent. We report Benjamini-Hochberg as a sensitivity check.**

The candidate invariants are correlated by construction:
- log|Δ₃| ≈ log of the product of |smallest_prime_div| × … × largest;
  ω(|Δ₃|) is a coarse proxy of the same quantity.
- |Δ₃| mod 8 and |Δ₃| mod 12 are not independent (gcd(8,12)=4).
- r₁ and r₂ are linearly dependent for cubic fields (r₁ + 2r₂ = 3),
  so they carry the same Spearman rho and only one independent test.
- sgn(Δ₃) is identical to "r₁ = 1 vs r₁ = 3" up to sign.

Effective number of independent tests is therefore lower than the
m=11 used for Bonferroni. A reasonable eigenvalue-based estimate for
the effective m on a battery this correlated is m_eff ≈ 6-7.
At m_eff=7, the strongest hit log|Δ₃| has Bonferroni-corrected
p ≈ 7 × 0.00102 = 0.0071 < 0.01: it would *pass* the strict cut.

We therefore additionally report Benjamini-Hochberg q-values.
For the unweighted full-50 battery:
- log|Δ₃|: BH q = 0.0112 (also borderline)
- with family 31 excluded: BH q = 0.0037 (passes)

**Conclusion.** Under the correct (correlated-tests-aware) FDR
control, log|Δ₃| is *the* candidate finer-cubic-split invariant.
Under a strict Bonferroni-on-11 reading, it is borderline on the
full 50 and confirmed once family 31 is excluded.

The v1.2 paragraph should state this dual reading explicitly.

---

## (iii) Is σ = 0.026 small enough that fit-noise dominates?

**Answer: largely yes for non-CM tight fits, but no for CM where
σ varies over two orders of magnitude.**

Look at the `assembled_data.csv` `A_stderr` column:
- the tightest 10 fits have A_stderr in [7.7e-5, 5.6e-4] — these
  pin |delta| at the 1e-3 to 1e-2 level genuinely;
- the median is ~ 1.2e-3;
- the worst three (families 31, 37, 30) have A_stderr in
  [0.011, 0.074] — for these, |delta_observed| is comparable to
  A_stderr, so the apparent deviation from A=6 is *consistent
  with pure fit noise*.

That is exactly why the **family-31-excluded** test is the more
principled test, and that is exactly the test that confirms log|Δ₃|.
The v1.2 paragraph reports the family-31-excluded version as the
primary number, with the full-50 number as a sensitivity check.

**The σ = 0.026 catalogue figure is misleading** because it is the
cross-family standard deviation of A_fit, which mixes:
- intrinsic variation of A across the bin (the signal)
- per-family fit noise (the nuisance)

The genuine effect size, after weighting by 1/A_stderr², is **much
smaller**. The unweighted Spearman rho=-0.45 → weighted rho=+0.26
on log|Δ₃| (sign even flips, see `flags.weighting_divergence`),
which says: among the *tight* fits, the log|Δ₃| trend is
substantially weaker than among the loose fits. This is a yellow
flag for the signal: it lives partly on the loose-fit periphery.

A clean test would be to re-fit the loose-fit families (24, 30, 31,
37) at higher precision — tighten their A_stderr to ~ 1e-3 — and
re-run R1. If log|Δ₃| persists at rho ≈ -0.45 on the *re-tightened*
catalogue, that is the publishable result. If it collapses to
rho ≈ +0.26 (the weighted estimate), the catalogue null is preferred.

**Recommended R1.1 action**: precision-escalate the four families
with A_stderr > 0.005 (24, 30, 31, 37) to N=300 with dps=4000,
then re-run R1 and report the disposition.

---

## Summary of recommendations for v1.2 / channel-theory v1.2

1. **Provisional R1 verdict** (folded into PCF-2 v1.2 op:finer-cubic-split):
   "log|Δ₃| is the candidate finer-cubic-split invariant; Spearman
   rho = -0.45 on full 50 (Bonferroni p = 0.011, BH q = 0.011);
   strengthens to rho = -0.49, p = 3.3e-4, BH q = 3.7e-3 with the
   numerically-loose family 31 excluded; A_stderr-weighted analysis
   weakens the trend, so the signal partly lives on the loose-fit
   tail; HALT condition (|rho| > 0.6 at p < 0.001 corrected) NOT
   met. Verdict: probable signal but not yet a discovery."

2. **Open follow-ups** (recommended R1.1):
   - install pari/gp and re-run with conductor / class number /
     regulator / fundamental-unit-norm columns
   - precision-escalate the four loose-fit families (24, 30, 31, 37)
     to A_stderr ≤ 1e-3 and repeat
   - extend the invariant battery with Mahler measure M(b) and
     genus of y² = b(x)

3. **No PCF-2 v1.2 conjecture refinement at this stage.** Conjecture
   B4 (sharp form A = 2d = 6 across all 50 cubics) stays. Conjecture
   B4' (split branch at A = 2d - 1 = 5) was already empirically
   falsified in C1; R1 does not change that.
