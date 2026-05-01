# Handoff — PCF2-SESSION-R1
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~45 minutes
**Status:** COMPLETE

## What was accomplished
Implemented Phases R1-1 .. R1-5 of the finer-cubic-split probe over the
50 cubic families from PCF2 Sessions B (37 CM) and C1 (13 non-CM).
Assembled per-family dataframe with 12 catalogue-level invariants,
ran Spearman + Pearson correlations against
$\delta = A_{\mathrm{fit}} - 6$ with Bonferroni and Benjamini-Hochberg
multiple-testing correction, weighted by $A_{\mathrm{stderr}}^{-2}$
re-test, and family-31-excluded sensitivity. Produced LaTeX
correlation table, AEAL claims, rubber-duck critique, and a draft
PCF-2 v1.2 paragraph insert.

## Key numerical findings
- log|Delta_3| Spearman rho = -0.451, raw p = 1.0e-3, Bonferroni p =
  1.1e-2, BH q = 1.1e-2 over full 50 (m=11 valid tests).
  *Source: r1_correlation_probe.py*
- Excluding family 31: log|Delta_3| Spearman rho = -0.492, raw p =
  3.3e-4, Bonferroni p = 3.7e-3, BH q = 3.7e-3. Family 31 MASKS
  the signal (it is the loosest fit, A_stderr = 0.074).
  *Source: r1_correlation_probe.py*
- A_stderr-weighted Spearman peak |rho| = 0.58 on Delta_3_sign / r1 /
  r2 (binary invariants, p unreliable via Fisher-z); on log|Delta_3|
  the weighted |rho| collapses to 0.26. Signal partly lives on the
  loose-fit tail. *Source: r1_correlation_probe.py*
- Best 2-invariant rank regression: (log|Delta_3|, CM_disc) with
  adjusted R^2 = 0.218; no pair clears the 0.4 threshold. No
  composite-invariant flag. *Source: r1_correlation_probe.py*
- HALT condition (|rho| > 0.6 at Bonferroni p < 0.001) NOT met on any
  invariant in any of the three test variants.

## Judgment calls made
1. **Class number / regulator / conductor / fundamental-unit norm /
   unit density / narrow class number** were recorded as NaN, since
   pari/gp is not installed in this venv (no `cypari2`, no `gp.exe`,
   no `sage`). The relay prompt explicitly permits NaN entries with
   a follow-up flag. This is the single largest interpretive caveat;
   I have flagged it in `rubber_duck_critique.md` (i) and recommended
   an R1.1 cycle once pari is installed. *Six tests remain open.*
2. **Weighted Spearman p-values** were computed via Fisher-z on
   effective sample size; for binary invariants (Delta_3_sign, r1,
   r2) this approximation is unreliable, so I did not promote the
   weighted |rho|=0.58 on those columns to a finding. The unweighted
   continuous invariant (log|Delta_3|) is the headline.
3. **BH q-values** added on top of the Bonferroni correction the
   prompt requested, because the candidate invariants are
   non-independent (smallest_prime_div, omega, log|Delta_3| all share
   information about |Delta_3|; r1+2*r2=3 forces a linear constraint).
   This matched the rubber-duck (ii) directive.

## Anomalies and open questions
- **Borderline-not-strict signal.** Bonferroni p = 1.1e-2 on full 50
  for log|Delta_3| narrowly *misses* the prompt's strict 0.01 cut but
  passes the 0.05 cut. With family 31 excluded (a defensible
  exclusion because A_stderr = 0.074 is two orders of magnitude above
  the median), the strict cut is passed. This is a real signal
  candidate, not a clean null and not a publishable finding. The
  rubber-duck recommends an R1.1 cycle before promoting it.
- **Weight inversion.** On log|Delta_3| the unweighted rho is -0.45
  and the weighted rho is +0.26. That is a sign flip; it is the
  cleanest evidence that the loose-fit families (24, 30, 31, 37)
  carry a disproportionate share of the apparent finer-cubic-split
  signal. **Recommended action:** precision-escalate those four
  families (N=300, dps=4000) and re-run R1; if the signal persists
  after re-tightening, that is the publishable result.
- **Catalogue-level vs field-level invariants.** All 12 tested
  invariants are computable from the polynomial coefficients alone.
  The six NaN columns are field-level (require number-field
  arithmetic). It is plausible -- and the rubber-duck flags this as
  a primary risk -- that the actual finer split lives on conductor
  or class number, neither of which we tested. **Strong
  recommendation:** install pari/gp before the next finer-split
  probe.
- **Family 31** (b = (n-1)^3 + 2 = n^3 - 3n^2 + 3n + 1) deserves a
  dedicated A-fit re-extraction. It has Delta_3 = -108 (same as
  families 3 and 30), bin = -_S3_CM, but anomalously loose A_stderr.
  It is plausibly under-converged at the standard N window.

## What would have been asked (if bidirectional)
- Should I install pari/gp into the venv (one-time apt-style step)
  and re-run R1 in this same session, or defer to R1.1? I deferred,
  but a single 15-minute install would close the largest open caveat.
- Should the v1.2 paragraph insert be staged into the actual PCF-2
  source (`tex/submitted/pcf2_main.tex`) or held as a draft until
  R1.1 closes? I held as a draft (`v12_paragraph_insert.tex`) since
  the verdict is provisional.
- Should I run R1 on the smaller +_C3_real sub-bin (n=2 families)
  alone? Sample size is too small for any statistical claim, so I
  did not, but Claude may have a structural reason to want it.

## Recommended next step
**R1.1 follow-up** (one prompt, ≈1 h): install pari/gp into the venv,
re-run the data assembly with conductor / class_number /
class_number_plus / regulator / fund_unit_norm / unit_density
populated, and *concurrently* precision-escalate families 24, 30, 31,
37 to A_stderr ≤ 1e-3. Re-run the correlation battery. The
disposition of log|Delta_3| under that combined upgrade is the
finer-cubic-split verdict for PCF-2 v1.2.

If R1.1 is not affordable, the alternative is a clean R2 prompt that
extends the battery only with the two coefficient-only invariants
flagged in the rubber-duck (Mahler measure and genus of y^2 = b(x)).
That is a pure-Python session and could ship in 30 minutes.

## Files committed
- r1_correlation_probe.py
- assembled_data.csv
- results.json
- correlation_table.tex
- claims.jsonl (5 entries)
- rubber_duck_critique.md
- v12_paragraph_insert.tex
- handoff.md
- halt_log.json (halt=false, reasons=[])
- discrepancy_log.json (empty)
- unexpected_finds.json (NaN columns flagged; weighting divergence on
  six binary/coarse invariants flagged; no family-31 single-driver
  flag; no pair adj R^2 > 0.4)
- run.log

## AEAL claim count
5 entries written to claims.jsonl this session.
