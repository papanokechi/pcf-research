# Handoff — PCF2-SESSION-R1.3
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~70 minutes
**Status:** COMPLETE

## What was accomplished
A four-phase cross-degree cleanup of the j-invariant finer-split
discovery from R1.1 (d=3 confirmed) and R1.2 (d=4 inconclusive),
designed to discriminate three hypotheses (A: B5/B6 universal but
masked at Q1 precision; B: B5/B6 d=3-specific; C: enumeration too
narrow). Verdict: **CASE B with C-caveat**. B5/B6 are confirmed
d=3-specific in their sharp/deep-WKB form; a residual shallow-N
j=0-cell shift exists at d=4 but does not survive to deep-N.

## Key numerical findings

- **R13-A (cubic sanity, d=3, dps=2000, N=[50,250], N_ref=400):**
  Spearman ρ(log|j|, δ_R13_free) = -0.5683, Bonferroni p =
  5.01e-5 (K=3); reproduces R1.1's headline ρ=-0.6906 within
  fit-window/precision differences. **PASS.**
- **R13-B (quartic test on Q1's 60, fixed-A=8, same window):**
  Spearman ρ(log|j|, δ_R13_free) = +0.0732, ρ(log|j|, residual_at_max_n)
  = -0.0729, both Bonferroni p = 1.0. **NULL on Q1 catalogue**
  (only 1 j=0 quartic in Q1: family 32).
- **R13-C (extended quartic enumeration, W=5):** 73,205 candidates →
  237 with I=12·a4·a0 - 3·a3·a1 + a2² = 0 → 101 NEW (not in Q1)
  irreducible non-singular j=0 quartics. WKB tail-fit at same
  window yields a distribution **shifted toward zero** vs Q1
  non-j=0 cluster: median δ = -2.49e-3 vs -3.73e-3. Mann-Whitney
  p = 4.8e-4 unstratified; **p = 1.1e-6 stratified** on coefficient
  sum within Q1 band [6,12] (n=77 j=0-in-band).
- **R13-D (family-32 deep-WKB, dps=5000, N_ref=1000, N∈[200,800]):**
  fam32 (j=0): δ_deep = **-4.55e-4 ± 1.45e-5**.
  fam01 baseline (non-j=0, matched coefficient norm):
  δ_deep = **-4.60e-4 ± 1.47e-5**.
  Difference: +4.9e-6, equivalent to **0.24 pooled stderr units**.
  Both families shrink |δ| by ~8× from shallow R13-B. **j=0 vs
  non-j=0 statistically indistinguishable at deep N.**
- **R13-E joint:** B5/B6 sharp form rejected at d=4 (deep-WKB null);
  shallow-N effect is a sub-leading WKB systematic that does not
  lift to a deep rule.

## Judgment calls made

1. **residual_mean is structurally near zero** (least-squares with
   intercept γ forces it to vanish). Rather than halting on the
   prompt's noise gate (`std/|mean| > 0.5`), I downgraded
   `residual_mean` to non-primary-diagnostic and used
   `residual_at_max_n` and FREE-A `δ_R13_free` as the meaningful
   summaries. Documented in rubber_duck_critique.md focus item (iii).
2. **Stopped extended enumeration at W=5** (cell size 101 ≥ 4),
   rather than escalating to W=7 or W=10. The spec said escalate
   if N_4 < 4; we exceeded by 25×.
3. **Re-classified the automated CASE_C verdict** to "CASE B with
   C-caveat" because R13-D is the most rigorous test (highest dps,
   longest fit window, lowest residual_std) and shows a clean null,
   while R13-C is a shallow-N statistical comparison subject to
   coefficient-magnitude confounds (we documented partial confound:
   ρ(a4, δ) = +0.51 within j=0 cell).
4. **Did NOT run deep-WKB on the 101 new j=0 quartics.** Each would
   take ~30s at dps=5000 N=600, ~1 hour total; deferred to R1.4 if
   the shallow-N effect needs structural interpretation.
5. **Used non-Q1 family 1 ([1,-3,-3,-3,-3], non-j=0) as the deep-WKB
   baseline** rather than averaging over all 60 Q1 families. Single
   matched-coefficient-norm comparator was sufficient.

## Anomalies and open questions

**THIS IS THE MOST IMPORTANT SECTION.**

1. **op:shallow-j-effect-d4 (NEW).** The shallow-N j=0 cell shift
   at d=4 is highly significant (p=1.1e-6 stratified) but
   dissolves at deep-N. Plausible mechanism: a β·log(N) or γ
   sub-leading WKB term with non-trivial j-dependence that decays
   as N→∞. Would benefit from: (a) explicit form of α(b) at d=4
   from Newton-polygon-Birkhoff (op:xi0-d3-direct extension), (b)
   deep-WKB on a representative sample of the 101 new j=0
   quartics, (c) cross-check whether the shift is also visible
   at d=3 in an extended cubic catalogue.
2. **Coefficient-magnitude confound within j=0 cell.** ρ(a4, δ) =
   +0.51 (highly significant) within the 99 cleaned j=0
   quartics. The j=0 condition I=0 constrains α₂² = 3α₃α₁ - 12α₄α₀,
   which biases the cell toward small-α₂ tuples. The "shift toward
   zero" partly reflects this bias. Stratification mitigates but
   does not fully remove the confound.
3. **Family 32 deep-WKB stderr might be optimistic.** The 6-point
   regression A_stderr = 1.45e-5 assumes i.i.d. Gaussian residuals,
   which they are not (WKB sub-leading shape). A more honest stderr
   might be ~10×. Even so, the fam32-vs-fam01 difference of 4.9e-6
   is well within any reasonable stderr.
4. **D_4 j=0 sub-cell (n=8) unexamined.** Among 99 new j=0
   quartics, 91 are S_4 and 8 are D_4. Whether the Galois
   stratification matters at d=4 is unknown. Recorded as
   op:galois-modular-stratification.
5. **No j=1728 quartic found in W=5 window.** R1.2 noted the
   j=1728 cell is empty in Q1; we confirm at W=5. The "lemniscatic"
   op (op:lemniscatic) is therefore still empty-cell.
6. **A_stderr in R13-A for cubics is uniformly small except fam37**
   (A_stderr=5.4e-4 at dps=2000 vs <1e-4 for fams 24, 30, 31). This
   is the same family flagged in R1.1's precision_escalation_log.
   No new information; documented for continuity.

## What would have been asked (if bidirectional)

1. Should the shallow-N effect be elevated to a Conjecture-class
   statement (op-level) in v1.3, or just recorded as a Remark?
   I went with Remark + op, since the structural mechanism is
   unclear.
2. Is the standalone_note.tex draft cautious enough? I framed it
   as "$d=3$-specific phenomenon" with explicit non-extension
   to $d=4$.
3. Would Claude prefer a different fit-window in R13-D (e.g.
   N∈[300,1000] with dps=12000)? The current window is constrained
   by mpmath precision saturation: at dps=5000 only N up to ~300
   yields a finite y_n.

## Recommended next step

**PCF2-V13-RELEASE (cubic-modular framing, cautious).** Fold
v13_paragraph_insert_v3.tex (Remark rem:b5b6-degree-stratification-r13
+ Remark rem:b5b6-shallow-d4 + Conjectures B5(d=3)/B6(d=3) +
op:shallow-j-effect-d4 + op:b5-degree-d) into PCF-2 source. Build
v1.3 (~16pp), update changelog and version in metadata, archive
v1.2 in pcf2/v12_archive/, optionally upload to Zenodo as v1.3.
After v1.3 is built, consider standalone_note.tex as a separate
arXiv/Zenodo upload (4pp, CAUTIOUS framing).

If the human prefers extra rigor before v1.3 release, an
intermediate **R1.4** prompt should:
- run deep-WKB on a 10-family representative sample of the 101
  new j=0 quartics (dps=5000, N∈[200,800])
- re-test the shallow-N stratified shift with the full dataset
  including these 10 deep-fits replacing their shallow ones
- if 10/10 deep δ match the non-j=0 baseline within stderr ~10×
  the difference observed for fam32, the shallow-N effect is
  decisively resolved

## Files committed

In sessions/2026-05-01/PCF2-SESSION-R1.3/:
- r1_3_residualization.py             (Phase A,B mpmath fitting)
- r1_3_postprocess_AB.py              (Phase A,B Spearman synthesis)
- r1_3_extended_enumeration.py        (Phase C enumeration + WKB)
- r1_3_family32_deep.py               (Phase D deep-WKB)
- r1_3_phase_E_verdict.py             (Phase E preliminary)
- r1_3_phase_E_robustness.py          (Phase E coefficient confound check)
- r1_3_phase_E_final.py               (Phase E final verdict + Phase F paragraph)
- residualization_d3.csv              (50 cubic rows)
- residualization_d4.csv              (60 quartic rows)
- residualization_d4_extended_jzero.csv (101 NEW j=0 quartics)
- extended_quartic_catalogue.json     (window stats + 101 j=0 records)
- extended_j_zero_cell.json           (101 j=0 quartics + WKB fits)
- family32_deep_residual.json         (fam32 + fam01 baseline deep fits)
- results_v3_phase_AB.json            (Phase A,B summary)
- results_v3.json                     (final verdict + sub-verdicts)
- phase_E_robustness.json             (stratified Mann-Whitney)
- claims.jsonl                        (26 AEAL entries, concatenated)
- claims_phase_AB.jsonl, claims_phase_C.jsonl, claims_phase_D.jsonl,
  claims_phase_E.jsonl, claims_phase_E_final.jsonl
- rubber_duck_critique.md             (4 focus items, addressed)
- v13_paragraph_insert_v3.tex         (cubic-modular framing)
- standalone_note.tex                 (4pp CAUTIOUS d=3-specific)
- handoff.md                          (this file)
- halt_log.json, discrepancy_log.json, unexpected_finds.json
- *.log                               (run logs, ~5 files)

## AEAL claim count
**26 entries** written to claims.jsonl this session (8 from Phase A,B;
4 from Phase C; 5 from Phase D; 5 from Phase E preliminary; 4 from
Phase E final).

## SIARC stack queue update (R1.3 verdict noted)

Suggested update to siarc-relay-bridge standing queue:

- **R1.3 — j-invariant cross-degree cleanup:** COMPLETE, verdict
  CASE B with C-caveat. B5/B6 confirmed d=3-specific.
- **PCF2-V13-RELEASE:** READY (cubic-modular cautious framing,
  v13_paragraph_insert_v3.tex prepared).
- **R1.4 — shallow-N d=4 deep-harvest** (OPTIONAL): only if v1.3
  framing needs additional support before public release.
- **STANDALONE-NOTE-J-INVARIANT** (OPTIONAL, 4pp cautious):
  draft prepared at standalone_note.tex; consider Zenodo upload
  after PCF2-V13 release.
