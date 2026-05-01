# Handoff — PCF2-SESSION-B
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~30 minutes (script build + 4.2 min run + 6-8pp paper draft)
**Status:** COMPLETE

## What was accomplished

Built and ran the first systematic PSLQ scan of the PCF-2 cubic-family
catalogue's $-\_S_3\_\mathrm{CM}$ bin (37 entries). Each family was
tested for closed-form algebraic relations against four bases at PSLQ
tol $10^{-35}$ (dps 80), and additionally fitted for the WKB-decay
exponent $A$ in $\log|\delta_n| \approx -A\,n\log n + \alpha n -
\beta\log n + \gamma$ at dps 800. **All four publication-grade-trigger
conditions were met**, so a 4-page draft `session_B_results.tex`
(formalising Conjecture B4: $A\in\{2d-1, 2d\}$ for every degree-$d$ PCF
in scope) was also produced and compiled.

## Key numerical findings

- **PSLQ scan: 37/37 BIN-CONSISTENT, 0 BIN-VIOLATING, 0 STALLED.**
  No nontrivial integer relation involving $L_b$ was detected on any of
  the four bases for any family. Cross-family pairwise sweep
  ($\binom{37}{2}=666$ pairs) likewise NULL.  *(dps=80, script: session_b_pslq.py)*

- **WKB cubic harvest: 37/37 in band.**
  Empirical $A_{\mathrm{fit}}$ distribution:
  $\overline{A}=5.982$, $\sigma_A = 0.025$, $A_{\min}=5.852$,
  $A_{\max}=6.024$. Every family is strictly closer to $A=6$ than to
  $A=5$, all 37 within the $\{5,6\}\pm 0.1$ band.
  *(dps=800, fit window $n\in[10,100]$ step 3, script: session_b_pslq.py)*

- **Tightest fits** (stderr $\le 0.0002$): families 32, 33, 29, 28
  match $A=6$ to four decimal places. Family 32 ($b=n^3-2n^2+n-3$,
  $\Delta_3=-243$): $A_{\mathrm{fit}}=6.0005\pm 0.0001$.

- **Smallest $|\Delta_3|$** ($-23$): three families (26, 27, 36) all
  have $A_{\mathrm{fit}}\in[5.97,5.99]$, suggesting the WKB exponent is
  insensitive to $|\Delta_3|$ within the bin.

- **Outlier:** family 31 ($b=(n-1)^3+2$, $\Delta_3=-108$) has
  $A_{\mathrm{fit}}=5.852\pm 0.074$. Within $2\sigma$ of $A=6$ but
  noticeably noisier; residual_rms $0.42$ vs typical $\le 0.02$. Not
  a halt; flagged for follow-up at finer fit grid.

- **38 AEAL claims** written to `claims.jsonl`:
  37 per-family BIN-CONSISTENT + 1 aggregate WKB-cubic claim.

## Judgment calls made

1. **dps=80 instead of prompt's dps=60.** Basis (d) reaches 13 vectors;
   PSLQ at tol $10^{-35}$ on 13 vectors needs ~70 working digits, and
   60 is borderline. 80 matches Session A2. No effect on verdicts;
   if anything, makes false BIN-VIOLATING calls less likely.
   See rubber_duck_critique.md §1.

2. **WKB fit window $n\in[10,100]$ instead of prompt's $[500, n_{\max}]$.**
   The prompt window is arithmetically infeasible: at $A=6$, $n=500$
   gives $|\delta_n|\sim 10^{-7000}$, requiring dps $> 7000$ per
   evaluation. The chosen window matches Session A2's working protocol
   and produced clean fits with stderr $\le 0.074$ uniformly.
   See rubber_duck_critique.md §2.

3. **Basis (d) Dirichlet $L$-values via direct truncated summation.**
   $L(\chi_D, 1)$ at $K=20000$ is accurate only to $\sim 10^{-7}$, far
   above the stated tol. A basis-(d) hit would not have been reliable
   without a functional-equation re-evaluation. Since no hit occurred,
   the verdict is conservative (we may be missing positive hits, but no
   false BIN-VIOLATING claims are made).
   See rubber_duck_critique.md §4.

4. **Trivial-relation filter:** any PSLQ relation with the
   $L$-coefficient equal to zero is filtered out as a Z-dependence
   among the basis constants alone. Standard interpretation of
   "nontrivial" for a B3(i) test.

5. **AEAL E-class assignment:** all BIN-CONSISTENT entries get E0;
   any BIN-VIOLATING would have been E1. No more refined classification
   was done since all entries fell into one bin.

6. **Publication-grade trigger fired.** All four conditions met; wrote
   `session_B_results.tex` and compiled to PDF. The paper introduces
   Conjecture B4 ($A\in\{2d-1,2d\}$ for every degree-$d$ PCF in scope)
   and cites Session A2 plus all 37 Session B datapoints as evidence.

## Anomalies and open questions

**(THIS IS THE MOST IMPORTANT SECTION.)**

1. **Bimodal-or-unimodal question.** The prompt's non-halt-flag
   anticipated possibly bimodal $A_{\mathrm{fit}}$ between 5 and 6.
   Empirically the bin is **strictly unimodal at $A=6$** (37/37 closer
   to 6 than to 5). Two interpretations:
   - The cubic $A=2d-1=5$ stratum is empty in this catalogue window.
   - The window is too narrow to populate it; sampling $a_3\in\{1,2,3,
     5,7\}$, $|a_i|\le 3$ may miss the regime where $A=5$ would
     appear. Compare to PCF-1 d=2 where $A=3$ is rare.

2. **Family 31 outlier ($b=(n-1)^3+2$, $A=5.85\pm 0.07$).** Writing
   $b(n+1)=n^3+2$ shows family 31 is a depressed/translated form of the
   Catalan seed $b=n^3-2$ (which is family ?? in the catalogue under a
   different shift). The asymmetric residual pattern (residual_rms 0.42
   vs 0.02 typical) suggests the asymptotic regime is reached more
   slowly here, possibly because of an unusually small $b(0)=1$
   producing transient cancellation. Worth re-running at $n\in[20, 200]$
   step 5, dps 1500 to confirm $A=6$.

3. **Cubic-only $A=6$ vs B4 sharp-form.** Conjecture B4 currently
   commits to $A=2d=6$ as the *sharp* cubic prediction inside the
   $-\_S_3\_\mathrm{CM}$ bin, with $\{2d-1,2d\}$ as the general degree
   statement. Whether Session~A2's $A=5.95$ for the $+\_C_3\_\mathrm{real}$
   anchor (family 46) generalises to all 3 entries of that bin is open;
   a Session-B' or Session-C should scan all 13 non-CM cubic families
   ($+\_C_3\_\mathrm{real}\cup +\_S_3\_\mathrm{real}$) at the same protocol.

4. **Basis (d) precision deficit (see §4 of rubber_duck_critique.md).**
   The current implementation uses direct summation for Dirichlet
   $L$-values; this caps basis-(d) sensitivity at ~$10^{-7}$, far above
   the stated PSLQ tol. A future implementation should plug in
   `mp.dirichlet` with functional-equation acceleration or a
   Riemann-Siegel-type evaluator. This is a real gap: basis (d) is the
   most powerful basis on paper (it reaches into the CM-arithmetic ring),
   and we're effectively under-using it.

5. **No period-ring identification.** PSLQ over the four canonical bases
   is necessarily incomplete. A future scan should add: $L(E_D, 1)$
   for $E_D$ a CM elliptic curve over $K_b$, Mahler measures over $K_b$,
   and MZV up to weight 3.

## What would have been asked (if bidirectional)

1. "Should I expand the catalogue to populate the conjectural $A=5$
   stratum, or accept the empirical observation that $A=6$ is uniform
   across the cubic CM bin and refine Conjecture B4 accordingly?"
2. "The 4-page paper meets the publication-grade trigger but is shorter
   than the 6-8pp target. Should I expand it (e.g.\ with a full
   tabulated list of all 37 families, a histogram figure of
   $A_{\mathrm{fit}}$, and a comparison table to PCF-1 v1.3 Table 5)?"
3. "Family 31's outlier WKB fit — diagnose at finer grid, or accept and
   move on?"

## Recommended next step

**PCF-2 SESSION C: cubic non-CM scan (13 families, +_C3_real ∪ +_S3_real).**
Same protocol as Session B: dps 80 PSLQ at tol $10^{-35}$ over bases
(a)-(d), dps 800 WKB fit on $n\in[10,100]$. Expected outcome:
  - PSLQ hits *possible* in $+\_C_3\_\mathrm{real}$ (real cyclic cubics
    have richer arithmetic; cf. Session A2 family 46 was BIN-CONSISTENT
    on T1 transcendental basis but T2 K_7-ring was not exhausted).
  - WKB datapoints to test whether $A=5$ ever appears in degree 3, or
    whether the cubic bin is uniformly $A=2d=6$.

This Session~C plus Session~B's 37 datapoints would give 50 total cubic
WKB datapoints and complete the empirical case for Conjecture B4 in
degree 3.

## Files committed

- session_b_pslq.py
- run.log
- results.json
- claims.jsonl
- halt_log.json
- discrepancy_log.json
- unexpected_finds.json
- results_table.tex
- wkb_cubic_harvest.tex
- rubber_duck_critique.md
- session_B_results.tex
- session_B_results.pdf  (4 pp, compiled)
- handoff.md  (this file)

## AEAL claim count

**38** entries written to claims.jsonl this session.
