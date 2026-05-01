# Handoff — PCF2-SESSION-Q1
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~25 minutes (script build + 16 s WKB harvest + Newton + paper draft + compile)
**Status:** COMPLETE

## What was accomplished

Tested Conjecture~B4 of PCF-2 v1.1 at degree $d=4$ (claim Q1-A) and the
op:xi0-degree-d conjecture from CHANNEL-THEORY-V11 at $d=4$ (claim Q1-B).
Built a 60-quartic lex-window catalogue, ran a WKB-exponent harvest at
$\mathrm{dps}=1200$, derived the universal Newton-polygon characteristic
root $\xi_0(b)=d/\beta_d^{1/d}$ from the $L_d=1-zB_d(\theta+1)-z^2$
operator, and verified it numerically across 8 representative quartics.
Both empirical claims fired, no halts. The publication-grade trigger
fired and a 3-page companion note `session_Q1_results.tex` was drafted
and compiled (`session_Q1_results.pdf`, 3 pp, 0 errors after fix).

## Key numerical findings

- **Quartic catalogue: 60 irreducible + 10 reducible-control.**
  Galois histogram $S_4:57,\,D_4:3$; signature $(2,1):59,\,(4,0):1$.
  No CM ($r_1=0$) families in window. Calibration anchors
  $x^4-2\to D_4$, $x^4+1\to V_4$, $x^4-x-1\to S_4$ all match.
  *(script: quartic_family_enumeration.py)*

- **WKB-quartic harvest: 60/60 b4_consistent_at_8 at dps=1200.**
  Mean $A_{\mathrm{fit}}=7.9544$, $\sigma_A=0.0037$, range
  $[7.9488, 7.9605]$. All 60 within $0.0512$ of $A=8$, all $\ge 0.95$
  from $A=7$. 0 STALLED, 0 split candidates, 0 violating, 0 halts.
  *(script: session_q1_wkb.py, fit window N in [10,60] step 2, N_ref=150)*

- **Newton polygon at d=4: $\xi_0(b)=4/a_4^{1/4}$ verified spread 0.**
  Derivation: edge $(0,0)\to(1,4)$ with char poly
  $\chi(c)=1-(a_4/256)c^4$. Empirical check on 8 representatives
  (4 catalogue, 3 anchors, 1 with $a_4=7$): $c(4):=\xi_0\cdot a_4^{1/4}=4$
  to spread 0 at dps=80.
  *(script: session_q1_newton.py)*

- **General-d formula derived (claim Q1-B SUPPORTED with c(d)=d):**
  $\xi_0(b)=d/\beta_d^{1/d}$ for any degree-$d$ PCF in standard class.
  Recovers $d=2$ Prop.~3.3.A and predicts $d=3$ value
  $\xi_0=3/\beta_3^{1/3}$.

- **5 AEAL claims** in `claims.jsonl` (2 Phase Q1-1, 1 Phase Q1-2,
  2 Phase Q1-3).

## Judgment calls made

1. **Catalogue size 60, not 30–80 evenly.** The lex window in the
   prompt is structurally heavy on $a_4=1,a_3=-3$ entries; we stopped
   at 60 to avoid wasting cycles on near-duplicate $S_4$-mixed entries.
   The prompt's lower bound (≥30) is exceeded.

2. **WKB fit window $[10,60]$ instead of $[10,100]$ (cubic protocol).**
   At $d=4$ the $A\,N\log N$ decay is faster, $|\delta_N|$ at $N=80$
   is below the $\mathrm{dps}=1200$ floor for many families, and the
   fit becomes ill-conditioned. $[10,60]$ at step 2 produced 26 usable
   points and stderr ≤ 0.0028 uniformly. Documented in
   wkb_quartic_harvest.tex and rubber_duck_critique.md §2.

3. **Reference $N_{\mathrm{ref}}=150$.** $L_b$ stable at 999+ digits
   (identical to working precision) at $N=150$ in 60/60 families,
   so $|L_N - L_{150}|$ is genuinely $|L_N - L|$ to working precision.
   No need for higher $N_{\mathrm{ref}}$.

4. **No PSLQ scan.** Prompt's two empirical claims are about WKB
   exponents and Newton-polygon roots, not algebraic identities. The
   catalogue contains no CM bin (the natural PSLQ target). Recommended
   as next session. See rubber_duck_critique.md §4.

5. **Newton polygon was derived analytically, not just per-family
   numerically.** The general-$d$ formula $\xi_0=d/\beta_d^{1/d}$ is
   the WKB characteristic root of the canonical operator $L_d$ and is
   bin-, signature-, and Galois-independent. The 8-representative
   numerical verification is therefore confirmation of an already-proven
   structural identity (operational sense), not an empirical
   regularity. See newton_polygon_d4.tex.

6. **Publication-grade trigger fired.** All four conditions met: ≥30
   families ($60$), ≥80% b4_consistent_at_8 (100%), Q1-B SUPPORTED,
   no halts. `session_Q1_results.tex` (3 pp) drafted and compiled.

## Anomalies and open questions

**(THIS IS THE MOST IMPORTANT SECTION.)**

1. **Systematic offset $\overline{A}\approx 2d - 0.05$ across $d=2,3,4$.**
   Mean $A_{\mathrm{fit}}$ at $d=4$ is $7.954$, $0.046$ below $A=8$;
   at $d=3$ was $5.973$ (Session B), $0.027$ below 6; at $d=2$ on
   the $A=4$ stratum was $\approx 3.97$. The $0.05$-per-degree
   pattern is consistent with a fit-window bias from the absent
   $-(\log\log N)$ correction term in the 4-parameter ansatz. None
   of the 60 quartic families is thereby pushed within $0.30$ of
   $A=7$, but a sharp ($\sigma=10^{-4}$) test of "$A=2d$" exactly
   would either need to extend the ansatz or move to a finer grid
   at larger $N$.

2. **No CM quartic in catalogue.** The lex window
   $a_4\in\{1,2,3,5,7\}$, $|a_{i<4}|\le 3$ produced zero totally
   imaginary quartic ($r_1=0$, signature $(0,2)$). The natural
   degree-4 analogue of the cubic $-\_S_3\_\mathrm{CM}$ bin (which is
   the only populated cubic bin and the one Session B harvested
   exclusively) is therefore unrepresented in this catalogue. The
   question whether a CM quartic could populate an $A=7$ branch is
   open.

3. **$\xi_0$ is operational, not Borel-radius.** The derived formula
   $\xi_0=d/\beta_d^{1/d}$ is the WKB characteristic root of the
   irregular singularity at $z=0$. It matches Prop.~3.3.A of
   CHANNEL-THEORY-V11 in the operational sense and equals the
   trans-series leading exponent. It is *not* established to be the
   Borel-plane convergence radius in the Pade-Borel sense; the
   Borel channel was flagged anomalous at $K=12$ for QL15/QL26 by
   BOREL-CHANNEL-K-SCAN. Establishing the identity in the Borel-radius
   sense at $d=4$ would require running the Newton-Birkhoff /
   CC-PIPELINE-G machinery on a representative quartic, which is the
   recommended P+1 follow-up.

4. **Claim Q1-A interpretation: "d=2 anomaly" hypothesis strengthened.**
   The result rules out the alternative interpretation
   "cubic is the anomaly". Both $d=3$ (Sessions B+C1, 50/50) and
   $d=4$ (Session Q1, 60/60) populate only the sharp branch $A=2d$.
   The split at $d=2$ remains uniquely populated by sgn($\Delta_2$).
   This refines op:d2-anomaly toward the working interpretation
   "$d=2$ is genuinely special."

5. **Quartic Galois discriminator (C_4 vs D_4) untested.** The
   `m=1`-branch of the resolvent-cubic algorithm uses
   `sympy.nsimplify` + numerator/denominator perfect-square test. No
   genuine $C_4$ calibration was run (anchors are $D_4, V_4, S_4$).
   The catalogue's 3 $D_4$ entries may include a misclassified
   $C_4$. This does not affect the B4 verdict (bin-insensitive at
   $\sigma=0.0037$), but should be pinned down before publication.
   Recommend $x^4-4x^2+2$ as a $C_4$ calibration anchor.

## What would have been asked (if bidirectional)

- Whether to extend the lex window to $a_4=1, |a_{i<4}|\le 5$ to pick
  up CM quartics (would have added ~50 more families and ~3 minutes
  compute, but probably worth the effort given §2 above).
- Whether to attempt a Borel-radius verification at $d=4$ on one
  representative (analogous to vquad_recovery.py at $d=2$). This is
  the natural P+1 session.
- Whether to drop session_Q1_results.tex into pcf2_program_statement.tex
  as a v1.2 patch in this session, or to defer to a separate
  PCF2-V12-RELEASE session.

## Recommended next step

**PCF2-SESSION-D1: CM-QUARTIC HARVEST.** Extend the lex window to
$a_4=1$, $|a_{0,1,2,3}|\le 5$ targeting totally imaginary quartics
(signature $(0,2)$); harvest WKB exponents and run PSLQ on a
$L(\chi_K, s)$-augmented basis for each CM family. Halt-on-positive
will fire if any CM quartic populates an $A=7$ branch; if none does,
B4 sharp form gains a CM-class bin and the conjecture is closer to
publishable.

Alternative P+1: **BOREL-CHANNEL-D4-PROBE.** Run newton_birkhoff.py
adapted for $d=4$ on family 35 ($+_$real$_D_4$, the unique totally
real quartic in the catalogue) to test whether the
operational-$\xi_0$ identity extends to a Borel-radius identity.

## Files committed

- quartic_family_enumeration.py
- quartic_family_catalogue.json
- session_q1_wkb.py
- results.json
- session_q1_newton.py
- newton_d4_results.json
- results_table.tex
- wkb_quartic_harvest.tex
- newton_polygon_d4.tex
- cross_degree_table.tex
- session_Q1_results.tex (publication-grade companion)
- session_Q1_results.pdf
- claims.jsonl  (5 entries)
- claims_phaseQ1_1.jsonl, claims_phaseQ1_2.jsonl, claims_phaseQ1_3.jsonl
- calibration_anchors.json
- galois_distribution_summary.json
- rubber_duck_critique.md
- handoff.md
- halt_log.json (no halts)
- discrepancy_log.json (empty)
- unexpected_finds.json (empty)
- enum.log, wkb_run.log, newton_run.log
- session_Q1_results.aux/log/out (LaTeX artefacts)

## AEAL claim count

5 entries written to `claims.jsonl` this session:
- 2 from Phase Q1-1 (catalogue enumeration + calibration anchors)
- 1 from Phase Q1-2 (WKB-quartic harvest aggregate)
- 2 from Phase Q1-3 (per-rep $\xi_0$ verification + general-$d$ formula)
