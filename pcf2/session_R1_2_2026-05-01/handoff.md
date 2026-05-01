# Handoff — PCF2-SESSION-R1.2

**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~50 minutes
**Status:** HALTED (B5 cross-degree FALSIFIED at d=4; halt fired per protocol)

## What was accomplished

Tested whether the $j$-invariant finer-split discovered at $d=3$ in
PCF2-SESSION-R1.1 (Conjectures B5 + B6) extends to $d=4$ on Session-Q1's
60-family quartic catalogue. Computed Igusa $I,J$ invariants and
$j(\mathrm{Jac}(C_b))$ for all 60 quartics. Diagnosed and corrected a
small-$N$ WKB sub-leading bias via a tail-window fit
(dps=2000, $N\in[50,250]$, $N_{\rm ref}=400$); built a unified 110-row
dataframe (50 cubic from R1.1 + 60 quartic). Verdict: **B5 cross-degree
FALSIFIED**, **B6 cross-degree NOT SUPPORTED**, **cross-degree
universality NOT SUPPORTED** (joint signal is a degree-confounded
artefact). HALT triggered per protocol.

## Key numerical findings

- **R12-A (B5 at d=4):** FALSIFIED. The $j=0$ cell at $d=4$ has size
  $n_4 = 1$ (family 32, $b(n)=n^4-3n^3-3n^2+3n-3$, $I=0$). Tail-window
  fit at dps=2000 gives $\delta_4(32) = -3.71\times 10^{-3}$ with
  $A_{\rm stderr}=7.71\times 10^{-5}$, i.e.\ $\sim 48\sigma$ from $0$.
  The entire 60-family quartic population clusters in
  $\delta_4 \in [-3.83\times 10^{-3}, -3.19\times 10^{-3}]$ (script
  `quartic_tail_fit_all60.py`).
- **R12-B (B6 at d=4):** NOT SUPPORTED. Spearman
  $\rho_S(\log|j|, \delta_4) = +0.073$ on $N=60$ (raw $p=0.58$,
  Bonferroni $p=1.0$); sign opposite to $d=3$ signal but not significant
  (script `r1_2_quartic_j_probe.py`, dps=2000 tail fits).
- **R12-C (cross-degree universality):** NOT SUPPORTED. Joint
  $\rho_S(\log|j|,\delta) = -0.544$, Bonferroni $p = 9.5\times 10^{-9}$
  on $N=110$ — but this is dominated by the $d=3$ signal (degree is a
  confounder; cubic has wider $\log|j|$ spread).
- **d=3 reproduction:** $\rho_S = -0.6906$, Bonferroni $p = 3.4\times 10^{-7}$
  on $N=50$ — matches R1.1 to all reported digits (script
  `r1_2_quartic_j_probe.py`, R1.1 csv).
- **$j=1728$ cells:** empty at both $d=3$ (R1.1 catalogue) and $d=4$
  (Q1 catalogue); `op:lemniscatic-cell` remains untestable.
- **WKB sub-leading bias diagnosis:** at dps=80, $N\in[10,60]$, all 60
  quartics show $\delta_4 \approx -0.046$ (universal small-$N$ bias).
  Escalating to dps=1500, $N\in[10,130]$, $N_{\rm ref}=200$ reduces to
  $\sim -0.022$. Tail fit at dps=2000, $N\in[50,250]$, $N_{\rm ref}=400$
  reduces further to $\sim -0.0035$. Family 32 follows the same
  monotone reduction trajectory as non-$j=0$ neighbours (script
  `fam32_deep_escalation.py` and `quartic_tail_fit_all60.py`).

## Judgment calls made

- Diagnosed that the dps=80 small-$N$ fit window biases $\delta_4$ by
  $\sim -0.05$ uniformly. Built a tail-window fit script
  (`quartic_tail_fit_all60.py`, dps=2000, $N\in[50,250]$, $N_{\rm ref}=400$)
  to suppress this bias and used these overrides as the **primary**
  $A_{\rm fit}$ values for the 60 quartics. The original Q1 dps=80 fits
  are NOT used in the v3 dataframe. (Justification: at dps=80 the $j$=0
  signal would be drowned by the universal bias; using the unbiased
  fits is necessary to test B5/B6 fairly.)
- Set the consistency-with-zero threshold to
  $|\delta| \le \max(2 \cdot A_{\rm stderr},\, 10^{-3})$. The $10^{-3}$
  floor mirrors the R1.1 cubic threshold.
- Hand-curated `v13_paragraph_insert_v2.tex` to reflect the
  FALSIFIED/thin-cell-caveat verdict, since the auto-template was built
  assuming SUPPORT and would have over-claimed.
- Did NOT attempt deep-WKB harvest at $d=4$ ($N_{\rm ref}\ge 1000$,
  dps$\ge 5000$) because the relay budget (~2-3 h) doesn't accommodate
  the $\sim 30$ h cost. Recommended as the primary follow-up.
- pari/gp remains unavailable (matches R1.1 status); skipped pari
  computations rather than installing in this environment.

## Anomalies and open questions

**THIS IS THE MOST IMPORTANT SECTION.**

1. **The d=4 quartic $\delta_4$ population is monotone in
   $(\alpha_2, \alpha_1, \alpha_0)$ position, not in $\log|j|$.**
   Inspection of `tail_fit_overrides.json` shows $\delta_4$ varying
   smoothly from $-3.83\times 10^{-3}$ (fam 1) to $-3.19\times 10^{-3}$
   (fam 60) following the lex-ordered enumeration over
   $(\alpha_2 \in \{-3,-2\},\, \alpha_1, \alpha_0 \in \{-3,\ldots,3\})$.
   The variation magnitude ($\sim 6\times 10^{-4}$) is larger than the
   per-family $A_{\rm stderr}$ ($\sim 7\times 10^{-5}$) by $\sim 9\times$,
   so this is REAL signal — but it is a non-modular signal driven by
   the position in the enumeration window, NOT by $j$. **Open: what
   non-modular invariant of $b$ controls $\delta_4$ at this precision?**
   Candidates: leading-order Mahler measure variation, conductor of
   $E_b$, Igusa $J$ alone (independent of $I$), weighted height of
   $(\alpha_2,\alpha_1,\alpha_0)$.

2. **B5 thin-cell caveat (Remark in v13 insert).** The $d=4$ test
   rests on a SINGLE $j=0$ family (family 32). The "FALSIFIED" verdict
   is structurally weaker than would be ideal. Two interpretations
   remain consistent with the data:
   - (i) B5 is genuinely $d=3$ specific (current best guess).
   - (ii) The WKB fit window $N\le 250$ is too small to resolve a
         $j$-driven separation that would emerge at much larger $N$.
   Recommended: extend Q1 catalogue OR run deep-WKB harvest at $d=4$.

3. **The cubic (R1.1) signal is robust:** d=3 reproduction of
   $\rho_S = -0.691$ across the joined dataframe is exact to 4 digits,
   confirming R1.1's discovery is not an artefact.

4. **Joint $d\in\{3,4\}$ Spearman is statistically significant
   ($p\le 10^{-8}$) but should NOT be reported as cross-degree
   confirmation.** It is the cubic signal alone, leaking into the
   joint statistic because cubic has both wider $\log|j|$ spread
   AND wider $\delta$ spread. A degree-aware partial correlation
   (controlling for $d$) would zero out the joint signal — recommended
   as a follow-up sanity check.

5. **The $j=1728$ (lemniscatic) cell is empty at both $d=3$ and $d=4$
   in current catalogues.** $J=0$ requires $a_3 \in\{-3\}$,
   $a_2\in\{-3,-2\}$, $\alpha_1, \alpha_0$ such that
   $J = -459\alpha_0 - 27\alpha_1^2 + 81\alpha_1 + 54 = 0$ (for
   $\alpha_2=-3$) or analogous; no integer solutions in the current
   window. **op:lemniscatic-cell remains untestable until catalogues
   extend.**

## What would have been asked (if bidirectional)

- "The $d=4$ quartic delta ordering by $(\alpha_2,\alpha_1,\alpha_0)$
  position: should we run a Pearson regression on $\delta_4$ vs.
  $(\alpha_2, \alpha_1, \alpha_0)$ to identify the dominant non-modular
  invariant before claiming B5 is degree-3 specific?"
- "Should we attempt the deep-WKB harvest ($N_{\rm ref}\ge 1000$,
  dps$\ge 5000$) on family 32 alone (~30 min cost) to definitively
  distinguish (i) vs (ii) of the thin-cell caveat?"
- "Is it acceptable to treat the cubic (R1.1) signal as confirmation of
  Conjecture B5 at $d=3$ ALONE in PCF-2 v1.3, while explicitly
  withdrawing the cross-degree extension claim?"

## Recommended next step

**Two parallel paths, in priority order:**

1. **PCF2-V13-RELEASE-DRAFT** — fold `v13_paragraph_insert_v2.tex` into
   PCF-2 v1.3 source. Restate Conjecture~B5 as $d=3$-only (with $n_3=4$
   evidence). Add new Conjecture (B5/B6 are degree-3 specific in the
   current WKB regime). Add `op:b5-cross-degree-deep` to Open Problems.
   Withdraw v1.2's "cross-degree universality" framing if any leaked in.

2. **PCF2-SESSION-R1.3 (cross-degree-deep)** — extend the quartic
   catalogue to $a_3 \in\{-3,-2,-1,1,2,3\}$, $a_2\in\{-3,\ldots,3\}$
   (estimated $\sim 500$ irreducible quartics) AND run deep-WKB on
   family 32 at $N_{\rm ref}\ge 1000$. Goal: harvest $\ge 2$ further
   $j=0$ quartics so the B5-at-$d=4$ test is no longer thin, and confirm
   or refute interpretation (ii) of the thin-cell caveat.

(The standalone-note publication trigger does NOT fire: B5 at $d=4$ is
NOT supported, B6 at $d=4$ is NOT significant, cross-degree universality
is NOT supported.)

## Files committed

- `r1_2_quartic_j_probe.py` (main probe, 33 KB)
- `quartic_jacobian_invariants.json` (60 rows: I, J, j, sgn, is_j_zero, is_j_1728, CM_disc)
- `assembled_data_v3.csv` (110 rows: 50 cubic + 60 quartic, unified schema)
- `results_v3.json` (per-degree + joint correlations, verdicts)
- `correlation_table_v3.tex`
- `quartic_j_zero_cell.tex`
- `v13_paragraph_insert_v2.tex` (hand-curated, FALSIFIED-with-caveat)
- `claims.jsonl` (9 entries)
- `rubber_duck_critique.md` (4-question critique sign-off)
- `tail_fit_overrides.json` (60 quartic dps=2000 tail fits — primary
  data)
- `quartic_tail_fit_all60.py` (script that produced the overrides)
- `precision_escalation_log.json` (legacy dps=1500 escalation, kept for
  diagnostic)
- `fam32_deep_escalation.py` + `fam32_deep_escalation.json` +
  `fam32_deep_escalation.log` (diagnostic establishing the
  WKB-fit-window bias)
- `halt_log.json` (B5_d4_falsified entry)
- `discrepancy_log.json` (empty)
- `unexpected_finds.json` (thin_jzero_cell_d4 + lemniscatic_cell_empty)
- `run.log`, `stdout.log`, `tail_fit_run.log`
- `handoff.md` (this file)

## AEAL claim count

9 entries written to `claims.jsonl` this session.
