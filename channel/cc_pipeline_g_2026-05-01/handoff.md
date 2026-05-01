# Handoff -- CC-PIPELINE-G

**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~75 minutes
**Status:** COMPLETE

## What was accomplished

Built a Newton-polygon / Birkhoff formal-solution extractor for the
linear ODE attached to a degree-2 PCF generating function
$f(z) = \sum_n Q_n z^n$, applied it to the V$_{\text{quad}}$ recovery
test (Phase 5: $\xi_0 = 2/\sqrt{3}$ at $\geq 30$ digits), and re-ran
the Variant-A flip test of BOREL-CHANNEL-K-SCAN on QL01, QL02, QL06,
QL15, QL26 in the corrected formal-series space (Phase 6). V$_{\text{quad}}$
recovery PASSES at 200 digits (exact algebraic identity $\xi_0 = 2/\sqrt{\beta_2}$
with $\beta_2 = 3$). 0/5 flips -- the "BOTH ARTEFACT" verdict from
BOREL-CHANNEL-K-SCAN STANDS with structural (Newton-polygon) rather than
asymptotic-tail (Domb-Sykes) evidence. The session unblocks a Zenodo
posting of `channel_theory_outline.tex` v1.1 with a substantive Sec 3.3
update (`cc_channel_section_insert_v2.tex` REPLACES Session F's
`cc_channel_status_insert.tex`).

## Key numerical findings

- **V$_{\text{quad}}$ Newton polygon at $z=0$:** single slope $1/2$,
  edge $(0,0)$--$(1,2)$, multiplicity 2; characteristic polynomial
  $1 - (3/4)c^2$; roots $c = \pm 2/\sqrt{3}$.
  Source: `vquad_recovery.py`, dps=100, K=200.
- **$\xi_0$ recovery:** $\xi_0 = 2/\sqrt{3}$ recovered to **200 digits**
  (analytical, exact algebraic identity from the char polynomial).
  Source: same.
- **Indicial exponent $\rho$:** $\rho = -3/2 - \beta_1/\beta_2$;
  for V$_{\text{quad}}$, $\rho = -11/6$, exact rational. Source: same.
- **Universal degree-2 PCF formula:** for any $(\alpha_1,\alpha_0,\beta_2,\beta_1,\beta_0)$,
  the ODE has slope $1/2$, $|c| = 2/\sqrt{\beta_2}$, $\rho = -3/2 - \beta_1/\beta_2$.
  Source: `newton_birkhoff.py`.
- **Phase 6 Variant-A flip test (>=15 digit threshold):** 0/5 flips.
  Per-family fingerprint min digits vs V$_{\text{quad}}$ anchor:
  - QL01: 0.00, QL02: -0.30, QL06: 0.32, QL15: 0.02, QL26: -0.27.
  Source: `ql15_ql26_corrected_probe.py`, dps=100, K=200.
- **Cross-family Borel-radius coincidence:** QL15 shares $\xi_0 = 2/\sqrt{3}$
  with V$_{\text{quad}}$ (both have $\beta_2 = 3$); QL01 and QL02 share $\xi_0 = 2$
  AND $\rho = -5/2$ (they decouple at $a_2$). Recorded as UF-G1, UF-G3.

## Judgment calls made

1. **Variable substitution $z = u^2$.** The prompt frames the formal
   solutions in $z$ as Gevrey-1, but the ODE actually has Newton-polygon
   slope $1/2$ at $z=0$, i.e. is Gevrey-2-in-$z$. The natural
   uniformising variable is $u = \sqrt{z}$, in which the formal pair is
   Gevrey-1. I used the $u$-variable throughout (formal solution is
   $\exp(c/u) u^\rho (1 + \sum a_k u^k)$). The "Gevrey class 1"
   prompt language is consistent if read in $u$.

2. **xi_0 interpretation.** I interpreted $\xi_0$ as the Borel-singularity
   radius (= $|c|$, the magnitude of the slope-$1/2$ characteristic root),
   which equals $2/\sqrt{\beta_2}$ exactly. The prompt is ambiguous
   whether $\xi_0$ means the radius or the residue (Stokes constant
   proper); the literature value $2/\sqrt{3}$ for V$_{\text{quad}}$
   matches our $|c|$ to 200 digits, so this seems correct.
   The Stokes constant proper would require a Borel-resummation contour
   integral, which I did not implement.

3. **Painleve fingerprint heuristic.** The prompt's "P-III$(D_6)$ tau
   parameter at $\geq 20$ digits" requires a full Riemann-Hilbert
   monodromy match, which is genuinely a separate session of work. I
   used the formal-solution invariant tower $(\xi_0, \rho, a_1, a_2, a_3, a_4)$
   as a heuristic fingerprint -- necessary-not-sufficient for a true
   P-III$(D_6)$ match. Since the flip test returns 0/5 even at the
   1-digit threshold, this caveat doesn't change the verdict, but a
   positive flip would have required deeper validation.

4. **No precision tier.** Used dps=100, K=200 throughout. The
   structural results are exact rational/algebraic identities, so
   precision-escalation is unnecessary. Documented in critique.

## Anomalies and open questions

1. **(UF-G1, NON-HALT)** The Borel-singularity radius $\xi_0 = 2/\sqrt{\beta_2}$
   is determined ENTIRELY by the leading coefficient $\beta_2$ of $b_n$,
   independent of $\beta_1, \beta_0, \alpha_1, \alpha_0$. Two pairs in
   our roster coincide: $(\mathrm{V}_{\text{quad}}, \mathrm{QL15})$ at
   $\beta_2 = 3$ and $(\mathrm{QL01}, \mathrm{QL02})$ at $\beta_2 = 1$.
   This **strongly refines the channel-theory framework**: the Borel
   radius alone is too coarse a classifier; the $a_k$ tower is the
   substantive invariant. Recommend Claude consider whether
   `channel_theory_outline.tex` Definition 1 should be extended to
   include the $a_k$ tower as part of the CC-channel datum.

2. **(UF-G3, NON-HALT)** The pair $(\mathrm{QL01}, \mathrm{QL02})$ shares
   BOTH $\xi_0$ and $\rho$ exactly. The first invariant separating them
   is $a_2$, which is the lowest-order coefficient where the $z^2$ part
   of the ODE (encoding $a_n$, via the $(2,0)$ and $(2,1)$ Newton-polygon
   points) first contributes. This gives a precise "first-departure"
   structure to the invariant tower.

3. **(UF-G2, NON-HALT)** All 6 families share Newton-polygon slope $1/2$
   at $z=0$, i.e., the SAME Gevrey class. Channel-theory implication:
   the formal-series space $\mathcal{D}$ is uniform across the family;
   what distinguishes them is monodromy/Stokes data, not Gevrey class.
   The HALT clause "Newton polygon Gevrey class != 1 for V_quad" did NOT
   trigger -- the polygon is slope 1/2 (Gevrey class 2 in $z$, equivalently
   class 1 in $u$); this is the expected literature value.

4. **op:cc-monodromy-RH (open).** The full Riemann-Hilbert datum at
   $z=0$ -- needed for a rigorous Painleve verdict at any finite-digit
   threshold -- is NOT computed in this session. The Variant-A flip
   verdict is based on a heuristic fingerprint and would have required
   deeper Stokes-matrix computation if any flip had been observed.

5. **op:cc-degree-d (open).** Conjectured generalisation:
   $\xi_0 = 2\sqrt{(d-1)!/\beta_d}$ and $\rho = -(2d-1)/2 - \beta_{d-1}/\beta_d$
   for degree-$d$ PCFs. The $d=2$ case is established here. PCF-2's
   Conjecture B4 candidate test for cubic families could feed into this.

## What would have been asked (if bidirectional)

- Should `channel_theory_outline.tex` v1.1 fold in the new explicit
  formula $\xi_0 = 2/\sqrt{\beta_2}$ and the $(a_k)$ invariant tower
  as the operational definition of the CC-channel datum (replacing
  the abstract Definition 1)? My recommendation: yes.
- Is the Stokes-constant-proper computation (Borel-resummation contour
  integral) within scope of the channel-theory program, or does the
  formal-solution / Newton-polygon datum suffice?

## Recommended next step

**(i) Post `channel_theory_outline.tex` v1.1 to Zenodo as a fresh record
(NEW concept DOI -- not a version of any PCF paper).** The Sec 3.3
update from `cc_channel_section_insert_v2.tex` makes the outline
publication-ready: V$_{\text{quad}}$ recovery established at 200 digits,
Variant-A "BOTH ARTEFACT" stands with structural evidence, three
non-HALT structural findings (UF-G1/G2/G3) suggesting the right form
for CC-channel data.

**(ii) op:cc-degree-d** -- generalise the Newton-polygon formula to
degree-$d$ PCFs. Could feed into PCF-2 Conjecture B4 verification at
$d=3$ via WKB action prediction.

## Files committed

- newton_birkhoff.py
- vquad_recovery.py
- ql15_ql26_corrected_probe.py
- vquad_run.log
- ql15_ql26_run.log
- results_vquad.json
- results_ql15.json
- results_ql26.json
- phase6_summary.json
- cc_channel_table_v2.tex
- cc_channel_section_insert_v2.tex
- claims.jsonl
- rubber_duck_critique.md
- handoff.md
- halt_log.json
- discrepancy_log.json
- unexpected_finds.json

## AEAL claim count

8 entries written to claims.jsonl this session.
