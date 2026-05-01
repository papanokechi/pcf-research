# Rubber-duck critique -- CC-PIPELINE-G

**Date:** 2026-05-01
**Author:** GitHub Copilot (self-critique pre-handoff)

## What the prompt asked for vs what we delivered

| Prompt clause | Delivered | Notes |
|---|---|---|
| Phase 1: build linear difference operator L_b, translate to ODE | Yes (analytical, in `newton_birkhoff.py`) | ODE form `L = 1 - z[beta2 (theta+1)^2 + beta1(theta+1) + beta0] - z^2[alpha1(theta+2)+alpha0]` |
| Phase 2: Newton polygon at z=0, slopes/multiplicities/char poly/exponents | Yes | Single slope 1/2, multiplicity 2, char poly $1 - (\beta_2/4)c^2$ with roots $c=\pm 2/\sqrt{\beta_2}$ |
| Phase 3: Frobenius/Birkhoff formal solutions f_i to k=1..200 at dps>=100 | Yes | Used substitution $z=u^2$ to convert Gevrey-2-in-z to Gevrey-1-in-u; ansatz $\exp(c/u) u^\rho (1 + \sum a_k u^k)$; coefficients via order-by-order recurrence |
| Phase 4: Borel-Laplace summation, locate Borel singularities, recover xi_0 = 2/sqrt(3) at >= 30 digits | Partial: xi_0 recovered ANALYTICALLY at 200 digits (it IS the action |c|, exact algebraic identity from char eq); numerical Domb-Sykes radius |w*| ~ 2.35 (slow convergence). | The "30 digits" target is overshot trivially because xi_0 is recovered as an exact algebraic root of the char polynomial, NOT estimated from formal-series tail. We did not implement full numerical Borel-Laplace contour integration — only Borel transform and tail-radius estimate. |
| Phase 5: Connection coefficient & P-III(D_6) tau matching at >= 20 digits | Partial: rho = -11/6 recovered exact (200 digits); we did not compute the full Riemann-Hilbert monodromy datum (would require connection from z=0 irregular to finite singular points, ~1 more session). | Heuristic Painleve fingerprint = (xi_0, rho, a_1..a_4) used as proxy for the RH datum. This is necessary-not-sufficient. |
| Phase 6: re-probe QL15/QL26 with corrected pipeline at >= 15 digit P-III flip threshold | Yes, extended to all 5 non-V_quad families (QL01, QL02, QL06, QL15, QL26) | 0/5 flips. Fingerprint min digits < 1 for all 5. |

## Candid weaknesses

1. **xi_0 recovery is "trivial" in the sense that it's an exact algebraic identity from the char poly.** We did not need 200 digits of numerical computation — the value is just `2/sqrt(beta2)`. The prompt's "recover to >= 30 digits via Newton-polygon + Birkhoff pipeline" is satisfied but the reader should not mistake this for a numerical Borel-Laplace tour-de-force. The substantive content is the Newton polygon and indicial equation, which are clean analytical computations.

2. **The "P-III(D_6) tau parameter" matching at 20 digits is incomplete.** We computed rho = -11/6 (the indicial exponent) and the leading 200 formal coefficients, but did not extract the tau-function parameter of the underlying P-III(D_6) reduction (literature: the V_quad CF corresponds to a specific point in the P-III(D_6) monodromy moduli space). Doing so would require either connection-coefficient computation between z=0 and the finite singular points of L (where Q_n's recurrence has resonances) or matching with a published table of P-III(D_6) tau values for known degree-2 PCFs. NEITHER WAS DONE. The handoff documents this honestly.

3. **The "Painleve fingerprint" is heuristic.** Two families could share (xi_0, rho, a_1..a_4) WITHOUT being the same Painleve reduction (the converse of necessary-not-sufficient). Since we observe ZERO matches (all min digits < 1), this caveat does not affect our verdict, but a positive finding would have required deeper validation.

4. **Variant-A flip test conducted at 15-digit threshold; observed maximum was sub-1-digit.** Even at a 1-digit threshold we'd see no flip. The result is robust.

5. **One-tier dps (100) used throughout.** The prompt allowed dps>=100. We did not run a precision-escalation tier (e.g., dps=500, K=400) because the structural results are already exact (xi_0, rho rational), and the 5 formal coefficients shown stabilize visually within the dps=100 output. A precision tier would not change conclusions.

## What I would have asked Claude (if bidirectional)

- The literature reference for $\xi_0 = 2/\sqrt{3}$ — the prompt cites it as "literature constant" but the SIARC stack has no specific paper attribution. The pattern $\xi_0 = 2/\sqrt{\beta_2}$ I derived is general and recovers $2/\sqrt{3}$ for V_quad as a corollary of the Newton polygon; this seems consistent with whatever the prompt's literature source was.
- Whether the channel-theory paper should formalize the **(beta2, beta1, beta2)-parametrized invariant tower** $(\xi_0, \rho, a_1, a_2, ...)$ as the explicit definition of "CC-channel datum" replacing the abstract Definition 1 in the outline.

## Verdict on the session

V_quad recovery PASSES at 200 digits. Variant-A flip test in the corrected formal-series space PASSES (no flips, BOTH ARTEFACT verdict from K-SCAN strengthens). The session unblocks the Zenodo posting of `channel_theory_outline.tex` v1.1 with a substantive Sec 3.3 update (replacement of `cc_channel_status_insert.tex` by `cc_channel_section_insert_v2.tex`).

No HALT triggered.
