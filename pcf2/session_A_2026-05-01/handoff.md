# Handoff — PCF2-SESSION-A
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~30 minutes
**Status:** COMPLETE

## What was accomplished
First cubic-family catalogue for PCF-2 built. Enumerated the first 50
$\mathbb{Z}$-primitive irreducible cubics $b(n)=\alpha_3 n^3+\alpha_2 n^2+\alpha_1 n+\alpha_0$
with $\alpha_3\in\{1,2,3,5,7\}$, $\alpha_2,\alpha_1,\alpha_0\in\{-3,\dots,3\}$
in lex order; computed exact $\Delta_3$, prime factorisation,
Galois group, splitting/CM field, Conjecture-B3(i) trichotomy bin, and
the partial CF limits $L_{500}, L_{1000}, L_{2000}$ at dps=200 with
$a(n)=1$. Recorded a 10-entry reducible-control list and verified all
three calibration anchors (Apéry, Catalan seed $n^3-2$, cyclic-$C_3$
probe $n^3-3n+1$).

## Key numerical findings
- Bin counts (50-family Session-A catalogue): `+_S3_real`=10, `+_C3_real`=3, `-_S3_CM`=37, `other`=0, `out`=0. (script: `cubic_family_enumeration.py`, dps=200)
- Reducible-control: 10 tuples in the same lex window factor over $\mathbb{Q}$ (e.g. $(1,-3,-3,1)\to(n+1)(n^2-4n+1)$, $(1,-3,-1,3)\to(n-3)(n-1)(n+1)$).
- Exotic small-discriminant cyclic cubics in the window: 3 candidates with $|\Delta_3|\le 100$, namely $(1,-3,0,1)$ and $(1,-3,0,3)$ both with $\Delta_3=81$ (splitting field $\mathbb{Q}(\zeta_9+\zeta_9^{-1})$), and $(1,-2,-1,1)$ with $\Delta_3=49$ (splitting field $\mathbb{Q}(\zeta_7+\zeta_7^{-1})$).
- CM-field histogram for the 37-element `-_S3_CM` bin spans 26 distinct imaginary quadratic fields; most populous are $\mathbb{Q}(\sqrt{-3})$ (5), $\mathbb{Q}(\sqrt{-23})$ (3), and $\mathbb{Q}(\sqrt{-7}), \mathbb{Q}(\sqrt{-19}), \mathbb{Q}(\sqrt{-26}), \mathbb{Q}(\sqrt{-59}), \mathbb{Q}(\sqrt{-239})$ each with multiplicity 2.
- All 50 families converge "fast" (per d23/d12 ratio test): the constant-$a$ continued fraction with cubic $b(n)$ is super-linearly convergent uniformly across the catalogue.

## Calibration anchors (mandatory cross-cuts) — all PASS
| Anchor | Expected | Computed | Match |
|---|---|---|---|
| Apéry $b=34n^3+51n^2+27n+5$ | reducible $(2n+1)(17n^2+17n+5)$, $\Delta_3=-459=-3^3\cdot 17$ | reducible $(2n+1)(17n^2+17n+5)$, $\Delta_3=-459$ | ✓ |
| Catalan seed $b=n^3-2$ | $\Delta_3=-108$, Gal=$S_3$, CM=$\mathbb{Q}(\sqrt{-3})$, bin `-_S3_CM` | $\Delta_3=-108$, Gal=$S_3$, CM=$\mathbb{Q}(\sqrt{-3})$, bin `-_S3_CM` | ✓ |
| Cyclic-$C_3$ $b=n^3-3n+1$ | $\Delta_3=81=9^2$, Gal=$C_3$, splitting field $\mathbb{Q}(\zeta_9+\zeta_9^{-1})$, bin `+_C3_real` | $\Delta_3=81$, square=True, Gal=$C_3$, bin `+_C3_real` | ✓ |

`anchors.all_match = True`. Halt rule (anchor disagreement) NOT triggered.

## Sample 5-row catalogue table
| id | $(\alpha_3,\alpha_2,\alpha_1,\alpha_0)$ | $\Delta_3$ | Gal | bin | CM field | $L_{2000}$ (40-digit) | rate |
|---|---|---:|---|---|---|---|---|
| 1 | $(1,-3,-3,-3)$ | $-864$ | $S_3$ | `-_S3_CM` | $\mathbb{Q}(\sqrt{-6})$ | $-3.123817690476524074720780117\ldots$ | fast |
| 2 | $(1,-3,-3,-2)$ | $-459$ | $S_3$ | `-_S3_CM` | $\mathbb{Q}(\sqrt{-51})$ | $-2.141189540536045156472943106\ldots$ | fast |
| 3 | $(1,-3,-3,-1)$ | $-108$ | $S_3$ | `-_S3_CM` | $\mathbb{Q}(\sqrt{-3})$ | $-1.164201930217769324417923652\ldots$ | fast |
| 4 | $(1,-3,-3, 2)$ | $+621$ | $S_3$ | `+_S3_real` | $\mathbb{Q}(\sqrt{69})$ | $+1.679769919258163267805355469\ldots$ | fast |
| 5 | $(1,-3,-3, 3)$ | $+756$ | $S_3$ | `+_S3_real` | $\mathbb{Q}(\sqrt{21})$ | $+2.532591481184963433840447529\ldots$ | fast |

(Family 19, $b=n^3-3n+1$, is the in-window cyclic-$C_3$ anchor: $\Delta_3=81$, $L_{2000}\approx 0.32727920\ldots$.)

## Bin distribution (for Sec 4 of an eventual PCF-2-results paper)
- `+_S3_real` (predicted: elementary closed form): 10
- `+_C3_real` (predicted: elementary closed form, possibly cyclic-cubic units): 3
- `-_S3_CM`  (predicted: Painlevé/Stokes obstruction, cubic analogue of PCF-1's $\Delta<0$ cell): 37
- `other`: 0
- `out` (reducible, in control list): 10 (separate)

The window is dominated by the $-_S3_CM$ cell (74%), which is the
analytically interesting one. Three cyclic-$C_3$ specimens are already
present in the 50-row sample; the cyclic cell will not need a window
widening.

## Judgment calls made
- **Window choice for the reducible control list**: kept it to the
  first 10 reducible Z-primitives interleaved with the irreducibles,
  rather than scanning the whole window, because the prompt caps the
  control list at 10. Apéry's polynomial is OUT of the enumeration
  window (coefficients exceed 3) and is reported separately under
  `calibration_anchors.json`, not in the reducible-control list.
- **CM-field convention** for $S_3$ cubics: I record
  $\mathbb{Q}(\sqrt{\Delta_3})$ reduced to its squarefree representative
  (e.g. $\Delta_3=-108\Rightarrow\mathbb{Q}(\sqrt{-3})$) since
  $\mathbb{Q}(\sqrt{D})$ depends only on the squarefree part of $D$.
- **Continued-fraction convention**: I evaluated $L_N=b(0)+a(1)/(b(1)+a(2)/(b(2)+\ldots))$ with $a(n)=1$, the natural "constant-$a$" PCF for Session A. Subsequent Sessions A2/A3 should sweep $a(n)$ over polynomial seeds.
- **"fast" / "exponential" / "linear" labelling**: heuristic ratio test
  $|L_{2000}-L_{1000}|/|L_{1000}-L_{500}|$, with thresholds $10^{-6}$ /
  $0.5$. All 50 families came out "fast"; this is consistent with the
  dominant cubic growth in $b(n)$ and is not a deep statement.

## Anomalies and open questions
- **The $-_S3_CM$ cell over-dominates** the small-coefficient window
  (37/50). This is morally correct (the discriminant of a generic cubic
  is negative), but it means Session-B / convergent-acceleration
  experiments on the $+_C3_real$ cell will need targeted enumeration
  beyond the small window — the three cyclic-$C_3$ specimens here are
  too few for a Painlevé-vs-elementary discrimination study.
- **Cyclic-$C_3$ anchor not unique**: the calibration probe
  $b=n^3-3n+1$ has $\Delta_3=81$. Family 20, $b=n^3-3n+3$, ALSO has
  $\Delta_3=81$. Both are in the same $\mathbb{Q}(\zeta_9+\zeta_9^{-1})$
  field (same discriminant of the splitting field). This is expected
  (the discriminant is invariant under translations of the variable
  giving the same Galois closure) but worth flagging — Session B should
  test whether the corresponding CF limits are arithmetically related.
- **Family 46, $b=n^3-2n^2-n+1$**, has $\Delta_3=49=7^2$, which is the
  smallest-disc real cyclic cubic in the window after $\Delta_3=81$.
  Its splitting field is $\mathbb{Q}(\zeta_7+\zeta_7^{-1})$ (the
  totally real subfield of $\mathbb{Q}(\zeta_7)$, conductor 7) — a
  well-studied cyclic cubic. **Recommended Session A2 anchor.**
- **No $\Delta_3=0$ degenerate cases** appeared (would have indicated
  a repeated root); good.
- **CF convergence label is "fast" everywhere**, which prevents the
  rate column from discriminating cells. A more refined rate
  diagnostic (e.g. $-\log|L_N-L_\infty|/\log N$ with a high-$N$
  reference) is needed for Session B if convergence-rate is to enter
  the trichotomy diagnostic.

## What would have been asked (if bidirectional)
- Should the reducible-control list be enlarged (10 is sparse for a
  statistical baseline) or kept tight per spec?
- For Session A2: should we sweep $a(n)\in\{1, n, n+1, n^2, n^2+n+1\}$
  uniformly across all 50 $b$-families, or focus on the 3+1 cyclic-$C_3$
  specimens (where the closed-form prediction is most testable)?
- Should the CM-field histogram be augmented with class-number data
  (a candidate cross-cut for Conjecture B3(ii))?

## Recommended next step
**PCF-2 Session B — cyclic-$C_3$ closed-form attack.** Take the four
cyclic-$C_3$ specimens identified here ($b=n^3-3n+1$, $b=n^3-3n+3$,
$b=n^3-2n^2-n+1$, plus a fourth $\Delta_3=49$ companion if one exists
in a slightly widened window) and attempt PSLQ identification of
$L_\infty$ in the integral closure of $\mathbb{Q}(\zeta_7+\zeta_7^{-1})$
and $\mathbb{Q}(\zeta_9+\zeta_9^{-1})$ respectively. This is the
analogue of the PCF-1 $\Delta>0$ closed-form sweep and should produce
the first cubic-cyclotomic continued-fraction identities.

## Files committed
- `cubic_family_enumeration.py`
- `cubic_family_catalogue.json`     (50 families + 10 reducible control)
- `galois_distribution_summary.json`
- `calibration_anchors.json`
- `claims.jsonl`                    (5 AEAL entries)
- `run.log`
- `halt_log.json`                   (empty `{}`)
- `discrepancy_log.json`            (empty `{}`)
- `unexpected_finds.json`           (one finding: see above re Family 46 and the duplicate $\Delta_3=81$ cyclic cubic)
- `handoff.md`

## AEAL claim count
5 entries written to `claims.jsonl` this session.
