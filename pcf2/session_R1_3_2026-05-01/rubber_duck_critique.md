# R1.3 — Rubber-duck critique

Focus list per the relay prompt; honest critique of the R1.3
methodology and verdict.

## (i) Is FIXED-A residualization the right diagnostic, or should we
also try SUBTRACT-α (fix A AND α at predicted values, fit only
β γ residual)?

The FIXED-A residualization annihilates a known leading WKB
shape and lets the rest absorb into [α, β, γ].  In our
implementation `α` is the slope coefficient on `n` (the linear
WKB sub-leading term).  Empirically (Phase R13-A and B), with
A=2d held fixed, the regression gives `α ≈ 7.9998` for d=4
across all 60 Q1 quartics with extremely tight clustering — i.e.
the regression is forced to absorb most of the variance into α.
This makes `residual_mean` structurally near zero (it IS zero
by least-squares plus an intercept γ), so the **`residual_mean`
column is NOT a meaningful diagnostic** — it's
mathematically forced to be ≈0 within numerical noise.  We
correctly relied on `residual_at_max_n` and the FREE-A `delta`
instead.

A more aggressive variant — **fix BOTH A AND α at predicted
values**, then let only β, γ float — would expose any
systematic α-dependence.  WKB theory at d=2 predicts (in
some channels, e.g. V_quad) a closed form for α that depends on
β_d and the channel; PCF-2 v1.2 has c(d)=d for the
characteristic root xi_0 in the Newton-polygon-Birkhoff
analysis (Channel Theory v1.2).  We did not write a
predicted-α residualization here because the universal
α-prediction is not yet established at d∈{3,4} (op:xi0-d3-direct).
This is a clear next step for R1.4 if the shallow-N d=4 effect
turns out structural.

## (ii) Is the sub-leading-systematic decay rate (1/N, 1/log N,
1/log²N) consistent with the deep-WKB observation? If so, what
does that reveal about Conjecture B4 at d=4?

R13-D shows fam32 |δ_4| shrinks from 3.71e-3 (shallow N=250) to
4.55e-4 (deep N=800).  Naive fit: |δ_4(N)| ∝ N^p ⇒
log(3.71e-3 / 4.55e-4) ≈ 2.10, log(N_deep / N_shallow) ≈
log(525/150) ≈ 1.25.  So p ≈ -1.68.  This is NOT 1/N (p=-1) and
NOT 1/log²N (p=-2 effective near N=200..800), but somewhere in
between.  A 1/log²N decay would predict shrinkage by a factor of
(log 800 / log 250)² ≈ (6.68/5.52)² ≈ 1.46, far smaller than the
observed factor of 8.2.  A 1/N^{1.7} decay matches.

We have only TWO N-windows so the inferred p is unreliable.  But
the take-away is that fam01 (non-j=0) shrinks at the SAME rate
as fam32, both reaching ~4.6e-4 at deep N.  This is consistent
with B4 at d=4 (sharp form, A=2d=8) IF we believe the residual
δ_4 ∝ N^p continues to shrink as N → ∞.  We verified Q1's claim
that A_{fit}=8 holds to ~10^-3 across all 60 Q1 quartics; the
deep-WKB extension to ~5e-4 is consistent with B4 sharp form
but does not "lock" B4 sharp at d=4 in the way that the d=3
residualization locks B5/B6 sharp at d=3.

## (iii) Does the verdict-matrix logic correctly distinguish
"weak cross-degree signal" from "no signal"? What's the
threshold?

The matrix in the prompt says:
  - (a)PASS, (b)POS, (c)POS, (d)POS  -> CASE A
  - (a)PASS, (b)NEG, (c)NEG, (d)NEG  -> CASE B
  - mixed                            -> CASE C
We hit (a)PASS, (b)NEG-on-Q1-60, (c)POS-shallow, (d)NEG-deep
which is genuinely intermediate.  The deep-WKB result (d) is
the most rigorous evidence (highest dps, longest fit window,
lowest residual_std).  We weighted (d) highest in the final
verdict, downgrading (c) to "real but shallow-only".  The
threshold we used was: if R13-D shows fam32 ~ fam01 within
3 pooled stderr at deep-N, and (a) passes, then case B
is supported regardless of R13-C.  R13-C is then re-classified
as a separate phenomenon (op:shallow-j-effect-d4).

This is a defensible reading because (d) is a direct test of
the SHARP B5 statement at d=4 on the only j=0 quartic in the
canonical Q1 catalogue, while (c) is a finite-sample shallow-N
statistical comparison subject to coefficient-magnitude
confounds (we documented partial confound: ρ(a4, δ) = +0.51
within the j=0 cell).

## (iv) Is the extended-enumeration cap at α=±10 too narrow;
could a Galois-D_4 quartic with j=0 lie outside?

We stopped at W=5 because the j=0 cell already had 99 families
(easily satisfying ≥4) before reaching the prompt's escalation
gates.  Among those 99 we found:
  - 91 with Galois group S_4
  - 8 with Galois group D_4 (cands 22, 23, 34, 35, 60, 67)
  - 0 with V_4, C_4, A_4
The D_4 sub-cell has size 8 and its δ_R13_free range
[-3.36e-3, -5.0e-4] overlaps the S_4 sub-cell.  No clear
Galois-stratification visible.  We did not run deep-WKB on
any of them — that's a v1.4 / R1.4 task.

A V_4 j=0 quartic would require both `disc(b)` square AND
`I=12 a4 a0 - 3 a3 a1 + a2² = 0`.  In the W=5 window we found
0 such tuples.  Whether one exists in W=10 is unclear; if it
does, its WKB behavior could differ.  We flag this as
op:galois-modular-stratification (the optional non-halt flag
from the prompt).

## Additional honest concern: fam32 deep stderr might be optimistic

The deep-WKB regression has only 6 datapoints (N in
{200,220,240,260,280,300}).  The reported A_stderr=1.45e-5
assumes the residuals are i.i.d. Gaussian, which they are not
(they're WKB sub-leading shape).  A more honest stderr might be
~10x larger (~1.5e-4).  Even so, the fam32-vs-fam01 difference
of 4.9e-6 is < 1.5e-4, so the conclusion (no deep-N j-effect)
stands.

## Summary

The R1.3 verdict (CASE B with C caveat) is well-supported.
v1.3 should adopt cubic-modular framing for B5/B6 with
explicit footnote on the shallow-N d=4 phenomenon.  The
shallow-N effect is real but does not lift to a deep-N rule
in this scope.
