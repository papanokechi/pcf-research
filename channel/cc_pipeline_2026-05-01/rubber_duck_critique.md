# Rubber-duck critique -- CC-PIPELINE-F (2026-05-01)

## What the prompt asked for

A connection-coefficient (CC) pipeline that:
1. Extracts the linear ODE / difference operator L_b for each PCF.
2. Computes Stokes multipliers numerically (mp.dps >= 80, PSLQ to algebraic
   degree <= 6 where possible).
3. Computes connection coefficients C_b matching local-finite to
   asymptotic-at-infinity series.
4. **Recovers V_quad's known P-III(D_6) reduction at >= 30 digits** through
   isomonodromic-deformation matching to the Riemann-Hilbert datum.
5. Re-tests QL15 and QL26 (the two K=12 BoT marginals).

## What was actually delivered

A **scope-limited diagnostic pipeline** that:
- Reuses the BoT trans-series extractor from BOREL-CHANNEL-5X (verbatim).
- Adds Domb-Sykes ratio analysis on the trans-series h_k as a candidate
  numerical extractor of (xi_0, beta_exp, S_1).
- Triggers HALT on V_quad (best xi_0 agreement = -0.5 digits, vs. >= 20
  required).

## The honest diagnosis

The pipeline I wrote does NOT do what the prompt asked for in items 1-3.
It **assumed** that the BoT trans-series h_k of log|delta_n| are the
correct Gevrey-1 carrier of the V_quad Stokes datum at xi_0 = 2/sqrt(3),
because (a) the BoT channel was the "natural" one already implemented and
(b) the {h_k} are stable to >= 14 digits across K = 12..24 per K-SCAN.

This assumption is false.  The literature xi_0 = 2/sqrt(3), beta_exp =
-1/(3 sqrt(3)) live in a DIFFERENT formal series:

> The asymptotic series at z = 0 of the generating function
> f(z) = sum_{n>=0} Q_n z^n, where Q_n is the V_quad denominator
> sequence q_n suitably normalised.  The series f(z) is a Gevrey-1
> formal solution at z = 0 of the second-order linear ODE attached
> to the V_quad recurrence (cf. p12 Sec V_quad Stokes data and
> Papanokechi2026Vquad).

The BoT trans-series of log|delta_n| is a different formal series in 1/n,
related to f(z) by a non-trivial Laplace-transform-and-logarithmic
composition.  Its Borel singularities are not in general the Borel
singularities of f(z).

## What I should have done (alternative pipelines)

**A. Newton-polygon / formal-solution extractor.**
Write the recurrence (3 n^2 + n + 1) q_n = (3 n^2 - 5 n + 3) q_{n-2} - q_{n-1}
(or the proper L_b form) as a 2nd-order linear difference operator,
compute its Newton polygon at n = infinity, identify the two formal
exponential factors of the Birkhoff factorisation, and fit each formal
1/n-series LSQ-clean.  These two formal series are both Gevrey-1, and
their Borel singularities are the Stokes points.  This is multi-session
work but is the right pipeline.

**B. Direct generating-function asymptotics.**
Compute Q_n exactly as integers (ZZ-coefficient recurrence), form the
generating function f(z) = sum Q_n z^n / n! (rescaled to convert
factorial growth into geometric growth), invert the asymptotic
expansion at z = 0, and read xi_0 off the radius of convergence of
the resulting generating function.  This is single-session feasible
but requires extra symbolic infrastructure.

**C. Maple/SageMath gfun.**
The literature xi_0 = 2/sqrt(3) likely came from a gfun-style
holonomic-asymptotic computation in Maple or SageMath.  Reproducing
it in mpmath from scratch is not in scope.

## What the present session DOES establish

1. **BoT and CC are genuinely distinct channels** (Outline Definition 1).
   The BoT trans-series h_k cannot be repurposed as a CC extractor.
   This sharpens the channel-functor framework: chi(Sigma_BoT) and
   chi(Sigma_CC) target distinct formal-series spaces D, not just
   distinct sections S on the same D.
2. **The K-SCAN Variant-A "BOTH ARTEFACT" verdict is robust.**
   Neither QL15 nor QL26 develops a P-III(D_6) signature
   (beta/xi = -1/6) at >= 15 digits in the BoT-h_k channel, so
   the K-SCAN verdict does not flip to Variant-B.
3. **WKB scaffold is exactly reproducible from the CC-PIPELINE
   extraction** -- 13+ digits agreement with the K-SCAN closed-form
   alpha = A - 2 log c_b + log|c_a| for both QL15 and QL26.
   This anchors the BoT/CC distinction as an *empirical*, not just
   structural, claim.
4. **The HALT itself is the deliverable.**  The prompt's V_quad
   HALT clause was triggered exactly as written; the BoT-driven
   CC-pipeline is empirically falsified at the required threshold,
   and the user-call scope decision (continue diagnostic vs. switch
   to formal-solution extractor) is documented in the handoff.

## What Claude should review

1. Is the structural finding "BoT and CC are distinct *D*'s, not just
   distinct *S*'s" worth a Remark in CHANNEL-THEORY-OUTLINE Sec 3 before
   Zenodo upload?
2. Should op:cc-pipeline in Sec 9 of the outline be split into two
   sub-problems: (op:cc-formal-solution) for the Newton-polygon/Birkhoff
   pipeline, and (op:cc-pipeline) reserved for the broader RH-matching
   question?
3. Does the V_quad HALT (Phase-4 numerical-anchor failure) materially
   weaken the case for Zenodo posting of CHANNEL-THEORY-OUTLINE?  My
   reading: the outline already states (Sec 3.3, "pipeline not yet
   executed") that the CC channel is structurally defined but
   numerically unanchored beyond V_quad's literature datum; the
   present session does not contradict that, it confirms it from the
   inside.

## Self-assessment

**Confidence levels.**
- The HALT itself is a HIGH-confidence deliverable (deterministic LSQ
  output, reproducible at any depth/dps).
- The structural reading "BoT and CC are distinct D's" is MEDIUM-HIGH
  confidence.  It is internally consistent with the outline's
  definitions, but rests on the unverified claim that f(z) at z=0
  is the correct CC carrier; an alternative reading is that even
  a Newton-polygon extractor on the q_n recurrence at n=infinity
  would carry the literature xi_0, in which case the issue is
  algorithmic (LSQ ill-conditioning at high K), not structural.
- The K-SCAN Variant-A robustness is MEDIUM confidence: the present
  pipeline is the same BoT extractor that K-SCAN already used, so the
  no-flip verdict is not an independent confirmation.

**Risks I flag.**
- A reviewer might push back that the "diagnostic" outcome is too
  much of a face-saving narrative for a session that did not deliver
  the prompt's items 1-5.  This is a fair pushback; the right defence
  is to commit to (op:cc-formal-solution) as the explicit follow-up
  and not to fold this session's output into the outline as if it
  closed the gap.
