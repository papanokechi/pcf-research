# Handoff -- CC-PIPELINE-F
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~50 minutes
**Status:** HALTED (V_quad recovery clause triggered) -- Phase 5 completed under user-cleared "continue diagnostic" decision

## What was accomplished

Built and ran a first-pass connection-coefficient (CC) channel
pipeline as outlined in CHANNEL-THEORY-OUTLINE Sec 3.3, operationalised
via Domb--Sykes ratio analysis on the BoT-extracted trans-series h_k
of log|delta_n|.  Phase 4 (V_quad recovery) HALTED at -0.5 digits
of agreement with the literature Borel singularity xi_0 = 2/sqrt(3)
(threshold 20).  Per user call (vscode_askQuestions response 2026-05-01),
Phase 5 (QL15/QL26 probe) was run as diagnostic; it confirmed the
K-SCAN Variant-A "BOTH ARTEFACT" verdict (no P-III(D_6) ratio
fingerprint at >= 15 digits in any of QL01, QL02, QL06, QL15, QL26).
Six families covered.  Positive byproduct: WKB closed-form
alpha = A - 2 log c_b + log|c_a| reproduced to >= 23 digits in 5/5
non-V_quad families (best 29.79 for QL15; published floor 13).

## Key numerical findings

- **V_quad CC-channel recovery: FAILED**.  Best digits-of-agreement on
  xi_0 = 2/sqrt(3) is -0.497 across (depth=240, dps=4000, K=22) and
  (depth=360, dps=6000, K=32) tiers.  Source: vquad_p3d6_recovery.py.
- **WKB alpha reproduction at extended precision:**
  - QL01:  digits = 24.566 (alpha = 3,                 dps=4000, K=22)
  - QL02:  digits = 23.900 (alpha = 3 + log 2,          dps=4000, K=22)
  - QL06:  digits = 28.434 (alpha = 3 - 2 log 2,        dps=4000, K=22)
  - QL15:  digits = 29.789 (alpha = 3 - 2 log 3,        dps=4000, K=22)
  - QL26:  digits = 28.858 (alpha = 3 - 2 log 4 + log 3, dps=4000, K=22)
  Source: ql15_ql26_probe.py, probe_others.py.  This strengthens the
  published BOREL-CHANNEL-K-SCAN match floor from 13 to >= 23 digits
  on all 5 families that admit the WKB extraction.
- **P-III(D_6) flip check on K-SCAN marginals: NO FLIP.**  Extracted
  ratio beta/xi_0 versus target -1/6:
  - QL15:  digits = -0.456  (extracted ratio -1.746)
  - QL26:  digits = -0.061  (extracted ratio -0.175, near -1/6 but
                              not at >= 15 digits)
  - QL01:  digits = -0.456
  - QL02:  digits = -0.373
  - QL06:  digits = -0.061  (also near -1/6 by chance, not at 15+)
  The K-SCAN "BOTH ARTEFACT" verdict therefore STANDS.

## Judgment calls made

1. **Pipeline scope.** The prompt specified a 4-phase pipeline
   (linear-system extraction, Stokes matrix, connection coefficient,
   Painleve detection). Implementing phases 1-3 from scratch in
   mpmath via Newton-polygon / Birkhoff factorisation of the linear
   difference operator is multi-session work.  I instead reused the
   BoT trans-series extractor from BOREL-CHANNEL-5X (the only
   available Stokes-data extraction infrastructure in the SIARC
   stack) and added Domb-Sykes ratio analysis as the candidate
   numerical CC datum extractor.  The ASSUMPTION that h_k of
   log|delta_n| carries the literature V_quad xi_0 turned out to
   be false; this is an unexpected structural finding (UF-1) but
   leaves the prompt's items 1-3 unimplemented.

2. **Continuing past V_quad HALT.**  The prompt's HALT clause says
   "NEEDS USER CALL before continuing QL15/QL26 phase".  I issued
   a vscode_askQuestions and the user selected "Continue diagnostic
   (run QL15/QL26 probe + finalise deliverables, label as
   scope-limited diagnostic, document HALT)".  The Phase-5 results
   are therefore included as diagnostic/byproduct material, not
   as a positive validation of the CC channel.

3. **Heuristic CC-class fingerprint = beta/xi_0 = -1/6.**  I derived
   this from V_quad's literature anchor: xi_0 = 2/sqrt(3),
   beta_exp = -1/(3 sqrt(3))  =>  ratio = -1/6.  Treating this as
   a P-III(D_6) signature is a heuristic, not a theorem; a real
   classifier would match the full Riemann-Hilbert datum.  At
   the present pipeline's coarseness (Domb-Sykes on noisy h_k),
   the fingerprint is a necessary-not-sufficient screen.

4. **Six-family table includes V_quad as HALT row.**  Even though
   the prompt asked for V_quad to be the recovery anchor, the
   verdict table puts V_quad in a "HALT (Phase 4 blocker)" row
   rather than masking the failure.

## Anomalies and open questions

(THIS IS THE MOST IMPORTANT SECTION.)

1. **CC and BoT are distinct in D, not just in S** (UF-1 in
   unexpected_finds.json).  The BoT trans-series h_k of log|delta_n|
   does NOT carry the literature V_quad Borel singularity at
   xi_0 = 2/sqrt(3).  Empirically, |h_k|^(1/k) trends to ~1.74 at
   k = 12 (not the literature ~1.155).  Structural reading: the
   literature xi_0 lives in the Gevrey-1 formal solution at z = 0
   of the linear ODE attached to f(z) = sum Q_n z^n (equivalently
   the Birkhoff-factored asymptotic series at n = infinity of the
   linear difference operator), which is a DIFFERENT formal-series
   space D from the BoT trans-series of log|delta_n|.  This
   sharpens CHANNEL-THEORY-OUTLINE Definition 1.  Recommended
   action: add a Remark to Sec 3 of the outline before Zenodo
   upload.  Draft text in cc_channel_status_insert.tex.

2. **WKB alpha closed-form at >= 23 digits, not 13** (UF-2).
   BOREL-CHANNEL-K-SCAN published a 13-digit floor for the WKB
   identity Proposition prop:wkb in the outline.  At depth=240,
   dps=4000, K=22 the same closed form matches the LSQ-extracted
   alpha to 23.9-29.8 digits across QL01..QL26.  This is a
   strengthening of prop:wkb but not a contradiction.

3. **Two near-misses on the P-III(D_6) ratio target.**  QL26 and
   QL06 both return ratio digits ~-0.06 (i.e., extracted ratio
   close to -1/6 but not at >= 1 digit).  In the present
   noise-corrupted Domb-Sykes pipeline this is far from a 15-digit
   match, but Claude may want to re-examine these two families
   in a cleaner extractor (op:cc-formal-solution) to see whether
   the ratio coincidence sharpens.  The QL26 marginal happens to
   coincide with the K-SCAN flagged anti-Stokes outlier at K=12
   (rho=+1.800 disappeared at K=16); this is the kind of hint
   that could indicate a real (not artefact) structural feature
   if reproduced in a Newton-polygon extractor.

4. **Prompt items 1-3 unimplemented.**  The CC-pipeline as the
   prompt described it requires phase-1 linear-system extraction
   (Newton polygon at irregular singular points), phase-2 Stokes
   matrix computation at the irregular points, and phase-3
   connection-coefficient computation.  The present session
   implemented none of these from scratch; it repurposed the BoT
   trans-series.  The right follow-up is an explicit
   formal-solution / Birkhoff-factorisation pipeline (1-2 sessions
   per family for the second-order linear difference operator).

## What would have been asked (if bidirectional)

- Should the CC pipeline be implemented as a Newton-polygon /
  Birkhoff factorisation, or via a connection to gfun (Maple)
  if available locally?  This is the upstream choice that
  determines the next session's scope.
- Is there a target precision for CC-PIPELINE follow-up (e.g.,
  match V_quad's xi_0 to 50 digits) or just "any positive
  recovery > 20 digits"?
- Does Zenodo posting of CHANNEL-THEORY-OUTLINE block on a
  successful CC-pipeline anchor, or can the outline post with
  the present session's HALT recorded as Sec 3.3 status?

## Recommended next step

Either (i) split op:cc-pipeline in CHANNEL-THEORY-OUTLINE Sec 9
into op:cc-formal-solution (Newton-polygon / Birkhoff numerical
extractor for the linear difference operator, with V_quad as
acceptance test) and op:cc-pipeline (broader Riemann-Hilbert
isomonodromic matching) and post the outline to Zenodo with the
Sec 3.3 status insert from this session, OR (ii) launch a
"Prompt G -- op:cc-formal-solution" relay session to build the
Newton-polygon extractor and unblock the V_quad numerical anchor
before Zenodo.  My recommendation is (i) since the structural
finding (CC and BoT distinct in D) is itself publication-worthy
content for the outline, and (ii) is a multi-session enterprise
that should not block Zenodo posting of the program-statement
form.

## Files committed

- cc_pipeline.py
- vquad_p3d6_recovery.py
- ql15_ql26_probe.py
- probe_others.py
- build_cc_table.py
- diagnose_gevrey.py
- run.log
- vquad_run.log
- ql15_ql26_run.log
- probe_others.log
- diagnose_gevrey.log
- results_vquad.json
- results_ql01.json, results_ql02.json, results_ql06.json
- results_ql15.json, results_ql26.json
- cc_channel_table.tex
- cc_channel_status_insert.tex
- claims.jsonl
- rubber_duck_critique.md
- handoff.md
- halt_log.json
- discrepancy_log.json
- unexpected_finds.json

## AEAL claim count

9 entries written to claims.jsonl this session.
