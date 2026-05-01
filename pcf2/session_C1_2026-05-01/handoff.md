# Handoff — PCF2-SESSION-C1
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~25 minutes (compute: 5.1 s)
**Status:** COMPLETE (with substantive empirical falsification — see below)

## What was accomplished
Completed the WKB-exponent harvest for the 13 non-CM cubic families
(bins `+_S3_real` and `+_C3_real`) in the PCF-2 cubic catalogue,
bringing degree-3 coverage to 50/50 (combined with Session B's 37 CM
families). Tested refined Conjecture B4'. **Outcome: B4' is empirically
falsified at d=3.** All 13 non-CM families have A_fit ≈ 6, not the
predicted A = 2d - 1 = 5. The surviving cubic conjecture is the
unsplit B4 sharp form A = 2d = 6 across all 50 cubic families. Updated
the Session-B paper draft to install B4' (per prompt) with a
falsification remark and a new §"WKB-exponent harvest, full degree-3
coverage (Session C1)".

## Key numerical findings
- 13/13 non-CM families pass the 30-digit L-stability gate at n=200
  (in fact all reach the dps=80 floor by n=150 → ≥80 stable digits).
  WKB fit at dps=800 over n∈[10,100] step 3, N_ref=300.
  Source: `session_c1_wkb.py` (this session).
- A_fit range: [5.884, 5.989], mean 5.973, stddev 0.026.
  Predicted by B4' for these bins: A = 5. Predicted by original B4: A∈{5,6}.
- Per-family A_fit (4 dp) and verdict:
  | id | bin | A_fit ± stderr | verdict (vs predicted=5) |
  |----|-----|----------------|--------------------------|
  | 4  | +_S3_real | 5.9719 ± 0.0018 | wrong-branch (≈6) |
  | 5  | +_S3_real | 5.9715 ± 0.0016 | wrong-branch (≈6) |
  | 9  | +_S3_real | 5.9774 ± 0.0017 | wrong-branch (≈6) |
  | 10 | +_S3_real | 5.9770 ± 0.0014 | wrong-branch (≈6) |
  | 14 | +_S3_real | 5.9813 ± 0.0012 | wrong-branch (≈6) |
  | 15 | +_S3_real | 5.9818 ± 0.0012 | wrong-branch (≈6) |
  | 19 | +_C3_real | 5.9860 ± 0.0009 | wrong-branch (≈6) |
  | 20 | +_C3_real | 5.9870 ± 0.0008 | wrong-branch (≈6) |
  | 37 | +_S3_real | 5.8836 ± 0.0479 | ambiguous (interior) |
  | 38 | +_S3_real | 5.9806 ± 0.0012 | wrong-branch (≈6) |
  | 39 | +_S3_real | 5.9810 ± 0.0011 | wrong-branch (≈6) |
  | 42 | +_S3_real | 5.9850 ± 0.0009 | wrong-branch (≈6) |
  | 46 | +_C3_real | 5.9889 ± 0.0007 | wrong-branch (≈6) |
- A_fit histogram (bin width 0.025, predicted=5):
    [5.875, 5.900): 1 (family 37)
    [5.950, 5.975): 4
    [5.975, 6.000): 8
  Zero families in [4.90, 5.10].
- Family 46 (the conductor-7 anchor cited in the C1 prompt as
  A_fit = 5.95) recomputes to A_fit = 5.989 ± 0.001 at dps=800. The
  prompt's anchor value 5.95 is itself closer to 6 than to 5; B4' was
  effectively falsified by its own anchor before C1 began.
- Combined 50-family cubic harvest:
    - 37 CM (Session B): mean A 5.982, sigma 0.025
    - 13 non-CM (C1):    mean A 5.973, sigma 0.026
    - 50/50 within 0.15 of A = 6, none within 0.10 of A = 5.

## Judgment calls made
- **dps for the WKB fit.** The C1 prompt specifies "mp.dps = 80 work"
  but in the same paragraph says "identical to Session~B step 1b".
  At dps = 80 the |δ_n| values for n ≳ 28 (with A ≈ 6) underflow and
  the fit collapses. I used Session B's dps = 800. Documented in
  `rubber_duck_critique.md` §3.
- **Off-spec verdict label.** The three verdicts in the prompt
  (consistent / violating / ambiguous) do not partition cases where
  A_fit is close to the *wrong* branch (within 0.10 of 6 when
  predicted is 5). I added a fourth label
  `b4_prime_consistent_wrong_branch` for these 12 families; they are
  also flagged in `unexpected_finds.json` per the prompt's NON-HALT
  FLAG clause. No halt logic affected. Documented in
  `rubber_duck_critique.md` §2.
- **post-condition installation in spite of falsification.** The
  prompt's POST-CONDITION says install B4' as Session B's stated
  conjecture "regardless of outcome". Done; immediately followed by a
  bold "Update (PCF-2 Session C1, ...)" remark inside the conjecture
  environment, plus a new §sec:C1-update spelling out the
  falsification and the surviving sharp form A = 2d = 6 on 50/50.
- **Unexpected-finds saturation.** Each non-CM family triggered both
  flags (closer-to-6 and ≥50-digit convergence at n=200), giving 26
  entries for 13 families. The "≥50 digits" flag is a dps-floor
  artefact (the 999.0 value caps `stable_digits` at the floor, which
  for dps=80 is ≥80 digits). Honest reading: gate was non-binding for
  all 13. Not corrected; documented in `rubber_duck_critique.md` §1.

## Anomalies and open questions
1. **B4' empirically falsified at d=3.** This was the explicit test
   of the session and is the headline result. The catalogue's
   trichotomy bin does **not** split A_fit at d=3, in contrast to
   PCF-1 v1.3 Theorem 5 where the discriminant sign **does** split A
   at d=2 (A ∈ {3, 4}). Open question: is the d=2 split a low-degree
   artefact (indicial-root resonance), or is there a degree-d threshold
   beyond which bin-dependence reappears? A Session C2 quartic probe is
   the natural next step.
2. **Family 37 outlier behaviour.** A_fit = 5.884 ± 0.048 stands out
   (stderr ~30× the other 12). The family is `n^3 - 2n^2 - 3n + 1`
   with Δ_3 = 257 (shared with family 9, which fits cleanly at
   5.9774). Worth a targeted re-fit at higher dps (say 1600) or a
   wider n-window to check whether the larger stderr is genuine
   (e.g. a slowly-varying β or γ correction) or an artefact of the
   particular n-grid.
3. **Self-inconsistency in the C1 prompt.** The prompt presented family
   46 as "calibration anchor at A_fit = 5.95" while predicting branch
   A = 5 for `+_C3_real`. 5.95 is closer to 6 than to 5; the anchor
   itself contradicts the conjecture being tested. The PCF2 Session A2
   commit message explicitly recorded "A_fit=5.95 -> 6=2d" — so the
   record was internally consistent and the prompt's framing of B4' was
   the bit out of date. Logged.
4. **L-stability gate calibration.** The 30-digit gate is wholly
   non-binding for non-CM cubic L-values at dps = 80 (all reach ≥80
   stable digits by n = 200). The gate would only bite at very small
   N or much larger d. The "striking convergence" non-halt flag
   triggered for all 13 families because of this; treat as
   instrumentation, not as a substantive flag.
5. **No PSLQ in C1 (per spec).** The 13 non-CM L-values are recorded
   to ≥80 digits at n = 200 in `results.json`. A future targeted PSLQ
   probe could test, e.g., whether `+_C3_real` L-values lie in the
   ring of integers of the cyclotomic field K_7 (since families 19/20
   and 46 all have Δ_3 ∈ {49, 81} = {7², 3⁴}).

## What would have been asked (if bidirectional)
- Should the prompt's POST-CONDITION still install B4' verbatim as
  the conjecture statement, given the data falsifies it?
  (I did install B4' but with an inline falsification remark and a
  new C1 update section — not just a silent install. If you'd prefer
  the conjecture environment to instead state the surviving "A = 2d
  unsplit" form with B4' demoted to a remark on a refuted hypothesis,
  one further pass on Session B's .tex would do it.)
- Should family 37's wider stderr be re-fit at higher dps before
  treating it as in-band?

## Recommended next step
**Session C2 — quartic catalogue probe.** Build a small quartic
catalogue (say 20–30 families along the same lex window) and run the
identical WKB pipeline. If A_fit ∈ {7, 8} unsplit, the d=2 bin-split
is confirmed as a low-degree artefact. If A_fit splits by some bin,
we have a new candidate for the right "channel-classifying" structural
invariant (likely tied to indicial-root resonance rather than to the
trichotomy bin). Either outcome is publishable.

## Files committed
- `session_c1_wkb.py`         (script, 22.4 KB)
- `run.log`                   (terminal trace)
- `stdout.log`                (mirror of run.log via Tee-Object)
- `results.json`              (13-family per-family record + summary)
- `wkb_noncm_table.tex`       (LaTeX table, 13 rows, for v1.1)
- `wkb_cubic_harvest_v2.tex`  (combined 50-row table; replaces Session B)
- `claims.jsonl`              (1 aggregate AEAL claim — see notes)
- `rubber_duck_critique.md`   (this session's self-critique)
- `handoff.md`                (this file)
- `halt_log.json`             (empty `{}` — no halt triggered)
- `discrepancy_log.json`      (empty `{}`)
- `unexpected_finds.json`     (26 entries for 13 families)

## AEAL claim count
**1** entry written to `claims.jsonl` this session (the aggregate
50-family Conjecture B4' status claim). No per-family
`b4_prime_consistent` claims, because zero families were consistent
with the predicted A = 5 branch — the falsification is the headline
finding and is recorded as the aggregate.
