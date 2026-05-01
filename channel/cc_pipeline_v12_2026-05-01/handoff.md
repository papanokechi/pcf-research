# Handoff -- CHANNEL-THEORY-V12
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~35 minutes
**Status:** COMPLETE

## What was accomplished

Surgical erratum revision of the v1.1 channel-theory outline driven by
PCF2-SESSION-Q1. Produced v1.2 (`channel_theory_outline.tex` +
compiled PDF, 14 pp) that corrects the empirically falsified v1.1
candidate $c(d) = 2\sqrt{(d-1)!}$ to the linear form $c(d) = d$
(i.e. $\xi_0(b) = d/\beta_d^{1/d}$), folds in the $d=4$ direct
Newton-polygon confirmation (spread 0 across 8 representatives at
$\mathrm{dps}=80$), and opens `op:xi0-d3-direct` to close the
remaining $d=3$ verification. Concept DOI 10.5281/zenodo.19941678
preserved; v1.1 (10.5281/zenodo.19941679) superseded. Bibliography,
abstract, "What is new in v1.2" footnote, §3.3 (new
Proposition 3.3.A', Conjecture 3.3.A*, Remark 3.3.E),
§3.4 WKB cross-reference, §10 SIARC-sessions list, and §9
open-problems block all updated. Submission_log patch fragment for
`tex/submitted/submission_log.txt` item 16 produced.

## Key numerical findings

- $c(4) = 4$ with spread $0$ at $\mathrm{dps}=80$ across $8$ quartic
  representatives (PCF2-SESSION-Q1, `session_q1_newton.py`,
  `newton_polygon_d4.tex`).
- v1.1 candidate $c(d) = 2\sqrt{(d-1)!}$ predicts
  $c(4) = 2\sqrt{6} \approx 4.899$; falsified by $\approx 22\%$.
- v1.2 corrected formula: $\xi_0(b) = d/\beta_d^{1/d}$
  (Conjecture 3.3.A*); proven $d=2$ (Proposition 3.3.A), verified
  $d=4$ (Proposition 3.3.A'), $d=3$ deferred to op:xi0-d3-direct.
- PDF compile: 3 pdflatex passes + bibtex, 14 pages, 0 undefined
  refs/citations.

## Judgment calls made

- **`op:b4-degree-d` cross-reference.** The prompt asked to "update
  op:b4-degree-d ... narrow the open problem to d ≥ 5", but
  `op:b4-degree-d` is a PCF-2 problem, not a channel-theory problem
  (no `op:b4-degree-d` block exists in this outline). Resolved by
  cross-referencing PCF-2's `op:b4-degree-d` in the §9 intro
  paragraph and in the WKB §3.4 cross-reference, noting
  PCF2-SESSION-Q1's $A=2d=8$ confirmation $60/60$ at $d=4$ and
  thus the narrowing to $d \ge 5$.
- **Page count = 14.** The non-halt window is $[12, 14]$, exactly
  hit. The "What is new in v1.2" content is a single paragraph
  inside §1 (per Phase V12-5(d) front-matter option) rather than a
  separate appendix; if v1.3 needs more space, splitting to
  appendix is the recommended action.
- **Placeholder bibliography entry.** Added
  `siarc_pcf2_v12_toappear` with the cite-all (concept) DOI of the
  PCF-2 stack (10.5281/zenodo.19936297) and `note = {Cite-all
  (concept) DOI; v1.2 to appear.}`. Not actually cited in the body
  to avoid an undefined-reference risk if PCF-2 v1.2 is delayed.
- **`prop:xi0-d4` proof.** Sketch covers $a_n \equiv 1$
  degree-$4$ PCF only; general-numerator quartic case is deferred,
  consistent with PCF2-SESSION-Q1's scope.
- **`\version` macro.** Defined locally for v1.2 only; consistent
  with v1.1's manual `(v1.1)` date string.

## Anomalies and open questions

- The exact rational $a_2$ in `siarc_cc_pipeline_g` carried in
  prop:noninj is "$Q_{L01}$: $a_2 = 17/512$; $Q_{L02}$:
  $a_2 = -239/512$"; this is a v1.1 statement carried unchanged
  through v1.2 and was not re-verified in this session. If
  `siarc_cc_pipeline_g` records different fractions, the v1.2
  Proposition 3.3.B has the same defect as v1.1 (out of scope of
  this erratum revision).
- The "Borel-Laplace caveat" in Theorem 3.3.D is unchanged; v1.2
  does not advance op:cc-formal-borel.
- v1.2 introduces *one* new conjecture (Conjecture 3.3.A*) and one
  new proposition (Proposition 3.3.A'); their cross-degree
  consistency with PCF-1 v1.3 Theorem 5.bis ($A = 2d$ at $d = 2$
  constant-$a$ row) and PCF-2 v1.1 Conjecture B4 ($A = 2d$ at
  $d = 3$, $50/50$) is verified by inspection but was not
  separately recomputed.
- Q1 measured $A_{\mathrm{mean}} = 7.954$ with $\sigma = 3.7
  \times 10^{-3}$ at $d=4$ (slight under-shoot of $A = 2d = 8$),
  flagged in PCF2-SESSION-Q1's own handoff but not a halt; v1.2
  does not re-flag.

## What would have been asked (if bidirectional)

- "Should op:xi0-d3-direct be staged as its own relay prompt now,
  or batched with PCF-2 v1.2's cubic-and-quartic catalogue
  follow-up?"
- "Is the cite-all DOI 10.5281/zenodo.19936297 stable as the PCF-2
  cite-all, or is the cite-all generated per-version on Zenodo?"

## Recommended next step

Issue relay prompt **CHANNEL-XI0-D3-DIRECT**: 1–2 hour Newton-polygon
test of $\xi_0 = 3/\beta_3^{1/3}$ at $d=3$. One cubic representative
per Galois bin from PCF-2 v1.1 (3 families: $+_{S_3}$ real,
$+_{C_3}$ real, $-_{S_3}$ CM). Adapt
`sessions/2026-05-01/PCF2-SESSION-Q1/session_q1_newton.py` to the
cubic operator $L_3 = 1 - z B_3(\theta+1) - z^2$, measure $\xi_0$ at
$\mathrm{dps} \ge 80$, compare with $c(3) = 3$. Resolves
op:xi0-d3-direct and fully closes Conjecture 3.3.A* at low degree.

## Files committed

- channel_theory_outline.tex
- channel_theory_outline.pdf (14 pages, 511,499 bytes)
- channel_theory_outline.aux / .bbl / .blg / .log / .out / .toc (build
  artefacts; OK to retain or strip per bridge convention)
- annotated_bibliography.bib
- archive/channel_theory_outline_v1.1.tex
- zenodo_description_v1.2.txt
- zenodo_notes_v1.2.txt
- submission_log_patch_item16.txt
- claims.jsonl (5 AEAL entries)
- rubber_duck_critique.md
- handoff.md
- halt_log.json (empty `{}`)
- discrepancy_log.json (empty `{}`)
- unexpected_finds.json (1 entry: v1.1 candidate falsified at d=4 --
  expected per session brief, no halt)

## AEAL claim count

5 entries written to `claims.jsonl` this session.

## Zenodo "New version" upload checklist

1. Open existing concept DOI 10.5281/zenodo.19941678 / record
   10.5281/zenodo.19941679 (v1.1).
2. Click **"New version"**. Title: keep verbatim.
3. Replace `channel_theory_outline.pdf` (and source ZIP if used)
   with v1.2 build.
4. Replace description with `zenodo_description_v1.2.txt` (verbatim
   body).
5. Replace version notes with `zenodo_notes_v1.2.txt` (verbatim
   body).
6. Bump version field: `1.2`. Publication date: 2026-05-01.
7. Related identifiers: keep all v1.1 links; add
   `Is documented by → https://github.com/papanokechi/siarc-relay-bridge/tree/main/sessions/2026-05-01/PCF2-SESSION-Q1/`.
8. Concept DOI must remain 10.5281/zenodo.19941678 (Zenodo enforces
   this on "New version"; verify on publish-confirm dialog).
9. Publish; record the new v1.2 record DOI in
   `tex/submitted/submission_log.txt` per
   `submission_log_patch_item16.txt`.

## pdflatex log diff vs v1.1

```
v1.1: Output written on channel_theory_outline.pdf (12 pages, 434914 bytes).
v1.2: Output written on channel_theory_outline.pdf (14 pages, 511499 bytes).
```
Both runs: 0 undefined refs, 0 undefined citations, 0 errors. The
+2 page growth is exactly the §3.3 v1.2 insertion (Proposition
3.3.A', proof, Conjecture 3.3.A*, Remark 3.3.E) plus the new
op:xi0-d3-direct in §9.
