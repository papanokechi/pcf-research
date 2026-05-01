# Handoff — PCF2-V12-RELEASE
**Date:** 2026-05-01
**Agent:** GitHub Copilot (VS Code)
**Session duration:** ~45 minutes
**Status:** COMPLETE

## What was accomplished
Produced PCF-2 v1.2 source and PDF, integrating Session Q1 (quartic
harvest at d=4) and Session R1 (finer-cubic-split correlation probe).
Conjecture B4 sharpened from "verified at d=3, 50/50" to
"verified at d∈{3,4}, 110/110 jointly". Universal Newton-polygon
characteristic root xi_0(b)=d/beta_d^(1/d) (c(d)=d) added as
Proposition (`prop:xi0-univ`), recovering the proven d=2 case and
falsifying the v1.1 Channel-Theory candidate c(d)=2*sqrt((d-1)!).
R1 borderline signal on log|Delta_3| recorded in a new subsection;
HALT not triggered; B4 retained. Open problems updated: op:b4-degree-d
narrowed to d>=5; op:finer-cubic-split reformulated as R1.1
follow-up. Channel Theory cite-all DOI 10.5281/zenodo.19941678 added
to references. v1.1 archived to
`tex/submitted/archive/pcf2_program_statement_v1.1.tex`.

## Key numerical findings
- **Conjecture B4 sharpened**: A=2d verified at d∈{3,4} on 110/110
  jointly harvested families. Cubic mean A_fit=5.978, σ=0.026
  (Sessions B+C1, dps=800). Quartic mean A_fit=7.954, σ=0.0037
  (Session Q1, dps=1200, fit window N∈[10,60]).
- **Universal xi_0-identity**: c(d)=d verified at d=4 to spread 0
  (relative <1e-50) across 8 representative quartics; falsifies
  v1.1 candidate c(4)=2*sqrt(6)~4.899.
- **R1 finer-cubic-split**: log|Delta_3| flagged as borderline
  candidate slow-trend invariant. Spearman ρ=-0.451 (full 50,
  Bonferroni p=1.1e-2 over m=11), ρ=-0.492 (excluding fam 31,
  Bonferroni p=3.7e-3); weighted |ρ|_w=0.26. R1 HALT threshold
  (|ρ|>0.6 at Bonferroni p<1e-3) NOT triggered.
- **PCF-2 v1.2 PDF**: 18 pages, compiled cleanly in 3 pdflatex
  passes, 0 errors, 0 undefined references. Within prompt's
  page-count band [14, 18].

## Judgment calls made
1. **Halt condition #1 not triggered.** Halt 1 reads "Channel Theory
   v1.2 has not yet landed on Zenodo". The v1.2 PDF exists locally
   and the prompt itself supplies cite-all DOI 10.5281/zenodo.19941678
   (= concept DOI, preserved across child versions). PCF-2 v1.2
   cites only the concept DOI, so the citation is live whether or
   not the v1.2 child upload has happened. Recorded the
   Zenodo-upload-of-CT-v1.2 as a follow-up user action in
   `unexpected_finds.json` rather than halting. The conservative
   alternative would have been to halt and demand explicit
   confirmation; I judged the cite-all-DOI semantics make the
   citation correct as written.
2. **Proposition env added to preamble.** v1.1 had no proposition
   environment (only conjecture, problem, definition, remark). The
   Q1-B universal xi_0 result is a *derivation*, properly framed as
   a Proposition rather than a Conjecture. Added one line:
   `\newtheorem{proposition}[conjecture]{Proposition}`. Counter
   shared with `conjecture` so the running-numbering is preserved.
3. **R1 paragraph absorbed as a `\subsection` of the B4 section,
   not as a top-level section.** The R1 content is logically
   commentary on Conjecture B4 (it is a finer-split investigation,
   not a separate result), so it sits inside the B4 section as
   `\subsection{Finer cubic split: log|Delta_3| as candidate
   slow-trend invariant (Session R1)}` (label `sec:R1-finer-split`).
   The Q1 quartic content, by contrast, is a top-level new section
   `sec:Q1-quartic` because it introduces a new degree.
4. **submission_log.txt patch staged but NOT applied in-place.**
   The submission_log.txt is the user's authoritative record;
   editing it requires the v1.2 child DOI which is minted only at
   Zenodo Publish-time. Produced `submission_log_patch_item15.txt`
   with explicit `<TBD-v1.2>` placeholder for the user to apply
   post-upload.
5. **No tex/submitted/ git push.** Per Standing Final Step, only
   the bridge and pcf-research mirrors are pushed. The
   `tex/submitted/pcf2_program_statement.{tex,pdf}` is updated
   in-place but not committed/pushed by this session, in keeping
   with the 'do not push the PDF until Zenodo upload is finalised'
   note from PCF2-V11-RELEASE handoff §3.
6. **Series-label "1.1" → "1.2" updated in both `\title` and
   `\date{...(v1.1)}` directives.** Verified the rendered PDF cover
   page reads "Version 1.2" and date "May 1, 2026 (v1.2)".

## Anomalies and open questions
1. **Channel Theory v1.2 Zenodo upload pending (user action).**
   Local v1.2 PDF exists at
   `sessions/2026-05-01/CHANNEL-THEORY-V12/channel_theory_outline.pdf`.
   Cite-all DOI 10.5281/zenodo.19941678 is preserved by Zenodo
   concept-DOI semantics, so the PCF-2 v1.2 reference is live, but
   the v1.2 *child* DOI of channel theory is not yet minted.
   **User action recommended.** Document does not block on this.
2. **Q1 lex window contains no CM quartic.** The 60-family d=4
   catalogue is dominated by mixed-S₄ extensions and contains 0
   V₄/C₄ totally imaginary quartics. The 60/60 b4_consistent_at_8
   verdict therefore does not yet rule out an A=7 branch on the
   imaginary-CM cell. Documented in §"Caveats" of the Q1 section
   and in op:b4-degree-d (which now suggests a CM-quartic cross-cut
   alongside the d=5 quintic harvest).
3. **R1.1 follow-up requires `pari/gp` install.** Six number-field
   invariants (h, h⁺, conductor, regulator, fund-unit-norm, unit
   density) deferred from R1 because pari/gp was not available in
   the agent's environment. Mahler measure and curve-genus
   invariants are also pending. Recorded in op:finer-cubic-split.
4. **No retraction of v1.1 conjectures.** v1.2 absorbs Q1 and R1
   without retracting any v1.1 statement. Conjecture B4
   (unsplit, A=2d) is the same conjecture, now verified at one
   more degree. B4' remains empirically falsified at d=3 and is
   now also empirically empty at d=4 (60/60 b4_consistent_at_8,
   no A=7 branch). No conjecture flip.

## What would have been asked (if bidirectional)
- Should we record v1.1 of channel theory (DOI 19941679) *and* v1.2
  in the bib, or only the cite-all concept DOI? (Chose: only
  concept DOI; preserves provenance and is the minimal correct
  reference.)
- Should the R1 paragraph be promoted to its own section
  ("Finer-cubic-split probe (R1)") rather than embedded as a
  subsection of B4? (Chose: subsection. R1 is commentary on B4.)
- Should the d=4 lex window be extended in v1.2 to include CM
  quartics, or left to a follow-up session? (Chose: follow-up.
  Extending the catalogue inside the v1.2 release window would
  require a second compute pass and would push session length
  beyond the prompt's 45–60-minute envelope.)

## Recommended next step
**Channel Theory v1.2 Zenodo upload.** Press Publish on a fresh
"New version" of concept DOI 10.5281/zenodo.19941678 with
sessions/2026-05-01/CHANNEL-THEORY-V12/channel_theory_outline.pdf
as the file. After the v1.2 child DOI is minted, update item 16
of submission_log.txt to record both the v1.1 (19941679) and v1.2
child DOIs. Concurrent: PCF-2 v1.2 Zenodo upload as a "New version"
of concept DOI 10.5281/zenodo.19936297, replacing item-15 v1.1 line
in submission_log.txt with the patch in
`submission_log_patch_item15.txt` after substituting the minted v1.2
child DOI.

After both uploads land:
1. **PCF-2 R1.1 follow-up** (op:finer-cubic-split): install
   pari/gp, add the six number-field invariants to the R1
   correlation battery, precision-escalate the four loose-fit
   cubic families.
2. **PCF-2 Session Q2 quartic-CM cross-cut**
   (op:b4-degree-d, complement): enumerate V₄/C₄ totally imaginary
   quartics (Z-primitive irreducible) in a wider lex window and
   run the WKB harvest, ruling out (or revealing) an A=7 branch
   on the imaginary-CM cell.

## Files committed (this session, sessions/2026-05-01/PCF2-V12-RELEASE/)
- `pcf2_program_statement.tex`              (v1.2 source)
- `pcf2_program_statement.pdf`              (v1.2 PDF, 18 pp, 0 errors)
- `archive/pcf2_program_statement_v1.1.tex` (preserved v1.1 source)
- `zenodo_description_v1.2.txt`             (v1.0+v1.1+v1.2 changelog)
- `zenodo_notes_v1.2.txt`                   (single-line Notes field)
- `submission_log_patch_item15.txt`         (item-15 replacement block,
                                             with <TBD-v1.2> placeholder)
- `claims.jsonl`                            (5 AEAL claims)
- `rubber_duck_critique.md`                 (self-critique)
- `handoff.md`                              (this file)
- `halt_log.json`, `discrepancy_log.json`   (both empty `{}`)
- `unexpected_finds.json`                   (3 non-blocking notes
                                             on Zenodo upload status)

## AEAL claim count
5 entries written to claims.jsonl this session (Conjecture B4
sharpened d∈{3,4}, c(d)=d cross-degree at d∈{2,4}, d=2 anomaly
status, log|Delta_3| candidate finer-cubic-split, v1.2 supersedes v1.1).
