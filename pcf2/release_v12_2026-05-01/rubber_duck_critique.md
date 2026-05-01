# Rubber-duck critique — PCF2-V12-RELEASE

Self-critique of the v1.2 release before staging. Standard SIARC rubber
duck format: "what would a hostile reviewer ask?"

## 1. Is Conjecture B4 actually sharper, or just larger N?

**Question.** v1.1 said "verified at d=3, 50/50". v1.2 says "verified at
d∈{3,4}, 110/110". The conjecture text is identical except for the
domain expansion. Is anything actually new?

**Response.** Two things are new. First, the conjecture is now tested
in a degree where the v1.1 candidate constant c(d) = 2√((d-1)!) (Channel
Theory v1.1) and the corrected c(d) = d (this v1.2) make different
predictions: c(4) = 2√6 ≈ 4.899 vs c(4) = 4. The d=4 verification
(spread 0 across 8 representatives) is therefore not just adding
data; it actively falsifies a previously published candidate. Second,
op:b4-degree-d is materially narrowed (from "d ≥ 4" to "d ≥ 5"), so
the conjecture's open-problem boundary moves.

## 2. The d=4 catalogue has no CM quartic. Is that a credibility hole?

**Question.** At d=3 the strongest evidence cell was −_S₃_CM (37/37).
At d=4 the lex window contains 0 V₄/C₄ totally imaginary quartics.
Could the A=8 result be an artefact of the missing CM cell?

**Response.** Acknowledged in §"Caveats" of the Q1 section verbatim,
and marked as a recommended cross-cut in op:b4-degree-d: "A
complementary cross-cut is a CM-quartic harvest at d=4 [...] to rule
out an A=7 branch on the imaginary-CM cell." The current evidence is
60/60 unsplit on a mixed-S₄-dominated catalogue, which is stronger
than v1.1's 13/13 +real bins at d=3 but weaker than 37/37 −CM at d=3.
We did not over-claim.

## 3. R1 borderline signal: is it being absorbed appropriately?

**Question.** Spearman ρ = −0.45, Bonferroni p = 1.1e−2 is in the
"interesting but not theorem-grade" zone. Is the v1.2 paragraph
strengthening or weakening Conjecture B4?

**Response.** The paragraph is presented as a candidate slow-trend
invariant, with the HALT-not-triggered status stated explicitly. The
weighted-by-stderr⁻² check weakens |ρ| to 0.26, which we report.
op:finer-cubic-split is reformulated as the R1.1 follow-up (pari/gp
six invariants + precision escalation on 4 loose-fit families + Mahler
measure + curve genus). Conjecture B4 is retained; it is not
strengthened by the borderline signal. A hostile reviewer cannot
extract a "you over-claimed R1" objection from the v1.2 text.

## 4. Channel Theory v1.2 not on Zenodo — does this break the citation?

**Question.** unexpected_finds.json notes that Channel Theory v1.2's
child DOI is not yet published. Citing 10.5281/zenodo.19941678 risks
a dead-link review.

**Response.** 10.5281/zenodo.19941678 is the *concept DOI* (cite-all),
which by Zenodo semantics is preserved across child versions. v1.1 of
channel theory is already published under that concept DOI (per
submission_log.txt item 16). The bibitem and Related-identifiers
reference the concept DOI only, so the citation is live regardless of
whether the v1.2 child upload has happened. Zenodo-upload-of-v1.2 is
recorded as a follow-up user action in handoff.md.

## 5. "What is new in v1.2" paragraph — too long for a program statement?

**Question.** Program statements traditionally have minimal front
matter. The "What is new in v1.2" italicised paragraph is 12 lines.

**Response.** Defensible: the SIARC convention (PCF-1 v1.3) carries a
similar paragraph after the abstract. The length matches v1.1's
abstract closing paragraph that this v1.2 paragraph replaces; net
length change is < 6 lines. Acceptable.

## 6. Did anything in the existing v1.1 prose silently get inconsistent?

**Question.** The Open Problems section was rewritten. Are there
forward-references in the rest of the document that now broken?

**Response.** Three refs were checked: (i) the Risk register table
references op:b4-degree-d in its mitigation column — kept consistent;
(ii) §"Why this matters" references conj:B3iv — unchanged; (iii) the
abstract's closing reference to op:d2-anomaly is updated to mention
both d=3 and d=4 unsplit. pdflatex pass 3 reports 0 undefined
references, which is the canonical check.

## 7. Is "Q1-A" / "Q1-B" terminology a violation of program-statement
   discipline?

**Question.** The v1.0/v1.1 PCF-2 was a program statement, not a
results paper. v1.2 now contains a Proposition (not just Conjectures)
plus phrases like "Claim Q1-B verdict: SUPPORTED". Is this a category
violation?

**Response.** Acknowledged tension. The Newton-polygon Proposition
prop:xi0-univ is mathematically a Lemma about the Newton polygon of
L_d at z=0 — it is a derivation, not a conjecture, and including it
inside a "program statement" follows the precedent of PCF-1 v1.3
which also carries a Theorem (Theorem 5, the d=2 WKB exponent
identity). The "Q1-A/B verdict" language is bridge-style and could be
demoted to softer prose ("we observe that...") in a future v1.3. Flagged
for follow-up but not blocking.

## 8. AEAL claims hash all to the v1.2 .tex? Is that meaningful?

**Question.** All five claims in claims.jsonl carry the same
output_hash (the .tex SHA256), except claim 5 which uses the .pdf
hash. Is that an evidence-class inflation?

**Response.** Yes, this is a known SIARC convention: when the source
of truth for a numerical claim is the integrated paper (and the
underlying scripts live in upstream session folders with their own
claims.jsonl), the v1.2 release records hash = release-tex-hash to
mark "this claim was absorbed by this release". The original
session-level claim hashes (Q1: claims_phaseQ1_*.jsonl; R1:
claims.jsonl) are preserved in the upstream folders. This is the same
pattern used in PCF2-V11-RELEASE. Acceptable.

## Summary

No HALT issues. Two non-blocking flags recorded in unexpected_finds.json
(Channel Theory v1.2 Zenodo upload pending; v1.2 child DOI placeholder
in submission log patch). Conjecture B4 is sharpened correctly; the d=4
result actively falsifies a previously published candidate constant
and is therefore non-trivially new. R1 is presented as borderline and
correctly does not promote to a theorem. Document compiles in 18 pp,
within the prompt's [14, 18] page-count band.
