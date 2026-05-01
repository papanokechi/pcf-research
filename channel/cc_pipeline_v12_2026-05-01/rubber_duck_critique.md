# Rubber-duck critique -- CHANNEL-THEORY-V12

**Focus question (from prompt):** did we leave any v1.1 token reference to
$c(d) = 2\sqrt{(d-1)!}$ anywhere it should not appear?

## Method

`pdftotext channel_theory_outline.pdf` followed by regex grep for
`(d - 1)!`, `2 sqrt`, `4.899`, etc. PDF rendering converts
`\sqrt{(d-1)!}` to `(d − 1)!` (with the radical glyph rendered as a
preceding char), so the test is robust against that rendering.

## Where the token still appears (PDF text)

| PDF line | Context | Verdict |
|----------|---------|---------|
| 37 | Abstract: "...falsifies the v1.1 candidate $c(d) = 2\sqrt{(d-1)!}$ (which agrees with $d=2$ by coincidence and disagrees with the measured $c(4)=4$ by $\approx 22\%$)..." | Expected — required falsification statement in abstract per Phase V12-5(c). |
| 319 | §3.3 sub-section heading paragraph: "...the universality identity is in fact linear in $d$, not $2\sqrt{(d-1)!}$ as preliminarily suggested in v1.1." | Expected — required prelude per Phase V12-2(b). |
| 355–365 | Remark 3.3.E: "Version~1.1 of this outline conjectured the candidate $c(d) = 2\sqrt{(d-1)!}$... falsified at $d=4$... prediction $c(4) = 2\sqrt{3!} = 2\sqrt{6} \approx 4.899$ -- a $\approx 22\%$ disagreement..." | Expected — falsification remark itself. |
| 749 | §9 op:xi0-degree-d body: "The earlier v1.1 candidate $c(d) = 2\sqrt{(d-1)!}$ is empirically falsified at $d=4$ (Theorem~\ref{rem:c-d-falsified})." | Expected — required by Phase V12-3(a) sentence replacement. |

All four occurrences are in falsification / erratum context. There are
**zero** occurrences in:
- the §1 (Position) body,
- §3.3.A or its proof,
- §3.3 worked example for $\Vquad$,
- §5 (no-go),
- §6 (bridge),
- §7 (related work),
- §8 (implications for Master Conjecture),
- the WKB §3.4 paragraph,
- any other open problem,
- any cross-reference to PCF-1 v1.3 / PCF-2 v1.1.

## Other v1.1-vs-v1.2 consistency checks

1. **Concept DOI 10.5281/zenodo.19941678**: present in abstract paragraph
   and §10 (related work, SIARC sessions); unchanged.
2. **v1.1 record DOI 10.5281/zenodo.19941679**: explicitly marked as
   superseded in the abstract "What is new in v1.2" footnote; unchanged
   elsewhere.
3. **`\version` macro**: bumped to `1.2`; date macro pulls
   `2026-05-01 (v\version)`.
4. **TOC**: reflects new subsubsection "Cross-degree extension at
   $d=4$ (v1.2)" inserted within §3.3 (verified by reading toc).
5. **Bibliography**: `siarc_pcf2_session_q1` and
   `siarc_pcf2_v12_toappear` added; old entries unchanged.
6. **PDF page count**: 14 (within [12, 14] non-halt window).
7. **0 undefined refs / citations** across 3 pdflatex passes + bibtex.
8. **Cross-reference labels**: `prop:xi0-d4`, `conj:xi0-univ`,
   `rem:c-d-falsified` are all defined and resolved (no
   `Reference '...' on page ... undefined` in the log).
9. **Consistency with PCF-1 v1.3**: Conjecture 3.3.A* generalises
   the proven $d=2$ case (PCF-1 Theorem 5.bis: $A=2d=4$ in degree-2
   constant-$a$ row) without contradicting it; $d=2$ branch of the
   formula gives $\xi_0 = 2/\sqrt{\beta_2}$, identical to the
   v1.1 Proposition 3.3.A.
10. **Consistency with PCF-2 v1.1**: Conjecture 3.3.A* sits over PCF-2
    Conjecture B4 ($A=2d$), with Q1 confirming both at $d=4$. No
    statement of v1.2 contradicts a standing PCF-2 v1.1 statement.

## Residual risks / suggestions

- The `\version` macro is locally defined in this file only; if a
  template-wide version macro is introduced later, this should be
  reconciled.
- The `siarc_pcf2_v12_toappear` placeholder uses the cite-all DOI
  (10.5281/zenodo.19936297) of the PCF-2 stack. If PCF-2 v1.2 ends up
  with a different cite-all DOI, this entry must be rewritten before
  the corresponding v1.3 of channel-theory.
- The proof sketch of Proposition 3.3.A' is correct only for the
  $a_n \equiv 1$ degree-4 PCF (numerator constant). General-numerator
  generalisations are not in scope of v1.2 and are deferred.

**Verdict.** No residual v1.1 token escapes the erratum
context. Halt conditions all clear.
