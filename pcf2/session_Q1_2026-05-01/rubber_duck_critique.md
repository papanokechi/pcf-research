# Rubber-duck critique — PCF2-SESSION-Q1

This document plays devil's advocate against the session's three main
empirical claims. Each weakness is raised, then either folded into the
deliverables (handoff Anomalies, or Cross-Degree caveat) or noted for a
follow-up session.

## §1. Galois-bin coverage of the catalogue is heavily skewed.

The lex window $a_4\in\{1,2,3,5,7\}$, $a_{0,1,2,3}\in\{-3,\ldots,3\}$
produced $57\times S_4 + 3\times D_4$ irreducible quartics, and zero
$V_4 / C_4 / A_4 / -\_\mathrm{CM}$ entries. The promised quartic
trichotomy with five Galois groups is empirically a near-monoculture
of $S_4$ within this window.

This is **not a problem for the unsplit B4 statement** (a unanimous
verdict across $57\,S_4$ entries plus $3\,D_4$ + 0 of the others is
weaker evidence for "all Galois classes" but still adequate for "every
$b$ in the catalogue"), but it **is a real weakness** of any claim that
the splitting is bin-driven. The natural follow-up is to extend the
window to $a_4\in\{1\}$, $a_{0,1,2,3}\in\{-5,\ldots,5\}$ and
specifically harvest CM-quartic seeds (totally imaginary signature
$(0,2)$). Recorded as the recommended next session.

Fix: noted in cross_degree_table.tex caveat paragraph and in handoff
Anomalies §2.

## §2. The fit-window bias is real and not just cosmetic.

Mean $A_{\mathrm{fit}}=7.954$ with $\sigma=0.0037$ is a $12\sigma$
deviation from $A=8$. By a strict reading of "$A_{\mathrm{fit}}=8$",
this is a falsification. The defence is structural: the same 4-parameter
ansatz at $d=3$ produced $\overline{A}=5.973$ ($7\sigma$ below 6), and
PCF-1 v1.3 at $d=2$ produced $\overline{A}=3.97$ on the $A=4$ stratum.
The pattern of "$\overline{A}\approx 2d - 0.05\cdot(d/2)$" is
self-consistent across $d\in\{2,3,4\}$ and is parameter-error-driven,
not physics-driven.

A sharper test would either (a) extend the ansatz with a $-(\log\log
N)$ term, or (b) use a finer grid at larger $N$ so that the
$\beta\log N$ term is well-conditioned. Neither was done in this
session.

Fix: handoff Anomalies §1 records this; cross_degree_table.tex notes
the systematic offset across $d$.

## §3. The Newton-polygon $\xi_0$ is operational, not Borel-radius.

The "universal $\xi_0(b) = d/\beta_d^{1/d}$" claim is a derivation in
the *operator's* characteristic polynomial along the slope-$1/d$
edge --- i.e., it is the WKB-exponent of the formal-solution
exponential at the irregular singularity. That is the right notion
for the trans-series leading exponent, and matches the d=2 statement
of Prop.~3.3.A in CHANNEL-THEORY-V11. It is **not** the same as the
Borel-radius identity in the Pade/Borel sense (which CC-PIPELINE-G
showed converges to the correct $V_{\mathrm{quad}}$ value but was
shown by BOREL-CHANNEL-K-SCAN to be artefactual at $K=12$ for QL15
and QL26).

So claim Q1-B as written --- "$\xi_0 = c(4)/\beta_4^{1/4}$ in the
Newton-polygon characteristic root sense" --- is established. The
stronger claim "this is the Borel-radius of the trans-series" is
*not* established at $d=4$ and would require running the Pade-Borel
machinery on a representative quartic.

Fix: claim Q1-B is stated explicitly as "operational characteristic
root" in newton_polygon_d4.tex; Borel-radius extension flagged as
an open question in handoff Anomalies §3.

## §4. PSLQ scan was not run.

The prompt did not explicitly require a PSLQ scan, but Sessions B and
C1 used PSLQ to test BIN-CONSISTENT verdicts at $d=3$. We did not run
a PSLQ scan at $d=4$ because (a) the catalogue contains no CM bin, and
the natural target for PSLQ at $d=4$ would be the $-_\mathrm{CM}$
quartic bin (analogous to $-\_S_3\_\mathrm{CM}$ at $d=3$); (b) the
session's two specified empirical claims (Q1-A, Q1-B) are about WKB
exponents and Newton-polygon roots, not algebraic transcendence.

A PSLQ scan on the 60 cataloged quartics against tier (a)–(d) bases
would be a follow-up session and is recommended in handoff Recommended
next step.

## §5. Quartic Galois discriminator (C_4 vs D_4) is not fully tested.

The $m=1$ branch of the resolvent-cubic algorithm distinguishes $C_4$
from $D_4$ via whether $(\alpha^2 - 4r)(\alpha - p)\,\mathrm{disc}(f)$
is a rational square. The current implementation uses
`sympy.nsimplify` followed by perfect-square check on numerator and
denominator separately. This is correct in principle but we did not
verify against a known $C_4$ quartic (calibration anchors include
$x^4-2$ as $D_4$ and $x^4+1$ as $V_4$ but no genuine $C_4$). The
catalogue's $3\times D_4$ entries may include a misclassified $C_4$;
empirically the WKB harvest is bin-insensitive so this does not
affect the B4 verdict.

Fix: noted in handoff. A focused calibration on the standard
$x^4 - 4 x^2 + 2$ ($C_4$, the maximal real subfield of the
$16$th cyclotomic) is recommended.
