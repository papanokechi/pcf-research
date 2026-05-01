# PCF2-SESSION-R1.2 — Rubber-duck critique

This critique is structured against the four focus questions in the
relay prompt.

## (i) Is the genus-1 reduction at d=4 the right analogue of d=3?

At d=3, the curve $C_b\colon y^2 = b(x)$ with $b$ a monic squarefree
cubic IS itself an elliptic curve $E_b$ (the affine point $(0, \sqrt{b(0)})$
gives a non-singular Weierstrass model after standard transformation).
At d=4, $C_b$ is a smooth genus-1 curve but is NOT canonically an
elliptic curve over $\Q$ (no marked rational point in general). The
Jacobian $\mathrm{Jac}(C_b)$ IS an elliptic curve over $\Q$, with
j-invariant $j = 6912 \cdot I^3 / (4 I^3 - J^2)$ computable from the
Igusa $I, J$ invariants of the binary quartic.

**Verdict:** the Jacobian is the correct analogue. There is no
conceptual concern with the $d=4$ identification.

**However**, a reasonable secondary worry: even though $E_b$ exists
abstractly, the *connection* between $E_b$ and the WKB asymptotic of
$\mathrm{CF}(b)$ may be more indirect at $d=4$ than at $d=3$. At $d=3$,
the period integrals on $E_b$ literally appear in the WKB phase via the
contour integral $\oint \sqrt{b(x)}\,dx$. At $d=4$, the analogous
integral $\oint dx / \sqrt{b(x)}$ is an *abelian integral on $C_b$*
which decomposes through $\mathrm{Jac}(C_b)$ but is not literally a
period of $E_b$ on its own. This may explain why the $j$-signal is
weaker (or invisible) at $d=4$.

**Open problem:** is the right $d=4$ invariant the *Mahler measure of
$E_b$* (i.e., $h(j)$, height of the j-invariant), or the *conductor*
(an arithmetic invariant rather than analytic), or some Igusa-J
analogue independent of $I$?

## (ii) Is the $j=0$ cell big enough at $d=4$ ($n\ge 2$) to credibly confirm B5?

NO. The $j=0$ cell at $d=4$ in the current Q1 catalogue contains
exactly $n_4 = 1$ family (family 32, $b(n) = n^4 - 3n^3 - 3n^2 + 3n - 3$).
This is below the $n \ge 2$ threshold the prompt set for "credible
cross-confirmation". The result of family 32 ($\delta_4 = -3.71\times 10^{-3}$,
i.e.\ $\sim 48\sigma$ from $0$, but **statistically indistinguishable from
the non-$j=0$ quartic cluster** $\delta_4 \in [-3.83, -3.19]\times 10^{-3}$)
is therefore reported as B5 FALSIFIED rather than confirmed, *and the
falsification itself is qualified by the thin-cell caveat in
Remark~\ref{rem:thin-jzero-d4}*.

**The strongest conclusion currently licensed by the data is:** the
$j$-invariant signal that produces a $\rho_S = -0.691$ correlation at
$d=3$ is **not visible at $d=4$ in the WKB regime $N \le N_{\rm ref}\le 400$**.

To strengthen to a sharper falsification (or to convert to confirmation),
two concrete extensions are recommended:

1. Extend the quartic catalogue to $a_3 \in \{-3,\ldots,-1,1,2,3\}$ and
   wider $a_2$ window so that $\ge 2$ further irreducible $j=0$ quartics
   are captured (estimated $\sim 5$--$8$ candidates from a window of
   $\sim 500$ irreducible quartics).
2. Run a deep-WKB harvest at $d=4$ with $N_{\rm ref} \ge 1000$,
   $N_{\rm fit} \in [200, 800]$, $\mathrm{dps}\ge 5000$. Estimated cost
   $\sim 30$ min for the $j=0$ family alone, $\sim 30$ h for all 60.

## (iii) For non-squarefree quartics, what is the right invariant?

Within the Q1 catalogue, all 60 families are *irreducible*
(filtered upstream), hence automatically squarefree. So the
non-squarefree case does not arise in the current dataset.

**Fall-back branch (recorded for future use):** for a reducible quartic
$b(x) = b_1(x) \cdot b_2(x)$ where one factor has degree $\ge 2$,
$C_b$ has lower geometric genus and the j-invariant of the Jacobian is
not well-defined as the j-invariant of an elliptic curve. The right
invariant in such degenerate cells is plausibly the j-invariant of the
*irreducible factor* of largest degree (degree-3 factor yields
j-invariant of the cubic; degree-4 reducible to two irreducible quadratics
yields a pair of CM-j-invariants on a product Jacobian). This branch
should be implemented if the quartic enumeration is ever extended to
include reducible $b(x)$.

## (iv) Is "cross-degree universality" really warranted given only $d \in \{3,4\}$?

NO. The current data does NOT support a "cross-degree universality"
claim. In fact:

- $d=3$ Spearman: $\rho_S = -0.6906$, Bonf $p = 3.4\times 10^{-7}$
- $d=4$ Spearman: $\rho_S = +0.073$, Bonf $p = 1.0$ (NOT significant; sign opposite)
- Joint $d \in \{3,4\}$ Spearman: $\rho_S = -0.544$, Bonf $p = 9.5\times 10^{-9}$
  — but this is **dominated by the $d=3$ signal** (degree is a confounder
  because cubic populations have wider $\log|j|$ spread).

The honest framing in the v1.3 paragraph insert is therefore:
> "The $j$-invariant finer-split is degree-$3$ specific in the WKB
>  regime $N \le 400$; cross-degree universality is FALSIFIED
>  (with thin-cell caveat at $d=4$)."

The earlier draft "universal across the genus-1 cell $d \in \{3,4\}$"
phrasing is **withdrawn**.

## Sign-off

The session correctly:
- Identifies the unique $j=0$ quartic in the Q1 catalogue.
- Diagnoses the small-$N$ WKB sub-leading bias (drops $\delta_4$
  from $\sim -0.046$ at dps=80 to $\sim -0.0035$ at dps=2000, $N\in[50,250]$).
- Verifies the $d=3$ signal is reproducible ($\rho = -0.6906$ matches
  R1.1 to all reported digits).
- Concludes B5/B6 do NOT extend to $d=4$ at current precision
  (HALT triggered).
- Recommends two concrete follow-ups (catalogue extension; deep WKB
  at $d=4$).

The HALT condition fired correctly: "$j_{\rm zero}$ cell at $d=4$ is
non-empty BUT all delta values inconsistent with $0$ at A_stderr"
$\Rightarrow$ B5 FALSIFIED at $d=4$ $\Rightarrow$ halt for user review.
