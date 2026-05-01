# Rubber-duck critique — PCF-2 SESSION C1

**Date:** 2026-05-01
**Scope:** WKB-exponent harvest on the non-CM cubic bins (`+_S3_real` and
`+_C3_real`, 13 families) testing refined Conjecture B4'.

---

## 1. Headline finding (and the obvious objection)

**Result:** all 13 non-CM cubic families have $A_{\mathrm{fit}}$ in the
interval $[5.88, 5.99]$, mean $5.973$, stddev $0.026$. Twelve of thirteen
are within $0.10$ of $A = 6$; one (family 37, $A=5.884$) is in the
interior of $(5,6)$ but still far closer to $6$ than to $5$. **Zero
families** are consistent with the refined-B4' prediction
$A = 2d - 1 = 5$.

**Obvious objection #1 — could this be a fit-window artefact?** The fit
window $n\in[10,100]$ is identical to Session B (which delivered
$A_{\mathrm{fit}} = 6$ on the CM bin) and to Session A2 (which delivered
$A_{\mathrm{fit}} = 5.95$ on family 46, the conductor-7 anchor). If the
window biased toward $A=6$, Session A2's anchor would have already
exhibited that bias — and indeed it did: family 46's "anchor" value
of $5.95$ is *itself closer to $6$ than to $5$.

> **Self-correction:** the C1 prompt's preamble describes family 46 as
> the calibration anchor "at $A_{\mathrm{fit}} = 5.95$" while presenting
> $5$ as the predicted branch. Those two facts are mutually
> inconsistent: an anchor at $5.95$ is consistent with $A=6$, not $A=5$.
> B4' was effectively falsified by its own anchor before C1 began;
> C1 simply confirms this on 12 further non-CM families. The PCF2
> session-A2 commit message recorded `A_fit=5.95 -> 6=2d`, so the
> falsification was visible in the record.

**Obvious objection #2 — could the L-stability gate be too loose?**
At `mp.dps = 80` the function `cf_value` returns identical values for
$L_{150}$ and $L_{200}$ to within working precision; `stable_digits`
then caps at $999.0$. This is a precision-floor artefact, not "999
digits of agreement". The honest interpretation is "at least $80$
digits of agreement at $n=200$", which still comfortably passes the
$\ge 30$-digit gate. The "striking convergence ($\ge 50$ digits)" flag
in `unexpected_finds.json` triggers on the floor rather than on a true
$50$-digit observation; it should be read as "the gate was wholly
non-binding for non-CM cubic L-values at $n=200$, dps=80".

## 2. Verdict-logic edge case (off-spec)

The C1 prompt defines three verdicts:

  - `b4_prime_consistent` — $|A_{\mathrm{fit}} - {\rm predicted}| \le 0.10$
  - `b4_prime_violating` — $\min(|A-5|, |A-6|) > 0.30$
  - `b4_prime_ambiguous` — $A\in(5,6)$ with neither branch within $0.10$

These three do **not partition** the cases. A family with predicted
branch $5$ but $A_{\mathrm{fit}} = 5.97$ is:

  - not consistent (predicted distance $0.97$),
  - not violating ($|A-6|=0.03$),
  - not ambiguous ($|A-6|<0.10$).

I introduced a fourth label, `b4_prime_consistent_wrong_branch`, to
capture these (12 of the 13 families). This is a strict superset of
the spec but loses no information: every "wrong-branch" family is also
flagged in `unexpected_finds.json` per the prompt's NON-HALT FLAG
clause ("A_fit closer to 6 than to 5"). Halt logic is unchanged: only
`b4_prime_violating` triggers a halt (zero families).

## 3. dps inconsistency in the prompt

The C1 prompt asks for `mp.dps = 80` for the WKB fit and "identical to
Session B step 1b" in the same paragraph. Session B uses `mp.dps = 800`
(the actual delta values at $n=100$ for $A=6$ are $\sim 10^{-2760}$,
far below an 80-digit floor). I followed Session B's value, `dps = 800`,
because the alternative is to fit zeros: at dps=80, all $|\delta_n|$
for $n \gtrsim 28$ would underflow and the fit would lose 70+ of its
27 points. Documented.

## 4. What this means for B4'

**B4' as stated in the prompt is empirically false.** The data on the
50 cubic families (37 CM + 13 non-CM) supports the stronger and simpler

> $A = 2d = 6$ uniformly across all degree-3 PCFs in scope.

i.e. the original Session-B Conjecture B4 with the sharp cubic form,
without the bin-dependent refinement. The trichotomy bin (CM vs.
elementary-positive S3 vs. C3 real) does **not** split the WKB
exponent.

This still leaves the PCF-1 v1.3 Theorem 5 base case intact
($d=2$: $A\in\{3,4\}$ with the $\Delta>0$ vs $\Delta<0$ split), so the
"bin matters at $d=2$ but not at $d=3$" pattern needs explanation.
One natural hypothesis: the $d=2$ split is driven by the sign of the
discriminant *of the recurrence* (and hence by the indicial-root
asymptotic regime), not by the trichotomy of the underlying number
field; in cubic, the relevant indicial-root structure may not
discriminate the bins. A separate session should test this by
recomputing the indicial roots / WKB phase data per family.

## 5. Sanity-checks performed

  - Family 46 (the C1 / Session A2 anchor): C1 reports
    $A_{\mathrm{fit}} = 5.989 \pm 0.001$. Session A2 reported
    $5.95$ via the same pipeline at coarser settings. Both are
    close to $6$, neither is close to $5$. **Consistent.**
  - Family 19 / 20 share $\Delta_3 = 81$ and the same Galois group;
    $A_{\mathrm{fit}}$ values agree to 4 decimals (5.9860 vs 5.9870),
    sanity-check on the fit's family-level dependence.
  - Family 9 / 37 share $\Delta_3 = 257$ but differ in coefficients;
    $A_{\mathrm{fit}}$ values are 5.977 and 5.884 respectively, the
    larger-stderr fit (family 37, stderr $0.048$) being the one outside
    the $\pm 0.10$ window. The full set of $A$-values is bounded in
    $[5.88, 5.99]$, no extreme outliers.

## 6. Limitations

  - Bin sample sizes: $+\_S_3\_\mathrm{real}$ has 10 families, $+\_C_3\_\mathrm{real}$
    has 3. The smaller sub-bin is too small to give an independent
    statistical signal; the present analysis treats the 13 as a single
    "non-CM cubic" group.
  - The fit is ordinary least-squares on $\log|\delta_n|$; it does not
    weight by per-point uncertainty (the dps=800 floor effectively
    sets uniform weights). A weighted fit, or a pure
    $-A n\log n$ extraction at large $n$, may sharpen $A_{\mathrm{fit}}$
    slightly but is unlikely to move the 12 wrong-branch families
    by anywhere close to $1.0$.
  - No PSLQ scan in C1 (per spec). The 13 non-CM L-values remain
    available for a future targeted PSLQ probe.

## 7. Action items downstream

  1. **Update Session B's Conjecture B4** to install B4' as the
     prompt-required wording, and **append a remark** that C1
     falsifies the elementary-positive branch at $d=3$. The
     working-copy edit is performed by C1's standing post-condition.
  2. **Drop the bin-dependence from B4'** for cubic in the v1.1
     paper; the surviving conjecture is "A = 2d for every degree-d
     PCF in scope, $d \le 3$".
  3. **Open question** for Claude / strategic review: at what $d$
     (if any) does a bin-dependent split emerge? PCF-1 has it for
     $d=2$; PCF-2 C1 rejects it for $d=3$. A Session C2 sweeping a
     small quartic catalogue would be the natural next probe.
