#!/usr/bin/env python3
"""
area2_stieltjes_test.py
Test Stieltjes integral representation for V_quad via Sturm-Liouville weight.

Hypothesis: for PCF with a(n)=1, b(n) in Z[n], the limit K satisfies some
integral relation involving w(x) = -c(x)/b(x), the SL weight.
"""

from mpmath import mp, mpf, quad, log, pi, sqrt, pslq, nstr, inf, fabs

mp.dps = 60  # extra guard digits beyond 50

# ── helpers ──────────────────────────────────────────────────────────────────

def pcf_limit(a_func, b_func, n_terms=300, dps=60):
    """Compute PCF limit K = b(0) + a(1)/(b(1) + a(2)/(b(2) + ...)) via
    backward recurrence (Lentz-modified Wallis)."""
    old_dps = mp.dps
    mp.dps = dps + 20
    # backward recurrence from n_terms down
    val = mpf(0)
    for n in range(n_terms, 0, -1):
        val = a_func(n) / (b_func(n) + val)
    val = b_func(0) + val
    mp.dps = old_dps
    return +val  # round to current dps


def compute_c(b_poly_coeffs):
    """For the exact ODE  d/dx[b(x)y'] + c(x)y = 0  with a(n)=1,
    the coefficient c(x) satisfies the relation from the continuum limit.

    For a(n)=1, b(n) polynomial of degree d, the ODE is:
      b(x)y'' + b'(x)y' + c(x)y = 0
    which is exact (self-adjoint): d/dx[b(x)y'] + c(x)y = 0.

    From the three-term recurrence y_{n+1} + y_{n-1} = b(n) y_n  (with a(n)=1),
    the continuum limit gives c(x) = -(shift operator - 2 + shift^{-1})(y)/y
    For polynomial b of degree d, c(x) is determined by:
      b(x) + 1/b(x) approx b(x) for large x, and the correction terms.

    Actually for a(n)=1 the recurrence is:
      b(n) y_n = y_{n+1} + y_n  ... no wait.
    
    Standard PCF: K = b(0) + 1/(b(1) + 1/(b(2) + ...))
    Three-term: y_{n+1} = b(n) y_n - a(n) y_{n-1}  with a(n)=1
    So: y_{n+1} = b(n) y_n - y_{n-1}
    
    Continuum: y''(x) + (2 - b(x)) y(x) ≈ 0  ... not quite.
    
    Actually from the T8 area2 work:
    The exact ODE is d/dx[b(x)y'] + c(x)y = 0.
    For the three-term recurrence y_{n+1} - b(n)y_n + y_{n-1} = 0,
    the WKB / continuum limit gives the ODE with:
      c(x) = 2 - b(x)  ... for the leading order.
    
    But more precisely, from the self-adjoint form:
    b(x)y'' + b'(x)y' + c(x)y = 0
    
    The characteristic equation of the recurrence at large n is:
    λ^2 - b(n)λ + 1 = 0  (since a(n)=1)
    
    For the SL weight, we use w(x) = -c(x)/b(x).
    
    From the T8 work: c(x) comes from matching the discrete and continuous.
    For degree-2 b(x), c(x) was identified as -x^2 for V_quad case.
    
    Let me just use the specific c(x) values given in the problem.
    """
    pass  # We'll use explicit c(x) for each case


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 & 2: V_quad — b(n) = 3n^2 + n + 1,  c(x) = -x^2
# ═══════════════════════════════════════════════════════════════════════════════

print("=" * 72)
print("V_quad: b(n) = 3n^2 + n + 1, a(n) = 1")
print("  w(x) = -c(x)/b(x) = x^2 / (3x^2 + x + 1)")
print("=" * 72)

V_quad = pcf_limit(lambda n: mpf(1), lambda n: 3*n**2 + n + 1)
print(f"\nV_quad = {nstr(V_quad, 50)}")

# w(x) = x^2 / (3x^2 + x + 1)
def w_vquad(x):
    return x**2 / (3*x**2 + x + 1)

print("\n── STEP 1: Numerical integrals ──")

# I1 = ∫_0^∞ w(x) dx
I1 = quad(w_vquad, [0, inf])
print(f"I1 = ∫_0^∞ w(x) dx          = {nstr(I1, 50)}")

# I2 = ∫_0^∞ w(x)/x dx
def f2(x):
    return x / (3*x**2 + x + 1)  # w(x)/x = x/(3x^2+x+1)
I2 = quad(f2, [0, inf])
print(f"I2 = ∫_0^∞ w(x)/x dx        = {nstr(I2, 50)}")

# I3 = ∫_0^∞ w(x)/(1+x) dx
def f3(x):
    return x**2 / ((3*x**2 + x + 1) * (1 + x))
I3 = quad(f3, [0, inf])
print(f"I3 = ∫_0^∞ w(x)/(1+x) dx    = {nstr(I3, 50)}")

# I4 = ∫_0^∞ w(x)/(x(1+x)) dx
def f4(x):
    return x / ((3*x**2 + x + 1) * (1 + x))
I4 = quad(f4, [0, inf])
print(f"I4 = ∫_0^∞ w(x)/(x(1+x)) dx = {nstr(I4, 50)}")

# I5 = ∫_0^1 w(x) dx
I5 = quad(w_vquad, [0, 1])
print(f"I5 = ∫_0^1 w(x) dx          = {nstr(I5, 50)}")

# I6 = ∫_0^∞ x·w(x) dx
def f6(x):
    return x**3 / (3*x**2 + x + 1)
# I6 diverges (degree 3 / degree 2 → ∞), use cutoff
try:
    I6 = quad(f6, [0, inf], error=True)
    if isinstance(I6, tuple):
        I6 = I6[0]
    print(f"I6 = ∫_0^∞ x·w(x) dx        = {nstr(I6, 50)}")
except:
    I6 = None
    print(f"I6 = ∫_0^∞ x·w(x) dx        = DIVERGENT (x·w ~ x^3/(3x^2) = x/3 → ∞)")

# I7 = -∫_0^∞ log(x)·w(x) dx
def f7(x):
    return -log(x) * x**2 / (3*x**2 + x + 1)
I7 = quad(f7, [0, inf])
print(f"I7 = -∫_0^∞ log(x)·w(x) dx  = {nstr(I7, 50)}")

# Additional: Stieltjes transform S(s) = ∫ w(x)/(s+x) dx at various s
print("\n── Stieltjes transforms S(s) = ∫_0^∞ w(x)/(s+x) dx ──")
for s_val in [mpf(1), mpf(2), mpf('0.5'), V_quad, 1/V_quad]:
    def fs(x, s=s_val):
        return x**2 / ((3*x**2 + x + 1) * (s + x))
    Ss = quad(fs, [0, inf])
    print(f"  S({nstr(s_val, 8)}) = {nstr(Ss, 40)}")

print("\n── STEP 2: Comparison to V_quad ──")
V2 = V_quad**2
Vinv = 1/V_quad

integrals = {"I1": I1, "I2": I2, "I3": I3, "I4": I4, "I5": I5, "I7": I7}
if I6 is not None:
    integrals["I6"] = I6

print(f"\n{'Name':>6} {'Value':>50} {'|I-V|':>14} {'|I-1/V|':>14} {'|I-V^2|':>14}")
print("-" * 100)
for name, val in integrals.items():
    dV = fabs(val - V_quad)
    dVinv = fabs(val - Vinv)
    dV2 = fabs(val - V2)
    print(f"{name:>6} {nstr(val, 40):>50} {nstr(dV, 4):>14} {nstr(dVinv, 4):>14} {nstr(dV2, 4):>14}")

# Find closest
best_name, best_diff, best_target = None, mpf(1e30), ""
for name, val in integrals.items():
    for tname, tval in [("V", V_quad), ("1/V", Vinv), ("V^2", V2)]:
        d = fabs(val - tval)
        if d < best_diff:
            best_name, best_diff, best_target = name, d, tname
print(f"\n→ CLOSEST: {best_name} to {best_target}, diff = {nstr(best_diff, 6)}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: PSLQ
# ═══════════════════════════════════════════════════════════════════════════════

print("\n── STEP 3: PSLQ ──")
mp.dps = 55
pslq_vec = [mpf(1), I1, I2, I3, I4, I5]
if I6 is not None:
    pslq_vec.append(I6)
pslq_vec.extend([I7, V_quad])

labels = ["1", "I1", "I2", "I3", "I4", "I5"]
if I6 is not None:
    labels.append("I6")
labels.extend(["I7", "V_quad"])

print(f"PSLQ input vector ({len(pslq_vec)} elements): {labels}")
rel = pslq(pslq_vec)
if rel is not None:
    # Check residual
    residual = sum(c * v for c, v in zip(rel, pslq_vec))
    print(f"PSLQ relation found: {rel}")
    print(f"Labels:               {labels}")
    terms = []
    for c, l in zip(rel, labels):
        if c != 0:
            terms.append(f"({c})*{l}")
    print(f"Relation: {' + '.join(terms)} = 0")
    print(f"Residual: {nstr(residual, 6)}")
    if fabs(residual) < mpf('1e-20'):
        print("→ PSLQ: HIT (residual < 1e-20)")
    else:
        print(f"→ PSLQ: WEAK (residual = {nstr(residual, 6)})")
else:
    print("→ PSLQ: NULL (no relation found)")

# Also try smaller PSLQ subsets
print("\n── PSLQ sub-tests ──")
for subset_labels, subset_vals in [
    (["1", "I1", "V_quad"], [mpf(1), I1, V_quad]),
    (["1", "I2", "V_quad"], [mpf(1), I2, V_quad]),
    (["1", "I3", "V_quad"], [mpf(1), I3, V_quad]),
    (["1", "I7", "V_quad"], [mpf(1), I7, V_quad]),
    (["1", "I1", "I2", "V_quad"], [mpf(1), I1, I2, V_quad]),
    (["1", "I2", "I3", "V_quad"], [mpf(1), I2, I3, V_quad]),
    (["1", "I3", "I4", "V_quad"], [mpf(1), I3, I4, V_quad]),
]:
    r = pslq(subset_vals)
    if r is not None:
        res = sum(c * v for c, v in zip(r, subset_vals))
        if fabs(res) < mpf('1e-20'):
            terms = [f"({c})*{l}" for c, l in zip(r, subset_labels) if c != 0]
            print(f"  HIT: {' + '.join(terms)} = 0  (res={nstr(res, 4)})")
        else:
            print(f"  WEAK {subset_labels}: {r}  (res={nstr(res, 4)})")
    else:
        print(f"  null {subset_labels}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: b(n) = n + 1  (simplest case)
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("STEP 4: b(n) = n + 1, a(n) = 1")
print("=" * 72)

mp.dps = 60
K_np1 = pcf_limit(lambda n: mpf(1), lambda n: n + 1)
print(f"PCF limit K = {nstr(K_np1, 50)}")

# Check known constants
from mpmath import e, euler
print(f"  1/e       = {nstr(1/e, 50)}")
print(f"  e - 1     = {nstr(e - 1, 50)}")
print(f"  e         = {nstr(e, 50)}")
print(f"  log(2)    = {nstr(log(2), 50)}")
print(f"  |K - 1|   = {nstr(fabs(K_np1 - 1), 6)}")
print(f"  |K - e|   = {nstr(fabs(K_np1 - e), 6)}")
print(f"  |K - 1/e| = {nstr(fabs(K_np1 - 1/e), 6)}")

# For b(n) = n+1: three-term is y_{n+1} = (n+1) y_n - y_{n-1}
# c(x) from the self-adjoint ODE: d/dx[(x+1)y'] + c(x)y = 0
# The characteristic eq: λ^2 - (x+1)λ + 1 = 0
# For the SL weight, c(x) needs to be determined.
# From the recurrence continuum limit: c(x) ≈ 2 - b(x) = 1 - x
# So w(x) = -c(x)/b(x) = (x-1)/(x+1)

# But w must be non-negative on the support, so support is [1, ∞)
print("\nFor c(x) = 1 - x (continuum approx):")
print("  w(x) = (x-1)/(x+1), support [1, ∞)")

def w_np1(x):
    return (x - 1) / (x + 1)

J1 = quad(lambda x: w_np1(x) / (1 + x), [1, inf])
print(f"  ∫_1^∞ (x-1)/(x+1)^2 dx = {nstr(J1, 40)}")

J2 = quad(lambda x: w_np1(x), [1, inf])
print(f"  ∫_1^∞ (x-1)/(x+1) dx   = DIVERGENT (→ ∞)")

# Actually c(x) might not be 1-x. Let's just test specific integrals.
# For b(x) = x+1, the simplest w is 1/(x+1)
print("\nAlternative: w(x) = 1/(x+1)")
J3 = quad(lambda x: 1 / (1 + x)**2, [0, inf])
print(f"  ∫_0^∞ 1/(1+x)^2 dx = {nstr(J3, 40)} (should be 1)")
print(f"  |K - 1| = {nstr(fabs(K_np1 - 1), 6)}")

# Try to identify K via PSLQ with common constants
pslq_check = pslq([mpf(1), K_np1, 1/e, e, log(2), pi])
if pslq_check:
    print(f"  PSLQ[1, K, 1/e, e, ln2, π] = {pslq_check}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: b(n) = 2n + 1  (Catalan-type)
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("STEP 5: b(n) = 2n + 1, a(n) = 1")
print("=" * 72)

K_2np1 = pcf_limit(lambda n: mpf(1), lambda n: 2*n + 1)
print(f"PCF limit K = {nstr(K_2np1, 50)}")

# Check known constants
print(f"  1/tanh(1) = {nstr(1/mp.tanh(1), 50)}")
print(f"  coth(1)   = {nstr(mp.coth(1), 50)}")
print(f"  tanh(1)   = {nstr(mp.tanh(1), 50)}")
print(f"  √e        = {nstr(sqrt(e), 50)}")
print(f"  |K-coth1| = {nstr(fabs(K_2np1 - mp.coth(1)), 6)}")

# For b(n) = 2n+1: this is the well-known CF for coth(1)
# Three-term: y_{n+1} = (2n+1)y_n - y_n-1
# c(x) from continuum: c(x) = 2 - (2x+1) = 1 - 2x
# w(x) = (2x-1)/(2x+1), support [1/2, ∞)

# Also: the standard formula is coth(1) = 1 + 1/(3 + 1/(5 + ...))
# which is exactly b(n) = 2n+1, a(n) = 1.

print(f"\n  K is coth(1)? diff = {nstr(fabs(K_2np1 - mp.coth(1)), 6)}")

# Test w(x) = 1/(2x+1)
def w_2np1(x):
    return 1 / (2*x + 1)

print("\nw(x) = 1/(2x+1):")
L1 = quad(w_2np1, [0, inf])
print(f"  ∫_0^∞ 1/(2x+1) dx        = DIVERGENT")

L2 = quad(lambda x: 1 / ((2*x+1)*(1+x)), [0, inf])
print(f"  ∫_0^∞ 1/((2x+1)(1+x)) dx = {nstr(L2, 40)}")
print(f"  2·log(2)                  = {nstr(2*log(2), 40)}")

# c(x) = 1 - 2x, w(x) = (2x-1)/(2x+1) on [1/2, ∞)
print("\nw(x) = (2x-1)/(2x+1) on [1/2, ∞):")
def w_cat(x):
    return (2*x - 1) / (2*x + 1)

L3 = quad(lambda x: w_cat(x) / (1 + x), [mpf(0.5), inf])
print(f"  ∫_{'{1/2}'}^∞ w(x)/(1+x) dx   = {nstr(L3, 40)}")
L4 = quad(lambda x: w_cat(x) / (x*(1+x)), [mpf(0.5), inf])
print(f"  ∫_{'{1/2}'}^∞ w(x)/(x(1+x)) dx = {nstr(L4, 40)}")

coth1 = mp.coth(1)
print(f"\n  |L2 - coth1| = {nstr(fabs(L2 - coth1), 6)}")
print(f"  |L3 - coth1| = {nstr(fabs(L3 - coth1), 6)}")
print(f"  |L4 - coth1| = {nstr(fabs(L4 - coth1), 6)}")

# PSLQ for coth(1) case
pslq_cat = pslq([mpf(1), L2, L3, L4, coth1])
if pslq_cat:
    print(f"  PSLQ[1, L2, L3, L4, coth1] = {pslq_cat}")


# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

print(f"\nCLOSEST INTEGRAL to V_quad: {best_name} → {best_target}, "
      f"|diff| = {nstr(best_diff, 6)}")

print(f"\nPCF limits:")
print(f"  b(n)=3n²+n+1 → K = {nstr(V_quad, 30)}")
print(f"  b(n)=n+1     → K = {nstr(K_np1, 30)}")
print(f"  b(n)=2n+1    → K = {nstr(K_2np1, 30)} ≈ coth(1)")

print(f"\nPATTERN: (see integral comparisons above)")
print(f"STIELTJES HYPOTHESIS: (determined by output)")
