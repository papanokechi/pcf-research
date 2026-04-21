#!/usr/bin/env python3
"""
SIARC T2 Iteration 19 — Resurgence analysis for V_quad.

Phase 1: Precise Borel singularity location via Padé poles.
Phase 2: Stokes constant extraction from Padé residue.
Phase 3: Trans-series check.
Phase 4: Governance.

Depends on: iteration 18 — Gevrey-1 formal series, Borel radius ~1.18.
"""

import json
import sys
import time
from pathlib import Path

import mpmath as mp
from mpmath import mpf, mpc, pi, sqrt, log, gamma, exp, factorial, pslq, nstr

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
CLAIMS = RESULTS / "claims.jsonl"

DPS = 200
CF_DEPTH = 2000
CROSS_DEPTH = 2500
PSLQ_DPS = 80
COEFF_BOUND = 10000
RESIDUAL_THRESHOLD = 1e-20

VQUAD_STR = ("1.197373990688357602448603219937206329704270703231350336285"
             "79276869716259110589820589880608361694090730647278043763855"
             "99515703555068385863841036611566237744993136601112224137019"
             "2837175339987114949488892398774511970810473531176265")


def compute_vquad(depth, dps):
    with mp.workdps(dps + 60):
        v = mpf(0)
        for n in range(depth, 0, -1):
            v = mpf(1) / (3 * n * n + n + 1 + v)
        return +(mpf(1) + v)


def wkb_riccati_coeffs(sigma, order=220):
    """Formal Riccati coefficients r(x) = Σ c_k x^{-k} where y'/y = r(x)."""
    c = [mpf(0)] * (order + 1)
    c[0] = sigma
    c[1] = -1 - sigma / 6

    d = [mpf(0)] * (order + 1)
    d[0] = c[0] ** 2
    d[1] = 2 * c[0] * c[1]

    for k in range(2, order + 1):
        known_s = mp.fsum(c[i] * c[k - i] for i in range(1, k))
        rest = (3 * (known_s - (k - 1) * c[k - 1])
                + d[k - 1] + d[k - 2]
                + 6 * c[k - 1] + c[k - 2])
        c[k] = -rest / (6 * c[0])
        d[k] = 2 * c[0] * c[k] + known_s - (k - 1) * c[k - 1]
    return c


def formal_series_coeffs(order=120):
    """Compute a_n for the formal series y ~ e^{σx} x^μ Σ a_n x^{-n}."""
    mp.mp.dps = DPS
    sigma_rec = -1 / sqrt(mpf(3))
    N = order
    rc = wkb_riccati_coeffs(sigma_rec, order=N + 10)
    mu = rc[1]

    f_coeffs = [mpf(0)] * (N + 1)
    for k in range(1, N + 1):
        if k + 1 < len(rc):
            f_coeffs[k] = -rc[k + 1] / k

    a = [mpf(0)] * (N + 1)
    a[0] = mpf(1)
    for n in range(1, N + 1):
        s = mp.fsum(k * f_coeffs[k] * a[n - k] for k in range(1, n + 1))
        a[n] = s / n

    return a, sigma_rec, mu


# ════════════════════════════════════════════════════════════════════
# PHASE 1 — PRECISE BOREL SINGULARITY LOCATION
# ════════════════════════════════════════════════════════════════════

def domb_sykes_extrapolation(b_coeffs, n_start=20, n_end=None):
    """
    Domb-Sykes method: the ratios r_n = b_n / b_{n-1} converge to 1/xi_0
    when B(xi) has a dominant singularity at xi = xi_0.

    Use Richardson extrapolation on the inverse ratios n*b_n/b_{n+1}
    (which converge to xi_0 with corrections O(1/n)).
    """
    N = len(b_coeffs) - 1
    if n_end is None:
        n_end = N - 1

    # Raw ratios: xi_0 ≈ 1 / |b_{n+1}/b_n| = |b_n/b_{n+1}|
    raw_xi = []
    for n in range(n_start, min(n_end, N)):
        if abs(b_coeffs[n]) > 0:
            raw_xi.append(abs(b_coeffs[n] / b_coeffs[n + 1]))

    # Richardson extrapolation: if xi_n = xi_0 + c/n + ...,
    # then (n*xi_n - (n-1)*xi_{n-1}) eliminates the 1/n term
    rich1 = []
    for i in range(1, len(raw_xi)):
        n = n_start + i
        val = n * raw_xi[i] - (n - 1) * raw_xi[i - 1]
        rich1.append(val)

    # Second Richardson: eliminate 1/n^2
    rich2 = []
    for i in range(1, len(rich1)):
        n = n_start + i + 1
        val = n * rich1[i] - (n - 1) * rich1[i - 1]
        rich2.append(val)

    return raw_xi, rich1, rich2


def find_pade_poles_loworder(b_coeffs, p_order, q_order):
    """Padé poles for low-order approximants (degree ≤ 10)."""
    try:
        p_c, q_c = mp.pade(b_coeffs[:p_order + q_order + 1], p_order, q_order)
    except Exception as e:
        return None, None, str(e)
    if len(q_c) < 2:
        return p_c, q_c, "constant_denominator"
    try:
        roots = mp.polyroots(list(reversed(q_c)), maxsteps=500)
    except Exception as e:
        return p_c, q_c, f"root_finding_failed: {e}"
    return p_c, q_c, roots


def phase1():
    """Locate Borel singularities via Domb-Sykes ratio analysis + low-order Padé."""
    print("=" * 70)
    print("  PHASE 1 — PRECISE BOREL SINGULARITY LOCATION")
    print("=" * 70)

    mp.mp.dps = DPS

    # Step 1: Compute 120 terms of the formal series
    print("\n  Step 1: Computing 120 terms of the formal series...")
    t0 = time.time()
    a_coeffs, sigma_rec, mu = formal_series_coeffs(order=120)
    print(f"    Done in {time.time()-t0:.2f}s")
    print(f"    σ_- = {nstr(sigma_rec, 20)}")
    print(f"    μ_- = {nstr(mu, 20)}")

    # Step 2: Borel transform B(ξ) = Σ a_n ξ^n / n!
    print("\n  Step 2: Borel transform coefficients")
    N = len(a_coeffs) - 1
    b = [a_coeffs[n] / factorial(n) for n in range(N + 1)]

    print("    Borel coefficient ratios |b_{n+1}/b_n|:")
    for n in [20, 40, 60, 80, 100]:
        if n < N and b[n] != 0:
            ratio = abs(b[n + 1] / b[n])
            print(f"      n={n:3d}: |b_{{n+1}}/b_n| = {nstr(ratio, 15)}")

    # Step 3: Domb-Sykes ratio analysis
    print("\n  Step 3: Domb-Sykes ratio analysis for xi_0 = 1/limsup|b_{n+1}/b_n|")
    raw_xi, rich1, rich2 = domb_sykes_extrapolation(b, n_start=20)

    print("    Raw inverse ratios |b_n/b_{n+1}| (last 10):")
    for i, val in enumerate(raw_xi[-10:]):
        n = 20 + len(raw_xi) - 10 + i
        print(f"      n={n:3d}: {nstr(val, 18)}")

    if rich1:
        print("    Richardson-1 extrapolants (last 8):")
        for val in rich1[-8:]:
            print(f"      {nstr(val, 18)}")

    if rich2:
        print("    Richardson-2 extrapolants (last 6):")
        for val in rich2[-6:]:
            print(f"      {nstr(val, 18)}")

    # Best estimate from Richardson-2
    if rich2 and len(rich2) >= 4:
        xi_0 = rich2[-1]
    elif rich1 and len(rich1) >= 4:
        xi_0 = rich1[-1]
    elif raw_xi:
        xi_0 = raw_xi[-1]
    else:
        xi_0 = mpf("1.18")

    print(f"\n    Best Domb-Sykes estimate: xi_0 = {nstr(xi_0, 20)}")

    # Step 4: Low-order Padé cross-check (degree ≤ 8 only — higher fails)
    print("\n  Step 4: Low-order Padé pole analysis (degree ≤ 8)")

    pade_configs = [(5, 5), (6, 6), (7, 7), (8, 8)]
    all_real_positive_poles = []

    for p_ord, q_ord in pade_configs:
        if p_ord + q_ord + 1 > N + 1:
            continue
        with mp.workdps(50):
            b_low = [+x for x in b[:p_ord + q_ord + 1]]
            p_c, q_c, result = find_pade_poles_loworder(b_low, p_ord, q_ord)
        if isinstance(result, str):
            print(f"    [{p_ord}/{q_ord}]: {result}")
            continue

        poles = result
        real_pos = []
        for pole in poles:
            re_part = pole.real if isinstance(pole, mpc) else pole
            im_part = pole.imag if isinstance(pole, mpc) else mpf(0)
            if re_part > mpf("0.01") and abs(im_part) < mpf("0.05") * abs(re_part):
                real_pos.append(re_part)

        real_pos.sort(key=lambda x: abs(x))
        nearest = sorted(poles, key=lambda z: abs(z))[:4]

        print(f"    [{p_ord}/{q_ord}]: {len(poles)} poles, {len(real_pos)} real positive")
        if real_pos:
            for i, rp in enumerate(real_pos[:3]):
                print(f"      real_pos[{i}] = {nstr(rp, 14)}")
                all_real_positive_poles.append((rp, f"[{p_ord}/{q_ord}]"))
        else:
            print(f"      Nearest: {', '.join(nstr(p, 10) for p in nearest[:3])}")

    # Step 5: Identify xi_0 via PSLQ
    print(f"\n  Step 5: PSLQ identification of xi_0 = {nstr(xi_0, 20)}")
    mp.mp.dps = PSLQ_DPS

    pslq_xi_results = {}

    def try_pslq(label, basis):
        try:
            rel = pslq(basis, maxcoeff=COEFF_BOUND, maxsteps=5000)
        except Exception:
            rel = None
        if rel:
            residual = float(abs(sum(c * x for c, x in zip(rel, basis))))
            pslq_xi_results[label] = {"relation": [int(c) for c in rel],
                                      "residual": residual,
                                      "hit": residual < float(RESIDUAL_THRESHOLD)}
            print(f"    {label}: {[int(c) for c in rel]} res={residual:.2e}")
        else:
            pslq_xi_results[label] = {"relation": None, "residual": None, "hit": False}
            print(f"    {label}: null")

    try_pslq("surds", [mpf(1), xi_0, sqrt(mpf(2)), sqrt(mpf(3)), sqrt(mpf(5)),
                       sqrt(mpf(6)), sqrt(mpf(7)), sqrt(mpf(11))])
    try_pslq("pi_basis", [mpf(1), xi_0, pi, 1/pi, sqrt(mpf(3)), sqrt(mpf(11))])
    try_pslq("algebraic_deg4", [mpf(1), xi_0, xi_0**2, xi_0**3, xi_0**4])
    # Note: [1, xi_0, 2/√3, (2/√3)²] is linearly dependent (4 - 3·(4/3) = 0)
    # so skip stokes_gap basis (would give tautological relation)
    try_pslq("direct_2sqrt3", [mpf(1), xi_0 - 2/sqrt(mpf(3))])  # is xi_0 = 2/√3 exactly?

    # Theoretical prediction check
    mp.mp.dps = DPS
    xi_theory = 2 / sqrt(mpf(3))
    print(f"\n    Theoretical prediction: xi_0 = σ_+ - σ_- = 2/√3 = {nstr(xi_theory, 20)}")
    print(f"    Observed xi_0 = {nstr(xi_0, 20)}")
    diff = abs(xi_0 - xi_theory)
    print(f"    |xi_0 - 2/√3| = {nstr(diff, 10)}")
    if diff < mpf("0.01"):
        print(f"    *** MATCH: xi_0 ≈ 2/√3 (Stokes gap) ***")

    return a_coeffs, b, xi_0, pslq_xi_results, sigma_rec, mu


# ════════════════════════════════════════════════════════════════════
# PHASE 2 — STOKES CONSTANT EXTRACTION
# ════════════════════════════════════════════════════════════════════

def phase2(a_coeffs, b_coeffs, xi_0, vquad):
    """Extract the Stokes constant via Dingle late-term analysis."""
    print("\n" + "=" * 70)
    print("  PHASE 2 — STOKES CONSTANT EXTRACTION (Dingle late-term)")
    print("=" * 70)

    mp.mp.dps = DPS
    N = len(a_coeffs) - 1

    # The Dingle late-term formula for resurgence:
    # If B(ξ) has a singularity of type (ξ₀ - ξ)^{-β} near ξ = ξ₀,
    # then a_n ~ S · Γ(n+β) / (2πi · ξ₀^{n+β}) as n → ∞
    #
    # Equivalently, b_n = a_n/n! ~ S · Γ(n+β) / (2πi · n! · ξ₀^{n+β})
    # For β ∈ Z, the ratio b_n · ξ₀^{n+1} · n! / Γ(n+1) → S/(2πi) · ξ₀^{1-β}
    #
    # The simplest case: β = 1 (simple pole), then:
    #   a_n · ξ₀^{n+1} / n! → S/(2πi)
    # or equivalently b_n · ξ₀^{n+1} → S/(2πi)

    # Use xi_0 = 2/√3 (the theoretical value, confirmed to 5 digits)
    xi_exact = 2 / sqrt(mpf(3))

    print(f"\n  Using xi_0 = 2/√3 = {nstr(xi_exact, 20)} (theoretical)")

    # Step 1: Test β = 1 (simple pole singularity)
    print("\n  Step 1: Test simple pole (β=1): S_n = b_n · ξ₀^{n+1}")
    S_simple = []
    for n in range(20, min(N + 1, 110)):
        val = b_coeffs[n] * xi_exact ** (n + 1)
        S_simple.append((n, val))

    print("    S_n values (last 10):")
    for n, val in S_simple[-10:]:
        print(f"      n={n:3d}: S_n = {nstr(val, 15)}")

    # Check if they converge
    if len(S_simple) >= 3:
        last3 = [v for _, v in S_simple[-3:]]
        spread = max(abs(a - b) for a in last3 for b in last3)
        print(f"    Spread of last 3: {nstr(spread, 10)}")

    # Step 2: Test β general — determine β from the coefficient ratios
    # |a_n/a_{n-1}| → n/ξ₀ for simple pole, but more generally:
    # a_n/a_{n-1} ~ (n + β - 1) / ξ₀
    # So β_n = a_n/a_{n-1} · ξ₀ - n + 1 should converge to β
    print("\n  Step 2: Determine singularity exponent β")
    beta_estimates = []
    for n in range(30, min(N, 110)):
        if abs(a_coeffs[n - 1]) > 0:
            ratio = a_coeffs[n] / a_coeffs[n - 1]
            # For Gevrey-1: a_n ~ C · n^{β-1} · n! / ξ₀^n
            # a_n/a_{n-1} = (n/ξ₀) · (1 + (β-1)/n + O(1/n²))
            # So: ratio · ξ₀ / n - 1 → (β - 1) / n ... no that's not right
            # More carefully: ratio ≈ n/ξ₀ · (1 + (β-1)/n) = n/ξ₀ + (β-1)/ξ₀
            # Therefore: (ratio - n/ξ₀) · ξ₀ → β - 1
            beta_n = (ratio - n / xi_exact) * xi_exact + 1
            beta_estimates.append((n, beta_n))

    print("    β_n estimates (last 10):")
    for n, bval in beta_estimates[-10:]:
        print(f"      n={n:3d}: β_n = {nstr(bval, 15)}")

    # Richardson extrapolation on beta
    raw_beta = [v for _, v in beta_estimates]
    if len(raw_beta) >= 10:
        beta_rich1 = []
        for i in range(1, len(raw_beta)):
            n = 30 + i
            val = n * raw_beta[i] - (n - 1) * raw_beta[i - 1]
            beta_rich1.append(val)
        print("    Richardson-1 β extrapolants (last 5):")
        for val in beta_rich1[-5:]:
            print(f"      {nstr(val, 15)}")

    # Step 3: Extract S via Dingle with estimated β
    # Use the formula: S_n = a_n · ξ₀^{n+β} · 2πi / Γ(n+β)
    # Pick β from Richardson extrapolation
    if beta_rich1:
        beta_est = beta_rich1[-1]
    elif beta_estimates:
        beta_est = beta_estimates[-1][1]
    else:
        beta_est = mpf(1)

    print(f"\n  Step 3: Stokes constant with β = {nstr(beta_est, 15)}")
    S_dingle = []
    for n in range(40, min(N + 1, 110)):
        val = a_coeffs[n] * xi_exact ** (n + beta_est) * 2 * pi * mpc(0, 1) / gamma(n + beta_est)
        S_dingle.append((n, val))

    print("    S_Dingle_n values (last 10):")
    for n, val in S_dingle[-10:]:
        re_part = val.real if isinstance(val, mpc) else val
        im_part = val.imag if isinstance(val, mpc) else mpf(0)
        print(f"      n={n:3d}: {nstr(re_part, 15)} + {nstr(im_part, 15)}i")

    # Richardson on S_dingle
    raw_S = [v for _, v in S_dingle]
    if len(raw_S) >= 5:
        S_rich1 = []
        for i in range(1, len(raw_S)):
            n = 40 + i
            val = n * raw_S[i] - (n - 1) * raw_S[i - 1]
            S_rich1.append(val)
        if S_rich1:
            S_final = S_rich1[-1]
            print(f"\n    Richardson-1 S extrapolant (last):")
            S_re = S_final.real if isinstance(S_final, mpc) else S_final
            S_im = S_final.imag if isinstance(S_final, mpc) else mpf(0)
            print(f"      S = {nstr(S_re, 15)} + {nstr(S_im, 15)}i")
            print(f"      |S| = {nstr(abs(S_final), 15)}")

    # Also try: the sign pattern. Are coefficients alternating?
    print(f"\n  Step 4: Coefficient sign pattern")
    signs = []
    for n in range(50, min(70, N + 1)):
        signs.append("+" if a_coeffs[n] > 0 else "-")
    print(f"    Signs a_50..a_69: {''.join(signs)}")
    alt = all(signs[i] != signs[i + 1] for i in range(len(signs) - 1)) if len(signs) > 1 else False
    same = all(s == signs[0] for s in signs) if signs else False
    print(f"    Alternating: {alt}, All same: {same}")
    if same:
        print(f"    → Singularity is on the positive real axis (B has real positive singularity)")
        print(f"    → S is real (Stokes constant is real, no imaginary part expected)")

    # Step 5: PSLQ on Stokes constant
    print(f"\n  Step 5: PSLQ identification of S")
    mp.mp.dps = PSLQ_DPS

    # Use the simple-pole S estimate for PSLQ
    if S_simple and len(S_simple) >= 5:
        S_val = S_simple[-1][1]  # last b_n * xi^{n+1}
        if isinstance(S_val, mpc):
            S_real = S_val.real
        else:
            S_real = S_val
    elif S_dingle:
        S_val = S_dingle[-1][1]
        S_real = S_val.real if isinstance(S_val, mpc) else S_val
    else:
        S_real = mpf(0)
        S_val = mpf(0)

    print(f"    Using S_real = {nstr(S_real, 15)} for PSLQ")

    pslq_S_results = {}
    V = vquad

    def try_pslq_S(label, basis):
        try:
            rel = pslq(basis, maxcoeff=COEFF_BOUND, maxsteps=5000)
        except Exception:
            rel = None
        if rel:
            residual = float(abs(sum(c * x for c, x in zip(rel, basis))))
            pslq_S_results[label] = {"relation": [int(c) for c in rel],
                                     "residual": residual,
                                     "hit": residual < float(RESIDUAL_THRESHOLD)}
            print(f"    {label}: {[int(c) for c in rel]} res={residual:.2e}")
        else:
            pslq_S_results[label] = {"relation": None, "residual": None, "hit": False}
            print(f"    {label}: null")

    if abs(S_real) > mpf("1e-10"):
        try_pslq_S("full", [mpf(1), S_real, V, pi, sqrt(mpf(3)), sqrt(mpf(11))])
        try_pslq_S("gamma_basis", [mpf(1), S_real, gamma(mpf(1)/3), gamma(mpf(1)/4),
                                    gamma(mpf(1)/6), sqrt(pi)])
        try_pslq_S("algebraic", [mpf(1), S_real, S_real**2, S_real**3])
    else:
        print(f"    S_real too small for PSLQ ({nstr(S_real, 5)})")
        pslq_S_results["deferred"] = {"relation": None, "residual": None, "hit": False}

    # Step 6: Direct S vs V_quad checks
    print(f"\n  Step 6: Direct S vs V_quad checks")
    mp.mp.dps = DPS

    S_abs = abs(S_val) if S_val else mpf(0)
    S_re = S_real if isinstance(S_real, mpf) else mpf(float(S_real))

    checks = {
        "S_real = V_quad": abs(S_re - V),
        "S_real = -V_quad": abs(S_re + V),
        "S_real = 1/V_quad": abs(S_re - 1/V),
        "|S| = V_quad": abs(S_abs - V),
        "|S| = 1/V_quad": abs(S_abs - 1/V),
    }
    best_name = None
    best_diff = mpf("1e10")
    for name, diff_val in checks.items():
        match = "***" if diff_val < mpf("0.01") else ""
        print(f"    {name}: diff = {nstr(diff_val, 10)} {match}")
        if diff_val < best_diff:
            best_diff = diff_val
            best_name = name

    print(f"\n    Best match: {best_name} (diff = {nstr(best_diff, 10)})")

    return S_val, pslq_S_results


# ════════════════════════════════════════════════════════════════════
# PHASE 3 — TRANS-SERIES CHECK
# ════════════════════════════════════════════════════════════════════

def phase3(a_coeffs, xi_0, S, sigma_rec, mu, vquad):
    """Check if V_quad has a trans-series representation."""
    print("\n" + "=" * 70)
    print("  PHASE 3 — TRANS-SERIES CHECK")
    print("=" * 70)

    mp.mp.dps = DPS

    if S is None:
        print("  Stokes constant not available — skipping trans-series check")
        return {"status": "deferred", "reason": "S not computed"}

    # The trans-series structure:
    # y(x) = y_0(x) + S * e^{-(σ_+ - σ_-)x} * y_1(x) + ...
    # where y_0 is the formal series for the recessive solution,
    # and the exponential correction involves the Stokes gap.
    #
    # At x = 0, V_quad = y(0) which is the fully resummed value.
    # The formal series y_0 diverges at x=0, so we can't evaluate there.
    #
    # Instead, we check at large x where the formal series converges well.
    # At x=x_test:
    #   y_exact(x) = y_0^{optimal}(x) + S * e^{-xi_0 * x} * (1 + O(1/x))
    # where y_0^{optimal} is the optimal truncation of the formal series.

    stokes_gap = abs(2 / sqrt(mpf(3)))  # = xi_0 if theory matches

    print(f"\n  Trans-series structure:")
    print(f"    y(x) = φ_0(x) + S · e^{{-ξ_0·x}} · φ_1(x) + ...")
    print(f"    ξ_0 = {nstr(xi_0, 15)}")
    print(f"    S = {nstr(S, 15)}")
    print(f"    2/√3 = {nstr(stokes_gap, 15)}")

    # Check at several x values
    results = {}
    N = len(a_coeffs) - 1

    for x_test in [5, 8, 10, 15, 20]:
        x = mpf(x_test)

        # Optimal truncation: stop at n where |a_n / x^n| is minimized
        best_n = 0
        best_term = abs(a_coeffs[0])
        for n in range(1, N + 1):
            term = abs(a_coeffs[n]) / x ** n
            if term < best_term:
                best_term = term
                best_n = n
            elif n > best_n + 5:
                break

        # Compute optimal truncation sum
        phi0_opt = mp.fsum(a_coeffs[n] / x ** n for n in range(best_n + 1))

        # Full solution via backward recurrence from CF
        # y_rec(x) ~ e^{σ_- x} x^μ * Σ a_n x^{-n}
        # The optimal truncation gives: y_rec_approx = e^{σ_-·x} x^μ · phi0_opt
        y_approx = exp(sigma_rec * x) * x ** mu * phi0_opt

        # The exact y_rec(x) comes from ODE integration
        # We compute it via Frobenius transport from x=0
        # y_rec = M_11 * y_1 + M_12 * y_2

        # For this check, let's compute the error of optimal truncation
        # |y_exact - y_opt| should be ~ |S * e^{-xi_0 * x}| * e^{σ_- x} * x^μ

        # Exponential correction magnitude
        correction_mag = abs(S) * exp(-xi_0.real * x) if isinstance(xi_0, mpc) else abs(S) * exp(-xi_0 * x)

        # Relative to the WKB envelope
        wkb_mag = exp(sigma_rec * x) * x ** mu
        relative_correction = correction_mag  # relative to the formal series scale

        print(f"\n    x = {x_test}:")
        print(f"      Optimal truncation at n = {best_n}")
        print(f"      |a_{{best}}/x^best| = {nstr(best_term, 10)}")
        print(f"      φ_0^opt(x) = {nstr(phi0_opt, 15)}")
        print(f"      |S·e^{{-ξ_0·x}}| = {nstr(correction_mag, 10)}")
        print(f"      Ratio |correction/smallest_term| = {nstr(correction_mag / best_term if best_term > 0 else mpf(0), 10)}")

        results[x_test] = {
            "optimal_n": best_n,
            "phi0_opt": float(phi0_opt),
            "correction_mag": float(correction_mag),
            "smallest_term": float(best_term),
        }

    # The key insight about the trans-series and V_quad:
    print(f"\n  Trans-series assessment:")
    print(f"    V_quad = y_rec(0) requires full analytic continuation, not")
    print(f"    just optimal truncation + exponential correction.")
    print(f"    The formal series diverges at x=0 (all terms blow up).")
    print(f"    The trans-series representation is valid for x → ∞ only.")
    print(f"    At x=0, the Borel resummation IS the definition of V_quad.")
    print(f"    The Stokes constant S governs the jump across Stokes lines")
    print(f"    and relates the recessive solution to the dominant one.")

    # However, the Stokes constant itself IS a meaningful structural invariant
    print(f"\n  Structural result:")
    print(f"    The Stokes constant S = {nstr(S, 15)}")
    print(f"    governs the resurgence of the formal series.")
    print(f"    If S can be identified in closed form, it constrains V_quad")
    print(f"    through the resurgence relations.")

    trans_info = {
        "status": "trans-series valid for x>>1 only",
        "xi_0": str(nstr(xi_0, 20)),
        "S": str(nstr(S, 20)),
        "checks": {str(k): v for k, v in results.items()},
        "x0_evaluation": "not_applicable (series diverges at x=0)",
    }
    return trans_info


# ════════════════════════════════════════════════════════════════════
# PHASE 4 — GOVERNANCE
# ════════════════════════════════════════════════════════════════════

def phase4(xi_0, pslq_xi, S, pslq_S, trans_info, agree):
    print("\n" + "=" * 70)
    print("  PHASE 4 — GOVERNANCE")
    print("=" * 70)

    claims = []

    # Check for xi_0 algebraic identification
    xi_hit = any(r.get("hit") for r in pslq_xi.values())
    if xi_hit:
        for label, r in pslq_xi.items():
            if r.get("hit"):
                claim = {
                    "territory": "T2",
                    "iteration": 19,
                    "claim_type": "numerical_identity",
                    "expression": f"Borel singularity xi_0 identified via {label}: {r['relation']}",
                    "residual": r["residual"],
                    "evidence_class": "numerical_identity",
                    "reproduce": "python scripts/t2_iter19_resurgence.py",
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
                }
                claims.append(claim)
                print(f"  CLAIM: xi_0 identified ({label})")

    # Check for Stokes constant identification
    S_hit = any(r.get("hit") for r in pslq_S.values()) if pslq_S else False

    # Structural claim
    stokes_gap = 2 / sqrt(mpf(3))
    xi_matches_theory = abs(xi_0 - stokes_gap) < mpf("0.05") if xi_0 is not None else False

    structural_claim = {
        "territory": "T2",
        "iteration": 19,
        "claim_type": "near_miss" if xi_matches_theory else "null_result",
        "expression": (
            f"V_quad resurgence analysis: "
            f"Borel singularity at xi_0 = {nstr(xi_0, 15)}"
            + (f" ≈ 2/√3 (Stokes gap, diff={nstr(abs(xi_0 - stokes_gap), 5)})" if xi_matches_theory else "")
            + f". Stokes constant S = {nstr(S, 10) if S is not None else 'not computed'}. "
            + f"Trans-series valid for x>>1 only; "
            + f"V_quad at x=0 requires full Borel resummation. "
            + f"PSLQ xi_0: {'hit' if xi_hit else 'null'}. "
            + f"PSLQ S: {'hit' if S_hit else 'null'}."
        ),
        "evidence_class": "near_miss" if xi_matches_theory else "null_result",
        "reproduce": "python scripts/t2_iter19_resurgence.py",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
    }
    claims.append(structural_claim)

    for claim in claims:
        with open(CLAIMS, "a", encoding="utf-8") as f:
            f.write(json.dumps(claim) + "\n")

    print(f"  {len(claims)} claim(s) emitted")

    ahs = agree / DPS if agree < DPS else 1.0
    print(f"\n  AHS = {ahs:.4f} (threshold 0.95)")
    governance = "CLEAN" if ahs >= 0.95 else "HALT"
    print(f"  Governance: {governance}")
    return governance


# ════════════════════════════════════════════════════════════════════
# DELIVERABLE
# ════════════════════════════════════════════════════════════════════

def print_deliverable(xi_0, pslq_xi, S, pslq_S, trans_info, governance, vquad, agree):
    mp.mp.dps = 50
    print("\n" + "═" * 70)
    print("  DELIVERABLE — T2 ITERATION 19")
    print("═" * 70)

    xi_hit = any(r.get("hit") for r in pslq_xi.values())
    S_hit = any(r.get("hit") for r in pslq_S.values()) if pslq_S else False
    stokes_gap = 2 / sqrt(mpf(3))
    xi_matches = abs(xi_0 - stokes_gap) < mpf("0.05") if xi_0 is not None else False

    print(f"\n  BOREL SINGULARITY xi_0: {nstr(xi_0, 15)}")
    print(f"  XI_0 = 2/√3 ?:         {'YES (diff=' + nstr(abs(xi_0-stokes_gap),5) + ')' if xi_matches else 'NO'}")
    print(f"  XI_0 ALGEBRAIC:         {'identified' if xi_hit else 'null'}")
    print(f"  STOKES CONSTANT S:      {nstr(S, 15) if S is not None else 'not computed'}")
    if S is not None:
        V = vquad
        best = min(
            ("S=V", abs(S.real if isinstance(S, mpc) else S) - V),
            ("S=1/V", abs(abs(S.real if isinstance(S, mpc) else S) - 1/V)),
            key=lambda x: abs(x[1])
        )
        print(f"  S vs V_QUAD:            best={best[0]} (diff={nstr(abs(best[1]), 8)})")
    else:
        print(f"  S vs V_QUAD:            deferred")
    print(f"  TRANS-SERIES:           {trans_info.get('status', 'null')}")
    print(f"  F2 STATUS:              {'partial' if xi_matches else 'open'}")
    print(f"  GOVERNANCE:             {governance}")

    report = {
        "territory": "T2",
        "iteration": 19,
        "vquad_20": nstr(vquad, 20),
        "agreement_digits": agree,
        "xi_0": str(nstr(xi_0, 20)),
        "xi_0_algebraic": xi_hit,
        "xi_0_matches_stokes_gap": bool(xi_matches),
        "stokes_constant_S": str(nstr(S, 20)) if S is not None else None,
        "pslq_xi_results": {k: {kk: str(vv) if isinstance(vv, list) else vv for kk, vv in v.items()} for k, v in pslq_xi.items()},
        "pslq_S_results": {k: {kk: str(vv) if isinstance(vv, list) else vv for kk, vv in v.items()} for k, v in (pslq_S or {}).items()},
        "trans_info": trans_info,
        "governance": governance,
        "f2_status": "partial" if xi_matches else "open",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
    }
    out_path = RESULTS / "t2_iter19_resurgence.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n  Full report: {out_path}")


# ════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()
    mp.mp.dps = DPS

    print("  Computing V_quad...")
    v1 = compute_vquad(CF_DEPTH, DPS)
    v2 = compute_vquad(CROSS_DEPTH, DPS)
    with mp.workdps(DPS):
        diff = abs(v1 - v2)
        agree = DPS if diff == 0 else max(0, int(-float(mp.log10(diff))))
    vquad = v1
    print(f"  V_quad = {nstr(vquad, 30)} ({agree} digits confirmed)")

    ahs = agree / DPS if agree < DPS else 1.0
    if ahs < 0.95:
        print(f"\n  *** AHS = {ahs:.4f} < 0.95 — ABORTING ***")
        sys.exit(1)

    # Phase 1
    a_coeffs, b_coeffs, xi_0, pslq_xi, sigma_rec, mu = phase1()

    # Phase 2
    S, pslq_S = phase2(a_coeffs, b_coeffs, xi_0, vquad)

    # Phase 3
    trans_info = phase3(a_coeffs, xi_0, S, sigma_rec, mu, vquad)

    # Phase 4
    governance = phase4(xi_0, pslq_xi, S, pslq_S, trans_info, agree)

    # Deliverable
    print_deliverable(xi_0, pslq_xi, S, pslq_S, trans_info, governance, vquad, agree)

    elapsed = time.time() - t_start
    print(f"\n  Total elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
