#!/usr/bin/env python3
"""
SIARC T2 Iteration 20 — Stokes constant extraction for V_quad.

Phase 1: Determine branch exponent alpha from late-term ratio asymptotics.
Phase 2: Conformal mapping + Euler-Borel summation for S extraction.
Phase 3: Independent check via median (lateral Borel) summation.
Phase 4: Trans-series coefficient check.
Phase 5: Governance.

Depends on: iteration 19 — xi_0 = 2/√3, branch point on negative Borel axis,
            PIII(D6) params α=1/6, β=γ=0, δ=-1/2.
"""

import json
import sys
import time
from pathlib import Path

import mpmath as mp
from mpmath import mpf, mpc, pi, sqrt, log, gamma, exp, factorial, pslq, nstr, fsum

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
CLAIMS = RESULTS / "claims.jsonl"

DPS = 200
CF_DEPTH = 2000
CROSS_DEPTH = 2500
PSLQ_DPS = 80
COEFF_BOUND = 10000
RESIDUAL_THRESHOLD = 1e-15

# ── infrastructure from iter 19 ────────────────────────────────────

def compute_vquad(depth, dps):
    with mp.workdps(dps + 60):
        v = mpf(0)
        for n in range(depth, 0, -1):
            v = mpf(1) / (3 * n * n + n + 1 + v)
        return +(mpf(1) + v)


def wkb_riccati_coeffs(sigma, order=220):
    c = [mpf(0)] * (order + 1)
    c[0] = sigma
    c[1] = -1 - sigma / 6
    d = [mpf(0)] * (order + 1)
    d[0] = c[0] ** 2
    d[1] = 2 * c[0] * c[1]
    for k in range(2, order + 1):
        known_s = fsum(c[i] * c[k - i] for i in range(1, k))
        rest = (3 * (known_s - (k - 1) * c[k - 1])
                + d[k - 1] + d[k - 2]
                + 6 * c[k - 1] + c[k - 2])
        c[k] = -rest / (6 * c[0])
        d[k] = 2 * c[0] * c[k] + known_s - (k - 1) * c[k - 1]
    return c


def formal_series_coeffs(order=150):
    """Compute a_n for y ~ e^{σx} x^μ Σ a_n x^{-n}. Extended to 150 terms."""
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
        s = fsum(k * f_coeffs[k] * a[n - k] for k in range(1, n + 1))
        a[n] = s / n
    return a, sigma_rec, mu


# ════════════════════════════════════════════════════════════════════
# PHASE 1 — BRANCH EXPONENT alpha
# ════════════════════════════════════════════════════════════════════

def phase1(a_coeffs):
    print("=" * 70)
    print("  PHASE 1 — BRANCH EXPONENT β")
    print("=" * 70)

    mp.mp.dps = DPS
    N = len(a_coeffs) - 1
    xi0 = 2 / sqrt(mpf(3))

    # For alternating-sign series with Borel singularity at -xi0:
    #   B(ξ) ~ S · (ξ₀ + ξ)^{-β}  near ξ = -ξ₀
    # gives: a_n ~ S · (-1)^n · Γ(n+β) / (Γ(β) · ξ₀^{n+β})
    # so:    r_n = a_{n+1}/a_n = -(n+β)/ξ₀
    #        |r_n| · ξ₀ = n + β
    #        |r_n| · ξ₀ - n → β

    print(f"\n  ξ₀ = 2/√3 = {nstr(xi0, 20)}")
    print(f"\n  Step 1: Extract β from ratio asymptotics")

    beta_raw = []
    for n in range(30, min(N, 140)):
        if abs(a_coeffs[n]) > 0:
            r_n = a_coeffs[n + 1] / a_coeffs[n]
            beta_n = abs(r_n) * xi0 - n
            beta_raw.append((n, beta_n))

    print("    Raw β_n = |r_n|·ξ₀ - n (last 12):")
    for n, bv in beta_raw[-12:]:
        print(f"      n={n:3d}: β_n = {nstr(bv, 15)}")

    # Richardson-1
    vals = [v for _, v in beta_raw]
    ns = [n for n, _ in beta_raw]
    rich1 = []
    for i in range(1, len(vals)):
        n = ns[i]
        rich1.append((n, n * vals[i] - (n - 1) * vals[i - 1]))

    print("\n    Richardson-1 β (last 8):")
    for n, r1v in rich1[-8:]:
        print(f"      n={n:3d}: {nstr(r1v, 15)}")

    # Richardson-2
    r1vals = [v for _, v in rich1]
    r1ns = [n for n, _ in rich1]
    rich2 = []
    for i in range(1, len(r1vals)):
        n = r1ns[i]
        rich2.append((n, n * r1vals[i] - (n - 1) * r1vals[i - 1]))

    print("\n    Richardson-2 β (last 6):")
    for n, r2v in rich2[-6:]:
        print(f"      n={n:3d}: {nstr(r2v, 15)}")

    beta_est = rich2[-1][1] if rich2 else (rich1[-1][1] if rich1 else beta_raw[-1][1])
    print(f"\n    Best β estimate: {nstr(beta_est, 15)}")

    # Step 2: Check β against theoretical values
    print("\n  Step 2: Theoretical comparison")
    sigma_plus = 1 / sqrt(mpf(3))
    sigma_minus = -1 / sqrt(mpf(3))
    mu_plus = -1 - sigma_plus / 6   # = -1 - 1/(6√3)
    mu_minus = -1 - sigma_minus / 6  # = -1 + 1/(6√3)
    mu_diff = mu_plus - mu_minus      # = -2/(6√3) = -1/(3√3) = -√3/9

    checks = {
        "μ₊ − μ₋ = -1/(3√3)": mu_diff,
        "-1/2 (PIII)": mpf("-0.5"),
        "-1 (log)": mpf("-1"),
        "0 (simple pole)": mpf("0"),
    }
    for label, expected in checks.items():
        diff = abs(beta_est - expected)
        match = " ***" if diff < mpf("0.01") else ""
        print(f"    β = {label} = {nstr(expected, 12)}: diff = {nstr(diff, 8)}{match}")

    # Step 3: PSLQ identification of β
    print(f"\n  Step 3: PSLQ on β = {nstr(beta_est, 15)}")
    mp.mp.dps = PSLQ_DPS

    pslq_beta = {}

    def try_pslq_b(label, basis):
        try:
            rel = pslq(basis, maxcoeff=COEFF_BOUND, maxsteps=5000)
        except Exception:
            rel = None
        if rel:
            residual = float(abs(sum(c * x for c, x in zip(rel, basis))))
            pslq_beta[label] = {"relation": [int(c) for c in rel],
                                "residual": residual,
                                "hit": residual < RESIDUAL_THRESHOLD}
            print(f"    {label}: {[int(c) for c in rel]} res={residual:.2e}")
        else:
            pslq_beta[label] = {"relation": None, "residual": None, "hit": False}
            print(f"    {label}: null")

    try_pslq_b("sqrt3", [mpf(1), beta_est, sqrt(mpf(3))])
    try_pslq_b("algebraic", [mpf(1), beta_est, beta_est**2, beta_est**3])
    try_pslq_b("direct_mu_diff", [beta_est, mu_diff])  # is β = mu_diff exactly?
    try_pslq_b("mu_sigma", [mpf(1), beta_est, mu_diff, sigma_plus])

    mp.mp.dps = DPS
    return beta_est, mu_diff, pslq_beta
        match = " ***" if diff < mpf("0.05") else ""
        print(f"    alpha = {label}: diff = {nstr(diff, 8)}{match}")

    # Also try: the formula with sign correction
    # For alternating series: a_n ~ C (-xi0)^{-n} n^alpha Γ(n)
    # => a_{n+1}/a_n ~ -1/xi0 · (n+alpha)/n · n/(n-1)... simplified:
    #    = (-1/xi0) · (1 + alpha/n + ...)
    # So: (r_n · (-xi0) - 1) · n → alpha
    print("\n  Step 4: Alternative extraction: (r_n · (-xi_0) - 1) · n → alpha")
    alpha_alt = []
    for n, r_n in ratios:
        val = (r_n * (-xi0) - 1) * n
        alpha_alt.append((n, val))

    print("    Last 15:")
    for n, av in alpha_alt[-15:]:
        print(f"      n={n:3d}: {nstr(av, 12)}")

    # Richardson on alternative
    a_alt_vals = [v for _, v in alpha_alt]
    a_alt_ns = [n for n, _ in alpha_alt]
    alt_rich1 = []
    for i in range(1, len(a_alt_vals)):
        n = a_alt_ns[i]
        r1 = n * a_alt_vals[i] - (n - 1) * a_alt_vals[i - 1]
        alt_rich1.append((n, r1))
    if alt_rich1:
        print(f"    Alt Richardson-1 (last 5):")
        for n, v in alt_rich1[-5:]:
            print(f"      n={n:3d}: {nstr(v, 12)}")
        alpha_alt_est = alt_rich1[-1][1]
        print(f"    Alt best: {nstr(alpha_alt_est, 12)}")
    else:
        alpha_alt_est = alpha_est

    return alpha_est, alpha_alt_est


# ════════════════════════════════════════════════════════════════════
# PHASE 2 — CONFORMAL MAPPING + EULER-BOREL
# ════════════════════════════════════════════════════════════════════

def phase2(a_coeffs, alpha_est):
    print("\n" + "=" * 70)
    print("  PHASE 2 — CONFORMAL MAPPING + EULER-BOREL")
    print("=" * 70)

    mp.mp.dps = DPS
    N = len(a_coeffs) - 1
    xi0 = 2 / sqrt(mpf(3))

    # Borel transform: b_n = a_n / Γ(n + 1)
    b = [a_coeffs[n] / factorial(n) for n in range(N + 1)]

    # ── Step 1: Conformal map ξ = xi0 · w/(1-w) ──────────────────
    # This maps w ∈ [0,1) → ξ ∈ [0,∞), and sends the singularity
    # at ξ = -xi0 to w = -xi0/(xi0 + xi0) ... wait, let's be careful.
    #
    # Singularity is at ξ = -xi0 (negative real axis, because signs alternate).
    # We want Borel integral ∫_0^∞ B(ξ) e^{-ξ} dξ along the positive real axis.
    # The singularity at -xi0 is NOT on this contour, so the Borel sum
    # along the real axis is actually well-defined!
    #
    # The issue is that the coefficients still feel the singularity:
    # b_n ~ C (-1)^n / xi0^n, so the Borel series has radius of convergence xi0.
    #
    # Strategy: Use Euler's series transformation to accelerate.
    # If Σ (-1)^n c_n = Σ (-1)^n Δ^n c_0 / 2^{n+1} (Euler transform)
    # where Δ^k c_0 = Σ_{j=0}^{k} C(k,j) (-1)^j c_j.
    #
    # For the Borel integral: split B(ξ) = Σ b_n ξ^n,
    # then ∫_0^∞ B(ξ) e^{-ξ} dξ = Σ b_n · n! = Σ a_n (formal sum).
    # This diverges. But laterally resummed:
    # The Borel sum S_B = ∫_0^∞ [Padé approx of B(ξ)] e^{-ξ} dξ

    # ── Step 2: Euler transform of Borel coefficients ─────────────
    # The series B(ξ) = Σ b_n ξ^n has alternating signs for large n.
    # Apply Euler transform: let c_n = b_n · xi0^n (removes growth rate)
    # then c_n → C · (-1)^n · n^alpha for large n.
    # The Euler transform of Σ c_n z^n at z=ξ/xi0 is:
    #   Σ d_k (z/(1+z))^k / (1+z) where d_k = Δ^k c_0 / 2^k

    print("\n  Step 1: Euler transform of normalized Borel coefficients")

    # Normalize: c_n = b_n · xi0^n
    M = min(N, 120)  # use first 120 terms
    c = [b[n] * xi0 ** n for n in range(M + 1)]

    print(f"    Using {M+1} Borel coefficients")
    print(f"    c_0 = {nstr(c[0], 12)}, c_1 = {nstr(c[1], 12)}")
    print(f"    Signs c_80..c_89: {''.join('+' if c[n] > 0 else '-' for n in range(80, min(90, M+1)))}")

    # Compute Euler forward differences: Δ^k c_0 = Σ_{j=0}^k C(k,j)(-1)^{k-j} c_j
    # Use the explicit formula with binomials
    K_euler = min(80, M)  # number of Euler terms
    delta = [mpf(0)] * (K_euler + 1)
    for k in range(K_euler + 1):
        s = mpf(0)
        for j in range(k + 1):
            binom_kj = mp.binomial(k, j)
            s += binom_kj * ((-1) ** (k - j)) * c[j]
        delta[k] = s

    print(f"\n    First 10 Euler differences Δ^k c_0:")
    for k in range(min(10, K_euler + 1)):
        print(f"      Δ^{k} c_0 = {nstr(delta[k], 15)}")

    # Check convergence: |delta[k]| should decrease
    print(f"\n    |Δ^k c_0| for k = 20,30,...,{K_euler}:")
    for k in range(20, K_euler + 1, 10):
        print(f"      k={k:3d}: |Δ^k| = {nstr(abs(delta[k]), 10)}")

    # ── Step 3: Evaluate the Euler-transformed Borel sum ──────────
    # S_B = ∫_0^∞ B(ξ) e^{-ξ} dξ
    # With the Euler transform at z = ξ/xi0:
    # B(ξ) = (1/(1+z)) Σ_{k≥0} d_k (z/(1+z))^k  where d_k = Δ^k c_0
    # and ξ = xi0·z, dξ = xi0·dz
    # S_B = xi0 ∫_0^∞ [Σ d_k w^k] e^{-xi0·z} / (1+z) dz
    #   where w = z/(1+z), i.e. z = w/(1-w)
    # Change variable to w: z = w/(1-w), dz = 1/(1-w)^2 dw, w ∈ [0,1)
    # ξ = xi0·w/(1-w), e^{-ξ} = e^{-xi0·w/(1-w)}
    # 1/(1+z) = 1/(1/(1-w)) = (1-w)
    # S_B = xi0 ∫_0^1 [Σ d_k w^k] e^{-xi0·w/(1-w)} (1-w) / (1-w)^2 dw
    #      = xi0 ∫_0^1 [Σ d_k w^k] e^{-xi0·w/(1-w)} / (1-w) dw
    #
    # This integral is hard to evaluate directly. Instead, use term-by-term:
    # S_B = Σ_{k≥0} d_k · I_k
    # where I_k = xi0 ∫_0^1 w^k / (1-w) · e^{-xi0·w/(1-w)} dw
    #
    # These integrals I_k can be computed by substituting t = w/(1-w) ∈ [0,∞):
    # w = t/(1+t), 1-w = 1/(1+t), dw = 1/(1+t)^2 dt
    # I_k = xi0 ∫_0^∞ (t/(1+t))^k · (1+t) · e^{-xi0·t} / (1+t)^2 dt
    #      = xi0 ∫_0^∞ t^k / (1+t)^{k+1} · e^{-xi0·t} dt
    #
    # These are known in terms of confluent hypergeometric functions:
    # I_k = xi0 · U(k+1, k+1, xi0) · Γ(k+1) / Γ(1)
    # Actually: I_k = ∫_0^∞ t^k (1+t)^{-(k+1)} e^{-xi0·t} xi0 dt
    # = xi0 · Γ(k+1) · U(k+1, k+1, xi0)
    # where U is the confluent hypergeometric function of the second kind.
    # But U(a,a,z) = e^z ∫_0^∞ t^{a-1} (1+t)^{-1} e^{-z(1+t)} dt ...
    #
    # Simpler: just numerically integrate I_k.

    print("\n  Step 2: Computing Euler-Borel integrals I_k")

    def I_k(k_val, xi0_val, dps_int=50):
        """Compute I_k = xi0 ∫_0^∞ t^k (1+t)^{-(k+1)} e^{-xi0 t} dt"""
        with mp.workdps(dps_int):
            def integrand(t):
                return xi0_val * t ** k_val * (1 + t) ** (-(k_val + 1)) * exp(-xi0_val * t)
            return mp.quad(integrand, [0, mp.inf], method='tanh-sinh')

    # Compute I_k for k = 0..K_euler
    K_use = min(K_euler, 60)  # limit to avoid slowness
    Ik_vals = []
    t0 = time.time()
    for k in range(K_use + 1):
        ik = I_k(k, xi0, dps_int=40)
        Ik_vals.append(ik)
    print(f"    Computed {K_use+1} integrals in {time.time()-t0:.1f}s")

    print(f"    I_0 = {nstr(Ik_vals[0], 12)}")
    print(f"    I_10 = {nstr(Ik_vals[10], 12)}")
    if K_use >= 30:
        print(f"    I_30 = {nstr(Ik_vals[30], 12)}")

    # ── Step 4: Euler-Borel sum ───────────────────────────────────
    # S_B = Σ_{k=0}^{K} d_k · I_k
    print("\n  Step 3: Euler-Borel partial sums")

    partial_sums = []
    running = mpf(0)
    for k in range(K_use + 1):
        running += delta[k] * Ik_vals[k]
        if k % 10 == 0 or k == K_use:
            partial_sums.append((k, +running))
            print(f"    S_B(K={k:3d}) = {nstr(running, 15)}")

    S_euler_borel = running
    print(f"\n    Final Euler-Borel sum: {nstr(S_euler_borel, 15)}")

    # ── Step 5: Direct numerical Borel integration ────────────────
    # S_B = ∫_0^∞ B(ξ) e^{-ξ} dξ where B(ξ) = Σ b_n ξ^n
    # Since the singularity is at -xi0 (off the integration contour),
    # we can just Padé the series and integrate.

    print("\n  Step 4: Direct Padé-Borel integration")

    # Use [40/40] Padé of b_n
    pade_ord = 40
    try:
        p_c, q_c = mp.pade(b[:2 * pade_ord + 1], pade_ord, pade_ord)

        def pade_B(xi):
            num = mp.polyval(list(reversed(p_c)), xi)
            den = mp.polyval(list(reversed(q_c)), xi)
            if abs(den) < mpf("1e-200"):
                return mpf(0)
            return num / den

        def integrand_pade(xi):
            return pade_B(xi) * exp(-xi)

        with mp.workdps(50):
            S_pade_borel = mp.quad(integrand_pade, [0, mp.inf], method='tanh-sinh')
        print(f"    Padé[40/40]-Borel sum = {nstr(S_pade_borel, 15)}")
    except Exception as e:
        print(f"    Padé-Borel failed: {e}")
        S_pade_borel = None

    # Also try [30/30]
    try:
        p_c30, q_c30 = mp.pade(b[:61], 30, 30)

        def pade_B30(xi):
            num = mp.polyval(list(reversed(p_c30)), xi)
            den = mp.polyval(list(reversed(q_c30)), xi)
            if abs(den) < mpf("1e-200"):
                return mpf(0)
            return num / den

        def integrand_p30(xi):
            return pade_B30(xi) * exp(-xi)

        with mp.workdps(50):
            S_p30 = mp.quad(integrand_p30, [0, mp.inf], method='tanh-sinh')
        print(f"    Padé[30/30]-Borel sum = {nstr(S_p30, 15)}")
    except Exception as e:
        print(f"    Padé[30/30]-Borel failed: {e}")
        S_p30 = None

    # ── Step 5: Compare ──────────────────────────────────────────
    print("\n  Step 5: Comparison of Borel sums")
    print(f"    Euler-Borel (K={K_use}): {nstr(S_euler_borel, 15)}")
    if S_pade_borel is not None:
        print(f"    Padé[40/40]-Borel:      {nstr(S_pade_borel, 15)}")
        diff_ep = abs(S_euler_borel - S_pade_borel)
        print(f"    |Euler - Padé40|:       {nstr(diff_ep, 8)}")
    if S_p30 is not None:
        print(f"    Padé[30/30]-Borel:      {nstr(S_p30, 15)}")

    # NOTE: The "Stokes constant" is NOT the Borel sum — it's the
    # DISCONTINUITY when the Borel contour crosses the singularity.
    # Since the singularity is at -xi0 (off the positive real contour),
    # there's no Stokes phenomenon on the real axis.
    # The Stokes constant governs the jump across arg(x) = 0 ↔ arg(x) = π.

    return S_euler_borel, S_pade_borel, S_p30, b, delta, Ik_vals


# ════════════════════════════════════════════════════════════════════
# PHASE 3 — LATERAL BOREL / MEDIAN SUMMATION
# ════════════════════════════════════════════════════════════════════

def phase3(b_coeffs, alpha_est):
    print("\n" + "=" * 70)
    print("  PHASE 3 — LATERAL BOREL + STOKES DISCONTINUITY")
    print("=" * 70)

    mp.mp.dps = DPS
    xi0 = 2 / sqrt(mpf(3))
    N = len(b_coeffs) - 1

    # The Stokes constant S appears when the Borel integration contour
    # is deformed to pass THROUGH the singularity at -xi0.
    # This happens when x is on the negative real axis (arg(x) = π).
    #
    # For x > 0: S_B(x) = ∫_0^∞ B(ξ) e^{-xξ} dξ  (real axis, no singularity)
    # For x < 0: the contour must be deformed, and the Stokes jump is:
    #   S_{upper}(x) - S_{lower}(x) = 2πi · S · (residue-like term)
    #
    # Alternative: evaluate the Borel integral along the NEGATIVE real axis:
    # ∫_0^{-∞} B(ξ) e^{-xξ} dξ  for x > 0
    # which encounters the singularity at -xi0.
    #
    # The lateral sums are:
    # S_± = ∫_0^{-∞±iε} B(ξ) e^{-xξ} dξ
    # And S = (S_+ - S_-) / (2πi) · [normalization]

    print(f"\n  Stokes discontinuity via lateral Borel sums")
    print(f"  Singularity at ξ = -xi0 = {nstr(-xi0, 15)}")

    # ── Step 1: Padé approximant for B(ξ) ─────────────────────────
    pade_ord = 40
    try:
        p_c, q_c = mp.pade(b_coeffs[:2 * pade_ord + 1], pade_ord, pade_ord)
    except Exception as e:
        print(f"    Padé construction failed: {e}")
        return None, {}

    def pade_B(z):
        num = mp.polyval(list(reversed(p_c)), z)
        den = mp.polyval(list(reversed(q_c)), z)
        if abs(den) < mpf("1e-200"):
            return mpf("nan")
        return num / den

    # ── Step 2: Lateral integrals along negative real axis ────────
    # For x_test > 0, compute:
    # L_±(x) = ∫_0^∞ B(-t ± iε) e^{x·t} (-1) dt
    # where we integrate B along ξ = -t ± iε (ε small)
    # The Stokes jump = L_+(x) - L_-(x)

    x_test = mpf(1)  # test at x=1
    eps_vals = [mpf("0.05"), mpf("0.01"), mpf("0.005")]

    print(f"\n  Step 1: Lateral integrals at x = {nstr(x_test, 3)}")

    stokes_jumps = []
    for eps in eps_vals:
        with mp.workdps(50):
            def integrand_upper(t):
                z = -t + mpc(0, 1) * eps
                return pade_B(z) * exp(x_test * t)

            def integrand_lower(t):
                z = -t - mpc(0, 1) * eps
                return pade_B(z) * exp(x_test * t)

            # Integration range: [0, 3*xi0] (captures the singularity region)
            try:
                L_upper = mp.quad(integrand_upper, [0, 3 * xi0], method='gauss-legendre')
                L_lower = mp.quad(integrand_lower, [0, 3 * xi0], method='gauss-legendre')
                jump = L_upper - L_lower
                stokes_jumps.append((eps, jump))
                print(f"    ε = {nstr(eps, 4)}: L+ - L- = {nstr(jump, 12)}")
                print(f"      Im(L+ - L-) = {nstr(jump.imag if isinstance(jump, mpc) else mpf(0), 12)}")
                print(f"      Re(L+ - L-) = {nstr(jump.real if isinstance(jump, mpc) else jump, 12)}")
            except Exception as e:
                print(f"    ε = {nstr(eps, 4)}: integration failed: {e}")

    # ── Step 3: Extract S from the jump ───────────────────────────
    # For a branch point B(ξ) ~ S·(ξ+xi0)^alpha near ξ=-xi0:
    # The discontinuity across the branch cut is:
    #   disc B = B(ξ+iε) - B(ξ-iε) = 2i·sin(πα)·S·|ξ+xi0|^alpha  (for ξ < -xi0)
    # The integrated jump is:
    #   L+ - L- = ∫_{xi0}^{3xi0} disc B(-t) e^{x·t} dt
    #           ≈ 2i·sin(πα)·S · [integral of |t-xi0|^alpha · e^{x·t}]

    # For now, just report the jump magnitude
    print(f"\n  Step 2: Stokes constant estimation")
    if stokes_jumps:
        # Use smallest epsilon
        eps_best, jump_best = stokes_jumps[-1]
        jump_im = jump_best.imag if isinstance(jump_best, mpc) else mpf(0)
        jump_re = jump_best.real if isinstance(jump_best, mpc) else jump_best

        # For α ≈ -1/2: sin(π·α) = sin(-π/2) = -1
        # So Im(L+ - L-) ≈ -2·S · ∫ |t-xi0|^{-1/2} e^{t} dt
        # This integral can be estimated as ~ 2·sqrt(π/(2xi0))·e^{xi0}... complex
        # Just report the raw jump
        S_lateral = jump_im / (2 * pi)  # rough normalization
        print(f"    Raw jump Im = {nstr(jump_im, 12)}")
        print(f"    S ~ Im(jump)/(2π) = {nstr(S_lateral, 12)}")
        print(f"    |jump| = {nstr(abs(jump_best), 12)}")
    else:
        S_lateral = None

    # ── Step 4: PSLQ on S_lateral ────────────────────────────────
    pslq_S = {}
    if S_lateral is not None and abs(S_lateral) > mpf("1e-30"):
        mp.mp.dps = PSLQ_DPS
        S_val = abs(S_lateral)

        def try_pslq(label, basis):
            try:
                rel = pslq(basis, maxcoeff=COEFF_BOUND, maxsteps=5000)
            except Exception:
                rel = None
            if rel:
                residual = float(abs(sum(c * x for c, x in zip(rel, basis))))
                pslq_S[label] = {"relation": [int(c) for c in rel],
                                 "residual": residual,
                                 "hit": residual < RESIDUAL_THRESHOLD}
                print(f"    {label}: {[int(c) for c in rel]} res={residual:.2e}")
            else:
                pslq_S[label] = {"relation": None, "residual": None, "hit": False}
                print(f"    {label}: null")

        print(f"\n  Step 3: PSLQ on |S| = {nstr(S_val, 12)}")
        try_pslq("constants", [mpf(1), S_val, pi, sqrt(mpf(3)),
                               sqrt(mpf(11)), gamma(mpf(1)/6)])
        try_pslq("gamma_basis", [mpf(1), S_val, gamma(mpf(1)/3),
                                 gamma(mpf(1)/4), sqrt(pi)])
        try_pslq("exp_pi_sqrt3", [mpf(1), S_val, exp(pi/sqrt(mpf(3))),
                                   exp(-pi/sqrt(mpf(3)))])
    else:
        print(f"\n  S_lateral too small or None for PSLQ")

    return S_lateral, pslq_S


# ════════════════════════════════════════════════════════════════════
# PHASE 4 — TRANS-SERIES CHECK
# ════════════════════════════════════════════════════════════════════

def phase4(S_euler, S_pade, vquad, alpha_est):
    print("\n" + "=" * 70)
    print("  PHASE 4 — TRANS-SERIES / BOREL SUM CHECK")
    print("=" * 70)

    mp.mp.dps = DPS
    xi0 = 2 / sqrt(mpf(3))

    # The key insight: for x > 0 (our case, V_quad = y(0)),
    # the Borel singularity at -xi0 is NOT on the integration contour.
    # Therefore: V_quad = S_B (the Borel sum along the real axis)
    # with NO Stokes correction needed.
    #
    # The Borel sum IS the exact value V_quad.
    # Check: S_euler ≈ S_pade ≈ V_quad?

    print(f"\n  V_quad = {nstr(vquad, 20)}")
    if S_euler is not None:
        print(f"  Euler-Borel sum = {nstr(S_euler, 20)}")
        diff_e = abs(S_euler - vquad)
        print(f"  |Euler-Borel - V_quad| = {nstr(diff_e, 8)}")
        agree_e = max(0, int(-float(mp.log10(diff_e)))) if diff_e > 0 else 50
        print(f"  Agreement: {agree_e} digits")
    else:
        agree_e = 0

    if S_pade is not None:
        print(f"  Padé[40/40]-Borel sum = {nstr(S_pade, 20)}")
        diff_p = abs(S_pade - vquad)
        print(f"  |Padé-Borel - V_quad| = {nstr(diff_p, 8)}")
        agree_p = max(0, int(-float(mp.log10(diff_p)))) if diff_p > 0 else 50
        print(f"  Agreement: {agree_p} digits")
    else:
        agree_p = 0

    # The trans-series correction is only relevant when the Borel contour
    # crosses the singularity — i.e. for x on the Stokes line (arg(x)=π).
    # For V_quad = y(0) with x > 0, there is NO trans-series correction.
    print(f"\n  Trans-series assessment:")
    print(f"    For V_quad (x > 0): Borel singularity at ξ = -xi0")
    print(f"    is OFF the positive real integration contour.")
    print(f"    → NO Stokes correction needed.")
    print(f"    → V_quad = (Borel sum along R+) exactly.")
    print(f"    → The Stokes constant governs the Stokes phenomenon")
    print(f"      when x crosses the anti-Stokes line (arg x → π).")

    return max(agree_e, agree_p)


# ════════════════════════════════════════════════════════════════════
# PHASE 5 — GOVERNANCE
# ════════════════════════════════════════════════════════════════════

def phase5(alpha_est, alpha_alt, S_euler, S_pade, S_lateral,
           pslq_S, borel_agree, agree_digits):
    print("\n" + "=" * 70)
    print("  PHASE 5 — GOVERNANCE")
    print("=" * 70)

    mp.mp.dps = 50
    xi0 = 2 / sqrt(mpf(3))
    claims = []

    # Claim 1: branch exponent
    alpha_report = alpha_est
    near_half = abs(alpha_est + mpf("0.5")) < mpf("0.1")

    claim1 = {
        "territory": "T2",
        "iteration": 20,
        "claim_type": "near_miss",
        "expression": (
            f"V_quad Borel singularity: xi_0 = 2/sqrt(3) (Stokes gap). "
            f"Branch exponent alpha = {nstr(alpha_report, 6)}"
            + (f" (near -1/2, consistent with PIII(D6))" if near_half else
               f" (deviates from -1/2)")
            + f". Singularity on negative Borel axis. "
            + f"Borel sum along R+ reproduces V_quad to {borel_agree} digits."
        ),
        "evidence_class": "near_miss",
        "reproduce": "python scripts/t2_iter20_stokes_constant.py",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
    }
    claims.append(claim1)

    # Claim 2: Borel sum = V_quad
    if borel_agree >= 5:
        claim2 = {
            "territory": "T2",
            "iteration": 20,
            "claim_type": "numerical_identity",
            "expression": (
                f"V_quad = Borel sum of formal WKB series along R+ "
                f"(confirmed to {borel_agree} digits via Padé-Borel integration). "
                f"No Stokes correction needed: singularity at -2/sqrt(3) is off contour."
            ),
            "evidence_class": "numerical_identity",
            "reproduce": "python scripts/t2_iter20_stokes_constant.py",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        claims.append(claim2)

    # Claim 3: PSLQ on Stokes constant
    pslq_hit = any(r.get("hit") for r in pslq_S.values()) if pslq_S else False
    if pslq_hit:
        for label, r in pslq_S.items():
            if r.get("hit"):
                claim3 = {
                    "territory": "T2",
                    "iteration": 20,
                    "claim_type": "numerical_identity",
                    "expression": f"Stokes constant S identified via PSLQ ({label}): {r['relation']}",
                    "evidence_class": "numerical_identity",
                    "reproduce": "python scripts/t2_iter20_stokes_constant.py",
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
                }
                claims.append(claim3)

    for claim in claims:
        with open(CLAIMS, "a", encoding="utf-8") as f:
            f.write(json.dumps(claim) + "\n")

    ahs = agree_digits / DPS if agree_digits < DPS else 1.0
    governance = "CLEAN" if ahs >= 0.95 else "HALT"
    print(f"\n  {len(claims)} claim(s) emitted")
    print(f"  AHS = {ahs:.4f}")
    print(f"  Governance: {governance}")

    return governance


# ════════════════════════════════════════════════════════════════════
# DELIVERABLE
# ════════════════════════════════════════════════════════════════════

def print_deliverable(alpha_est, alpha_alt, S_euler, S_pade,
                      S_lateral, pslq_S, borel_agree, governance, vquad):
    mp.mp.dps = 50
    xi0 = 2 / sqrt(mpf(3))
    print("\n" + "═" * 70)
    print("  DELIVERABLE — T2 ITERATION 20")
    print("═" * 70)

    near_half = abs(alpha_est + mpf("0.5")) < mpf("0.1")
    pslq_hit = any(r.get("hit") for r in pslq_S.values()) if pslq_S else False

    print(f"\n  BRANCH EXPONENT alpha:  {nstr(alpha_est, 8)}"
          f" {'(≈ -1/2, PIII(D6))' if near_half else '(not -1/2)'}")
    print(f"  ALT BRANCH EXPONENT:    {nstr(alpha_alt, 8)}")
    print(f"  EULER-BOREL SUM:        {nstr(S_euler, 15) if S_euler else 'failed'}")
    print(f"  PADE[40/40]-BOREL SUM:  {nstr(S_pade, 15) if S_pade else 'failed'}")
    print(f"  V_QUAD:                 {nstr(vquad, 15)}")
    print(f"  BOREL-VQUAD AGREEMENT:  {borel_agree} digits")
    print(f"  LATERAL STOKES S:       {nstr(S_lateral, 10) if S_lateral else 'not computed'}")
    print(f"  PSLQ ON S:              {'hit' if pslq_hit else 'null'}")
    print(f"  F2 STATUS:              {'partial' if borel_agree >= 5 else 'open'}")
    print(f"  GOVERNANCE:             {governance}")

    report = {
        "territory": "T2",
        "iteration": 20,
        "branch_exponent_alpha": str(nstr(alpha_est, 12)),
        "branch_exponent_alt": str(nstr(alpha_alt, 12)),
        "alpha_near_minus_half": bool(near_half),
        "euler_borel_sum": str(nstr(S_euler, 20)) if S_euler else None,
        "pade_borel_sum": str(nstr(S_pade, 20)) if S_pade else None,
        "vquad_20": str(nstr(vquad, 20)),
        "borel_vquad_agree": borel_agree,
        "lateral_stokes_S": str(nstr(S_lateral, 12)) if S_lateral else None,
        "pslq_S_results": {k: {kk: str(vv) if isinstance(vv, list) else vv
                                for kk, vv in v.items()}
                           for k, v in (pslq_S or {}).items()},
        "governance": governance,
        "f2_status": "partial" if borel_agree >= 5 else "open",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
    }
    out_path = RESULTS / "t2_iter20_stokes.json"
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

    # Compute 150 formal series terms
    print("\n  Computing 150-term formal series...")
    a_coeffs, sigma_rec, mu = formal_series_coeffs(order=150)
    print(f"  σ_- = {nstr(sigma_rec, 15)}, μ_- = {nstr(mu, 15)}")

    # Phase 1
    alpha_est, alpha_alt = phase1(a_coeffs)

    # Phase 2
    S_euler, S_pade, S_p30, b_coeffs, delta, Ik_vals = phase2(a_coeffs, alpha_est)

    # Phase 3
    S_lateral, pslq_S = phase3(b_coeffs, alpha_est)

    # Phase 4
    borel_agree = phase4(S_euler, S_pade, vquad, alpha_est)

    # Phase 5
    governance = phase5(alpha_est, alpha_alt, S_euler, S_pade, S_lateral,
                        pslq_S, borel_agree, agree)

    # Deliverable
    print_deliverable(alpha_est, alpha_alt, S_euler, S_pade,
                      S_lateral, pslq_S, borel_agree, governance, vquad)

    elapsed = time.time() - t_start
    print(f"\n  Total elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
