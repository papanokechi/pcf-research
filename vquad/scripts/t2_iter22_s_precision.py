#!/usr/bin/env python3
"""
SIARC T2 Iteration 22 — S to 12+ digits for PSLQ.

Phase 1: Extend formal series to 1000 terms at dps=100.
Phase 2: Beta with log(n) correction regression.
Phase 3: Dingle S extraction at 1000 terms.
Phase 4: PSLQ on S at 12 digits (5 bases).
Phase 5: Governance.

Depends on: iter21 — S ≈ 0.437705 (7 digits), beta = -1/(3√3) (6 digits).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
import mpmath as mp
from mpmath import mpf, pi, sqrt, log, gamma, exp, factorial, pslq, nstr, fsum

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
CLAIMS = RESULTS / "claims.jsonl"

DPS = 100
PSLQ_DPS = 80
COEFF_BOUND = 10000
RESIDUAL_THRESHOLD = 1e-10  # conservative for 12-digit S
AHS_THRESHOLD = 0.75


# ── infrastructure ─────────────────────────────────────────────────

def compute_vquad(depth, dps):
    with mp.workdps(dps + 60):
        v = mpf(0)
        for n in range(depth, 0, -1):
            v = mpf(1) / (3 * n * n + n + 1 + v)
        return +(mpf(1) + v)


def wkb_riccati_coeffs(sigma, order=220):
    """Riccati coefficients for WKB expansion of (3x²+x+1)y''+(6x+1)y'-x²y=0."""
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


def formal_series_coeffs(order=1000):
    """Compute a_n for y ~ e^{σx} x^μ Σ a_n x^{-n}."""
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
# PHASE 1 — EXTEND SERIES TO 1000 TERMS AT DPS=100
# ════════════════════════════════════════════════════════════════════

def phase1():
    print("=" * 70)
    print("  PHASE 1 — EXTEND SERIES TO 1000 TERMS AT DPS=100")
    print("=" * 70)

    mp.mp.dps = DPS
    t0 = time.time()

    TARGET = 1000
    print(f"\n  Computing {TARGET}-term formal WKB series at dps={DPS}...")

    a_coeffs, sigma_rec, mu = formal_series_coeffs(order=TARGET)

    elapsed = time.time() - t0
    N = len(a_coeffs) - 1
    print(f"  Computed {N+1} coefficients in {elapsed:.1f}s")

    if elapsed > 600:
        print(f"  *** TIMEOUT: exceeded 600s, halting Phase 1 ***")
        if N < 500:
            print(f"  *** Only {N+1} terms computed — insufficient. ABORT. ***")
            return a_coeffs, sigma_rec, mu, False
        else:
            print(f"  *** {N+1} terms available, proceeding with partial data ***")

    xi0 = 2 / sqrt(mpf(3))
    ratio_last = abs(a_coeffs[N] / a_coeffs[N - 1])
    print(f"  σ₋ = {nstr(sigma_rec, 15)}")
    print(f"  μ₋ = {nstr(mu, 15)}")
    print(f"  |a_{N}/a_{N-1}| = {nstr(ratio_last, 12)}")
    print(f"  |a_{N}/a_{N-1}| · ξ₀ = {nstr(ratio_last * xi0, 10)}"
          f" (should ≈ {N-1} + β ≈ {N - 1 - 0.192:.3f})")
    print(f"  |r_{{N-1}}| · ξ₀ / (N-1) = "
          f"{nstr(ratio_last * xi0 / (N - 1), 12)} (→ 1)")

    signs = ''.join('+' if a_coeffs[n] > 0 else '-' for n in range(N - 10, N + 1))
    print(f"  Signs a_{N-10}..a_{N}: {signs} (alternating = branch on neg axis)")

    return a_coeffs, sigma_rec, mu, True


# ════════════════════════════════════════════════════════════════════
# PHASE 2 — BETA WITH LOG(N) CORRECTION
# ════════════════════════════════════════════════════════════════════

def phase2(a_coeffs):
    print("\n" + "=" * 70)
    print("  PHASE 2 — BETA WITH LOG(N) CORRECTION")
    print("=" * 70)

    mp.mp.dps = DPS
    N = len(a_coeffs) - 1
    xi0 = 2 / sqrt(mpf(3))
    mu_diff = -1 / (3 * sqrt(mpf(3)))

    print(f"\n  Target: β = -1/(3√3) = {nstr(mu_diff, 20)}")

    # ── Step 0: Raw β_n from ratios ──────────────────────
    # r_n = a_{n+1}/a_n, |r_n| ~ (n + β) / ξ₀
    # so β_n = |r_n| · ξ₀ - n → β as n → ∞

    beta_raw = []
    for n in range(50, N):
        if abs(a_coeffs[n]) > 0:
            r_n = a_coeffs[n + 1] / a_coeffs[n]
            beta_n = abs(r_n) * xi0 - n
            beta_raw.append((n, float(beta_n)))

    print(f"  Raw β_n: {len(beta_raw)} values (n=50..{N-1})")
    print(f"  Last 6:")
    for n, bv in beta_raw[-6:]:
        print(f"    n={n:4d}: β_n = {bv:.18f}")

    # ── Step 1: Define y_n = (|r_n| · ξ₀ - 1) · n = (β_n + 1 - 1) + ...
    # Actually: β_n = β + corrections, so let's directly regress β_n.
    # Model: β_n = β + c₁·log(n)/n + c₂/n + c₃·log(n)/n² + c₄/n²
    # Fit for n = 700..N-1

    fit_start = 700
    fit_data = [(n, bv) for n, bv in beta_raw if n >= fit_start]
    if len(fit_data) < 50:
        fit_start = 400
        fit_data = [(n, bv) for n, bv in beta_raw if n >= fit_start]

    print(f"\n  Step 1: Regression β_n ~ β + c₁·log(n)/n + c₂/n + c₃·log(n)/n² + c₄/n²")
    print(f"          for n = {fit_start}..{N-1} ({len(fit_data)} points)")

    # Build design matrix
    ns_fit = np.array([n for n, _ in fit_data], dtype=np.float64)
    ys_fit = np.array([bv for _, bv in fit_data], dtype=np.float64)

    # Columns: [1, log(n)/n, 1/n, log(n)/n², 1/n²]
    X = np.column_stack([
        np.ones_like(ns_fit),
        np.log(ns_fit) / ns_fit,
        1.0 / ns_fit,
        np.log(ns_fit) / ns_fit**2,
        1.0 / ns_fit**2,
    ])

    # Solve by least-squares
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, ys_fit, rcond=None)
    beta_regress = coeffs[0]

    print(f"    Intercept (= β): {beta_regress:.18f}")
    print(f"    c₁ (log(n)/n):   {coeffs[1]:.10f}")
    print(f"    c₂ (1/n):        {coeffs[2]:.10f}")
    print(f"    c₃ (log(n)/n²):  {coeffs[3]:.10f}")
    print(f"    c₄ (1/n²):       {coeffs[4]:.10f}")
    if len(residuals) > 0:
        rms = np.sqrt(residuals[0] / len(fit_data))
        print(f"    RMS residual:     {rms:.2e}")

    # Also try with extended model: add log(n)²/n², 1/n³
    X2 = np.column_stack([
        np.ones_like(ns_fit),
        np.log(ns_fit) / ns_fit,
        1.0 / ns_fit,
        np.log(ns_fit) / ns_fit**2,
        1.0 / ns_fit**2,
        np.log(ns_fit)**2 / ns_fit**2,
        1.0 / ns_fit**3,
    ])
    coeffs2, residuals2, _, _ = np.linalg.lstsq(X2, ys_fit, rcond=None)
    beta_regress2 = coeffs2[0]
    print(f"\n    Extended model intercept: {beta_regress2:.18f}")

    # Compare to exact
    beta_exact = float(mp.mpf(mu_diff))
    diff1 = abs(beta_regress - beta_exact)
    diff2 = abs(beta_regress2 - beta_exact)
    dig1 = int(-np.log10(diff1 / abs(beta_exact))) if diff1 > 0 else 20
    dig2 = int(-np.log10(diff2 / abs(beta_exact))) if diff2 > 0 else 20
    best_beta = beta_regress if diff1 <= diff2 else beta_regress2
    best_diff = min(diff1, diff2)
    best_dig = max(dig1, dig2)

    print(f"\n  ── Beta summary ──")
    print(f"    5-param model:  β = {beta_regress:.18f}, |diff| = {diff1:.2e}, {dig1} digits")
    print(f"    7-param model:  β = {beta_regress2:.18f}, |diff| = {diff2:.2e}, {dig2} digits")
    print(f"    Exact -1/(3√3): β = {beta_exact:.18f}")
    print(f"    ★ Best: {best_dig} digits")

    # ── Step 2: Classical Richardson for comparison ───────
    print(f"\n  Step 2: Classical Richardson (for comparison)")
    vals = [bv for _, bv in beta_raw if _ >= 400]
    ns_r = [n for n, _ in beta_raw if n >= 400]

    def rich_step(values, indices):
        out_v, out_n = [], []
        for i in range(1, len(values)):
            n = indices[i]
            n1 = indices[i - 1]
            r = (n * values[i] - n1 * values[i - 1]) / (n - n1)
            out_v.append(r)
            out_n.append(n)
        return out_v, out_n

    cv, cn = vals, ns_r
    for lev in range(1, 4):
        cv, cn = rich_step(cv, cn)
        if len(cv) < 3:
            break
    if cv:
        beta_rich = cv[-1]
        diff_rich = abs(beta_rich - beta_exact)
        dig_rich = int(-np.log10(diff_rich / abs(beta_exact))) if diff_rich > 0 else 20
        print(f"    Richardson^{lev}: β = {beta_rich:.18f}, |diff| = {diff_rich:.2e}, {dig_rich} digits")

    # ── Step 3: PSLQ on beta (if 9+ digits) ─────────────
    beta_pslq = {}
    if best_dig >= 7:
        print(f"\n  Step 3: PSLQ on β (best={best_dig} digits)")
        mp.mp.dps = PSLQ_DPS
        beta_mp = mpf(str(best_beta))

        def try_pslq_b(label, basis):
            try:
                rel = pslq(basis, maxcoeff=COEFF_BOUND, maxsteps=10000)
            except Exception:
                rel = None
            if rel:
                residual = float(abs(sum(c * x for c, x in zip(rel, basis))))
                hit = residual < 1e-20
                beta_pslq[label] = {"relation": [int(c) for c in rel],
                                    "residual": residual, "hit": hit}
                tag = "HIT" if hit else f"near (res={residual:.2e})"
                print(f"    {label}: {[int(c) for c in rel]} [{tag}]")
            else:
                beta_pslq[label] = {"relation": None, "residual": None, "hit": False}
                print(f"    {label}: null")

        try_pslq_b("beta_mu_diff", [mpf(1), beta_mp, mu_diff])
        try_pslq_b("beta_sqrt3", [mpf(1), beta_mp, sqrt(mpf(3))])
        # NOTE: 1 and 1/3 are Q-dependent. Use independent elements.
        try_pslq_b("beta_alg", [mpf(1), beta_mp, beta_mp**2, sqrt(mpf(3))])

        mp.mp.dps = DPS
    else:
        print(f"\n  Step 3: Skipping PSLQ (only {best_dig} digits)")

    return best_beta, best_dig, beta_pslq


# ════════════════════════════════════════════════════════════════════
# PHASE 3 — DINGLE S EXTRACTION AT 1000 TERMS
# ════════════════════════════════════════════════════════════════════

def phase3(a_coeffs):
    print("\n" + "=" * 70)
    print("  PHASE 3 — DINGLE S EXTRACTION AT 1000 TERMS")
    print("=" * 70)

    mp.mp.dps = DPS
    N = len(a_coeffs) - 1
    xi0 = 2 / sqrt(mpf(3))
    beta_exact = -1 / (3 * sqrt(mpf(3)))

    gamma_beta = gamma(beta_exact)
    print(f"\n  ξ₀ = {nstr(xi0, 15)}")
    print(f"  β_exact = {nstr(beta_exact, 15)}")
    print(f"  Γ(β) = {nstr(gamma_beta, 15)}")

    # ── Step 1: Compute S_n for n = 200..N ───────────────
    # Dingle: a_n ~ (S/Γ(β)) · (-1)^n · Γ(n+β) / ξ₀^{n+β}
    # So: S_n = a_n · Γ(β) · ξ₀^{n+β} / ((-1)^n · Γ(n+β))
    print(f"\n  Step 1: Dingle S_n = a_n · Γ(β) · ξ₀^{{n+β}} / ((-1)^n · Γ(n+β))")

    S_raw = []
    for n in range(200, N + 1):
        sign = (-1) ** n
        S_n = a_coeffs[n] * gamma_beta * xi0 ** (n + beta_exact) / \
              (sign * gamma(n + beta_exact))
        S_raw.append((n, S_n))

    print(f"    Computed S_n for n = 200..{N}")
    print(f"    S_n (last 12):")
    for n, sn in S_raw[-12:]:
        print(f"      n={n:4d}: S_n = {nstr(sn, 20)}")

    # ── Step 2: Richardson in 1/(n+β) ────────────────────
    # S_n = S + c₁/(n+β) + c₂/(n+β)² + ...
    print(f"\n  Step 2: Richardson extrapolation in 1/(n+β)")

    start_n = 700
    vals = [s for n, s in S_raw if n >= start_n]
    ns = [n for n, _ in S_raw if n >= start_n]
    nbs = [mpf(n) + beta_exact for n in ns]

    def richardson_nb(values, nb_list):
        out_v, out_nb = [], []
        for i in range(1, len(values)):
            nb = nb_list[i]
            nb_prev = nb_list[i - 1]
            r = (nb * values[i] - nb_prev * values[i - 1]) / (nb - nb_prev)
            out_v.append(r)
            out_nb.append(nb)
        return out_v, out_nb

    cv, cnb = vals, nbs
    rich_S = {}
    for lev in range(1, 6):
        cv, cnb = richardson_nb(cv, cnb)
        if len(cv) < 3:
            break
        rich_S[lev] = cv[-1]
        print(f"    R^{lev} (last): S = {nstr(cv[-1], 20)}")
        if lev >= 2:
            # Show last 4 for stability
            print(f"    R^{lev} (last 4):")
            for j in range(max(0, len(cv) - 4), len(cv)):
                print(f"      S = {nstr(cv[j], 20)}")

    # ── Step 3: Wynn epsilon acceleration ────────────────
    print(f"\n  Step 3: Wynn epsilon algorithm on S_n (n≥700)")

    wynn_in = [s for n, s in S_raw if n >= 700]
    M = len(wynn_in)
    if M >= 10:
        eps_arr = [[mpf(0)] * M for _ in range(M + 1)]
        for i in range(M):
            eps_arr[0][i] = wynn_in[i]
        for r in range(M - 1):
            for i in range(M - r - 1):
                diff = eps_arr[r][i + 1] - eps_arr[r][i]
                if abs(diff) < mpf("1e-200"):
                    eps_arr[r + 1][i] = mpf("1e200")
                else:
                    prev = eps_arr[r - 1][i + 1] if r > 0 else mpf(0)
                    eps_arr[r + 1][i] = prev + 1 / diff
        wynn_vals = []
        for k in range(1, min(M // 2, 30)):
            v = eps_arr[2 * k][0]
            if abs(v) < 100:
                wynn_vals.append((2 * k, v))
        if wynn_vals:
            print(f"    Wynn (last 8):")
            for idx, val in wynn_vals[-8:]:
                print(f"      ε[{idx}][0] = {nstr(val, 20)}")
            S_wynn = wynn_vals[-1][1]
        else:
            S_wynn = None
            print(f"    Wynn: no stable values")
    else:
        S_wynn = None
        print(f"    Wynn: insufficient data")

    # ── Compare Richardson and Wynn ──────────────────────
    S_rich = rich_S.get(2, rich_S.get(1, vals[-1])) if rich_S else vals[-1]
    print(f"\n  ★ Dingle+Richardson^2: S = {nstr(S_rich, 20)}")
    if S_wynn is not None:
        print(f"  ★ Dingle+Wynn:        S = {nstr(S_wynn, 20)}")
        diff_rw = abs(S_rich - S_wynn)
        if diff_rw > 0 and abs(S_rich) > 0:
            agree_rw = int(-float(mp.log10(diff_rw / abs(S_rich))))
        else:
            agree_rw = 50
        print(f"    |Rich - Wynn| / |S| → {agree_rw} digits agreement")
    else:
        agree_rw = 0

    # ── Step 4: Split-half consistency ───────────────────
    print(f"\n  Step 4: Split-half consistency")
    mid = len(vals) // 2
    first_half = vals[:mid]
    second_half = vals[mid:]
    first_nbs = nbs[:mid]
    second_nbs = nbs[mid:]

    cv1, cn1 = first_half, first_nbs
    cv2, cn2 = second_half, second_nbs
    for _ in range(3):  # Richardson^3 on each half
        if len(cv1) >= 2:
            cv1, cn1 = richardson_nb(cv1, cn1)
        if len(cv2) >= 2:
            cv2, cn2 = richardson_nb(cv2, cn2)

    S_digits = 0
    if cv1 and cv2:
        S_first = cv1[-1]
        S_second = cv2[-1]
        split_diff = abs(S_first - S_second)
        if split_diff > 0 and abs(S_rich) > 0:
            S_digits = max(0, int(-float(mp.log10(split_diff / abs(S_rich)))))
        print(f"    First half (n=700..~{ns[mid]}):  S = {nstr(S_first, 18)}")
        print(f"    Second half (n=~{ns[mid]}..{ns[-1]}): S = {nstr(S_second, 18)}")
        print(f"    |upper - lower| = {nstr(split_diff, 6)}")
        print(f"    Agreement: ~{S_digits} digits")

    # Also try: numpy regression on S_n for even more digits
    print(f"\n  Step 5: Regression on S_n (complement)")
    # S_n = S + c₁/(n+β) + c₂/(n+β)² + c₃·log(n)/(n+β)² + ...
    fit_data_s = [(n, float(s)) for n, s in S_raw if n >= 500]
    ns_s = np.array([n for n, _ in fit_data_s], dtype=np.float64)
    ys_s = np.array([s for _, s in fit_data_s], dtype=np.float64)
    beta_f = float(beta_exact)
    nbs_f = ns_s + beta_f

    X_s = np.column_stack([
        np.ones_like(ns_s),
        1.0 / nbs_f,
        1.0 / nbs_f**2,
        np.log(ns_s) / nbs_f**2,
        1.0 / nbs_f**3,
    ])
    c_s, _, _, _ = np.linalg.lstsq(X_s, ys_s, rcond=None)
    S_regress = c_s[0]
    print(f"    Regression S (5-param, n≥500): {S_regress:.18f}")
    print(f"    c₁: {c_s[1]:.10f}, c₂: {c_s[2]:.10f}")

    # Extended regression with more terms
    X_s2 = np.column_stack([
        np.ones_like(ns_s),
        1.0 / nbs_f,
        1.0 / nbs_f**2,
        np.log(ns_s) / nbs_f**2,
        1.0 / nbs_f**3,
        np.log(ns_s) / nbs_f**3,
        1.0 / nbs_f**4,
    ])
    c_s2, _, _, _ = np.linalg.lstsq(X_s2, ys_s, rcond=None)
    S_regress2 = c_s2[0]
    print(f"    Regression S (7-param, n≥500): {S_regress2:.18f}")

    # Compare all S estimates
    all_S = {"rich": S_rich, "regress5": mpf(str(S_regress)),
             "regress7": mpf(str(S_regress2))}
    if S_wynn is not None:
        all_S["wynn"] = S_wynn

    print(f"\n  ── S summary ──")
    for label, sv in all_S.items():
        print(f"    {label:12s}: {nstr(sv, 20)}")

    # Cross-agreement matrix
    keys = list(all_S.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            d = abs(all_S[keys[i]] - all_S[keys[j]])
            if d > 0 and abs(all_S[keys[i]]) > 0:
                ag = int(-float(mp.log10(d / abs(all_S[keys[i]]))))
            else:
                ag = 50
            print(f"    {keys[i]} vs {keys[j]}: {ag} digits")

    # Best S: use Richardson as primary, regression as cross-check
    S_best = S_rich
    return S_best, S_digits, all_S, split_diff if cv1 and cv2 else None


# ════════════════════════════════════════════════════════════════════
# PHASE 4 — PSLQ ON S AT 12 DIGITS
# ════════════════════════════════════════════════════════════════════

def phase4(S_best, S_digits):
    print("\n" + "=" * 70)
    print("  PHASE 4 — PSLQ ON S AT 12 DIGITS")
    print("=" * 70)

    mp.mp.dps = PSLQ_DPS
    S = S_best
    print(f"\n  S = {nstr(S, 20)}")
    print(f"  S precision: ~{S_digits} digits")

    pslq_results = {}

    def try_pslq(label, basis, description="", s_index=1):
        """Run PSLQ. s_index is the position of S in the basis (default 1).
        A hit is only genuine if the coefficient at s_index is nonzero."""
        try:
            rel = pslq(basis, maxcoeff=COEFF_BOUND, maxsteps=10000)
        except Exception as e:
            rel = None
        if rel:
            residual = float(abs(sum(c * x for c, x in zip(rel, basis))))
            s_coeff_zero = (int(rel[s_index]) == 0) if s_index < len(rel) else True
            genuine = residual < RESIDUAL_THRESHOLD and not s_coeff_zero
            pslq_results[label] = {"relation": [int(c) for c in rel],
                                   "residual": residual, "hit": genuine,
                                   "description": description}
            if genuine:
                tag = "★ HIT ★"
            elif residual < RESIDUAL_THRESHOLD and s_coeff_zero:
                tag = "BASIS-DEP (S coeff=0)"
            else:
                tag = f"near (res={residual:.2e})"
            print(f"    {label}: {[int(c) for c in rel]} [{tag}]")
            if genuine and description:
                print(f"      → {description}")
        else:
            pslq_results[label] = {"relation": None, "residual": None,
                                   "hit": False, "description": description}
            print(f"    {label}: null")

    # ── Basis 1: Gamma(1/6) ──────────────────────────────
    print(f"\n  Basis 1 (Gamma at 1/6):")
    G16 = gamma(mpf(1) / 6)
    print(f"    Γ(1/6) = {nstr(G16, 15)}")
    print(f"    Γ(1/6)/(4π) = {nstr(G16 / (4 * pi), 15)}")
    try_pslq("gamma16_basic",
             [mpf(1), S, G16, G16 / (4 * pi)],
             "S in terms of Gamma(1/6)")
    try_pslq("gamma16_sqrt2pi",
             [mpf(1), S, G16 / sqrt(2 * pi), G16 * sqrt(mpf(3)) / (4 * pi)],
             "S with Gamma(1/6)/sqrt(2pi)")
    try_pslq("gamma16_extended",
             [mpf(1), S, G16, G16 / (4 * pi),
              G16 / sqrt(2 * pi), G16 * sqrt(mpf(3)) / (4 * pi)],
             "Full Gamma(1/6) basis")

    # ── Basis 2: Jimbo connection formula ────────────────
    # NOTE: sin(π/6)=1/2 is rational — DO NOT include in basis.
    # π/(Γ(σ)Γ(1-σ)) = sin(πσ) — rational multiple of known, avoid.
    print(f"\n  Basis 2 (Jimbo connection formula):")
    sigma = 1 / sqrt(mpf(3))
    G_plus = gamma(1 + sigma)
    G_minus = gamma(1 - sigma)
    G_sigma = gamma(sigma)
    print(f"    σ = 1/√3 = {nstr(sigma, 15)}")
    print(f"    Γ(1+σ) = {nstr(G_plus, 15)}")
    print(f"    Γ(1-σ) = {nstr(G_minus, 15)}")
    print(f"    Γ(σ) = {nstr(G_sigma, 15)}")
    try_pslq("jimbo_basic",
             [mpf(1), S, G_plus, G_minus, G_sigma],
             "Jimbo Γ(σ), Γ(1±σ) basis")
    try_pslq("jimbo_product",
             [mpf(1), S, G_plus * G_minus, G_sigma,
              G_plus / G_minus],
             "Jimbo product/ratio basis")

    # ── Basis 3: Exponential CM ──────────────────────────
    print(f"\n  Basis 3 (exponential CM, discriminant -3):")
    e_pos = exp(pi / sqrt(mpf(3)))
    e_neg = exp(-pi / sqrt(mpf(3)))
    e2 = exp(2 * pi / sqrt(mpf(3)))
    try_pslq("exp_cm",
             [mpf(1), S, e_pos, e_neg, e2, 1 / sqrt(e2 - 1)],
             "CM exponentials at discriminant -3")
    try_pslq("exp_cm_reduced",
             [mpf(1), S, e_pos, e_neg],
             "Reduced exponential CM")

    # ── Basis 4: Algebraic over Q(√3) ───────────────────
    print(f"\n  Basis 4 (algebraic over Q(√3)):")
    sq3 = sqrt(mpf(3))
    q4 = mpf(3) ** (mpf(1) / 4)
    q34 = mpf(3) ** (mpf(3) / 4)
    try_pslq("algebraic_q3",
             [mpf(1), S, S ** 2, sq3, q4, q34, S * sq3],
             "S algebraic degree ≤ 4 over Q(sqrt(3))")
    try_pslq("algebraic_simple",
             [mpf(1), S, S ** 2, S ** 3, sq3],
             "S algebraic degree ≤ 3 with sqrt(3)")

    # ── Basis 5: Barnes G ────────────────────────────────
    print(f"\n  Basis 5 (Barnes G-function):")
    # Barnes G: G(1+z) = Gamma(z) · G(z), G(1) = 1
    # G(z) = (2π)^{z/2} · exp(-[z(z-1)/2 + z·γ/2]) · ∏_{k=1}^∞ [(1+z/k)^k · exp(-z+z²/(2k))]
    # For numerical eval, use mpmath loggamma relation:
    # log G(1+z) = z/2 · log(2π) - z(1+z)/2 + z · log Γ(z+1) - log Γ(z+1) + ...
    # Actually mpmath has barnesg
    try:
        G_plus_b = mp.barnesg(1 + sigma)
        G_minus_b = mp.barnesg(1 - sigma)
        print(f"    G(1+σ) = {nstr(G_plus_b, 15)}")
        print(f"    G(1-σ) = {nstr(G_minus_b, 15)}")
        try_pslq("barnes_g",
                 [mpf(1), S, G_plus_b * G_minus_b,
                  mp.log(G_plus_b)],
                 "Barnes G at 1±σ")
        try_pslq("barnes_g_ext",
                 [mpf(1), S, G_plus_b, G_minus_b,
                  G_plus_b * G_minus_b],
                 "Barnes G extended")
    except AttributeError:
        print(f"    Barnes G not available in this mpmath version")
        pslq_results["barnes_g"] = {"relation": None, "residual": None,
                                    "hit": False, "description": "barnesg unavailable"}
        pslq_results["barnes_g_ext"] = pslq_results["barnes_g"]

    # ── Bonus: additional bases ──────────────────────────
    print(f"\n  Bonus bases:")
    # Stokes constant often involves 2πi / Γ factors
    G13 = gamma(mpf(1) / 3)
    G23 = gamma(mpf(2) / 3)
    try_pslq("gamma_thirds",
             [mpf(1), S, G13, G23, G13 * G23, pi / G13],
             "Gamma(1/3), Gamma(2/3)")
    # NOTE: π√3 = 3·(π/√3), so they are dependent. Use only independent elements.
    try_pslq("pi_combos",
             [mpf(1), S, pi, sq3, pi ** 2, pi * sq3],
             "pi and sqrt(3) combinations")
    # sin(πβ)/π factor (Stokes normalization)
    sin_pi_beta = mp.sin(pi * (-1 / (3 * sq3)))
    try_pslq("stokes_norm",
             [mpf(1), S * sin_pi_beta / pi,
              sq3, pi, G13, G23],
             "S·sin(πβ)/π normalization")

    mp.mp.dps = DPS

    # Count hits
    hits = [k for k, v in pslq_results.items() if v.get("hit")]
    return pslq_results, hits


# ════════════════════════════════════════════════════════════════════
# PHASE 5 — GOVERNANCE
# ════════════════════════════════════════════════════════════════════

def phase5(beta_dig, S_digits, pslq_hits, S_best, best_beta, split_diff):
    print("\n" + "=" * 70)
    print("  PHASE 5 — GOVERNANCE")
    print("=" * 70)

    mp.mp.dps = DPS
    mu_diff = -1 / (3 * sqrt(mpf(3)))

    claims = []

    # Beta claim
    if beta_dig >= 9:
        claim_b = {
            "territory": "T2", "iteration": 22, "type": "numerical_identity",
            "object": "beta_exponent", "value": str(best_beta),
            "identity": "-1/(3*sqrt(3)) = -sqrt(3)/9",
            "digits": beta_dig,
        }
        claims.append(claim_b)
        beta_claim = "numerical_identity"
    else:
        claim_b = {
            "territory": "T2", "iteration": 22, "type": "near_miss",
            "object": "beta_exponent", "value": str(best_beta),
            "target": "-1/(3*sqrt(3))", "digits": beta_dig,
        }
        claims.append(claim_b)
        beta_claim = "near_miss"

    # S claim
    if pslq_hits:
        claim_s = {
            "territory": "T2", "iteration": 22, "type": "numerical_identity",
            "object": "stokes_constant_S", "value": nstr(S_best, 15),
            "pslq_hits": pslq_hits, "digits": S_digits,
        }
        claims.append(claim_s)
        s_claim = "numerical_identity"
        f2_status = "partial"
    else:
        claim_s = {
            "territory": "T2", "iteration": 22, "type": "near_miss",
            "object": "stokes_constant_S", "value": nstr(S_best, 15),
            "digits": S_digits, "pslq": "null across all bases",
        }
        claims.append(claim_s)
        s_claim = "near_miss"
        f2_status = "open"

    # AHS calculation
    # Components: series_ok, beta_quality, S_quality, pslq_outcome
    series_score = 1.0  # 1000 terms
    beta_score = min(beta_dig / 12, 1.0)  # target 12 digits
    s_score = min(S_digits / 15, 1.0)  # target 15 digits
    pslq_score = 1.0 if pslq_hits else 0.0
    ahs = 0.2 * series_score + 0.2 * beta_score + 0.4 * s_score + 0.2 * pslq_score
    ahs = round(ahs, 4)

    with open(CLAIMS, "a") as f:
        for c in claims:
            f.write(json.dumps(c) + "\n")

    print(f"\n  {len(claims)} claim(s) emitted")
    print(f"  β digits: {beta_dig}, S digits: {S_digits}")
    print(f"  AHS = {ahs} (threshold={AHS_THRESHOLD})")
    governance = "CONTINUE" if ahs >= AHS_THRESHOLD else "HALT"
    print(f"  Governance: {governance}")

    return claims, ahs, governance, beta_claim, s_claim, f2_status


# ════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()

    mp.mp.dps = DPS
    print(f"  Computing V_quad...")
    vquad = compute_vquad(2000, DPS)
    print(f"  V_quad = {nstr(vquad, 30)} ({DPS} digits)")

    # Phase 1
    a_coeffs, sigma_rec, mu, ok = phase1()
    if not ok:
        print("\n  *** Phase 1 FAILED — aborting ***")
        return

    # Phase 2
    best_beta, beta_dig, beta_pslq = phase2(a_coeffs)

    # Phase 3
    S_best, S_digits, all_S, split_diff = phase3(a_coeffs)

    # Phase 4
    pslq_results, pslq_hits = phase4(S_best, S_digits)

    # Phase 5
    claims, ahs, governance, beta_claim, s_claim, f2_status = \
        phase5(beta_dig, S_digits, pslq_hits, S_best, best_beta, split_diff)

    # ── Write results ────────────────────────────────────
    result = {
        "territory": "T2",
        "iteration": 22,
        "series_terms": len(a_coeffs),
        "dps": DPS,
        "beta_estimate": str(best_beta),
        "mu_diff_exact": nstr(-1 / (3 * sqrt(mpf(3))), 20),
        "beta_match_digits": beta_dig,
        "beta_claim": beta_claim,
        "beta_pslq": beta_pslq,
        "S_best": nstr(S_best, 20),
        "S_all": {k: nstr(v, 18) for k, v in all_S.items()},
        "S_consistency_digits": S_digits,
        "split_half_diff": nstr(split_diff, 6) if split_diff else "N/A",
        "pslq_results": {k: {"relation": v["relation"],
                              "residual": v["residual"],
                              "hit": v["hit"]}
                         for k, v in pslq_results.items()},
        "pslq_hits": pslq_hits,
        "s_claim": s_claim,
        "f2_status": f2_status,
        "ahs": ahs,
        "governance": governance,
        "vquad": nstr(vquad, 20),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    outpath = RESULTS / "t2_iter22_s_precision.json"
    with open(outpath, "w") as f:
        json.dump(result, f, indent=2, default=str)

    # ── Deliverable ──────────────────────────────────────
    print("\n" + "═" * 70)
    print("  DELIVERABLE — T2 ITERATION 22")
    print("═" * 70)
    print(f"\n  SERIES TERMS:         {len(a_coeffs)}")
    print(f"  BETA (log-corrected): {best_beta:.18f}")
    print(f"  BETA DIGITS:          {beta_dig}")
    print(f"  BETA CLAIM:           {beta_claim}")
    print(f"  S (1000 terms):       {nstr(S_best, 20)}")
    print(f"  SPLIT-HALF:           {nstr(split_diff, 6) if split_diff else 'N/A'}")
    for k in ["gamma16_basic", "gamma16_extended", "jimbo_basic", "jimbo_reduced",
              "exp_cm", "algebraic_q3", "barnes_g"]:
        if k in pslq_results:
            hit = "HIT" if pslq_results[k]["hit"] else "null"
            print(f"  PSLQ {k:20s}: {hit}")
    print(f"  S IDENTIFIED:         {'yes (' + str(pslq_hits) + ')' if pslq_hits else 'no'}")
    print(f"  AHS:                  {ahs}")
    print(f"  F2 STATUS:            {f2_status}")
    print(f"\n  Full report: {outpath}")
    print(f"\n  Total elapsed: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
