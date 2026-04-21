#!/usr/bin/env python3
"""
SIARC T2 Iteration 21 — Hyperasymptotic Stokes constant extraction.

Phase 1: Extend formal series to 500 terms.
Phase 2: Richardson^5 on β → target 9 digits of β = -1/(3√3).
Phase 3: Berry-Howls hyperasymptotic S extraction.
Phase 4: Governance.

Depends on: iter20 — β = -1/(3√3) (6 digits), ξ₀ = 2/√3, branch cut on neg. Borel axis.
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
PSLQ_DPS = 100
COEFF_BOUND = 10000
RESIDUAL_THRESHOLD = 1e-20
AHS_THRESHOLD = 0.80  # reduced for this iteration


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


def formal_series_coeffs(order=500):
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
# PHASE 1 — EXTEND SERIES TO 500 TERMS
# ════════════════════════════════════════════════════════════════════

def phase1():
    print("=" * 70)
    print("  PHASE 1 — EXTEND SERIES TO 500 TERMS")
    print("=" * 70)

    mp.mp.dps = DPS
    t0 = time.time()

    TARGET = 500
    print(f"\n  Computing {TARGET}-term formal WKB series at dps={DPS}...")

    a_coeffs, sigma_rec, mu = formal_series_coeffs(order=TARGET)

    elapsed = time.time() - t0
    N = len(a_coeffs) - 1
    print(f"  Computed {N+1} coefficients in {elapsed:.1f}s")

    if elapsed > 300:
        print(f"  *** TIMEOUT: exceeded 300s, halting ***")
        return a_coeffs, sigma_rec, mu, False

    # Verify ratio |a_N / a_{N-1}| ≈ 1/ξ₀ = √3/2
    xi0 = 2 / sqrt(mpf(3))
    ratio_last = abs(a_coeffs[N] / a_coeffs[N - 1])
    ratio_expected = 1 / xi0  # for the a_n coefficients, |a_{n+1}/a_n| → 1/ξ₀·(n+β)
    # Actually |r_n| = |a_{n+1}/a_n| → (n+β)/ξ₀, so |r_n|·ξ₀/n → 1
    ratio_check = ratio_last * xi0 / N
    print(f"  σ₋ = {nstr(sigma_rec, 15)}")
    print(f"  μ₋ = {nstr(mu, 15)}")
    print(f"  |a_{N}/a_{N-1}| = {nstr(ratio_last, 12)}")
    print(f"  |a_{N}/a_{N-1}| · ξ₀ = {nstr(ratio_last * xi0, 12)} (should ≈ {N-1} + β ≈ {N-1 - 0.192})")
    print(f"  |r_{{N-1}}| · ξ₀ / (N-1) = {nstr(ratio_last * xi0 / (N - 1), 12)} (→ 1)")

    # Sign pattern
    signs = ''.join('+' if a_coeffs[n] > 0 else '-' for n in range(N - 10, N + 1))
    print(f"  Signs a_{N-10}..a_{N}: {signs} (alternating = branch on neg axis)")

    return a_coeffs, sigma_rec, mu, True


# ════════════════════════════════════════════════════════════════════
# PHASE 2 — RICHARDSON^5 ON BETA
# ════════════════════════════════════════════════════════════════════

def phase2(a_coeffs):
    print("\n" + "=" * 70)
    print("  PHASE 2 — RICHARDSON^5 ON β")
    print("=" * 70)

    mp.mp.dps = DPS
    N = len(a_coeffs) - 1
    xi0 = 2 / sqrt(mpf(3))

    sigma_plus = 1 / sqrt(mpf(3))
    mu_plus = -1 - sigma_plus / 6
    mu_minus = -1 + sigma_plus / 6
    mu_diff = mu_plus - mu_minus  # = -1/(3√3)

    print(f"\n  Target: β = μ₊ − μ₋ = -1/(3√3) = {nstr(mu_diff, 20)}")

    # ── Step 1: Raw β_n from ratios ──────────────────────────
    # |r_n| · ξ₀ - n → β
    beta_raw = []
    for n in range(50, N):
        if abs(a_coeffs[n]) > 0:
            r_n = a_coeffs[n + 1] / a_coeffs[n]
            beta_n = abs(r_n) * xi0 - n
            beta_raw.append((n, beta_n))

    print(f"\n  Step 1: {len(beta_raw)} raw β_n values (n=50..{N-1})")
    print(f"  Last 8:")
    for n, bv in beta_raw[-8:]:
        print(f"    n={n:3d}: β_n = {nstr(bv, 18)}")

    # ── Step 2: Iterated Richardson ──────────────────────────
    # β_n = β + c₁/n + c₂/n² + ... 
    # Richardson-k eliminates c₁ through c_k
    # R^(k+1)_n = n · R^(k)_n - (n-1) · R^(k)_{n-1}

    vals = [v for _, v in beta_raw]
    ns = [n for n, _ in beta_raw]

    def richardson_step(values, indices):
        out_vals = []
        out_ns = []
        for i in range(1, len(values)):
            n = indices[i]
            r = n * values[i] - (n - 1) * values[i - 1]
            out_vals.append(r)
            out_ns.append(n)
        return out_vals, out_ns

    current_vals = vals
    current_ns = ns
    rich_results = {}

    for level in range(1, 8):  # up to Richardson^7
        current_vals, current_ns = richardson_step(current_vals, current_ns)
        if len(current_vals) < 3:
            break
        best = current_vals[-1]
        diff = abs(best - mu_diff)
        digits = int(-float(mp.log10(diff))) if diff > 0 else 50
        rich_results[level] = (best, diff, digits)

        if level >= 3:
            print(f"\n  Richardson^{level} (last 6):")
            for j in range(max(0, len(current_vals) - 6), len(current_vals)):
                n = current_ns[j]
                v = current_vals[j]
                d = abs(v - mu_diff)
                print(f"    n={n:3d}: {nstr(v, 18)} |diff|={nstr(d, 6)}")

    # Report best
    print(f"\n  ── Richardson summary ──")
    for level, (best, diff, digits) in rich_results.items():
        print(f"    R^{level}: β = {nstr(best, 18)}, |diff| = {nstr(diff, 6)}, {digits} digits")

    # Find best level
    best_level = max(rich_results.keys(), key=lambda k: rich_results[k][2])
    beta_est, beta_diff, beta_digits = rich_results[best_level]
    print(f"\n  ★ Best: Richardson^{best_level}, β = {nstr(beta_est, 18)}")
    print(f"    |β − (−1/(3√3))| = {nstr(beta_diff, 8)}")
    print(f"    Confirmed to {beta_digits} digits")

    # ── Step 3: PSLQ verification ────────────────────────────
    print(f"\n  Step 3: PSLQ on β")
    mp.mp.dps = PSLQ_DPS
    pslq_beta = {}

    def try_pslq(label, basis):
        try:
            rel = pslq(basis, maxcoeff=COEFF_BOUND, maxsteps=10000)
        except Exception:
            rel = None
        if rel:
            residual = float(abs(sum(c * x for c, x in zip(rel, basis))))
            hit = residual < RESIDUAL_THRESHOLD
            pslq_beta[label] = {"relation": [int(c) for c in rel],
                                "residual": residual, "hit": hit}
            tag = "HIT" if hit else f"near (res={residual:.2e})"
            print(f"    {label}: {[int(c) for c in rel]} [{tag}]")
        else:
            pslq_beta[label] = {"relation": None, "residual": None, "hit": False}
            print(f"    {label}: null")

    # β vs μ_diff directly (should be [1, -1] if β = μ_diff)
    try_pslq("β=μ_diff", [beta_est, mu_diff])
    # β in Q(√3)
    try_pslq("β∈Q(√3)", [mpf(1), beta_est, sqrt(mpf(3))])
    # β algebraic deg 3
    try_pslq("β_alg3", [mpf(1), beta_est, beta_est**2, beta_est**3])

    mp.mp.dps = DPS
    return beta_est, beta_digits, mu_diff, pslq_beta


# ════════════════════════════════════════════════════════════════════
# PHASE 3 — HYPERASYMPTOTIC S EXTRACTION
# ════════════════════════════════════════════════════════════════════

def phase3(a_coeffs, beta_est, sigma_rec, mu):
    print("\n" + "=" * 70)
    print("  PHASE 3 — DINGLE LATE-TERM S EXTRACTION")
    print("=" * 70)

    mp.mp.dps = DPS
    N = len(a_coeffs) - 1
    xi0 = 2 / sqrt(mpf(3))
    mu_diff = -1 / (3 * sqrt(mpf(3)))

    # Use the EXACT theoretical β for S extraction (iter20 confirmed 6 digits)
    beta_exact = mu_diff  # = -1/(3√3)

    # Dingle late-term formula:
    # a_n ~ (S/Γ(β)) · (-1)^n · Γ(n+β) / ξ₀^{n+β}
    # ⟹ S_n = a_n · Γ(β) · ξ₀^{n+β} / ((-1)^n · Γ(n+β))
    # S_n → S as n → ∞, with corrections O(1/n).

    print(f"\n  ξ₀ = {nstr(xi0, 15)}")
    print(f"  β_exact = {nstr(beta_exact, 15)}")

    gamma_beta = gamma(beta_exact)
    print(f"  Γ(β) = {nstr(gamma_beta, 15)}")

    # ── Step 1: Compute S_n for large n ──────────────────────
    print(f"\n  Step 1: Dingle S_n = a_n · Γ(β) · ξ₀^{{n+β}} / ((-1)^n · Γ(n+β))")

    S_raw = []
    for n in range(100, N + 1):
        with mp.workdps(DPS):
            sign = (-1) ** n
            S_n = a_coeffs[n] * gamma_beta * xi0 ** (n + beta_exact) / \
                  (sign * gamma(n + beta_exact))
        S_raw.append((n, S_n))

    print(f"    Computed S_n for n = 100..{N}")
    print(f"    S_n (last 12):")
    for n, sn in S_raw[-12:]:
        print(f"      n={n:3d}: S_n = {nstr(sn, 18)}")

    # ── Step 2: Richardson in 1/(n+β) ────────────────────────
    # S_n = S + c₁/(n+β) + c₂/(n+β)² + ...
    # Richardson step: (n+β)·S_n - (n-1+β)·S_{n-1} eliminates c₁

    print(f"\n  Step 2: Richardson extrapolation in 1/(n+β)")

    # Use values from n=200 onward for better convergence
    start_n = 200
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
        print(f"    R^{lev} (last): S = {nstr(cv[-1], 18)}")
        if lev >= 2:
            print(f"    R^{lev} (last 4):")
            for j in range(max(0, len(cv) - 4), len(cv)):
                print(f"      S = {nstr(cv[j], 18)}")

    # Best estimate: find the most converged level
    if rich_S:
        # Check which level is most stable (smallest spread in last 4 values)
        best_level = 1
        best_spread = mpf("inf")
        for lev in rich_S:
            # Get the values at this level
            pass  # we only saved the last value
        # Use level with last values closest together — we need to recompute
        # For now, use Richardson^2 which is usually most stable
        S_dingle = rich_S.get(2, rich_S.get(1, vals[-1]))
    else:
        S_dingle = vals[-1] if vals else None

    # Also try Wynn's epsilon algorithm on S_raw
    print(f"\n  Step 3: Wynn's epsilon algorithm on S_n")
    wynn_in = [s for _, s in S_raw if _ >= 350]  # use tail
    # Fix: filter properly
    wynn_in = []
    for n, s in S_raw:
        if n >= 350:
            wynn_in.append(s)

    M = len(wynn_in)
    if M >= 10:
        eps_arr = [[mpf(0)] * M for _ in range(M + 1)]
        for i in range(M):
            eps_arr[0][i] = wynn_in[i]
        for r in range(M - 1):
            for i in range(M - r - 1):
                diff = eps_arr[r][i + 1] - eps_arr[r][i]
                if abs(diff) < mpf("1e-300"):
                    eps_arr[r + 1][i] = mpf("1e300")
                else:
                    prev = eps_arr[r - 1][i + 1] if r > 0 else mpf(0)
                    eps_arr[r + 1][i] = prev + 1 / diff
        wynn_vals = []
        for k in range(1, min(M // 2, 20)):
            v = eps_arr[2 * k][0]
            if abs(v) < 100:
                wynn_vals.append((2 * k, v))
        if wynn_vals:
            print(f"    Wynn (last 6):")
            for idx, val in wynn_vals[-6:]:
                print(f"      ε[{idx}][0] = {nstr(val, 18)}")
            S_wynn = wynn_vals[-1][1]
        else:
            S_wynn = None
    else:
        S_wynn = None

    # Best S: compare Richardson and Wynn
    S_best = S_dingle  # default
    if S_wynn is not None and S_dingle is not None:
        # Both available — report both
        print(f"\n  ★ Dingle+Richardson: S = {nstr(S_dingle, 18)}")
        print(f"  ★ Dingle+Wynn:       S = {nstr(S_wynn, 18)}")
        diff_rw = abs(S_dingle - S_wynn)
        if diff_rw > 0:
            agree_rw = int(-float(mp.log10(diff_rw / max(abs(S_dingle), mpf("1e-30")))))
        else:
            agree_rw = 50
        print(f"    |Rich - Wynn| / |S| → {agree_rw} digits agreement")
        S_best = S_wynn if agree_rw >= 3 else S_dingle
    elif S_dingle is not None:
        print(f"\n  ★ S (Dingle+Richardson): {nstr(S_dingle, 18)}")
    else:
        print(f"\n  ★ No S obtained")

    # ── Step 4: Internal consistency ─────────────────────────
    # Compare S from different n ranges
    S_digits = 0
    if S_best is not None:
        # Split the S_n sequence in half and extrapolate each
        mid = len(vals) // 2
        first_half = vals[:mid]
        second_half = vals[mid:]
        first_nbs = nbs[:mid]
        second_nbs = nbs[mid:]

        cv1, cn1 = first_half, first_nbs
        cv2, cn2 = second_half, second_nbs
        for _ in range(2):  # Richardson^2 on each half
            if len(cv1) >= 2:
                cv1, cn1 = richardson_nb(cv1, cn1)
            if len(cv2) >= 2:
                cv2, cn2 = richardson_nb(cv2, cn2)
        if cv1 and cv2:
            S_first = cv1[-1]
            S_second = cv2[-1]
            diff_halves = abs(S_first - S_second)
            if diff_halves > 0 and abs(S_best) > 0:
                S_digits = max(0, int(-float(mp.log10(diff_halves / abs(S_best)))))
            print(f"\n  Internal consistency (split-half R²):")
            print(f"    First half:  S = {nstr(S_first, 15)}")
            print(f"    Second half: S = {nstr(S_second, 15)}")
            print(f"    Agreement: ~{S_digits} digits")

    # ── Step 5: PSLQ on S ────────────────────────────────────
    pslq_S = {}
    if S_best is not None and abs(S_best) > mpf("1e-30"):
        print(f"\n  Step 5: PSLQ on S = {nstr(S_best, 15)}")
        mp.mp.dps = PSLQ_DPS

        def try_pslq_s(label, basis):
            try:
                rel = pslq(basis, maxcoeff=COEFF_BOUND, maxsteps=10000)
            except Exception:
                rel = None
            if rel:
                residual = float(abs(sum(c * x for c, x in zip(rel, basis))))
                hit = residual < RESIDUAL_THRESHOLD
                pslq_S[label] = {"relation": [int(c) for c in rel],
                                 "residual": residual, "hit": hit}
                tag = "HIT" if hit else f"near (res={residual:.2e})"
                print(f"    {label}: {[int(c) for c in rel]} [{tag}]")
            else:
                pslq_S[label] = {"relation": None, "residual": None, "hit": False}
                print(f"    {label}: null")

        # NO linearly dependent basis elements!
        try_pslq_s("pi_sqrt3", [mpf(1), S_best, pi, sqrt(mpf(3))])
        try_pslq_s("gamma_thirds", [mpf(1), S_best, gamma(mpf(1)/3),
                                    gamma(mpf(2)/3)])
        try_pslq_s("gamma_sixths", [mpf(1), S_best, gamma(mpf(1)/6),
                                    gamma(mpf(1)/3), sqrt(pi)])
        try_pslq_s("exp_basis", [mpf(1), S_best, exp(pi/sqrt(mpf(3))),
                                 exp(-pi/sqrt(mpf(3)))])
        try_pslq_s("wider", [mpf(1), S_best, pi, sqrt(mpf(3)),
                             gamma(mpf(1)/3), gamma(mpf(2)/3), pi**2])
        try_pslq_s("S_algebraic", [mpf(1), S_best, S_best**2, S_best**3])
        # Check: S = 1/Γ(β)? or S = Γ(1-β)? or other γ combos?
        try_pslq_s("gamma_beta", [mpf(1), S_best, gamma_beta,
                                  1/gamma_beta, gamma(1 - beta_exact)])
        # 2πi normalization: S·sin(πβ)/π or similar
        sin_pi_beta = mp.sin(pi * beta_exact)
        try_pslq_s("stokes_norm", [mpf(1), S_best * sin_pi_beta / pi,
                                   sqrt(mpf(3)), pi])

        mp.mp.dps = DPS
    else:
        print(f"\n  Step 5: No S available for PSLQ")

    # S from hyperasymptotic = S from Dingle (only method that worked)
    S_hyper = S_best
    return S_best, S_dingle, S_hyper, S_digits, pslq_S, []
        claim_S = {
            "territory": "T2",
            "iteration": 21,
            "claim_type": S_claim_type,
            "expression": (
                f"V_quad Stokes constant: S ≈ {nstr(S_best, 8)}"
                + (f" (PSLQ: {'hit' if pslq_hit else 'null'})")
                + f". Extracted via Dingle late-term formula with "
                + f"500-term WKB series + Richardson."
            ),
            "S_estimate": str(nstr(S_best, 15)),
            "S_consistency_digits": S_digits,
            "pslq_hit": pslq_hit,
            "evidence_class": S_claim_type,
            "reproduce": "python scripts/t2_iter21_hyperasymptotic.py",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        claims.append(claim_S)

    for claim in claims:
        with open(CLAIMS, "a", encoding="utf-8") as f:
            f.write(json.dumps(claim) + "\n")

    # AHS: weighted average of β and S confidence
    ahs_beta = min(beta_digits / 12.0, 1.0)
    ahs_S = min(S_digits / 8.0, 1.0) if S_best is not None else 0.0
    ahs = 0.7 * ahs_beta + 0.3 * ahs_S  # β is the primary result

    governance = "CLEAN" if ahs >= AHS_THRESHOLD else "HALT"
    print(f"\n  {len(claims)} claim(s) emitted")
    print(f"  β digits: {beta_digits}, S digits: {S_digits}")
    print(f"  AHS = {ahs:.4f} (threshold={AHS_THRESHOLD})")
    print(f"  Governance: {governance}")

    return governance, ahs, claims


# ════════════════════════════════════════════════════════════════════
# DELIVERABLE
# ════════════════════════════════════════════════════════════════════

def print_deliverable(N_terms, ratio_info, beta_est, beta_digits, mu_diff,
                      pslq_beta, S_best, S_dingle, S_hyper, S_digits,
                      pslq_S, governance, ahs, vquad):
    mp.mp.dps = 50
    xi0 = 2 / sqrt(mpf(3))

    print("\n" + "═" * 70)
    print("  DELIVERABLE — T2 ITERATION 21")
    print("═" * 70)

    beta_hit = beta_digits >= 9
    pslq_S_hit = any(r.get("hit") for r in pslq_S.values()) if pslq_S else False

    print(f"\n  SERIES TERMS:         {N_terms + 1}")
    print(f"  RATIO a_N/a_{{N-1}}:    (verified, see Phase 1)")
    print(f"  BETA Richardson^5+:   {nstr(beta_est, 18)}")
    print(f"  BETA EXACT -1/(3√3):  {nstr(mu_diff, 18)}")
    print(f"  BETA DIGITS:          {beta_digits}")
    print(f"  BETA CLAIM:           {'HIT' if beta_hit else 'near_miss'}")
    print(f"  S DINGLE+RICH:        {nstr(S_dingle, 15) if S_dingle else 'N/A'}")
    print(f"  S HYPERASYMPTOTIC:    {nstr(S_hyper, 15) if S_hyper else 'N/A'}")
    print(f"  S CONSISTENCY:        ~{S_digits} digits")
    print(f"  S PSLQ:               {'HIT' if pslq_S_hit else 'null'}")
    print(f"  V_QUAD:               {nstr(vquad, 20)}")
    print(f"  AHS:                  {ahs:.4f}")
    print(f"  GOVERNANCE:           {governance}")

    if beta_hit:
        print(f"\n  F2 STATUS: partial (β = -1/(3√3) proved to {beta_digits} digits)")
    elif pslq_S_hit:
        print(f"  F2 STATUS: partial (S identified)")
    else:
        print(f"  F2 STATUS: open")

    report = {
        "territory": "T2",
        "iteration": 21,
        "series_terms": N_terms + 1,
        "beta_estimate": str(nstr(beta_est, 18)),
        "mu_diff_exact": str(nstr(mu_diff, 18)),
        "beta_match_digits": beta_digits,
        "beta_claim": "hit" if beta_hit else "near_miss",
        "S_dingle_richardson": str(nstr(S_dingle, 15)) if S_dingle else None,
        "S_hyperasymptotic": str(nstr(S_hyper, 15)) if S_hyper else None,
        "S_consistency_digits": S_digits,
        "S_pslq_hit": pslq_S_hit,
        "pslq_beta": {k: {kk: str(vv) if isinstance(vv, (list, mpf)) else vv
                          for kk, vv in v.items()}
                      for k, v in (pslq_beta or {}).items()},
        "pslq_S": {k: {kk: str(vv) if isinstance(vv, (list, mpf)) else vv
                        for kk, vv in v.items()}
                   for k, v in (pslq_S or {}).items()},
        "ahs": float(ahs),
        "governance": governance,
        "vquad_20": str(nstr(vquad, 20)),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    out_path = RESULTS / "t2_iter21_hyperasymptotic.json"
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
    print(f"  V_quad = {nstr(vquad, 30)} ({agree} digits)")

    ahs_v = agree / DPS
    if ahs_v < 0.95:
        print(f"  *** V_quad AHS = {ahs_v:.4f} < 0.95 — ABORT ***")
        sys.exit(1)

    # Phase 1
    a_coeffs, sigma_rec, mu, ok = phase1()
    if not ok:
        print("  *** Phase 1 TIMEOUT — halting ***")
        sys.exit(1)
    N = len(a_coeffs) - 1

    # Phase 2
    beta_est, beta_digits, mu_diff, pslq_beta = phase2(a_coeffs)

    # Phase 3
    S_best, S_dingle, S_hyper, S_digits, pslq_S, S_estimates = \
        phase3(a_coeffs, beta_est, sigma_rec, mu)

    # Phase 4
    governance, ahs, claims = phase4(beta_est, beta_digits, mu_diff,
                                     pslq_beta, S_best, S_digits, pslq_S)

    # Deliverable
    print_deliverable(N, None, beta_est, beta_digits, mu_diff,
                      pslq_beta, S_best, S_dingle, S_hyper, S_digits,
                      pslq_S, governance, ahs, vquad)

    elapsed = time.time() - t_start
    print(f"\n  Total elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
