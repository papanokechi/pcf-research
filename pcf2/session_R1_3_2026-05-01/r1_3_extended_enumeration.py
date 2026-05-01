"""PCF-2 SESSION R1.3 -- Phase R13-C: extended quartic enumeration

Window:
    a4 in {1, 2, 3, 5, 7}
    a3, a2, a1, a0 in {-W, ..., W} for W = 5 (initial),
                                escalate to W = 7 then W = 10 if j=0 cell
                                size < 4.

For each candidate b(n) = a4 n^4 + a3 n^3 + a2 n^2 + a1 n + a0 we check:
    - gcd(a4,...,a0) == 1                              (Z-primitive)
    - I_quartic(coeffs) = 12*a4*a0 - 3*a3*a1 + a2^2 == 0
                                                       (j(Jac(C_b)) = 0)
    - irreducibility over Q
    - non-singular (Delta_quartic = 4 I^3 - J^2 != 0)
    - b is squarefree as a polynomial in n
The j=0 cell is harvested.  We also count (without expensive irred test)
the size of the broader irreducible catalogue at each window for
context.

For each NEW j=0 quartic family (not in Q1's 60), we run the
fixed-A=8 tail-window WKB at dps=2000, N=[50,250], N_ref=400 (matching
Phase R13-B), and record A_R13_free, delta_R13_free, residual_mean,
residual_at_max_n, residual_std.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
SES_Q1 = ROOT / "sessions" / "2026-05-01" / "PCF2-SESSION-Q1"

DPS = 2000
N_GRID = list(range(50, 252, 2))
N_REF = 400
LOG = HERE / "r1_3_extended_enumeration.log"

n_sym = sp.symbols("n")


def L(s):
    line = f"[{time.strftime('%H:%M:%S')}] {s}"
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


# --- enumeration helpers (mirrors Q1 conventions) -----------------

def igusa_I(a4, a3, a2, a1, a0):
    return 12 * a4 * a0 - 3 * a3 * a1 + a2 ** 2


def igusa_J(a4, a3, a2, a1, a0):
    return (72 * a4 * a2 * a0 - 27 * a4 * a1 ** 2 - 27 * a3 ** 2 * a0
            + 9 * a3 * a2 * a1 - 2 * a2 ** 3)


def primitive(coeffs):
    g = 0
    for c in coeffs:
        g = math.gcd(g, abs(c))
    return g == 1


def is_irreducible_quartic(coeffs):
    a4, a3, a2, a1, a0 = coeffs
    poly = a4 * n_sym ** 4 + a3 * n_sym ** 3 + a2 * n_sym ** 2 + a1 * n_sym + a0
    facs = sp.factor_list(poly, n_sym, domain="QQ")[1]
    return len(facs) == 1 and facs[0][1] == 1 and sp.degree(facs[0][0], n_sym) == 4


def is_squarefree_poly(coeffs):
    a4, a3, a2, a1, a0 = coeffs
    p = sp.Poly(a4 * n_sym ** 4 + a3 * n_sym ** 3 + a2 * n_sym ** 2
                + a1 * n_sym + a0, n_sym)
    return sp.gcd(p, p.diff(n_sym)).degree() == 0


def Delta_4(coeffs):
    a4, a3, a2, a1, a0 = coeffs
    return int(sp.discriminant(a4 * n_sym ** 4 + a3 * n_sym ** 3
                               + a2 * n_sym ** 2 + a1 * n_sym + a0, n_sym))


def galois_quartic_quick(coeffs):
    """Light-weight Galois group identification (cf. Q1 enum script)."""
    a4, a3, a2, a1, a0 = coeffs
    A4 = sp.Rational(a4); A3 = sp.Rational(a3); A2 = sp.Rational(a2)
    A1 = sp.Rational(a1); A0 = sp.Rational(a0)
    c3 = A3 / A4; c2 = A2 / A4; c1 = A1 / A4; c0 = A0 / A4
    s = -c3 / 4
    m = sp.symbols("_m")
    expr = ((m + s) ** 4 + c3 * (m + s) ** 3 + c2 * (m + s) ** 2
            + c1 * (m + s) + c0)
    poly = sp.Poly(sp.expand(expr), m)
    coeffs5 = poly.all_coeffs()
    while len(coeffs5) < 5:
        coeffs5 = [sp.Integer(0)] + coeffs5
    p, q, r = coeffs5[2], coeffs5[3], coeffs5[4]
    y = sp.symbols("_y")
    R = y ** 3 - 2 * p * y ** 2 + (p ** 2 - 4 * r) * y + q ** 2
    rat = [s for s in sp.solve(R, y) if s.is_rational]
    delta4 = Delta_4(coeffs)
    delta_is_square = (delta4 >= 0
                       and math.isqrt(delta4) ** 2 == delta4)
    if len(rat) >= 3:
        return "V_4"
    if len(rat) == 1:
        alpha = rat[0]
        df = (16 * p ** 4 * r - 4 * p ** 3 * q ** 2 - 128 * p ** 2 * r ** 2
              + 144 * p * q ** 2 * r - 27 * q ** 4 + 256 * r ** 3)
        test_val = sp.simplify((alpha ** 2 - 4 * r) * (alpha - p) * df)
        if test_val == 0:
            return "C_4_or_D_4"
        try:
            tv = sp.nsimplify(test_val, rational=True)
            num, den = sp.fraction(tv)
            num = int(num); den = int(den)
            if num * den < 0:
                return "D_4"
            r1 = math.isqrt(abs(num)); r2 = math.isqrt(abs(den))
            if r1 * r1 == abs(num) and r2 * r2 == abs(den):
                return "C_4"
            return "D_4"
        except Exception:
            return "D_4"
    return "A_4" if delta_is_square else "S_4"


# --- Q1 catalogue lookup ------------------------------------------

def load_q1_set():
    cat = json.load(open(SES_Q1 / "quartic_family_catalogue.json"))["families"]
    return {(f["alpha_4"], f["alpha_3"], f["alpha_2"],
             f["alpha_1"], f["alpha_0"]): f["family_id"] for f in cat}


# --- enumeration --------------------------------------------------

def enumerate_window(W, q1_set, max_log_every=2000):
    """Yield candidate j=0 (I=0) primitive quartics in the window
    (a4 in {1,2,3,5,7}; a3..a0 in [-W..W]).  Skips Q1's 60.  Skips
    duplicates from a smaller window if known."""
    a4_vals = [1, 2, 3, 5, 7]
    rng = list(range(-W, W + 1))
    n_total = 0
    n_primitive = 0
    n_jzero = 0
    n_irred_jzero = 0
    out = []
    for a4 in a4_vals:
        for a3 in rng:
            for a2 in rng:
                for a1 in rng:
                    for a0 in rng:
                        n_total += 1
                        coeffs = (a4, a3, a2, a1, a0)
                        if not primitive(coeffs):
                            continue
                        n_primitive += 1
                        if igusa_I(*coeffs) != 0:
                            continue
                        n_jzero += 1
                        # in Q1 already?
                        if coeffs in q1_set:
                            continue
                        if not is_irreducible_quartic(coeffs):
                            continue
                        if not is_squarefree_poly(coeffs):
                            continue
                        J = igusa_J(*coeffs)
                        if J == 0:
                            # Delta = -J^2 = 0 if J=0, so degenerate
                            continue
                        n_irred_jzero += 1
                        # full record
                        out.append({
                            "coeffs": list(coeffs),
                            "I_quartic": 0,
                            "J_quartic": int(J),
                            "Delta_4_exact": Delta_4(coeffs),
                            "in_Q1": False,
                        })
    return out, {"window": W, "n_candidates_total": n_total,
                 "n_primitive": n_primitive,
                 "n_I_zero": n_jzero,
                 "n_I_zero_new_irred_nondegenerate": n_irred_jzero}


# --- WKB tail fit (mirror of Phase A/B) ---------------------------

def cf_value(coeffs_leading_first, N, dps):
    with mp.workdps(dps):
        ms = [mp.mpf(c) for c in coeffs_leading_first]

        def b(k):
            v = ms[0]
            kk = mp.mpf(k)
            for c in ms[1:]:
                v = v * kk + c
            return v
        x = b(N)
        for k in range(N - 1, -1, -1):
            x = b(k) + mp.mpf(1) / x
        return +x


def compute_y(coeffs, n_grid, n_ref, dps):
    with mp.workdps(dps):
        L_ref = cf_value(coeffs, n_ref, dps)
        ns, ys = [], []
        for N in n_grid:
            d = abs(cf_value(coeffs, N, dps) - L_ref)
            if d == 0:
                continue
            ns.append(N)
            ys.append(float(mp.log(d)))
    return np.array(ns, dtype=float), np.array(ys, dtype=float)


def fit_pair(ns, ys, A_fix):
    z = ys + A_fix * ns * np.log(ns)
    X = np.column_stack([ns, -np.log(ns), np.ones_like(ns)])
    sol, *_ = np.linalg.lstsq(X, z, rcond=None)
    alpha, beta, gamma = (float(c) for c in sol)
    resid = z - X @ sol
    fixed = {"alpha": alpha, "beta": beta, "gamma": gamma,
             "n_points": int(len(ys)),
             "residual_mean": float(resid.mean()),
             "residual_std": float(resid.std(ddof=1)) if len(resid) > 1 else 0.0,
             "residual_at_max_n": float(resid[-1])}
    Xf = np.column_stack([-ns * np.log(ns), ns, -np.log(ns), np.ones_like(ns)])
    sol2, *_ = np.linalg.lstsq(Xf, ys, rcond=None)
    A, alphaf, betaf, gammaf = (float(c) for c in sol2)
    residf = ys - Xf @ sol2
    s2 = float(np.sum(residf ** 2) / max(len(ys) - 4, 1))
    cov = s2 * np.linalg.inv(Xf.T @ Xf)
    se = np.sqrt(np.diag(cov))
    free = {"A": A, "alpha": alphaf, "beta": betaf, "gamma": gammaf,
            "A_stderr": float(se[0]),
            "n_points": int(len(ys))}
    return fixed, free


# --- main ---------------------------------------------------------

def main():
    if LOG.exists():
        LOG.unlink()
    L(f"R13-C extended quartic enumeration: dps={DPS}, fit window N=[50,250], N_ref=400")
    q1_set = load_q1_set()
    L(f"Q1 catalogue: {len(q1_set)} known families")

    summary = {}
    final_window = 5
    j0_records = []
    for W in (5, 7, 10):
        L(f"--- Enumerating window W={W} (a3..a0 in [-{W}..{W}]) ---")
        recs, stats = enumerate_window(W, q1_set)
        L(f"  totals: candidates={stats['n_candidates_total']}, "
          f"primitive={stats['n_primitive']}, I=0={stats['n_I_zero']}, "
          f"I=0 new+irred+nondeg={stats['n_I_zero_new_irred_nondegenerate']}")
        summary[f"W={W}"] = stats
        if len(recs) >= 4:
            j0_records = recs
            final_window = W
            L(f"  cell size {len(recs)} >= 4 at W={W}, stopping escalation.")
            break
        else:
            j0_records = recs
            final_window = W
            if W < 10:
                L(f"  cell size {len(recs)} < 4, escalating window.")

    L(f"Final j=0 cell (extended, NEW): {len(j0_records)} families at W={final_window}")
    for r in j0_records:
        c = r["coeffs"]
        try:
            r["Galois_group"] = galois_quartic_quick(c)
        except Exception as e:
            r["Galois_group"] = f"err:{e}"

    # write extended catalogue (j=0 cell + their stats)
    cat_out = {
        "task_id": "PCF2-SESSION-R1.3",
        "phase": "C",
        "windows_tried": list(summary.keys()),
        "window_stats": summary,
        "final_window": final_window,
        "j_zero_new_cell": j0_records,
    }
    with open(HERE / "extended_quartic_catalogue.json", "w") as fh:
        json.dump(cat_out, fh, indent=2)
    L(f"  wrote extended_quartic_catalogue.json")

    # WKB on each new j=0 family
    L(f"--- WKB tail fit on {len(j0_records)} new j=0 quartics ---")
    A_FIX = 8.0
    j0_results = []
    for i, r in enumerate(j0_records):
        coeffs = r["coeffs"]
        t0 = time.time()
        ns, ys = compute_y(coeffs, N_GRID, N_REF, DPS)
        fixed, free = fit_pair(ns, ys, A_FIX)
        dt = time.time() - t0
        delta = free["A"] - 8.0
        r2 = dict(r)
        r2.update({
            "fixed_A8_residual_mean": fixed["residual_mean"],
            "fixed_A8_residual_std": fixed["residual_std"],
            "fixed_A8_residual_at_max_n": fixed["residual_at_max_n"],
            "fixed_A8_n_points": fixed["n_points"],
            "free_A": free["A"],
            "free_A_stderr": free["A_stderr"],
            "delta_R13_free": delta,
            "dt_seconds": dt,
        })
        j0_results.append(r2)
        L(f"  cand {i+1:>2d}/{len(j0_records)} coeffs={coeffs} Gal={r2['Galois_group']:<12s} "
          f"A_free={free['A']:.5f} (delta={delta:+.4e}) "
          f"resid@maxN={fixed['residual_at_max_n']:+.4e} "
          f"std={fixed['residual_std']:.4e} dt={dt:.1f}s")

    # save
    with open(HERE / "extended_j_zero_cell.json", "w") as fh:
        json.dump({
            "task_id": "PCF2-SESSION-R1.3",
            "phase": "C",
            "config": {"dps": DPS, "N_grid": [50, 250, 2], "N_ref": N_REF,
                       "A_fixed": A_FIX, "final_window": final_window},
            "results": j0_results,
        }, fh, indent=2)
    L("wrote extended_j_zero_cell.json")

    # CSV append to residualization_d4.csv: emit a separate file
    # `residualization_d4_extended.csv` (Q1 60 + new j=0).  We keep
    # the original residualization_d4.csv intact.
    if j0_results:
        with open(HERE / "residualization_d4_extended_jzero.csv", "w",
                  newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["coeffs", "Galois_group", "A_R13_free",
                        "A_stderr_R13_free", "delta_R13_free",
                        "residual_mean", "residual_std",
                        "residual_at_max_n", "n_points"])
            for r in j0_results:
                w.writerow([r["coeffs"], r["Galois_group"],
                            r["free_A"], r["free_A_stderr"],
                            r["delta_R13_free"],
                            r["fixed_A8_residual_mean"],
                            r["fixed_A8_residual_std"],
                            r["fixed_A8_residual_at_max_n"],
                            r["fixed_A8_n_points"]])

    # claims
    h = hashlib.sha256(
        (HERE / "extended_j_zero_cell.json").read_bytes()).hexdigest()
    claims = []
    def cl(text):
        claims.append({"claim": text, "evidence_type": "computation",
                       "dps": DPS, "reproducible": True,
                       "script": "r1_3_extended_enumeration.py",
                       "output_hash": h})
    cl(f"R13-C: extended-window quartic enumeration final_window=W={final_window}; "
       f"window stats: {json.dumps(summary)}")
    cl(f"R13-C: NEW j=0 cell (I=12*a4*a0-3*a3*a1+a2^2=0, primitive, "
       f"irreducible, non-singular, not in Q1's 60): "
       f"{len(j0_records)} family/families")
    if j0_results:
        deltas = [r["delta_R13_free"] for r in j0_results]
        resids = [r["fixed_A8_residual_at_max_n"] for r in j0_results]
        cl(f"R13-C: NEW j=0 quartic delta_R13_free range: "
           f"[{min(deltas):.4e}, {max(deltas):.4e}]; "
           f"residual_at_max_n range: [{min(resids):.4e}, {max(resids):.4e}]")
        # comparison with Q1 non-j=0 cluster
        with open(HERE / "residualization_d4.csv") as fh:
            q1rows = list(csv.DictReader(fh))
        cluster_deltas = [float(r["delta_R13_free"]) for r in q1rows
                          if float(r["log_abs_j"]) != 0.0]
        cl(f"R13-C: Q1 non-j=0 quartic cluster delta_R13_free range: "
           f"[{min(cluster_deltas):.4e}, {max(cluster_deltas):.4e}], "
           f"median={np.median(cluster_deltas):.4e}; new j=0 cell median "
           f"delta = {np.median(deltas):.4e}")

    with open(HERE / "claims_phase_C.jsonl", "w") as fh:
        for c in claims:
            fh.write(json.dumps(c) + "\n")
    L(f"wrote {len(claims)} claims to claims_phase_C.jsonl")


if __name__ == "__main__":
    main()
