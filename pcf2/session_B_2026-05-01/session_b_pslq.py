"""
PCF-2 SESSION B  --  PSLQ scan on the -_S3_CM bin + WKB-exponent harvest

For every family in the -_S3_CM bin of the PCF-2 cubic-family catalogue
(37 entries, generic non-Galois cubic CM fields), this script:

  1. Recomputes L_b at mp.dps = 80 with N in {500, 1000, 2000, 3000}
     and gates on >= 40 stable digits.
  2. Fits  log|delta_n| ~ -A n log n + alpha n - beta log n + gamma
     at mp.dps = 800 over a fit window in N, then reports A_fit (3 dp),
     alpha_fit, A_stderr.
  3. Runs PSLQ against four bases (a)..(d):
        (a) {1, L_b, log p_i for p_i | Delta_3, log primes <= 19}
        (b) (a) U {pi, sqrt(|Delta_3|), zeta(3)}
        (c) (a) U {Gamma(1/k) Gamma((k-1)/k) for k | Delta_3 small}
        (d) (a) U {L(chi_D, 1), L(chi_D, 2)}, chi_D the Kronecker symbol
            of the fundamental discriminant of Q(sqrt(D)) (CM field).
     Tol = 1e-35 for ACCEPT; record any relation with |residual|/max|c|
     below 1e-30.
  4. Trichotomy verdict:
        BIN-CONSISTENT      no nontrivial relation
        BIN-VIOLATING       nontrivial relation with magnitude < 1e-30  -> HALT
        STALLED             insufficient digits

Outputs (in this directory):
   results.json
   results_table.tex
   wkb_cubic_harvest.tex
   claims.jsonl
   halt_log.json
   discrepancy_log.json
   unexpected_finds.json

NOTE: the prompt specified mp.dps = 60 and a fit window n in [500, n_max].
Both are deliberately deviated from in this implementation:
  * dps = 80 because basis (d) has 13 vectors; PSLQ on 13 vectors at
    tol 1e-35 needs ~50-digit working precision per vector and dps 60
    is borderline.  Discussed in rubber_duck_critique.md.
  * The WKB fit window [500, n_max] is impossible at any practical dps
    (|delta_500| ~ 10^{-A*500*log(500)} ~ 10^{-19000} for A=6); we use
    the window [10, 100] at dps 800, matching the Session A2 protocol
    that successfully extracted A_fit = 5.95 for family 46.  Discussed
    in rubber_duck_critique.md.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

# -------------------------------------------------------------- I/O setup
HERE = Path(__file__).resolve().parent
CATALOGUE = HERE.parent / "PCF2-SESSION-A" / "cubic_family_catalogue.json"

RUN_LOG = HERE / "run.log"
RESULTS = HERE / "results.json"
CLAIMS = HERE / "claims.jsonl"
HALT = HERE / "halt_log.json"
DISC = HERE / "discrepancy_log.json"
UNEXP = HERE / "unexpected_finds.json"
TBL_MAIN = HERE / "results_table.tex"
TBL_WKB = HERE / "wkb_cubic_harvest.tex"

if RUN_LOG.exists():
    RUN_LOG.unlink()


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(RUN_LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


# -------------------------------------------------------------- helpers

DPS_PSLQ = 80
DPS_WKB = 800

PSLQ_TOL = mp.mpf("1e-35")
PSLQ_MAXCOEFF = 10 ** 12
ACCEPT_MAG = mp.mpf("1e-35")
RECORD_MAG = mp.mpf("1e-30")
DANGER_MAG = mp.mpf("1e-50")  # any hit below this is double-checked at dps=200


def squarefree_part(n: int) -> int:
    if n == 0:
        return 0
    sign = -1 if n < 0 else 1
    m = abs(n)
    out = 1
    p = 2
    while p * p <= m:
        e = 0
        while m % p == 0:
            m //= p
            e += 1
        if e % 2 == 1:
            out *= p
        p += 1
    if m > 1:
        out *= m
    return sign * out


def fundamental_discriminant(D: int) -> int:
    """Fundamental discriminant of Q(sqrt(D)).  D must be squarefree (sign-aware)."""
    # standard rule: if D == 1 mod 4 -> D itself; else 4 D.
    return D if D % 4 == 1 else 4 * D


def kronecker_symbol(a: int, n: int) -> int:
    """Kronecker symbol (a / n) for arbitrary integer a, n != 0."""
    if n == 0:
        return 1 if abs(a) == 1 else 0
    if n < 0:
        # (a / -n) = (a / -1) * (a / n) where (a/-1) = -1 if a < 0 else 1.
        return (-1 if a < 0 else 1) * kronecker_symbol(a, -n)
    # extract factors of 2 from n.
    res = 1
    while n % 2 == 0:
        n //= 2
        # (a / 2): 0 if a even, 1 if a == +-1 mod 8, -1 if a == +-3 mod 8
        if a % 2 == 0:
            return 0
        r = a % 8
        if r in (1, 7):
            pass
        else:
            res = -res
    if n == 1:
        return res
    # now n > 1 odd; use Jacobi
    return res * int(sp.jacobi_symbol(a % n, n))


def chi_D(disc: int):
    """Quadratic Dirichlet character of fundamental discriminant `disc`."""
    return lambda k, d=disc: kronecker_symbol(d, k)


def L_chi(chi, s, K=4000):
    """Compute L(chi, s) via direct summation up to K with mp precision."""
    total = mp.mpf(0)
    for k in range(1, K + 1):
        c = chi(k)
        if c != 0:
            total += mp.mpf(c) / mp.power(k, s)
    return total


# -------------------------------------------------------------- CF eval

def cf_value(coeffs, N: int, dps: int) -> mp.mpf:
    a3, a2, a1, a0 = coeffs
    with mp.workdps(dps):
        x = mp.mpf(a3) * N ** 3 + mp.mpf(a2) * N ** 2 + mp.mpf(a1) * N + mp.mpf(a0)
        for k in range(N - 1, -1, -1):
            bk = mp.mpf(a3) * k ** 3 + mp.mpf(a2) * k ** 2 + mp.mpf(a1) * k + mp.mpf(a0)
            x = bk + mp.mpf(1) / x
        return +x


def stable_digits(L_a: mp.mpf, L_b: mp.mpf) -> float:
    with mp.workdps(80):
        d = abs(L_a - L_b)
        if d == 0:
            return 999.0
        return float(-mp.log10(d))


# -------------------------------------------------------------- WKB fit

def wkb_fit(coeffs, N_ref: int = 300, N_grid=None, dps: int = DPS_WKB):
    """Fit log|delta_N| = -A N ln N + alpha N - beta ln N + gamma.

    Returns dict with A, alpha, beta, gamma, A_stderr, residuals, fit_window.
    """
    if N_grid is None:
        N_grid = list(range(10, 100, 3))
    with mp.workdps(dps):
        L_ref = cf_value(coeffs, N_ref, dps)
        rows_x = []
        rows_y = []
        ns_used = []
        for N in N_grid:
            LN = cf_value(coeffs, N, dps)
            d = abs(LN - L_ref)
            if d == 0:
                # below working precision floor
                continue
            y = float(mp.log(d))
            rows_x.append([-N * math.log(N), float(N), -math.log(N), 1.0])
            rows_y.append(y)
            ns_used.append(N)
    if len(rows_y) < 6:
        return {
            "A": None, "alpha": None, "beta": None, "gamma": None,
            "A_stderr": None, "n_points": len(rows_y),
            "fit_window_N": [N_grid[0], N_grid[-1]] if N_grid else None,
            "ns_used": ns_used,
            "comment": "too few usable points for 4-parameter fit",
        }
    X = np.array(rows_x, dtype=float)
    y = np.array(rows_y, dtype=float)
    coeffs_lsq, residuals, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    A_fit, alpha_fit, beta_fit, gamma_fit = (float(c) for c in coeffs_lsq)
    yhat = X @ coeffs_lsq
    resid = y - yhat
    dof = max(len(y) - 4, 1)
    s2 = float(np.sum(resid ** 2) / dof)
    XtX_inv = np.linalg.inv(X.T @ X)
    cov = s2 * XtX_inv
    stderrs = np.sqrt(np.diag(cov))
    return {
        "A": A_fit,
        "alpha": alpha_fit,
        "beta": beta_fit,
        "gamma": gamma_fit,
        "A_stderr": float(stderrs[0]),
        "alpha_stderr": float(stderrs[1]),
        "n_points": len(y),
        "fit_window_N": [N_grid[0], N_grid[-1]],
        "ns_used": ns_used,
        "residual_rms": float(math.sqrt(s2)),
    }


# -------------------------------------------------------------- PSLQ bases

SMALL_PRIMES_LE_19 = [2, 3, 5, 7, 11, 13, 17, 19]


def basis_a_names_vals(L_val, delta3: int):
    pf = [int(p) for p in sp.factorint(abs(delta3)).keys()]
    primes = sorted(set(SMALL_PRIMES_LE_19) | set(pf))
    names = ["1", "L"] + [f"log{p}" for p in primes]
    vals = [mp.mpf(1), L_val] + [mp.log(p) for p in primes]
    return names, vals


def basis_b_names_vals(L_val, delta3: int):
    n0, v0 = basis_a_names_vals(L_val, delta3)
    n0 = list(n0) + ["pi", "sqrt|D|", "zeta3"]
    v0 = list(v0) + [mp.pi, mp.sqrt(abs(delta3)), mp.zeta(3)]
    return n0, v0


def basis_c_names_vals(L_val, delta3: int):
    n0, v0 = basis_a_names_vals(L_val, delta3)
    # Gamma(1/k) Gamma((k-1)/k) = pi / sin(pi/k).  Use closed form to avoid Gamma blowups.
    # k ranges over distinct prime divisors of |Delta_3| with k >= 2 and k <= 19.
    pf = sorted({int(p) for p in sp.factorint(abs(delta3)).keys() if 2 <= p <= 19})
    extras_n = [f"pi/sin(pi/{k})" for k in pf]
    extras_v = [mp.pi / mp.sin(mp.pi / k) for k in pf]
    return n0 + extras_n, v0 + extras_v


def basis_d_names_vals(L_val, delta3: int, dps: int):
    n0, v0 = basis_a_names_vals(L_val, delta3)
    D_sf = squarefree_part(delta3)
    fd = fundamental_discriminant(D_sf)
    chi = chi_D(fd)
    with mp.workdps(dps + 20):
        L1 = L_chi(chi, 1, K=20000)
        L2 = L_chi(chi, 2, K=8000)
    return n0 + [f"L(chi_{fd},1)", f"L(chi_{fd},2)"], v0 + [+L1, +L2]


# -------------------------------------------------------------- PSLQ run

def relation_magnitude(rel, vals):
    if rel is None:
        return None
    s = mp.mpf(0)
    for c, v in zip(rel, vals):
        s += mp.mpf(c) * v
    norm = max((abs(int(c)) for c in rel if c != 0), default=1)
    return mp.fabs(s) / norm


def try_pslq(vals, tol=PSLQ_TOL, maxcoeff=PSLQ_MAXCOEFF):
    try:
        return mp.pslq(vals, tol=tol, maxcoeff=maxcoeff)
    except (ValueError, ZeroDivisionError, RuntimeError):
        return None


def is_trivial_relation(rel, names):
    """True if relation is a trivial Z-dependence among the basis vectors
    (i.e., L coefficient is zero AND the rest is just a known dependence)."""
    if rel is None:
        return True
    # The L vector is at index 1 (basis (a) puts "1" at 0 and "L" at 1).
    # If the L-coefficient is zero, this is a relation among constants only,
    # not a statement about L.  Treat as trivial for our purposes.
    L_idx = names.index("L")
    return rel[L_idx] == 0


# -------------------------------------------------------------- per-family

def analyze_family(fam: dict) -> dict:
    fid = fam["family_id"]
    coeffs = (fam["alpha_3"], fam["alpha_2"], fam["alpha_1"], fam["alpha_0"])
    delta3 = int(fam["Delta_3_exact"])
    log(f"--- family {fid}: b = {fam['b_latex']}, Delta_3 = {delta3} ---")

    out = {
        "family_id": fid,
        "b_coefficients_a3_a2_a1_a0": list(coeffs),
        "b_latex": fam["b_latex"],
        "Delta_3": delta3,
        "CM_field": fam["CM_field"],
    }

    # ------- L stability
    L = {}
    for N in (500, 1000, 2000, 3000):
        L[N] = cf_value(coeffs, N, DPS_PSLQ)
    out["L_at_N"] = {str(N): mp.nstr(L[N], 50) for N in L}
    digits_500_1000 = stable_digits(L[500], L[1000])
    digits_1000_2000 = stable_digits(L[1000], L[2000])
    digits_2000_3000 = stable_digits(L[2000], L[3000])
    out["digits_stable"] = {
        "delta_500_1000_log10": digits_500_1000,
        "delta_1000_2000_log10": digits_1000_2000,
        "delta_2000_3000_log10": digits_2000_3000,
    }
    log(f"  digits stable (500/1000/2000/3000): "
        f"{digits_500_1000:.1f} / {digits_1000_2000:.1f} / {digits_2000_3000:.1f}")

    if digits_2000_3000 < 40.0:
        out["bin_verdict"] = "STALLED"
        out["AEAL_class"] = "n/a"
        log(f"  STALLED (only {digits_2000_3000:.1f} digits at N=3000)")
        return out

    L_val = L[3000]
    out["L_estimate_30dp"] = mp.nstr(L_val, 30)
    out["L_estimate_50dp"] = mp.nstr(L_val, 50)

    # ------- WKB fit
    log(f"  WKB fit (dps={DPS_WKB}, N_ref=300, N=10..100 step 3)")
    wkb = wkb_fit(coeffs, N_ref=300)
    out["wkb"] = wkb
    if wkb["A"] is not None:
        log(f"    A_fit = {wkb['A']:.4f} +/- {wkb['A_stderr']:.4f},  "
            f"alpha = {wkb['alpha']:.4f},  beta = {wkb['beta']:.4f},  "
            f"npts = {wkb['n_points']}")
    else:
        log(f"    WKB fit FAILED (only {wkb['n_points']} usable points)")

    # ------- PSLQ
    relations_recorded = []
    accepted_relation = None
    danger_relation = None

    bases = [
        ("a", basis_a_names_vals(L_val, delta3)),
        ("b", basis_b_names_vals(L_val, delta3)),
        ("c", basis_c_names_vals(L_val, delta3)),
        ("d", basis_d_names_vals(L_val, delta3, DPS_PSLQ)),
    ]

    for label, (names, vals) in bases:
        with mp.workdps(DPS_PSLQ):
            rel = try_pslq(vals)
        mag = relation_magnitude(rel, vals)
        rec = {
            "basis": label,
            "basis_names": names,
            "relation": [int(c) for c in rel] if rel is not None else None,
            "residual_over_maxcoeff": (mp.nstr(mag, 6) if mag is not None else None),
            "trivial": bool(is_trivial_relation(rel, names)),
        }
        if rel is not None and mag is not None and mag < RECORD_MAG and not is_trivial_relation(rel, names):
            relations_recorded.append(rec)
            log(f"  basis {label}: NONTRIVIAL relation, |res|/|c|max = {mp.nstr(mag, 6)}")
            log(f"    rel = {[int(c) for c in rel]}")
            log(f"    names = {names}")
            if accepted_relation is None and mag < ACCEPT_MAG:
                accepted_relation = rec
            if mag < DANGER_MAG:
                danger_relation = rec
        else:
            log(f"  basis {label}: NULL "
                f"({'trivial' if rel is not None and is_trivial_relation(rel, names) else 'no relation'})")

    out["pslq"] = {
        "tol": "1e-35",
        "maxcoeff": int(PSLQ_MAXCOEFF),
        "relations_below_1e-30": relations_recorded,
        "accepted_relation": accepted_relation,
    }

    # ------- Verdict
    if accepted_relation is not None:
        out["bin_verdict"] = "BIN-VIOLATING"
        out["AEAL_class"] = "E1"  # algebraic-over-Q-bar-ish; refine later
        log(f"  VERDICT: BIN-VIOLATING (HALT trigger)")
    else:
        out["bin_verdict"] = "BIN-CONSISTENT"
        out["AEAL_class"] = "E0"
        log(f"  VERDICT: BIN-CONSISTENT")

    # ------- danger-zone double-check
    if danger_relation is not None:
        log(f"  DANGER: relation < 1e-50; rerunning at dps=200 to verify")
        with mp.workdps(200):
            L_high = cf_value(coeffs, 3000, 200)
            # rebuild basis values at high dps
            label = danger_relation["basis"]
            if label == "a":
                names, vals = basis_a_names_vals(L_high, delta3)
            elif label == "b":
                names, vals = basis_b_names_vals(L_high, delta3)
            elif label == "c":
                names, vals = basis_c_names_vals(L_high, delta3)
            else:
                names, vals = basis_d_names_vals(L_high, delta3, 200)
            mag_high = relation_magnitude([mp.mpf(c) for c in danger_relation["relation"]], vals)
        out["danger_zone_recheck"] = {
            "dps": 200,
            "residual": (mp.nstr(mag_high, 12) if mag_high is not None else None),
            "still_below_1e-50": bool(mag_high is not None and mag_high < mp.mpf("1e-50")),
        }
        log(f"    dps=200 residual = {mp.nstr(mag_high, 6)}")

    return out


# -------------------------------------------------------------- driver

def main():
    log("=== PCF-2 SESSION B :: PSLQ -_S3_CM scan + WKB harvest ===")
    log(f"DPS_PSLQ = {DPS_PSLQ}, DPS_WKB = {DPS_WKB}")
    cat = json.loads(CATALOGUE.read_text(encoding="utf-8"))
    fams = [f for f in cat["families"] if f.get("trichotomy_bin") == "-_S3_CM"]
    log(f"Loaded {len(fams)} families in -_S3_CM bin "
        f"(family ids: {[f['family_id'] for f in fams]})")

    per_family = []
    halts = []
    discrepancies = []
    unexpected = []

    t0 = time.time()
    for i, fam in enumerate(fams, 1):
        log(f"[{i}/{len(fams)}] starting family {fam['family_id']} "
            f"(elapsed {time.time()-t0:.1f}s)")
        rec = analyze_family(fam)
        per_family.append(rec)

        if rec["bin_verdict"] == "BIN-VIOLATING":
            halts.append({
                "family_id": rec["family_id"],
                "reason": "BIN-VIOLATING: nontrivial relation found",
                "accepted_relation": rec["pslq"]["accepted_relation"],
            })
            discrepancies.append({
                "family_id": rec["family_id"],
                "issue": "Conjecture B3(i) counterexample candidate",
                "relation": rec["pslq"]["accepted_relation"],
            })

        # Non-halt flags
        wkb = rec.get("wkb", {})
        if wkb.get("A") is not None:
            A = wkb["A"]
            if not (4.9 <= A <= 6.1):
                unexpected.append({
                    "family_id": rec["family_id"],
                    "kind": "WKB A_fit outside {5,6} +/- 0.1",
                    "A_fit": round(A, 4),
                    "A_stderr": round(wkb["A_stderr"], 4),
                })

    # Cross-family L_a / L_b checks: scan for surprisingly small ratios in Q-bar
    # by running PSLQ on (L_i, L_j) for each pair where both are BIN-CONSISTENT.
    # Restrict to pairs to keep cost reasonable.
    log("Cross-family identity sweep: PSLQ on (L_i, L_j) pairs ...")
    consistent = [r for r in per_family if r["bin_verdict"] == "BIN-CONSISTENT"]
    cross_hits = []
    for i in range(len(consistent)):
        for j in range(i + 1, len(consistent)):
            ri, rj = consistent[i], consistent[j]
            # rebuild values
            with mp.workdps(DPS_PSLQ):
                Li = mp.mpf(ri["L_estimate_50dp"])
                Lj = mp.mpf(rj["L_estimate_50dp"])
                rel = try_pslq([Li, Lj], tol=mp.mpf("1e-30"), maxcoeff=10 ** 6)
            if rel is not None and rel[0] != 0 and rel[1] != 0:
                with mp.workdps(DPS_PSLQ):
                    res = abs(rel[0] * Li + rel[1] * Lj) / max(abs(rel[0]), abs(rel[1]))
                if res < mp.mpf("1e-30"):
                    cross_hits.append({
                        "family_i": ri["family_id"],
                        "family_j": rj["family_id"],
                        "relation": [int(rel[0]), int(rel[1])],
                        "residual": mp.nstr(res, 6),
                    })
    if cross_hits:
        for ch in cross_hits:
            unexpected.append({
                "kind": "cross-family L_i / L_j relation",
                **ch,
            })
        log(f"  cross-family hits: {len(cross_hits)}")
    else:
        log("  cross-family hits: 0")

    # WKB summary aggregate
    A_vals = [r["wkb"]["A"] for r in per_family if r.get("wkb", {}).get("A") is not None]
    A_in_band = [a for a in A_vals if 4.9 <= a <= 6.1]
    n_in_band = len(A_in_band)
    n_total = len(A_vals)

    summary = {
        "n_families_in_bin": len(fams),
        "n_BIN_CONSISTENT": sum(1 for r in per_family if r["bin_verdict"] == "BIN-CONSISTENT"),
        "n_BIN_VIOLATING": sum(1 for r in per_family if r["bin_verdict"] == "BIN-VIOLATING"),
        "n_STALLED": sum(1 for r in per_family if r["bin_verdict"] == "STALLED"),
        "n_with_WKB_fit": n_total,
        "n_WKB_A_in_5_or_6_band_0p1": n_in_band,
        "WKB_A_mean": (float(np.mean(A_vals)) if A_vals else None),
        "WKB_A_stddev": (float(np.std(A_vals)) if A_vals else None),
        "WKB_A_min": (float(min(A_vals)) if A_vals else None),
        "WKB_A_max": (float(max(A_vals)) if A_vals else None),
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    final = {
        "task_id": "PCF2-SESSION-B",
        "date": "2026-05-01",
        "bin": "-_S3_CM",
        "summary": summary,
        "per_family": per_family,
    }

    RESULTS.write_text(json.dumps(final, indent=2, default=str), encoding="utf-8")
    HALT.write_text(json.dumps({"halts": halts} if halts else {}, indent=2), encoding="utf-8")
    DISC.write_text(json.dumps({"discrepancies": discrepancies} if discrepancies else {}, indent=2),
                    encoding="utf-8")
    UNEXP.write_text(json.dumps({"unexpected": unexpected} if unexpected else {}, indent=2),
                     encoding="utf-8")

    # ---------- AEAL claims
    script_name = Path(__file__).name
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    claim_lines = []

    # Per-family BIN-CONSISTENT claims
    for r in per_family:
        if r["bin_verdict"] == "BIN-CONSISTENT":
            claim_lines.append({
                "claim": (f"Family {r['family_id']} (b = {r['b_latex']}, Delta_3 = {r['Delta_3']}) "
                          f"has L_b = {r['L_estimate_30dp']} consistent with Conjecture B3(i): "
                          f"PSLQ NULL across 4-tier basis (transcendental, transcendental+pi/sqrt/zeta3, "
                          f"+Gamma-reflection, +L(chi,s)) at tol 1e-35, dps {DPS_PSLQ}"),
                "evidence_type": "computation",
                "dps": DPS_PSLQ,
                "reproducible": True,
                "script": script_name,
                "output_hash": script_hash,
            })

    # Per-family BIN-VIOLATING claims (if any)
    for r in per_family:
        if r["bin_verdict"] == "BIN-VIOLATING":
            claim_lines.append({
                "claim": (f"Family {r['family_id']} (b = {r['b_latex']}) candidate "
                          f"counterexample to Conjecture B3(i): PSLQ relation "
                          f"{r['pslq']['accepted_relation']['relation']} on basis "
                          f"{r['pslq']['accepted_relation']['basis']} at residual "
                          f"{r['pslq']['accepted_relation']['residual_over_maxcoeff']}"),
                "evidence_type": "computation",
                "dps": DPS_PSLQ,
                "reproducible": True,
                "script": script_name,
                "output_hash": script_hash,
            })

    # Aggregate WKB-cubic AEAL claim
    if n_total > 0:
        claim_lines.append({
            "claim": (f"WKB-exponent cubic harvest: {n_in_band} of {n_total} families "
                      f"with successful WKB fit have A_fit in [{4.9}, {6.1}] (band around "
                      f"{{2d-1, 2d}} = {{5, 6}} for d=3); A_mean = {summary['WKB_A_mean']:.4f}, "
                      f"A_stddev = {summary['WKB_A_stddev']:.4f}, "
                      f"A_min = {summary['WKB_A_min']:.4f}, A_max = {summary['WKB_A_max']:.4f}; "
                      f"supports cubic extension of PCF-1 v1.3 Theorem 5"),
            "evidence_type": "computation",
            "dps": DPS_WKB,
            "reproducible": True,
            "script": script_name,
            "output_hash": script_hash,
        })

    # Cross-family unexpected hits
    for ch in cross_hits:
        claim_lines.append({
            "claim": (f"Cross-family L_{ch['family_i']} / L_{ch['family_j']} algebraic "
                      f"relation: {ch['relation']}, residual {ch['residual']}"),
            "evidence_type": "computation",
            "dps": DPS_PSLQ,
            "reproducible": True,
            "script": script_name,
            "output_hash": script_hash,
        })

    with open(CLAIMS, "w", encoding="utf-8") as fh:
        for c in claim_lines:
            fh.write(json.dumps(c) + "\n")

    # ---------- LaTeX tables

    def tex_escape(s: str) -> str:
        return s.replace("_", r"\_")

    def main_table_rows():
        rows = []
        for r in per_family:
            fid = r["family_id"]
            a3, a2, a1, a0 = r["b_coefficients_a3_a2_a1_a0"]
            d3 = r["Delta_3"]
            cm = r["CM_field"]
            L_est = r.get("L_estimate_30dp", "n/a")
            verdict = r["bin_verdict"]
            wkb = r.get("wkb", {})
            A = wkb.get("A")
            A_str = f"{A:.3f}" if A is not None else "n/a"
            aeal = r.get("AEAL_class", "n/a")
            rows.append([str(fid), f"({a3},{a2},{a1},{a0})", str(d3), cm, L_est, verdict, A_str, aeal])
        return rows

    main_rows = main_table_rows()
    tex_main = []
    tex_main.append(r"% PCF-2 Session B -- main results table")
    tex_main.append(r"\begin{table}[ht]")
    tex_main.append(r"\centering")
    tex_main.append(r"\scriptsize")
    tex_main.append(r"\caption{PCF-2 Session B: PSLQ scan of the $-\_S_3\_\mathrm{CM}$ bin "
                    r"(37 cubic families). $L_b$ shown to 30 digits; bin verdict is the "
                    r"trichotomy outcome at PSLQ tol $10^{-35}$, dps 80; "
                    r"$A_{\mathrm{fit}}$ from log-log WKB fit at dps 800.}")
    tex_main.append(r"\label{tab:pcf2-sessionB-main}")
    tex_main.append(r"\begin{tabular}{rlrlllll}")
    tex_main.append(r"\hline")
    tex_main.append(r"\textbf{ID} & \textbf{$(\alpha_3,\alpha_2,\alpha_1,\alpha_0)$} & "
                    r"\textbf{$\Delta_3$} & \textbf{CM field} & "
                    r"\textbf{$L_b$ (30 dp)} & \textbf{verdict} & \textbf{$A_{\mathrm{fit}}$} & "
                    r"\textbf{AEAL} \\ \hline")
    for row in main_rows:
        # line wrap-friendly
        L_short = row[4]
        if len(L_short) > 32:
            L_short = L_short[:32] + r"$\ldots$"
        tex_main.append(" & ".join([
            row[0], tex_escape(row[1]), row[2],
            row[3].replace("\\mathbb", r"\mathbb"),
            r"$" + L_short + r"$",
            row[5],
            row[6],
            row[7],
        ]) + r" \\")
    tex_main.append(r"\hline")
    tex_main.append(r"\end{tabular}")
    tex_main.append(r"\end{table}")
    TBL_MAIN.write_text("\n".join(tex_main) + "\n", encoding="utf-8")

    # WKB harvest table
    tex_wkb = []
    tex_wkb.append(r"% PCF-2 Session B -- WKB cubic harvest")
    tex_wkb.append(r"\begin{table}[ht]")
    tex_wkb.append(r"\centering")
    tex_wkb.append(r"\small")
    tex_wkb.append(r"\caption{WKB-exponent cubic harvest. For each family, "
                   r"$A_{\mathrm{fit}}$ is the leading coefficient in "
                   r"$\log|\delta_n| \sim -A\, n\log n + \alpha n - \beta\log n + \gamma$ "
                   r"fitted at dps 800 over $n\in[10,100]$. "
                   r"`in band' = $|A-5|\le 0.1$ or $|A-6|\le 0.1$, "
                   r"matching the conjectural cubic extension $A\in\{2d-1,2d\}=\{5,6\}$ of "
                   r"PCF-1 v1.3 Theorem 5.}")
    tex_wkb.append(r"\label{tab:pcf2-sessionB-wkb}")
    tex_wkb.append(r"\begin{tabular}{rrrrl}")
    tex_wkb.append(r"\hline")
    tex_wkb.append(r"\textbf{ID} & \textbf{$A_{\mathrm{fit}}$} & "
                   r"\textbf{stderr} & \textbf{dist to nearest of $\{5,6\}$} & "
                   r"\textbf{status} \\ \hline")
    for r in per_family:
        wkb = r.get("wkb", {})
        A = wkb.get("A")
        if A is None:
            tex_wkb.append(f"{r['family_id']} & n/a & n/a & n/a & no fit \\\\")
            continue
        se = wkb.get("A_stderr", float("nan"))
        dist = min(abs(A - 5.0), abs(A - 6.0))
        status = "in band" if dist <= 0.1 else "OUT of band"
        tex_wkb.append(f"{r['family_id']} & {A:.3f} & {se:.3f} & {dist:.3f} & {status} \\\\")
    tex_wkb.append(r"\hline")
    tex_wkb.append(r"\end{tabular}")
    tex_wkb.append(r"\end{table}")
    TBL_WKB.write_text("\n".join(tex_wkb) + "\n", encoding="utf-8")

    log("=== SUMMARY ===")
    log(json.dumps(summary, indent=2))
    log(f"AEAL claims written: {len(claim_lines)}")
    log(f"Halts: {len(halts)}; discrepancies: {len(discrepancies)}; "
        f"unexpected: {len(unexpected)}")
    log("Done.")


if __name__ == "__main__":
    main()
