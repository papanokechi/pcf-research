"""
PCF-2 SESSION C1  --  WKB-exponent harvest, non-CM cubic bins (+_S3_real, +_C3_real)

Completes degree-3 coverage of the cubic-family catalogue (50/50 with Session B)
and tests refined Conjecture B4':

    For a degree-d PCF (1, b) of generic type in scope:
        A_fit = 2d - 1   if trichotomy bin is elementary-positive
                          (+_S_d real or +_C_d real),
        A_fit = 2d       otherwise.

For d = 3 the elementary-positive prediction is A = 5; Session B established
A_fit = 6 unimodally on the -_S3_CM bin.

Pipeline (per family):
  1. Compute L_b at mp.dps = 80, n in {50, 100, 150, 200}; gate on >= 30
     stable digits between n=150 and n=200.
  2. WKB fit (identical to Session B step 1b):
        log|delta_n| ~ -A n log n + alpha n - beta log n + gamma
     fitted at mp.dps = 800 over n in [10, 100] step 3 (ordinary LSQ).
  3. B4'-verdict per family:
        b4_prime_consistent  iff |A_fit - predicted_branch| <= 0.10
        b4_prime_violating   iff min(|A_fit - 5|, |A_fit - 6|) > 0.30
        b4_prime_ambiguous   else
     predicted_branch = 5 for both +_S3_real and +_C3_real (elem-positive).

NO PSLQ scan in C1.  Striking convergence (>= 50 stable digits at n=200)
flagged in unexpected_finds.json.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
CATALOGUE = HERE.parent / "PCF2-SESSION-A" / "cubic_family_catalogue.json"

RUN_LOG = HERE / "run.log"
RESULTS = HERE / "results.json"
CLAIMS = HERE / "claims.jsonl"
HALT = HERE / "halt_log.json"
DISC = HERE / "discrepancy_log.json"
UNEXP = HERE / "unexpected_finds.json"
TBL_NONCM = HERE / "wkb_noncm_table.tex"
TBL_FULL = HERE / "wkb_cubic_harvest_v2.tex"

if RUN_LOG.exists():
    RUN_LOG.unlink()


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(RUN_LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


# ------------------------------------------------------------ params

DPS_L = 80          # for L_b stability check at small N
DPS_WKB = 800       # for WKB fit (matches Session B)
N_GRID_L = (50, 100, 150, 200)
WKB_NS = list(range(10, 100, 3))
N_REF = 300

PREDICTED_BRANCH = {"+_S3_real": 5, "+_C3_real": 5}
NON_CM_BINS = ("+_S3_real", "+_C3_real")


# ------------------------------------------------------------ CF eval

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


# ------------------------------------------------------------ WKB fit

def wkb_fit(coeffs, N_ref: int = N_REF, N_grid=None, dps: int = DPS_WKB):
    if N_grid is None:
        N_grid = list(WKB_NS)
    with mp.workdps(dps):
        L_ref = cf_value(coeffs, N_ref, dps)
        rows_x, rows_y, ns_used = [], [], []
        for N in N_grid:
            LN = cf_value(coeffs, N, dps)
            d = abs(LN - L_ref)
            if d == 0:
                continue
            y = float(mp.log(d))
            rows_x.append([-N * math.log(N), float(N), -math.log(N), 1.0])
            rows_y.append(y)
            ns_used.append(N)
    if len(rows_y) < 6:
        return {
            "A": None, "alpha": None, "beta": None, "gamma": None,
            "A_stderr": None, "alpha_stderr": None,
            "n_points": len(rows_y),
            "fit_window_N": [N_grid[0], N_grid[-1]] if N_grid else None,
            "ns_used": ns_used,
            "comment": "too few usable points for 4-parameter fit",
        }
    X = np.array(rows_x, dtype=float)
    y = np.array(rows_y, dtype=float)
    coeffs_lsq, *_ = np.linalg.lstsq(X, y, rcond=None)
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


# ------------------------------------------------------------ verdict

def b4_prime_verdict(A_fit: float, predicted: int):
    d5 = abs(A_fit - 5.0)
    d6 = abs(A_fit - 6.0)
    d_pred = abs(A_fit - float(predicted))
    if d_pred <= 0.10:
        return "b4_prime_consistent"
    if min(d5, d6) > 0.30:
        return "b4_prime_violating"
    if 5.0 <= A_fit <= 6.0 and d5 > 0.10 and d6 > 0.10:
        return "b4_prime_ambiguous"
    # closer to the wrong branch (within 0.10) -- still a B4' falsifier of the
    # bin assignment but not absolute B4' falsifier (A in {5,6})
    if d5 <= 0.10 or d6 <= 0.10:
        return "b4_prime_consistent_wrong_branch"
    return "b4_prime_ambiguous"


# ------------------------------------------------------------ per-family

def analyze_family(fam: dict) -> dict:
    fid = fam["family_id"]
    bin_ = fam["trichotomy_bin"]
    coeffs = (fam["alpha_3"], fam["alpha_2"], fam["alpha_1"], fam["alpha_0"])
    delta3 = int(fam["Delta_3_exact"])
    log(f"--- family {fid} [{bin_}]: b = {fam['b_latex']}, Delta_3 = {delta3} ---")

    out = {
        "family_id": fid,
        "b_coefficients_a3_a2_a1_a0": list(coeffs),
        "b_latex": fam["b_latex"],
        "Delta_3": delta3,
        "trichotomy_bin": bin_,
    }

    # ---- L stability gate
    L = {N: cf_value(coeffs, N, DPS_L) for N in N_GRID_L}
    out["L_at_N"] = {str(N): mp.nstr(L[N], 50) for N in L}
    d_50_100 = stable_digits(L[50], L[100])
    d_100_150 = stable_digits(L[100], L[150])
    d_150_200 = stable_digits(L[150], L[200])
    out["digits_stable"] = {
        "delta_50_100_log10": d_50_100,
        "delta_100_150_log10": d_100_150,
        "delta_150_200_log10": d_150_200,
    }
    digits_stable_at_200 = d_150_200
    out["digits_stable_at_200"] = digits_stable_at_200
    log(f"  digits stable (50/100/150/200 pairs): "
        f"{d_50_100:.1f} / {d_100_150:.1f} / {d_150_200:.1f}")

    if digits_stable_at_200 < 30.0:
        out["wkb_fit_run"] = False
        out["status"] = "STALLED"
        out["b4_prime_predicted_branch"] = PREDICTED_BRANCH[bin_]
        out["b4_prime_verdict"] = "STALLED"
        log(f"  STALLED ({digits_stable_at_200:.1f} digits at n=200, gate=30)")
        return out

    out["L_estimate_30dp"] = mp.nstr(L[200], 30)
    out["L_estimate_50dp"] = mp.nstr(L[200], 50)

    # ---- WKB fit
    log(f"  WKB fit (dps={DPS_WKB}, N_ref={N_REF}, n in {WKB_NS[0]}..{WKB_NS[-1]} step 3)")
    wkb = wkb_fit(coeffs)
    out["wkb"] = wkb
    out["wkb_fit_run"] = True

    if wkb["A"] is None:
        out["status"] = "WKB_FAILED"
        out["b4_prime_predicted_branch"] = PREDICTED_BRANCH[bin_]
        out["b4_prime_verdict"] = "WKB_FAILED"
        log(f"    WKB fit FAILED ({wkb['n_points']} usable points)")
        return out

    A_fit = wkb["A"]
    A_se = wkb["A_stderr"]
    alpha_fit = wkb["alpha"]
    log(f"    A_fit = {A_fit:.4f} +/- {A_se:.4f}, alpha = {alpha_fit:.4f}, "
        f"beta = {wkb['beta']:.4f}, npts = {wkb['n_points']}")

    out["A_fit"] = round(A_fit, 4)
    out["A_stderr"] = round(A_se, 4)
    out["alpha_fit"] = round(alpha_fit, 4)
    pred = PREDICTED_BRANCH[bin_]
    out["b4_prime_predicted_branch"] = pred
    verdict = b4_prime_verdict(A_fit, pred)
    out["b4_prime_verdict"] = verdict
    out["status"] = "OK"
    log(f"    predicted branch = {pred}; verdict = {verdict}")
    return out


# ------------------------------------------------------------ main

def main():
    log("=== PCF-2 SESSION C1 :: non-CM cubic WKB harvest ===")
    log(f"DPS_L = {DPS_L}, DPS_WKB = {DPS_WKB}, "
        f"N_GRID_L = {N_GRID_L}, WKB_NS = [{WKB_NS[0]}..{WKB_NS[-1]} step 3]")

    cat = json.loads(CATALOGUE.read_text(encoding="utf-8"))
    fams = [f for f in cat["families"] if f.get("trichotomy_bin") in NON_CM_BINS]
    log(f"Loaded {len(fams)} families in non-CM bins "
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

        # HALT: violating
        if rec.get("b4_prime_verdict") == "b4_prime_violating":
            halts.append({
                "family_id": rec["family_id"],
                "bin": rec["trichotomy_bin"],
                "A_fit": rec.get("A_fit"),
                "A_stderr": rec.get("A_stderr"),
                "reason": ("b4_prime_violating: A_fit > 0.30 from BOTH 5 and 6, "
                           "falsifies Conjecture B4' (A in {2d-1, 2d}) for d=3"),
            })
            discrepancies.append({
                "family_id": rec["family_id"],
                "issue": "Conjecture B4' falsifier: A_fit outside {5, 6} band",
                "A_fit": rec.get("A_fit"),
            })

        # NON-HALT FLAGS
        if rec.get("status") == "OK":
            A = rec["A_fit"]
            d5 = abs(A - 5.0)
            d6 = abs(A - 6.0)
            if d6 < d5:
                # closer to 6 than to 5 in an elementary-positive bin
                unexpected.append({
                    "family_id": rec["family_id"],
                    "bin": rec["trichotomy_bin"],
                    "kind": ("non-CM family with A_fit closer to 6 than to 5: "
                             "would falsify B4'"),
                    "A_fit": A,
                    "A_stderr": rec.get("A_stderr"),
                    "predicted_branch": rec["b4_prime_predicted_branch"],
                })
            if rec.get("digits_stable_at_200", 0) >= 50.0:
                unexpected.append({
                    "family_id": rec["family_id"],
                    "bin": rec["trichotomy_bin"],
                    "kind": "striking convergence (>= 50 stable digits at n=200)",
                    "digits_stable_at_200": rec["digits_stable_at_200"],
                })

    # ---- branch-disagreement check
    branches = {bin_: set() for bin_ in NON_CM_BINS}
    for r in per_family:
        if r.get("status") == "OK":
            d5 = abs(r["A_fit"] - 5.0)
            d6 = abs(r["A_fit"] - 6.0)
            nearest = 5 if d5 <= d6 else 6
            branches[r["trichotomy_bin"]].add(nearest)
    if len({tuple(sorted(s)) for s in branches.values() if s}) > 1:
        unexpected.append({
            "kind": ("+_S3_real and +_C3_real give DIFFERENT empirical branches "
                     "-- B4' would need bin-level refinement"),
            "branches_per_bin": {k: sorted(list(v)) for k, v in branches.items()},
        })

    # ---- stalled-rate halt
    n_stalled = sum(1 for r in per_family if r.get("status") == "STALLED")
    if n_stalled >= max(1, int(0.25 * len(per_family))):
        halts.append({
            "reason": (f"stalled-rate halt: {n_stalled} of {len(per_family)} families "
                       f"have digits_stable_at_200 < 30 (>= 25% threshold)"),
            "n_stalled": n_stalled,
            "n_total": len(per_family),
        })

    # ---- aggregate
    A_vals = [r["A_fit"] for r in per_family if r.get("status") == "OK"]
    A_consistent_5 = [a for a in A_vals if abs(a - 5.0) <= 0.10]
    A_in_band = [a for a in A_vals if abs(a - 5.0) <= 0.10 or abs(a - 6.0) <= 0.10]
    summary = {
        "n_families": len(fams),
        "n_OK": sum(1 for r in per_family if r.get("status") == "OK"),
        "n_STALLED": n_stalled,
        "n_WKB_FAILED": sum(1 for r in per_family if r.get("status") == "WKB_FAILED"),
        "n_b4_prime_consistent": sum(1 for r in per_family if r.get("b4_prime_verdict") == "b4_prime_consistent"),
        "n_b4_prime_consistent_wrong_branch": sum(1 for r in per_family if r.get("b4_prime_verdict") == "b4_prime_consistent_wrong_branch"),
        "n_b4_prime_ambiguous": sum(1 for r in per_family if r.get("b4_prime_verdict") == "b4_prime_ambiguous"),
        "n_b4_prime_violating": sum(1 for r in per_family if r.get("b4_prime_verdict") == "b4_prime_violating"),
        "n_in_band_5_or_6_0p1": len(A_in_band),
        "n_consistent_with_predicted_5": len(A_consistent_5),
        "A_fit_mean": (float(np.mean(A_vals)) if A_vals else None),
        "A_fit_stddev": (float(np.std(A_vals)) if A_vals else None),
        "A_fit_min": (float(min(A_vals)) if A_vals else None),
        "A_fit_max": (float(max(A_vals)) if A_vals else None),
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    final = {
        "task_id": "PCF2-SESSION-C1",
        "date": "2026-05-01",
        "bins": list(NON_CM_BINS),
        "predicted_branch": "A = 2d - 1 = 5 (elementary-positive bins)",
        "summary": summary,
        "per_family": per_family,
    }

    RESULTS.write_text(json.dumps(final, indent=2, default=str), encoding="utf-8")
    HALT.write_text(json.dumps({"halts": halts} if halts else {}, indent=2), encoding="utf-8")
    DISC.write_text(json.dumps({"discrepancies": discrepancies} if discrepancies else {}, indent=2),
                    encoding="utf-8")
    UNEXP.write_text(json.dumps({"unexpected": unexpected} if unexpected else {}, indent=2),
                     encoding="utf-8")

    # ---- AEAL claims
    script_name = Path(__file__).name
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    claim_lines = []

    for r in per_family:
        if r.get("b4_prime_verdict") == "b4_prime_consistent":
            claim_lines.append({
                "claim": (f"Family {r['family_id']} (b = {r['b_latex']}, bin = "
                          f"{r['trichotomy_bin']}, Delta_3 = {r['Delta_3']}) has "
                          f"WKB exponent A_fit = {r['A_fit']:.4f} +/- {r['A_stderr']:.4f}, "
                          f"consistent with refined Conjecture B4' branch "
                          f"A = 2d-1 = {r['b4_prime_predicted_branch']} "
                          f"(elementary-positive cubic); fit at dps {DPS_WKB} over "
                          f"n in [{WKB_NS[0]}, {WKB_NS[-1]}] step 3"),
                "evidence_type": "computation",
                "dps": DPS_WKB,
                "reproducible": True,
                "script": script_name,
                "output_hash": script_hash,
            })

    n_consistent = summary["n_b4_prime_consistent"]
    claim_lines.append({
        "claim": (f"Conjecture B4' verified on {n_consistent + 37} of 50 cubic families "
                  f"(C1 contributes {n_consistent} of {len(fams)} non-CM, "
                  f"Session B contributed 37 of 37 -_S3_CM at A=2d=6 with mean 5.982 "
                  f"sigma 0.025); C1 mean A_fit = "
                  f"{summary['A_fit_mean']:.4f} (stddev "
                  f"{summary['A_fit_stddev']:.4f}, min {summary['A_fit_min']:.4f}, "
                  f"max {summary['A_fit_max']:.4f}) on the +_S3_real and +_C3_real bins, "
                  f"predicted branch A = 2d-1 = 5"),
        "evidence_type": "computation",
        "dps": DPS_WKB,
        "reproducible": True,
        "script": script_name,
        "output_hash": script_hash,
    })

    with open(CLAIMS, "w", encoding="utf-8") as fh:
        for c in claim_lines:
            fh.write(json.dumps(c) + "\n")

    # ---- LaTeX tables
    write_noncm_table(per_family)
    write_full_50_table(per_family, cat)

    log("=== SUMMARY ===")
    log(json.dumps(summary, indent=2))
    log(f"AEAL claims written: {len(claim_lines)}")
    log(f"Halts: {len(halts)}; discrepancies: {len(discrepancies)}; "
        f"unexpected: {len(unexpected)}")
    log("Done.")


# ------------------------------------------------------------ tables

def _bin_tex(bin_: str) -> str:
    return (bin_
            .replace("+_S3_real", r"$+\_S_3\_\mathrm{real}$")
            .replace("+_C3_real", r"$+\_C_3\_\mathrm{real}$")
            .replace("-_S3_CM", r"$-\_S_3\_\mathrm{CM}$"))


def _coeff_tex(c) -> str:
    return f"({c[0]},{c[1]},{c[2]},{c[3]})"


def write_noncm_table(per_family):
    lines = []
    lines.append(r"% PCF-2 Session C1 -- non-CM cubic WKB harvest")
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\small")
    lines.append(r"\caption{WKB-exponent harvest on the elementary-positive "
                 r"cubic bins ($+\_S_3\_\mathrm{real}$ and $+\_C_3\_\mathrm{real}$, "
                 r"13 families). $A_{\mathrm{fit}}$ from "
                 r"$\log|\delta_n| \sim -A\, n\log n + \alpha n - \beta\log n + \gamma$ "
                 r"fit at dps 800 over $n\in[10,100]$ step 3, $N_{\mathrm{ref}} = 300$. "
                 r"`pred.\ branch' is Conjecture~B4$'$ prediction $A = 2d - 1 = 5$.}")
    lines.append(r"\label{tab:pcf2-sessionC1-wkb}")
    lines.append(r"\begin{tabular}{rlrlllll}")
    lines.append(r"\hline")
    lines.append(r"\textbf{ID} & \textbf{$(\alpha_3,\alpha_2,\alpha_1,\alpha_0)$} & "
                 r"\textbf{$\Delta_3$} & \textbf{bin} & "
                 r"\textbf{$A_{\mathrm{fit}}$} & \textbf{stderr} & "
                 r"\textbf{pred.\ branch} & \textbf{verdict} \\ \hline")
    for r in per_family:
        fid = r["family_id"]
        c = r["b_coefficients_a3_a2_a1_a0"]
        d3 = r["Delta_3"]
        b_tex = _bin_tex(r["trichotomy_bin"])
        if r.get("status") == "OK":
            A = f"{r['A_fit']:.4f}"
            se = f"{r['A_stderr']:.4f}"
        else:
            A = "n/a"
            se = "n/a"
        pred = r.get("b4_prime_predicted_branch", "?")
        verdict = r.get("b4_prime_verdict", "?").replace("_", r"\_")
        lines.append(f"{fid} & {_coeff_tex(c)} & {d3} & {b_tex} & "
                     f"{A} & {se} & {pred} & {verdict} \\\\")
    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    TBL_NONCM.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_full_50_table(per_family, cat):
    """Combined 50-row WKB table: Session B (-_S3_CM, predicted A=6) + C1 non-CM (predicted A=5)."""
    sessB_results_path = HERE.parent / "PCF2-SESSION-B" / "results.json"
    sessB = json.loads(sessB_results_path.read_text(encoding="utf-8"))
    sessB_by_id = {r["family_id"]: r for r in sessB["per_family"]}
    sessC1_by_id = {r["family_id"]: r for r in per_family}

    rows = []
    for fam in cat["families"]:
        fid = fam["family_id"]
        bin_ = fam.get("trichotomy_bin", "?")
        coeffs = (fam.get("alpha_3"), fam.get("alpha_2"),
                  fam.get("alpha_1"), fam.get("alpha_0"))
        d3 = fam.get("Delta_3_exact")
        if bin_ == "-_S3_CM":
            r = sessB_by_id.get(fid)
            wkb = (r or {}).get("wkb", {}) if r else {}
            A = wkb.get("A")
            se = wkb.get("A_stderr")
            pred = 6
        elif bin_ in NON_CM_BINS:
            r = sessC1_by_id.get(fid)
            A = (r or {}).get("A_fit") if r else None
            se = (r or {}).get("A_stderr") if r else None
            pred = 5
        else:
            # reducible-control or other bins -- no scope
            A = None
            se = None
            pred = None
        if A is None:
            A_str, se_str = "n/a", "n/a"
            verdict = "out-of-scope" if pred is None else "n/a"
        else:
            A_str = f"{A:.4f}"
            se_str = f"{se:.4f}" if se is not None else "n/a"
            d_pred = abs(A - float(pred))
            d5 = abs(A - 5.0)
            d6 = abs(A - 6.0)
            if d_pred <= 0.10:
                verdict = r"$\checkmark$"
            elif min(d5, d6) > 0.30:
                verdict = "VIOLATE"
            else:
                verdict = "ambiguous"
        rows.append((fid, coeffs, bin_, d3, A_str, se_str, pred, verdict))

    lines = []
    lines.append(r"% PCF-2 Session C1 -- combined 50-row WKB cubic harvest")
    lines.append(r"% Replaces sessions/2026-05-01/PCF2-SESSION-B/wkb_cubic_harvest.tex")
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\scriptsize")
    lines.append(r"\caption{Combined 50-row WKB-exponent cubic harvest "
                 r"(Sessions B + C1). $A_{\mathrm{fit}}$ fitted from "
                 r"$\log|\delta_n| \sim -A\, n\log n + \alpha n - \beta\log n + \gamma$ "
                 r"at dps 800. `pred.' is the Conjecture~B4$'$ branch: "
                 r"$2d - 1 = 5$ for the elementary-positive bins "
                 r"$+\_S_3\_\mathrm{real}$ and $+\_C_3\_\mathrm{real}$, "
                 r"$2d = 6$ for $-\_S_3\_\mathrm{CM}$. "
                 r"Out-of-scope rows are reducible controls or non-cubic.}")
    lines.append(r"\label{tab:pcf2-cubic-harvest-v2}")
    lines.append(r"\begin{tabular}{rlllrlll}")
    lines.append(r"\hline")
    lines.append(r"\textbf{ID} & \textbf{$(\alpha_3,\alpha_2,\alpha_1,\alpha_0)$} & "
                 r"\textbf{bin} & \textbf{$\Delta_3$} & "
                 r"\textbf{pred.} & \textbf{$A_{\mathrm{fit}}$} & "
                 r"\textbf{stderr} & \textbf{verdict} \\ \hline")
    for fid, coeffs, bin_, d3, A_str, se_str, pred, verdict in rows:
        c = coeffs if coeffs[0] is not None else (None, None, None, None)
        ctex = _coeff_tex(c) if c[0] is not None else "n/a"
        bt = _bin_tex(bin_)
        pred_str = str(pred) if pred is not None else "--"
        lines.append(f"{fid} & {ctex} & {bt} & {d3 if d3 is not None else 'n/a'} & "
                     f"{pred_str} & {A_str} & {se_str} & {verdict} \\\\")
    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    TBL_FULL.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
