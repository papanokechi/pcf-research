"""
PCF-2 SESSION Q1 -- WKB harvest at d = 4 (Phase Q1-2)

For every family in the quartic catalogue, compute L_b stability and fit
    log|delta_n| ~ -A n log n + alpha n - beta log n + gamma
to test Conjecture B4 at d = 4 (claim Q1-A): predicted A = 2d = 8.

Outputs (in this directory):
    results.json
    discrepancy_log.json
    unexpected_finds.json
    halt_log.json (only on halt)
    claims_phaseQ1_2.jsonl
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
CATALOGUE = HERE / "quartic_family_catalogue.json"
RESULTS = HERE / "results.json"
DISC = HERE / "discrepancy_log.json"
UNEXP = HERE / "unexpected_finds.json"
HALT = HERE / "halt_log.json"
CLAIMS = HERE / "claims_phaseQ1_2.jsonl"
RUN_LOG = HERE / "wkb_run.log"

if RUN_LOG.exists():
    RUN_LOG.unlink()


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(RUN_LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


# ---------- params

DPS_STAB = 80
DPS_WKB = 1200          # working precision for WKB delta extraction
N_REF = 150             # reference N for delta_n := |L_N - L_ref|
N_GRID_LO, N_GRID_HI = 10, 60
N_GRID_STEP = 2

# B4 verdict thresholds
B4_TARGET = 8.0
B4_BAND_TIGHT = 0.10
B4_BAND_LOOSE = 0.30

# Halt threshold for publishable counterexample:
HALT_AFIT_DIST = 0.30
HALT_STDERR = 0.05


# ---------- CF eval (quartic)

def cf_value(coeffs, N: int, dps: int) -> mp.mpf:
    a4, a3, a2, a1, a0 = coeffs
    with mp.workdps(dps):
        b_N = (mp.mpf(a4) * N ** 4 + mp.mpf(a3) * N ** 3
               + mp.mpf(a2) * N ** 2 + mp.mpf(a1) * N + mp.mpf(a0))
        x = b_N
        for k in range(N - 1, -1, -1):
            bk = (mp.mpf(a4) * k ** 4 + mp.mpf(a3) * k ** 3
                  + mp.mpf(a2) * k ** 2 + mp.mpf(a1) * k + mp.mpf(a0))
            x = bk + mp.mpf(1) / x
        return +x


def stable_digits(L_a: mp.mpf, L_b: mp.mpf) -> float:
    with mp.workdps(80):
        d = abs(L_a - L_b)
        if d == 0:
            return 999.0
        return float(-mp.log10(d))


# ---------- WKB fit

def wkb_fit(coeffs, n_ref: int = N_REF, dps: int = DPS_WKB):
    n_grid = list(range(N_GRID_LO, N_GRID_HI + 1, N_GRID_STEP))
    with mp.workdps(dps):
        L_ref = cf_value(coeffs, n_ref, dps)
        rows_x = []
        rows_y = []
        ns_used = []
        for N in n_grid:
            LN = cf_value(coeffs, N, dps)
            d = abs(LN - L_ref)
            if d == 0:
                continue
            y = float(mp.log(d))
            rows_x.append([-N * math.log(N), float(N), -math.log(N), 1.0])
            rows_y.append(y)
            ns_used.append(N)
    if len(rows_y) < 6:
        return {"A": None, "alpha": None, "beta": None, "gamma": None,
                "A_stderr": None, "n_points": len(rows_y),
                "ns_used": ns_used,
                "comment": "too few usable points"}
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
        "fit_window_N": [N_GRID_LO, N_GRID_HI],
        "ns_used": ns_used,
        "residual_rms": float(math.sqrt(s2)),
    }


def b4_verdict(A_fit, A_stderr) -> str:
    if A_fit is None:
        return "STALLED"
    d8 = abs(A_fit - 8.0)
    d7 = abs(A_fit - 7.0)
    d6 = abs(A_fit - 6.0)
    d5 = abs(A_fit - 5.0)
    if d8 <= B4_BAND_TIGHT:
        return "b4_consistent_at_8"
    if d7 <= B4_BAND_TIGHT and d8 > B4_BAND_LOOSE:
        return "b4_split_candidate_at_7"
    # outside both, check far from {5,6,7,8}
    if (d8 > B4_BAND_LOOSE and d7 > B4_BAND_LOOSE
            and d6 > B4_BAND_LOOSE and d5 > B4_BAND_LOOSE):
        return "b4_violating"
    if d8 > B4_BAND_TIGHT and d7 > B4_BAND_TIGHT:
        return "b4_ambiguous"
    return "b4_borderline"


# ---------- per-family

def analyze_family(fam: dict) -> dict:
    fid = fam["family_id"]
    coeffs = (fam["alpha_4"], fam["alpha_3"], fam["alpha_2"],
              fam["alpha_1"], fam["alpha_0"])
    log(f"--- family {fid}: b = {fam['b_latex']}, D4 = {fam['Delta_4_exact']} "
        f"Gal = {fam['Galois_group']} bin = {fam['trichotomy_bin']} ---")

    out = {
        "family_id": fid,
        "b_coefficients": list(coeffs),
        "b_latex": fam["b_latex"],
        "Delta_4": fam["Delta_4_exact"],
        "Galois_group": fam["Galois_group"],
        "trichotomy_bin": fam["trichotomy_bin"],
    }

    # L stability at low dps
    L = {}
    for N in (50, 100, 150):
        L[N] = cf_value(coeffs, N, DPS_STAB)
    digs_50_100 = stable_digits(L[50], L[100])
    digs_100_150 = stable_digits(L[100], L[150])
    out["L_at_N"] = {str(N): mp.nstr(L[N], 30) for N in L}
    out["digits_stable"] = {
        "delta_50_100_log10": digs_50_100,
        "delta_100_150_log10": digs_100_150,
    }
    out["L_estimate_30dp"] = mp.nstr(L[150], 30)
    log(f"  L stability: 50/100={digs_50_100:.1f}  100/150={digs_100_150:.1f}")

    # WKB fit
    log(f"  WKB fit (dps={DPS_WKB}, N_ref={N_REF}, N in [{N_GRID_LO},{N_GRID_HI}] step {N_GRID_STEP})")
    t1 = time.time()
    wkb = wkb_fit(coeffs)
    out["wkb"] = wkb
    if wkb["A"] is not None:
        log(f"    A_fit = {wkb['A']:.4f} +/- {wkb['A_stderr']:.4f}, "
            f"alpha = {wkb['alpha']:.3f}, beta = {wkb['beta']:.3f}, "
            f"npts = {wkb['n_points']}, t = {time.time()-t1:.1f}s")
    else:
        log(f"    WKB FAILED ({wkb.get('comment', '')})")

    A = wkb.get("A")
    se = wkb.get("A_stderr")
    out["b4_verdict"] = b4_verdict(A, se)
    log(f"    verdict: {out['b4_verdict']}")
    return out


# ---------- driver

def main():
    log("=== PCF-2 SESSION Q1 :: Phase Q1-2 WKB harvest at d=4 ===")
    cat = json.loads(CATALOGUE.read_text(encoding="utf-8"))
    fams = cat["families"]
    log(f"Loaded {len(fams)} families")

    per_family = []
    halts = []
    discrepancies = []
    unexpected = []
    stalls = 0

    t0 = time.time()
    for i, fam in enumerate(fams, 1):
        log(f"[{i}/{len(fams)}] family {fam['family_id']} elapsed={time.time()-t0:.1f}s")
        try:
            rec = analyze_family(fam)
        except Exception as exc:
            log(f"  EXCEPTION: {exc}")
            rec = {
                "family_id": fam["family_id"],
                "b_coefficients": [fam.get(k) for k in
                                   ("alpha_4", "alpha_3", "alpha_2", "alpha_1", "alpha_0")],
                "b_latex": fam.get("b_latex"),
                "Delta_4": fam.get("Delta_4_exact"),
                "Galois_group": fam.get("Galois_group"),
                "trichotomy_bin": fam.get("trichotomy_bin"),
                "wkb": {"A": None, "comment": f"exception: {exc}"},
                "b4_verdict": "STALLED",
            }
        per_family.append(rec)
        if rec["b4_verdict"] == "STALLED":
            stalls += 1
        if rec["b4_verdict"] == "b4_violating":
            wkb = rec.get("wkb", {})
            if (wkb.get("A_stderr") is not None
                    and wkb["A_stderr"] < HALT_STDERR
                    and abs(wkb["A"] - 7.0) > HALT_AFIT_DIST
                    and abs(wkb["A"] - 8.0) > HALT_AFIT_DIST):
                halts.append({
                    "family_id": rec["family_id"],
                    "reason": "publishable B4 counterexample candidate",
                    "A_fit": wkb["A"], "A_stderr": wkb["A_stderr"],
                })
                discrepancies.append({
                    "family_id": rec["family_id"],
                    "issue": "Conjecture B4 counterexample candidate at d=4",
                    "A_fit": wkb["A"], "A_stderr": wkb["A_stderr"],
                    "b_latex": rec["b_latex"],
                })

        # non-halt flags
        wkb = rec.get("wkb", {})
        if wkb.get("A") is not None:
            A = wkb["A"]
            if rec["b4_verdict"] == "b4_split_candidate_at_7":
                unexpected.append({
                    "family_id": rec["family_id"],
                    "kind": "b4_split_candidate at A_fit ~ 7",
                    "A_fit": round(A, 4),
                    "Galois": rec["Galois_group"],
                    "bin": rec["trichotomy_bin"],
                })
            elif rec["b4_verdict"] not in ("b4_consistent_at_8", "b4_violating"):
                unexpected.append({
                    "family_id": rec["family_id"],
                    "kind": f"verdict={rec['b4_verdict']}",
                    "A_fit": round(A, 4),
                    "A_stderr": round(wkb["A_stderr"], 4),
                })

    # WKB summary
    A_vals = [r["wkb"]["A"] for r in per_family if r.get("wkb", {}).get("A") is not None]
    A_at_8 = [a for a in A_vals if abs(a - 8.0) <= B4_BAND_TIGHT]
    A_at_7 = [a for a in A_vals if abs(a - 7.0) <= B4_BAND_TIGHT]

    verdict_counts: dict[str, int] = {}
    for r in per_family:
        verdict_counts[r["b4_verdict"]] = verdict_counts.get(r["b4_verdict"], 0) + 1

    summary = {
        "n_families": len(fams),
        "verdict_counts": verdict_counts,
        "n_with_WKB_fit": len(A_vals),
        "n_within_0p1_of_8": len(A_at_8),
        "n_within_0p1_of_7": len(A_at_7),
        "WKB_A_mean": (float(np.mean(A_vals)) if A_vals else None),
        "WKB_A_stddev": (float(np.std(A_vals)) if A_vals else None),
        "WKB_A_min": (float(min(A_vals)) if A_vals else None),
        "WKB_A_max": (float(max(A_vals)) if A_vals else None),
        "stall_fraction": stalls / max(len(fams), 1),
        "elapsed_seconds": round(time.time() - t0, 1),
    }

    final = {
        "task_id": "PCF2-SESSION-Q1",
        "phase": "Q1-2",
        "date": "2026-05-01",
        "summary": summary,
        "per_family": per_family,
    }

    RESULTS.write_text(json.dumps(final, indent=2, default=str), encoding="utf-8")
    DISC.write_text(json.dumps({"discrepancies": discrepancies} if discrepancies else {},
                               indent=2), encoding="utf-8")
    UNEXP.write_text(json.dumps({"unexpected": unexpected} if unexpected else {},
                                indent=2), encoding="utf-8")

    if halts or summary["stall_fraction"] >= 0.25:
        halt_payload = {
            "halts": halts,
            "stall_fraction": summary["stall_fraction"],
            "stall_threshold": 0.25,
        }
        HALT.write_text(json.dumps(halt_payload, indent=2), encoding="utf-8")
        log(f"  HALT WRITTEN: {halt_payload}")

    # claims
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    claims_lines = []
    claims_lines.append({
        "claim": (f"WKB-quartic harvest (60 cataloged families): {len(A_at_8)}/{len(A_vals)} "
                  f"within 0.10 of A=2d=8; mean={summary['WKB_A_mean']:.4f}, "
                  f"sigma={summary['WKB_A_stddev']:.4f}, "
                  f"min={summary['WKB_A_min']:.4f}, max={summary['WKB_A_max']:.4f}; "
                  f"verdict counts: {verdict_counts}"),
        "evidence_type": "computation", "dps": DPS_WKB, "reproducible": True,
        "script": "session_q1_wkb.py", "output_hash": script_hash,
    })
    if halts:
        for h in halts:
            claims_lines.append({
                "claim": (f"Family {h['family_id']} candidate B4 counterexample at d=4: "
                          f"A_fit={h['A_fit']:.4f} +/- {h['A_stderr']:.4f}"),
                "evidence_type": "computation", "dps": DPS_WKB, "reproducible": True,
                "script": "session_q1_wkb.py", "output_hash": script_hash,
            })
    with open(CLAIMS, "w", encoding="utf-8") as f:
        for c in claims_lines:
            f.write(json.dumps(c, ensure_ascii=False) + "\n")

    log("=" * 72)
    log(f"Verdict counts: {verdict_counts}")
    log(f"WKB stats: mean={summary['WKB_A_mean']:.4f}, sigma={summary['WKB_A_stddev']:.4f}, "
        f"min={summary['WKB_A_min']:.4f}, max={summary['WKB_A_max']:.4f}")
    log(f"Stall fraction: {summary['stall_fraction']:.2%}")
    log(f"Wrote: {RESULTS.name}, {DISC.name}, {UNEXP.name}, {CLAIMS.name}")


if __name__ == "__main__":
    main()
