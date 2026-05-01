"""PCF-2 SESSION R1.3 -- Phase R13-D: family-32 deep WKB

Family 32: b(n) = n^4 - 3 n^3 - 3 n^2 + 3 n - 3 (the unique j=0 quartic
in the Q1 catalogue).

Run a single-family deep-WKB with:
    dps    = 5000
    N_ref  = 1000
    fit window N in [200, 800] step 20  (=> 31 points)
Compute:
    - FIXED-A=8 fit          -> residual_mean_deep, residual_at_max_n_deep,
                                 residual_std_deep
    - FREE-A 4-parameter fit -> A_fit_deep, delta_4_deep, A_stderr_deep
For comparison we also load the corresponding tail-fit of Q1 family 1
(non-j=0 representative) at the same deep settings, to give a deep-N
non-j=0 cluster baseline.
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
LOG = HERE / "r1_3_family32_deep.log"

DPS = 5000
N_GRID = list(range(200, 801, 20))   # 200, 220, ..., 800  -> 31 pts
N_REF = 1000

FAM32_COEFFS = [1, -3, -3, 3, -3]   # j(Jac) = 0
# A non-j=0 deep baseline: Q1 family 1 = [1, -3, -3, -3, -3]
FAM_BASELINE_COEFFS = [1, -3, -3, -3, -3]


def L(s):
    line = f"[{time.strftime('%H:%M:%S')}] {s}"
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def cf_value(coeffs, N, dps):
    with mp.workdps(dps):
        ms = [mp.mpf(c) for c in coeffs]

        def b(k):
            v = ms[0]; kk = mp.mpf(k)
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


def deep_fit(coeffs, label):
    L(f"  computing y_n for {label} (coeffs={coeffs}, dps={DPS}, "
      f"N_ref={N_REF}, N=[{N_GRID[0]}..{N_GRID[-1]}] step 20)")
    t0 = time.time()
    ns, ys = compute_y(coeffs, N_GRID, N_REF, DPS)
    dt = time.time() - t0
    L(f"  y_n done in {dt:.1f}s, n_points={len(ys)}")

    # FIXED-A=8 fit
    A_fix = 8.0
    z = ys + A_fix * ns * np.log(ns)
    X = np.column_stack([ns, -np.log(ns), np.ones_like(ns)])
    sol, *_ = np.linalg.lstsq(X, z, rcond=None)
    alpha_fix, beta_fix, gamma_fix = (float(c) for c in sol)
    resid_fix = z - X @ sol
    fixed = {
        "A_fixed": A_fix,
        "alpha": alpha_fix, "beta": beta_fix, "gamma": gamma_fix,
        "n_points": int(len(ys)),
        "residual_mean": float(resid_fix.mean()),
        "residual_std": float(resid_fix.std(ddof=1)) if len(resid_fix) > 1 else 0.0,
        "residual_at_max_n": float(resid_fix[-1]),
        "residual_at_n_ref": float(resid_fix[-1]),  # closest to N_ref=1000 we have
        "residuals": resid_fix.tolist(),
        "ns": ns.tolist(),
    }

    # FREE-A 4-param fit
    Xf = np.column_stack([-ns * np.log(ns), ns, -np.log(ns), np.ones_like(ns)])
    sol2, *_ = np.linalg.lstsq(Xf, ys, rcond=None)
    A, alphaf, betaf, gammaf = (float(c) for c in sol2)
    residf = ys - Xf @ sol2
    s2 = float(np.sum(residf ** 2) / max(len(ys) - 4, 1))
    cov = s2 * np.linalg.inv(Xf.T @ Xf)
    se = np.sqrt(np.diag(cov))
    free = {"A": A, "alpha": alphaf, "beta": betaf, "gamma": gammaf,
            "A_stderr": float(se[0]),
            "n_points": int(len(ys)),
            "residual_rms": float(math.sqrt(s2))}

    L(f"  FIXED A=8: alpha={alpha_fix:+.6f}, beta={beta_fix:+.6f}, "
      f"gamma={gamma_fix:+.6f}")
    L(f"  FIXED resid: mean={fixed['residual_mean']:+.4e}, "
      f"std={fixed['residual_std']:.4e}, at_max_n(N=800)={fixed['residual_at_max_n']:+.4e}")
    L(f"  FREE 4-param: A={A:.7f} +/- {se[0]:.2e}, delta=A-8={A-8.0:+.4e}")

    return {
        "label": label, "coeffs": coeffs,
        "fixed": fixed, "free": free,
        "delta_deep": A - 8.0,
        "n_points": int(len(ys)),
        "dt_seconds": dt,
    }


def main():
    if LOG.exists():
        LOG.unlink()
    L(f"R13-D family-32 deep-WKB: dps={DPS}, N_ref={N_REF}, "
      f"fit window N=[{N_GRID[0]}..{N_GRID[-1]}] step 20")

    fam32 = deep_fit(FAM32_COEFFS, "fam32_jZero")
    base = deep_fit(FAM_BASELINE_COEFFS, "fam01_baseline_nonjZero")

    out = {
        "task_id": "PCF2-SESSION-R1.3",
        "phase": "D",
        "config": {"dps": DPS, "N_grid": [N_GRID[0], N_GRID[-1], 20],
                   "N_ref": N_REF},
        "fam32": fam32,
        "fam01_baseline": base,
        "shallow_R13B_fam32": {
            "delta_R13_free": -3.7081e-3,
            "residual_at_max_n_N250": +2.2844e-3,
            "residual_std_N250": 1.2774e-3,
            "source": "Phase R13-B fixed-A=8 fit, dps=2000, N=[50,250]",
        },
        "shallow_R13B_fam01": {
            "delta_R13_free": -3.8265e-3,
            "residual_at_max_n_N250": +4.6489e-3,  # not sure; placeholder
            "source": "R1.2 tail_fit_overrides fam 1, dps=2000, N=[50,250]",
        },
    }

    # save residuals as separate JSON
    json_path = HERE / "family32_deep_residual.json"
    with open(json_path, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    L(f"wrote {json_path.name}")

    # claims
    h = hashlib.sha256(json_path.read_bytes()).hexdigest()
    claims = []
    def cl(text):
        claims.append({"claim": text, "evidence_type": "computation",
                       "dps": DPS, "reproducible": True,
                       "script": "r1_3_family32_deep.py",
                       "output_hash": h})

    cl(f"R13-D: family-32 deep WKB (dps={DPS}, N_ref={N_REF}, "
       f"window=[{N_GRID[0]},{N_GRID[-1]}]) FREE 4-param A_fit = "
       f"{fam32['free']['A']:.7f} +/- {fam32['free']['A_stderr']:.2e}, "
       f"delta_deep = {fam32['delta_deep']:+.4e}")
    cl(f"R13-D: family-32 deep WKB FIXED A=8 residual_at_max_n(N=800) = "
       f"{fam32['fixed']['residual_at_max_n']:+.4e}, residual_std = "
       f"{fam32['fixed']['residual_std']:.4e}, residual_mean = "
       f"{fam32['fixed']['residual_mean']:+.4e}")
    cl(f"R13-D: non-j=0 baseline (Q1 family 1, b=n^4-3n^3-3n^2-3n-3) deep "
       f"WKB FREE A_fit = {base['free']['A']:.7f} +/- "
       f"{base['free']['A_stderr']:.2e}, delta_deep = "
       f"{base['delta_deep']:+.4e}")
    cl(f"R13-D: comparison of fam32 vs fam01 deep delta_deep: "
       f"{fam32['delta_deep']:+.4e} vs {base['delta_deep']:+.4e}; "
       f"diff = {fam32['delta_deep'] - base['delta_deep']:+.4e}")
    cl(f"R13-D: deep-WKB shrinks fam32 |delta| from "
       f"{abs(-3.7081e-3):.4e} (R13-B shallow) to "
       f"{abs(fam32['delta_deep']):.4e} (deep)")

    out_jsonl = HERE / "claims_phase_D.jsonl"
    with open(out_jsonl, "w") as fh:
        for c in claims:
            fh.write(json.dumps(c) + "\n")
    L(f"wrote {len(claims)} claims to {out_jsonl.name}")


if __name__ == "__main__":
    main()
