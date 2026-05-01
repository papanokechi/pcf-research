"""Tail-window WKB fit for all 60 Q1 quartic families, to mitigate
the fit-window bias diagnosed in fam32_deep_escalation.

Uses dps=2000, N=range(50, 252, 2), N_ref=400. ~3s/family => ~3 min total.
"""
from __future__ import annotations
import json, math, time
from pathlib import Path
import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
Q1_DIR = ROOT / "sessions" / "2026-05-01" / "PCF2-SESSION-Q1"

DPS = 2000
N_GRID = list(range(50, 252, 2))
N_REF = 400


def cf_value(coeffs, N, dps):
    a4, a3, a2, a1, a0 = coeffs
    with mp.workdps(dps):
        a4m = mp.mpf(a4); a3m = mp.mpf(a3); a2m = mp.mpf(a2)
        a1m = mp.mpf(a1); a0m = mp.mpf(a0)
        x = a4m * N**4 + a3m * N**3 + a2m * N**2 + a1m * N + a0m
        for k in range(N - 1, -1, -1):
            bk = a4m * k**4 + a3m * k**3 + a2m * k**2 + a1m * k + a0m
            x = bk + mp.mpf(1) / x
        return +x


def wkb_tail_fit(coeffs):
    with mp.workdps(DPS):
        L_ref = cf_value(coeffs, N_REF, DPS)
        rows_x, rows_y = [], []
        for N in N_GRID:
            LN = cf_value(coeffs, N, DPS)
            d = abs(LN - L_ref)
            if d == 0:
                continue
            y = float(mp.log(d))
            rows_x.append([-N * math.log(N), float(N), -math.log(N), 1.0])
            rows_y.append(y)
    X = np.array(rows_x); y = np.array(rows_y)
    coeffs_lsq, *_ = np.linalg.lstsq(X, y, rcond=None)
    A, alpha, beta, gamma = (float(c) for c in coeffs_lsq)
    yhat = X @ coeffs_lsq
    resid = y - yhat
    dof = max(len(y) - 4, 1)
    s2 = float(np.sum(resid**2) / dof)
    XtX_inv = np.linalg.inv(X.T @ X)
    cov = s2 * XtX_inv
    stderrs = np.sqrt(np.diag(cov))
    return {
        "A": A, "alpha": alpha, "beta": beta, "gamma": gamma,
        "A_stderr": float(stderrs[0]),
        "n_points": len(y),
        "fit_window_N": [N_GRID[0], N_GRID[-1]],
        "residual_rms": float(math.sqrt(s2)),
        "dps": DPS, "N_ref": N_REF,
    }


def main():
    log_file = HERE / "tail_fit_run.log"
    log_file.write_text("", encoding="utf-8")
    def L(s):
        line = f"[{time.strftime('%H:%M:%S')}] {s}"
        print(line, flush=True)
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(line + "\n")

    L(f"Tail-window WKB fit for 60 Q1 quartics: dps={DPS}, N=[{N_GRID[0]}..{N_GRID[-1]}], N_ref={N_REF}")
    q1 = json.load(open(Q1_DIR / "results.json"))
    out = {}
    t_start = time.time()
    for ff in q1["per_family"]:
        fid = ff["family_id"]
        coeffs = ff["b_coefficients"]
        t0 = time.time()
        res = wkb_tail_fit(coeffs)
        dt = time.time() - t0
        delta = res["A"] - 8.0
        L(f"  fam {fid:>2}: A={res['A']:.6f} +/- {res['A_stderr']:.2e}, "
          f"delta={delta:+.4e}, dt={dt:.1f}s")
        res["delta"] = delta
        res["dt_seconds"] = dt
        out[fid] = res
    L(f"Total elapsed: {time.time()-t_start:.1f}s")

    out_path = HERE / "tail_fit_overrides.json"
    with open(out_path, "w") as f:
        json.dump({
            "config": {"dps": DPS, "N_grid": N_GRID, "N_ref": N_REF},
            "overrides": {str(k): v for k, v in out.items()},
        }, f, indent=2)
    L(f"  wrote {out_path.name}")


if __name__ == "__main__":
    main()
