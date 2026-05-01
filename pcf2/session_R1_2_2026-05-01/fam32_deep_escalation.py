"""Deep WKB escalation for the j=0 quartic candidate (family 32) at d=4.
Pushes N_ref and fit window much further to distinguish "slow convergence
to 0" from "B5 genuinely does not extend to d=4".

Family 32: b(n) = n^4 - 3n^3 - 3n^2 + 3n - 3, with I=0 (j(Jac)=0).
"""
from __future__ import annotations
import json, math, time
from pathlib import Path
import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
COEFFS = [1, -3, -3, 3, -3]  # [a4,a3,a2,a1,a0]


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


def wkb_fit(coeffs, dps, N_grid, N_ref):
    with mp.workdps(dps):
        L_ref = cf_value(coeffs, N_ref, dps)
        rows_x, rows_y = [], []
        for N in N_grid:
            LN = cf_value(coeffs, N, dps)
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
        "fit_window_N": [N_grid[0], N_grid[-1]],
        "residual_rms": float(math.sqrt(s2)),
        "dps": dps, "N_ref": N_ref,
    }


def main():
    log = []
    def L(s):
        line = f"[{time.strftime('%H:%M:%S')}] {s}"
        print(line, flush=True); log.append(line)

    L("Deep WKB escalation for family 32 (j=0 quartic).")
    runs = []
    for dps, Nmax, Nref in [
        (1500, 130, 200),  # baseline
        (1500, 200, 300),
        (2000, 250, 400),
    ]:
        N_grid = list(range(10, Nmax + 1, 2))
        t0 = time.time()
        res = wkb_fit(COEFFS, dps, N_grid, Nref)
        dt = time.time() - t0
        delta = res["A"] - 8.0
        L(f"  dps={dps}, N=[{N_grid[0]}..{N_grid[-1]}], N_ref={Nref}: "
          f"A={res['A']:.6f} +/- {res['A_stderr']:.2e}, "
          f"delta={delta:+.4e}, dt={dt:.1f}s")
        res["delta"] = delta
        res["dt_seconds"] = dt
        runs.append(res)

    # Also: tail-window fit (drop small N to suppress sub-leading)
    L("Tail-window-only fit (drop N < 50)")
    for dps, Nmax, Nref in [(2000, 250, 400)]:
        N_grid = list(range(50, Nmax + 1, 2))
        t0 = time.time()
        res = wkb_fit(COEFFS, dps, N_grid, Nref)
        dt = time.time() - t0
        delta = res["A"] - 8.0
        L(f"  TAIL dps={dps}, N=[{N_grid[0]}..{N_grid[-1]}], N_ref={Nref}: "
          f"A={res['A']:.6f} +/- {res['A_stderr']:.2e}, "
          f"delta={delta:+.4e}, dt={dt:.1f}s")
        res["delta"] = delta
        res["tail_window"] = True
        runs.append(res)

    out = HERE / "fam32_deep_escalation.json"
    with open(out, "w") as f:
        json.dump({"family_id": 32, "coeffs": COEFFS, "runs": runs,
                   "log": log}, f, indent=2)
    L(f"  wrote {out.name}")


if __name__ == "__main__":
    main()
