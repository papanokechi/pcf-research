"""PCF-2 SESSION R1.3 -- Phase R13-A and R13-B: fixed-A residualization

For every cubic family (R1.1 dataset, 50 rows) and every Q1 quartic family
(60 rows), we recompute the WKB tail data
    y_n := log|L_N - L_{N_ref}|
at dps=2000, N in N_GRID, N_ref=400, and then fit the FIXED-A model
    y_n = -A * n * log(n) + alpha * n - beta * log(n) + gamma + residual(n)
with A := 2 * d held fixed.

The fitted [alpha, beta, gamma] absorb the leading sub-leading WKB
shape; the residual sequence is what survives.  Per family we report
    residual_mean := mean residual over the fit window
    residual_std  := std of residuals
    residual_at_n_ref := residual at the largest n in the window
                         (used as a tail-extrapolation proxy; we cannot
                          evaluate residual at exactly N=N_ref because
                          y_{N_ref} is undefined when L_ref = L_{N_ref})

We also do an UNCONSTRAINED 4-parameter fit (A free) on the same y data
to extract A_fit_R13 and delta_d := A_fit_R13 - 2d, for cross-comparison
with R1.1 (cubic) and R1.2 tail_fit_overrides (quartic).

For Spearman correlation tests we use scipy.stats.spearmanr.

OUTPUTS:
    residualization_d3.csv
    residualization_d4.csv
    results_v3_phase_AB.json   (statistics, written; merged into results_v3 later)
    claims_phase_AB.jsonl
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
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
SES_A = ROOT / "sessions" / "2026-05-01" / "PCF2-SESSION-A"
SES_Q1 = ROOT / "sessions" / "2026-05-01" / "PCF2-SESSION-Q1"
SES_R11 = ROOT / "sessions" / "2026-05-01" / "PCF2-SESSION-R1.1"
SES_R12 = ROOT / "sessions" / "2026-05-01" / "PCF2-SESSION-R1.2"

DPS = 2000
N_GRID = list(range(50, 252, 2))
N_REF = 400
LOG = HERE / "r1_3_residualization.log"


def L(s: str):
    line = f"[{time.strftime('%H:%M:%S')}] {s}"
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def cf_value(coeffs_leading_first, N, dps):
    """Compute the truncated continued fraction L_N = b_0 + 1/(b_1 + 1/(...
    + 1/b_N)).  coeffs_leading_first are integer/rational [a_d, ..., a_0]."""
    with mp.workdps(dps):
        ms = [mp.mpf(c) for c in coeffs_leading_first]
        deg = len(ms) - 1

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


def compute_y_series(coeffs_leading_first, n_grid, n_ref, dps):
    with mp.workdps(dps):
        L_ref = cf_value(coeffs_leading_first, n_ref, dps)
        ys = []
        ns = []
        for N in n_grid:
            LN = cf_value(coeffs_leading_first, N, dps)
            d = abs(LN - L_ref)
            if d == 0:
                continue
            ys.append(float(mp.log(d)))
            ns.append(N)
    return np.array(ns, dtype=float), np.array(ys, dtype=float)


def fixed_A_fit(ns, ys, A_fixed):
    """Fit y = -A n log n + alpha n - beta log n + gamma with A fixed.
    Returns dict with alpha, beta, gamma, residuals, residual_mean, std, etc."""
    z = ys + A_fixed * ns * np.log(ns)  # subtract -A n log n  =>  add  A n log n
    X = np.column_stack([ns, -np.log(ns), np.ones_like(ns)])
    sol, *_ = np.linalg.lstsq(X, z, rcond=None)
    alpha, beta, gamma = (float(c) for c in sol)
    zhat = X @ sol
    resid = z - zhat  # equivalently y - yhat
    dof = max(len(ys) - 3, 1)
    s2 = float(np.sum(resid ** 2) / dof)
    XtX_inv = np.linalg.inv(X.T @ X)
    cov = s2 * XtX_inv
    se = np.sqrt(np.diag(cov))
    return {
        "A_fixed": float(A_fixed),
        "alpha": alpha, "beta": beta, "gamma": gamma,
        "alpha_stderr": float(se[0]),
        "beta_stderr": float(se[1]),
        "gamma_stderr": float(se[2]),
        "n_points": int(len(ys)),
        "residual_rms": float(math.sqrt(s2)),
        "residual_mean": float(resid.mean()),
        "residual_std": float(resid.std(ddof=1)) if len(resid) > 1 else 0.0,
        "residual_at_max_n": float(resid[-1]),
        "residuals": resid.tolist(),
        "ns": ns.tolist(),
    }


def free_A_fit(ns, ys):
    X = np.column_stack([-ns * np.log(ns), ns, -np.log(ns), np.ones_like(ns)])
    sol, *_ = np.linalg.lstsq(X, ys, rcond=None)
    A, alpha, beta, gamma = (float(c) for c in sol)
    yhat = X @ sol
    resid = ys - yhat
    dof = max(len(ys) - 4, 1)
    s2 = float(np.sum(resid ** 2) / dof)
    XtX_inv = np.linalg.inv(X.T @ X)
    cov = s2 * XtX_inv
    se = np.sqrt(np.diag(cov))
    return {
        "A": A, "alpha": alpha, "beta": beta, "gamma": gamma,
        "A_stderr": float(se[0]),
        "residual_rms": float(math.sqrt(s2)),
        "n_points": int(len(ys)),
    }


# ---------------------------------------------------------------- inputs

def load_cubics():
    cat = json.load(open(SES_A / "cubic_family_catalogue.json"))["families"]
    # log|j| from R1.1 csv
    j_map = {}
    with open(SES_R11 / "assembled_data_v2.csv") as fh:
        rd = csv.DictReader(fh)
        for r in rd:
            j_map[int(r["family_id"])] = {
                "log_abs_j": float(r["log_abs_j"]),
                "j_invariant": float(r["j_invariant"]),
                "A_fit_R11": float(r["A_fit"]),
                "delta_R11": float(r["delta"]),
                "Galois_group": r["Galois_group"],
                "bin": r["bin"],
            }
    rows = []
    for f in cat:
        fid = f["family_id"]
        rows.append({
            "family_id": fid,
            "coeffs": [f["alpha_3"], f["alpha_2"], f["alpha_1"], f["alpha_0"]],
            **j_map.get(fid, {}),
        })
    return rows


def load_quartics():
    cat = json.load(open(SES_Q1 / "quartic_family_catalogue.json"))["families"]
    jac = {r["family_id"]: r for r in json.load(
        open(SES_R12 / "quartic_jacobian_invariants.json"))["rows"]}
    # tail_fit_overrides (R1.2) free-A baseline for delta comparison
    tfo = json.load(open(SES_R12 / "tail_fit_overrides.json"))["overrides"]
    rows = []
    for f in cat:
        fid = f["family_id"]
        coeffs = [f["alpha_4"], f["alpha_3"], f["alpha_2"],
                  f["alpha_1"], f["alpha_0"]]
        jr = jac.get(fid, {})
        ov = tfo.get(str(fid), {})
        # R1.2 used log_abs_j = log(|j| + 1). For consistency we use the
        # SAME definition (so that ties at j=0 land at log_abs_j = 0).
        rows.append({
            "family_id": fid,
            "coeffs": coeffs,
            "log_abs_j": float(jr.get("log_abs_j", float("nan"))),
            "j_invariant": float(jr.get("j_invariant", float("nan"))),
            "is_j_zero": bool(jr.get("is_j_zero", False)),
            "is_j_1728": bool(jr.get("is_j_1728", False)),
            "Galois_group": f.get("Galois_group", ""),
            "bin": f.get("trichotomy_bin", ""),
            "A_fit_R12": float(ov.get("A", float("nan"))),
            "delta_R12": float(ov.get("delta", float("nan"))),
        })
    return rows


# ---------------------------------------------------------------- main

def fit_population(rows, d, label):
    L(f"=== Fitting {label}: {len(rows)} families, d={d}, A_fixed={2*d} ===")
    A_fixed = float(2 * d)
    out = []
    for i, r in enumerate(rows):
        t0 = time.time()
        ns, ys = compute_y_series(r["coeffs"], N_GRID, N_REF, DPS)
        fix = fixed_A_fit(ns, ys, A_fixed)
        free = free_A_fit(ns, ys)
        dt = time.time() - t0
        rec = dict(r)
        rec.update({
            "fixed": {k: v for k, v in fix.items() if k not in ("residuals", "ns")},
            "free": free,
            "residuals": fix["residuals"],
            "ns_used": fix["ns"],
            "delta_R13_free": free["A"] - 2 * d,
            "dt_seconds": dt,
        })
        out.append(rec)
        L(f"  fam {r['family_id']:>3d}: A_free={free['A']:.5f} (delta={free['A']-2*d:+.4e}), "
          f"resid_mean={fix['residual_mean']:+.4e}, std={fix['residual_std']:.4e}, "
          f"@maxN={fix['residual_at_max_n']:+.4e}, dt={dt:.1f}s")
    return out


def write_csv(rows, d, path):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "family_id", "log_abs_j", "j_invariant",
            "A_fit_R13_free", "A_stderr_R13_free", "delta_R13_free",
            "residual_mean", "residual_std", "residual_at_max_n",
            "alpha_fixedA", "beta_fixedA", "gamma_fixedA",
            "n_points",
        ])
        for r in rows:
            w.writerow([
                r["family_id"],
                r.get("log_abs_j"),
                r.get("j_invariant"),
                r["free"]["A"],
                r["free"]["A_stderr"],
                r["delta_R13_free"],
                r["fixed"]["residual_mean"],
                r["fixed"]["residual_std"],
                r["fixed"]["residual_at_max_n"],
                r["fixed"]["alpha"],
                r["fixed"]["beta"],
                r["fixed"]["gamma"],
                r["fixed"]["n_points"],
            ])


def correlations(rows, key_x, key_ys, label, K_bonf):
    out = {}
    for k in key_ys:
        x = np.array([r[key_x] for r in rows], dtype=float)
        y = np.array([_get_nested(r, k) for r in rows], dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if len(x) < 5:
            out[k] = {"label": label, "n": int(len(x)), "rho": None, "raw_p": None,
                      "bonf_p": None, "K": K_bonf}
            continue
        rho, p = stats.spearmanr(x, y)
        out[k] = {
            "label": label, "n": int(len(x)),
            "rho": float(rho), "raw_p": float(p),
            "bonf_p": float(min(1.0, p * K_bonf)),
            "K": K_bonf,
        }
    return out


def _get_nested(r, k):
    if "." in k:
        a, b = k.split(".", 1)
        return r.get(a, {}).get(b, float("nan"))
    return r.get(k, float("nan"))


def main():
    if LOG.exists():
        LOG.unlink()
    L(f"R1.3 residualization: dps={DPS}, N=[{N_GRID[0]}..{N_GRID[-1]} step 2], N_ref={N_REF}")

    cubics = load_cubics()
    quartics = load_quartics()

    cubic_rows = fit_population(cubics, d=3, label="cubic (R13-A sanity)")
    quartic_rows = fit_population(quartics, d=4, label="quartic (R13-B test)")

    write_csv(cubic_rows, 3, HERE / "residualization_d3.csv")
    write_csv(quartic_rows, 4, HERE / "residualization_d4.csv")
    L(f"  wrote residualization_d3.csv ({len(cubic_rows)} rows)")
    L(f"  wrote residualization_d4.csv ({len(quartic_rows)} rows)")

    # ------- Spearman: residual_mean, residual_at_max_n, delta_R13_free
    K = 3  # 3 residual flavours
    corr_d3 = correlations(
        cubic_rows, "log_abs_j",
        ["fixed.residual_mean", "fixed.residual_at_max_n", "delta_R13_free"],
        "d3_full50", K)
    corr_d4 = correlations(
        quartic_rows, "log_abs_j",
        ["fixed.residual_mean", "fixed.residual_at_max_n", "delta_R13_free"],
        "d4_full60", K)

    # also the d=3 reproduction of R1.1 headline (delta vs log|j|) for sanity
    # (delta_R13_free ~ 4-param free A, mirrors R1.1's 4-param fit but at
    # higher precision and on tail window)
    L("--- Spearman rho(log|j|, *) at d=3 ---")
    for k, v in corr_d3.items():
        L(f"  {k:<35s} n={v['n']}, rho={v['rho']:+.4f}, raw_p={v['raw_p']:.3g}, bonf_p={v['bonf_p']:.3g}")
    L("--- Spearman rho(log|j|, *) at d=4 ---")
    for k, v in corr_d4.items():
        L(f"  {k:<35s} n={v['n']}, rho={v['rho']:+.4f}, raw_p={v['raw_p']:.3g}, bonf_p={v['bonf_p']:.3g}")

    # save
    out = {
        "task_id": "PCF2-SESSION-R1.3",
        "phase": "AB",
        "config": {"dps": DPS, "N_grid": [N_GRID[0], N_GRID[-1], 2],
                   "N_ref": N_REF},
        "n_cubic": len(cubic_rows),
        "n_quartic": len(quartic_rows),
        "corr_d3": corr_d3,
        "corr_d4": corr_d4,
        "cubic": [{"family_id": r["family_id"],
                   "log_abs_j": r["log_abs_j"],
                   "A_R13_free": r["free"]["A"],
                   "A_R13_stderr": r["free"]["A_stderr"],
                   "delta_R13_free": r["delta_R13_free"],
                   "residual_mean": r["fixed"]["residual_mean"],
                   "residual_std": r["fixed"]["residual_std"],
                   "residual_at_max_n": r["fixed"]["residual_at_max_n"]}
                  for r in cubic_rows],
        "quartic": [{"family_id": r["family_id"],
                     "log_abs_j": r["log_abs_j"],
                     "is_j_zero": r["is_j_zero"],
                     "A_R13_free": r["free"]["A"],
                     "A_R13_stderr": r["free"]["A_stderr"],
                     "delta_R13_free": r["delta_R13_free"],
                     "residual_mean": r["fixed"]["residual_mean"],
                     "residual_std": r["fixed"]["residual_std"],
                     "residual_at_max_n": r["fixed"]["residual_at_max_n"]}
                    for r in quartic_rows],
    }
    with open(HERE / "results_v3_phase_AB.json", "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    L("wrote results_v3_phase_AB.json")

    # claims
    claims = []
    def claim(text, **k):
        rec = {"claim": text, "evidence_type": "computation",
               "dps": DPS, "reproducible": True,
               "script": "r1_3_residualization.py",
               "output_hash": k.pop("output_hash", "TBD"),
               **k}
        claims.append(rec)

    # output hash = sha256 of results json
    h = hashlib.sha256((HERE / "results_v3_phase_AB.json").read_bytes()).hexdigest()

    rho_d3_rmean = corr_d3["fixed.residual_mean"]["rho"]
    rho_d3_atmaxn = corr_d3["fixed.residual_at_max_n"]["rho"]
    rho_d3_delta = corr_d3["delta_R13_free"]["rho"]
    rho_d4_rmean = corr_d4["fixed.residual_mean"]["rho"]
    rho_d4_atmaxn = corr_d4["fixed.residual_at_max_n"]["rho"]
    rho_d4_delta = corr_d4["delta_R13_free"]["rho"]

    claim(f"R13-A: cubic residualization (FIXED A=6, dps=2000, N=[50,250], "
          f"N_ref=400) yields Spearman rho(log|j|, residual_mean) = {rho_d3_rmean:.4f} "
          f"on n={corr_d3['fixed.residual_mean']['n']} cubic families "
          f"(Bonferroni p = {corr_d3['fixed.residual_mean']['bonf_p']:.3g}, K=3)",
          output_hash=h)
    claim(f"R13-A: cubic Spearman rho(log|j|, residual_at_max_n=N=250) = "
          f"{rho_d3_atmaxn:.4f} (Bonferroni p = "
          f"{corr_d3['fixed.residual_at_max_n']['bonf_p']:.3g})", output_hash=h)
    claim(f"R13-A: cubic Spearman rho(log|j|, delta_R13_free) = "
          f"{rho_d3_delta:.4f} (Bonferroni p = "
          f"{corr_d3['delta_R13_free']['bonf_p']:.3g}); compare R1.1 headline "
          f"rho = -0.6906", output_hash=h)
    claim(f"R13-B: quartic residualization (FIXED A=8, same window) yields "
          f"Spearman rho(log|j|, residual_mean) = {rho_d4_rmean:.4f} on "
          f"n={corr_d4['fixed.residual_mean']['n']} quartic families "
          f"(Bonferroni p = {corr_d4['fixed.residual_mean']['bonf_p']:.3g}, K=3)",
          output_hash=h)
    claim(f"R13-B: quartic Spearman rho(log|j|, residual_at_max_n=N=250) = "
          f"{rho_d4_atmaxn:.4f} (Bonferroni p = "
          f"{corr_d4['fixed.residual_at_max_n']['bonf_p']:.3g})", output_hash=h)
    claim(f"R13-B: quartic Spearman rho(log|j|, delta_R13_free) = "
          f"{rho_d4_delta:.4f} (Bonferroni p = "
          f"{corr_d4['delta_R13_free']['bonf_p']:.3g}); compare R1.2 headline "
          f"rho = +0.073", output_hash=h)

    with open(HERE / "claims_phase_AB.jsonl", "w") as fh:
        for c in claims:
            fh.write(json.dumps(c) + "\n")
    L(f"wrote {len(claims)} claims to claims_phase_AB.jsonl")


if __name__ == "__main__":
    main()
