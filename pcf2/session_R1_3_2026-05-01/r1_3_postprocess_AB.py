"""Post-process: read residualization_{d3,d4}.csv to compute the
Phase R13-A and R13-B correlations and emit results_v3_phase_AB.json
and claims_phase_AB.jsonl.  Avoids re-running the 5-min mpmath fit."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent


def load_csv(path):
    with open(path) as fh:
        return [{k: (float(v) if v not in ("", "nan", "None") and k != "family_id" else (v))
                 for k, v in row.items()} for row in csv.DictReader(fh)]


def spearman_safe(x, y):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 5:
        return float("nan"), float("nan"), int(len(x))
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p), int(len(x))


def corr_block(rows, key_x, key_ys, label, K):
    out = {}
    for k in key_ys:
        xs = [r[key_x] for r in rows]
        ys = [r[k] for r in rows]
        rho, p, n = spearman_safe(xs, ys)
        bonf = min(1.0, p * K) if math.isfinite(p) else float("nan")
        out[k] = {"label": label, "n": n, "rho": rho, "raw_p": p,
                  "bonf_p": bonf, "K": K}
    return out


def main():
    cubic = load_csv(HERE / "residualization_d3.csv")
    quartic = load_csv(HERE / "residualization_d4.csv")
    print(f"loaded {len(cubic)} cubic, {len(quartic)} quartic rows")

    K = 3  # 3 residual flavours per degree
    keys = ["residual_mean", "residual_at_max_n", "delta_R13_free"]
    corr_d3 = corr_block(cubic, "log_abs_j", keys, "d3_full50", K)
    corr_d4 = corr_block(quartic, "log_abs_j", keys, "d4_full60", K)

    print("=== Spearman rho(log|j|, *) at d=3 ===")
    for k, v in corr_d3.items():
        print(f"  {k:<25s} n={v['n']}, rho={v['rho']:+.4f}, raw_p={v['raw_p']:.3g}, bonf_p={v['bonf_p']:.3g}")
    print("=== Spearman rho(log|j|, *) at d=4 ===")
    for k, v in corr_d4.items():
        print(f"  {k:<25s} n={v['n']}, rho={v['rho']:+.4f}, raw_p={v['raw_p']:.3g}, bonf_p={v['bonf_p']:.3g}")

    # additional diagnostic: residual_std distribution
    rstd_d3 = [r["residual_std"] for r in cubic]
    rstd_d4 = [r["residual_std"] for r in quartic]
    rmean_d3 = [r["residual_mean"] for r in cubic]
    rmean_d4 = [r["residual_mean"] for r in quartic]
    print(f"d3 residual_std: min={min(rstd_d3):.3e}, max={max(rstd_d3):.3e}, "
          f"median={np.median(rstd_d3):.3e}")
    print(f"d4 residual_std: min={min(rstd_d4):.3e}, max={max(rstd_d4):.3e}, "
          f"median={np.median(rstd_d4):.3e}")
    print(f"d3 residual_mean range: [{min(rmean_d3):.3e}, {max(rmean_d3):.3e}]")
    print(f"d4 residual_mean range: [{min(rmean_d4):.3e}, {max(rmean_d4):.3e}]")

    # is residual_std/|residual_mean| > 0.5 most of the time?
    bad_d3 = sum(1 for r in cubic if abs(r["residual_mean"]) > 0
                 and r["residual_std"] / abs(r["residual_mean"]) > 0.5)
    bad_d4 = sum(1 for r in quartic if abs(r["residual_mean"]) > 0
                 and r["residual_std"] / abs(r["residual_mean"]) > 0.5)
    print(f"families with std/|mean| > 0.5: d3={bad_d3}/{len(cubic)}, "
          f"d4={bad_d4}/{len(quartic)}")

    # j=0 cell at d=4: family 32 only
    j0 = [r for r in quartic if r["log_abs_j"] == 0.0]
    print(f"d=4 j=0 cell size: {len(j0)} (family ids: "
          f"{[r['family_id'] for r in j0]})")
    nonj0 = [r for r in quartic if r["log_abs_j"] != 0.0]
    if j0 and nonj0:
        print(f"  j=0 residual_mean: {[r['residual_mean'] for r in j0]}")
        print(f"  non-j=0 residual_mean: min={min(r['residual_mean'] for r in nonj0):.3e}, "
              f"max={max(r['residual_mean'] for r in nonj0):.3e}, "
              f"median={np.median([r['residual_mean'] for r in nonj0]):.3e}")
        print(f"  j=0 delta_R13_free: {[r['delta_R13_free'] for r in j0]}")
        print(f"  non-j=0 delta_R13_free: min={min(r['delta_R13_free'] for r in nonj0):.3e}, "
              f"max={max(r['delta_R13_free'] for r in nonj0):.3e}")

    # save
    out = {
        "task_id": "PCF2-SESSION-R1.3",
        "phase": "AB",
        "config": {"dps": 2000, "N_grid": [50, 250, 2], "N_ref": 400,
                   "fixed_A_d3": 6.0, "fixed_A_d4": 8.0},
        "n_cubic": len(cubic), "n_quartic": len(quartic),
        "corr_d3": corr_d3, "corr_d4": corr_d4,
        "residual_std_d3_summary": {"min": float(min(rstd_d3)),
                                     "max": float(max(rstd_d3)),
                                     "median": float(np.median(rstd_d3))},
        "residual_std_d4_summary": {"min": float(min(rstd_d4)),
                                     "max": float(max(rstd_d4)),
                                     "median": float(np.median(rstd_d4))},
        "residual_mean_d3_range": [float(min(rmean_d3)), float(max(rmean_d3))],
        "residual_mean_d4_range": [float(min(rmean_d4)), float(max(rmean_d4))],
        "noisy_families_d3": int(bad_d3),
        "noisy_families_d4": int(bad_d4),
        "d4_j_zero_cell": [r["family_id"] for r in j0],
    }
    out_path = HERE / "results_v3_phase_AB.json"
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\nwrote {out_path.name}")

    # claims
    h = hashlib.sha256(out_path.read_bytes()).hexdigest()
    claims = []
    def cl(text):
        claims.append({"claim": text, "evidence_type": "computation",
                       "dps": 2000, "reproducible": True,
                       "script": "r1_3_residualization.py + r1_3_postprocess_AB.py",
                       "output_hash": h})

    cl(f"R13-A: cubic residualization (FIXED A=6, dps=2000, N=[50,250], N_ref=400) "
       f"yields Spearman rho(log|j|, residual_mean) = "
       f"{corr_d3['residual_mean']['rho']:+.4f} on n={corr_d3['residual_mean']['n']} "
       f"cubic families (Bonferroni p = {corr_d3['residual_mean']['bonf_p']:.3g}, K=3)")
    cl(f"R13-A: cubic Spearman rho(log|j|, residual_at_max_n=N=250) = "
       f"{corr_d3['residual_at_max_n']['rho']:+.4f} (Bonferroni p = "
       f"{corr_d3['residual_at_max_n']['bonf_p']:.3g})")
    cl(f"R13-A: cubic Spearman rho(log|j|, delta_R13_free) = "
       f"{corr_d3['delta_R13_free']['rho']:+.4f} (Bonferroni p = "
       f"{corr_d3['delta_R13_free']['bonf_p']:.3g}); compare R1.1 headline "
       f"rho = -0.6906 (n=50, R1.1 reported 4-param fit at smaller window).")
    cl(f"R13-B: quartic residualization (FIXED A=8, same window) yields "
       f"Spearman rho(log|j|, residual_mean) = "
       f"{corr_d4['residual_mean']['rho']:+.4f} on n={corr_d4['residual_mean']['n']} "
       f"quartic families (Bonferroni p = {corr_d4['residual_mean']['bonf_p']:.3g}, K=3)")
    cl(f"R13-B: quartic Spearman rho(log|j|, residual_at_max_n=N=250) = "
       f"{corr_d4['residual_at_max_n']['rho']:+.4f} (Bonferroni p = "
       f"{corr_d4['residual_at_max_n']['bonf_p']:.3g})")
    cl(f"R13-B: quartic Spearman rho(log|j|, delta_R13_free) = "
       f"{corr_d4['delta_R13_free']['rho']:+.4f} (Bonferroni p = "
       f"{corr_d4['delta_R13_free']['bonf_p']:.3g}); compare R1.2 headline "
       f"rho = +0.073")
    cl(f"R13-A: residual_std on cubic catalogue: median = "
       f"{np.median(rstd_d3):.3e}, max = {max(rstd_d3):.3e}; ratio std/|mean| > 0.5 "
       f"in {bad_d3}/{len(cubic)} families (noise gate)")
    cl(f"R13-B: residual_std on quartic catalogue: median = "
       f"{np.median(rstd_d4):.3e}, max = {max(rstd_d4):.3e}; ratio std/|mean| > 0.5 "
       f"in {bad_d4}/{len(quartic)} families (noise gate)")

    out_jsonl = HERE / "claims_phase_AB.jsonl"
    with open(out_jsonl, "w") as fh:
        for c in claims:
            fh.write(json.dumps(c) + "\n")
    print(f"wrote {len(claims)} claims to {out_jsonl.name}")


if __name__ == "__main__":
    main()
