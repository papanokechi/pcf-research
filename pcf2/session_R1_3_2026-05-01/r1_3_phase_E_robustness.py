"""R13-E robustness: control phase-C extended j=0 cell shift for
coefficient magnitude.  If the apparent "shift toward zero" of j=0
deltas is driven by larger coefficient norms (which mechanically
inflate sub-leading WKB shape), then the C-signal is a confound."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent

C = json.load(open(HERE / "extended_j_zero_cell.json"))
results = C["results"]
# clean: |delta| < 0.05
clean = [r for r in results if abs(r["delta_R13_free"]) < 0.05]
print(f"Phase R13-E robustness: {len(clean)} clean j=0 quartics")

deltas = np.array([r["delta_R13_free"] for r in clean])
coeffs = [r["coeffs"] for r in clean]
# coefficient magnitude
abs_sum = np.array([abs(c[1]) + abs(c[2]) + abs(c[3]) + abs(c[4]) for c in coeffs])
abs_max = np.array([max(abs(c[1]), abs(c[2]), abs(c[3]), abs(c[4])) for c in coeffs])
log_disc = np.array([math.log(abs(r["Delta_4_exact"]) + 1) for r in clean])
a4_arr = np.array([c[0] for c in coeffs])

print("\n=== Spearman: delta vs coefficient invariants (j=0 cell only) ===")
for name, arr in [("|a3|+|a2|+|a1|+|a0|", abs_sum),
                  ("max(|a3|,...,|a0|)", abs_max),
                  ("log(|Delta_4|+1)", log_disc),
                  ("a4", a4_arr)]:
    rho, p = stats.spearmanr(arr, deltas)
    print(f"  {name:<25s} rho={rho:+.4f}, raw_p={p:.3g}")

# also Pearson on |delta| vs coefficient size
print("\n=== Pearson: |delta| vs coefficient invariants ===")
abs_delta = np.abs(deltas)
for name, arr in [("|a3|+|a2|+|a1|+|a0|", abs_sum),
                  ("max(|a3|,...,|a0|)", abs_max),
                  ("a4", a4_arr)]:
    r, p = stats.pearsonr(arr, abs_delta)
    print(f"  {name:<25s} pearson={r:+.4f}, p={p:.3g}")

# -------- baseline non-j=0 cluster from Q1 60 (residualization_d4.csv)
import csv
with open(HERE / "residualization_d4.csv") as fh:
    d4rows = list(csv.DictReader(fh))
d4nonj0 = [r for r in d4rows if float(r["log_abs_j"]) != 0.0]
# Q1 abs_sum
q1_coeffs = []
import json
qjac = json.load(open(HERE.parent / "PCF2-SESSION-R1.2" / "quartic_jacobian_invariants.json"))["rows"]
qjac_map = {r["family_id"]: r["coeffs"] for r in qjac}
q1_clean = [(int(r["family_id"]), float(r["delta_R13_free"]),
             qjac_map[int(r["family_id"])])
            for r in d4nonj0]
q1_abs_sum = np.array([abs(c[1]) + abs(c[2]) + abs(c[3]) + abs(c[4])
                        for _, _, c in q1_clean])
q1_deltas = np.array([d for _, d, _ in q1_clean])

print(f"\n=== Coefficient-magnitude comparison ===")
print(f"j=0 cell (n={len(deltas)})    abs_sum range: [{abs_sum.min()}, {abs_sum.max()}], median={np.median(abs_sum):.1f}")
print(f"Q1 non-j=0 cluster (n={len(q1_deltas)}) abs_sum range: "
      f"[{q1_abs_sum.min()}, {q1_abs_sum.max()}], median={np.median(q1_abs_sum):.1f}")

# Stratified Mann-Whitney: compare j=0 vs non-j=0 within MATCHED
# coefficient bands.  Take j=0 quartics with abs_sum in the Q1 range,
# and compare those.
band = (q1_abs_sum.min(), q1_abs_sum.max())
j0_in_band = deltas[(abs_sum >= band[0]) & (abs_sum <= band[1])]
print(f"\n=== Stratified test: j=0 quartics whose abs_sum is within Q1 cluster band ===")
print(f"j=0 in band: n={len(j0_in_band)}, range [{j0_in_band.min():.4e}, {j0_in_band.max():.4e}], median={np.median(j0_in_band):.4e}")
print(f"Q1 non-j=0:  n={len(q1_deltas)}, range [{q1_deltas.min():.4e}, {q1_deltas.max():.4e}], median={np.median(q1_deltas):.4e}")
if len(j0_in_band) >= 5:
    U, p_mw = stats.mannwhitneyu(j0_in_band, q1_deltas, alternative="two-sided")
    print(f"Mann-Whitney U={U:.1f}, p={p_mw:.3g}")
else:
    p_mw = float("nan")
    U = float("nan")

# write supplemental json
out = {
    "task_id": "PCF2-SESSION-R1.3",
    "phase": "E_robustness",
    "n_j0_clean": int(len(deltas)),
    "j0_delta_summary": {
        "min": float(deltas.min()), "max": float(deltas.max()),
        "median": float(np.median(deltas))
    },
    "j0_abs_sum_summary": {
        "min": int(abs_sum.min()), "max": int(abs_sum.max()),
        "median": float(np.median(abs_sum))
    },
    "q1_nonj0_delta_summary": {
        "min": float(q1_deltas.min()), "max": float(q1_deltas.max()),
        "median": float(np.median(q1_deltas))
    },
    "q1_abs_sum_summary": {
        "min": int(q1_abs_sum.min()), "max": int(q1_abs_sum.max()),
        "median": float(np.median(q1_abs_sum))
    },
    "stratified": {
        "band": [float(band[0]), float(band[1])],
        "n_j0_in_band": int(len(j0_in_band)),
        "j0_in_band_median": float(np.median(j0_in_band)) if len(j0_in_band) else None,
        "j0_in_band_min": float(j0_in_band.min()) if len(j0_in_band) else None,
        "j0_in_band_max": float(j0_in_band.max()) if len(j0_in_band) else None,
        "q1_nonj0_median": float(np.median(q1_deltas)),
        "mannwhitney_U": float(U),
        "mannwhitney_p": float(p_mw),
    },
}
out_path = HERE / "phase_E_robustness.json"
with open(out_path, "w") as fh:
    json.dump(out, fh, indent=2)
print(f"\nwrote {out_path.name}")
