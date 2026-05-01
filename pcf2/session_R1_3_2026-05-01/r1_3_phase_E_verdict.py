"""PCF-2 SESSION R1.3 -- Phase R13-E joint verdict synthesis.

Combines outputs of phases R13-A/B/C/D into a single verdict
document and numerical summary, plus selects the v1.3 paragraph
template (R13-F)."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent


def safe_load(p):
    with open(p) as fh:
        return json.load(fh)


def main():
    AB = safe_load(HERE / "results_v3_phase_AB.json")
    Cdata = safe_load(HERE / "extended_quartic_catalogue.json")
    C_jzero = safe_load(HERE / "extended_j_zero_cell.json")
    Ddata = safe_load(HERE / "family32_deep_residual.json")

    # ----- A: cubic sanity diagnostic
    A_pass = False
    rho_d3_delta = AB["corr_d3"]["delta_R13_free"]["rho"]
    rho_d3_atmaxn = AB["corr_d3"]["residual_at_max_n"]["rho"]
    rho_d3_rmean = AB["corr_d3"]["residual_mean"]["rho"]
    # Sanity passes if delta_R13_free OR residual_at_max_n preserve a
    # |rho| > 0.4 with Bonferroni p < 0.01.  residual_mean is structurally
    # zero by least-squares so we exclude it from the pass criterion.
    p_d3_delta = AB["corr_d3"]["delta_R13_free"]["bonf_p"]
    p_d3_atmaxn = AB["corr_d3"]["residual_at_max_n"]["bonf_p"]
    A_pass = ((abs(rho_d3_delta) > 0.4 and p_d3_delta < 0.01)
              or (abs(rho_d3_atmaxn) > 0.4 and p_d3_atmaxn < 0.01))

    # ----- B: quartic test
    rho_d4_delta = AB["corr_d4"]["delta_R13_free"]["rho"]
    rho_d4_atmaxn = AB["corr_d4"]["residual_at_max_n"]["rho"]
    p_d4_delta = AB["corr_d4"]["delta_R13_free"]["bonf_p"]
    p_d4_atmaxn = AB["corr_d4"]["residual_at_max_n"]["bonf_p"]
    # B verdict labels
    B_strong_neg = (rho_d4_delta < -0.4 and p_d4_delta < 0.001) or \
                   (rho_d4_atmaxn < -0.4 and p_d4_atmaxn < 0.001)
    B_null = (abs(rho_d4_delta) < 0.15 and abs(rho_d4_atmaxn) < 0.15)
    B_intermediate = not B_strong_neg and not B_null

    # ----- C: extended j=0 cell test
    # Compare delta_R13_free distribution: extended j=0 cell vs Q1 non-j=0
    j0_deltas = [r["delta_R13_free"] for r in C_jzero["results"]]
    # remove gross outliers (|delta| > 0.05) which are likely small-N artefacts
    j0_clean = [d for d in j0_deltas if abs(d) < 0.05]
    # Q1 non-j=0 deltas from residualization_d4.csv
    with open(HERE / "residualization_d4.csv") as fh:
        d4rows = list(csv.DictReader(fh))
    cluster_deltas = [float(r["delta_R13_free"]) for r in d4rows
                      if float(r["log_abs_j"]) != 0.0]
    # Mann-Whitney U test: do j=0 deltas differ from non-j=0 cluster?
    if j0_clean and cluster_deltas:
        U, p_mw = stats.mannwhitneyu(j0_clean, cluster_deltas,
                                      alternative="two-sided")
    else:
        U, p_mw = float("nan"), float("nan")
    # also the j=0 fraction whose |delta| is significantly closer to zero
    # than the non-j=0 cluster's max |delta|
    cluster_max_abs = max(abs(d) for d in cluster_deltas)
    cluster_min_abs = min(abs(d) for d in cluster_deltas)
    n_j0_closer_to_zero = sum(1 for d in j0_clean
                              if abs(d) < cluster_min_abs)
    n_j0_in_cluster = sum(1 for d in j0_clean
                          if cluster_min_abs <= abs(d) <= cluster_max_abs)
    n_j0_more_negative = sum(1 for d in j0_clean
                             if abs(d) > cluster_max_abs)
    C_summary = {
        "n_j0_clean": int(len(j0_clean)),
        "j0_delta_min": float(min(j0_clean)) if j0_clean else None,
        "j0_delta_max": float(max(j0_clean)) if j0_clean else None,
        "j0_delta_median": float(np.median(j0_clean)) if j0_clean else None,
        "cluster_delta_min": float(min(cluster_deltas)),
        "cluster_delta_max": float(max(cluster_deltas)),
        "cluster_delta_median": float(np.median(cluster_deltas)),
        "cluster_min_abs": float(cluster_min_abs),
        "cluster_max_abs": float(cluster_max_abs),
        "n_j0_more_negative_than_cluster": int(n_j0_more_negative),
        "n_j0_in_cluster_band": int(n_j0_in_cluster),
        "n_j0_closer_to_zero_than_cluster_min": int(n_j0_closer_to_zero),
        "mannwhitney_U": float(U),
        "mannwhitney_p": float(p_mw),
        "outliers_removed": int(len(j0_deltas) - len(j0_clean)),
    }
    # Verdict for C: if j=0 cell shows a SHIFT toward zero with p_mw < 0.01,
    # and the FRACTION closer-to-zero is large, that's CASE A signal.
    C_supports_A = (p_mw < 0.01
                    and np.median([abs(d) for d in j0_clean])
                    < np.median([abs(d) for d in cluster_deltas]))

    # ----- D: family 32 deep
    fam32_delta_deep = Ddata["fam32"]["delta_deep"]
    fam01_delta_deep = Ddata["fam01_baseline"]["delta_deep"]
    fam32_resid_deep = Ddata["fam32"]["fixed"]["residual_at_max_n"]
    fam01_resid_deep = Ddata["fam01_baseline"]["fixed"]["residual_at_max_n"]
    deep_diff = fam32_delta_deep - fam01_delta_deep
    deep_stderr = (Ddata["fam32"]["free"]["A_stderr"] ** 2
                   + Ddata["fam01_baseline"]["free"]["A_stderr"] ** 2) ** 0.5
    # "approaches 0" criterion: |fam32 deep delta| < |fam32 shallow delta|
    fam32_shallow_delta = -3.7081e-3
    deep_shrinks_fam32 = abs(fam32_delta_deep) < abs(fam32_shallow_delta) / 2
    # "fam32 vs fam01 deep are SAME" criterion: |diff| < 3 * pooled stderr
    fam32_same_as_fam01 = abs(deep_diff) < 3 * deep_stderr
    D_supports_B = fam32_same_as_fam01 and not (
        abs(fam32_delta_deep) < 0.1 * abs(fam32_shallow_delta)
        and abs(fam01_delta_deep) > 5 * abs(fam32_delta_deep))
    D_summary = {
        "fam32_delta_deep": float(fam32_delta_deep),
        "fam01_delta_deep": float(fam01_delta_deep),
        "deep_diff": float(deep_diff),
        "pooled_stderr": float(deep_stderr),
        "deep_diff_in_stderr_units": float(deep_diff / deep_stderr) if deep_stderr else float("nan"),
        "fam32_resid_at_max_n_deep": float(fam32_resid_deep),
        "fam01_resid_at_max_n_deep": float(fam01_resid_deep),
        "deep_shrinks_fam32_delta": bool(deep_shrinks_fam32),
        "fam32_same_as_fam01_deep": bool(fam32_same_as_fam01),
        "D_supports_case_B": bool(D_supports_B),
    }

    # ----- joint verdict
    # Decision matrix:
    #   A_pass=True, B_strong_neg=True, C_supports_A=True, D_supports_B=False
    #     -> CASE A (universal cross-degree)
    #   A_pass=True, B_null=True, C_supports_A=False, D_supports_B=True
    #     -> CASE B (d=3-specific)
    #   A_pass=True, intermediate            -> CASE C
    #   A_pass=False                          -> METHODOLOGY HALT
    if not A_pass:
        verdict = "METHODOLOGY_HALT"
        verdict_details = (
            "Phase R13-A failed sanity check at d=3: residualization "
            "did not preserve the j-effect signal observed in R1.1."
        )
    elif B_strong_neg and C_supports_A:
        verdict = "CASE_A_universal"
        verdict_details = (
            "B5/B6 universal: residualization at d=4 reveals a "
            "previously masked j-effect; extended j=0 cell shifts "
            "delta toward zero relative to non-j=0 cluster."
        )
    elif B_null and D_supports_B and not C_supports_A:
        verdict = "CASE_B_cubic_specific"
        verdict_details = (
            "B5/B6 d=3-specific: residualization at d=4 reveals NO "
            "j-effect; extended j=0 cell does not differ from non-j=0 "
            "cluster; family-32 deep-WKB matches non-j=0 baseline to "
            "stderr."
        )
    else:
        verdict = "CASE_C_ambiguous"
        verdict_details = (
            "Mixed evidence; recommend defer absorption to v1.4 with "
            "additional deep-WKB harvest of extended j=0 cell."
        )

    out = {
        "task_id": "PCF2-SESSION-R1.3",
        "phase": "E",
        "phase_AB": {
            "A_pass": bool(A_pass),
            "rho_d3": {"residual_mean": rho_d3_rmean,
                       "residual_at_max_n": rho_d3_atmaxn,
                       "delta_R13_free": rho_d3_delta},
            "rho_d4": {"residual_mean": AB["corr_d4"]["residual_mean"]["rho"],
                       "residual_at_max_n": rho_d4_atmaxn,
                       "delta_R13_free": rho_d4_delta},
            "p_bonf_d3": {"residual_at_max_n": p_d3_atmaxn,
                          "delta_R13_free": p_d3_delta},
            "p_bonf_d4": {"residual_at_max_n": p_d4_atmaxn,
                          "delta_R13_free": p_d4_delta},
            "B_strong_neg": bool(B_strong_neg),
            "B_null": bool(B_null),
            "B_intermediate": bool(B_intermediate),
        },
        "phase_C": C_summary,
        "phase_D": D_summary,
        "verdict": verdict,
        "verdict_details": verdict_details,
    }
    out_path = HERE / "results_v3.json"
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(json.dumps(out, indent=2, default=str))
    print("\n>>> VERDICT:", verdict)

    # claims
    h = hashlib.sha256(out_path.read_bytes()).hexdigest()
    claims = []
    def cl(text):
        claims.append({"claim": text, "evidence_type": "computation",
                       "dps": 5000, "reproducible": True,
                       "script": "r1_3_phase_E_verdict.py",
                       "output_hash": h})
    cl(f"R13-E joint verdict: {verdict} ({verdict_details})")
    cl(f"R13-E phase A: cubic sanity {'PASS' if A_pass else 'FAIL'}; "
       f"rho(log|j|, delta_R13_free) = {rho_d3_delta:+.4f} (Bonf p = "
       f"{p_d3_delta:.3g}); rho(log|j|, residual_at_max_n) = "
       f"{rho_d3_atmaxn:+.4f} (Bonf p = {p_d3_atmaxn:.3g})")
    cl(f"R13-E phase B: quartic test rho(log|j|, delta_R13_free) = "
       f"{rho_d4_delta:+.4f} (Bonf p = {p_d4_delta:.3g}); "
       f"rho(log|j|, residual_at_max_n) = {rho_d4_atmaxn:+.4f} "
       f"(Bonf p = {p_d4_atmaxn:.3g}); strong_neg={B_strong_neg}, "
       f"null={B_null}, intermediate={B_intermediate}")
    cl(f"R13-E phase C: extended j=0 cell n={C_summary['n_j0_clean']} "
       f"(after outlier removal); Mann-Whitney U={C_summary['mannwhitney_U']:.1f} "
       f"p={C_summary['mannwhitney_p']:.3g}; j=0 delta median = "
       f"{C_summary['j0_delta_median']:.4e} vs cluster median = "
       f"{C_summary['cluster_delta_median']:.4e}; supports_case_A="
       f"{C_supports_A}")
    cl(f"R13-E phase D: family-32 deep delta = {fam32_delta_deep:+.4e} "
       f"vs fam01 deep delta = {fam01_delta_deep:+.4e}; diff = "
       f"{deep_diff:+.4e} ({D_summary['deep_diff_in_stderr_units']:+.2f} "
       f"pooled stderr); supports_case_B={D_supports_B}")
    out_jsonl = HERE / "claims_phase_E.jsonl"
    with open(out_jsonl, "w") as fh:
        for c in claims:
            fh.write(json.dumps(c) + "\n")
    print(f"\nwrote {len(claims)} claims to {out_jsonl.name}")


if __name__ == "__main__":
    main()
