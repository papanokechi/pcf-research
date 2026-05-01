"""ql15_ql26_corrected_probe.py -- Phase 6 of CC-PIPELINE-G.

Apply the Newton-polygon + Birkhoff formal-solution pipeline to all
five non-V_quad Delta<0 degree-2 PCFs (QL01, QL02, QL06, QL15, QL26)
and run the Variant-A flip test in the corrected formal-series space.

Per-family:
  - Newton polygon slope and char equation -> action c, |xi_0_family|
  - Indicial exponent rho = -3/2 - beta1/beta2
  - Formal series coefficients a_1, ..., a_K (in u = sqrt(z))
  - Borel radius |w*| (numerical, Gevrey-1-in-u Domb-Sykes tail)
  - Painleve fingerprint vs V_quad reference (Phase 5 verdict at >= 15
    digits): match (=> P-III(D_6)) or no-match (=> trans-numerical
    or different reduction).

The HALT clause "QL15 or QL26 detects Painleve at >= 15 digits in the
corrected space" is checked here.
"""

from __future__ import annotations
import json
from pathlib import Path
import sys
import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from newton_birkhoff import (
    build_op, newton, formal_solve, borel_radius, painleve_fingerprint
)

DPS = 100
K_FORMAL = 200
THRESHOLD = 15  # >= digits to trigger Variant-A flip

# Family roster: V_quad reference + 5 others
FAMILIES = [
    {"name": "V_quad", "alpha1": 0,  "alpha0": 1, "beta2": 3, "beta1": 1,  "beta0": 1, "Delta": -11, "lit_xi0": "2/sqrt(3)"},
    {"name": "QL01",   "alpha1": 1,  "alpha0": 0, "beta2": 1, "beta1": 1,  "beta0": 1, "Delta": -3,  "lit_xi0": "2"},
    {"name": "QL02",   "alpha1": 2,  "alpha0": 1, "beta2": 1, "beta1": 1,  "beta0": 1, "Delta": -4,  "lit_xi0": "2"},
    {"name": "QL06",   "alpha1": 1,  "alpha0": 0, "beta2": 2, "beta1": -1, "beta0": 1, "Delta": -7,  "lit_xi0": "sqrt(2)"},
    {"name": "QL15",   "alpha1": 1,  "alpha0": 0, "beta2": 3, "beta1": -2, "beta0": 2, "Delta": -20, "lit_xi0": "2/sqrt(3)"},
    {"name": "QL26",   "alpha1": -3, "alpha0": 1, "beta2": 4, "beta1": -2, "beta0": 2, "Delta": -28, "lit_xi0": "1"},
]


def process(spec):
    op = build_op(family=spec["name"],
                  alpha1=spec["alpha1"], alpha0=spec["alpha0"],
                  beta2=spec["beta2"], beta1=spec["beta1"], beta0=spec["beta0"])
    np_data = newton(op)
    fp = formal_solve(op, K=K_FORMAL, dps=DPS, sign=+1)
    bp = borel_radius(fp, n_use=K_FORMAL // 2)

    xi0 = np_data["slopes"][0]["abs_action_xi0"]
    rho = fp["rho"]
    return {
        "family": spec["name"],
        "Delta": spec["Delta"],
        "recurrence": {
            "a_n": f"{spec['alpha1']} n + {spec['alpha0']}",
            "b_n": f"{spec['beta2']} n^2 + {spec['beta1']} n + {spec['beta0']}",
        },
        "lit_xi0_str": spec["lit_xi0"],
        "xi0_extracted": xi0,
        "xi0_str": mp.nstr(xi0, 30),
        "rho": rho,
        "rho_str": mp.nstr(rho, 30),
        "a_first10": [fp["a"][k] for k in range(1, 11)],
        "borel_radius_numeric": bp["w_star_abs_estimate"],
        "newton_polygon": np_data["slopes"][0],
        "formal": fp,
    }


def main():
    mp.mp.dps = DPS
    results = []
    for spec in FAMILIES:
        print(f"\n=== {spec['name']} (Delta = {spec['Delta']}) ===")
        r = process(spec)
        print(f"  beta2 = {spec['beta2']}, xi_0 = 2/sqrt(beta2) = {r['xi0_str'][:40]}")
        print(f"  rho = -3/2 - beta1/beta2 = {r['rho_str'][:30]}")
        print(f"  a_1..a_5 = {[mp.nstr(a, 10) for a in r['a_first10'][:5]]}")
        print(f"  Borel radius |w*| (numerical tail) = {mp.nstr(r['borel_radius_numeric'], 8)}")
        results.append(r)

    # V_quad reference for fingerprint comparison
    vquad = next(r for r in results if r["family"] == "V_quad")

    print("\n" + "=" * 70)
    print("VARIANT-A FLIP TEST (Phase 6): per-family Painleve fingerprint vs V_quad")
    print("=" * 70)
    flips = []
    table_rows = []
    for r in results:
        fp_compare = painleve_fingerprint(r["formal"], vquad["formal"], threshold=THRESHOLD)
        digits_xi0 = fp_compare["digits"]["c"]
        digits_rho = fp_compare["digits"]["rho"]
        a_digits = [fp_compare["digits"][f"a_{k}"] for k in range(1, 5)]
        min_d = fp_compare["min_digits"]
        verdict = fp_compare["verdict"]
        if r["family"] != "V_quad" and verdict == "match":
            flips.append(r["family"])
            painleve_class = "P-III(D_6) [FLIP candidate]"
        elif r["family"] == "V_quad":
            painleve_class = "P-III(D_6) [reference anchor]"
        else:
            painleve_class = "trans-numerical (no Painleve match in corrected space)"
        print(f"  {r['family']:<8} digits(xi0)={digits_xi0:6.2f} digits(rho)={digits_rho:6.2f} "
              f"a_1..a_4 digits={[f'{d:5.2f}' for d in a_digits]} -> {verdict}")
        table_rows.append({
            "family": r["family"],
            "Delta": r["Delta"],
            "newton_slope": "1/2",
            "xi0_str": r["xi0_str"],
            "xi0_lit_str": next(s["lit_xi0"] for s in FAMILIES if s["name"] == r["family"]),
            "rho_str": r["rho_str"],
            "borel_radius_num": mp.nstr(r["borel_radius_numeric"], 6),
            "fingerprint_min_digits": float(min_d),
            "verdict": verdict,
            "painleve_class": painleve_class,
        })

    print(f"\nFLIPS detected (>={THRESHOLD} digit P-III match in corrected space): {flips}")
    halt_flip = len(flips) > 0

    # Per-family JSON dump
    for r in results:
        family_safe = r["family"].lower().replace("_", "")
        out = {
            "family": r["family"],
            "Delta": r["Delta"],
            "recurrence": r["recurrence"],
            "newton_polygon": {
                "slope": "1/2",
                "char_polynomial": str(r["newton_polygon"]["char_poly"]),
                "abs_action_xi0": str(r["xi0_extracted"]),
                "lit_xi0": r["lit_xi0_str"],
                "gevrey_class_in_z": 2,
                "gevrey_class_in_u": 1,
            },
            "rho": str(r["rho"]),
            "formal_a_first10": [str(a) for a in r["a_first10"]],
            "borel_radius_numeric": str(r["borel_radius_numeric"]),
            "K_formal": K_FORMAL,
            "dps": DPS,
        }
        if r["family"] in ("QL15", "QL26"):
            Path(HERE / f"results_{family_safe}.json").write_text(
                json.dumps(out, indent=2, default=str), encoding="utf-8"
            )
            print(f"Wrote results_{family_safe}.json")

    # Summary
    summary = {
        "families": [
            {**row, "fingerprint_min_digits": row["fingerprint_min_digits"]}
            for row in table_rows
        ],
        "flips_at_threshold": flips,
        "threshold_digits": THRESHOLD,
        "variant_A_verdict": "FLIPPED" if halt_flip else "BOTH ARTEFACT (stands)",
    }
    Path(HERE / "phase6_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )

    # cc_channel_table_v2.tex
    write_cc_table(table_rows)

    return summary, halt_flip


def write_cc_table(rows):
    lines = []
    lines.append("\\begin{table}[ht]")
    lines.append("\\centering")
    lines.append("\\small")
    lines.append("\\caption{CC-channel via Newton-polygon / Birkhoff: 6 families, "
                 "$\\Delta<0$ degree-2 PCFs.  All families share Newton-polygon slope "
                 "$1/2$ at $z=0$ (Gevrey-2 in $z$, Gevrey-1 in $u=\\sqrt{z}$). "
                 "$\\xi_0 = |c|$ where $c$ is the slope-$1/2$ characteristic root "
                 "$c^2 = 4/\\beta_2$. Painleve verdict at $\\geq 15$-digit fingerprint "
                 "agreement vs the V$_{\\text{quad}}$ P-III($D_6$) anchor "
                 "(formal-solution $(\\xi_0,\\rho,a_1,a_2,a_3,a_4)$).}")
    lines.append("\\label{tab:cc-channel-v2}")
    lines.append("\\begin{tabular}{l c c c c c}")
    lines.append("\\toprule")
    lines.append("Family & $\\Delta$ & $\\xi_0$ (extr.) vs lit. & $\\rho$ & "
                 "Fingerprint min digits & Painlev\\'e verdict \\\\")
    lines.append("\\midrule")
    for r in rows:
        line = (f"{r['family']} & {r['Delta']} & "
                f"$={r['xi0_lit_str']}$ ({r['fingerprint_min_digits']:.1f}d) & "
                f"${r['rho_str'][:8]}$ & "
                f"{r['fingerprint_min_digits']:.2f} & "
                f"{r['painleve_class']} \\\\")
        lines.append(line)
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    lines.append("\\end{table}")
    Path(HERE / "cc_channel_table_v2.tex").write_text("\n".join(lines), encoding="utf-8")
    print("Wrote cc_channel_table_v2.tex")


if __name__ == "__main__":
    summary, halt_flip = main()
    print(f"\n=== FINAL VERDICT ===")
    print(f"Variant-A flip test (corrected space): {summary['variant_A_verdict']}")
    print(f"Flips: {summary['flips_at_threshold']}")
