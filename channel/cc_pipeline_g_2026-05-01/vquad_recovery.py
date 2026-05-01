"""V_quad CC-channel recovery via Newton polygon + Birkhoff.

Phases 1-5: build the linear ODE for f(z) = sum Q_n z^n with V_quad
recurrence (a_n=1, b_n=3 n^2 + n + 1), compute the Newton polygon at
z=0 (slope 1/2, char eq c^2 = 4/beta2 = 4/3, c = +/- 2/sqrt(3)),
construct the Frobenius/Birkhoff formal solutions in u = sqrt(z),
extract the Gevrey-1 Borel-plane radius, and verify the literature
constant xi_0 = 2/sqrt(3) to >= 30 digits.
"""

from __future__ import annotations
import json
import hashlib
from pathlib import Path
import mpmath as mp

import sys
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from newton_birkhoff import (
    build_op, newton, formal_solve, borel_radius
)

DPS = 100
K_FORMAL = 200


def file_hash(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()[:16]


def main():
    mp.mp.dps = DPS

    # V_quad: a_n = 0 n + 1,  b_n = 3 n^2 + n + 1
    op = build_op(family="V_quad",
                  alpha1=0, alpha0=1,
                  beta2=3, beta1=1, beta0=1)

    # Phase 2: Newton polygon
    np_data = newton(op)
    xi0_lit = mp.mpf(2) / mp.sqrt(mp.mpf(3))
    xi0_extracted = np_data["slopes"][0]["abs_action_xi0"]
    xi0_err = abs(xi0_extracted - xi0_lit)
    if xi0_err == 0:
        digits_xi0 = mp.mpf(2 * DPS)
    else:
        digits_xi0 = -mp.log10(xi0_err / xi0_lit)
    print(f"[Phase 2] Newton polygon slope=1/2, char eq 1 - 3 c^2/4")
    print(f"          c (extracted)   = {xi0_extracted}")
    print(f"          xi_0 (literature) = 2/sqrt(3) = {xi0_lit}")
    print(f"          Digits agreement: {digits_xi0}")

    # Phase 3-4: formal solutions
    fp = formal_solve(op, K=K_FORMAL, dps=DPS, sign=+1)
    fm = formal_solve(op, K=K_FORMAL, dps=DPS, sign=-1)
    print(f"[Phase 3] Formal solution f_+ : c = {fp['c']}, rho = {fp['rho']}")
    print(f"          a_1..a_5 sample:")
    for k in range(1, 6):
        print(f"            a_{k} = {mp.nstr(fp['a'][k], 25)}")

    # Phase 4: Borel radius (in u, Gevrey-1)
    bp = borel_radius(fp, n_use=K_FORMAL // 2)
    bm = borel_radius(fm, n_use=K_FORMAL // 2)
    print(f"[Phase 4] Borel radius |w*| (sign +) = {mp.nstr(bp['w_star_abs_estimate'], 20)}")
    print(f"          Borel radius |w*| (sign -) = {mp.nstr(bm['w_star_abs_estimate'], 20)}")

    # The Borel-plane singularity location is the OTHER characteristic root
    # (resurgence: alien derivative connects f_+ with f_-).  Distance
    # |w*| = |c_+ - c_-| / 2 = |c|.  So expected w* = 2 c0 = 4/sqrt(3)?
    # Actually the standard resurgence result for a slope-1 Gevrey-1
    # series at u=0 with action c/u gives a Borel singularity at
    # w = c (the action of the partner solution), so |w*| = |c| = 2/sqrt(3).
    expected_borel_radius = xi0_lit
    err_b = abs(bp["w_star_abs_estimate"] - expected_borel_radius) / expected_borel_radius
    digits_b = -mp.log10(err_b) if err_b > 0 else mp.mpf(2 * DPS)
    print(f"          |w*| expected = 2/sqrt(3) = {expected_borel_radius}")
    print(f"          Digits agreement (numerical Borel radius): {digits_b}")

    # Phase 5: P-III(D_6) tau parameter
    # The Riemann-Hilbert datum of P-III(D_6) at the V_quad reduction is
    # encoded by (rho, xi_0).  We report rho.
    rho = fp["rho"]
    print(f"[Phase 5] P-III(D_6) tau-related invariant rho = {rho}")

    out = {
        "family": "V_quad",
        "recurrence": {"a_n": "1", "b_n": "3 n^2 + n + 1"},
        "newton_polygon": {
            "slope": "1/2",
            "edge": "(0,0) -- (1,2)",
            "multiplicity": 2,
            "char_polynomial": "1 - (3/4) c^2",
            "gevrey_class_in_z": 2,
            "gevrey_class_in_u": 1,
        },
        "xi_0_literature": str(xi0_lit),
        "xi_0_extracted_analytic": str(xi0_extracted),
        "xi_0_digits_agreement_analytic": float(digits_xi0),
        "rho": str(rho),
        "formal_a_coefficients_sample": {
            f"a_{k}": str(fp["a"][k]) for k in range(1, 11)
        },
        "borel_radius_numeric_plus": str(bp["w_star_abs_estimate"]),
        "borel_radius_numeric_minus": str(bm["w_star_abs_estimate"]),
        "borel_radius_digits_agreement": float(digits_b),
        "borel_sample_radii": [(k, r) for (k, r) in bp["radii_sample"]],
        "K_formal": K_FORMAL,
        "dps": DPS,
        "verdict": "RECOVERY_PASS" if float(digits_xi0) >= 30 else "RECOVERY_FAIL",
    }
    Path(HERE / "results_vquad.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8"
    )
    print(f"\nVerdict: {out['verdict']}  (digits xi_0 = {digits_xi0})")
    print("Wrote results_vquad.json")


if __name__ == "__main__":
    main()
