"""
PCF-2 Session A -- CUBIC-FAMILY-ENUMERATION

Builds the first cubic-family catalogue for the PCF-2 program.
For each candidate b(n) = a3 n^3 + a2 n^2 + a1 n + a0 (Z-primitive,
irreducible over Q), computes:
    - Delta_3 = disc(b, n)  (exact integer)
    - prime factorization of Delta_3
    - Galois group of splitting field over Q  (C_3 vs S_3)
    - splitting field / CM field invariant Q(sqrt(Delta_3))
    - Conjecture B3(i) trichotomy bin
    - numerical limit estimate L_N for N in {500, 1000, 2000} at dps=200,
      with a(n) = 1, using the standard CF
        L = b(0) + a(1)/( b(1) + a(2)/( b(2) + ... ))
      (the recurrence-free continued-fraction tail).

Family enumeration order:
    a3 in {1, 2, 3, 5, 7}
    a2, a1, a0 in {-3, -2, -1, 0, 1, 2, 3}  (lex order: a3, a2, a1, a0)
    keep only:  gcd(a3, a2, a1, a0) == 1  AND  b irreducible over Q
    take first 50 such tuples.

Reducible-control list: first 10 (Z-primitive but reducible) tuples
in the same lex order.

Coefficient ordering convention: this script enumerates by ASCENDING
degree internally for sympy poly(b, n), i.e. b = a3*n^3 + a2*n^2 +
a1*n + a0. The catalogue records alpha_3, alpha_2, alpha_1, alpha_0
explicitly so there is no ambiguity in downstream consumers.

Apery anchor: b = 34 n^3 + 51 n^2 + 27 n + 5 is OUT of the enumeration
ranges (coefficients > 3) but is computed separately as a calibration.
Same for the Catalan seed b = n^3 - 2 (in range, should land in -_S3_CM)
and the cyclic-C_3 probe b = n^3 - 3 n + 1 (in range, should land in
+_C3_real).

Outputs written next to this script:
    cubic_family_catalogue.json
    galois_distribution_summary.json
    calibration_anchors.json
    claims.jsonl
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
from typing import Any

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------- utils

def fmt_field(disc: int) -> str:
    """LaTeX for Q(sqrt(D)), reduced to squarefree part."""
    if disc == 0:
        return r"\mathbb{Q}"
    sign = -1 if disc < 0 else 1
    d = abs(disc)
    sq = 1
    p = 2
    while p * p <= d:
        while d % (p * p) == 0:
            d //= p * p
            sq *= p
        p += 1
    sf = sign * d
    if sf == 1:
        return r"\mathbb{Q}"
    return rf"\mathbb{{Q}}(\sqrt{{{sf}}})"


def squarefree_part(n: int) -> int:
    if n == 0:
        return 0
    sign = -1 if n < 0 else 1
    m = abs(n)
    out = 1
    p = 2
    while p * p <= m:
        e = 0
        while m % p == 0:
            m //= p
            e += 1
        if e % 2 == 1:
            out *= p
        p += 1
    if m > 1:
        out *= m
    return sign * out


def factorize_int(n: int) -> dict[str, int]:
    if n == 0:
        return {"0": 1}
    fac = sp.factorint(abs(n))
    out = {str(int(p)): int(e) for p, e in fac.items()}
    if n < 0:
        out["-1"] = 1
    return out


def is_perfect_square(n: int) -> bool:
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def galois_cubic(b_poly: sp.Poly, delta3: int) -> str:
    """Galois group of splitting field of an irreducible Q-cubic."""
    if delta3 == 0:
        # would mean a repeated root => reducible over Q-bar but the
        # poly itself may still be irreducible over Q only if disc=0 and
        # gcd(b, b') is degree 1; but disc=0 implies common root in Q-bar
        # which for monic-rational case forces a rational root, hence
        # reducible. So delta3=0 + irreducible should not occur.
        return "degenerate"
    return "C_3" if is_perfect_square(delta3) else "S_3"


def trichotomy_bin(irreducible: bool, delta3: int, gal: str) -> str:
    if not irreducible:
        return "out"
    if delta3 > 0 and gal == "C_3":
        return "+_C3_real"
    if delta3 > 0 and gal == "S_3":
        return "+_S3_real"
    if delta3 < 0 and gal == "S_3":
        return "-_S3_CM"
    return "other"


def splitting_field_latex(delta3: int, gal: str) -> str:
    if gal == "C_3":
        return r"\text{real cyclic cubic, totally real}"
    if gal == "S_3":
        sub = fmt_field(delta3)
        if delta3 > 0:
            return rf"K_b\;[\mathbb{{Q}}\!:\!K_b]=6,\ \text{{quadratic resolvent }}{sub}\ \text{{real}}"
        return rf"K_b\;[\mathbb{{Q}}\!:\!K_b]=6,\ \text{{quadratic resolvent }}{sub}\ \text{{imag.}}"
    return "n/a"


# --------------------------------------------------------------------- CF eval

def evaluate_cf(b_coeffs: tuple[int, int, int, int],
                a_func=lambda n: mp.mpf(1),
                N: int = 500,
                dps: int = 200) -> mp.mpf:
    """
    Evaluate the partial continued fraction
        L_N = b(0) + a(1)/( b(1) + a(2)/( b(2) + ... + a(N)/b(N) ))
    where b(n) = a3 n^3 + a2 n^2 + a1 n + a0 and a(n) is supplied by
    the caller (default: a(n) = 1 for Session A).
    """
    a3, a2, a1, a0 = b_coeffs

    def b(n: int) -> mp.mpf:
        return mp.mpf(a3) * n ** 3 + mp.mpf(a2) * n ** 2 + mp.mpf(a1) * n + mp.mpf(a0)

    with mp.workdps(dps):
        x = b(N)
        for n in range(N - 1, -1, -1):
            an1 = a_func(n + 1)
            if x == 0:
                # avoid div-by-zero; nudge with tiny epsilon (degenerate)
                x = mp.mpf("1e-{}".format(dps // 2))
            x = b(n) + an1 / x
        return x


def cf_estimate(b_coeffs, dps: int = 200) -> dict[str, str]:
    """Compute L_N for N in {500, 1000, 2000} and a convergence label."""
    out: dict[str, Any] = {}
    L = {}
    for N in (500, 1000, 2000):
        L[N] = evaluate_cf(b_coeffs, N=N, dps=dps)
        out[f"L_{N}"] = mp.nstr(L[N], 40)
    # convergence rate proxy
    with mp.workdps(dps):
        d12 = abs(L[1000] - L[500])
        d23 = abs(L[2000] - L[1000])
        if d12 == 0 and d23 == 0:
            rate = "fast"
        elif d23 == 0:
            rate = "exponential"
        else:
            ratio = d23 / d12 if d12 > 0 else mp.mpf("inf")
            r = float(ratio)
            if r < 1e-6:
                rate = "exponential"
            elif r < 0.5:
                rate = "fast"
            else:
                rate = "linear"
            out["ratio_d23_over_d12"] = mp.nstr(ratio, 6)
    out["L_estimate"] = out["L_2000"]
    out["convergence_rate"] = rate
    return out


# --------------------------------------------------------------------- core

n_sym = sp.symbols("n")


def analyze_family(family_id: int,
                   coeffs: tuple[int, int, int, int],
                   include_cf: bool = True) -> dict[str, Any]:
    a3, a2, a1, a0 = coeffs
    b_expr = a3 * n_sym ** 3 + a2 * n_sym ** 2 + a1 * n_sym + a0
    b_poly = sp.Poly(b_expr, n_sym)

    factors = sp.factor_list(b_expr, n_sym, domain="QQ")
    # factors = (lead_coeff, [(factor_expr, multiplicity), ...])
    is_irred = (len(factors[1]) == 1 and factors[1][0][1] == 1
                and sp.degree(factors[1][0][0], n_sym) == 3)

    delta3 = int(sp.discriminant(b_expr, n_sym))
    delta3_fac = factorize_int(delta3)

    if not is_irred:
        gal = "reducible"
        split = "n/a"
        cm_field = "n/a"
        bin_ = "out"
        factorization_str = sp.latex(sp.factor(b_expr, n_sym))
    else:
        gal = galois_cubic(b_poly, delta3)
        split = splitting_field_latex(delta3, gal)
        if gal == "C_3":
            cm_field = "none (totally real)"
        else:
            # S_3 case: CM-field invariant per B3(ii)
            cm_field = fmt_field(delta3)
        bin_ = trichotomy_bin(True, delta3, gal)
        factorization_str = sp.latex(b_expr)

    rec: dict[str, Any] = {
        "family_id": family_id,
        "alpha_3": a3,
        "alpha_2": a2,
        "alpha_1": a1,
        "alpha_0": a0,
        "b_latex": sp.latex(b_expr),
        "irreducible": bool(is_irred),
        "factorization": factorization_str,
        "Delta_3_exact": delta3,
        "Delta_3_factorization": delta3_fac,
        "Delta_3_sign": int((delta3 > 0) - (delta3 < 0)),
        "Delta_3_is_square": is_perfect_square(delta3),
        "Galois_group": gal,
        "splitting_field": split,
        "CM_field": cm_field,
        "trichotomy_bin": bin_,
    }

    if include_cf:
        rec.update(cf_estimate(coeffs, dps=200))

    return rec


def enumerate_candidates():
    """Yield (a3, a2, a1, a0) in lex order with gcd == 1."""
    a3_vals = [1, 2, 3, 5, 7]
    small = [-3, -2, -1, 0, 1, 2, 3]
    for a3 in a3_vals:
        for a2 in small:
            for a1 in small:
                for a0 in small:
                    g = math.gcd(math.gcd(a3, abs(a2)), math.gcd(abs(a1), abs(a0)))
                    if g != 1:
                        continue
                    yield (a3, a2, a1, a0)


def is_irreducible_quick(coeffs):
    a3, a2, a1, a0 = coeffs
    b_expr = a3 * n_sym ** 3 + a2 * n_sym ** 2 + a1 * n_sym + a0
    facs = sp.factor_list(b_expr, n_sym, domain="QQ")[1]
    return len(facs) == 1 and facs[0][1] == 1 and sp.degree(facs[0][0], n_sym) == 3


# --------------------------------------------------------------------- main

def main():
    catalogue = []
    reducible_control = []
    fid = 0
    rid = 0
    for coeffs in enumerate_candidates():
        if is_irreducible_quick(coeffs):
            if fid >= 50:
                continue
            fid += 1
            rec = analyze_family(fid, coeffs, include_cf=True)
            catalogue.append(rec)
            print(f"  irred #{fid:2d}: a=({coeffs[0]},{coeffs[1]},{coeffs[2]},{coeffs[3]})  "
                  f"Delta3={rec['Delta_3_exact']:>9d}  Gal={rec['Galois_group']:<3s}  "
                  f"bin={rec['trichotomy_bin']}")
        else:
            if rid >= 10:
                continue
            rid += 1
            rec = analyze_family(1000 + rid, coeffs, include_cf=False)
            rec["control_id"] = rid
            reducible_control.append(rec)
            print(f"  reduc #{rid:2d}: a=({coeffs[0]},{coeffs[1]},{coeffs[2]},{coeffs[3]})  "
                  f"factorization={rec['factorization']}")
        if fid >= 50 and rid >= 10:
            break

    # ---- bin distribution
    bins = {"+_S3_real": 0, "+_C3_real": 0, "-_S3_CM": 0, "other": 0, "out": 0}
    cm_hist: dict[str, int] = {}
    for rec in catalogue:
        bins[rec["trichotomy_bin"]] = bins.get(rec["trichotomy_bin"], 0) + 1
        if rec["trichotomy_bin"] == "-_S3_CM":
            cm = rec["CM_field"]
            cm_hist[cm] = cm_hist.get(cm, 0) + 1

    # exotic: cyclic cubics with |Delta_3| <= 100
    exotic = [
        {"family_id": r["family_id"], "coeffs": [r["alpha_3"], r["alpha_2"], r["alpha_1"], r["alpha_0"]],
         "Delta_3": r["Delta_3_exact"]}
        for r in catalogue
        if r["Galois_group"] == "C_3" and abs(r["Delta_3_exact"]) <= 100
    ]

    summary = {
        "n_families": len(catalogue),
        "n_reducible_control": len(reducible_control),
        "bin_counts": bins,
        "CM_field_histogram_for_minus_S3_CM_bin": cm_hist,
        "exotic_small_disc_cyclic_cubics": exotic,
    }

    # ---- calibration anchors
    anchors = {}
    # (1) Apery: b = 34 n^3 + 51 n^2 + 27 n + 5  (out of enumeration window)
    anchor_apery = analyze_family(-1, (34, 51, 27, 5), include_cf=False)
    anchors["apery_b_n"] = {
        "expected": {"reducible": True, "factorization_q": "(2n+1)(17n^2+17n+5)",
                     "Delta_3": -459, "Galois": "C_2 (quadratic) for the irreducible cubic factor only"},
        "computed": {
            "irreducible": anchor_apery["irreducible"],
            "factorization": anchor_apery["factorization"],
            "Delta_3": anchor_apery["Delta_3_exact"],
            "Delta_3_factorization": anchor_apery["Delta_3_factorization"],
        },
        "match": (anchor_apery["irreducible"] is False
                  and anchor_apery["Delta_3_exact"] == -459),
    }
    # (2) Catalan-cubic: b = n^3 - 2  (in window: a3=1, a2=0, a1=0, a0=-2)
    anchor_cat = analyze_family(-2, (1, 0, 0, -2), include_cf=False)
    # disc(n^3-2) = -27 * (-2)^2 = -108
    anchors["catalan_seed_n3_minus_2"] = {
        "expected": {"Delta_3": -108, "Galois": "S_3", "CM_field_sqfree_disc": -3},
        "computed": {
            "Delta_3": anchor_cat["Delta_3_exact"],
            "Galois": anchor_cat["Galois_group"],
            "CM_field": anchor_cat["CM_field"],
            "trichotomy_bin": anchor_cat["trichotomy_bin"],
        },
        "match": (anchor_cat["Delta_3_exact"] == -108
                  and anchor_cat["Galois_group"] == "S_3"
                  and squarefree_part(-108) == -3
                  and anchor_cat["trichotomy_bin"] == "-_S3_CM"),
    }
    # (3) cyclic-C_3 probe: b = n^3 - 3 n + 1  (in window)
    anchor_c3 = analyze_family(-3, (1, 0, -3, 1), include_cf=False)
    # disc(n^3 + p n + q) = -4 p^3 - 27 q^2 = -4(-3)^3 - 27 = 108 - 27 = 81
    anchors["cyclic_C3_n3_minus_3n_plus_1"] = {
        "expected": {"Delta_3": 81, "is_square": True, "Galois": "C_3",
                     "splitting_field": r"\mathbb{Q}(\zeta_9 + \zeta_9^{-1})"},
        "computed": {
            "Delta_3": anchor_c3["Delta_3_exact"],
            "is_square": anchor_c3["Delta_3_is_square"],
            "Galois": anchor_c3["Galois_group"],
            "trichotomy_bin": anchor_c3["trichotomy_bin"],
            "splitting_field": anchor_c3["splitting_field"],
        },
        "match": (anchor_c3["Delta_3_exact"] == 81
                  and anchor_c3["Galois_group"] == "C_3"
                  and anchor_c3["trichotomy_bin"] == "+_C3_real"),
    }
    anchors["all_match"] = all(a["match"] for a in anchors.values() if isinstance(a, dict) and "match" in a)

    # ---- write outputs
    cat_path = os.path.join(HERE, "cubic_family_catalogue.json")
    sum_path = os.path.join(HERE, "galois_distribution_summary.json")
    anc_path = os.path.join(HERE, "calibration_anchors.json")
    claims_path = os.path.join(HERE, "claims.jsonl")

    with open(cat_path, "w", encoding="utf-8") as f:
        json.dump({"families": catalogue, "reducible_control": reducible_control},
                  f, indent=2, ensure_ascii=False, default=str)
    with open(sum_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=str)
    with open(anc_path, "w", encoding="utf-8") as f:
        json.dump(anchors, f, indent=2, ensure_ascii=False, default=str)

    # ---- AEAL claims
    def file_hash(p):
        return hashlib.sha256(open(p, "rb").read()).hexdigest()

    cat_hash = file_hash(cat_path)
    sum_hash = file_hash(sum_path)
    anc_hash = file_hash(anc_path)

    claims = [
        {"claim": f"Enumerated 50 Z-primitive irreducible cubics b(n) with a3 in {{1,2,3,5,7}}, a2,a1,a0 in {{-3..3}}, lex order; bin counts: {bins}",
         "evidence_type": "computation", "dps": 200, "reproducible": True,
         "script": "cubic_family_enumeration.py", "output_hash": cat_hash},
        {"claim": f"Galois distribution and CM-field histogram across the 50-family Session-A catalogue: {summary['CM_field_histogram_for_minus_S3_CM_bin']}",
         "evidence_type": "computation", "dps": 200, "reproducible": True,
         "script": "cubic_family_enumeration.py", "output_hash": sum_hash},
        {"claim": "Calibration: b=34n^3+51n^2+27n+5 is REDUCIBLE over Q with disc=-459=-3^3*17 (Apery, out-of-scope of B3(i)).",
         "evidence_type": "computation", "dps": 200, "reproducible": True,
         "script": "cubic_family_enumeration.py", "output_hash": anc_hash},
        {"claim": "Calibration: b=n^3-2 has disc=-108, Gal=S_3, CM-field=Q(sqrt(-3)), bin=-_S3_CM (Catalan seed).",
         "evidence_type": "computation", "dps": 200, "reproducible": True,
         "script": "cubic_family_enumeration.py", "output_hash": anc_hash},
        {"claim": "Calibration: b=n^3-3n+1 has disc=81=9^2, Gal=C_3, totally real cyclic cubic in Q(zeta_9+zeta_9^{-1}), bin=+_C3_real.",
         "evidence_type": "computation", "dps": 200, "reproducible": True,
         "script": "cubic_family_enumeration.py", "output_hash": anc_hash},
    ]
    with open(claims_path, "w", encoding="utf-8") as f:
        for c in claims:
            f.write(json.dumps(c, ensure_ascii=False) + "\n")

    print()
    print("=" * 72)
    print("Bin distribution:")
    for k, v in bins.items():
        print(f"  {k:<12s} : {v}")
    print(f"Reducible-control count: {len(reducible_control)}")
    print(f"CM-field histogram (-_S3_CM bin): {cm_hist}")
    print(f"Exotic small-disc cyclic cubics: {len(exotic)}")
    print()
    print("Calibration anchors:")
    for k in ("apery_b_n", "catalan_seed_n3_minus_2", "cyclic_C3_n3_minus_3n_plus_1"):
        a = anchors[k]
        print(f"  {k:<35s} match={a['match']}")
    print(f"  all_match = {anchors['all_match']}")
    print()
    print(f"Wrote: {cat_path}")
    print(f"Wrote: {sum_path}")
    print(f"Wrote: {anc_path}")
    print(f"Wrote: {claims_path}")

    if not anchors["all_match"]:
        print()
        print("HALT: at least one calibration anchor mismatched")
        with open(os.path.join(HERE, "halt_log.json"), "w", encoding="utf-8") as f:
            json.dump({"reason": "calibration anchor mismatch", "anchors": anchors},
                      f, indent=2, default=str)


if __name__ == "__main__":
    main()
