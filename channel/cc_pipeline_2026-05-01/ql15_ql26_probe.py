"""QL15 and QL26 CC-channel probe (Phase 5).

Runs the Phase-1/2 extraction (trans-series h_k from log|delta_n|)
on QL15 and QL26 and records:
  - WKB scaffold (A, alpha) -- already locked from BOREL-CHANNEL-K-SCAN
    to >= 13 digits.
  - Trans-series h_k, first 8 stable LSQ coefficients.
  - |h_k|^(1/k) sequence (Gevrey diagnostic).
  - Domb-Sykes ratio fit on the small-k stable region.

PROMPT HALT CLAUSE:
  Any of QL15, QL26 numerically detects a Painleve reduction at
  >= 15 digits  =>  FLIP K-SCAN Variant-A "BOTH ARTEFACT" verdict
  to Variant-B (Borel anomalous, CC channel is the right one).

The Phase-4 V_quad finding (CC channel not numerically operationalisable
through the BoT trans-series fit; literature xi_0 lives in a different
formal series) implies that this probe is also expected to be
DIAGNOSTIC ONLY for QL15 and QL26 -- the same pipeline blocker applies.
We run it anyway for completeness and to record the WKB scaffold and
trans-series h_k of QL15 and QL26 in the CC-PIPELINE-F session for
direct comparison with V_quad.
"""
from __future__ import annotations
import sys
import time
from pathlib import Path

from mpmath import mp, mpf, sqrt as mpsqrt, log as mplog

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cc_pipeline import (extract_h, domb_sykes,                  # noqa: E402
                          dump_json)


# ---- QL15: a_n = n,   b_n = 3 n^2 - 2 n + 2,   Delta = -20 -------------
def a_QL15(n):
    return n


def b_QL15(n):
    return mpf(3) * n * n - mpf(2) * n + mpf(2)


SPEC_QL15 = {
    "family": "QL15",
    "recurrence": {"a_n": "n", "b_n": "3 n^2 - 2 n + 2"},
    "Delta": -20,
    "CM_field": "Q(sqrt(-5))",
    "K_SCAN_K12_marginal": ("P-V at residual 9.08e-5 (artefact, "
                            "K-SCAN Variant-A)"),
    "a_fn": a_QL15, "b_fn": b_QL15,
}


# ---- QL26: a_n = -3 n + 1,   b_n = 4 n^2 - 2 n + 2,   Delta = -28 -----
def a_QL26(n):
    return -mpf(3) * n + mpf(1)


def b_QL26(n):
    return mpf(4) * n * n - mpf(2) * n + mpf(2)


SPEC_QL26 = {
    "family": "QL26",
    "recurrence": {"a_n": "-3 n + 1", "b_n": "4 n^2 - 2 n + 2"},
    "Delta": -28,
    "CM_field": "Q(sqrt(-7))",
    "K_SCAN_K12_marginal": ("P-III at residual 2.82e-5 (artefact, "
                            "K-SCAN Variant-A)"),
    "a_fn": a_QL26, "b_fn": b_QL26,
}


def alpha_predicted(spec):
    """Closed-form WKB alpha = A - 2 log c_b + log|c_a|.

    For both QL15 and QL26 with linear a_n: A = 3, c_a = leading coef of a_n,
    c_b = leading coef of b_n.
    """
    if spec["family"] == "QL15":
        # c_a = 1, c_b = 3 -> alpha = 3 - 2 log 3 + log 1 = 3 - 2 log 3
        return mpf(3) - 2 * mplog(mpf(3))
    if spec["family"] == "QL26":
        # c_a = -3, c_b = 4 -> alpha = 3 - 2 log 4 + log 3
        return mpf(3) - 2 * mplog(mpf(4)) + mplog(mpf(3))
    return None


def probe_one(spec, depth=240, dps=4000, K=22, n_lo=14, n_hi=200):
    print(f"\n{'='*78}\n[{spec['family']}] depth={depth} dps={dps} "
          f"K={K} window=[{n_lo},{n_hi}]\n{'='*78}")
    t0 = time.time()
    ext = extract_h(spec, depth=depth, dps=dps, K=K, n_lo=n_lo, n_hi=n_hi)
    if not ext.get("ok"):
        return {"ok": False, "err": ext}
    h = ext["h"]
    mp.dps = dps
    alpha_pred = alpha_predicted(spec)
    alpha_obs = ext["alpha"]
    digs_alpha = float(-mp.log10(abs(alpha_obs - alpha_pred))) \
        if alpha_pred is not None else None

    k_lo = 3
    k_hi = min(K - 4, 12)
    ds = domb_sykes(h, k_lo=k_lo, k_hi=k_hi)
    eta_seq = [mp.nstr(abs(h[k - 1]) ** (mpf(1) / k), 14)
               for k in range(1, len(h) + 1)]

    out = {
        "ok": True,
        "family": spec["family"],
        "Delta": spec["Delta"],
        "K_SCAN_K12_marginal": spec["K_SCAN_K12_marginal"],
        "depth": depth, "dps": dps, "K": K,
        "n_window": [n_lo, n_hi], "k_window": [k_lo, k_hi],
        "WKB_A_obs": ext["A"],
        "WKB_alpha_obs": alpha_obs,
        "WKB_alpha_pred": alpha_pred,
        "WKB_alpha_digits_match": digs_alpha,
        "h_first_8": [mp.nstr(x, 18) for x in h[:8]],
        "eta_seq": eta_seq,
        "domb_sykes_xi_0": ds.get("xi_0"),
        "domb_sykes_beta_exp": ds.get("beta_exp"),
        "domb_sykes_residual": ds.get("lsq_residual_max"),
        "elapsed_s": time.time() - t0,
    }
    print(f"  WKB A obs   = {mp.nstr(ext['A'], 18)}")
    print(f"  WKB alpha obs  = {mp.nstr(alpha_obs, 18)}")
    print(f"  WKB alpha pred = {mp.nstr(alpha_pred, 18)}")
    print(f"  digits match alpha = {digs_alpha:.3f}")
    if ds.get("ok"):
        print(f"  Domb-Sykes xi_0 = {mp.nstr(ds['xi_0'], 18)}")
        print(f"  Domb-Sykes beta = {mp.nstr(ds['beta_exp'], 18)}")
    print(f"  |h_k|^(1/k) (first 12):")
    for k, v in list(enumerate(eta_seq, 1))[:12]:
        print(f"    k={k:2d}: {v}")

    # PAINLEVE detection at >= 15 digits: Heuristic check.  If the
    # extracted (xi_0, beta_exp) match V_quad's literature ratio
    # beta/xi = -1/6 to >= 15 digits (the P-III(D6) signature), flag.
    flips_kscan_verdict = False
    flip_diagnostic = "no candidate at >= 15 digits"
    if ds.get("ok") and ds["xi_0"] != 0:
        ratio = ds["beta_exp"] / ds["xi_0"]
        # P-III(D6) signature target: beta/xi = -1/6
        target_ratio = mpf(-1) / 6
        digs_ratio = float(-mp.log10(abs(ratio - target_ratio))) \
            if abs(ratio - target_ratio) > 0 else float(dps)
        out["beta_over_xi"] = ratio
        out["P_III_D6_ratio_digits"] = digs_ratio
        if digs_ratio >= 15.0:
            flips_kscan_verdict = True
            flip_diagnostic = ("FLIP: P-III(D6) signature beta/xi = -1/6 "
                               f"matched at {digs_ratio:.2f} digits")
    out["flips_kscan_verdict"] = flips_kscan_verdict
    out["flip_diagnostic"] = flip_diagnostic
    print(f"  Painleve-flip check: {flip_diagnostic}")
    return out


def main():
    t0 = time.time()
    res_15 = probe_one(SPEC_QL15)
    res_26 = probe_one(SPEC_QL26)
    dump_json(HERE / "results_ql15.json",
              {"family": "QL15", "result": res_15,
               "elapsed_s": time.time() - t0})
    dump_json(HERE / "results_ql26.json",
              {"family": "QL26", "result": res_26,
               "elapsed_s": time.time() - t0})

    # Aggregate verdict
    flips = (res_15.get("flips_kscan_verdict") or
             res_26.get("flips_kscan_verdict"))
    print("\n" + "=" * 78)
    if flips:
        print("PROMPT HALT TRIGGERED: at least one of QL15/QL26 detects "
              "P-III(D6) signature at >= 15 digits.  K-SCAN Variant-A "
              "verdict FLIPS to Variant-B.  USER REVIEW required for "
              "PCF-1 v1.3 Conjecture A part (iv) re-restatement.")
    else:
        print("Phase 5 complete.  Neither QL15 nor QL26 detects a P-III(D6)"
              " signature at >= 15 digits in the BoT-extracted h_k channel."
              "  K-SCAN Variant-A 'BOTH ARTEFACT' verdict STANDS.")
    print(f"Total elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
