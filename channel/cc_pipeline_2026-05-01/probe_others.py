"""Extend the Phase-5 probe to QL01, QL02, QL06 so all 6 Delta<0
degree-2 families appear in cc_channel_table.tex."""
from __future__ import annotations
import sys
import time
from pathlib import Path

from mpmath import mp, mpf, log as mplog

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cc_pipeline import extract_h, domb_sykes, dump_json  # noqa: E402


# QL01: a_n = n,  b_n = n^2 + n + 1,  Delta = -3
def a_QL01(n):
    return n


def b_QL01(n):
    return n * n + n + mpf(1)


# QL02: a_n = 2n+1,  b_n = n^2 + n + 1,  Delta = -4
def a_QL02(n):
    return mpf(2) * n + mpf(1)


def b_QL02(n):
    return n * n + n + mpf(1)


# QL06: a_n = n,  b_n = 2 n^2 - n + 1,  Delta = -7
def a_QL06(n):
    return n


def b_QL06(n):
    return mpf(2) * n * n - n + mpf(1)


SPECS = [
    {"family": "QL01", "Delta": -3, "CM_field": "Q(sqrt(-3))",
     "K_SCAN_K12_marginal": "n/a (not in K-scan)",
     "a_fn": a_QL01, "b_fn": b_QL01,
     "alpha_pred": "3.0",       # A=3, c_a=1, c_b=1
     "alpha_pred_fn": lambda: mpf(3)},
    {"family": "QL02", "Delta": -4, "CM_field": "Q(i)",
     "K_SCAN_K12_marginal": "n/a (not in K-scan)",
     "a_fn": a_QL02, "b_fn": b_QL02,
     "alpha_pred": "3 + log 2",  # A=3, c_a=2, c_b=1 -> 3 - 0 + log 2
     "alpha_pred_fn": lambda: mpf(3) + mplog(mpf(2))},
    {"family": "QL06", "Delta": -7, "CM_field": "Q(sqrt(-7))",
     "K_SCAN_K12_marginal": "n/a (not in K-scan)",
     "a_fn": a_QL06, "b_fn": b_QL06,
     "alpha_pred": "3 - 2 log 2",  # A=3, c_a=1, c_b=2
     "alpha_pred_fn": lambda: mpf(3) - 2 * mplog(mpf(2))},
]


def probe(spec, depth=240, dps=4000, K=22, n_lo=14, n_hi=200):
    print(f"\n[{spec['family']}]")
    t0 = time.time()
    ext = extract_h(spec, depth=depth, dps=dps, K=K, n_lo=n_lo, n_hi=n_hi)
    if not ext.get("ok"):
        return {"ok": False, "err": ext}
    mp.dps = dps
    alpha_pred = spec["alpha_pred_fn"]()
    digs_alpha = float(-mp.log10(abs(ext["alpha"] - alpha_pred)))
    h = ext["h"]
    ds = domb_sykes(h, k_lo=3, k_hi=min(K - 4, 12))
    ratio = ds["beta_exp"] / ds["xi_0"] if ds.get("ok") and ds["xi_0"] != 0 else None
    target_ratio = mpf(-1) / 6
    digs_ratio = (float(-mp.log10(abs(ratio - target_ratio)))
                  if ratio is not None and abs(ratio - target_ratio) > 0
                  else (float(dps) if ratio is not None else None))
    flips = bool(digs_ratio is not None and digs_ratio >= 15.0)
    out = {
        "ok": True,
        "family": spec["family"],
        "Delta": spec["Delta"],
        "CM_field": spec["CM_field"],
        "K_SCAN_K12_marginal": spec["K_SCAN_K12_marginal"],
        "depth": depth, "dps": dps, "K": K,
        "WKB_A_obs": ext["A"],
        "WKB_alpha_obs": ext["alpha"],
        "WKB_alpha_pred": alpha_pred,
        "WKB_alpha_pred_str": spec["alpha_pred"],
        "WKB_alpha_digits_match": digs_alpha,
        "h_first_8": [mp.nstr(x, 18) for x in h[:8]],
        "domb_sykes_xi_0": ds.get("xi_0"),
        "domb_sykes_beta_exp": ds.get("beta_exp"),
        "beta_over_xi": ratio,
        "P_III_D6_ratio_digits": digs_ratio,
        "flips_kscan_verdict": flips,
        "elapsed_s": time.time() - t0,
    }
    print(f"  alpha digits = {digs_alpha:.3f}; ratio digits to -1/6 = "
          f"{digs_ratio if digs_ratio is None else f'{digs_ratio:.3f}'}; "
          f"flip = {flips}")
    return out


def main():
    t0 = time.time()
    for spec in SPECS:
        res = probe(spec)
        path = HERE / f"results_{spec['family'].lower()}.json"
        dump_json(path, {"family": spec["family"], "result": res,
                          "elapsed_s": time.time() - t0})
        print(f"Wrote {path.name}")
    print(f"\nTotal elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
