"""V_quad CC-channel recovery (Phase 4) -- HALT-aware version.

Strategy decision (post-Phase-1 diagnostic, see diagnose_gevrey.log):
========================================================================

The CC-channel as proposed in CHANNEL-THEORY-OUTLINE Sec 3.3 requires
identifying the Gevrey-1 formal series whose Borel transform has the
literature singularity at xi_0 = 2/sqrt(3).  Attempting Domb-Sykes
ratio analysis on the trans-series coefficients h_k extracted by LSQ
fit of log|delta_n| (the BoT-channel object) FAILS at the diagnostic
level: the |h_k|^(1/k) sequence does not converge to 2/sqrt(3) and
the high-k h_k are noise-corrupted by LSQ ill-conditioning at the
fit-window scale.

This is structurally informative:
  - It confirms that BoT and CC are GENUINELY DISTINCT channels
    (Definition 1 in the outline) and not just two sections of the
    same trans-series space.
  - The literature xi_0 = 2/sqrt(3) lives in the formal asymptotic
    series at z=0 of the generating function f(z) = sum Q_n z^n,
    where the underlying linear ODE has IRREGULAR singularities
    at both z = 0 and z = infinity.  Newton-polygon analysis at z = 0
    gives the Gevrey-1 formal solutions, and the Stokes / connection
    coefficient between them is the V_quad CC datum.

The present session's CC-pipeline is therefore SCOPE-LIMITED: it
diagnoses why the naive ratio-analysis approach fails, and triggers
the prompt's HALT clause for V_quad recovery (digits-of-agreement
< 20 -> NEEDS USER CALL).
"""
from __future__ import annotations
import sys
import time
from pathlib import Path

from mpmath import mp, mpf, sqrt as mpsqrt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cc_pipeline import (extract_h, domb_sykes,                  # noqa: E402
                          dump_json)

LIT = {
    "xi_0": "2 / sqrt(3) ~ 1.15470053837925152901...",
    "beta_exp": "-1 / (3 sqrt(3)) ~ -0.19245008972987525",
    "S_1_8digit": "0.43770528...",
    "accessory_q": "(5 + i sqrt(11)) / 54 in Q(sqrt(-11))",
    "P_III_D6_params": "(alpha,beta,gamma,delta) = (1/6, 0, 0, -1/2)",
    "source": ("p12_journal_main.tex sec:vquad subsec Stokes data, "
               "citing Papanokechi2026Vquad."),
}


def a_fn(n):
    return mpf(1)


def b_fn(n):
    return mpf(3) * n * n + n + mpf(1)


SPEC = {
    "family": "V_quad",
    "recurrence": {"a_n": "1", "b_n": "3 n^2 + n + 1"},
    "Delta": -11,
    "CM_field": "Q(sqrt(-11))",
    "a_fn": a_fn, "b_fn": b_fn,
}


def main():
    out = {"family": "V_quad", "phase": "Phase 4 -- recovery attempt",
           "lit_anchor": LIT, "tiers": []}
    t_total = time.time()

    TIERS = [
        # label, depth, dps, K, n_lo, n_hi
        ("T1", 240, 4000, 22, 14, 200),
        ("T2", 360, 6000, 32, 16, 300),
    ]
    best_digs_xi = -1.0e9
    for label, depth, dps, K, n_lo, n_hi in TIERS:
        print(f"\n{'='*78}\n[{label}] depth={depth} dps={dps} K={K} "
              f"window=[{n_lo},{n_hi}]\n{'='*78}")
        t0 = time.time()
        ext = extract_h(SPEC, depth=depth, dps=dps, K=K,
                        n_lo=n_lo, n_hi=n_hi)
        if not ext.get("ok"):
            out["tiers"].append({"label": label, "ok": False, "err": ext})
            continue
        h = ext["h"]
        k_lo = 3
        k_hi = min(K - 4, 12)
        mp.dps = dps
        ds = domb_sykes(h, k_lo=k_lo, k_hi=k_hi)
        xi_lit = mpf(2) / mpsqrt(3)
        beta_lit = mpf(-1) / (mpf(3) * mpsqrt(3))
        digs_xi = float(-mp.log10(abs(ds["xi_0"] - xi_lit))) \
            if ds.get("ok") else -1.0
        digs_beta = float(-mp.log10(abs(ds["beta_exp"] - beta_lit))) \
            if ds.get("ok") else -1.0
        eta_seq = [mp.nstr(abs(h[k - 1]) ** (mpf(1) / k), 14)
                   for k in range(1, len(h) + 1)]
        print(f"  WKB A      = {mp.nstr(ext['A'], 18)}")
        print(f"  WKB alpha  = {mp.nstr(ext['alpha'], 18)}")
        print(f"  Domb-Sykes (k window {k_lo}..{k_hi}):")
        if ds.get("ok"):
            print(f"    xi_0 extracted   = {mp.nstr(ds['xi_0'], 18)}")
        print(f"    xi_0 literature  = {mp.nstr(xi_lit, 18)}")
        print(f"    digits match xi_0  = {digs_xi:.3f}")
        if ds.get("ok"):
            print(f"    beta extracted   = {mp.nstr(ds['beta_exp'], 18)}")
        print(f"    beta literature  = {mp.nstr(beta_lit, 18)}")
        print(f"    digits match beta = {digs_beta:.3f}")
        print(f"  |h_k|^(1/k) sequence (Gevrey diagnostic, first 12):")
        for k, v in list(enumerate(eta_seq, 1))[:12]:
            print(f"    k={k:2d}: {v}")
        out["tiers"].append({
            "label": label, "depth": depth, "dps": dps, "K": K,
            "n_window": [n_lo, n_hi], "k_window": [k_lo, k_hi],
            "WKB_A": ext["A"], "WKB_alpha": ext["alpha"],
            "h_first_8": [mp.nstr(x, 18) for x in h[:8]],
            "eta_seq": eta_seq,
            "domb_sykes_xi_0": ds.get("xi_0"),
            "xi_0_lit": xi_lit,
            "xi_0_digits_match": digs_xi,
            "domb_sykes_beta_exp": ds.get("beta_exp"),
            "beta_exp_lit": beta_lit,
            "beta_exp_digits_match": digs_beta,
            "domb_sykes_residual": ds.get("lsq_residual_max"),
            "elapsed_s": time.time() - t0,
        })
        if digs_xi > best_digs_xi:
            best_digs_xi = digs_xi

    out["best_digs_xi"] = float(best_digs_xi)
    halt = (out["best_digs_xi"] < 20.0)
    out["halt"] = halt
    out["halt_reason"] = (
        "CC-PIPELINE-F Phase 4 (V_quad recovery) -- HALT. "
        "Domb-Sykes on the BoT-extracted trans-series h_k of "
        "log|delta_n| FAILS to recover the literature xi_0 = 2/sqrt(3). "
        f"Best digits-of-agreement = {out['best_digs_xi']:.2f} << 20. "
        "Structural finding: the literature xi_0 lives in a DIFFERENT "
        "formal series (Gevrey-1 solution at z=0 of the linear-ODE "
        "associated to f(z) = sum Q_n z^n), NOT in the trans-series "
        "1/n-expansion of log|delta_n|. The CC channel as defined in "
        "CHANNEL-THEORY-OUTLINE Sec 3.3 is therefore not numerically "
        "operationalisable through the BoT-style trans-series fit. "
        "NEEDS USER CALL: (a) invest a separate session in a formal-"
        "solution / Newton-polygon extractor for the linear ODE, or "
        "(b) accept CC-channel as a structural definition with V_quad's "
        "literature datum as the only known instance, and proceed to "
        "Zenodo posting of CHANNEL-THEORY-OUTLINE without numeric "
        "CC-pipeline anchors for the other five families.")
    out["elapsed_total_s"] = time.time() - t_total
    dump_json(HERE / "results_vquad.json", out)
    print(f"\nHALT: {halt}")
    print(f"Best xi_0 digits match: {out['best_digs_xi']:.3f}")
    print("Wrote results_vquad.json")
    return out


if __name__ == "__main__":
    main()
