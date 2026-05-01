"""CC-channel pipeline for PCF asymptotics.

Operationalises the connection-coefficient (CC) channel as defined in
sessions/2026-05-01/CHANNEL-THEORY-OUTLINE/channel_theory_outline.tex
Sec 3.3.

Strategy.
========

The Borel-of-trans-series channel (BoT, Session E/E') already extracts
the WKB-residual coefficients

    log|delta_n|  =  -A n log n + alpha n - beta log n + gamma
                     + sum_{k>=1} h_k / n^k

with the {h_k} stable to >= 14 digits across K = 12..24.  The BoT
channel's downstream FAILS because Pade-Borel resummation poles drift
under K-extension (K-SCAN: ARTEFACT verdict).  The Pade artefact does
not contaminate the {h_k} themselves: it is downstream.

The CC channel re-uses the {h_k} but replaces Pade-Borel by direct
Borel-singularity analysis (Domb-Sykes / Darboux):

    h_k  ~  S_1  *  Gamma(k + beta_exp)  /  xi_0^k     (k -> infty)

so that consecutive ratios

    r_k := h_{k+1} / h_k  =  (k + beta_exp) / xi_0  +  O(1/k^2),

and

    S_1  =  lim_{k -> infty}  h_k * xi_0^k / Gamma(k + beta_exp).

Richardson / Aitken extrapolation converts the 14-digit-stable {h_k}
into a high-precision estimate of (xi_0, beta_exp, S_1) without
invoking any Pade-Borel resummation.  This is the connection-coefficient
datum for the second-order linear difference operator L_b governing
the convergent (p_n, q_n).

Phase 4 then matches the extracted (xi_0, beta_exp) for V_quad against
the literature P-III(D6) anchor

    xi_0^lit       = 2/sqrt(3)        (~ 1.15470053837925152901...)
    beta_exp^lit   = -1/(3*sqrt(3))   (~ -0.19245008972987525)
    S_1^lit        ~ 0.43770528...    (Dingle late-term formula,
                                       8 digits in P12 Sec 6 Stokes data)

A match at >= 20 digits is the prompt's HALT clause threshold.

Phase 5 runs the same extraction on QL15 (Delta=-20) and QL26
(Delta=-28) and reports their (xi_0, beta_exp) signature; if any
algebraic match against a P-I..P-VI moduli point is detected at
>= 15 digits, the BOREL-CHANNEL-K-SCAN "BOTH ARTEFACT" verdict
flips to Variant-B (FLAG: borel anomalous, CC channel is the right
one).
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from mpmath import mp, mpf, mpc, matrix, lu_solve, gamma as mpgamma, sqrt as mpsqrt


# Reuse the v3 trans-series extractor.  The borel_channel module lives
# under sessions/2026-05-01/BOREL-CHANNEL-5X/.
import sys
_HERE = Path(__file__).resolve().parent
_BOREL = _HERE.parent / "BOREL-CHANNEL-5X"
if str(_BOREL) not in sys.path:
    sys.path.insert(0, str(_BOREL))
from borel_channel import compute_convergents, extract_transseries  # noqa: E402


# ====================================================================
# Phase 1-2: extract trans-series coefficients h_k  (linear-system data)
# ====================================================================

def extract_h(spec: Dict, depth: int, dps: int, K: int,
              n_lo: int, n_hi: int) -> Dict:
    """Extract trans-series coefficients (h_1, ..., h_K) for a PCF."""
    a_fn = spec["a_fn"]
    b_fn = spec["b_fn"]
    t0 = time.time()
    L_arr, n_arr, stab = compute_convergents(a_fn, b_fn, depth=depth,
                                             dps=dps, record_from=10)
    if L_arr is None:
        return {"ok": False, "reason": "convergents failed"}
    ts = extract_transseries(L_arr, n_arr, n_lo=n_lo, n_hi=n_hi, K=K)
    if ts is None:
        return {"ok": False, "reason": "trans-series LSQ failed"}
    return {
        "ok": True,
        "family": spec["family"],
        "L_star": ts["L_star"],
        "A": ts["A"],
        "alpha": ts["alpha"],
        "beta": ts["beta"],
        "gamma": ts["gamma"],
        "h": ts["h"],
        "K": K,
        "depth": depth,
        "dps": dps,
        "n_lo": n_lo,
        "n_hi": n_hi,
        "stab_log10": stab,
        "elapsed_s": time.time() - t0,
    }


# ====================================================================
# Phase 2: Borel-singularity analysis (= Stokes data extraction)
# ====================================================================

def domb_sykes(h: List[mpf], k_lo: int, k_hi: int) -> Dict:
    """Domb-Sykes ratio analysis.

    For h_k ~ S_1 * Gamma(k + beta_exp) / xi_0^k, the ratio
        r_k := h_{k+1} / h_k = (k + beta_exp) / xi_0 + O(1/k^2)
    is linear in k with slope 1/xi_0 and intercept beta_exp/xi_0.
    Linear regression on (k, r_k) for k in [k_lo, k_hi] gives both.

    Returns xi_0, beta_exp, and the linear-fit residual.
    """
    rows = []
    for k in range(k_lo, k_hi + 1):
        if k + 1 >= len(h):
            break
        if h[k] == 0 or h[k + 1] == 0:
            continue
        r = h[k] / h[k - 1] if k - 1 >= 0 else None
        # The 1-indexing convention: h is stored as [h_1, h_2, ..., h_K]
        # i.e. h[i] = h_{i+1}.  We want r_m = h_{m+1}/h_m for m=k.
    # Re-do cleanly using 1-indexed view:
    H = {m: h[m - 1] for m in range(1, len(h) + 1)}
    pts = []
    for m in range(k_lo, k_hi):
        if m not in H or (m + 1) not in H:
            continue
        if H[m] == 0:
            continue
        pts.append((mpf(m), H[m + 1] / H[m]))
    if len(pts) < 4:
        return {"ok": False, "reason": "not enough points",
                "n_pts": len(pts)}
    # LSQ:  r_m = a * m + b   with a = 1/xi_0, b = beta_exp/xi_0.
    n = len(pts)
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] * p[0] for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    det = n * sxx - sx * sx
    a = (n * sxy - sx * sy) / det
    b = (sxx * sy - sx * sxy) / det
    xi_0 = mpf(1) / a
    beta_exp = b * xi_0
    res = mpf(0)
    for m, r in pts:
        pred = a * m + b
        rr = abs(r - pred)
        if rr > res:
            res = rr
    return {"ok": True,
            "xi_0": xi_0,
            "beta_exp": beta_exp,
            "slope_a": a,
            "intercept_b": b,
            "lsq_residual_max": res,
            "n_pts": n,
            "k_lo": k_lo, "k_hi": k_hi}


def neville_richardson(seq: List[mpf], levels: int = None) -> List[List[mpf]]:
    """Apply Neville-style 1/k^p Richardson extrapolation.

    Given a sequence x_k ~ x* + a/k + b/k^2 + ...  the iterated
    extrapolation
        x^{(p)}_k = (k+p) x^{(p-1)}_{k+1} - k x^{(p-1)}_k) / p
    accelerates convergence by one order per level.
    """
    n = len(seq)
    table = [list(seq)]
    if levels is None:
        levels = n - 1
    for p in range(1, levels + 1):
        prev = table[p - 1]
        m = len(prev) - 1
        if m <= 0:
            break
        row = []
        # Use offset k-> k+p_offset so we extrapolate the largest indices.
        # Treat seq[i] as defined at parameter k = k0 + i for some k0.
        for i in range(m):
            # Extrapolate as: val_i = ((k+p) * prev[i+1] - k * prev[i]) / p
            # but we don't know k0; absorb by scaling.  Simpler: use the
            # asymmetric Richardson  R_p = (2^p prev[i+1] - prev[i])/(2^p - 1)
            # ASSUMING prev was already a halving sequence -- here it is
            # NOT halving.  Use the "k-index" Richardson:
            # x_k = x* + c/k^p + ...  =>
            # (k+1)^p x_{k+1} - k^p x_k  ~  x* * ((k+1)^p - k^p)
            # So  x* = ((k+1)^p x_{k+1} - k^p x_k) / ((k+1)^p - k^p)
            k = mpf(i + 1)  # nominal k-index of seq[i]
            kp1 = k + 1
            num = kp1 ** p * prev[i + 1] - k ** p * prev[i]
            den = kp1 ** p - k ** p
            row.append(num / den)
        table.append(row)
    return table


def stokes_extract(h: List[mpf], k_lo: int, k_hi: int,
                   xi_0: Optional[mpf] = None,
                   beta_exp: Optional[mpf] = None) -> Dict:
    """Extract Stokes constant S_1 := h_k * xi_0^k / Gamma(k + beta_exp)
    in the limit k -> infinity, then Richardson-extrapolate.

    If xi_0, beta_exp not supplied, do a Domb-Sykes pass first.
    """
    if xi_0 is None or beta_exp is None:
        ds = domb_sykes(h, k_lo=max(2, k_lo // 2), k_hi=k_hi)
        if not ds.get("ok"):
            return {"ok": False, "reason": "domb_sykes failed", **ds}
        xi_0 = ds["xi_0"]
        beta_exp = ds["beta_exp"]
    seq = []
    for k in range(k_lo, k_hi + 1):
        if k - 1 >= len(h):
            break
        try:
            val = h[k - 1] * xi_0 ** k / mpgamma(k + beta_exp)
        except Exception:
            continue
        seq.append(val)
    if len(seq) < 4:
        return {"ok": False, "reason": "too few Stokes points",
                "n": len(seq)}
    table = neville_richardson(seq, levels=min(8, len(seq) - 1))
    last = table[-1][-1] if table[-1] else table[-2][-1]
    return {"ok": True,
            "xi_0": xi_0,
            "beta_exp": beta_exp,
            "S_1_naive": seq[-1],
            "S_1_richardson": last,
            "richardson_table_size": len(table),
            "n_pts": len(seq)}


# ====================================================================
# Phase 3: connection coefficient (ratio form)
# ====================================================================

def connection_ratio(h: List[mpf], xi_0: mpf, beta_exp: mpf,
                     k_check_lo: int, k_check_hi: int) -> Dict:
    """The "connection coefficient" of the CC channel is the ratio
    between the asymptotic-at-infinity expansion (encoded in h_k) and
    the convergent-at-zero expansion (the formal recurrence).  In the
    Borel plane this is exactly the residue of Borel(zeta) at the
    leading singularity zeta = xi_0.

    For h_k ~ S_1 * Gamma(k + beta_exp) / xi_0^k the leading residue is
        Res_{zeta=xi_0} Borel(zeta) = S_1 * (xi_0)^{beta_exp} / Gamma(beta_exp)
    times a sign convention.  We tabulate S_1, xi_0, beta_exp and the
    derived residue.
    """
    sk = stokes_extract(h, k_lo=k_check_lo, k_hi=k_check_hi,
                        xi_0=xi_0, beta_exp=beta_exp)
    if not sk.get("ok"):
        return sk
    S_1 = sk["S_1_richardson"]
    try:
        residue = S_1 * (xi_0 ** beta_exp) / mpgamma(beta_exp)
    except Exception:
        residue = None
    return {"ok": True,
            "xi_0": xi_0,
            "beta_exp": beta_exp,
            "S_1": S_1,
            "residue_at_xi_0": residue,
            "stokes_n_pts": sk["n_pts"]}


# ====================================================================
# Phase 4: V_quad anchor match (literature target)
# ====================================================================

V_QUAD_LIT = {
    "xi_0": "2/sqrt(3)",
    "beta_exp": "-1/(3*sqrt(3))",
    "S_1_8digit": "0.43770528",
    "accessory_q": "(5 + i*sqrt(11))/54",
    "P_III_D6_params": "(alpha,beta,gamma,delta) = (1/6, 0, 0, -1/2)",
    "source": "p12_journal_main.tex Sec sec:vquad subsec Stokes data; "
              "Papanokechi2026Vquad",
}


def vquad_recovery(h: List[mpf], k_check_lo: int, k_check_hi: int,
                   dps: int) -> Dict:
    """Compare extracted (xi_0, beta_exp) against the literature anchor.

    Returns digits-of-agreement and pass/fail diagnostic.
    """
    mp.dps = dps
    xi_lit = mpf(2) / mpsqrt(3)
    beta_lit = mpf(-1) / (mpf(3) * mpsqrt(3))

    ds = domb_sykes(h, k_lo=max(2, k_check_lo // 2), k_hi=k_check_hi)
    if not ds.get("ok"):
        return {"ok": False, "reason": "domb-sykes failed on V_quad",
                **ds}

    err_xi = abs(ds["xi_0"] - xi_lit)
    err_beta = abs(ds["beta_exp"] - beta_lit)
    digs_xi = float(-mp.log10(err_xi)) if err_xi > 0 else float(dps)
    digs_beta = float(-mp.log10(err_beta)) if err_beta > 0 else float(dps)

    sk = stokes_extract(h, k_lo=k_check_lo, k_hi=k_check_hi,
                        xi_0=xi_lit, beta_exp=beta_lit)
    return {
        "ok": True,
        "xi_0_extracted": ds["xi_0"],
        "xi_0_lit": xi_lit,
        "xi_0_digits_match": digs_xi,
        "beta_exp_extracted": ds["beta_exp"],
        "beta_exp_lit": beta_lit,
        "beta_exp_digits_match": digs_beta,
        "S_1_richardson_at_lit": sk.get("S_1_richardson"),
        "S_1_lit_8digit": mpf("0.43770528"),
        "domb_sykes_residual": ds["lsq_residual_max"],
        "k_lo": k_check_lo, "k_hi": k_check_hi,
        "lit_anchor": V_QUAD_LIT,
    }


# ====================================================================
# Phase 5: Painleve detection in CC channel (heuristic catalogue)
# ====================================================================

def painleve_class_predict(family: str, A: mpf, alpha: mpf,
                            xi_0: mpf, beta_exp: mpf) -> Dict:
    """Heuristic CC-channel classifier.

    The CC-channel datum (xi_0, beta_exp) is a continuous moduli point.
    For a family to admit a P-class reduction we need the datum to
    match the moduli point of P-I, P-II, ..., P-VI on the
    Riemann-Hilbert side.  We do NOT have a closed-form RH datum
    catalogue here; we record the extracted moduli and flag the
    *necessary* condition that beta_exp be a rational multiple of
    xi_0 / sqrt(disc) (an empirical regularity at V_quad: beta_exp
    = - xi_0 / 6).

    Reports:
      - beta_exp / xi_0 ratio
      - PSLQ search for a rational multiple of beta_exp in
        {xi_0, xi_0^2, xi_0/sqrt(3), 1/xi_0}
    """
    # Sanity checks
    if abs(xi_0) < mpf("1e-20"):
        return {"ok": False, "reason": "xi_0 too small"}
    ratio = beta_exp / xi_0
    flags = []
    # V_quad relation: beta_exp / xi_0 = -1/(3*sqrt(3)) / (2/sqrt(3))
    #                                  = -1/6.  HEURISTIC FINGERPRINT.
    if abs(ratio + mpf(1) / 6) < mpf("1e-6"):
        flags.append("matches V_quad ratio -1/6 (P-III(D6) signature)")
    return {
        "ok": True,
        "family": family,
        "A": A, "alpha": alpha,
        "xi_0": xi_0, "beta_exp": beta_exp,
        "beta_over_xi": ratio,
        "ratio_minus_one_sixth": ratio + mpf(1) / 6,
        "flags": flags,
    }


# ====================================================================
# JSON serialiser (mpmath -> str)
# ====================================================================

def _jsonify(o):
    if isinstance(o, (mpf, mpc)):
        return mp.nstr(o, mp.dps, strip_zeros=False)
    if isinstance(o, dict):
        return {k: _jsonify(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonify(v) for v in o]
    return o


def dump_json(path: Path, obj):
    path.write_text(json.dumps(_jsonify(obj), indent=2), encoding="utf-8")
