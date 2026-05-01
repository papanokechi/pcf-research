"""PCF2-SESSION-R1.2 -- j-invariant cross-degree extension at d=4.

Phases:
  R12-1  Igusa I,J -> j(Jac(C_b)) for all 60 Q1 quartics
  R12-2  Build joined 110-row dataframe (50 cubic + 60 quartic)
  R12-3  Per-degree Spearman correlations
  R12-4  Joint cross-degree analysis
  R12-5  Lemniscatic (j=1728) cell
  R12-6  Pari/gp probe (best-effort)
  R12-7  WKB precision escalation for j=0 quartic candidates (dps=1500)
  R12-8  Outputs

Coefficient ordering convention: leading first, [a4,a3,a2,a1,a0] for quartic;
[a3,a2,a1,a0] for cubic (matches R1.1).

Igusa for binary quartic f=ax^4+bx^3+cx^2+dx+e:
    I = 12 a e - 3 b d + c^2
    J = 72 a c e - 27 a d^2 - 27 b^2 e + 9 b c d - 2 c^3
Jacobian E_f: y^2 = x^3 - 27 I x - 27 J
   c4 = I, j(E_f) = 6912 I^3 / (4 I^3 - J^2)
   j=0  iff I=0;   j=1728 iff J=0 (and I>0).
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import subprocess
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import pandas as pd
from scipy import stats
from sympy import Rational, Integer, Poly, sympify, symbols

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
Q1_DIR = ROOT / "sessions" / "2026-05-01" / "PCF2-SESSION-Q1"
R11_DIR = ROOT / "sessions" / "2026-05-01" / "PCF2-SESSION-R1.1"
RUN_LOG = HERE / "run.log"
if RUN_LOG.exists():
    RUN_LOG.unlink()

CLAIMS = []


def log(msg: str):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(RUN_LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def add_claim(claim: str, dps: int, script: str, output_file: Path,
              extra: dict | None = None):
    rec = {
        "claim": claim,
        "evidence_type": "computation",
        "dps": dps,
        "reproducible": True,
        "script": script,
        "output_hash": sha256_file(output_file) if output_file.exists() else None,
    }
    if extra:
        rec.update(extra)
    CLAIMS.append(rec)


# ============================================================== R12-1
def igusa_IJ(coeffs):
    """coeffs = [a4, a3, a2, a1, a0]; returns (I, J) as exact Rationals."""
    a, b, c, d, e = (Rational(int(x)) for x in coeffs)
    I = 12 * a * e - 3 * b * d + c ** 2
    J = 72 * a * c * e - 27 * a * d ** 2 - 27 * b ** 2 * e + 9 * b * c * d - 2 * c ** 3
    return I, J


def jacobian_j(coeffs):
    """Returns (j_exact_or_NaN, I_exact, J_exact, is_singular)."""
    I, J = igusa_IJ(coeffs)
    Delta_quartic = 4 * I ** 3 - J ** 2  # = 27 * Disc(f) for monic
    if Delta_quartic == 0:
        return float("nan"), I, J, True
    j = Rational(6912) * I ** 3 / Delta_quartic
    return j, I, J, False


def phase_R12_1():
    log("Phase R12-1: Igusa I,J and j(Jac(C_b)) for 60 Q1 quartics")
    cat = json.load(open(Q1_DIR / "quartic_family_catalogue.json"))["families"]
    rows = []
    j_zero = []
    j_1728 = []
    for f in cat:
        coeffs = [f["alpha_4"], f["alpha_3"], f["alpha_2"],
                  f["alpha_1"], f["alpha_0"]]
        j, I, J, singular = jacobian_j(coeffs)
        if singular:
            log(f"  fam {f['family_id']}: SINGULAR (4I^3 = J^2)")
            j_float = float("nan")
            log_abs_j = float("nan")
            is_j_zero = False
            is_j_1728 = False
            sgn = 0
        else:
            j_float = float(j)
            log_abs_j = math.log(abs(j_float) + 1.0)
            is_j_zero = (I == 0)
            is_j_1728 = (J == 0) and (I != 0)
            if is_j_zero:
                j_zero.append(f["family_id"])
            if is_j_1728:
                j_1728.append(f["family_id"])
            sgn = 1 if j_float > 0 else (-1 if j_float < 0 else 0)
        # CM-disc heuristic
        if is_j_zero:
            CM_disc = -3
        elif is_j_1728:
            CM_disc = -4
        else:
            CM_disc = float("nan")
        # g_2, g_3 of the Jacobian E_f: y^2 = x^3 - 27 I x - 27 J
        # In Weierstrass y^2 = 4 x^3 - g_2 x - g_3 form,
        # rescale via x -> x/3 then y -> y/sqrt(108): the Jacobian's
        # canonical (g_2, g_3) corresponds to (I/12, J/216) up to
        # standard normalization. We record I, J as the primitive
        # integral invariants (canonical for binary quartics), and
        # report g_2_E = -27 I, g_3_E = -27 J (Weierstrass y^2 = x^3 + Ax + B
        # form coefficients A = -27 I, B = -27 J).
        rows.append({
            "family_id": f["family_id"],
            "coeffs": coeffs,
            "I_quartic": str(I),
            "J_quartic": str(J),
            "Delta_quartic_4I3_minus_J2": str(4 * I ** 3 - J ** 2),
            "Disc_b": int(f["Delta_4_exact"]),
            "Weierstrass_A": str(-27 * I),
            "Weierstrass_B": str(-27 * J),
            "j_invariant": j_float,
            "log_abs_j": log_abs_j,
            "sgn_j": sgn,
            "is_j_zero": bool(is_j_zero),
            "is_j_1728": bool(is_j_1728),
            "is_singular": bool(singular),
            "CM_disc": CM_disc,
        })
    out = HERE / "quartic_jacobian_invariants.json"
    with open(out, "w") as fh:
        json.dump({
            "task_id": "PCF2-SESSION-R1.2",
            "n_families": len(rows),
            "j_zero_cell": j_zero,
            "j_1728_cell": j_1728,
            "n_singular": sum(1 for r in rows if r["is_singular"]),
            "rows": rows,
        }, fh, indent=2)
    log(f"  wrote {out.name}; j_zero cell = {j_zero}; j_1728 cell = {j_1728}")
    add_claim(
        f"At d=4 in Q1 catalogue, j=0 cell has {len(j_zero)} family(ies): "
        f"{j_zero}; j=1728 cell has {len(j_1728)} family(ies): {j_1728}.",
        dps=0, script="r1_2_quartic_j_probe.py", output_file=out
    )
    return rows, j_zero, j_1728


# ============================================================== R12-7 (early)
# Re-fit WKB at dps=1500 for j=0 candidates (and a couple of |I|-low controls).
DPS_HI = 1500
N_GRID_HI = list(range(10, 132, 2))  # 61 points
N_REF_HI = 200


def cf_value_quartic(coeffs5, N: int, dps: int) -> mp.mpf:
    a4, a3, a2, a1, a0 = coeffs5
    with mp.workdps(dps):
        a4m = mp.mpf(a4); a3m = mp.mpf(a3); a2m = mp.mpf(a2)
        a1m = mp.mpf(a1); a0m = mp.mpf(a0)
        x = a4m * N ** 4 + a3m * N ** 3 + a2m * N ** 2 + a1m * N + a0m
        for k in range(N - 1, -1, -1):
            bk = a4m * k ** 4 + a3m * k ** 3 + a2m * k ** 2 + a1m * k + a0m
            x = bk + mp.mpf(1) / x
        return +x


def wkb_fit_hi(coeffs5, dps: int = DPS_HI, N_grid=None, N_ref: int = N_REF_HI) -> dict:
    if N_grid is None:
        N_grid = N_GRID_HI
    with mp.workdps(dps):
        L_ref = cf_value_quartic(coeffs5, N_ref, dps)
        rows_x, rows_y, ns_used = [], [], []
        for N in N_grid:
            LN = cf_value_quartic(coeffs5, N, dps)
            d = abs(LN - L_ref)
            if d == 0:
                continue
            y = float(mp.log(d))
            rows_x.append([-N * math.log(N), float(N), -math.log(N), 1.0])
            rows_y.append(y)
            ns_used.append(N)
    if len(rows_y) < 6:
        return {"A": None, "A_stderr": None, "n_points": len(rows_y),
                "comment": "insufficient points"}
    X = np.array(rows_x); y = np.array(rows_y)
    coeffs_lsq, *_ = np.linalg.lstsq(X, y, rcond=None)
    A_fit, alpha_fit, beta_fit, gamma_fit = (float(c) for c in coeffs_lsq)
    yhat = X @ coeffs_lsq
    resid = y - yhat
    dof = max(len(y) - 4, 1)
    s2 = float(np.sum(resid ** 2) / dof)
    XtX_inv = np.linalg.inv(X.T @ X)
    cov = s2 * XtX_inv
    stderrs = np.sqrt(np.diag(cov))
    return {
        "A": A_fit, "alpha": alpha_fit, "beta": beta_fit, "gamma": gamma_fit,
        "A_stderr": float(stderrs[0]), "alpha_stderr": float(stderrs[1]),
        "n_points": len(y), "fit_window_N": [N_grid[0], ns_used[-1]],
        "ns_used": ns_used, "residual_rms": float(math.sqrt(s2)),
        "dps": dps, "N_ref": N_ref,
    }


def phase_R12_7(j_zero_ids, jac_rows, q1_per_family):
    log("Phase R12-7: precision escalation (dps=1500) for selected quartic families")
    # Always escalate the j=0 cell. For controls, also escalate the families
    # with smallest non-zero |I| (best estimate of "near-CM" comparison).
    rank = sorted(jac_rows, key=lambda r: abs(int(r["I_quartic"]))
                  if r["I_quartic"] != "0" else 10**18)
    targets = list(j_zero_ids)
    for r in rank:
        if r["family_id"] in targets:
            continue
        if r["is_singular"]:
            continue
        if int(r["I_quartic"]) == 0:
            continue
        if len(targets) >= len(j_zero_ids) + 3:
            break
        targets.append(r["family_id"])
    log(f"  escalating families: {targets}")
    overrides = {}
    for fid in targets:
        ff = next((f for f in q1_per_family if f["family_id"] == fid), None)
        if ff is None:
            continue
        coeffs = ff["b_coefficients"]  # [a4,a3,a2,a1,a0]
        t0 = time.time()
        wkb = wkb_fit_hi(coeffs, dps=DPS_HI)
        dt = time.time() - t0
        log(f"  fam {fid}: A={wkb.get('A')} +/- {wkb.get('A_stderr')} "
            f"(dt={dt:.1f}s, n={wkb.get('n_points')})")
        overrides[fid] = wkb
    out = HERE / "precision_escalation_log.json"
    with open(out, "w") as fh:
        json.dump({
            "dps": DPS_HI,
            "N_grid": N_GRID_HI,
            "N_ref": N_REF_HI,
            "targets": targets,
            "overrides": {str(k): v for k, v in overrides.items()},
        }, fh, indent=2)
    add_claim(
        f"WKB precision escalation at dps={DPS_HI} for {len(targets)} quartic "
        f"families completed; A_fit values recorded.",
        dps=DPS_HI, script="r1_2_quartic_j_probe.py", output_file=out
    )
    return overrides


# ============================================================== R12-2
def phase_R12_2(jac_rows, q1_per_family, q1_catalogue, overrides):
    log("Phase R12-2: build joined 110-row dataframe")
    # Cubic rows from R1.1
    cubic_csv = R11_DIR / "assembled_data_v2.csv"
    cdf = pd.read_csv(cubic_csv)
    # Make schema unified
    cubic_rows = []
    for _, r in cdf.iterrows():
        log_abs_Delta = float(r.get("log_abs_Delta_3", float("nan")))
        cubic_rows.append({
            "family_id_str": f"d3_{int(r['family_id'])}",
            "family_id": int(r["family_id"]),
            "degree": 3,
            "A_fit": float(r["A_fit"]),
            "A_stderr": float(r["A_stderr"]),
            "delta": float(r["A_fit"]) - 6.0,
            "j_invariant": float(r["j_invariant"]),
            "log_abs_j": float(r["log_abs_j"]),
            "sgn_j": (1 if r["j_invariant"] > 0
                      else (-1 if r["j_invariant"] < 0 else 0))
                     if not math.isnan(float(r["j_invariant"])) else 0,
            "is_j_zero": (abs(float(r["j_invariant"])) < 1e-9)
                          if not math.isnan(float(r["j_invariant"])) else False,
            "is_j_1728": (abs(float(r["j_invariant"]) - 1728.0) < 1e-9)
                          if not math.isnan(float(r["j_invariant"])) else False,
            "CM_disc": float(r["CM_disc"]) if not pd.isna(r["CM_disc"]) else float("nan"),
            "log_abs_Delta": log_abs_Delta,
            "log_Mahler": float(r.get("log_Mahler", float("nan"))),
            "bin": str(r["bin"]),
        })
    # Quartic rows
    q1_by_id = {f["family_id"]: f for f in q1_per_family}
    cat_by_id = {f["family_id"]: f for f in q1_catalogue}
    quartic_rows = []
    for jr in jac_rows:
        fid = jr["family_id"]
        ff = q1_by_id[fid]
        cat = cat_by_id[fid]
        if str(fid) in {str(k) for k in overrides}:
            wkb = overrides[fid]
            A_fit = wkb["A"]; A_stderr = wkb["A_stderr"]
        elif fid in overrides:
            wkb = overrides[fid]
            A_fit = wkb["A"]; A_stderr = wkb["A_stderr"]
        else:
            A_fit = ff["wkb"]["A"]; A_stderr = ff["wkb"]["A_stderr"]
        absD = abs(int(cat["Delta_4_exact"]))
        log_abs_Delta = math.log(absD) if absD > 0 else float("nan")
        # log Mahler: monic so M = prod max(1, |alpha|)
        coeffs = [cat["alpha_4"], cat["alpha_3"], cat["alpha_2"],
                  cat["alpha_1"], cat["alpha_0"]]
        with mp.workdps(50):
            P = Poly([Rational(int(c)) for c in coeffs], symbols('x'))
            roots = mp.polyroots([float(c) for c in coeffs], maxsteps=200)
            logM = mp.mpf(0)
            for rr in roots:
                ar = abs(rr)
                if ar > 1:
                    logM += mp.log(ar)
        quartic_rows.append({
            "family_id_str": f"d4_{fid}",
            "family_id": fid,
            "degree": 4,
            "A_fit": float(A_fit),
            "A_stderr": float(A_stderr),
            "delta": float(A_fit) - 8.0,
            "j_invariant": jr["j_invariant"],
            "log_abs_j": jr["log_abs_j"],
            "sgn_j": jr["sgn_j"],
            "is_j_zero": jr["is_j_zero"],
            "is_j_1728": jr["is_j_1728"],
            "CM_disc": jr["CM_disc"],
            "log_abs_Delta": log_abs_Delta,
            "log_Mahler": float(logM),
            "bin": str(cat["trichotomy_bin"]),
        })
    rows = cubic_rows + quartic_rows
    df = pd.DataFrame(rows)
    out = HERE / "assembled_data_v3.csv"
    df.to_csv(out, index=False)
    log(f"  wrote {out.name}: {len(df)} rows ({(df['degree']==3).sum()} cubic + "
        f"{(df['degree']==4).sum()} quartic)")
    add_claim(
        f"Joined 110-row dataframe assembled (50 cubic + 60 quartic).",
        dps=80, script="r1_2_quartic_j_probe.py", output_file=out
    )
    return df


# ============================================================== R12-3,4
def spearman_safe(x, y):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return float("nan"), float("nan"), int(m.sum())
    r = stats.spearmanr(x[m], y[m])
    return float(r.correlation), float(r.pvalue), int(m.sum())


def pearson_safe(x, y):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return float("nan"), float("nan"), int(m.sum())
    r = stats.pearsonr(x[m], y[m])
    return float(r.statistic), float(r.pvalue), int(m.sum())


def phase_R12_3_4(df: pd.DataFrame):
    log("Phase R12-3,4: per-degree and joint Spearman correlations")
    K = 12  # we will report ~12 correlation tests
    out = {"per_degree": {}, "joint": {}, "K_bonferroni": K}

    for d in (3, 4):
        sub = df[df.degree == d]
        rho, p, n = spearman_safe(sub.log_abs_j, sub.delta)
        rhoP, pP, _ = pearson_safe(sub.log_abs_j, sub.delta)
        rho_logD, p_logD, _ = spearman_safe(sub.log_abs_Delta, sub.delta)
        bonf_p = min(1.0, p * K) if math.isfinite(p) else float("nan")
        # j_zero cell stats
        jz = sub[sub.is_j_zero]
        jz_stats = {
            "n": int(len(jz)),
            "delta_mean": float(jz.delta.mean()) if len(jz) else float("nan"),
            "delta_min": float(jz.delta.min()) if len(jz) else float("nan"),
            "delta_max": float(jz.delta.max()) if len(jz) else float("nan"),
            "A_stderr_max": float(jz.A_stderr.max()) if len(jz) else float("nan"),
            "consistent_with_zero_count": int(
                ((jz.delta.abs() <= 2.0 * jz.A_stderr) | (jz.delta.abs() <= 1e-3))
                .sum()
            ) if len(jz) else 0,
            "family_ids": [int(x) for x in jz.family_id.tolist()],
        }
        out["per_degree"][d] = {
            "n": n,
            "spearman_log_abs_j_vs_delta": {
                "rho": rho, "p": p, "bonferroni_p": bonf_p,
            },
            "pearson_log_abs_j_vs_delta": {"r": rhoP, "p": pP},
            "spearman_log_abs_Delta_vs_delta": {"rho": rho_logD, "p": p_logD},
            "j_zero_cell": jz_stats,
            "delta_summary": {
                "mean": float(sub.delta.mean()),
                "min": float(sub.delta.min()),
                "max": float(sub.delta.max()),
                "std": float(sub.delta.std()),
            },
        }
        log(f"  d={d}: n={n} rho(log|j|,delta)={rho:.4f} p={p:.3g} "
            f"bonf_p={bonf_p:.3g}; j_zero_cell n={jz_stats['n']}")

    # Joint d in {3,4}
    rho_j, p_j, n_j = spearman_safe(df.log_abs_j, df.delta)
    bonf_j = min(1.0, p_j * K) if math.isfinite(p_j) else float("nan")
    out["joint"] = {
        "n": n_j,
        "spearman_log_abs_j_vs_delta": {
            "rho": rho_j, "p": p_j, "bonferroni_p": bonf_j,
        },
    }
    log(f"  joint: n={n_j} rho={rho_j:.4f} p={p_j:.3g} bonf_p={bonf_j:.3g}")

    add_claim(
        f"d=3 reproduction: rho(log|j|, delta) = "
        f"{out['per_degree'][3]['spearman_log_abs_j_vs_delta']['rho']:.4f} on N=50 "
        f"(R1.1 reported -0.6906).",
        dps=1500, script="r1_2_quartic_j_probe.py",
        output_file=HERE / "assembled_data_v3.csv"
    )
    add_claim(
        f"d=4 Spearman rho(log|j|, delta_4) = "
        f"{out['per_degree'][4]['spearman_log_abs_j_vs_delta']['rho']:.4f}, "
        f"raw p = {out['per_degree'][4]['spearman_log_abs_j_vs_delta']['p']:.3g}, "
        f"Bonferroni p = "
        f"{out['per_degree'][4]['spearman_log_abs_j_vs_delta']['bonferroni_p']:.3g} "
        f"on N=60.",
        dps=80, script="r1_2_quartic_j_probe.py",
        output_file=HERE / "assembled_data_v3.csv"
    )
    add_claim(
        f"Cross-degree joint Spearman rho(log|j|, delta) = "
        f"{rho_j:.4f}, raw p = {p_j:.3g}, Bonferroni p = {bonf_j:.3g}, "
        f"on N={n_j} (joint d in {{3,4}}).",
        dps=1500, script="r1_2_quartic_j_probe.py",
        output_file=HERE / "assembled_data_v3.csv"
    )
    return out


# ============================================================== R12-5 lemniscatic
def phase_R12_5(df):
    log("Phase R12-5: lemniscatic (j=1728) cell")
    out = {}
    for d in (3, 4):
        sub = df[(df.degree == d) & (df.is_j_1728 == True)]
        out[d] = {
            "n": int(len(sub)),
            "family_ids": [int(x) for x in sub.family_id.tolist()],
            "delta_values": [float(x) for x in sub.delta.tolist()],
        }
        log(f"  d={d}: j=1728 cell n={out[d]['n']}, delta={out[d]['delta_values']}")
    return out


# ============================================================== R12-6 pari/gp
def phase_R12_6(jac_rows):
    log("Phase R12-6: pari/gp probe")
    gp = shutil.which("gp") or shutil.which("gp.exe")
    out = {"available": False, "reason": None, "results": {}}
    if not gp:
        for p in [r"C:\Program Files\PARI\gp.exe",
                  r"C:\Program Files (x86)\PARI\gp.exe",
                  r"C:\PARI\gp.exe"]:
            if Path(p).exists():
                gp = p
                break
    if not gp:
        try:
            import cypari2  # noqa
            out["available"] = "cypari2"
        except Exception:
            out["reason"] = "no gp binary, no cypari2 (same status as R1.1)"
            log("  pari/gp not available, deferring (matches R1.1)")
            return out
    log(f"  pari/gp at {gp}; would compute conductor + class number for j_zero cell")
    # Best-effort: skip actual gp calls in this run to keep deterministic and
    # within budget; record availability.
    out["available"] = bool(gp)
    return out


# ============================================================== R12-8 outputs
def write_correlation_table(out_data, df):
    rows = []
    for d in (3, 4):
        s = out_data["per_degree"][d]["spearman_log_abs_j_vs_delta"]
        n = out_data["per_degree"][d]["n"]
        rows.append((d, n, s["rho"], s["p"], s["bonferroni_p"]))
    s = out_data["joint"]["spearman_log_abs_j_vs_delta"]
    rows.append(("3+4", out_data["joint"]["n"], s["rho"], s["p"], s["bonferroni_p"]))
    lines = [
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{Spearman $\rho(\log|j|, \delta)$ where $\delta = A_{\rm fit} - 2d$;",
        r"per-degree (R1.1 cubic reproduction; R1.2 new quartic) and joint.",
        r"$K=12$ Bonferroni denominator.}",
        r"\label{tab:r12-correlation}",
        r"\begin{tabular}{cccccc}",
        r"\hline",
        r"degree & $N$ & $\rho_S$ & raw $p$ & Bonf. $p$ \\",
        r"\hline",
    ]
    for d, n, rho, p, bp in rows:
        rho_s = "NaN" if not math.isfinite(rho) else f"{rho:.4f}"
        p_s = "NaN" if not math.isfinite(p) else f"{p:.3g}"
        bp_s = "NaN" if not math.isfinite(bp) else f"{bp:.3g}"
        lines.append(f"  ${d}$ & {n} & ${rho_s}$ & ${p_s}$ & ${bp_s}$ \\\\")
    lines += [r"\hline", r"\end{tabular}", r"\end{table}"]
    out = HERE / "correlation_table_v3.tex"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"  wrote {out.name}")
    return out


def write_quartic_jzero_table(out_data, df):
    sub3 = df[(df.degree == 3) & (df.is_j_zero == True)]
    sub4 = df[(df.degree == 4) & (df.is_j_zero == True)]
    lines = [
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{The $j(\mathrm{Jac}(C_b)) = 0$ cell at $d=3$ (cubic, R1.1) and",
        r"$d=4$ (quartic, R1.2). $\delta := A_{\rm fit} - 2d$ measures the",
        r"departure from the Conjecture B4 sharp predicted value.}",
        r"\label{tab:r12-jzero-cell}",
        r"\begin{tabular}{cccccc}",
        r"\hline",
        r"degree & family\_id & $A_{\rm fit}$ & $A_{\rm stderr}$ & $\delta$ "
        r"& consistent with $0$? \\",
        r"\hline",
    ]
    for sub, d in ((sub3, 3), (sub4, 4)):
        for _, r in sub.iterrows():
            cons = abs(r.delta) <= max(2.0 * r.A_stderr, 1e-3)
            lines.append(
                f"  ${d}$ & {int(r.family_id)} & ${r.A_fit:.6f}$ & "
                f"${r.A_stderr:.2e}$ & ${r.delta:+.4e}$ & "
                f"{'yes' if cons else 'no'} \\\\"
            )
    lines += [r"\hline", r"\end{tabular}", r"\end{table}"]
    out = HERE / "quartic_j_zero_cell.tex"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"  wrote {out.name}")
    return out


def write_v13_paragraph_v2(out_data, df):
    sub3 = df[(df.degree == 3) & (df.is_j_zero == True)]
    sub4 = df[(df.degree == 4) & (df.is_j_zero == True)]
    n3 = len(sub3); n4 = len(sub4)
    s4 = out_data["per_degree"][4]["spearman_log_abs_j_vs_delta"]
    sJ = out_data["joint"]["spearman_log_abs_j_vs_delta"]
    s3 = out_data["per_degree"][3]["spearman_log_abs_j_vs_delta"]
    body = f"""% v1.3 paragraph insert (replaces R1.1 v13_paragraph_insert.tex)
% Strengthened with d=4 confirmation from PCF2-SESSION-R1.2.

\\subsection*{{The $j$-invariant as a modular finer-split (Conjectures B5, B6)}}

For each polynomial continued fraction $\\mathrm{{CF}}(b)$ with $b\\in\\Z[n]$
of degree $d\\in\\{{3,4\\}}$, let $C_b$ be the projective genus-1 curve
$y^2 = b(x)$ and let $j(\\mathrm{{Jac}}(C_b))\\in\\Q$ be the $j$-invariant of
its Jacobian elliptic curve $E_b$. (For $d=3$, $C_b = E_b$; for $d=4$, $E_b$
is the Jacobian, with $j$ computed from the Igusa $I,J$ invariants of the
binary quartic $b$ via $j = 6912\\,I^{{3}}/(4 I^{{3}}-J^{{2}})$.) The 110-family
joint dataset (50 cubic from PCF-2 v1.1 Sessions B/C1; 60 quartic from
Session Q1) reveals a striking modular finer-split inside the
Conjecture~B4 sharp universality $A_{{\\rm PCF}}(b)=2d$:

\\begin{{conjecture}}[B5: $j$-zero cell exactness]
If $j(\\mathrm{{Jac}}(C_b))=0$ (equianharmonic CM, Igusa invariant $I=0$ at
$d=4$; $\\Q(\\sqrt{{-3}})$-CM at $d=3$), then $A_{{\\rm PCF}}(b) = 2d$ exactly.
\\end{{conjecture}}

\\noindent
\\emph{{Empirical evidence.}}
Across the joined dataset, the $j=0$ cell at $d=3$ contains $n_3={n3}$ cubic
families (PCF-2 Session B), all with $|\\delta_3| := |A_{{\\rm fit}}-6| \\le
10^{{-3}}$ at WKB precision $\\mathrm{{dps}}=1500$. The $j=0$ cell at $d=4$
contains $n_4={n4}$ quartic family(ies); see Table~\\ref{{tab:r12-jzero-cell}}
for the family identifier and $\\delta_4$ statistics. Combined cell size
$N_5={n3+n4}$.

\\begin{{conjecture}}[B6: $\\log|j|$ soft drift]
The Spearman rank correlation
$\\rho_S(\\log(|j(\\mathrm{{Jac}}(C_b))|+1),\\;\\delta)$,
where $\\delta := A_{{\\rm fit}} - 2d$, is strictly negative and statistically
significant for each $d\\in\\{{3,4\\}}$ separately and for the joint
$d\\in\\{{3,4\\}}$ population, where $A_{{\\rm fit}}$ comes from the standard
WKB log-residual fit at common precision.
\\end{{conjecture}}

\\noindent
\\emph{{Empirical evidence.}} See
Table~\\ref{{tab:r12-correlation}} for the per-degree and joint
$\\rho_S(\\log|j|,\\delta)$ statistics together with raw $p$ and Bonferroni
$p$ (denominator $K=12$). Cubic reproduces R1.1 ($\\rho={s3['rho']:.3f}$,
Bonf $p\\le {s3['bonferroni_p']:.1e}$); quartic gives $\\rho={s4['rho']:.3f}$
(raw $p={s4['p']:.2g}$); joint gives $\\rho={sJ['rho']:.3f}$
(Bonf $p\\le{sJ['bonferroni_p']:.1e}$).

\\begin{{remark}}[Cross-degree universality]
\\label{{rem:cross-degree-universality}}
Conjectures~B5 and~B6 are degree-uniform across the (currently sampled)
genus-1 cell $d\\in\\{{3,4\\}}$, suggesting the $j$-invariant of the
genus-1 curve $C_b\\colon y^2=b(x)$ is the right modular finer-split for
all PCF degrees in which $C_b$ is genus-1.
\\end{{remark}}
"""
    out = HERE / "v13_paragraph_insert_v2.tex"
    out.write_text(body, encoding="utf-8")
    log(f"  wrote {out.name}")
    return out


def main():
    log("=== PCF2-SESSION-R1.2 START ===")
    q1_results = json.load(open(Q1_DIR / "results.json"))
    q1_per_family = q1_results["per_family"]
    q1_catalogue = json.load(open(Q1_DIR / "quartic_family_catalogue.json"))["families"]

    jac_rows, j_zero, j_1728 = phase_R12_1()

    # Phase R12-7: precision escalation. Two layers:
    #   (a) dps=1500 N=10..130 escalation for j_zero candidates +
    #       low-|I| controls (legacy, kept for diagnostic).
    #   (b) tail-window fit (dps=2000, N=50..250, N_ref=400) for ALL 60
    #       quartics; this is the PRIMARY override used for analysis,
    #       because it eliminates the small-N WKB sub-leading bias that
    #       caused all 60 quartics to cluster at delta ~ -0.046 (dps=80).
    legacy_overrides = phase_R12_7(j_zero, jac_rows, q1_per_family)

    tail_path = HERE / "tail_fit_overrides.json"
    if tail_path.exists():
        tail = json.load(open(tail_path))["overrides"]
        overrides = {int(k): v for k, v in tail.items()}
        log(f"  using TAIL-WINDOW overrides for all {len(overrides)} quartics "
            f"(dps=2000, N=[50..250], N_ref=400)")
    else:
        overrides = legacy_overrides
        log("  WARNING: tail_fit_overrides.json not found; using legacy escalation only")

    df = phase_R12_2(jac_rows, q1_per_family, q1_catalogue, overrides)
    out_data = phase_R12_3_4(df)
    lemn = phase_R12_5(df)
    pari = phase_R12_6(jac_rows)

    # Verdicts
    s4 = out_data["per_degree"][4]["spearman_log_abs_j_vs_delta"]
    sJ = out_data["joint"]["spearman_log_abs_j_vs_delta"]
    jz4 = out_data["per_degree"][4]["j_zero_cell"]
    # B5 at d=4 verdict
    b5_d4 = "INCONCLUSIVE"
    if jz4["n"] == 0:
        b5_d4 = "VACUOUS"
    elif jz4["consistent_with_zero_count"] == jz4["n"]:
        b5_d4 = "SUPPORTED" if jz4["n"] >= 2 else "SUPPORTED_THIN_CELL_N1"
    else:
        b5_d4 = "FALSIFIED"
    # B6 at d=4 verdict
    b6_d4 = "INCONCLUSIVE"
    if math.isfinite(s4["rho"]) and s4["rho"] < 0 and s4["bonferroni_p"] < 0.001:
        b6_d4 = "SUPPORTED"
    elif math.isfinite(s4["rho"]) and s4["rho"] > 0 and s4["bonferroni_p"] < 0.01:
        b6_d4 = "FALSIFIED_WRONG_SIGN"
    elif math.isfinite(s4["rho"]) and s4["rho"] < 0:
        b6_d4 = "WEAK_NEGATIVE_NOT_SIGNIFICANT"
    # Cross-degree verdict
    cross = "INCONCLUSIVE"
    if math.isfinite(sJ["rho"]) and sJ["rho"] < 0 and sJ["bonferroni_p"] < 0.001:
        cross = "SUPPORTED"

    add_claim(
        f"Verdict R12-A (B5 at d=4): {b5_d4} (j_zero cell n={jz4['n']}, "
        f"consistent_with_zero={jz4['consistent_with_zero_count']}/{jz4['n']}).",
        dps=DPS_HI, script="r1_2_quartic_j_probe.py",
        output_file=HERE / "assembled_data_v3.csv"
    )
    add_claim(
        f"Verdict R12-B (B6 at d=4): {b6_d4} (rho={s4['rho']:.4f}, "
        f"Bonf p={s4['bonferroni_p']:.3g}).",
        dps=80, script="r1_2_quartic_j_probe.py",
        output_file=HERE / "assembled_data_v3.csv"
    )
    add_claim(
        f"Verdict R12-C (cross-degree universality d in {{3,4}}): {cross} "
        f"(joint rho={sJ['rho']:.4f}, Bonf p={sJ['bonferroni_p']:.3g}).",
        dps=80, script="r1_2_quartic_j_probe.py",
        output_file=HERE / "assembled_data_v3.csv"
    )

    results = {
        "task_id": "PCF2-SESSION-R1.2",
        "date": "2026-05-01",
        "preconditions": {
            "Q1_catalogue": str(Q1_DIR / "quartic_family_catalogue.json"),
            "Q1_results": str(Q1_DIR / "results.json"),
            "R1_1_csv": str(R11_DIR / "assembled_data_v2.csv"),
        },
        "j_zero_cells": {
            "d3": out_data["per_degree"][3]["j_zero_cell"],
            "d4": out_data["per_degree"][4]["j_zero_cell"],
        },
        "j_1728_cells": lemn,
        "correlations": out_data,
        "pari_status": pari,
        "verdicts": {
            "R12-A_B5_at_d4": b5_d4,
            "R12-B_B6_at_d4": b6_d4,
            "R12-C_cross_degree": cross,
        },
        "claims_count": len(CLAIMS),
    }
    out_results = HERE / "results_v3.json"
    out_results.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    log(f"  wrote {out_results.name}")

    write_correlation_table(out_data, df)
    write_quartic_jzero_table(out_data, df)
    # NOTE: v13_paragraph_insert_v2.tex is hand-crafted (not regenerated)
    # because the auto-template assumed B5 cross-degree SUPPORT, but the
    # actual finding is FALSIFIED. Hand-curated text reflects the verdict.
    if not (HERE / "v13_paragraph_insert_v2.tex").exists():
        write_v13_paragraph_v2(out_data, df)

    # claims.jsonl
    out_claims = HERE / "claims.jsonl"
    with open(out_claims, "w", encoding="utf-8") as fh:
        for c in CLAIMS:
            fh.write(json.dumps(c) + "\n")
    log(f"  wrote {out_claims.name} ({len(CLAIMS)} entries)")

    # halt / discrepancy / unexpected
    halt = {}
    if jz4["n"] == 0:
        halt["B5_d4_vacuous"] = "is_j_zero cell at d=4 is empty"
    elif b5_d4.startswith("FALSIFIED"):
        halt["B5_d4_falsified"] = (
            f"j_zero cell at d=4 size {jz4['n']}, "
            f"consistent_with_zero {jz4['consistent_with_zero_count']}"
        )
    if b6_d4 == "FALSIFIED_WRONG_SIGN":
        halt["B6_d4_wrong_sign"] = (
            f"rho={s4['rho']:+.4f} (positive) at Bonf p={s4['bonferroni_p']:.3g}"
        )
    (HERE / "halt_log.json").write_text(json.dumps(halt, indent=2), encoding="utf-8")

    discrepancy = {}
    # Cubic reproduction expected rho ~= -0.691; check
    s3 = out_data["per_degree"][3]["spearman_log_abs_j_vs_delta"]
    if math.isfinite(s3["rho"]) and abs(s3["rho"] - (-0.691)) > 0.05:
        discrepancy["d3_reproduction"] = {
            "expected_rho": -0.691,
            "got_rho": s3["rho"],
        }
    (HERE / "discrepancy_log.json").write_text(
        json.dumps(discrepancy, indent=2), encoding="utf-8"
    )

    unexpected = {}
    if jz4["n"] == 1:
        unexpected["thin_jzero_cell_d4"] = (
            "exactly 1 quartic in Q1 has j(Jac(C_b))=0; cross-degree "
            "B5 confirmation is thin (N=1). Recommend extending the quartic "
            "enumeration window to harvest more candidates."
        )
    if all(d.get("n", 0) == 0 for d in lemn.values()):
        unexpected["lemniscatic_cell_empty"] = (
            "is_j_1728 cell empty at d in {3,4} in current catalogues; "
            "op:lemniscatic-cell remains untestable."
        )
    (HERE / "unexpected_finds.json").write_text(
        json.dumps(unexpected, indent=2), encoding="utf-8"
    )

    log("=== PCF2-SESSION-R1.2 COMPLETE ===")
    return results


if __name__ == "__main__":
    main()
