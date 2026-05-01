"""
PCF2-SESSION-R1 -- finer-cubic-split probe.

Re-analyse the 50 cubic A_fit datapoints from PCF2 Sessions B and C1 for
catalogue-level invariants that correlate with delta = A_fit - 6.

Pipeline phases R1-1 .. R1-5 per relay prompt.

Outputs (alongside this script):
  assembled_data.csv
  results.json
  correlation_table.tex
  claims.jsonl
  rubber_duck_critique.md
  handoff.md
  halt_log.json, discrepancy_log.json, unexpected_finds.json
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import sys
import shutil
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
from sympy import Integer, Rational, factorint, isprime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
SESSION_DIR = os.path.dirname(os.path.abspath(__file__))
CATALOGUE = os.path.join(ROOT, "sessions", "2026-05-01", "PCF2-SESSION-A",
                         "cubic_family_catalogue.json")
SESSION_B = os.path.join(ROOT, "sessions", "2026-05-01", "PCF2-SESSION-B",
                         "results.json")
SESSION_C1 = os.path.join(ROOT, "sessions", "2026-05-01", "PCF2-SESSION-C1",
                          "results.json")

# --- Phase R1-1: data assembly --------------------------------------------- #


def parse_cm_disc(cm_field_latex: str) -> int | None:
    """Extract d from \\mathbb{Q}(\\sqrt{d}) (signed integer)."""
    if not cm_field_latex:
        return None
    m = re.search(r"sqrt\{(-?\d+)\}", cm_field_latex)
    if m:
        return int(m.group(1))
    return None


def assemble_dataframe() -> pd.DataFrame:
    cat = json.load(open(CATALOGUE))["families"]
    catB = json.load(open(SESSION_B))["per_family"]
    catC = json.load(open(SESSION_C1))["per_family"]

    by_id_cat = {f["family_id"]: f for f in cat}
    fits = {}
    for f in catB:
        fits[f["family_id"]] = f
    for f in catC:
        fits[f["family_id"]] = f

    rows = []
    for fid, c in by_id_cat.items():
        if fid not in fits:
            continue  # only families with a WKB fit
        ff = fits[fid]
        wkb = ff.get("wkb", {}) or {}
        A_fit = wkb.get("A")
        A_stderr = wkb.get("A_stderr")
        alpha_fit = wkb.get("alpha")
        alpha_stderr = wkb.get("alpha_stderr")
        if A_fit is None:
            continue
        Delta = int(c["Delta_3_exact"])
        absD = abs(Delta)
        Dfact = factorint(absD)
        prime_divs = sorted(int(p) for p in Dfact.keys())
        omega = len(prime_divs)
        smallest_p = prime_divs[0] if prime_divs else 0
        # Embedding signature for cubic: Delta>0 -> 3 real / 0 cplx pairs;
        # Delta<0 -> 1 real / 1 complex pair.
        if Delta > 0:
            r1, r2 = 3, 0
        else:
            r1, r2 = 1, 1
        cm_d = parse_cm_disc(c.get("CM_field", ""))
        rows.append(dict(
            family_id=fid,
            alpha_3=int(c["alpha_3"]),
            alpha_2=int(c["alpha_2"]),
            alpha_1=int(c["alpha_1"]),
            alpha_0=int(c["alpha_0"]),
            Delta_3=Delta,
            Delta_3_sign=1 if Delta > 0 else (-1 if Delta < 0 else 0),
            abs_Delta_3=absD,
            log_abs_Delta_3=math.log(absD),
            absD_mod_8=absD % 8,
            absD_mod_12=absD % 12,
            absD_parity=absD % 2,
            smallest_prime_div=smallest_p,
            omega_Delta_3=omega,
            Galois_group=str(c["Galois_group"]),
            splitting_field=str(c.get("splitting_field", "")),
            CM_field=str(c.get("CM_field", "")),
            CM_disc=cm_d if cm_d is not None else float("nan"),
            r1=r1,
            r2=r2,
            bin=str(c["trichotomy_bin"]),
            # Class number / regulator / fund-unit-norm: pari/gp not available
            # in this environment (cypari2 not installed, gp.exe not on PATH).
            # Recorded as NaN; flagged in unexpected_finds.json.
            conductor=float("nan"),
            class_number=float("nan"),
            class_number_plus=float("nan"),
            regulator=float("nan"),
            fund_unit_norm=float("nan"),
            unit_density=float("nan"),
            A_fit=float(A_fit),
            A_stderr=float(A_stderr) if A_stderr is not None else float("nan"),
            alpha_fit=float(alpha_fit) if alpha_fit is not None else float("nan"),
            alpha_stderr=float(alpha_stderr) if alpha_stderr is not None else float("nan"),
            delta=float(A_fit) - 6.0,
        ))
    df = pd.DataFrame(rows).sort_values("family_id").reset_index(drop=True)
    return df


# --- Phase R1-2..R1-5 helpers --------------------------------------------- #


CANDIDATE_INVARIANTS = [
    # (column, kind: "ordinal"/"categorical")
    ("alpha_3", "ordinal"),
    ("Delta_3_sign", "ordinal"),
    ("absD_parity", "ordinal"),
    ("absD_mod_8", "ordinal"),
    ("absD_mod_12", "ordinal"),
    ("log_abs_Delta_3", "ordinal"),
    ("smallest_prime_div", "ordinal"),
    ("omega_Delta_3", "ordinal"),
    ("Galois_group", "categorical"),
    ("r1", "ordinal"),  # # real embeddings
    ("r2", "ordinal"),  # # complex pairs
    ("CM_disc", "ordinal"),  # NaN for non-CM bins -> only over CM rows
    ("conductor", "ordinal"),
    ("class_number", "ordinal"),
    ("class_number_plus", "ordinal"),
    ("regulator", "ordinal"),
    ("fund_unit_norm", "ordinal"),
    ("unit_density", "ordinal"),
]


def _spearman(x: np.ndarray, y: np.ndarray):
    if len(x) < 4:
        return float("nan"), float("nan")
    if np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan"), float("nan")
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def _pearson(x: np.ndarray, y: np.ndarray):
    if len(x) < 4 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan"), float("nan")
    r, p = stats.pearsonr(x, y)
    return float(r), float(p)


def _weighted_spearman(x: np.ndarray, y: np.ndarray, w: np.ndarray):
    """Weighted Spearman: rank-correlate then compute weighted Pearson on ranks."""
    if len(x) < 4 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan"), float("nan")
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    w = w / w.sum()
    mx = np.sum(w * rx)
    my = np.sum(w * ry)
    cov = np.sum(w * (rx - mx) * (ry - my))
    vx = np.sum(w * (rx - mx) ** 2)
    vy = np.sum(w * (ry - my) ** 2)
    if vx <= 0 or vy <= 0:
        return float("nan"), float("nan")
    rho = cov / math.sqrt(vx * vy)
    # Approximate p via Fisher z (treat effective n = sum(w)^2 / sum(w^2))
    n_eff = (w.sum() ** 2) / np.sum(w ** 2)
    if n_eff <= 3:
        return float(rho), float("nan")
    z = math.atanh(max(min(rho, 0.9999999), -0.9999999)) * math.sqrt((n_eff - 3) / 1.06)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return float(rho), float(p)


def _kruskal_categorical(cat: np.ndarray, y: np.ndarray):
    """Kruskal-Wallis as a categorical-vs-continuous test; returns (eta_proxy, p)."""
    groups = defaultdict(list)
    for c, v in zip(cat, y):
        groups[c].append(v)
    if len(groups) < 2 or any(len(g) < 2 for g in groups.values()):
        return float("nan"), float("nan")
    H, p = stats.kruskal(*groups.values())
    n = len(y)
    eta = H / (n - 1)  # crude effect-size proxy in [0,1)
    return float(eta), float(p)


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    """Return BH-adjusted q-values (NaN-safe)."""
    arr = np.array(pvals, dtype=float)
    mask = ~np.isnan(arr)
    n = mask.sum()
    out = np.full_like(arr, np.nan)
    if n == 0:
        return out.tolist()
    idx = np.where(mask)[0]
    sub = arr[mask]
    order = np.argsort(sub)
    ranked = sub[order]
    q = np.empty_like(ranked)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        val = ranked[i] * n / (i + 1)
        prev = min(prev, val)
        q[i] = prev
    q_unsorted = np.empty_like(ranked)
    q_unsorted[order] = q
    out[idx] = q_unsorted
    return out.tolist()


def run_correlations(df: pd.DataFrame, label: str, weighted: bool):
    """Apply Phase R1-2/R1-4 single-invariant battery to df, return list of dicts."""
    delta = df["delta"].to_numpy(dtype=float)
    if weighted:
        std = df["A_stderr"].to_numpy(dtype=float)
        std = np.where(np.isnan(std) | (std <= 0), np.nanmedian(std), std)
        weights = 1.0 / (std ** 2)
    else:
        weights = None

    out = []
    for col, kind in CANDIDATE_INVARIANTS:
        if col not in df.columns:
            continue
        x_raw = df[col]
        if kind == "categorical":
            x_arr = x_raw.astype(str).to_numpy()
            mask = x_arr != "nan"
            if mask.sum() < 4:
                eta = float("nan"); p_kw = float("nan")
            else:
                eta, p_kw = _kruskal_categorical(x_arr[mask], delta[mask])
            out.append(dict(invariant=col, kind="categorical", n=int(mask.sum()),
                             rho_or_eta=eta, p=p_kw, pearson_r=float("nan"),
                             pearson_p=float("nan"), label=label))
        else:
            x_num = pd.to_numeric(x_raw, errors="coerce").to_numpy(dtype=float)
            mask = ~np.isnan(x_num) & ~np.isnan(delta)
            if mask.sum() < 4:
                rho = p = pr = pp = float("nan")
            else:
                if weighted:
                    rho, p = _weighted_spearman(x_num[mask], delta[mask], weights[mask])
                else:
                    rho, p = _spearman(x_num[mask], delta[mask])
                pr, pp = _pearson(x_num[mask], delta[mask])
            out.append(dict(invariant=col, kind="ordinal", n=int(mask.sum()),
                             rho_or_eta=rho, p=p, pearson_r=pr, pearson_p=pp,
                             label=label))
    # Bonferroni and BH adjustment over the *valid* (non-NaN p) tests.
    valid_idx = [i for i, r in enumerate(out) if not math.isnan(r["p"])]
    m = len(valid_idx)
    pvals = [out[i]["p"] for i in valid_idx]
    bh = benjamini_hochberg(pvals)
    for k, i in enumerate(valid_idx):
        out[i]["bonferroni_p"] = min(out[i]["p"] * m, 1.0)
        out[i]["bh_q"] = bh[k]
    for r in out:
        r.setdefault("bonferroni_p", float("nan"))
        r.setdefault("bh_q", float("nan"))
        r["m_tests"] = m
    return out


# --- Phase R1-3: cheap multi-invariant search ----------------------------- #


def best_pair_regression(df: pd.DataFrame):
    """Best 2-invariant linear regression (rank-transformed) for delta."""
    y = stats.rankdata(df["delta"].to_numpy(dtype=float))
    candidates = []
    for col, kind in CANDIDATE_INVARIANTS:
        if kind == "categorical":
            continue
        x = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        if np.isnan(x).all():
            continue
        # Replace NaN with column median (only for pair-search descriptive purposes)
        med = np.nanmedian(x) if not np.isnan(x).all() else 0.0
        x = np.where(np.isnan(x), med, x)
        if np.all(x == x[0]):
            continue
        candidates.append((col, stats.rankdata(x)))

    best = []
    for (c1, x1), (c2, x2) in combinations(candidates, 2):
        X = np.column_stack([x1, x2, np.ones_like(x1)])
        try:
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        except np.linalg.LinAlgError:
            continue
        yhat = X @ beta
        ss_res = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        if ss_tot <= 0:
            continue
        r2 = 1 - ss_res / ss_tot
        n = len(y); k = 2
        adj = 1 - (1 - r2) * (n - 1) / max(n - k - 1, 1)
        best.append((c1, c2, float(r2), float(adj)))
    best.sort(key=lambda t: -t[3])
    return best[:10]


# --- Phase R1-5: family-31 sensitivity ------------------------------------- #


def family_31_sensitivity(df: pd.DataFrame):
    sub = df[df["family_id"] != 31].reset_index(drop=True)
    return run_correlations(sub, label="exclude_fam31", weighted=False)


# --- Output formatting ---------------------------------------------------- #


def write_correlation_tex(combined: list[dict], path: str) -> None:
    rows = sorted(combined, key=lambda r: (
        float("inf") if math.isnan(r.get("bonferroni_p", float("nan"))) else r["bonferroni_p"]
    ))
    with open(path, "w", encoding="utf-8") as f:
        f.write("% PCF2-SESSION-R1 single-invariant correlation table\n")
        f.write("\\begin{tabular}{llrrrrrr}\n")
        f.write("\\hline\n")
        f.write("Invariant & Variant & $n$ & $\\rho$/$\\eta$ & $p$ & $p_{\\text{Bonf}}$ & $q_{\\text{BH}}$ & $m$ \\\\\n")
        f.write("\\hline\n")
        for r in rows:
            def fmt(x, d=3):
                if x is None or (isinstance(x, float) and math.isnan(x)):
                    return "--"
                if isinstance(x, float):
                    return f"{x:.{d}g}"
                return str(x)
            f.write(" & ".join([
                r["invariant"].replace("_", "\\_"),
                r["label"].replace("_", "\\_"),
                str(r["n"]),
                fmt(r["rho_or_eta"]),
                fmt(r["p"]),
                fmt(r.get("bonferroni_p", float("nan"))),
                fmt(r.get("bh_q", float("nan"))),
                str(r["m_tests"]),
            ]) + " \\\\\n")
        f.write("\\hline\n\\end{tabular}\n")


def main():
    df = assemble_dataframe()
    csv_path = os.path.join(SESSION_DIR, "assembled_data.csv")
    df.to_csv(csv_path, index=False)
    print(f"Assembled {len(df)} rows -> {csv_path}")
    print(df[["family_id", "bin", "Delta_3", "A_fit", "delta", "A_stderr"]].to_string(index=False))

    # Phase R1-2: unweighted Spearman battery
    res_unw = run_correlations(df, label="unweighted_full50", weighted=False)
    # Phase R1-4: weighted
    res_w = run_correlations(df, label="weighted_full50", weighted=True)
    # Phase R1-5: family-31 excluded
    res_no31 = family_31_sensitivity(df)

    # Phase R1-3: best 2-invariant pair regression
    pairs = best_pair_regression(df)

    # Halt-condition checks
    halt_triggered = False
    halt_reasons = []
    for r in res_unw + res_w + res_no31:
        rho = r["rho_or_eta"]
        bp = r.get("bonferroni_p", float("nan"))
        if not math.isnan(rho) and not math.isnan(bp):
            if abs(rho) > 0.6 and bp < 0.001:
                halt_triggered = True
                halt_reasons.append({
                    "invariant": r["invariant"], "label": r["label"],
                    "rho_or_eta": rho, "bonferroni_p": bp, "raw_p": r["p"]
                })

    # NaN columns flag (class number etc) - record in unexpected_finds, not halt
    nan_cols = [c for c in ["conductor", "class_number", "class_number_plus",
                             "regulator", "fund_unit_norm", "unit_density"]
                if df[c].isna().all()]

    # Significant unweighted vs weighted divergence
    flags_weight_diff = []
    for u, w in zip(res_unw, res_w):
        if u["invariant"] != w["invariant"]:
            continue
        if math.isnan(u["rho_or_eta"]) or math.isnan(w["rho_or_eta"]):
            continue
        if abs(u["rho_or_eta"]) > 1e-9 and (
            abs(w["rho_or_eta"]) > 2 * abs(u["rho_or_eta"])
            or abs(u["rho_or_eta"]) > 2 * abs(w["rho_or_eta"])
        ):
            flags_weight_diff.append({"invariant": u["invariant"],
                                       "rho_unw": u["rho_or_eta"],
                                       "rho_w": w["rho_or_eta"]})

    # Family-31 driven correlations
    flags_fam31 = []
    for full, no31 in zip(res_unw, res_no31):
        if full["invariant"] != no31["invariant"]:
            continue
        if math.isnan(full["rho_or_eta"]) or math.isnan(no31["rho_or_eta"]):
            continue
        if abs(full["rho_or_eta"]) > 0.3 and abs(no31["rho_or_eta"]) < 0.5 * abs(full["rho_or_eta"]):
            flags_fam31.append({"invariant": full["invariant"],
                                 "rho_full": full["rho_or_eta"],
                                 "rho_no31": no31["rho_or_eta"]})

    # Multi-invariant pair flag
    flags_pair = [dict(c1=p[0], c2=p[1], r2=p[2], adj_r2=p[3]) for p in pairs if p[3] > 0.4]

    # Write results.json
    summary = dict(
        task_id="PCF2-SESSION-R1",
        date="2026-05-01",
        n_families=len(df),
        n_with_class_number=int((~df["class_number"].isna()).sum()),
        delta_summary=dict(mean=float(df["delta"].mean()),
                           std=float(df["delta"].std(ddof=1)),
                           min=float(df["delta"].min()),
                           max=float(df["delta"].max())),
        halt_triggered=halt_triggered,
        halt_reasons=halt_reasons,
        unweighted_full50=res_unw,
        weighted_full50=res_w,
        unweighted_exclude_fam31=res_no31,
        best_pair_regressions=[
            dict(c1=p[0], c2=p[1], r2=p[2], adj_r2=p[3]) for p in pairs
        ],
        flags=dict(
            nan_columns=nan_cols,
            weighting_divergence=flags_weight_diff,
            family_31_driven=flags_fam31,
            pair_with_adj_r2_gt_0p4=flags_pair,
        ),
        environment_notes=(
            "pari/gp not available; cypari2 not installed; sage not available. "
            "Class numbers, conductors, regulators, and fundamental-unit norms "
            "recorded as NaN. The remaining 12 invariants (alpha_3, Delta_3_sign, "
            "absD_parity, absD_mod_8/12, log|D|, smallest_prime_div, omega(D), "
            "Galois_group, r1, r2, CM_disc) form the actual test family."
        ),
    )
    with open(os.path.join(SESSION_DIR, "results.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # correlation_table.tex (combined)
    combined = res_unw + res_w + res_no31
    write_correlation_tex(combined, os.path.join(SESSION_DIR, "correlation_table.tex"))

    # halt_log / discrepancy / unexpected
    with open(os.path.join(SESSION_DIR, "halt_log.json"), "w", encoding="utf-8") as f:
        json.dump(dict(halt=halt_triggered, reasons=halt_reasons), f, indent=2)
    with open(os.path.join(SESSION_DIR, "discrepancy_log.json"), "w", encoding="utf-8") as f:
        json.dump({}, f, indent=2)
    with open(os.path.join(SESSION_DIR, "unexpected_finds.json"), "w", encoding="utf-8") as f:
        json.dump(dict(
            nan_columns=nan_cols,
            weighting_divergence=flags_weight_diff,
            family_31_driven=flags_fam31,
            pair_with_adj_r2_gt_0p4=flags_pair,
        ), f, indent=2)

    print("\nHALT triggered:", halt_triggered)
    print("Best pair (adj R^2):",
          pairs[0] if pairs else "none")

    # AEAL claims (>= 3)
    self_path = os.path.abspath(__file__)
    self_hash = hashlib.sha256(open(self_path, "rb").read()).hexdigest()
    out_hash = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

    claims = []

    # Pick the strongest unweighted finding (max |rho_or_eta| with valid p).
    unw_valid = [r for r in res_unw if not math.isnan(r["rho_or_eta"]) and not math.isnan(r["p"])]
    no31_valid = [r for r in res_no31 if not math.isnan(r["rho_or_eta"]) and not math.isnan(r["p"])]
    if unw_valid:
        top = max(unw_valid, key=lambda r: abs(r["rho_or_eta"]))
        passes_strict = abs(top["rho_or_eta"]) > 0.4 and top.get("bonferroni_p", 1.0) < 0.01
        passes_loose = abs(top["rho_or_eta"]) > 0.4 and top.get("bonferroni_p", 1.0) < 0.05
        if passes_strict:
            claim_text = (
                f"single-invariant signal CONFIRMED on full 50: {top['invariant']} "
                f"correlates with delta = A_fit - 6 at Spearman rho={top['rho_or_eta']:.3f}, "
                f"raw p={top['p']:.3g}, Bonferroni-corrected p={top['bonferroni_p']:.3g} "
                f"(m={top['m_tests']}); passes |rho|>0.4 AND Bonferroni p<0.01"
            )
        elif passes_loose:
            claim_text = (
                f"borderline single-invariant signal on full 50: {top['invariant']} "
                f"with delta = A_fit - 6 has Spearman rho={top['rho_or_eta']:.3f}, "
                f"raw p={top['p']:.3g}, Bonferroni-corrected p={top['bonferroni_p']:.3g} "
                f"(m={top['m_tests']}) -- passes |rho|>0.4 and Bonferroni p<0.05 but "
                "narrowly misses Bonferroni p<0.01 strict cut"
            )
        else:
            claim_text = (
                "no single catalogue invariant correlates with delta = A_fit - 6 at "
                "|rho|>0.4, Bonferroni-corrected p<0.05 over the full 50; m="
                f"{top['m_tests']} valid tests; strongest hit: {top['invariant']} "
                f"rho={top['rho_or_eta']:.3f}, raw p={top['p']:.3g}, "
                f"Bonferroni p={top['bonferroni_p']:.3g}"
            )
    else:
        claim_text = "no valid Spearman test produced a finite p-value; pipeline degenerate"

    claims.append(dict(
        claim=claim_text,
        evidence_type="computation",
        dps=15,
        reproducible=True,
        script="r1_correlation_probe.py",
        output_hash=out_hash,
    ))

    # If family-31-excluded passes strict cut (|rho|>0.4, Bonferroni p<0.01), add a claim.
    if no31_valid:
        top31 = max(no31_valid, key=lambda r: abs(r["rho_or_eta"]))
        if abs(top31["rho_or_eta"]) > 0.4 and top31.get("bonferroni_p", 1.0) < 0.01:
            claims.append(dict(
                claim=(
                    f"with the noisy outlier family 31 (stderr 0.074) excluded, "
                    f"{top31['invariant']} passes the strict cut: "
                    f"Spearman rho={top31['rho_or_eta']:.3f}, raw p={top31['p']:.3g}, "
                    f"Bonferroni-corrected p={top31['bonferroni_p']:.3g}, "
                    f"BH q={top31.get('bh_q', float('nan')):.3g} (m={top31['m_tests']}); "
                    "family 31 MASKS rather than drives the signal"
                ),
                evidence_type="computation",
                dps=15,
                reproducible=True,
                script="r1_correlation_probe.py",
                output_hash=out_hash,
            ))

    # Weighted-vs-unweighted comparison claim
    w_valid = [r for r in res_w if not math.isnan(r["rho_or_eta"])]
    if unw_valid and w_valid:
        max_unw = max(abs(r["rho_or_eta"]) for r in unw_valid)
        max_w = max(abs(r["rho_or_eta"]) for r in w_valid)
        claims.append(dict(
            claim=(
                f"A_stderr-weighted Spearman peak |rho|={max_w:.3f} vs unweighted "
                f"peak |rho|={max_unw:.3f}; ratio={max_w/max(max_unw,1e-9):.2f} "
                "(does/does not exceed factor 2 -- see flags.weighting_divergence)"
            ),
            evidence_type="computation",
            dps=15,
            reproducible=True,
            script="r1_correlation_probe.py",
            output_hash=out_hash,
        ))

    # Family-31 sensitivity claim
    if unw_valid and no31_valid:
        max_unw = max(abs(r["rho_or_eta"]) for r in unw_valid)
        max_no31 = max(abs(r["rho_or_eta"]) for r in no31_valid)
        claims.append(dict(
            claim=(
                f"family-31-excluded Spearman peak |rho|={max_no31:.3f} vs full-50 "
                f"peak |rho|={max_unw:.3f}; "
                + ("family 31 drives apparent signal" if max_no31 < 0.5 * max_unw
                   else "signal (or null) persists without family 31")
            ),
            evidence_type="computation",
            dps=15,
            reproducible=True,
            script="r1_correlation_probe.py",
            output_hash=out_hash,
        ))

    # Multi-invariant pair claim (descriptive)
    if pairs:
        c1, c2, r2, adj = pairs[0]
        claims.append(dict(
            claim=(
                f"best 2-invariant rank regression: ({c1}, {c2}) with "
                f"adjusted R^2={adj:.3f}, raw R^2={r2:.3f} (descriptive only)"
            ),
            evidence_type="computation",
            dps=15,
            reproducible=True,
            script="r1_correlation_probe.py",
            output_hash=out_hash,
        ))

    with open(os.path.join(SESSION_DIR, "claims.jsonl"), "w", encoding="utf-8") as f:
        for c in claims:
            f.write(json.dumps(c) + "\n")

    print(f"Wrote {len(claims)} AEAL claims")


if __name__ == "__main__":
    main()
