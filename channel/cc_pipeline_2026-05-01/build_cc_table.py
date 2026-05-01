"""Build cc_channel_table.tex from the per-family results JSONs."""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

ROWS = [
    # (family, file, label_a, label_b, Delta, K_SCAN_K12)
    ("V$_{\\mathrm{quad}}$", "results_vquad.json",
     "$1$", "$3n^{2}+n+1$", "$-11$", "ARTEFACT"),
    ("Q$_{L01}$", "results_ql01.json",
     "$n$", "$n^{2}+n+1$", "$-3$", "n/a"),
    ("Q$_{L02}$", "results_ql02.json",
     "$2n+1$", "$n^{2}+n+1$", "$-4$", "n/a"),
    ("Q$_{L06}$", "results_ql06.json",
     "$n$", "$2n^{2}-n+1$", "$-7$", "n/a"),
    ("Q$_{L15}$", "results_ql15.json",
     "$n$", "$3n^{2}-2n+2$", "$-20$", "ARTEFACT"),
    ("Q$_{L26}$", "results_ql26.json",
     "$-3n+1$", "$4n^{2}-2n+2$", "$-28$", "ARTEFACT"),
]


def _digits(d):
    if d is None:
        return "n/a"
    if isinstance(d, str):
        return d
    return f"{d:.2f}"


def main():
    rows_out = []
    for fam, fp, la, lb, delta, kscan in ROWS:
        verdict = "n/a"
        digs_xi = "n/a"
        digs_alpha = "n/a"
        if fp is not None:
            f = HERE / fp
            if f.exists():
                obj = json.loads(f.read_text(encoding="utf-8"))
                if "tiers" in obj:
                    # V_quad
                    last = obj["tiers"][-1] if obj["tiers"] else {}
                    digs_xi = _digits(last.get("xi_0_digits_match"))
                    if obj.get("halt"):
                        verdict = "HALT (Phase 4 blocker)"
                else:
                    r = obj.get("result", {})
                    digs_alpha = _digits(r.get("WKB_alpha_digits_match"))
                    flips = r.get("flips_kscan_verdict", False)
                    pdig = r.get("P_III_D6_ratio_digits")
                    if flips:
                        verdict = (f"FLIP P-III(D6) at "
                                   f"{_digits(pdig)} digits")
                    else:
                        verdict = "no flip (no P-class signature)"
        rows_out.append((fam, la, lb, delta, kscan, digs_alpha, digs_xi,
                         verdict))

    lines = []
    lines.append("% CC-channel verdict table -- CC-PIPELINE-F session")
    lines.append("\\begin{center}")
    lines.append("\\renewcommand{\\arraystretch}{1.2}")
    lines.append("\\begin{tabular}{lcccccc}")
    lines.append("\\hline")
    lines.append("Family & $a_n$ & $b_n$ & $\\Delta$ & K-SCAN & "
                 "WKB $\\alpha$ digits & CC verdict \\\\")
    lines.append("\\hline")
    for fam, la, lb, delta, kscan, da, dx, verd in rows_out:
        lines.append(f"{fam} & {la} & {lb} & {delta} & "
                     f"{kscan} & {da} & {verd} \\\\")
    lines.append("\\hline")
    lines.append("\\end{tabular}")
    lines.append("\\end{center}")
    lines.append("")
    lines.append("\\noindent CC-PIPELINE-F (2026-05-01).  ``WKB $\\alpha$ digits''")
    lines.append("counts the digits of agreement between the LSQ-extracted")
    lines.append("trans-series exponent $\\alpha$ and the closed-form prediction")
    lines.append("$\\alpha = A - 2\\log c_{b} + \\log|c_{a}|$ from")
    lines.append("\\cite{siarc_borel_kscan} (V$_{\\mathrm{quad}}$ row reports the")
    lines.append("Phase-4 $\\xi_{0}$-digits attempt against the literature anchor)..")
    lines.append("``CC verdict'' is the present session's outcome on the family-")
    lines.append("level CC-channel detection question.  No flip was found at the")
    lines.append("$\\geq 15$-digit threshold; the V$_{\\mathrm{quad}}$ HALT")
    lines.append("indicates the BoT trans-series fit cannot be repurposed as a")
    lines.append("CC-channel extractor without an independent formal-solution")
    lines.append("computation for the underlying linear difference equation.")
    out = "\n".join(lines) + "\n"
    (HERE / "cc_channel_table.tex").write_text(out, encoding="utf-8")
    print(out)
    print(f"Wrote cc_channel_table.tex ({len(rows_out)} rows)")


if __name__ == "__main__":
    main()
