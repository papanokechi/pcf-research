"""Diagnose Gevrey class of the trans-series h_k for V_quad."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cc_pipeline import extract_h
from mpmath import mp, mpf, sqrt as mpsqrt, log as mplog


def a_fn(n):
    return mpf(1)


def b_fn(n):
    return mpf(3) * n * n + n + mpf(1)


SPEC = {"family": "V_quad", "a_fn": a_fn, "b_fn": b_fn}

mp.dps = 800
ext = extract_h(SPEC, depth=180, dps=800, K=18, n_lo=12, n_hi=140)
h = ext["h"]

print("|h_k|^(1/k) sequence (potential geometric inverse-radius):")
for k, v in enumerate(h, 1):
    eta = abs(v) ** (mpf(1) / k)
    if v != 0:
        slog = mplog(abs(v)) / k
    else:
        slog = "-inf"
    print(f"  k={k:2d}: |h_k|^(1/k) = {mp.nstr(eta, 14)}, "
          f"log|h_k|/k = {mp.nstr(slog, 14) if v != 0 else slog}")

print()
print("ratio h_{k+1}/h_k:")
for k in range(1, len(h)):
    if h[k - 1] != 0:
        r = h[k] / h[k - 1]
        print(f"  k={k:2d}: r_k = {mp.nstr(r, 14)}")

print()
print("ratio (h_{k+1}/h_k) / (k):  -- if Gevrey-1, this -> 1/xi_0")
for k in range(1, len(h)):
    if h[k - 1] != 0:
        r = h[k] / h[k - 1]
        print(f"  k={k:2d}: r_k/k = {mp.nstr(r / k, 14)}")

print()
print("2/sqrt(3) =", mp.nstr(mpf(2) / mpsqrt(3), 14))
print("sqrt(3)/2 =", mp.nstr(mpsqrt(3) / 2, 14))
print("alpha (V_quad WKB) =", mp.nstr(ext["alpha"], 14))
print("4 - 2 log 3 =", mp.nstr(mpf(4) - 2 * mplog(3), 14))
