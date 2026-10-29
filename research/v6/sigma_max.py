"""Largest sigma meeting DFR <= 2^-128 (exact convolution), q=3329, n=768, and its security."""
import math, sys
sys.path.insert(0, "../v6"); sys.path.insert(0, "../v2")
from sigma_sweep import dfr, lg, bits
n, q = 768, 3329
for code in (0, 1):
    lo, hi = 1.0, 1.7
    for _ in range(7):
        mid = (lo + hi) / 2
        d = dfr(n, q, mid)[code]
        lo, hi = (mid, hi) if lg(d) <= -128 else (lo, mid)
    name = "standard code" if code == 0 else "E8 code     "
    print(f"{name}: sigma_max = {lo:.3f}  -> security std noise {bits(n,q,lo):.0f} bits, "
          f"E8 noise {bits(n,q,lo/math.sqrt(2)):.0f} bits", flush=True)
