"""
AFA-12 v6 — Is there a sigma that makes the v2.0 design (plain LWE, n=768) work?
================================================================================
Design of the v2.0 document: plain LWE, A uniform n x n, s, e, r, e1, e2 all
from the same noise distribution (std sigma), 256-bit message, no compression.
For each sigma we compute, with validated tools:
  - security: primal estimate (validated vs Kyber), for standard noise and for
    integral-E8 noise (equivalent to standard noise with sigma/sqrt2, see v2)
  - decryption failure: EXACT distribution of w = e.r - s.e1 + e2 (convolution),
    standard 1-bit encoding and E8 message encoding (union bound, as v3)
Two moduli: q = 3329 (as in v2.0 / ML-KEM) and q = 2^15 (FrodoKEM-style).
Feasible = security >= 128 bits AND DFR <= 2^-128.
"""
import math
import sys
import numpy as np
sys.path.insert(0, "../v2"); sys.path.insert(0, "../v3")
from estimator import primal_beta
from e8_message_coding import to_arr, conv, power, tail, mul_dist, add_dist


def dgauss(s):
    R = int(math.ceil(12 * s)) + 1
    x = np.arange(-R, R + 1)
    p = np.exp(-x.astype(float) ** 2 / (2 * s * s)); p /= p.sum()
    return dict(zip(x.tolist(), p.tolist()))


def dfr(n, q, s):
    D = dgauss(s)
    prod = to_arr(mul_dist(D, D))
    w = conv(power(prod, 2 * n), to_arr(D))         # e.r - s.e1 (2n products) + e2
    p1 = 2 * tail(w, math.ceil(q / 4))
    p2 = tail(conv(w, w), q / 2 - 1)
    p8 = tail(power(w, 8), q - 4)
    std = 256 * p1
    e8 = 32 * (112 * p2 + 128 * p8)
    return std, e8


def lg(p):
    return math.log2(p) if p > 0 else -9999


def bits(n, q, s):
    b, _ = primal_beta(n, q, s)
    return 0.292 * b if b else float("inf")


if __name__ == "__main__":
    n = 768
    for q, sigmas in [(3329, [0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 12.97]),
                      (2 ** 15, [1.0, 2.0, 2.8, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0])]:
        print(f"\nq = {q}, n = {n}  (v2.0 sigma=32.5 in rho-convention = std 12.97)")
        print(f"{'std':>6}{'sec std':>9}{'sec E8noise':>12}{'DFR std':>10}{'DFR E8code':>12}   feasible")
        for s in sigmas:
            bs, be = bits(n, q, s), bits(n, q, s / math.sqrt(2))
            ds, de = dfr(n, q, s)
            f = []
            if bs >= 128 and lg(ds) <= -128: f.append("std noise+std code")
            if bs >= 128 and lg(de) <= -128: f.append("std noise+E8 code")
            if be >= 128 and lg(de) <= -128: f.append("E8 noise+E8 code")
            print(f"{s:>6.2f}{bs:>9.0f}{be:>12.0f}{lg(ds):>10.0f}{lg(de):>12.0f}   {', '.join(f) or '-'}",
                  flush=True)
