"""
AFA-12 v3 — Reconcile our E8 DFR (Kyber1024, du=11, dv=5: 2^-335) with
Liu & Sakzad 2026 (2^-286).
Hypothesis A: model difference — exact convolution (ours) vs Gaussian/CLT
approximation of the per-coefficient noise (common in coding papers).
We compute both for standard and E8 encodings with the same union bounds.
"""
import math
import sys
import numpy as np
from math import erfc, sqrt
sys.path.insert(0, "../v3")
from e8_message_coding import kyber_noise, failure_standard, failure_e8, log2


def log2_gauss_tail(t, sd):
    """log2 P(N(0,sd^2) >= t), stable for large t."""
    z = t / sd
    if z < 30:
        return math.log2(0.5 * erfc(z / sqrt(2)))
    # asymptotic: phi(z)/z
    return (-z * z / 2 - math.log(z * sqrt(2 * math.pi))) / math.log(2)


def gauss_standard(var, q, n=256):
    return math.log2(n * 2) + log2_gauss_tail(math.ceil(q / 4), sqrt(var))


def gauss_e8(var, q, n=256):
    a = math.log2(112) + log2_gauss_tail(q / 2 - 1, sqrt(2 * var))
    b = math.log2(128) + log2_gauss_tail(q - 4, sqrt(8 * var))
    m = max(a, b)
    return math.log2(n // 8) + m + math.log2(2 ** (a - m) + 2 ** (b - m)), a, b


if __name__ == "__main__":
    q = 3329
    cases = [("ML-KEM-512  (10,4)", 2, 3, 2, 10, 4),
             ("ML-KEM-768  (10,4)", 3, 2, 2, 10, 4),
             ("ML-KEM-1024 (11,5)", 4, 2, 2, 11, 5),
             ("ML-KEM-1024 (10,5)", 4, 2, 2, 10, 5)]
    print(f"{'case':<20}{'var(w)':>9}{'std exact':>11}{'std gauss':>11}"
          f"{'E8 exact':>10}{'E8 gauss':>10}")
    for name, k, e1, e2, du, dv in cases:
        w = kyber_noise(k, e1, e2, du, dv, q)
        lo, a = w
        xs = np.arange(lo, lo + len(a))
        var = float((a * xs.astype(float) ** 2).sum())
        se = log2(failure_standard(w, q))
        ee = log2(failure_e8(w, q)[0])
        sg = gauss_standard(var, q)
        eg, pa, pb = gauss_e8(var, q)
        print(f"{name:<20}{var:>9.1f}{se:>11.1f}{sg:>11.1f}{ee:>10.1f}{eg:>10.1f}")
    print("\nLiu & Sakzad 2026 (Kyber1024, l=8 GTD): (11,5) -> 2^-286, (10,5) -> 2^-239")
