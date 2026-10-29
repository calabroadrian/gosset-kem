"""
AFA-12 v5 — Exact DFR: standard vs E8 vs Leech message coding in ML-KEM
========================================================================
Same noise (ML-KEM CBD + compression), exact per-coefficient distribution
(Kyber-spec method), union bound over minimal vectors, coefficients treated
independent (validated empirically for E8 in v3).

Codes (torus Z_q^l, lattice L with qZ^l in L):
  standard: 1 bit/coef, threshold q/4 per coefficient
  E8      : L = (q/2)E8,              2^8 points / 8 coefs  (1.0 bit/coef), d_min = 0.707 q
  Leech   : L = (q/8)Lambda24 (integral coords, min norm^2 32, contains 8Z^24),
            2^36 points / 24 coefs (1.5 bit/coef), d_min = 0.707 q
            256-bit message -> 8 Leech blocks = 192 coefficients of v
            (the other 64 coefficients of v are not transmitted)
Leech minimal vectors (196560):
  1104  of shape (+-4,+-4,0^22)   -> +-w_i +- w_j          >= q/2 - 1
  97152 of shape (+-2^8,0^16)     -> sum_8 +-w             >= q - 4
  98304 of shape (-+3,+-1^23)     -> 3w + sum_23 +-w       >= 2q - 13
(slack terms absorb rounding of the non-integer codeword coordinates)
"""
import math
import sys
import numpy as np
sys.path.insert(0, "../v3")
from e8_message_coding import kyber_noise, tail, conv, power

Q = 3329


def dilate(d, c):
    lo, a = d
    out = np.zeros((len(a) - 1) * c + 1)
    out[::c] = a
    return lo * c, out


def lg(p):
    return math.log2(p) if p > 0 else float("-inf")


def dfr_all(w):
    p1 = 2 * tail(w, math.ceil(Q / 4))
    p2 = tail(conv(w, w), Q / 2 - 1)
    w8 = power(w, 8)
    p8 = tail(w8, Q - 4)
    w23 = conv(power(w, 16), conv(power(w, 4), conv(power(w, 2), w)))
    p24 = tail(conv(dilate(w, 3), w23), 2 * Q - 13)
    std = 256 * p1
    e8 = 32 * (112 * p2 + 128 * p8)
    leech = 8 * (1104 * p2 + 97152 * p8 + 98304 * p24)
    return std, e8, leech, (p1, p2, p8, p24)


def ct_bytes(k, du, dv, nv):
    return (k * 256 * du + nv * dv) // 8


if __name__ == "__main__":
    sets = [("ML-KEM-512", 2, 3, 2, 10, 4, -139),
            ("ML-KEM-768", 3, 2, 2, 10, 4, -164),
            ("ML-KEM-1024", 4, 2, 2, 11, 5, -174)]
    print("A) SAME PARAMETERS (du, dv as in FIPS 203)")
    print(f"{'scheme':<12}{'target':>7}{'std':>9}{'E8':>9}{'Leech':>9}"
          f"{'ct std/E8':>11}{'ct Leech':>10}")
    for name, k, e1, e2, du, dv, tgt in sets:
        w = kyber_noise(k, e1, e2, du, dv)
        s, e, l, _ = dfr_all(w)
        print(f"{name:<12}{tgt:>7}{lg(s):>9.1f}{lg(e):>9.1f}{lg(l):>9.1f}"
              f"{ct_bytes(k, du, dv, 256):>11}{ct_bytes(k, du, dv, 192):>10}", flush=True)

    print("\nB) SMALLEST CIPHERTEXT meeting the published DFR target (search du, dv)")
    for name, k, e1, e2, du0, dv0, tgt in sets:
        best = {}
        for du in range(du0 - 3, du0 + 1):
            for dv in range(2, 7):
                w = kyber_noise(k, e1, e2, du, dv)
                s, e, l, _ = dfr_all(w)
                for code, dfr, nv in (("std", s, 256), ("E8", e, 256), ("Leech", l, 192)):
                    if lg(dfr) <= tgt:
                        c = ct_bytes(k, du, dv, nv)
                        if code not in best or c < best[code][0]:
                            best[code] = (c, du, dv, lg(dfr))
        base = ct_bytes(k, du0, dv0, 256)
        print(f"{name} (FIPS 203 ciphertext {base} B, target 2^{tgt}):")
        for code in ("std", "E8", "Leech"):
            if code in best:
                c, du, dv, d = best[code]
                print(f"   {code:<6} {c:>5} B  ({(c-base)/base*100:+.1f}%)  du={du} dv={dv}  DFR 2^{d:.1f}")
        sys.stdout.flush()
