"""
AFA-12 v3 — E8 in the right place: message encoding ("Binary Embedding" layer)
==============================================================================
Noise stays standard (CBD, as in ML-KEM). Only the message encoding changes:
  standard : 1 bit per coefficient, m_i -> round(q/2)*m_i, decode |w_i| < q/4
  E8 code  : 8 bits per 8 coefficients, message = point of (q/2)*E8 mod q*Z^8
             (E8 contains 2Z^8, so (q/2)E8 / qZ^8 has exactly 256 points)
Same rate (1 bit/coef). Minimum distance grows by sqrt(2).

Per-coefficient decryption noise is computed EXACTLY (direct convolution of
positive distributions, as in the Kyber reference failure scripts):
  w = e.r + e2 + dv - s.(e1 + du)
Failure (union bounds, coefficients treated independent, as in Kyber spec):
  standard: per 8 coefs  <= sum_i P(|w_i| > q/4)
  E8      : per 8 coefs  <= 112 * P(w_i + w_j >= q/2 - 1)
                           + 128 * P(sum_{8} +-w_i >= q - 4)
            (240 minimal vectors of E8; -1/-4 absorb rounding of q/2, q/4)
"""
import math
import numpy as np


def cbd(eta):
    vals = np.arange(-eta, eta + 1)
    p = np.array([math.comb(2 * eta, eta + v) for v in vals], float) / 4 ** eta
    return dict(zip(vals.tolist(), p.tolist()))


def compress_error(q, d):
    """Distribution of Decompress(Compress(x)) - x (centered) for x uniform mod q."""
    out = {}
    for x in range(q):
        c = round(x * (2 ** d) / q) % (2 ** d)
        y = round(c * q / (2 ** d))
        e = (y - x) % q
        if e > q // 2:
            e -= q
        out[e] = out.get(e, 0) + 1 / q
    return out


def add_dist(a, b):
    out = {}
    for x, px in a.items():
        for y, py in b.items():
            out[x + y] = out.get(x + y, 0) + px * py
    return out


def mul_dist(a, b):
    out = {}
    for x, px in a.items():
        for y, py in b.items():
            out[x * y] = out.get(x * y, 0) + px * py
    return out


def to_arr(d):
    lo, hi = min(d), max(d)
    arr = np.zeros(hi - lo + 1)
    for k, v in d.items():
        arr[k - lo] += v
    return lo, arr


def conv(a, b):
    (la, aa), (lb, bb) = a, b
    c = np.convolve(aa, bb)
    c[c < 1e-300] = 0.0
    nz = np.nonzero(c)[0]
    return la + lb + nz[0], c[nz[0]:nz[-1] + 1]


def power(a, n):
    result, base = (0, np.array([1.0])), a
    while n:
        if n & 1:
            result = conv(result, base)
        base = conv(base, base)
        n >>= 1
    return result


def kyber_noise(k, eta1, eta2, du, dv, q=3329, n=256):
    s = cbd(eta1)
    e1_du = add_dist(cbd(eta2), compress_error(q, du))
    e2_dv = add_dist(cbd(eta2), compress_error(q, dv))
    t1 = to_arr(mul_dist(s, e1_du))            # s * (e1 + du)
    t2 = to_arr(mul_dist(cbd(eta1), cbd(eta1)))  # e * r
    w = conv(power(t1, k * n), power(t2, k * n))
    w = conv(w, to_arr(e2_dv))
    return w


def tail(dist, t):
    lo, a = dist
    xs = np.arange(lo, lo + len(a))
    return float(a[xs >= t].sum())


def log2(p):
    return math.log2(p) if p > 0 else float("-inf")


def failure_standard(w, q, n=256):
    p = tail(w, math.ceil(q / 4)) * 2          # symmetric
    return n * p


def failure_e8(w, q, n=256):
    pair = conv(w, w)
    p2 = tail(pair, q / 2 - 1)
    eight = power(w, 8)
    p8 = tail(eight, q - 4)
    per_block = 112 * p2 + 128 * p8
    return (n // 8) * per_block, p2, p8


if __name__ == "__main__":
    q = 3329
    configs = [
        # name, k, eta1, eta2, du, dv, published log2 failure
        ("ML-KEM-512 ", 2, 3, 2, 10, 4, -139),
        ("ML-KEM-768 ", 3, 2, 2, 10, 4, -164),
        ("ML-KEM-1024", 4, 2, 2, 11, 5, -174),
    ]
    print("VALIDATION + E8 MESSAGE CODING (same noise, same rate)")
    print(f"{'scheme':<12}{'published':>10}{'ours std':>10}{'E8 code':>10}")
    res = {}
    for name, k, e1, e2, du, dv, pub in configs:
        w = kyber_noise(k, e1, e2, du, dv, q)
        fs = failure_standard(w, q)
        fe, _, _ = failure_e8(w, q)
        res[name] = (w, fs, fe)
        print(f"{name:<12}{pub:>10}{log2(fs):>10.1f}{log2(fe):>10.1f}")

    print("\nSPENDING THE GAIN on ML-KEM-512 (keep failure <= published 2^-139):")
    print(f"{'du':>4}{'dv':>4}{'ct bytes':>10}{'std log2 d':>12}{'E8 log2 d':>11}")
    for du, dv in [(10, 4), (10, 3), (9, 4), (9, 3), (8, 4), (8, 3)]:
        w = kyber_noise(2, 3, 2, du, dv, q)
        fs, (fe, _, _) = failure_standard(w, q), failure_e8(w, q)
        ct = (2 * 256 * du + 256 * dv) // 8
        print(f"{du:>4}{dv:>4}{ct:>10}{log2(fs):>12.1f}{log2(fe):>11.1f}")
