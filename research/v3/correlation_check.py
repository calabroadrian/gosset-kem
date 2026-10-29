"""
AFA-12 v3 — Does the independence assumption hold for E8 message coding?
=========================================================================
Real ML-KEM-style PKE (ring Z_q[X]/(X^256+1), module rank k, CBD noise,
du/dv compression), with inflated noise (eta) so that failures are
observable. For each setting we measure the per-8-coefficient-block failure
rate three ways, for both encodings:
  real : full encryption/decryption (all natural correlations)
  iid  : coefficients drawn independently from the exact marginal of w
  bound: union bound used for the cryptographic-size figures
If real == iid (within statistical error), the independence assumption used
by us and by Liu & Sakzad (2026) is validated.
"""
import math
import sys
import time
import numpy as np
sys.path.insert(0, "../v3")
from e8_message_coding import kyber_noise, tail, conv, power

Q, N = 3329, 256
HALF = Q / 2


# ------------------------------------------------------------ ring arithmetic
def polymul(a, b):
    """Batched negacyclic product mod X^256+1. a,b: (..., 256) int arrays."""
    L = 2 * N
    fa = np.fft.rfft(a.astype(np.float64), L)
    fb = np.fft.rfft(b.astype(np.float64), L)
    c = np.rint(np.fft.irfft(fa * fb, L)).astype(np.int64)
    return (c[..., :N] - c[..., N:]) % Q


def cbd(rng, eta, shape):
    a = rng.integers(0, 2, size=shape + (eta,)).sum(-1)
    b = rng.integers(0, 2, size=shape + (eta,)).sum(-1)
    return a - b


def compress(x, d):
    return np.rint(x * (2 ** d) / Q).astype(np.int64) % (2 ** d)


def decompress(c, d):
    return np.rint(c * Q / (2 ** d)).astype(np.int64)


def center(x):
    x = x % Q
    return np.where(x > Q // 2, x - Q, x)


# ------------------------------------------------------------ encodings
def e8_reps():
    reps = []
    for m in range(128):  # D8 mod 2: even-weight binary vectors
        bits = [(m >> i) & 1 for i in range(7)]
        bits.append(sum(bits) % 2)
        reps.append(np.array(bits, float))
    for m in range(128):  # D8+1/2 mod 2: entries in {1/2, 3/2}, even # of 3/2
        bits = [(m >> i) & 1 for i in range(7)]
        bits.append(sum(bits) % 2)
        reps.append(0.5 + np.array(bits, float))
    R = np.array(reps)
    assert len({tuple(r) for r in R}) == 256
    return R


REPS = e8_reps()


def nearest_d8(z):
    f = np.rint(z)
    odd = (f.sum(-1) % 2) != 0
    err = z - f
    idx = np.argmax(np.abs(err), axis=-1)
    rows = np.nonzero(odd)
    sel = (rows[0], rows[1], idx[odd]) if z.ndim == 3 else (rows[0], idx[odd])
    f[sel] += np.where(err[sel] > 0, 1.0, -1.0)
    return f


def nearest_e8(z):
    a = nearest_d8(z)
    b = nearest_d8(z - 0.5) + 0.5
    da = ((z - a) ** 2).sum(-1)
    db = ((z - b) ** 2).sum(-1)
    return np.where((da <= db)[..., None], a, b)


def e8_decode_rep(y):
    """y: (..., 8) centered values -> E8 coset representative mod 2."""
    p = nearest_e8(y / HALF)
    return np.mod(p, 2.0)


# ------------------------------------------------------------ real scheme
def run_real(k, eta, du, dv, B, rng):
    A = rng.integers(0, Q, size=(k, k, N))
    s = cbd(rng, eta, (k, N))
    e = cbd(rng, eta, (k, N))
    t = sum(polymul(A[:, j], s[j][None]) for j in range(k)) + e  # (k,N): t_i = sum_j A_ij s_j
    t %= Q
    r = cbd(rng, eta, (B, k, N))
    e1 = cbd(rng, eta, (B, k, N))
    e2 = cbd(rng, eta, (B, N))
    # u_j = sum_i A_ij r_i + e1_j   (A^T r)
    u = np.stack([sum(polymul(A[i, j][None], r[:, i]) for i in range(k)) for j in range(k)], 1)
    u = (u + e1) % Q
    tr = sum(polymul(t[i][None], r[:, i]) for i in range(k))
    # messages
    msg_std = rng.integers(0, 2, size=(B, N))
    msg_e8 = rng.integers(0, 256, size=(B, N // 8))
    enc_std = (msg_std * int(round(HALF)))
    enc_e8 = np.rint(REPS[msg_e8] * HALF).astype(np.int64).reshape(B, N)
    out = {}
    for name, enc in (("std", enc_std), ("e8", enc_e8)):
        v = (tr + e2 + enc) % Q
        uc = decompress(compress(u, du), du)
        vc = decompress(compress(v, dv), dv)
        su = sum(polymul(s[i][None], uc[:, i]) for i in range(k))
        y = center(vc - su)
        if name == "std":
            dec = (np.abs(center(y - int(round(HALF)))) < Q / 4).astype(int)
            bad = (dec != msg_std).reshape(B, N // 8, 8).any(-1)
        else:
            rep = e8_decode_rep(center(y).reshape(B, N // 8, 8).astype(float))
            bad = np.abs(rep - REPS[msg_e8]).sum(-1) > 1e-9
        out[name] = bad
    return out


# ------------------------------------------------------------ iid model
def sample_marginal(w, size, rng):
    lo, a = w
    cdf = np.cumsum(a)
    cdf /= cdf[-1]
    idx = np.searchsorted(cdf, rng.random(size), side="right")
    return lo + np.minimum(idx, len(a) - 1)


def run_iid(w, blocks, rng):
    x = sample_marginal(w, (blocks, 8), rng)
    bad_std = (np.abs(x) >= Q / 4).any(-1)
    # E8: decode noise around the zero codeword (lattice decoder is translation-invariant)
    rep = e8_decode_rep(x.astype(float))
    bad_e8 = np.abs(rep).sum(-1) > 1e-9
    return bad_std, bad_e8


def bounds(w):
    p_std = 1 - (1 - 2 * tail(w, math.ceil(Q / 4))) ** 8
    p2 = tail(conv(w, w), Q / 2 - 1)
    p8 = tail(power(w, 8), Q - 4)
    return p_std, min(1.0, 112 * p2 + 128 * p8)


def ci(p, n):
    if n == 0:
        return 0.0
    return 1.96 * math.sqrt(max(p * (1 - p), 1e-300) / n)


if __name__ == "__main__":
    settings = [(2, 14, 10, 4), (2, 16, 10, 4), (2, 18, 10, 4), (2, 22, 10, 4), (3, 16, 10, 4), (3, 18, 10, 4)]
    n_keys = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    rng = np.random.default_rng(2026)
    print(f"{'k':>2}{'eta':>4}{'enc':>5}{'blocks':>10}{'real':>12}{'iid':>12}{'bound':>12}"
          f"{'real/iid':>10}")
    rows = []
    for k, eta, du, dv in settings:
        t0 = time.time()
        w = kyber_noise(k, eta, eta, du, dv)
        b_std, b_e8 = bounds(w)
        real = {"std": [], "e8": []}
        for _ in range(n_keys):  # fresh key every batch (correlation lives within a key)
            o = run_real(k, eta, du, dv, B, rng)
            for kk in real:
                real[kk].append(o[kk].ravel())
        real = {kk: np.concatenate(v) for kk, v in real.items()}
        nb = len(real["std"])
        i_std, i_e8 = run_iid(w, nb, rng)
        for name, rr, ii, bb in (("std", real["std"], i_std, b_std), ("e8", real["e8"], i_e8, b_e8)):
            pr, pi = rr.mean(), ii.mean()
            ratio = pr / pi if pi > 0 else float("nan")
            print(f"{k:>2}{eta:>4}{name:>5}{nb:>10}{pr:>12.3e}{pi:>12.3e}{bb:>12.3e}{ratio:>10.2f}"
                  f"   (±{ci(pr, nb):.1e} / ±{ci(pi, nb):.1e})", flush=True)
            rows.append([k, eta, 0 if name == "std" else 1, nb, pr, pi, bb])
        print(f"   [{time.time()-t0:.0f}s]", flush=True)
    np.save("../v3/correlation_check.npy", np.array(rows))
