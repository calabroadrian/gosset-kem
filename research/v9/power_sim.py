"""
AFA-12 v9 — Simulated power side-channel assessment of message decoding
=======================================================================
Pre-silicon methodology (standard in side-channel evaluation): every intermediate
value of the decoder leaks HW(value as 32-bit word) + Gaussian measurement noise.
Target: message decoding in decapsulation (the step exploited by published
single-trace / deep-learning attacks on Kyber implementations).
Compared under identical conditions, 8 coefficients -> 8 message bits per block:
  STD : ML-KEM poly_tomsg (official formula), intermediates per coefficient
  E8  : constant-time E8 decoder (e8fast scalar logic), intermediates per coordinate
Noise w: exact decryption-noise marginal of each scheme (ML-KEM-768 / AFA12-768-MaxSec).
Tests: (1) TVLA fixed-vs-random Welch t; (2) profiled LDA attack (Gaussian templates,
pooled covariance) recovering each message bit from ONE trace; byte recovery rate.
This is a simulation, not a physical measurement.
"""
import sys, math
import numpy as np
sys.path.insert(0, "../v3")
from e8_message_coding import kyber_noise
import correlation_check as cc

Q = 3329
rng = np.random.default_rng(42)


def hw32(x):
    x = np.asarray(x, dtype=np.int64) & 0xFFFFFFFF
    x = x.astype(np.uint64)
    c = np.zeros(x.shape, dtype=np.int64)
    for s in range(32):
        c += ((x >> np.uint64(s)) & np.uint64(1)).astype(np.int64)
    return c


def std_intermediates(w):                          # w: (N, 8) values in [0,q)
    t1 = (w << 1) + 1665
    t2 = (t1 * 80635) & 0xFFFFFFFF
    t3 = t2 >> 28
    bits = t3 & 1
    acc = np.cumsum(bits << np.arange(8), axis=1)   # running msg byte
    pts = [w, t1, t2, t3, bits, acc]
    return np.concatenate([hw32(p) for p in pts], axis=1), (bits << np.arange(8)).sum(1)


def e8_intermediates(w):                           # w: (N, 8) in [0,q) ; mirrors e8fast.c cand()
    y = np.where(w > Q // 2, w - Q, w)
    Y = 4 * y
    U = 2 * Q
    pts = [Y]
    res = {}
    for off in (0, Q):
        X = Y - off
        f = (X >= U // 2).astype(np.int64) - (X < -U // 2) - (X < -3 * U // 2)
        err = X - U * f
        a = np.abs(err)
        idx = np.argmax(a, axis=1)
        par = f.sum(1) & 1
        sg = np.where(err[np.arange(len(err)), idx] >= 0, 1, -1)
        f2 = f.copy(); f2[np.arange(len(f)), idx] += par * sg
        e2 = Y - off - U * f2
        d = (e2 * e2).sum(1)
        pts += [f, err, a, f2]
        res[off] = (f2, d, idx, par)
    useB = (res[Q][1] < res[0][1]).astype(np.int64)
    fsel = np.where(useB[:, None] == 1, res[Q][0], res[0][0])
    msg = ((fsel[:, :7] & 1) << np.arange(7)).sum(1) | (useB << 7)
    extra = np.stack([res[0][1] & 0xFFFFFFFF, res[Q][1] & 0xFFFFFFFF, useB, res[0][2], res[Q][2], msg], axis=1)
    return np.concatenate([hw32(p) for p in pts] + [hw32(extra)], axis=1), msg


def e8_encode_block(m):                            # m: (N,) bytes -> (N,8) codeword values
    bits = (m[:, None] >> np.arange(7)) & 1
    par = bits.sum(1) & 1
    b = np.concatenate([bits, par[:, None]], axis=1)
    coset = (m >> 7) & 1
    return np.where(coset[:, None] == 1, 832 + 1665 * b, 1665 * b)


def make_traces(kind, N, sigma, fixed=None):
    if kind == "std":
        w_noise = cc.sample_marginal(kyber_noise(3, 2, 2, 10, 4), (N, 8), rng)
        m = rng.integers(0, 256, N) if fixed is None else np.full(N, fixed)
        cw = ((m[:, None] >> np.arange(8)) & 1) * 1665
        L, dec = std_intermediates((cw + w_noise) % Q)
    else:
        w_noise = cc.sample_marginal(kyber_noise(3, 3, 2, 10, 4), (N, 8), rng)
        m = rng.integers(0, 256, N) if fixed is None else np.full(N, fixed)
        L, dec = e8_intermediates((e8_encode_block(m) + w_noise) % Q)
    T = L + rng.normal(0, sigma, L.shape)
    return T, m, dec


def tvla(kind, sigma, N=20000):
    Tf, _, _ = make_traces(kind, N, sigma, fixed=0x5A)
    Tr, _, _ = make_traces(kind, N, sigma)
    t = (Tf.mean(0) - Tr.mean(0)) / np.sqrt(Tf.var(0) / N + Tr.var(0) / N + 1e-12)
    return float(np.max(np.abs(t)))


def lda_attack(kind, sigma, Ntr=30000, Nte=10000):
    Ttr, mtr, _ = make_traces(kind, Ntr, sigma)
    Tte, mte, _ = make_traces(kind, Nte, sigma)
    mu, sd = Ttr.mean(0), Ttr.std(0) + 1e-9
    Ttr, Tte = (Ttr - mu) / sd, (Tte - mu) / sd
    preds = np.zeros((Nte, 8), dtype=np.int64)
    for b in range(8):
        y = (mtr >> b) & 1
        m0, m1 = Ttr[y == 0].mean(0), Ttr[y == 1].mean(0)
        C = np.cov(np.concatenate([Ttr[y == 0] - m0, Ttr[y == 1] - m1]).T) + 1e-6 * np.eye(Ttr.shape[1])
        wv = np.linalg.solve(C, m1 - m0)
        thr = wv @ (m0 + m1) / 2
        preds[:, b] = (Tte @ wv > thr).astype(np.int64)
    truth = (mte[:, None] >> np.arange(8)) & 1
    bit_acc = float((preds == truth).mean())
    byte_acc = float((preds == truth).all(1).mean())
    return bit_acc, byte_acc


if __name__ == "__main__":
    print(f"{'decoder':<6}{'noise sigma':>12}{'TVLA max|t|':>13}{'bit accuracy':>14}{'byte (8 bits) from 1 trace':>28}")
    for sigma in (0.5, 1.0, 2.0, 4.0, 8.0):
        for kind in ("std", "e8"):
            t = tvla(kind, sigma)
            ba, by = lda_attack(kind, sigma)
            print(f"{kind:<6}{sigma:>12.1f}{t:>13.1f}{ba*100:>13.1f}%{by*100:>27.2f}%", flush=True)
    print("(random guessing: bit 50%, byte 0.39%)")
