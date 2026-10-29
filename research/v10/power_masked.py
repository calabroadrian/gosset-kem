"""
AFA-12 v10 — simulated power leakage of the MASKED E8 decoder vs unmasked
=========================================================================
Same model as v9/power_sim.py: each share leaks HW(share) + N(0, sigma^2).
Fresh masks and fresh gadget randomness in every trace.
Tests on one 8-coefficient block (realistic ML-KEM-768 decryption noise):
  1. first-order TVLA, fixed vs random input, every trace point
  2. single-trace profiled LDA attack on message bits (same attack as v9)
  3. second-order TVLA: centered product of the two shares of each recorded value
     (expected to leak: first-order masking only defeats first-order attacks)
"""
import ctypes, sys
import numpy as np
sys.path.insert(0, "../v3"); sys.path.insert(0, "../v9")
import correlation_check as cc
from e8_message_coding import kyber_noise
from power_sim import e8_encode_block, lda_attack as lda_unmasked, tvla as tvla_unmasked

Q = 3329
M = ctypes.CDLL("../v10/libe8masked.so")
U16 = np.ctypeslib.ndpointer(np.uint16, flags="C_CONTIGUOUS")
U8P = ctypes.POINTER(ctypes.c_uint8)
M.e8_block_masked.argtypes = [U16, U16, U8P, U8P]
M.masked_trace.argtypes = [U16, ctypes.c_int]; M.masked_seed.argtypes = [ctypes.c_uint64]
rng = np.random.default_rng(7)
noise_marg = kyber_noise(3, 3, 2, 10, 4)


def traces(N, sigma, fixed_w=None):
    buf = np.zeros(20000, np.uint16)
    m = rng.integers(0, 256, N)
    if fixed_w is None:
        W = (e8_encode_block(m) + cc.sample_marginal(noise_marg, (N, 8), rng)) % Q
    else:
        W = np.tile(fixed_w, (N, 1))
    T = None
    a, b = ctypes.c_uint8(), ctypes.c_uint8()
    M.masked_seed(int(rng.integers(1, 2**63)))
    for i in range(N):
        w0 = rng.integers(0, Q, 8); w1 = (W[i] - w0) % Q
        M.masked_trace(buf, len(buf))
        M.e8_block_masked(w0.astype(np.uint16), w1.astype(np.uint16), ctypes.byref(a), ctypes.byref(b))
        n = M.masked_trace_len()
        if T is None: T = np.zeros((N, n), np.float32)
        T[i] = buf[:n]
    T += rng.normal(0, sigma, T.shape).astype(np.float32)
    return T, m


def welch(A, B):
    return (A.mean(0) - B.mean(0)) / np.sqrt(A.var(0) / len(A) + B.var(0) / len(B) + 1e-12)


def lda(Ttr, mtr, Tte, mte):
    mu, sd = Ttr.mean(0), Ttr.std(0) + 1e-9
    Ttr, Tte = (Ttr - mu) / sd, (Tte - mu) / sd
    # keep the 300 points with highest variance ratio to bound covariance size
    preds = np.zeros((len(Tte), 8), np.int64)
    for bit in range(8):
        y = (mtr >> bit) & 1
        m0, m1 = Ttr[y == 0].mean(0), Ttr[y == 1].mean(0)
        sel = np.argsort(-np.abs(m1 - m0))[:300]
        A, B = Ttr[:, sel], Tte[:, sel]
        C = np.cov(np.concatenate([A[y == 0] - m0[sel], A[y == 1] - m1[sel]]).T) + 1e-3 * np.eye(len(sel))
        wv = np.linalg.solve(C, m1[sel] - m0[sel])
        preds[:, bit] = (B @ wv > wv @ (m0[sel] + m1[sel]) / 2).astype(np.int64)
    truth = (mte[:, None] >> np.arange(8)) & 1
    return float((preds == truth).mean()), float((preds == truth).all(1).mean())


if __name__ == "__main__":
    sigma = 1.0
    fixed = (e8_encode_block(np.array([0x5A])) + cc.sample_marginal(noise_marg, (1, 8), rng))[0] % Q
    Tf, _ = traces(20000, sigma, fixed)
    Tr, _ = traces(20000, sigma)
    print(f"masked block: {Tf.shape[1]} leakage points per trace (both shares of every intermediate)")
    t1 = welch(Tf, Tr)
    print(f"[1] first-order TVLA, sigma={sigma}, 20k+20k traces: max|t| = {np.abs(t1).max():.2f}  "
          f"({'no first-order leak' if np.abs(t1).max() < 4.5 else 'LEAK'});  unmasked decoder: max|t| = {tvla_unmasked('e8', sigma):.1f}")
    Ttr, mtr = traces(30000, sigma); Tte, mte = traces(10000, sigma)
    ba, by = lda(Ttr, mtr, Tte, mte)
    ub, uy = lda_unmasked("e8", sigma)
    print(f"[2] single-trace LDA attack, sigma={sigma}: masked bits {ba*100:.1f}% bytes {by*100:.2f}%   |   "
          f"unmasked bits {ub*100:.1f}% bytes {uy*100:.2f}%   (random: 50% / 0.39%)")
    P = lambda T: (T[:, 0::2] - T[:, 0::2].mean(0)) * (T[:, 1::2] - T[:, 1::2].mean(0))
    t2 = welch(P(Tf), P(Tr))
    print(f"[3] second-order TVLA (product of the 2 shares of each value): max|t| = {np.abs(t2).max():.1f}  "
          f"({'second-order leak, as expected for 1st-order masking' if np.abs(t2).max() > 4.5 else 'none'})")
