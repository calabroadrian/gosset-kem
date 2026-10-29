"""
AFA-12 v5 — independence check for 24-coefficient Leech blocks.
The dominant Leech failure term is the linear form 3*w_a + sum_{23 others} +-w_i
over a 24-coefficient block. We measure P(form >= T) on REAL ML-KEM-style
decryption noise (full ring arithmetic, compression, fresh keys) and on
IID coefficients from the exact marginal, at inflated noise so events are
observable. Same for the 8-coefficient octad form (sum_8 +-w).
"""
import sys, time, math
import numpy as np
sys.path.insert(0, "../v3")
import correlation_check as cc
from e8_message_coding import kyber_noise

Q = 3329


def real_noise(k, eta, du, dv, B, rng):
    A = rng.integers(0, Q, size=(k, k, 256))
    s = cc.cbd(rng, eta, (k, 256)); e = cc.cbd(rng, eta, (k, 256))
    t = (sum(cc.polymul(A[:, j], s[j][None]) for j in range(k)) + e) % Q
    r = cc.cbd(rng, eta, (B, k, 256)); e1 = cc.cbd(rng, eta, (B, k, 256)); e2 = cc.cbd(rng, eta, (B, 256))
    u = (np.stack([sum(cc.polymul(A[i, j][None], r[:, i]) for i in range(k)) for j in range(k)], 1) + e1) % Q
    v = (sum(cc.polymul(t[i][None], r[:, i]) for i in range(k)) + e2) % Q      # message = 0
    uc = cc.decompress(cc.compress(u, du), du); vc = cc.decompress(cc.compress(v, dv), dv)
    return cc.center(vc - sum(cc.polymul(s[i][None], uc[:, i]) for i in range(k)))  # (B,256)


def forms(w, rng):
    """w: (N,256) noise -> leech form values (N*8,) and octad form values (N*8,)"""
    blk = w[:, :192].reshape(-1, 24).astype(np.int64)
    sg = rng.choice([-1, 1], size=blk.shape)
    wt = np.ones(24, np.int64); wt[0] = 3
    leech = (blk * sg * wt).sum(1)
    oct_ = w[:, :64].reshape(-1, 8).astype(np.int64)
    sg8 = rng.choice([-1, 1], size=oct_.shape)
    octad = (oct_ * sg8).sum(1)
    return leech, octad


if __name__ == "__main__":
    rng = np.random.default_rng(31)
    print(f"{'k':>2}{'eta':>4}{'form':>7}{'T':>7}{'samples':>11}{'real':>12}{'iid':>12}{'ratio':>8}")
    for k, eta, nkeys in [(2, 13, 3000), (2, 16, 1500)]:
        t0 = time.time()
        ws = [real_noise(k, eta, 10, 4, 100, rng) for _ in range(nkeys)]
        W = np.concatenate(ws)
        L_real, O_real = forms(W, rng)
        marg = kyber_noise(k, eta, eta, 10, 4)
        Wi = cc.sample_marginal(marg, W.shape, rng)
        L_iid, O_iid = forms(Wi, rng)
        for name, xr, xi in (("leech", L_real, L_iid), ("octad", O_real, O_iid)):
            for qtl in (1e-3, 1e-4):
                T = float(np.quantile(xi, 1 - qtl))
                pr, pi = (xr >= T).mean(), (xi >= T).mean()
                n = len(xr)
                print(f"{k:>2}{eta:>4}{name:>7}{T:>7.0f}{n:>11}{pr:>12.3e}{pi:>12.3e}{pr/pi:>8.2f}"
                      f"   (±{1.96*math.sqrt(pr*(1-pr)/n):.1e} / ±{1.96*math.sqrt(pi*(1-pi)/n):.1e})", flush=True)
        print(f"   [{time.time()-t0:.0f}s]", flush=True)
