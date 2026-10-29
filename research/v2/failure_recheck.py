"""
AFA-12 v2 — Phase 2 re-check with the correct (exact, unrounded) samplers.
Regev/Kyber-style decryption noise per coefficient:
    w = <e, r> - <s, e1> + e2     (n-dim inner products, all terms ~ noise)
Compares Z-Gaussian noise vs integral-E8 noise at EQUAL per-coordinate
variance. Question: does the E8 shape change the tail of w (failure rate)?
"""
import math
import sys
import numpy as np
sys.path.insert(0, "../v2")
from e8_construction_a import IntegralE8Sampler, ZGaussianSampler, matched_s_for_L


def calibrate_Z(var):
    lo, hi = 0.05, 20
    for _ in range(60):
        mid = (lo + hi) / 2
        if ZGaussianSampler(mid).coord_variance() < var:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def sample_w(smp, n, N, chunk=2000):
    out = []
    for _ in range(N // chunk):
        e = smp.sample(chunk * n).reshape(chunk, n)
        r = smp.sample(chunk * n).reshape(chunk, n)
        s = smp.sample(chunk * n).reshape(chunk, n)
        e1 = smp.sample(chunk * n).reshape(chunk, n)
        e2 = smp.sample(chunk * 8).reshape(chunk, 8)[:, 0]
        out.append((e * r).sum(1) - (s * e1).sum(1) + e2)
    return np.concatenate(out)


if __name__ == "__main__":
    n, N, var = 256, 200_000, 1.0
    sZ, sL = calibrate_Z(var), matched_s_for_L(var)
    Z = ZGaussianSampler(sZ, seed=11)
    L = IntegralE8Sampler(sL, seed=12)
    print(f"per-coord variance: Z={Z.coord_variance():.4f}  E8={L.coord_variance():.4f}")
    wZ, wL = sample_w(Z, n, N), sample_w(L, n, N)
    theory = 2 * n * var * var + var
    print(f"var(w): Z={wZ.var():.2f}  E8={wL.var():.2f}  theory={theory:.2f}")
    k = lambda w: float(((w - w.mean()) ** 4).mean() / w.var() ** 2 - 3)
    print(f"excess kurtosis of w: Z={k(wZ):+.4f}  E8={k(wL):+.4f}")
    sd = math.sqrt(theory)
    print(f"{'threshold':>12}{'P(|w|>t) Z':>14}{'P(|w|>t) E8':>14}")
    rows = []
    for mult in [2, 3, 4, 4.5]:
        t = mult * sd
        pz, pl = float((np.abs(wZ) > t).mean()), float((np.abs(wL) > t).mean())
        rows.append((mult, pz, pl))
        print(f"{mult:>9.1f} sd{pz:>14.2e}{pl:>14.2e}")
    # two-sample KS
    a, b = np.sort(wZ), np.sort(wL)
    grid = np.union1d(a, b)
    ks = float(np.max(np.abs(np.searchsorted(a, grid, side='right') / len(a)
                             - np.searchsorted(b, grid, side='right') / len(b))))
    crit = 1.36 * math.sqrt(2 / N)
    print(f"KS statistic: {ks:.5f}  (5% critical value {crit:.5f}) -> "
          f"{'indistinguishable' if ks < crit else 'different'}")
    np.save("../v2/failure_recheck.npy",
            np.array([[wZ.var(), wL.var(), k(wZ), k(wL), ks, crit]]))
