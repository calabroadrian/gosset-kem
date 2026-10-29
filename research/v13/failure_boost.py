"""
AFA-12 v13 — decryption-failure (failure boosting) attack cost, D'Anvers et al. (PKC 2019) methodology.

Decrypted error per coefficient:  delta = e.r - s.(e1 + c_u) + e2 + c_v  =  <S, C_i> + G_i
  S = (e, -s)  secret, each coordinate CBD(eta1), variance sigma_S^2 = eta1/2
  C = (r, e1 + c_u) attacker-known (the attacker grinds messages m; r, e1, e2, c_u, c_v follow from m)
  G = e2 + c_v attacker-known offset.
Weak-ciphertext criterion: ||C||^2 >= T. Given T, delta_i ~ N(0, sigma_S^2 T) + G_i.
  alpha(T) = P(||C||^2 >= T)   (Chernoff upper bound: favours the attacker)
  beta(T)  = P(decryption failure | ||C||^2 = T)
     standard: 256 * P(|delta| > q/4)
     E8:       32 blocks * [112 P(+-d_a +-d_b > q/2) + 128 P(1/2 sum +-d_i > q/2)]   (Voronoi relevant vectors)
Attack: precompute 1/alpha candidates per weak ciphertext, query 1/beta weak ciphertexts to see one failure.
  work W = 1/(alpha*beta), queries Q = 1/beta, with the NIST limit Q <= 2^64.
"""
import math, numpy as np
from scipy.special import log_ndtr
from scipy.optimize import minimize_scalar

q = 3329
LN2 = math.log(2)

def cbd(eta):
    from math import comb
    v = np.arange(-eta, eta + 1); p = np.array([comb(2 * eta, eta + x) for x in v], float) / 4 ** eta
    return v, p

def rounding_err(d):
    x = np.arange(q)
    c = ((x.astype(np.int64) << d) + q // 2) // q % (1 << d)
    y = (c * q + (1 << (d - 1))) >> d
    e = (y - x) % q; e = np.where(e > q // 2, e - q, e)
    vals, cnt = np.unique(e, return_counts=True)
    return vals, cnt / q

def conv(a, b):
    va, pa = a; vb, pb = b
    lo = va[0] + vb[0]; p = np.convolve(pa, np.concatenate([pb]) if False else pb)
    # vectors must be contiguous integer supports
    full_a = np.zeros(va[-1] - va[0] + 1); full_a[va - va[0]] = pa
    full_b = np.zeros(vb[-1] - vb[0] + 1); full_b[vb - vb[0]] = pb
    p = np.convolve(full_a, full_b)
    return np.arange(lo, lo + len(p)), p

def negate(a): return (-a[0][::-1], a[1][::-1])

def log_Q(x):                               # log of Gaussian upper tail
    return log_ndtr(-x)

def logsumexp(a):
    m = np.max(a); return m + math.log(np.sum(np.exp(a - m))) if np.isfinite(m) else m

class Scheme:
    def __init__(self, name, k, eta1, eta2, du, dv, e8):
        self.name, self.k, self.e8 = name, k, e8
        self.nk = 256 * k; self.s2 = eta1 / 2
        self.r = cbd(eta1)
        self.y = conv(cbd(eta2), rounding_err(du))                 # e1 + c_u
        self.G = conv(cbd(eta2), rounding_err(dv))                 # e2 + c_v
        G2 = conv(self.G, self.G); self.GA = G2                    # +-G_a +-G_b (symmetric)
        G8 = self.G
        for _ in range(7): G8 = conv(G8, self.G)
        self.GB = G8                                               # sum of 8 (+-G_i), symmetric
        self.ET = self.nk * (np.sum(self.r[1] * self.r[0] ** 2) + np.sum(self.y[1] * self.y[0] ** 2))

    def log_alpha(self, T):                 # Chernoff bound on P(||C||^2 >= T), natural log
        if T <= self.ET: return 0.0
        def f(lam):
            mr = np.log(np.sum(self.r[1] * np.exp(lam * self.r[0] ** 2)))
            my = np.log(np.sum(self.y[1] * np.exp(lam * self.y[0] ** 2)))
            return self.nk * (mr + my) - lam * T
        res = minimize_scalar(f, bounds=(1e-9, 0.45), method='bounded')
        return min(0.0, res.fun)

    def log_beta(self, T, ratio=1.0):
        s = math.sqrt(self.s2 * T)
        if not self.e8:
            g, p = self.G
            lp = np.log(p[p > 0]); g = g[p > 0]
            a = np.concatenate([lp + log_Q((q / 4 - g) / s), lp + log_Q((q / 4 + g) / s)])
            return math.log(256) + logsumexp(a)
        sd = s * math.sqrt(2)
        g, p = self.GA; m = p > 0
        lA = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m]) / sd))
        g, p = self.GB; m = p > 0
        lB = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m] / 2) / sd))
        base = logsumexp(np.array([math.log(112) + lA, math.log(128) + lB]))
        if ratio == 1.0:
            return math.log(32) + base
        # one attacker-chosen direction (and its negation) with variance multiplied by `ratio`
        sdr = sd * math.sqrt(ratio)
        g, p = self.GB; m = p > 0
        lBr = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m] / 2) / sdr))
        g, p = self.GA; m = p > 0
        lAr = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m]) / sdr))
        lr = max(lAr, lBr)
        return math.log(32) + logsumexp(np.array([base, math.log(2) + lr]))

    def attack(self, qmax_log2=64, log_alpha_extra=None, ratio=1.0):
        best_c = (1e9, None, None, None); best_u = (1e9, None, None, None)
        for f in np.linspace(1.0, 3.5, 251):
            T = f * self.ET
            la = self.log_alpha(T) + (log_alpha_extra or 0.0); lb = self.log_beta(T, ratio)
            W = -(la + lb) / LN2; Qq = -lb / LN2
            if W < best_u[0]: best_u = (W, Qq, f, -la / LN2)
            if Qq <= qmax_log2 and W < best_c[0]: best_c = (W, Qq, f, -la / LN2)
        return best_c, best_u

SETS = [
    Scheme("ML-KEM-512", 2, 3, 2, 10, 4, False), Scheme("AFA12-512-Compact", 2, 3, 2, 9, 4, True), Scheme("AFA12-512-MaxSec", 2, 4, 3, 10, 4, True),
    Scheme("ML-KEM-768", 3, 2, 2, 10, 4, False), Scheme("AFA12-768-Compact", 3, 2, 2, 9, 5, True), Scheme("AFA12-768-MaxSec", 3, 3, 2, 10, 4, True),
    Scheme("ML-KEM-1024", 4, 2, 2, 11, 5, False), Scheme("AFA12-1024-Compact", 4, 2, 2, 10, 4, True), Scheme("AFA12-1024-MaxSec", 4, 3, 2, 11, 5, True),
]

if __name__ == "__main__":
    print("FAILURE BOOSTING (norm criterion). log2 values. beta at average ciphertext ~ DFR scale check.")
    print(f"{'scheme':20s} {'log2 beta(E[T])':>16s} | {'with Q<=2^64: work':>20s} {'queries':>8s} {'T/E[T]':>7s} | {'no query limit: work':>21s} {'queries':>8s}")
    for S in SETS:
        b0 = S.log_beta(S.ET) / LN2
        c, u = S.attack()
        cs = f"{c[0]:8.1f} {c[1]:8.1f} {c[2]:7.2f}" if c[2] is not None else f"{'infeasible':>20s} {'-':>8s} {'-':>7s}"
        print(f"{S.name:20s} {b0:16.1f} | {cs:>38s} | {u[0]:21.1f} {u[1]:8.1f}", flush=True)
