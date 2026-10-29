"""
AFA-12 v13 — E8-specific boosting: correlated decryption errors inside a block.

Cov(delta_a, delta_b) = sigma_S^2 * A(b-a), A(d) = negacyclic autocorrelation of C = (r, e1+c_u) summed over the 2k
polynomials. Standard decoding looks at one coefficient (correlation irrelevant); the E8 decoder fails along
Voronoi directions u, whose variance is sigma_S^2 * u^T M u = 2 sigma_S^2 T * R(u), with
  type A (e_a +- e_b):         R = 1 +- rho_d,                          rho_d = A(d)/T
  type B ((+-1/2)^8, even):    R = 1 + sum_{a<b} s_a s_b rho_{b-a} / 4
An attacker grinds messages for ciphertexts with large max_u R (in addition to large ||C||).
1) Monte Carlo: distribution of rho_d and of R_max for random ciphertexts; check the Gaussian model of the tail.
2) Attack cost with both criteria (T and R), the boosted direction counted 14 times per block (conservative).
"""
import math, sys, itertools, numpy as np
from failure_boost import SETS, LN2, Scheme
from scipy.special import log_ndtr

rng = np.random.default_rng(12)

def sample_C(S, n):
    """n random ciphertext vectors C (shape n x 2k x 256): r ~ CBD(eta1), y = e1 + c_u."""
    def draw(dist, size):
        v, p = dist; return rng.choice(v, size=size, p=p)
    r = draw(S.r, (n, S.k, 256)); y = draw(S.y, (n, S.k, 256))
    return np.concatenate([r, y], axis=1).astype(np.float64)

def negacyclic_autocorr(C, dmax=7):
    """A(d) = sum_poly sum_j c_j * c_{j+d} with negacyclic wrap (x^256 = -1)."""
    out = []
    for d in range(1, dmax + 1):
        sh = np.concatenate([C[..., d:], -C[..., :d]], axis=-1)
        out.append(np.sum(C * sh, axis=(-1, -2)))
    return np.stack(out, axis=-1)

# type-B sign patterns (even number of minus signs), modulo global sign: 64 patterns
PATS = []
for bits in range(128):
    s = [1] + [(-1 if (bits >> i) & 1 else 1) for i in range(7)]
    if sum(1 for x in s if x < 0) % 2 == 0: PATS.append(np.array(s))
PATS = np.array(PATS)                                   # 64 x 8
# coefficient of rho_d in R for each pattern: sum over a<b with b-a=d of s_a s_b / 4
WB = np.zeros((len(PATS), 7))
for a in range(8):
    for b in range(a + 1, 8):
        WB[:, b - a - 1] += PATS[:, a] * PATS[:, b] / 4

def rmax(rho):
    RA = 1 + np.max(np.abs(rho), axis=-1)
    RB = 1 + np.max(np.abs(rho @ WB.T), axis=-1)
    return np.maximum(RA, RB)

def log_alpha_R(x, v):
    """Gaussian union bound on P(R_max >= 1 + x), rho_d ~ N(0, v) iid."""
    terms = [math.log(14) + log_ndtr(-x / math.sqrt(v))]
    for w in WB:
        terms.append(math.log(2) + log_ndtr(-x / math.sqrt(v * np.sum(w ** 2))))
    m = max(terms); return m + math.log(sum(math.exp(t - m) for t in terms))

def log_beta_dir(S, T, R, mult=14):
    """beta with mult boosted directions per block at variance ratio R (union with the unboosted rest)."""
    from failure_boost import log_Q, logsumexp, q
    s = math.sqrt(S.s2 * T); sd = s * math.sqrt(2)
    g, p = S.GA; m = p > 0; lA = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m]) / sd))
    g, p = S.GB; m = p > 0; lB = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m] / 2) / sd))
    sdr = sd * math.sqrt(R)
    g, p = S.GB; m = p > 0; lBr = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m] / 2) / sdr))
    g, p = S.GA; m = p > 0; lAr = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m]) / sdr))
    base = logsumexp(np.array([math.log(112) + lA, math.log(128) + lB, math.log(mult) + max(lAr, lBr)]))
    return math.log(32) + base

if __name__ == "__main__":
    nmc = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    print(f"MONTE CARLO of in-block correlation ({nmc} random ciphertexts per scheme)")
    VAR = {}
    for S in SETS:
        if not S.e8: continue
        rhos, Rs = [], []
        for b in range(0, nmc, 20000):
            C = sample_C(S, min(20000, nmc - b))
            T = np.sum(C ** 2, axis=(-1, -2)); rho = negacyclic_autocorr(C) / T[:, None]
            rhos.append(rho); Rs.append(rmax(rho))
        rho = np.concatenate(rhos); R = np.concatenate(Rs); v = float(np.var(rho)); VAR[S.name] = v
        # tail check: empirical vs Gaussian-union model
        chk = []
        for x in (3, 4, 4.5):
            xx = x * math.sqrt(v); emp = float(np.mean(R >= 1 + xx)); mod = math.exp(log_alpha_R(xx, v))
            chk.append(f"P(R>=1+{x:.1f}sd) emp {emp:.2e} model {mod:.2e}")
        print(f"  {S.name:20s} sd(rho)={math.sqrt(v):.4f}  R_max: mean {R.mean():.3f}, p99.99 {np.quantile(R, 0.9999):.3f}, max {R.max():.3f}")
        print("      " + " | ".join(chk), flush=True)
    print("\nATTACK COST with norm + correlation grinding (log2). Q <= 2^64 and unlimited queries.")
    for S in SETS:
        if not S.e8:
            c, u = S.attack(); print(f"  {S.name:20s} (standard decoding: correlation irrelevant)   Q<=2^64 work {c[0]:7.1f} | unlimited work {u[0]:6.1f}"); continue
        v = VAR[S.name]; bc = (1e9,); bu = (1e9,)
        for f in np.linspace(1.0, 3.5, 126):
            T = f * S.ET; la_T = S.log_alpha(T)
            for xs in np.linspace(0, 12, 49):
                x = xs * math.sqrt(v); la = la_T + (log_alpha_R(x, v) if xs > 0 else 0.0)
                la = min(la, 0.0); lb = log_beta_dir(S, T, 1 + x)
                W = -(la + lb) / LN2; Qq = -lb / LN2
                if W < bu[0]: bu = (W, Qq, f, xs)
                if Qq <= 64 and W < bc[0]: bc = (W, Qq, f, xs)
        c0, u0 = S.attack()
        cs = f"{bc[0]:7.1f} (T {bc[2]:.2f}x, R +{bc[3]:.1f}sd)" if len(bc) > 1 else "infeasible"
        print(f"  {S.name:20s} norm only: Q<=2^64 {c0[0]:7.1f}, unlimited {u0[0]:6.1f} | norm+corr: Q<=2^64 {cs}, unlimited {bu[0]:6.1f} (T {bu[2]:.2f}x, R +{bu[3]:.1f}sd)  delta {bu[0]-u0[0]:+.1f}", flush=True)
