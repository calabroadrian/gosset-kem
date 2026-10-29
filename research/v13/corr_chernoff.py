"""
AFA-12 v13 — rigorous-style tail of the in-block correlation statistic.
For a Voronoi direction u the statistic is Z = C^T Q C / C^T C, with Q = sum_l w_l (P^l + P^-l)/2 per polynomial
(P = negacyclic shift). Negacyclic shifts are diagonalised by the odd Fourier basis: the eigenvalues of Q are
mu_j = sum_l w_l cos(pi (2j+1) l / 256), j = 0..255, for each of the 2k polynomials.
P(Z >= x) = P(C^T (Q - x I) C >= 0) <= min_lambda prod (1 - 2 lambda sigma_p^2 (mu_j - x))^(-1/2)   (Gaussian C)
Checked against Monte Carlo (CBD + rounding entries) where MC has data, then used for the adversarial DFR.
"""
import json, math, sys, itertools, numpy as np
from scipy.optimize import minimize_scalar
from failure_boost import SETS, LN2
from corr_boost import sample_C, PATS
from corr_layout import anc_all, layout_struct, CONSEC, log_beta

J = np.arange(256)
WPAT = np.array([[p[i] * p[j] / 4 for (i, j) in itertools.combinations(range(8), 2)] for p in PATS])   # 64 x 28

def eig(lags, w):
    return np.cos(np.pi * np.outer(2 * J + 1, np.array(lags)) / 256) @ np.array(w)

def log_tail(mu, x, sig2):
    """Chernoff bound on P(Z >= x); sig2 = list of (variance, multiplicity) per polynomial type."""
    d = mu - x
    if d.max() <= 0: return -np.inf
    lmax = min(1 / (2 * s * d.max()) for s, _ in sig2) * 0.999999
    def f(lam):
        return sum(m * -0.5 * np.sum(np.log1p(-2 * lam * s * d)) for s, m in sig2)
    r = minimize_scalar(f, bounds=(0, lmax), method='bounded')
    return min(0.0, r.fun)

def directions(LS, kind):
    """list of (eigenvalue vector, count): type A (one lag, both signs) and type B (64 patterns) per block."""
    out = []; blocks = LS if kind == "distinct" else LS[:1]
    for pairs in blocks:
        lags = [d for (_, _, d) in pairs]
        for l in sorted(set(lags)):
            mu = eig([l], [1.0]); out.append((mu, 1)); out.append((-mu, 1))
        for w in WPAT:
            out.append((eig(lags, w), 1))
    return out

def log_alpha_R(x, dirs, sig2):
    terms = [math.log(c) + log_tail(mu, x, sig2) for mu, c in dirs]
    terms = [t for t in terms if np.isfinite(t)]
    if not terms: return -np.inf
    m = max(terms); return min(0.0, m + math.log(sum(math.exp(t - m) for t in terms)))

if __name__ == "__main__":
    nmc = int(sys.argv[1]) if len(sys.argv) > 1 else 60000
    NEW = json.load(open("layout.json"))["blocks"]
    LAY = {"consecutive": layout_struct(CONSEC), "distinct": layout_struct(NEW)}
    from corr_layout import rmax_layout, wsum
    print("TAIL MODEL CHECK: P(R_max >= r), Monte Carlo vs Chernoff (Gaussian quadratic form)")
    table = {}
    for S in SETS:
        if not S.e8: continue
        sr = S.s2 if False else float(np.sum(S.r[1] * S.r[0] ** 2)); sy = float(np.sum(S.y[1] * S.y[0] ** 2))
        sig2 = [(sr, S.k), (sy, S.k)]
        # the ratio Z uses T = C^T C: normalise variances so that E[C^T C] matches (scale-free statistic)
        Rs = {k_: [] for k_ in LAY}
        for b in range(0, nmc, 10000):
            C = sample_C(S, min(10000, nmc - b)); A = anc_all(C); rho = A / A[:, :1]
            for k_, LS in LAY.items(): Rs[k_].append(rmax_layout(rho, LS))
        line = f"{S.name:20s}"
        DIRS = {k_: directions(LS, k_) for k_, LS in LAY.items()}
        for k_ in LAY:
            R = np.concatenate(Rs[k_])
            chk = []
            for qq in (0.99, 0.999):
                r = float(np.quantile(R, qq)); chk.append(f"P>={r:.3f}: MC {1 - qq:.0e} bound {math.exp(log_alpha_R(r - 1, DIRS[k_], sig2)):.1e}")
            line += f" | {k_}: " + ", ".join(chk)
        print(line, flush=True)
        table[S.name] = (sig2, DIRS)
    print("\nADVERSARIAL DFR (log2 best per-ciphertext failure probability after grinding 2^X messages), Chernoff tail")
    print(f"{'scheme':20s} {'X':>4s} {'consecutive':>12s} {'distinct-lag':>13s} {'ML-KEM (norm only)':>20s}")
    std = {2: SETS[0], 3: SETS[3], 4: SETS[6]}
    for S in SETS:
        if not S.e8: continue
        sig2, DIRS = table[S.name]
        for X in (64, 128):
            out = []
            for k_, LS in LAY.items():
                mult = 14 * 32 if k_ == "consecutive" else 2 * wsum(LS)[1]
                xs = np.linspace(0, 0.9, 61); laR = [0.0] + [log_alpha_R(x, DIRS[k_], sig2) for x in xs[1:]]
                best = -1e9
                for f in np.linspace(1.0, 3.0, 81):
                    T = f * S.ET; laT = S.log_alpha(T)
                    for x, lr in zip(xs, laR):
                        if -(laT + lr) / LN2 <= X: best = max(best, log_beta(S, T, 1 + x, mult) / LN2)
                out.append(best)
            Sd = std[S.k]; bn = -1e9
            for f in np.linspace(1.0, 3.5, 251):
                T = f * Sd.ET
                if -Sd.log_alpha(T) / LN2 <= X: bn = max(bn, Sd.log_beta(T) / LN2)
            print(f"{S.name:20s} {X:4d} {out[0]:12.1f} {out[1]:13.1f} {bn:20.1f}", flush=True)
