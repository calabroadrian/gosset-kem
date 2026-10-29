"""
AFA-12 v14 — tight tail for the in-block correlation statistic (Lugannani-Rice saddlepoint) instead of the
Chernoff bound used in v13, then the adversarial DFR and the total first-failure cost again.
Y = C^T (Q - xI) C = sum_p sum_j sigma_p^2 (mu_j - x) g^2   (g standard normal, eigenvalues of v13)
K(l) = -1/2 sum log(1 - 2 l a),  a = sigma^2 (mu - x).  P(Y >= 0) ~ 1 - Phi(w) + phi(w) (1/u - 1/w),
w = sqrt(-2 K(l*)), u = l* sqrt(K''(l*)), K'(l*) = 0.
Validation: single-direction tails against Monte Carlo with the true CBD/rounding entries, and P(R_max >= r).
"""
import sys, json, math, itertools, numpy as np
sys.path.insert(0, '../v13')
from scipy.optimize import brentq
from scipy.stats import norm
from failure_boost import SETS, LN2
from corr_boost import sample_C, PATS
from corr_layout import anc_all, layout_struct, CONSEC, log_beta, wsum, rmax_layout
from corr_chernoff import directions, eig, log_tail as log_tail_chernoff

def log_tail_saddle(mu, x, sig2):
    a = np.concatenate([np.repeat(s * (mu - x), m) for s, m in sig2])
    if a.max() <= 0: return -np.inf
    hi = 1 / (2 * a.max())
    Kp = lambda l: np.sum(a / (1 - 2 * l * a))
    if Kp(0) >= 0: return 0.0
    l = brentq(Kp, 0, hi * (1 - 1e-12))
    K = -0.5 * np.sum(np.log1p(-2 * l * a)); K2 = np.sum(2 * a * a / (1 - 2 * l * a) ** 2)
    if -2 * K <= 0: return 0.0
    w = math.sqrt(-2 * K); u = l * math.sqrt(K2)
    # log of 1-Phi(w) + phi(w)(1/u - 1/w) = log phi(w) + log(Mills(w) + 1/u - 1/w), Mills(w) = Q(w)/phi(w)
    from scipy.special import erfcx
    mills = math.sqrt(math.pi / 2) * erfcx(w / math.sqrt(2))
    inner = mills + 1 / u - 1 / w
    if inner <= 0: inner = 1 / u
    return min(0.0, norm.logpdf(w) + math.log(inner))

def log_alpha_R(x, dirs, sig2, tail):
    terms = [math.log(c) + tail(mu, x, sig2) for mu, c in dirs]
    terms = [t for t in terms if np.isfinite(t)]
    if not terms: return -np.inf
    m = max(terms); return min(0.0, m + math.log(sum(math.exp(t - m) for t in terms)))

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    NEW = json.load(open("../v13/layout.json"))["blocks"]
    LAY = {"consecutive": layout_struct(CONSEC), "distinct": layout_struct(NEW)}
    names = sys.argv[2:] or [S.name for S in SETS if S.e8]
    for S in SETS:
        if S.name not in names: continue
        sr = float(np.sum(S.r[1] * S.r[0] ** 2)); sy = float(np.sum(S.y[1] * S.y[0] ** 2)); sig2 = [(sr, S.k), (sy, S.k)]
        if mode == "validate":
            n = 200000; Z1 = []; Z2 = []; Rs = {k_: [] for k_ in LAY}
            W = np.array([[p[i] * p[j] / 4 for (i, j) in itertools.combinations(range(8), 2)] for p in PATS])
            lagsB = [d for (_, _, d) in LAY["distinct"][0]]
            for b in range(0, n, 10000):
                C = sample_C(S, 10000); A = anc_all(C); rho = A / A[:, :1]
                Z1.append(rho[:, 1]); Z2.append(rho[:, lagsB] @ W[0])
                for k_, LS in LAY.items(): Rs[k_].append(rmax_layout(rho, LS))
            Z1 = np.concatenate(Z1); Z2 = np.concatenate(Z2)
            print(f"{S.name}: single-direction tails, Monte Carlo vs saddlepoint vs Chernoff")
            for nm, Z, mu in (("type A, lag 1", Z1, eig([1], [1.0])), ("type B, distinct block 0, pattern 0", Z2, eig(lagsB, W[0]))):
                for qq in (1e-2, 1e-3, 1e-4):
                    x = float(np.quantile(Z, 1 - qq))
                    print(f"   {nm:36s} P(Z>={x:.4f}): MC {qq:.0e}  saddle {math.exp(log_tail_saddle(mu, x, sig2)):.2e}  chernoff {math.exp(log_tail_chernoff(mu, x, sig2)):.2e}")
            for k_, LS in LAY.items():
                R = np.concatenate(Rs[k_]); D = directions(LS, k_)
                for qq in (1e-2, 1e-3, 1e-4):
                    r = float(np.quantile(R, 1 - qq))
                    print(f"   R_max {k_:11s} P(R>={r:.3f}): MC {qq:.0e}  saddle-union {math.exp(log_alpha_R(r - 1, D, sig2, log_tail_saddle)):.2e}")
            sys.stdout.flush()
        else:   # costs
            out = []
            for k_, LS in LAY.items():
                D = directions(LS, k_); mult = 14 * 32 if k_ == "consecutive" else 2 * wsum(LS)[1]
                xs = np.linspace(0, 1.2, 49); laR = [0.0] + [log_alpha_R(x, D, sig2, log_tail_saddle) for x in xs[1:]]
                adv = {}; best = (1e9, 0)
                for f in np.linspace(1.0, 3.5, 101):
                    T = f * S.ET; laT = S.log_alpha(T)
                    for x, lr in zip(xs, laR):
                        lb = log_beta(S, T, 1 + x, mult); la = laT + lr
                        for X in (64, 128):
                            if -la / LN2 <= X: adv[X] = max(adv.get(X, -1e9), lb / LN2)
                        W = -(la + lb) / LN2; Qq = -lb / LN2
                        if Qq <= 64 and W < best[0]: best = (W, Qq)
                out.append((adv[64], adv[128], best[0]))
            print(f"{S.name:20s} consecutive: DFR@2^64 {out[0][0]:7.1f} @2^128 {out[0][1]:7.1f} cost {out[0][2]:6.1f} | "
                  f"distinct: DFR@2^64 {out[1][0]:7.1f} @2^128 {out[1][1]:7.1f} cost {out[1][2]:6.1f}", flush=True)
