"""
AFA-12 v13 — correlation boosting for an arbitrary block layout (consecutive vs distinct-lag layout.json).
Cov(delta_x, delta_y) = sigma^2 * Anc(y - x), Anc = negacyclic autocorrelation of C summed over the 2k polys.
Monte Carlo of R_max and of the per-direction tails, then the adversarial DFR with grinding 2^32/2^64/2^128.
"""
import json, math, sys, itertools, numpy as np
from scipy.special import log_ndtr
from failure_boost import SETS, LN2, log_Q, logsumexp, q
from corr_boost import sample_C, PATS

def anc_all(C):
    """Anc(d), d = 0..255, for C of shape (n, 2k, 256)."""
    F = np.fft.rfft(C, 512, axis=-1); L = np.fft.irfft(F * np.conj(F), 512, axis=-1)   # L[d] = sum c_i c_{i+d}, L[512-d] = L(-d)
    L = L.sum(axis=1)                                                                  # sum over polys
    A = np.empty((C.shape[0], 256)); A[:, 0] = L[:, 0]
    for d in range(1, 256): A[:, d] = L[:, d] - L[:, 512 - (256 - d)]                  # L(d) - L(d-256)
    return A

def layout_struct(blocks):
    """per block: list of (a_idx, b_idx, lag d=y-x) for sorted positions."""
    out = []
    for blk in blocks:
        blk = sorted(blk); out.append([(i, j, blk[j] - blk[i]) for i, j in itertools.combinations(range(8), 2)])
    return out

def rmax_layout(rho_all, LS):
    """rho_all: (n, 256) = Anc/T. Returns max over blocks and Voronoi directions of the variance ratio R."""
    n = rho_all.shape[0]; best = np.ones(n)
    # weights of type-B patterns per pair index (64 x 28)
    W = np.array([[p[i] * p[j] / 4 for (i, j) in itertools.combinations(range(8), 2)] for p in PATS])
    for pairs in LS:
        lags = [d for (_, _, d) in pairs]; r = rho_all[:, lags]                       # (n, 28)
        RA = 1 + np.abs(r).max(-1); RB = 1 + np.abs(r @ W.T).max(-1)
        best = np.maximum(best, np.maximum(RA, RB))
    return best

CONSEC = [list(range(8 * b, 8 * b + 8)) for b in range(32)]

def wsum(LS):
    """max over blocks/patterns of sqrt(sum of squared weights per distinct lag) (type-B sd multiplier) and
    max number of (block,pair) entries sharing one lag (type-A multiplicity)."""
    W = np.array([[p[i] * p[j] / 4 for (i, j) in itertools.combinations(range(8), 2)] for p in PATS])
    mx = 0
    for pairs in LS:
        lags = np.array([d for (_, _, d) in pairs])
        for w in W:
            agg = {}
            for l, x in zip(lags, w): agg[l] = agg.get(l, 0) + x
            mx = max(mx, math.sqrt(sum(v * v for v in agg.values())))
    cnt = {}
    for pairs in LS:
        for (_, _, d) in pairs: cnt[d] = cnt.get(d, 0) + 1
    return mx, max(cnt.values())

def log_alpha_R_layout(x, v, LS, kind):
    """Gaussian union bound on P(R_max >= 1 + x)."""
    W = np.array([[p[i] * p[j] / 4 for (i, j) in itertools.combinations(range(8), 2)] for p in PATS])
    terms = []
    blocks = LS if kind == "distinct" else LS[:1]          # consecutive: all blocks share the same lags
    for pairs in blocks:
        lags = np.array([d for (_, _, d) in pairs]); ul = sorted(set(lags.tolist()))
        terms.append(math.log(2 * len(ul)) + log_ndtr(-x / math.sqrt(v)))
        for w in W:
            agg = {}
            for l, c in zip(lags, w): agg[l] = agg.get(l, 0) + c
            terms.append(math.log(2) + log_ndtr(-x / math.sqrt(v * sum(c * c for c in agg.values()))))
    m = max(terms); return m + math.log(sum(math.exp(t - m) for t in terms))

def log_beta(S, T, R, mult):
    s = math.sqrt(S.s2 * T); sd = s * math.sqrt(2)
    g, p = S.GA; m = p > 0; lA = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m]) / sd))
    g, p = S.GB; m = p > 0; lB = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m] / 2) / sd))
    sdr = sd * math.sqrt(R)
    g, p = S.GB; m = p > 0; lBr = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m] / 2) / sdr))
    g, p = S.GA; m = p > 0; lAr = logsumexp(np.log(p[m]) + log_Q((q / 2 - g[m]) / sdr))
    base = math.log(32) + logsumexp(np.array([math.log(112) + lA, math.log(128) + lB]))
    return logsumexp(np.array([base, math.log(mult) + max(lAr, lBr)]))

if __name__ == "__main__":
    nmc = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    NEW = json.load(open("layout.json"))["blocks"]
    LAY = {"consecutive": layout_struct(CONSEC), "distinct": layout_struct(NEW)}
    for k_, LS in LAY.items():
        mx, ma = wsum(LS); print(f"layout {k_:11s}: max type-B sd multiplier {mx:.3f}, max blocks-pairs sharing one lag {ma}")
    res = {}
    for S in SETS:
        if not S.e8: continue
        Rs = {k_: [] for k_ in LAY}; rhos = []
        for b in range(0, nmc, 10000):
            C = sample_C(S, min(10000, nmc - b)); A = anc_all(C); rho = A / A[:, :1]
            rhos.append(rho[:, 1:8].ravel())
            for k_, LS in LAY.items(): Rs[k_].append(rmax_layout(rho, LS))
        v = float(np.var(np.concatenate(rhos)))
        line = f"{S.name:20s} sd(rho) {math.sqrt(v):.4f}"
        for k_ in LAY:
            R = np.concatenate(Rs[k_]); res[(S.name, k_)] = (v, R)
            x = np.quantile(R, 0.999) - 1
            line += f" | {k_}: R p99.9 {1 + x:.3f} (model P {math.exp(log_alpha_R_layout(x, v, LAY[k_], k_)):.1e} vs 1e-3), max {R.max():.3f}"
        print(line, flush=True)
    print("\nADVERSARIAL DFR (log2 of the best per-ciphertext failure probability after grinding 2^X messages)")
    print(f"{'scheme':20s} {'X':>4s} {'consecutive':>12s} {'distinct-lag':>13s} {'ML-KEM same level (norm only)':>30s}")
    std = {2: SETS[0], 3: SETS[3], 4: SETS[6]}
    for S in SETS:
        if not S.e8: continue
        for X in (64, 128):
            out = []
            for k_, LS in LAY.items():
                v = res[(S.name, k_)][0]; mult = 14 * 32 if k_ == "consecutive" else 2 * wsum(LS)[1]
                best = -1e9
                for f in np.linspace(1.0, 3.0, 101):
                    T = f * S.ET; laT = S.log_alpha(T)
                    for xs in np.linspace(0, 16, 65):
                        la = laT + (log_alpha_R_layout(xs * math.sqrt(v), v, LS, k_) if xs > 0 else 0.0)
                        if -min(la, 0.0) / LN2 <= X: best = max(best, log_beta(S, T, 1 + xs * math.sqrt(v), mult) / LN2)
                out.append(best)
            Sd = std[S.k]; bn = -1e9
            for f in np.linspace(1.0, 3.5, 251):
                T = f * Sd.ET
                if -Sd.log_alpha(T) / LN2 <= X: bn = max(bn, Sd.log_beta(T) / LN2)
            print(f"{S.name:20s} {X:4d} {out[0]:12.1f} {out[1]:13.1f} {bn:30.1f}", flush=True)
