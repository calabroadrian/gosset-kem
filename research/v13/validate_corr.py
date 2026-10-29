"""
AFA-12 v13 — validate the correlation-boosting model with REAL ring encryptions (inflated noise so failures
are measurable). For every ciphertext we know C = (r, e1 + c_u); we compute T = ||C||^2 and the in-block
autocorrelations rho_d, predict its per-block failure probability with the Gaussian model (all 240 Voronoi
directions, variance 2 sigma^2 T R_u), and compare predicted vs observed failures per R-bucket.
"""
import sys, math, numpy as np
sys.path.insert(0, '../v3')
from correlation_check import polymul, cbd, compress, decompress, center, e8_decode_rep, REPS, HALF, Q, N
from scipy.special import log_ndtr

def pairs():
    P = []
    for a in range(8):
        for b in range(a + 1, 8):
            P.append((a, b))
    return P
PAIRS = pairs()
PATS = []
for bits in range(128):
    s = [1] + [(-1 if (bits >> i) & 1 else 1) for i in range(7)]
    if sum(1 for x in s if x < 0) % 2 == 0: PATS.append(np.array(s))
PATS = np.array(PATS)
WB = np.zeros((len(PATS), 7))
for a in range(8):
    for b in range(a + 1, 8):
        WB[:, b - a - 1] += PATS[:, a] * PATS[:, b] / 4

def offset_dist(eta2, dv):
    from math import comb
    v = np.arange(-eta2, eta2 + 1); p = np.array([comb(2 * eta2, eta2 + x) for x in v], float) / 4 ** eta2
    x = np.arange(Q); y = decompress(compress(x, dv), dv); e = center(y - x)
    ev, ec = np.unique(e, return_counts=True)
    full = np.zeros(ev[-1] - ev[0] + 1); full[ev - ev[0]] = ec / Q
    g = np.convolve(p, full); lo = v[0] + ev[0]
    return np.arange(lo, lo + len(g)), g

def self_conv(d, n):
    v, p = d; out_v, out_p = v, p
    for _ in range(n - 1):
        out_p = np.convolve(out_p, p); out_v = np.arange(out_v[0] + v[0], out_v[0] + v[0] + len(out_p))
    return out_v, out_p

def run(k, eta, du, dv, B, rng):
    A = rng.integers(0, Q, size=(k, k, N))
    s = cbd(rng, eta, (k, N)); e = cbd(rng, eta, (k, N))
    t = (sum(polymul(A[:, j], s[j][None]) for j in range(k)) + e) % Q
    r = cbd(rng, eta, (B, k, N)); e1 = cbd(rng, eta, (B, k, N)); e2 = cbd(rng, eta, (B, N))
    u = np.stack([sum(polymul(A[i, j][None], r[:, i]) for i in range(k)) for j in range(k)], 1)
    u = (u + e1) % Q
    tr = sum(polymul(t[i][None], r[:, i]) for i in range(k))
    msg = rng.integers(0, 256, size=(B, N // 8))
    enc = np.rint(REPS[msg] * HALF).astype(np.int64).reshape(B, N)
    v = (tr + e2 + enc) % Q
    uc = decompress(compress(u, du), du); vc = decompress(compress(v, dv), dv)
    su = sum(polymul(s[i][None], uc[:, i]) for i in range(k))
    y = center(vc - su)
    rep = e8_decode_rep(center(y).reshape(B, N // 8, 8).astype(float))
    bad = np.abs(rep - REPS[msg]).sum(-1) > 1e-9                     # (B, 32)
    cu = center(uc - u)
    C = np.concatenate([r, e1 + cu], axis=1).astype(np.float64)    # (B, 2k, 256)
    return bad, C

def autocorr(C):
    out = []
    for d in range(1, 8):
        sh = np.concatenate([C[..., d:], -C[..., :d]], axis=-1)
        out.append(np.sum(C * sh, axis=(-1, -2)))
    return np.stack(out, -1)

if __name__ == "__main__":
    k, eta, du, dv = 2, int(sys.argv[1]) if len(sys.argv) > 1 else 12, 10, 4
    nbatch = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    rng = np.random.default_rng(7)
    s2 = eta / 2
    G = offset_dist(eta, dv)
    GA = self_conv(G, 2); GB = self_conv(G, 8)
    BAD, PRED, RM, TT = [], [], [], []
    for it in range(nbatch):
        bad, C = run(k, eta, du, dv, 2000, rng)
        T = np.sum(C ** 2, axis=(-1, -2)); rho = autocorr(C) / T[:, None]
        # ratios for all 240 directions: 112 type A (28 pairs x 4 signs -> +-(e_a +- e_b): ratio 1 +- rho), 128 type B
        RA = np.concatenate([1 + rho[:, [b - a - 1 for a, b in PAIRS]], 1 - rho[:, [b - a - 1 for a, b in PAIRS]]], 1)  # 56 distinct, x2 for sign
        RB = 1 + rho @ WB.T                                                                                            # 64 distinct, x2 for sign
        sd = np.sqrt(2 * s2 * T)[:, None]
        def pfail(Rs, Gd, scale, sd):
            gv, gp = Gd; m = gp > 1e-18; gv = gv[m]; lp = np.log(gp[m])
            z = (Q / 2 - gv[None, None, :] * scale) / (sd[:, :, None] * np.sqrt(Rs)[:, :, None])
            return np.exp(log_ndtr(-z) + lp).sum(-1)
        pblock = np.concatenate([2 * pfail(RA[c:c+100], GA, 1.0, sd[c:c+100]).sum(-1) + 2 * pfail(RB[c:c+100], GB, 0.5, sd[c:c+100]).sum(-1)
                                 for c in range(0, len(T), 100)])
        BAD.append(bad.sum(-1)); PRED.append(32 * pblock); RM.append(np.maximum(RA.max(-1), RB.max(-1))); TT.append(T)
        if it % 10 == 9: print(f"  batch {it+1}/{nbatch}", flush=True)
    BAD = np.concatenate(BAD); PRED = np.concatenate(PRED); RM = np.concatenate(RM); TT = np.concatenate(TT)
    n = len(BAD)
    print(f"\nREAL ENCRYPTIONS k={k}, eta={eta}, du={du}, dv={dv}: {n} ciphertexts, {BAD.sum()} block failures, predicted {PRED.sum():.1f}")
    qs = np.quantile(RM, [0, 0.5, 0.9, 0.99, 1.0])
    for lo, hi, lab in ((qs[0], qs[1], "R_max bottom 50%"), (qs[1], qs[2], "R_max 50-90%"), (qs[2], qs[3], "R_max 90-99%"), (qs[3], qs[4] + 1, "R_max top 1%")):
        m = (RM >= lo) & (RM < hi)
        ob, pr = BAD[m].sum(), PRED[m].sum()
        print(f"  {lab:18s} (R {lo:.3f}-{min(hi, qs[4]):.3f}): ciphertexts {m.sum():6d}  failures observed {ob:5d}  predicted {pr:8.1f}  "
              f"rate per ct {ob / m.sum():.2e} vs {pr / m.sum():.2e}")
    # norm buckets as control
    qt = np.quantile(TT, [0, 0.5, 0.99, 1.0])
    for lo, hi, lab in ((qt[0], qt[1], "T bottom 50%"), (qt[2], qt[3] + 1, "T top 1%")):
        m = (TT >= lo) & (TT < hi)
        print(f"  {lab:18s}: ciphertexts {m.sum():6d}  failures observed {BAD[m].sum():5d}  predicted {PRED[m].sum():8.1f}")
