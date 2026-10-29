"""
AFA-12 v14 — correlation exposure of the Leech code (24 coordinates per block, 10 blocks in 256 coefficients)
and of E8, consecutive vs optimised layouts.
Exposure metric = sd multiplier of the variance ratio R(u) - 1 for the worst minimal vector u, in units of
sd(rho): sqrt(sum over lag classes of (sum of |weights| of the pairs at that lag)^2), weight of pair (a,b) =
2 |u_a u_b| / |u|^2 (worst-case signs). Leech minimal vectors (scaled by sqrt8, norm^2 32):
  (+-4, +-4, 0^22)            1 pair,   weight 1
  (+-2^8, 0^16) on octads    28 pairs, weight 1/4
  (-+3, +-1^23)              276 pairs: 23 with weight 3/16, 253 with weight 1/16
A random 24-set cannot have 276 distinct lag classes (only 127 exist), so Leech is inherently more exposed.
"""
import itertools, numpy as np, sys
sys.path.insert(0, '../v5')
N = 256
def lagc(a, b): d = (b - a) % N; return min(d, N - d)

def golay_octads():
    H = np.load('../v5/turyn_H.npy') if False else None
    # extended Golay code from the standard generator [I | A] (MacWilliams-Sloane), octads = weight-8 words
    A = np.array([[1,1,0,1,1,1,0,0,0,1,0,1],[1,0,1,1,1,0,0,0,1,0,1,1],[0,1,1,1,0,0,0,1,0,1,1,1],[1,1,1,0,0,0,1,0,1,1,0,1],
                  [1,1,0,0,0,1,0,1,1,0,1,1],[1,0,0,0,1,0,1,1,0,1,1,1],[0,0,0,1,0,1,1,0,1,1,1,1],[0,0,1,0,1,1,0,1,1,1,0,1],
                  [0,1,0,1,1,0,1,1,1,0,0,1],[1,0,1,1,0,1,1,1,0,0,0,1],[0,1,1,0,1,1,1,0,0,0,1,1],[1,1,1,1,1,1,1,1,1,1,1,0]])
    G = np.concatenate([np.eye(12, dtype=int), A], 1)
    octs = []
    for m in range(1, 4096):
        v = np.array([(m >> i) & 1 for i in range(12)]) @ G % 2
        if v.sum() == 8: octs.append(np.nonzero(v)[0])
    return octs

def exposure(pos, octads):
    """pos: 24 (or 8) coefficient positions of one block. Returns worst sd multiplier per vector type."""
    n = len(pos); res = {}
    def agg(pairs_w):
        a = {}
        for (i, j), w in pairs_w: l = lagc(pos[i], pos[j]); a[l] = a.get(l, 0) + w
        return np.sqrt(sum(v * v for k, v in a.items() if k != 128))
    res['pair (4,4)'] = 1.0
    if n == 24:
        res['octad (2^8)'] = max(agg([((i, j), 0.25) for i, j in itertools.combinations(o, 2)]) for o in octads)
        best = 0
        for c in range(24):
            pw = [((i, j), (3 / 16 if c in (i, j) else 1 / 16)) for i, j in itertools.combinations(range(24), 2)]
            best = max(best, agg(pw))
        res['(3,1^23)'] = best
    else:
        res['E8 half (1/2^8)'] = agg([((i, j), 0.25) for i, j in itertools.combinations(range(8), 2)])
    return res

def search_leech(seed, iters=4000):
    rng = np.random.default_rng(seed); octs = golay_octads()
    perm = rng.permutation(N)[:240].reshape(10, 24)
    def cost(P): return max(max(exposure(list(b), octs).values()) for b in P)
    cur = cost(perm)
    for it in range(iters):
        b1, b2 = rng.integers(10, size=2); i, j = rng.integers(24, size=2)
        P = perm.copy(); P[b1, i], P[b2, j] = P[b2, j], P[b1, i]
        c = cost(P[[b1, b2]]) if b1 != b2 else cost(P[[b1]])
        old = cost(perm[[b1, b2]]) if b1 != b2 else cost(perm[[b1]])
        if c <= old: perm = P
    return perm, cost(perm)

if __name__ == "__main__":
    octs = golay_octads(); print(f"Golay octads: {len(octs)} (expected 759)")
    cons = [list(range(24 * b, 24 * b + 24)) for b in range(10)]
    e = exposure(cons[0], octs); print("Leech consecutive block:", {k: round(v, 3) for k, v in e.items()})
    import json
    E8 = json.load(open('../v13/layout.json'))['blocks']
    print("E8 consecutive block:", {k: round(v, 3) for k, v in exposure(list(range(8)), octs).items()})
    print("E8 v13 layout (worst block):", round(max(exposure(b, octs)['E8 half (1/2^8)'] for b in E8), 3))
    best = None
    for seed in range(3):
        P, c = search_leech(seed, int(sys.argv[1]) if len(sys.argv) > 1 else 1500)
        print(f"Leech optimised layout, seed {seed}: worst exposure {c:.3f}", flush=True)
        if best is None or c < best[1]: best = (P, c)
    P = best[0]
    print("Leech best layout per type (worst block):", {k: round(max(exposure(list(b), octs)[k] for b in P), 3) for k in exposure(list(P[0]), octs)})
    json.dump({"blocks": [[int(x) for x in b] for b in P]}, open('leech_layout.json', 'w'))
