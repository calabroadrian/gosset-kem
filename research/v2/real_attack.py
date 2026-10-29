"""
AFA-12 v2 — Real lattice attacks on small LWE instances
========================================================
Primal (Kannan embedding) attack with our own LLL/BKZ, comparing:
  Z      : error ~ D_Z (coordinatewise), variance v
  E8-naive: error ~ D_L (L = ConstructionA(H8) = sqrt2*E8), same variance v,
            attacker ignores the structure
  E8-aware: same E8 error, attacker intersects the embedding lattice with
            L^(m/8) x Z  (lattice = 2*L1 + q*L2, valid since gcd(2,q)=1)
  Z-1/sqrt2: Z error with variance v/2 (the predicted equivalent of E8-aware)
Secret uniform mod q (so only the error carries structure).
Success = the reduced basis contains +-(e, 1).
"""
import math
import sys
import time
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, "../v2")
from e8_construction_a import IntegralE8Sampler, ZGaussianSampler, H8, matched_s_for_L


# ---------------------------------------------------------------- integer HNF
def hnf_basis(gens, M):
    """Row-style HNF of an integer generating set, for a lattice containing
    M*Z^d. Entries of non-pivot rows are reduced mod M (keeps coefficients
    small; valid because M*e_i are lattice vectors and are appended)."""
    ncols = len(gens[0])
    A = [[int(v) % M for v in r] for r in gens]
    for i in range(ncols):
        r = [0] * ncols
        r[i] = M
        A.append(r)
    basis = []
    rows = A
    for col in range(ncols):
        rows = [[v % M for v in r] for r in rows]
        rows = [r for r in rows if any(r)]
        # ensure M*e_col remains available below
        r = [0] * ncols
        r[col] = M
        rows.append(r)
        # gcd-eliminate column `col` among rows
        piv = None
        rest = []
        for r in rows:
            if r[col] == 0:
                rest.append(r)
                continue
            if piv is None:
                piv = r
                continue
            # extended gcd combine piv and r
            a, b = piv[col], r[col]
            g, x, y = _egcd(a, b)
            new_piv = [x * p + y * s for p, s in zip(piv, r)]
            new_r = [(b // g) * p - (a // g) * s for p, s in zip(piv, r)]
            piv = new_piv
            if any(new_r):
                rest.append(new_r)
        if piv is not None:
            if piv[col] < 0:
                piv = [-v for v in piv]
            basis.append(piv)
        rows = rest
    return basis


def _egcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        qq, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - qq * x1
        y0, y1 = y1, y0 - qq * y1
    return a, x0, y0


# ---------------------------------------------------------------- LLL / BKZ
def lll(B, delta=0.99):
    B = np.array(B, dtype=np.float64)
    n = B.shape[0]

    def gso(B):
        Q, R = np.linalg.qr(B.T)
        diag = np.diag(R)
        mu = (R / diag[:, None]).T
        return mu, diag ** 2

    mu, Bn = gso(B)
    k = 1
    it = 0
    while k < n:
        it += 1
        for j in range(k - 1, -1, -1):
            r = round(mu[k, j])
            if r:
                B[k] -= r * B[j]
                mu[k, :j + 1] -= r * mu[j, :j + 1]
        if Bn[k] >= (delta - mu[k, k - 1] ** 2) * Bn[k - 1]:
            k += 1
        else:
            B[[k, k - 1]] = B[[k - 1, k]]
            mu, Bn = gso(B)
            k = max(k - 1, 1)
        if it > 200000:
            break
    return np.round(B).astype(np.int64)


def bkz_lite(B, beta, tours=2):
    """Very small BKZ: LLL + for each block, run LLL with deep 'pressure'
    by projecting via local enumeration substitute (pairwise insertion).
    For our small dims we use iterated LLL on rotated blocks (heuristic)."""
    B = lll(B)
    n = B.shape[0]
    for _ in range(tours):
        for i in range(0, n - 1):
            j = min(n, i + beta)
            blk = B[i:j].copy()
            rng = np.random.default_rng(i)
            U = np.eye(len(blk), dtype=np.int64)
            # random unimodular shake then LLL inside block
            for a in range(len(blk) - 1):
                U[a + 1] += rng.integers(-1, 2) * U[a]
            blk = lll(U @ blk)
            B[i:j] = blk
        B = lll(B)
    return B


# ---------------------------------------------------------------- instance
def make_instance(n, m, q, sampler, rng):
    A = rng.integers(0, q, size=(m, n))
    s = rng.integers(0, q, size=n)
    e = sampler.sample(m)
    b = (A @ s + e) % q
    return A, b, e


def embedding_generators(A, b, q):
    """Generators of L1 = {(x,t): x == t*b - A s mod q}."""
    m, n = A.shape
    gens = []
    for i in range(n):
        gens.append(list(A[:, i] % q) + [0])
    for j in range(m):
        r = [0] * (m + 1)
        r[j] = q
        gens.append(r)
    gens.append(list(b % q) + [1])
    r = [0] * (m + 1)
    r[m] = q
    gens.append(r)
    return gens


def L_blocks_generators(m):
    """Generators of L^(m/8) x Z."""
    gens = []
    G = np.array([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 1, 0, 0],
                  [0, 0, 0, 0, 1, 1, 1, 1], [0, 1, 0, 1, 0, 1, 0, 1]])
    for blk in range(m // 8):
        for g in G:
            r = [0] * (m + 1)
            r[8 * blk:8 * blk + 8] = list(g)
            gens.append(r)
        for i in range(8):
            r = [0] * (m + 1)
            r[8 * blk + i] = 2
            gens.append(r)
    r = [0] * (m + 1)
    r[m] = 1
    gens.append(r)
    return gens


def attack(A, b, e, q, aware, beta):
    m = A.shape[0]
    g1 = embedding_generators(A, b, q)
    if aware:
        g2 = L_blocks_generators(m)
        gens = [[2 * v for v in r] for r in g1] + [[q * v for v in r] for r in g2]
        B = hnf_basis(gens, 2 * q)
    else:
        B = hnf_basis(g1, q)
    B = bkz_lite(B, beta) if beta > 2 else lll(B)
    target = np.array(list(e) + [1])
    for row in B:
        if np.array_equal(row, target) or np.array_equal(row, -target):
            return True
    return False


def trial(args):
    kind, n, m, q, var, beta, seed = args
    rng = np.random.default_rng(seed)
    if kind == "Z":
        smp = ZGaussianSampler(math.sqrt(var), seed=seed)
    elif kind == "Z_half":
        smp = ZGaussianSampler(math.sqrt(var / 2), seed=seed)
    else:
        smp = IntegralE8Sampler(S_L_CACHE[var], seed=seed)
    A, b, e = make_instance(n, m, q, smp, rng)
    return attack(A, b, e, q, aware=(kind == "E8_aware"), beta=beta)


S_L_CACHE = {}


def calibrate_Z_var(var):
    """s for Z-sampler such that coord variance == var."""
    lo, hi = 0.05, 20
    for _ in range(60):
        mid = (lo + hi) / 2
        if ZGaussianSampler(mid).coord_variance() < var:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


if __name__ == "__main__":
    n, m, q = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    trials = int(sys.argv[4])
    beta = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    vars_ = [float(x) for x in sys.argv[6].split(",")] if len(sys.argv) > 6 else [1, 2, 4, 8]
    for v in vars_:
        S_L_CACHE[v] = matched_s_for_L(v)
    kinds = ["Z", "E8_naive", "E8_aware", "Z_half"]
    print(f"n={n} m={m} q={q} trials={trials} reduction={'LLL' if beta<=2 else 'BKZ-lite'+str(beta)}")
    print(f"{'var':>6}" + "".join(f"{k:>11}" for k in kinds))
    results = {}
    t0 = time.time()
    with Pool(2) as pool:
        for v in vars_:
            line = f"{v:>6.2f}"
            for k in kinds:
                args = [(k, n, m, q, v, beta, 1000 * hash((k, v)) % 10**6 + i) for i in range(trials)]
                # Z variance: use calibrated s so var is exact
                ok = pool.map(trial, args)
                rate = sum(ok) / trials
                results[(v, k)] = rate
                line += f"{rate*100:>10.0f}%"
            print(line, flush=True)
    print(f"elapsed {time.time()-t0:.0f}s")
    np.save(f"../v2/attack_n{n}_m{m}_q{q}_b{beta}.npy",
            np.array([[v] + [results[(v, k)] for k in kinds] for v in vars_]))
