"""
AFA-12 v2 — Exact discrete Gaussian over the integral E8 lattice
=================================================================
L = Construction A of the extended Hamming code H8 [8,4,4]:
    L = { x in Z^8 : x mod 2 in H8 }  ==  sqrt(2) * E8   (det 16, min norm^2 4)

This is the correct way to get *integer* E8-structured noise (no rounding).
The original harness rounded E8 points with numpy's round-half-to-even,
which mapped every half-integer E8 point to an all-even vector: the noise
used in phases 2-4 was not E8 noise.

Exact sampler: D_{L,s} is a mixture over the 16 cosets c in H8 of
D_{2Z^8 + c, s}; given the coset, coordinates are independent, each from
D_{2Z + c_i, s}. Coset weight = prod_i rho_s(2Z + c_i).
Gaussian convention: rho_s(x) = exp(-x^2 / (2 s^2)).
"""
import itertools
import math
import numpy as np


def hamming_h8():
    """Extended Hamming [8,4,4] code: 16 codewords, weights {0,4,8}."""
    G = np.array([
        [1, 1, 1, 1, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 1, 1, 1],
        [0, 1, 0, 1, 0, 1, 0, 1],
    ], dtype=np.int64)
    words = []
    for bits in itertools.product([0, 1], repeat=4):
        words.append((np.array(bits) @ G) % 2)
    W = np.array(words, dtype=np.int64)
    assert len({tuple(w) for w in W}) == 16
    assert set(W.sum(1).tolist()) == {0, 4, 8}
    return W


H8 = hamming_h8()
TAIL = 14  # support radius in units of s


def _support(s):
    R = int(math.ceil(TAIL * s)) + 2
    return np.arange(-R, R + 1)


def coord_tables(s):
    """Per-parity 1D distributions D_{2Z+b, s} for b in {0,1}."""
    x = _support(s)
    rho = np.exp(-(x.astype(float) ** 2) / (2 * s * s))
    tabs = {}
    for b in (0, 1):
        mask = (x % 2) == b
        xs, ps = x[mask], rho[mask]
        tabs[b] = (xs, ps / ps.sum(), ps.sum())
    return tabs


class IntegralE8Sampler:
    """Exact D_{L, s} sampler for L = ConstructionA(H8) ~ sqrt(2) E8."""

    def __init__(self, s, seed=0):
        self.s = s
        self.rng = np.random.default_rng(seed)
        self.tabs = coord_tables(s)
        m0, m1 = self.tabs[0][2], self.tabs[1][2]
        w = H8.sum(1)
        logw = (8 - w) * math.log(m0) + w * math.log(m1)
        p = np.exp(logw - logw.max())
        self.coset_p = p / p.sum()
        self.cdf = {b: np.cumsum(self.tabs[b][1]) for b in (0, 1)}

    def sample_blocks(self, k):
        cos = self.rng.choice(16, size=k, p=self.coset_p)
        bits = H8[cos]                       # (k, 8)
        out = np.empty((k, 8), dtype=np.int64)
        u = self.rng.random((k, 8))
        for b in (0, 1):
            xs = self.tabs[b][0]
            idx = np.searchsorted(self.cdf[b], u[bits == b], side="right")
            idx = np.minimum(idx, len(xs) - 1)
            out[bits == b] = xs[idx]
        return out

    def sample(self, n):
        assert n % 8 == 0
        return self.sample_blocks(n // 8).reshape(-1)

    # ---- exact analytic quantities ----
    def coord_variance(self):
        var = 0.0
        for c, pc in zip(H8, self.coset_p):
            for b in c:
                xs, ps, _ = self.tabs[int(b)]
                var += pc * float((ps * xs.astype(float) ** 2).sum())
        return var / 8

    def block_entropy_bits(self):
        H = -float((self.coset_p * np.log2(self.coset_p)).sum())
        hb = {}
        for b in (0, 1):
            _, ps, _ = self.tabs[b]
            ps = ps[ps > 0]
            hb[b] = -float((ps * np.log2(ps)).sum())
        for c, pc in zip(H8, self.coset_p):
            H += pc * sum(hb[int(b)] for b in c)
        return H


class ZGaussianSampler:
    """Exact D_{Z, s} coordinatewise sampler (reference: standard LWE)."""

    def __init__(self, s, seed=0):
        self.s = s
        self.rng = np.random.default_rng(seed)
        x = _support(s)
        rho = np.exp(-(x.astype(float) ** 2) / (2 * s * s))
        self.xs, self.ps = x, rho / rho.sum()
        self.cdf = np.cumsum(self.ps)

    def sample(self, n):
        idx = np.searchsorted(self.cdf, self.rng.random(n), side="right")
        return self.xs[np.minimum(idx, len(self.xs) - 1)]

    def coord_variance(self):
        return float((self.ps * self.xs.astype(float) ** 2).sum())

    def block_entropy_bits(self):
        p = self.ps[self.ps > 0]
        return 8 * -float((p * np.log2(p)).sum())


def matched_s_for_L(target_var):
    """Find s such that Var(D_{L,s}) per coordinate == target_var (bisection)."""
    lo, hi = 0.2, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if IntegralE8Sampler(mid).coord_variance() < target_var:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def membership_ok(blocks):
    """x in L  <=>  x mod 2 is an H8 codeword."""
    codes = {tuple(w) for w in H8}
    return np.array([tuple(r) in codes for r in (blocks % 2)])


if __name__ == "__main__":
    print("=" * 70)
    print("Integral E8 (Construction A of Hamming [8,4,4]) — exact analysis")
    print("=" * 70)

    # Lattice facts
    B = np.vstack([H8[[1, 2, 4, 8]], 2 * np.eye(8, dtype=np.int64)])
    print(f"H8 weights: {sorted(set(H8.sum(1).tolist()))} (doubly-even)")

    smp = IntegralE8Sampler(1.5, seed=1)
    blocks = smp.sample_blocks(20000)
    print(f"Membership check (20000 blocks): {membership_ok(blocks).mean()*100:.2f}% in L")
    print(f"Sum parity even: {(blocks.sum(1) % 2 == 0).mean()*100:.2f}% "
          f"(structural: Hamming words have even weight)")
    print(f"Coordinate mod 2 pattern known: always one of 16 codewords "
          f"-> 4 of 8 parity bits determined per block")

    print("\nEqual-variance comparison, per 8D block:")
    print(f"{'Var/coord':>10} {'s_Z':>7} {'s_L':>7} {'H_Z (bits)':>11} "
          f"{'H_L (bits)':>11} {'loss':>7}")
    rows = []
    for target in [0.5, 1.0, 1.5, 2.0, 4.0, 9.0]:
        # Z sampler: find s with variance target
        lo, hi = 0.1, 20.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if ZGaussianSampler(mid).coord_variance() < target:
                lo = mid
            else:
                hi = mid
        sZ = 0.5 * (lo + hi)
        sL = matched_s_for_L(target)
        HZ = ZGaussianSampler(sZ).block_entropy_bits()
        HL = IntegralE8Sampler(sL).block_entropy_bits()
        rows.append((target, sZ, sL, HZ, HL, HZ - HL))
        print(f"{target:>10.2f} {sZ:>7.3f} {sL:>7.3f} {HZ:>11.3f} {HL:>11.3f} {HZ-HL:>7.3f}")

    print("\nTheory: for s >> 1, entropy loss -> log2(det L) = log2(16) = 4 bits/block,")
    print("i.e. 0.5 bit per coordinate; n=256 -> ~128 bits less entropy than D_Z.")
    np.save("../v2/entropy_table.npy", np.array(rows))
