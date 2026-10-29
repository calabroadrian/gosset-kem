"""
AFA-12 v2 — Primal uSVP estimator ("2016 estimate", ADPS16 / Kyber spec)
=========================================================================
Embedding (Bai–Galbraith style, secret and error same std sigma):
  d = n + m + 1, vol = q^m
  BKZ-beta recovers the target iff
        sigma * sqrt(beta) <= delta(beta)^(2*beta - d - 1) * vol^(1/d)
Core-SVP cost: 0.292*beta (classical), 0.265*beta (quantum).

Structured noise: if secret and error lie in L^(blocks) with index D per
8-dim block (L = sqrt(2)E8 integral: D = 16), the attacker intersects the
embedding lattice with that sublattice (q odd, so the mod-2 and mod-q
conditions are independent by CRT): volume is multiplied by D^((n+m)/8)
while the target vector is unchanged.

Validation targets (Kyber round-3 spec, primal core-SVP beta):
  Kyber-512  ~ 406 (eta1=3, sigma=sqrt(1.5)),  Kyber-768 ~ 623,  Kyber-1024 ~ 873
"""
import math


def delta(beta):
    return ((math.pi * beta) ** (1.0 / beta) * beta / (2 * math.pi * math.e)) ** (1.0 / (2 * (beta - 1)))


def primal_beta(n, q, sigma, m_max=None, log_index_per_coord=0.0):
    """Smallest beta for which some m <= m_max succeeds.
    log_index_per_coord: ln(D)/8 per structured coordinate (0 for Z noise)."""
    m_max = m_max or n
    for beta in range(40, 1500):
        lhs = math.log(sigma) + 0.5 * math.log(beta)
        ld = math.log(delta(beta))
        for m in range(max(1, beta - n), m_max + 1):
            d = n + m + 1
            if beta > d:
                continue
            log_vol = m * math.log(q) + (n + m) * log_index_per_coord
            rhs = (2 * beta - d - 1) * ld + log_vol / d
            if lhs <= rhs:
                return beta, m
    return None, None


if __name__ == "__main__":
    print("VALIDATION vs Kyber round-3 published primal beta")
    print(f"{'scheme':<12}{'n':>6}{'sigma':>8}{'beta':>7}{'published':>11}{'bits':>7}")
    refs = [("Kyber-512", 512, math.sqrt(1.5), 406),
            ("Kyber-768", 768, 1.0, 623),
            ("Kyber-1024", 1024, 1.0, 873)]
    for name, n, s, pub in refs:
        b, m = primal_beta(n, 3329, s)
        print(f"{name:<12}{n:>6}{s:>8.3f}{b:>7}{pub:>11}{0.292*b:>7.1f}")

    print("\nEQUAL-VARIANCE COMPARISON  (q=3329, per-coord variance 1.0)")
    lnD = math.log(16) / 8
    print(f"{'n':>6}{'beta Z':>9}{'beta E8':>9}{'delta':>7}{'bits Z':>8}{'bits E8':>9}"
          f"{'  sigma_Z equiv':>16}")
    for n in [256, 512, 768, 1024]:
        bz, _ = primal_beta(n, 3329, 1.0)
        be, _ = primal_beta(n, 3329, 1.0, log_index_per_coord=lnD)
        # find sigma for Z-noise that matches E8 beta
        lo, hi = 0.1, 1.0
        for _ in range(40):
            mid = (lo + hi) / 2
            bm, _ = primal_beta(n, 3329, mid)
            if bm < be:
                lo = mid
            else:
                hi = mid
        print(f"{n:>6}{bz:>9}{be:>9}{be-bz:>7}{0.292*bz:>8.1f}{0.292*be:>9.1f}{hi:>16.3f}")

    print("\nTo keep Kyber-768-level beta with E8 noise you must raise variance:")
    target, _ = primal_beta(768, 3329, 1.0)
    for v in [1.0, 1.5, 2.0, 2.5]:
        b, _ = primal_beta(768, 3329, math.sqrt(v), log_index_per_coord=lnD)
        print(f"  var={v:.1f}: beta={b} (target {target})")
