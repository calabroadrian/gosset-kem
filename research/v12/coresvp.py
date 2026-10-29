"""
AFA-12 v12 — core-SVP security estimates (primal uSVP and dual), Kyber/NIST methodology.
Cost model: classical 0.292*beta, quantum 0.265*beta (core-SVP, as in the Kyber spec).
Two LWE instances per parameter set:
  KEY : secret s ~ CBD(eta1), error e ~ CBD(eta1), n = k*256, m <= k*256 samples
  CT  : secret r ~ CBD(eta1), error e1 ~ CBD(eta2) (+ optional du rounding), m <= k*256
Validation: Kyber Round-3 spec, Table 4 (core-SVP classical primal/dual: 118/117, 182/181, 256/253).
"""
import math, sys

def delta(b):
    return ((math.pi * b) ** (1.0 / b) * b / (2 * math.pi * math.e)) ** (1.0 / (2 * (b - 1)))

def primal(n, q, ss, se, mmax):
    """min beta s.t. se*sqrt(beta) <= delta^(2beta-d) * vol^(1/d), d=n+m+1, secret rescaled to se."""
    best = None
    for b in range(60, 1600):
        lg_d = math.log(delta(b))
        ok = False
        for m in range(max(1, n // 4), mmax + 1, 4):
            d = n + m + 1
            lvol = m * math.log(q) + n * math.log(se / ss)
            if math.log(se) + 0.5 * math.log(b) <= (2 * b - d) * lg_d + lvol / d:
                ok = True; mm = m; break
        if ok:
            return b, mm
    return None

def dual(n, q, ss, se, mmax):
    """min over (beta, m) of 0.292 beta + max(0, log2(1/eps^2) - 0.2075 beta)  (Kyber spec dual)."""
    best = (1e9, None, None)
    for b in range(60, 1600, 1):
        lg_d = math.log(delta(b))
        for m in range(max(1, n // 4), mmax + 1, 8):
            d = n + m
            # dual lattice {(x,y): x^T A = y^T mod q}, vol q^n; y rescaled by ss/se so that
            # <x,e> + <y,s> has std se*|(x,y')|; shortest vector found ~ delta^(d-1) * vol^(1/d)
            l = math.exp((d - 1) * lg_d + n * math.log(q * ss / se) / d)
            tau = l * se / q
            eps = 4 * math.exp(-2 * math.pi ** 2 * tau ** 2)
            if eps <= 0: continue
            cost = 0.292 * b + max(0.0, -2 * math.log2(eps) - 0.2075 * b) if eps < 1 else 0.292 * b
            if cost < best[0]:
                best = (cost, b, m)
        if 0.292 * b > best[0]:
            break
    return best

CBD = lambda eta: math.sqrt(eta / 2.0)

def evaluate(name, k, eta1, eta2, du, dv, q=3329, rounding=False):
    n = 256 * k
    s1 = CBD(eta1)
    out = {}
    pb, _ = primal(n, q, s1, s1, n)
    dc, db, _ = dual(n, q, s1, s1, n)
    out['key'] = (pb, dc, db)
    se = CBD(eta2)
    if rounding:
        se = math.sqrt(eta2 / 2.0 + (q / 2 ** du) ** 2 / 12.0)
    pb2, _ = primal(n, q, s1, se, n)
    dc2, db2, _ = dual(n, q, s1, se, n)
    out['ct'] = (pb2, dc2, db2)
    return out

if __name__ == '__main__':
    sets = [
        ("ML-KEM-512", 2, 3, 2, 10, 4), ("ML-KEM-768", 3, 2, 2, 10, 4), ("ML-KEM-1024", 4, 2, 2, 11, 5),
        ("AFA12-512-Compact", 2, 3, 2, 9, 4), ("AFA12-768-Compact", 3, 2, 2, 9, 5), ("AFA12-1024-Compact", 4, 2, 2, 10, 4),
        ("AFA12-512-MaxSec", 2, 4, 3, 10, 4), ("AFA12-768-MaxSec", 3, 3, 2, 10, 4), ("AFA12-1024-MaxSec", 4, 3, 2, 11, 5),
    ]
    only = sys.argv[1:]
    for s in sets:
        if only and s[0] not in only: continue
        r = evaluate(*s)
        rr = evaluate(*s, rounding=True)
        kp, kd, kdb = r['key']; cp, cd, cdb = r['ct']; rp, rd, _ = rr['ct']
        print(f"{s[0]:20s} KEY primal b={kp} ({0.292*kp:.1f}/{0.265*kp:.1f})  dual b={kdb} ({kd:.1f})  | "
              f"CT primal b={cp} ({0.292*cp:.1f}) dual ({cd:.1f}) | CT+rounding primal b={rp} ({0.292*rp:.1f}) dual ({rd:.1f})", flush=True)
