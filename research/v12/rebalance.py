"""
AFA-12 v12 — re-balance MaxSec: maximize the weakest attack over BOTH LWE instances
(key: eta1/eta1; ciphertext: eta1/eta2), with ciphertext size <= FIPS and DFR (E8 code) <= FIPS DFR.
"""
import sys, math
sys.path.insert(0, '../v3')
from e8_message_coding import kyber_noise, failure_standard, failure_e8, log2
from coresvp import primal, dual, CBD

q = 3329
FIPS = {2: (3, 2, 10, 4), 3: (2, 2, 10, 4), 4: (2, 2, 11, 5)}

def sec(k, e1, e2, du, rounding):
    n = 256 * k; s1 = CBD(e1)
    se = math.sqrt(e2 / 2 + (q / 2 ** du) ** 2 / 12) if rounding else CBD(e2)
    kp = 0.292 * primal(n, q, s1, s1, n)[0]; kd = dual(n, q, s1, s1, n)[0]
    cp = 0.292 * primal(n, q, s1, se, n)[0]; cd = dual(n, q, s1, se, n)[0]
    return min(kp, kd, cp, cd)

for k in (2, 3, 4):
    f = FIPS[k]
    ctmax = 32 * (k * f[2] + f[3])
    dfr_std = log2(failure_standard(kyber_noise(k, *f), q))
    base_c, base_r = sec(k, f[0], f[1], f[2], False), sec(k, f[0], f[1], f[2], True)
    print(f"k={k} ML-KEM: ct {ctmax} B, DFR 2^{dfr_std:.0f}, security conservative {base_c:.1f}, with rounding {base_r:.1f}", flush=True)
    rows = []
    for e1 in (2, 3, 4, 5):
        for e2 in (2, 3, 4):
            for du in (9, 10, 11):
                for dv in (3, 4, 5):
                    ct = 32 * (k * du + dv)
                    if ct > ctmax: continue
                    d = log2(failure_e8(kyber_noise(k, e1, e2, du, dv), q)[0])
                    if d > dfr_std: continue
                    c, r = sec(k, e1, e2, du, False), sec(k, e1, e2, du, True)
                    rows.append((c, r, e1, e2, du, dv, ct, d))
    rows.sort(reverse=True)
    for c, r, e1, e2, du, dv, ct, d in rows[:5]:
        print(f"   eta1={e1} eta2={e2} du={du} dv={dv}  ct {ct} B  DFR 2^{d:.0f}  conservative {c:.1f} ({c-base_c:+.1f})  with rounding {r:.1f} ({r-base_r:+.1f})", flush=True)
