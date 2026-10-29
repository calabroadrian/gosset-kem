"""
AFA-12 v14 — directional failure boosting (D'Anvers, Rossi, Virdia, EUROCRYPT 2020), geometric model,
for the standard decoder and the E8 decoder.
After n failures the attacker holds an estimate E of the secret direction with cos(E, S) = rho_n.
He then grinds ciphertexts whose decryption-error vector for some decoding direction is aligned with E.

Per decoding direction:  standard: <S, C_i> + G_i > q/4      (norm of C_i = c, offset variance v_G, M = 256)
                         E8:       <S, C_u> + <G,u> > q/2    (norm of C_u = c*sqrt2, offset variance 2 v_G, M = 32*240)
Grinding: alpha_y = min(1, M * Q(y sqrt(d)))  (cos of a random direction with E ~ N(0, 1/d)), times the norm
grinding alpha_T. Given alignment y: <S, C_dir> = s c k [y rho + sqrt((1-y^2)(1-rho^2)) Z / sqrt d].
Each failure yields a vector with cos t = thr / (s c k) to S; n failures: rho_n = sqrt(n) t / sqrt(1 + (n-1) t^2).
Simplifications favour the attacker: the failing coefficient / direction is assumed known (no rotation search),
and every failure adds full information.
"""
import sys, math, numpy as np
sys.path.insert(0, '../v13')
from scipy.special import log_ndtr
from failure_boost import SETS, LN2, q

def gvar(S):
    g, p = S.G; return float(np.sum(p * g ** 2))

def next_failure_cost(S, rho, qlog=64):
    d = 2 * S.nk; s = math.sqrt(d * S.s2)
    if S.e8: thr, kap, M = q / 2, math.sqrt(2), 32 * 240
    else: thr, kap, M = q / 4, 1.0, 256
    vG = gvar(S) * kap ** 2
    best = (1e9, 0, 0, 0)
    for f in np.linspace(1.0, 3.0, 81):
        T = f * S.ET; c = math.sqrt(T); laT = S.log_alpha(T)
        for y in np.linspace(0, 0.9, 181):
            la_y = 0.0 if y == 0 else min(0.0, math.log(M) + log_ndtr(-y * math.sqrt(d)))
            mean = s * c * kap * y * rho
            sd = math.sqrt(s * s * c * c * kap * kap * (1 - y * y) * (1 - rho * rho) / d + vG)
            lb = log_ndtr(-(thr - mean) / sd)
            W = -(laT + la_y + lb) / LN2; Qq = -lb / LN2
            if Qq <= qlog and W < best[0]: best = (W, Qq, f, y)
    return best

def rho_after(S, n):
    d = 2 * S.nk; s = math.sqrt(d * S.s2); c = math.sqrt(S.ET)
    thr, kap = (q / 2, math.sqrt(2)) if S.e8 else (q / 4, 1.0)
    t = min(0.99, thr / (s * c * kap))
    return math.sqrt(n) * t / math.sqrt(1 + (n - 1) * t * t) if n > 0 else 0.0, t

if __name__ == "__main__":
    print("DIRECTIONAL FAILURE BOOSTING: log2 work for failure number n+1 given n failures (each with <= 2^64 queries)")
    print(f"{'scheme':20s} {'t (cos per failure)':>19s} | " + " | ".join(f"n={n} (rho)" for n in range(0, 6)))
    for S in SETS:
        cells = []; _, t = rho_after(S, 1)
        for n in range(0, 6):
            rho, _ = rho_after(S, n); W = next_failure_cost(S, rho)
            cells.append(f"{W[0]:6.1f} ({rho:.2f})")
        print(f"{S.name:20s} {t:19.3f} | " + " | ".join(cells), flush=True)
