"""AFA-12 v13 — DFR under adversarial selection: best per-ciphertext failure probability an attacker reaches
when grinding 2^X messages (alpha >= 2^-X), using norm (all schemes) and in-block correlation (E8)."""
import math, numpy as np
from failure_boost import SETS, LN2
from corr_boost import log_alpha_R, log_beta_dir
SD = {"AFA12-512-Compact": 0.0350, "AFA12-512-MaxSec": 0.0314, "AFA12-768-Compact": 0.0303,
      "AFA12-768-MaxSec": 0.0257, "AFA12-1024-Compact": 0.0232, "AFA12-1024-MaxSec": 0.0221}   # from corr_boost MC
print(f"{'scheme':20s} {'avg-ct beta':>11s} | best log2 beta with grinding 2^32 / 2^64 / 2^128  (norm only)  | with correlation (E8)")
for S in SETS:
    row_n, row_c = [], []
    for X in (32, 64, 128):
        bn = -1e9; bc = -1e9
        for f in np.linspace(1.0, 3.5, 251):
            T = f * S.ET; la_T = S.log_alpha(T)
            if -la_T / LN2 <= X: bn = max(bn, S.log_beta(T) / LN2)
            if S.e8:
                v = SD[S.name] ** 2
                for xs in np.linspace(0, 16, 65):
                    la = la_T + (log_alpha_R(xs * math.sqrt(v), v) if xs > 0 else 0)
                    if -min(la, 0) / LN2 <= X: bc = max(bc, log_beta_dir(S, T, 1 + xs * math.sqrt(v)) / LN2)
        row_n.append(bn); row_c.append(bc)
    b0 = S.log_beta(S.ET) / LN2
    cs = " / ".join(f"{x:7.1f}" for x in row_c) if S.e8 else "      (n/a)"
    print(f"{S.name:20s} {b0:11.1f} | " + " / ".join(f"{x:7.1f}" for x in row_n) + f"            | {cs}", flush=True)
