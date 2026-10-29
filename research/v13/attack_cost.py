"""AFA-12 v13 — total failure-boosting cost with the Chernoff correlation tail: min work W = 1/(alpha*beta)
subject to Q = 1/beta <= 2^64, for consecutive and distinct-lag layouts (and ML-KEM, norm only)."""
import json, math, sys, numpy as np
from failure_boost import SETS, LN2
from corr_layout import layout_struct, CONSEC, log_beta, wsum
from corr_chernoff import directions, log_alpha_R
names = sys.argv[1:] or [S.name for S in SETS if S.e8]
NEW = json.load(open("layout.json"))["blocks"]
LAY = {"consecutive": layout_struct(CONSEC), "distinct": layout_struct(NEW)}
for S in SETS:
    if S.name not in names: continue
    sr = float(np.sum(S.r[1] * S.r[0] ** 2)); sy = float(np.sum(S.y[1] * S.y[0] ** 2)); sig2 = [(sr, S.k), (sy, S.k)]
    out = []
    for k_, LS in LAY.items():
        D = directions(LS, k_); mult = 14 * 32 if k_ == "consecutive" else 2 * wsum(LS)[1]
        xs = np.linspace(0, 1.2, 49); laR = [0.0] + [log_alpha_R(x, D, sig2) for x in xs[1:]]
        best = (1e9, 0)
        for f in np.linspace(1.0, 3.5, 101):
            T = f * S.ET; laT = S.log_alpha(T)
            for x, lr in zip(xs, laR):
                lb = log_beta(S, T, 1 + x, mult); W = -(laT + lr + lb) / LN2; Qq = -lb / LN2
                if Qq <= 64 and W < best[0]: best = (W, Qq)
        out.append(best)
    c0, _ = S.attack()
    print(f"{S.name:20s} norm only {c0[0]:7.1f} | consecutive {out[0][0]:7.1f} (Q 2^{out[0][1]:.1f}) | distinct {out[1][0]:7.1f} (Q 2^{out[1][1]:.1f})", flush=True)
