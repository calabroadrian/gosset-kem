"""Deeper-tail check of the independence assumption (more keys, lower noise)."""
import sys, time, numpy as np
sys.path.insert(0, "../v3")
import correlation_check as cc
from e8_message_coding import kyber_noise
rng = np.random.default_rng(777)
rows = []
for k, eta, nkeys in [(2, 13, 4000), (2, 12, 6000)]:
    t0 = time.time()
    w = kyber_noise(k, eta, eta, 10, 4)
    b_std, b_e8 = cc.bounds(w)
    cnt = {"std": 0, "e8": 0}; nb = 0
    for _ in range(nkeys):
        o = cc.run_real(k, eta, 10, 4, 100, rng)
        cnt["std"] += int(o["std"].sum()); cnt["e8"] += int(o["e8"].sum()); nb += o["std"].size
    i_std, i_e8 = 0, 0
    for _ in range(nb // 1_000_000 + 1):
        a, b = cc.run_iid(w, 1_000_000, rng); i_std += int(a.sum()); i_e8 += int(b.sum())
    ni = (nb // 1_000_000 + 1) * 1_000_000
    for name, cr, ci_, bb in (("std", cnt["std"], i_std, b_std), ("e8", cnt["e8"], i_e8, b_e8)):
        pr, pi = cr / nb, ci_ / ni
        print(f"k={k} eta={eta} {name}: real {pr:.3e} ({cr} ev / {nb})  iid {pi:.3e} ({ci_} ev / {ni})  "
              f"bound {bb:.3e}  real/iid {pr/pi if pi else float('nan'):.2f}", flush=True)
        rows.append([k, eta, 0 if name == "std" else 1, nb, cr, ni, ci_, bb])
    print(f"  [{time.time()-t0:.0f}s]", flush=True)
np.save("../v3/tail_check.npy", np.array(rows, dtype=float))
