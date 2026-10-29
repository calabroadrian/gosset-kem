"""AFA-12 v2 — results figure (reads saved experiment outputs)."""
import math
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "../v2")
from estimator import primal_beta

C = {"Z": "#2563eb", "E8_naive": "#9ca3af", "E8_aware": "#dc2626", "Z_half": "#f59e0b",
     "ink": "#111827", "muted": "#6b7280", "grid": "#e5e7eb"}

fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
fig.suptitle("AFA-12 v2 — E8-structured noise at equal variance: what the experiments show",
             fontsize=14, fontweight="bold", color=C["ink"])

# (1) real attacks
ax = axes[0]
d = np.load("../v2/attack_n24_m48_q257_b2.npy")
v = d[:, 0]
labels = {"Z": "Standard noise (D_Z)", "E8_naive": "E8 noise, attacker ignores structure",
          "E8_aware": "E8 noise, structure-aware attacker", "Z_half": "Standard noise, variance / 2"}
for i, k in enumerate(["Z", "E8_naive", "E8_aware", "Z_half"]):
    ls = "--" if k == "Z_half" else "-"
    ax.plot(v, d[:, i + 1] * 100, ls, marker="o", lw=2, color=C[k], label=labels[k])
ax.set_xlabel("per-coordinate noise variance")
ax.set_ylabel("primal attack success (%)")
ax.set_title("Real LLL attacks (n=24, m=48, q=257, 40 trials/point)", fontsize=10)
ax.legend(fontsize=8, frameon=False)
ax.grid(color=C["grid"])

# (2) estimator at cryptographic size
ax = axes[1]
ns = [256, 512, 768, 1024]
lnD = math.log(16) / 8
bz = [primal_beta(n, 3329, 1.0)[0] for n in ns]
be = [primal_beta(n, 3329, 1.0, log_index_per_coord=lnD)[0] for n in ns]
x = np.arange(len(ns))
ax.bar(x - 0.18, [0.292 * b for b in bz], 0.36, color=C["Z"], label="standard noise")
ax.bar(x + 0.18, [0.292 * b for b in be], 0.36, color=C["E8_aware"], label="E8 noise")
for i in range(len(ns)):
    ax.text(x[i] + 0.18, 0.292 * be[i] + 3, f"-{0.292*(bz[i]-be[i]):.0f}", ha="center",
            fontsize=9, color=C["E8_aware"])
ax.set_xticks(x, [f"n={n}" for n in ns])
ax.set_ylabel("core-SVP security (bits, classical)")
ax.set_title("Primal estimate, q=3329, variance 1.0\n(estimator validated vs Kyber: 406/624/874)",
             fontsize=10)
ax.legend(fontsize=8, frameon=False)
ax.grid(axis="y", color=C["grid"])

# (3) entropy
ax = axes[2]
t = np.load("../v2/entropy_table.npy")
ax.plot(t[:, 0], t[:, 3], "-o", color=C["Z"], lw=2, label="standard noise, 8 coords")
ax.plot(t[:, 0], t[:, 4], "-o", color=C["E8_aware"], lw=2, label="E8 noise, 1 block")
ax.set_xscale("log")
ax.set_xlabel("per-coordinate variance")
ax.set_ylabel("entropy per 8-dim block (bits)")
ax.set_title("Entropy at equal variance: exactly log2(det L) = 4 bits lost", fontsize=10)
ax.legend(fontsize=8, frameon=False)
ax.grid(color=C["grid"])

for a in axes:
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("../v2/afa12_v2_results.png", dpi=150)
print("saved afa12_v2_results.png")
