"""AFA-12 v13 — the new layout must not change the average failure rate: real ring encryptions at inflated
noise (k=2, eta=13), same ciphertexts decoded with blocks taken consecutively vs with E8_LAYOUT."""
import sys, json, numpy as np
sys.path.insert(0, '../v3')
from correlation_check import polymul, cbd, compress, decompress, center, e8_decode_rep, REPS, HALF, Q, N
LAY = np.array([x for b in json.load(open('layout.json'))['blocks'] for x in sorted(b)])
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
k, eta, du, dv = 2, 13, 10, 4
fails = {"consecutive": 0, "layout": 0}; blocks = 0
for it in range(int(sys.argv[1]) if len(sys.argv) > 1 else 50):
    B = 2000
    A = rng.integers(0, Q, size=(k, k, N)); s = cbd(rng, eta, (k, N)); e = cbd(rng, eta, (k, N))
    t = (sum(polymul(A[:, j], s[j][None]) for j in range(k)) + e) % Q
    r = cbd(rng, eta, (B, k, N)); e1 = cbd(rng, eta, (B, k, N)); e2 = cbd(rng, eta, (B, N))
    u = (np.stack([sum(polymul(A[i, j][None], r[:, i]) for i in range(k)) for j in range(k)], 1) + e1) % Q
    tr = sum(polymul(t[i][None], r[:, i]) for i in range(k))
    msg = rng.integers(0, 256, size=(B, N // 8)); pts = np.rint(REPS[msg] * HALF).astype(np.int64).reshape(B, N)
    for name in fails:
        enc = pts if name == "consecutive" else np.empty_like(pts)
        if name == "layout": enc[:, LAY] = pts
        v = (tr + e2 + enc) % Q
        uc = decompress(compress(u, du), du); vc = decompress(compress(v, dv), dv)
        y = center(vc - sum(polymul(s[i][None], uc[:, i]) for i in range(k)))
        yy = y if name == "consecutive" else y[:, LAY]
        rep = e8_decode_rep(center(yy).reshape(B, N // 8, 8).astype(float))
        fails[name] += int((np.abs(rep - REPS[msg]).sum(-1) > 1e-9).sum())
    blocks += B * 32
for nm, f in fails.items():
    p = f / blocks; print(f"{nm:12s} block failure {p:.3e} ({f}/{blocks}), 95% CI +-{1.96 * np.sqrt(p * (1 - p) / blocks):.1e}")
print("v3 prediction (independence model, eta=13): 1.705e-05")
