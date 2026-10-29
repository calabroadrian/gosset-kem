"""
AFA-12 v4 — correctness of the C codec against the Python reference.
1. multiply-shift rounding exact over the whole input range
2. encode/decode round trip for all 256 messages per block
3. noisy decoding: C decoder == Python nearest-E8 decoder (millions of blocks),
   any disagreement must be an exact distance tie
4. failure rate of C decoder on the real-scheme noise == Python figures
"""
import ctypes
import sys
import numpy as np

sys.path.insert(0, "../v3")
import correlation_check as cc

lib = ctypes.CDLL("../v4/libe8codec.so")
lib.e8_round_div_U_test.restype = ctypes.c_int64
lib.e8_round_div_U_test.argtypes = [ctypes.c_int64]
Q, U = 3329, 2 * 3329
I16P = np.ctypeslib.ndpointer(np.int16, flags="C_CONTIGUOUS")
U8P = np.ctypeslib.ndpointer(np.uint8, flags="C_CONTIGUOUS")
lib.e8_encode.argtypes = [U8P, I16P]
lib.e8_decode.argtypes = [I16P, U8P]
lib.e8_decode_leaky.argtypes = [I16P, U8P]

# 1. rounding
bad = 0
for Y in range(-8 * U, 8 * U + 1):
    exact = (2 * Y + U) // (2 * U)          # round half up
    if lib.e8_round_div_U_test(Y) != exact:
        bad += 1
print(f"[1] multiply-shift rounding: {bad} mismatches over {16*U+1} inputs")

# 2. round trip
msgs = np.arange(256, dtype=np.uint8).reshape(8, 32)
ok = True
for row in msgs:
    enc = np.zeros(256, np.int16)
    lib.e8_encode(np.ascontiguousarray(row), enc)
    for dec_fn in (lib.e8_decode, lib.e8_decode_leaky):
        out = np.zeros(32, np.uint8)
        dec_fn(enc, out)
        ok &= np.array_equal(out, row)
print(f"[2] encode/decode round trip, all 256 messages: {'OK' if ok else 'FAIL'}")


def c_index_to_rep(m):
    """C message index -> E8 coset representative mod 2 (8 coords)."""
    bits = [(m >> i) & 1 for i in range(7)]
    bits.append(sum(bits) % 2)
    return np.array(bits, float) + (0.5 if m >> 7 else 0.0)


REP_OF = np.array([c_index_to_rep(m) for m in range(256)])


def dist2_to_coset(y, rep):
    """squared distance from y (8,) to nearest point of (q/2)*(rep + 2Z^8)."""
    z = y / (Q / 2) - rep
    z = z - 2 * np.round(z / 2)
    return float((z * z).sum())


# 3. noisy decoding vs Python
rng = np.random.default_rng(5)
blocks_total, mism, ties = 0, 0, 0
for trial in range(200):
    msg = rng.integers(0, 256, 32).astype(np.uint8)
    enc = np.zeros(256, np.int16)
    lib.e8_encode(msg, enc)
    scale = rng.choice([50, 200, 400, 600, 900])
    noise = np.rint(rng.normal(0, scale, 256)).astype(np.int64)
    y = ((enc.astype(np.int64) + noise) % Q).astype(np.int16)
    for _ in range(50):
        pass
    out = np.zeros(32, np.uint8)
    lib.e8_decode(y, out)
    py = cc.e8_decode_rep(cc.center(y.astype(np.int64)).reshape(32, 8).astype(float))
    for b in range(32):
        blocks_total += 1
        crep = REP_OF[out[b]]
        if not np.allclose(np.mod(crep, 2), py[b]):
            mism += 1
            yb = cc.center(y.astype(np.int64)).reshape(32, 8)[b].astype(float)
            if abs(dist2_to_coset(yb, crep) - dist2_to_coset(yb, py[b])) < 1e-9:
                ties += 1
print(f"[3] noisy decode: {blocks_total} blocks, {mism} disagreements, {ties} of them exact ties")

# larger randomized sweep (vectorized)
N = 2_000_000
y = rng.integers(0, Q, size=(N // 32, 256)).astype(np.int16)
outs = np.zeros((N // 32, 32), np.uint8)
for i in range(len(y)):
    lib.e8_decode(np.ascontiguousarray(y[i]), outs[i])
py = cc.e8_decode_rep(cc.center(y.astype(np.int64)).reshape(-1, 8).astype(float))
crep = np.mod(REP_OF[outs.reshape(-1)], 2)
dis = np.nonzero(np.abs(crep - py).sum(-1) > 1e-9)[0]
yy = cc.center(y.astype(np.int64)).reshape(-1, 8).astype(float)
nt = sum(abs(dist2_to_coset(yy[j], REP_OF[outs.reshape(-1)[j]]) - dist2_to_coset(yy[j], py[j])) < 1e-9
         for j in dis)
print(f"[3b] uniform random inputs: {N//8*8//8} blocks... "
      f"{len(py)} blocks, {len(dis)} disagreements, {nt} exact ties")

# 4. failure rate on real-scheme noise (k=2, eta=13) with the C decoder
rng = np.random.default_rng(2027)
fails, blocks = 0, 0
for _ in range(1500):
    k, eta, du, dv, B = 2, 13, 10, 4, 100
    A = rng.integers(0, Q, size=(k, k, 256))
    s = cc.cbd(rng, eta, (k, 256)); e = cc.cbd(rng, eta, (k, 256))
    t = (sum(cc.polymul(A[:, j], s[j][None]) for j in range(k)) + e) % Q
    r = cc.cbd(rng, eta, (B, k, 256)); e1 = cc.cbd(rng, eta, (B, k, 256)); e2 = cc.cbd(rng, eta, (B, 256))
    u = (np.stack([sum(cc.polymul(A[i, j][None], r[:, i]) for i in range(k)) for j in range(k)], 1) + e1) % Q
    tr = sum(cc.polymul(t[i][None], r[:, i]) for i in range(k))
    msg = rng.integers(0, 256, size=(B, 32)).astype(np.uint8)
    enc = np.zeros((B, 256), np.int16)
    for i in range(B):
        lib.e8_encode(msg[i], enc[i])
    v = (tr + e2 + enc) % Q
    uc = cc.decompress(cc.compress(u, du), du); vc = cc.decompress(cc.compress(v, dv), dv)
    yv = ((vc - sum(cc.polymul(s[i][None], uc[:, i]) for i in range(k))) % Q).astype(np.int16)
    out = np.zeros(32, np.uint8)
    for i in range(B):
        lib.e8_decode(np.ascontiguousarray(yv[i]), out)
        fails += int((out != msg[i]).sum()); blocks += 32
print(f"[4] C decoder on real ML-KEM-style noise (k=2, eta=13): block failure "
      f"{fails/blocks:.3e} ({fails} / {blocks});  v3 Python measurement 1.711e-05, bound 1.705e-05")
