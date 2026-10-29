"""
AFA-12 v5 — correctness of the C Leech codec.
1. round trip for random 36-bit messages
2. C decoder vs independent brute-force nearest point (numpy, all 8192 cosets)
   on noisy and uniform inputs: same distance (ties allowed to differ in label)
3. lattice sanity: min distance / kissing structure through the encoder
4. block failure on REAL ML-KEM-style noise vs exact union-bound prediction
"""
import ctypes, math, sys, time
import numpy as np
sys.path.insert(0, "../v3"); sys.path.insert(0, "../v5")
import correlation_check as cc
from e8_message_coding import kyber_noise, tail, conv, power
from leech_eval import dilate

Q = 3329
lib = ctypes.CDLL("../v5/libleech.so")
I16P = np.ctypeslib.ndpointer(np.int16, flags="C_CONTIGUOUS")
lib.leech_encode_block.argtypes = [ctypes.c_uint64, I16P]
lib.leech_decode_block.argtypes = [I16P]
lib.leech_decode_block.restype = ctypes.c_uint64
H, HP = np.load("../v5/turyn_H.npy")


def golay():
    cw = []
    for x in HP:
        for a in H:
            for b in H:
                w = (int(a ^ x) << 16) | (int(b ^ x) << 8) | int(a ^ b ^ x)
                cw.append([(w >> (23 - i)) & 1 for i in range(24)])
    return np.array(cw)


C = golay()                                        # (4096, 24)


def enc(m):
    out = np.zeros(24, np.int16); lib.leech_encode_block(ctypes.c_uint64(int(m)), out); return out


def point_of_msg(m):
    """lattice point p (0..7 per coord, mod 8) of a message, via the encoder rounding inverse"""
    y = enc(m).astype(np.int64)
    return np.rint(8 * y / Q).astype(np.int64) % 8


def brute_nearest_dist(y):
    """min over Lambda24 points of ||8y/q - p||^2 (p in lattice, modulo 8Z^24)."""
    Y = 8.0 * cc.center(y.astype(np.int64)) / Q            # (24,)
    best = np.inf
    for h in (0, 1):
        T = h + 2 * C                                       # (4096,24) residues mod 4
        t = (Y[None, :] - T) / 4.0
        k = np.floor(t + 0.5)
        e = (t - k) * 4.0
        d = (e ** 2).sum(1)
        mis = ((k.sum(1) - h) % 2) != 0
        s = np.where(e >= 0, 1.0, -1.0)
        fl = ((e - 4 * s) ** 2 - e ** 2).min(1)
        d = d + np.where(mis, fl, 0.0)
        best = min(best, d.min())
    return best


def dist_to_msg(y, m):
    Y = 8.0 * cc.center(y.astype(np.int64)) / Q
    p = point_of_msg(m)
    diff = Y - p
    diff -= 8 * np.round(diff / 8)
    return float((diff ** 2).sum())


rng = np.random.default_rng(9)

# 1. round trip
ok = all(lib.leech_decode_block(enc(m)) == m for m in rng.integers(0, 2 ** 36, 3000, dtype=np.uint64))
print(f"[1] round trip, 3000 random 36-bit messages: {'OK' if ok else 'FAIL'}")

# 3. lattice sanity: minimal squared distance between distinct points (sample) = 32
pts = np.array([point_of_msg(m) for m in rng.integers(0, 2 ** 36, 400, dtype=np.uint64)])
p0 = point_of_msg(0)
dmin = min(float(((((p - p0 + 4) % 8) - 4) ** 2).sum()) for p in pts if not np.array_equal(p, p0))
print(f"[3] sampled min squared distance to origin point: {dmin} (Leech min norm 32 expected >= 32)")

# 2. C vs brute force
t0 = time.time(); n = 0; bad = 0; ties = 0
for trial in range(1500):
    if trial % 2:
        y = rng.integers(0, Q, 24).astype(np.int16)
    else:
        m = int(rng.integers(0, 2 ** 36, dtype=np.uint64))
        y = ((enc(m).astype(np.int64) + np.rint(rng.normal(0, rng.choice([150, 300, 450]), 24)).astype(np.int64)) % Q).astype(np.int16)
    mc = lib.leech_decode_block(np.ascontiguousarray(y))
    dc, db = dist_to_msg(y, mc), brute_nearest_dist(y)
    n += 1
    if abs(dc - db) > 1e-6:
        bad += 1
print(f"[2] C vs brute-force nearest point: {n} blocks, {bad} with larger distance than optimum "
      f"[{time.time()-t0:.0f}s]")

# 4. real-noise block failure
k, eta = 2, 11
w = kyber_noise(k, eta, eta, 10, 4)
p2 = tail(conv(w, w), Q / 2 - 1); p8 = tail(power(w, 8), Q - 4)
w23 = conv(power(w, 16), conv(power(w, 4), conv(power(w, 2), w)))
p24 = tail(conv(dilate(w, 3), w23), 2 * Q - 13)
pred = 1104 * p2 + 97152 * p8 + 98304 * p24
fails = blocks = 0
t0 = time.time()
for _ in range(1500):
    A = rng.integers(0, Q, size=(k, k, 256))
    s = cc.cbd(rng, eta, (k, 256)); e = cc.cbd(rng, eta, (k, 256))
    t = (sum(cc.polymul(A[:, j], s[j][None]) for j in range(k)) + e) % Q
    B = 50
    r = cc.cbd(rng, eta, (B, k, 256)); e1 = cc.cbd(rng, eta, (B, k, 256)); e2 = cc.cbd(rng, eta, (B, 256))
    u = (np.stack([sum(cc.polymul(A[i, j][None], r[:, i]) for i in range(k)) for j in range(k)], 1) + e1) % Q
    msgs = rng.integers(0, 2 ** 36, size=(B, 8), dtype=np.uint64)
    encv = np.zeros((B, 256), np.int64)
    for bi in range(B):
        for j in range(8):
            encv[bi, 24 * j:24 * j + 24] = enc(msgs[bi, j])
    v = (sum(cc.polymul(t[i][None], r[:, i]) for i in range(k)) + e2 + encv) % Q
    uc = cc.decompress(cc.compress(u, 10), 10); vc = cc.decompress(cc.compress(v, 4), 4)
    yv = ((vc - sum(cc.polymul(s[i][None], uc[:, i]) for i in range(k))) % Q).astype(np.int16)
    for bi in range(B):
        for j in range(8):
            got = lib.leech_decode_block(np.ascontiguousarray(yv[bi, 24 * j:24 * j + 24]))
            fails += int(got != msgs[bi, j]); blocks += 1
print(f"[4] real ML-KEM noise (k=2, eta={eta}): Leech block failure {fails/blocks:.3e} "
      f"({fails}/{blocks}), union-bound prediction {pred:.3e} [{time.time()-t0:.0f}s]")
