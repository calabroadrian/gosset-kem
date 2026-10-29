"""Check: multiplication-free reformulation of the E8 decoder == verified C decoder (e8_decode_fast)."""
import ctypes, numpy as np
Q = 3329
L = ctypes.CDLL("../v4/libe8fast_chk.so")
I16P = np.ctypeslib.ndpointer(np.int16, flags="C_CONTIGUOUS"); U8P = np.ctypeslib.ndpointer(np.uint8, flags="C_CONTIGUOUS")
L.e8_decode_fast.argtypes = [I16P, U8P]
def linear_decode(w):                       # w: (N,8) in [0,q)
    y = np.where(w > Q // 2, w - Q, w); Y = 4 * y
    s = (Y < 0).astype(np.int64); A = np.abs(Y)
    t = (A > Q).astype(np.int64)            # no ties: 4y != +-q (q odd)
    a = np.where(t == 1, 2 * Q - A, A)      # |error| of candidate D8
    parA = t.sum(1) & 1; parB = s.sum(1) & 1
    amax = a.max(1); amin = a.min(1)
    iA = np.argmax(a, 1); iB = np.argmin(a, 1)   # first index on ties
    D = a.sum(1) - 4 * Q + parA * (2 * Q - 2 * amax) - parB * 2 * amin
    useB = (D > 0).astype(np.int64)
    hA = np.zeros_like(a); hA[np.arange(len(a)), iA] = 1
    hB = np.zeros_like(a); hB[np.arange(len(a)), iB] = 1
    bA = t ^ (parA[:, None] & hA); bB = s ^ (parB[:, None] & hB)
    bits = np.where(useB[:, None] == 1, bB, bA)
    return (bits[:, :7] << np.arange(7)).sum(1) | (useB << 7)
rng = np.random.default_rng(3); N = 400000
w = rng.integers(0, Q, size=(N // 32, 256)).astype(np.int16)
ref = np.zeros((N // 32, 32), np.uint8)
for i in range(len(w)): L.e8_decode_fast(np.ascontiguousarray(w[i]), ref[i])
lin = linear_decode(w.reshape(-1, 8).astype(np.int64)).reshape(-1, 32)
print(f"linear reformulation vs verified decoder: {int((lin == ref).sum())}/{ref.size} blocks identical")
