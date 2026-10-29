"""Correctness of the masked E8 decoder vs the verified unmasked decoder."""
import ctypes, numpy as np
Q = 3329
M = ctypes.CDLL("../v10/libe8masked.so")
R = ctypes.CDLL("../v4/libe8fast_chk.so")
U16 = np.ctypeslib.ndpointer(np.uint16, flags="C_CONTIGUOUS"); U8 = np.ctypeslib.ndpointer(np.uint8, flags="C_CONTIGUOUS")
I16 = np.ctypeslib.ndpointer(np.int16, flags="C_CONTIGUOUS")
M.e8_decode_masked.argtypes = [U16, U16, U8, U8, ctypes.c_int]; R.e8_decode_fast.argtypes = [I16, U8]
M.masked_seed.argtypes = [ctypes.c_uint64]
rng = np.random.default_rng(11); ok = tot = 0; shares_uniform = []
for it in range(3000):
    w = rng.integers(0, Q, 256).astype(np.int64)
    if it % 3 == 0:                                   # also exact codewords and boundary values
        w[:16] = rng.choice([0, 1, 832, 833, 1664, 1665, 1666, 2496, 2497, 3328], 16)
    w0 = rng.integers(0, Q, 256); w1 = (w - w0) % Q
    m0, m1, ref = np.zeros(32, np.uint8), np.zeros(32, np.uint8), np.zeros(32, np.uint8)
    M.masked_seed(int(rng.integers(1, 2**63)))
    M.e8_decode_masked(w0.astype(np.uint16), w1.astype(np.uint16), m0, m1, it % 2)
    R.e8_decode_fast(w.astype(np.int16), ref)
    ok += int(((m0 ^ m1) == ref).sum()); tot += 32
    shares_uniform.append(m0.copy())
s = np.concatenate(shares_uniform)
print(f"masked decoder (unshared output) vs verified decoder: {ok}/{tot} blocks identical (with and without shuffling)")
print(f"output share m0 alone: mean {s.mean():.1f} (uniform 127.5), distinct values {len(np.unique(s))}/256")
