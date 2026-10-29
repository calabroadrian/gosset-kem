"""Regenerate kat/GOSSET_KEM_KAT.json (deterministic). Run from the repository root after `make`."""
import ctypes, hashlib, json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = ctypes.CDLL(os.path.join(ROOT, "build", "libgosset.so"))
for f in ("hyb_pk_bytes", "hyb_sk_bytes", "hyb_ct_bytes"): getattr(L, f).restype = ctypes.c_size_t
K = [2, 3, 4] * 3; DUV = [(10, 4), (10, 4), (11, 5), (9, 4), (9, 5), (10, 4), (10, 4), (10, 4), (11, 5)]
NAMES = ["ML-KEM-512", "ML-KEM-768", "ML-KEM-1024", "Gosset-KEM-512-Compact", "Gosset-KEM-768-Compact", "Gosset-KEM-1024-Compact",
         "Gosset-KEM-512-MaxSec", "Gosset-KEM-768-MaxSec", "Gosset-KEM-1024-MaxSec"]
def seed(label, i, n=32): return hashlib.shake_256(f"GOSSET-KEM-KAT-v1|{label}|{i}".encode()).digest(n)
out = {"description": "Gosset-KEM known-answer tests (specification 1.0). Seeds = SHAKE256('GOSSET-KEM-KAT-v1|<label>|<i>'). "
       "Standard ML-KEM entries included as cross-check (must equal FIPS 203).", "kem": [], "hybrid": []}
for vid in range(9):
    k = K[vid]; ekl, dkl, ctl = 384 * k + 32, 768 * k + 96, 32 * (DUV[vid][0] * k + DUV[vid][1])
    for i in range(5):
        d, z, m = seed(f"d{vid}", i), seed(f"z{vid}", i), seed(f"m{vid}", i)
        ek, dk, c, Kk, K2 = (ctypes.create_string_buffer(n) for n in (ekl, dkl, ctl, 32, 32))
        L.mlkem_keygen_internal(vid, d, z, ek, dk); L.mlkem_encaps_internal(vid, ek.raw, m, Kk, c); L.mlkem_decaps(vid, dk.raw, c.raw, K2)
        assert Kk.raw == K2.raw
        out["kem"].append({"variant": NAMES[vid], "tc": i, "d": d.hex(), "z": z.hex(), "m": m.hex(),
                           "ek": ek.raw.hex(), "dk": dk.raw.hex(), "c": c.raw.hex(), "K": Kk.raw.hex()})
for mode, nm in ((1, "Hybrid-C (ML-KEM-768+X25519)"), (7, "Hybrid-S (Gosset-KEM-768-MaxSec+X25519)")):
    for i in range(5):
        d, z, skx, m, ex = (seed(f"h{mode}{t}", i) for t in "dzxme")
        pk, sk, ct, Kk, K2 = (ctypes.create_string_buffer(n) for n in (L.hyb_pk_bytes(mode), L.hyb_sk_bytes(mode), L.hyb_ct_bytes(mode), 32, 32))
        L.hyb_keygen_internal(mode, d, z, skx, pk, sk); L.hyb_encaps_internal(mode, pk.raw, m, ex, ct, Kk); L.hyb_decaps(mode, sk.raw, ct.raw, K2)
        assert Kk.raw == K2.raw
        out["hybrid"].append({"mode": nm, "tc": i, "d": d.hex(), "z": z.hex(), "sk_x25519": skx.hex(), "m": m.hex(), "eph_x25519": ex.hex(),
                              "pk": pk.raw.hex(), "sk": sk.raw.hex(), "ct": ct.raw.hex(), "K": Kk.raw.hex()})
s = json.dumps(out, indent=1)
open(os.path.join(ROOT, "kat", "GOSSET_KEM_KAT.json"), "w").write(s)
print(f"{len(out['kem'])} KEM + {len(out['hybrid'])} hybrid vectors, SHA3-256 {hashlib.sha3_256(s.encode()).hexdigest()}")
