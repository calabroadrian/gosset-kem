"""AFA-12 v9 — deterministic known-answer tests for the non-FIPS variants and hybrid modes (audit package)."""
import ctypes, hashlib, json
L = ctypes.CDLL("../v9/libhybrid.so")
for f in ("hyb_pk_bytes", "hyb_sk_bytes", "hyb_ct_bytes"): getattr(L, f).restype = ctypes.c_size_t
K = [2, 3, 4] * 3; DUV = [(10, 4), (10, 4), (11, 5), (9, 4), (9, 5), (10, 4), (10, 4), (10, 4), (11, 5)]
NAMES = ["ML-KEM-512", "ML-KEM-768", "ML-KEM-1024", "AFA12-512-Compact", "AFA12-768-Compact", "AFA12-1024-Compact",
         "AFA12-512-MaxSec", "AFA12-768-MaxSec", "AFA12-1024-MaxSec"]
def seed(label, i, n=32): return hashlib.shake_256(f"AFA12-KAT-v1|{label}|{i}".encode()).digest(n)
out = {"description": "AFA-12 known-answer tests. Seeds = SHAKE256('AFA12-KAT-v1|<label>|<i>'). "
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
for mode, nm in ((1, "Hybrid-C (ML-KEM-768+X25519)"), (7, "Hybrid-S (AFA12-768-MaxSec+X25519)")):
    for i in range(5):
        d, z, skx, m, ex = (seed(f"h{mode}{t}", i) for t in "dzxme")
        pk, sk, ct, Kk, K2 = (ctypes.create_string_buffer(n) for n in (L.hyb_pk_bytes(mode), L.hyb_sk_bytes(mode), L.hyb_ct_bytes(mode), 32, 32))
        L.hyb_keygen_internal(mode, d, z, skx, pk, sk); L.hyb_encaps_internal(mode, pk.raw, m, ex, ct, Kk); L.hyb_decaps(mode, sk.raw, ct.raw, K2)
        assert Kk.raw == K2.raw
        out["hybrid"].append({"mode": nm, "tc": i, "d": d.hex(), "z": z.hex(), "sk_x25519": skx.hex(), "m": m.hex(), "eph_x25519": ex.hex(),
                              "pk": pk.raw.hex(), "sk": sk.raw.hex(), "ct": ct.raw.hex(), "K": Kk.raw.hex()})
s = json.dumps(out, indent=1)
open("../../kat/AFA12_KAT.json", "w").write(s)
print(f"KAT file: {len(out['kem'])} KEM vectors + {len(out['hybrid'])} hybrid vectors, {len(s)//1024} KiB, SHA3-256 {hashlib.sha3_256(s.encode()).hexdigest()}")
# cross-check: standard entries must equal FIPS-validated reference build (libmlkem.so, NIST-validated)
R = ctypes.CDLL("../v7/libmlkem.so"); ok = 0; n = 0
for e in out["kem"][:15]:
    vid = NAMES.index(e["variant"]); k = K[vid]
    ek, dk = ctypes.create_string_buffer(384*k+32), ctypes.create_string_buffer(768*k+96)
    R.mlkem_keygen_internal(vid, bytes.fromhex(e["d"]), bytes.fromhex(e["z"]), ek, dk); n += 1; ok += ek.raw.hex() == e["ek"]
print(f"standard-variant entries match the NIST-validated reference build: {ok}/{n}")
