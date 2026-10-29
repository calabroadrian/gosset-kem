"""AFA-12 v9 — hybrid KEM tests: independent recomputation of the combiner, round trips,
tamper resistance of each component, FIPS 203 input checks through the public API."""
import ctypes, hashlib, os
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey, X25519PublicKey
from cryptography.hazmat.primitives import serialization as ser
L = ctypes.CDLL("../v9/libhybrid.so")
for f in ("hyb_pk_bytes", "hyb_sk_bytes", "hyb_ct_bytes"): getattr(L, f).restype = ctypes.c_size_t
LABEL = b"AFA12-HYBRID-v1"
EK = {1: 1184, 7: 1184}; CT = {1: 1088, 7: 1088}
for mode, name in ((1, "Mode C: ML-KEM-768 + X25519"), (7, "Mode S: AFA12-768-MaxSec + X25519")):
    pkl, skl, ctl = L.hyb_pk_bytes(mode), L.hyb_sk_bytes(mode), L.hyb_ct_bytes(mode)
    comb = rt = 0; N = 300
    for _ in range(N):
        d, z, skx, m, ex = (os.urandom(32) for _ in range(5))
        pk, sk = ctypes.create_string_buffer(pkl), ctypes.create_string_buffer(skl)
        L.hyb_keygen_internal(mode, d, z, skx, pk, sk)
        ct, K, K2 = ctypes.create_string_buffer(ctl), ctypes.create_string_buffer(32), ctypes.create_string_buffer(32)
        L.hyb_encaps_internal(mode, pk.raw, m, ex, ct, K)
        L.hyb_decaps(mode, sk.raw, ct.raw, K2); rt += K.raw == K2.raw
        # independent recomputation
        Kpq, cpq = ctypes.create_string_buffer(32), ctypes.create_string_buffer(CT[mode])
        L.mlkem_encaps_internal(mode, pk.raw, m, Kpq, cpq)
        pkx = X25519PrivateKey.from_private_bytes(skx).public_key().public_bytes(ser.Encoding.Raw, ser.PublicFormat.Raw)
        ctx = X25519PrivateKey.from_private_bytes(ex).public_key().public_bytes(ser.Encoding.Raw, ser.PublicFormat.Raw)
        ssx = X25519PrivateKey.from_private_bytes(ex).exchange(X25519PublicKey.from_public_bytes(pkx))
        Kref = hashlib.sha3_256(LABEL + Kpq.raw + ssx + ctx + pkx).digest()
        comb += (Kref == K.raw and ct.raw == cpq.raw + ctx and pk.raw[EK[mode]:] == pkx)
    # tampering of each component, through the public API (OS randomness)
    pk, sk, ct, K, K2 = (ctypes.create_string_buffer(n) for n in (pkl, skl, ctl, 32, 32))
    assert L.hyb_keygen(mode, pk, sk) == 0 and L.hyb_encaps(mode, pk.raw, ct, K) == 0
    t_pq = t_x = 0
    for i in range(100):
        c1 = bytearray(ct.raw); c1[i * 7 % CT[mode]] ^= 1
        L.hyb_decaps_checked(mode, sk.raw, bytes(c1), K2); t_pq += K2.raw != K.raw
        c2 = bytearray(ct.raw); c2[CT[mode] + i % 32] ^= 1 << (i % 8)
        L.hyb_decaps_checked(mode, sk.raw, bytes(c2), K2); t_x += K2.raw != K.raw
    L.hyb_decaps_checked(mode, sk.raw, ct.raw, K2); ok_api = K2.raw == K.raw
    badpk = bytearray(pk.raw); badpk[0] = 0xFF; badpk[1] |= 0x0F            # coefficient 4095 >= q
    rej_pk = L.hyb_encaps(mode, bytes(badpk), ct, K) == -1
    badsk = bytearray(sk.raw); badsk[1200] ^= 1                             # breaks H(ek) inside dk
    rej_sk = L.hyb_decaps_checked(mode, bytes(badsk), ct.raw, K2) == -1
    print(f"{name}  (pk {pkl} B, sk {skl} B, ct {ctl} B)")
    print(f"   combiner recomputed independently (OpenSSL X25519 + hashlib SHA3): {comb}/{N}")
    print(f"   round trips (internal): {rt}/{N};  public API with OS randomness: {'OK' if ok_api else 'FAIL'}")
    print(f"   tampered lattice part -> different key: {t_pq}/100;  tampered X25519 part -> different key: {t_x}/100")
    print(f"   invalid public key rejected: {rej_pk};  corrupted secret key rejected: {rej_sk}")
