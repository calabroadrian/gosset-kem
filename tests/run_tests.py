"""Gosset-KEM test suite. Run from the repository root after `make` (or simply `make test`)."""
import ctypes, hashlib, json, os, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B = os.path.join(ROOT, "build")
L = ctypes.CDLL(os.path.join(B, "libgosset.so"))
R = ctypes.CDLL(os.path.join(B, "libgosset_ref.so"))
M = ctypes.CDLL(os.path.join(B, "libgosset_masked.so"))
for f in ("hyb_pk_bytes", "hyb_sk_bytes", "hyb_ct_bytes"): getattr(L, f).restype = ctypes.c_size_t
M.mk_seed.argtypes = [ctypes.c_uint64]
H = bytes.fromhex
NAMES = ["ML-KEM-512", "ML-KEM-768", "ML-KEM-1024", "Gosset-KEM-512-Compact", "Gosset-KEM-768-Compact", "Gosset-KEM-1024-Compact",
         "Gosset-KEM-512-MaxSec", "Gosset-KEM-768-MaxSec", "Gosset-KEM-1024-MaxSec"]
K = [2, 3, 4] * 3
DUV = [(10, 4), (10, 4), (11, 5), (9, 4), (9, 5), (10, 4), (10, 4), (10, 4), (11, 5)]
XPLAT_DIGEST = "80254431656df5bfd7fa9f150bc04577"
sizes = lambda v: (384 * K[v] + 32, 768 * K[v] + 96, 32 * (DUV[v][0] * K[v] + DUV[v][1]))
failures = 0

def report(name, ok, n):
    global failures
    failures += ok != n
    print(f"  {'PASS' if ok == n else 'FAIL'}  {name}: {ok}/{n}")

print("[1] NIST ACVP FIPS 203 vectors (ML-KEM-512/768/1024)")
IDS = {"ML-KEM-512": 0, "ML-KEM-768": 1, "ML-KEM-1024": 2}
kg = json.load(open(os.path.join(ROOT, "kat", "keygen.json")))
ok = n = 0
for g in kg["testGroups"]:
    i = IDS[g["parameterSet"]]; ekl, dkl, _ = sizes(i)
    for t in g["tests"]:
        ek, dk = ctypes.create_string_buffer(ekl), ctypes.create_string_buffer(dkl)
        L.mlkem_keygen_internal(i, H(t["d"]), H(t["z"]), ek, dk); n += 1; ok += ek.raw == H(t["ek"]) and dk.raw == H(t["dk"])
report("keyGen", ok, n)
ed = json.load(open(os.path.join(ROOT, "kat", "encdec.json")))
res = {}
for g in ed["testGroups"]:
    fn = g["function"]; i = IDS[g["parameterSet"]]; _, _, ctl = sizes(i)
    for t in g["tests"]:
        if fn == "encapsulation":
            k, c = ctypes.create_string_buffer(32), ctypes.create_string_buffer(ctl)
            L.mlkem_encaps_internal(i, H(t["ek"]), H(t["m"]), k, c); good = k.raw == H(t["k"]) and c.raw == H(t["c"])
        elif fn == "decapsulation":
            k = ctypes.create_string_buffer(32); L.mlkem_decaps(i, H(t["dk"]), H(t["c"]), k); good = k.raw == H(t["k"])
        elif fn == "encapsulationKeyCheck":
            good = bool(L.mlkem_check_ek(i, H(t["ek"]))) == t["testPassed"]
        elif fn == "decapsulationKeyCheck":
            good = bool(L.mlkem_check_dk(i, H(t["dk"]))) == t["testPassed"]
        else:
            continue
        a, b = res.get(fn, (0, 0)); res[fn] = (a + good, b + 1)
for fn, (a, b) in res.items(): report(fn, a, b)

print("[2] Gosset-KEM known-answer tests (kat/GOSSET_KEM_KAT.json)")
kat = json.load(open(os.path.join(ROOT, "kat", "GOSSET_KEM_KAT.json")))
ok = 0
for e in kat["kem"]:
    v = NAMES.index(e["variant"]); ekl, dkl, ctl = sizes(v)
    ek, dk, c, k = (ctypes.create_string_buffer(x) for x in (ekl, dkl, ctl, 32))
    L.mlkem_keygen_internal(v, H(e["d"]), H(e["z"]), ek, dk); L.mlkem_encaps_internal(v, ek.raw, H(e["m"]), k, c)
    k2 = ctypes.create_string_buffer(32); L.mlkem_decaps(v, dk.raw, c.raw, k2)
    ok += ek.raw.hex() == e["ek"] and dk.raw.hex() == e["dk"] and c.raw.hex() == e["c"] and k.raw.hex() == e["K"] == k2.raw.hex()
report("KEM vectors (9 parameter sets)", ok, len(kat["kem"]))
ok = 0
for e in kat["hybrid"]:
    mode = 1 if e["mode"].startswith("Hybrid-C") else 7
    pk, sk, ct, k, k2 = (ctypes.create_string_buffer(x) for x in (L.hyb_pk_bytes(mode), L.hyb_sk_bytes(mode), L.hyb_ct_bytes(mode), 32, 32))
    L.hyb_keygen_internal(mode, H(e["d"]), H(e["z"]), H(e["sk_x25519"]), pk, sk)
    L.hyb_encaps_internal(mode, pk.raw, H(e["m"]), H(e["eph_x25519"]), ct, k); L.hyb_decaps(mode, sk.raw, ct.raw, k2)
    ok += pk.raw.hex() == e["pk"] and ct.raw.hex() == e["ct"] and k.raw.hex() == e["K"] == k2.raw.hex()
report("hybrid vectors (modes C and S)", ok, len(kat["hybrid"]))

print("[3] reference build == optimized build, implicit rejection")
same = rej = 0; n = 0
for v in range(9):
    ekl, dkl, ctl = sizes(v)
    for _ in range(40):
        d, z, m = os.urandom(32), os.urandom(32), os.urandom(32); outs = []
        for lib in (R, L):
            ek, dk, c, k, k2 = (ctypes.create_string_buffer(x) for x in (ekl, dkl, ctl, 32, 32))
            lib.mlkem_keygen_internal(v, d, z, ek, dk); lib.mlkem_encaps_internal(v, ek.raw, m, k, c); lib.mlkem_decaps(v, dk.raw, c.raw, k2)
            outs.append((ek.raw, dk.raw, c.raw, k.raw, k2.raw))
        n += 1; same += outs[0] == outs[1] and outs[1][3] == outs[1][4]
        cb = bytearray(outs[1][2]); cb[n % ctl] ^= 1 << (n % 8); k3 = ctypes.create_string_buffer(32)
        L.mlkem_decaps(v, outs[1][1], bytes(cb), k3); rej += k3.raw != outs[1][3]
report("identical outputs, 9 sets x 40", same, n)
report("tampered ciphertext -> unrelated key", rej, n)

print("[4] first-order masked decapsulation == unmasked (E8 sets, valid and tampered)")
ok = n = 0
for v in range(3, 9):
    ekl, dkl, ctl = sizes(v)
    for it in range(20):
        ek, dk, c, ke, ku, km = (ctypes.create_string_buffer(x) for x in (ekl, dkl, ctl, 32, 32, 32))
        L.mlkem_keygen_internal(v, os.urandom(32), os.urandom(32), ek, dk); L.mlkem_encaps_internal(v, ek.raw, os.urandom(32), ke, c)
        M.mk_seed(int.from_bytes(os.urandom(8), "little")); L.mlkem_decaps(v, dk.raw, c.raw, ku)
        f = M.mlkem_decaps_masked(v, dk.raw, c.raw, km); n += 1; ok += km.raw == ku.raw == ke.raw and f == 1
        cb = bytearray(c.raw); cb[it * 13 % ctl] ^= 1 << (it % 8)
        L.mlkem_decaps(v, dk.raw, bytes(cb), ku); f = M.mlkem_decaps_masked(v, dk.raw, bytes(cb), km); n += 1; ok += km.raw == ku.raw and f == 0
report("masked == unmasked", ok, n)

print("[5] X25519 (RFC 7748)")
def X(k, u): o = ctypes.create_string_buffer(32); L.x25519(o, k, u); return o.raw
r = X(H("a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4"), H("e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c")).hex()
k = u = (9).to_bytes(32, "little")
for i in range(1, 1001):
    k, u = X(k, u), k
    if i == 1: r1 = k.hex()
report("RFC 7748 5.2 vector, iterated 1 and 1000",
       (r == "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552")
       + (r1 == "422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079")
       + (k.hex() == "684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51"), 3)

print("[6] cross-platform digest (all schemes, hybrid, DEM, E8 decoders)")
out = subprocess.run([os.path.join(B, "xplat_kat")], capture_output=True, text=True).stdout
report(f"GLOBAL digest == {XPLAT_DIGEST}", int(XPLAT_DIGEST in out and "internal failures: 0" in out), 1)

print("\nRESULT:", "ALL TESTS PASSED" if failures == 0 else f"{failures} TEST GROUP(S) FAILED")
sys.exit(1 if failures else 0)
