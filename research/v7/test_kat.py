"""AFA-12 v7 — validate our ML-KEM against NIST ACVP FIPS 203 vectors (keyGen, encap, decap)."""
import ctypes, json, os, sys
L = ctypes.CDLL("../v7/libmlkem.so")
IDS = {"ML-KEM-512": 0, "ML-KEM-768": 1, "ML-KEM-1024": 2}
K = {0: 2, 1: 3, 2: 4}
DU = {0: (10, 4), 1: (10, 4), 2: (11, 5)}
ekb = lambda i: 384 * K[i] + 32
dkb = lambda i: 768 * K[i] + 96
ctb = lambda i: 32 * (DU[i][0] * K[i] + DU[i][1])
H = bytes.fromhex
res = {}

kg = json.load(open("../../kat/keygen.json"))
for g in kg["testGroups"]:
    i = IDS[g["parameterSet"]]; ok = 0
    for t in g["tests"]:
        ek = ctypes.create_string_buffer(ekb(i)); dk = ctypes.create_string_buffer(dkb(i))
        L.mlkem_keygen_internal(i, H(t["d"]), H(t["z"]), ek, dk)
        ok += (ek.raw == H(t["ek"]) and dk.raw == H(t["dk"]))
    res[(g["parameterSet"], "keyGen")] = (ok, len(g["tests"]))

ed = json.load(open("../../kat/encdec.json"))
for g in ed["testGroups"]:
    if g["function"] not in ("encapsulation", "decapsulation"):
        continue
    i = IDS[g["parameterSet"]]; ok = 0
    for t in g["tests"]:
        if g["function"] == "encapsulation":
            k = ctypes.create_string_buffer(32); c = ctypes.create_string_buffer(ctb(i))
            L.mlkem_encaps_internal(i, H(t["ek"]), H(t["m"]), k, c)
            ok += (k.raw == H(t["k"]) and c.raw == H(t["c"]))
        else:
            k = ctypes.create_string_buffer(32)
            L.mlkem_decaps(i, H(t["dk"]), H(t["c"]), k)
            ok += (k.raw == H(t["k"]))
    res[(g["parameterSet"], g["function"])] = (ok, len(g["tests"]))

allok = True
for (ps, fn), (ok, n) in sorted(res.items()):
    print(f"  {ps:<12} {fn:<14} {ok:>3}/{n}  {'OK' if ok == n else 'FAIL'}")
    allok &= ok == n
print("NIST ACVP FIPS 203 vectors:", "ALL PASS" if allok else "FAILURES")
