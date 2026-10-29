"""AFA-12 v7 — (1) FIPS 203 input checks vs NIST ACVP keyCheck vectors, both builds;
(2) reference vs optimized build produce identical outputs for all 6 variants (incl. E8)."""
import ctypes, json, os
H = bytes.fromhex
IDS = {"ML-KEM-512": 0, "ML-KEM-768": 1, "ML-KEM-1024": 2}
libs = {n: ctypes.CDLL(f"../v7/{n}") for n in ("libmlkem.so", "libmlkem_opt.so")}
ed = json.load(open("../../kat/encdec.json"))
for name, L in libs.items():
    tot = ok = 0
    for g in ed["testGroups"]:
        if g["function"] not in ("encapsulationKeyCheck", "decapsulationKeyCheck"):
            continue
        i = IDS[g["parameterSet"]]
        for t in g["tests"]:
            r = L.mlkem_check_ek(i, H(t["ek"])) if g["function"] == "encapsulationKeyCheck" else L.mlkem_check_dk(i, H(t["dk"]))
            tot += 1; ok += (bool(r) == t["testPassed"])
    print(f"[1] {name:<17} ACVP key checks: {ok}/{tot} {'OK' if ok == tot else 'FAIL'}")
A, B = libs["libmlkem.so"], libs["libmlkem_opt.so"]
K = [2, 3, 4, 2, 3, 4]; DUV = [(10, 4), (10, 4), (11, 5), (9, 4), (9, 5), (10, 4)]
same = 0; n = 0
for vid in range(6):
    k = K[vid]; ekl, dkl, ctl = 384 * k + 32, 768 * k + 96, 32 * (DUV[vid][0] * k + DUV[vid][1])
    for _ in range(200):
        d, z, m = os.urandom(32), os.urandom(32), os.urandom(32)
        outs = []
        for L in (A, B):
            ek, dk = ctypes.create_string_buffer(ekl), ctypes.create_string_buffer(dkl)
            L.mlkem_keygen_internal(vid, d, z, ek, dk)
            kk, c, k2 = ctypes.create_string_buffer(32), ctypes.create_string_buffer(ctl), ctypes.create_string_buffer(32)
            L.mlkem_encaps_internal(vid, ek, m, kk, c); L.mlkem_decaps(vid, dk, c, k2)
            outs.append((ek.raw, dk.raw, c.raw, kk.raw, k2.raw))
        n += 1; same += outs[0] == outs[1] and outs[0][3] == outs[0][4]
print(f"[2] reference vs optimized build, all 6 variants: {same}/{n} identical (ek, dk, c, K) and K==K'")
