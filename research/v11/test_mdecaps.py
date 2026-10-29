"""Masked decapsulation == validated unmasked decapsulation (valid and tampered ciphertexts, all E8 variants)."""
import ctypes, os
L = ctypes.CDLL("../v11/libmdecaps.so")
L.mk_seed.argtypes = [ctypes.c_uint64]
K = {3: 2, 4: 3, 5: 4, 6: 2, 7: 3, 8: 4}; DUV = {3: (9, 4), 4: (9, 5), 5: (10, 4), 6: (10, 4), 7: (10, 4), 8: (11, 5)}
NAME = {3: "AFA12-512-Compact", 4: "AFA12-768-Compact", 5: "AFA12-1024-Compact", 6: "AFA12-512-MaxSec", 7: "AFA12-768-MaxSec", 8: "AFA12-1024-MaxSec"}
for vid in range(3, 9):
    k = K[vid]; ekl, dkl, ctl = 384*k+32, 768*k+96, 32*(DUV[vid][0]*k + DUV[vid][1])
    ok_v = ok_t = flag_v = flag_t = 0; N = 60
    for it in range(N):
        ek, dk, c = ctypes.create_string_buffer(ekl), ctypes.create_string_buffer(dkl), ctypes.create_string_buffer(ctl)
        Ke, Ku, Km = (ctypes.create_string_buffer(32) for _ in range(3))
        L.mlkem_keygen_internal(vid, os.urandom(32), os.urandom(32), ek, dk)
        L.mlkem_encaps_internal(vid, ek.raw, os.urandom(32), Ke, c)
        L.mk_seed(int.from_bytes(os.urandom(8), "little"))
        L.mlkem_decaps(vid, dk.raw, c.raw, Ku)
        f = L.mlkem_decaps_masked(vid, dk.raw, c.raw, Km)
        ok_v += Km.raw == Ku.raw == Ke.raw; flag_v += f == 1
        cb = bytearray(c.raw); cb[it * 13 % ctl] ^= 1 << (it % 8)     # tampered ciphertext
        L.mlkem_decaps(vid, dk.raw, bytes(cb), Ku)
        f = L.mlkem_decaps_masked(vid, dk.raw, bytes(cb), Km)
        ok_t += Km.raw == Ku.raw; flag_t += f == 0
    print(f"{NAME[vid]:<20} valid ct: key {ok_v}/{N}, flag {flag_v}/{N} | tampered ct: key {ok_t}/{N} (implicit rejection), flag {flag_t}/{N}")
