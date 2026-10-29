"""AFA-12 v9 — X25519 validation: RFC 7748 vectors + 2000 random cases vs OpenSSL (python cryptography)."""
import ctypes
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey, X25519PublicKey
from cryptography.hazmat.primitives import serialization as ser
L = ctypes.CDLL("../v9/libx25519.so")
def X(k, u): o = ctypes.create_string_buffer(32); L.x25519(o, k, u); return o.raw
# RFC 7748 5.2 inputs; expected outputs cross-checked with OpenSSL (vector 2 value confirmed by OpenSSL)
vec = [("a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4", "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c"),
       ("4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d", "e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493")]
ok = True
for k, u in vec:
    ref = X25519PrivateKey.from_private_bytes(bytes.fromhex(k)).exchange(X25519PublicKey.from_public_bytes(bytes.fromhex(u)))
    ok &= X(bytes.fromhex(k), bytes.fromhex(u)) == ref
ok &= X(bytes.fromhex(vec[0][0]), bytes.fromhex(vec[0][1])).hex() == "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552"
print(f"[1] RFC 7748 sec. 5.2 vectors (vs RFC value and OpenSSL): {'OK' if ok else 'FAIL'}")
k = u = (9).to_bytes(32, "little")
for i in range(1, 1001):
    k, u = X(k, u), k
    if i == 1: r1 = k.hex()
print(f"[2] RFC 7748 iterated 1 / 1000: {'OK' if r1 == '422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079' else 'FAIL'} / "
      f"{'OK' if k.hex() == '684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51' else 'FAIL'}")
good = 0
for _ in range(2000):
    a, b = X25519PrivateKey.generate(), X25519PrivateKey.generate()
    ar = a.private_bytes(ser.Encoding.Raw, ser.PrivateFormat.Raw, ser.NoEncryption())
    bpub = b.public_key().public_bytes(ser.Encoding.Raw, ser.PublicFormat.Raw)
    good += X(ar, (9).to_bytes(32, "little")) == a.public_key().public_bytes(ser.Encoding.Raw, ser.PublicFormat.Raw) and X(ar, bpub) == a.exchange(b.public_key())
print(f"[3] vs OpenSSL, 2000 random key pairs + shared secrets: {good}/2000")
