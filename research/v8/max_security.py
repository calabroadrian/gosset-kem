"""
AFA-12 v8 — spend the E8 gain on SECURITY instead of size.
Same ciphertext/key sizes as FIPS 203 (same k, du, dv). Raise the noise (eta1, eta2)
as far as the decryption-failure target allows, with standard vs E8 message coding.
Security: primal estimate validated vs Kyber (406/624/874); sigma = sqrt(eta1/2).
Failure: exact distribution (Kyber-spec method), union bound for E8 (validated in v3/v5).
"""
import math, sys
sys.path.insert(0, "../v2"); sys.path.insert(0, "../v3")
from estimator import primal_beta
from e8_message_coding import kyber_noise, failure_standard, failure_e8, log2

levels = [("ML-KEM-512", 2, 10, 4, -139, 3, 2), ("ML-KEM-768", 3, 10, 4, -164, 2, 2), ("ML-KEM-1024", 4, 11, 5, -174, 2, 2)]
print(f"{'level':<12}{'code':<6}{'eta1':>5}{'eta2':>5}{'log2 DFR':>10}{'beta':>6}{'classical':>10}{'quantum':>9}")
for name, k, du, dv, tgt, e1_0, e2_0 in levels:
    for code in ("std", "E8"):
        best = None
        for e1 in range(2, 9):
            for e2 in range(2, 9):
                w = kyber_noise(k, e1, e2, du, dv)
                d = log2(failure_standard(w, 3329) if code == "std" else failure_e8(w, 3329)[0])
                if d > tgt:
                    continue
                b, _ = primal_beta(256 * k, 3329, math.sqrt(e1 / 2))
                key = (b, e2)
                if best is None or key > best[0]:
                    best = (key, e1, e2, d, b)
        (_, e1, e2, d, b) = best
        print(f"{name:<12}{code:<6}{e1:>5}{e2:>5}{d:>10.1f}{b:>6}{0.292*b:>10.1f}{0.265*b:>9.1f}", flush=True)
    b0, _ = primal_beta(256 * k, 3329, math.sqrt(e1_0 / 2))
    print(f"{'  FIPS 203':<12}{'std':<6}{e1_0:>5}{e2_0:>5}{'':>10}{b0:>6}{0.292*b0:>10.1f}{0.265*b0:>9.1f}")
