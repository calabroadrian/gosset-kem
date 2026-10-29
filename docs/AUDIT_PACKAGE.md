# Gosset-KEM — External Audit Package

*Formerly developed under the working name AFA-12. Paths of the form `vN/…` refer to `research/vN/`; the production code is in `src/`.*

**Version:** 1.0 (2026-10-24) · **Author:** Adrian Calabro · **Status:** research prototype, not yet audited

This package gives an independent reviewer everything needed to evaluate Gosset-KEM: the exact specification, the claims and their evidence, reproducible tests, known limitations, and the questions we want answered.

---

## 1. What Gosset-KEM is

Gosset-KEM is a **variant of ML-KEM (FIPS 203)**. Every component is identical to FIPS 203 except:

| Element | FIPS 203 | Gosset-KEM |
|---|---|---|
| Message encoding (K-PKE.Encrypt, line `μ ← Decompress₁(ByteDecode₁(m))`) | 1 bit per coefficient, μᵢ = mᵢ·1665 | **E8 lattice code**: byte b (0..255) of m → 8 coefficients (see §2.1) |
| Message decoding (K-PKE.Decrypt, line `m ← ByteEncode₁(Compress₁(w))`) | per-coefficient threshold | **nearest point of (q/2)·E8** per 8 coefficients (§2.2) |
| Parameters | Table 2 of FIPS 203 | §2.3 (Compact or MaxSec) |
| CBD | η ∈ {2,3} | adds **η = 4** (same definition as Alg. 8 with η=4) |

Everything else (G, H, J, PRF, XOF, SampleNTT, NTT, Fujisaki-Okamoto with implicit rejection, key formats, input checks of §7.2/§7.3) is unchanged.

### 2.1 Encoding (per byte j of m, coefficients E8_LAYOUT[8j..8j+7])
- coset c = bit 7 of b; x₀..x₆ = bits 0..6 of b; x₇ = x₀ ⊕ … ⊕ x₆ (even weight).
- If c = 0: μᵢ = 1665·xᵢ. If c = 1: μᵢ = 832 + 1665·xᵢ.
- (Equivalently: point (q/2)·(x + c/2) of (q/2)·E8, rounded; E8 = D8 ∪ (D8 + ½).)
- **Block layout (v13):** μ value i of block j is placed at coefficient E8_LAYOUT[8j+i] (`v7/e8_layout.h`, a fixed public permutation of 0..255). In every block the 28 pairwise index distances (classes min(d, 256−d)) are all distinct. This defeats correlation-based failure boosting; see C9. The decoder gathers the same 8 coefficients.

### 2.2 Decoding (per block: the 8 coefficients E8_LAYOUT[8j..8j+7], w in [0, q))
- yᵢ = w centered in (−q/2, q/2]; Yᵢ = 4yᵢ; U = 2q.
- For off ∈ {0, q}: fᵢ = round_half_up((Yᵢ − off)/U); if Σfᵢ odd, move the coordinate with largest |Yᵢ − off − U·fᵢ| one step toward Yᵢ (first index on ties); dist_off = Σ(Yᵢ − off − U·fᵢ)².
- Choose off = q iff dist_q < dist_0. Output byte = (bit 7 = [off = q]) | Σ_{i<7} (fᵢ mod 2)·2ⁱ.
- Implementation: constant time, value barriers against compiler branch insertion (see §5).

### 2.3 Parameter sets

| Name | k | η₁ | η₂ | d_u | d_v | ek | dk | ct | DFR (exact, union bound) | Core-SVP classical / quantum |
|---|---|---|---|---|---|---|---|---|---|---|
| Gosset-KEM-512-Compact | 2 | 3 | 2 | 9 | 4 | 800 | 1632 | **704** | 2⁻¹⁴²·⁰ | 118.6 / 107.6 (= ML-KEM-512) |
| Gosset-KEM-768-Compact | 3 | 2 | 2 | 9 | 5 | 1184 | 2400 | **1024** | 2⁻¹⁹⁴·¹ | 182.2 / 165.4 (= ML-KEM-768) |
| Gosset-KEM-1024-Compact | 4 | 2 | 2 | 10 | 4 | 1568 | 3168 | **1408** | 2⁻²¹⁵·⁰ | 255.2 / 231.6 (= ML-KEM-1024) |
| Gosset-KEM-512-MaxSec | 2 | **4** | **3** | 10 | 4 | 800 | 1632 | 768 | 2⁻¹⁴⁶·⁶ | key primal 123.8 / 112.4; weakest link 122.3 / 111.0 |
| Gosset-KEM-768-MaxSec | 3 | **3** | 2 | 10 | 4 | 1184 | 2400 | 1088 | 2⁻¹⁶⁶·⁸ | key primal 192.7 / 174.9; weakest link 190.7 / 173.0 |
| Gosset-KEM-1024-MaxSec | 4 | **3** | 2 | 11 | 5 | 1568 | 3168 | 1568 | 2⁻¹⁸⁹·⁴ | key primal 269.5 / 244.6; weakest link 263.4 / 239.0 |

The weakest-link column is the minimum over {primal, dual} × {key, ciphertext} with d_u rounding credited (v12/). Without the rounding credit the classical figures are 119.5 / 185.5 / 259.9.

### 2.4 Hybrid KEM (X-Wing-style)
pk = ek ‖ pkₓ, sk = dk ‖ skₓ ‖ pkₓ, ct = c ‖ ctₓ; **K = SHA3-256("GOSSET-KEM-HYBRID-v1" ‖ K_pq ‖ ssₓ ‖ ctₓ ‖ pkₓ)**.
Mode C = ML-KEM-768 + X25519 (FIPS 203 component). Mode S = Gosset-KEM-768-MaxSec + X25519.

### 2.5 DEM
(K_enc ‖ K_mac) = SHAKE256("GOSSET-KEM-DEM-KDF-v1" ‖ K, 64); M_pad = ISO/IEC 7816-4 to a multiple of 64 bytes; C = M_pad ⊕ SHAKE256("GOSSET-KEM-DEM-ENC-v1" ‖ K_enc ‖ nonce₁₆); tag = HMAC-SHA3-256(K_mac, "GOSSET-KEM-DEM-MAC-v1" ‖ nonce ‖ C); verify before decrypt.

---

## 3. Security claims and evidence

| # | Claim | Evidence (reproducible) | Strength |
|---|---|---|---|
| C1 | Implementation of standard ML-KEM is correct | NIST ACVP: 180/180 keyGen/encap/decap + 60/60 key checks; ref vs optimized builds identical (1,200/1,200) | Strong |
| C2 | Gosset-KEM hardness = Module-LWE with the stated noise; MaxSec gain over ML-KEM | Core-SVP primal + dual on both key and ciphertext LWE instances (v12/coresvp.py), reproducing Kyber R3 primal 118.6/181.9/254.9 and dual 117.0/180.2/252.6. Weakest-link gain: **+5.3/+10.5/+10.8 bits** classical (ciphertext rounding credited) or **+6.1/+5.3/+7.3** (not credited); quantum +4.8/+9.5/+9.8. The ciphertext instance is the bottleneck at 768/1024. A search over (η₁, η₂, d_u, d_v) under FIPS ct size and DFR constraints confirms the MaxSec sets are optimal | **Confirmed with the official lattice-estimator** (commit f31b86f, all attacks except arora-gb/bkw, MATZOV cost model; v12/official_estimate.py). The strongest attack is dual-hybrid throughout. Weakest-link log₂ rop for ML-KEM → MaxSec: 512: 136.0→141.6 / 139.7→144.9; 768: 196.4→201.2 / 196.4→207.2; 1024: 262.3→267.6 / 262.3→271.7 (each pair: rounding not credited / credited). **Gain +4.8 to +10.8 bits**, matching core-SVP within 1 bit |
| C3 | DFR as stated | Exact convolution (reproduces FIPS 203 DFRs), union bound; real-encryption validation at inflated noise matches the bound (E8: 1.711e-5 vs 1.705e-5) | Strong in measurable range; extrapolation is the same assumption FIPS 203 uses |
| C4 | E8/Leech noise *in the error* would weaken LWE | Entropy theorem (−log₂ det L), validated estimator, real LLL attacks | Strong (negative result; not used in Gosset-KEM) |
| C5 | Decoder and decapsulation are constant time | Assembly inspection; dudect Welch-t ≤ 3.1 at 10⁵–10⁷ measurements; positive control detected (|t| ≈ 26,000) | One CPU, one compiler (gcc 13, x86-64) |
| C6 | Hybrid key is secure if either component is | X-Wing-style combiner (binds ctₓ and pkₓ); recomputed independently with OpenSSL + hashlib (600/600) | Combiner argument by analogy with X-Wing; **no formal proof for this exact construction** |
| C7 | X25519 correct | RFC 7748 vectors + 2,000 random cases vs OpenSSL | Strong |
| C9 | Failure-boosting attacks are infeasible (≤ 2⁶⁴ queries) | v13/: norm criterion (D'Anvers 2019) and an **E8-specific in-block correlation criterion**. The per-ciphertext model was validated on 80k real encryptions (299 failures observed vs 303 predicted). The correlation tail uses a Chernoff bound on the negacyclic quadratic form, about 10× above the Monte Carlo, so it favours the attacker. **First-failure cost, Q ≤ 2⁶⁴: 2¹⁸⁵ / 2²²⁵ (512 Compact/MaxSec), 2³⁵⁸ / 2³⁷⁹ (768), 2⁷⁸³ / 2⁵⁴⁵ (1024)** | **v14:** a tight Lugannani-Rice tail (matches Monte Carlo within 15% down to 10⁻⁴) gives first-failure costs of 2¹⁹² / 2²³² (512), 2³⁶⁵ / 2³⁸⁷ (768) and 2⁷⁹¹ / 2⁵⁵² (1024). Directional boosting (D'Anvers et al. 2020 geometric model, attacker-favourable) makes later failures cheaper in every scheme, but they stay costlier for Gosset-KEM than for ML-KEM (e.g. 2¹⁴⁵ vs 2⁸⁵ for the second failure at level 512), so the first failure dominates. Multi-target attacks are not counted, as for ML-KEM. With the earlier consecutive layout the bound gave **2¹²² / 2¹⁴³ at level 512, at or below the claimed level**, so the layout is mandatory |
| C8 | DEM: IND-CCA style authenticated encryption | Encrypt-then-MAC, key separation; 528/528 tamperings rejected; HMAC-SHA3 = Python hmac | Standard construction; SHAKE256 used as stream cipher (non-standard choice, see Q5) |

---

## 4. Known limitations and risks (please verify)

| Risk | Severity | Status |
|---|---|---|
| **Power/EM side channels of the E8 decoder.** Simulation (Hamming-weight model + profiled LDA attack, single trace): at measurement noise σ=8 the standard decoder leaks ~1.1% of message bytes, **the E8 decoder ~15.2%**. More secret-dependent operations → more leakage. | **High for devices exposed to physical access** (smart cards, IoT, HSM without shielding) | **First-order masked full decapsulation prototype** (v11/): masked u·s/NTT, bitsliced masked E8 decoder (34k cycles, 12.3× faster than v10), masked Keccak (SHA3-512/SHAKE256), masked CBD+B2A re-encryption, masked ciphertext comparison, key selection on shares. Bit-exact 720/720 (valid + tampered, 6 parameter sets). Cost ×4.7–6.4 vs unmasked. Simulated TVLA (HW model, σ=1, 43,956 points, up to 30k+30k traces, 4 runs): first-order max|t| 4.5–5.04, Bonferroni threshold 4.87, no point repeats across independent runs (ISO 17825 dual-run) → no first-order leak detected; second order |t|≈170 (expected). Pending: production CSPRNG (prototype: xorshift128+), key stored masked with refresh, hardware validation, higher order. Until hardware-validated, deploy Mode S only where the attacker has no physical access. |
| Compiler-inserted branches | Medium | Found twice: Leech decoder (gcc -O3) and the reference E8 decoder `v4/e8codec.c` under **clang 18 -O3** (dudect |t| ≈ 1,700; gcc was clean). Both fixed with value barriers. After the fix, every decoder and the full decapsulation pass under gcc 13 and clang 18 (max |t| ≤ 2.5). The KEM always used `e8_decode_avx2`, which was never affected. Windows builds (MinGW gcc and clang) give bit-identical outputs to Linux (local/xplat_kat.c). **AMD Zen 3 (Ryzen 5 5600X):** all decoders are clean (≤ 1.8 at 10M), and clang decaps is clean for all 12 sets. The gcc binary shows a reproducible **fixed-input** effect: |t| 10–19 at 1M, but only 5–10 cycles of 85k (< 0.012%). The repeated identical input runs faster, and ML-KEM-512 shows the same effect. Valid vs invalid ciphertext stays indistinguishable (|t| ≤ 1.9), so this is classified as a low-impact microarchitectural artifact on public input. Recommendation: clang builds, plus dudect per target CPU |
| η = 4 is outside FIPS 203 | Low (cryptographic), High (compliance) | CBD definition unchanged; validated vs the bit-level definition |
| Not FIPS 203 → not certifiable as is | High (compliance) | See §6 |
| **Correlation-based failure boosting (fixed in v13)** | Was high for level 512 | Consecutive E8 blocks let an attacker who grinds ciphertexts for one large autocorrelation lag raise many pair-variances at once. The adversarial DFR after 2⁶⁴ grinding was 2⁻⁶⁶ (512-Compact), vs 2⁻¹⁰⁴ for ML-KEM-512. With the distinct-lag layout it is 2⁻⁸¹, and the attack cost rises from 2¹²² to 2¹⁸⁵. The average DFR is unchanged (verified on real encryptions). The KATs were regenerated (SHA3-256 4b1c0e37…; the v12 KATs with consecutive blocks remain in the research history) |
| Prior art: E8 coding for Kyber (Liu & Sakzad, arXiv:2601.08452) | — | Gosset-KEM contributions: exact DFRs for all levels, real-encryption validation, MaxSec mode, constant-time C implementation, hybrid + DEM integration |

---

## 5. Reproduction

| Folder | Content | Command |
|---|---|---|
| v2/ | negative result (structured noise) | `python3 estimator.py`, `python3 real_attack.py 24 48 257 40 2 5,6,7,8,10,12,14` |
| v3/ | exact DFR, correlation validation | `python3 e8_message_coding.py`, `python3 correlation_check.py 400 100` |
| v4/ | E8 codec in C, constant time | `gcc -O3 -mavx2 bench_fast.c e8fast.c e8codec.c -lm && ./a.out 5000000` |
| v5/ | Leech evaluation and decoder | `python3 leech_eval.py`, `python3 test_leech.py` |
| v7/ | full ML-KEM + Gosset-KEM, DEM, demo | `python3 test_kat_opt.py`, `python3 test_checks_equiv.py`, `./bench_maxsec`, `./demo` |
| v8/ | MaxSec parameter search | `python3 max_security.py` |
| v9/ | X25519, hybrid, power simulation, KATs | `python3 test_x25519.py`, `python3 test_hybrid.py`, `./bench_hybrid`, `python3 power_sim.py`, `python3 gen_kat.py` |
| v14/ | saddlepoint correlation tail, directional boosting, Leech layout exposure | `python3 corr_saddle.py validate`, `python3 corr_saddle.py costs`, `python3 directional.py`, `python3 leech_layout.py` |
| v13/ | failure boosting (norm + E8 correlation), layout search, real-encryption validation | `python3 failure_boost.py`, `python3 corr_chernoff.py`, `python3 attack_cost.py`, `python3 validate_corr.py 14 40`, `python3 dfr_layout_check.py 60` |
| v12/ | core-SVP primal + dual, key + ciphertext instances, MaxSec rebalance search, official lattice-estimator | `python3 coresvp.py`, `python3 rebalance.py`, `git clone https://github.com/malb/lattice-estimator` then `python3 official_estimate.py rough` / `full` (Sage required, e.g. `pip install passagemath-standard`) |
| v10/ | masked E8 decoder (per block) | `python3 test_masked.py`, `python3 power_masked.py` |
| v11/ | bitsliced masked decoder, full masked decapsulation, TVLA | `python3 test_mdecaps.py`, `./bench_mdecaps <reps> <traces> <control 0/1> <seed>` (build line in the file header) |

Known-answer tests: `audit/GOSSET_KEM_KAT.json` (45 KEM + 10 hybrid vectors; seeds derived from SHAKE256 labels; the 15 standard ML-KEM entries match the NIST-validated build).

---

## 6. Certification path (FIPS 140-3 environments)

1. **Gosset-KEM itself is not an approved algorithm** and cannot be FIPS-validated as is.
2. **Compliant deployment today: Hybrid Mode C** (ML-KEM-768 FIPS 203 + X25519), combined with an approved hash-based derivation where the approved shared secret is included. This is the same pattern as X25519MLKEM768 in TLS. NIST guidance for combining a shared secret from an approved scheme with other secrets: SP 800-56C Rev. 2; KEM combiners are discussed in SP 800-227. *The exact combiner format must be checked by the certification lab.*
3. **Compliance + Gosset-KEM margin:** a three-component hybrid (ML-KEM-768 + Gosset-KEM-768-MaxSec + X25519). Approved component present; Gosset-KEM adds margin. Cost: about 2.1 KB more per exchange.
4. **Long term:** submit Gosset-KEM encodings/parameters as an IETF CFRG draft and as public comment for a future FIPS 203 revision; obtain independent cryptanalysis first.

---

## 7. Questions for the auditor

1. Is the union-bound + independence model for the E8 DFR acceptable at 2⁻¹⁴⁰…2⁻²¹⁵ (validated empirically to ~10⁻⁶)? The correlation effect on failure boosting is analysed in v13 (C9). Please check the Chernoff tail model and whether directional failure boosting changes the conclusion.
2. Is the MaxSec security estimate (primal, key noise η₁ only) adequate, or should dual/hybrid attacks and the ciphertext LWE instance be estimated separately with the lattice-estimator?
3. Does increasing η change any assumption in the FO/IND-CCA2 proof of ML-KEM beyond the DFR term?
4. Is the X-Wing-style combiner with Gosset-KEM (not ML-KEM) as the PQ component covered by the X-Wing security argument? What extra properties of Gosset-KEM are required (e.g., ciphertext collision resistance)?
5. DEM: is SHAKE256 as a keystream acceptable, or should AES-256-GCM / ChaCha20-Poly1305 be mandated for deployment?
6. Side channels: confirm or refute the simulated power leakage on real hardware (e.g., ChipWhisperer on Cortex-M4); recommend a masking strategy for the E8 decoder within a masked FO.
7. Review constant-time properties on other compilers/architectures (clang, MSVC, ARM).
