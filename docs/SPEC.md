# Gosset-KEM Specification, version 1.0

**Status:** research specification, not a standard. Not approved by NIST. Do not use in production without independent review.
**Author:** Adrian Calabro · **Date:** 2026-10-29 · **Reference implementation:** https://github.com/calabroadrian/gosset-kem (`src/`).
**Name:** after Thorold Gosset, discoverer of the 4₂₁ polytope whose 240 vertices are the minimal vectors of E8 ("Gosset lattice"). Working name during research: AFA-12.

Gosset-KEM is a family of key-encapsulation mechanisms derived from **ML-KEM (FIPS 203)**. The only algorithmic changes are:

1. the message is encoded with the **E8 lattice code** (8 message bits per 8 coefficients) instead of 1 bit per coefficient;
2. the 32 E8 blocks are placed at the coefficient positions given by a fixed public table, **E8_LAYOUT**;
3. new parameter sets, including **η = 4** for the centered binomial distribution.

Everything not stated here is **exactly FIPS 203**: G, H, J, PRF, XOF, SampleNTT, SamplePolyCBD, NTT and NTT⁻¹, ByteEncode/ByteDecode, Compress/Decompress, the Fujisaki–Okamoto transform with implicit rejection (ML-KEM.KeyGen / Encaps / Decaps), the key and ciphertext formats, and the input checks of FIPS 203 §7.2–7.3. Section numbers such as "Alg. 14" refer to FIPS 203.

Notation: q = 3329, n = 256. Coefficients are integers mod q. For x ∈ ℤ_q, **cent(x)** is the representative in (−q/2, q/2], i.e. in [−1664, 1664].

---

## 1. Parameter sets

| Name | id | k | η₁ | η₂ | d_u | d_v | Message code | ek (B) | dk (B) | ct (B) |
|---|---|---|---|---|---|---|---|---|---|---|
| ML-KEM-512 | 0 | 2 | 3 | 2 | 10 | 4 | FIPS 203 | 800 | 1632 | 768 |
| ML-KEM-768 | 1 | 3 | 2 | 2 | 10 | 4 | FIPS 203 | 1184 | 2400 | 1088 |
| ML-KEM-1024 | 2 | 4 | 2 | 2 | 11 | 5 | FIPS 203 | 1568 | 3168 | 1568 |
| Gosset-KEM-512-Compact | 3 | 2 | 3 | 2 | **9** | 4 | E8 | 800 | 1632 | **704** |
| Gosset-KEM-768-Compact | 4 | 3 | 2 | 2 | **9** | **5** | E8 | 1184 | 2400 | **1024** |
| Gosset-KEM-1024-Compact | 5 | 4 | 2 | 2 | **10** | **4** | E8 | 1568 | 3168 | **1408** |
| Gosset-KEM-512-MaxSec | 6 | 2 | **4** | **3** | 10 | 4 | E8 | 800 | 1632 | 768 |
| Gosset-KEM-768-MaxSec | 7 | 3 | **3** | 2 | 10 | 4 | E8 | 1184 | 2400 | 1088 |
| Gosset-KEM-1024-MaxSec | 8 | 4 | **3** | 2 | 11 | 5 | E8 | 1568 | 3168 | 1568 |

Ids 0–2 are unmodified ML-KEM, included for comparison and cross-checks. Sizes are 384k + 32, 768k + 96 and 32(k·d_u + d_v) bytes, as in FIPS 203.

- **Compact** sets keep the ML-KEM noise, so they keep its lattice security, and spend the E8 decoding margin on smaller ciphertexts.
- **MaxSec** sets keep the ML-KEM sizes and spend the margin on more noise, which buys more lattice security.

## 2. Centered binomial distribution with η = 4

SamplePolyCBD_η (Alg. 8) is used unchanged with η = 4. PRF_η(s, b) returns 64·η = 256 bytes, and each coefficient is (Σ_{j<4} b_{8i+j}) − (Σ_{j<4} b_{8i+4+j}) over the bit array b. The reference implementation checks this sampler bit-exactly against a direct implementation of Alg. 8.

## 3. E8 message code

E8 = D8 ∪ (D8 + ½·𝟙), where D8 = {x ∈ ℤ⁸ : Σxᵢ even}. A message byte selects one of 256 coset representatives of E8 / 2E8; the code points are (q/2)·E8 reduced mod q.

### 3.1 Block layout

```
E8_LAYOUT[256] = {
      0,  18,  58, 159, 162, 197, 227, 228,  27,  42,  78, 112, 133, 145, 209, 235,
     43,  47,  54,  74,  97, 190, 231, 246,  96, 142, 194, 204, 217, 218, 226, 243,
     21,  37,  61,  84, 161, 166, 203, 244,  32,  71,  77,  79, 138, 157, 223, 252,
     30,  36,  73, 124, 168, 184, 195, 216,   1,  45,  59, 101, 104, 129, 186, 255,
      9,  10,  29,  44,  72, 105, 193, 248,  16,  57,  60,  70,  82, 180, 189, 245,
     39,  65, 103, 114, 153, 155, 242, 249,  11,  86, 128, 134, 232, 237, 241, 253,
     13, 122, 143, 151, 158, 182, 196, 251,   3,  20,  38,  48, 126, 139, 222, 254,
     12,  25,  31,  34,  81,  91, 111, 115,  26,  55,  89, 106, 167, 225, 229, 234,
     40,  56,  63,  99, 107, 132, 137, 207,  17,  19,  98, 108, 201, 202, 214, 233,
      8,  15,  35,  52, 120, 169, 221, 224,  53,  64, 117, 118, 152, 177, 179, 191,
     46,  92, 136, 172, 187, 192, 220, 250,  14,  24,  66, 100, 123, 131, 149, 208,
      6,  22,  51, 119, 154, 156, 160, 212,  85, 140, 164, 176, 183, 205, 230, 247,
      2,  62,  83,  94, 181, 185, 200, 238,   5,  50, 113, 121, 174, 175, 210, 213,
      7,  76, 102, 116, 148, 206, 215, 236,  23,  69,  87,  93, 109, 110, 165, 198,
     67, 125, 135, 144, 146, 170, 188, 219,   4,  49,  80, 150, 163, 173, 178, 211,
     68,  75,  90,  95, 130, 141, 239, 240,  28,  33,  41,  88, 127, 147, 171, 199 }
```

Block j (0 ≤ j < 32) uses the coefficients P_j = (E8_LAYOUT[8j], …, E8_LAYOUT[8j+7]), slot i ↦ coefficient E8_LAYOUT[8j+i].

E8_LAYOUT is a permutation of 0..255 with this property: **inside every block, the 28 pairwise index distances, taken as classes min(d, 256 − d), are all distinct.** This property is a security requirement (§6.3). An implementation MUST use exactly this table, because interoperability depends on it.

### 3.2 E8Encode(m) → μ ∈ ℤ_q²⁵⁶

Input: a 32-byte message m. For each block j = 0..31, with b = m[j]:

1. c ← bit 7 of b (the coset bit). xᵢ ← bit i of b for i = 0..6. x₇ ← x₀ ⊕ x₁ ⊕ … ⊕ x₆.
2. For i = 0..7: μ[E8_LAYOUT[8j+i]] ← 1665·xᵢ if c = 0, or 832 + 1665·xᵢ if c = 1.

(These are round((q/2)·(xᵢ + c/2)) with halves rounded up: 832.25 → 832, 1664.5 → 1665, 2496.75 → 2497.)

### 3.3 E8Decode(w) → m

Input: w ∈ ℤ_q²⁵⁶. For each block j = 0..31:

1. For i = 0..7: Yᵢ ← 4·cent(w[E8_LAYOUT[8j+i]]). Let U = 2q = 6658.
2. For each offset off ∈ {0, q} compute a candidate:
   1. fᵢ ← ⌊(Yᵢ − off)/U + ½⌋ (round half up) and eᵢ ← Yᵢ − off − U·fᵢ.
   2. If Σfᵢ is odd: let t be the **smallest** index with |e_t| = maxᵢ|eᵢ|; set f_t ← f_t + 1 if e_t ≥ 0, else f_t ← f_t − 1; recompute e_t.
   3. dist_off ← Σᵢ eᵢ².
3. If dist_q < dist_0 (strictly): output bit 7 = 1 and use the f of off = q. Otherwise output bit 7 = 0 and use the f of off = 0.
4. Output bit i (i = 0..6) = fᵢ mod 2.

This is nearest-point decoding in (q/2)·E8 (Conway–Sloane), using the scaled coordinates Y = 4y so that all quantities are integers. Ties are broken as stated, so the result is fully deterministic.

### 3.4 Changes inside K-PKE

- **K-PKE.Encrypt (Alg. 14):** replace μ ← Decompress₁(ByteDecode₁(m)) by **μ ← E8Encode(m)** for the E8 parameter sets.
- **K-PKE.Decrypt (Alg. 15):** replace m ← ByteEncode₁(Compress₁(w)) by **m ← E8Decode(w)**.

Everything else, including how d_u and d_v are applied, is unchanged.

## 4. Implementation requirements

1. **Constant time.** E8Decode depends on secret data. It MUST run without secret-dependent branches, memory indices or variable-latency instructions.
   - Compilers can insert branches into mask logic: gcc -O3 did so in a Leech decoder, and clang 18 -O3 in a reference E8 decoder.
   - Use value barriers and verify each build with a dudect-style test.
   - E8_LAYOUT gathers and scatters use public indices only.
2. **Masking.** Devices exposed to power or EM measurement need masking. A first-order masked decapsulation, including a bitsliced masked E8 decoder, is provided in `src/masked/` as a prototype.
   - The masked comparison selects the key on shares; only the selected key is unmasked.
   - The prototype uses a xorshift PRNG for masks. A production version MUST replace it with a CSPRNG.
3. **Key separation.** For the same seed d, ML-KEM-512 and Gosset-KEM-512-Compact produce **identical keys**: same k and η₁, and G(d‖k) does not encode the parameter set. The same holds for the other levels.
   - Applications MUST NOT reuse a key pair across parameter sets.
   - A future version may add a parameter-set byte to G's input.

## 5. Hybrid mode (optional)

pk = ek ‖ pk_X, sk = dk ‖ sk_X ‖ pk_X, ct = c ‖ ct_X, and

K = SHA3-256("GOSSET-KEM-HYBRID-v1" ‖ K_pq ‖ ss_X ‖ ct_X ‖ pk_X),

where X is X25519 (RFC 7748) and K_pq is the Gosset-KEM or ML-KEM shared key. The structure follows X-Wing.
- Mode C = ML-KEM-768 + X25519.
- Mode S = Gosset-KEM-768-MaxSec + X25519.

## 6. Security considerations

### 6.1 Lattice security
The underlying problem is Module-LWE with the stated noise. Estimates with the lattice-estimator (commit f31b86f, all attacks, MATZOV cost model) give the weakest link over the key and ciphertext instances, in log₂ of bit operations:

| Level | ML-KEM | Gosset-KEM MaxSec | Gain |
|---|---|---|---|
| 512 | 139.7 | 144.9 | +5.2 |
| 768 | 196.4 | 207.2 | +10.8 |
| 1024 | 262.3 | 271.7 | +9.4 |

These figures credit the ciphertext rounding noise. Without that credit the gain is +4.8 to +5.6 bits. Compact sets equal ML-KEM.

### 6.2 Decryption failure rate (average ciphertext)
Exact convolution with a union bound over the 240 Voronoi-relevant vectors of E8, validated on real encryptions at inflated noise (measured/predicted ratio 1.00 ± 0.05):

| Set | 512 | 768 | 1024 |
|---|---|---|---|
| ML-KEM | 2⁻¹³⁹ | 2⁻¹⁶⁵ | 2⁻¹⁷⁵ |
| Gosset-KEM Compact | 2⁻¹⁴² | 2⁻¹⁹⁴ | 2⁻²¹⁵ |
| Gosset-KEM MaxSec | 2⁻¹⁴⁷ | 2⁻¹⁶⁷ | 2⁻¹⁸⁹ |

### 6.3 Failure-boosting attacks and why E8_LAYOUT exists
A multi-coordinate decoder fails along directions u that combine several coefficients. The decryption errors δ_a, δ_b of two coefficients have covariance σ²·A(b − a), where A is the negacyclic autocorrelation of the attacker-known ciphertext randomness. An attacker who grinds messages can select ciphertexts with large A at chosen distances, which raises their failure probability along E8 directions. The standard 1-bit decoder is immune to this effect.

- With **consecutive** blocks, one distance appears in up to 7 pairs per block. The first-failure cost (Q ≤ 2⁶⁴ queries) then falls to **2¹²⁷** for Gosset-KEM-512-Compact, below its claimed level.
- E8_LAYOUT makes the 28 distances in each block distinct. The first-failure cost becomes **2¹⁹² / 2²³² (512 Compact / MaxSec), 2³⁶⁵ / 2³⁸⁷ (768) and 2⁷⁹¹ / 2⁵⁵² (1024)**. These come from a saddlepoint tail that matches Monte Carlo within 15% down to 10⁻⁴.
- After a first failure, directional boosting makes later failures cheaper in every scheme, but they remain costlier for Gosset-KEM than for ML-KEM (2nd failure: 2¹⁴⁵ vs 2⁸⁵ at level 512). The first failure dominates.
- **Codes of higher dimension are more exposed.** For the Leech lattice (24 coordinates), no layout can make all 276 distances distinct. The best layout found leaves an exposure of 2.14, against 1.32 for E8 with E8_LAYOUT.

### 6.4 Not covered
- Formal proofs beyond the FIPS 203 Fujisaki–Okamoto argument with the stated failure rates.
- Hardware side-channel validation.
- Multi-target attacks, which are excluded as for ML-KEM.

## 7. Test vectors

- `kat/GOSSET_KEM_KAT.json`: 45 KEM vectors (5 per parameter set) and 10 hybrid vectors. Seeds are SHAKE256("GOSSET-KEM-KAT-v1|label|i"). SHA3-256 of the file: `4b1c0e3761c8b2415fd2fe5ecf75bb10fc4f67a3787807fad88ea6cfa963de0f`.
- The 15 ML-KEM entries equal FIPS 203. `kat/keygen.json` and `kat/encdec.json` hold the NIST ACVP vectors, which pass 180/180, plus 60/60 key checks.
- `tests/xplat_kat.c`: a cross-platform digest over all schemes, the hybrid modes, the DEM and the decoders, `80254431656df5bfd7fa9f150bc04577`.

## 8. Version history

| Version | Change |
|---|---|
| 1.0 (2026-10-29) | First specification, published as Gosset-KEM (research builds used the working name AFA-12 and different domain-separation labels). It includes E8_LAYOUT. Earlier research builds (≤ v12) used consecutive blocks; they are not compatible with 1.0 and are insecure against correlation-based failure boosting at level 512. |
