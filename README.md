# Gosset-KEM

**E8 message coding for ML-KEM (FIPS 203), with a correlation-based failure-boosting attack on multi-dimensional lattice codes and its countermeasure.**

> Research prototype. Not standardised, not audited. Do not use in production.

Named after Thorold Gosset, who discovered the 4₂₁ polytope whose 240 vertices are the minimal vectors of E8, the "Gosset lattice". Developed under the working name **AFA-12**, which still appears in `research/`.

Gosset-KEM replaces the 1-bit-per-coefficient message encoding of ML-KEM with the 8-dimensional Gosset lattice **E8**. The decoding margin can be spent in two ways:

| Mode | What you get | Cost |
|---|---|---|
| **Compact** | Ciphertexts 6–10 % smaller, same lattice security as ML-KEM | ≤ 0.9 % decapsulation time |
| **MaxSec** | +5 to +11 bits against the best lattice attacks (lattice-estimator, MATZOV), same sizes as ML-KEM | +2 to +7 % time |

While analysing failure-boosting attacks we found that multi-dimensional codes such as E8 and Leech open a new attack surface. An attacker who chooses ciphertexts with correlated decryption errors inside a block raises their failure probability, while the standard decoder is immune. With consecutive blocks this puts the level-1 Compact set below its claimed security (2^127). A fixed public coefficient layout with distinct in-block index distances (**E8_LAYOUT**) restores a large margin (2^192) at no cost. Details are in the paper.

## Contents

| Path | What |
|---|---|
| `docs/paper/gosset-kem.pdf` | Paper (source: `gosset-kem.tex`) |
| `docs/SPEC.md` | **Specification v1.0**: exact algorithms, parameter sets, E8_LAYOUT, test vectors |
| `docs/AUDIT_PACKAGE.md` | Claims, evidence, known risks, questions for reviewers |
| `docs/STATUS_es.md` | Full research log (Spanish) |
| `src/` | C implementation: ML-KEM + Gosset-KEM (AVX2, constant time), FIPS 202, X25519, hybrid KEM, DEM |
| `src/masked/` | First-order masked decapsulation (bitsliced masked E8 decoder, masked Keccak) |
| `src/ref/` | Independent reference build (bit-level CBD), used for cross-checks |
| `tests/` | Test suite |
| `bench/` | Benchmarks and dudect-style constant-time tests |
| `kat/` | NIST ACVP vectors for ML-KEM and Gosset-KEM known-answer tests |
| `research/` | All analysis scripts, v2 to v14, to reproduce every number |

## Build and test

Requirements: Linux or macOS on x86-64 with AVX2, gcc or clang, Python 3.

```sh
make test           # NIST ACVP, Gosset-KEM KATs, reference == optimized, masked == unmasked, X25519, cross-platform digest
make bench          # cycles and constant-time tests
make CC=clang test  # same with clang
```

Expected final line: `RESULT: ALL TESTS PASSED`.

## Parameter sets

| Name | k | η₁ | η₂ | d_u | d_v | ct (bytes) | DFR | Lattice security (log2, weakest link) |
|---|---|---|---|---|---|---|---|---|
| Gosset-KEM-512-Compact | 2 | 3 | 2 | 9 | 4 | 704 | 2^-142 | = ML-KEM-512 |
| Gosset-KEM-768-Compact | 3 | 2 | 2 | 9 | 5 | 1024 | 2^-194 | = ML-KEM-768 |
| Gosset-KEM-1024-Compact | 4 | 2 | 2 | 10 | 4 | 1408 | 2^-215 | = ML-KEM-1024 |
| Gosset-KEM-512-MaxSec | 2 | 4 | 3 | 10 | 4 | 768 | 2^-147 | 144.9 (ML-KEM 139.7) |
| Gosset-KEM-768-MaxSec | 3 | 3 | 2 | 10 | 4 | 1088 | 2^-167 | 207.2 (ML-KEM 196.4) |
| Gosset-KEM-1024-MaxSec | 4 | 3 | 2 | 11 | 5 | 1568 | 2^-189 | 271.7 (ML-KEM 262.3) |

## Citation

Repository: https://github.com/calabroadrian/gosset-kem

See `CITATION.cff`. Author: Adrian Calabro.

## License

Apache License 2.0 (`LICENSE`). The NIST ACVP test vectors in `kat/` are public-domain works of the U.S. Government.

---

### Resumen en español

Gosset-KEM cambia la forma en que ML-KEM, el estándar post-cuántico de EE.UU., empaqueta el mensaje. En lugar de 1 bit por coeficiente usa la geometría E8, que tolera más ruido. Ese margen se usa para dos cosas: mensajes cifrados más chicos (modo Compact) o más seguridad con los mismos tamaños (modo MaxSec, entre +5 y +11 bits).

El aporte principal del proyecto es haber encontrado un ataque nuevo contra este tipo de códigos y su defensa. El atacante elige mensajes cifrados cuyos errores están correlacionados dentro de un bloque. Con bloques consecutivos, el nivel 512 queda por debajo de su seguridad declarada. Con la distribución pública E8_LAYOUT, el margen vuelve a ser amplio y sin costo extra.
