# Research history (v2 – v14)

The project was developed under the working name **AFA-12**; scripts and outputs here keep that name and the research-era domain labels ("AFA12-…"). The published scheme is **Gosset-KEM** (`../docs/SPEC.md`).

Every analysis behind the paper and the specification, in the order it was done. The production code is in `../src`; the folders here are the historical working copies, kept so that every number can be reproduced.

Build the shared libraries the scripts load with `make` (in this folder), then run each script **from its own folder**, e.g. `cd v13 && python3 failure_boost.py`. Python 3 with numpy and scipy is required; `v12/official_estimate.py` also needs Sage (e.g. `pip install passagemath-standard`) and a clone of https://github.com/malb/lattice-estimator inside `v12/`.

| Folder | Topic | Main result |
|---|---|---|
| v2 | E8-structured **noise** (original idea) | Negative: loses log2(index) bits; confirmed with real LLL attacks |
| v3 | E8 **message** coding, exact DFR, prior work, validation with real encryptions | DFR 2^-241 / -279 / -335 at ML-KEM parameters |
| v4 | Constant-time E8 decoder in C (scalar, AVX2), dudect | ~1,300 cycles AVX2 |
| v5 | Leech lattice code | Smaller ciphertexts, but ~850k cycles |
| v6 | Noise sweep for the "v2.0" design | sigma = 32.5 never decrypts; valid window found |
| v7 | Full ML-KEM + AFA-12 implementation, NIST ACVP, DEM | 180/180 NIST vectors |
| v8 | MaxSec parameter search | + security at FIPS sizes |
| v9 | X25519, hybrid KEM, power simulation, KATs, audit package | Mode C / Mode S |
| v10 | First-order masked E8 decoder | TVLA first order clean (simulated) |
| v11 | Bitsliced masked decoder, full masked decapsulation | 12.3x faster, 720/720 exact |
| v12 | Security with core-SVP and the official lattice-estimator | +5 to +11 bits for MaxSec |
| v13 | Failure boosting, **correlation attack**, layout countermeasure | 2^127 -> 2^192 at level 512 |
| v14 | Saddlepoint tail, directional boosting, Leech exposure | Conclusions confirmed |
| local | Windows / AMD Zen 3 test kit and timing diagnostics | Cross-platform identical; Zen 3 artefact classified |

Note: folders up to v12 use consecutive E8 blocks; the E8_LAYOUT of specification 1.0 was introduced in v13.
