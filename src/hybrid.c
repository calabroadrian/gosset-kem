/*
 * Gosset-KEM v9 — Hybrid KEM: (ML-KEM or Gosset-KEM) + X25519, X-Wing-style combiner
 * ===========================================================================
 *   pk = ek_pq || pk_x            sk = dk_pq || sk_x || pk_x            ct = c_pq || ct_x
 *   K  = SHA3-256("GOSSET-KEM-HYBRID-v1" || K_pq || ss_x || ct_x || pk_x)
 * The shared key stays secret if EITHER component is unbroken (lattice OR elliptic curve).
 * Modes:  id 1 = ML-KEM-768 (FIPS 203)  -> "Mode C", compliance-friendly
 *         id 7 = Gosset-KEM-768-MaxSec        -> "Mode S", maximum security
 * Public API draws randomness from the OS (getrandom) and wipes secrets.
 */
#include <stdint.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#include <bcrypt.h>
#else
#include <sys/random.h>
#endif
#include "fips202.h"
#include "mlkem.h"
#include "x25519.h"
#include "hybrid.h"

static const char LABEL[] = "GOSSET-KEM-HYBRID-v1";

void secure_zero(void *p, size_t n) { volatile uint8_t *v = p; while (n--) *v++ = 0; }
static int os_random(uint8_t *b, size_t n) {
#ifdef _WIN32
    return BCryptGenRandom(NULL, b, (ULONG)n, BCRYPT_USE_SYSTEM_PREFERRED_RNG) == 0 ? 0 : -1;
#else
    while (n) { ssize_t r = getrandom(b, n, 0); if (r <= 0) return -1; b += r; n -= (size_t)r; }
    return 0;
#endif
}

size_t hyb_pk_bytes(int id) { return mlkem_ek_bytes(mlkem_get(id)) + 32; }
size_t hyb_sk_bytes(int id) { return mlkem_dk_bytes(mlkem_get(id)) + 64; }
size_t hyb_ct_bytes(int id) { return mlkem_ct_bytes(mlkem_get(id)) + 32; }

static void combine(uint8_t K[32], const uint8_t Kpq[32], const uint8_t ssx[32], const uint8_t ctx[32], const uint8_t pkx[32]) {
    uint8_t in[sizeof LABEL - 1 + 128];
    memcpy(in, LABEL, sizeof LABEL - 1);
    memcpy(in + sizeof LABEL - 1, Kpq, 32); memcpy(in + sizeof LABEL - 1 + 32, ssx, 32);
    memcpy(in + sizeof LABEL - 1 + 64, ctx, 32); memcpy(in + sizeof LABEL - 1 + 96, pkx, 32);
    sha3_256(K, in, sizeof in);
    secure_zero(in, sizeof in);
}

/* deterministic cores (for test vectors) */
void hyb_keygen_internal(int id, const uint8_t d[32], const uint8_t z[32], const uint8_t skx[32], uint8_t *pk, uint8_t *sk) {
    size_t ekl = mlkem_ek_bytes(mlkem_get(id)), dkl = mlkem_dk_bytes(mlkem_get(id));
    uint8_t pkx[32];
    mlkem_keygen_internal(id, d, z, pk, sk);
    x25519_base(pkx, skx);
    memcpy(pk + ekl, pkx, 32);
    memcpy(sk + dkl, skx, 32); memcpy(sk + dkl + 32, pkx, 32);
}
void hyb_encaps_internal(int id, const uint8_t *pk, const uint8_t m[32], const uint8_t ex[32], uint8_t *ct, uint8_t K[32]) {
    size_t ekl = mlkem_ek_bytes(mlkem_get(id)), cl = mlkem_ct_bytes(mlkem_get(id));
    const uint8_t *pkx = pk + ekl;
    uint8_t Kpq[32], ssx[32];
    mlkem_encaps_internal(id, pk, m, Kpq, ct);
    x25519_base(ct + cl, ex);
    x25519(ssx, ex, pkx);
    combine(K, Kpq, ssx, ct + cl, pkx);
    secure_zero(Kpq, 32); secure_zero(ssx, 32);
}
void hyb_decaps(int id, const uint8_t *sk, const uint8_t *ct, uint8_t K[32]) {
    size_t dkl = mlkem_dk_bytes(mlkem_get(id)), cl = mlkem_ct_bytes(mlkem_get(id));
    uint8_t Kpq[32], ssx[32];
    mlkem_decaps(id, sk, ct, Kpq);
    x25519(ssx, sk + dkl, ct + cl);
    combine(K, Kpq, ssx, ct + cl, sk + dkl + 32);
    secure_zero(Kpq, 32); secure_zero(ssx, 32);
}

/* public API with OS randomness and input checks */
int hyb_keygen(int id, uint8_t *pk, uint8_t *sk) {
    uint8_t seed[96];
    if (os_random(seed, 96)) return -1;
    hyb_keygen_internal(id, seed, seed + 32, seed + 64, pk, sk);
    secure_zero(seed, 96);
    return 0;
}
int hyb_encaps(int id, const uint8_t *pk, uint8_t *ct, uint8_t K[32]) {
    if (!mlkem_check_ek(id, pk)) return -1;                  /* FIPS 203 input check */
    uint8_t r[64];
    if (os_random(r, 64)) return -1;
    hyb_encaps_internal(id, pk, r, r + 32, ct, K);
    secure_zero(r, 64);
    return 0;
}
int hyb_decaps_checked(int id, const uint8_t *sk, const uint8_t *ct, uint8_t K[32]) {
    if (!mlkem_check_dk(id, sk)) return -1;                  /* FIPS 203 input check */
    hyb_decaps(id, sk, ct, K);
    return 0;
}
