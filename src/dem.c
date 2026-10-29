/*
 * Gosset-KEM v7 — DEM layer (Ley 4 + capa DEM), corrected
 * ===================================================
 * Fixes the V6 design, which used the same KEM key K for the keystream and
 * for the MAC. Now:
 *   (K_enc || K_mac) = SHAKE256("GOSSET-KEM-DEM-KDF-v1" || K, 64)        key separation
 *   M_pad            = M || 0x80 || 0x00...  (ISO/IEC 7816-4, multiple of 64 B)  Ley 4
 *   C                = M_pad XOR SHAKE256("GOSSET-KEM-DEM-ENC-v1" || K_enc || nonce, |M_pad|)
 *   tag              = HMAC-SHA3-256(K_mac, "GOSSET-KEM-DEM-MAC-v1" || nonce || C)
 * Decryption verifies the tag in constant time BEFORE decrypting/unpadding
 * (encrypt-then-MAC: no padding or decryption oracle). Nonce: 16 bytes, unique per K.
 */
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "fips202.h"
#include "dem.h"

#define RATE 136   /* SHA3-256 block size */

void hmac_sha3_256(const uint8_t *key, size_t klen, const uint8_t *msg, size_t mlen, uint8_t out[32]) {
    uint8_t k0[RATE] = {0}, ipad[RATE], opad[RATE], inner[32];
    if (klen > RATE) sha3_256(k0, key, klen); else memcpy(k0, key, klen);
    for (int i = 0; i < RATE; i++) { ipad[i] = k0[i] ^ 0x36; opad[i] = k0[i] ^ 0x5c; }
    keccak_ctx c;
    keccak_init(&c, RATE, 0x06); keccak_absorb(&c, ipad, RATE); keccak_absorb(&c, msg, mlen);
    keccak_finalize(&c); keccak_squeeze(&c, inner, 32);
    keccak_init(&c, RATE, 0x06); keccak_absorb(&c, opad, RATE); keccak_absorb(&c, inner, 32);
    keccak_finalize(&c); keccak_squeeze(&c, out, 32);
}

static void derive(const uint8_t K[32], uint8_t kenc[32], uint8_t kmac[32]) {
    static const char lbl[] = "GOSSET-KEM-DEM-KDF-v1";
    uint8_t in[sizeof lbl - 1 + 32], out[64];
    memcpy(in, lbl, sizeof lbl - 1); memcpy(in + sizeof lbl - 1, K, 32);
    shake256(out, 64, in, sizeof in);
    memcpy(kenc, out, 32); memcpy(kmac, out + 32, 32);
}

static void keystream_xor(const uint8_t kenc[32], const uint8_t nonce[16], uint8_t *buf, size_t len) {
    static const char lbl[] = "GOSSET-KEM-DEM-ENC-v1";
    keccak_ctx c; uint8_t ks[RATE];
    keccak_init(&c, RATE, 0x1F);
    keccak_absorb(&c, (const uint8_t *)lbl, sizeof lbl - 1); keccak_absorb(&c, kenc, 32); keccak_absorb(&c, nonce, 16);
    keccak_finalize(&c);
    for (size_t off = 0; off < len; off += RATE) {
        size_t n = len - off < RATE ? len - off : RATE;
        keccak_squeeze(&c, ks, n);
        for (size_t i = 0; i < n; i++) buf[off + i] ^= ks[i];
    }
}

static void mac(const uint8_t kmac[32], const uint8_t nonce[16], const uint8_t *C, size_t clen, uint8_t tag[32]) {
    static const char lbl[] = "GOSSET-KEM-DEM-MAC-v1";
    size_t L = sizeof lbl - 1 + 16 + clen;
    uint8_t *m = malloc(L);
    memcpy(m, lbl, sizeof lbl - 1); memcpy(m + sizeof lbl - 1, nonce, 16); memcpy(m + sizeof lbl - 1 + 16, C, clen);
    hmac_sha3_256(kmac, 32, m, L, tag);
    free(m);
}

size_t dem_padded_len(size_t mlen) { return (mlen / 64 + 1) * 64; }

void dem_encrypt(const uint8_t K[32], const uint8_t nonce[16], const uint8_t *M, size_t mlen, uint8_t *C, uint8_t tag[32]) {
    uint8_t kenc[32], kmac[32]; derive(K, kenc, kmac);
    size_t plen = dem_padded_len(mlen);
    memcpy(C, M, mlen); C[mlen] = 0x80; memset(C + mlen + 1, 0, plen - mlen - 1);   /* ISO/IEC 7816-4 */
    keystream_xor(kenc, nonce, C, plen);
    mac(kmac, nonce, C, plen, tag);
}

/* returns plaintext length, or -1 if authentication fails (nothing is decrypted) */
long dem_decrypt(const uint8_t K[32], const uint8_t nonce[16], const uint8_t *C, size_t clen, const uint8_t tag[32], uint8_t *M) {
    if (clen == 0 || clen % 64) return -1;
    uint8_t kenc[32], kmac[32], t2[32]; derive(K, kenc, kmac);
    mac(kmac, nonce, C, clen, t2);
    uint8_t d = 0; for (int i = 0; i < 32; i++) d |= (uint8_t)(t2[i] ^ tag[i]);
    if (d) return -1;                                        /* reject before decrypting */
    memcpy(M, C, clen); keystream_xor(kenc, nonce, M, clen);
    long i = (long)clen - 1;
    while (i >= 0 && M[i] == 0x00) i--;                      /* authenticated data: no padding oracle */
    if (i < 0 || M[i] != 0x80) return -1;
    return i;
}
