/* AFA-12 cross-platform equivalence test: deterministic KEM / hybrid / DEM / E8 runs hashed per scheme.
   The digests must be identical on every platform and compiler. */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "../v7/fips202.h"
#include "../v7/mlkem.h"
#include "../v7/dem.h"
#include "../v9/hybrid.h"
void e8_encode(const uint8_t msg[32], int16_t out[256]);
void e8_decode_avx2(const int16_t in[256], uint8_t msg[32]);
void e8_decode_fast(const int16_t in[256], uint8_t msg[32]);

static void seed(uint8_t *o, size_t n, const char *lab, int id, int i) {
    uint8_t in[64]; int l = snprintf((char *)in, sizeof in, "%s/%d/%d", lab, id, i); shake256(o, n, in, (size_t)l);
}
static void hex(const uint8_t *h) { for (int i = 0; i < 16; i++) printf("%02x", h[i]); }

int main(void) {
    static uint8_t ek[1600], dk[3200], c[1600], pk[1700], sk[3300], ct[1700];
    uint8_t d[32], z[32], m[32], K[32], K2[32], e[32], x[32], all[32];
    keccak_ctx A; keccak_init(&A, 136, 0x06);
    int fails = 0;
    for (int id = 0; id < 9; id++) {
        const mlkem_params *P = mlkem_get(id);
        keccak_ctx H; keccak_init(&H, 136, 0x06);
        for (int i = 0; i < 25; i++) {
            seed(d, 32, "d", id, i); seed(z, 32, "z", id, i); seed(m, 32, "m", id, i);
            mlkem_keygen_internal(id, d, z, ek, dk);
            mlkem_encaps_internal(id, ek, m, K, c);
            mlkem_decaps(id, dk, c, K2);
            if (memcmp(K, K2, 32) || !mlkem_check_ek(id, ek) || !mlkem_check_dk(id, dk)) fails++;
            keccak_absorb(&H, ek, mlkem_ek_bytes(P)); keccak_absorb(&H, dk, mlkem_dk_bytes(P));
            keccak_absorb(&H, c, mlkem_ct_bytes(P)); keccak_absorb(&H, K, 32);
            c[i % mlkem_ct_bytes(P)] ^= 1; mlkem_decaps(id, dk, c, K2);            /* implicit rejection */
            if (!memcmp(K, K2, 32)) fails++;
            keccak_absorb(&H, K2, 32);
        }
        uint8_t h[32]; keccak_finalize(&H); keccak_squeeze(&H, h, 32); keccak_absorb(&A, h, 32);
        printf("%-20s ", P->name); hex(h); printf("\n");
    }
    { /* hybrid modes C and S */
        int ids[2] = {1, 7};
        for (int t = 0; t < 2; t++) {
            keccak_ctx H; keccak_init(&H, 136, 0x06); int id = ids[t];
            for (int i = 0; i < 25; i++) {
                seed(d, 32, "hd", id, i); seed(z, 32, "hz", id, i); seed(x, 32, "hx", id, i);
                seed(m, 32, "hm", id, i); seed(e, 32, "he", id, i);
                hyb_keygen_internal(id, d, z, x, pk, sk); hyb_encaps_internal(id, pk, m, e, ct, K); hyb_decaps(id, sk, ct, K2);
                if (memcmp(K, K2, 32)) fails++;
                keccak_absorb(&H, pk, hyb_pk_bytes(id)); keccak_absorb(&H, ct, hyb_ct_bytes(id)); keccak_absorb(&H, K, 32);
            }
            uint8_t h[32]; keccak_finalize(&H); keccak_squeeze(&H, h, 32); keccak_absorb(&A, h, 32);
            printf("%-20s ", t ? "Hybrid-S" : "Hybrid-C"); hex(h); printf("\n");
        }
    }
    { /* DEM and E8 decoders on noisy inputs */
        keccak_ctx H; keccak_init(&H, 136, 0x06);
        static uint8_t M[1000], C[1100], M2[1100]; uint8_t nonce[16], tag[32];
        for (int i = 0; i < 40; i++) {
            size_t ml = (size_t)(i * 25); seed(K, 32, "dk", 0, i); seed(nonce, 16, "dn", 0, i); seed(M, ml ? ml : 1, "dm", 0, i);
            dem_encrypt(K, nonce, M, ml, C, tag); size_t cl = dem_padded_len(ml);
            if (dem_decrypt(K, nonce, C, cl, tag, M2) != (long)ml || memcmp(M, M2, ml)) fails++;
            C[0] ^= 1; if (dem_decrypt(K, nonce, C, cl, tag, M2) >= 0) fails++; C[0] ^= 1;
            keccak_absorb(&H, C, cl); keccak_absorb(&H, tag, 32);
            int16_t w[256]; uint8_t r1[32], r2[32], noise[512];
            seed(m, 32, "em", 0, i); seed(noise, 512, "en", 0, i); e8_encode(m, w);
            for (int j = 0; j < 256; j++) { int v = w[j] + ((int)(noise[2*j] | noise[2*j+1] << 8) % 1201) - 600; w[j] = (int16_t)(((v % 3329) + 3329) % 3329); }
            e8_decode_avx2(w, r1); e8_decode_fast(w, r2); if (memcmp(r1, r2, 32)) fails++;
            keccak_absorb(&H, r1, 32);
        }
        uint8_t h[32]; keccak_finalize(&H); keccak_squeeze(&H, h, 32); keccak_absorb(&A, h, 32);
        printf("%-20s ", "DEM+E8"); hex(h); printf("\n");
    }
    keccak_finalize(&A); keccak_squeeze(&A, all, 32);
    printf("GLOBAL               "); hex(all); printf("\ninternal failures: %d\n", fails);
    return fails != 0;
}
