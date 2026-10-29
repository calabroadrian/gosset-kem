/* AFA-12 v7 — end-to-end demo: AFA12-768-MaxSec KEM + corrected DEM (diagram V6 flow) */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mlkem.h"
#include "dem.h"

static void rnd(uint8_t *b, size_t n) { FILE *f = fopen("/dev/urandom", "rb"); if (!f || fread(b, 1, n, f) != n) exit(1); fclose(f); }
static void hex(const char *l, const uint8_t *b, size_t n) { printf("  %-22s", l); for (size_t i = 0; i < n && i < 16; i++) printf("%02x", b[i]); printf("%s\n", n > 16 ? "..." : ""); }

int main(void) {
    const int id = 7;                                   /* AFA12-768-MaxSec */
    const mlkem_params *P = mlkem_get(id);
    uint8_t d[32], z[32], m[32], ek[1600], dk[3200], c[1600], K[32], K2[32], nonce[16], tag[32];
    const char *msg = "Orden de Operacion / Datos Sensibles";
    uint8_t C[256], M[256];
    rnd(d, 32); rnd(z, 32); rnd(m, 32); rnd(nonce, 16);

    printf("== %s  (ek %zu B, dk %zu B, ct %zu B)\n", P->name, mlkem_ek_bytes(P), mlkem_dk_bytes(P), mlkem_ct_bytes(P));
    mlkem_keygen_internal(id, d, z, ek, dk);
    printf("  receiver key checks: ek %s, dk %s\n", mlkem_check_ek(id, ek) ? "valid" : "INVALID", mlkem_check_dk(id, dk) ? "valid" : "INVALID");
    mlkem_encaps_internal(id, ek, m, K, c);             /* sender: Leyes 1,2,3,5 */
    size_t mlen = strlen(msg), plen = dem_padded_len(mlen);
    dem_encrypt(K, nonce, (const uint8_t *)msg, mlen, C, tag);   /* Ley 4 + DEM */
    hex("KEM ciphertext c", c, mlkem_ct_bytes(P)); hex("DEM ciphertext C", C, plen); hex("tag", tag, 32);
    printf("  message %zu B -> padded/encrypted %zu B (length hidden to 64-byte blocks)\n", mlen, plen);

    mlkem_decaps(id, dk, c, K2);                        /* receiver: Ley 6 + FO */
    long r = dem_decrypt(K2, nonce, C, plen, tag, M);
    printf("  [legit]          K == K' : %s | DEM: %s\n", memcmp(K, K2, 32) ? "NO" : "yes",
           r < 0 ? "REJECTED" : "authenticated");
    if (r >= 0) printf("  recovered: \"%.*s\"\n", (int)r, M);

    uint8_t cbad[1600]; memcpy(cbad, c, mlkem_ct_bytes(P)); cbad[100] ^= 1;
    mlkem_decaps(id, dk, cbad, K2);
    r = dem_decrypt(K2, nonce, C, plen, tag, M);
    printf("  [tampered KEM ct] K == K' : %s | DEM: %s  (implicit rejection -> unrelated key)\n",
           memcmp(K, K2, 32) ? "no" : "YES", r < 0 ? "rejected" : "ACCEPTED");

    mlkem_decaps(id, dk, c, K2); C[5] ^= 0x20;
    r = dem_decrypt(K2, nonce, C, plen, tag, M);
    printf("  [tampered DEM C]  DEM: %s\n", r < 0 ? "rejected before decryption" : "ACCEPTED");
    return 0;
}
