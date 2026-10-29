/* AFA-12 — localize a decapsulation timing difference (Welch t per stage and per input class).
 * Tests on one parameter set (default AFA12-512-MaxSec, id 6):
 *  A  full decaps: fixed valid ct  vs random ct               (original test)
 *  B  full decaps: fixed valid ct  vs random VALID ct         (accept path in both classes)
 *  C  full decaps: random valid ct vs random ct               (valid vs invalid, no fixed class)
 *  D  full decaps: fixed valid ct #1 vs fixed valid ct #2     (two fixed inputs)
 *  S1..S6 stages on fixed vs random inputs: pke_decrypt, G=sha3_512, J=shake256, pke_encrypt, compare+select, cbd(eta1)
 * Output: max|t| and the mean cycle difference per class. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <x86intrin.h>
#include "../v7/mlkem_opt.c"

static uint64_t rs = 0x853C49E6748FEA9BULL;
static inline uint64_t xr(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static void rnd(uint8_t *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = (uint8_t)xr(); }
static inline uint64_t cyc(void) { unsigned a; _mm_lfence(); uint64_t t = __rdtscp(&a); _mm_lfence(); return t; }
static int cmpu(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y); }

static double welch(const uint64_t *v, const int *c, int n, double thr, double *dm) {
    double s[2] = {0}, q[2] = {0}, k[2] = {0};
    for (int i = 0; i < n; i++) { double x = (double)v[i]; if (x > thr) continue; s[c[i]] += x; q[c[i]] += x * x; k[c[i]]++; }
    double m0 = s[0] / k[0], m1 = s[1] / k[1]; *dm = m0 - m1;
    return (m0 - m1) / sqrt((q[0] / k[0] - m0 * m0) / k[0] + (q[1] / k[1] - m1 * m1) / k[1]);
}
static void report(const char *name, uint64_t *t, int *cl, int n) {
    uint64_t *so = malloc(8 * (size_t)n); memcpy(so, t, 8 * (size_t)n); qsort(so, n, 8, cmpu);
    double p[] = {0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0}, w = 0, dmw = 0;
    for (int k = 0; k < 7; k++) { double dm, tv = welch(t, cl, n, p[k] >= 1 ? 1e300 : (double)so[(int)(p[k] * n)], &dm); if (fabs(tv) > fabs(w)) { w = tv; dmw = dm; } }
    printf("  %-52s max|t| = %6.2f  mean diff %+8.1f cycles (median %llu)  %s\n", name, fabs(w), dmw, (unsigned long long)so[n / 2], fabs(w) > 4.5 ? "LEAK" : "ok");
    free(so); fflush(stdout);
}

#define NB 256
int main(int argc, char **argv) {
    int id = argc > 1 ? atoi(argv[1]) : 6, n = argc > 2 ? atoi(argv[2]) : 500000;
    const mlkem_params *P = &PARAMS[id]; size_t L = mlkem_ct_bytes(P); int k = P->k;
    static uint8_t ek[1600], dk[3200], f1[1600], f2[1600], K[32], d[32], z[32], m[32];
    rnd(d, 32); rnd(z, 32); mlkem_keygen_internal(id, d, z, ek, dk);
    rnd(m, 32); mlkem_encaps_internal(id, ek, m, K, f1);
    rnd(m, 32); mlkem_encaps_internal(id, ek, m, K, f2);
    uint64_t *t = malloc(8 * (size_t)n); int *cl = malloc(4 * (size_t)n);
    static uint8_t cb[NB][1600], mb[NB][32], rb[NB][32];
    printf("DIAGNOSTIC %s, n=%d per test\n", P->name, n);
    for (int test = 0; test < 4; test++) {
        for (int b0 = 0; b0 < n; b0 += NB) {
            int mm = n - b0 < NB ? n - b0 : NB;
            for (int j = 0; j < mm; j++) {
                int c = (int)(xr() & 1); cl[b0 + j] = c;
                if (test == 0) { if (!c) memcpy(cb[j], f1, L); else rnd(cb[j], L); }
                if (test == 1) { if (!c) memcpy(cb[j], f1, L); else { rnd(m, 32); mlkem_encaps_internal(id, ek, m, K, cb[j]); } }
                if (test == 2) { if (!c) { rnd(m, 32); mlkem_encaps_internal(id, ek, m, K, cb[j]); } else rnd(cb[j], L); }
                if (test == 3) memcpy(cb[j], c ? f2 : f1, L);
            }
            for (int j = 0; j < mm; j++) { uint64_t a = cyc(); mlkem_decaps(id, dk, cb[j], K); uint64_t e = cyc(); t[b0 + j] = e - a; }
        }
        const char *nm[] = {"A decaps: fixed valid vs random", "B decaps: fixed valid vs random valid",
                            "C decaps: random valid vs random invalid", "D decaps: fixed valid #1 vs fixed valid #2"};
        report(nm[test], t, cl, n);
    }
    /* stage tests: class 0 fixed input, class 1 random input */
    uint8_t mfix[32], rfix[32], g[64], in[64], cp[1600], jin[32 + 1600], kb[32];
    pke_decrypt(P, dk, f1, mfix); { uint8_t i2[64]; memcpy(i2, mfix, 32); memcpy(i2 + 32, dk + 768 * k + 32, 32); sha3_512(g, i2, 64); memcpy(rfix, g + 32, 32); }
    for (int st = 0; st < 6; st++) {
        for (int b0 = 0; b0 < n; b0 += NB) {
            int mm = n - b0 < NB ? n - b0 : NB;
            for (int j = 0; j < mm; j++) {
                int c = (int)(xr() & 1); cl[b0 + j] = c;
                if (!c) { memcpy(cb[j], f1, L); memcpy(mb[j], mfix, 32); memcpy(rb[j], rfix, 32); }
                else { rnd(cb[j], L); rnd(mb[j], 32); rnd(rb[j], 32); }
            }
            for (int j = 0; j < mm; j++) {
                uint64_t a = cyc();
                switch (st) {
                case 0: pke_decrypt(P, dk, cb[j], m); break;
                case 1: memcpy(in, mb[j], 32); memcpy(in + 32, dk + 768 * k + 32, 32); sha3_512(g, in, 64); break;
                case 2: memcpy(jin, z, 32); memcpy(jin + 32, cb[j], L); shake256(kb, 32, jin, 32 + L); break;
                case 3: pke_encrypt(P, ek, mb[j], rb[j], cp); break;
                case 4: { uint8_t diff = 0; for (size_t i = 0; i < L; i++) diff |= (uint8_t)(cb[j][i] ^ f2[i]);
                          uint8_t mask = (uint8_t)(-(int)((-(uint32_t)diff) >> 31)); __asm__ volatile("" : "+r"(mask));
                          for (int i = 0; i < 32; i++) K[i] = (uint8_t)(g[i] ^ (mask & (g[i] ^ kb[i]))); } break;
                case 5: { poly y; for (int i = 0; i < k; i++) cbd(&y, rb[j], (uint8_t)i, P->eta1); __asm__ volatile("" :: "r"(&y) : "memory"); } break;
                }
                uint64_t e = cyc(); t[b0 + j] = e - a;
            }
        }
        const char *nm[] = {"S1 pke_decrypt (ct)", "S2 G = SHA3-512(m||h)", "S3 J = SHAKE256(z||c)",
                            "S4 pke_encrypt (m, r)  [re-encryption]", "S5 compare + select", "S6 cbd(eta1) x k"};
        report(nm[st], t, cl, n);
    }
    return 0;
}
