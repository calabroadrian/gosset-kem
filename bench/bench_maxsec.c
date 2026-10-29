/* Gosset-KEM v7 — ML-KEM vs ML-KEM+E8: sizes, round trips, implicit rejection, cycles, constant time */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>
#include "mlkem.h"

static uint64_t st = 0x853C49E6748FEA9BULL;
static inline uint64_t xr(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static void rnd(uint8_t *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = (uint8_t)xr(); }
static inline uint64_t cyc(void) { unsigned a; _mm_lfence(); uint64_t t = __rdtscp(&a); _mm_lfence(); return t; }
static int cmpu(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y); }

typedef struct { uint8_t ek[1600], dk[3200], c[1600], K[32], K2[32]; } bufs;

static double welch(const uint64_t *v, const int *c, int n, double thr) {
    double s[2] = {0}, q[2] = {0}, k[2] = {0};
    for (int i = 0; i < n; i++) { double x = (double)v[i]; if (x > thr) continue; s[c[i]] += x; q[c[i]] += x * x; k[c[i]]++; }
    double m0 = s[0] / k[0], m1 = s[1] / k[1];
    return (m0 - m1) / sqrt((q[0] / k[0] - m0 * m0) / k[0] + (q[1] / k[1] - m1 * m1) / k[1]);
}

int main(int argc, char **argv) {
    int trips = argc > 1 ? atoi(argv[1]) : 20000, reps = argc > 2 ? atoi(argv[2]) : 2000, ctn = argc > 3 ? atoi(argv[3]) : 100000;
    bufs B; uint8_t d[32], z[32], m[32];
    printf("%-16s %6s %6s %6s %9s %s\n", "scheme", "ek B", "dk B", "ct B", "roundtrip", "implicit-reject");
    for (int id = 0; id < 9; id++) {
        if (id >= 3 && id < 6) continue;   /* compact mode measured in bench_kem */
        const mlkem_params *P = mlkem_get(id);
        /* round trips */
        int bad = 0;
        for (int t = 0; t < trips; t++) {
            rnd(d, 32); rnd(z, 32); rnd(m, 32);
            if (t % 50 == 0) mlkem_keygen_internal(id, d, z, B.ek, B.dk);
            mlkem_encaps_internal(id, B.ek, m, B.K, B.c);
            mlkem_decaps(id, B.dk, B.c, B.K2);
            bad += memcmp(B.K, B.K2, 32) != 0;
        }
        /* implicit rejection: tampered ciphertext must give K != real K */
        int rej_ok = 0;
        for (int t = 0; t < 200; t++) {
            rnd(m, 32); mlkem_encaps_internal(id, B.ek, m, B.K, B.c);
            B.c[xr() % mlkem_ct_bytes(P)] ^= (uint8_t)(1u << (xr() % 8));
            mlkem_decaps(id, B.dk, B.c, B.K2);
            rej_ok += memcmp(B.K, B.K2, 32) != 0;
        }
        /* cycles: measured later, interleaved with the standard variant */
        printf("%-16s %6zu %6zu %6zu %9s %d/200\n", P->name, mlkem_ek_bytes(P), mlkem_dk_bytes(P),
               mlkem_ct_bytes(P), bad ? "FAIL" : "OK", rej_ok);
        fflush(stdout);
    }
    printf("\nCYCLES (median, standard and E8 interleaved in the same loop, %d reps each)\n", reps);
    printf("%-14s %12s %12s %12s %12s %12s %12s\n", "level", "keygen std", "keygen Max", "encaps std", "encaps Max", "decaps std", "decaps Max");
    for (int lv = 0; lv < 3; lv++) {
        uint64_t *T[6]; for (int q = 0; q < 6; q++) T[q] = malloc(8 * (size_t)reps);
        for (int r = 0; r < reps; r++)
            for (int v = 0; v < 2; v++) {
                int id = lv + 6 * ((r + v) & 1);          /* alternate order every rep */
                int slot = (id >= 6);
                rnd(d, 32); rnd(z, 32); rnd(m, 32);
                uint64_t a = cyc(); mlkem_keygen_internal(id, d, z, B.ek, B.dk); uint64_t b = cyc();
                mlkem_encaps_internal(id, B.ek, m, B.K, B.c); uint64_t c = cyc();
                mlkem_decaps(id, B.dk, B.c, B.K2); uint64_t e = cyc();
                T[0 + slot][r] = b - a; T[2 + slot][r] = c - b; T[4 + slot][r] = e - c;
            }
        double md[6];
        for (int q = 0; q < 6; q++) { qsort(T[q], reps, 8, cmpu); md[q] = (double)T[q][reps / 2]; free(T[q]); }
        printf("%-14s %12.0f %12.0f %12.0f %12.0f %12.0f %12.0f\n", mlkem_get(lv)->name, md[0], md[1], md[2], md[3], md[4], md[5]);
        printf("%-14s %12s %+11.1f%% %12s %+11.1f%% %12s %+11.1f%%   (decaps delta %+.0f cycles)\n", "  MaxSec vs std", "",
               100 * (md[1] - md[0]) / md[0], "", 100 * (md[3] - md[2]) / md[2], "", 100 * (md[5] - md[4]) / md[4], md[5] - md[4]);
        fflush(stdout);
    }
    /* constant-time decaps: class 0 valid ciphertext (fixed), class 1 random ciphertext */
    printf("\nCONSTANT-TIME DECAPS (Welch t, fixed valid ct vs random ct, n=%d)\n", ctn);
    for (int id = 0; id < 9; id++) {
        if (id >= 3 && id < 6) continue;
        const mlkem_params *P = mlkem_get(id); size_t L = mlkem_ct_bytes(P);
        rnd(d, 32); rnd(z, 32); rnd(m, 32);
        mlkem_keygen_internal(id, d, z, B.ek, B.dk);
        uint8_t fixed[1600]; mlkem_encaps_internal(id, B.ek, m, B.K, fixed);
        uint64_t *t = malloc(8 * (size_t)ctn), *so = malloc(8 * (size_t)ctn); int *cl = malloc(4 * (size_t)ctn);
        const int BB = 256; uint8_t (*cb)[1600] = malloc(1600 * BB);
        for (int base0 = 0; base0 < ctn; base0 += BB) {
            int mm = ctn - base0 < BB ? ctn - base0 : BB;
            for (int j = 0; j < mm; j++) { cl[base0 + j] = (int)(xr() & 1); if (!cl[base0 + j]) memcpy(cb[j], fixed, L); else rnd(cb[j], L); }
            for (int j = 0; j < mm; j++) { uint64_t a = cyc(); mlkem_decaps(id, B.dk, cb[j], B.K2); uint64_t b = cyc(); t[base0 + j] = b - a; }
        }
        memcpy(so, t, 8 * (size_t)ctn); qsort(so, ctn, 8, cmpu);
        double p[] = {0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0}, w = 0;
        for (int k = 0; k < 7; k++) { double tv = welch(t, cl, ctn, p[k] >= 1 ? 1e300 : (double)so[(int)(p[k] * ctn)]); if (fabs(tv) > fabs(w)) w = tv; }
        printf("  %-16s max|t| = %6.2f  %s\n", P->name, fabs(w), fabs(w) > 4.5 ? "LEAK" : "no leak detected");
        free(t); free(so); free(cl); free(cb); fflush(stdout);
    }
    return 0;
}
