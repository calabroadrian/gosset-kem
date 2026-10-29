/* AFA-12 v5 — Leech decoder benchmark + dudect-style constant-time test */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>
#include "leech.h"

void e8_decode_fast(const int16_t in[256], uint8_t msg[32]);
void std_decode_ref(const int16_t in[256], uint8_t msg[32]);

static uint64_t st = 0x2545F4914F6CDD1DULL;
static inline uint64_t xr(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static inline uint64_t cyc(void) { unsigned a; _mm_lfence(); uint64_t t = __rdtscp(&a); _mm_lfence(); return t; }
static int cmpu(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y); }

static void leech_msg(const int16_t *in, uint8_t *out) { leech_decode(in, (uint64_t *)out); }
typedef void (*fn)(const int16_t *, uint8_t *);

static double med(fn f, int reps) {
    int16_t in[256]; uint8_t out[64];
    uint64_t *t = malloc(8 * (size_t)reps);
    for (int r = 0; r < reps; r++) {
        for (int i = 0; i < 256; i++) in[i] = (int16_t)(xr() % 3329);
        uint64_t a = cyc(); f(in, out); uint64_t b = cyc(); t[r] = b - a;
    }
    qsort(t, reps, 8, cmpu); double m = (double)t[reps / 2]; free(t); return m;
}

static double welch(const uint64_t *v, const int *c, int n, double thr) {
    double s[2] = {0}, q[2] = {0}, k[2] = {0};
    for (int i = 0; i < n; i++) { double x = (double)v[i]; if (x > thr) continue; s[c[i]] += x; q[c[i]] += x * x; k[c[i]]++; }
    double m0 = s[0] / k[0], m1 = s[1] / k[1];
    return (m0 - m1) / sqrt((q[0] / k[0] - m0 * m0) / k[0] + (q[1] / k[1] - m1 * m1) / k[1]);
}

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 300000;
    med(leech_msg, 200);
    double s = med(std_decode_ref, 200000), e = med(e8_decode_fast, 200000), l = med(leech_msg, 5000);
    printf("BENCHMARK (median cycles per decapsulation-message decode)\n");
    printf("  ML-KEM poly_tomsg (256 coefs)      : %9.0f\n", s);
    printf("  E8 constant-time (32 blocks)       : %9.0f  (%.1fx)\n", e, e / s);
    printf("  Leech constant-time ML (8 blocks)  : %9.0f  (%.0fx)\n", l, l / s);

    /* dudect: fixed (encoded message, no noise) vs random */
    int16_t fixed[256] = {0};
    for (int b = 0; b < 8; b++) leech_encode_block(xr() & ((1ULL << 36) - 1), fixed + 24 * b);
    uint64_t *t = malloc(8 * (size_t)n), *so = malloc(8 * (size_t)n); int *c = malloc(4 * (size_t)n);
    const int B = 1024; int16_t (*inb)[256] = malloc(sizeof(int16_t) * 256 * B); uint8_t out[64];
    for (int base = 0; base < n; base += B) {
        int m = n - base < B ? n - base : B;
        for (int j = 0; j < m; j++) { c[base + j] = (int)(xr() & 1);
            if (!c[base + j]) memcpy(inb[j], fixed, sizeof fixed); else for (int i = 0; i < 256; i++) inb[j][i] = (int16_t)(xr() % 3329); }
        for (int j = 0; j < m; j++) { uint64_t a = cyc(); leech_msg(inb[j], out); uint64_t b = cyc(); t[base + j] = b - a; }
    }
    memcpy(so, t, 8 * (size_t)n); qsort(so, n, 8, cmpu);
    double p[] = {0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0}, w = 0;
    for (int k = 0; k < 7; k++) { double tv = welch(t, c, n, p[k] >= 1 ? 1e300 : (double)so[(int)(p[k] * n)]); if (fabs(tv) > fabs(w)) w = tv; }
    printf("\nCONSTANT-TIME TEST (Welch t, n=%d): leech_decode max|t| = %.2f  %s\n", n, fabs(w), fabs(w) > 4.5 ? "LEAK" : "no leak detected");
    return 0;
}
