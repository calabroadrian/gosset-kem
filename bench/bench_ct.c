/*
 * Gosset-KEM v4 — benchmark + dudect-style constant-time test
 * Timing test (Reparaz, Balasch, Verbauwhede 2017 methodology, own code):
 *   class 0: fixed input  = an exact D8 codeword, no noise
 *   class 1: random input = uniform coefficients mod q
 * Classes interleaved at random; Welch t-test on cycle counts, also on
 * measurements cropped at several percentiles. |t| > 4.5 => timing leak.
 * Positive control: e8_decode_leaky must FAIL this test.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>
#include "e8codec.h"

static inline uint64_t cycles(void) {
    unsigned aux; _mm_lfence(); uint64_t t = __rdtscp(&aux); _mm_lfence(); return t;
}

static uint64_t rng_state = 88172645463325252ULL;
static inline uint64_t xr(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y);
}

typedef void (*dec_fn)(const int16_t *, uint8_t *);

static double median_cycles(dec_fn f, int reps) {
    int16_t in[256]; uint8_t out[32];
    uint64_t *t = malloc(sizeof(uint64_t) * reps);
    for (int r = 0; r < reps; r++) {
        for (int i = 0; i < 256; i++) in[i] = (int16_t)(xr() % 3329);
        uint64_t a = cycles(); f(in, out); uint64_t b = cycles();
        t[r] = b - a;
    }
    qsort(t, reps, sizeof(uint64_t), cmp_u64);
    double m = (double)t[reps / 2]; free(t); return m;
}

static double enc_median(int reps) {
    uint8_t msg[32]; int16_t out[256];
    uint64_t *t = malloc(sizeof(uint64_t) * reps);
    for (int r = 0; r < reps; r++) {
        for (int i = 0; i < 32; i++) msg[i] = (uint8_t)xr();
        uint64_t a = cycles(); e8_encode(msg, out); uint64_t b = cycles();
        t[r] = b - a;
    }
    qsort(t, reps, sizeof(uint64_t), cmp_u64);
    double m = (double)t[reps / 2]; free(t); return m;
}

/* Welch t on samples x (class0) y (class1) keeping only values <= thr */
static double welch(const uint64_t *v, const int *cls, int n, double thr) {
    double s[2] = {0, 0}, s2[2] = {0, 0}; double c[2] = {0, 0};
    for (int i = 0; i < n; i++) {
        double x = (double)v[i]; if (x > thr) continue;
        s[cls[i]] += x; s2[cls[i]] += x * x; c[cls[i]] += 1;
    }
    double m0 = s[0] / c[0], m1 = s[1] / c[1];
    double v0 = s2[0] / c[0] - m0 * m0, v1 = s2[1] / c[1] - m1 * m1;
    return (m0 - m1) / sqrt(v0 / c[0] + v1 / c[1]);
}

static double dudect(dec_fn f, int n, const char *name) {
    uint64_t *t = malloc(sizeof(uint64_t) * n);
    int *cls = malloc(sizeof(int) * n);
    int16_t fixed[256]; uint8_t out[32], msg[32];
    for (int i = 0; i < 32; i++) msg[i] = (uint8_t)(xr() & 0x7f);   /* D8 coset */
    e8_encode(msg, fixed);
    /* dudect practice: prepare a whole batch of inputs BEFORE measuring, so
       both classes see identical memory/cache handling */
    const int BATCH = 4096;
    int16_t (*inb)[256] = malloc(sizeof(int16_t) * 256 * BATCH);
    for (int base = 0; base < n; base += BATCH) {
        int m = (n - base < BATCH) ? n - base : BATCH;
        for (int j = 0; j < m; j++) {
            cls[base + j] = (int)(xr() & 1);
            if (cls[base + j] == 0) memcpy(inb[j], fixed, sizeof fixed);
            else for (int i = 0; i < 256; i++) inb[j][i] = (int16_t)(xr() % 3329);
        }
        for (int j = 0; j < m; j++) {
            uint64_t a = cycles(); f(inb[j], out); uint64_t b = cycles();
            t[base + j] = b - a;
        }
    }
    free(inb);
    uint64_t *sorted = malloc(sizeof(uint64_t) * n);
    memcpy(sorted, t, sizeof(uint64_t) * n);
    qsort(sorted, n, sizeof(uint64_t), cmp_u64);
    double pct[] = {0.10, 0.25, 0.50, 0.75, 0.90, 0.99, 1.0};
    double worst = 0;
    printf("  %-16s", name);
    for (int k = 0; k < 7; k++) {
        double thr = pct[k] >= 1.0 ? 1e300 : (double)sorted[(int)(pct[k] * n)];
        double tv = welch(t, cls, n, thr);
        printf(" %7.2f", tv);
        if (fabs(tv) > fabs(worst)) worst = tv;
    }
    printf("   max|t|=%6.2f  %s\n", fabs(worst), fabs(worst) > 4.5 ? "LEAK" : "no leak detected");
    free(t); free(cls); free(sorted);
    return fabs(worst);
}

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 2000000;
    /* warm up */
    median_cycles(e8_decode, 20000);
    printf("BENCHMARK (median cycles per 256-coefficient message, %s)\n", "rdtscp");
    double sd = median_cycles(std_decode, 200000);
    double ed = median_cycles(e8_decode, 200000);
    double ee = enc_median(200000);
    printf("  std_decode (ML-KEM poly_tomsg) : %8.0f\n", sd);
    printf("  e8_decode  (constant-time)     : %8.0f   (%.1fx)\n", ed, ed / sd);
    printf("  e8_encode                      : %8.0f\n", ee);
    printf("\nCONSTANT-TIME TEST (Welch t; columns = crop at p10 p25 p50 p75 p90 p99 all), n=%d\n", n);
    dudect(e8_decode, n, "e8_decode");
    dudect(std_decode, n, "std_decode");
    dudect(e8_decode_leaky, n, "e8_decode_LEAKY");
    return 0;
}
