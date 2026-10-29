/* AFA-12 v4 — correctness + benchmark + constant-time test for optimized decoders */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>
#include "e8codec.h"

void e8_decode_fast(const int16_t in[256], uint8_t msg[32]);
void e8_decode_avx2(const int16_t in[256], uint8_t msg[32]);
void std_decode_ref(const int16_t in[256], uint8_t msg[32]);
typedef void (*dec_fn)(const int16_t *, uint8_t *);

static uint64_t st = 0x9E3779B97F4A7C15ULL;
static inline uint64_t xr(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static inline uint64_t cyc(void) { unsigned a; _mm_lfence(); uint64_t t = __rdtscp(&a); _mm_lfence(); return t; }
static int cmpu(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y); }

/* squared distance (scaled by 16) from y to nearest point of coset given by message index */
static long long coset_dist(const int16_t *y, uint8_t m) {
    double best = 0; int par = 0; double rep[8];
    for (int i = 0; i < 8; i++) { int bit = i < 7 ? (m >> i) & 1 : par; par ^= bit; rep[i] = bit + ((m >> 7) ? 0.5 : 0.0); }
    for (int i = 0; i < 8; i++) {
        double yc = y[i] > 3329 / 2 ? y[i] - 3329 : y[i];
        double z = yc / 1664.5 - rep[i];
        z -= 2 * floor(z / 2 + 0.5);
        best += z * z;
    }
    return llround(best * 1e9);
}

static double med(dec_fn f, int reps) {
    int16_t (*in)[256] = malloc(sizeof(int16_t) * 256 * 1024); uint8_t out[32];
    uint64_t *t = malloc(8 * (size_t)reps);
    for (int r = 0; r < reps; r++) {
        if (r % 1024 == 0) for (int j = 0; j < 1024; j++) for (int i = 0; i < 256; i++) in[j][i] = (int16_t)(xr() % 3329);
        uint64_t a = cyc(); f(in[r % 1024], out); uint64_t b = cyc(); t[r] = b - a;
    }
    qsort(t, reps, 8, cmpu); double m = (double)t[reps / 2]; free(t); free(in); return m;
}

static double welch(const uint64_t *v, const int *c, int n, double thr) {
    double s[2] = {0}, q[2] = {0}, k[2] = {0};
    for (int i = 0; i < n; i++) { double x = (double)v[i]; if (x > thr) continue; s[c[i]] += x; q[c[i]] += x * x; k[c[i]]++; }
    double m0 = s[0] / k[0], m1 = s[1] / k[1];
    return (m0 - m1) / sqrt((q[0] / k[0] - m0 * m0) / k[0] + (q[1] / k[1] - m1 * m1) / k[1]);
}

static void dudect(dec_fn f, int n, const char *name) {
    uint64_t *t = malloc(8 * (size_t)n), *s = malloc(8 * (size_t)n); int *c = malloc(4 * (size_t)n);
    int16_t fixed[256]; uint8_t msg[32], out[32];
    for (int i = 0; i < 32; i++) msg[i] = (uint8_t)(xr() & 0x7f);
    e8_encode(msg, fixed);
    const int B = 4096; int16_t (*inb)[256] = malloc(sizeof(int16_t) * 256 * B);
    for (int base = 0; base < n; base += B) {
        int m = n - base < B ? n - base : B;
        for (int j = 0; j < m; j++) { c[base + j] = (int)(xr() & 1);
            if (!c[base + j]) memcpy(inb[j], fixed, sizeof fixed); else for (int i = 0; i < 256; i++) inb[j][i] = (int16_t)(xr() % 3329); }
        for (int j = 0; j < m; j++) { uint64_t a = cyc(); f(inb[j], out); uint64_t b = cyc(); t[base + j] = b - a; }
    }
    memcpy(s, t, 8 * (size_t)n); qsort(s, n, 8, cmpu);
    double p[] = {0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0}, w = 0;
    for (int k = 0; k < 7; k++) { double tv = welch(t, c, n, p[k] >= 1 ? 1e300 : (double)s[(int)(p[k] * n)]); if (fabs(tv) > fabs(w)) w = tv; }
    printf("  %-16s max|t| = %6.2f  %s\n", name, fabs(w), fabs(w) > 4.5 ? "LEAK" : "no leak detected");
    free(t); free(s); free(c); free(inb);
}

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 5000000;
    /* correctness vs the verified v4 decoder: identical output, or an exact tie */
    long blocks = 0, diff_fast = 0, diff_avx = 0, nontie = 0;
    for (int r = 0; r < 200000; r++) {
        int16_t in[256]; uint8_t a[32], b[32], c[32];
        if (r & 1) for (int i = 0; i < 256; i++) in[i] = (int16_t)(xr() % 3329);
        else { uint8_t m[32]; for (int i = 0; i < 32; i++) m[i] = (uint8_t)xr(); e8_encode(m, in);
               for (int i = 0; i < 256; i++) { int v = in[i] + (int)(xr() % 1201) - 600; in[i] = (int16_t)(((v % 3329) + 3329) % 3329); } }
        e8_decode(in, a); e8_decode_fast(in, b); e8_decode_avx2(in, c);
        for (int k = 0; k < 32; k++) {
            blocks++;
            if (a[k] != b[k]) { diff_fast++; if (coset_dist(in + 8 * k, a[k]) != coset_dist(in + 8 * k, b[k])) nontie++; }
            if (a[k] != c[k]) { diff_avx++;  if (coset_dist(in + 8 * k, a[k]) != coset_dist(in + 8 * k, c[k])) nontie++; }
        }
    }
    printf("CORRECTNESS vs verified decoder: %ld blocks, fast diffs %ld, avx2 diffs %ld, non-tie diffs %ld\n",
           blocks, diff_fast, diff_avx, nontie);
    med(e8_decode_avx2, 20000);
    double s = med(std_decode_ref, 300000), e0 = med(e8_decode, 300000), e1 = med(e8_decode_fast, 300000), e2 = med(e8_decode_avx2, 300000);
    printf("\nBENCHMARK (median cycles / 256 coefficients)\n");
    printf("  ML-KEM poly_tomsg (reference)  : %7.0f\n", s);
    printf("  E8 decoder v4.0 (first version): %7.0f  (%.1fx)\n", e0, e0 / s);
    printf("  E8 decoder fast (scalar int32) : %7.0f  (%.1fx)\n", e1, e1 / s);
    printf("  E8 decoder AVX2                : %7.0f  (%.1fx)\n", e2, e2 / s);
    printf("\nCONSTANT-TIME TEST, n=%d\n", n);
    dudect(e8_decode_fast, n, "e8_decode_fast");
    dudect(e8_decode_avx2, n, "e8_decode_avx2");
    dudect(std_decode_ref, n, "std_decode_ref");
    dudect(e8_decode_leaky, n, "LEAKY control");
    return 0;
}
