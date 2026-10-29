/* AFA-12 v9 — hybrid KEM cycles and constant-time decapsulation */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>
#include "hybrid.h"
#include "x25519.h"
static uint64_t st = 0x9E3779B97F4A7C15ULL;
static inline uint64_t xr(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static void rnd(uint8_t *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = (uint8_t)xr(); }
static inline uint64_t cyc(void) { unsigned a; _mm_lfence(); uint64_t t = __rdtscp(&a); _mm_lfence(); return t; }
static int cmpu(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y); }
static double welch(const uint64_t *v, const int *c, int n, double thr) {
    double s[2] = {0}, q[2] = {0}, k[2] = {0};
    for (int i = 0; i < n; i++) { double x = (double)v[i]; if (x > thr) continue; s[c[i]] += x; q[c[i]] += x * x; k[c[i]]++; }
    double m0 = s[0] / k[0], m1 = s[1] / k[1];
    return (m0 - m1) / sqrt((q[0] / k[0] - m0 * m0) / k[0] + (q[1] / k[1] - m1 * m1) / k[1]);
}
int main(int argc, char **argv) {
    int reps = argc > 1 ? atoi(argv[1]) : 2000, n = argc > 2 ? atoi(argv[2]) : 50000;
    uint8_t pk[1300], sk[2600], ct[1200], K[32], K2[32], d[32], z[32], x[32], m[32], e[32];
    { uint64_t *t = malloc(8 * (size_t)reps); uint8_t o[32];
      for (int r = 0; r < reps; r++) { rnd(x, 32); uint64_t a = cyc(); x25519_base(o, x); t[r] = cyc() - a; }
      qsort(t, reps, 8, cmpu); printf("X25519 scalar multiplication: %llu cycles (median)\n\n", (unsigned long long)t[reps / 2]); free(t); }
    printf("%-34s %8s %8s %8s %10s %10s %10s\n", "mode", "pk B", "sk B", "ct B", "keygen", "encaps", "decaps");
    int modes[2] = {1, 7}; const char *names[2] = {"C: ML-KEM-768 + X25519", "S: AFA12-768-MaxSec + X25519"};
    uint64_t *T[2][3]; for (int mi = 0; mi < 2; mi++) for (int q = 0; q < 3; q++) T[mi][q] = malloc(8 * (size_t)reps);
    for (int r = 0; r < reps; r++)
        for (int v = 0; v < 2; v++) {
            int mi = (r + v) & 1, id = modes[mi];                 /* interleaved, alternating order */
            rnd(d, 32); rnd(z, 32); rnd(x, 32); rnd(m, 32); rnd(e, 32);
            uint64_t a = cyc(); hyb_keygen_internal(id, d, z, x, pk, sk); uint64_t b = cyc();
            hyb_encaps_internal(id, pk, m, e, ct, K); uint64_t c = cyc();
            hyb_decaps(id, sk, ct, K2); uint64_t f = cyc();
            T[mi][0][r] = b - a; T[mi][1][r] = c - b; T[mi][2][r] = f - c;
        }
    for (int mi = 0; mi < 2; mi++) {
        for (int q = 0; q < 3; q++) qsort(T[mi][q], reps, 8, cmpu);
        printf("%-34s %8zu %8zu %8zu %10llu %10llu %10llu\n", names[mi], hyb_pk_bytes(modes[mi]), hyb_sk_bytes(modes[mi]), hyb_ct_bytes(modes[mi]),
               (unsigned long long)T[mi][0][reps / 2], (unsigned long long)T[mi][1][reps / 2], (unsigned long long)T[mi][2][reps / 2]);
    }
    printf("\nCONSTANT-TIME HYBRID DECAPS (Welch t, fixed valid ct vs random ct, n=%d)\n", n);
    for (int mi = 0; mi < 2; mi++) {
        int id = modes[mi]; size_t L = hyb_ct_bytes(id);
        rnd(d, 32); rnd(z, 32); rnd(x, 32); rnd(m, 32); rnd(e, 32);
        hyb_keygen_internal(id, d, z, x, pk, sk); uint8_t fixed[1200]; hyb_encaps_internal(id, pk, m, e, fixed, K);
        uint64_t *t = malloc(8 * (size_t)n), *so = malloc(8 * (size_t)n); int *cl = malloc(4 * (size_t)n);
        const int B = 256; uint8_t (*cb)[1200] = malloc(1200 * B);
        for (int base = 0; base < n; base += B) {
            int mm = n - base < B ? n - base : B;
            for (int j = 0; j < mm; j++) { cl[base + j] = (int)(xr() & 1); if (!cl[base + j]) memcpy(cb[j], fixed, L); else rnd(cb[j], L); }
            for (int j = 0; j < mm; j++) { uint64_t a = cyc(); hyb_decaps(id, sk, cb[j], K2); t[base + j] = cyc() - a; }
        }
        memcpy(so, t, 8 * (size_t)n); qsort(so, n, 8, cmpu);
        double p[] = {0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0}, w = 0;
        for (int k = 0; k < 7; k++) { double tv = welch(t, cl, n, p[k] >= 1 ? 1e300 : (double)so[(int)(p[k] * n)]); if (fabs(tv) > fabs(w)) w = tv; }
        printf("  %-32s max|t| = %5.2f  %s\n", names[mi], fabs(w), fabs(w) > 4.5 ? "LEAK" : "no leak detected");
        free(t); free(so); free(cl); free(cb);
    }
    return 0;
}
