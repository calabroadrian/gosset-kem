/* build: gcc -O3 -mavx2 -o bench_mdecaps bench_mdecaps.c masked_decaps.c masked_core.c masked_keccak.c masked_e8bs.c ../v7/fips202.c ../v4/e8fast.c ../v4/e8codec.c -lm */
/* masked vs unmasked decapsulation: cycles (interleaved) and first/second-order TVLA (simulated HW leakage) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>
#include "../v7/mlkem.h"
int mlkem_decaps_masked(int id, const uint8_t *dk, const uint8_t *c, uint8_t K[32]);
void mk_seed(uint64_t s); void mk_tvla_init(int max_points, double noise); void mk_trace_begin(int c); int mk_trace_end(void);
double mk_tvla_first(int n); double mk_tvla_second(int n); void mk_tvla_diag(int n, double thr, int *count, int *argmax, double *maxv);
static uint64_t st = 0x51ED2701ULL ^ 0xA5A5F00DULL;
static inline uint64_t xr(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static void rnd(uint8_t *b, size_t n) { for (size_t i = 0; i < n; i++) b[i] = (uint8_t)xr(); }
static inline uint64_t cyc(void){unsigned a;_mm_lfence();uint64_t t=__rdtscp(&a);_mm_lfence();return t;}
static int cmpu(const void*a,const void*b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return (x>y)-(x<y);}
int main(int argc, char **argv) {
    uint64_t seedarg = argc > 4 ? strtoull(argv[4], 0, 0) : 0; int reps = argc > 1 ? atoi(argv[1]) : 300, ntr = argc > 2 ? atoi(argv[2]) : 5000, control = argc > 3 ? atoi(argv[3]) : 0;
    uint8_t d[32], z[32], m[32], ek[1600], dk[3200], c[1600], K[32];
    if (seedarg) st ^= seedarg;
    printf("CYCLES (median, masked and unmasked interleaved, %d reps)\n", reps);
    for (int id = 6; id <= 8 && reps > 0; id++) {
        const mlkem_params *P = mlkem_get(id);
        uint64_t *tu = malloc(8 * reps), *tm = malloc(8 * reps);
        for (int r = 0; r < reps; r++) {
            rnd(d, 32); rnd(z, 32); rnd(m, 32); mk_seed(xr());
            mlkem_keygen_internal(id, d, z, ek, dk); mlkem_encaps_internal(id, ek, m, K, c);
            uint64_t a = cyc(); mlkem_decaps(id, dk, c, K); uint64_t b = cyc();
            mlkem_decaps_masked(id, dk, c, K); uint64_t e = cyc();
            tu[r] = b - a; tm[r] = e - b;
        }
        qsort(tu, reps, 8, cmpu); qsort(tm, reps, 8, cmpu);
        printf("  %-20s unmasked %9llu   masked %10llu   (x%.1f)\n", P->name, (unsigned long long)tu[reps/2], (unsigned long long)tm[reps/2], (double)tm[reps/2]/tu[reps/2]);
        free(tu); free(tm);
    }
    /* TVLA: AFA12-768-MaxSec, fixed valid ciphertext vs random ciphertexts, same key, fresh masks */
    int id = 7; const mlkem_params *P = mlkem_get(id); size_t L = mlkem_ct_bytes(P);
    rnd(d, 32); rnd(z, 32); rnd(m, 32);
    mlkem_keygen_internal(id, d, z, ek, dk); uint8_t fixed[1600]; mlkem_encaps_internal(id, ek, m, K, fixed);
    mk_tvla_init(400000, 1.0);
    int npts = 0;
    for (int t = 0; t < 2 * ntr; t++) {
        int cl = t & 1;
        if (cl == 0 && !control) memcpy(c, fixed, L); else rnd(c, L);
        mk_seed(xr());
        mk_trace_begin(cl); mlkem_decaps_masked(id, dk, c, K); npts = mk_trace_end();
    }
    printf("\nTVLA on full masked decapsulation (%s), %s, sigma=1, %d+%d traces, %d leakage points per trace\n", P->name, control ? "CONTROL random-vs-random" : "fixed-vs-random", ntr, ntr, npts);
    printf("  first order : max|t| = %.2f  %s\n", mk_tvla_first(npts), mk_tvla_first(npts) > 4.5 ? "LEAK" : "no first-order leak detected");
    { int c45, c40, am; double mv; mk_tvla_diag(npts, 4.5, &c45, &am, &mv); mk_tvla_diag(npts, 4.0, &c40, &am, &mv);
      printf("  first order diagnostics: points |t|>4.5: %d (null expectation %.2f), |t|>4.0: %d (null expectation %.2f), argmax point %d\n",
             c45, npts * 6.8e-6, c40, npts * 6.33e-5, am); }
    printf("  second order: max|t| = %.2f  %s\n", mk_tvla_second(npts), mk_tvla_second(npts) > 4.5 ? "leak (expected for 1st-order masking)" : "none detected");
    return 0;
}
