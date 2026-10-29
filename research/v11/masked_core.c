/* AFA-12 v11 — masking core: PRNG, leakage accumulators, bitsliced masked arithmetic */
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "masked_core.h"

uint64_t mk_rs0 = 0x9E3779B97F4A7C15ULL, mk_rs1 = 0xD1B54A32D192ED03ULL;
void mk_seed(uint64_t s) { mk_rs0 = s ^ 0x9E3779B97F4A7C15ULL; mk_rs1 = (s * 0xBF58476D1CE4E5B9ULL) | 1; for (int i = 0; i < 8; i++) mk_rnd(); }

/* ---- leakage: online Welch accumulators per trace point (two classes) ---- */
int mk_leak_on = 0;
static double *A_sum[2], *A_sq[2], *P_sum[2], *P_sq[2];   /* first order and second-order (product of share pair) */
static long A_n[2]; static int cls = 0, pt = 0, maxpts = 0; static double sigma = 1.0;
static double prev_hw; static int have_prev = 0;
static uint64_t gs = 0x1234567ULL;
static double gauss(void) {                                   /* Box-Muller on a separate PRNG */
    gs ^= gs << 13; gs ^= gs >> 7; gs ^= gs << 17; double u1 = ((gs >> 11) + 1.0) / 9007199254740993.0;
    gs ^= gs << 13; gs ^= gs >> 7; gs ^= gs << 17; double u2 = (gs >> 11) / 9007199254740992.0;
    return sqrt(-2 * log(u1)) * cos(6.283185307179586 * u2);
}
void mk_tvla_init(int max_points, double noise) {
    maxpts = max_points; sigma = noise;
    for (int c = 0; c < 2; c++) {
        A_sum[c] = calloc(max_points, sizeof(double)); A_sq[c] = calloc(max_points, sizeof(double));
        P_sum[c] = calloc(max_points, sizeof(double)); P_sq[c] = calloc(max_points, sizeof(double)); A_n[c] = 0;
    }
}
void mk_trace_begin(int c) { cls = c; pt = 0; have_prev = 0; mk_leak_on = 1; }
int mk_trace_end(void) { mk_leak_on = 0; A_n[cls]++; return pt; }
void mk_leak_word(uint64_t v) {
    if (pt >= maxpts) return;
    double h = (double)__builtin_popcountll(v) + sigma * gauss();
    A_sum[cls][pt] += h; A_sq[cls][pt] += h * h;
    if (have_prev) {                                          /* pairs (s0, s1) are leaked consecutively */
        double p = (prev_hw - 32.0) * (h - 32.0);
        P_sum[cls][pt] += p; P_sq[cls][pt] += p * p; have_prev = 0;
    } else { prev_hw = h; have_prev = 1; }
    pt++;
}
static double tmax(double **S, double **SQ, int n) {
    double m = 0;
    for (int i = 0; i < n; i++) {
        double m0 = S[0][i] / A_n[0], m1 = S[1][i] / A_n[1];
        double v0 = SQ[0][i] / A_n[0] - m0 * m0, v1 = SQ[1][i] / A_n[1] - m1 * m1;
        double den = sqrt(v0 / A_n[0] + v1 / A_n[1]);
        if (den > 0) { double t = fabs((m0 - m1) / den); if (t > m) m = t; }
    }
    return m;
}
double mk_tvla_first(int npts) { return tmax(A_sum, A_sq, npts); }
double mk_tvla_second(int npts) { return tmax(P_sum, P_sq, npts); }

/* ---- bitsliced masked arithmetic ---- */
void bs_add(mw *z, const mw *x, const mw *y, int n) {
    mw c = {0, 0};
    for (int j = 0; j < n; j++) {
        mw t = mx(x[j], y[j]);
        mw s = mx(t, c);
        if (j < n - 1) { mw g = secand(x[j], y[j]); mw p = secand(c, t); c = mx(g, p); }
        z[j] = s;
    }
}
void bs_addpub(mw *z, const mw *x, const uint64_t *k, int n) {
    mw c = {0, 0};
    for (int j = 0; j < n; j++) {
        mw t = mx(x[j], mw_pub(k[j]));
        mw s = mx(t, c);
        if (j < n - 1) { mw g = mandp(x[j], k[j]); mw p = secand(c, t); c = mx(g, p); }
        z[j] = s;
    }
}
void bs_addconst(mw *z, const mw *x, uint32_t c, int n) {
    uint64_t k[32]; for (int j = 0; j < n; j++) k[j] = ((c >> j) & 1) ? ~0ULL : 0; bs_addpub(z, x, k, n);
}
void bs_neg(mw *z, const mw *x, int n) {
    mw t[32]; for (int j = 0; j < n; j++) t[j] = mnot(x[j]); bs_addconst(z, t, 1, n);
}
void bs_sel(mw *z, mw m, const mw *a, const mw *b, int n) { for (int j = 0; j < n; j++) z[j] = secsel(m, a[j], b[j]); }

void bs_transpose_in(uint64_t *bits, const uint32_t *v, int lanes, int n) {
    for (int j = 0; j < n; j++) { uint64_t w = 0; for (int l = 0; l < lanes; l++) w |= (uint64_t)((v[l] >> j) & 1) << l; bits[j] = w; }
}
void bs_transpose_out(uint32_t *v, const uint64_t *bits, int lanes, int n) {
    for (int l = 0; l < lanes; l++) { uint32_t x = 0; for (int j = 0; j < n; j++) x |= (uint32_t)((bits[j] >> l) & 1) << j; v[l] = x; }
}
void bs_from_pub_lanes(mw *z, const uint32_t *v, int lanes, int n) {
    uint64_t b[32]; bs_transpose_in(b, v, lanes, n);
    for (int j = 0; j < n; j++) z[j] = mw_share(b[j]);
}

/* diagnostics: exceedance counts and arg-max of first-order |t| */
void mk_tvla_diag(int n, double thr, int *count, int *argmax, double *maxv) {
    int c = 0, am = -1; double m = 0;
    for (int i = 0; i < n; i++) {
        double m0 = A_sum[0][i] / A_n[0], m1 = A_sum[1][i] / A_n[1];
        double v0 = A_sq[0][i] / A_n[0] - m0 * m0, v1 = A_sq[1][i] / A_n[1] - m1 * m1;
        double den = sqrt(v0 / A_n[0] + v1 / A_n[1]);
        if (den <= 0) continue;
        double t = fabs((m0 - m1) / den);
        if (t > thr) c++;
        if (t > m) { m = t; am = i; }
    }
    *count = c; *argmax = am; *maxv = m;
}
