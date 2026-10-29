/*
 * AFA-12 v11 — first-order Boolean masking core
 * A masked 64-bit word is a pair of shares (x = s0 ^ s1). In bitsliced mode each of
 * the 64 bit positions is an independent lane (64 coefficients processed at once).
 */
#ifndef MASKED_CORE_H
#define MASKED_CORE_H
#include <stdint.h>

typedef struct { uint64_t s0, s1; } mw;

/* randomness: xorshift128+ for the prototype (production: CSPRNG / TRNG) */
extern uint64_t mk_rs0, mk_rs1;
static inline uint64_t mk_rnd(void) {
    uint64_t x = mk_rs0, y = mk_rs1; mk_rs0 = y; x ^= x << 23; mk_rs1 = x ^ y ^ (x >> 17) ^ (y >> 26); return mk_rs1 + y;
}
void mk_seed(uint64_t s);

/* leakage instrumentation (Hamming weight of each share; online TVLA accumulators) */
extern int mk_leak_on;
void mk_leak_word(uint64_t v);
#define MK_LEAK(m) do { if (mk_leak_on) { mk_leak_word((m).s0); mk_leak_word((m).s1); } } while (0)

static inline mw mw_pub(uint64_t c) { mw r = {c, 0}; return r; }
static inline mw mw_share(uint64_t v) { uint64_t r = mk_rnd(); mw m = {v ^ r, r}; return m; }
static inline mw mx(mw a, mw b) { mw r = {a.s0 ^ b.s0, a.s1 ^ b.s1}; return r; }
static inline mw mnot(mw a) { mw r = {~a.s0, a.s1}; return r; }
static inline mw mandp(mw a, uint64_t c) { mw r = {a.s0 & c, a.s1 & c}; return r; }
static inline mw mshl(mw a, int k) { mw r = {a.s0 << k, a.s1 << k}; return r; }
static inline mw mshr(mw a, int k) { mw r = {a.s0 >> k, a.s1 >> k}; return r; }
static inline mw mrol(mw a, int k) { if (!k) return a; mw r = {(a.s0 << k) | (a.s0 >> (64 - k)), (a.s1 << k) | (a.s1 >> (64 - k))}; return r; }
static inline uint64_t munmask(mw a) { return a.s0 ^ a.s1; }

/* ISW AND with one fresh random word */
static inline mw secand(mw a, mw b) {
    uint64_t r = mk_rnd();
    mw z;
    z.s0 = (a.s0 & b.s0) ^ r;
    uint64_t t = (r ^ (a.s0 & b.s1)) ^ (a.s1 & b.s0);
    z.s1 = (a.s1 & b.s1) ^ t;
    MK_LEAK(z);
    return z;
}
static inline mw secsel(mw m, mw a, mw b) { return mx(b, secand(m, mx(a, b))); }   /* per lane: m ? a : b */
static inline mw mrefresh(mw a) { uint64_t r = mk_rnd(); mw z = {a.s0 ^ r, a.s1 ^ r}; return z; }

/* ---------------- bitsliced masked integers: x[j] = bit j of every lane ---------------- */
void bs_add(mw *z, const mw *x, const mw *y, int n);              /* z = x + y mod 2^n */
void bs_addpub(mw *z, const mw *x, const uint64_t *k, int n);     /* z = x + k (k public, per-lane bitsliced) */
void bs_addconst(mw *z, const mw *x, uint32_t c, int n);          /* same constant in every lane */
void bs_neg(mw *z, const mw *x, int n);                           /* z = -x mod 2^n */
void bs_sel(mw *z, mw m, const mw *a, const mw *b, int n);        /* per lane m ? a : b */
void bs_from_pub_lanes(mw *z, const uint32_t *v, int lanes, int n); /* public/uniform values -> fresh Boolean shares */
void bs_transpose_in(uint64_t *bits, const uint32_t *v, int lanes, int n);
void bs_transpose_out(uint32_t *v, const uint64_t *bits, int lanes, int n);
#endif
