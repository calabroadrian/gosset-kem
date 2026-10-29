/*
 * AFA-12 v10 — First-order masked E8 decoder (prototype)
 * ======================================================
 * Input : arithmetic shares of each coefficient, w = (w0 + w1) mod q   (w0 uniform)
 * Output: Boolean shares of the message byte,     m = m0 XOR m1
 * The secret never appears unshared. Every value is a pair of 32-bit Boolean
 * shares (x = s0 ^ s1); arithmetic is done modulo 2^NB (NB = 20 bits, values < 2^18).
 *
 * Gadgets
 *   SecAND : ISW multiplication with one fresh random word per call
 *   SecAdd : Boolean-masked addition, carry chain built from SecAND (NB iterations)
 *   A2B    : each arithmetic share (individually uniform) is Boolean-shared, then
 *            SecAdd + masked conditional subtraction of q
 * Decoder: multiplication-free reformulation (verified identical to the reference):
 *   y centered, Y = 4y, s_i = [Y<0], |Y|, t_i = [|Y| > q], a_i = t ? 2q-|Y| : |Y|
 *   parA = xor t, parB = xor s, one-hot argmax/argmin of a (first index on ties)
 *   D = sum a - 4q + parA (2q - 2 amax) - parB (2 amin);  useB = [D > 0]
 *   bits = useB ? s ^ (parB & hB) : t ^ (parA & hA);  byte = bits[0..6] | useB<<7
 * Randomness: xorshift PRNG for the prototype (production: CSPRNG/TRNG).
 * Leakage hook: LEAK() records the Hamming weight of each share separately
 * (first-order value model) for the simulation in power_masked.py.
 */
#include <stdint.h>
#include <string.h>

#define Q 3329
#define NB 20
#define NBMASK ((1u << NB) - 1)

typedef struct { uint32_t s0, s1; } m32;

/* ---------------- randomness and leakage trace ---------------- */
static uint64_t rs = 0x243F6A8885A308D3ULL;
static inline uint32_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (uint32_t)rs; }
static uint16_t *TR = 0; static int TRN = 0, TRMAX = 0;
static inline void leak1(uint32_t v) { if (TR && TRN < TRMAX) TR[TRN++] = (uint16_t)__builtin_popcount(v & NBMASK); }
#define LEAK(m) do { leak1((m).s0); leak1((m).s1); } while (0)

/* ---------------- linear operations (share-wise) ---------------- */
static inline m32 mx(m32 a, m32 b) { m32 r = {a.s0 ^ b.s0, a.s1 ^ b.s1}; return r; }
static inline m32 mnot(m32 a) { m32 r = {~a.s0, a.s1}; return r; }
static inline m32 mshl(m32 a, int k) { m32 r = {a.s0 << k, a.s1 << k}; return r; }
static inline m32 mandc(m32 a, uint32_t c) { m32 r = {a.s0 & c, a.s1 & c}; return r; }
static inline m32 mconst(uint32_t c) { uint32_t r = rnd(); m32 m = {c ^ r, r}; return m; }
static inline m32 msignmask(m32 a) {                 /* all-ones iff bit NB-1 set */
    m32 r = {(uint32_t)-(int32_t)((a.s0 >> (NB - 1)) & 1), (uint32_t)-(int32_t)((a.s1 >> (NB - 1)) & 1)}; return r;
}
static inline m32 mbit0(m32 a) { return mandc(a, 1); }
static inline m32 mbroadcast(m32 bit) { m32 r = {(uint32_t)-(int32_t)(bit.s0 & 1), (uint32_t)-(int32_t)(bit.s1 & 1)}; return r; }

/* ---------------- nonlinear gadgets ---------------- */
static inline m32 secand(m32 a, m32 b) {
    uint32_t r = rnd();
    m32 z;
    z.s0 = (a.s0 & b.s0) ^ r;
    uint32_t t = (r ^ (a.s0 & b.s1)) ^ (a.s1 & b.s0);
    z.s1 = (a.s1 & b.s1) ^ t;
    LEAK(z);
    return z;
}
static m32 secadd(m32 a, m32 b) {                    /* (a + b) mod 2^NB, Boolean-masked */
    m32 w = secand(a, b), x = mx(a, b), u = {0, 0};
    for (int k = 0; k < NB - 1; k++) { u = mshl(mx(secand(u, x), w), 1); }
    m32 z = mx(x, u); LEAK(z);
    return z;
}
static inline m32 secsel(m32 mask, m32 a, m32 b) {  /* mask ? a : b  */
    return mx(b, secand(mask, mx(a, b)));
}
static inline m32 secneg(m32 a) { return secadd(mnot(a), mconst(1)); }

/* arithmetic shares mod q  ->  Boolean shares of w in [0,q) */
static m32 a2b_modq(uint16_t w0, uint16_t w1) {
    m32 b0 = mconst(w0), b1 = mconst(w1);            /* each share alone is uniform: safe to share */
    m32 z = secadd(b0, b1);                          /* in [0, 2q) */
    m32 zq = secadd(z, mconst((uint32_t)(-Q) & NBMASK));
    return secsel(msignmask(zq), z, zq);             /* z - q < 0 ? z : z - q */
}

/* masked decode of one 8-coefficient block; returns Boolean shares of the byte */
void e8_block_masked(const uint16_t w0[8], const uint16_t w1[8], uint8_t *m0, uint8_t *m1) {
    m32 a[8], sumA = {0, 0}, tvec = {0, 0}, svec = {0, 0}, parA = {0, 0}, parB = {0, 0};
    for (int i = 0; i < 8; i++) {
        m32 w = a2b_modq(w0[i], w1[i]);
        m32 c = secadd(w, mconst((uint32_t)(-(Q / 2 + 1)) & NBMASK));       /* w - 1665 */
        m32 y = secsel(msignmask(c), w, secadd(w, mconst((uint32_t)(-Q) & NBMASK)));  /* centered */
        m32 Y = mshl(y, 2); LEAK(Y);
        m32 sm = msignmask(Y);                                              /* s_i mask */
        m32 absY = secadd(mx(Y, sm), mbit0(sm));                            /* |Y| */
        m32 tq = secadd(mconst((uint32_t)Q), secneg(absY));                 /* q - |Y| */
        m32 tm = msignmask(tq);                                             /* t_i = [|Y| > q] */
        m32 alt = secadd(mconst(2u * Q), secneg(absY));                     /* 2q - |Y| */
        a[i] = secsel(tm, alt, absY); LEAK(a[i]);
        sumA = secadd(sumA, a[i]);
        tvec = mx(tvec, mshl(mbit0(tm), i)); svec = mx(svec, mshl(mbit0(sm), i));
        parA = mx(parA, mbit0(tm)); parB = mx(parB, mbit0(sm));
    }
    /* first-index argmax / argmin of a, as masked one-hot bytes */
    m32 amax = a[0], amin = a[0], hA = mconst(1), hB = mconst(1);
    for (int i = 1; i < 8; i++) {
        m32 gt = msignmask(secadd(amax, secneg(a[i])));                     /* amax - a_i < 0  <=> a_i > amax */
        amax = secsel(gt, a[i], amax); hA = secsel(gt, mconst(1u << i), hA);
        m32 lt = msignmask(secadd(a[i], secneg(amin)));                     /* a_i - amin < 0 */
        amin = secsel(lt, a[i], amin); hB = secsel(lt, mconst(1u << i), hB);
    }
    m32 pA = mbroadcast(parA), pB = mbroadcast(parB);
    m32 T1 = secadd(mconst(2u * Q), secneg(mshl(amax, 1)));                 /* 2q - 2 amax */
    m32 T2 = mshl(amin, 1);                                                 /* 2 amin */
    m32 D = secadd(sumA, mconst((uint32_t)(-4 * Q) & NBMASK));
    D = secadd(D, secand(pA, T1));
    D = secadd(D, secneg(secand(pB, T2)));
    m32 useB = msignmask(secneg(D));                                        /* -D < 0  <=> D > 0 */
    m32 bA = mx(tvec, secand(pA, hA)), bB = mx(svec, secand(pB, hB));
    m32 bits = secsel(useB, bB, bA);
    m32 byte = mx(mandc(bits, 0x7F), mandc(useB, 0x80)); LEAK(byte);
    *m0 = (uint8_t)byte.s0; *m1 = (uint8_t)byte.s1;
}

/* full message: 32 blocks; optional shuffling of block order (hiding countermeasure) */
void e8_decode_masked(const uint16_t w0[256], const uint16_t w1[256], uint8_t m0[32], uint8_t m1[32], int shuffle) {
    uint8_t order[32]; for (int i = 0; i < 32; i++) order[i] = (uint8_t)i;
    if (shuffle) for (int i = 31; i > 0; i--) { int j = (int)(rnd() % (uint32_t)(i + 1)); uint8_t t = order[i]; order[i] = order[j]; order[j] = t; }
    for (int k = 0; k < 32; k++) { int b = order[k]; e8_block_masked(w0 + 8 * b, w1 + 8 * b, &m0[b], &m1[b]); }
}

/* test hooks */
void masked_seed(uint64_t s) { rs = s ? s : 1; }
void masked_trace(uint16_t *buf, int max) { TR = buf; TRN = 0; TRMAX = max; }
int masked_trace_len(void) { return TRN; }
