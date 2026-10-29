/*
 * AFA-12 v9 — X25519 (RFC 7748), own implementation
 * Field GF(2^255-19) in radix 2^51 (5 x 64-bit limbs, 128-bit products).
 * Montgomery ladder with constant-time conditional swap; no secret-dependent
 * branches or memory indices.
 */
#include <stdint.h>
#include <string.h>
#include "x25519.h"

typedef uint64_t fe[5];
typedef unsigned __int128 u128;
#define M51 ((1ULL << 51) - 1)

static uint64_t ld64(const uint8_t *p) { uint64_t r = 0; for (int i = 7; i >= 0; i--) r = (r << 8) | p[i]; return r; }

static void fe_frombytes(fe h, const uint8_t s[32]) {
    uint64_t a = ld64(s), b = ld64(s + 8), c = ld64(s + 16), d = ld64(s + 24);
    h[0] = a & M51;
    h[1] = ((a >> 51) | (b << 13)) & M51;
    h[2] = ((b >> 38) | (c << 26)) & M51;
    h[3] = ((c >> 25) | (d << 39)) & M51;
    h[4] = (d >> 12) & M51;                         /* bit 255 ignored (RFC 7748) */
}

static void fe_carry(fe h) {
    uint64_t c;
    c = h[0] >> 51; h[0] &= M51; h[1] += c;
    c = h[1] >> 51; h[1] &= M51; h[2] += c;
    c = h[2] >> 51; h[2] &= M51; h[3] += c;
    c = h[3] >> 51; h[3] &= M51; h[4] += c;
    c = h[4] >> 51; h[4] &= M51; h[0] += 19 * c;
    c = h[0] >> 51; h[0] &= M51; h[1] += c;
}

static void fe_tobytes(uint8_t s[32], const fe f) {
    fe h; memcpy(h, f, sizeof h);
    fe_carry(h); fe_carry(h);
    uint64_t q = (h[0] + 19) >> 51;                 /* q = 1 iff h >= p */
    q = (h[1] + q) >> 51; q = (h[2] + q) >> 51; q = (h[3] + q) >> 51; q = (h[4] + q) >> 51;
    h[0] += 19 * q;
    uint64_t c;
    c = h[0] >> 51; h[0] &= M51; h[1] += c;
    c = h[1] >> 51; h[1] &= M51; h[2] += c;
    c = h[2] >> 51; h[2] &= M51; h[3] += c;
    c = h[3] >> 51; h[3] &= M51; h[4] += c;
    h[4] &= M51;
    uint64_t w0 = h[0] | (h[1] << 51), w1 = (h[1] >> 13) | (h[2] << 38), w2 = (h[2] >> 26) | (h[3] << 25), w3 = (h[3] >> 39) | (h[4] << 12);
    uint64_t w[4] = {w0, w1, w2, w3};
    for (int i = 0; i < 4; i++) for (int j = 0; j < 8; j++) s[8 * i + j] = (uint8_t)(w[i] >> (8 * j));
}

static void fe_add(fe h, const fe f, const fe g) { for (int i = 0; i < 5; i++) h[i] = f[i] + g[i]; }
static void fe_sub(fe h, const fe f, const fe g) {       /* f + 4p - g, limbs stay positive */
    h[0] = f[0] + 0x1FFFFFFFFFFFB4ULL - g[0];
    for (int i = 1; i < 5; i++) h[i] = f[i] + 0x1FFFFFFFFFFFFCULL - g[i];
    fe_carry(h);
}
static void fe_mul(fe h, const fe f, const fe g) {
    u128 f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4];
    uint64_t g0 = g[0], g1 = g[1], g2 = g[2], g3 = g[3], g4 = g[4];
    uint64_t g1_19 = 19 * g1, g2_19 = 19 * g2, g3_19 = 19 * g3, g4_19 = 19 * g4;
    u128 r0 = f0 * g0 + f1 * g4_19 + f2 * g3_19 + f3 * g2_19 + f4 * g1_19;
    u128 r1 = f0 * g1 + f1 * g0 + f2 * g4_19 + f3 * g3_19 + f4 * g2_19;
    u128 r2 = f0 * g2 + f1 * g1 + f2 * g0 + f3 * g4_19 + f4 * g3_19;
    u128 r3 = f0 * g3 + f1 * g2 + f2 * g1 + f3 * g0 + f4 * g4_19;
    u128 r4 = f0 * g4 + f1 * g3 + f2 * g2 + f3 * g1 + f4 * g0;
    uint64_t c;
    r1 += (uint64_t)(r0 >> 51); h[0] = (uint64_t)r0 & M51;
    r2 += (uint64_t)(r1 >> 51); h[1] = (uint64_t)r1 & M51;
    r3 += (uint64_t)(r2 >> 51); h[2] = (uint64_t)r2 & M51;
    r4 += (uint64_t)(r3 >> 51); h[3] = (uint64_t)r3 & M51;
    c = (uint64_t)(r4 >> 51); h[4] = (uint64_t)r4 & M51;
    h[0] += 19 * c;
    c = h[0] >> 51; h[0] &= M51; h[1] += c;
}
static void fe_sq(fe h, const fe f) { fe_mul(h, f, f); }
static void fe_mul_small(fe h, const fe f, uint64_t k) {
    u128 c = 0;
    for (int i = 0; i < 5; i++) { c += (u128)f[i] * k; h[i] = (uint64_t)c & M51; c >>= 51; }
    h[0] += 19 * (uint64_t)c;
    fe_carry(h);
}
static void fe_sqn(fe h, const fe f, int n) { fe_sq(h, f); for (int i = 1; i < n; i++) fe_sq(h, h); }
static void fe_invert(fe out, const fe z) {               /* z^(p-2), ref10 addition chain */
    fe z2, z9, z11, z2_5_0, z2_10_0, z2_20_0, z2_50_0, z2_100_0, t;
    fe_sq(z2, z); fe_sqn(t, z2, 2); fe_mul(z9, t, z); fe_mul(z11, z9, z2);
    fe_sq(t, z11); fe_mul(z2_5_0, t, z9);
    fe_sqn(t, z2_5_0, 5); fe_mul(z2_10_0, t, z2_5_0);
    fe_sqn(t, z2_10_0, 10); fe_mul(z2_20_0, t, z2_10_0);
    fe_sqn(t, z2_20_0, 20); fe_mul(t, t, z2_20_0);
    fe_sqn(t, t, 10); fe_mul(z2_50_0, t, z2_10_0);
    fe_sqn(t, z2_50_0, 50); fe_mul(z2_100_0, t, z2_50_0);
    fe_sqn(t, z2_100_0, 100); fe_mul(t, t, z2_100_0);
    fe_sqn(t, t, 50); fe_mul(t, t, z2_50_0);
    fe_sqn(t, t, 5); fe_mul(out, t, z11);
}
static void fe_cswap(fe a, fe b, uint64_t swap) {
    uint64_t mask = 0 - swap;
    __asm__ volatile("" : "+r"(mask));
    for (int i = 0; i < 5; i++) { uint64_t x = mask & (a[i] ^ b[i]); a[i] ^= x; b[i] ^= x; }
}

void x25519(uint8_t out[32], const uint8_t scalar[32], const uint8_t point[32]) {
    uint8_t k[32]; memcpy(k, scalar, 32);
    k[0] &= 248; k[31] &= 127; k[31] |= 64;
    fe x1, x2 = {1}, z2 = {0}, x3, z3 = {1}, A, AA, B, BB, E, C, D, DA, CB, t;
    fe_frombytes(x1, point); memcpy(x3, x1, sizeof x3);
    uint64_t swap = 0;
    for (int pos = 254; pos >= 0; pos--) {
        uint64_t bit = (k[pos >> 3] >> (pos & 7)) & 1;
        swap ^= bit; fe_cswap(x2, x3, swap); fe_cswap(z2, z3, swap); swap = bit;
        fe_add(A, x2, z2); fe_sq(AA, A); fe_sub(B, x2, z2); fe_sq(BB, B); fe_sub(E, AA, BB);
        fe_add(C, x3, z3); fe_sub(D, x3, z3); fe_mul(DA, D, A); fe_mul(CB, C, B);
        fe_add(t, DA, CB); fe_sq(x3, t);
        fe_sub(t, DA, CB); fe_sq(t, t); fe_mul(z3, x1, t);
        fe_mul(x2, AA, BB);
        fe_mul_small(t, E, 121665); fe_add(t, AA, t); fe_mul(z2, E, t);
    }
    fe_cswap(x2, x3, swap); fe_cswap(z2, z3, swap);
    fe_invert(z2, z2); fe_mul(x2, x2, z2);
    fe_tobytes(out, x2);
    volatile uint8_t *vk = k; for (int i = 0; i < 32; i++) vk[i] = 0;
}

void x25519_base(uint8_t out[32], const uint8_t scalar[32]) {
    uint8_t nine[32] = {9}; x25519(out, scalar, nine);
}
