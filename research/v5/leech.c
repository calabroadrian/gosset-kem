/*
 * AFA-12 v5 — Constant-time maximum-likelihood Leech lattice codec for ML-KEM
 * ===========================================================================
 * Code: L = (q/8) * Lambda24 (integral coordinates, min norm^2 32, 8Z^24 in L).
 * Lambda24 = U_{c in Golay} (2c + 4 D24)  U  (1 + 2c + 4 (D24 + e1))
 *   -> point p = h + 2c + 4z, h in {0,1}, c Golay codeword, sum(z) == h (mod 2)
 * Golay via Turyn: c = (a^x, b^x, a^b^x), a,b in H, x in H' (two Hamming [8,4,4]).
 * Message per block (36 bits): h (1) | x,a,b indices (12) | z_0..z_22 mod 2 (23).
 *
 * Decoder: exact nearest point. For every half h and every coordinate i and
 * codeword bit b, the nearest point of h + 2b + 4Z is found with comparisons;
 * per-octet tables (cost sum, z-parity, cheapest parity flip) for all 256
 * octet patterns; then all 2 x 4096 (h, codeword) candidates are scored with
 * three table reads indexed by PUBLIC loop counters. Selection, final lookups
 * and parity repair use masks only: no secret-dependent branches or indices.
 */
#include <stdint.h>
#include "leech.h"

#define Q 3329

static const uint8_t HH[16] = {0, 85, 15, 90, 60, 105, 51, 102, 240, 165, 255, 170, 204, 153, 195, 150};
static const uint8_t HP[16] = {0, 108, 58, 86, 75, 39, 113, 29, 197, 169, 255, 147, 142, 226, 180, 216};

/* value barrier: hides mask values from the optimizer so it cannot turn
   masked selections back into secret-dependent branches (found by dudect:
   gcc -O3 emitted je/cmove on the selection masks without it) */
static inline int64_t vb(int64_t x) { __asm__ volatile("" : "+r"(x)); return x; }
static inline int64_t m64(int64_t c) { return vb(-c); }             /* 0/1 -> 0/-1 */
static inline int64_t lt64(int64_t a, int64_t b) { return vb((a - b) >> 63); }
static inline int64_t sel64(int64_t m, int64_t a, int64_t b) { return b ^ (vb(m) & (a ^ b)); }

static inline int bit_of(uint32_t w24, int i) {                    /* coordinate i of 24-bit word */
    return (int)((w24 >> (23 - i)) & 1);
}

static uint8_t ct_lookup16(const uint8_t *t, int64_t idx) {
    uint8_t r = 0;
    for (int j = 0; j < 16; j++) {
        int64_t eq = ~(((j - idx) | (idx - j)) >> 63);
        r |= (uint8_t)(eq & t[j]);
    }
    return r;
}

void leech_encode_block(uint64_t m, int16_t out[24]) {
    int h = (int)(m & 1);
    int xi = (int)((m >> 1) & 15), ai = (int)((m >> 5) & 15), bi = (int)((m >> 9) & 15);
    uint32_t c = ((uint32_t)(HH[ai] ^ HP[xi]) << 16) | ((uint32_t)(HH[bi] ^ HP[xi]) << 8) |
                 (uint32_t)(HH[ai] ^ HH[bi] ^ HP[xi]);
    int zs = 0;
    for (int i = 0; i < 24; i++) {
        int z = (i < 23) ? (int)((m >> (13 + i)) & 1) : ((h - zs) & 1);
        zs += z;
        int p = h + 2 * bit_of(c, i) + 4 * z;                   /* 0..7 */
        out[i] = (int16_t)((Q * p + 4) / 8 % Q);                  /* round(q p / 8) */
    }
}

uint64_t leech_decode_block(const int16_t in[24]) {
    int64_t cost[2][24][2], flip[2][24][2], kk[2][24][2], sg[2][24][2];
    int64_t par[2][24][2];
    for (int i = 0; i < 24; i++) {
        int64_t y = in[i];
        y -= Q & lt64(Q / 2, y);
        int64_t Y = 8 * y;
        for (int h = 0; h < 2; h++)
            for (int b = 0; b < 2; b++) {
                int64_t t = Y - (int64_t)Q * (2 * b + h);
                /* k = round((t/q)/4) in {-2..1}: [t>=2q] - [t<-2q] - [t<-6q] */
                int64_t k = (1 + lt64(t, 2 * Q)) + lt64(t, -2 * Q) + lt64(t, -6 * Q);
                int64_t e = t - 4 * (int64_t)Q * k;
                int64_t s = (e >> 63) | 1;                          /* +1 if e>=0 else -1 */
                int64_t e2 = e - 4 * (int64_t)Q * s;
                cost[h][i][b] = e * e;
                flip[h][i][b] = e2 * e2 - e * e;
                kk[h][i][b] = k;
                sg[h][i][b] = s;
                par[h][i][b] = k & 1;
            }
    }
    /* per-octet tables over all 256 patterns */
    static __thread int64_t Tc[2][3][256], Tf[2][3][256];
    static __thread int64_t Tp[2][3][256];
    for (int h = 0; h < 2; h++)
        for (int o = 0; o < 3; o++)
            for (int w = 0; w < 256; w++) {
                int64_t c = 0, p = 0, f = INT64_MAX;
                for (int j = 0; j < 8; j++) {
                    int b = (w >> (7 - j)) & 1, i = 8 * o + j;
                    c += cost[h][i][b];
                    p ^= par[h][i][b];
                    int64_t fl = flip[h][i][b];
                    f = sel64(lt64(fl, f), fl, f);
                }
                Tc[h][o][w] = c; Tp[h][o][w] = p; Tf[h][o][w] = f;
            }
    /* score all 2 x 4096 candidates (public loop indices) */
    int64_t best = INT64_MAX, bh = 0, bx = 0, ba = 0, bb = 0;
    for (int h = 0; h < 2; h++)
        for (int x = 0; x < 16; x++)
            for (int a = 0; a < 16; a++)
                for (int b = 0; b < 16; b++) {
                    int w1 = HH[a] ^ HP[x], w2 = HH[b] ^ HP[x], w3 = HH[a] ^ HH[b] ^ HP[x];
                    int64_t d = Tc[h][0][w1] + Tc[h][1][w2] + Tc[h][2][w3];
                    int64_t mis = (Tp[h][0][w1] ^ Tp[h][1][w2] ^ Tp[h][2][w3] ^ h) & 1;
                    int64_t f = Tf[h][0][w1];
                    f = sel64(lt64(Tf[h][1][w2], f), Tf[h][1][w2], f);
                    f = sel64(lt64(Tf[h][2][w3], f), Tf[h][2][w3], f);
                    d += m64(mis) & f;
                    int64_t m = lt64(d, best);
                    best = sel64(m, d, best);
                    bh = sel64(m, h, bh); bx = sel64(m, x, bx); ba = sel64(m, a, ba); bb = sel64(m, b, bb);
                }
    /* reconstruct winning codeword (masked lookups) */
    uint8_t hx = ct_lookup16(HP, bx), ha = ct_lookup16(HH, ba), hb = ct_lookup16(HH, bb);
    uint32_t c = ((uint32_t)(ha ^ hx) << 16) | ((uint32_t)(hb ^ hx) << 8) | (uint32_t)(ha ^ hb ^ hx);
    /* per-coordinate z with masked selection over (h, b), then parity repair */
    int64_t z[24], fl[24], s[24], ptot = 0;
    for (int i = 0; i < 24; i++) {
        int64_t b = bit_of(c, i), mh = m64(bh), mb = m64(b);
        #define PICK(A) sel64(mh, sel64(mb, A[1][i][1], A[1][i][0]), sel64(mb, A[0][i][1], A[0][i][0]))
        z[i] = PICK(kk); fl[i] = PICK(flip); s[i] = PICK(sg);
        #undef PICK
        ptot ^= z[i];
    }
    int64_t mis = m64((ptot ^ bh) & 1), bf = INT64_MAX, bi = 0;
    for (int i = 0; i < 24; i++) {
        int64_t m = lt64(fl[i], bf);
        bf = sel64(m, fl[i], bf); bi = sel64(m, i, bi);
    }
    for (int i = 0; i < 24; i++) {
        int64_t hit = mis & ~((((int64_t)i - bi) | (bi - (int64_t)i)) >> 63);
        z[i] += hit & s[i];
    }
    uint64_t msg = (uint64_t)(bh & 1) | ((uint64_t)bx << 1) | ((uint64_t)ba << 5) | ((uint64_t)bb << 9);
    for (int i = 0; i < 23; i++) msg |= (uint64_t)(z[i] & 1) << (13 + i);
    return msg;
}

/* full message: 8 blocks over coefficients 0..191 (256-bit key uses 32 of 36 bits/block) */
void leech_decode(const int16_t in[192], uint64_t out[8]) {
    for (int b = 0; b < 8; b++) out[b] = leech_decode_block(in + 24 * b);
}
