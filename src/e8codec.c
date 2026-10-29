/*
 * Gosset-KEM v4 — Constant-time E8 message codec for ML-KEM (q = 3329)
 * =================================================================
 * Message: 8 bits per block of 8 coefficients, encoded as a point of
 * (q/2)*E8 mod q*Z^8. Decoding = nearest E8 point (Conway–Sloane), done
 * in integer arithmetic with no secret-dependent branches or memory indices.
 *
 * Scaling: work with Y = 4*y (y = centered coefficient), lattice unit
 * U = 2q (= 4 * q/2). D8 candidate: U*f, f in Z^8, sum f even.
 * Half coset candidate: U*g + q (= 4 * q/2 * (g + 1/2)), sum g even.
 *
 * Message index m (0..255): bit 7 selects the coset (0: D8, 1: D8+1/2),
 * bits 0..6 are coordinate bits 0..6 mod 2; bit for coordinate 7 is the
 * parity of bits 0..6 (forced by the even-sum condition).
 */
#include <stdint.h>
#include "e8codec.h"

#define Q 3329
#define U (2 * Q)

/* constant-time helpers (arithmetic shifts on int64 are used as masks).
 * vb64 is a value barrier: clang -O3 otherwise turns the mask logic below into
 * secret-dependent branches (found by dudect on clang 18, |t| > 1000). */
static inline int64_t vb64(int64_t x) { __asm__ volatile("" : "+r"(x)); return x; }
static inline int64_t ct_abs(int64_t x) { int64_t s = vb64(x >> 63); return (x ^ s) - s; }
static inline int64_t ct_lt(int64_t a, int64_t b) { return vb64((a - b) >> 63); }      /* -1 if a<b else 0 */
static inline int64_t ct_select(int64_t mask, int64_t a, int64_t b) { return b ^ (vb64(mask) & (a ^ b)); }

/* round(Y / U) for |Y| <= 8U, no division instruction on secret data:
 * floor(Z / U) for Z = Y + U/2 + 8U in [0, 17U) via multiply-shift.
 * M = ceil(2^32 / U); exactness verified exhaustively by the test suite. */
static inline int64_t ct_round_div_U(int64_t Y) {
    const uint64_t M = (((uint64_t)1 << 32) + U - 1) / U;
    uint64_t Z = (uint64_t)(Y + U / 2 + 8 * U);
    return (int64_t)((Z * M) >> 32) - 8;
}

/* Nearest point of the D8-type candidate: f = round((Y - off)/U), fix parity
 * by moving the worst-rounded coordinate the other way. Returns squared
 * distance (scaled) and writes f. */
static inline int64_t d8_candidate(const int64_t Y[8], int64_t off, int64_t f[8]) {
    int64_t err[8], parity = 0, best = -1, bidx = 0, bsgn = 0, dist = 0;
    for (int i = 0; i < 8; i++) {
        f[i] = ct_round_div_U(Y[i] - off);
        err[i] = Y[i] - off - U * f[i];
        parity ^= f[i];
        int64_t a = ct_abs(err[i]);
        int64_t m = ct_lt(best, a);                      /* a > best */
        best = ct_select(m, a, best);
        bidx = ct_select(m, i, bidx);
        bsgn = ct_select(m, (err[i] >> 63) | 1, bsgn);   /* sign(err): +1 / -1 */
    }
    int64_t odd = vb64(-(parity & 1));                         /* all-ones if sum odd */
    for (int i = 0; i < 8; i++) {
        int64_t hit = odd & ~vb64((((int64_t)i - bidx) | (bidx - (int64_t)i)) >> 63); /* i == bidx */
        f[i] += hit & bsgn;
        int64_t e = Y[i] - off - U * f[i];
        dist += e * e;
    }
    return dist;
}

void e8_encode(const uint8_t msg[32], int16_t out[256]) {
    for (int b = 0; b < 32; b++) {
        uint8_t m = msg[b];
        int16_t coset = (int16_t)(m >> 7);
        int16_t par = 0;
        for (int i = 0; i < 8; i++) {
            int16_t bit = (i < 7) ? (int16_t)((m >> i) & 1) : par;
            par ^= bit;
            /* value = round(q/2 * x), x = bit + coset/2 in {0, 1/2, 1, 3/2} */
            int16_t v0 = (int16_t)(bit * 1665);               /* x in {0, 1}     */
            int16_t v1 = (int16_t)(832 + bit * 1665);         /* x in {1/2, 3/2} */
            int16_t mask = (int16_t)-coset;
            out[8 * b + i] = (int16_t)(v0 ^ (mask & (v0 ^ v1)));
        }
    }
}

void e8_decode(const int16_t in[256], uint8_t msg[32]) {
    for (int b = 0; b < 32; b++) {
        int64_t Y[8], f[8], g[8];
        for (int i = 0; i < 8; i++) {
            int64_t y = in[8 * b + i];
            y = y % Q; y += Q & (y >> 63);                    /* [0, q) */
            y -= Q & ct_lt(Q / 2, y);                         /* centered */
            Y[i] = 4 * y;
        }
        int64_t dA = d8_candidate(Y, 0, f);
        int64_t dB = d8_candidate(Y, Q, g);
        int64_t useB = ct_lt(dB, dA);
        uint8_t m = 0;
        for (int i = 0; i < 7; i++) {
            int64_t bit = ct_select(useB, g[i], f[i]) & 1;
            m |= (uint8_t)(bit << i);
        }
        m |= (uint8_t)((useB & 1) << 7);
        msg[b] = m;
    }
}

/* Reference ML-KEM 1-bit decoding (poly_tomsg), for benchmarking. */
void std_decode(const int16_t in[256], uint8_t msg[32]) {
    for (int i = 0; i < 32; i++) {
        msg[i] = 0;
        for (int j = 0; j < 8; j++) {
            uint32_t t = (uint32_t)((in[8 * i + j] % Q + Q) % Q);
            t <<= 1; t += 1665; t *= 80635; t >>= 28; t &= 1;   /* as in the Kyber reference */
            msg[i] |= (uint8_t)(t << j);
        }
    }
}

/* Deliberately LEAKY decoder (positive control for the timing test):
 * branches on secret data and exits early. Same output as e8_decode. */
void e8_decode_leaky(const int16_t in[256], uint8_t msg[32]) {
    for (int b = 0; b < 32; b++) {
        int64_t Y[8], f[8], g[8];
        int quiet = 1;
        for (int i = 0; i < 8; i++) {
            int64_t y = in[8 * b + i] % Q; if (y < 0) y += Q; if (y > Q / 2) y -= Q;
            Y[i] = 4 * y;
            if (Y[i] != 0 && Y[i] != -4 * 1664) quiet = 0;   /* exact D8 codeword */
        }
        if (quiet) {                                  /* secret-dependent shortcut */
            uint8_t m = 0;
            for (int i = 0; i < 7; i++) m |= (uint8_t)((Y[i] != 0) << i);
            msg[b] = m;
            continue;
        }
        int64_t dA = d8_candidate(Y, 0, f), dB = d8_candidate(Y, Q, g);
        uint8_t m = 0;
        if (dB < dA) { for (int i = 0; i < 7; i++) m |= (uint8_t)((g[i] & 1) << i); m |= 0x80; }
        else         { for (int i = 0; i < 7; i++) m |= (uint8_t)((f[i] & 1) << i); }
        msg[b] = m;
    }
}

int64_t e8_round_div_U_test(int64_t Y) { return ct_round_div_U(Y); }
