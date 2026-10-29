/*
 * Gosset-KEM v4 — optimized constant-time E8 decoder
 * Input contract (same as ML-KEM poly_tomsg after reduction): coefficients
 * in [0, q). Observations used:
 *  - Y = 4*y_centered lies in (-U, U), so round(Y/U) is in {-1,0,1} and
 *    round((Y-q)/U) in {-1,0}: rounding = 2 constant-time comparisons, no
 *    multiplication or division.
 *  - int32 arithmetic everywhere (max squared distance < 2^31).
 *  - scalar version (portable) + AVX2 version: one 8-coordinate block =
 *    one 256-bit register (8 x int32).
 */
#include <stdint.h>
#include <immintrin.h>

#define Q 3329
#define U (2 * Q)

/* value barrier: keeps the optimizer from turning masks into secret-dependent branches */
static inline int32_t vb32(int32_t x) { __asm__ volatile("" : "+r"(x)); return x; }
static inline int32_t lt(int32_t a, int32_t b) { return vb32((a - b) >> 31); }   /* -1 if a<b */

static inline int32_t cand(const int32_t Y[8], int32_t off, int32_t f[8]) {
    int32_t err[8], par = 0, best = -1, bidx = 0, bsg = 0, d = 0;
    for (int i = 0; i < 8; i++) {
        int32_t X = Y[i] - off;
        /* round half up of X/U for X in (-1.5U, U) */
        f[i] = -lt(U / 2 - 1, X) + lt(X, -U / 2) + lt(X, -3 * U / 2);
        err[i] = X - U * f[i];
        par ^= f[i];
        int32_t s = err[i] >> 31, a = (err[i] ^ s) - s;   /* branch-free abs */
        int32_t m = lt(best, a);
        best ^= m & (best ^ a);
        bidx ^= m & (bidx ^ i);
        bsg ^= m & (bsg ^ ((err[i] >> 31) | 1));
    }
    int32_t odd = vb32(-(par & 1));
    for (int i = 0; i < 8; i++) {
        int32_t hit = odd & ~(((i - bidx) | (bidx - i)) >> 31);
        f[i] += hit & bsg;
        int32_t e = Y[i] - off - U * f[i];
        d += e * e;
    }
    return d;
}

void e8_decode_fast(const int16_t in[256], uint8_t msg[32]) {
    for (int b = 0; b < 32; b++) {
        int32_t Y[8], f[8], g[8];
        for (int i = 0; i < 8; i++) {
            int32_t y = in[8 * b + i];
            y -= Q & lt(Q / 2, y);
            Y[i] = 4 * y;
        }
        int32_t dA = cand(Y, 0, f), dB = cand(Y, Q, g);
        int32_t useB = lt(dB, dA);
        uint32_t m = 0;
        for (int i = 0; i < 7; i++) m |= (uint32_t)(((f[i] ^ (useB & (f[i] ^ g[i]))) & 1) << i);
        msg[b] = (uint8_t)(m | ((uint32_t)(useB & 1) << 7));
    }
}

/* ---------------- AVX2: one block per register ---------------- */
static inline __m256i cand_avx(__m256i Y, int32_t off, __m256i *fo, int32_t *dist) {
    const __m256i vU = _mm256_set1_epi32(U), one = _mm256_set1_epi32(1);
    __m256i X = _mm256_sub_epi32(Y, _mm256_set1_epi32(off));
    __m256i up  = _mm256_cmpgt_epi32(X, _mm256_set1_epi32(U / 2 - 1));     /* X >= U/2  -> -1 */
    __m256i dn1 = _mm256_cmpgt_epi32(_mm256_set1_epi32(-U / 2), X);       /* X < -U/2  -> -1 */
    __m256i dn2 = _mm256_cmpgt_epi32(_mm256_set1_epi32(-3 * U / 2), X);   /* X < -3U/2 -> -1 */
    __m256i f = _mm256_add_epi32(_mm256_sub_epi32(dn1, up), dn2);         /* f = [X>=U/2] - [X<-U/2] - [X<-3U/2] */
    __m256i err = _mm256_sub_epi32(X, _mm256_mullo_epi32(vU, f));
    __m256i a = _mm256_abs_epi32(err);
    /* horizontal max of a, then first lane equal to max */
    __m256i mx = _mm256_max_epi32(a, _mm256_permutevar8x32_epi32(a, _mm256_setr_epi32(4,5,6,7,0,1,2,3)));
    mx = _mm256_max_epi32(mx, _mm256_shuffle_epi32(mx, 0x4E));
    mx = _mm256_max_epi32(mx, _mm256_shuffle_epi32(mx, 0xB1));
    __m256i eq = _mm256_cmpeq_epi32(a, mx);
    uint32_t eqbits = (uint32_t)_mm256_movemask_ps(_mm256_castsi256_ps(eq));
    uint32_t first = eqbits & (~eqbits + 1);                               /* lowest set bit */
    __m256i lanebit = _mm256_setr_epi32(1, 2, 4, 8, 16, 32, 64, 128);
    __m256i sel = _mm256_cmpeq_epi32(_mm256_and_si256(lanebit, _mm256_set1_epi32((int)first)), lanebit);
    /* parity of sum f */
    __m256i s = _mm256_add_epi32(f, _mm256_permutevar8x32_epi32(f, _mm256_setr_epi32(4,5,6,7,0,1,2,3)));
    s = _mm256_add_epi32(s, _mm256_shuffle_epi32(s, 0x4E));
    s = _mm256_add_epi32(s, _mm256_shuffle_epi32(s, 0xB1));
    __m256i odd = _mm256_sub_epi32(_mm256_setzero_si256(), _mm256_and_si256(s, one)); /* all-ones if odd */
    __m256i sgn = _mm256_or_si256(_mm256_srai_epi32(err, 31), one);       /* +1 / -1 */
    f = _mm256_add_epi32(f, _mm256_and_si256(_mm256_and_si256(odd, sel), sgn));
    err = _mm256_sub_epi32(X, _mm256_mullo_epi32(vU, f));
    __m256i sq = _mm256_mullo_epi32(err, err);
    sq = _mm256_add_epi32(sq, _mm256_permutevar8x32_epi32(sq, _mm256_setr_epi32(4,5,6,7,0,1,2,3)));
    sq = _mm256_add_epi32(sq, _mm256_shuffle_epi32(sq, 0x4E));
    sq = _mm256_add_epi32(sq, _mm256_shuffle_epi32(sq, 0xB1));
    *dist = _mm256_cvtsi256_si32(sq);
    *fo = f;
    return f;
}

void e8_decode_avx2(const int16_t in[256], uint8_t msg[32]) {
    const __m256i vQ = _mm256_set1_epi32(Q), half = _mm256_set1_epi32(Q / 2);
    for (int b = 0; b < 32; b++) {
        __m256i y = _mm256_cvtepi16_epi32(_mm_loadu_si128((const __m128i *)(in + 8 * b)));
        y = _mm256_sub_epi32(y, _mm256_and_si256(vQ, _mm256_cmpgt_epi32(y, half)));
        __m256i Y = _mm256_slli_epi32(y, 2);
        __m256i f, g; int32_t dA, dB;
        cand_avx(Y, 0, &f, &dA);
        cand_avx(Y, Q, &g, &dB);
        int32_t useB = lt(dB, dA);
        __m256i mB = _mm256_set1_epi32(useB);
        __m256i h = _mm256_xor_si256(f, _mm256_and_si256(mB, _mm256_xor_si256(f, g)));
        h = _mm256_slli_epi32(_mm256_and_si256(h, _mm256_set1_epi32(1)), 31);
        uint32_t bits = (uint32_t)_mm256_movemask_ps(_mm256_castsi256_ps(h)) & 0x7f;
        msg[b] = (uint8_t)(bits | ((uint32_t)(useB & 1) << 7));
    }
}

/* Kyber reference poly_tomsg (coefficients in [0,q)) — fair baseline */
void std_decode_ref(const int16_t in[256], uint8_t msg[32]) {
    for (int i = 0; i < 32; i++) {
        msg[i] = 0;
        for (int j = 0; j < 8; j++) {
            uint32_t t = (uint32_t)in[8 * i + j];
            t <<= 1; t += 1665; t *= 80635; t >>= 28; t &= 1;
            msg[i] |= (uint8_t)(t << j);
        }
    }
}
