/* AFA-12 v11 — first-order masked Keccak-f[1600], SHA3-512 and SHAKE256 (masked prefix || public suffix) */
#include <string.h>
#include "masked_core.h"

static const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808AULL, 0x8000000080008000ULL,
    0x000000000000808BULL, 0x0000000080000001ULL, 0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008AULL, 0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000AULL,
    0x000000008000808BULL, 0x800000000000008BULL, 0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL, 0x000000000000800AULL, 0x800000008000000AULL,
    0x8000000080008081ULL, 0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL};
static const int R[25] = {0, 1, 62, 28, 27, 36, 44, 6, 55, 20, 3, 10, 43, 25, 39, 41, 45, 15, 21, 8, 18, 2, 61, 56, 14};

void mkeccakf(mw a[25]) {
    for (int r = 0; r < 24; r++) {
        mw c[5], d[5], b[25];
        for (int x = 0; x < 5; x++) c[x] = mx(mx(mx(a[x], a[x + 5]), mx(a[x + 10], a[x + 15])), a[x + 20]);
        for (int x = 0; x < 5; x++) d[x] = mx(c[(x + 4) % 5], mrol(c[(x + 1) % 5], 1));
        for (int y = 0; y < 5; y++)
            for (int x = 0; x < 5; x++) {
                int X = y, Y = (2 * x + 3 * y) % 5;
                b[X + 5 * Y] = mrol(mx(a[x + 5 * y], d[x]), R[x + 5 * y]);
            }
        for (int y = 0; y < 5; y++)
            for (int x = 0; x < 5; x++)
                a[x + 5 * y] = mx(b[x + 5 * y], secand(mnot(b[(x + 1) % 5 + 5 * y]), b[(x + 2) % 5 + 5 * y]));
        a[0].s0 ^= RC[r];
    }
}

static void xor_byte(mw st[25], int pos, uint8_t v0, uint8_t v1) {
    st[pos / 8].s0 ^= (uint64_t)v0 << (8 * (pos % 8));
    st[pos / 8].s1 ^= (uint64_t)v1 << (8 * (pos % 8));
}
/* sponge: masked prefix (m0 ^ m1, mlen bytes) || public suffix (p, plen bytes) */
static void msponge(int rate, uint8_t pad, const uint8_t *m0, const uint8_t *m1, int mlen,
                    const uint8_t *p, int plen, uint8_t *o0, uint8_t *o1, int olen) {
    mw st[25]; memset(st, 0, sizeof st);
    int pos = 0;
    for (int i = 0; i < mlen + plen; i++) {
        if (i < mlen) xor_byte(st, pos, m0[i], m1[i]); else xor_byte(st, pos, p[i - mlen], 0);
        if (++pos == rate) { mkeccakf(st); pos = 0; }
    }
    xor_byte(st, pos, pad, 0); xor_byte(st, rate - 1, 0x80, 0);
    mkeccakf(st);
    pos = 0;
    for (int i = 0; i < olen; i++) {
        if (pos == rate) { mkeccakf(st); pos = 0; }
        o0[i] = (uint8_t)(st[pos / 8].s0 >> (8 * (pos % 8)));
        o1[i] = (uint8_t)(st[pos / 8].s1 >> (8 * (pos % 8)));
        pos++;
    }
}
void msha3_512(uint8_t o0[64], uint8_t o1[64], const uint8_t *m0, const uint8_t *m1, int mlen, const uint8_t *p, int plen) {
    msponge(72, 0x06, m0, m1, mlen, p, plen, o0, o1, 64);
}
void mshake256(uint8_t *o0, uint8_t *o1, int olen, const uint8_t *m0, const uint8_t *m1, int mlen, const uint8_t *p, int plen) {
    msponge(136, 0x1F, m0, m1, mlen, p, plen, o0, o1, olen);
}
