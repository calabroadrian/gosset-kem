/*
 * AFA-12 v11 — first-order masked decapsulation (E8 variants: Compact and MaxSec)
 * ==============================================================================
 * Secret data handled only as shares:
 *   s (secret key)       arithmetic shares mod q
 *   w = v - s^T u        arithmetic shares (linear, share-wise NTT/basemul)
 *   m' (message)         Boolean shares        <- bitsliced masked E8 decoder
 *   (K', r') = G(m'||h)  Boolean shares        <- masked SHA3-512
 *   PRF(r', N)           Boolean shares        <- masked SHAKE256
 *   y, e1, e2            arithmetic shares     <- masked CBD (bitsliced adder + B2A mod q)
 *   mu = Encode_E8(m')   arithmetic shares     <- 1-bit B2A mod q, linear encoding
 *   u', v'               arithmetic shares     (linear)
 *   c' == c ?            masked interval check per coefficient (A2B + comparisons with
 *                        PUBLIC decompression intervals of c), masked AND of all checks;
 *                        only the final 1-bit result is unmasked
 *   Kbar = J(z||c)       Boolean shares        <- masked SHAKE256 (z masked)
 * Output: the selected key is unmasked at the end (it is handed to the application).
 * Prototype notes: the secret key is masked at load time (production: store it masked
 * and refresh per call); PRNG is xorshift128+ (production: CSPRNG).
 */
#include "../v7/mlkem_opt.c"
#include "masked_core.h"

void msha3_512(uint8_t o0[64], uint8_t o1[64], const uint8_t *m0, const uint8_t *m1, int mlen, const uint8_t *p, int plen);
void mshake256(uint8_t *o0, uint8_t *o1, int olen, const uint8_t *m0, const uint8_t *m1, int mlen, const uint8_t *p, int plen);
void e8_decode_masked_bs(const uint16_t w0[256], const uint16_t w1[256], uint8_t m0[32], uint8_t m1[32]);
void a2b_modq_bs(mw *x, const uint32_t *v0, const uint32_t *v1, int lanes, int nbw);

typedef struct { poly s0, s1; } mpoly;            /* arithmetic shares, each in [0,q) */

static inline uint32_t rnd_modq(void) { return (uint32_t)(((mk_rnd() & 0xFFFFFFFFULL) * Q) >> 32); }

/* Boolean-masked bitsliced small integer (nv bits) -> arithmetic shares mod q */
static void b2a_modq_bs(int16_t *a0, int16_t *a1, const mw *v, int nv, int is_signed, int lanes) {
    mw x[14], qs[14], rb[14], s[14], sq[14], z[14];
    for (int j = 0; j < 14; j++) x[j] = j < nv ? v[j] : (is_signed ? v[nv - 1] : mw_pub(0));
    if (is_signed) {                                         /* negative -> add q (masked select by sign) */
        mw neg = x[13];
        for (int j = 0; j < 14; j++) qs[j] = ((Q >> j) & 1) ? neg : mw_pub(0);
        bs_add(x, x, qs, 14);
    }
    uint32_t r[64];
    for (int l = 0; l < lanes; l++) r[l] = rnd_modq();
    bs_from_pub_lanes(rb, r, lanes, 14);
    bs_add(s, x, rb, 14);                                    /* x + r in [0, 2q) */
    bs_addconst(sq, s, (uint32_t)(-Q) & 0x3FFF, 14);
    bs_sel(z, sq[13], s, sq, 14);                            /* (x + r) mod q, still masked */
    uint64_t zb[14]; uint32_t zv[64];
    for (int j = 0; j < 14; j++) zb[j] = munmask(z[j]);      /* uniform: reveals nothing */
    bs_transpose_out(zv, zb, lanes, 14);
    for (int l = 0; l < lanes; l++) { a0[l] = (int16_t)zv[l]; a1[l] = (int16_t)((Q - r[l]) % Q); }
}

/* masked CBD_eta from Boolean-masked PRF output */
static void cbd_masked(mpoly *f, const uint8_t *b0, const uint8_t *b1, int eta) {
    for (int g = 0; g < 4; g++) {
        mw A[4], B[4];
        for (int j = 0; j < eta; j++) {
            uint64_t a0 = 0, a1 = 0, c0 = 0, c1 = 0;
            for (int l = 0; l < 64; l++) {
                int i = 64 * g + l, pa = 2 * i * eta + j, pb = 2 * i * eta + eta + j;
                a0 |= (uint64_t)((b0[pa / 8] >> (pa % 8)) & 1) << l; a1 |= (uint64_t)((b1[pa / 8] >> (pa % 8)) & 1) << l;
                c0 |= (uint64_t)((b0[pb / 8] >> (pb % 8)) & 1) << l; c1 |= (uint64_t)((b1[pb / 8] >> (pb % 8)) & 1) << l;
            }
            A[j].s0 = a0; A[j].s1 = a1; B[j].s0 = c0; B[j].s1 = c1;
        }
        mw S[4], T[4], t[4], nT[4], V[4];
        for (int j = 0; j < 4; j++) { S[j] = mw_pub(0); T[j] = mw_pub(0); }
        for (int j = 0; j < eta; j++) {                      /* popcounts with bitsliced adders */
            t[0] = A[j]; t[1] = t[2] = t[3] = mw_pub(0); bs_add(S, S, t, 4);
            t[0] = B[j]; bs_add(T, T, t, 4);
        }
        bs_neg(nT, T, 4); bs_add(V, S, nT, 4);               /* S - T in [-4, 4], 4-bit two's complement */
        b2a_modq_bs(f->s0.c + 64 * g, f->s1.c + 64 * g, V, 4, 1, 64);
    }
}

/* masked E8 encoding of a Boolean-masked message -> arithmetic shares */
static void encode_masked(mpoly *mu, const uint8_t m0[32], const uint8_t m1[32]) {
    mw bit[9];
    for (int i = 0; i < 8; i++) {
        uint64_t x0 = 0, x1 = 0;
        for (int b = 0; b < 32; b++) { x0 |= (uint64_t)((m0[b] >> i) & 1) << b; x1 |= (uint64_t)((m1[b] >> i) & 1) << b; }
        bit[i].s0 = x0; bit[i].s1 = x1;
    }
    mw coset = bit[7];
    bit[7] = mw_pub(0); for (int i = 0; i < 7; i++) bit[7] = mx(bit[7], bit[i]);   /* x7 = parity of x0..x6 */
    int16_t c0[32], c1[32], x0[32], x1[32];
    b2a_modq_bs(c0, c1, &coset, 1, 0, 32);
    for (int i = 0; i < 8; i++) {
        b2a_modq_bs(x0, x1, &bit[i], 1, 0, 32);
        for (int b = 0; b < 32; b++) {
            mu->s0.c[E8_LAYOUT[8 * b + i]] = (int16_t)((1665 * x0[b] + 832 * c0[b]) % Q);
            mu->s1.c[E8_LAYOUT[8 * b + i]] = (int16_t)((1665 * x1[b] + 832 * c1[b]) % Q);
        }
    }
}

/* public decompression intervals: compress_d(x) == c  <=>  (x - LO[c]) mod q < WID[c] */
static int16_t LO[12][2048], WID[12][2048]; static int iv_done[12];
static void intervals(int d) {
    if (iv_done[d]) return;
    int n = 1 << d; int first[2048], cnt[2048];
    for (int c = 0; c < n; c++) { first[c] = -1; cnt[c] = 0; }
    for (int x = 0; x < Q; x++) { int c = compress((int16_t)x, d); cnt[c]++; }
    for (int c = 0; c < n; c++) {                            /* start = x with compress(x)=c and compress(x-1)!=c */
        for (int x = 0; x < Q; x++) if (compress((int16_t)x, d) == c && compress((int16_t)((x + Q - 1) % Q), d) != c) { first[c] = x; break; }
        LO[d][c] = (int16_t)(first[c] < 0 ? 0 : first[c]); WID[d][c] = (int16_t)cnt[c];
    }
    iv_done[d] = 1;
}
static mw check_poly(const mpoly *p, const poly *cc, int d, mw acc) {
    intervals(d);
    for (int g = 0; g < 4; g++) {
        uint32_t v0[64], v1[64], off[64], wid[64];
        for (int l = 0; l < 64; l++) {
            v0[l] = (uint32_t)p->s0.c[64 * g + l]; v1[l] = (uint32_t)p->s1.c[64 * g + l];
            int c = cc->c[64 * g + l];
            off[l] = (uint32_t)((Q - LO[d][c]) % Q); wid[l] = (uint32_t)(-WID[d][c]) & 0x3FFF;
        }
        mw x[14], s[14], sq[14], xp[14], df[14]; uint64_t ko[14], kw[14];
        a2b_modq_bs(x, v0, v1, 64, 14);
        bs_transpose_in(ko, off, 64, 14); bs_transpose_in(kw, wid, 64, 14);
        bs_addpub(s, x, ko, 14);                             /* x + (q - lo) in [0, 2q) */
        bs_addconst(sq, s, (uint32_t)(-Q) & 0x3FFF, 14);
        bs_sel(xp, sq[13], s, sq, 14);                       /* (x - lo) mod q */
        bs_addpub(df, xp, kw, 14);                           /* (x - lo) - width < 0  <=> inside */
        acc = secand(acc, df[13]);
    }
    return acc;
}

int mlkem_decaps_masked(int id, const uint8_t *dk, const uint8_t *c, uint8_t K[32]) {
    const mlkem_params *P = &PARAMS[id]; init_zetas();
    if (!P->e8) return -1;
    int k = P->k; size_t ctlen = mlkem_ct_bytes(P);
    const uint8_t *ek = dk + 384 * k, *h = dk + 768 * k + 32, *z = dk + 768 * k + 64;
    /* ---- load and mask the secret key (prototype) ---- */
    mpoly s[4]; uint8_t z0[32], z1[32];
    for (int i = 0; i < k; i++) {
        poly t; byte_decode(&t, dk + 384 * i, 12);
        for (int n = 0; n < N; n++) { int16_t r = (int16_t)rnd_modq(); s[i].s0.c[n] = r; s[i].s1.c[n] = (int16_t)((t.c[n] - r + Q) % Q); }
    }
    for (int i = 0; i < 32; i++) { uint8_t r = (uint8_t)mk_rnd(); z0[i] = z[i] ^ r; z1[i] = r; }
    /* ---- decrypt (linear, share-wise) ---- */
    poly u[4], v, cu[4], cv; int off = 0;
    for (int i = 0; i < k; i++) {
        byte_decode(&cu[i], c + off, P->du); off += 32 * P->du;
        for (int n = 0; n < N; n++) u[i].c[n] = decompress(cu[i].c[n], P->du);
        ntt(&u[i]);
    }
    byte_decode(&cv, c + off, P->dv);
    for (int n = 0; n < N; n++) v.c[n] = decompress(cv.c[n], P->dv);
    uint16_t w0[256], w1[256];
    for (int sh = 0; sh < 2; sh++) {
        poly acc; memset(&acc, 0, sizeof acc);
        for (int i = 0; i < k; i++) basemul_acc(&acc, sh ? &s[i].s1 : &s[i].s0, &u[i]);
        poly_reduce(&acc); invntt_tomont(&acc);
        for (int n = 0; n < N; n++) {
            int16_t x = sh ? (int16_t)(-acc.c[n]) : (int16_t)(v.c[n] - acc.c[n]);
            (sh ? w1 : w0)[n] = (uint16_t)canon(x);
        }
    }
    /* ---- masked decode, G, re-encryption ---- */
    uint8_t m0[32], m1[32], g0[64], g1[64];
    if (mk_leak_on) for (int n = 0; n < N; n++) { mk_leak_word(w0[n]); mk_leak_word(w1[n]); }
    { uint16_t t0[256], t1[256]; for (int j = 0; j < N; j++) { t0[j] = w0[E8_LAYOUT[j]]; t1[j] = w1[E8_LAYOUT[j]]; }
      e8_decode_masked_bs(t0, t1, m0, m1); }   /* v13 layout */
    msha3_512(g0, g1, m0, m1, 32, h, 32);
    poly A[4][4], t[4];
    const uint8_t *rho = ek + 384 * k;
    for (int i = 0; i < k; i++) byte_decode(&t[i], ek + 384 * i, 12);
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) sample_ntt(&A[i][j], rho, (uint8_t)j, (uint8_t)i);
    mpoly y[4], e1[4], e2, mu; uint8_t nonce = 0, b0[256], b1[256];
    for (int i = 0; i < k; i++) { mshake256(b0, b1, 64 * P->eta1, g0 + 32, g1 + 32, 32, &nonce, 1); nonce++; cbd_masked(&y[i], b0, b1, P->eta1); }
    for (int i = 0; i < k; i++) { mshake256(b0, b1, 64 * P->eta2, g0 + 32, g1 + 32, 32, &nonce, 1); nonce++; cbd_masked(&e1[i], b0, b1, P->eta2); }
    mshake256(b0, b1, 64 * P->eta2, g0 + 32, g1 + 32, 32, &nonce, 1); cbd_masked(&e2, b0, b1, P->eta2);
    encode_masked(&mu, m0, m1);
    mpoly up[4], vp;
    for (int sh = 0; sh < 2; sh++) {
        poly yh[4];
        for (int j = 0; j < k; j++) { yh[j] = sh ? y[j].s1 : y[j].s0; ntt(&yh[j]); }
        for (int i = 0; i < k; i++) {
            poly acc; memset(&acc, 0, sizeof acc);
            for (int j = 0; j < k; j++) basemul_acc(&acc, &A[j][i], &yh[j]);
            poly_reduce(&acc); invntt_tomont(&acc);
            poly *e = sh ? &e1[i].s1 : &e1[i].s0, *o = sh ? &up[i].s1 : &up[i].s0;
            for (int n = 0; n < N; n++) o->c[n] = canon((int16_t)(acc.c[n] + e->c[n]));
        }
        poly acc; memset(&acc, 0, sizeof acc);
        for (int j = 0; j < k; j++) basemul_acc(&acc, &t[j], &yh[j]);
        poly_reduce(&acc); invntt_tomont(&acc);
        poly *ee = sh ? &e2.s1 : &e2.s0, *mm = sh ? &mu.s1 : &mu.s0, *o = sh ? &vp.s1 : &vp.s0;
        for (int n = 0; n < N; n++) o->c[n] = canon((int16_t)(barrett((int16_t)(acc.c[n] + ee->c[n])) + mm->c[n]));
    }
    if (mk_leak_on) {
        for (int i = 0; i < k; i++) for (int n = 0; n < N; n++) { mk_leak_word((uint16_t)up[i].s0.c[n]); mk_leak_word((uint16_t)up[i].s1.c[n]); }
        for (int n = 0; n < N; n++) { mk_leak_word((uint16_t)vp.s0.c[n]); mk_leak_word((uint16_t)vp.s1.c[n]); }
    }
    /* ---- masked comparison with the received ciphertext ---- */
    mw acc = mw_pub(~0ULL);
    for (int i = 0; i < k; i++) acc = check_poly(&up[i], &cu[i], P->du, acc);
    acc = check_poly(&vp, &cv, P->dv, acc);
    for (int sh = 32; sh >= 1; sh >>= 1) acc = secand(acc, mshr(acc, sh));
    int valid = (int)(munmask(acc) & 1);
    /* ---- implicit rejection key, masked ---- */
    uint8_t kb0[32], kb1[32];
    mshake256(kb0, kb1, 32, z0, z1, 32, c, (int)ctlen);
    /* select on SHARES first, then unmask only the selected key (never unmask K' of a
       rejected ciphertext: that would be a plaintext-checking oracle) */
    uint8_t msk = (uint8_t)(-valid), k0[32], k1[32];
    for (int i = 0; i < 32; i++) {
        k0[i] = (uint8_t)((g0[i] & msk) | (kb0[i] & (uint8_t)~msk));
        k1[i] = (uint8_t)((g1[i] & msk) | (kb1[i] & (uint8_t)~msk));
    }
    for (int i = 0; i < 32; i++) K[i] = k0[i] ^ k1[i];
    return valid;
}
