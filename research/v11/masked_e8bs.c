/*
 * AFA-12 v11 — bitsliced first-order masked E8 decoder (32 blocks in parallel)
 * Same multiplication-free algorithm as v10 (verified identical to the reference
 * decoder), but every masked word carries 32 lanes (one per 8-coefficient block),
 * so each SecAND processes the whole message at once.
 */
#include <string.h>
#include "masked_core.h"

#define Q 3329
#define NBW 18
#define MASKW ((1u << NBW) - 1)

/* arithmetic shares mod q (v0 + v1 mod q, each in [0,q)) -> bitsliced Boolean shares, NBW bits */
void a2b_modq_bs(mw *x, const uint32_t *v0, const uint32_t *v1, int lanes, int nbw) {
    mw a[14], b[14], z[14], zq[14];
    bs_from_pub_lanes(a, v0, lanes, 14);                 /* each arithmetic share alone is uniform */
    bs_from_pub_lanes(b, v1, lanes, 14);
    bs_add(z, a, b, 14);                                 /* [0, 2q) */
    bs_addconst(zq, z, (uint32_t)(-Q) & 0x3FFF, 14);
    bs_sel(x, zq[13], z, zq, 14);                        /* z - q < 0 ? z : z - q */
    for (int j = 14; j < nbw; j++) x[j] = mw_pub(0);
}

void e8_decode_masked_bs(const uint16_t w0[256], const uint16_t w1[256], uint8_t m0[32], uint8_t m1[32]) {
    mw a[8][NBW], sumA[NBW], tb[8], sb[8], parA = {0, 0}, parB = {0, 0};
    for (int j = 0; j < NBW; j++) sumA[j] = mw_pub(0);
    for (int i = 0; i < 8; i++) {
        uint32_t v0[32], v1[32];
        for (int b = 0; b < 32; b++) { v0[b] = w0[8 * b + i]; v1[b] = w1[8 * b + i]; }
        mw w[NBW], c[NBW], wq[NBW], y[NBW], Y[NBW], t1[NBW], one[NBW], absY[NBW], nabs[NBW], tq[NBW], alt[NBW];
        a2b_modq_bs(w, v0, v1, 32, NBW);
        bs_addconst(c, w, (uint32_t)(-(Q / 2 + 1)) & MASKW, NBW);   /* w - 1665 */
        bs_addconst(wq, w, (uint32_t)(-Q) & MASKW, NBW);            /* w - q */
        bs_sel(y, c[NBW - 1], w, wq, NBW);                          /* centered */
        Y[0] = Y[1] = mw_pub(0); for (int j = 2; j < NBW; j++) Y[j] = y[j - 2];   /* Y = 4y */
        mw sm = Y[NBW - 1];
        for (int j = 0; j < NBW; j++) { t1[j] = mx(Y[j], sm); one[j] = mw_pub(0); }
        one[0] = sm;
        bs_add(absY, t1, one, NBW);                                 /* |Y| */
        bs_neg(nabs, absY, NBW);
        bs_addconst(tq, nabs, Q, NBW);                              /* q - |Y| */
        bs_addconst(alt, nabs, 2 * Q, NBW);                         /* 2q - |Y| */
        mw tm = tq[NBW - 1];                                        /* t = [|Y| > q] */
        bs_sel(a[i], tm, alt, absY, NBW);
        bs_add(sumA, sumA, a[i], NBW);
        tb[i] = tm; sb[i] = sm; parA = mx(parA, tm); parB = mx(parB, sm);
    }
    mw amax[NBW], amin[NBW], hA[8], hB[8], d[NBW], na[NBW];
    memcpy(amax, a[0], sizeof amax); memcpy(amin, a[0], sizeof amin);
    for (int k = 0; k < 8; k++) { hA[k] = mw_pub(k == 0 ? ~0ULL : 0); hB[k] = hA[k]; }
    for (int i = 1; i < 8; i++) {
        bs_neg(na, a[i], NBW); bs_add(d, amax, na, NBW); mw gt = d[NBW - 1];     /* a_i > amax */
        bs_sel(amax, gt, a[i], amax, NBW);
        for (int k = 0; k < 8; k++) hA[k] = secsel(gt, mw_pub(k == i ? ~0ULL : 0), hA[k]);
        bs_neg(na, amin, NBW); bs_add(d, a[i], na, NBW); mw lt = d[NBW - 1];     /* a_i < amin */
        bs_sel(amin, lt, a[i], amin, NBW);
        for (int k = 0; k < 8; k++) hB[k] = secsel(lt, mw_pub(k == i ? ~0ULL : 0), hB[k]);
    }
    mw am2[NBW], T1[NBW], P[NBW], D[NBW], nD[NBW];
    am2[0] = mw_pub(0); for (int j = 1; j < NBW; j++) am2[j] = amax[j - 1];
    bs_neg(T1, am2, NBW); bs_addconst(T1, T1, 2 * Q, NBW);                     /* 2q - 2 amax */
    bs_addconst(D, sumA, (uint32_t)(-4 * Q) & MASKW, NBW);
    for (int j = 0; j < NBW; j++) P[j] = secand(parA, T1[j]);
    bs_add(D, D, P, NBW);
    P[0] = mw_pub(0); for (int j = 1; j < NBW; j++) P[j] = secand(parB, amin[j - 1]);   /* parB * 2 amin */
    bs_neg(P, P, NBW); bs_add(D, D, P, NBW);
    bs_neg(nD, D, NBW);
    mw useB = nD[NBW - 1];                                                     /* D > 0 */
    mw out[8];
    for (int i = 0; i < 7; i++) {
        mw bA = mx(tb[i], secand(parA, hA[i])), bB = mx(sb[i], secand(parB, hB[i]));
        out[i] = secsel(useB, bB, bA);
    }
    out[7] = useB;
    for (int b = 0; b < 32; b++) {
        uint8_t x0 = 0, x1 = 0;
        for (int i = 0; i < 8; i++) { x0 |= (uint8_t)(((out[i].s0 >> b) & 1) << i); x1 |= (uint8_t)(((out[i].s1 >> b) & 1) << i); }
        m0[b] = x0; m1[b] = x1;
    }
}
