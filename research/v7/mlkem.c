/*
 * AFA-12 v7 — ML-KEM (FIPS 203) reference-style implementation + E8 message-coding variant
 * Own implementation written from the FIPS 203 algorithms (KeyGen/Encaps/Decaps internal).
 * Both variants share every line except the message encode/decode and du/dv:
 *   standard : Decompress_1(ByteDecode_1(m)) / ByteEncode_1(Compress_1(w))   (FIPS 203)
 *   E8       : e8_encode(m) / e8_decode_avx2(w), and smaller du/dv (compact mode)
 * Arithmetic: canonical representatives in [0,q), reductions by constant modulus
 * (compiled to multiply-shift, no division instructions). Not optimized like the
 * official AVX2 code: absolute cycle counts are higher, relative cost is what we compare.
 */
#include <stdint.h>
#include <string.h>
#include "fips202.h"
#include "e8_layout.h"
#include "mlkem.h"

#define N 256
#define Q 3329

void e8_encode(const uint8_t msg[32], int16_t out[256]);
void e8_decode_avx2(const int16_t in[256], uint8_t msg[32]);

static const mlkem_params PARAMS[9] = {
    {2, 3, 2, 10, 4, 0, "ML-KEM-512"},  {3, 2, 2, 10, 4, 0, "ML-KEM-768"},  {4, 2, 2, 11, 5, 0, "ML-KEM-1024"},
    {2, 3, 2, 9, 4, 1, "ML-KEM-512+E8"}, {3, 2, 2, 9, 5, 1, "ML-KEM-768+E8"}, {4, 2, 2, 10, 4, 1, "ML-KEM-1024+E8"},
    {2, 4, 3, 10, 4, 1, "AFA12-512-MaxSec"}, {3, 3, 2, 10, 4, 1, "AFA12-768-MaxSec"}, {4, 3, 2, 11, 5, 1, "AFA12-1024-MaxSec"}};
const mlkem_params *mlkem_get(int id) { return &PARAMS[id]; }

typedef struct { int16_t c[N]; } poly;
static int16_t ZETAS[128];

static inline int16_t red(int32_t x) { int32_t r = x % Q; r += Q & (r >> 31); return (int16_t)r; }

static int bitrev7(int i) { int r = 0; for (int b = 0; b < 7; b++) r |= ((i >> b) & 1) << (6 - b); return r; }
static void init_zetas(void) {
    static int done = 0; if (done) return;
    for (int i = 0; i < 128; i++) {
        int32_t z = 1; for (int e = 0; e < bitrev7(i); e++) z = (z * 17) % Q;
        ZETAS[i] = (int16_t)z;
    }
    done = 1;
}

static void ntt(poly *f) {
    int i = 1;
    for (int len = 128; len >= 2; len >>= 1)
        for (int start = 0; start < N; start += 2 * len) {
            int32_t z = ZETAS[i++];
            for (int j = start; j < start + len; j++) {
                int16_t t = red(z * f->c[j + len]);
                f->c[j + len] = red(f->c[j] - t);
                f->c[j] = red(f->c[j] + t);
            }
        }
}

static void invntt(poly *f) {
    int i = 127;
    for (int len = 2; len <= 128; len <<= 1)
        for (int start = 0; start < N; start += 2 * len) {
            int32_t z = ZETAS[i--];
            for (int j = start; j < start + len; j++) {
                int16_t t = f->c[j];
                f->c[j] = red(t + f->c[j + len]);
                f->c[j + len] = red(z * (f->c[j + len] - t));
            }
        }
    for (int j = 0; j < N; j++) f->c[j] = red(3303 * (int32_t)f->c[j]);
}

static void basemul_acc(poly *r, const poly *a, const poly *b) {    /* r += a o b in NTT domain */
    for (int i = 0; i < 64; i++) {
        int32_t g = ZETAS[64 + i];
        for (int s = 0; s < 2; s++) {
            int32_t gg = s ? Q - g : g;
            int k = 4 * i + 2 * s;
            int32_t a0 = a->c[k], a1 = a->c[k + 1], b0 = b->c[k], b1 = b->c[k + 1];
            int32_t c0 = red(a0 * b0 + red(a1 * b1) * gg);
            int32_t c1 = red(a0 * b1 + a1 * b0);
            r->c[k] = red(r->c[k] + c0);
            r->c[k + 1] = red(r->c[k + 1] + c1);
        }
    }
}

static void sample_ntt(poly *a, const uint8_t rho[32], uint8_t j, uint8_t i) {
    uint8_t seed[34]; memcpy(seed, rho, 32); seed[32] = j; seed[33] = i;
    keccak_ctx c; keccak_init(&c, 168, 0x1F); keccak_absorb(&c, seed, 34); keccak_finalize(&c);
    int n = 0; uint8_t b[168];
    while (n < N) {
        keccak_squeeze(&c, b, 168);
        for (int p = 0; p + 3 <= 168 && n < N; p += 3) {
            int d1 = b[p] + 256 * (b[p + 1] & 15), d2 = (b[p + 1] >> 4) + 16 * b[p + 2];
            if (d1 < Q) a->c[n++] = (int16_t)d1;
            if (d2 < Q && n < N) a->c[n++] = (int16_t)d2;
        }
    }
}

static void cbd(poly *f, const uint8_t sigma[32], uint8_t nonce, int eta) {
    uint8_t in[33], buf[64 * 4]; memcpy(in, sigma, 32); in[32] = nonce;
    shake256(buf, 64 * eta, in, 33);
    for (int i = 0; i < N; i++) {
        int x = 0, y = 0;
        for (int j = 0; j < eta; j++) {
            int bx = 2 * i * eta + j, by = 2 * i * eta + eta + j;
            x += (buf[bx / 8] >> (bx % 8)) & 1;
            y += (buf[by / 8] >> (by % 8)) & 1;
        }
        f->c[i] = red(x - y);
    }
}

/* word-based bit packing (as in optimized implementations) */
static void byte_encode(uint8_t *out, const poly *f, int d) {
    uint64_t acc = 0; int bits = 0, pos = 0;
    for (int i = 0; i < N; i++) {
        acc |= (uint64_t)(uint16_t)f->c[i] << bits; bits += d;
        while (bits >= 8) { out[pos++] = (uint8_t)acc; acc >>= 8; bits -= 8; }
    }
}
static void byte_decode(poly *f, const uint8_t *in, int d) {
    uint64_t acc = 0; int bits = 0, pos = 0; uint32_t mask = (1u << d) - 1;
    for (int i = 0; i < N; i++) {
        while (bits < d) { acc |= (uint64_t)in[pos++] << bits; bits += 8; }
        int32_t v = (int32_t)(acc & mask); acc >>= d; bits -= d;
        f->c[i] = (int16_t)(d == 12 ? v % Q : v);
    }
}
static inline int16_t compress(int16_t x, int d) { return (int16_t)((((uint32_t)x << d) + Q / 2) / Q & ((1u << d) - 1)); }
static inline int16_t decompress(int16_t y, int d) { return (int16_t)(((uint32_t)y * Q + (1u << (d - 1))) >> d); }

/* ---------------- K-PKE ---------------- */
static void pke_keygen(const mlkem_params *P, const uint8_t d[32], uint8_t *ek, uint8_t *dk) {
    int k = P->k; uint8_t in[33], g[64]; memcpy(in, d, 32); in[32] = (uint8_t)k;
    sha3_512(g, in, 33);
    const uint8_t *rho = g, *sigma = g + 32; uint8_t nonce = 0;
    poly A[4][4], s[4], e[4], t[4];
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) sample_ntt(&A[i][j], rho, (uint8_t)j, (uint8_t)i);
    for (int i = 0; i < k; i++) cbd(&s[i], sigma, nonce++, P->eta1);
    for (int i = 0; i < k; i++) cbd(&e[i], sigma, nonce++, P->eta1);
    for (int i = 0; i < k; i++) { ntt(&s[i]); ntt(&e[i]); }
    for (int i = 0; i < k; i++) {
        t[i] = e[i];
        for (int j = 0; j < k; j++) basemul_acc(&t[i], &A[i][j], &s[j]);
    }
    for (int i = 0; i < k; i++) byte_encode(ek + 384 * i, &t[i], 12);
    memcpy(ek + 384 * k, rho, 32);
    for (int i = 0; i < k; i++) byte_encode(dk + 384 * i, &s[i], 12);
}

static void pke_encrypt(const mlkem_params *P, const uint8_t *ek, const uint8_t m[32], const uint8_t r[32], uint8_t *c) {
    int k = P->k; uint8_t nonce = 0; const uint8_t *rho = ek + 384 * k;
    poly A[4][4], t[4], y[4], e1[4], e2, u[4], v, mu;
    for (int i = 0; i < k; i++) byte_decode(&t[i], ek + 384 * i, 12);
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) sample_ntt(&A[i][j], rho, (uint8_t)j, (uint8_t)i);
    for (int i = 0; i < k; i++) cbd(&y[i], r, nonce++, P->eta1);
    for (int i = 0; i < k; i++) cbd(&e1[i], r, nonce++, P->eta2);
    cbd(&e2, r, nonce, P->eta2);
    for (int i = 0; i < k; i++) ntt(&y[i]);
    for (int i = 0; i < k; i++) {
        memset(&u[i], 0, sizeof(poly));
        for (int j = 0; j < k; j++) basemul_acc(&u[i], &A[j][i], &y[j]);
        invntt(&u[i]);
        for (int n = 0; n < N; n++) u[i].c[n] = red(u[i].c[n] + e1[i].c[n]);
    }
    memset(&v, 0, sizeof v);
    for (int j = 0; j < k; j++) basemul_acc(&v, &t[j], &y[j]);
    invntt(&v);
    if (P->e8) { int16_t tq[N]; e8_encode(m, tq); for (int j = 0; j < N; j++) mu.c[E8_LAYOUT[j]] = tq[j]; }   /* v13 layout */
    else for (int n = 0; n < N; n++) mu.c[n] = (int16_t)(((m[n / 8] >> (n % 8)) & 1) * 1665);
    for (int n = 0; n < N; n++) v.c[n] = red(v.c[n] + e2.c[n] + mu.c[n]);
    int off = 0; poly tmp;
    for (int i = 0; i < k; i++) {
        for (int n = 0; n < N; n++) tmp.c[n] = compress(u[i].c[n], P->du);
        byte_encode(c + off, &tmp, P->du); off += 32 * P->du;
    }
    for (int n = 0; n < N; n++) tmp.c[n] = compress(v.c[n], P->dv);
    byte_encode(c + off, &tmp, P->dv);
}

static void pke_decrypt(const mlkem_params *P, const uint8_t *dk, const uint8_t *c, uint8_t m[32]) {
    int k = P->k; poly u[4], s[4], v, w;
    int off = 0;
    for (int i = 0; i < k; i++) {
        byte_decode(&u[i], c + off, P->du); off += 32 * P->du;
        for (int n = 0; n < N; n++) u[i].c[n] = decompress(u[i].c[n], P->du);
        ntt(&u[i]);
        byte_decode(&s[i], dk + 384 * i, 12);
    }
    byte_decode(&v, c + off, P->dv);
    for (int n = 0; n < N; n++) v.c[n] = decompress(v.c[n], P->dv);
    memset(&w, 0, sizeof w);
    for (int i = 0; i < k; i++) basemul_acc(&w, &s[i], &u[i]);
    invntt(&w);
    for (int n = 0; n < N; n++) w.c[n] = red(v.c[n] - w.c[n]);
    if (P->e8) { int16_t tq[N]; for (int j = 0; j < N; j++) tq[j] = w.c[E8_LAYOUT[j]]; e8_decode_avx2(tq, m); }   /* v13 layout */
    else { memset(m, 0, 32); for (int n = 0; n < N; n++) m[n / 8] |= (uint8_t)(compress(w.c[n], 1) << (n % 8)); }
}

/* ---------------- ML-KEM ---------------- */
size_t mlkem_ek_bytes(const mlkem_params *P) { return 384 * P->k + 32; }
size_t mlkem_dk_bytes(const mlkem_params *P) { return 768 * P->k + 96; }
size_t mlkem_ct_bytes(const mlkem_params *P) { return 32 * (P->du * P->k + P->dv); }

void mlkem_keygen_internal(int id, const uint8_t d[32], const uint8_t z[32], uint8_t *ek, uint8_t *dk) {
    const mlkem_params *P = &PARAMS[id]; init_zetas();
    int k = P->k;
    pke_keygen(P, d, ek, dk);
    memcpy(dk + 384 * k, ek, 384 * k + 32);
    sha3_256(dk + 768 * k + 32, ek, 384 * k + 32);
    memcpy(dk + 768 * k + 64, z, 32);
}

void mlkem_encaps_internal(int id, const uint8_t *ek, const uint8_t m[32], uint8_t K[32], uint8_t *c) {
    const mlkem_params *P = &PARAMS[id]; init_zetas();
    uint8_t in[64], g[64];
    memcpy(in, m, 32); sha3_256(in + 32, ek, mlkem_ek_bytes(P));
    sha3_512(g, in, 64);
    memcpy(K, g, 32);
    pke_encrypt(P, ek, m, g + 32, c);
}

void mlkem_decaps(int id, const uint8_t *dk, const uint8_t *c, uint8_t K[32]) {
    const mlkem_params *P = &PARAMS[id]; init_zetas();
    int k = P->k; size_t ctlen = mlkem_ct_bytes(P);
    const uint8_t *ek = dk + 384 * k, *h = dk + 768 * k + 32, *z = dk + 768 * k + 64;
    uint8_t m[32], in[64], g[64], cp[1600], kbar[32], jin[32 + 1600];
    pke_decrypt(P, dk, c, m);
    memcpy(in, m, 32); memcpy(in + 32, h, 32);
    sha3_512(g, in, 64);
    memcpy(jin, z, 32); memcpy(jin + 32, c, ctlen);
    shake256(kbar, 32, jin, 32 + ctlen);
    pke_encrypt(P, ek, m, g + 32, cp);
    uint8_t diff = 0;
    for (size_t i = 0; i < ctlen; i++) diff |= (uint8_t)(c[i] ^ cp[i]);
    uint8_t mask = (uint8_t)(-(int)((-(uint32_t)diff) >> 31));   /* 0xFF if different */
    __asm__ volatile("" : "+r"(mask));
    for (int i = 0; i < 32; i++) K[i] = (uint8_t)(g[i] ^ (mask & (g[i] ^ kbar[i])));
}

/* failure-test hook: decrypt only, returns recovered message */
void mlkem_pke_decrypt(int id, const uint8_t *dk, const uint8_t *c, uint8_t m[32]) {
    init_zetas(); pke_decrypt(&PARAMS[id], dk, c, m);
}

/* ---------------- FIPS 203 input checks (Sec. 7.2 / 7.3) ---------------- */
/* Encapsulation key check: every 12-bit coefficient of ek must be < q
   (equivalently ByteEncode12(ByteDecode12(ek)) == ek). Returns 1 if valid. */
int mlkem_check_ek(int id, const uint8_t *ek) {
    const mlkem_params *P = &PARAMS[id]; int bad = 0;
    for (int i = 0; i < 128 * P->k; i++) {
        uint32_t a = ek[3 * i] | ((uint32_t)(ek[3 * i + 1] & 15) << 8);
        uint32_t b = (ek[3 * i + 1] >> 4) | ((uint32_t)ek[3 * i + 2] << 4);
        bad |= (a >= Q) | (b >= Q);
    }
    return !bad;
}
/* Decapsulation key check: H(ek embedded in dk) == stored hash. Returns 1 if valid. */
int mlkem_check_dk(int id, const uint8_t *dk) {
    const mlkem_params *P = &PARAMS[id]; int k = P->k; uint8_t h[32];
    sha3_256(h, dk + 384 * k, 384 * k + 32);
    uint8_t d = 0; for (int i = 0; i < 32; i++) d |= (uint8_t)(h[i] ^ dk[768 * k + 32 + i]);
    return d == 0;
}
