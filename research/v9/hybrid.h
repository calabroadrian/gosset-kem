#ifndef HYBRID_H
#define HYBRID_H
#include <stdint.h>
#include <stddef.h>
void secure_zero(void *p, size_t n);
size_t hyb_pk_bytes(int id);
size_t hyb_sk_bytes(int id);
size_t hyb_ct_bytes(int id);
void hyb_keygen_internal(int id, const uint8_t d[32], const uint8_t z[32], const uint8_t skx[32], uint8_t *pk, uint8_t *sk);
void hyb_encaps_internal(int id, const uint8_t *pk, const uint8_t m[32], const uint8_t ex[32], uint8_t *ct, uint8_t K[32]);
void hyb_decaps(int id, const uint8_t *sk, const uint8_t *ct, uint8_t K[32]);
int hyb_keygen(int id, uint8_t *pk, uint8_t *sk);
int hyb_encaps(int id, const uint8_t *pk, uint8_t *ct, uint8_t K[32]);
int hyb_decaps_checked(int id, const uint8_t *sk, const uint8_t *ct, uint8_t K[32]);
#endif
