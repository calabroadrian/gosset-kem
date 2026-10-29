#ifndef MLKEM_H
#define MLKEM_H
#include <stdint.h>
#include <stddef.h>
typedef struct { int k, eta1, eta2, du, dv, e8; const char *name; } mlkem_params;
const mlkem_params *mlkem_get(int id);   /* 0..2 standard, 3..5 E8 compact, 6..8 AFA-12 MaxSec (mlkem_opt.c) */
size_t mlkem_ek_bytes(const mlkem_params *P);
size_t mlkem_dk_bytes(const mlkem_params *P);
size_t mlkem_ct_bytes(const mlkem_params *P);
void mlkem_keygen_internal(int id, const uint8_t d[32], const uint8_t z[32], uint8_t *ek, uint8_t *dk);
void mlkem_encaps_internal(int id, const uint8_t *ek, const uint8_t m[32], uint8_t K[32], uint8_t *c);
void mlkem_decaps(int id, const uint8_t *dk, const uint8_t *c, uint8_t K[32]);
void mlkem_pke_decrypt(int id, const uint8_t *dk, const uint8_t *c, uint8_t m[32]);
int mlkem_check_ek(int id, const uint8_t *ek);
int mlkem_check_dk(int id, const uint8_t *dk);
#endif
