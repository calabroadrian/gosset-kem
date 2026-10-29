#ifndef DEM_H
#define DEM_H
#include <stdint.h>
#include <stddef.h>
void hmac_sha3_256(const uint8_t *key, size_t klen, const uint8_t *msg, size_t mlen, uint8_t out[32]);
size_t dem_padded_len(size_t mlen);
void dem_encrypt(const uint8_t K[32], const uint8_t nonce[16], const uint8_t *M, size_t mlen, uint8_t *C, uint8_t tag[32]);
long dem_decrypt(const uint8_t K[32], const uint8_t nonce[16], const uint8_t *C, size_t clen, const uint8_t tag[32], uint8_t *M);
#endif
