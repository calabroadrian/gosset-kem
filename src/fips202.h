#ifndef FIPS202_H
#define FIPS202_H
#include <stdint.h>
#include <stddef.h>
typedef struct { uint64_t st[25]; unsigned rate, pos; uint8_t pad; } keccak_ctx;
void keccakf(uint64_t st[25]);
void keccak_init(keccak_ctx *c, unsigned rate, uint8_t pad);
void keccak_absorb(keccak_ctx *c, const uint8_t *in, size_t len);
void keccak_finalize(keccak_ctx *c);
void keccak_squeeze(keccak_ctx *c, uint8_t *out, size_t len);
void sha3_256(uint8_t out[32], const uint8_t *in, size_t len);
void sha3_512(uint8_t out[64], const uint8_t *in, size_t len);
void shake128(uint8_t *out, size_t outlen, const uint8_t *in, size_t len);
void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t len);
#endif
