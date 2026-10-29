#ifndef LEECH_H
#define LEECH_H
#include <stdint.h>
void leech_encode_block(uint64_t m, int16_t out[24]);
uint64_t leech_decode_block(const int16_t in[24]);
void leech_decode(const int16_t in[192], uint64_t out[8]);
#endif
