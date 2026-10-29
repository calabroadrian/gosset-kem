#ifndef E8CODEC_H
#define E8CODEC_H
#include <stdint.h>
void e8_encode(const uint8_t msg[32], int16_t out[256]);
void e8_decode(const int16_t in[256], uint8_t msg[32]);
void std_decode(const int16_t in[256], uint8_t msg[32]);
void e8_decode_leaky(const int16_t in[256], uint8_t msg[32]);
int64_t e8_round_div_U_test(int64_t Y);
#endif
