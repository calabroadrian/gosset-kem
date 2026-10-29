#include <stdio.h>
#include <x86intrin.h>
#include "mlkem.c"
static inline uint64_t cyc(void){unsigned a;_mm_lfence();uint64_t t=__rdtscp(&a);_mm_lfence();return t;}
#define BENCH(name, expr) { uint64_t best=~0ULL; for(int r=0;r<2000;r++){uint64_t a=cyc(); expr; uint64_t b=cyc(); if(b-a<best)best=b-a;} printf("  %-28s %8llu cycles\n", name, (unsigned long long)best); }
int main(void){
  init_zetas(); poly p, q2, r; uint8_t seed[32]={1}, buf[1600]={0}, out[64]; uint64_t st[25]={0};
  for(int i=0;i<256;i++){p.c[i]=i*7%3329;q2.c[i]=i*13%3329;}
  BENCH("keccakf (1 permutation)", keccakf(st));
  BENCH("sample_ntt (1 poly of A)", sample_ntt(&r, seed, 0, 1));
  BENCH("cbd eta=2 (1 poly)", cbd(&r, seed, 3, 2));
  BENCH("ntt", ntt(&p));
  BENCH("invntt", invntt(&p));
  BENCH("basemul_acc", basemul_acc(&r,&p,&q2));
  BENCH("byte_encode 12", byte_encode(buf,&p,12));
  BENCH("sha3_512(64B)", sha3_512(out, buf, 64));
  BENCH("shake256 J(800B)", shake256(out,32,buf,800));
  BENCH("e8_decode_avx2", e8_decode_avx2(p.c, out));
  return 0;}
