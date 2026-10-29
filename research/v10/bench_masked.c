#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <x86intrin.h>
void e8_decode_masked(const uint16_t w0[256], const uint16_t w1[256], uint8_t m0[32], uint8_t m1[32], int shuffle);
void e8_decode_avx2(const int16_t in[256], uint8_t msg[32]);
void std_decode_ref(const int16_t in[256], uint8_t msg[32]);
static inline uint64_t cyc(void){unsigned a;_mm_lfence();uint64_t t=__rdtscp(&a);_mm_lfence();return t;}
static int cmpu(const void*a,const void*b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return (x>y)-(x<y);}
int main(void){
  uint16_t w0[256],w1[256]; int16_t w[256]; uint8_t a[32],b[32]; uint64_t t[4][2001];
  for(int r=0;r<2001;r++){
    for(int i=0;i<256;i++){w[i]=(int16_t)(rand()%3329); w0[i]=(uint16_t)(rand()%3329); w1[i]=(uint16_t)((w[i]-w0[i]+3329)%3329);}
    uint64_t c0=cyc(); std_decode_ref(w,a); uint64_t c1=cyc(); e8_decode_avx2(w,a); uint64_t c2=cyc();
    e8_decode_masked(w0,w1,a,b,0); uint64_t c3=cyc(); e8_decode_masked(w0,w1,a,b,1); uint64_t c4=cyc();
    t[0][r]=c1-c0; t[1][r]=c2-c1; t[2][r]=c3-c2; t[3][r]=c4-c3;
  }
  const char*n[4]={"ML-KEM poly_tomsg (unmasked)","E8 AVX2 (unmasked)","E8 masked, 1st order","E8 masked + shuffling"};
  for(int k=0;k<4;k++){qsort(t[k],2001,8,cmpu); printf("  %-30s %9llu cycles\n",n[k],(unsigned long long)t[k][1000]);}
  return 0;}
