#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <x86intrin.h>
void e8_decode_masked(const uint16_t w0[256], const uint16_t w1[256], uint8_t m0[32], uint8_t m1[32], int shuffle);
void e8_decode_masked_bs(const uint16_t w0[256], const uint16_t w1[256], uint8_t m0[32], uint8_t m1[32]);
void e8_decode_avx2(const int16_t in[256], uint8_t msg[32]);
static inline uint64_t cyc(void){unsigned a;_mm_lfence();uint64_t t=__rdtscp(&a);_mm_lfence();return t;}
static int cmpu(const void*a,const void*b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return (x>y)-(x<y);}
int main(void){ uint16_t w0[256],w1[256]; int16_t w[256]; uint8_t a[32],b[32]; static uint64_t t[3][3001];
  for(int r=0;r<3001;r++){ for(int i=0;i<256;i++){w[i]=(int16_t)(rand()%3329); w0[i]=(uint16_t)(rand()%3329); w1[i]=(uint16_t)((w[i]-w0[i]+3329)%3329);}
    for(int v=0;v<3;v++){ int k=(r+v)%3; uint64_t c0=cyc();
      if(k==0) e8_decode_avx2(w,a); else if(k==1) e8_decode_masked(w0,w1,a,b,0); else e8_decode_masked_bs(w0,w1,a,b);
      t[k][r]=cyc()-c0; } }
  const char*n[3]={"E8 AVX2 (unmasked)","E8 masked v10 (per block)","E8 masked v11 (bitsliced)"};
  for(int k=0;k<3;k++){qsort(t[k],3001,8,cmpu); printf("  %-30s %9llu cycles\n",n[k],(unsigned long long)t[k][1500]);}
  printf("  speed-up bitsliced vs v10: %.1fx\n",(double)t[1][1500]/t[2][1500]); return 0; }
