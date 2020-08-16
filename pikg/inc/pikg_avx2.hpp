#ifndef H_PIKG_AVX2
#define H_PIKG_AVX2
#include <immintrin.h>
struct __m256dx2{
  __m256d v0,v1;
};
struct __m256dx3{
  __m256d v0,v1,v2;
};
struct __m256dx4{
  __m256d v0,v1,v2,v3;
};
struct __m256x2{
  __m256 v0,v1;
};
struct __m256x3{
  __m256 v0,v1,v2;
};
struct __m256x4{
  __m256 v0,v1,v2,v3;
};
static inline __m256dx2 _mm256_set1_pdx2(const PIKG::F64vec2 v){
  __m256dx2 ret;
  ret.v0 = _mm256_set1_pd(v.x);
  ret.v1 = _mm256_set1_pd(v.y);
  return ret;
}
static inline __m256x2  _mm256_set1_psx2(const PIKG::F32vec2 v){
  __m256x2 ret;
  ret.v0 = _mm256_set1_ps(v.x);
  ret.v1 = _mm256_set1_ps(v.y);
  return ret;
}
static inline __m256dx3 _mm256_set1_pdx3(const PIKG::F64vec v){
  __m256dx3 ret;
  ret.v0 = _mm256_set1_pd(v.x);
  ret.v1 = _mm256_set1_pd(v.y);
  ret.v2 = _mm256_set1_pd(v.z);
  return ret;
}
static inline __m256x3  _mm256_set1_psx3(const PIKG::F32vec v){
  __m256x3 ret;
  ret.v0 = _mm256_set1_ps(v.x);
  ret.v1 = _mm256_set1_ps(v.y);
  ret.v2 = _mm256_set1_ps(v.z);
  return ret;
}
#endif
