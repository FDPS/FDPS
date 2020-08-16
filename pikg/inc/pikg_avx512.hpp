#ifndef H_PIKG_AVX_512
#define H_PIKG_AVX_512
#include <immintrin.h>
struct __m512dx2{
  __m512d v0,v1;
};
struct __m512dx3{
  __m512d v0,v1,v2;
};
struct __m512dx4{
  __m512d v0,v1,v2,v3;
};
struct __m512x2{
  __m512 v0,v1;
};
struct __m512x3{
  __m512 v0,v1,v2;
};
struct __m512x4{
  __m512 v0,v1,v2,v3;
};
static inline __m512dx2 _mm512_set1_pdx2(const PIKG::F64vec2 v){
  __m512dx2 ret;
  ret.v0 = _mm512_set1_pd(v.x);
  ret.v1 = _mm512_set1_pd(v.y);
  return ret;
}
static inline __m512x2  _mm512_set1_psx2(const PIKG::F32vec2 v){
  __m512x2 ret;
  ret.v0 = _mm512_set1_ps(v.x);
  ret.v1 = _mm512_set1_ps(v.y);
  return ret;
}
static inline __m512dx3 _mm512_set1_pdx3(const PIKG::F64vec v){
  __m512dx3 ret;
  ret.v0 = _mm512_set1_pd(v.x);
  ret.v1 = _mm512_set1_pd(v.y);
  ret.v2 = _mm512_set1_pd(v.z);
  return ret;
}
static inline __m512x3  _mm512_set1_psx3(const PIKG::F32vec v){
  __m512x3 ret;
  ret.v0 = _mm512_set1_ps(v.x);
  ret.v1 = _mm512_set1_ps(v.y);
  ret.v2 = _mm512_set1_ps(v.z);
  return ret;
}
#endif
