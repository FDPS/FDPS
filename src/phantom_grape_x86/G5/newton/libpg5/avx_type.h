#ifndef __AVX_TYPE__
#define __AVX_TYPE__

typedef double v4df __attribute__((vector_size(32)));
typedef float  v8sf __attribute__((vector_size(32)));

#define ALIGN64 __attribute__ ((aligned(64)))
#define NUNROLL 4
//#define NUNROLL 2

typedef struct ipdata{
  float x[4];
  float y[4];
  float z[4];
  float eps2[4];
} Ipdata, *pIpdata;

typedef struct jpdata{
  float x, y, z, m;
} Jpdata, *pJpdata;

typedef struct jpdata0{
  float xm[2][4], ep[2][4];
} Jpdata0, *pJpdata0;

typedef struct fodata{
  float ax[4];
  float ay[4];
  float az[4];
  float phi[4];
} Fodata, *pFodata;
#endif /* __AVX_TYPE__ */
