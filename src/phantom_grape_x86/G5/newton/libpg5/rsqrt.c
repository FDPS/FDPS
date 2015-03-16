#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#define LEN 1000000
#define SCALE (1)

//static double rsqrt_correct = 1.0;

#if 1
double rsqrt(double x){
  __m128 xfv, xv;
  float xf;

  xf = (float)x;
  xfv = _mm_load_ss(&xf);
  xv = _mm_rsqrt_ss(xfv);
  _mm_store_ss(&xf, xv);
  //asm("rsqrtss %1, %0":"=x"(xf):"x"(xf));
  return (double)xf; //  * rsqrt_correct;
}
#else
static double rsqrt_correct = 1.0;
double rsqrt(double x){
  float xf = x;
  float yf;
  asm("rsqrtss %1, %0":"=x"(yf):"x"(xf));
  xf *= yf;
  xf *= yf;
  xf -= 3.0f;
  xf *= yf;
  xf *= -0.5;
  return (double)xf  * rsqrt_correct;
}
#endif

double drsqrt(double x){
  return 1./sqrt(x);
}

double rsqrt_bias(){
  int i;
  union {int i; float f;} buf[2];
  /*
    int array[2];
    float volatile *xfp = (float *)array;
    int volatile *varray = array;
  */
  double x0, x1, y0, y1;
  double weight;
  double err0, err1;
  double SumWeight, SumErr, SumSqrErr;
  double bias, mse;
  
  SumWeight = SumErr = SumSqrErr = 0.;
  for(i=0;i<(1<<23);i++){

    buf[0].i = (127<<23) + i;
    buf[1].i = (128<<23) + i;

    x0 = (double)buf[0].f;
    x1 = (double)buf[1].f;
    y0 = rsqrt(x0);
    y1 = rsqrt(x1);
    err0 = y0*y0*x0 - 1.;
    err1 = y1*y1*x1 - 1.;
    // printf("%f %f %e, %f %f %e\n", x0, y0, err0, x1, y1, err1);

    weight = y0*y0; /* 1/x */
    SumWeight += weight;
    SumErr += weight*(err0+err1);
    SumSqrErr += weight*(err0*err0 + err1*err1);
  }
  SumWeight *= 2.0;
  bias = 0.5 * SumErr / SumWeight;
  mse  = 0.5 * sqrt(SumSqrErr / SumWeight);
  fprintf(stderr, "rsqrt: MSE = %e,  Bias = %e\n", mse, bias);
	
  return bias;
}

#if 0
int main(){
  /*
    int i, exp;
    double frac;
    double x, y;
    double delta;
    double sum = 0.0;
	  
    for(i=0;i<LEN;i++){
    x = (i+1.) * 0.11;
    y = rsqrt(x);
    delta = x*y*y - 1.0;
    sum += delta;
    frac = frexp(x, &exp);
    if(exp&1) frac*=2;
    printf("%f %f %e %e\n", x, frac, delta*SCALE, sum);
    }
  */
  // fprintf(stderr, "n: %d, sum: %e\n", LEN, sum);
  rsqrt_correct = 1.0- rsqrt_bias();
  rsqrt_bias();
  
  return 0;
}
#endif
