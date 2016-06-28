#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifdef __FFTW2__
#include "sfftw.h"
#include "srfftw.h"
#elif __FFTW3__
#include "fftw3.h"
#define cmplx_re(c) ((c)[0])
#define cmplx_im(c) ((c)[1])
#else
#error "Neither __FFTW2__ nor __FFTW3__ defined"
#endif

#define TINY (1.0e-30)

#include "constants.h"
#include "run_param.h"

//#define __SAVITZKY_GOLAY_FILTERING__

#ifdef __SAVITZKY_GOLAY_FILTERING__
#define NFACT_PKBIN (5)
#else
#define NFACT_PKBIN (1)
#endif

#define NBIN_PK    (nmesh*NFACT_PKBIN)

#define MESH(ix, iy, iz) (mesh[(iz)+nmesh_p2*((iy)+nmesh*(ix))])
#define MESH_HAT(ix, iy, iz) (mesh_hat[(iz)+(nmesh/2+1)*((iy)+nmesh*(ix))])

float correction(float k, float k_N){
  float ss = sin(0.5*PI*k/k_N);
  return (1.0-SQR(ss)+0.13333333*QUAD(ss));
}

void calc_power(float* mesh, struct run_param *this_run, int nmesh)
{
  int ix, iy, iz;
  int ik;
  int nmesh_p2;

  nmesh_p2 = nmesh+2;

  static int forward_plan_created = 0;

#ifdef __FFTW2__
  static rfftwnd_plan forward_plan;
#elif __FFTW3__
  static fftwf_plan forward_plan;
#endif

  float length;
  float pk_norm;

  float *pk, *weight;

  float dk,dk_base;

  length = 1.0;
  pk_norm = SQR((float)nmesh)*SQR((float)nmesh)*SQR((float)nmesh);

#ifdef __SAVITZKY_GOLAY_FILTERING__
  dk_base = 2.0*PI/length;
  dk = dk_base/(float)(NFACT_PKBIN-1);
  if (NBIN_PK < nmesh*dk_base/dk) {
    fprintf(stderr, "NBIN_PK too small\n");
    exit(EXIT_FAILURE);
  }
#else
  dk = 2.0*PI/length;
#endif

  pk = (float *) malloc(sizeof(float)*NBIN_PK);
  weight = (float *) malloc(sizeof(float)*NBIN_PK);

  for(ik=0;ik<NBIN_PK;ik++) {
    pk[ik] = 0.0;
    weight[ik] = 0.0;
  }

#ifdef __FFTW2__
  if(forward_plan_created == 0) {
    forward_plan = 
      rfftw3d_create_plan(nmesh, nmesh, nmesh,
                          FFTW_REAL_TO_COMPLEX,
                          FFTW_ESTIMATE | FFTW_IN_PLACE);
    forward_plan_created = 1;
  }
  rfftwnd_one_real_to_complex(forward_plan, (fftw_real *)mesh,  NULL);

  fftw_complex *mesh_hat;
  mesh_hat = (fftw_complex *) mesh;
#elif __FFTW3__
  if(forward_plan_created == 0) {
    forward_plan = 
      fftwf_plan_dft_r2c_3d(nmesh, nmesh, nmesh,
			    mesh, (fftwf_complex *)mesh, FFTW_ESTIMATE);
    forward_plan_created = 1;
  }
  fftwf_execute(forward_plan);

  fftwf_complex *mesh_hat;
  mesh_hat = (fftwf_complex *) mesh;
#endif


  for(ix=0;ix<nmesh;ix++) {
    float kx, ky, kz;

    if(ix<=nmesh/2) {
      kx = (float)ix;
    }else{
      kx = (float)(nmesh-ix);
    }

    for(iy=0;iy<nmesh;iy++) {

      if(iy<=nmesh/2) {
	ky = (float)iy;
      }else{
	ky = (float)(nmesh-iy);
      }

      for(iz=0;iz<nmesh_p2/2;iz++){
	kz = (float)iz;

#ifdef __SAVITZKY_GOLAY_FILTERING__
	ik = (int) (sqrt(SQR(kx)+SQR(ky)+SQR(kz))*dk_base/dk);
#else
	ik = (int) (sqrt(SQR(kx)+SQR(ky)+SQR(kz)));
#endif

#ifdef __FFTW2__
	pk[ik] += SQR(MESH_HAT(ix,iy,iz).re)
        	 +SQR(MESH_HAT(ix,iy,iz).im);
#elif  __FFTW3__
	pk[ik] += SQR(cmplx_re(MESH_HAT(ix, iy, iz)))
	         +SQR(cmplx_im(MESH_HAT(ix, iy, iz)));
#endif
	weight[ik] += 1.0;
      }

    }

  }

  for(ik=0;ik<NBIN_PK;ik++) {
    pk[ik] /= ((weight[ik]+TINY)*pk_norm);
    pk[ik] /= correction(ik*dk, NBIN_PK/2*dk);
  }

#ifdef __SAVITZKY_GOLAY_FILTERING__

  // savitzky + golay filtering coefficient for smoothing
  coeff[0] = 4.195809E-02;
  coeff[1] =-1.048952E-01;
  coeff[2] =-2.331010E-02;
  coeff[3] = 1.398601E-01;
  coeff[4] = 2.797203E-01;
  coeff[5] = 3.333334E-01;
  coeff[6] = 2.797203E-01;
  coeff[7] = 1.398601E-01;
  coeff[8] =-2.331010E-02;
  coeff[9] =-1.048952E-01;
  coeff[10]= 4.195809E-02;

  // savitzky + golay filtering coefficient for the 1st derivative
  dcoeff[0] =  5.827501E-02;
  dcoeff[1] = -5.710955E-02;
  dcoeff[2] = -1.033411E-01;
  dcoeff[3] = -9.770782E-02;
  dcoeff[4] = -5.749804E-02;
  dcoeff[5] =  0.000000E+00;
  dcoeff[6] =  5.749804E-02;
  dcoeff[7] =  9.770782E-02;
  dcoeff[8] =  1.033411E-01;
  dcoeff[9] =  5.710955E-02;
  dcoeff[10]= -5.827501E-02;

  int ikk,k;
  for(ikk=5;ikk<NBIN_PK/2-5;ikk++) {
    float sum, dsum, kmid;
    
    sum = dsum = 0.0;
    for(k=0;k<11;k++) {
      sum += pk[i-5+k]*coeff[k];
      dsum += log(pk[i-5+k]*dcoeff[k]);
    }
    kmid = dk*((float)ikk+0.5);
    printf("%12.4e %12.4e %12.4e %12.4e\n", 
	   kmid, power[i], sum, kmid*dsum/dk);
  }
#else
  int ikk;
  for(ikk=1;ikk<NBIN_PK/2;ikk++){
    printf("%12.4e\t%12.4e\t%12.4e\n",dk*((float)ikk+0.5),pk[ikk]
           ,(log(pk[ikk])-log(pk[ikk-1]))/(log(dk*ikk)-log(dk*(ikk-1))));
  }
#endif

  free(pk);
  free(weight);
}

