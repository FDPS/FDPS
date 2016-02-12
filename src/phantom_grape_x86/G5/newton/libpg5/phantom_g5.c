#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <immintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#include "avx_type.h"
#include "gp5util.h"

#define NUM_PIPE (4)

#ifndef MAXDEV
#define MAXDEV (128)
#endif /* MAXDEV */

static double Eps;

static struct Ptcl_Mem {
  Ipdata iptcl;
  Fodata fout;
  Jpdata jptcl[JMEMSIZE];
  Jpdata0 jptcl0[JMEMSIZE/2];
  int nbody, pad[15];
} ptcl_mem[MAXDEV] ALIGN64;

static float Acc_correct = 1.0;
static float Pot_correct = -1.0;
static __m128 Acc_correctV;
static __m128 Pot_correctV;

int g5_get_number_of_pipelines(void) 
{
  return NUM_PIPE;
}

int g5_get_jmemsize(void) 
{
  return JMEMSIZE;
}

void g5_open(void)
{
  static int init_call = 1;
  if(init_call) {
    double rsqrt_bias();
    double bias = rsqrt_bias();
    float acc_corr = 1.0 - 3.0*bias;
    float pot_corr = -(1.0-bias);
    Acc_correct = acc_corr;
    Pot_correct = pot_corr;
    Acc_correctV = _mm_set_ps(acc_corr,acc_corr,acc_corr,acc_corr);
    Pot_correctV = _mm_set_ps(pot_corr,pot_corr,pot_corr,pot_corr);
    init_call = 0;
  }
  return;
}

void g5_close() 
{
  return;
}

void g5_set_eps_to_all(double eps) 
{
  Eps = eps;
}

void g5_set_range(double xmin, double xmax, double mmin)
{
  return;
}

void g5_set_nMC(int devid, int n)
{
    assert(devid < MAXDEV);
  struct Ptcl_Mem *pm = ptcl_mem + devid;
  pm->nbody = n;
}

void g5_set_n(int n) 
{
  g5_set_nMC(0, n);
}

void g5_set_xiMC(int devid, int ni, double (*xi)[3])
{
  int i;
    assert(devid < MAXDEV);
  struct Ptcl_Mem *pm = ptcl_mem + devid;

  assert(ni <= NUM_PIPE);
  for(i=0;i<ni;i++) {
    float eps2 = Eps*Eps;
    pm->iptcl.x[i] = (float)xi[i][0];
    pm->iptcl.y[i] = (float)xi[i][1];
    pm->iptcl.z[i] = (float)xi[i][2];
    pm->iptcl.eps2[i] = eps2;
  }
}

void g5_set_xiMC0(int devid, int ni, double (*xi)[3],
		  double *eps2)
{
  int i;
    assert(devid < MAXDEV);
  struct Ptcl_Mem *pm = ptcl_mem + devid;

  assert(ni <= NUM_PIPE);
  for(i=0;i<ni;i++) {
    pm->iptcl.x[i] = (float)xi[i][0];
    pm->iptcl.y[i] = (float)xi[i][1];
    pm->iptcl.z[i] = (float)xi[i][2];
    pm->iptcl.eps2[i] = eps2[i];
  }
}

void g5_set_xi(int ni, double (*xi)[3]) 
{
  g5_set_xiMC(0, ni, xi);
}

void g5_set_xmjMC(int devid, int adr, int nj, double (*xj)[3], double *mj) 
{
  int j;
    assert(devid < MAXDEV);
  struct Ptcl_Mem *pm = ptcl_mem + devid;

  for(j=adr;j<adr+nj;j++) {
    __m256d pd = {xj[j][0], xj[j][1], xj[j][2], mj[j]};
    __m128  ps = _mm256_cvtpd_ps(pd);
    *(__m128 *)(pm->jptcl+j) = ps;
  }

  int rsdl = (NUNROLL - (nj % NUNROLL)) % NUNROLL;
  for(j=nj;j<nj+rsdl;j++){
    __m256d pd = {0.0, 0.0, 0.0, 0.0};
    __m128  ps = _mm256_cvtpd_ps(pd);
    *(__m128 *)(pm->jptcl+j) = ps;
  }
}

void g5_set_xmjMC0(int devid, int adr, int nj, double (*xj)[3],
		   double *mj, double *epsj2) 
{
  int j;
    assert(devid < MAXDEV);
  struct Ptcl_Mem *pm = ptcl_mem + devid;

  assert(adr % 2 == 0);
  for(j=adr;j<adr+nj;j+=2) {
    int jadr = j / 2;
    pm->jptcl0[jadr].xm[0][0] = (float)xj[j][0];
    pm->jptcl0[jadr].xm[0][1] = (float)xj[j][1];
    pm->jptcl0[jadr].xm[0][2] = (float)xj[j][2];
    pm->jptcl0[jadr].xm[0][3] = (float)mj[j];
    pm->jptcl0[jadr].xm[1][0] = (float)xj[j+1][0];
    pm->jptcl0[jadr].xm[1][1] = (float)xj[j+1][1];
    pm->jptcl0[jadr].xm[1][2] = (float)xj[j+1][2];
    pm->jptcl0[jadr].xm[1][3] = (float)mj[j+1];
    pm->jptcl0[jadr].ep[0][0] = (float)epsj2[j];
    pm->jptcl0[jadr].ep[0][1] = (float)epsj2[j];
    pm->jptcl0[jadr].ep[0][2] = (float)epsj2[j];
    pm->jptcl0[jadr].ep[0][3] = (float)epsj2[j];
    pm->jptcl0[jadr].ep[1][0] = (float)epsj2[j+1];
    pm->jptcl0[jadr].ep[1][1] = (float)epsj2[j+1];
    pm->jptcl0[jadr].ep[1][2] = (float)epsj2[j+1];
    pm->jptcl0[jadr].ep[1][3] = (float)epsj2[j+1];
  }

  int rsdl = (NUNROLL - (nj % NUNROLL)) % NUNROLL;
  for(j=nj;j<nj+rsdl;j+=2){
      int jj, jadr = j / 2;
      for(jj = 0; jj < 2; jj++){
          int jp = jadr * 2 + jj;
          if(jp < nj)
              continue;
          pm->jptcl0[jadr].xm[jj][0] = 0.0f;
          pm->jptcl0[jadr].xm[jj][1] = 0.0f;
          pm->jptcl0[jadr].xm[jj][2] = 0.0f;
          pm->jptcl0[jadr].xm[jj][3] = 0.0f;
          pm->jptcl0[jadr].ep[jj][0] = 1.0f;
          pm->jptcl0[jadr].ep[jj][1] = 1.0f;
          pm->jptcl0[jadr].ep[jj][2] = 1.0f;
          pm->jptcl0[jadr].ep[jj][3] = 1.0f;
      }
  }

}

void g5_set_xmj(int adr, int nj, double (*xj)[3], double *mj) 
{
  g5_set_xmjMC(0, adr, nj, xj, mj);
}

void g5_set_xmj0(int adr, int nj, double (*xj)[3],
		 double *mj, double *epsj2) 
{
  g5_set_xmjMC0(0, adr, nj, xj, mj, epsj2);
}

void g5_runMC(int devid) 
{
    assert(devid < MAXDEV);
  struct Ptcl_Mem *pm = ptcl_mem + devid;
  void GravityKernel(pIpdata, pFodata, pJpdata, int);
  GravityKernel(&(pm->iptcl), &(pm->fout), pm->jptcl, pm->nbody);
}

void g5_runMC0(int devid) 
{
    assert(devid < MAXDEV);
  struct Ptcl_Mem *pm = ptcl_mem + devid;
  void GravityKernel0(pIpdata, pFodata, pJpdata0, int);
  GravityKernel0(&(pm->iptcl), &(pm->fout), pm->jptcl0, pm->nbody);
}

void g5_run(void) 
{
  g5_runMC(0);
}

void g5_get_forceMC(int devid, int ni, double (*a)[3], double *pot) 
{
    assert(devid < MAXDEV);
  assert(ni <= NUM_PIPE);

  struct Ptcl_Mem *pm = ptcl_mem + devid;
  int i;

#if 1
  *(__m128 *)pm->fout.ax = 
    _mm_mul_ps(*(__m128 *)(pm->fout.ax),Acc_correctV);
  *(__m128 *)pm->fout.ay = 
    _mm_mul_ps(*(__m128 *)(pm->fout.ay),Acc_correctV);
  *(__m128 *)pm->fout.az = 
    _mm_mul_ps(*(__m128 *)(pm->fout.az),Acc_correctV);
  *(__m128 *)pm->fout.phi = 
    _mm_mul_ps(*(__m128 *)(pm->fout.phi),Pot_correctV);
#endif

  for(i=0;i<ni;i++) {
    a[i][0] = (double)(pm->fout.ax[i]);
    a[i][1] = (double)(pm->fout.ay[i]);
    a[i][2] = (double)(pm->fout.az[i]);
    pot[i]  = (double)(pm->fout.phi[i]);
  }

}

void g5_get_force(int ni, double (*a)[3], double *pot) 
{
  g5_get_forceMC(0, ni, a, pot);
}

void g5_calculate_force_on_xMC(int devid, double (*x)[3], double (*a)[3], 
                               double *p, int ni)
{
    assert(devid < MAXDEV);
  int off;
  int np = g5_get_number_of_pipelines();
  for(off=0;off<ni;off+=np) {
    int nii = np < ni-off ? np : ni-off;
    g5_set_xiMC(devid, nii, x+off);
    g5_runMC(devid);
    g5_get_forceMC(devid, nii, a+off, p+off);
  }
}

void g5_calculate_force_on_xMC0(int devid, double (*x)[3],
				double (*a)[3], double *p,
				int ni, double *eps2)
{
    assert(devid < MAXDEV);
  int off;
  int np = g5_get_number_of_pipelines();
  for(off=0;off<ni;off+=np) {
    int nii = np < ni-off ? np : ni-off;
    g5_set_xiMC0(devid, nii, x+off, eps2+off);
    g5_runMC0(devid);
    g5_get_forceMC(devid, nii, a+off, p+off);
  }
}

#ifndef ENABLE_OPENMP
void g5_calculate_force_on_x(double (*x)[3], double (*a)[3], double *p, int ni)
{
  g5_calculate_force_on_xMC(0, x, a, p, ni);
}

void g5_calculate_force_on_x0(double (*x)[3], double (*a)[3],
			      double *p, int ni, double *eps2)
{
  g5_calculate_force_on_xMC0(0, x, a, p, ni, eps2);
}
#else
#include <omp.h>
void g5_calculate_force_on_x(double (*x)[3], double (*a)[3], double *p, 
                             int nitot)
{
  int off;
  const int np = g5_get_number_of_pipelines();
#pragma omp parallel for
  for(off=0; off<nitot; off+=np) {
    int tid = omp_get_thread_num();
    int ni = np < nitot-off ? np : nitot-off;
    g5_set_xiMC(tid, ni, x+off);
    {
      void GravityKernel(pIpdata, pFodata, pJpdata, int);
      pIpdata ip = &ptcl_mem[tid].iptcl;
      pFodata fo = &ptcl_mem[tid].fout;
      pJpdata jp = ptcl_mem[0].jptcl;
      int nbody  = ptcl_mem[0].nbody;
      GravityKernel(ip, fo, jp, nbody);
    }
    g5_get_forceMC(tid, ni, a+off, p+off);
  }
}

void g5_calculate_force_on_x0(double (*x)[3], double (*a)[3], double *p, 
			      int nitot, double *eps2)
{
  int off;
  const int np = g5_get_number_of_pipelines();
#pragma omp parallel for
  for(off=0; off<nitot; off+=np) {
    int tid = omp_get_thread_num();
    int ni = np < nitot-off ? np : nitot-off;
    g5_set_xiMC0(tid, ni, x+off, eps2+off);
    {
      void GravityKernel0(pIpdata, pFodata, pJpdata0, int);
      pIpdata ip = &ptcl_mem[tid].iptcl;
      pFodata fo = &ptcl_mem[tid].fout;
      pJpdata0 jp = ptcl_mem[0].jptcl0;
      int nbody  = ptcl_mem[0].nbody;
      GravityKernel0(ip, fo, jp, nbody);
    }
    g5_get_forceMC(tid, ni, a+off, p+off);
  }
}
#endif
