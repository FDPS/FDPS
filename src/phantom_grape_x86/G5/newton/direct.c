#define ERRORTEST 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <sys/time.h>
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif
#include "gp5util.h"
#define NJMAX (JMEMSIZE)
// #define NJMAX 65536

#define rdtscll(val) do { \
	unsigned int a,d; \
	asm volatile("rdtsc" : "=a"(a), "=d"(d)); \
	(val) = ((unsigned long)a) | (((unsigned long)d)<<32); \
} while(0)

#if 1
#define GIGAHELTZ 3.8
double get_dtime(void){
  struct timeval tv;

  gettimeofday(&tv, NULL);
  return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 0.001 * 0.001);
}
#endif

void
get_cputime(double *laptime, double *sprittime);

void
readnbody(int *nj, double *mj, double (*xj)[3], double (*vj)[3], char *fname)
{
  int i, dummy, fi;
  double dummyd;
  FILE *fp;
  
  fp = fopen(fname, "r");
  if (fp == NULL)
    {
      perror("readnbody");
      exit(1);
    }
  fi = fscanf(fp, "%d\n", nj);
  fi = fscanf(fp, "%d\n", &dummy);
  fi = fscanf(fp, "%lf\n", &dummyd);
  fi = fprintf(stderr, "nj: %d\n", *nj);
  for (i = 0; i < *nj; i++)
    {
      fi = fscanf(fp, "%lf\n", mj+i);
    }
  for (i = 0; i < *nj; i++)
    {
      fi = fscanf(fp, "%lf %lf %lf\n",
		  xj[i]+0, xj[i]+1, xj[i]+2);
    }
  for (i = 0; i < *nj; i++)
    {
      fi = fscanf(fp, "%lf %lf %lf\n",
		  vj[i]+0, vj[i]+1, vj[i]+2);
    }
}

void
writenbody(int nj, double *mj, double (*xj)[3], double (*vj)[3], char *fname)
{
     int i;
     FILE *fp;

     fp = fopen(fname, "w");
     fprintf(fp, "%d\n", nj);
     fprintf(fp, "%d\n", 3);
     fprintf(fp, "%e\n", 0.0);
     for (i = 0; i < nj; i++)
     {
	  fprintf(fp, "%e\n", mj[i]);
     }
     for (i = 0; i < nj; i++)
     {
	  fprintf(fp, "%e %e %e\n",
		  xj[i][0], xj[i][1], xj[i][2]);
     }
     for (i = 0; i < nj; i++)
     {
	  fprintf(fp, "%e %e %e\n",
		  vj[i][0], vj[i][1], vj[i][2]);
     }
}

void
calc_gravity(double *mj, double (*xj)[3], double (*vj)[3],
	     double eps, double (*a)[3], double *p, int nj)
{
     double epsinv;
     int i;
     double cycle;

     // g5_set_xj(0, nj, xj);
     // g5_set_mj(0, nj, mj);
     g5_set_xmj(0, nj, xj, mj);
     g5_set_eps_to_all(eps);
     g5_set_n(nj);
     double st1 = get_dtime();
     g5_calculate_force_on_x(xj, a, p, nj);

     double st2 = get_dtime();
     cycle = (double)(st2 - st1) * GIGAHELTZ * 1e9 / ((double)nj*(double)nj/4);
#ifdef ENABLE_OPENMP
#pragma omp parallel
#if 1
	 {
		 if(omp_get_thread_num() == 0) cycle *= omp_get_num_threads();
	 }
#else
	 cycle *= omp_get_num_threads();
#endif
#endif
	 printf("gravity %f cycle per loop\n", cycle);

     for (i = 0; i < nj; i++)
     {
	  p[i] = -p[i];
     }
     if (eps != 0.0)
     {
	  epsinv = 1.0/eps;
	  for (i = 0; i < nj; i++)
	  {
	       p[i] = p[i] + mj[i] * epsinv;
	  }
     }
}

#ifdef SYMMETRIC
void
calc_gravity0(double *mj, double (*xj)[3], double (*vj)[3],
	      double *epsj2, double (*a)[3], double *p, int nj)
{
     double epsinv;
     int i;
     double cycle;

     g5_set_xmj0(0, nj, xj, mj, epsj2);
     g5_set_n(nj);
     double st1 = get_dtime();
     g5_calculate_force_on_x0(xj, a, p, nj, epsj2);

     double st2 = get_dtime();
     cycle = (double)(st2 - st1) * GIGAHELTZ * 1e9 / ((double)nj*(double)nj/4);
#ifdef ENABLE_OPENMP
#pragma omp parallel
#if 1
	 {
		 if(omp_get_thread_num() == 0) cycle *= omp_get_num_threads();
	 }
#else
	 cycle *= omp_get_num_threads();
#endif
#endif
	 printf("gravity %f cycle per loop\n", cycle);

     for (i = 0; i < nj; i++)
     {
	  p[i] = -p[i];
	  if (epsj2[i] != 0.0)
	  {
              epsinv = 1.0 / (sqrt(2.0) * sqrt(epsj2[i]));
	      p[i] = p[i] + mj[i] * epsinv;
	  }
     }
}
#endif

void
push_velocity(double (*vj)[3], double (*a)[3], double dt, int nj)
{
     int j, k;

     for (j = 0; j < nj; j++)
     {
	  for (k = 0; k < 3; k++)
	  {
	       vj[j][k] += dt * a[j][k];
	  }
     }
}

void
push_position(double (*xj)[3], double (*vj)[3], double (*a)[3],
	      double dt, int nj)
{
     int j, k;

     for (j = 0; j < nj; j++)
     {
	  for (k = 0; k < 3; k++)
	  {
	       xj[j][k] += dt * vj[j][k];
	  }
     }
}

void
energy(double *mj, double (*vj)[3], double *p, int nj, double *ke, double *pe)
{
     int i, k;
     
     *pe = 0;
     *ke = 0;
     for (i = 0; i < nj; i++)
     {
	  *pe += mj[i] * p[i];
	  for (k = 0; k < 3; k++)
	  {
	       *ke += 0.5 * mj[i] * vj[i][k] * vj[i][k];
	  }
     }
     *pe /= 2.0;
}

int
main(int argc, char **argv)
{
     static double mj[NJMAX], xj[NJMAX][3], vj[NJMAX][3], epsj2[NJMAX];
     static double a[NJMAX][3], p[NJMAX];
     double xmax, xmin, mmin;
     double time;
//     double eps, dt, endt;
     double dt, endt;
     double e, e0, ke, pe;
     double LapTime, SpritTime, IntPerSec, Gflops;

     int nj;
     int nstep, step;
     dt = 0.01;
     endt = 10.0;
     time = 0.0;
     nstep = endt/dt;
     xmax = 10.0;
     xmin = -10.0;

     if (argc < 3)
     {
	  fprintf(stderr, "usage: %s <infile> <outfile>\n",  argv[0]);
	  exit(1);
     }
  
     readnbody(&nj, mj, xj, vj, argv[1]);
     mmin = mj[0];
#if ERRORTEST == 1
     double eps;
     eps = 4.0 / (double)nj;
#else
#ifdef SYMMETRIC
     int i;
     for(i = 0; i < nj; i++)
       epsj2[i] = (0.01 + 0.01 * (double)i / (double)nj) * (0.01 + 0.01 * (double)i / (double)nj);
     //     mj[1021] = 1.0;
     //     mj[1022] = 1.0;
     //     mj[1023] = 1.0;
#else
     double eps;
     eps = 0.02;
#endif
#endif
     g5_open();
     g5_set_range(xmin, xmax, mmin);
#ifdef SYMMETRIC
     calc_gravity0(mj, xj, vj, epsj2, a, p, nj);
#else
     calc_gravity(mj, xj, vj, eps, a, p, nj);
#endif
     energy(mj, vj, p, nj, &ke, &pe);
     e0 = ke+pe;

#if ERRORTEST == 1
     int i;
     char out[1024];
     FILE *fp;
     sprintf(out, "pl%03dk_eps4n_avx.ap", nj / 1024);
     fp = fopen(out, "w");
     for(i = 0; i < nj; i++)
       fprintf(fp, "%5d %+.16e %+.16e\n", i, sqrt(a[i][0]*a[i][0]+a[i][1]*a[i][1]+a[i][2]*a[i][2]), p[i]);
     fclose(fp);
     exit(0);
#endif

	 // TimeStart = (double)clock() / CLOCKS_PER_SEC;
	 get_cputime(&LapTime, &SpritTime);
     for (step = 1; step < nstep; step++)
     {
	  push_velocity(vj, a, 0.5*dt, nj);
	  push_position(xj, vj, a, dt, nj);
	  time = time + dt;
#ifdef SYMMETRIC
	  calc_gravity0(mj, xj, vj, epsj2, a, p, nj);
#else
	  calc_gravity(mj, xj, vj, eps, a, p, nj);
#endif
	  push_velocity(vj, a, 0.5*dt, nj);
#ifdef ANIM
	  plot_star(xj, nj, time, 0.3, mj, mj[0]);
#endif /* ANIM */

	  if (step % (nstep/10) == 0) {
	      energy(mj, vj, p, nj, &ke, &pe);
	      e = ke+pe;
		  // TimeEnd = (double)clock() / CLOCKS_PER_SEC;
		 get_cputime(&LapTime, &SpritTime);
		 IntPerSec = ((double)nj * (double)nj * (long)(nstep/10)) / LapTime;
		  Gflops = IntPerSec * 38. * 1.e-9; 
	      printf("step: %d time: %e\n", step, time);
	      printf("e: %e de: %e\n", e, e-e0);
	      printf("ke: %e pe: %e\n", ke, pe);
	      printf("ke/pe: %e\n\n", ke/pe);
	      printf("%e interaction per sec, %f Gflops \n", IntPerSec, Gflops);
		  // TimeStart = TimeEnd;
	  }
     }
     g5_close();
     writenbody(nj, mj, xj, vj, argv[2]);

     return 0;

}

