#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
#include "avx_type.h"
#include "gravity.h"
#include "timeprof.h"
#include "gp6util.h"

static int flag_predict_j_particle;
static PrdPosVel prdposvel[NPIPES];
static NewAccJrk newaccjrk[NPIPES];
static double acccorr = - 0.125;
static double potcorr =   0.5;

// ----------------------------------------------------

void g6_open(int gpid)
{
  int nthread;

  nthread = omp_get_max_threads();
  avx_open(nthread);

  return;
}

void g6_open_(int *gpid)
{
  g6_open(*gpid);
  return;
}

// ----------------------------------------------------

void g6_set_tunit(int tunit)
{
  // do nothing
  return;
}

void g6_set_tunit_(int *tunit)
{
  g6_set_tunit(*tunit);
  return;
}

// ----------------------------------------------------

void g6_set_xunit(int xunit)
{
  // do nothing
  return;
}

void g6_set_xunit_(int *xunit)
{
  g6_set_xunit(*xunit);
  return;
}

// ----------------------------------------------------

void g6_reset(int gpid)
{
  // do nothing
  return;
}

void g6_reset_(int *gpid)
{
  g6_reset(*gpid);
  return;
}

// ----------------------------------------------------

void g6_reset_fofpga(int gpid)
{
  // do nothing
  return;
}

void g6_reset_fofpga_(int *gpid)
{
  g6_reset_fofpga(*gpid);
  return;
}

// ----------------------------------------------------

void g6_close(int gpid)
{
  avx_close();
  return;
}

void g6_close_(int *gpid)
{
  g6_close(*gpid);
  return;
}

// ----------------------------------------------------

int g6_npipes(void)
{
  return NPIPES;
}

int g6_npipes_(void)
{
  return g6_npipes();
}

// ----------------------------------------------------

void g6_set_j_particle(int gpid, int padr, int pidx,
		       double t, double dt, double mss, 
		       double *snp, double *jrk, double *acc,
		       double *vel, double *pos)
{
  avx_set_j_particle(padr, pidx, t, mss, pos, vel, acc, jrk);
  return;
}

void g6_set_j_particle_(int *gpid, int *padr, int *pidx,
		       double *t, double *dt, double *mss, 
		       double *snp, double *jrk, double *acc,
		       double *vel, double *pos)
{
  g6_set_j_particle(*gpid, *padr, *pidx, *t, *dt, *mss, snp, jrk, acc, vel, pos);
  return;
}

// ----------------------------------------------------

void g6_set_ti(int gpid, double time)
{
  avx_set_ti(time);
  flag_predict_j_particle = 1;
  return;
}

void g6_set_ti_(int *gpid, double *time)
{
  g6_set_ti(*gpid, *time);
  return;
}

// ----------------------------------------------------

void g6calc_firsthalf(int gpid, int nj, int ni, int *pidx,
		      double (*pos)[3], double (*vel)[3],
		      double (*acc)[3], double (*jrk)[3], double *pot,
		      double eps2, double *h2)
{
  if(flag_predict_j_particle == 1){
    avx_predict_j_particle(nj);
    flag_predict_j_particle = 0;
  }
  avx_initialize_neighbourlist();
  return;
}

void g6calc_firsthalf_(int *gpid, int *nj, int *ni, int *pidx,
		      double (*pos)[3], double (*vel)[3],
		      double (*acc)[3], double (*jrk)[3], double *pot,
		      double *eps2, double *h2)
{
  g6calc_firsthalf(*gpid, *nj, *ni, pidx, pos, vel, acc, jrk, pot, *eps2, h2);
  return;
}

// ----------------------------------------------------

int g6calc_lasthalf(int gpid, int nj, int ni, int *pidx,
		    double (*pos)[3], double (*vel)[3], double eps2, double *h2,
		    double (*acc)[3], double (*jrk)[3], double *pot)
{
  int i;

  for(i = 0; i < ni; i++){
    prdposvel[i].xpos = pos[i][0];
    prdposvel[i].ypos = pos[i][1];
    prdposvel[i].zpos = pos[i][2];
    prdposvel[i].xvel = vel[i][0];
    prdposvel[i].yvel = vel[i][1];
    prdposvel[i].zvel = vel[i][2];
    prdposvel[i].id   = (float)pidx[i];
    prdposvel[i].eps2 = (float)eps2;
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < ni; i += 2)
    gravity_kernel(nj, &prdposvel[i], &newaccjrk[i]);

  for(i = 0; i < ni; i++){
    acc[i][0] = acccorr * newaccjrk[i].xacc;
    acc[i][1] = acccorr * newaccjrk[i].yacc;
    acc[i][2] = acccorr * newaccjrk[i].zacc;
    jrk[i][0] = acccorr * newaccjrk[i].xjrk;
    jrk[i][1] = acccorr * newaccjrk[i].yjrk;
    jrk[i][2] = acccorr * newaccjrk[i].zjrk;    
    pot[i]    = potcorr * newaccjrk[i].pot;
  }  

  return 0;
}

int g6calc_lasthalf_(int *gpid, int *nj, int *ni, int *pidx,
		    double (*pos)[3], double (*vel)[3], double *eps2, double *h2,
		    double (*acc)[3], double (*jrk)[3], double *pot)
{
  
  return g6calc_lasthalf(*gpid, *nj, *ni, pidx, pos, vel, *eps2, h2, acc, jrk, pot);
}

// ----------------------------------------------------

int g6calc_lasthalf2(int gpid, int nj, int ni, int *pidx,
		     double (*pos)[3], double (*vel)[3], double eps2, double *h2,
		     double (*acc)[3], double (*jrk)[3], double *pot, int *nnb)
{
  int i;

  for(i = 0; i < ni; i++){
    prdposvel[i].xpos = pos[i][0];
    prdposvel[i].ypos = pos[i][1];
    prdposvel[i].zpos = pos[i][2];
    prdposvel[i].xvel = vel[i][0];
    prdposvel[i].yvel = vel[i][1];
    prdposvel[i].zvel = vel[i][2];
    prdposvel[i].id   = (float)pidx[i];
    prdposvel[i].eps2 = (float)eps2;
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < ni; i += 2)
    gravity_kernel2(nj, &prdposvel[i], &newaccjrk[i]);

  for(i = 0; i < ni; i++){
    acc[i][0] = acccorr * newaccjrk[i].xacc;
    acc[i][1] = acccorr * newaccjrk[i].yacc;
    acc[i][2] = acccorr * newaccjrk[i].zacc;
    jrk[i][0] = acccorr * newaccjrk[i].xjrk;
    jrk[i][1] = acccorr * newaccjrk[i].yjrk;
    jrk[i][2] = acccorr * newaccjrk[i].zjrk;    
    pot[i]    = potcorr * newaccjrk[i].pot;
    nnb[i]    =           newaccjrk[i].nnb;
  }  

  return 0;
}

int g6calc_lasthalf2_(int *gpid, int *nj, int *ni, int *pidx,
		    double (*pos)[3], double (*vel)[3], double *eps2, double *h2,
		      double (*acc)[3], double (*jrk)[3], double *pot, int *nnb)
{
  
  return g6calc_lasthalf2(*gpid, *nj, *ni, pidx, pos, vel, *eps2, h2, acc, jrk, pot, nnb);
}

// ----------------------------------------------------

int g6calc_lasthalf2p(int gpid, int nj, int ni, int *pidx,
		      double (*pos)[3], double (*vel)[3], double eps2, double *h2,
		      double (*acc)[3], double (*jrk)[3], double *pot, int *nnb, float *rnnb)
{
  int i;

  for(i = 0; i < ni; i++){
    prdposvel[i].xpos = pos[i][0];
    prdposvel[i].ypos = pos[i][1];
    prdposvel[i].zpos = pos[i][2];
    prdposvel[i].xvel = vel[i][0];
    prdposvel[i].yvel = vel[i][1];
    prdposvel[i].zvel = vel[i][2];
    prdposvel[i].id   = (float)pidx[i];
    prdposvel[i].eps2 = (float)eps2;
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < ni; i += 2)
    gravity_kernel2(nj, &prdposvel[i], &newaccjrk[i]);

  for(i = 0; i < ni; i++){
    acc[i][0] = acccorr * newaccjrk[i].xacc;
    acc[i][1] = acccorr * newaccjrk[i].yacc;
    acc[i][2] = acccorr * newaccjrk[i].zacc;
    jrk[i][0] = acccorr * newaccjrk[i].xjrk;
    jrk[i][1] = acccorr * newaccjrk[i].yjrk;
    jrk[i][2] = acccorr * newaccjrk[i].zjrk;    
    pot[i]    = potcorr * newaccjrk[i].pot;
    nnb[i]    =           newaccjrk[i].nnb;
    rnnb[i]   =           newaccjrk[i].rnnb;
  }  

  return 0;
}

int g6calc_lasthalf2p_(int *gpid, int *nj, int *ni, int *pidx,
		       double (*pos)[3], double (*vel)[3], double *eps2, double *h2,
		       double (*acc)[3], double (*jrk)[3], double *pot, int *nnb, float *rnnb)
{
  
  return g6calc_lasthalf2p(*gpid, *nj, *ni, pidx, pos, vel, *eps2, h2, acc, jrk, pot, nnb, rnnb);
}

// ----------------------------------------------------

int g6calc_lasthalfn(int gpid, int nj, int ni, int *pidx,
		    double (*pos)[3], double (*vel)[3], double eps2, double *h2,
		    double (*acc)[3], double (*jrk)[3], double *pot)
{
  int i;

  for(i = 0; i < ni; i++){
    prdposvel[i].xpos = pos[i][0];
    prdposvel[i].ypos = pos[i][1];
    prdposvel[i].zpos = pos[i][2];
    prdposvel[i].xvel = vel[i][0];
    prdposvel[i].yvel = vel[i][1];
    prdposvel[i].zvel = vel[i][2];
    prdposvel[i].id   = (float)pidx[i];
    prdposvel[i].eps2 = (float)eps2;
    prdposvel[i].h2   = (float)h2[i];
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < ni; i += 2){
    int ithread = omp_get_thread_num();
    gravity_kerneln(nj, &prdposvel[i], &newaccjrk[i], i, ithread);
  }

  for(i = 0; i < ni; i++){
    acc[i][0] = acccorr * newaccjrk[i].xacc;
    acc[i][1] = acccorr * newaccjrk[i].yacc;
    acc[i][2] = acccorr * newaccjrk[i].zacc;
    jrk[i][0] = acccorr * newaccjrk[i].xjrk;
    jrk[i][1] = acccorr * newaccjrk[i].yjrk;
    jrk[i][2] = acccorr * newaccjrk[i].zjrk;    
    pot[i]    = potcorr * newaccjrk[i].pot;
  }  

  return 0;
}

int g6calc_lasthalfn_(int *gpid, int *nj, int *ni, int *pidx,
		    double (*pos)[3], double (*vel)[3], double *eps2, double *h2,
		    double (*acc)[3], double (*jrk)[3], double *pot)
{
  
  return g6calc_lasthalfn(*gpid, *nj, *ni, pidx, pos, vel, *eps2, h2, acc, jrk, pot);
}

// ----------------------------------------------------

int g6calc_lasthalf2n(int gpid, int nj, int ni, int *pidx,
		     double (*pos)[3], double (*vel)[3], double eps2, double *h2,
		     double (*acc)[3], double (*jrk)[3], double *pot, int *nnb)
{
  int i;

  for(i = 0; i < ni; i++){
    prdposvel[i].xpos = pos[i][0];
    prdposvel[i].ypos = pos[i][1];
    prdposvel[i].zpos = pos[i][2];
    prdposvel[i].xvel = vel[i][0];
    prdposvel[i].yvel = vel[i][1];
    prdposvel[i].zvel = vel[i][2];
    prdposvel[i].id   = (float)pidx[i];
    prdposvel[i].eps2 = (float)eps2;
    prdposvel[i].h2   = (float)h2[i];
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < ni; i += 2){
    int ithread = omp_get_thread_num();
    gravity_kernel2n(nj, &prdposvel[i], &newaccjrk[i], i, ithread);
  }

  for(i = 0; i < ni; i++){
    acc[i][0] = acccorr * newaccjrk[i].xacc;
    acc[i][1] = acccorr * newaccjrk[i].yacc;
    acc[i][2] = acccorr * newaccjrk[i].zacc;
    jrk[i][0] = acccorr * newaccjrk[i].xjrk;
    jrk[i][1] = acccorr * newaccjrk[i].yjrk;
    jrk[i][2] = acccorr * newaccjrk[i].zjrk;    
    pot[i]    = potcorr * newaccjrk[i].pot;
    nnb[i]    =           newaccjrk[i].nnb;
  }  

  return 0;
}

int g6calc_lasthalf2n_(int *gpid, int *nj, int *ni, int *pidx,
		    double (*pos)[3], double (*vel)[3], double *eps2, double *h2,
		      double (*acc)[3], double (*jrk)[3], double *pot, int *nnb)
{
  
  return g6calc_lasthalf2n(*gpid, *nj, *ni, pidx, pos, vel, *eps2, h2, acc, jrk, pot, nnb);
}

// ----------------------------------------------------

int g6calc_lasthalf2np(int gpid, int nj, int ni, int *pidx,
		       double (*pos)[3], double (*vel)[3], double eps2, double *h2,
		       double (*acc)[3], double (*jrk)[3], double *pot, int *nnb, float *rnnb)
{
  int i;

  for(i = 0; i < ni; i++){
    prdposvel[i].xpos = pos[i][0];
    prdposvel[i].ypos = pos[i][1];
    prdposvel[i].zpos = pos[i][2];
    prdposvel[i].xvel = vel[i][0];
    prdposvel[i].yvel = vel[i][1];
    prdposvel[i].zvel = vel[i][2];
    prdposvel[i].id   = (float)pidx[i];
    prdposvel[i].eps2 = (float)eps2;
    prdposvel[i].h2   = (float)h2[i];
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < ni; i += 2){
    int ithread = omp_get_thread_num();
    gravity_kernel2n(nj, &prdposvel[i], &newaccjrk[i], i, ithread);
  }

  for(i = 0; i < ni; i++){
    acc[i][0] = acccorr * newaccjrk[i].xacc;
    acc[i][1] = acccorr * newaccjrk[i].yacc;
    acc[i][2] = acccorr * newaccjrk[i].zacc;
    jrk[i][0] = acccorr * newaccjrk[i].xjrk;
    jrk[i][1] = acccorr * newaccjrk[i].yjrk;
    jrk[i][2] = acccorr * newaccjrk[i].zjrk;    
    pot[i]    = potcorr * newaccjrk[i].pot;
    nnb[i]    =           newaccjrk[i].nnb;
    rnnb[i]   =           newaccjrk[i].rnnb;
  }  

  return 0;
}

int g6calc_lasthalf2np_(int *gpid, int *nj, int *ni, int *pidx,
			double (*pos)[3], double (*vel)[3], double *eps2, double *h2,
			double (*acc)[3], double (*jrk)[3], double *pot, int *nnb, float *rnnb)
{
  
  return g6calc_lasthalf2np(*gpid, *nj, *ni, pidx, pos, vel, *eps2, h2, acc, jrk, pot, nnb, rnnb);
}

// ----------------------------------------------------

void g6_setup_njdata(int gpid, int nj)
{
  // do nothing
  return;
}

void g6_setup_njdata_(int *gpid, int *nj)
{
  g6_setup_njdata(*gpid, *nj);
  return;
}

// ----------------------------------------------------

int g6_read_neighbour_list(int gpid)
{
  return avx_get_neighbourlist_error();
  // 0: success, -1: hard err, 1: overflow of hardware memory
}

int g6_read_neighbour_list_(int *gpid)
{
  return g6_read_neighbour_list(*gpid);;
}

// ----------------------------------------------------

int g6_get_neighbour_list(int gpid, int ipipe, int maxlen, int *nblen, int *nbl)
{
  return avx_get_neighbourlist(ipipe, maxlen, nblen, nbl);
  // 0: success, 1: nb > maxlen
}

int g6_get_neighbour_list_(int *gpid, int *ipipe, int *maxlen, int *nblen, int *nbl)
{
  return g6_get_neighbour_list(*gpid, *ipipe, *maxlen, nblen, nbl);;
}

// ----------------------------------------------------

void g6_debugfunc(int tmp)
{
  //  avx_debugfunc();
  FILE *fp;

  fp = fopen("dump.dat", "a");
  fprintf(fp, "%d\n", tmp);
  fclose(fp);
  return;
}

void g6_debugfunc_(int *tmp)
{
  g6_debugfunc(*tmp);
  return;
}

// ----------------------------------------------------

void g6_debugfunc_double(double tmp)
{
  FILE *fp;

  fp = fopen("dump.dat", "a");
  fprintf(fp, "%+e\n", tmp);
  fclose(fp);
  return;
}

void g6_debugfunc_double_(double *tmp)
{
  g6_debugfunc_double(*tmp);
  return;
}

// ----------------------------------------------------

void g6_dump(double tim, int n, double *m, double (*x)[3], double (*v)[3])
{
  int i;
  char out[1024];
  FILE *fp;

  sprintf(out, "dump.dat.%05d", (int)tim);
  fp = fopen(out, "w");
  for(i = 0; i < n; i++){
    fprintf(fp, "%+.16e", m[i]);
    fprintf(fp, " %+.16e %+.16e %+.16e", x[i][0], x[i][1], x[i][2]);
    fprintf(fp, " %+.16e %+.16e %+.16e", v[i][0], v[i][1], v[i][2]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return;
}

void g6_dump_(double *tim, int *n, double *m, double (*x)[3], double (*v)[3])
{
  g6_dump(*tim, *n, m, x, v);
  return;
}
