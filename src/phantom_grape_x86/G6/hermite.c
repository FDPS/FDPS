#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "banana01.h"
#include "particle.h"
#include "globaltime.h"
#include "calc_energy.h"
#include "timestep.h"
#include "timeprof.h"
#include "hermite.h"
#include "gp6util.h"

//#define NPIPES 48
#define NPIPES 128

#define SECONDORDER 5e-1
#define THIRDORDER  1.666666666666667e-1
#define FOURTHORDER 4.166666666666667e-2
#define FIFTHORDER  8.333333333333333e-3

static REAL etas, etac, etas_g, etac_g;
static struct Predposvel  posvel[NMAX],  posvel_in[NMAX];
static struct NewAccJerk accjerk[NMAX], accjerk_in[NMAX];
static int neighbour[NMAX][1024];

void initialize_accuracy(REAL etastmp, REAL etactmp)
{
  int i;

  etas   = etastmp;
  etac   = etactmp;
  etas_g = etas / 4.0;
  etac_g = etac / 4.0;
  for(i = 0; i < NMAX; i++)
    accjerk[i].nptr = neighbour[i];
  return;
}

static void predict_exmotion(REAL t_gl, int ni, int *index, struct Particle *particle, struct Address *address)
{
  int i, iadr, padr, id;
  REAL dt, dt2, dt3;

  for(i = 0; i < ni; i++){
    iadr = index[i];
    id   = particle[iadr].id;
    dt   = t_gl - particle[iadr].t;
    dt2 = dt * dt;
    dt3 = dt * dt2;

    posvel[i].xpos = particle[iadr].xpos + particle[iadr].xvel * dt + particle[iadr].xacc * dt2 * SECONDORDER + particle[iadr].xjrk * dt3 * THIRDORDER;
    posvel[i].ypos = particle[iadr].ypos + particle[iadr].yvel * dt + particle[iadr].yacc * dt2 * SECONDORDER + particle[iadr].yjrk * dt3 * THIRDORDER;
    posvel[i].zpos = particle[iadr].zpos + particle[iadr].zvel * dt + particle[iadr].zacc * dt2 * SECONDORDER + particle[iadr].zjrk * dt3 * THIRDORDER;
    posvel[i].xvel = particle[iadr].xvel + particle[iadr].xacc * dt + particle[iadr].xjrk * dt2 * SECONDORDER;
    posvel[i].yvel = particle[iadr].yvel + particle[iadr].yacc * dt + particle[iadr].yjrk * dt2 * SECONDORDER;
    posvel[i].zvel = particle[iadr].zvel + particle[iadr].zacc * dt + particle[iadr].zjrk * dt2 * SECONDORDER;
    posvel[i].id   = (float)id;
  }
  return;
}

// ---------------------------------------------------------------------------

static void correct(struct Particle *pi0, struct Particle *pi1, REAL s0[DIM], REAL c0[DIM])
{
  int k;
  REAL dts, dtsinv;
  REAL a0[DIM], j0[DIM];
  REAL x1[DIM], v1[DIM], a1[DIM], j1[DIM];

  dts    = pi0->dt;
  dtsinv = 1.0 / dts;

  a0[0] = pi0->xacc;
  a0[1] = pi0->yacc;
  a0[2] = pi0->zacc;
  j0[0] = pi0->xjrk;
  j0[1] = pi0->yjrk;
  j0[2] = pi0->zjrk;

  x1[0] = pi1->xpos;
  x1[1] = pi1->ypos;
  x1[2] = pi1->zpos;
  v1[0] = pi1->xvel;
  v1[1] = pi1->yvel;
  v1[2] = pi1->zvel;
  a1[0] = pi1->xacc;
  a1[1] = pi1->yacc;
  a1[2] = pi1->zacc;
  j1[0] = pi1->xjrk;
  j1[1] = pi1->yjrk;
  j1[2] = pi1->zjrk;

  for(k = 0; k < DIM; k++){
    s0[k]  = 2.0 * (- 3.0 * (a0[k] - a1[k]) - (2.0 * j0[k] + j1[k]) * dts) * dtsinv * dtsinv;
    c0[k]  = 6.0 * (2.0 * (a0[k] - a1[k]) + (j0[k] + j1[k]) * dts) * dtsinv * dtsinv * dtsinv;
    x1[k] += s0[k] * dts * dts * dts * dts * FOURTHORDER + c0[k] * dts * dts * dts * dts * dts * FIFTHORDER;
    v1[k] += s0[k] * dts * dts * dts * THIRDORDER + c0[k] * dts * dts * dts * dts * FOURTHORDER;
  }

  pi1->xpos = x1[0];
  pi1->ypos = x1[1];
  pi1->zpos = x1[2];
  pi1->xvel = v1[0];
  pi1->yvel = v1[1];
  pi1->zvel = v1[2];

  return;
}

// ---------------------------------------------------------------------------

void startup_sgmotion(int ni, int *index, struct Particle *particle)
{
  int i, iadr;
  struct Particle p0;

  for(i = 0; i < ni; i++){
    iadr = index[i];

    p0.xacc = accjerk[i].xacc;
    p0.yacc = accjerk[i].yacc;
    p0.zacc = accjerk[i].zacc;
    p0.xjrk = accjerk[i].xjrk;
    p0.yjrk = accjerk[i].yjrk;
    p0.zjrk = accjerk[i].zjrk;
    p0.pot  = accjerk[i].pot;
    p0.in   = accjerk[i].in;

    p0.dt   = dt_startup(etas, particle[iadr].t, p0);

    particle[iadr].dt   = p0.dt;
    particle[iadr].xacc = p0.xacc;
    particle[iadr].yacc = p0.yacc;
    particle[iadr].zacc = p0.zacc;
    particle[iadr].xjrk = p0.xjrk;
    particle[iadr].yjrk = p0.yjrk;
    particle[iadr].zjrk = p0.zjrk;
    particle[iadr].pot  = p0.pot;
    particle[iadr].in   = p0.in;

  }
}

// ---------------------------------------------------------------------------

void correct_sgmotion(int ni, int *index, struct Particle *particle)
{
  int i, iadr;
  REAL s0[DIM], c0[DIM];
  struct Particle p0, p1;
  
  for(i = 0; i < ni; i++){
    iadr = index[i];

    p0.t    = particle[iadr].t + particle[iadr].dt;
    p0.dt   = particle[iadr].dt;
    p0.xpos = particle[iadr].xpos;
    p0.ypos = particle[iadr].ypos;
    p0.zpos = particle[iadr].zpos;
    p0.xvel = particle[iadr].xvel;
    p0.yvel = particle[iadr].yvel;
    p0.zvel = particle[iadr].zvel;
    p0.xacc = particle[iadr].xacc;
    p0.yacc = particle[iadr].yacc;
    p0.zacc = particle[iadr].zacc;
    p0.xjrk = particle[iadr].xjrk;
    p0.yjrk = particle[iadr].yjrk;
    p0.zjrk = particle[iadr].zjrk;

    p1.t    = p0.t;
    p1.xpos = posvel[i].xpos;
    p1.ypos = posvel[i].ypos;
    p1.zpos = posvel[i].zpos;
    p1.xvel = posvel[i].xvel;
    p1.yvel = posvel[i].yvel;
    p1.zvel = posvel[i].zvel;
    p1.xacc = accjerk[i].xacc;
    p1.yacc = accjerk[i].yacc;
    p1.zacc = accjerk[i].zacc;
    p1.xjrk = accjerk[i].xjrk;
    p1.yjrk = accjerk[i].yjrk;
    p1.zjrk = accjerk[i].zjrk;
    p1.pot  = accjerk[i].pot;
    p1.in   = accjerk[i].in;

    correct(&p0, &p1, s0, c0);

    p1.dt = dt_criterion(etac, &p0, &p1, s0, c0);

    particle[iadr].t    = p1.t;
    particle[iadr].dt   = p1.dt;
    particle[iadr].xpos = p1.xpos;
    particle[iadr].ypos = p1.ypos;
    particle[iadr].zpos = p1.zpos;
    particle[iadr].xvel = p1.xvel;
    particle[iadr].yvel = p1.yvel;
    particle[iadr].zvel = p1.zvel;
    particle[iadr].xacc = p1.xacc;
    particle[iadr].yacc = p1.yacc;
    particle[iadr].zacc = p1.zacc;
    particle[iadr].xjrk = p1.xjrk;
    particle[iadr].yjrk = p1.yjrk;
    particle[iadr].zjrk = p1.zjrk;
    particle[iadr].pot  = p1.pot;
    particle[iadr].in   = p1.in;
  }
}

// ---------------------------------------------------------------------------
// **************************** time integration *****************************
// ---------------------------------------------------------------------------

void integrate_exmotion_g6avx(int nj, double eps2, int ni, int *index, REAL t_gl, struct Particle *particle, struct Address *address, void (*corrector)(int, int *, struct Particle *))
{
  int i, ii, nii;
  int npipe;
  static int pidx[NPIPES];
  static double pos[NPIPES][3], vel[NPIPES][3], acc[NPIPES][3], jrk[NPIPES][3], pot[NPIPES], h2[NPIPES];
  
  npipe = g6_npipes();
  g6_set_ti(0, t_gl);
  predict_exmotion(t_gl, ni, index, particle, address);
  for(i = 0; i < ni; i += npipe){
    nii = npipe;
    if(i + npipe > ni) nii = ni - i;
    for(ii = 0; ii < nii; ii++){
      pidx[ii]   = (int)posvel[i+ii].id;
      pos[ii][0] = posvel[i+ii].xpos;
      pos[ii][1] = posvel[i+ii].ypos;
      pos[ii][2] = posvel[i+ii].zpos;
      vel[ii][0] = posvel[i+ii].xvel;
      vel[ii][1] = posvel[i+ii].yvel;
      vel[ii][2] = posvel[i+ii].zvel;
    }
    g6calc_firsthalf(0, nj, nii, pidx, pos, vel, acc, jrk, pot, eps2, h2);
    g6calc_lasthalf(0, nj, nii, pidx, pos, vel, eps2, h2, acc, jrk, pot);
    for(ii = 0; ii < nii; ii++){
      accjerk[i+ii].xacc = acc[ii][0];
      accjerk[i+ii].yacc = acc[ii][1];
      accjerk[i+ii].zacc = acc[ii][2];
      accjerk[i+ii].pot  = pot[ii];
      accjerk[i+ii].xjrk = jrk[ii][0];
      accjerk[i+ii].yjrk = jrk[ii][1];
      accjerk[i+ii].zjrk = jrk[ii][2];
    }
  }
  corrector(ni, index, particle);

  int iadr;
  REAL gpsnp[3], gpjrk[3], gpacc[3], gpvel[3], gppos[3];
  REAL over6, over2;
  over6 = 1.0 / 6.0;
  over2 = 0.5;
  for(i = 0; i < ni; i++){
    iadr = index[i];

    gpsnp[0] = gpsnp[1] = gpsnp[2] = 0.0;

    gpjrk[0] = particle[iadr].xjrk * over6;
    gpjrk[1] = particle[iadr].yjrk * over6;
    gpjrk[2] = particle[iadr].zjrk * over6;

    gpacc[0] = particle[iadr].xacc * over2;
    gpacc[1] = particle[iadr].yacc * over2;
    gpacc[2] = particle[iadr].zacc * over2;

    gpvel[0] = particle[iadr].xvel;
    gpvel[1] = particle[iadr].yvel;
    gpvel[2] = particle[iadr].zvel;

    gppos[0] = particle[iadr].xpos;
    gppos[1] = particle[iadr].ypos;
    gppos[2] = particle[iadr].zpos;

    g6_set_j_particle(0, iadr, iadr, particle[iadr].t, particle[iadr].dt, particle[iadr].mass, gpsnp, gpjrk, gpacc, gpvel, gppos);
  }

  return;
}

void integrate_exmotion_g6avx2(int nj, double eps2, int ni, int *index, REAL t_gl, struct Particle *particle, struct Address *address, void (*corrector)(int, int *, struct Particle *))
{
  int i, ii, nii;
  int npipe;
  static int pidx[NPIPES], nnb[NPIPES];
  static double pos[NPIPES][3], vel[NPIPES][3], acc[NPIPES][3], jrk[NPIPES][3], pot[NPIPES], h2[NPIPES];
  
  npipe = g6_npipes();
  g6_set_ti(0, t_gl);
  predict_exmotion(t_gl, ni, index, particle, address);
  for(i = 0; i < ni; i += npipe){
    nii = npipe;
    if(i + npipe > ni) nii = ni - i;
    for(ii = 0; ii < nii; ii++){
      pidx[ii]   = (int)posvel[i+ii].id;
      pos[ii][0] = posvel[i+ii].xpos;
      pos[ii][1] = posvel[i+ii].ypos;
      pos[ii][2] = posvel[i+ii].zpos;
      vel[ii][0] = posvel[i+ii].xvel;
      vel[ii][1] = posvel[i+ii].yvel;
      vel[ii][2] = posvel[i+ii].zvel;
    }
    g6calc_firsthalf(0, nj, nii, pidx, pos, vel, acc, jrk, pot, eps2, h2);
    g6calc_lasthalf2(0, nj, nii, pidx, pos, vel, eps2, h2, acc, jrk, pot, nnb);
    for(ii = 0; ii < nii; ii++){
      accjerk[i+ii].xacc = acc[ii][0];
      accjerk[i+ii].yacc = acc[ii][1];
      accjerk[i+ii].zacc = acc[ii][2];
      accjerk[i+ii].pot  = pot[ii];
      accjerk[i+ii].xjrk = jrk[ii][0];
      accjerk[i+ii].yjrk = jrk[ii][1];
      accjerk[i+ii].zjrk = jrk[ii][2];
      accjerk[i+ii].in   = nnb[ii];
    }
  }
  corrector(ni, index, particle);

  int iadr;
  REAL gpsnp[3], gpjrk[3], gpacc[3], gpvel[3], gppos[3];
  REAL over6, over2;
  over6 = 1.0 / 6.0;
  over2 = 0.5;
  for(i = 0; i < ni; i++){
    iadr = index[i];

    gpsnp[0] = gpsnp[1] = gpsnp[2] = 0.0;

    gpjrk[0] = particle[iadr].xjrk * over6;
    gpjrk[1] = particle[iadr].yjrk * over6;
    gpjrk[2] = particle[iadr].zjrk * over6;

    gpacc[0] = particle[iadr].xacc * over2;
    gpacc[1] = particle[iadr].yacc * over2;
    gpacc[2] = particle[iadr].zacc * over2;

    gpvel[0] = particle[iadr].xvel;
    gpvel[1] = particle[iadr].yvel;
    gpvel[2] = particle[iadr].zvel;

    gppos[0] = particle[iadr].xpos;
    gppos[1] = particle[iadr].ypos;
    gppos[2] = particle[iadr].zpos;

    g6_set_j_particle(0, iadr, iadr, particle[iadr].t, particle[iadr].dt, particle[iadr].mass, gpsnp, gpjrk, gpacc, gpvel, gppos);
  }

  return;
}

void integrate_exmotion_g6avxn(int nj, double eps2, int ni, int *index, REAL t_gl, struct Particle *particle, struct Address *address, void (*corrector)(int, int *, struct Particle *))
{
  int i, ii, nii;
  int npipe, dummy;
  static int pidx[NPIPES];
  static double pos[NPIPES][3], vel[NPIPES][3], acc[NPIPES][3], jrk[NPIPES][3], pot[NPIPES], h2[NPIPES];
  
  npipe = g6_npipes();
  g6_set_ti(0, t_gl);
  predict_exmotion(t_gl, ni, index, particle, address);
  for(i = 0; i < ni; i += npipe){
    nii = npipe;
    if(i + npipe > ni) nii = ni - i;
    for(ii = 0; ii < nii; ii++){
      pidx[ii]   = (int)posvel[i+ii].id;
      pos[ii][0] = posvel[i+ii].xpos;
      pos[ii][1] = posvel[i+ii].ypos;
      pos[ii][2] = posvel[i+ii].zpos;
      vel[ii][0] = posvel[i+ii].xvel;
      vel[ii][1] = posvel[i+ii].yvel;
      vel[ii][2] = posvel[i+ii].zvel;
      h2[ii]     = 0.02;
    }
    g6calc_firsthalf(0, nj, nii, pidx, pos, vel, acc, jrk, pot, eps2, h2);
    g6calc_lasthalfn(0, nj, nii, pidx, pos, vel, eps2, h2, acc, jrk, pot);
    dummy = g6_read_neighbour_list(0);
    for(ii = 0; ii < nii; ii++){
      accjerk[i+ii].xacc = acc[ii][0];
      accjerk[i+ii].yacc = acc[ii][1];
      accjerk[i+ii].zacc = acc[ii][2];
      accjerk[i+ii].pot  = pot[ii];
      accjerk[i+ii].xjrk = jrk[ii][0];
      accjerk[i+ii].yjrk = jrk[ii][1];
      accjerk[i+ii].zjrk = jrk[ii][2];
      g6_get_neighbour_list(0, ii, 1024, &accjerk[i+ii].nn, accjerk[i+ii].nptr);
    }
  }
  corrector(ni, index, particle);

  int iadr;
  int j;
  for(i = 0; i < ni; i++){
    iadr = index[i];
    particle[iadr].nn = accjerk[i].nn;
    for(j = 0; j < accjerk[i].nn; j++)
      particle[iadr].nptr[j] = accjerk[i].nptr[j];
  }

  REAL gpsnp[3], gpjrk[3], gpacc[3], gpvel[3], gppos[3];
  REAL over6, over2;
  over6 = 1.0 / 6.0;
  over2 = 0.5;
  for(i = 0; i < ni; i++){
    iadr = index[i];

    gpsnp[0] = gpsnp[1] = gpsnp[2] = 0.0;

    gpjrk[0] = particle[iadr].xjrk * over6;
    gpjrk[1] = particle[iadr].yjrk * over6;
    gpjrk[2] = particle[iadr].zjrk * over6;

    gpacc[0] = particle[iadr].xacc * over2;
    gpacc[1] = particle[iadr].yacc * over2;
    gpacc[2] = particle[iadr].zacc * over2;

    gpvel[0] = particle[iadr].xvel;
    gpvel[1] = particle[iadr].yvel;
    gpvel[2] = particle[iadr].zvel;

    gppos[0] = particle[iadr].xpos;
    gppos[1] = particle[iadr].ypos;
    gppos[2] = particle[iadr].zpos;

    g6_set_j_particle(0, iadr, iadr, particle[iadr].t, particle[iadr].dt, particle[iadr].mass, gpsnp, gpjrk, gpacc, gpvel, gppos);
  }

  return;
}

void integrate_exmotion_g6avx2n(int nj, double eps2, int ni, int *index, REAL t_gl, struct Particle *particle, struct Address *address, void (*corrector)(int, int *, struct Particle *))
{
  int i, ii, nii;
  int npipe, dummy;
  static int pidx[NPIPES], nnb[NPIPES];
  static double pos[NPIPES][3], vel[NPIPES][3], acc[NPIPES][3], jrk[NPIPES][3], pot[NPIPES], h2[NPIPES];
  
  npipe = g6_npipes();
  g6_set_ti(0, t_gl);
  predict_exmotion(t_gl, ni, index, particle, address);
  for(i = 0; i < ni; i += npipe){
    nii = npipe;
    if(i + npipe > ni) nii = ni - i;
    for(ii = 0; ii < nii; ii++){
      pidx[ii]   = (int)posvel[i+ii].id;
      pos[ii][0] = posvel[i+ii].xpos;
      pos[ii][1] = posvel[i+ii].ypos;
      pos[ii][2] = posvel[i+ii].zpos;
      vel[ii][0] = posvel[i+ii].xvel;
      vel[ii][1] = posvel[i+ii].yvel;
      vel[ii][2] = posvel[i+ii].zvel;
      h2[ii]     = 0.02;
    }
    g6calc_firsthalf(0, nj, nii, pidx, pos, vel, acc, jrk, pot, eps2, h2);
    g6calc_lasthalf2n(0, nj, nii, pidx, pos, vel, eps2, h2, acc, jrk, pot, nnb);
    dummy = g6_read_neighbour_list(0);
    for(ii = 0; ii < nii; ii++){
      accjerk[i+ii].xacc = acc[ii][0];
      accjerk[i+ii].yacc = acc[ii][1];
      accjerk[i+ii].zacc = acc[ii][2];
      accjerk[i+ii].pot  = pot[ii];
      accjerk[i+ii].xjrk = jrk[ii][0];
      accjerk[i+ii].yjrk = jrk[ii][1];
      accjerk[i+ii].zjrk = jrk[ii][2];
      accjerk[i+ii].in   = nnb[ii];
      g6_get_neighbour_list(0, ii, 1024, &accjerk[i+ii].nn, accjerk[i+ii].nptr);
    }
  }
  corrector(ni, index, particle);

  int iadr;
  int j;
  for(i = 0; i < ni; i++){
    iadr = index[i];
    particle[iadr].nn = accjerk[i].nn;
    for(j = 0; j < accjerk[i].nn; j++)
      particle[iadr].nptr[j] = accjerk[i].nptr[j];
  }

  REAL gpsnp[3], gpjrk[3], gpacc[3], gpvel[3], gppos[3];
  REAL over6, over2;
  over6 = 1.0 / 6.0;
  over2 = 0.5;
  for(i = 0; i < ni; i++){
    iadr = index[i];

    gpsnp[0] = gpsnp[1] = gpsnp[2] = 0.0;

    gpjrk[0] = particle[iadr].xjrk * over6;
    gpjrk[1] = particle[iadr].yjrk * over6;
    gpjrk[2] = particle[iadr].zjrk * over6;

    gpacc[0] = particle[iadr].xacc * over2;
    gpacc[1] = particle[iadr].yacc * over2;
    gpacc[2] = particle[iadr].zacc * over2;

    gpvel[0] = particle[iadr].xvel;
    gpvel[1] = particle[iadr].yvel;
    gpvel[2] = particle[iadr].zvel;

    gppos[0] = particle[iadr].xpos;
    gppos[1] = particle[iadr].ypos;
    gppos[2] = particle[iadr].zpos;

    g6_set_j_particle(0, iadr, iadr, particle[iadr].t, particle[iadr].dt, particle[iadr].mass, gpsnp, gpjrk, gpacc, gpvel, gppos);
  }

  return;
}
