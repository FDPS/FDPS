#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "banana01.h"
#include "particle.h"
#include "globaltime.h"
#include "timestep.h"

#define BIT 51

static REAL dtmax, dtmin;

static REAL calc_modulus(REAL, REAL, REAL);
static void sorting(int, int, int *, struct Particle *);
static void exchange(int, int, int []);

void initialize_timestep(REAL dtmaxtmp)
{
  dtmax = dtmaxtmp;
  dtmin = dtmax * 1.0 / pow(2.0, (REAL)(BIT));
  return;
}

REAL dt_startup(REAL eta_s, REAL t_gl, struct Particle particle)
{
  int k;
  REAL dtcrit, dtblock;
  REAL acc, jrk;

  acc = particle.xacc * particle.xacc + particle.yacc * particle.yacc + particle.zacc * particle.zacc;
  jrk = particle.xjrk * particle.xjrk + particle.yjrk * particle.yjrk + particle.zjrk * particle.zjrk;
  dtcrit = eta_s * sqrt(acc / jrk);
  dtblock = dtmax;
  while(dtblock > dtcrit && dtblock > dtmin)
    dtblock = dtblock * 0.5;
  while(calc_modulus(t_gl, dtblock, dtmax) != 0.0){
    dtblock = dtblock * 0.5;
    if(dtblock < dtmin){
      fprintf(stderr, "Err: timestep criterion dtcrit: %e time: %e dtcand: %e\n", dtcrit, t_gl, dtblock);
      exit(1);
    }
  }
#ifdef TWOBODYTEST
  if(acc == 0.0 && jrk == 0.0)
    dtblock = 0.125;
#endif
  return dtblock;
}

REAL dt_criterion(REAL eta, struct Particle *pi0, struct Particle *pi1, REAL *s0, REAL *c0)
{
  int k;
  REAL dti, dts, ts;
  REAL s1[DIM];
  REAL acc, jrk, snp, crc;

  ts  = pi0->t;
  dts = pi0->dt;
  for(k = 0; k < DIM; k++)
    s1[k] = s0[k] + c0[k] * dts;
  acc = pi1->xacc * pi1->xacc + pi1->yacc * pi1->yacc + pi1->zacc * pi1->zacc;
  jrk = pi1->xjrk * pi1->xjrk + pi1->yjrk * pi1->yjrk + pi1->zjrk * pi1->zjrk;
  snp = s1[0] * s1[0] + s1[1] * s1[1] + s1[2] * s1[2];
  crc = c0[0] * c0[0] + c0[1] * c0[1] + c0[2] * c0[2];

  acc = sqrt(acc);
  jrk = sqrt(jrk);
  snp = sqrt(snp);
  crc = sqrt(crc);
  dti = eta * sqrt((acc * snp + jrk * jrk) / (jrk * crc + snp * snp));
  if(dti < dts){
    while(dts > dti && dts > dtmin)
      dts = dts * 0.5;
  }else if(dti > 2.0 * dts)
    if(calc_modulus(ts, 2.0 * dts, dtmax) == 0.0)
      dts = 2.0 * dts;
  return dts;
}

void global_time(REAL *t_gl, int n, int *ni, int *index, struct Particle *particle)
{
  REAL nextt;

  *t_gl = 1e30;
  sorting(0, (*ni)-1, index, particle);

  *t_gl = particle[index[0]].t + particle[index[0]].dt;

  for((*ni) = 0, nextt = particle[index[(*ni)]].t + particle[index[(*ni)]].dt; *t_gl == nextt && (*ni) < n; (*ni)++, nextt = particle[index[(*ni)]].t + particle[index[(*ni)]].dt)
    ;
  return;
}

void distribute_index(struct Actives *index, struct Particle *particle, struct Address *address)
{
  int i, iadr;

  index->nisg = 0;
  index->nicm = 0;
  index->nirl = 0;  
  for(i = 0; i < index->ni; i++){
    iadr = index->index[i];
    index->indexsg[index->nisg] = iadr;
    ++(index->nisg);
  }
  return;
}

void change_actives(int insertadr, struct Actives *index)
{
  int i, flag, iadr, tmpadr;

  for(flag = 0, iadr = insertadr, i = index->ni; ; i++){
    if(index->index[i] == insertadr)
      flag = 1;
    tmpadr          = index->index[i];
    index->index[i] = iadr;
    iadr            = tmpadr;
    if(flag == 1)
      break;
  }
  ++(index->ni);
  return;
}

void reset_time(int n, struct Particle *particle, struct Globaltime *tim)
{
  int i;

  tim->t_gl = 0.0;
  (tim->t_mod)++;
  for(i = 0; i < n; i++)
    particle[i].t = 0.0;
}

static REAL calc_modulus(REAL num, REAL mod, REAL modmax)
{
  REAL tmpmod;
  
  tmpmod = modmax;
  while(num > 0.0){
    while(tmpmod > num)
      tmpmod *= 0.5;
    if(tmpmod < mod)
      break;
    num -= tmpmod;
  }
  return num;
}

static void sorting(int left, int right, int *index, struct Particle *particle)
{
  int i, j, axis_index;
  REAL axis;

  axis_index = (left + right) / 2;
  axis       = particle[index[axis_index]].dt;
  i = left;
  j = right;
  for( ; ; ){
    while(particle[index[i]].dt < axis)
      i++;
    while(particle[index[j]].dt > axis)
      j--;
    if(i >= j) break;
    exchange(i, j, index);
    i++;
    j--;
  }
  if(left < i - 1)
    sorting(left, i-1, index, particle);
  if(j + 1 < right)
    sorting(j+1, right, index, particle);
}

static void exchange(int i, int j, int index[])
{
  int tmpi;
  
  tmpi = index[i];
  index[i] = index[j];
  index[j] = tmpi;
}
