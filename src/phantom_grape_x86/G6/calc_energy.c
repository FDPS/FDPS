#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "banana01.h"
#include "particle.h"
#include "calc_energy.h"

REAL calc_total_kinetic_energy(int n, struct Particle *particle)
{
  int i, k;
  REAL v2, ek;

  for(ek = 0.0, i = 0; i < n; i++){
    v2  = particle[i].xvel * particle[i].xvel;
    v2 += particle[i].yvel * particle[i].yvel;
    v2 += particle[i].zvel * particle[i].zvel;
    ek += particle[i].mass * v2;
  }
  return 0.5 * ek;
}

#ifndef ACCURACYTEST
REAL calc_total_potential_energy(int n, struct Particle *particle)
{
  int i;
  REAL ew;

  for(ew = 0.0, i = 0; i < n; i++)
    ew += particle[i].mass * particle[i].pot;
  
  return 0.5 * ew;
}
#else
REAL calc_total_potential_energy(int n, struct Particle *particle)
{
  int i, j;
  REAL dx, dy, dz, r2, rinv;
  struct Particle *iptr, *jptr;
  REAL ew = 0.0;
  REAL eps2 = 16.0 / ((REAL)n * (REAL)n);
  
  for(i = 0, iptr = particle; i < n; i++, iptr++){
    for(j = i + 1, jptr = iptr + 1; j < n; j++, jptr++){
      dx    = jptr->xpos - iptr->xpos;
      dy    = jptr->ypos - iptr->ypos;
      dz    = jptr->zpos - iptr->zpos;
      r2    = dx * dx + dy * dy + dz * dz + eps2;
      rinv  = 1.0 / sqrt(r2);
      rinv *= iptr->mass * jptr->mass;
      ew   -= rinv;
    }
  }
  return ew;
}
#endif
