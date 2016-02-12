#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "banana01.h"
#include "particle.h"
#include "energy.h"
#include "step.h"
#include "globaltime.h"
#include "io.h"
#include "calc_energy.h"
#include "hermite.h"
#include "timestep.h"
#include "timeprof.h"
#include "gp6util.h"

void initialize_iparticle(struct Particle **);
void calc_total_energy(int, struct Particle *, struct Energy *);
void set_particle_on_grape(int, struct Particle *);
void initialtest(int, int *, struct Particle *);
void neighbourlisttest(int, int *, int, double, struct Particle *);

int main(int argc, char **argv)
{
  /* star data */
  int n;
  REAL epsinv, eps2;
  static struct Particle particle[NMAX];
  static struct Address  address[NMAX];
  /* timestep */
  struct Globaltime tim;
  REAL dtstar, dtstarinv;
  REAL eta, eta_s;
  static struct Actives index;
  static struct Step step;
  /* energy */
  static struct Energy energies;
  /* io */
  FILE *fp;
  int first;
  REAL dtbin, dtbininv, tend;
  char outputfile[1000], inputfile[1000];
  /* tmp */
  int i;
  int fi;

  fp = fopen("input.list", "r");
  fi = fscanf(fp, "%d%lf%lf", &first, &tend, &dtbin);
  fi = fscanf(fp, "%lf", &epsinv);
  fi = fscanf(fp, "%lf%lf%lf", &eta, &eta_s, &dtstarinv);
  fi = fscanf(fp, "%s%s", inputfile, outputfile);
  fclose(fp);

  dtstar   = 1.0 / dtstarinv;
  dtbininv = 1.0 / dtbin;
  if(epsinv == 0.0)
    eps2 = 0.0;
  else
    eps2     = 1.0 / (epsinv * epsinv);
  initialize_accuracy(eta_s, eta);
  initialize_timestep(dtstar);

  g6_open(0);
  g6_set_tunit(45);
  g6_set_xunit(45);

  if(first == 0){
    fp = fopen(inputfile, "r");
    input_startup(fp, dtstar, &n, particle, address);
    fclose(fp);
    tim.t_gl  = 0.0;
    tim.t_mod = 0;
    index.ni = n;
    for(i = 0; i < index.ni; i++)
      index.index[i] = i;
    distribute_index(&index, particle, address);
    step.tb = get_dtime();
    set_particle_on_grape(n, particle);
#ifdef NORMAL
    integrate_exmotion_g6avx(n, eps2, index.nisg, index.indexsg, tim.t_gl, particle, address, &startup_sgmotion);
#else
#ifdef ONLYNNB
    integrate_exmotion_g6avx2(n, eps2, index.nisg, index.indexsg, tim.t_gl, particle, address, &startup_sgmotion);
#else
#ifdef ONLYNEIGHBOUR
    integrate_exmotion_g6avxn(n, eps2, index.nisg, index.indexsg, tim.t_gl, particle, address, &startup_sgmotion);
#else
    integrate_exmotion_g6avx2n(n, eps2, index.nisg, index.indexsg, tim.t_gl, particle, address, &startup_sgmotion);
#endif
#endif
#endif
    step.step += index.ni;
    step.blck += 1;
    calc_total_energy(n, particle, &energies);
    energies.ei = energies.e;
    energies.ep = energies.e;
  }else{
    index.ni = n;
    for(i = 0; i < index.ni; i++)
      index.index[i] = i;
    fprintf(stderr, "error in a variable \"first\".\n");
    exit(0);
  }

#ifdef INITIALTEST
  //  neighbourlisttest(index.nisg, index.indexsg, n, eps2, particle);
  initialtest(index.nisg, index.indexsg, particle);
  exit(0);
#endif

#ifdef ACCURACYTEST
  REAL t1, t2;
  REAL de_past = 0.0;
  t1 = get_dtime();
  tend = 7.5;
#else
  output_parameter(stdout, inputfile, epsinv, tend);
  output_data(stdout, n, tim.t_mod * dtstar + tim.t_gl, &energies, &step);
#endif

  while(tim.t_mod * dtstar + tim.t_gl < tend){
    global_time(&tim.t_gl, n, &index.ni, index.index, particle);
    distribute_index(&index, particle, address);

#ifdef NORMAL
    integrate_exmotion_g6avx(n, eps2, index.nisg, index.indexsg, tim.t_gl, particle, address, &correct_sgmotion);
#else
#ifdef ONLYNNB
    integrate_exmotion_g6avx2(n, eps2, index.nisg, index.indexsg, tim.t_gl, particle, address, &correct_sgmotion);
#else
#ifdef ONLYNEIGHBOUR
    integrate_exmotion_g6avxn(n, eps2, index.nisg, index.indexsg, tim.t_gl, particle, address, &correct_sgmotion);
#else
    integrate_exmotion_g6avx2n(n, eps2, index.nisg, index.indexsg, tim.t_gl, particle, address, &correct_sgmotion);
#endif
#endif
#endif
    step.step += index.ni;
    step.blck += 1;

    if(tim.t_gl - (int)(tim.t_gl * dtstarinv) * dtstar == 0.0){
      calc_total_energy(n, particle, &energies);
#ifdef ACCURACYTEST
      if(tim.t_mod * dtstar + tim.t_gl - (int)((tim.t_mod * dtstar + tim.t_gl) / 0.375) * 0.375 == 0.0){
	de_past    += fabs((energies.e - energies.ep) / energies.ep);
	energies.ep = energies.e;
      }
#else
      output_data(stdout, n, tim.t_mod * dtstar + tim.t_gl, &energies, &step);
#endif

      if(fabs((energies.e - energies.ei) / (- energies.ei)) > 1e-3){
	fprintf(stderr, "Energy error is more than 0.1 %% !\n");
	exit(EXIT_SUCCESS);
      }

      reset_time(n, particle, &tim);

      set_particle_on_grape(n, particle);
    }

    if(tim.t_mod * dtstar + tim.t_gl - (int)((tim.t_mod * dtstar + tim.t_gl) * dtbininv) * dtbin == 0){
      char outsnap[1000];
      sprintf(outsnap, "%s.bin.%04d", outputfile, (int)(tim.t_mod * dtstar + tim.t_gl));
      fp = fopen(outsnap, "wb");
      output_snapshot(fp, &tim, n, particle, &energies);
      fclose(fp);
    }

  }

#ifdef ACCURACYTEST
  t2 = get_dtime();
  REAL tnow = tim.t_mod * dtstar + tim.t_gl;
  printf("%d %e %e %.5f %.6f %e %e\n", n, step.step / (double)(n) / (tnow / (2.0 * sqrt(2.0))), de_past / (tnow / 0.375), eta, eta_s, tend, t2 - t1);
#endif

  g6_close(0);

  exit(EXIT_SUCCESS);
}

void calc_total_energy(int n, struct Particle *particle, struct Energy *energies)
{
  int i, k;

  energies->ek = calc_total_kinetic_energy(n, particle);
  energies->ew = calc_total_potential_energy(n, particle);
  energies->e  = energies->ek + energies->ew;

  return;
}

void set_particle_on_grape(int n, struct Particle *particle)
{
  int i;
  double gpsnp[3], gpjrk[3], gpacc[3], gpvel[3], gppos[3];
  double over2, over6;

  over2 = 0.5;
  over6 = 1.0 / 6.0;
  for(i = 0; i < n; i++){
    gpsnp[0] = 0.0;
    gpsnp[1] = 0.0;
    gpsnp[2] = 0.0;
    
    gpjrk[0] = particle[i].xjrk * over6;
    gpjrk[1] = particle[i].yjrk * over6;
    gpjrk[2] = particle[i].zjrk * over6;

    gpacc[0] = particle[i].xacc * over2;
    gpacc[1] = particle[i].yacc * over2;
    gpacc[2] = particle[i].zacc * over2;

    gpvel[0] = particle[i].xvel;
    gpvel[1] = particle[i].yvel;
    gpvel[2] = particle[i].zvel;

    gppos[0] = particle[i].xpos;
    gppos[1] = particle[i].ypos;
    gppos[2] = particle[i].zpos;

    g6_set_j_particle(0, i, i, particle[i].t, particle[i].dt, particle[i].mass, gpsnp, gpjrk, gpacc, gpvel, gppos);
  }

  return;
}

void initialtest(int ni, int *index, struct Particle *particle)
{
  int i, j, iadr;
  REAL acc, jrk;
  for(i = 0; i < ni; i++){
    iadr = index[i];
    acc  = particle[iadr].xacc * particle[iadr].xacc;
    acc += particle[iadr].yacc * particle[iadr].yacc;
    acc += particle[iadr].zacc * particle[iadr].zacc;
    jrk  = particle[iadr].xjrk * particle[iadr].xjrk;
    jrk += particle[iadr].yjrk * particle[iadr].yjrk;
    jrk += particle[iadr].zjrk * particle[iadr].zjrk;
    printf("%4d %+.16e %+.16e %+.16e %4d %4d\n", iadr, sqrt(acc), sqrt(jrk), particle[iadr].pot, particle[iadr].in, particle[iadr].nn);
#if 0
    printf("%4d %4d", iadr, particle[iadr].nn);
    for(j = 0; j < particle[iadr].nn; j++)
      printf(" %4d", particle[iadr].nptr[j]);
    printf("\n");
#endif
  }
}

void neighbourlisttest(int ni, int *index, int n, double eps2, struct Particle *particle)
{
  int i, j, iadr;
  double dx, dy, dz, r2;
  
  for(i = 0; i < ni; i++){
    iadr = index[i];
    particle[iadr].nn = 0;
    for(j = 0; j < n; j++){
      if(iadr == j)
	continue;
      dx = particle[j].xpos - particle[iadr].xpos;
      dy = particle[j].ypos - particle[iadr].ypos;
      dz = particle[j].zpos - particle[iadr].zpos;
      r2 = dx * dx + dy * dy + dz * dz + eps2;
      if(r2 < 0.02){
	particle[iadr].nptr[particle[iadr].nn] = particle[j].id;
	++(particle[iadr].nn);
      }
    }
  }
  return;
}
