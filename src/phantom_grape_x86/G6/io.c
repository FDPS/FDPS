#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "banana01.h"
#include "particle.h"
#include "energy.h"
#include "globaltime.h"
#include "step.h"
#include "io.h"
#include "timeprof.h"

void input_startup(FILE *fp, REAL dtmax, int *n, struct Particle *particle, struct Address *address)
{
  int i, k, fi;
  int ntmp;
  REAL ttmp;
  
  fi = fscanf(fp, "%lf", &ttmp);
  fi = fscanf(fp, "%d", n);
  
  for(i = 0; i < *n; i++){
    fi = fscanf(fp, "%d", &ntmp);
    fi = fscanf(fp, "%f", &(particle[i].mass));
    fi = fscanf(fp, "%lf%lf%lf", &(particle[i].xpos), &(particle[i].ypos), &(particle[i].zpos));
    fi = fscanf(fp, "%lf%lf%lf", &(particle[i].xvel), &(particle[i].yvel), &(particle[i].zvel));

    particle[i].t  = 0.0;
    particle[i].dt = dtmax;

    particle[i].id = i;
    particle[i].nn = 0;
    particle[i].nptr = (int*)malloc(sizeof(int) * 1024);

    address[i].jgp = 0;
    address[i].adr = i;
  }
  return;
}

void input_snapshot(FILE *fp, struct Globaltime *tim, int *n, struct Particle *particle, struct Energy *energies)
{
  int fi;

  fi = fread(tim, sizeof(struct Globaltime), 1, fp);
  fi = fread(n, sizeof(int), 1, fp);
  fi = fread(particle, sizeof(struct Particle), (*n), fp);
  fi = fread(energies, sizeof(struct Energy), 1, fp);
  return;
}

void output_snapshot(FILE *fp, struct Globaltime *tim, int n, struct Particle *particle, struct Energy *energies)
{
  fwrite(tim, sizeof(struct Globaltime), 1, fp);
  fwrite(&n, sizeof(int), 1, fp);
  fwrite(particle, sizeof(struct Particle), n, fp);
  fwrite(energies, sizeof(struct Energy), 1, fp);
  return;
}

void output_parameter(FILE *fp, char *inputfile, REAL epsinv, REAL tend)
{
  fprintf(fp, "initial: %s\n", inputfile);
  if(epsinv == 0.0)
    fprintf(fp, "eps: %e\n", 0.0);
  else
    fprintf(fp, "eps: %e\n", 1.0 / epsinv);
  fprintf(fp, "tend: %.f\n", tend);
  fprintf(fp, "\n");
  return;
}

void output_data(FILE *fp, int n, REAL tim, struct Energy *energies, struct Step *step)
{
  step->tf = get_dtime();
  fprintf(stderr, "Time: %9.3f\n", tim);
  fprintf(fp, "Time: %9.3f\n", tim);
  fprintf(fp, "---Energy----------------------------------------------------------------------\n");
  fprintf(fp, "tot: %.8f ", energies->e);
  fprintf(fp, "w: %.8f k: %.8f ", energies->ew, energies->ek);
  fprintf(fp, "de: %+.3e de/|e|: %+.3e\n", energies->e - energies->ei, (energies->e - energies->ei) / (- energies->ei));
  fprintf(fp, "---Time------------------------------------------------------------------------\n");
  fprintf(fp, "tcalc: %.2e", step->tf - step->tb);
  fprintf(fp, " %.2e Gflops", (REAL)step->step * (REAL)n * OPERATION / (step->tf - step->tb) * 1e-9);
#if 0
  fprintf(fp, " niave: %.2e\n", (REAL)step->step / (REAL)step->blck);
#else
  fprintf(fp, " niave: %.2e", (REAL)step->step / (REAL)step->blck);
  fprintf(fp, " %e\n", (REAL)step->step);
#endif
  fprintf(fp, "\n");
  fflush(fp);
  step->tb   = get_dtime();
  step->step = 0;
  step->blck = 0;
  return;
}
