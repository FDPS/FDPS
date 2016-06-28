#include <stdio.h>
#include <stdlib.h>

#include "cosmology.h"
#include "run_param.h"
#include "particle.h"

void output_header(struct run_param *this_run, FILE *output_fp)
{
  int cnt;

  cnt = 0;

  cnt += fwrite(&(this_run->mpi_nproc), sizeof(int), 1, output_fp);
  cnt += fwrite(&(this_run->mpi_rank), sizeof(int), 1, output_fp);
  cnt += fwrite(&(this_run->npart_local), sizeof(int64_t), 1, output_fp);
  cnt += fwrite(&(this_run->npart_total), sizeof(int64_t), 1, output_fp);
  cnt += fwrite(&(this_run->cosm.omega_m), sizeof(COSM_REAL), 1, output_fp);
  cnt += fwrite(&(this_run->cosm.omega_b), sizeof(COSM_REAL), 1, output_fp);
  cnt += fwrite(&(this_run->cosm.omega_v), sizeof(COSM_REAL), 1, output_fp);
  cnt += fwrite(&(this_run->cosm.omega_nu), sizeof(COSM_REAL), 1, output_fp);
  cnt += fwrite(&(this_run->cosm.hubble), sizeof(COSM_REAL), 1, output_fp);
  cnt += fwrite(&(this_run->anow), sizeof(COSM_REAL), 1, output_fp);
  cnt += fwrite(&(this_run->tnow), sizeof(COSM_REAL), 1, output_fp);

  cnt += fwrite(&(this_run->lunit), sizeof(double), 1, output_fp);
  cnt += fwrite(&(this_run->munit), sizeof(double), 1, output_fp);
  cnt += fwrite(&(this_run->tunit), sizeof(double), 1, output_fp);

}

void output_data_all(struct particle *ptcl, struct run_param *this_run,
		     char *filename)
{
  FILE *fp;

  fp = fopen(filename, "w");

  /* This is for a concatenated file.*/
  this_run->mpi_rank = this_run->mpi_nproc;
  this_run->npart_local = this_run->npart_total;

  output_header(this_run, fp);

  for(int64_t i=0;i<this_run->npart_total;i++) {
    fwrite(&(ptcl[i].mass), sizeof(float), 1, fp);
    fwrite(&(ptcl[i].eps),  sizeof(float), 1, fp);
    fwrite(&(ptcl[i].xpos), sizeof(double), 1, fp);
    fwrite(&(ptcl[i].ypos), sizeof(double), 1, fp);
    fwrite(&(ptcl[i].zpos), sizeof(double), 1, fp);
    fwrite(&(ptcl[i].xvel), sizeof(double), 1, fp);
    fwrite(&(ptcl[i].yvel), sizeof(double), 1, fp);
    fwrite(&(ptcl[i].zvel), sizeof(double), 1, fp);
  }

  fclose(fp);
}
