#include <stdio.h>
#include <stdlib.h>

#include "cosmology.h"
#include "run_param.h"
#include "particle.h"

void input_header_file(struct run_param *this_run, char *prefix)
{

  FILE *input_fp;
  static char filename[128];

  sprintf(filename, "%s-0", prefix);

  input_fp = fopen(filename,"r");

  int cnt;

  cnt = 0;

  cnt += fread(&(this_run->mpi_nproc), sizeof(int), 1, input_fp);
  cnt += fread(&(this_run->mpi_rank), sizeof(int), 1, input_fp);
  cnt += fread(&(this_run->npart_local), sizeof(int64_t), 1, input_fp);
  cnt += fread(&(this_run->npart_total), sizeof(int64_t), 1, input_fp);
  cnt += fread(&(this_run->cosm.omega_m), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->cosm.omega_b), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->cosm.omega_v), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->cosm.omega_nu), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->cosm.hubble), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->anow), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->tnow), sizeof(COSM_REAL), 1, input_fp);

  cnt += fread(&(this_run->lunit), sizeof(double), 1, input_fp);
  cnt += fread(&(this_run->munit), sizeof(double), 1, input_fp);
  cnt += fread(&(this_run->tunit), sizeof(double), 1, input_fp);

  fclose(input_fp);
}

void input_header(struct run_param *this_run, FILE *input_fp)
{
  int cnt;

  cnt = 0;

  cnt += fread(&(this_run->mpi_nproc), sizeof(int), 1, input_fp);
  cnt += fread(&(this_run->mpi_rank), sizeof(int), 1, input_fp);
  cnt += fread(&(this_run->npart_local), sizeof(int64_t), 1, input_fp);
  cnt += fread(&(this_run->npart_total), sizeof(int64_t), 1, input_fp);
  cnt += fread(&(this_run->cosm.omega_m), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->cosm.omega_b), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->cosm.omega_v), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->cosm.omega_nu), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->cosm.hubble), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->anow), sizeof(COSM_REAL), 1, input_fp);
  cnt += fread(&(this_run->tnow), sizeof(COSM_REAL), 1, input_fp);

  cnt += fread(&(this_run->lunit), sizeof(double), 1, input_fp);
  cnt += fread(&(this_run->munit), sizeof(double), 1, input_fp);
  cnt += fread(&(this_run->tunit), sizeof(double), 1, input_fp);

}

void input_data(struct particle *ptcl, struct run_param *this_run, 
		char *filename)
{
  FILE *input_fp;

  input_fp = fopen(filename,"r");

  input_header(this_run, input_fp);

  for(int64_t i=0;i<this_run->npart_local;i++) {
    fread(&(ptcl[i].mass), sizeof(float), 1, input_fp);
    fread(&(ptcl[i].eps),  sizeof(float), 1, input_fp);
    fread(&(ptcl[i].xpos), sizeof(double), 1, input_fp);
    fread(&(ptcl[i].ypos), sizeof(double), 1, input_fp);
    fread(&(ptcl[i].zpos), sizeof(double), 1, input_fp);
    fread(&(ptcl[i].xvel), sizeof(double), 1, input_fp);
    fread(&(ptcl[i].yvel), sizeof(double), 1, input_fp);
    fread(&(ptcl[i].zvel), sizeof(double), 1, input_fp);
  }

  fclose(input_fp);
}

void input_data_all(struct particle *ptcl, struct run_param *this_run,
		    char *prefix)
{

  static char filename[128];
  int64_t npart_read;
  int nproc,npart_total;
  FILE *fp;

  npart_read = 0;

  /* Initially, read the header part of the first file */
  sprintf(filename,"%s-0", prefix);

  fp = fopen(filename,"r");
  input_header(this_run, fp);
  nproc = this_run->mpi_nproc;
  npart_total = this_run->npart_total;
  fclose(fp);

  for(int ip = 0;ip < nproc;ip++) {
    sprintf(filename,"%s-%d", prefix, ip);

    input_data(ptcl+npart_read, this_run, filename);

    npart_read += this_run->npart_local;

  }

  if(npart_read != this_run->npart_total) {
    fprintf(stderr, "Inconsistent Number of Particles !\n");
    exit(EXIT_FAILURE);
  }

  
}
