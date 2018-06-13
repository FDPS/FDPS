#ifndef __RUN_PARAM_H__
#define __RUN_PARAM_H__

#include <stdio.h>
#include <stdint.h>

#include "cosmology.h"

struct run_param {
  int step;
  int mpi_nproc, mpi_rank;

  int64_t npart_total;
  int64_t npart_local;

  COSM_REAL tnow;
  COSM_REAL anow;
  COSM_REAL znow;
  COSM_REAL hnow;

  COSM_REAL tend;

  double lunit,munit,tunit;

  struct cosmology cosm;

  char model_name[64];

  int noutput, output_indx;
  COSM_REAL *output_timing;

  FILE *diag_fp, *param_fp;
};

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) (SQR(x)*SQR(x))

#endif /* __RUN_PARAM_H__ */
