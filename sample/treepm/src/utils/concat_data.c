#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#ifdef __ISOLATED__
#error Isolated boundary condition is not appropriate for this exectuable.
#endif /* __ISOLATED__ */

#include "run_param.h"
#include "cosmology.h"
#include "particle.h"

void input_header_file(struct run_param*, char*);
void input_data_all(struct particle *, struct run_param*, char*);
void output_data_all(struct particle *, struct run_param*, char*);

int main(int argc, char **argv) 
{

  struct particle *ptcl;
  struct run_param this_run;
  float *mesh;

  if(argc != 3) {
    fprintf(stderr, "Usage :: %s <input_prefix> <output_data>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  input_header_file(&this_run, argv[1]);

  ptcl = (struct particle *)malloc(sizeof(struct particle)*this_run.npart_total);

  int nmesh = atoi(argv[2]);
  int nmesh_p2 = nmesh+2;

  mesh = (float *) malloc(sizeof(float)*nmesh*nmesh*nmesh_p2);
  
  input_data_all(ptcl, &this_run, argv[1]);

  output_data_all(ptcl, &this_run, argv[2]);

  free(ptcl);
}

