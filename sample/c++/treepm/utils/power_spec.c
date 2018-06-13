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
void calc_power(float *, struct run_param*, int);

#define RHO(ix,iy,iz) mesh[(iz)+nmesh_p2*((iy)+nmesh*(ix))]

void calc_mesh_density(struct particle *ptcl, float *mesh, 
		       int nmesh, struct run_param *this_run)
{
  int nmesh_p2, nmesh_total;

  nmesh_p2 = nmesh+2;
  nmesh_total = nmesh*nmesh*nmesh_p2;

  /* zero out */
  for(int i=0;i<nmesh_total;i++) mesh[i] = 0.0;

  for(int p=0;p<this_run->npart_total;p++) {
    float xpos, ypos, zpos;  /* normalized coordinate 0.0 < {x,y,z}pos < 1.0 */

    float xt1, dx1;
    float wi11,wi21,wi31,wi12,wi22,wi32,wi13,wi23,wi33;
    int   iw11,iw21,iw31,iw12,iw22,iw32,iw13,iw23,iw33;

    xpos = ptcl[p].xpos;
    ypos = ptcl[p].ypos;
    zpos = ptcl[p].zpos;

    xt1  = xpos*(float)nmesh-0.5;
    iw21 = (int)(xt1 + 0.5);
    dx1  = xt1 - (float)iw21;
    wi11 = 0.5*(0.5-dx1)*(0.5-dx1);
    wi21 = 0.75-dx1*dx1;
    wi31 = 0.5*(0.5+dx1)*(0.5+dx1);
    iw11 = iw21-1;
    iw31 = iw21+1;
    if(iw21 == 0){
      iw11 = nmesh-1;
    }else if(iw21 == nmesh-1){
      iw31 = 0;
    }else if(iw21 == nmesh){
      iw21 = 0;
      iw31 = 1;
    }
    
    xt1  = ypos*(float)nmesh-0.5;
    iw22 = (int)(xt1 + 0.5);
    dx1  = xt1 - (float)iw22;
    wi12 = 0.5*(0.5-dx1)*(0.5-dx1);
    wi22 = 0.75-dx1*dx1;
    wi32 = 0.5*(0.5+dx1)*(0.5+dx1);
    iw12 = iw22-1;
    iw32 = iw22+1;
    if(iw22 == 0){
      iw12 = nmesh-1;
    }else if(iw22 == nmesh-1){
      iw32 = 0;
    }else if(iw22 == nmesh){
      iw22 = 0;
      iw32 = 1;
    }

    xt1  = zpos*(float)nmesh-0.5;
    iw23 = (int)(xt1 + 0.5);
    dx1  = xt1 - (float)iw23;
    wi13 = 0.5*(0.5-dx1)*(0.5-dx1);
    wi23 = 0.75-dx1*dx1;
    wi33 = 0.5*(0.5+dx1)*(0.5+dx1);
    iw13 = iw23-1;
    iw33 = iw23+1;
    if(iw23 == 0){
      iw13 = nmesh-1;
    }else if(iw23 == nmesh-1){
      iw33 = 0;
    }else if(iw23 == nmesh){
      iw23 = 0;
      iw33 = 1;
    }

    wi11 *= ptcl[p].mass;
    wi21 *= ptcl[p].mass;
    wi31 *= ptcl[p].mass;

    RHO(iw11,iw12,iw13)=RHO(iw11,iw12,iw13)+wi11*wi12*wi13;
    RHO(iw21,iw12,iw13)=RHO(iw21,iw12,iw13)+wi21*wi12*wi13;
    RHO(iw31,iw12,iw13)=RHO(iw31,iw12,iw13)+wi31*wi12*wi13;
    RHO(iw11,iw22,iw13)=RHO(iw11,iw22,iw13)+wi11*wi22*wi13;
    RHO(iw21,iw22,iw13)=RHO(iw21,iw22,iw13)+wi21*wi22*wi13;
    RHO(iw31,iw22,iw13)=RHO(iw31,iw22,iw13)+wi31*wi22*wi13;
    RHO(iw11,iw32,iw13)=RHO(iw11,iw32,iw13)+wi11*wi32*wi13;
    RHO(iw21,iw32,iw13)=RHO(iw21,iw32,iw13)+wi21*wi32*wi13;
    RHO(iw31,iw32,iw13)=RHO(iw31,iw32,iw13)+wi31*wi32*wi13;
    
    RHO(iw11,iw12,iw23)=RHO(iw11,iw12,iw23)+wi11*wi12*wi23;
    RHO(iw21,iw12,iw23)=RHO(iw21,iw12,iw23)+wi21*wi12*wi23;
    RHO(iw31,iw12,iw23)=RHO(iw31,iw12,iw23)+wi31*wi12*wi23;
    RHO(iw11,iw22,iw23)=RHO(iw11,iw22,iw23)+wi11*wi22*wi23;
    RHO(iw21,iw22,iw23)=RHO(iw21,iw22,iw23)+wi21*wi22*wi23;
    RHO(iw31,iw22,iw23)=RHO(iw31,iw22,iw23)+wi31*wi22*wi23;
    RHO(iw11,iw32,iw23)=RHO(iw11,iw32,iw23)+wi11*wi32*wi23;
    RHO(iw21,iw32,iw23)=RHO(iw21,iw32,iw23)+wi21*wi32*wi23;
    RHO(iw31,iw32,iw23)=RHO(iw31,iw32,iw23)+wi31*wi32*wi23;

    RHO(iw11,iw12,iw33)=RHO(iw11,iw12,iw33)+wi11*wi12*wi33;
    RHO(iw21,iw12,iw33)=RHO(iw21,iw12,iw33)+wi21*wi12*wi33;
    RHO(iw31,iw12,iw33)=RHO(iw31,iw12,iw33)+wi31*wi12*wi33;
    RHO(iw11,iw22,iw33)=RHO(iw11,iw22,iw33)+wi11*wi22*wi33;
    RHO(iw21,iw22,iw33)=RHO(iw21,iw22,iw33)+wi21*wi22*wi33;
    RHO(iw31,iw22,iw33)=RHO(iw31,iw22,iw33)+wi31*wi22*wi33;
    RHO(iw11,iw32,iw33)=RHO(iw11,iw32,iw33)+wi11*wi32*wi33;
    RHO(iw21,iw32,iw33)=RHO(iw21,iw32,iw33)+wi21*wi32*wi33;
    RHO(iw31,iw32,iw33)=RHO(iw31,iw32,iw33)+wi31*wi32*wi33;

  }

  return;
}

int main(int argc, char **argv) 
{

  struct particle *ptcl;
  struct run_param this_run;
  float *mesh;

  if(argc != 3) {
    fprintf(stderr, "Usage :: %s <input_prefix> <nmesh>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  input_header_file(&this_run, argv[1]);

  ptcl = (struct particle *)malloc(sizeof(struct particle)*this_run.npart_total);

  int nmesh = atoi(argv[2]);
  int nmesh_p2 = nmesh+2;

  mesh = (float *) malloc(sizeof(float)*nmesh*nmesh*nmesh_p2);
    
  input_data_all(ptcl, &this_run, argv[1]);

  calc_mesh_density(ptcl, mesh, nmesh, &this_run);
  calc_power(mesh, &this_run, nmesh);

  free(mesh);
  free(ptcl);
}

