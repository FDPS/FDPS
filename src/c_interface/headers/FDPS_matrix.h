#pragma once
#include "FDPS_basic.h"

//**** PS::F32mat
typedef struct {
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
   fdps_f32 xx,yy,zz,xy,xz,yz;
#else
   fdps_f32 xx,yy,xy;
#endif
} fdps_f32mat;

//**** PS::F64mat
typedef struct {
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
   fdps_f64 xx,yy,zz,xy,xz,yz;
#else
   fdps_f64 xx,yy,xy;
#endif
} fdps_f64mat;
