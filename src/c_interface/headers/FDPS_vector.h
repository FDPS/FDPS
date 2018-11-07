#pragma once
#include "FDPS_basic.h"

//**** PS::F32vec
typedef struct  {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
   fdps_f32 x,y;
#else
   fdps_f32 x,y,z;
#endif
} fdps_f32vec;

//**** PS::F64vec
typedef struct  {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
   fdps_f64 x,y;
#else
   fdps_f64 x,y,z;
#endif
} fdps_f64vec;
