#pragma once
#include "FDPS_basic.h"
#include "FDPS_vector.h"
#include "FDPS_matrix.h"

#ifdef PARTICLE_SIMULATOR_SPMOM_F32
typedef fdps_s32    fdps_sSP;
typedef fdps_f32    fdps_fSP;
typedef fdps_f32vec fdps_fSPvec;
typedef fdps_f32mat fdps_fSPmat;
#else
typedef fdps_s64    fdps_sSP;
typedef fdps_f64    fdps_fSP;
typedef fdps_f64vec fdps_fSPvec;
typedef fdps_f64mat fdps_fSPmat;
#endif

//**** PS::SPJMonopole
typedef struct {
   fdps_fSP mass;
   fdps_fSPvec pos;
} fdps_spj_monopole;

//**** PS::SPJQuadrupole
typedef struct {
   fdps_fSP mass;
   fdps_fSPvec pos;
   fdps_fSPmat quad;
} fdps_spj_quadrupole;

//**** PS::SPJMonopoleGeometricCenter
typedef struct {
   fdps_sSP n_ptcl;
   fdps_fSP charge;
   fdps_fSPvec pos;
} fdps_spj_monopole_geomcen;

//**** PS::SPJDipoleGeometricCenter
typedef struct {
   fdps_sSP n_ptcl;
   fdps_fSP charge;
   fdps_fSPvec pos;
   fdps_fSPvec dipole;
} fdps_spj_dipole_geomcen;

//**** PS::SPJQuadrupoleGeometricCenter
typedef struct {
   fdps_sSP n_ptcl;
   fdps_fSP charge;
   fdps_fSPvec pos;
   fdps_fSPvec dipole;
   fdps_fSPmat quadrupole;
} fdps_spj_quadrupole_geomcen;

//**** PS::SPJMonopoleScatter
typedef struct {
   fdps_fSP mass;
   fdps_fSPvec pos;
} fdps_spj_monopole_scatter;

//**** PS::SPJQuadrupoleScatter
typedef struct {
   fdps_fSP mass;
   fdps_fSPvec pos;
   fdps_fSPmat quad;
} fdps_spj_quadrupole_scatter;

//**** PS::SPJMonopoleCutoff
typedef struct {
   fdps_fSP mass;
   fdps_fSPvec pos;
} fdps_spj_monopole_cutoff;
