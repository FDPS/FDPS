!==================================
!   MODULE: FDPS super particle
!==================================
module fdps_super_particle
   use, intrinsic :: iso_c_binding
   use fdps_vector
   use fdps_matrix
   implicit none

#ifdef PARTICLE_SIMULATOR_SPMOM_F32
   !**** PS::SPJMonopole
   type, public, bind(c) :: fdps_spj_monopole
      real(kind=c_float) :: mass
      type(fdps_f32vec) :: pos
   end type fdps_spj_monopole

   !**** PS::SPJQuadrupole
   type, public, bind(c) :: fdps_spj_quadrupole
      real(kind=c_float) :: mass
      type(fdps_f32vec)  :: pos
      type(fdps_f32mat)  :: quad
   end type fdps_spj_quadrupole

   !**** PS::SPJMonopoleGeometricCenter
   type, public, bind(c) :: fdps_spj_monopole_geomcen
      integer(kind=c_int) :: n_ptcl
      real(kind=c_float) :: charge
      type(fdps_f32vec) :: pos
   end type fdps_spj_monopole_geomcen

   !**** PS::SPJDipoleGeometricCenter
   type, public, bind(c) :: fdps_spj_dipole_geomcen
      integer(kind=c_int) :: n_ptcl
      real(kind=c_float) :: charge
      type(fdps_f32vec) :: pos
      type(fdps_f32vec) :: dipole
   end type fdps_spj_dipole_geomcen

   !**** PS::SPJQuadrupoleGeometricCenter
   type, public, bind(c) :: fdps_spj_quadrupole_geomcen
      integer(kind=c_int) :: n_ptcl
      real(kind=c_float) :: charge
      type(fdps_f32vec) :: pos
      type(fdps_f32vec) :: dipole
      type(fdps_f32mat) :: quadrupole
   end type fdps_spj_quadrupole_geomcen

   !**** PS::SPJMonopoleScatter
   type, public, bind(c) :: fdps_spj_monopole_scatter
      real(kind=c_float) :: mass
      type(fdps_f32vec) :: pos
   end type fdps_spj_monopole_scatter

   !**** PS::SPJQuadrupoleScatter
   type, public, bind(c) :: fdps_spj_quadrupole_scatter
      real(kind=c_float) :: mass
      type(fdps_f32vec)  :: pos
      type(fdps_f32mat)  :: quad
   end type fdps_spj_quadrupole_scatter

   !**** PS::SPJMonopoleSymmetry
   type, public, bind(c) :: fdps_spj_monopole_symmetry
      real(kind=c_float) :: mass
      type(fdps_f32vec) :: pos
   end type fdps_spj_monopole_symmetry

   !**** PS::SPJQuadrupoleSymmetry
   type, public, bind(c) :: fdps_spj_quadrupole_symmetry
      real(kind=c_float) :: mass
      type(fdps_f32vec) :: pos
      type(fdps_f32mat) :: quad
   end type fdps_spj_quadrupole_symmetry

   !**** PS::SPJMonopoleCutoff
   type, public, bind(c) :: fdps_spj_monopole_cutoff
      real(kind=c_float) :: mass
      type(fdps_f32vec) :: pos
   end type fdps_spj_monopole_cutoff
#else
   !**** PS::SPJMonopole
   type, public, bind(c) :: fdps_spj_monopole
      real(kind=c_double) :: mass
      type(fdps_f64vec) :: pos
   end type fdps_spj_monopole

   !**** PS::SPJQuadrupole
   type, public, bind(c) :: fdps_spj_quadrupole
      real(kind=c_double) :: mass
      type(fdps_f64vec)  :: pos
      type(fdps_f64mat)  :: quad
   end type fdps_spj_quadrupole

   !**** PS::SPJMonopoleGeometricCenter
   type, public, bind(c) :: fdps_spj_monopole_geomcen
      integer(kind=c_long_long) :: n_ptcl
      real(kind=c_double) :: charge
      type(fdps_f64vec) :: pos
   end type fdps_spj_monopole_geomcen

   !**** PS::SPJDipoleGeometricCenter
   type, public, bind(c) :: fdps_spj_dipole_geomcen
      integer(kind=c_long_long) :: n_ptcl
      real(kind=c_double) :: charge
      type(fdps_f64vec) :: pos
      type(fdps_f64vec) :: dipole
   end type fdps_spj_dipole_geomcen

   !**** PS::SPJQuadrupoleGeometricCenter
   type, public, bind(c) :: fdps_spj_quadrupole_geomcen
      integer(kind=c_long_long) :: n_ptcl
      real(kind=c_double) :: charge
      type(fdps_f64vec) :: pos
      type(fdps_f64vec) :: dipole
      type(fdps_f64mat) :: quadrupole
   end type fdps_spj_quadrupole_geomcen

   !**** PS::SPJMonopoleScatter
   type, public, bind(c) :: fdps_spj_monopole_scatter
      real(kind=c_double) :: mass
      type(fdps_f64vec) :: pos
   end type fdps_spj_monopole_scatter

   !**** PS::SPJQuadrupoleScatter
   type, public, bind(c) :: fdps_spj_quadrupole_scatter
      real(kind=c_double) :: mass
      type(fdps_f64vec)  :: pos
      type(fdps_f64mat)  :: quad
   end type fdps_spj_quadrupole_scatter

   !**** PS::SPJMonopoleSymmetry
   type, public, bind(c) :: fdps_spj_monopole_symmetry
      real(kind=c_double) :: mass
      type(fdps_f64vec) :: pos
   end type fdps_spj_monopole_symmetry

   !**** PS::SPJQuadrupoleSymmetry
   type, public, bind(c) :: fdps_spj_quadrupole_symmetry
      real(kind=c_double) :: mass
      type(fdps_f64vec) :: pos
      type(fdps_f64mat) :: quad
   end type fdps_spj_quadrupole_symmetry

   !**** PS::SPJMonopoleCutoff
   type, public, bind(c) :: fdps_spj_monopole_cutoff
      real(kind=c_double) :: mass
      type(fdps_f64vec) :: pos
   end type fdps_spj_monopole_cutoff
#endif

   ! [TODO]
   !    PS::SPJMonopolePeriodic
   !    PS::SPJMonopoleCutoffScatter

end module fdps_super_particle
