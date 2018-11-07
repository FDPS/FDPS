#include "macro_defs.h"
!================================
!   MODULE: User defined types
!================================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   use fdps_super_particle
   use mathematical_constants
   implicit none

   !* Public parameters & variables
   real(kind=c_double), parameter, public :: specific_heat_ratio = 5.0d0/3.0d0
   real(kind=c_double), parameter, public :: CFL_dyn = 0.3d0
   real(kind=c_double), parameter, public :: CFL_hydro = 0.3d0
   real(kind=c_double), parameter, public :: scf_smth = 1.25d0
   integer(kind=c_int), parameter, public :: N_neighbor = 50
   real(kind=c_double), public :: eps_grav
   real(kind=c_double), public :: mass_avg
   real(kind=c_double), public :: dt_max

   !**** Force types
   type, public, bind(c) :: force_grav !$fdps Force
      !$fdps clear
      type(fdps_f64vec) :: acc
      real(kind=c_double) :: pot
   end type force_grav

   type, public, bind(c) :: force_dens !$fdps Force
      !$fdps clear smth=keep
      integer(kind=c_int) :: flag
      real(kind=c_double) :: dens
      real(kind=c_double) :: smth
      real(kind=c_double) :: gradh
      real(kind=c_double) :: divv
      type(fdps_f64vec) :: rotv
   end type force_dens

   type, public, bind(c) :: force_hydro !$fdps Force
      !$fdps clear 
      type(fdps_f64vec) :: acc
      real(kind=c_double) :: eng_dot
      real(kind=c_double) :: ent_dot
      real(kind=c_double) :: dt
   end type force_hydro

   !**** Full particle type
   type, public, bind(c) :: fp_nbody !$fdps FP
      !$fdps copyFromForce force_grav (acc,acc) (pot,pot)
      integer(kind=c_long_long) :: id !$fdps id
      real(kind=c_double) :: mass !$fdps charge
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel
      type(fdps_f64vec) :: acc
      real(kind=c_double) :: pot
   end type fp_nbody

   type, public, bind(c) :: fp_sph !$fdps FP
      !$fdps copyFromForce force_grav (acc,acc_grav) (pot,pot_grav)
      !$fdps copyFromForce force_dens (flag,flag) (dens,dens) (smth,smth) (gradh,gradh) (divv,divv) (rotv,rotv)
      !$fdps copyFromForce force_hydro (acc,acc_hydro) (eng_dot,eng_dot) (ent_dot,ent_dot) (dt,dt)
      integer(kind=c_long_long) :: id !$fdps id
      real(kind=c_double) :: mass !$fdps charge
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel 
      type(fdps_f64vec) :: acc_grav
      real(kind=c_double) :: pot_grav
      type(fdps_f64vec)   :: acc_hydro
      integer(kind=c_int) :: flag
      real(kind=c_double) :: dens
      real(kind=c_double) :: eng
      real(kind=c_double) :: ent 
      real(kind=c_double) :: pres
      real(kind=c_double) :: smth
      real(kind=c_double) :: gradh
      real(kind=c_double) :: divv
      type(fdps_f64vec)   :: rotv
      real(kind=c_double) :: balsw
      real(kind=c_double) :: snds
      real(kind=c_double) :: eng_dot
      real(kind=c_double) :: ent_dot
      real(kind=c_double) :: dt
      type(fdps_f64vec)   :: vel_half
      real(kind=c_double) :: eng_half
      real(kind=c_double) :: ent_half
   end type fp_sph

   !**** Essential particle types
   type, public, bind(c) :: ep_grav !$fdps EPI,EPJ
      !$fdps copyFromFP fp_nbody (id,id) (mass,mass) (pos,pos)
      !$fdps copyFromFP fp_sph (id,id) (mass,mass) (pos,pos)
      integer(kind=c_long_long) :: id !$fdps id
      real(kind=c_double) :: mass !$fdps charge
      type(fdps_f64vec) :: pos !$fdps position
   end type ep_grav

   type, public, bind(c) :: ep_hydro !$fdps EPI,EPJ
      !$fdps copyFromFP fp_sph (id,id) (pos,pos) (vel,vel) (mass,mass) (smth,smth) (dens,dens) (pres,pres) (gradh,gradh) (snds,snds) (balsw,balsw)
      integer(kind=c_long_long) :: id !$fdps id
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel
      real(kind=c_double) :: mass !$fdps charge
      real(kind=c_double) :: smth !$fdps rsearch
      real(kind=c_double) :: dens
      real(kind=c_double) :: pres
      real(kind=c_double) :: gradh
      real(kind=c_double) :: snds
      real(kind=c_double) :: balsw
   end type ep_hydro

   !* Public routines
   public :: W
   public :: gradW
   public :: dWdh
   public :: calc_gravity_ep_ep
   public :: calc_gravity_ep_sp
   public :: calc_density
   public :: calc_hydro_force

   contains

   !-------------------------------------------------------------------
   pure function W(r,h)
      ! M4 Cubic spline kernel
      ! (see Eq. (4) in Springel (2005)[MNRAS,364,1105])
      implicit none
      real(kind=c_double) :: W
      real(kind=c_double), intent(in) :: r,h
      !* Local variables
      real(kind=c_double) :: u,cc,u2,s
     
      u = r/h
      cc=8.0d0/(pi*h*h*h)
      if (u <= 0.5d0) then
          u2 = u*u
          W = cc*(1.0d0+u2*6.0d0*(-1.0d0+u))
      else if ((0.5d0 < u) .and. (u <= 1.0d0)) then
          s = 1.0d0-u
          W = cc*2.0d0*s*s*s
      else
          W = 0.0d0
      end if

   end function W

   !-------------------------------------------------------------------
   pure function gradW(dr,h)
      ! This subroutine gives \nabla W(r,h), i.e.,
      ! \dfrac{\partial W(r,h)}{\partial r}\dfrac{dr}{r}.
      implicit none
      type(fdps_f64vec) :: gradW
      type(fdps_f64vec), intent(in) :: dr
      real(kind=c_double), intent(in) :: h
      !* Local variables
      real(kind=c_double) :: r,u,cc

      r = dsqrt(dr%x * dr%x &
               +dr%y * dr%y &
               +dr%z * dr%z)
      u=r/h;
      cc = -48.0d0/(pi*h*h*h*h);
#if defined(USE_PRESCR_OF_THOMAS_COUCHMAN_1992)
      if (u <= 1.0d0/3.0d0) then
         cc = cc*(1.0d0/3.0d0)/(r)
      else if ((1.0d0/3.0d0 < u) .and. (u <= 0.5d0)) then
         cc = cc*u*(2.0d0-3.0d0*u)/(r)
      else if ((0.5d0 < u) .and. (u < 1.0d0)) then
         cc = cc*(1.0d0-u)*(1.0d0-u)/(r)
      else 
         cc = 0.0d0
      end if
#else
      if ((0.0d0 < u) .and. (u <= 0.5d0)) then
         cc = cc*u*(2.0d0-3.0d0*u)/(r)
      else if ((0.5d0 < u) .and. (u < 1.0d0)) then
         cc = cc*(1.0d0-u)*(1.0d0-u)/(r)
      else 
         ! r=0 case is included in this branch
         cc = 0.0d0
      end if
#endif
      gradW%x = dr%x * cc
      gradW%y = dr%y * cc
      gradW%z = dr%z * cc

   end function gradW

   !-------------------------------------------------------------------
   pure function dWdh(r, h)
      ! This subroutine gives dW(r,h)/dh, i.e.,
      ! \dfrac{\partial W(r,h)}{\partial h}.
      implicit none
      real(kind=c_double) :: dWdh
      real(kind=c_double), intent(in) :: r,h
      !* Local variables
      real(kind=c_double) :: u,u2,s,cc

      u=r/h
      cc=-24.0d0/(pi*h*h*h*h)
      if (u <= 0.5d0) then
         u2 = u*u
         dWdh = cc*(1.0d0            &
                   +u2*(-10.0d0      &
                        +12.0d0*u))
      else if ((0.5d0 < u) .and. (u < 1.0d0)) then
         s = 1.0d0-u
         dWdh = cc*2.0d0*s*s*(1.0d0-2.0d0*u)
      else
         dWdh = 0.0d0
      end if

   end function dWdh

   !**** Interaction functions
#if defined(ENABLE_PHANTOM_GRAPE_X86)
   subroutine calc_gravity_ep_ep(ep_i,n_ip,ep_j,n_jp,f) bind(c)
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
      use omp_lib
#endif
      use phantom_grape_g5_x86
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(ep_grav), dimension(n_ip), intent(in) :: ep_i
      type(ep_grav), dimension(n_jp), intent(in) :: ep_j
      type(force_grav), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      integer(c_int) :: nipipe,njpipe,devid
      real(c_double), dimension(3,n_ip) :: xi,ai
      real(c_double), dimension(n_ip) :: pi
      real(c_double), dimension(3,n_jp) :: xj
      real(c_double), dimension(n_jp) :: mj

      nipipe = n_ip
      njpipe = n_jp
      do i=1,n_ip
         xi(1,i) = ep_i(i)%pos%x
         xi(2,i) = ep_i(i)%pos%y
         xi(3,i) = ep_i(i)%pos%z
         ai(1,i) = 0.0d0
         ai(2,i) = 0.0d0
         ai(3,i) = 0.0d0
         pi(i)   = 0.0d0
      end do
      do j=1,n_jp
         xj(1,j) = ep_j(j)%pos%x
         xj(2,j) = ep_j(j)%pos%y
         xj(3,j) = ep_j(j)%pos%z
         mj(j)   = ep_j(j)%mass
      end do
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
      devid = omp_get_thread_num()
      ! [IMPORTANT NOTE]
      !   The subroutine calc_gravity_ep_ep is called by a OpenMP thread
      !   in the FDPS. This means that here is already in the parallel region.
      !   So, you can use omp_get_thread_num() without !$OMP parallel directives.
      !   If you use them, a nested parallel resions is made and the gravity
      !   calculation will not be performed correctly.
#else
      devid = 0
#endif
      call g5_set_xmjMC(devid, 0, n_jp, xj, mj)
      call g5_set_nMC(devid, n_jp)
      call g5_calculate_force_on_xMC(devid, xi, ai, pi, n_ip)
      do i=1,n_ip
         f(i)%acc%x = f(i)%acc%x + ai(1,i)
         f(i)%acc%y = f(i)%acc%y + ai(2,i)
         f(i)%acc%z = f(i)%acc%z + ai(3,i)
         f(i)%pot   = f(i)%pot   - pi(i)
      end do
   end subroutine calc_gravity_ep_ep

   subroutine calc_gravity_ep_sp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
      use omp_lib
#endif
      use phantom_grape_g5_x86
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(ep_grav), dimension(n_ip), intent(in) :: ep_i
      type(fdps_spj_monopole), dimension(n_jp), intent(in) :: ep_j
      type(force_grav), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      integer(c_int) :: nipipe,njpipe,devid
      real(c_double), dimension(3,n_ip) :: xi,ai
      real(c_double), dimension(n_ip) :: pi
      real(c_double), dimension(3,n_jp) :: xj
      real(c_double), dimension(n_jp) :: mj

      nipipe = n_ip
      njpipe = n_jp
      do i=1,n_ip
         xi(1,i) = ep_i(i)%pos%x
         xi(2,i) = ep_i(i)%pos%y
         xi(3,i) = ep_i(i)%pos%z
         ai(1,i) = 0.0d0
         ai(2,i) = 0.0d0
         ai(3,i) = 0.0d0
         pi(i)   = 0.0d0
      end do
      do j=1,n_jp
         xj(1,j) = ep_j(j)%pos%x
         xj(2,j) = ep_j(j)%pos%y
         xj(3,j) = ep_j(j)%pos%z
         mj(j)   = ep_j(j)%mass
      end do
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
      devid = omp_get_thread_num()
      ! [IMPORTANT NOTE]
      !   The subroutine calc_gravity_ep_sp is called by a OpenMP thread
      !   in the FDPS. This means that here is already in the parallel region.
      !   So, you can use omp_get_thread_num() without !$OMP parallel directives.
      !   If you use them, a nested parallel resions is made and the gravity
      !   calculation will not be performed correctly.
#else
      devid = 0
#endif
      call g5_set_xmjMC(devid, 0, n_jp, xj, mj)
      call g5_set_nMC(devid, n_jp)
      call g5_calculate_force_on_xMC(devid, xi, ai, pi, n_ip)
      do i=1,n_ip
         f(i)%acc%x = f(i)%acc%x + ai(1,i)
         f(i)%acc%y = f(i)%acc%y + ai(2,i)
         f(i)%acc%z = f(i)%acc%z + ai(3,i)
         f(i)%pot   = f(i)%pot   - pi(i)
      end do
   end subroutine calc_gravity_ep_sp
#else
   subroutine calc_gravity_ep_ep(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(kind=c_int), intent(in), value :: n_ip,n_jp
      type(ep_grav), dimension(n_ip), intent(in) :: ep_i
      type(ep_grav), dimension(n_jp), intent(in) :: ep_j
      type(force_grav), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(kind=c_int) :: i,j
      real(kind=c_double) :: eps2,poti,r3_inv,r_inv
      type(fdps_f64vec) :: xi,ai,rij
      !* Compute force
      eps2 = eps_grav * eps_grav
      do i=1,n_ip
         xi%x = ep_i(i)%pos%x
         xi%y = ep_i(i)%pos%y
         xi%z = ep_i(i)%pos%z
         ai%x = 0.0d0
         ai%y = 0.0d0
         ai%z = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            rij%x  = xi%x - ep_j(j)%pos%x
            rij%y  = xi%y - ep_j(j)%pos%y
            rij%z  = xi%z - ep_j(j)%pos%z
            r3_inv = rij%x*rij%x &
                   + rij%y*rij%y &
                   + rij%z*rij%z &
                   + eps2
            r_inv  = 1.0d0/dsqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * ep_j(j)%mass
            r3_inv = r3_inv * r_inv
            ai%x   = ai%x - r3_inv * rij%x
            ai%y   = ai%y - r3_inv * rij%y
            ai%z   = ai%z - r3_inv * rij%z
            poti   = poti - r_inv
         end do
         f(i)%acc%x = f(i)%acc%x + ai%x
         f(i)%acc%y = f(i)%acc%y + ai%y
         f(i)%acc%z = f(i)%acc%z + ai%z
         f(i)%pot   = f(i)%pot   + poti
      end do
   end subroutine calc_gravity_ep_ep

   subroutine calc_gravity_ep_sp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(kind=c_int), intent(in), value :: n_ip,n_jp
      type(ep_grav), dimension(n_ip), intent(in) :: ep_i
      type(fdps_spj_monopole), dimension(n_jp), intent(in) :: ep_j
      type(force_grav), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(kind=c_int) :: i,j
      real(kind=c_double) :: eps2,poti,r3_inv,r_inv
      type(fdps_f64vec) :: xi,ai,rij
      !* Compute force
      eps2 = eps_grav * eps_grav
      do i=1,n_ip
         xi%x = ep_i(i)%pos%x
         xi%y = ep_i(i)%pos%y
         xi%z = ep_i(i)%pos%z
         ai%x = 0.0d0
         ai%y = 0.0d0
         ai%z = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            rij%x  = xi%x - ep_j(j)%pos%x
            rij%y  = xi%y - ep_j(j)%pos%y
            rij%z  = xi%z - ep_j(j)%pos%z
            r3_inv = rij%x*rij%x &
                   + rij%y*rij%y &
                   + rij%z*rij%z &
                   + eps2
            r_inv  = 1.0d0/dsqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * ep_j(j)%mass
            r3_inv = r3_inv * r_inv
            ai%x   = ai%x - r3_inv * rij%x
            ai%y   = ai%y - r3_inv * rij%y
            ai%z   = ai%z - r3_inv * rij%z
            poti   = poti - r_inv
         end do
         f(i)%acc%x = f(i)%acc%x + ai%x
         f(i)%acc%y = f(i)%acc%y + ai%y
         f(i)%acc%z = f(i)%acc%z + ai%z
         f(i)%pot   = f(i)%pot   + poti
      end do
   end subroutine calc_gravity_ep_sp
#endif

   subroutine calc_density(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(kind=c_int), intent(in), value :: n_ip,n_jp
      type(ep_hydro), dimension(n_ip), intent(in) :: ep_i
      type(ep_hydro), dimension(n_jp), intent(in) :: ep_j
      type(force_dens), dimension(n_ip), intent(inout) :: f
      !* Local parameters
      real(kind=c_double), parameter :: eps=1.0d-6
      !* Local variables
      integer(kind=c_int) :: i,j
      integer(kind=c_int) :: n_unchanged
      real(kind=c_double) :: M,M_trgt
      real(kind=c_double) :: dens,drho_dh
      real(kind=c_double) :: h,h_max_alw,h_L,h_U,dh,dh_prev
      type(fdps_f64vec) :: dr,dv,gradW_i

#if defined(ENABLE_VARIABLE_SMOOTHING_LENGTH)
      real(kind=c_double), dimension(n_jp) :: mj,rij
      M_trgt = mass_avg * N_neighbor
      do i=1,n_ip
          dens = 0.0d0
          h_max_alw = ep_i(i)%smth ! maximum allowance
          h = h_max_alw / SCF_smth
          ! Note that we increase smth by a factor of scf_smth
          ! before calling calc_density().
          h_L = 0.0d0
          h_U = h_max_alw
          dh_prev = 0.0d0
          n_unchanged = 0
          ! Software cache
          do j=1,n_jp
             mj(j) = ep_j(j)%mass
             dr%x = ep_i(i)%pos%x - ep_j(j)%pos%x
             dr%y = ep_i(i)%pos%y - ep_j(j)%pos%y
             dr%z = ep_i(i)%pos%z - ep_j(j)%pos%z
             rij(j) = dsqrt(dr%x * dr%x &
                           +dr%y * dr%y &
                           +dr%z * dr%z)
          end do
          iteration_loop: do
              ! Calculate density
              dens = 0.0d0
              do j=1,n_jp
                 dens = dens + mj(j) * W(rij(j), h)
              end do
              ! Check if the current value of the smoohting length satisfies 
              ! Eq.(5) in Springel (2005).
              M = 4.0d0 * pi * h * h * h * dens / 3.0d0
              if ((h < h_max_alw) .and. (dabs(M/M_trgt - 1.0d0) < eps)) then
                  ! In this case, Eq.(5) holds within a specified accuracy.
                  f(i)%flag = 1
                  f(i)%dens = dens
                  f(i)%smth = h
                  exit iteration_loop
              end if
              if (((h == h_max_alw) .and. (M < M_trgt)) .or. (n_unchanged == 4)) then
                  ! In this case, we skip this particle forcibly.
                  ! In order to determine consistently the density
                  ! and the smoohting length for this particle,
                  ! we must re-perform calcForceAllAndWriteBack().
                  f(i)%flag = 0
                  f(i)%dens = dens
                  f(i)%smth = h_max_alw
                  exit iteration_loop
              end if
              ! Update h_L & h_U
              if (M < M_trgt) then
                 if (h_L < h) h_L = h
              else if (M_trgt < M) then
                 if (h < h_U) h_U = h
              end if
              dh = h_U - h_L
              if (dh == dh_prev) then
                 n_unchanged = n_unchanged + 1
              else 
                 dh_prev = dh
                 n_unchanged = 0
              end if
              ! Update smoothing length
              h = ((3.0d0 * M_trgt)/(4.0d0 * pi * dens))**(1.0d0/3.0d0)
              if ((h <= h_L) .or. (h == h_U)) then
                 ! In this case, we switch to the bisection search.
                 ! The inclusion of '=' in the if statement is very
                 ! important to escape a limit cycle.
                 h = 0.5d0 * (h_L + h_U)
              else if (h_U < h) then
                 h = h_U
              end if
          end do iteration_loop
          ! Calculate grad-h term
          if (f(i)%flag == 1) then
              drho_dh = 0.0d0
              do j=1,n_jp
                 drho_dh = drho_dh + mj(j) * dWdh(rij(j), h)
              end do
              f(i)%gradh = 1.0d0 / (1.0d0 + (h * drho_dh) / (3.0d0 * dens))
          else 
              f(i)%gradh = 1.0d0 ! dummy value
          end if
          ! Compute \div v & \rot v for Balsara switch
#if defined(USE_BALSARA_SWITCH)
          do j=1,n_jp
             dr%x = ep_i(i)%pos%x - ep_j(j)%pos%x
             dr%y = ep_i(i)%pos%y - ep_j(j)%pos%y
             dr%z = ep_i(i)%pos%z - ep_j(j)%pos%z
             dv%x = ep_i(i)%vel%x - ep_j(j)%vel%x
             dv%y = ep_i(i)%vel%y - ep_j(j)%vel%y
             dv%z = ep_i(i)%vel%z - ep_j(j)%vel%z
             gradW_i = gradW(dr, f(i)%smth)
             f(i)%divv = f(i)%divv - mj(j) * (dv%x * gradW_i%x &
                                             +dv%y * gradW_i%y &
                                             +dv%z * gradW_i%z)
             f(i)%rotv%x = f(i)%rotv%x - mj(j) * (dv%y * gradW_i%z - dv%z * gradW_i%y)
             f(i)%rotv%y = f(i)%rotv%y - mj(j) * (dv%z * gradW_i%x - dv%x * gradW_i%z)
             f(i)%rotv%z = f(i)%rotv%z - mj(j) * (dv%x * gradW_i%y - dv%y * gradW_i%x)
          end do
          f(i)%divv   = f(i)%divv   / f(i)%dens
          f(i)%rotv%x = f(i)%rotv%x / f(i)%dens
          f(i)%rotv%y = f(i)%rotv%y / f(i)%dens
          f(i)%rotv%z = f(i)%rotv%z / f(i)%dens
#endif
      end do
#else
      double precision :: mj,rij
      do i=1,n_ip
         f(i)%dens = 0.0d0
         do j=1,n_jp
            dr%x = ep_j(j)%pos%x - ep_i(i)%pos%x
            dr%y = ep_j(j)%pos%y - ep_i(i)%pos%y
            dr%z = ep_j(j)%pos%z - ep_i(i)%pos%z
            rij = dsqrt(dr%x * dr%x &
                       +dr%y * dr%y &
                       +dr%z * dr%z)
            f(i)%dens = f(i)%dens &
                      + ep_j(j)%mass * W(rij,ep_i(i)%smth)
         end do
         f(i)%smth = ep_i(i)%smth
         f(i)%gradh = 1.0d0
         ! Compute \div v & \rot v for Balsara switch
#if defined(USE_BALSARA_SWITCH)
         do j=1,n_jp
            mj = ep_j(j)%mass
            dr%x = ep_i(i)%pos%x - ep_j(j)%pos%x
            dr%y = ep_i(i)%pos%y - ep_j(j)%pos%y
            dr%z = ep_i(i)%pos%z - ep_j(j)%pos%z
            dv%x = ep_i(i)%vel%x - ep_j(j)%vel%x
            dv%y = ep_i(i)%vel%y - ep_j(j)%vel%y
            dv%z = ep_i(i)%vel%z - ep_j(j)%vel%z
            gradW_i = gradW(dr, f(i)%smth)
            f(i)%divv = f(i)%divv - mj * (dv%x * gradW_i%x &
                                         +dv%y * gradW_i%y &
                                         +dv%z * gradW_i%z)
            f(i)%rotv%x = f(i)%rotv%x - mj * (dv%y * gradW_i%z - dv%z * gradW_i%y)
            f(i)%rotv%y = f(i)%rotv%y - mj * (dv%z * gradW_i%x - dv%x * gradW_i%z)
            f(i)%rotv%z = f(i)%rotv%z - mj * (dv%x * gradW_i%y - dv%y * gradW_i%x)
         end do
         f(i)%divv   = f(i)%divv   / f(i)%dens
         f(i)%rotv%x = f(i)%rotv%x / f(i)%dens
         f(i)%rotv%y = f(i)%rotv%y / f(i)%dens
         f(i)%rotv%z = f(i)%rotv%z / f(i)%dens
#endif
      end do
#endif

   end subroutine calc_density

   !**** Interaction function
   subroutine calc_hydro_force(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(kind=c_int), intent(in), value :: n_ip,n_jp
      type(ep_hydro), dimension(n_ip), intent(in) :: ep_i
      type(ep_hydro), dimension(n_jp), intent(in) :: ep_j
      type(force_hydro), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(kind=c_int) :: i,j
      real(kind=c_double) :: mass_i,mass_j,smth_i,smth_j, &
                             dens_i,dens_j,pres_i,pres_j, &
                             gradh_i,gradh_j,balsw_i,balsw_j, &
                             snds_i,snds_j
      real(kind=c_double) :: povrho2_i,povrho2_j, &
                             v_sig_max,dr_dv,w_ij,v_sig,AV
      type(fdps_f64vec) :: pos_i,pos_j,vel_i,vel_j, &
                           dr,dv,gradW_i,gradW_j,gradW_ij
      do i=1,n_ip
         !* Zero-clear
         v_sig_max = 0.0d0
         !* Extract i-particle info.
         pos_i = ep_i(i)%pos
         vel_i = ep_i(i)%vel
         mass_i  = ep_i(i)%mass
         smth_i  = ep_i(i)%smth
         dens_i  = ep_i(i)%dens
         pres_i  = ep_i(i)%pres
         gradh_i = ep_i(i)%gradh
         balsw_i = ep_i(i)%balsw
         snds_i  = ep_i(i)%snds
         povrho2_i = pres_i/(dens_i*dens_i)
         do j=1,n_jp
            !* Extract j-particle info.
            pos_j%x = ep_j(j)%pos%x
            pos_j%y = ep_j(j)%pos%y
            pos_j%z = ep_j(j)%pos%z
            vel_j%x = ep_j(j)%vel%x
            vel_j%y = ep_j(j)%vel%y
            vel_j%z = ep_j(j)%vel%z
            mass_j  = ep_j(j)%mass
            smth_j  = ep_j(j)%smth
            dens_j  = ep_j(j)%dens
            pres_j  = ep_j(j)%pres
            gradh_j = ep_j(j)%gradh
            balsw_j = ep_j(j)%balsw
            snds_j  = ep_j(j)%snds
            povrho2_j = pres_j/(dens_j*dens_j)
            !* Compute dr & dv
            dr%x = pos_i%x - pos_j%x
            dr%y = pos_i%y - pos_j%y
            dr%z = pos_i%z - pos_j%z
            dv%x = vel_i%x - vel_j%x
            dv%y = vel_i%y - vel_j%y
            dv%z = vel_i%z - vel_j%z
            !* Compute the signal velocity
            dr_dv = dr%x * dv%x + dr%y * dv%y + dr%z * dv%z
            if (dr_dv < 0.0d0) then
               w_ij = dr_dv / sqrt(dr%x * dr%x + dr%y * dr%y + dr%z * dr%z)
            else
               w_ij = 0.0d0
            end if
            v_sig = snds_i + snds_j - 3.0d0 * w_ij
            v_sig_max = max(v_sig_max, v_sig)
            !* Compute the artificial viscosity
            AV = - 0.5d0*v_sig*w_ij / (0.5d0*(dens_i+dens_j)) * 0.5d0*(balsw_i+balsw_j)
            !* Compute the average of the gradients of kernel
            gradW_i  = gradW(dr,smth_i)
            gradW_j  = gradW(dr,smth_j)
            gradW_ij%x = 0.5d0 * (gradW_i%x + gradW_j%x)
            gradW_ij%y = 0.5d0 * (gradW_i%y + gradW_j%y)
            gradW_ij%z = 0.5d0 * (gradW_i%z + gradW_j%z)
            !* Compute the acceleration and the heating rate
            f(i)%acc%x = f(i)%acc%x - mass_j*(gradh_i * povrho2_i * gradW_i%x &
                                             +gradh_j * povrho2_j * gradW_j%x &
                                             +AV * gradW_ij%x)
            f(i)%acc%y = f(i)%acc%y - mass_j*(gradh_i * povrho2_i * gradW_i%y &
                                             +gradh_j * povrho2_j * gradW_j%y &
                                             +AV * gradW_ij%y)
            f(i)%acc%z = f(i)%acc%z - mass_j*(gradh_i * povrho2_i * gradW_i%z &
                                             +gradh_j * povrho2_j * gradW_j%z &
                                             +AV * gradW_ij%z)
            f(i)%eng_dot = f(i)%eng_dot                                      &
                         + mass_j * gradh_i * povrho2_i * (dv%x * gradW_i%x  &
                                                          +dv%y * gradW_i%y  &
                                                          +dv%z * gradW_i%z) &
                         + mass_j * 0.5d0 * AV * (dv%x * gradW_ij%x          &
                                                 +dv%y * gradW_ij%y          &
                                                 +dv%z * gradW_ij%z)
            f(i)%ent_dot = f(i)%ent_dot                           &
                         + 0.5 * mass_j * AV * (dv%x * gradW_ij%x &
                                               +dv%y * gradW_ij%y &
                                               +dv%z * gradW_ij%z)
         end do
         f(i)%ent_dot = f(i)%ent_dot                  &
                      * (specific_heat_ratio - 1.0d0) &
                      / dens_i**(specific_heat_ratio - 1.0d0)
         f(i)%dt = CFL_hydro*2.0d0*smth_i/v_sig_max
      end do
   end subroutine calc_hydro_force

end module user_defined_types

