!===============================
!   MODULE: User defined types
!===============================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   use fdps_super_particle
   implicit none

   !* Public variables
   real(kind=c_double), public :: eps_grav ! gravitational softening

   !**** Full particle type
   type, public, bind(c) :: full_particle !$fdps FP,EPI,EPJ,Force
      !$fdps copyFromForce full_particle (pot,pot) (acc,acc)
      !$fdps copyFromFP full_particle (id,id) (mass,mass) (pos,pos) 
      !$fdps clear id=keep, mass=keep, pos=keep, vel=keep
      integer(kind=c_long_long) :: id
      real(kind=c_double)  mass !$fdps charge
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel !$fdps velocity
      real(kind=c_double) :: pot
      type(fdps_f64vec) :: acc
   end type full_particle

   !* The following types are used in PIKG-generated kenrels
   type, public, bind(c) :: epi_grav
      type(fdps_f32vec) :: pos
   end type epi_grav

   type, public, bind(c) :: epj_grav
      type(fdps_f32vec) :: pos
      real(kind=c_float) :: mass
   end type epj_grav

   type, public, bind(c) :: force_grav
      type(fdps_f32vec) :: acc
      real(kind=c_float) :: pot
   end type force_grav
  
   contains

   !**** Interaction function (particle-particle)
#if defined(ENABLE_PHANTOM_GRAPE_X86)
   subroutine calc_gravity_ep_ep(ep_i,n_ip,ep_j,n_jp,f) bind(c)
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
      use omp_lib
#endif
      use phantom_grape_g5_x86
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(full_particle), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
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
      !   The subroutine calc_gravity_pp is called by a OpenMP thread
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
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(fdps_spj_monopole), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
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
      !   The subroutine calc_gravity_psp is called by a OpenMP thread
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
#elif defined(ENABLE_PIKG_KERNEL_X86)
   subroutine calc_gravity_ep_ep(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      use, intrinsic :: iso_c_binding
      use pikg_module_ep_ep
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(full_particle), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      type(epi_grav), dimension(n_ip), target :: ep_i_tmp
      type(epj_grav), dimension(n_jp), target :: ep_j_tmp
      type(force_grav), dimension(n_ip), target :: f_tmp

      if (n_ip > 0) then
         do i=1,n_ip
            ep_i_tmp(i)%pos%x = ep_i(i)%pos%x - ep_i(1)%pos%x
            ep_i_tmp(i)%pos%y = ep_i(i)%pos%y - ep_i(1)%pos%y
            ep_i_tmp(i)%pos%z = ep_i(i)%pos%z - ep_i(1)%pos%z
            f_tmp(i)%acc%x = 0.0
            f_tmp(i)%acc%y = 0.0
            f_tmp(i)%acc%z = 0.0
            f_tmp(i)%pot   = 0.0
         end do
         do j=1,n_jp
            ep_j_tmp(j)%pos%x = ep_j(j)%pos%x - ep_i(1)%pos%x
            ep_j_tmp(j)%pos%y = ep_j(j)%pos%y - ep_i(1)%pos%y
            ep_j_tmp(j)%pos%z = ep_j(j)%pos%z - ep_i(1)%pos%z
            ep_j_tmp(j)%mass  = ep_j(j)%mass
         end do
         call pikg_calc_grav_ep_ep(c_loc(ep_i_tmp), n_ip, &
                                   c_loc(ep_j_tmp), n_jp, &
                                   c_loc(f_tmp))
         do i=1,n_ip
            f(i)%acc%x = f(i)%acc%x + f_tmp(i)%acc%x
            f(i)%acc%y = f(i)%acc%y + f_tmp(i)%acc%y
            f(i)%acc%z = f(i)%acc%z + f_tmp(i)%acc%z
            f(i)%pot   = f(i)%pot   + f_tmp(i)%pot
         end do
      end if

   end subroutine calc_gravity_ep_ep

   subroutine calc_gravity_ep_sp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      use, intrinsic :: iso_c_binding
      use pikg_module_ep_ep
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(fdps_spj_monopole), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      type(epi_grav), dimension(n_ip), target :: ep_i_tmp
      type(epj_grav), dimension(n_jp), target :: ep_j_tmp
      type(force_grav), dimension(n_ip), target :: f_tmp

      if (n_ip > 0) then
         do i=1,n_ip
            ep_i_tmp(i)%pos%x = ep_i(i)%pos%x - ep_i(1)%pos%x
            ep_i_tmp(i)%pos%y = ep_i(i)%pos%y - ep_i(1)%pos%y
            ep_i_tmp(i)%pos%z = ep_i(i)%pos%z - ep_i(1)%pos%z
            f_tmp(i)%acc%x = 0.0
            f_tmp(i)%acc%y = 0.0
            f_tmp(i)%acc%z = 0.0
            f_tmp(i)%pot   = 0.0
         end do
         do j=1,n_jp
            ep_j_tmp(j)%pos%x = ep_j(j)%pos%x - ep_i(1)%pos%x
            ep_j_tmp(j)%pos%y = ep_j(j)%pos%y - ep_i(1)%pos%y
            ep_j_tmp(j)%pos%z = ep_j(j)%pos%z - ep_i(1)%pos%z
            ep_j_tmp(j)%mass  = ep_j(j)%mass
         end do
         call pikg_calc_grav_ep_ep(c_loc(ep_i_tmp), n_ip, &
                                   c_loc(ep_j_tmp), n_jp, &
                                   c_loc(f_tmp))
         do i=1,n_ip
            f(i)%acc%x = f(i)%acc%x + f_tmp(i)%acc%x
            f(i)%acc%y = f(i)%acc%y + f_tmp(i)%acc%y
            f(i)%acc%z = f(i)%acc%z + f_tmp(i)%acc%z
            f(i)%pot   = f(i)%pot   + f_tmp(i)%pot
         end do
      end if

   end subroutine calc_gravity_ep_sp
#else
   subroutine calc_gravity_ep_ep(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(full_particle), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      real(c_double) :: eps2,poti,r3_inv,r_inv
      type(fdps_f64vec) :: xi,ai,rij

      !* Compute force
      eps2 = eps_grav * eps_grav
      do i=1,n_ip
         xi = ep_i(i)%pos
         ai = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            rij%x  = xi%x - ep_j(j)%pos%x
            rij%y  = xi%y - ep_j(j)%pos%y
            rij%z  = xi%z - ep_j(j)%pos%z
            r3_inv = rij%x*rij%x &
                   + rij%y*rij%y &
                   + rij%z*rij%z &
                   + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * ep_j(j)%mass
            r3_inv = r3_inv * r_inv
            ai%x   = ai%x - r3_inv * rij%x
            ai%y   = ai%y - r3_inv * rij%y
            ai%z   = ai%z - r3_inv * rij%z
            poti   = poti - r_inv
            ! [IMPORTANT NOTE]
            !   In the innermost loop, we use the components of vectors
            !   directly for vector operations because of the following
            !   reasion. Except for intel compilers with `-ipo` option,
            !   most of Fortran compilers use function calls to perform
            !   vector operations like rij = x - ep_j(j)%pos.
            !   This significantly slow downs the speed of the code.
            !   By using the components of vector directly, we can avoid 
            !   these function calls.
         end do
         f(i)%pot = f(i)%pot + poti
         f(i)%acc = f(i)%acc + ai
      end do

   end subroutine calc_gravity_ep_ep

   !**** Interaction function (particle-super particle)
   subroutine calc_gravity_ep_sp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      implicit none
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(fdps_spj_monopole), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      real(c_double) :: eps2,poti,r3_inv,r_inv
      type(fdps_f64vec) :: xi,ai,rij

      eps2 = eps_grav * eps_grav
      do i=1,n_ip
         xi = ep_i(i)%pos
         ai = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            rij%x  = xi%x - ep_j(j)%pos%x
            rij%y  = xi%y - ep_j(j)%pos%y
            rij%z  = xi%z - ep_j(j)%pos%z
            r3_inv = rij%x*rij%x &
                   + rij%y*rij%y &
                   + rij%z*rij%z &
                   + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * ep_j(j)%mass
            r3_inv = r3_inv * r_inv
            ai%x   = ai%x - r3_inv * rij%x
            ai%y   = ai%y - r3_inv * rij%y
            ai%z   = ai%z - r3_inv * rij%z
            poti   = poti - r_inv
         end do
         f(i)%pot = f(i)%pot + poti
         f(i)%acc = f(i)%acc + ai
      end do

   end subroutine calc_gravity_ep_sp
#endif

end module user_defined_types
