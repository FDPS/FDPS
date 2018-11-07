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
