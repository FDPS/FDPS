!================================
!   MODULE: User defined types
!================================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   implicit none

   !* Private parameters
   real(kind=c_double), parameter, private :: pi=datan(1.0d0)*4.0d0
   !* Public parameters
   real(kind=c_double), parameter, public :: kernel_support_radius=2.5d0

   !**** Force types
   type, public, bind(c) :: force_dens !$fdps Force
      !$fdps clear smth=keep
      real(kind=c_double) :: dens
      real(kind=c_double) :: smth
   end type force_dens

   type, public, bind(c) :: force_hydro !$fdps Force
      !$fdps clear 
      type(fdps_f64vec) :: acc
      real(kind=c_double) :: eng_dot
      real(kind=c_double) :: dt
   end type force_hydro

   !**** Full particle type
   type, public, bind(c) :: full_particle !$fdps FP
      !$fdps copyFromForce force_dens (dens,dens)
      !$fdps copyFromForce force_hydro (acc,acc) (eng_dot,eng_dot) (dt,dt)
      real(kind=c_double) :: mass !$fdps charge
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel 
      type(fdps_f64vec) :: acc
      real(kind=c_double) :: dens
      real(kind=c_double) :: eng
      real(kind=c_double) :: pres
      real(kind=c_double) :: smth !$fdps rsearch
      real(kind=c_double) :: snds
      real(kind=c_double) :: eng_dot
      real(kind=c_double) :: dt
      integer(kind=c_long_long) :: id
      type(fdps_f64vec) :: vel_half
      real(kind=c_double) :: eng_half
   end type full_particle

   !**** Essential particle type
   type, public, bind(c) :: essential_particle !$fdps EPI,EPJ
      !$fdps copyFromFP full_particle (id,id) (pos,pos) (vel,vel) (mass,mass) (smth,smth) (dens,dens) (pres,pres) (snds,snds)
      integer(kind=c_long_long) :: id !$fdps id
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel
      real(kind=c_double) :: mass !$fdps charge
      real(kind=c_double) :: smth !$fdps rsearch
      real(kind=c_double) :: dens
      real(kind=c_double) :: pres
      real(kind=c_double) :: snds
   end type essential_particle

   !* Public routines
   public :: W
   public :: gradW
   public :: calc_density
   public :: calc_hydro_force

   contains

   !-------------------------------------------------------------------
   pure function W(dr,h)
      implicit none
      real(kind=c_double) :: W
      type(fdps_f64vec), intent(in) :: dr
      real(kind=c_double), intent(in) :: h
      !* Local variables
      real(kind=c_double) :: s,s1,s2
      
      s = dsqrt(dr%x*dr%x &
               +dr%y*dr%y &
               +dr%z*dr%z)/h
      s1 = 1.0d0 - s
      if (s1 < 0.0d0) s1 = 0.0d0
      s2 = 0.5d0 - s
      if (s2 < 0.0d0) s2 = 0.0d0
      W = (s1*s1*s1) - 4.0d0*(s2*s2*s2)
      W = W * 16.0d0/(pi*h*h*h)

   end function W

   !-------------------------------------------------------------------
   pure function gradW(dr,h)
      implicit none
      type(fdps_f64vec) :: gradW
      type(fdps_f64vec), intent(in) :: dr
      real(kind=c_double), intent(in) :: h
      !* Local variables
      real(kind=c_double) :: dr_abs,s,s1,s2,coef

      dr_abs = dsqrt(dr%x*dr%x &
                    +dr%y*dr%y &
                    +dr%z*dr%z)
      s = dr_abs/h
      s1 = 1.0d0 - s
      if (s1 < 0.0d0) s1 = 0.0d0
      s2 = 0.5d0 - s
      if (s2 < 0.0d0) s2 = 0.0d0
      coef = - 3.0d0*(s1*s1) + 12.0d0*(s2*s2)
      coef = coef * 16.0d0/(pi*h*h*h)
      coef = coef / (dr_abs*h + 1.0d-6*h)
      gradW%x = dr%x * coef
      gradW%y = dr%y * coef
      gradW%z = dr%z * coef

   end function gradW

   !**** Interaction function
   subroutine calc_density(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(kind=c_int), intent(in), value :: n_ip,n_jp
      type(essential_particle), dimension(n_ip), intent(in) :: ep_i
      type(essential_particle), dimension(n_jp), intent(in) :: ep_j
      type(force_dens), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(kind=c_int) :: i,j
      type(fdps_f64vec) :: dr

      do i=1,n_ip
         f(i)%dens = 0.0d0
         do j=1,n_jp
            dr%x = ep_j(j)%pos%x - ep_i(i)%pos%x
            dr%y = ep_j(j)%pos%y - ep_i(i)%pos%y
            dr%z = ep_j(j)%pos%z - ep_i(i)%pos%z
            f(i)%dens = f(i)%dens &
                      + ep_j(j)%mass * W(dr,ep_i(i)%smth)
         end do
      end do

   end subroutine calc_density

   !**** Interaction function
   subroutine calc_hydro_force(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(kind=c_int), intent(in), value :: n_ip,n_jp
      type(essential_particle), dimension(n_ip), intent(in) :: ep_i
      type(essential_particle), dimension(n_jp), intent(in) :: ep_j
      type(force_hydro), dimension(n_ip), intent(inout) :: f
      !* Local parameters
      real(kind=c_double), parameter :: C_CFL=0.3d0
      !* Local variables
      integer(kind=c_int) :: i,j
      real(kind=c_double) :: mass_i,mass_j,smth_i,smth_j, &
                             dens_i,dens_j,pres_i,pres_j, &
                             snds_i,snds_j
      real(kind=c_double) :: povrho2_i,povrho2_j, &
                             v_sig_max,dr_dv,w_ij,v_sig,AV
      type(fdps_f64vec) :: pos_i,pos_j,vel_i,vel_j, &
                           dr,dv,gradW_ij

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
            AV = - 0.5d0*v_sig*w_ij / (0.5d0*(dens_i+dens_j))
            !* Compute the average of the gradients of kernel
            gradW_ij = 0.5d0 * (gradW(dr,smth_i) + gradW(dr,smth_j))
            !* Compute the acceleration and the heating rate
            f(i)%acc%x = f(i)%acc%x - mass_j*(povrho2_i+povrho2_j+AV)*gradW_ij%x
            f(i)%acc%y = f(i)%acc%y - mass_j*(povrho2_i+povrho2_j+AV)*gradW_ij%y
            f(i)%acc%z = f(i)%acc%z - mass_j*(povrho2_i+povrho2_j+AV)*gradW_ij%z
            f(i)%eng_dot = f(i)%eng_dot &
                         + mass_j * (povrho2_i + 0.5d0*AV) &
                          *(dv%x * gradW_ij%x &
                           +dv%y * gradW_ij%y &
                           +dv%z * gradW_ij%z)
         end do
         f(i)%dt = C_CFL*2.0d0*smth_i/(v_sig_max*kernel_support_radius)
      end do
      ! [IMPORTANT NOTE]
      !   In the innermost loop, we use the components of vectors
      !   directly for vector operations because of the following
      !   reasion. Except for intel compilers with `-ipo` option,
      !   most of Fortran compilers use function calls to perform
      !   vector operations like rij = x - ep_j(j)%pos.
      !   This significantly slow downs the speed of the code.
      !   By using the components of vector directly, we can avoid 
      !   these function calls.

   end subroutine calc_hydro_force

end module user_defined_types

