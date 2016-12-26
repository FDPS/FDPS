!=================================
!   MODULE: User defined types 
!=================================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   use fdps_super_particle
   implicit none

   !**** Force class
   type, public, bind(c) :: nbody_pp_results !$fdps Force
      !$fdps clear
      real(kind=c_double) :: pot
      type(fdps_f64vec) :: agrv
   end type nbody_pp_results

   !**** Full Particle Class
   type, public, bind(c) :: nbody_fp !$fdps FP
      !$fdps copyFromForce nbody_pp_results (pot,pot) (agrv,agrv)
      !$fdps copyFromForcePM agrv_pm
      integer(kind=c_long_long) :: id
      real(kind=c_double) :: m !$fdps charge
      real(kind=c_double) :: rc !$fdps rsearch
      type(fdps_f64vec) :: x !$fdps position
      type(fdps_f64vec) :: v,v_half
      type(fdps_f64vec) :: agrv
      real(kind=c_double) :: pot
      type(fdps_f32vec) :: agrv_pm
      real(kind=c_float) :: pot_pm
   end type nbody_fp

   !**** Essential Particle Class
   type, public, bind(c) :: nbody_ep !$fdps EPI,EPJ
      !$fdps copyFromFP nbody_fp (id,id) (m,m) (rc,rc) (x,x)
      integer(kind=c_long_long) :: id
      real(kind=c_double) :: m !$fdps charge
      real(kind=c_double) :: rc !$fdps rsearch
      type(fdps_f64vec) :: x !$fdps position
   end type nbody_ep

   !**** Crystal Parameters class
   type, public, bind(c) :: crystal_parameters
      integer(kind=c_int) :: nptcl_per_side
      type(fdps_f64vec) :: pos_vertex
   end type crystal_parameters

   !* Public routines
   public :: S2_pcut
   public :: S2_fcut
   public :: calc_force_ep_ep
   public :: calc_force_ep_sp

   contains

   !-------------------------------------------------------------------- 
   pure function S2_pcut(xi)
      ! This is the potential cutoff function where we used Eq.(8.75)
      ! in Hockney & Eastwood (1987).
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double) :: S2_pcut
      real(kind=c_double), intent(in) :: xi
      if (xi <= 1.0d0) then
         S2_pcut = 1.0d0             &
                 - xi*(208.0d0       &
                 + (xi*xi)*(-112.0d0 &
                 + (xi*xi)*(56.0d0   &
                 + xi*(-14.0d0       &
                 + xi*(-8.0d0        &
                 + 3.0d0*xi)))))/140.0d0
      else if ((1.0d0 < xi) .and. (xi < 2.0d0)) then
         S2_pcut = 1.0d0        &
                 - (12.0d0      &
                 + xi*(128.0d0  &
                 + xi*(224.0d0  &
                 + xi*(-448.0d0 &
                 + xi*(280.0d0  &
                 + xi*(-56.0d0  &
                 + xi*(-14.0d0  &
                 + xi*(8.0d0    &
                 - xi))))))))/140.0d0
      else
         S2_pcut = 0.0d0
      end if
   end function S2_pcut

   !-------------------------------------------------------------------- 
   pure function S2_fcut(xi)
      ! This function returns 1 - R(\xi), where \xi is r/(a/2), a is the
      ! scale length of the cutoff function, and R(\xi) is almost the same
      ! as the function defined as Eq.(8-72) in Hockney & Eastwood (1987).
      ! The only difference is that [1/(r/(a/2))]^2 is factored out 
      ! in this function from Eq.(8-72).
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double) :: S2_fcut
      real(kind=c_double), intent(in) :: xi
      if (xi <= 1.0d0) then
         S2_fcut = 1.0d0               &
                 - (xi*xi*xi)*(224.0d0 &
                 + (xi*xi)*(-224.0d0   &
                 + xi*(70.0d0          &
                 + xi*(48.0d0          &
                 - 21.0d0*xi))))/140.0d0
      else if ((1.0d0 < xi) .and. (xi < 2.0d0)) then
         S2_fcut = 1.0d0             &
                 - (12.0d0           &
                 + (xi*xi)*(-224.0d0 &
                 + xi*(896.0d0       &
                 + xi*(-840.0d0      &
                 + xi*(224.0d0       &
                 + xi*(70.0d0        &
                 + xi*(-48.0d0       &
                 + 7.0d0*xi)))))))/140.0d0
      else
         S2_fcut = 0.0d0
      end if
   end function S2_fcut
  
   !-------------------------------------------------------------------- 
   subroutine calc_force_ep_ep(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(nbody_ep), dimension(n_ip), intent(in) :: ep_i
      type(nbody_ep), dimension(n_jp), intent(in) :: ep_j
      type(nbody_pp_results), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      real(c_double) :: rij,rinv,rinv3,xi
      type(fdps_f64vec) :: dx

      do i=1,n_ip
         do j=1,n_jp
            dx%x = ep_i(i)%x%x - ep_j(j)%x%x
            dx%y = ep_i(i)%x%y - ep_j(j)%x%y
            dx%z = ep_i(i)%x%z - ep_j(j)%x%z
            rij  = dsqrt(dx%x * dx%x &
                        +dx%y * dx%y &
                        +dx%z * dx%z)
            if ((ep_i(i)%id == ep_j(j)%id) .and. (rij == 0.0d0)) cycle
            rinv = 1.0d0/rij
            rinv3 = rinv*rinv*rinv
            xi = 2.0d0*rij/ep_i(i)%rc
            f(i)%pot    = f(i)%pot    + ep_j(j)%m * S2_pcut(xi) * rinv
            f(i)%agrv%x = f(i)%agrv%x + ep_j(j)%m * S2_fcut(xi) * rinv3 * dx%x
            f(i)%agrv%y = f(i)%agrv%y + ep_j(j)%m * S2_fcut(xi) * rinv3 * dx%y
            f(i)%agrv%z = f(i)%agrv%z + ep_j(j)%m * S2_fcut(xi) * rinv3 * dx%z
         end do
         !* Self-interaction term
         f(i)%pot = f(i)%pot - ep_i(i)%m * (208.0d0/(70.0d0*ep_i(i)%rc))
      end do

   end subroutine calc_force_ep_ep

   !-------------------------------------------------------------------- 
   subroutine calc_force_ep_sp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(nbody_ep), dimension(n_ip), intent(in) :: ep_i
      type(fdps_spj_monopole_cutoff), dimension(n_jp), intent(in) :: ep_j
      type(nbody_pp_results), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      real(c_double) :: rij,rinv,rinv3,xi
      type(fdps_f64vec) :: dx

      do i=1,n_ip
         do j=1,n_jp
            dx%x = ep_i(i)%x%x - ep_j(j)%pos%x
            dx%y = ep_i(i)%x%y - ep_j(j)%pos%y
            dx%z = ep_i(i)%x%z - ep_j(j)%pos%z
            rij  = dsqrt(dx%x * dx%x &
                        +dx%y * dx%y &
                        +dx%z * dx%z)
            rinv = 1.0d0/rij
            rinv3 = rinv*rinv*rinv
            xi = 2.0d0*rij/ep_i(i)%rc
            f(i)%pot    = f(i)%pot    + ep_j(j)%mass * S2_pcut(xi) * rinv
            f(i)%agrv%x = f(i)%agrv%x + ep_j(j)%mass * S2_fcut(xi) * rinv3 * dx%x
            f(i)%agrv%y = f(i)%agrv%y + ep_j(j)%mass * S2_fcut(xi) * rinv3 * dx%y
            f(i)%agrv%z = f(i)%agrv%z + ep_j(j)%mass * S2_fcut(xi) * rinv3 * dx%z
         end do
      end do

   end subroutine calc_force_ep_sp

end module user_defined_types
