#include "macro_defs.h"
!-----------------------------------------------------------------------
!/////////////////////     S U B R O U T I N E     /////////////////////
!///////////////////// < I N I T I A L _ K I C K > /////////////////////
!-----------------------------------------------------------------------
subroutine initial_kick(psys_num_nbody,psys_num_sph,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   double precision, intent(in) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(fdps_f64vec) :: acc_tot
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph

   !* Update N-body system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   do i=1,nptcl_loc
      ptcl_nbody(i)%vel%x = ptcl_nbody(i)%vel%x + 0.5d0 * dt * ptcl_nbody(i)%acc%x
      ptcl_nbody(i)%vel%y = ptcl_nbody(i)%vel%y + 0.5d0 * dt * ptcl_nbody(i)%acc%y
      ptcl_nbody(i)%vel%z = ptcl_nbody(i)%vel%z + 0.5d0 * dt * ptcl_nbody(i)%acc%z
   end do
   nullify(ptcl_nbody)

   !* Update SPH system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   do i=1,nptcl_loc
      acc_tot%x = ptcl_sph(i)%acc_grav%x + ptcl_sph(i)%acc_hydro%x
      acc_tot%y = ptcl_sph(i)%acc_grav%y + ptcl_sph(i)%acc_hydro%y
      acc_tot%z = ptcl_sph(i)%acc_grav%z + ptcl_sph(i)%acc_hydro%z
      ptcl_sph(i)%vel_half%x = ptcl_sph(i)%vel%x + 0.5d0 * dt * acc_tot%x
      ptcl_sph(i)%vel_half%y = ptcl_sph(i)%vel%y + 0.5d0 * dt * acc_tot%y
      ptcl_sph(i)%vel_half%z = ptcl_sph(i)%vel%z + 0.5d0 * dt * acc_tot%z
#if !defined(ISOTHERMAL_EOS)
      ptcl_sph(i)%eng_half = ptcl_sph(i)%eng + 0.5d0 * dt * ptcl_sph(i)%eng_dot
      ptcl_sph(i)%ent_half = ptcl_sph(i)%ent + 0.5d0 * dt * ptcl_sph(i)%ent_dot
#endif
   end do
   nullify(ptcl_sph)

end subroutine initial_kick

!-----------------------------------------------------------------------
!///////////////////////   S U B R O U T I N E   ///////////////////////
!/////////////////////// < F U L L _ D R I F T > ///////////////////////
!-----------------------------------------------------------------------
subroutine full_drift(psys_num_nbody,psys_num_sph,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   double precision, intent(in) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph

   !* Update N-body system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   do i=1,nptcl_loc
      ptcl_nbody(i)%pos%x = ptcl_nbody(i)%pos%x + dt * ptcl_nbody(i)%vel%x
      ptcl_nbody(i)%pos%y = ptcl_nbody(i)%pos%y + dt * ptcl_nbody(i)%vel%y
      ptcl_nbody(i)%pos%z = ptcl_nbody(i)%pos%z + dt * ptcl_nbody(i)%vel%z
   end do
   nullify(ptcl_nbody)

   !* Update SPH system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   do i=1,nptcl_loc
      ptcl_sph(i)%pos%x = ptcl_sph(i)%pos%x + dt * ptcl_sph(i)%vel_half%x
      ptcl_sph(i)%pos%y = ptcl_sph(i)%pos%y + dt * ptcl_sph(i)%vel_half%y
      ptcl_sph(i)%pos%z = ptcl_sph(i)%pos%z + dt * ptcl_sph(i)%vel_half%z
   end do
   nullify(ptcl_sph)

end subroutine full_drift

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////  < P R E D I C T >  /////////////////////////
!-----------------------------------------------------------------------
subroutine predict(psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num
   double precision, intent(in) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(fdps_f64vec) :: acc_tot
   type(fdps_controller) :: fdps_ctrl
   type(fp_sph), dimension(:), pointer :: ptcl

   !* Update SPH system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      acc_tot%x = ptcl(i)%acc_grav%x + ptcl(i)%acc_hydro%x
      acc_tot%y = ptcl(i)%acc_grav%y + ptcl(i)%acc_hydro%y
      acc_tot%z = ptcl(i)%acc_grav%z + ptcl(i)%acc_hydro%z
      ptcl(i)%vel%x = ptcl(i)%vel%x + dt * acc_tot%x
      ptcl(i)%vel%y = ptcl(i)%vel%y + dt * acc_tot%y
      ptcl(i)%vel%z = ptcl(i)%vel%z + dt * acc_tot%z
#if !defined(ISOTHERMAL_EOS)
      ptcl(i)%eng = ptcl(i)%eng + dt * ptcl(i)%eng_dot
      ptcl(i)%ent = ptcl(i)%ent + dt * ptcl(i)%ent_dot
#endif
   end do
   nullify(ptcl)

end subroutine predict

!-----------------------------------------------------------------------
!///////////////////////   S U B R O U T I N E   ///////////////////////
!/////////////////////// < F I N A L _ K I C K > ///////////////////////
!-----------------------------------------------------------------------
subroutine final_kick(psys_num_nbody,psys_num_sph,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   double precision, intent(in) :: dt
   !* Local parameters
   double precision, parameter :: coeff = 0.1d0
   !* Local variables
   integer :: i,nptcl_loc
   double precision :: frac_dump
   type(fdps_f64vec) :: acc_tot
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph

   !* Update N-body system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   do i=1,nptcl_loc
      ptcl_nbody(i)%vel%x = ptcl_nbody(i)%vel%x + 0.5d0 * dt * ptcl_nbody(i)%acc%x
      ptcl_nbody(i)%vel%y = ptcl_nbody(i)%vel%y + 0.5d0 * dt * ptcl_nbody(i)%acc%y
      ptcl_nbody(i)%vel%z = ptcl_nbody(i)%vel%z + 0.5d0 * dt * ptcl_nbody(i)%acc%z
   end do
   nullify(ptcl_nbody)

   !* Update SPH system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   do i=1,nptcl_loc
      acc_tot%x = ptcl_sph(i)%acc_grav%x + ptcl_sph(i)%acc_hydro%x
      acc_tot%y = ptcl_sph(i)%acc_grav%y + ptcl_sph(i)%acc_hydro%y
      acc_tot%z = ptcl_sph(i)%acc_grav%z + ptcl_sph(i)%acc_hydro%z
      ptcl_sph(i)%vel%x = ptcl_sph(i)%vel_half%x + 0.5d0 * dt * acc_tot%x
      ptcl_sph(i)%vel%y = ptcl_sph(i)%vel_half%y + 0.5d0 * dt * acc_tot%y
      ptcl_sph(i)%vel%z = ptcl_sph(i)%vel_half%z + 0.5d0 * dt * acc_tot%z
#if !defined(ISOTHERMAL_EOS)
      ptcl_sph(i)%eng = ptcl_sph(i)%eng_half + 0.5d0 * dt * ptcl_sph(i)%eng_dot
      ptcl_sph(i)%ent = ptcl_sph(i)%ent_half + 0.5d0 * dt * ptcl_sph(i)%ent_dot
#endif
   end do
#if defined(DUMP_VELOCITY_OF_SPH_PARTICLE)
   do i=1,nptcl_loc
      frac_dump = exp(- coeff * (CFL_hydro/0.1d0) * ptcl_sph(i)%snds * dt / ptcl_sph(i)%smth)
      ptcl_sph(i)%vel%x = ptcl_sph(i)%vel%x * frac_dump
      ptcl_sph(i)%vel%y = ptcl_sph(i)%vel%y * frac_dump
      ptcl_sph(i)%vel%z = ptcl_sph(i)%vel%z * frac_dump
   end do
#endif

   nullify(ptcl_sph)

end subroutine final_kick

