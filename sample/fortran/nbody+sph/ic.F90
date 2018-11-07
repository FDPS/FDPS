#include "macro_defs.h"

!=============================
! MODULE: Tipsy file reader
!=============================
module tipsy_file_reader
   use, intrinsic :: iso_c_binding
   implicit none

   ! The following two classes are used in function readTipsyFile,
   ! which is used to read particle data created by MAGI.
   type, public, bind(c):: magi_tipsy_header
      real(kind=c_double) :: time
      integer(kind=c_int) :: nbodies
      integer(kind=c_int) :: ndim
      integer(kind=c_int) :: nsph
      integer(kind=c_int) :: ndark
      integer(kind=c_int) :: nstar
   end type magi_tipsy_header
  
   type, public, bind(c) :: magi_tipsy_particle 
      real(kind=c_float) :: mass 
      real(kind=c_float) :: pos(3)
      real(kind=c_float) :: vel(3)
      real(kind=c_float) :: eps
      integer(kind=c_int) :: idx
   end type magi_tipsy_particle

   ! Public routine
   public :: read_tipsy_file

   contains

   subroutine read_tipsy_file(file_name, psys_num)
      ! This function is used to read particle data created by
      ! MAGI (https://bitbucket.org/ymiki/magi). The particle
      ! data must be in the TIPSY format.
      use fdps_module
      use user_defined_types
      implicit none
      character(len=*), intent(in) :: file_name
      integer, intent(in) :: psys_num
      !* Local variables
      integer :: i
      type(magi_tipsy_header) :: header
      type(magi_tipsy_particle) :: ptcl_tipsy
      type(fdps_controller) :: fdps_ctrl
      type(fp_nbody), dimension(:), pointer :: ptcl

      ! Read file
      open(unit=9,file=trim(file_name),action='read', &
           form='unformatted',access='stream',status='old')
         read(9)header
         write(*,*)'nbodies = ',header%nbodies
         call fdps_ctrl%set_nptcl_loc(psys_num,header%nbodies)
         call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
         do i=1,header%nbodies
            read(9)ptcl_tipsy
            ptcl(i)%mass  = ptcl_tipsy%mass
            ptcl(i)%pos%x = ptcl_tipsy%pos(1)
            ptcl(i)%pos%y = ptcl_tipsy%pos(2)
            ptcl(i)%pos%z = ptcl_tipsy%pos(3)
            ptcl(i)%vel%x = ptcl_tipsy%vel(1)
            ptcl(i)%vel%y = ptcl_tipsy%vel(2)
            ptcl(i)%vel%z = ptcl_tipsy%vel(3)
         end do
      close(unit=9)

   end subroutine read_tipsy_file

end module tipsy_file_reader

!-----------------------------------------------------------------------
!///////////////////////   S U B R O U T I N E   ///////////////////////
!///////////////////////  < G A L A X Y _ I C >  ///////////////////////
!-----------------------------------------------------------------------
subroutine galaxy_IC(psys_num_nbody,psys_num_sph, &
                     bc,                          &
                     pos_root_domain_low,         &
                     pos_root_domain_high,        &
                     time_dump,dt_dump,time_end)
   use mathematical_constants
   use physical_constants
   use fdps_vector
   use fdps_module
   use user_defined_types
   use tipsy_file_reader
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   integer(kind=c_int), intent(inout) :: bc
   type(fdps_f64vec), intent(inout) :: pos_root_domain_low, &
                                       pos_root_domain_high
   double precision, intent(inout) :: time_dump,dt_dump,time_end
   !* Local parameters
   !- Definitions of the code units of MAGI
   !    [Important]
   !    (1) The values MUST BE consistent with "the computational units"
   !        written in the file ./magi_data/doc/unit.txt, which is output
   !        by MAGI when we create a particle data with MAGI.
   !    (2) The MAGI's code units are DIFFERENT for unit systems
   !        a user choose in the file ./magi_data/cfg/Galaxy.tipsy.
   !        For detail, read Section "Unit systems in inc/constants.[c h]"
   !        in https://bitbucket.org/ymiki/magi.
   !    (3) In this sample code, "Galactic scale" unit is adopted.
   !        It is consistent with ./magi_data/cfg/Galaxy.cfg.
   double precision, parameter :: magi_unit_mass = 1.0d8 * Msolar
   double precision, parameter :: magi_unit_leng = kpc
   double precision, parameter :: magi_unit_time = 1.0d2 * Myr
   double precision, parameter :: magi_unit_velc = magi_unit_leng/magi_unit_time
   !- Definitions of the model parameters for a gaseous exponential disk
   integer, parameter :: nptcl_sph = 2**18
   double precision, parameter :: m_gas = 1.0d10 * Msolar
   double precision, parameter :: rs = 7.0d0 * kpc ! scale radius
   double precision, parameter :: rt = 12.5d0 * kpc ! truncation radius
   double precision, parameter :: zd = 4.0d2 * pc ! scale height
   double precision, parameter :: zt = 1.0d0 * kpc ! truncation height
   double precision, parameter :: temp = 1.0d4 ! gas temperature 
   double precision, parameter :: mu = 0.5d0 ! mean molecular weight relative to the mass of hydrogen
   double precision, parameter :: eps = 1.0d-6 ! iteraction accuracy 
   !- Definitions of parameters for output
   integer, parameter :: nbin=64
   double precision, parameter :: safety=1.001d0
   !* Local variables
   integer :: i,indx
   integer :: nptcl_nbody
   integer :: mtts_num
   double precision :: x,y,z,r2,r
   double precision :: r_low,r_high,r_new
   double precision :: val,val_trgt
   double precision :: reldiff
   double precision :: z_new
   double precision :: unit_mass,unit_leng,unit_time,unit_velc,unit_eng
   double precision :: rmax,dr,zmax,dz
   double precision, dimension(nbin) :: sigma
   integer, dimension(nbin) :: dist_func
   character(len=64) :: file_name
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph

   ! Initialize pseudorandom number generator
   call fdps_ctrl%create_mtts(mtts_num)
   call fdps_ctrl%mtts_init_genrand(mtts_num,0)
   ! Place Nbody particles
   file_name = "./magi_data/dat/Galaxy.tipsy"
   call read_tipsy_file(file_name, psys_num_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   nptcl_nbody = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   do i=1,nptcl_nbody
      ptcl_nbody(i)%id  = i
      ptcl_nbody(i)%acc = 0.0d0
      ptcl_nbody(i)%pot = 0.0d0
   end do
   ! Place SPH particles
   call fdps_ctrl%set_nptcl_loc(psys_num_sph,nptcl_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   do i=1,nptcl_sph
      ! First make a uniform disk with a finite thickness
      do 
         x = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * rt
         y = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * rt
         z = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * zt
         r2 = (x * x + y * y)
         if ((r2 < rt * rt) .and. (dabs(z) < zt)) exit
      end do
      r = dsqrt(r2)
      ! Then re-scale particle position to generate an exponential disk
      r_low = 0.0d0
      r_high = rt
      val_trgt = (r/rt) * (r/rt)
      do
         r_new = 0.5d0 * (r_low + r_high)
         val = (1.0d0 - (r_new/rs + 1.0d0)*exp(-r_new/rs)) & 
             / (1.0d0 - (   rt/rs + 1.0d0)*exp(-   rt/rs))
         if (val < val_trgt) r_low  = r_new
         if (val > val_trgt) r_high = r_new
         reldiff = 2.0d0 * dabs(r_low - r_high)/(r_low + r_high)
         if (reldiff < eps) then
            r_new = 0.5d0 * (r_low + r_high)
            exit
         end if
      end do
      if (z >= 0.0) then
         z_new = - zd * dlog(1.0d0 - (z/zt) * (1.0d0 - exp(-zt/zd)))
      else
         z_new = zd * dlog(1.0d0 + (z/zt) * (1.0d0 - exp(-zt/zd)))
      end if
      ! Set  
      ptcl_sph(i)%id   = i + nptcl_nbody
      ptcl_sph(i)%mass = m_gas / nptcl_sph
      ptcl_sph(i)%pos%x = (r_new / r) * x
      ptcl_sph(i)%pos%y = (r_new / r) * y
      ptcl_sph(i)%pos%z = z_new
      ptcl_sph(i)%vel       = 0.0d0
      ptcl_sph(i)%acc_grav  = 0.0d0
      ptcl_sph(i)%pot_grav  = 0.0d0
      ptcl_sph(i)%acc_hydro = 0.0d0
      ptcl_sph(i)%eng = (kBoltz * temp)/((specific_heat_ratio - 1.0) * mu * Mhydrogen)
      ptcl_sph(i)%smth = (Rt*Rt*zt)**(1.0d0/3.0d0) * (dble(N_neighbor)/dble(nptcl_sph))**(1.0d0/3.0d0)
   end do
   ! Unit convertion (MAGI unit, CGS unit -> G=M=R=1 system)
   unit_mass = 0.0d0
   unit_leng = 0.0d0
   do i=1,nptcl_nbody
      ptcl_nbody(i)%mass = ptcl_nbody(i)%mass * magi_unit_mass
      ptcl_nbody(i)%pos  = ptcl_nbody(i)%pos  * magi_unit_leng
      ptcl_nbody(i)%vel  = ptcl_nbody(i)%vel  * magi_unit_velc
      unit_mass = unit_mass + ptcl_nbody(i)%mass
      r = dsqrt(ptcl_nbody(i)%pos * ptcl_nbody(i)%pos)
      if (r > unit_leng) unit_leng = r
   end do
   write(*,10)"Total mass in N-body particles = ",unit_mass/Msolar," [Msolar]"
   10 format(a,1es25.16e3,a)
   do i=1,nptcl_sph
      unit_mass = unit_mass + ptcl_sph(i)%mass
      r = dsqrt(ptcl_sph(i)%pos * ptcl_sph(i)%pos)
      if (r > unit_leng) unit_leng = r
   end do
   unit_time = dsqrt(unit_leng**3.0d0/(Ggrav * unit_mass))
   unit_velc = unit_leng / unit_time
   unit_eng  = unit_mass * unit_velc * unit_velc
   write(*,20)"unit_mass = ",unit_mass," [g]    = ",unit_mass/Msolar," [Msolar]" 
   write(*,20)"unit_leng = ",unit_leng," [cm]   = ",unit_leng/kpc," [kpc]"
   write(*,20)"unit_time = ",unit_time," [s]    = ",unit_time/Gyr," [Gyr]"
   write(*,20)"unit_velc = ",unit_velc," [cm/s] = ",unit_velc/km," [km/s]"
   20 format(a,1es25.16e3,a,1es25.16e3,a)
   do i=1,nptcl_nbody
      ptcl_nbody(i)%mass = ptcl_nbody(i)%mass / unit_mass
      ptcl_nbody(i)%pos  = ptcl_nbody(i)%pos  / unit_leng
      ptcl_nbody(i)%vel  = ptcl_nbody(i)%vel  / unit_velc
   end do
   do i=1,nptcl_sph
      ptcl_sph(i)%mass = ptcl_sph(i)%mass / unit_mass
      ptcl_sph(i)%pos  = ptcl_sph(i)%pos  / unit_leng
      ptcl_sph(i)%vel  = ptcl_sph(i)%vel  / unit_velc
      ptcl_sph(i)%smth = ptcl_sph(i)%smth / unit_leng
      ptcl_sph(i)%eng  = ptcl_sph(i)%eng  / (unit_eng/unit_mass)
   end do
   ! Set boundary condition
   bc = fdps_bc_open
   ! Set gravitational softening
   eps_grav = 1.0d-3
   ! Set I/O intervals
   dt_dump = 0.01d0
   time_dump = dt_dump
   time_end = 1.0d0
   ! Set maximum timestep
   dt_max = 1.0d-4
   ! Output
   write(*,*)"An initial condition for isolated galaxy simulation is made."
   write(*,*)"N_nbody = ",nptcl_nbody,", N_sph = ",nptcl_sph

   ! In the following, we compute the surface gas density and
   ! the distribution function of particle's z coordinates.
   ! Compute the maxium clyndrical radius of gas particles
   rmax = 0.0
   do i=1,nptcl_sph
      r = dsqrt(ptcl_sph(i)%pos%x * ptcl_sph(i)%pos%x &
               +ptcl_sph(i)%pos%y * ptcl_sph(i)%pos%y)
      if (r > rmax) rmax = r
   end do
   ! Compute the surface density
   dr = (safety * rmax)/nbin
   do i=1,nptcl_sph
      r = dsqrt(ptcl_sph(i)%pos%x * ptcl_sph(i)%pos%x &
               +ptcl_sph(i)%pos%y * ptcl_sph(i)%pos%y)
      indx = r/dr + 1
      sigma(indx) = sigma(indx) + ptcl_sph(i)%mass
   end do
   do i=1,nbin
      sigma(i) = sigma(i) / (2.0d0 * pi * (2*i - 1) * dr*dr)
   end do
   ! Output the surface density
   file_name = "Sigma.txt"
   open(unit=9,file=trim(file_name),action='write',status='replace')
      do i=1,nbin
         write(9,100)(i - 0.5d0) * dr, sigma(i)
      end do
      100 format(1es25.16e3,1x,1es25.16e3)
   close(unit=9)
   ! Compute the distribution function 
   zmax = zt/unit_leng
   dz = (safety * 2.0d0 * zmax)/nbin
   dist_func = 0
   do i=1,nptcl_sph
      indx = (ptcl_sph(i)%pos%z + safety * zmax)/dz + 1
      dist_func(indx) = dist_func(indx) + 1
   end do
   ! Output the distribution function
   file_name = "dist_func.txt"
   open(unit=9,file=trim(file_name),action='write',status='replace')
      do i=1,nbin
         z = (i - 0.5) * dz - safety * zmax
         write(9,200)z,dist_func(i)
      end do
      200 format(1es25.16e3,1x,i10)
   close(unit=9)
   
end subroutine galaxy_IC

!-----------------------------------------------------------------------
!////////////              S U B R O U T I N E              ////////////
!//////////// < C O L D _ C O L L A P S E _ T E S T _ I C > ////////////
!-----------------------------------------------------------------------
subroutine cold_collapse_test_IC(psys_num_nbody,psys_num_sph, &
                                 bc,                          &
                                 pos_root_domain_low,         &
                                 pos_root_domain_high,        &
                                 time_dump,dt_dump,time_end)
   use mathematical_constants
   use physical_constants
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   integer(kind=c_int), intent(inout) :: bc
   type(fdps_f64vec), intent(inout) :: pos_root_domain_low, &
                                       pos_root_domain_high
   double precision, intent(inout) :: time_dump,dt_dump,time_end
   !* Local parameters
   double precision, parameter :: M_nbody = 0.5d0
   double precision, parameter :: M_sph   = 0.5d0
   integer, parameter :: N_nbody = 512
   integer, parameter :: N_sph   = 512
   double precision, parameter :: R_nbody = 3.0d0
   double precision, parameter :: R_sph   = 3.0d0
   double precision, parameter :: E_bind_nbody = - 3.0d0 * M_nbody * M_nbody / 5.0d0 / R_nbody
   double precision, parameter :: E_bind_sph   = - 3.0d0 * M_sph * M_sph / 5.0d0 / R_sph
   double precision, parameter :: virial_ratio = 0.5d0
   !* Local variables
   integer :: i,nptcl_loc
   integer :: mtts_num
   double precision :: r2
   double precision :: dens,p
   double precision :: dens_nbody,dens_sph,unit_dens,unit_time
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph

   ! Initialize pseudorandom number generator
   call fdps_ctrl%create_mtts(mtts_num)
   call fdps_ctrl%mtts_init_genrand(mtts_num,0)
   ! Place Nbody particles
   call fdps_ctrl%set_nptcl_loc(psys_num_nbody,N_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   do i=1,N_nbody
       ptcl_nbody(i)%id = i
       ptcl_nbody(i)%mass = M_nbody / N_nbody
       do
           ptcl_nbody(i)%pos%x = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * R_nbody
           ptcl_nbody(i)%pos%y = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * R_nbody
           ptcl_nbody(i)%pos%z = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * R_nbody
           r2 = ptcl_nbody(i)%pos * ptcl_nbody(i)%pos
           if (r2 <= R_nbody * R_nbody) exit
       end do
       ptcl_nbody(i)%vel = 0.0d0
       ptcl_nbody(i)%acc = 0.0d0
       ptcl_nbody(i)%pot = 0.0d0
   end do
   ! Place SPH particles
   call fdps_ctrl%set_nptcl_loc(psys_num_sph,N_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   do i=1,N_sph
       ptcl_sph(i)%id = i + N_nbody
       ptcl_sph(i)%mass = M_sph / N_sph
       do
           ptcl_sph(i)%pos%x = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * R_sph
           ptcl_sph(i)%pos%y = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * R_sph
           ptcl_sph(i)%pos%z = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) * R_sph
           r2 = ptcl_sph(i)%pos * ptcl_sph(i)%pos
           if (r2 <= R_sph * R_sph) exit
       end do
       ptcl_sph(i)%vel = 0.0d0
       ptcl_sph(i)%acc_grav  = 0.0d0
       ptcl_sph(i)%pot_grav  = 0.0d0
       ptcl_sph(i)%acc_hydro = 0.0d0
       ptcl_sph(i)%eng = virial_ratio * abs(E_bind_sph) / M_sph  ! specific thermal energy
       dens = M_sph / (4.0d0 * pi * R_sph * R_sph * R_sph / 3.0d0) 
       p = specific_heat_ratio - 1.0d0
       ptcl_sph(i)%ent = p * ptcl_sph(i)%eng / dens**(p)
       ptcl_sph(i)%smth = R_sph * (dble(N_neighbor)/dble(N_sph))**(1.0d0/3.0d0) ! smoothing length
   end do
   ! Set boundary condition
   bc = fdps_bc_open
   ! Set other parameters
   eps_grav = 0.01d0 * R_nbody
   ! Set I/O intervals
   dens_nbody = 3.0d0 * M_nbody / (4.0d0 * pi * R_nbody * R_nbody * R_nbody)
   dens_sph   = 3.0d0 * M_sph / (4.0d0 * pi * R_sph * R_sph * R_sph)
   unit_dens = dens_nbody + dens_sph
   unit_time = dsqrt(3.0d0 * pi/ (32.0d0 * unit_dens))
   write(*,*)"unit_dens = ",unit_dens
   write(*,*)"unit_time = ",unit_time
   dt_dump   = 0.1d0 * unit_time
   time_dump = dt_dump
   time_end  = 2.0 * unit_time
   ! Set the maximum timestep
   dt_max = 1.0d-3

end subroutine cold_collapse_test_IC

!-----------------------------------------------------------------------
!///////////////////       S U B R O U T I N E       ///////////////////
!/////////////////// < E V R A R D _ T E S T _ I C > ///////////////////
!-----------------------------------------------------------------------
subroutine Evrard_test_IC(psys_num_nbody,psys_num_sph, &
                          bc,                          &
                          pos_root_domain_low,         &
                          pos_root_domain_high,        &
                          time_dump,dt_dump,time_end,  &
                          gen_mode)
   use mathematical_constants
   use physical_constants
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   integer(kind=c_int), intent(inout) :: bc
   type(fdps_f64vec), intent(inout) :: pos_root_domain_low, &
                                       pos_root_domain_high
   double precision, intent(inout) :: time_dump,dt_dump,time_end
   integer, intent(in) :: gen_mode
   !* Local parameters
   integer, parameter :: N_nbody=1
   integer, parameter :: N_ptcl_per_side=64
   double precision, parameter :: radius_of_sphere = 1.0d0
   double precision, parameter :: M_sph = 1.0d0
   !* Local variables
   integer :: i,j,k,id
   integer :: nptcl_loc
   integer :: nptcl_in
   integer :: N_sph
   double precision :: x,y,z,r,r2,r_new
   double precision :: dens,p
   double precision :: unit_dens,unit_time
   type(fdps_controller) :: fdps_ctrl
   type(fdps_f64vec), dimension(:), allocatable :: pos
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph
   character(len=64) :: file_name
   !* External routines
   double precision, external :: get_pos_cell_center

   !* Place Nbody particles
   !  here, we place only one dummy, mass-less particle.
   call fdps_ctrl%set_nptcl_loc(psys_num_nbody,N_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   do i=1,N_nbody
       ptcl_nbody(i)%id = i
       ptcl_nbody(i)%mass = 0.0d0
       ptcl_nbody(i)%pos  = 0.0d0
       ptcl_nbody(i)%vel  = 0.0d0
       ptcl_nbody(i)%acc  = 0.0d0
       ptcl_nbody(i)%pot  = 0.0d0
   end do
   !* Place SPH particles
   if (gen_mode == 0) then
      ! In this mode, we create an initial distribution of particles
      ! by rescaling the positions of particles which are placed in a grid.
      !* (1) Count # of particles in a sphere of radius 1
      n_sph = 0
      ! x-loop
      do i=1,n_ptcl_per_side
          x = get_pos_cell_center(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, i)
          ! y-loop
          do j=1,n_ptcl_per_side
              y = get_pos_cell_center(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, j)
              ! z-loop
              do k=1,n_ptcl_per_side
                  z = get_pos_cell_center(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, k)
                  r = dsqrt(x*x + y*y + z*z)
                  if (r <= radius_of_sphere) n_sph = n_sph + 1
              end do
          end do
      end do
      write(*,*)'n_sph = ',n_sph
      if (n_sph == 0) then
         call fdps_ctrl%PS_abort()
         stop 1
      end if
      call fdps_ctrl%set_nptcl_loc(psys_num_sph,n_sph)
      call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
      ! (2) Actually place particles
      id = 1
      ! x-loop
      do i=1,n_ptcl_per_side
          x = get_pos_cell_center(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, i)
          ! y-loop
          do j=1,n_ptcl_per_side
              y = get_pos_cell_center(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, j)
              ! z-loop
              do k=1,n_ptcl_per_side
                  z = get_pos_cell_center(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, k)
                  r = dsqrt(x*x + y*y + z*z)
                  if (r <= radius_of_sphere) then
                      r_new = radius_of_sphere * (r/radius_of_sphere)**(1.5d0)
                      ptcl_sph(id)%id    = id + N_nbody
                      ptcl_sph(id)%mass  = M_sph / N_sph
                      ptcl_sph(id)%pos%x = (r_new / r) * x
                      ptcl_sph(id)%pos%y = (r_new / r) * y
                      ptcl_sph(id)%pos%z = (r_new / r) * z
                      ptcl_sph(id)%vel = 0.0d0
                      ptcl_sph(id)%acc_grav  = 0.0d0
                      ptcl_sph(id)%pot_grav  = 0.0d0
                      ptcl_sph(id)%acc_hydro = 0.0d0
                      ptcl_sph(id)%eng = 0.05d0 * M_sph / radius_of_sphere
                      dens = M_sph / (2.0d0 * pi * radius_of_sphere * radius_of_sphere * r_new)
                      p = specific_heat_ratio - 1.0d0
                      ptcl_sph(id)%ent = p * ptcl_sph(id)%eng / dens**(p)
                      ptcl_sph(id)%smth = radius_of_sphere * (dble(N_neighbor)/dble(N_sph))**(1.0d0/3.0d0)
                      id = id + 1
                  end if
              end do
          end do
      end do
   else
      ! In this mode, we set an initial distribution of particles by reading a file.
      ! (1) Read particle data
      file_name = "result/glass_data_header.dat"
      open(unit=9,file=trim(file_name),action='read', &
           form='unformatted',access='stream',status='old')
         read(9)nptcl_in
      close(unit=9)
      allocate( pos(nptcl_in) )
      file_name = "result/glass_data_body.dat"
      open(unit=9,file=trim(file_name),action='read', &
           form='unformatted',access='stream',status='old')
         do i=1,nptcl_in
            read(9)pos(i)%x,pos(i)%y,pos(i)%z
         end do
      close(unit=9)
      write(*,*)nptcl_in," particles are read."
      ! For debug
      file_name = 'glass_data.txt'
      open(unit=9,file=trim(file_name),action='write',status='replace')
         do i=1,nptcl_in
            write(9,'(3es25.16e3)')pos(i)%x,pos(i)%y,pos(i)%z
         end do
      close(unit=9)
      ! (2) Count # of particles in a sphere of radius 1
      n_sph = 0
      do i=1,nptcl_in
         r2 = pos(i) * pos(i)
         if (r2 < 1.0d0) n_sph = n_sph + 1
      end do
      write(*,*)n_sph," particles of them will be used to make a Evrard sphere."
      ! (3) Place SPH particles
      call fdps_ctrl%set_nptcl_loc(psys_num_sph,n_sph)
      call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
      j = 0
      do i=1,n_sph
         ptcl_sph(i)%id = i 
         ptcl_sph(i)%mass = M_sph / n_sph
         do
            j = j + 1
            r2 = pos(j) * pos(j)
            if (r2 < 1.0d0) exit
         end do
         r = dsqrt(r2)
         if (r > 0.0d0) then
             r_new = radius_of_sphere * r**(1.5d0)
             ptcl_sph(i)%pos = (r_new / r) * pos(j)
         else
             ptcl_sph(i)%pos = pos(j)
         end if
         ptcl_sph(i)%vel       = 0.0d0
         ptcl_sph(i)%acc_grav  = 0.0d0
         ptcl_sph(i)%pot_grav  = 0.0d0
         ptcl_sph(i)%acc_hydro = 0.0d0
         ptcl_sph(i)%eng = 0.05d0 * M_sph / radius_of_sphere
         ptcl_sph(i)%smth = radius_of_sphere * (dble(N_neighbor)/dble(N_sph))**(1.0d0/3.0d0)
         ! Note that entropy is determined in setEntropy().
      end do
   end if
   ! Set boundary condition
   bc = fdps_bc_open
   ! Set the other parameters
   eps_grav = 1.0e-4 * radius_of_sphere
   if (fdps_ctrl%get_rank() == 0) then
       write(*,*)"An initial condition for the Evrard test is made."
       write(*,*)"(N_sph = ",N_sph,")"
   end if
   ! Set I/O intervals
   unit_dens = 3.0d0 * M_sph / (4.0d0 * pi * radius_of_sphere**(3.0d0))
   unit_time = dsqrt(pi * pi / 8.0d0) * radius_of_sphere**(1.5d0) / dsqrt(M_sph)
   write(*,*)"unit_dens = ",unit_dens
   write(*,*)"unit_time = ",unit_time
   dt_dump   = 0.1d0 * unit_time
   time_dump = dt_dump
   time_end  = unit_time
   ! Set the maximum timestep
   dt_max = 1.0d-3

end subroutine Evrard_test_IC

!-----------------------------------------------------------------------
!//////////////              F U N C T I O N              //////////////
!////////////// < G E T _ P O S _ C E L L _ C E N T E R > //////////////
!-----------------------------------------------------------------------
function get_pos_cell_center(pos_left_bound,  &
                             pos_right_bound, &
                             number_of_cells, &
                             i)
    implicit none
    double precision :: get_pos_cell_center
    double precision, intent(in) :: pos_left_bound,pos_right_bound
    integer, intent(in) :: number_of_cells,i
    !* Local variables
    double precision :: dx,center

    !* Error check
    if ((i <= 0) .or. (number_of_cells < i) .or. &
        (pos_left_bound > pos_right_bound)) then
       stop 1
    end if

    dx = (pos_right_bound - pos_left_bound) / number_of_cells
    if ( mod(number_of_cells,2) == 0) then
        ! # of cells is even.
        if (i <= number_of_cells/2) then
            get_pos_cell_center = pos_left_bound + dx * (i - 0.5d0)
        else 
            get_pos_cell_center = pos_right_bound - dx * (number_of_cells - i + 0.5d0)
        end if
    else
        ! # of cells is odd.
        center = 0.5d0 * (pos_left_bound + pos_right_bound)
        if (i <= number_of_cells/2) then
            get_pos_cell_center = center - dx * (number_of_cells/2 - i)
        else 
            get_pos_cell_center = center + dx * (i - number_of_cells/2)
        end if
    end if
end function get_pos_cell_center

!-----------------------------------------------------------------------
!////////////////////      S U B R O U T I N E      ////////////////////
!//////////////////// < M A K E _ G L A S S _ I C > ////////////////////
!-----------------------------------------------------------------------
subroutine make_glass_IC(psys_num_nbody,psys_num_sph, &
                         bc,                          &
                         pos_root_domain_low,         &
                         pos_root_domain_high,        &
                         time_dump,dt_dump,time_end)
   use mathematical_constants
   use physical_constants
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   integer(kind=c_int), intent(inout) :: bc
   type(fdps_f64vec), intent(inout) :: pos_root_domain_low, &
                                       pos_root_domain_high
   double precision, intent(inout) :: time_dump,dt_dump,time_end
   !* Local parameters
   integer, parameter :: N_nbody = 1 ! dummy value
   integer, parameter :: N_sph   = 2**18
   !* Local variables
   integer :: i,nptcl_loc
   integer :: mtts_num
   double precision :: x,y,z
   double precision :: tcross
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph

   ! Initialize pseudorandom number generator
   call fdps_ctrl%create_mtts(mtts_num)
   call fdps_ctrl%mtts_init_genrand(mtts_num,0)
   ! Place Nbody particles
   call fdps_ctrl%set_nptcl_loc(psys_num_nbody,N_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   do i=1,N_nbody
       ptcl_nbody(i)%id   = i
       ptcl_nbody(i)%mass = 0.0d0
       ptcl_nbody(i)%pos  = 0.0d0
       ptcl_nbody(i)%vel  = 0.0d0
       ptcl_nbody(i)%acc  = 0.0d0
       ptcl_nbody(i)%pot  = 0.0d0
   end do
   ! Place SPH particles
   call fdps_ctrl%set_nptcl_loc(psys_num_sph,N_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   do i=1,N_sph
       ptcl_sph(i)%id = i + N_nbody
       ptcl_sph(i)%mass = 8.0d0 / N_sph
       do
           x = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) 
           y = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) 
           z = (2.0d0 * fdps_ctrl%mtts_genrand_res53(mtts_num) - 1.0d0) 
           if ((-1.0d0 <= x) .and. (x < 1.0d0) .and. &
               (-1.0d0 <= y) .and. (y < 1.0d0) .and. &
               (-1.0d0 <= z) .and. (z < 1.0d0)) then 
               exit
           end if    
       end do
       ptcl_sph(i)%pos%x     = x
       ptcl_sph(i)%pos%y     = y
       ptcl_sph(i)%pos%z     = z
       ptcl_sph(i)%vel       = 0.0d0
       ptcl_sph(i)%acc_grav  = 0.0d0
       ptcl_sph(i)%pot_grav  = 0.0d0
       ptcl_sph(i)%acc_hydro = 0.0d0
       ptcl_sph(i)%eng       = 1.0d0
       ptcl_sph(i)%smth      = 2.0d0 * (dble(N_neighbor)/dble(N_sph))**(1.0d0/3.0d0) 
       ! [Notes]
       !   (1) The value of the specific thermal energy is chosen 
       !       so that the sound speed is nearly equal to 1.
       !   (2) The value of the entropy is determined in set_entropy().
   end do
   ! Set boundary condition
   bc = fdps_bc_periodic_xyz 
   pos_root_domain_low%x  = - 1.0d0
   pos_root_domain_low%y  = - 1.0d0
   pos_root_domain_low%z  = - 1.0d0
   pos_root_domain_high%x = 1.0d0 
   pos_root_domain_high%y = 1.0d0 
   pos_root_domain_high%z = 1.0d0 
   ! Set gravitational softening
   eps_grav = 0.01d0 
   ! Set I/O intervals
   tcross = 2.0d0 * dsqrt(3.0d0)
   write(*,*)"The sound crossing time = ",tcross
   dt_dump   = 4.0d0 * tcross
   time_dump = dt_dump
   time_end  = 64.0 * tcross
   ! Set maximum timestep
   dt_max = 1.0d-2

end subroutine make_glass_IC
