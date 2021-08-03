!-----------------------------------------------------------------------
!///////////////////// < M A I N   R O U T I N E > /////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   !* Local parameters
   !-(force parameters)
   real, parameter :: theta = 0.5
   integer, parameter :: n_leaf_limit = 8
   integer, parameter :: n_group_limit = 64
   !-(domain decomposition)
   real, parameter :: coef_ema=0.3
   !-(IO)
   integer, parameter :: output_interval=10
   !* Local variables
   integer :: i,j,k,ierr
   integer :: nstep
   integer :: psys_num,dinfo_num
   integer :: tree_num_dens,tree_num_hydro
   integer :: n_tot,n_loc
   logical :: clear
   double precision :: time,dt,end_time
   type(fdps_f64vec) :: pos_ll,pos_ul
   type(fdps_controller) :: fdps_ctrl
   type(full_particle), dimension(:), pointer :: ptcl
   type(c_funptr) :: pfunc_ep_ep
   !-(IO)
   character(len=64) :: filename
   !* External routines
   double precision, external :: get_timestep

   !* Initialize FDPS
   call fdps_ctrl%PS_Initialize()

   !* Make an instance of ParticleSystem and initialize it
   call fdps_ctrl%create_psys(psys_num,'full_particle')
   call fdps_ctrl%init_psys(psys_num)

   !* Make an initial condition and initialize the particle system
   call setup_IC(fdps_ctrl,psys_num,end_time,pos_ll,pos_ul)

   !* Make an instance of DomainInfo and initialize it
   call fdps_ctrl%create_dinfo(dinfo_num)
   call fdps_ctrl%init_dinfo(dinfo_num,coef_ema)
   call fdps_ctrl%set_boundary_condition(dinfo_num,fdps_bc_periodic_xyz)
   call fdps_ctrl%set_pos_root_domain(dinfo_num,pos_ll,pos_ul)

   !* Perform domain decomposition and exchange particles
   call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
   call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

   !* Make two tree structures
   n_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   !** dens_tree (used for the density calculation)
   call fdps_ctrl%create_tree(tree_num_dens, &
                              "Short,force_dens,essential_particle,essential_particle,Gather")
   call fdps_ctrl%init_tree(tree_num_dens,n_loc,theta, &
                            n_leaf_limit,n_group_limit)

   !** hydro_tree (used for the force calculation)
   call fdps_ctrl%create_tree(tree_num_hydro, &
                              "Short,force_hydro,essential_particle,essential_particle,Symmetry")
   call fdps_ctrl%init_tree(tree_num_hydro,n_loc,theta, &
                            n_leaf_limit,n_group_limit)

   !* Compute density, pressure, acceleration due to pressure gradient
   pfunc_ep_ep = c_funloc(calc_density)
   call fdps_ctrl%calc_force_all_and_write_back(tree_num_dens, &
                                                pfunc_ep_ep,   &
                                                psys_num,      &
                                                dinfo_num)
   call set_pressure(fdps_ctrl,psys_num)
   pfunc_ep_ep = c_funloc(calc_hydro_force)
   call fdps_ctrl%calc_force_all_and_write_back(tree_num_hydro, &
                                                pfunc_ep_ep,   &
                                                psys_num,      &
                                                dinfo_num)
   !* Get timestep
   dt = get_timestep(fdps_ctrl,psys_num)

   !* Main loop for time integration
   nstep = 0; time = 0.0d0
   do 
      !* Leap frog: Initial Kick & Full Drift
      call initial_kick(fdps_ctrl,psys_num,dt)
      call full_drift(fdps_ctrl,psys_num,dt)

      !* Adjust the positions of the SPH particles that run over
      !  the computational boundaries.
      call fdps_ctrl%adjust_pos_into_root_domain(psys_num,dinfo_num)

      !* Leap frog: Predict
      call predict(fdps_ctrl,psys_num,dt)

      !* Perform domain decomposition and exchange particles again
      call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
      call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

      !* Compute density, pressure, acceleration due to pressure gradient
      pfunc_ep_ep = c_funloc(calc_density)
      call fdps_ctrl%calc_force_all_and_write_back(tree_num_dens, &
                                                   pfunc_ep_ep,   &
                                                   psys_num,      &
                                                   dinfo_num)
      call set_pressure(fdps_ctrl,psys_num)
      pfunc_ep_ep = c_funloc(calc_hydro_force)
      call fdps_ctrl%calc_force_all_and_write_back(tree_num_hydro, &
                                                   pfunc_ep_ep,    &
                                                   psys_num,       &
                                                   dinfo_num)

      !* Get a new timestep
      dt = get_timestep(fdps_ctrl,psys_num)

      !* Leap frog: Final Kick
      call final_kick(fdps_ctrl,psys_num,dt)

      !* Output result files
      if (mod(nstep,output_interval) == 0) then
         call output(fdps_ctrl,psys_num,nstep)
         call check_cnsrvd_vars(fdps_ctrl,psys_num)
      end if

      !* Output information to STDOUT
      if (fdps_ctrl%get_rank() == 0) then
         write(*,200)time,nstep
         200 format("================================"/ &
                    "time  = ",1es25.16e3/              &
                    "nstep = ",i6/                      &
                    "================================")
      end if

      !* Termination condition
      if (time >= end_time) exit

      !* Update time & step
      time  = time + dt
      nstep = nstep + 1

   end do
   call fdps_ctrl%ps_finalize()
   stop 0

   !* Finalize FDPS
   call fdps_ctrl%PS_Finalize()

end subroutine f_main

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!///////////////////////// < S E T U P _ I C > /////////////////////////
!-----------------------------------------------------------------------
subroutine setup_IC(fdps_ctrl,psys_num,end_time,pos_ll,pos_ul)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   double precision, intent(inout) :: end_time
   type(fdps_f64vec) :: pos_ll,pos_ul
   !* Local variables
   integer :: i,ii,irank
   integer :: nprocs, myrank
   integer :: nx_left, ny_left, nz_left
   integer :: nx_right, ny_right, nz_right
   integer :: nptcl_loc, nptcl_glb
   integer :: nptcl_quot, nptcl_rem
   integer :: i_first, i_last
   double precision :: dens_left, dens_right
   double precision :: eng_left, eng_right
   double precision :: sv, mass
   double precision :: x,y,z,dx,dy,dz
   type(full_particle), dimension(:), pointer :: ptcl
   character(len=5)  :: proc_num
   character(len=64) :: fname

   !* Get # of MPI processes and rank number
   nprocs = fdps_ctrl%get_num_procs()
   myrank = fdps_ctrl%get_rank()

   !* Set the box size
   pos_ll%x = 0.0d0
   pos_ll%y = 0.0d0
   pos_ll%z = 0.0d0
   pos_ul%x = 1.0d0
   pos_ul%y = pos_ul%x / 8.0d0
   pos_ul%z = pos_ul%x / 8.0d0

   !* Set the left and right states
   dens_left  = 1.0d0
   eng_left   = 2.5d0
   dens_right = 0.5d0
   eng_right  = 2.5d0

   !* Set the separation of particle of the left state
   dx = 1.0d0 / 128.0d0
   dy = dx
   dz = dx

   !* Count # of particles in the computational domain
   !** (1) Left-half
   nx_left = 0
   x = 0.0d0
   do
      nx_left = nx_left + 1
      x = x + dx
      if (x >= 0.5d0*pos_ul%x) exit
   end do 
   ny_left = 0
   y = 0.0d0
   do 
      ny_left = ny_left + 1
      y = y + dy
      if (y >= pos_ul%y) exit
   end do
   nz_left = 0
   z = 0.0d0
   do
      nz_left = nz_left + 1
      z = z + dz
      if (z >= pos_ul%z) exit
   end do
   !** (2) Right-half
   nx_right = 0
   x = 0.5d0*pos_ul%x
   do
      nx_right = nx_right + 1
      x = x + (dens_left/dens_right)*dx
      if (x >= pos_ul%x) exit
   end do
   ny_right = 0
   y = 0.0d0
   do 
      ny_right = ny_right + 1
      y = y + dy
      if (y >= pos_ul%y) exit
   end do
   nz_right = 0
   z = 0.0d0
   do 
      nz_right = nz_right + 1
      z = z + dz
      if (z >= pos_ul%z) exit
   end do
   !** (3) calculate # of particles
   nptcl_glb = (nx_left * ny_left * nz_left) &
             + (nx_right * ny_right * nz_right)
   if (myrank == 0) then
      write(*,*)'nptcl_glb = ',nptcl_glb
   end if

   !* Set # of local particles
   nptcl_quot = nptcl_glb / nprocs
   nptcl_rem  = mod(nptcl_glb, nprocs)
   if (myrank == 0) then
      write(*,*)'nptcl_quot = ',nptcl_quot
      write(*,*)'nptcl_rem  = ',nptcl_rem
   end if
   nptcl_loc = nptcl_quot
   if (myrank < nptcl_rem) nptcl_loc = nptcl_loc + 1
   call fdps_ctrl%set_nptcl_loc(psys_num, nptcl_loc)
   i_first = 1 + nptcl_quot * myrank
   if (myrank > (nptcl_rem-1)) then
      i_first = i_first + nptcl_rem
   else
      i_first = i_first + myrank
   end if
   i_last = i_first + nptcl_loc - 1
   do irank=0, nprocs-1
      if (irank == myrank) then
         write(*,"(4(a,1x,i8))")'myrank =',myrank, &
                   ', nptcl_loc =',nptcl_loc, &
                   ', i_first =',i_first, &
                   ', i_last =',i_last
      end if
      call fdps_ctrl%barrier()
   end do

   !* Set local particles
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   sv = (pos_ul%x * pos_ul%y * pos_ul%z) / nptcl_glb ! specific volume
   mass = 0.5d0*(dens_left+dens_right) * sv
   i = 1
   !** (1) Left-half
   x = 0.0d0
   do 
      y = 0.0d0
      do 
         z = 0.0d0
         do 
            if ((i_first <= i) .and. (i <= i_last)) then
               ii = i - i_first + 1
               ptcl(ii)%mass  = mass
               ptcl(ii)%pos%x = x
               ptcl(ii)%pos%y = y
               ptcl(ii)%pos%z = z
               ptcl(ii)%dens  = dens_left
               ptcl(ii)%eng   = eng_left
               ptcl(ii)%smth  = kernel_support_radius * 0.012d0
               ptcl(ii)%id    = i
            end if
            i = i + 1
            z = z + dz
            if (z >= pos_ul%z) exit
         end do
         y = y + dy
         if (y >= pos_ul%y) exit
      end do
      x = x + dx
      if (x >= 0.5d0*pos_ul%x) exit
   end do
   !** (2) Right-half
   x = 0.5d0*pos_ul%x
   do
      y = 0.0d0
      do 
         z = 0.0d0
         do 
            if ((i_first <= i) .and. (i <= i_last)) then
               ii = i - i_first + 1
               ptcl(ii)%mass  = mass
               ptcl(ii)%pos%x = x
               ptcl(ii)%pos%y = y
               ptcl(ii)%pos%z = z
               ptcl(ii)%dens  = dens_right
               ptcl(ii)%eng   = eng_right
               ptcl(ii)%smth  = kernel_support_radius * 0.012d0
               ptcl(ii)%id    = i
            end if
            i = i + 1
            z = z + dz
            if (z >= pos_ul%z) exit
         end do
         y = y + dy
         if (y >= pos_ul%y) exit
      end do
      x = x + (dens_left/dens_right)*dx
      if (x >= pos_ul%x) exit
   end do

   !* Check the initial distribution
   !write(proc_num,"(i5.5)")myrank
   !fname = "initial" // proc_num // ".dat"
   !open(unit=9,file=trim(fname),action='write',status='replace')
   !   do i=1,nptcl_loc
   !      write(9,'(3es25.16e3)')ptcl(i)%pos%x, &
   !                             ptcl(i)%pos%y, &
   !                             ptcl(i)%pos%z
   !   end do
   !close(unit=9)

   !* Set the end time
   end_time = 0.12d0

   !* Inform to STDOUT
   if (fdps_ctrl%get_rank() == 0) then
      write(*,*)"setup..."
   end if
  !call fdps_ctrl%ps_finalize()
  !stop 0

end subroutine setup_IC

!-----------------------------------------------------------------------
!/////////////////////     S U B R O U T I N E     /////////////////////
!///////////////////// < G E T _ T I M E S T E P > /////////////////////
!-----------------------------------------------------------------------
function get_timestep(fdps_ctrl,psys_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   real(kind=c_double) :: get_timestep
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer, intent(in) :: psys_num
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl
   real(kind=c_double) :: dt_loc

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   dt_loc = 1.0d30
   do i=1,nptcl_loc
      dt_loc = min(dt_loc, ptcl(i)%dt)
   end do
   nullify(ptcl)

   !* Reduction
   call fdps_ctrl%get_min_value(dt_loc,get_timestep)

end function get_timestep

!-----------------------------------------------------------------------
!/////////////////////     S U B R O U T I N E     /////////////////////
!///////////////////// < I N I T I A L _ K I C K > /////////////////////
!-----------------------------------------------------------------------
subroutine initial_kick(fdps_ctrl,psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer, intent(in) :: psys_num
   double precision, intent(in) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%vel_half = ptcl(i)%vel + 0.5d0 * dt * ptcl(i)%acc
      ptcl(i)%eng_half = ptcl(i)%eng + 0.5d0 * dt * ptcl(i)%eng_dot
   end do
   nullify(ptcl)

end subroutine initial_kick

!-----------------------------------------------------------------------
!///////////////////////   S U B R O U T I N E   ///////////////////////
!/////////////////////// < F U L L _ D R I F T > ///////////////////////
!-----------------------------------------------------------------------
subroutine full_drift(fdps_ctrl,psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer, intent(in) :: psys_num
   double precision, intent(in) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%pos = ptcl(i)%pos + dt * ptcl(i)%vel_half
   end do
   nullify(ptcl)

end subroutine full_drift

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////  < P R E D I C T >  /////////////////////////
!-----------------------------------------------------------------------
subroutine predict(fdps_ctrl,psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer, intent(in) :: psys_num
   double precision, intent(in) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%vel = ptcl(i)%vel + dt * ptcl(i)%acc
      ptcl(i)%eng = ptcl(i)%eng + dt * ptcl(i)%eng_dot
   end do
   nullify(ptcl)

end subroutine predict

!-----------------------------------------------------------------------
!///////////////////////   S U B R O U T I N E   ///////////////////////
!/////////////////////// < F I N A L _ K I C K > ///////////////////////
!-----------------------------------------------------------------------
subroutine final_kick(fdps_ctrl,psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer, intent(in) :: psys_num
   double precision, intent(in) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%vel = ptcl(i)%vel_half + 0.5d0 * dt * ptcl(i)%acc
      ptcl(i)%eng = ptcl(i)%eng_half + 0.5d0 * dt * ptcl(i)%eng_dot
   end do
   nullify(ptcl)

end subroutine final_kick

!-----------------------------------------------------------------------
!/////////////////////     S U B R O U T I N E     /////////////////////
!///////////////////// < S E T _ P R E S S U R E > /////////////////////
!-----------------------------------------------------------------------
subroutine set_pressure(fdps_ctrl,psys_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer, intent(in) :: psys_num
   !* Local parameters
   double precision, parameter :: hcr=1.4d0
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%pres = (hcr - 1.0d0) * ptcl(i)%dens * ptcl(i)%eng
      ptcl(i)%snds = dsqrt(hcr * ptcl(i)%pres / ptcl(i)%dens)
   end do
   nullify(ptcl)

end subroutine set_pressure

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////   < O U T P U T >   /////////////////////////
!-----------------------------------------------------------------------
subroutine output(fdps_ctrl,psys_num,nstep)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   integer, intent(IN) :: nstep
   !* Local parameters
   character(len=16), parameter :: root_dir="result"
   character(len=16), parameter :: file_prefix_1st="snap"
   character(len=16), parameter :: file_prefix_2nd="proc"
   !* Local variables
   integer :: i,nptcl_loc
   integer :: myrank
   character(len=5) :: file_num,proc_num
   character(len=64) :: cmd,sub_dir,fname
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get the rank number
   myrank = fdps_ctrl%get_rank()

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data 
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)

   !* Output
   write(file_num,"(i5.5)")nstep
   write(proc_num,"(i5.5)")myrank
   fname =  trim(root_dir) // "/" &
         // trim(file_prefix_1st) // file_num // "-" &
         // trim(file_prefix_2nd) // proc_num // ".dat"
   open(unit=9,file=trim(fname),action='write',status='replace')
      do i=1,nptcl_loc
         write(9,100)ptcl(i)%id,ptcl(i)%mass, &
                     ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z, &
                     ptcl(i)%vel%x,ptcl(i)%vel%y,ptcl(i)%vel%z, &
                     ptcl(i)%dens,ptcl(i)%eng,ptcl(i)%pres
         100 format(i8,1x,10e25.16e3)
      end do
   close(unit=9)
   nullify(ptcl)

end subroutine output

!-----------------------------------------------------------------------
!////////////////          S U B R O U T I N E          ////////////////
!//////////////// < C H E C K _ C N S R V D _ V A R S > ////////////////
!-----------------------------------------------------------------------
subroutine check_cnsrvd_vars(fdps_ctrl,psys_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer, intent(in) :: psys_num
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl
   type(fdps_f64vec) :: mom_loc,mom
   real(kind=c_double) :: eng_loc,eng

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   mom_loc = 0.0d0; eng_loc = 0.0d0
   do i=1,nptcl_loc
      mom_loc = mom_loc + ptcl(i)%vel * ptcl(i)%mass
      eng_loc = eng_loc + ptcl(i)%mass &
                         *(ptcl(i)%eng &
                          +0.5d0*ptcl(i)%vel*ptcl(i)%vel)
   end do
   nullify(ptcl)

   !* Reduction & output
   call fdps_ctrl%get_sum(eng_loc,eng)
   call fdps_ctrl%get_sum(mom_loc%x,mom%x)
   call fdps_ctrl%get_sum(mom_loc%y,mom%y)
   call fdps_ctrl%get_sum(mom_loc%z,mom%z)
   if (fdps_ctrl%get_rank() == 0) then
      write(*,100)eng
      write(*,100)mom%x
      write(*,100)mom%y
      write(*,100)mom%z
      100 format(1es25.16e3)
   end if

end subroutine check_cnsrvd_vars
