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
   integer :: dens_tree_num,hydro_tree_num
   integer :: ntot,nloc
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
   ntot = fdps_ctrl%get_nptcl_glb(psys_num)
   !** dens_tree (used for the density calculation)
   call fdps_ctrl%create_tree(dens_tree_num, &
                              "Short,dens_force,essential_particle,essential_particle,Gather")
   call fdps_ctrl%init_tree(dens_tree_num,3*ntot,theta, &
                            n_leaf_limit,n_group_limit)

   !** hydro_tree (used for the force calculation)
   call fdps_ctrl%create_tree(hydro_tree_num, &
                              "Short,hydro_force,essential_particle,essential_particle,Symmetry")
   call fdps_ctrl%init_tree(hydro_tree_num,3*ntot,theta, &
                            n_leaf_limit,n_group_limit)

   !* Compute density, pressure, acceleration due to pressure gradient
   pfunc_ep_ep = c_funloc(calc_density)
   call fdps_ctrl%calc_force_all_and_write_back(dens_tree_num, &
                                                pfunc_ep_ep,   &
                                                psys_num,      &
                                                dinfo_num)
   call set_pressure(fdps_ctrl,psys_num)
   pfunc_ep_ep = c_funloc(calc_hydro_force)
   call fdps_ctrl%calc_force_all_and_write_back(hydro_tree_num, &
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
      call fdps_ctrl%calc_force_all_and_write_back(dens_tree_num, &
                                                   pfunc_ep_ep,   &
                                                   psys_num,      &
                                                   dinfo_num)
      call set_pressure(fdps_ctrl,psys_num)
      pfunc_ep_ep = c_funloc(calc_hydro_force)
      call fdps_ctrl%calc_force_all_and_write_back(hydro_tree_num, &
                                                   pfunc_ep_ep,   &
                                                   psys_num,      &
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
   integer :: i
   integer :: nprocs,myrank
   integer :: nptcl_glb
   double precision :: dens_L,dens_R,eng_L,eng_R
   double precision :: x,y,z,dx,dy,dz
   double precision :: dx_tgt,dy_tgt,dz_tgt
   type(full_particle), dimension(:), pointer :: ptcl
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

   !* Make an initial condition at RANK 0
   if (myrank == 0) then
      !* Set the left and right states
      dens_L = 1.0d0
      eng_L  = 2.5d0
      dens_R = 0.5d0
      eng_R  = 2.5d0
      !* Set the separation of particle of the left state
      dx = 1.0d0 / 128.0d0
      dy = dx
      dz = dx
      !* Set the number of local particles
      nptcl_glb = 0
      !** (1) Left-half
      x = 0.0d0
      do 
         y = 0.0d0
         do 
            z = 0.0d0
            do 
               nptcl_glb = nptcl_glb + 1
               z = z + dz
               if (z >= pos_ul%z) exit
            end do
            y = y + dy
            if (y >= pos_ul%y) exit
         end do
         x = x + dx
         if (x >= 0.5d0*pos_ul%x) exit
      end do
      write(*,*)'nptcl_glb(L)   = ',nptcl_glb
      !** (2) Right-half
      x = 0.5d0*pos_ul%x
      do
         y = 0.0d0
         do 
            z = 0.0d0
            do 
               nptcl_glb = nptcl_glb + 1
               z = z + dz
               if (z >= pos_ul%z) exit
            end do
            y = y + dy
            if (y >= pos_ul%y) exit
         end do
         x = x + (dens_L/dens_R)*dx
         if (x >= pos_ul%x) exit
      end do
      write(*,*)'nptcl_glb(L+R) = ',nptcl_glb
      !* Place SPH particles
      call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_glb)
      call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
      i = 0
      !** (1) Left-half
      x = 0.0d0
      do 
         y = 0.0d0
         do 
            z = 0.0d0
            do 
               i = i + 1
               ptcl(i)%id    = i
               ptcl(i)%pos%x = x
               ptcl(i)%pos%y = y
               ptcl(i)%pos%z = z
               ptcl(i)%dens  = dens_L
               ptcl(i)%eng   = eng_L
               z = z + dz
               if (z >= pos_ul%z) exit
            end do
            y = y + dy
            if (y >= pos_ul%y) exit
         end do
         x = x + dx
         if (x >= 0.5d0*pos_ul%x) exit
      end do
      write(*,*)'nptcl(L)   = ',i
      !** (2) Right-half
      x = 0.5d0*pos_ul%x
      do
         y = 0.0d0
         do 
            z = 0.0d0
            do 
               i = i + 1
               ptcl(i)%id    = i
               ptcl(i)%pos%x = x
               ptcl(i)%pos%y = y
               ptcl(i)%pos%z = z
               ptcl(i)%dens  = dens_R
               ptcl(i)%eng   = eng_R
               z = z + dz
               if (z >= pos_ul%z) exit
            end do
            y = y + dy
            if (y >= pos_ul%y) exit
         end do
         x = x + (dens_L/dens_R)*dx
         if (x >= pos_ul%x) exit
      end do
      write(*,*)'nptcl(L+R) = ',i
      !* Set particle mass and smoothing length
      do i=1,nptcl_glb
         ptcl(i)%mass = 0.5d0*(dens_L+dens_R)        &
                      * (pos_ul%x*pos_ul%y*pos_ul%z) &
                      / nptcl_glb
         ptcl(i)%smth = kernel_support_radius * 0.012d0
      end do

      !* Check the initial distribution
     !fname = "initial.dat"
     !open(unit=9,file=trim(fname),action='write',status='replace')
     !   do i=1,nptcl_glb
     !      write(9,'(3es25.16e3)')ptcl(i)%pos%x, &
     !                             ptcl(i)%pos%y, &
     !                             ptcl(i)%pos%z
     !   end do
     !close(unit=9)

   else
      call fdps_ctrl%set_nptcl_loc(psys_num,0)
   end if

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
