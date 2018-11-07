!-----------------------------------------------------------------------
!/////////////////////// < M A I N  R O U T I N E > ////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use fdps_module
#if defined(ENABLE_PHANTOM_GRAPE_X86)
   use phantom_grape_g5_x86
#endif
   use user_defined_types
   implicit none
   !* Local parameters
   integer, parameter :: ntot=2**10
   !-(force parameters)
   real, parameter :: theta = 0.5
   integer, parameter :: n_leaf_limit = 8
   integer, parameter :: n_group_limit = 64
   !-(domain decomposition)
   real, parameter :: coef_ema=0.3
   !-(timing parameters)
   double precision, parameter :: time_end = 10.0d0
   double precision, parameter :: dt = 1.0d0/128.0d0
   double precision, parameter :: dt_diag = 1.0d0/8.0d0
   double precision, parameter :: dt_snap = 1.0d0
   !* Local variables
   integer :: i,j,k,num_loop,ierr
   integer :: psys_num,dinfo_num,tree_num
   integer :: nloc
   logical :: clear
   double precision :: ekin0,epot0,etot0
   double precision :: ekin1,epot1,etot1
   double precision :: time_diag,time_snap,time_sys
   double precision :: r,acc
   type(fdps_controller) :: fdps_ctrl
   type(full_particle), dimension(:), pointer :: ptcl
   type(c_funptr) :: pfunc_ep_ep,pfunc_ep_sp
   !-(IO)
   character(len=64) :: fname
   integer(c_int) :: np

   !* Initialize FDPS
   call fdps_ctrl%PS_Initialize()

   !* Create domain info object
   call fdps_ctrl%create_dinfo(dinfo_num)
   call fdps_ctrl%init_dinfo(dinfo_num,coef_ema)

   !* Create particle system object
   call fdps_ctrl%create_psys(psys_num,'full_particle')
   call fdps_ctrl%init_psys(psys_num)

   !* Create tree object
   call fdps_ctrl%create_tree(tree_num, &
                              "Long,full_particle,full_particle,full_particle,Monopole")
   call fdps_ctrl%init_tree(tree_num,ntot,theta, &
                            n_leaf_limit,n_group_limit)

   !* Make an initial condition
   call setup_IC(fdps_ctrl,psys_num,ntot)

   !* Domain decomposition and exchange particle
   call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
   call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

#if defined(ENABLE_PHANTOM_GRAPE_X86)
    call g5_open()
    call g5_set_eps_to_all(eps_grav);
#endif

   !* Compute force at the initial time
   pfunc_ep_ep = c_funloc(calc_gravity_ep_ep)
   pfunc_ep_sp = c_funloc(calc_gravity_ep_sp)
   call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                pfunc_ep_ep, &
                                                pfunc_ep_sp, &
                                                psys_num,    &
                                                dinfo_num)
   !* Compute energies at the initial time
   clear = .true.
   call calc_energy(fdps_ctrl,psys_num,etot0,ekin0,epot0,clear)

   !* Time integration
   time_diag = 0.0d0
   time_snap = 0.0d0
   time_sys  = 0.0d0
   num_loop = 0
   do 
      !* Output
     !if (fdps_ctrl%get_rank() == 0) then
     !   write(*,50)num_loop,time_sys
     !   50 format('(num_loop, time_sys) = ',i5,1x,1es25.16e3)
     !end if
      if ( (time_sys >= time_snap) .or. &
           (((time_sys + dt) - time_snap) > (time_snap - time_sys)) ) then
         call output(fdps_ctrl,psys_num)
         time_snap = time_snap + dt_snap
      end if
      
      !* Compute energies and output the results
      clear = .true.
      call calc_energy(fdps_ctrl,psys_num,etot1,ekin1,epot1,clear)
      if (fdps_ctrl%get_rank() == 0) then
         if ( (time_sys >= time_diag) .or. &
              (((time_sys + dt) - time_diag) > (time_diag - time_sys)) ) then
            write(*,100)time_sys,(etot1-etot0)/etot0
            100 format("time: ",1es20.10e3,", energy error: ",1es20.10e3)
            time_diag = time_diag + dt_diag
         end if
      end if

      !* Leapfrog: Kick-Drift
      call kick(fdps_ctrl,psys_num,0.5d0*dt)
      time_sys = time_sys + dt
      call drift(fdps_ctrl,psys_num,dt)

      !* Domain decomposition & exchange particle
      if (mod(num_loop,4) == 0) then
         call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
      end if
      call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

      !* Force calculation
      pfunc_ep_ep = c_funloc(calc_gravity_ep_ep)
      pfunc_ep_sp = c_funloc(calc_gravity_ep_sp)
      call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                   pfunc_ep_ep, &
                                                   pfunc_ep_sp, &
                                                   psys_num,    &
                                                   dinfo_num)
      !* Leapfrog: Kick
      call kick(fdps_ctrl,psys_num,0.5d0*dt)

      !* Update num_loop
      num_loop = num_loop + 1

      !* Termination
      if (time_sys >= time_end) then
         exit
      end if
   end do

#if defined(ENABLE_PHANTOM_GRAPE_X86)
   call g5_close()
#endif

   !* Finalize FDPS
   call fdps_ctrl%PS_Finalize()

end subroutine f_main

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!///////////////////////// < S E T U P _ I C > /////////////////////////
!-----------------------------------------------------------------------
subroutine setup_IC(fdps_ctrl,psys_num,nptcl_glb)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num,nptcl_glb
   !* Local parameters
   double precision, parameter :: m_tot=1.0d0
   double precision, parameter :: rmax=3.0d0,r2max=rmax*rmax
   !* Local variables
   integer :: i,j,k,ierr
   integer :: nprocs,myrank
   double precision :: r2,cm_mass
   type(fdps_f64vec) :: cm_pos,cm_vel,pos
   type(full_particle), dimension(:), pointer :: ptcl
   character(len=64) :: fname

   !* Get # of MPI processes and rank number
   nprocs = fdps_ctrl%get_num_procs()
   myrank = fdps_ctrl%get_rank()

   !* Make an initial condition at RANK 0
   if (myrank == 0) then
      !* Set # of local particles
      call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_glb)

      !* Create an uniform sphere of particles
      !** get the pointer to full particle data
      call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
      !** initialize Mersenne twister
      call fdps_ctrl%MT_init_genrand(0)
      do i=1,nptcl_glb
         ptcl(i)%id   = i 
         ptcl(i)%mass = m_tot/nptcl_glb
         do 
            ptcl(i)%pos%x = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * rmax
            ptcl(i)%pos%y = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * rmax
            ptcl(i)%pos%z = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * rmax
            r2 = ptcl(i)%pos*ptcl(i)%pos
            if ( r2 < r2max ) exit
         end do
         ptcl(i)%vel = 0.0d0
      end do

      !* Correction
      cm_pos  = 0.0d0
      cm_vel  = 0.0d0
      cm_mass = 0.0d0
      do i=1,nptcl_glb
         cm_pos  = cm_pos   + ptcl(i)%mass * ptcl(i)%pos
         cm_vel  = cm_vel   + ptcl(i)%mass * ptcl(i)%vel
         cm_mass = cm_mass  + ptcl(i)%mass
      end do
      cm_pos = cm_pos/cm_mass
      cm_vel = cm_vel/cm_mass
      do i=1,nptcl_glb
         ptcl(i)%pos = ptcl(i)%pos - cm_pos
         ptcl(i)%vel = ptcl(i)%vel - cm_vel
      end do
      
      !* Output
     !fname = 'initial.dat'
     !open(unit=9,file=trim(fname),action='write',status='replace', &
     !     form='unformatted',access='stream')
     !open(unit=9,file=trim(fname),action='write',status='replace')
     !   do i=1,nptcl_glb
     !     !write(9)ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z
     !      write(9,'(3es25.16e3)')ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z
     !   end do
     !close(unit=9)

      !* Release the pointer
      nullify( ptcl )

   else
      call fdps_ctrl%set_nptcl_loc(psys_num,0)
   end if

   !* Set the gravitational softening
   eps_grav = 1.0d0/32.0d0

end subroutine setup_IC

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////     < K I C K >     /////////////////////////
!-----------------------------------------------------------------------
subroutine kick(fdps_ctrl,psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   double precision, intent(IN) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%vel = ptcl(i)%vel + ptcl(i)%acc * dt
   end do
   nullify(ptcl)

end subroutine kick

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////    < D R I F T >    /////////////////////////
!-----------------------------------------------------------------------
subroutine drift(fdps_ctrl,psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   double precision, intent(IN) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%pos = ptcl(i)%pos + ptcl(i)%vel * dt
   end do
   nullify(ptcl)

end subroutine drift

!-----------------------------------------------------------------------
!//////////////////////    S U B R O U T I N E    //////////////////////
!////////////////////// < C A L C _ E N E R G Y > //////////////////////
!-----------------------------------------------------------------------
subroutine calc_energy(fdps_ctrl,psys_num,etot,ekin,epot,clear)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   double precision, intent(INOUT) :: etot,ekin,epot
   logical, intent(IN) :: clear
   !* Local variables
   integer :: i,nptcl_loc
   double precision :: etot_loc,ekin_loc,epot_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Clear energies
   if (clear .eqv. .true.) then
      etot = 0.0d0
      ekin = 0.0d0
      epot = 0.0d0
   end if

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)

   !* Compute energies
   ekin_loc = 0.0d0
   epot_loc = 0.0d0 
   do i=1,nptcl_loc
      ekin_loc = ekin_loc + ptcl(i)%mass * ptcl(i)%vel * ptcl(i)%vel
      epot_loc = epot_loc + ptcl(i)%mass * (ptcl(i)%pot + ptcl(i)%mass/eps_grav)
   end do
   ekin_loc = ekin_loc * 0.5d0
   epot_loc = epot_loc * 0.5d0
   etot_loc = ekin_loc + epot_loc
   call fdps_ctrl%get_sum(ekin_loc,ekin)
   call fdps_ctrl%get_sum(epot_loc,epot)
   call fdps_ctrl%get_sum(etot_loc,etot)

   !* Release the pointer
   nullify(ptcl)

end subroutine calc_energy

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////   < O U T P U T >   /////////////////////////
!-----------------------------------------------------------------------
subroutine output(fdps_ctrl,psys_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
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
   !* Static variables
   integer, save :: snap_num=0

   !* Get the rank number
   myrank = fdps_ctrl%get_rank()

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data 
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)

   !* Output
   write(file_num,"(i5.5)")snap_num
   write(proc_num,"(i5.5)")myrank
   fname =  trim(root_dir) // "/" &
         // trim(file_prefix_1st) // file_num // "-" &
         // trim(file_prefix_2nd) // proc_num // ".dat"
   open(unit=9,file=trim(fname),action='write',status='replace')
      do i=1,nptcl_loc
         write(9,100)ptcl(i)%id,ptcl(i)%mass, &
                     ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z, &
                     ptcl(i)%vel%x,ptcl(i)%vel%y,ptcl(i)%vel%z
         100 format(i8,1x,7e25.16e3)
      end do
   close(unit=9)
   nullify(ptcl)

   !* Update snap_num
   snap_num = snap_num + 1

end subroutine output

