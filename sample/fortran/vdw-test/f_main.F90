!-----------------------------------------------------------------------
!/////////////////////// < M A I N  R O U T I N E > ////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use fdps_module
   use user_defined_types
   implicit none
   !* Local parameters
   integer, parameter :: ntot = 512
   !-(force parameters)
   real, parameter :: theta = 0.5
   integer, parameter :: n_leaf_limit = 8
   integer, parameter :: n_group_limit = 64
   !-(timing parameters)
   double precision, parameter :: time_end = 10.0d0
   double precision, parameter :: dt = 1.0d0/128.0d0
   double precision, parameter :: dt_diag = 1.0d0/16.0d0
   double precision, parameter :: dt_snap = 1.0d0
   !-(domain parameters)
   real, parameter :: coef_ema = 0.3
   double precision, parameter :: boxdh=5.0d0
   !* Local variables
   integer :: i,j,k,num_loop
   integer :: psys_num,dinfo_num,tree_num
   integer :: nloc
   logical :: clear
   double precision :: ekin0,epot0,etot0
   double precision :: ekin1,epot1,etot1
   double precision :: time_diag,time_snap,time_sys
   double precision :: r,acc
   type(fdps_f64vec) :: pos_ll,pos_ul
   type(fdps_controller) :: fdps_ctrl
   type(fplj), dimension(:), pointer :: ptcl
   type(c_funptr) :: pfunc_ep_ep
   !-(IO)
   character(len=64) :: fname

   !* Initialize FDPS
   call fdps_ctrl%PS_Initialize()

   !* Create particle system object
   call fdps_ctrl%create_psys(psys_num,'fplj')
   call fdps_ctrl%init_psys(psys_num)

   !* Make an initial condition
   call setup_IC(fdps_ctrl,psys_num,ntot,boxdh)

   !* Create domain info object
   call fdps_ctrl%create_dinfo(dinfo_num)
   call fdps_ctrl%init_dinfo(dinfo_num,coef_ema)
   call fdps_ctrl%set_boundary_condition(dinfo_num,fdps_bc_periodic_xyz)
   pos_ll%x = -boxdh; pos_ll%y = -boxdh; pos_ll%z = -boxdh
   pos_ul%x =  boxdh; pos_ul%y =  boxdh; pos_ul%z =  boxdh
   call fdps_ctrl%set_pos_root_domain(dinfo_num,pos_ll,pos_ul)

   !* Domain decomposition and exchange particle
   call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
   call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

   !* Make a tree structure
   call fdps_ctrl%create_tree(tree_num, &
                              "Short,fplj,fplj,fplj,Scatter")
   call fdps_ctrl%init_tree(tree_num,3*ntot,theta, &
                            n_leaf_limit,n_group_limit)

   !* Compute force
   pfunc_ep_ep = c_funloc(calc_force_fpfp) 
   call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                pfunc_ep_ep, &
                                                psys_num,    &
                                                dinfo_num)

   !* Compute energies at the initial time
   clear = .true.
   call calc_energy(fdps_ctrl,psys_num,etot0,ekin0,epot0,clear)

   !* Initial half-kick
   call kick(fdps_ctrl,psys_num,0.5d0*dt);

   !* Time integration
   time_diag = 0.0d0
   time_snap = 0.0d0
   time_sys  = 0.0d0
   num_loop = 0
   do
      !* Output
      if ( (time_sys >= time_snap) .or. &
           (((time_sys + dt) - time_snap) > (time_snap - time_sys)) ) then
         call output(fdps_ctrl,psys_num)
         time_snap = time_snap + dt_snap
      end if

      !* Leapfrog: Full-drift
      call drift(fdps_ctrl,psys_num,dt,boxdh);

      !* Domain decomposition & exchange particle
      if (mod(num_loop,4) == 0) then
         call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
      end if
      call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

      !* Force calculation
      pfunc_ep_ep = c_funloc(calc_force_fpfp)
      call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                   pfunc_ep_ep, &
                                                   psys_num,    &
                                                   dinfo_num)

      !* Half-Kick 
      call kick(fdps_ctrl,psys_num,0.5d0*dt);
      time_sys = time_sys + dt

      !* Compute energies
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

      !* Leapfrog: Last-kick
      call kick(fdps_ctrl,psys_num,0.5d0*dt)

      !* Update num_loop
      num_loop = num_loop + 1

      !* Termination
      if (time_sys >= time_end) then
         exit
      end if
   end do

   !* Finalize FDPS
   call fdps_ctrl%PS_Finalize()

end subroutine f_main

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!///////////////////////// < S E T U P _ I C > /////////////////////////
!-----------------------------------------------------------------------
subroutine setup_IC(fdps_ctrl,psys_num,nptcl_glb,boxdh)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num,nptcl_glb
   double precision, intent(IN) :: boxdh
   !* Local parameters
   double precision, parameter :: vmax=0.1d0
   !* Local variables
   integer :: i,j,k,l,ierr
   integer :: nprocs,myrank
   integer :: nptcl_1d
   double precision :: dx,cm_mass
   type(fdps_f64vec) :: cm_pos,cm_vel,pos
   type(fplj), dimension(:), pointer :: ptcl
   character(len=64) :: fname
   !* External routines
   double precision, external :: cbrt

   !* Get # of MPI processes and rank number
   nprocs = fdps_ctrl%get_num_procs()
   myrank = fdps_ctrl%get_rank()

   !* Check if the number of particles is in the form of n^{3}
  !nptcl_1d = int((dble(nptcl_glb))**(1.0d0/3.0d0))
   nptcl_1d = int(cbrt(dble(nptcl_glb)))
   if (nptcl_1d*nptcl_1d*nptcl_1d /= nptcl_glb) then
      if (myrank == 0) then
         write(*,*)'nptcl_glb does not have a cubit root!'
         write(*,*)'nptcl_glb = ',nptcl_glb
         write(*,*)'nptcl_1d  = ',nptcl_1d
      end if
      call fdps_ctrl%PS_abort()
      stop 1
   else
      if (myrank == 0) write(*,*)'nptcl_1d = ',nptcl_1d
   end if
   

   !* Make an initial condition at RANK 0
   if (myrank == 0) then
      !* Set # of local particles
      call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_glb)

      !* Create an uniform grid of particles
      !** get the pointer to full particle data
      call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
      !** set the positions
      dx = (2.0d0*boxdh)/nptcl_1d
      l=0
      do k=1,nptcl_1d
      do j=1,nptcl_1d
      do i=1,nptcl_1d
         l = l + 1
         ptcl(l)%pos%x = - boxdh + dx * (i - 0.5)
         ptcl(l)%pos%y = - boxdh + dx * (j - 0.5)
         ptcl(l)%pos%z = - boxdh + dx * (k - 0.5)
      end do
      end do
      end do
      !** set velocities 
      call fdps_ctrl%MT_init_genrand(0) 
      do i=1,nptcl_glb
         ptcl(i)%id   = i 
         ! Charge
         ptcl(i)%mass = 1.0d0
         ! velocity
         ptcl(i)%vel%x = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * vmax
         ptcl(i)%vel%y = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * vmax
         ptcl(i)%vel%z = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * vmax
         !* search radius
         ptcl(i)%search_radius = 4.0d0
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
      fname = 'initial.dat'
     !open(unit=9,file=trim(fname),action='write',status='replace', &
     !     form='unformatted',access='stream')
      open(unit=9,file=trim(fname),action='write',status='replace')
         do i=1,nptcl_glb
           !write(9)ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z
            write(9,'(3es25.16e3)')ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z
         end do
      close(unit=9)

      !* Release the pointer
      nullify( ptcl )

   else
      call fdps_ctrl%set_nptcl_loc(psys_num,0)
   end if

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
   type(fplj), dimension(:), pointer :: ptcl

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
subroutine drift(fdps_ctrl,psys_num,dt,boxdh)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   double precision, intent(IN) :: dt,boxdh
   !* Local variables
   integer :: i,nptcl_loc
   type(fplj), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%pos = ptcl(i)%pos + ptcl(i)%vel * dt
      ! Adjust the position so that the particle is in the box
      if (ptcl(i)%pos%x < -boxdh) ptcl(i)%pos%x = ptcl(i)%pos%x + 2.0d0*boxdh
      if (ptcl(i)%pos%x >  boxdh) ptcl(i)%pos%x = ptcl(i)%pos%x - 2.0d0*boxdh
      if (ptcl(i)%pos%y < -boxdh) ptcl(i)%pos%y = ptcl(i)%pos%y + 2.0d0*boxdh
      if (ptcl(i)%pos%y >  boxdh) ptcl(i)%pos%y = ptcl(i)%pos%y - 2.0d0*boxdh
      if (ptcl(i)%pos%z < -boxdh) ptcl(i)%pos%z = ptcl(i)%pos%z + 2.0d0*boxdh
      if (ptcl(i)%pos%z >  boxdh) ptcl(i)%pos%z = ptcl(i)%pos%z - 2.0d0*boxdh
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
   type(fplj), dimension(:), pointer :: ptcl

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
      epot_loc = epot_loc + ptcl(i)%mass * ptcl(i)%pot 
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
   type(fplj), dimension(:), pointer :: ptcl
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

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////     < C B R T >     /////////////////////////
!-----------------------------------------------------------------------
pure double precision function cbrt(a) result(ret)
   implicit none
   double precision, intent(in) :: a
   !* Local variables
   double precision :: x1,x2

   x1 = a
   x2 = a
   do
      x1 = x2
      x2 = (a - x1*x1*x1) / (3*x1*x1) + x1;
      if (abs((x2 - x1) / x1) < 1d-15) exit
   end do
   ret = x2

end function cbrt
