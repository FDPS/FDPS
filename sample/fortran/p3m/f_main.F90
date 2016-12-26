!-----------------------------------------------------------------------
!//////////////////////// M A I N   F U N C T I O N ////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use fdps_module
   use user_defined_types
   implicit none
   !* Local parameters
   !-(# of error measurements)
   integer(kind=c_int), parameter :: num_trials = 512
  !integer(kind=c_int), parameter :: num_trials = 1
   !-(tree)
   real(kind=c_float), parameter :: theta = 0.0
   integer(kind=c_int), parameter :: n_leaf_limit = 8
   integer(kind=c_int), parameter :: n_group_limit = 64
   !-(domain decomposition)
   real(kind=c_float), parameter :: coef_ema = 0.3
   !* Local variables
   integer(kind=c_int) :: i,nstep
   integer(kind=c_int) :: nptcl_1d,nptcl_loc
   integer(kind=c_int) :: psys_num,dinfo_num,tree_num,pm_num
   character(len=64) :: fname
   real(kind=c_double) :: relerr,relerr1
   type(fdps_f32vec) :: pos32
   type(fdps_controller) :: fdps_ctrl
   type(crystal_parameters) :: NaCl_params
   type(nbody_fp), dimension(:), pointer :: ptcl
   type(c_funptr) :: pfunc_ep_ep,pfunc_ep_sp
   !* Static variables
   logical, save :: is_tree_initialized=.false.
   logical, save :: first_output=.true.
   
   !* Initialize FDPS 
   call fdps_ctrl%ps_initialize()

   !* Make an instance of particle system and initialize it
   call fdps_ctrl%create_psys(psys_num,'nbody_fp')
   call fdps_ctrl%init_psys(psys_num)

   !* Make an instance of domain info. and initialize it
   call fdps_ctrl%create_dinfo(dinfo_num)
   call fdps_ctrl%init_dinfo(dinfo_num,coef_ema)

   !* Make an instance of ParticleMesh object
   call fdps_ctrl%create_pm(pm_num)

   !* Initialize Mersenne twister pseudo-random number generator
   call fdps_ctrl%MT_init_genrand(0)

   !==================================================================
   !* Compute relative energy errors of the Madelung energy
   !  due to the P^{3}M method for different # of particles
   !  and for different configurations.
   !==================================================================
   do nptcl_1d=4,32,2
      !* Information to STDOUT
      if (fdps_ctrl%get_rank() == 0) then
         write(*,100)nptcl_1d
         100 format("Processing ",i4," case...")
      end if

      NaCl_params%nptcl_per_side = nptcl_1d
      relerr = 0.0d0
      do nstep=1,num_trials
         !* [1] Randomly choose a configuration of the grid
         NaCl_params%pos_vertex%x = fdps_ctrl%MT_genrand_res53()
         NaCl_params%pos_vertex%y = fdps_ctrl%MT_genrand_res53()
         NaCl_params%pos_vertex%z = fdps_ctrl%MT_genrand_res53()

         !* [2] Make a NaCl crystal
         call setup_NaCl_crystal(fdps_ctrl, &
                                 psys_num,  &
                                 dinfo_num, &
                                 NaCl_params)

         !* [3] Initialize tree if needed
         if (is_tree_initialized .eqv. .false.) then
            nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
            call fdps_ctrl%create_tree(tree_num, &
                                       "Long,nbody_pp_results,nbody_ep,nbody_ep,MonopoleWithCutoff")
            call fdps_ctrl%init_tree(tree_num,3*nptcl_loc,theta, &
                                     n_leaf_limit,n_group_limit)
            is_tree_initialized = .true.
         end if

         !* [4] Compute force and potential with P^{3}M method
         !* [4-1] Get the pointer to FP and # of local particles
         nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
         call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
         !* [4-2] PP part
         pfunc_ep_ep = c_funloc(calc_force_ep_ep)
         pfunc_ep_sp = c_funloc(calc_force_ep_sp)
         call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                      pfunc_ep_ep, &
                                                      pfunc_ep_sp, &
                                                      psys_num,    &
                                                      dinfo_num)
         !* [4-3] PM part
         call fdps_ctrl%calc_pm_force_all_and_write_back(pm_num,   &
                                                         psys_num, &
                                                         dinfo_num)
         do i=1,nptcl_loc
            pos32 = ptcl(i)%x
            call fdps_ctrl%get_pm_potential(pm_num,pos32,ptcl(i)%pot_pm)
         end do
         !* [4-4] Compute the total acceleration and potential
         do i=1,nptcl_loc
            ptcl(i)%pot  = ptcl(i)%pot  - ptcl(i)%pot_pm
            ptcl(i)%agrv = ptcl(i)%agrv - ptcl(i)%agrv_pm
         end do

         !* [5] Compare with the result of the Ewald summation
         call calc_energy_error(fdps_ctrl,psys_num,relerr1)
         relerr = relerr + relerr1

      end do
      relerr = relerr/num_trials

      !* Output relative error
      if (fdps_ctrl%get_rank() == 0) then
         !* Output to file
         fname = "EnergyError.dat"
         if (first_output .eqv. .true.) then
            open(unit=9,file=trim(fname),action='write',status='replace')
            first_output = .false.
         else
            open(unit=9,file=trim(fname),action='write',status='old',position='append')
         end if
            write(9,'(i6,1x,1es25.16e3)')nptcl_1d,relerr
         close(unit=9)
         !* Output to STDOUT
         write(*,200)num_trials,nptcl_1d,relerr
         200 format("********** Result of this experiment **********"/ &
                    "   num_trials     = ",i6/                         &
                    "   nptcl_per_side = ",i6/                         &
                    "   Relative Error = ",1es25.16e3/                 &
                    "***********************************************")
      end if

   end do

   !* Finalize FDPS
   call fdps_ctrl%ps_finalize()

end subroutine f_main

!-----------------------------------------------------------------------
!///////////////           S U B R O U T I N E           ///////////////
!/////////////// < S E T U P _ N A C L _ C R Y S T A L > ///////////////
!-----------------------------------------------------------------------
subroutine setup_NaCl_crystal(fdps_ctrl,psys_num,dinfo_num,params)
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   use mpi
#endif
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer(c_int), intent(in) :: psys_num,dinfo_num
   type(crystal_parameters), intent(in) :: params
   !* Local variables
   integer(c_int) :: i,j,k,ii,irank,ierr
   integer(c_int) :: nprocs,myrank,id
   integer(c_int) :: nptcl,nptcl_loc,nptcl_rem;
   integer(c_int) :: i_start,i_end
   integer(c_int) :: size_of_mesh
   integer(c_int), dimension(:), allocatable :: nptcl_list,sendbuf
   real(c_double) :: cutoff_radius
   real(c_double) :: xmin,xmax,ymin,ymax,zmin,zmax
   real(c_double) :: m,x,y,z
   type(nbody_fp), dimension(:), pointer :: ptcl
   type(fdps_f64vec) :: pos_ll,pos_ul
   !-(IO)
   character(len=64) :: fname

   !* Get # of MPI processes and Rank number
   nprocs = fdps_ctrl%get_num_procs()
   myrank = fdps_ctrl%get_rank()

   !* Define the parameters
   nptcl = params%nptcl_per_side &
         * params%nptcl_per_side &
         * params%nptcl_per_side
   if (mod(params%nptcl_per_side,2) /= 0) then
      if (myrank == 0) then
         write(*,100)params%nptcl_per_side
         100 format("numPtcl_per_side is an invalid value: ",i6)
      end if
      flush(6)
      call fdps_ctrl%ps_finalize()
      stop 1
   end if

   !* Compute nptcl_list
   allocate( nptcl_list(0:nprocs-1) )
   nptcl_list = 0
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   !* MPI-parallel
   nptcl_loc = nptcl/nprocs
   nptcl_rem = mod(nptcl,nprocs)
   if (nptcl_rem /= 0) then
      if ((myrank+1) <= nptcl_rem) then
         nptcl_loc = nptcl_loc + 1
      end if
   end if
   allocate( sendbuf(0:nprocs-1) )
   sendbuf = nptcl_loc
   call MPI_AlltoAll(sendbuf,1,MPI_Integer,    &
                     nptcl_list,1,MPI_Integer, &
                     MPI_COMM_WORLD,ierr)
   deallocate( sendbuf )
   i_start = 1
   if (myrank > 0) then
      do irank=0,myrank-1
         i_start = i_start + nptcl_list(irank)
      end do
   end if
   i_end = i_start + nptcl_loc - 1
#else
   !* Serial
   nptcl_loc = nptcl
   i_start = 1
   i_end   = i_start + nptcl_loc - 1
#endif
   call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_loc)

   !* Make a particle distribution
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   size_of_mesh = fdps_ctrl%get_pm_mesh_num()
   cutoff_radius = fdps_ctrl%get_pm_cutoff_radius()
   id = 1
   do i=1,params%nptcl_per_side
      do j=1,params%nptcl_per_side
         do k=1,params%nptcl_per_side
            !* Substitute
            if ((i_start <= id) .and. (id <= i_end)) then
               ii = id - i_start + 1
               ptcl(ii)%id = id
               if (mod(i+j+k-3,2) == 0) then
                  ptcl(ii)%m =  1.0d0                 &
                             / (params%nptcl_per_side &
                               *params%nptcl_per_side)
               else
                  ptcl(ii)%m = -1.0d0                 &
                             / (params%nptcl_per_side &
                               *params%nptcl_per_side)
               end if
               ptcl(ii)%rc = cutoff_radius/size_of_mesh
               ptcl(ii)%x%x = (params%pos_vertex%x + i-1) * (1.0d0/params%nptcl_per_side)
               ptcl(ii)%x%y = (params%pos_vertex%y + j-1) * (1.0d0/params%nptcl_per_side)
               ptcl(ii)%x%z = (params%pos_vertex%z + k-1) * (1.0d0/params%nptcl_per_side)
               ptcl(ii)%v   = 0.0d0
            end if

            !* Update id
            id = id + 1
         end do
      end do
   end do

   !* Domain decomposition & exchange particle
   call fdps_ctrl%set_boundary_condition(dinfo_num,fdps_bc_periodic_xyz)
   pos_ll%x = 0.0d0; pos_ll%y = 0.0d0; pos_ll%z = 0.0d0
   pos_ul%x = 1.0d0; pos_ul%y = 1.0d0; pos_ul%z = 1.0d0
   call fdps_ctrl%set_pos_root_domain(dinfo_num,pos_ll,pos_ul)
   call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
   call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

end subroutine setup_NaCl_crystal

!-----------------------------------------------------------------------
!////////////////          S U B R O U T I N E          ////////////////
!//////////////// < C A L C _ E N E R G Y _ E R R O R > ////////////////
!-----------------------------------------------------------------------
subroutine calc_energy_error(fdps_ctrl,psys_num,relerr)
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(in) :: fdps_ctrl
   integer(c_int), intent(in) :: psys_num
   real(c_double), intent(inout) :: relerr
   !* Local parameters
   real(c_double), parameter :: e0=-1.7475645946332
   ! where e0 is the analytical solution for \sum q_{i}*\phi_{i}
   ! obtained by the PM^{3} method (computed by K.Nitadori).
   !* Local variables
   integer(c_int) :: i,nptcl_loc
   real(c_double) :: ebuf,esum
   type(nbody_fp), dimension(:), pointer :: ptcl

   !* Get the pointer to FP
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)

   !* Compute total energy
   ebuf=0.0d0; esum=0.0d0
   do i=1,nptcl_loc
      ebuf = ebuf + ptcl(i)%m * ptcl(i)%pot
   end do
   call fdps_ctrl%get_sum(ebuf,esum)
   
   !* Release the pointer
   nullify(ptcl)

   relerr = dabs((esum-e0)/e0)
   
end subroutine calc_energy_error
