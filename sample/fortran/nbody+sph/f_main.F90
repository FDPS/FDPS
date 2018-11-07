#include "macro_defs.h"
!-----------------------------------------------------------------------
!///////////////////// < M A I N   R O U T I N E > /////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use fdps_vector
   use fdps_module
#if defined(ENABLE_PHANTOM_GRAPE_X86)
   use phantom_grape_g5_x86
#endif
   use user_defined_types
   implicit none
   !* Local parameters
   !-(force parameters)
   real, parameter :: theta = 0.5
   integer, parameter :: n_leaf_limit = 8
   integer, parameter :: n_group_limit = 64
   !-(domain decomposition)
   real, parameter :: coef_ema=0.3
   !-(others)
   logical(kind=c_bool), parameter :: clear=.true.
   logical(kind=c_bool), parameter :: unclear=.false.
   !* Local variables
   integer :: i,j,k,ierr
   integer :: offset
   integer :: nstep
   integer :: psys_num_nbody, psys_num_sph
   integer :: dinfo_num
   integer :: tree_num_grav, tree_num_dens, tree_num_hydro
   integer :: nptcl_loc_nbody,nptcl_loc_sph,nptcl_loc_all
   double precision :: time,dt,time_dump,dt_dump,time_end
   double precision :: t_start,t_grav,t_hydro
   double precision :: t_start_1step,t_1step,t_1step_max
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph
   type(c_funptr) :: pfunc_ep_ep, pfunc_ep_sp
   type(force_grav) :: f_grav
   character(len=5) :: proc_num
   character(len=64) :: fname
   !* External routines
   double precision, external :: get_timestep

   !* Initialize FDPS
   call fdps_ctrl%PS_Initialize()

   !* Make instances of ParticleSystems and initialize them
   call fdps_ctrl%create_psys(psys_num_nbody,'fp_nbody')
   call fdps_ctrl%init_psys(psys_num_nbody)
   call fdps_ctrl%create_psys(psys_num_sph,'fp_sph')
   call fdps_ctrl%init_psys(psys_num_sph)

   !* Make an instance of DomainInfo and initialize it
   call fdps_ctrl%create_dinfo(dinfo_num)
   call fdps_ctrl%init_dinfo(dinfo_num,coef_ema)

   !* Make an initial condition and initialize the particle system
   call setup_IC(psys_num_nbody, psys_num_sph, dinfo_num, &
                 time_dump, dt_dump, time_end)

   !* Perform domain decomposition
   call fdps_ctrl%collect_sample_particle(dinfo_num, psys_num_nbody, clear)
   call fdps_ctrl%collect_sample_particle(dinfo_num, psys_num_sph, unclear)
   call fdps_ctrl%decompose_domain(dinfo_num)

   !* Perform particle exchange
   call fdps_ctrl%exchange_particle(psys_num_nbody,dinfo_num)
   call fdps_ctrl%exchange_particle(psys_num_sph,dinfo_num)

   !* Make three tree structures
   nptcl_loc_sph   = max(fdps_ctrl%get_nptcl_loc(psys_num_sph),1)
   nptcl_loc_nbody = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   nptcl_loc_all   = nptcl_loc_nbody + nptcl_loc_sph
   !** tree for gravity calculation
   call fdps_ctrl%create_tree(tree_num_grav, &
                              "Long,force_grav,ep_grav,ep_grav,Monopole")
   call fdps_ctrl%init_tree(tree_num_grav, 3*nptcl_loc_all, theta, &
                            n_leaf_limit, n_group_limit)
   !** tree for the density calculation
   call fdps_ctrl%create_tree(tree_num_dens, &
                              "Short,force_dens,ep_hydro,ep_hydro,Gather")
   call fdps_ctrl%init_tree(tree_num_dens, 3*nptcl_loc_sph, theta, &
                            n_leaf_limit, n_group_limit)

   !** tree for the hydrodynamic force calculation
   call fdps_ctrl%create_tree(tree_num_hydro, &
                              "Short,force_hydro,ep_hydro,ep_hydro,Symmetry")
   call fdps_ctrl%init_tree(tree_num_hydro, 3*nptcl_loc_sph, theta, &
                            n_leaf_limit, n_group_limit)

   !* Perform force calculations
#if defined(ENABLE_PHANTOM_GRAPE_X86)
    call g5_open()
    call g5_set_eps_to_all(eps_grav);
#endif
   !** Gravity calculation
   t_start = fdps_ctrl%get_wtime()
#if defined(ENABLE_GRAVITY_INTERACT)
   call fdps_ctrl%set_particle_local_tree(tree_num_grav, psys_num_nbody)
   call fdps_ctrl%set_particle_local_tree(tree_num_grav, psys_num_sph, unclear)
   pfunc_ep_ep = c_funloc(calc_gravity_ep_ep)
   pfunc_ep_sp = c_funloc(calc_gravity_ep_sp)
   call fdps_ctrl%calc_force_making_tree(tree_num_grav, &
                                         pfunc_ep_ep,   &
                                         pfunc_ep_sp,   &
                                         dinfo_num)
   nptcl_loc_nbody = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody, ptcl_nbody)
   do i=1,nptcl_loc_nbody
       call fdps_ctrl%get_force(tree_num_grav, i, f_grav)
       ptcl_nbody(i)%acc%x = f_grav%acc%x
       ptcl_nbody(i)%acc%y = f_grav%acc%y
       ptcl_nbody(i)%acc%z = f_grav%acc%z
       ptcl_nbody(i)%pot   = f_grav%pot
   end do
   offset = nptcl_loc_nbody
   nptcl_loc_sph = fdps_ctrl%get_nptcl_loc(psys_num_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph, ptcl_sph)
   do i=1,nptcl_loc_sph
       call fdps_ctrl%get_force(tree_num_grav, i + offset, f_grav)
       ptcl_sph(i)%acc_grav%x = f_grav%acc%x
       ptcl_sph(i)%acc_grav%y = f_grav%acc%y
       ptcl_sph(i)%acc_grav%z = f_grav%acc%z
       ptcl_sph(i)%pot_grav   = f_grav%pot
   end do
#endif
   t_grav = fdps_ctrl%get_wtime() - t_start
   !** SPH calculations
   t_start = fdps_ctrl%get_wtime()
#if defined(ENABLE_HYDRO_INTERACT)
   call calc_density_wrapper(psys_num_sph, dinfo_num, tree_num_dens)
   call set_entropy(psys_num_sph)
   call set_pressure(psys_num_sph)
   pfunc_ep_ep = c_funloc(calc_hydro_force)
   call fdps_ctrl%calc_force_all_and_write_back(tree_num_hydro, &
                                                pfunc_ep_ep,    &
                                                psys_num_sph,   &
                                                dinfo_num)
#endif
   t_hydro = fdps_ctrl%get_wtime() - t_start

   !* Set the initial velocity of gas particle
#if defined(SET_CIRCULAR_VELOCITY) 
   call set_circular_velocity(psys_num_sph)
#endif

   !* Get timestep
   dt = get_timestep(psys_num_nbody, psys_num_sph)

   !* Check the conserved variables
   call check_cnsrvd_vars(psys_num_nbody,psys_num_sph,0.0d0)

   !* Main loop for time integration
   nstep = 1; time = 0.0d0
   do 
      if (fdps_ctrl%get_rank() == 0) then
          write(*,200)"nstep = ",nstep, &
                      " dt = ",dt,      &
                      " time = ",time,  &
                      " time_end = ",time_end
          200 format(a,i8,3(a,1es25.16e3))
      end if
      t_start_1step = fdps_ctrl%get_wtime()

      !* Leap frog: Initial Kick & Full Drift
      call initial_kick(psys_num_nbody,psys_num_sph,dt)
      call full_drift(psys_num_nbody,psys_num_sph,dt)
      if (fdps_ctrl%get_boundary_condition(dinfo_num) /= fdps_bc_open) then
         call fdps_ctrl%adjust_pos_into_root_domain(psys_num_nbody,dinfo_num)
         call fdps_ctrl%adjust_pos_into_root_domain(psys_num_sph,dinfo_num)
      end if

      !* Leap frog: Predict
      call predict(psys_num_sph,dt)

      !* Perform domain decomposition again
      call fdps_ctrl%collect_sample_particle(dinfo_num, psys_num_nbody, clear)
      call fdps_ctrl%collect_sample_particle(dinfo_num, psys_num_sph, unclear)
      call fdps_ctrl%decompose_domain(dinfo_num)

      !* Exchange the particles between the (MPI) processes
      call fdps_ctrl%exchange_particle(psys_num_nbody,dinfo_num)
      call fdps_ctrl%exchange_particle(psys_num_sph,dinfo_num)

      !* Perform force calculations
      !** Gravity calculation
      t_start = fdps_ctrl%get_wtime()
#if defined(ENABLE_GRAVITY_INTERACT)
      call fdps_ctrl%set_particle_local_tree(tree_num_grav, psys_num_nbody)
      call fdps_ctrl%set_particle_local_tree(tree_num_grav, psys_num_sph, unclear)
      pfunc_ep_ep = c_funloc(calc_gravity_ep_ep)
      pfunc_ep_sp = c_funloc(calc_gravity_ep_sp)
      call fdps_ctrl%calc_force_making_tree(tree_num_grav, &
                                            pfunc_ep_ep,   &
                                            pfunc_ep_sp,   &
                                            dinfo_num)
      nptcl_loc_nbody = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
      call fdps_ctrl%get_psys_fptr(psys_num_nbody, ptcl_nbody)
      do i=1,nptcl_loc_nbody
          call fdps_ctrl%get_force(tree_num_grav, i, f_grav)
          ptcl_nbody(i)%acc%x = f_grav%acc%x
          ptcl_nbody(i)%acc%y = f_grav%acc%y
          ptcl_nbody(i)%acc%z = f_grav%acc%z
          ptcl_nbody(i)%pot   = f_grav%pot
      end do
      offset = nptcl_loc_nbody
      nptcl_loc_sph = fdps_ctrl%get_nptcl_loc(psys_num_sph)
      call fdps_ctrl%get_psys_fptr(psys_num_sph, ptcl_sph)
      do i=1,nptcl_loc_sph
          call fdps_ctrl%get_force(tree_num_grav, i + offset, f_grav)
          ptcl_sph(i)%acc_grav%x = f_grav%acc%x
          ptcl_sph(i)%acc_grav%y = f_grav%acc%y
          ptcl_sph(i)%acc_grav%z = f_grav%acc%z
          ptcl_sph(i)%pot_grav   = f_grav%pot
      end do
#endif
      t_grav = fdps_ctrl%get_wtime() - t_start
      !** SPH calculations
      t_start = fdps_ctrl%get_wtime()
#if defined(ENABLE_HYDRO_INTERACT)
      call calc_density_wrapper(psys_num_sph, dinfo_num, tree_num_dens)
      call set_pressure(psys_num_sph)
      pfunc_ep_ep = c_funloc(calc_hydro_force)
      call fdps_ctrl%calc_force_all_and_write_back(tree_num_hydro, &
                                                   pfunc_ep_ep,    &
                                                   psys_num_sph,   &
                                                   dinfo_num)
#endif
      t_hydro = fdps_ctrl%get_wtime() - t_start

      !* Get a new timestep
      dt = get_timestep(psys_num_nbody,psys_num_sph)

      !* Leap frog: Final Kick
      call final_kick(psys_num_nbody,psys_num_sph,dt)

      !* Get the elapsed time for this step
      t_1step = fdps_ctrl%get_wtime() - t_start_1step
      call fdps_ctrl%get_max_value(t_1step,t_1step_max)
      if (fdps_ctrl%get_rank() == 0) then
         write(*,*)'t_1step_max = ',t_1step_max
      end if

      !* Output result files
      if (time > time_dump) then
         call output(psys_num_nbody,psys_num_sph)
         if (fdps_ctrl%get_rank() == 0) then
            write(*,*)"============================================"
            write(*,*)"output a file at time = ",time 
            write(*,*)"============================================"
         end if
         time_dump = time_dump + dt_dump
      end if

      !* Check the conserved variables
      call check_cnsrvd_vars(psys_num_nbody,psys_num_sph,time)

      !* Cehck the amplitude of density fluctuation
#if defined(CHECK_DENSITY_FLUCTUATION)
      if (mod(nstep,100) == 0) then
         call check_dens_fluc(psys_num_sph)
      end if
#endif

      !* Termination condition
      if (time >= time_end) exit

      !* Update time & step
      time  = time + dt
      nstep = nstep + 1

      !* For debug
      if (nstep == 100) exit
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
subroutine setup_IC(psys_num_nbody,psys_num_sph,dinfo_num, &
                    time_dump,dt_dump,time_end)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph,dinfo_num
   double precision, intent(inout) :: time_dump,dt_dump,time_end
   !* Local variables
   integer :: i,nptcl_loc,nptcl_glb
   integer(kind=c_int) :: bc
   real(kind=c_double) :: m_sum_loc,m_sum
   type(fdps_f64vec) :: pos_root_domain_low,pos_root_domain_high
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph
   character(len=64) :: fname

   if (fdps_ctrl%get_rank() == 0) then
#if (INITIAL_CONDITION == 0)
     call galaxy_IC(psys_num_nbody,psys_num_sph, &
                    bc,                          &
                    pos_root_domain_low,         &
                    pos_root_domain_high,        &
                    time_dump,dt_dump,time_end)
#elif (INITIAL_CONDITION == 1)
     call cold_collapse_test_IC(psys_num_nbody,psys_num_sph, &
                                bc,                          &
                                pos_root_domain_low,         &
                                pos_root_domain_high,        &
                                time_dump,dt_dump,time_end)
#elif (INITIAL_CONDITION == 2)
     call Evrard_test_IC(psys_num_nbody,psys_num_sph, &
                         bc,                          &
                         pos_root_domain_low,         &
                         pos_root_domain_high,        &
                         time_dump,dt_dump,time_end,1)
#elif (INITIAL_CONDITION == 3)
     call make_glass_IC(psys_num_nbody,psys_num_sph, &
                        bc,                          &
                        pos_root_domain_low,         &
                        pos_root_domain_high,        &
                        time_dump,dt_dump,time_end)
#else
#error Invalid IC number is specified.
#endif
      !* Check the initial distribution
      fname = "IC.txt"
      open(unit=9,file=trim(fname),action='write',status='replace')
         nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
         call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
         do i=1,nptcl_loc
            write(9,'(3es25.16e3)')ptcl_nbody(i)%pos%x, &
                                   ptcl_nbody(i)%pos%y, &
                                   ptcl_nbody(i)%pos%z
         end do
         nullify(ptcl_nbody)
         write(9,*)''
         write(9,*)''
         nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_sph)
         call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
         do i=1,nptcl_loc
            write(9,'(3es25.16e3)')ptcl_sph(i)%pos%x, &
                                   ptcl_sph(i)%pos%y, &
                                   ptcl_sph(i)%pos%z
         end do
         nullify(ptcl_sph)
      close(unit=9)
   else 
      call fdps_ctrl%set_nptcl_loc(psys_num_nbody,0)
      call fdps_ctrl%set_nptcl_loc(psys_num_sph,0)
   end if

   !* Broadcast from RANK 0 
   call fdps_ctrl%broadcast(bc, 1, 0)
   call fdps_ctrl%broadcast(pos_root_domain_low%x, 1, 0)
   call fdps_ctrl%broadcast(pos_root_domain_low%y, 1, 0)
   call fdps_ctrl%broadcast(pos_root_domain_low%z, 1, 0)
   call fdps_ctrl%broadcast(pos_root_domain_high%x, 1, 0)
   call fdps_ctrl%broadcast(pos_root_domain_high%y, 1, 0)
   call fdps_ctrl%broadcast(pos_root_domain_high%z, 1, 0)
   call fdps_ctrl%broadcast(eps_grav, 1, 0)
   call fdps_ctrl%broadcast(time_dump, 1, 0)
   call fdps_ctrl%broadcast(dt_dump, 1, 0)
   call fdps_ctrl%broadcast(time_end, 1, 0)
   call fdps_ctrl%broadcast(dt_max, 1, 0)
   
   !* Set the boundary condition and the size of the computational domain if needed.
   call fdps_ctrl%set_boundary_condition(dinfo_num,bc)
   if (bc /= fdps_bc_open) then
       call fdps_ctrl%set_pos_root_domain(dinfo_num,           &
                                          pos_root_domain_low, &
                                          pos_root_domain_high)
   end if

   !* Compute the average mass of SPH particles
   nptcl_glb = fdps_ctrl%get_nptcl_glb(psys_num_sph)
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   m_sum_loc = 0.0d0; m_sum = 0.0d0
   do i=1,nptcl_loc
      m_sum_loc = m_sum_loc + ptcl_sph(i)%mass
   end do
   call fdps_ctrl%get_sum(m_sum_loc, m_sum)
   mass_avg = m_sum / nptcl_glb
   nullify(ptcl_sph)

   !* Inform to STDOUT
   if (fdps_ctrl%get_rank() == 0) then
      write(*,*)"setup_IC() is completed."
   end if
   !call fdps_ctrl%ps_finalize()
   !stop 0

end subroutine setup_IC

!-----------------------------------------------------------------------
!/////////////////////     S U B R O U T I N E     /////////////////////
!///////////////////// < G E T _ T I M E S T E P > /////////////////////
!-----------------------------------------------------------------------
function get_timestep(psys_num_nbody,psys_num_sph)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   real(kind=c_double) :: get_timestep
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   !* Local variables
   integer :: i,nptcl_loc
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph
   type(fdps_f64vec) :: acc_tot
   real(kind=c_double) :: dt_loc,acc

   dt_loc = huge(0.0d0)
   if (dt_max > 0.0d0) dt_loc = dt_max

   !* timescale for N-body system
#if defined(ENABLE_GRAVITY_INTERACT)
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   do i=1,nptcl_loc
      acc = dsqrt(ptcl_nbody(i)%acc%x * ptcl_nbody(i)%acc%x &
                 +ptcl_nbody(i)%acc%y * ptcl_nbody(i)%acc%y &
                 +ptcl_nbody(i)%acc%z * ptcl_nbody(i)%acc%z)
      if (acc > 0.0) then
         dt_loc = min(dt_loc, CFL_dyn * dsqrt(eps_grav / acc))
      end if
   end do
#endif

   !* Timescale for SPH system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   do i=1,nptcl_loc
#if defined(ENABLE_GRAVITY_INTERACT)
      acc_tot%x = ptcl_sph(i)%acc_grav%x + ptcl_sph(i)%acc_hydro%x
      acc_tot%y = ptcl_sph(i)%acc_grav%y + ptcl_sph(i)%acc_hydro%y
      acc_tot%z = ptcl_sph(i)%acc_grav%z + ptcl_sph(i)%acc_hydro%z
      acc = dsqrt(acc_tot%x * acc_tot%x &
                 +acc_tot%y * acc_tot%y &
                 +acc_tot%z * acc_tot%z)
      if (acc > 0.0) then
         dt_loc = min(dt_loc, CFL_dyn * dsqrt(eps_grav / acc))
      end if
#endif
#if defined(ENABLE_HYDRO_INTERACT)
      dt_loc = min(dt_loc, ptcl_sph(i)%dt)
#endif
   end do
 
   !* Reduction
   call fdps_ctrl%get_min_value(dt_loc,get_timestep)

   !* Release the pointers
   nullify(ptcl_nbody)
   nullify(ptcl_sph)

end function get_timestep

!-----------------------------------------------------------------------
!/////////////             S U B R O U T I N E             /////////////
!///////////// < C A L C _ D E N S I T Y _ W R A P P E R > /////////////
!-----------------------------------------------------------------------
subroutine calc_density_wrapper(psys_num,dinfo_num,tree_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num,dinfo_num,tree_num
   !* Local variables
   integer :: i,nptcl_loc,nptcl_glb
   integer :: n_compl_loc,n_compl
   type(fdps_controller) :: fdps_ctrl
   type(fp_sph), dimension(:), pointer :: ptcl
   type(c_funptr) :: pfunc_ep_ep

#if defined(ENABLE_VARIABLE_SMOOTHING_LENGTH)
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   nptcl_glb = fdps_ctrl%get_nptcl_glb(psys_num)
   call fdps_ctrl%get_psys_fptr(psys_num, ptcl)
   pfunc_ep_ep = c_funloc(calc_density)
   ! Determine the density and the smoothing length
   ! so that Eq.(6) in Springel (2005) holds within a specified accuracy.
   do 
       ! Increase smoothing length 
       do i=1,nptcl_loc
           ptcl(i)%smth = scf_smth * ptcl(i)%smth 
       end do
       ! Compute density, etc.
       call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    psys_num,    &
                                                    dinfo_num)
       ! Check convergence
       n_compl_loc = 0; n_compl = 0
       do i=1,nptcl_loc
           if (ptcl(i)%flag == 1) n_compl_loc = n_compl_loc + 1
       end do
       call fdps_ctrl%get_sum(n_compl_loc, n_compl)
       if (n_compl == nptcl_glb) exit
   end do
   !* Release the pointer
   nullify(ptcl)
#else
   pfunc_ep_ep = c_funloc(calc_density)
   call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                pfunc_ep_ep, &
                                                psys_num,    &
                                                dinfo_num)
#endif

end subroutine calc_density_wrapper

!-----------------------------------------------------------------------
!/////////////////////     S U B R O U T I N E     /////////////////////
!///////////////////// < S E T _ P R E S S U R E > /////////////////////
!-----------------------------------------------------------------------
subroutine set_pressure(psys_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num
   !* Local variables
   integer :: i,nptcl_loc
   type(fdps_controller) :: fdps_ctrl
   type(fp_sph), dimension(:), pointer :: ptcl

   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
#if defined(ISOTHERMAL_EOS)
      ! In this case, eng = const.
      ptcl(i)%pres = (specific_heat_ratio - 1.0d0) * ptcl(i)%dens * ptcl(i)%eng
      ptcl(i)%ent  = ptcl(i)%pres / ptcl(i)%dens**(specific_heat_ratio)
#else
#if defined(USE_ENTROPY)
      ptcl(i)%pres = ptcl(i)%ent * ptcl(i)%dens**(specific_heat_ratio)
      ptcl(i)%eng  = ptcl(i)%pres / ((specific_heat_ratio - 1.0d0) * ptcl(i)%dens)
#else
      ptcl(i)%pres = (specific_heat_ratio - 1.0d0) * ptcl(i)%dens * ptcl(i)%eng
      ptcl(i)%ent  = ptcl(i)%pres / ptcl(i)%dens**(specific_heat_ratio)
#endif
#endif
      ptcl(i)%snds = dsqrt(specific_heat_ratio * ptcl(i)%pres / ptcl(i)%dens)
#if defined(USE_BALSARA_SWITCH)
      ptcl(i)%balsw = abs(ptcl(i)%divv) / (abs(ptcl(i)%divv)                  &
                                          +dsqrt(ptcl(i)%rotv * ptcl(i)%rotv) &
                                          +1.0d-4 * ptcl(i)%snds / ptcl(i)%smth)
#else
      ptcl(i)%balsw = 1.0d0
#endif
   end do
   nullify(ptcl)

end subroutine set_pressure

!-----------------------------------------------------------------------
!//////////////////////    S U B R O U T I N E    //////////////////////
!////////////////////// < S E T _ E N T R O P Y > //////////////////////
!-----------------------------------------------------------------------
subroutine set_entropy(psys_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num
   !* Local variables
   integer :: i,nptcl_loc
   type(fdps_controller) :: fdps_ctrl
   type(fp_sph), dimension(:), pointer :: ptcl

   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%ent = (specific_heat_ratio - 1.0d0) * ptcl(i)%eng &
                  / ptcl(i)%dens**(specific_heat_ratio - 1.0d0)
   end do
   nullify(ptcl)

end subroutine set_entropy

!-----------------------------------------------------------------------
!//////////////////////    S U B R O U T I N E    //////////////////////
!////////////////////// < S E T _ E N T R O P Y > //////////////////////
!-----------------------------------------------------------------------
subroutine set_circular_velocity(psys_num)
   use fdps_vector
   use fdps_module
   use mathematical_constants
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num
   !* Local variables
   integer :: i,nptcl_loc
   double precision :: r,rcyl,phi,theta
   double precision :: vel_circ,T_rot
   character(len=5) :: file_num
   character(len=64) :: file_name
   type(fdps_f64vec) :: acc,base_vect_r
   type(fdps_controller) :: fdps_ctrl
   type(fp_sph), dimension(:), pointer :: ptcl

   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      acc%x = ptcl(i)%acc_grav%x + ptcl(i)%acc_hydro%x
      acc%y = ptcl(i)%acc_grav%y + ptcl(i)%acc_hydro%y
      acc%z = ptcl(i)%acc_grav%z + ptcl(i)%acc_hydro%z
      r    = dsqrt(ptcl(i)%pos%x * ptcl(i)%pos%x &
                  +ptcl(i)%pos%y * ptcl(i)%pos%y &
                  +ptcl(i)%pos%z * ptcl(i)%pos%z)
      rcyl = dsqrt(ptcl(i)%pos%x * ptcl(i)%pos%x &
                  +ptcl(i)%pos%y * ptcl(i)%pos%y)
      phi = datan2(ptcl(i)%pos%y, ptcl(i)%pos%x)
      theta = datan2(rcyl, ptcl(i)%pos%z)
      base_vect_r%x = sin(theta) * cos(phi)
      base_vect_r%y = sin(theta) * sin(phi)
      base_vect_r%z = cos(theta)
      vel_circ = dsqrt(r * abs(acc%x * base_vect_r%x &
                              +acc%y * base_vect_r%y &
                              +acc%z * base_vect_r%z))
      ptcl(i)%vel%x = - vel_circ * sin(phi)
      ptcl(i)%vel%y =   vel_circ * cos(phi)
      ptcl(i)%vel%z = 0.0
   end do

#if 0
   ! [for debug] 
   ! Output the velocity field on the x-y plane
   write(file_num,"(i5.5)")fdps_ctrl%get_rank()
   file_name = "velc_fld" // file_num // ".txt"
   open(unit=9,file=trim(file_name),action='write',status='replace')
      do i=1,nptcl_loc
         write(9,'(4es25.16)')ptcl(i)%pos%x,ptcl(i)%pos%y, &
                              ptcl(i)%vel%x,ptcl(i)%pos%y
      end do
   close(unit=9)
   ! Output the rotation curve
   file_name = "rot_curve" // file_num // ".txt"
   open(unit=9,file=trim(file_name),action='write',status='replace')
      do i=1,nptcl_loc
         r    = dsqrt(ptcl(i)%pos * ptcl(i)%pos)
         rcyl = dsqrt(ptcl(i)%pos%x * ptcl(i)%pos%x &
                     +ptcl(i)%pos%y * ptcl(i)%pos%y)
         vel_circ = dsqrt(ptcl(i)%vel * ptcl(i)%vel)
         T_rot = 2.0 * pi * r / vel_circ
         write(9,'(3es25.16e3)')rcyl,vel_circ,T_rot
      end do
   close(unit=9)
   !call fdps_ctrl%ps_finalize()
   !stop 0
#endif

   nullify(ptcl)

end subroutine set_circular_velocity

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////   < O U T P U T >   /////////////////////////
!-----------------------------------------------------------------------
subroutine output(psys_num_nbody,psys_num_sph)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(IN) :: psys_num_nbody,psys_num_sph
   !* Local parameters
   character(len=16), parameter :: root_dir="result"
   character(len=16), parameter :: nbody_file_prefix_1st="nbody"
   character(len=16), parameter :: sph_file_prefix_1st="sph"
   character(len=16), parameter :: file_prefix_2nd="proc"
   !* Local variables
   integer :: i,nptcl_loc
   integer :: myrank
   character(len=5) :: file_num,proc_num
   character(len=64) :: cmd,sub_dir,fname
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph
   !* Static variables
   integer, save :: ndump=1

   myrank = fdps_ctrl%get_rank()

   !* Output N-body data
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   write(file_num,"(i5.5)")ndump
   write(proc_num,"(i5.5)")myrank
   fname =  trim(root_dir) // "/" &
         // trim(nbody_file_prefix_1st) // file_num // "-" &
         // trim(file_prefix_2nd) // proc_num // ".dat"
   open(unit=9,file=trim(fname),action='write',status='replace')
      do i=1,nptcl_loc
         write(9,100)ptcl_nbody(i)%id,ptcl_nbody(i)%mass, &
                     ptcl_nbody(i)%pos%x,ptcl_nbody(i)%pos%y,ptcl_nbody(i)%pos%z, &
                     ptcl_nbody(i)%vel%x,ptcl_nbody(i)%vel%y,ptcl_nbody(i)%vel%z
         100 format(i8,1x,7e25.16e3)
      end do
   close(unit=9)
   nullify(ptcl_nbody)

   !* Output SPH data
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   write(file_num,"(i5.5)")ndump
   write(proc_num,"(i5.5)")myrank
   fname =  trim(root_dir) // "/" &
         // trim(sph_file_prefix_1st) // file_num // "-" &
         // trim(file_prefix_2nd) // proc_num // ".dat"
   open(unit=9,file=trim(fname),action='write',status='replace')
      do i=1,nptcl_loc
         write(9,200)ptcl_sph(i)%id,ptcl_sph(i)%mass, &
                     ptcl_sph(i)%pos%x,ptcl_sph(i)%pos%y,ptcl_sph(i)%pos%z, &
                     ptcl_sph(i)%vel%x,ptcl_sph(i)%vel%y,ptcl_sph(i)%vel%z, &
                     ptcl_sph(i)%dens,ptcl_sph(i)%eng,ptcl_sph(i)%pres
         200 format(i8,1x,10e25.16e3)
      end do
   close(unit=9)
   nullify(ptcl_sph)

   !* Update ndump
   ndump = ndump + 1

end subroutine output

!-----------------------------------------------------------------------
!////////////////          S U B R O U T I N E          ////////////////
!//////////////// < C H E C K _ C N S R V D _ V A R S > ////////////////
!-----------------------------------------------------------------------
subroutine check_cnsrvd_vars(psys_num_nbody,psys_num_sph,time)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   integer, intent(in) :: psys_num_nbody,psys_num_sph
   double precision, intent(in) :: time
   !* Local variables
   integer :: i,nptcl_loc
   type(fdps_controller) :: fdps_ctrl
   type(fp_nbody), dimension(:), pointer :: ptcl_nbody
   type(fp_sph), dimension(:), pointer :: ptcl_sph
   type(fdps_f64vec) :: mom_loc,mom
   real(kind=c_double) :: ekin_loc,epot_loc,eth_loc
   real(kind=c_double) :: ekin,epot,eth
   real(kind=c_double) :: emech,etot
   real(kind=c_double) :: relerr_mech,relerr_tot
   !* Static variables
   logical, save :: is_initialized = .false.
   real(kind=c_double), save :: emech_ini,etot_ini

   ekin_loc = 0.0d0
   epot_loc = 0.0d0
   eth_loc  = 0.0d0
   mom_loc  = 0.0d0
   !* Sum over N-body system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_nbody)
   call fdps_ctrl%get_psys_fptr(psys_num_nbody,ptcl_nbody)
   do i=1,nptcl_loc
      ekin_loc = ekin_loc &
               + 0.5d0 * ptcl_nbody(i)%mass * (ptcl_nbody(i)%vel%x * ptcl_nbody(i)%vel%x &
                                              +ptcl_nbody(i)%vel%y * ptcl_nbody(i)%vel%y &
                                              +ptcl_nbody(i)%vel%z * ptcl_nbody(i)%vel%z)
      epot_loc = epot_loc + 0.5d0 * ptcl_nbody(i)%mass * (ptcl_nbody(i)%pot &
                                                         +ptcl_nbody(i)%mass / eps_grav)
      mom_loc%x = mom_loc%x + ptcl_nbody(i)%mass * ptcl_nbody(i)%vel%x
      mom_loc%y = mom_loc%y + ptcl_nbody(i)%mass * ptcl_nbody(i)%vel%y
      mom_loc%z = mom_loc%z + ptcl_nbody(i)%mass * ptcl_nbody(i)%vel%z
   end do
   nullify(ptcl_nbody)
   !* Sum over SPH system
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num_sph)
   call fdps_ctrl%get_psys_fptr(psys_num_sph,ptcl_sph)
   do i=1,nptcl_loc
      ekin_loc = ekin_loc &
               + 0.5d0 * ptcl_sph(i)%mass * (ptcl_sph(i)%vel%x * ptcl_sph(i)%vel%x &
                                            +ptcl_sph(i)%vel%y * ptcl_sph(i)%vel%y &
                                            +ptcl_sph(i)%vel%z * ptcl_sph(i)%vel%z)
      epot_loc = epot_loc + 0.5d0 * ptcl_sph(i)%mass * (ptcl_sph(i)%pot_grav &
                                                       +ptcl_sph(i)%mass / eps_grav)
      eth_loc  = eth_loc  + ptcl_sph(i)%mass * ptcl_sph(i)%eng
      mom_loc%x = mom_loc%x + ptcl_sph(i)%mass * ptcl_sph(i)%vel%x
      mom_loc%y = mom_loc%y + ptcl_sph(i)%mass * ptcl_sph(i)%vel%y
      mom_loc%z = mom_loc%z + ptcl_sph(i)%mass * ptcl_sph(i)%vel%z
   end do
   nullify(ptcl_sph)

   !* Reduction 
   call fdps_ctrl%get_sum(ekin_loc,ekin)
   call fdps_ctrl%get_sum(epot_loc,epot)
   call fdps_ctrl%get_sum(eth_loc,eth)
   call fdps_ctrl%get_sum(mom_loc%x,mom%x)
   call fdps_ctrl%get_sum(mom_loc%y,mom%y)
   call fdps_ctrl%get_sum(mom_loc%z,mom%z)

   if (is_initialized .eqv. .false.) then
      emech_ini = ekin + epot
      etot_ini  = ekin + epot + eth
      is_initialized = .true.
   end if

   !* Output
   if (fdps_ctrl%get_rank() == 0) then
      emech = ekin + epot
      etot  = ekin + epot + eth
      relerr_mech = abs((emech - emech_ini)/emech_ini)
      relerr_tot  = abs((etot  - etot_ini)/etot_ini)
      write(*,100)"-------------------------"
      write(*,200)"E_kin  =",ekin
      write(*,200)"E_pot  =",epot
      write(*,200)"E_th   =",eth
      write(*,300)"E_mech =",emech,"(time =",time,"rel.err. =",relerr_mech,")"
      write(*,300)"E_tot  =",etot,"(time =",time,"rel.err. =",relerr_tot,")"
      write(*,200)"Mom (x) =",mom%x
      write(*,200)"Mom (y) =",mom%y
      write(*,200)"Mom (z) =",mom%z
      write(*,100)"-------------------------"
      100 format(a)
      200 format(a,1x,1es25.16e3)
      300 format(3(a,1x,1es25.16e3,1x),a)
   end if

end subroutine check_cnsrvd_vars

!-----------------------------------------------------------------------
!//////////////////        S U B R O U T I N E        //////////////////
!////////////////// < C H E C K _ D E N S _ F L U C > //////////////////
!-----------------------------------------------------------------------
subroutine check_dens_fluc(psys_num)
    use fdps_vector
    use fdps_module
    use user_defined_types 
    implicit none
    integer, intent(in) :: psys_num
    !* Local parameters
    character(len=16), parameter :: root_dir="result"
    double precision, parameter :: eps = 3.05d-3
    !* Local variables
    integer :: i
    integer :: nptcl_loc,nptcl_glb
    double precision :: tmp
    double precision :: dens_avg,dens_disp
    double precision :: fluc_str
    character(len=5) :: file_num
    character(len=64) :: file_name
    type(fdps_controller) :: fdps_ctrl 
    type(fp_sph), dimension(:), pointer :: ptcl

    nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
    nptcl_glb = fdps_ctrl%get_nptcl_glb(psys_num)
    call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
    ! Compute the average of density
    tmp = 0.0
    do i=1,nptcl_loc
        tmp = tmp + ptcl(i)%dens
    end do
    call fdps_ctrl%get_sum(tmp,dens_avg)
    dens_avg = dens_avg / nptcl_glb
    ! Compute the dispersion of density
    tmp = 0.0
    do i=1,nptcl_loc
        tmp = tmp + (ptcl(i)%dens - dens_avg)**(2.0d0)
    end do
    call fdps_ctrl%get_sum(tmp,dens_disp)
    dens_disp = dsqrt(dens_disp/nptcl_glb)
    ! Output the status 
    fluc_str = dens_disp / dens_avg
    if (fdps_ctrl%get_rank() == 0) then
        write(*,*)"###########################"
        write(*,*)"avg.       = ",dens_avg
        write(*,*)"disp.      = ",dens_disp
        write(*,*)"disp./avg. = ",fluc_str
        write(*,*)"###########################"
    end if
    ! Output data and end the simulation if the fluctuation is small
    if (fluc_str < eps) then
        if (fdps_ctrl%get_rank() == 0) then
           file_name = trim(root_dir) // "/" &
                     // "glass_data_header.dat"
           open(unit=9,file=trim(file_name),action='write', &
                form='unformatted',access='stream',status='replace')
              write(9)nptcl_glb
           close(unit=9)
        end if
        write(file_num,"(i5.5)")fdps_ctrl%get_rank()
        file_name = trim(root_dir) // "/" &
                  // "glass_data" // file_num // ".dat"
        open(unit=9,file=trim(file_name),action='write', &
            form='unformatted',access='stream',status='replace')
            do i=1,nptcl_loc
                write(9)ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z
            end do
        close(unit=9)
        if (fdps_ctrl%get_rank() == 0) then
            write(*,*)"A glass-like distribution is obtained."
            write(*,*)"The particle position data is output as files ", &
                      trim(file_name),", etc."
        end if
        call fdps_ctrl%ps_finalize()
        stop 0
    end if

end subroutine check_dens_fluc
