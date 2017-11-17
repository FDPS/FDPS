!==================================
!    MODULE: FDPS module
!==================================
module FDPS_module
   use, intrinsic :: iso_c_binding
   implicit none

   !* Private parameters
   !**** space dimension
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
   integer, parameter, private :: space_dim=3
#else
   integer, parameter, private :: space_dim=2
#endif

   !* Enum types
   !**** PS::BOUNDARY_CONDITION
   enum, bind(c)
      enumerator :: fdps_bc_open
      enumerator :: fdps_bc_periodic_x
      enumerator :: fdps_bc_periodic_y
      enumerator :: fdps_bc_periodic_z
      enumerator :: fdps_bc_periodic_xy
      enumerator :: fdps_bc_periodic_xz
      enumerator :: fdps_bc_periodic_yz
      enumerator :: fdps_bc_periodic_xyz
      enumerator :: fdps_bc_shearing_box
      enumerator :: fdps_bc_user_defined
   end enum
   !**** PS::INTERACTION_LIST_MODE
   enum, bind(c)
      enumerator :: fdps_make_list
      enumerator :: fdps_make_list_for_reuse
      enumerator :: fdps_reuse_list
   end enum

   !**** FDPS controller
   type, public :: FDPS_controller
      ! Note that this structure is introduced to avoid that
      ! users implement `external` statements in their codes.
   contains
      ! Methods
      !-(basic functions)
      procedure :: PS_initialize
      procedure :: PS_finalize
      procedure :: PS_abort
      !-(ParticleSystem functions)
      procedure :: create_psys
      procedure :: delete_psys
      procedure :: init_psys
      procedure :: get_psys_info
      procedure :: get_psys_memsize
      procedure :: get_psys_time_prof
      procedure :: clear_psys_time_prof
      procedure :: set_nptcl_smpl
      procedure :: set_nptcl_loc
      procedure :: get_nptcl_loc
      procedure :: get_nptcl_glb
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of get_psys_fptr*() 
      !           are generated.
      ! fdps-autogen:get_psys_fptr:method;
      !-----------------------------------------------
      procedure :: exchange_particle
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of add_particle*() 
      !           are generated.
      ! fdps-autogen:add_particle:method;
      !-----------------------------------------------
      procedure :: remove_particle
      procedure :: adjust_pos_into_root_domain
      procedure :: sort_particle
      !-(DomainInfo functions)
      procedure :: create_dinfo
      procedure :: delete_dinfo
      procedure :: init_dinfo
      procedure :: get_dinfo_time_prof
      procedure :: clear_dinfo_time_prof
      procedure :: set_nums_domain
      procedure :: set_boundary_condition
      procedure, private :: set_pos_root_domain_a32
      procedure, private :: set_pos_root_domain_a64
      procedure, private :: set_pos_root_domain_v32
      procedure, private :: set_pos_root_domain_v64
      generic :: set_pos_root_domain => set_pos_root_domain_a32, &
                                        set_pos_root_domain_a64, &
                                        set_pos_root_domain_v32, &
                                        set_pos_root_domain_v64
      procedure :: collect_sample_particle
      procedure :: decompose_domain
      procedure :: decompose_domain_all
      !-(Tree functions)
      procedure :: create_tree
      procedure :: delete_tree
      procedure :: init_tree
      procedure :: get_tree_info
      procedure :: get_tree_memsize
      procedure :: get_tree_time_prof
      procedure :: clear_tree_time_prof
      procedure :: get_num_interact_ep_ep_loc
      procedure :: get_num_interact_ep_sp_loc
      procedure :: get_num_interact_ep_ep_glb
      procedure :: get_num_interact_ep_sp_glb
      procedure :: clear_num_interact
      procedure :: get_num_tree_walk_loc
      procedure :: get_num_tree_walk_glb
      procedure, private :: calc_force_all_and_write_back_s
      procedure, private :: calc_force_all_and_write_back_l
      generic   :: calc_force_all_and_write_back => calc_force_all_and_write_back_s, &
                                                    calc_force_all_and_write_back_l
      procedure, private :: calc_force_all_s
      procedure, private :: calc_force_all_l
      generic   :: calc_force_all => calc_force_all_s, &
                                     calc_force_all_l
      procedure, private :: calc_force_making_tree_s
      procedure, private :: calc_force_making_tree_l
      generic   :: calc_force_making_tree => calc_force_making_tree_s, &
                                             calc_force_making_tree_l
      procedure, private :: calc_force_and_write_back_s
      procedure, private :: calc_force_and_write_back_l
      generic   :: calc_force_and_write_back => calc_force_and_write_back_s, &
                                                calc_force_and_write_back_l
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of get_neighbor_list*() 
      !           are generated.
      ! fdps-autogen:get_neighbor_list:method;
      !-----------------------------------------------
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of get_epj_from_id*()
      !           are generated.
      ! fdps-autogen:get_epj_from_id:method;
      !-----------------------------------------------
      !-(MPI comm. functions)
      procedure :: get_rank 
      procedure :: get_rank_multi_dim
      procedure :: get_num_procs
      procedure :: get_num_procs_multi_dim
      procedure :: get_logical_and
      procedure :: get_logical_or
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of reduction
      !           routines such as get_min_value*
      !           are generated.
      ! fdps-autogen:get_min_value:method;
      ! fdps-autogen:get_max_value:method;
      ! fdps-autogen:get_sum:method;
      ! fdps-autogen:broadcast:method;
      !-----------------------------------------------
      procedure :: get_wtime

      !-(Utility functions)
      procedure :: MT_init_genrand
      procedure :: MT_genrand_int31
      procedure :: MT_genrand_real1
      procedure :: MT_genrand_real2
      procedure :: MT_genrand_res53

      !-(Particle Mesh functions)
#if PARTICLE_SIMULATOR_USE_PM_MODULE
      procedure :: create_pm
      procedure :: delete_pm
      procedure :: get_pm_mesh_num
      procedure :: get_pm_cutoff_radius
      procedure :: set_dinfo_of_pm 
      procedure :: set_psys_of_pm  
      procedure, private :: get_pm_force_a32
      procedure, private :: get_pm_force_a64
      procedure, private :: get_pm_force_v32
      procedure, private :: get_pm_force_v64
      generic :: get_pm_force => get_pm_force_a32, &
                                 get_pm_force_a64, &
                                 get_pm_force_v32, &
                                 get_pm_force_v64
      procedure, private :: get_pm_potential_a32
      procedure, private :: get_pm_potential_a64
      procedure, private :: get_pm_potential_v32
      procedure, private :: get_pm_potential_v64
      generic :: get_pm_potential => get_pm_potential_a32, &
                                     get_pm_potential_a64, &
                                     get_pm_potential_v32, &
                                     get_pm_potential_v64
      procedure :: calc_pm_force_only
      procedure :: calc_pm_force_all_and_write_back
#endif
   end type FDPS_controller

   !* Private routines
   private :: PS_initialize
   private :: PS_finalize
   private :: PS_abort
   !-(ParticleSystem functions)
   private :: create_psys
   private :: delete_psys
   private :: init_psys
   private :: get_psys_info
   private :: get_psys_memsize
   private :: get_psys_time_prof
   private :: clear_psys_time_prof
   private :: set_nptcl_smpl
   private :: set_nptcl_loc
   private :: get_nptcl_loc
   private :: get_nptcl_glb
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           get_psys_fptr*() are declared.
   ! fdps-autogen:get_psys_fptr:decl;
   !-------------------------------------------
   private :: exchange_particle
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           add_particle*() are declared.
   ! fdps-autogen:add_particle:decl;
   !-------------------------------------------
   private :: remove_particle
   private :: adjust_pos_into_root_domain
   private :: sort_particle
   !-(DomainInfo functions)
   private :: create_dinfo
   private :: delete_dinfo
   private :: init_dinfo
   private :: get_dinfo_time_prof
   private :: clear_dinfo_time_prof
   private :: set_nums_domain
   private :: set_boundary_condition
   private :: set_pos_root_domain_a32
   private :: set_pos_root_domain_a64
   private :: set_pos_root_domain_v32
   private :: set_pos_root_domain_v64
   private :: collect_sample_particle
   private :: decompose_domain
   private :: decompose_domain_all
   !-(Tree functions)
   private :: create_tree
   private :: delete_tree
   private :: init_tree
   private :: get_tree_info
   private :: get_tree_memsize
   private :: get_tree_time_prof
   private :: clear_tree_time_prof
   private :: get_num_interact_ep_ep_loc
   private :: get_num_interact_ep_sp_loc
   private :: get_num_interact_ep_ep_glb
   private :: get_num_interact_ep_sp_glb
   private :: clear_num_interact
   private :: get_num_tree_walk_loc
   private :: get_num_tree_walk_glb
   private :: calc_force_all_and_write_back_s
   private :: calc_force_all_and_write_back_l
   private :: calc_force_all_s
   private :: calc_force_all_l
   private :: calc_force_making_tree_s
   private :: calc_force_making_tree_l
   private :: calc_force_and_write_back_s
   private :: calc_force_and_write_back_l
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           get_neighbor_list*() are declared.
   ! fdps-autogen:get_neighbor_list:decl;
   !-------------------------------------------
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           get_epj_from_id*() are declared.
   ! fdps-autogen:get_epj_from_id:decl;
   !-------------------------------------------
   !-(MPI comm. functions)
   private :: get_rank 
   private :: get_rank_multi_dim
   private :: get_num_procs
   private :: get_num_procs_multi_dim
   private :: get_logical_and
   private :: get_logical_or
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           for reduction operations are declared.
   ! fdps-autogen:get_min_value:decl;
   ! fdps-autogen:get_max_value:decl;
   ! fdps-autogen:get_sum:decl;
   ! fdps-autogen:broadcast:decl;
   !-------------------------------------------
   private :: get_wtime
   !-(Utility functions)
   private :: MT_init_genrand
   private :: MT_genrand_int31
   private :: MT_genrand_real1
   private :: MT_genrand_real2
   private :: MT_genrand_res53
   !-(ParticleMesh functions)
#ifdef PARTICLE_SIMULATOR_USE_PM_MODULE
   private :: create_pm
   private :: delete_pm
   private :: get_pm_mesh_num
   private :: get_pm_cutoff_radius
   private :: set_dinfo_of_pm 
   private :: set_psys_of_pm  
   private :: get_pm_force_a32
   private :: get_pm_force_a64
   private :: get_pm_force_v32
   private :: get_pm_force_v64
   private :: get_pm_potential_a32
   private :: get_pm_potential_a64
   private :: get_pm_potential_v32
   private :: get_pm_potential_v64
   private :: calc_pm_force_only
   private :: calc_pm_force_all_and_write_back
#endif
   !-(Error)
   private :: print_errmsg

   !* C++ function interfaces
   interface
      !--------------------------
      !  Initializer/Finalizer
      !--------------------------
      subroutine fdps_initialize() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
      end subroutine fdps_initialize

      subroutine fdps_finalize() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
      end subroutine fdps_finalize

      subroutine fdps_abort(err_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: err_num
      end subroutine fdps_abort

      !----------------------
      !  Particle System 
      !----------------------
      subroutine fdps_create_psys(psys_num,psys_info) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(inout) :: psys_num
         character(kind=c_char), dimension(*), intent(in) :: psys_info
      end subroutine fdps_create_psys

      subroutine fdps_delete_psys(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
      end subroutine fdps_delete_psys

      subroutine fdps_init_psys(psys_num)  bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
      end subroutine fdps_init_psys

      subroutine fdps_get_psys_info(psys_num,psys_info,charlen)  bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
         character(kind=c_char), dimension(*), intent(inout) :: psys_info
         integer(kind=c_size_t), intent(inout) :: charlen
      end subroutine fdps_get_psys_info

      function fdps_get_psys_memsize(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_psys_memsize
         integer(kind=c_int), value, intent(in) :: psys_num
      end function fdps_get_psys_memsize

      subroutine fdps_get_psys_time_prof(psys_num,prof) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_time_profile
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
         type(fdps_time_prof), intent(inout) :: prof
      end subroutine fdps_get_psys_time_prof

      subroutine fdps_clear_psys_time_prof(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
      end subroutine fdps_clear_psys_time_prof

      subroutine fdps_set_nptcl_smpl(psys_num,nptcl) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,nptcl
      end subroutine fdps_set_nptcl_smpl

      subroutine fdps_set_nptcl_loc(psys_num,nptcl) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,nptcl
      end subroutine fdps_set_nptcl_loc

      function fdps_get_nptcl_loc(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_nptcl_loc
         integer(kind=c_int), value, intent(in) :: psys_num
      end function fdps_get_nptcl_loc

      function fdps_get_nptcl_glb(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_nptcl_glb
         integer(kind=c_int), value, intent(in) :: psys_num
      end function fdps_get_nptcl_glb

      subroutine fdps_get_psys_cptr(psys_num,cptr) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
         type(c_ptr), intent(inout) :: cptr
      end subroutine fdps_get_psys_cptr

      subroutine fdps_exchange_particle(psys_num,dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,dinfo_num
      end subroutine fdps_exchange_particle

      !----------------------------------------------------------
      ! [Comment] A place where the interface statements for 
      !           fdps_add_particle*() are generated.
      ! fdps-autogen:add_particle:if;
      !----------------------------------------------------------

      subroutine fdps_remove_particle(psys_num,nptcl,ptcl_indx) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,nptcl
         integer(kind=c_int), dimension(*), intent(in) :: ptcl_indx
      end subroutine fdps_remove_particle

      subroutine fdps_adjust_pos_into_root_domain(psys_num,dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,dinfo_num
      end subroutine fdps_adjust_pos_into_root_domain

      !----------------------------------------------------------
      ! [Comment] A place where the interface statements for 
      !           fdps_sort_particle*() are generated.
      ! fdps-autogen:sort_particle:if;
      !----------------------------------------------------------

      !----------------------
      !  Domain Info
      !----------------------
      subroutine fdps_create_dinfo(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(inout) :: dinfo_num 
      end subroutine fdps_create_dinfo

      subroutine fdps_delete_dinfo(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num 
      end subroutine fdps_delete_dinfo

      subroutine fdps_init_dinfo(dinfo_num,coef_ema) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num 
         real(kind=c_float), value, intent(in) :: coef_ema
      end subroutine fdps_init_dinfo

      subroutine fdps_get_dinfo_time_prof(dinfo_num,prof) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_time_profile
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num
         type(fdps_time_prof), intent(inout) :: prof
      end subroutine fdps_get_dinfo_time_prof

      subroutine fdps_clear_dinfo_time_prof(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num
      end subroutine fdps_clear_dinfo_time_prof

      subroutine fdps_set_nums_domain(dinfo_num,nx,ny,nz) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num,nx,ny,nz
      end subroutine fdps_set_nums_domain

      subroutine fdps_set_boundary_condition(dinfo_num,bc) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num
         integer(kind=c_int), value, intent(in) :: bc
      end subroutine fdps_set_boundary_condition

      subroutine fdps_set_pos_root_domain(dinfo_num,low,high) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_vector
         implicit none
         integer(kind=c_int), value, intent(in):: dinfo_num
         type(fdps_f32vec), intent(in) :: low,high
      end subroutine fdps_set_pos_root_domain

      subroutine fdps_collect_sample_particle(dinfo_num,psys_num, &
                                              clear,weight) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num,psys_num
         logical(kind=c_bool), value, intent(in) :: clear
         real(kind=c_float), value, intent(in) :: weight
      end subroutine fdps_collect_sample_particle

      subroutine fdps_decompose_domain(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num
      end subroutine fdps_decompose_domain

      subroutine fdps_decompose_domain_all(dinfo_num,psys_num, &
                                           weight) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num,psys_num
         real(kind=c_float), value, intent(in) :: weight
      end subroutine fdps_decompose_domain_all

      !----------------------
      !  Tree 
      !----------------------
      subroutine fdps_create_tree(tree_num,tree_info) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(inout) :: tree_num
         character(kind=c_char), dimension(*), intent(in) :: tree_info
      end subroutine fdps_create_tree

      subroutine fdps_delete_tree(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
      end subroutine fdps_delete_tree

      subroutine fdps_init_tree(tree_num,nptcl,theta, &
                                n_leaf_limit,n_group_limit) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,nptcl
         real(kind=c_float), value, intent(in) :: theta
         integer(kind=c_int), value, intent(in) :: n_leaf_limit,n_group_limit
      end subroutine fdps_init_tree

      subroutine fdps_get_tree_info(tree_num,tree_info,charlen) bind(c) 
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
         character(kind=c_char), dimension(*), intent(inout) :: tree_info
         integer(kind=c_size_t), intent(inout) :: charlen
      end subroutine fdps_get_tree_info

      function fdps_get_tree_memsize(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_tree_memsize
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_tree_memsize

      subroutine fdps_get_tree_time_prof(tree_num,prof) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_time_profile
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
         type(fdps_time_prof), intent(inout) :: prof
      end subroutine fdps_get_tree_time_prof

      subroutine fdps_clear_tree_time_prof(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
      end subroutine fdps_clear_tree_time_prof

      function fdps_get_num_interact_ep_ep_loc(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_interact_ep_ep_loc
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_interact_ep_ep_loc

      function fdps_get_num_interact_ep_sp_loc(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_interact_ep_sp_loc
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_interact_ep_sp_loc

      function fdps_get_num_interact_ep_ep_glb(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_interact_ep_ep_glb
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_interact_ep_ep_glb

      function fdps_get_num_interact_ep_sp_glb(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_interact_ep_sp_glb
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_interact_ep_sp_glb

      subroutine fdps_clear_num_interact(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
      end subroutine fdps_clear_num_interact

      function fdps_get_num_tree_walk_loc(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_tree_walk_loc
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_tree_walk_loc

      function fdps_get_num_tree_walk_glb(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_tree_walk_glb
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_tree_walk_glb

      !----------------------------------------------------------
      ! [Comment] A place where the interface statements for 
      !           fdps_calc_force_*() are generated.
      ! fdps-autogen:calc_force_all_and_write_back:if;
      ! fdps-autogen:calc_force_all:if;
      ! fdps-autogen:calc_force_making_tree:if;
      ! fdps-autogen:calc_force_and_write_back:if;
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! [Comment] A place where the interface statements for
      !           fdps_get_neighbor_list*() are generated.
      ! fdps-autogen:get_neighbor_list:if;
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! [Comment] A place where the interface statements for
      !           fdps_get_epj_from_id*() are generated.
      ! fdps-autogen:get_epj_from_id:if;
      !----------------------------------------------------------

      !----------------------
      !  MPI comm. 
      !----------------------
      function fdps_get_rank() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_rank
      end function fdps_get_rank

      function fdps_get_rank_multi_dim(id) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_rank_multi_dim
         integer(kind=c_int), value, intent(in) :: id
      end function fdps_get_rank_multi_dim

      function fdps_get_num_procs() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_num_procs
      end function fdps_get_num_procs

      function fdps_get_num_procs_multi_dim(id) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_num_procs_multi_dim
         integer(kind=c_int), value, intent(in) :: id
      end function fdps_get_num_procs_multi_dim

      subroutine fdps_get_logical_and(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         logical(kind=c_bool), value, intent(in) :: f_in
         logical(kind=c_bool), intent(inout) :: f_out
      end subroutine fdps_get_logical_and

      subroutine fdps_get_logical_or(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         logical(kind=c_bool), value, intent(in) :: f_in
         logical(kind=c_bool), intent(inout) :: f_out
      end subroutine fdps_get_logical_or

      !-------------------------------------------
      ! [Comment] A place where the interface statements for 
      !           reduction functions are generated.
      ! fdps-autogen:get_min_value:ftn_if;
      ! fdps-autogen:get_max_value:ftn_if;
      ! fdps-autogen:get_sum:ftn_if;
      ! fdps-autogen:broadcast:ftn_if;
      !-------------------------------------------

      function fdps_get_wtime() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_get_wtime
      end function fdps_get_wtime

      !----------------------
      !  Utility
      !----------------------
      subroutine fdps_mt_init_genrand(s) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: s 
      end subroutine fdps_mt_init_genrand

      function fdps_mt_genrand_int31() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_mt_genrand_int31
      end function fdps_mt_genrand_int31

      function fdps_mt_genrand_real1() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_mt_genrand_real1
      end function fdps_mt_genrand_real1

      function fdps_mt_genrand_real2() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_mt_genrand_real2
      end function fdps_mt_genrand_real2

      function fdps_mt_genrand_real3() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_mt_genrand_real3
      end function fdps_mt_genrand_real3

      function fdps_mt_genrand_res53() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_mt_genrand_res53
      end function fdps_mt_genrand_res53

      !----------------------
      !  Particle Mesh
      !----------------------
#if PARTICLE_SIMULATOR_USE_PM_MODULE
      subroutine fdps_create_pm(pm_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(inout) :: pm_num
      end subroutine fdps_create_pm

      subroutine fdps_delete_pm(pm_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: pm_num
      end subroutine fdps_delete_pm

      function fdps_get_pm_mesh_num() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_pm_mesh_num
      end function fdps_get_pm_mesh_num

      function fdps_get_pm_cutoff_radius() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_get_pm_cutoff_radius
      end function fdps_get_pm_cutoff_radius

      subroutine fdps_set_dinfo_of_pm(pm_num,dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: pm_num,dinfo_num
      end subroutine fdps_set_dinfo_of_pm

      subroutine fdps_set_psys_of_pm(pm_num,psys_num,clear) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: pm_num,psys_num
         logical(kind=c_bool), value, intent(in) :: clear
      end subroutine fdps_set_psys_of_pm

      subroutine fdps_get_pm_force(pm_num,pos,f) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_vector
         implicit none
         integer(kind=c_int), value, intent(in) :: pm_num
         type(fdps_f32vec), intent(in) :: pos
         type(fdps_f32vec), intent(inout) :: f
      end subroutine fdps_get_pm_force

      subroutine fdps_get_pm_potential(pm_num,pos,pot) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_vector
         implicit none
         integer(kind=c_int), value, intent(in) :: pm_num
         type(fdps_f32vec), intent(in) :: pos
         real(kind=c_float), intent(inout) :: pot
      end subroutine fdps_get_pm_potential

      subroutine fdps_calc_pm_force_only(pm_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: pm_num
      end subroutine fdps_calc_pm_force_only

      subroutine fdps_calc_pm_force_all_and_write_back(pm_num,   &
                                                       psys_num, &
                                                       dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: pm_num,psys_num,dinfo_num
      end subroutine fdps_calc_pm_force_all_and_write_back
#endif

   end interface

   contains

   !----------------------------------------------------------
   subroutine PS_initialize(this)
      implicit none
      class(FDPS_controller) :: this

      call fdps_initialize()

   end subroutine PS_initialize 

   !----------------------------------------------------------
   subroutine PS_finalize(this)
      implicit none
      class(FDPS_controller) :: this

      call fdps_finalize()

   end subroutine PS_finalize 

   !----------------------------------------------------------
   subroutine PS_abort(this,err_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN), optional :: err_num

      if (present(err_num)) then
         call fdps_abort(err_num)
      else
         call fdps_abort(-1)
      end if

   end subroutine PS_abort

   !##########################################################
   subroutine create_psys(this,psys_num,psys_info_in)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(INOUT) :: psys_num
      character(len=*,kind=c_char), intent(IN) :: psys_info_in
     !character(len=:,kind=c_char), allocatable :: psys_info
      ! Note that Fortran compiler in K-computer does not support
      ! auto-reallocation introduced in Fortran 2003.
      ! (see http://www.nag-j.co.jp/fortran/fortran2003/Fortran2003_3_7.html)
      character(len=len(psys_info_in)+1,kind=c_char) :: psys_info

      psys_info = trim(psys_info_in) // c_null_char
      call fdps_create_psys(psys_num,psys_info)

   end subroutine create_psys

   !----------------------------------------------------------
   subroutine delete_psys(this,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      call fdps_delete_psys(psys_num)

   end subroutine delete_psys

   !----------------------------------------------------------
   subroutine init_psys(this,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      call fdps_init_psys(psys_num)

   end subroutine init_psys

   !----------------------------------------------------------
   subroutine get_psys_info(this,psys_num,psys_info)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num
      character(len=*,kind=c_char), intent(INOUT) :: psys_info
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: buf
      integer(kind=c_size_t) :: charlen

      call fdps_get_psys_info(psys_num,buf,charlen)
      psys_info = buf(1:charlen) ! copy
      ! [** Important **]
       !    You should use the intrinsic function trim() when
       !    you compare psys_info with an immediate value of
       !    a string of characters.
       !    cf.) 
       !    write(*,100)psys_info
       !    write(*,100)trim(psys_info)
       !    100 format(a,'|')

   end subroutine get_psys_info

   !----------------------------------------------------------
   function get_psys_memsize(this,psys_num)
      implicit none
      integer(kind=c_long_long) :: get_psys_memsize
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      get_psys_memsize =  fdps_get_psys_memsize(psys_num)
      
   end function get_psys_memsize

   !----------------------------------------------------------
   subroutine get_psys_time_prof(this,psys_num,prof)
      use fdps_time_profile
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num
      type(fdps_time_prof), intent(INOUT) :: prof

      call fdps_get_psys_time_prof(psys_num,prof)
      
   end subroutine get_psys_time_prof

   !----------------------------------------------------------
   subroutine clear_psys_time_prof(this,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      call fdps_clear_psys_time_prof(psys_num)
      
   end subroutine clear_psys_time_prof

   !----------------------------------------------------------
   subroutine set_nptcl_smpl(this,psys_num,nptcl)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num,nptcl

      call fdps_set_nptcl_smpl(psys_num,nptcl)

   end subroutine set_nptcl_smpl

   !----------------------------------------------------------
   subroutine set_nptcl_loc(this,psys_num,nptcl)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num,nptcl

      call fdps_set_nptcl_loc(psys_num,nptcl)

   end subroutine set_nptcl_loc

   !----------------------------------------------------------
   function get_nptcl_loc(this,psys_num)
      implicit none
      integer :: get_nptcl_loc
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      get_nptcl_loc =  fdps_get_nptcl_loc(psys_num)
      
   end function get_nptcl_loc

   !----------------------------------------------------------
   function get_nptcl_glb(this,psys_num)
      implicit none
      integer :: get_nptcl_glb
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      get_nptcl_glb =  fdps_get_nptcl_glb(psys_num)
      
   end function get_nptcl_glb

   !----------------------------------------------------------
   ! [Comment] A place where the definitions or implementations
   !           get_psys_fptr*() are generated.
   ! fdps-autogen:get_psys_fptr:impl;

   !----------------------------------------------------------
   subroutine exchange_particle(this,psys_num,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num,dinfo_num

      call fdps_exchange_particle(psys_num,dinfo_num)

   end subroutine exchange_particle

   !----------------------------------------------------------
   ! [Comment] A place where the definitions or implementations
   !           add_particle*() are generated.
   ! fdps-autogen:add_particle:impl;

   !----------------------------------------------------------
   subroutine remove_particle(this,psys_num,nptcl,ptcl_indx)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num
      integer(kind=c_int), intent(IN) :: nptcl
      integer(kind=c_int), dimension(nptcl), intent(IN) :: ptcl_indx

      call fdps_remove_particle(psys_num,nptcl,ptcl_indx)

   end subroutine remove_particle

   !----------------------------------------------------------
   subroutine adjust_pos_into_root_domain(this,psys_num,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num,dinfo_num

      call fdps_adjust_pos_into_root_domain(psys_num,dinfo_num)

   end subroutine adjust_pos_into_root_domain

   !----------------------------------------------------------
   subroutine sort_particle(this,psys_num,pfunc_comp)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num
      type(c_funptr), intent(in) :: pfunc_comp
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: psys_info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name

      call get_psys_info(this,psys_num,psys_info)

      select case (trim(psys_info)) 
      !--------------
      ! fdps-autogen:sort_particle:impl;
      !--------------
      case default
         errmsg = "Unknow psys_num is specified"
         func_name = "sort_particle"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine sort_particle

   !##########################################################
   subroutine create_dinfo(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(INOUT) :: dinfo_num

      call fdps_create_dinfo(dinfo_num)

   end subroutine create_dinfo

   !----------------------------------------------------------
   subroutine delete_dinfo(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num

      call fdps_delete_dinfo(dinfo_num)

   end subroutine delete_dinfo

   !----------------------------------------------------------
   subroutine init_dinfo(this,dinfo_num,coef_ema)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      real(kind=c_float), intent(IN), optional :: coef_ema

      if (present(coef_ema)) then
         call fdps_init_dinfo(dinfo_num,coef_ema)
      else
         call fdps_init_dinfo(dinfo_num,1.0_c_float)
      end if

   end subroutine init_dinfo

   !----------------------------------------------------------
   subroutine get_dinfo_time_prof(this,dinfo_num,prof)
      use fdps_time_profile
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      type(fdps_time_prof), intent(INOUT) :: prof

      call fdps_get_dinfo_time_prof(dinfo_num,prof)
      
   end subroutine get_dinfo_time_prof

   !----------------------------------------------------------
   subroutine clear_dinfo_time_prof(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num

      call fdps_clear_dinfo_time_prof(dinfo_num)
      
   end subroutine clear_dinfo_time_prof

   !----------------------------------------------------------
   subroutine set_nums_domain(this,dinfo_num,nx,ny,nz)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num,nx,ny
      integer(kind=c_int), intent(IN), optional :: nz

      if (present(nz)) then
         call fdps_set_nums_domain(dinfo_num,nx,ny,nz)
      else
         call fdps_set_nums_domain(dinfo_num,nx,ny,1_c_int)
      end if

   end subroutine set_nums_domain

   !----------------------------------------------------------
   subroutine set_boundary_condition(this,dinfo_num,bc)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num,bc
      
      call fdps_set_boundary_condition(dinfo_num,bc)

   end subroutine set_boundary_condition

   !----------------------------------------------------------
   subroutine set_pos_root_domain_a32(this,dinfo_num,low,high)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      real(kind=c_float), dimension(space_dim), intent(in) :: low,high
      !* Local variables
      type(fdps_f32vec) :: low_,high_
      
      low_ = low; high_ = high
      call fdps_set_pos_root_domain(dinfo_num,low_,high_)

   end subroutine set_pos_root_domain_a32

   !----------------------------------------------------------
   subroutine set_pos_root_domain_a64(this,dinfo_num,low,high)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      real(kind=c_double), dimension(space_dim), intent(in) :: low,high
      !* Local variables
      type(fdps_f32vec) :: low_,high_
      
      low_ = low; high_ = high
      call fdps_set_pos_root_domain(dinfo_num,low_,high_)

   end subroutine set_pos_root_domain_a64

   !----------------------------------------------------------
   subroutine set_pos_root_domain_v32(this,dinfo_num,low,high)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      type(fdps_f32vec), intent(in) :: low,high 

      call fdps_set_pos_root_domain(dinfo_num,low,high)

   end subroutine set_pos_root_domain_v32

   !----------------------------------------------------------
   subroutine set_pos_root_domain_v64(this,dinfo_num,low,high)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      type(fdps_f64vec), intent(in) :: low,high 
      !* Local variables
      type(fdps_f32vec) :: low_,high_
      
      low_ = low; high_ = high
      call fdps_set_pos_root_domain(dinfo_num,low_,high_)

   end subroutine set_pos_root_domain_v64

   !----------------------------------------------------------
   subroutine collect_sample_particle(this,dinfo_num,psys_num, &
                                      clear,weight)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num,psys_num
      logical(kind=c_bool), intent(IN), optional :: clear
      real(kind=c_float), intent(IN), optional :: weight
      !* Local variables
      logical(kind=c_bool) :: clear_
      real(kind=c_float) :: weight_

      if (present(clear)) then
         clear_ = clear
      else
         clear_ = .true. ! default value
      end if
      if (present(weight)) then
         weight_ = weight
      else
         weight_ = 1.0_c_float ! default value
      end if
      call fdps_collect_sample_particle(dinfo_num,psys_num, &
                                        clear_,weight_)

   end subroutine collect_sample_particle

   !----------------------------------------------------------
   subroutine decompose_domain(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num

      call fdps_decompose_domain(dinfo_num)

   end subroutine decompose_domain

   !----------------------------------------------------------
   subroutine decompose_domain_all(this,dinfo_num,psys_num, &
                                   weight)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num,psys_num
      real(kind=c_float), intent(IN), optional :: weight
      !* Local variables
      real(kind=c_float) :: wgh

      if (present(weight)) then
         call fdps_decompose_domain_all(dinfo_num,psys_num,weight)
      else
         wgh = real(fdps_get_nptcl_loc(psys_num),kind=c_float)
         call fdps_decompose_domain_all(dinfo_num,psys_num,wgh)
      end if

   end subroutine decompose_domain_all

   !##########################################################
   subroutine create_tree(this,tree_num,tree_info_in)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(INOUT) :: tree_num
      character(len=*,kind=c_char), intent(IN) :: tree_info_in
     !character(len=:,kind=c_char), allocatable :: tree_info
      ! See the comment in API create_psys for this comment-out.
      character(len=len(tree_info_in)+1,kind=c_char) :: tree_info

      tree_info = trim(tree_info_in) // c_null_char
      call fdps_create_tree(tree_num,tree_info)

   end subroutine create_tree

   !----------------------------------------------------------
   subroutine delete_tree(this,tree_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      call fdps_delete_tree(tree_num)

   end subroutine delete_tree

   !----------------------------------------------------------
   subroutine init_tree(this,tree_num,nptcl, &
                        theta,n_leaf_limit,n_group_limit)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,nptcl
      real(kind=c_float), intent(IN), optional :: theta
      integer(kind=c_int), intent(IN), optional :: n_leaf_limit,n_group_limit
      !* Local variables
      real(kind=c_float) :: theta_
      integer(kind=c_int) :: n_leaf_limit_,n_group_limit_

      if (present(theta)) then
         theta_ = theta
      else
         theta_ = 0.7_c_float
      end if
      if (present(n_leaf_limit)) then
         n_leaf_limit_ = n_leaf_limit
      else
         n_leaf_limit_ = 8_c_int
      end if
      if (present(n_group_limit)) then
         n_group_limit_ = n_group_limit
      else
         n_group_limit_ = 64_c_int
      end if

      call fdps_init_tree(tree_num,      &
                          nptcl,         &
                          theta_,        &
                          n_leaf_limit_, &
                          n_group_limit_)

   end subroutine init_tree

   !----------------------------------------------------------
   subroutine get_tree_info(this,tree_num,tree_info)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num
      character(len=*,kind=c_char), intent(INOUT) :: tree_info
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: buf
      integer(kind=c_size_t) :: charlen

      call fdps_get_tree_info(tree_num,buf,charlen)
      tree_info = buf(1:charlen) ! copy
      ! [** Important **]
      !    You should use the intrinsic function trim() when
      !    you compare tree_info with an immediate value of
      !    a string of characters.
      !    cf.) 
      !    write(*,100)tree_info
      !    write(*,100)trim(tree_info)
      !    100 format(a,'|')

   end subroutine get_tree_info

   !----------------------------------------------------------
   function get_tree_memsize(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_tree_memsize
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_tree_memsize =  fdps_get_tree_memsize(tree_num)
      
   end function get_tree_memsize

   !----------------------------------------------------------
   subroutine get_tree_time_prof(this,tree_num,prof)
      use fdps_time_profile
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num
      type(fdps_time_prof), intent(INOUT) :: prof

      call fdps_get_tree_time_prof(tree_num,prof)
      
   end subroutine get_tree_time_prof

   !----------------------------------------------------------
   subroutine clear_tree_time_prof(this,tree_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      call fdps_clear_tree_time_prof(tree_num)
      
   end subroutine clear_tree_time_prof

   !----------------------------------------------------------
   function get_num_interact_ep_ep_loc(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_interact_ep_ep_loc
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_interact_ep_ep_loc = fdps_get_num_interact_ep_ep_loc(tree_num)
      
   end function get_num_interact_ep_ep_loc

   !----------------------------------------------------------
   function get_num_interact_ep_sp_loc(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_interact_ep_sp_loc
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_interact_ep_sp_loc = fdps_get_num_interact_ep_sp_loc(tree_num)
      
   end function get_num_interact_ep_sp_loc

   !----------------------------------------------------------
   function get_num_interact_ep_ep_glb(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_interact_ep_ep_glb
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_interact_ep_ep_glb = fdps_get_num_interact_ep_ep_glb(tree_num)
      
   end function get_num_interact_ep_ep_glb

   !----------------------------------------------------------
   function get_num_interact_ep_sp_glb(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_interact_ep_sp_glb
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_interact_ep_sp_glb = fdps_get_num_interact_ep_sp_glb(tree_num)
      
   end function get_num_interact_ep_sp_glb

   !----------------------------------------------------------
   subroutine clear_num_interact(this,tree_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      call fdps_clear_num_interact(tree_num)

   end subroutine clear_num_interact

   !----------------------------------------------------------
   function get_num_tree_walk_loc(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_tree_walk_loc
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_tree_walk_loc = fdps_get_num_tree_walk_loc(tree_num)
      
   end function get_num_tree_walk_loc

   !----------------------------------------------------------
   function get_num_tree_walk_glb(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_tree_walk_glb
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_tree_walk_glb = fdps_get_num_tree_walk_glb(tree_num)
      
   end function get_num_tree_walk_glb

   !----------------------------------------------------------
   subroutine calc_force_all_and_write_back_s(this,        &
                                              tree_num,    &
                                              pfunc_ep_ep, &
                                              psys_num,    &
                                              dinfo_num,   &
                                              list_mode)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      integer(kind=c_int), optional, intent(IN) :: list_mode
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      logical(kind=c_bool) :: clear_
      integer(kind=c_int) :: list_mode_
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name

      !* Process Optional arguments
      clear_ = .true.
      if (present(list_mode)) then
         list_mode_ = list_mode
      else
         list_mode_ = fdps_make_list ! default value
      end if

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info)
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_all_and_write_back:impl:short;
      !--------------
      case default
         errmsg = "The combination psys_num and tree_num is invalid"
         func_name = "calc_force_all_and_write_back_s"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine calc_force_all_and_write_back_s

   !----------------------------------------------------------
   subroutine calc_force_all_and_write_back_l(this,        &
                                              tree_num,    &
                                              pfunc_ep_ep, &
                                              pfunc_ep_sp, &
                                              psys_num,    &
                                              dinfo_num,   &
                                              list_mode)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      integer(kind=c_int), optional, intent(IN) :: list_mode
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      logical(kind=c_bool) :: clear_
      integer(kind=c_int) :: list_mode_
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name

      !* Process optional arguments
      clear_ = .true.
      if (present(list_mode)) then
         list_mode_ = list_mode
      else
         list_mode_ = fdps_make_list ! default value
      end if

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info) 
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_all_and_write_back:impl:long;
      !--------------
      case default
         errmsg = "The combination psys_num and tree_num is invalid"
         func_name = "calc_force_all_and_write_back_l"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine calc_force_all_and_write_back_l

   !----------------------------------------------------------
   subroutine calc_force_all_s(this,        &
                               tree_num,    &
                               pfunc_ep_ep, &
                               psys_num,    &
                               dinfo_num,   &
                               list_mode)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      integer(kind=c_int), optional, intent(IN) :: list_mode
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      logical(kind=c_bool) :: clear_
      integer(kind=c_int) :: list_mode_
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name

      !* Process optional arguments
      clear_ = .true.
      if (present(list_mode)) then
         list_mode_ = list_mode
      else
         list_mode_ = fdps_make_list ! default value
      end if

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info)
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_all:impl:short;
      !--------------
      case default
         errmsg = "The combination psys_num and tree_num is invalid"
         func_name = "calc_force_all_s"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine calc_force_all_s

   !----------------------------------------------------------
   subroutine calc_force_all_l(this,        &
                               tree_num,    &
                               pfunc_ep_ep, &
                               pfunc_ep_sp, &
                               psys_num,    &
                               dinfo_num,   &
                               list_mode)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      integer(kind=c_int), optional, intent(IN) :: list_mode
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      logical(kind=c_bool) :: clear_
      integer(kind=c_int) :: list_mode_
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name

      !* Process optional arguments
      clear_ = .true.
      if (present(list_mode)) then
         list_mode_ = list_mode
      else
         list_mode_ = fdps_make_list ! default value
      end if

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info) 
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_all:impl:long;
      !--------------
      case default
         errmsg = "The combination psys_num and tree_num is invalid"
         func_name = "calc_force_all_l"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine calc_force_all_l

   !----------------------------------------------------------
   subroutine calc_force_making_tree_s(this,        &
                                       tree_num,    &
                                       pfunc_ep_ep, &
                                       dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      logical(kind=c_bool) :: clear_
      character(len=bufsize,kind=c_char) :: info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name
      
      !* Process optioanl arguments
      clear_ = .true.

      call get_tree_info(this,tree_num,info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_making_tree:impl:short;
      !--------------
      case default
         errmsg = "tree_num passed is invalid"
         func_name = "calc_force_making_tree_s"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine calc_force_making_tree_s

   !----------------------------------------------------------
   subroutine calc_force_making_tree_l(this,        &
                                       tree_num,    &
                                       pfunc_ep_ep, &
                                       pfunc_ep_sp, &
                                       dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      logical(kind=c_bool) :: clear_
      character(len=bufsize,kind=c_char) :: info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name

      !* Process optional arguments
      clear_ = .true.

      call get_tree_info(this,tree_num,info) 

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_making_tree:impl:long;
      !--------------
      case default
         errmsg = "tree_num passed is invalid"
         func_name = "calc_force_making_tree_l"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine calc_force_making_tree_l

   !----------------------------------------------------------
   subroutine calc_force_and_write_back_s(this,        &
                                          tree_num,    &
                                          pfunc_ep_ep, &
                                          psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      logical(kind=c_bool) :: clear_
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name

      !* Process optional arguments
      clear_ = .true.

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info)
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_and_write_back:impl:short;
      !--------------
      case default
         errmsg = "The combination psys_num and tree_num is invalid"
         func_name = "calc_force_and_write_back_s"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine calc_force_and_write_back_s

   !----------------------------------------------------------
   subroutine calc_force_and_write_back_l(this,        &
                                          tree_num,    &
                                          pfunc_ep_ep, &
                                          pfunc_ep_sp, &
                                          psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      logical(kind=c_bool) :: clear_
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info
      !-(To throw errors)
      character(len=64) :: errmsg,func_name

      !* Process optional arguments
      clear_ = .true.

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info) 
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_and_write_back:impl:long;
      !--------------
      case default
         errmsg = "The combination psys_num and tree_num is invalid"
         func_name = "calc_force_and_write_back_l"
         call print_errmsg(errmsg,func_name)
         call PS_abort(this)
         stop 1
      end select

   end subroutine calc_force_and_write_back_l

   !----------------------------------------------------------
   ! [Comment] A place where the implementations of
   !           get_neighbor_list*() are generated.
   ! fdps-autogen:get_neighbor_list:impl;

   !----------------------------------------------------------
   ! [Comment] A place where the implementations of
   !           get_epj_from_id*() are generated.
   ! fdps-autogen:get_epj_from_id:impl;

   !##########################################################
   function get_rank(this)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: get_rank

      get_rank = fdps_get_rank()
      
   end function get_rank

   !----------------------------------------------------------
   function get_rank_multi_dim(this,id)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: get_rank_multi_dim
      integer(kind=c_int), intent(IN) :: id

      get_rank_multi_dim = fdps_get_rank_multi_dim(id)
      
   end function get_rank_multi_dim

   !----------------------------------------------------------
   function get_num_procs(this)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: get_num_procs
      
      get_num_procs = fdps_get_num_procs()

   end function get_num_procs

   !----------------------------------------------------------
   function get_num_procs_multi_dim(this,id)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: get_num_procs_multi_dim
      integer(kind=c_int), intent(IN) :: id
      
      get_num_procs_multi_dim = fdps_get_num_procs_multi_dim(id)

   end function get_num_procs_multi_dim

   !----------------------------------------------------------
   subroutine get_logical_and(this,f_in,f_out)
      implicit none
      class(FDPS_controller) :: this
      logical(kind=c_bool), intent(IN) :: f_in
      logical(kind=c_bool), intent(INOUT) :: f_out

      call fdps_get_logical_and(f_in,f_out)

   end subroutine get_logical_and

   !----------------------------------------------------------
   subroutine get_logical_or(this,f_in,f_out)
      implicit none
      class(FDPS_controller) :: this
      logical(kind=c_bool), intent(IN) :: f_in
      logical(kind=c_bool), intent(INOUT) :: f_out

      call fdps_get_logical_or(f_in,f_out)

   end subroutine get_logical_or

   !----------------------------------------------------------
   ! [Comment] A place where the implementations of
   !           reduction functions are generated.
   ! fdps-autogen:get_min_value:ftn_impl;
   ! fdps-autogen:get_max_value:ftn_impl;
   ! fdps-autogen:get_sum:ftn_impl;
   ! fdps-autogen:broadcast:ftn_impl;
   !----------------------------------------------------------

   !----------------------------------------------------------
   function get_wtime(this)
      implicit none
      class(FDPS_controller) :: this 
      real(kind=c_double) :: get_wtime

      get_wtime = fdps_get_wtime()

   end function get_wtime

   !##########################################################
   subroutine MT_init_genrand(this,s)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: s

      call fdps_mt_init_genrand(s)
      
   end subroutine MT_init_genrand

   !----------------------------------------------------------
   function MT_genrand_int31(this)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: MT_genrand_int31

      MT_genrand_int31 = fdps_mt_genrand_int31()

   end function MT_genrand_int31

   !----------------------------------------------------------
   function MT_genrand_real1(this)
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double) :: MT_genrand_real1

      MT_genrand_real1 = fdps_mt_genrand_real1()

   end function MT_genrand_real1

   !----------------------------------------------------------
   function MT_genrand_real2(this)
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double) :: MT_genrand_real2

      MT_genrand_real2 = fdps_mt_genrand_real2()

   end function MT_genrand_real2

   !----------------------------------------------------------
   function MT_genrand_real3(this)
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double) :: MT_genrand_real3

      MT_genrand_real3 = fdps_mt_genrand_real3()

   end function MT_genrand_real3

   !----------------------------------------------------------
   function MT_genrand_res53(this)
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double) :: MT_genrand_res53

      MT_genrand_res53 = fdps_mt_genrand_res53()

   end function MT_genrand_res53

   !##########################################################
#if PARTICLE_SIMULATOR_USE_PM_MODULE
   subroutine create_pm(this,pm_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(INOUT) :: pm_num

      call fdps_create_pm(pm_num)

   end subroutine create_pm

   !----------------------------------------------------------
   subroutine delete_pm(this,pm_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: pm_num

      call fdps_delete_pm(pm_num)

   end subroutine delete_pm

   !----------------------------------------------------------
   function get_pm_mesh_num(this)
      implicit none
      integer(kind=c_int) :: get_pm_mesh_num
      class(FDPS_controller) :: this

      get_pm_mesh_num = fdps_get_pm_mesh_num()
   
   end function get_pm_mesh_num

   !----------------------------------------------------------
   function get_pm_cutoff_radius(this)
      implicit none
      real(kind=c_double) :: get_pm_cutoff_radius
      class(FDPS_controller) :: this

      get_pm_cutoff_radius = fdps_get_pm_cutoff_radius()

   end function get_pm_cutoff_radius

   !----------------------------------------------------------
   subroutine set_dinfo_of_pm(this,pm_num,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num,dinfo_num

      call fdps_set_dinfo_of_pm(pm_num,dinfo_num)

   end subroutine set_dinfo_of_pm

   !----------------------------------------------------------
   subroutine set_psys_of_pm(this,pm_num,psys_num,clear)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num,psys_num
      logical(kind=c_bool), intent(in), optional :: clear
      !* Local variables
      logical(kind=c_bool) :: clear_

      if (present(clear)) then
         clear_ = clear
      else
         clear_ = .true. ! default value
      end if
      call fdps_set_psys_of_pm(pm_num,psys_num,clear_)

   end subroutine set_psys_of_pm

   !----------------------------------------------------------
   subroutine get_pm_force_a32(this,pm_num,pos,f)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num
      real(kind=c_float), dimension(space_dim), intent(in) :: pos
      real(kind=c_float), dimension(space_dim), intent(inout) :: f
      !* Local variables
      integer :: i
      type(fdps_f32vec) :: pos_,f_

      pos_ = pos
      call fdps_get_pm_force(pm_num,pos_,f_)
      f(1) = f_%x
      f(2) = f_%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      f(3) = f_%z
#endif

   end subroutine get_pm_force_a32

   !----------------------------------------------------------
   subroutine get_pm_force_a64(this,pm_num,pos,f)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num
      real(kind=c_double), dimension(space_dim), intent(in) :: pos
      real(kind=c_double), dimension(space_dim), intent(inout) :: f
      !* Local variables
      type(fdps_f32vec) :: pos_,f_

      pos_ = pos
      call fdps_get_pm_force(pm_num,pos_,f_)
      f(1) = f_%x
      f(2) = f_%y
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
      f(3) = f_%z
#endif

   end subroutine get_pm_force_a64

   !----------------------------------------------------------
   subroutine get_pm_force_v32(this,pm_num,pos,f)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num
      type(fdps_f32vec), intent(in) :: pos
      type(fdps_f32vec), intent(inout) :: f

      call fdps_get_pm_force(pm_num,pos,f)

   end subroutine get_pm_force_v32

   !----------------------------------------------------------
   subroutine get_pm_force_v64(this,pm_num,pos,f)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num
      type(fdps_f64vec), intent(in) :: pos
      type(fdps_f64vec), intent(inout) :: f
      !* Local variables
      type(fdps_f32vec) :: pos_,f_

      pos_ = pos
      call fdps_get_pm_force(pm_num,pos_,f_)
      f = f_

   end subroutine get_pm_force_v64

   !----------------------------------------------------------
   subroutine get_pm_potential_a32(this,pm_num,pos,pot)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num
      real(kind=c_float), dimension(space_dim), intent(in) :: pos
      real(kind=c_float), intent(inout) :: pot
      !* Local variables
      type(fdps_f32vec) :: pos_

      pos_ = pos
      call fdps_get_pm_potential(pm_num,pos_,pot)

   end subroutine get_pm_potential_a32

   !----------------------------------------------------------
   subroutine get_pm_potential_a64(this,pm_num,pos,pot)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num
      real(kind=c_double), dimension(space_dim), intent(in) :: pos
      real(kind=c_float), intent(inout) :: pot
      !* Local variables
      type(fdps_f32vec) :: pos_

      pos_ = pos
      call fdps_get_pm_potential(pm_num,pos_,pot)

   end subroutine get_pm_potential_a64

   !----------------------------------------------------------
   subroutine get_pm_potential_v32(this,pm_num,pos,pot)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num
      type(fdps_f32vec), intent(in) :: pos
      real(kind=c_float), intent(inout) :: pot

      call fdps_get_pm_potential(pm_num,pos,pot)

   end subroutine get_pm_potential_v32

   !----------------------------------------------------------
   subroutine get_pm_potential_v64(this,pm_num,pos,pot)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num
      type(fdps_f64vec), intent(in) :: pos
      real(kind=c_float), intent(inout) :: pot
      !* Local variables
      type(fdps_f32vec) :: pos_

      pos_ = pos
      call fdps_get_pm_potential(pm_num,pos_,pot)

   end subroutine get_pm_potential_v64

   !----------------------------------------------------------
   subroutine calc_pm_force_only(this,pm_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num

      call fdps_calc_pm_force_only(pm_num)

   end subroutine calc_pm_force_only

   !----------------------------------------------------------
   subroutine calc_pm_force_all_and_write_back(this,pm_num, &
                                               psys_num,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: pm_num,psys_num,dinfo_num

      call fdps_calc_pm_force_all_and_write_back(pm_num,   &
                                                 psys_num, &
                                                 dinfo_num)
      
   end subroutine calc_pm_force_all_and_write_back
#endif

   !##########################################################
   subroutine print_errmsg(msg,func_name)
      implicit none
      character(len=*), intent(in) :: msg,func_name
      
      write(*,100)trim(msg),trim(func_name)
      100 format("*** PS_FTN_IF_ERROR ***"/ &
                 "message:  ",a/            &
                 "function: ",a/            &
                 "file    : FDPS_module.F90"/)
      flush(6)

   end subroutine print_errmsg

end module FDPS_module
