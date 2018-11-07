!================================
!   MODULE: FDPS time profile
!================================
module fdps_time_profile
   use, intrinsic :: iso_c_binding
   implicit none

   !**** PS::TimeProfile
   type, public, bind(c) :: fdps_time_prof
      real(kind=c_double) :: collect_sample_particle
      real(kind=c_double) :: decompose_domain
      real(kind=c_double) :: exchange_particle
      real(kind=c_double) :: set_particle_local_tree
      real(kind=c_double) :: set_particle_global_tree
      real(kind=c_double) :: make_local_tree
      real(kind=c_double) :: make_global_tree
      real(kind=c_double) :: set_root_cell
      real(kind=c_double) :: calc_force
      real(kind=c_double) :: calc_moment_local_tree
      real(kind=c_double) :: calc_moment_global_tree
      real(kind=c_double) :: make_LET_1st
      real(kind=c_double) :: make_LET_2nd
      real(kind=c_double) :: exchange_LET_1st
      real(kind=c_double) :: exchange_LET_2nd
      real(kind=c_double) :: write_back

      real(kind=c_double) :: morton_sort_local_tree
      real(kind=c_double) :: link_cell_local_tree
      real(kind=c_double) :: morton_sort_global_tree
      real(kind=c_double) :: link_cell_global_tree

      real(kind=c_double) :: make_local_tree_tot
      ! = make_local_tree + calc_moment_local_tree
      real(kind=c_double) :: make_global_tree_tot
      real(kind=c_double) :: exchange_LET_tot
      ! = make_LET_1st + make_LET_2nd + exchange_LET_1st + exchange_LET_2nd

      real(kind=c_double) :: calc_force__core__walk_tree
      real(kind=c_double) :: calc_force__core__keep_list
      real(kind=c_double) :: calc_force__core__copy_ep
      real(kind=c_double) :: calc_force__core__dispatch
      real(kind=c_double) :: calc_force__core__retrieve

      real(kind=c_double) :: calc_force__make_ipgroup
      real(kind=c_double) :: calc_force__core
      real(kind=c_double) :: calc_force__copy_original_order

      real(kind=c_double) :: exchange_particle__find_particle
      real(kind=c_double) :: exchange_particle__exchange_particle

      real(kind=c_double) :: decompose_domain__sort_particle_1st
      real(kind=c_double) :: decompose_domain__sort_particle_2nd
      real(kind=c_double) :: decompose_domain__sort_particle_3rd
      real(kind=c_double) :: decompose_domain__gather_particle

      real(kind=c_double) :: decompose_domain__setup
      real(kind=c_double) :: decompose_domain__determine_coord_1st
      real(kind=c_double) :: decompose_domain__migrae_particle_1st
      real(kind=c_double) :: decompose_domain__determine_coord_2nd
      real(kind=c_double) :: decompose_domain__determine_coord_3rd
      real(kind=c_double) :: decompose_domain__exchange_pos_domain

      real(kind=c_double) :: exchange_LET_1st__a2a_n
      real(kind=c_double) :: exchange_LET_1st__icomm_sp
      real(kind=c_double) :: exchange_LET_1st__a2a_sp
      real(kind=c_double) :: exchange_LET_1st__icomm_ep
      real(kind=c_double) :: exchange_LET_1st__a2a_ep

      real(kind=c_double) :: add_moment_as_sp_local
      real(kind=c_double) :: add_moment_as_sp_global
   end type fdps_time_prof

end module fdps_time_profile

