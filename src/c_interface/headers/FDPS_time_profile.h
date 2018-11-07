
//**** PS::TimeProfile
typedef struct {
   double collect_sample_particle;
   double decompose_domain;
   double exchange_particle;
   double set_particle_local_tree;
   double set_particle_global_tree;
   double make_local_tree;
   double make_global_tree;
   double set_root_cell;
   double calc_force;
   double calc_moment_local_tree;
   double calc_moment_global_tree;
   double make_LET_1st;
   double make_LET_2nd;
   double exchange_LET_1st;
   double exchange_LET_2nd;
   double write_back;

   double morton_sort_local_tree;
   double link_cell_local_tree;
   double morton_sort_global_tree;
   double link_cell_global_tree;

   double make_local_tree_tot; // = make_local_tree + calc_moment_local_tree
   double make_global_tree_tot;
   double exchange_LET_tot; // = make_LET_1st + make_LET_2nd + exchange_LET_1st + exchange_LET_2nd

   double calc_force__core__walk_tree;
   double calc_force__core__keep_list;
   double calc_force__core__copy_ep;
   double calc_force__core__dispatch;
   double calc_force__core__retrieve;

   double calc_force__make_ipgroup;
   double calc_force__core;
   double calc_force__copy_original_order;

   double exchange_particle__find_particle;
   double exchange_particle__exchange_particle;

   double decompose_domain__sort_particle_1st;
   double decompose_domain__sort_particle_2nd;
   double decompose_domain__sort_particle_3rd;
   double decompose_domain__gather_particle;

   double decompose_domain__setup;
   double decompose_domain__determine_coord_1st;
   double decompose_domain__migrae_particle_1st;
   double decompose_domain__determine_coord_2nd;
   double decompose_domain__determine_coord_3rd;
   double decompose_domain__exchange_pos_domain;

   double exchange_LET_1st__a2a_n;
   double exchange_LET_1st__icomm_sp;
   double exchange_LET_1st__a2a_sp;
   double exchange_LET_1st__icomm_ep;
   double exchange_LET_1st__a2a_ep;

   double add_moment_as_sp_local;
   double add_moment_as_sp_global;
} fdps_time_prof;
