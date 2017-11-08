/* Standard headers */
#include <stdio.h>
#include <iostream>
#include <cmath>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "FDPS_Manipulators.h"
using namespace FDPS_Manipulators;

extern "C" {

//----------------------
//  Basic
//----------------------
void fdps_initialize() {
   PS_Initialize();
}
void fdps_finalize() {
   PS_Finalize();
}
void fdps_abort(const int err_num) {
   PS_Abort(err_num);
}

//----------------------
//  Particle System 
//----------------------
void fdps_create_psys(int *psys_num,
                      char *psys_info) {
   create_psys(psys_num,psys_info);
}
void fdps_delete_psys(const int psys_num) {
   delete_psys(psys_num);
}
void fdps_init_psys(const int psys_num) {
   init_psys(psys_num);
}
void fdps_get_psys_info(const int psys_num,
                        char *psys_info,
                        size_t *charlen) {
   get_psys_info(psys_num,psys_info,charlen);
}
long long int fdps_get_psys_memsize(const int psys_num) {
   return get_psys_memsize(psys_num);
}
void fdps_get_psys_time_prof(const int psys_num,
                             PS::TimeProfile *prof) {
   get_psys_time_prof(psys_num,prof);
}
void fdps_clear_psys_time_prof(const int psys_num) {
   clear_psys_time_prof(psys_num);
}
void fdps_set_nptcl_smpl(const int psys_num,
                         const int nptcl) {
   set_nptcl_smpl(psys_num,nptcl);
}
void fdps_set_nptcl_loc(const int psys_num,
                        const int nptcl) {
   set_nptcl_loc(psys_num,nptcl);
}
int fdps_get_nptcl_loc(const int psys_num) {
   return get_nptcl_loc(psys_num);
}
int fdps_get_nptcl_glb(const int psys_num) {
   return get_nptcl_glb(psys_num);
}
void fdps_get_psys_cptr(const int psys_num,
                        void **cptr) {
   get_psys_cptr(psys_num,cptr);
}
void fdps_exchange_particle(const int psys_num,
                            const int dinfo_num) {
   exchange_particle(psys_num,dinfo_num);
}
//------------------------------------------------------------------
// [Comment] A place where fdps_add_particle_*() are generated.
// fdps-autogen:add_particle;
//------------------------------------------------------------------
void fdps_remove_particle(const int psys_num, 
                          const int nptcl,
                          int *ptcl_indx) {
   remove_particle(psys_num,nptcl,ptcl_indx);
}
void fdps_adjust_pos_into_root_domain(const int psys_num,
                                      const int dinfo_num) {
   adjust_pos_into_root_domain(psys_num,dinfo_num);
}
//------------------------------------------------------------------
// [Comment] A place where fdps_sort_particle_*() are generated.
// fdps-autogen:sort_particle;
//------------------------------------------------------------------

//----------------------
//  Domain Info
//----------------------
void fdps_create_dinfo(int *dinfo_num) {
   create_dinfo(dinfo_num);
}
void fdps_delete_dinfo(const int dinfo_num) {
   delete_dinfo(dinfo_num);
}
void fdps_init_dinfo(const int dinfo_num,
                     const float coef_ema) {
   init_dinfo(dinfo_num,coef_ema);
}
void fdps_get_dinfo_time_prof(const int dinfo_num, 
                              PS::TimeProfile *prof) {
   get_dinfo_time_prof(dinfo_num,prof);
}
void fdps_clear_dinfo_time_prof(const int dinfo_num) {
   clear_dinfo_time_prof(dinfo_num);
}
void fdps_set_nums_domain(const int dinfo_num,
                          const int nx,
                          const int ny,
                          const int nz) {
   set_nums_domain(dinfo_num,nx,ny,nz);
}
void fdps_set_boundary_condition(const int dinfo_num,
                                 const enum PS_BOUNDARY_CONDITION bc) {
   set_boundary_condition(dinfo_num,bc);
}
void fdps_set_pos_root_domain(const int dinfo_num,
                              const PS::F32vec *low,
                              const PS::F32vec *high) {
   set_pos_root_domain(dinfo_num,low,high);
}
void fdps_collect_sample_particle(const int dinfo_num,
                                  const int psys_num, 
                                  const bool clear,
                                  const float weight) {
   collect_sample_particle(dinfo_num,psys_num,clear,weight);
}
void fdps_decompose_domain(const int dinfo_num) {
   decompose_domain(dinfo_num);
}
void fdps_decompose_domain_all(const int dinfo_num,
                               const int psys_num, 
                               const float weight) {
   decompose_domain_all(dinfo_num,psys_num,weight);
}

//----------------------
//  Tree 
//----------------------
void fdps_create_tree(int *tree_num,
                      char *tree_info) {
   create_tree(tree_num,tree_info);
}
void fdps_delete_tree(const int tree_num) {
   delete_tree(tree_num);
}
void fdps_init_tree(const int tree_num,
                    const int nptcl,
                    const float theta,
                    const int n_leaf_limit,
                    const int n_group_limit) {
   init_tree(tree_num,nptcl,theta,
             n_leaf_limit,n_group_limit);
}
void fdps_get_tree_info(const int tree_num,
                        char *tree_info,
                        size_t *charlen) {
   get_tree_info(tree_num,tree_info,charlen);
}
long long int fdps_get_tree_memsize(const int tree_num) {
   return get_tree_memsize(tree_num);
}
void fdps_get_tree_time_prof(const int tree_num,
                             PS::TimeProfile *prof) {
   get_tree_time_prof(tree_num,prof);
}
void fdps_clear_tree_time_prof(const int tree_num) {
   clear_tree_time_prof(tree_num);
}
long long int fdps_get_num_interact_ep_ep_loc(const int tree_num) {
   return get_num_interact_ep_ep_loc(tree_num);
}
long long int fdps_get_num_interact_ep_sp_loc(const int tree_num) {
   return get_num_interact_ep_sp_loc(tree_num);
}
long long int fdps_get_num_interact_ep_ep_glb(const int tree_num) {
   return get_num_interact_ep_ep_glb(tree_num);
}
long long int fdps_get_num_interact_ep_sp_glb(const int tree_num) {
   return get_num_interact_ep_sp_glb(tree_num);
}
void fdps_clear_num_interact(const int tree_num) {
   return clear_num_interact(tree_num);
}
long long int fdps_get_num_tree_walk_loc(const int tree_num) {
   return get_num_tree_walk_loc(tree_num);
}
long long int fdps_get_num_tree_walk_glb(const int tree_num) {
   return get_num_tree_walk_glb(tree_num);
}

//------------------------------------------------------------------
// [Comment] A place where fdps_calc_force_*() are generated.
// fdps-autogen:calc_force_all_and_write_back;
// fdps-autogen:calc_force_all;
// fdps-autogen:calc_force_making_tree;
// fdps-autogen:calc_force_and_write_back;
//------------------------------------------------------------------

//------------------------------------------------------------------
// [Comment] A place where fdps_get_neighbor_list_*() are generated.
// fdps-autogen:get_neighbor_list;
//------------------------------------------------------------------

//------------------------------------------------------------------
// [Comment] A place where fdps_get_epj_from_id_*() are generated.
// fdps-autogen:get_epj_from_id;
//------------------------------------------------------------------

//----------------------
//  MPI comm. 
//----------------------
int fdps_get_rank() {
   return PS::Comm::getRank();
}
int fdps_get_num_procs() {
   return PS::Comm::getNumberOfProc();
}
int fdps_get_rank_multi_dim(const int id) {
   return PS::Comm::getRankMultiDim(id);
}
int fdps_get_num_procs_multi_dim(const int id) {
   return PS::Comm::getNumberOfProcMultiDim(id);
}
void fdps_get_logical_and(const bool in,
                          bool *out) {
   *out = PS::Comm::synchronizeConditionalBranchAND(in);
}
void fdps_get_logical_or(const bool in,
                         bool *out) {
   *out = PS::Comm::synchronizeConditionalBranchOR(in);
}

//------------------------------------------------------------------
// [Comment] A place where fdps_get_min_value*(), etc. are generated.
// fdps-autogen:get_min_value;
// fdps-autogen:get_max_value;
// fdps-autogen:get_sum;
// fdps-autogen:broadcast;
//------------------------------------------------------------------

double fdps_get_wtime() {
   return PS::GetWtime();
}

//----------------------
//  Utility
//----------------------
void fdps_mt_init_genrand(const int s) {
   mt_init_genrand(s);
}
int fdps_mt_genrand_int31() {
   return mt_genrand_int31();
}
double fdps_mt_genrand_real1() {
   return mt_genrand_real1();
}
double fdps_mt_genrand_real2() {
   return mt_genrand_real2();
}
double fdps_mt_genrand_real3() {
   return mt_genrand_real3();
}
double fdps_mt_genrand_res53() {
   return mt_genrand_res53();
}

//----------------------
//  Particle Mesh
//----------------------
#ifdef PARTICLE_SIMULATOR_USE_PM_MODULE
void fdps_create_pm(int *pm_num) {
   create_pm(pm_num);
}
void fdps_delete_pm(const int pm_num) {
   delete_pm(pm_num);
}
int fdps_get_pm_mesh_num() {
   return get_pm_mesh_num();
}
double fdps_get_pm_cutoff_radius() {
   return get_pm_cutoff_radius();
}
void fdps_set_dinfo_of_pm(const int pm_num, 
                          const int dinfo_num) {
   set_dinfo_of_pm(pm_num,dinfo_num);
}
void fdps_set_psys_of_pm(const int pm_num, 
                         const int psys_num,
                         const bool clear) {
   set_psys_of_pm(pm_num,psys_num,clear);
}
void fdps_get_pm_force(const int pm_num,
                       const PS::F32vec *pos,
                       PS::F32vec *force) {
   get_pm_force(pm_num,pos,force);
}
void fdps_get_pm_potential(const int pm_num, 
                           const PS::F32vec *pos,
                           PS::F32 *pot) {
   get_pm_potential(pm_num,pos,pot);
}
void fdps_calc_pm_force_only(const int pm_num) {
   calc_pm_force_only(pm_num);
}
void fdps_calc_pm_force_all_and_write_back(const int pm_num,
                                           const int psys_num,
                                           const int dinfo_num) {
   calc_pm_force_all_and_write_back(pm_num,psys_num,dinfo_num);
}
#endif

} // END of extern "C"
