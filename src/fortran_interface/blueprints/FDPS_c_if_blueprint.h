#pragma once
/* Standard headers */
#include <stdio.h>
/* FDPS headers */
#include "FDPS_enum.h"
#include "FDPS_vector.h"
#include "FDPS_matrix.h"
#include "FDPS_super_particle.h"
#include "FDPS_time_profile.h"

//----------------------
//  Basic
//----------------------
void fdps_initialize();
void fdps_finalize();
void fdps_abort(const int err_num);

//----------------------
//  Particle System 
//----------------------
void fdps_create_psys(int *psys_num,
                      char *psys_info);
void fdps_delete_psys(const int psys_num);
void fdps_init_psys(const int psys_num);
void fdps_get_psys_info(const int psys_num,
                        char *psys_info,
                        size_t *charlen);
long long int fdps_get_psys_memsize(const int psys_num);
void fdps_get_psys_time_prof(const int psys_num,
                             fdps_time_prof *prof);
void fdps_clear_psys_time_prof(const int psys_num);
void fdps_set_nptcl_smpl(const int psys_num,
                         const int nptcl);
void fdps_set_nptcl_loc(const int psys_num,
                        const int nptcl);
int fdps_get_nptcl_loc(const int psys_num);
int fdps_get_nptcl_glb(const int psys_num);
void * fdps_get_psys_cptr(const int psys_num);
void fdps_exchange_particle(const int psys_num,
                            const int dinfo_num);
void fdps_add_particle(const int psys_num,
                       const void *cptr_to_fp);
void fdps_remove_particle(const int psys_num, 
                          const int nptcl,
                          int *ptcl_indx);
void fdps_adjust_pos_into_root_domain(const int psys_num,
                                      const int dinfo_num);
void fdps_sort_particle(const int psys_num,
                        _Bool (*pfunc_comp)(const void *, const void *));

//----------------------
//  Domain Info
//----------------------
void fdps_create_dinfo(int *dinfo_num);
void fdps_delete_dinfo(const int dinfo_num);
void fdps_init_dinfo(const int dinfo_num,
                     const float coef_ema);
void fdps_get_dinfo_time_prof(const int dinfo_num, 
                              fdps_time_prof *prof);
void fdps_clear_dinfo_time_prof(const int dinfo_num);
void fdps_set_nums_domain(const int dinfo_num,
                          const int nx,
                          const int ny,
                          const int nz);
void fdps_set_boundary_condition(const int dinfo_num,
                                 const int bc);
int fdps_get_boundary_condition(const int dinfo_num);
void fdps_set_pos_root_domain(const int dinfo_num,
                              const fdps_f32vec *low,
                              const fdps_f32vec *high);
void fdps_collect_sample_particle(const int dinfo_num,
                                  const int psys_num, 
                                  const _Bool clear,
                                  const float weight);
void fdps_decompose_domain(const int dinfo_num);
void fdps_decompose_domain_all(const int dinfo_num,
                               const int psys_num, 
                               const float weight);

//----------------------
//  Tree 
//----------------------
void fdps_create_tree(int *tree_num,
                      char *tree_info);
void fdps_delete_tree(const int tree_num);
void fdps_init_tree(const int tree_num,
                    const int nptcl,
                    const float theta,
                    const int n_leaf_limit,
                    const int n_group_limit);
void fdps_get_tree_info(const int tree_num,
                        char *tree_info,
                        size_t *charlen);
long long int fdps_get_tree_memsize(const int tree_num);
void fdps_get_tree_time_prof(const int tree_num,
                             fdps_time_prof *prof);
void fdps_clear_tree_time_prof(const int tree_num);
long long int fdps_get_num_interact_ep_ep_loc(const int tree_num);
long long int fdps_get_num_interact_ep_sp_loc(const int tree_num);
long long int fdps_get_num_interact_ep_ep_glb(const int tree_num);
long long int fdps_get_num_interact_ep_sp_glb(const int tree_num);
void fdps_clear_num_interact(const int tree_num);
long long int fdps_get_num_tree_walk_loc(const int tree_num);
long long int fdps_get_num_tree_walk_glb(const int tree_num);
void fdps_set_particle_local_tree(const int tree_num,
                                  const int psys_num,
                                  const _Bool clear);
void fdps_get_force(const int tree_num, 
                    const fdps_s32 i,
                    const void *cptr_to_force);
void fdps_calc_force_all_and_write_back(const int tree_num,
                                        void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                        void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                        const int psys_num,
                                        const int dinfo_num,
                                        const _Bool clear,
                                        const int list_mode);
void fdps_calc_force_all(const int tree_num,
                         void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                         void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                         const int psys_num,
                         const int dinfo_num,
                         const _Bool clear,
                         const int list_mode);
void fdps_calc_force_making_tree(const int tree_num,
                                 void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                 void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                 const int dinfo_num,
                                 const _Bool clear);
void fdps_calc_force_and_write_back(const int tree_num,
                                    void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                    void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                    const int psys_num,
                                    const _Bool clear);
void fdps_get_neighbor_list(const int tree_num,
                            const fdps_f64vec *pos,
                            const fdps_f64 r_search,
                            int *num_epj,
                            void **cptr_to_epj);
void * fdps_get_epj_from_id(const int tree_num,
                            const fdps_s64 id);

//----------------------
//  MPI comm. 
//----------------------
int fdps_get_rank();
int fdps_get_num_procs();
int fdps_get_rank_multi_dim(const int id);
int fdps_get_num_procs_multi_dim(const int id);
_Bool fdps_get_logical_and(const _Bool in);
_Bool fdps_get_logical_or(const _Bool in);

fdps_s32 fdps_get_min_value_s32(const fdps_s32 f_in);
fdps_s64 fdps_get_min_value_s64(const fdps_s64 f_in);
fdps_u32 fdps_get_min_value_u32(const fdps_u32 f_in);
fdps_u64 fdps_get_min_value_u64(const fdps_u64 f_in);
fdps_f32 fdps_get_min_value_f32(const fdps_f32 f_in);
fdps_f64 fdps_get_min_value_f64(const fdps_f64 f_in);
void fdps_get_min_value_w_id_f32(const fdps_f32 f_in,
                                 const int i_in,
                                 fdps_f32 *f_out,
                                 int *i_out);
void fdps_get_min_value_w_id_f64(const fdps_f64 f_in,
                                 const int i_in,
                                 fdps_f64 *f_out,
                                 int *i_out);

fdps_s32 fdps_get_max_value_s32(const fdps_s32 f_in);
fdps_s64 fdps_get_max_value_s64(const fdps_s64 f_in);
fdps_u32 fdps_get_max_value_u32(const fdps_u32 f_in);
fdps_u64 fdps_get_max_value_u64(const fdps_u64 f_in);
fdps_f32 fdps_get_max_value_f32(const fdps_f32 f_in);
fdps_f64 fdps_get_max_value_f64(const fdps_f64 f_in);
void fdps_get_max_value_w_id_f32(const fdps_f32 f_in,
                                 const int i_in,
                                 fdps_f32 *f_out,
                                 int *i_out);
void fdps_get_max_value_w_id_f64(const fdps_f64 f_in,
                                 const int i_in,
                                 fdps_f64 *f_out,
                                 int *i_out);

fdps_s32 fdps_get_sum_s32(const fdps_s32 f_in);
fdps_s64 fdps_get_sum_s64(const fdps_s64 f_in);
fdps_u32 fdps_get_sum_u32(const fdps_u32 f_in);
fdps_u64 fdps_get_sum_u64(const fdps_u64 f_in);
fdps_f32 fdps_get_sum_f32(const fdps_f32 f_in);
fdps_f64 fdps_get_sum_f64(const fdps_f64 f_in);

void fdps_broadcast_s32(fdps_s32 *val, int n, int src);
void fdps_broadcast_s64(fdps_s64 *val, int n, int src);
void fdps_broadcast_u32(fdps_u32 *val, int n, int src);
void fdps_broadcast_u64(fdps_u64 *val, int n, int src);
void fdps_broadcast_f32(fdps_f32 *val, int n, int src);
void fdps_broadcast_f64(fdps_f64 *val, int n, int src);

double fdps_get_wtime();
void fdps_barrier();

//----------------------
//  Utility
//----------------------
void fdps_mt_init_genrand(const int s);
int fdps_mt_genrand_int31();
double fdps_mt_genrand_real1();
double fdps_mt_genrand_real2();
double fdps_mt_genrand_real3();
double fdps_mt_genrand_res53();

void fdps_create_mtts(int * mtts_num);
void fdps_delete_mtts(const int mtts_num);
void fdps_mtts_init_genrand(const int mtts_num,
                            const int s);
int fdps_mtts_genrand_int31(const int mtts_num);
double fdps_mtts_genrand_real1(const int mtts_num);
double fdps_mtts_genrand_real2(const int mtts_num);
double fdps_mtts_genrand_real3(const int mtts_num);
double fdps_mtts_genrand_res53(const int mtts_num);

//----------------------
//  Particle Mesh
//----------------------
#ifdef PARTICLE_SIMULATOR_USE_PM_MODULE
void fdps_create_pm(int *pm_num);
void fdps_delete_pm(const int pm_num);
int fdps_get_pm_mesh_num();
double fdps_get_pm_cutoff_radius();
void fdps_set_dinfo_of_pm(const int pm_num, 
                          const int dinfo_num);
void fdps_set_psys_of_pm(const int pm_num, 
                         const int psys_num,
                         const _Bool clear);
void fdps_get_pm_force(const int pm_num,
                       const fdps_f32vec *pos,
                       fdps_f32vec *force);
void fdps_get_pm_potential(const int pm_num, 
                           const fdps_f32vec *pos,
                           fdps_f32 *pot);
void fdps_calc_pm_force_only(const int pm_num);
void fdps_calc_pm_force_all_and_write_back(const int pm_num,
                                           const int psys_num,
                                           const int dinfo_num);
#endif
