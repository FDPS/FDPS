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

// CommInfo static tables
const int fdps_comm_table_size = 32;
static int fdps_comm_table_initialized = 0;
static fdps_comm_info  comm_table[fdps_comm_table_size];
static int  comm_table_used[fdps_comm_table_size];

    

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
void * fdps_get_psys_cptr(const int psys_num) {
   return get_psys_cptr(psys_num);
}
void fdps_exchange_particle(const int psys_num,
                            const int dinfo_num) {
   exchange_particle(psys_num,dinfo_num);
}
void fdps_add_particle(const int psys_num,
                       const void *cptr_to_fp) {
   add_particle(psys_num,cptr_to_fp);
}
void fdps_remove_particle(const int psys_num, 
                          const int nptcl,
                          int *ptcl_indx) {
   remove_particle(psys_num,nptcl,ptcl_indx);
}
void fdps_adjust_pos_into_root_domain(const int psys_num,
                                      const int dinfo_num) {
   adjust_pos_into_root_domain(psys_num,dinfo_num);
}
void fdps_sort_particle(const int psys_num,
                        bool (*pfunc_comp)(const void *, const void *)) {
   sort_particle(psys_num,pfunc_comp);
}
void fdps_set_psys_comm_info(int psys_num, int ci){
   set_psys_comm_info(psys_num, comm_table+ci);
}

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
   float coef_ema_;
   if ((0.0 <= coef_ema) && (coef_ema <= 1.0)) {
      coef_ema_ = coef_ema;
   } else {
      coef_ema_ = FDPS_DFLT_VAL_COEF_EMA;
   }
   init_dinfo(dinfo_num,coef_ema_);
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
int fdps_get_boundary_condition(const int dinfo_num) {
   return get_boundary_condition(dinfo_num);
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
   float weight_;
   if (weight >= 0.0) {
      weight_ = weight;
   } else {
      weight_ = (float) get_nptcl_loc(psys_num);
   }
   collect_sample_particle(dinfo_num,psys_num,clear,weight_);
}
void fdps_decompose_domain(const int dinfo_num) {
   decompose_domain(dinfo_num);
}
void fdps_decompose_domain_all(const int dinfo_num,
                               const int psys_num, 
                               const float weight) {
   float weight_;
   if (weight >= 0.0) {
      weight_ = weight;
   } else {
      weight_ = (float) get_nptcl_loc(psys_num);
   }
   decompose_domain_all(dinfo_num,psys_num,weight_);
}
void fdps_set_dinfo_comm_info(int dinfo_num, int ci){
    
   set_dinfo_comm_info(dinfo_num, comm_table+ci);
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
   float theta_;
   if (theta >= 0.0) {
      theta_ = theta;
   } else {
      theta_ = FDPS_DFLT_VAL_THETA;
   }
   int n_leaf_limit_;
   if (n_leaf_limit > 0) {
      n_leaf_limit_ = n_leaf_limit;
   } else {
      n_leaf_limit_ = FDPS_DFLT_VAL_N_LEAF_LIMIT;
   }
   int n_group_limit_;
   if (n_group_limit > 0) {
      n_group_limit_ = n_group_limit;
   } else {
      n_group_limit_ = FDPS_DFLT_VAL_N_GROUP_LIMIT;
   }
   init_tree(tree_num,nptcl,theta_,
             n_leaf_limit_,n_group_limit_);
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
void fdps_set_particle_local_tree(const int tree_num,
                                  const int psys_num,
                                  const bool clear) {
   set_particle_local_tree(tree_num,psys_num,clear);
}
void fdps_get_force(const int tree_num, 
                    const PS::S32 i,
                    const void *cptr_to_force) {
   get_force(tree_num,i,cptr_to_force);
}
void fdps_calc_force_all_and_write_back(const int tree_num,
                                        void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                        void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                        const int psys_num,
                                        const int dinfo_num,
                                        const bool clear,
                                        const enum PS_INTERACTION_LIST_MODE list_mode) {
   enum PS_INTERACTION_LIST_MODE list_mode_;
   if (list_mode >= 0) {
      list_mode_ = list_mode;
   } else {
      list_mode_ = FDPS_DFLT_VAL_LIST_MODE;
   }
   calc_force_all_and_write_back(tree_num,
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 psys_num,
                                 dinfo_num,
                                 clear,
                                 list_mode_);
}
void fdps_calc_force_all(const int tree_num,
                         void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                         void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                         const int psys_num,
                         const int dinfo_num,
                         const bool clear,
                         const enum PS_INTERACTION_LIST_MODE list_mode) {
   enum PS_INTERACTION_LIST_MODE list_mode_;
   if (list_mode >= 0) {
      list_mode_ = list_mode;
   } else {
      list_mode_ = FDPS_DFLT_VAL_LIST_MODE;
   }
   calc_force_all(tree_num,
                  pfunc_ep_ep,
                  pfunc_ep_sp,
                  psys_num,
                  dinfo_num,
                  clear,
                  list_mode_);
}
void fdps_calc_force_making_tree(const int tree_num,
                                 void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                 void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                 const int dinfo_num,
                                 const bool clear) {
   calc_force_making_tree(tree_num,
                          pfunc_ep_ep,
                          pfunc_ep_sp,
                          dinfo_num,
                          clear);
}
void fdps_calc_force_and_write_back(const int tree_num,
                                    void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                    void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                    const int psys_num,
                                    const bool clear) {
   calc_force_and_write_back(tree_num,
                             pfunc_ep_ep,
                             pfunc_ep_sp,
                             psys_num,
                             clear);
}
void fdps_get_neighbor_list(const int tree_num,
                            const PS::F64vec *pos,
                            const PS::F64 r_search,
                            int *num_epj,
                            void **cptr_to_epj) {
   get_neighbor_list(tree_num,
                     pos,
                     r_search,
                     num_epj,
                     cptr_to_epj);
}
void * fdps_get_epj_from_id(const int tree_num,
                            const PS::S64 id) {
   return get_epj_from_id(tree_num, id);
}
void fdps_set_tree_comm_info(int tree_num, int ci){
   set_tree_comm_info(tree_num, comm_table+ci);
}
void fdps_set_exchange_let_mode(int tree_num, int mode){
   set_exchange_let_mode(tree_num, mode);
}

//----------------------
//  MPI comm. 
//----------------------
#include "FDPS_comm_info.h"


void fdps_ci_raw_delete(fdps_comm_info * ci)
{
    ((PS::CommInfo*)ci)->free();
}

void fdps_ci_raw_initialized_comm_table_if_not()
{
   if (fdps_comm_table_initialized ==0){        
       for(int i=0;i<fdps_comm_table_size; i++){
           comm_table_used[i]=0;
       }
       fdps_comm_table_initialized = 1;
    }
}

int  fdps_ci_raw_register_comm(fdps_comm_info  p)
{
   fdps_ci_raw_initialized_comm_table_if_not();
   for(int i=0;i<fdps_comm_table_size; i++){
       if (comm_table_used[i]==0){
           comm_table[i]=p;
	   comm_table_used[i]=1;
           return i;
       }
   }
   return -1;
}

#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
#define MPI_Comm int
#define MPI_Fint int

#endif
fdps_comm_info  fdps_ci_raw_initialize(MPI_Comm* comm)
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    PS::CommInfo CI = PS::CommInfo( *comm);
#else
    PS::CommInfo CI = PS::CommInfo();
#endif
    return *((fdps_comm_info *) (&CI));
}
void fdps_ci_raw_set_communicator(fdps_comm_info* ci,MPI_Comm * comm)
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    ((PS::CommInfo*)(ci))->setCommunicator(*comm);
#else
    ((PS::CommInfo*)(ci))->setCommunicator();
#endif

}


int fdps_ci_initialize(MPI_Fint comm)
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Comm ccomm = MPI_Comm_f2c(comm);
#else
    MPI_Comm ccomm = comm;
#endif
    return fdps_ci_raw_register_comm(fdps_ci_raw_initialize(&ccomm));
}
void fdps_ci_set_communicator(int ci, MPI_Fint comm)
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Comm ccomm = MPI_Comm_f2c(comm);
#else
    MPI_Comm ccomm = comm;
#endif
    fdps_ci_raw_set_communicator(comm_table+ci, &ccomm);
}

void fdps_ci_delete(int ci)
{
    fdps_ci_raw_delete(comm_table+ci);
    comm_table_used[ci] = 0;
}



fdps_comm_info  fdps_ci_raw_create(fdps_comm_info  ci, int n, int rank[])
{
    PS::CommInfo CI =((PS::CommInfo*)(&ci))->create(n, rank);
    return *((fdps_comm_info *) (&CI));
}
fdps_comm_info  fdps_ci_raw_split(fdps_comm_info  ci,int color, int key)
{
    PS::CommInfo CI =((PS::CommInfo*)(&ci))->split(color, key);
    return *((fdps_comm_info *) (&CI));
}
int  fdps_ci_create(int ci, int n, int rank[])
{
    return fdps_ci_raw_register_comm(fdps_ci_raw_create(comm_table[ci], n, rank));
}
int fdps_ci_split(int ci,int color, int key)
{
    return fdps_ci_raw_register_comm(fdps_ci_raw_split(comm_table[ci], color, key));
}



int fdps_ci_raw_get_rank(fdps_comm_info * ci) {
   return ((PS::CommInfo*)(ci))->getRank();
}
int fdps_ci_get_rank(int ci) {
   return ((PS::CommInfo*)(comm_table+ci))->getRank();
}
int fdps_get_rank() {
   return PS::Comm::getRank();
}
int fdps_ci_raw_get_num_procs(fdps_comm_info * ci) {
   return ((PS::CommInfo*)(ci))->getNumberOfProc();
}
int fdps_ci_get_num_procs(int ci) {
   return ((PS::CommInfo*)(comm_table+ci))->getNumberOfProc();
}
int fdps_get_num_procs() {
   return PS::Comm::getNumberOfProc();
}
int fdps_ci_raw_get_rank_multi_dim(fdps_comm_info * ci, const int id) {
   return ((PS::CommInfo*)(ci))->getRankMultiDim(id);
}
int fdps_ci_get_rank_multi_dim(int ci, const int id) {
   return ((PS::CommInfo*)(comm_table+ci))->getRankMultiDim(id);
}
int fdps_get_rank_multi_dim(const int id) {
   return PS::Comm::getRankMultiDim(id);
}
int fdps_ci_raw_get_num_procs_multi_dim(fdps_comm_info * ci, const int id) {
   return ((PS::CommInfo*)(ci))->getNumberOfProcMultiDim(id);
}
int fdps_ci_get_num_procs_multi_dim(int ci, const int id) {
   return ((PS::CommInfo*)(comm_table+ci))->getNumberOfProcMultiDim(id);
}
int fdps_get_num_procs_multi_dim(const int id) {
   return PS::Comm::getNumberOfProcMultiDim(id);
}
bool fdps_ci_raw_get_logical_and(fdps_comm_info * ci, const bool in) {
   return ((PS::CommInfo*)(ci))->synchronizeConditionalBranchAND(in);
}
bool fdps_ci_get_logical_and(int ci, const bool in) {
   return ((PS::CommInfo*)(comm_table+ci))->synchronizeConditionalBranchAND(in);
}
bool fdps_get_logical_and(const bool in) {
   return PS::Comm::synchronizeConditionalBranchAND(in);
}
bool fdps_ci_raw_get_logical_or(fdps_comm_info * ci, const bool in) {
   return ((PS::CommInfo*)(ci))->synchronizeConditionalBranchOR(in);
}
bool fdps_ci_get_logical_or(int ci, const bool in) {
   return ((PS::CommInfo*)(comm_table+ci))->synchronizeConditionalBranchOR(in);
}
bool fdps_get_logical_or(const bool in) {
   return PS::Comm::synchronizeConditionalBranchOR(in);
}

PS::S32 fdps_ci_raw_get_min_value_s32(fdps_comm_info * ci, const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMinValue(tmp);
}
PS::S32 fdps_ci_get_min_value_s32(int ci, const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMinValue(tmp);
}
PS::S32 fdps_get_min_value_s32(const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return PS::Comm::getMinValue(tmp);
}
PS::S64 fdps_ci_raw_get_min_value_s64(fdps_comm_info * ci, const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMinValue(tmp);
}
PS::S64 fdps_ci_get_min_value_s64(int ci, const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMinValue(tmp);
}
PS::S64 fdps_get_min_value_s64(const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return PS::Comm::getMinValue(tmp);
}
PS::U32 fdps_ci_raw_get_min_value_u32(fdps_comm_info * ci, const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMinValue(tmp);
}
PS::U32 fdps_ci_get_min_value_u32(int ci, const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMinValue(tmp);
}
PS::U32 fdps_get_min_value_u32(const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return PS::Comm::getMinValue(tmp);
}
PS::U64 fdps_ci_raw_get_min_value_u64(fdps_comm_info * ci, const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMinValue(tmp);
}
PS::U64 fdps_ci_get_min_value_u64(int ci, const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMinValue(tmp);
}
PS::U64 fdps_get_min_value_u64(const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return PS::Comm::getMinValue(tmp);
}
PS::F32 fdps_ci_raw_get_min_value_f32(fdps_comm_info * ci, const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMinValue(tmp);
}
PS::F32 fdps_ci_get_min_value_f32(int ci, const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMinValue(tmp);
}
PS::F32 fdps_get_min_value_f32(const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return PS::Comm::getMinValue(tmp);
}
PS::F64 fdps_ci_raw_get_min_value_f64(fdps_comm_info * ci, const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMinValue(tmp);
}
PS::F64 fdps_ci_get_min_value_f64(int ci, const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMinValue(tmp);
}
PS::F64 fdps_get_min_value_f64(const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return PS::Comm::getMinValue(tmp);
}
void fdps_ci_raw_get_min_value_w_id_f32(fdps_comm_info * ci, const PS::F32 f_in,
                                 const int i_in,
                                 PS::F32 *f_out,
                                 int *i_out) {
    ((PS::CommInfo*)(ci))->getMinValue(f_in,i_in,*f_out,*i_out);
}
void fdps_ci_get_min_value_w_id_f32(int ci, const PS::F32 f_in,
                                 const int i_in,
                                 PS::F32 *f_out,
                                 int *i_out) {
    ((PS::CommInfo*)(comm_table+ci))->getMinValue(f_in,i_in,*f_out,*i_out);
}
void fdps_get_min_value_w_id_f32(const PS::F32 f_in,
                                 const int i_in,
                                 PS::F32 *f_out,
                                 int *i_out) {
    PS::Comm::getMinValue(f_in,i_in,*f_out,*i_out);
}
void fdps_ci_raw_get_min_value_w_id_f64(fdps_comm_info * ci, const PS::F64 f_in,
                                 const int i_in,
                                 PS::F64 *f_out,
                                 int *i_out) {
    ((PS::CommInfo*)(ci))->getMinValue(f_in,i_in,*f_out,*i_out);
}
void fdps_ci_get_min_value_w_id_f64(int ci, const PS::F64 f_in,
                                 const int i_in,
                                 PS::F64 *f_out,
                                 int *i_out) {
    ((PS::CommInfo*)(comm_table+ci))->getMinValue(f_in,i_in,*f_out,*i_out);
}
void fdps_get_min_value_w_id_f64(const PS::F64 f_in,
                                 const int i_in,
                                 PS::F64 *f_out,
                                 int *i_out) {
    PS::Comm::getMinValue(f_in,i_in,*f_out,*i_out);
}

PS::S32 fdps_ci_raw_get_max_value_s32(fdps_comm_info * ci, const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMaxValue(tmp);
}
PS::S32 fdps_ci_get_max_value_s32(int ci, const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMaxValue(tmp);
}
PS::S32 fdps_get_max_value_s32(const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return PS::Comm::getMaxValue(tmp);
}
PS::S64 fdps_ci_raw_get_max_value_s64(fdps_comm_info * ci, const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMaxValue(tmp);
}
PS::S64 fdps_ci_get_max_value_s64(int ci, const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMaxValue(tmp);
}
PS::S64 fdps_get_max_value_s64(const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return PS::Comm::getMaxValue(tmp);
}
PS::U32 fdps_ci_raw_get_max_value_u32(fdps_comm_info * ci, const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMaxValue(tmp);
}
PS::U32 fdps_ci_get_max_value_u32(int ci, const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMaxValue(tmp);
}
PS::U32 fdps_get_max_value_u32(const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return PS::Comm::getMaxValue(tmp);
}
PS::U64 fdps_ci_raw_get_max_value_u64(fdps_comm_info * ci, const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMaxValue(tmp);
}
PS::U64 fdps_ci_get_max_value_u64(int ci, const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMaxValue(tmp);
}
PS::U64 fdps_get_max_value_u64(const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return PS::Comm::getMaxValue(tmp);
}
PS::F32 fdps_ci_raw_get_max_value_f32(fdps_comm_info * ci, const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMaxValue(tmp);
}
PS::F32 fdps_ci_get_max_value_f32(int ci, const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMaxValue(tmp);
}
PS::F32 fdps_get_max_value_f32(const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return PS::Comm::getMaxValue(tmp);
}
PS::F64 fdps_ci_raw_get_max_value_f64(fdps_comm_info * ci, const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getMaxValue(tmp);
}
PS::F64 fdps_ci_get_max_value_f64(int ci, const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getMaxValue(tmp);
}
PS::F64 fdps_get_max_value_f64(const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return PS::Comm::getMaxValue(tmp);
}
void fdps_ci_raw_get_max_value_w_id_f32(fdps_comm_info * ci, const PS::F32 f_in,
                                 const int i_in,
                                 PS::F32 *f_out,
                                 int *i_out) {
    ((PS::CommInfo*)(ci))->getMaxValue(f_in,i_in,*f_out,*i_out);
}
void fdps_ci_get_max_value_w_id_f32(int ci, const PS::F32 f_in,
                                 const int i_in,
                                 PS::F32 *f_out,
                                 int *i_out) {
    ((PS::CommInfo*)(comm_table+ci))->getMaxValue(f_in,i_in,*f_out,*i_out);
}
void fdps_get_max_value_w_id_f32(const PS::F32 f_in,
                                 const int i_in,
                                 PS::F32 *f_out,
                                 int *i_out) {
    PS::Comm::getMaxValue(f_in,i_in,*f_out,*i_out);
}
void fdps_ci_raw_get_max_value_w_id_f64(fdps_comm_info * ci, const PS::F64 f_in,
                                 const int i_in,
                                 PS::F64 *f_out,
                                 int *i_out) {
    ((PS::CommInfo*)(ci))->getMaxValue(f_in,i_in,*f_out,*i_out);
}
void fdps_ci_get_max_value_w_id_f64(int ci, const PS::F64 f_in,
                                 const int i_in,
                                 PS::F64 *f_out,
                                 int *i_out) {
    ((PS::CommInfo*)(comm_table+ci))->getMaxValue(f_in,i_in,*f_out,*i_out);
}
void fdps_get_max_value_w_id_f64(const PS::F64 f_in,
                                 const int i_in,
                                 PS::F64 *f_out,
                                 int *i_out) {
    PS::Comm::getMaxValue(f_in,i_in,*f_out,*i_out);
}

PS::S32 fdps_ci_raw_get_sum_s32(fdps_comm_info * ci, const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getSum(tmp);
}
PS::S32 fdps_ci_get_sum_s32(int ci, const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getSum(tmp);
}
PS::S32 fdps_get_sum_s32(const PS::S32 f_in) {
    PS::S32 tmp = f_in;
    return PS::Comm::getSum(tmp);
}
PS::S64 fdps_ci_raw_get_sum_s64(fdps_comm_info * ci, const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getSum(tmp);
}
PS::S64 fdps_ci_get_sum_s64(int ci, const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getSum(tmp);
}
PS::S64 fdps_get_sum_s64(const PS::S64 f_in) {
    PS::S64 tmp = f_in;
    return PS::Comm::getSum(tmp);
}
PS::U32 fdps_ci_raw_get_sum_u32(fdps_comm_info * ci, const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getSum(tmp);
}
PS::U32 fdps_ci_get_sum_u32(int ci, const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getSum(tmp);
}
PS::U32 fdps_get_sum_u32(const PS::U32 f_in) {
    PS::U32 tmp = f_in;
    return PS::Comm::getSum(tmp);
}
PS::U64 fdps_ci_raw_get_sum_u64(fdps_comm_info * ci, const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getSum(tmp);
}
PS::U64 fdps_ci_get_sum_u64(int ci, const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getSum(tmp);
}
PS::U64 fdps_get_sum_u64(const PS::U64 f_in) {
    PS::U64 tmp = f_in;
    return PS::Comm::getSum(tmp);
}
PS::F32 fdps_ci_raw_get_sum_f32(fdps_comm_info * ci, const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getSum(tmp);
}
PS::F32 fdps_ci_get_sum_f32(int ci, const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getSum(tmp);
}
PS::F32 fdps_get_sum_f32(const PS::F32 f_in) {
    PS::F32 tmp = f_in;
    return PS::Comm::getSum(tmp);
}
PS::F64 fdps_ci_raw_get_sum_f64(fdps_comm_info * ci, const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return ((PS::CommInfo*)(ci))->getSum(tmp);
}
PS::F64 fdps_ci_get_sum_f64(int ci, const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return ((PS::CommInfo*)(comm_table+ci))->getSum(tmp);
}
PS::F64 fdps_get_sum_f64(const PS::F64 f_in) {
    PS::F64 tmp = f_in;
    return PS::Comm::getSum(tmp);
}

void fdps_ci_raw_broadcast_s32(fdps_comm_info * ci, PS::S32 *val, int n, int src) {
    ((PS::CommInfo*)(ci))->broadcast(val,n,src);
}
void fdps_ci_broadcast_s32(int ci, PS::S32 *val, int n, int src) {
    ((PS::CommInfo*)(comm_table+ci))->broadcast(val,n,src);
}
void fdps_broadcast_s32(PS::S32 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}
void fdps_ci_raw_broadcast_s64(fdps_comm_info * ci, PS::S64 *val, int n, int src) {
    ((PS::CommInfo*)(ci))->broadcast(val,n,src);
}
void fdps_ci_broadcast_s64(int ci, PS::S64 *val, int n, int src) {
    ((PS::CommInfo*)(comm_table+ci))->broadcast(val,n,src);
}
void fdps_broadcast_s64(PS::S64 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}
void fdps_ci_raw_broadcast_u32(fdps_comm_info * ci, PS::U32 *val, int n, int src) {
    ((PS::CommInfo*)(ci))->broadcast(val,n,src);
}
void fdps_ci_broadcast_u32(int ci, PS::U32 *val, int n, int src) {
    ((PS::CommInfo*)(comm_table+ci))->broadcast(val,n,src);
}
void fdps_broadcast_u32(PS::U32 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}
void fdps_ci_raw_broadcast_u64(fdps_comm_info * ci, PS::U64 *val, int n, int src) {
    ((PS::CommInfo*)(ci))->broadcast(val,n,src);
}
void fdps_ci_broadcast_u64(int ci, PS::U64 *val, int n, int src) {
    ((PS::CommInfo*)(comm_table+ci))->broadcast(val,n,src);
}
void fdps_broadcast_u64(PS::U64 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}
void fdps_ci_raw_broadcast_f32(fdps_comm_info * ci, PS::F32 *val, int n, int src) {
    ((PS::CommInfo*)(ci))->broadcast(val,n,src);
}
void fdps_ci_broadcast_f32(int ci, PS::F32 *val, int n, int src) {
    ((PS::CommInfo*)(comm_table+ci))->broadcast(val,n,src);
}
void fdps_broadcast_f32(PS::F32 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}
void fdps_ci_raw_broadcast_f64(fdps_comm_info * ci, PS::F64 *val, int n, int src) {
    ((PS::CommInfo*)(ci))->broadcast(val,n,src);
}
void fdps_ci_broadcast_f64(int ci, PS::F64 *val, int n, int src) {
    ((PS::CommInfo*)(comm_table+ci))->broadcast(val,n,src);
}
void fdps_broadcast_f64(PS::F64 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

double fdps_ci_raw_get_wtime(fdps_comm_info * ci) {
   return PS::GetWtime();
}
double fdps_ci_get_wtime(int ci) {
   return PS::GetWtime();
}
double fdps_get_wtime() {
   return PS::GetWtime();
}

void fdps_ci_raw_barrier(fdps_comm_info * ci) {
   ((PS::CommInfo*)(ci))->barrier();
}
void fdps_ci_barrier(int ci) {
   ((PS::CommInfo*)(comm_table+ci))->barrier();
}
void fdps_barrier() {
   PS::Comm::barrier();
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

void fdps_create_mtts(int * mtts_num) {
   create_mtts(mtts_num);
}
void fdps_delete_mtts(const int mtts_num) {
   delete_mtts(mtts_num);
}
void fdps_mtts_init_genrand(const int mtts_num,
                            const int s) {
   mtts_init_genrand(mtts_num,s);
}
int fdps_mtts_genrand_int31(const int mtts_num) {
   return mtts_genrand_int31(mtts_num);
}
double fdps_mtts_genrand_real1(const int mtts_num) {
   return mtts_genrand_real1(mtts_num);
}
double fdps_mtts_genrand_real2(const int mtts_num) {
   return mtts_genrand_real2(mtts_num);
}
double fdps_mtts_genrand_real3(const int mtts_num) {
   return mtts_genrand_real3(mtts_num);
}
double fdps_mtts_genrand_res53(const int mtts_num) {
   return mtts_genrand_res53(mtts_num);
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
   set_dinfo_of_pm(pom_num,dinfo_num);
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
