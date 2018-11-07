#pragma once
/* Standard headers */
#include <string>
#include <vector>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.hpp"

namespace FDPS_Manipulators {
   //--------------------------------
   // Parameters
   //--------------------------------
   enum PS_BOUNDARY_CONDITION {
      BC_OPEN,
      BC_PERIODIC_X,
      BC_PERIODIC_Y,
      BC_PERIODIC_Z,
      BC_PERIODIC_XY,
      BC_PERIODIC_XZ,
      BC_PERIODIC_YZ,
      BC_PERIODIC_XYZ,
      BC_SHEARING_BOX,
      BC_USER_DEFINED
   };
   enum PS_INTERACTION_LIST_MODE {
      MAKE_LIST,
      MAKE_LIST_FOR_REUSE,
      REUSE_LIST
   };

   //--------------------------------
   // Initializer of FDPS_Manipulators
   //--------------------------------
   extern void Initialize(int argc, char *argv[]);
   //--------------------------------
   // FDPS Initializer/Finalizer
   //--------------------------------
   extern void PS_Initialize();
   extern void PS_Finalize();
   extern void PS_Abort(const int err_num);
   //--------------------------------
   // ParticleSystem manipulators
   //--------------------------------
   extern void create_psys(int *psys_num,
                           char *psys_info);
   extern void delete_psys(const int psys_num);
   extern void init_psys(const int psys_num);
   extern void get_psys_info(const int psys_num,
                             char *psys_info,
                             size_t *charlen);
   extern long long int get_psys_memsize(const int psys_num);
   extern void get_psys_time_prof(const int psys_num, 
                                  PS::TimeProfile *prof);
   extern void clear_psys_time_prof(const int psys_num);
   extern void set_nptcl_smpl(const int psys_num,
                              const int numPtcl);
   extern void set_nptcl_loc(const int psys_num,
                             const int numPtcl);
   extern int get_nptcl_loc(const int psys_num);
   extern int get_nptcl_glb(const int psys_num);
   extern void * get_psys_cptr(const int psys_num);
   extern void exchange_particle(const int psys_num,
                                 const int dinfo_num);
   extern void add_particle(const int psys_num,
                            const void *cptr_to_fp);
   extern void remove_particle(const int psys_num, 
                               const int numPtcl,
                               int *ptcl_indx);
   extern void adjust_pos_into_root_domain(const int psys_num,
                                           const int dinfo_num);
   extern void sort_particle(const int psys_num,
                             bool (*pfunc_comp)(const void *, const void *));

   //----------------------------
   // DomainInfo manipulators
   //----------------------------
   extern void create_dinfo(int *dinfo_num);
   extern void delete_dinfo(const int dinfo_num);
   extern void init_dinfo(const int dinfo_num,
                          const float coef_ema);
   extern void get_dinfo_time_prof(const int dinfo_num, 
                                   PS::TimeProfile *prof);
   extern void clear_dinfo_time_prof(const int dinfo_num);
   extern void set_nums_domain(const int dinfo_num,
                               const int nx,
                               const int ny,
                               const int nz);
   extern void set_boundary_condition(const int dinfo_num, 
                                      const enum PS_BOUNDARY_CONDITION bc);
   extern int get_boundary_condition(const int dinfo_num);
   extern void set_pos_root_domain(const int dinfo_num,
                                   const PS::F32vec *low,
                                   const PS::F32vec *high);
   extern void collect_sample_particle(const int dinfo_num,
                                       const int psys_num, 
                                       const bool clear,
                                       const float weight);
   extern void decompose_domain(const int dinfo_num);
   extern void decompose_domain_all(const int dinfo_num,
                                    const int psys_num,
                                    const float weight);

   //----------------------------
   // TreeForForce manipulators
   //----------------------------
   extern void create_tree(int *tree_num,
                           char *tree_info);
   extern void delete_tree(const int tree_num);
   extern void init_tree(const int tree_num,
                         const int numPtcl,
                         const float theta, 
                         const int n_leaf_limit,
                         const int n_group_limit);
   extern void get_tree_info(const int tree_num,
                             char *tree_info,
                             size_t *charlen);
   extern long long int get_tree_memsize(const int tree_num);
   extern void get_tree_time_prof(const int tree_num, 
                                  PS::TimeProfile *prof);
   extern void clear_tree_time_prof(const int tree_num);
   extern long long int get_num_interact_ep_ep_loc(const int tree_num);
   extern long long int get_num_interact_ep_sp_loc(const int tree_num);
   extern long long int get_num_interact_ep_ep_glb(const int tree_num);
   extern long long int get_num_interact_ep_sp_glb(const int tree_num);
   extern void clear_num_interact(const int tree_num);
   extern long long int get_num_tree_walk_loc(const int tree_num);
   extern long long int get_num_tree_walk_glb(const int tree_num);
   extern void set_particle_local_tree(const int tree_num, 
                                       const int psys_num,
                                       const bool clear);
   extern void get_force(const int tree_num, 
                         const PS::S32 i,
                         const void *cptr_to_force);
   extern void calc_force_all_and_write_back(const int tree_num,
                                             void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                             void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                             const int psys_num,
                                             const int dinfo_num,
                                             const bool clear,
                                             const enum PS_INTERACTION_LIST_MODE list_mode);
   extern void calc_force_all(const int tree_num,
                              void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                              void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                              const int psys_num,
                              const int dinfo_num,
                              const bool clear,
                              const enum PS_INTERACTION_LIST_MODE list_mode);
   extern void calc_force_making_tree(const int tree_num,
                                      void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                      void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                      const int dinfo_num,
                                      const bool clear);
   extern void calc_force_and_write_back(const int tree_num,
                                         void *(pfunc_ep_ep)(void *, int, void *, int, void *),
                                         void *(pfunc_ep_sp)(void *, int, void *, int, void *),
                                         const int psys_num,
                                         const bool clear);
   extern void get_neighbor_list(const int tree_num,
                                 const PS::F64vec *pos,
                                 const PS::F64 r_search,
                                 int *num_epj,
                                 void **cptr_to_epj);
   extern void * get_epj_from_id(const int tree_num,
                                 const PS::S64 id);


   //----------------------------
   // Utility functions
   //----------------------------
   extern void mt_init_genrand(int s);
   extern int mt_genrand_int31(void);
   extern double mt_genrand_real1(void);
   extern double mt_genrand_real2(void);
   extern double mt_genrand_real3(void);
   extern double mt_genrand_res53(void);

   extern void create_mtts(int* mtts_num);
   extern void delete_mtts(const int mtts_num);
   extern void mtts_init_genrand(const int mtts_num, const int s);
   extern int mtts_genrand_int31(const int mtts_num);
   extern double mtts_genrand_real1(const int mtts_num);
   extern double mtts_genrand_real2(const int mtts_num);
   extern double mtts_genrand_real3(const int mtts_num);
   extern double mtts_genrand_res53(const int mtts_num);

   //----------------------------
   // Particle Mesh
   //----------------------------
#ifdef PARTICLE_SIMULATOR_USE_PM_MODULE
   extern void create_pm(int *pm_num);
   extern void delete_pm(const int pm_num);
   extern int get_pm_mesh_num();
   extern double get_pm_cutoff_radius();
   extern void set_dinfo_of_pm(const int pm_num,
                               const int dinfo_num);
   extern void set_psys_of_pm(const int pm_num,
                              const int psys_num,
                              const bool clear);
   extern void get_pm_force(const int pm_num,
                            const PS::F32vec *pos,
                            PS::F32vec *force);
   extern void get_pm_potential(const int pm_num,
                                const PS::F32vec *pos,
                                PS::F32 *pot);
   extern void calc_pm_force_only(const int pm_num);
   extern void calc_pm_force_all_and_write_back(const int pm_num,
                                                const int psys_num,
                                                const int dinfo_num);
#endif
}
