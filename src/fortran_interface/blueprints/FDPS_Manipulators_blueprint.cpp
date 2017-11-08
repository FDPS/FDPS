/* Standard headers */
#include <string.h>
#include <string>
#include <cstring>
#include <vector>
#include <typeinfo>
/* FDPS headers */
#include <particle_simulator.hpp>
#ifdef PARTICLE_SIMULATOR_USE_PM_MODULE
#include <particle_mesh.hpp>
#include <particle_mesh_class.hpp>
#include <param.h>
#include <param_fdps.h>
#endif
/* User-defined headers */
#include "FDPS_Manipulators.h"

namespace FDPS_Manipulators {
   //===================
   // Private classes 
   //===================
   class Psys_Data {
   public:
      int id;
      void * ptr;
      std::string info;
   };

   class Dinfo_Data {
   public:
      int id;
      void * ptr;
   };

   class Tree_Data {
   public:
      int id;
      void * ptr;
      std::string info;
      char * info_char;
      std::string type;
      std::string Force;
      std::string EPI;
      std::string EPJ;
      std::string mode;
   };

   class PM_Data {
   public:
      int id;
      void * ptr;
   };

   class NBS_Target_Data {
   public:
      PS::F64vec pos;
      PS::F64 r_search;

      PS::F64vec getPos() const {
         return this->pos;
      }

      PS::F64 getRSearch() const {
         return this->r_search;
      }
   };

   //===================
   // Private variables 
   //===================
   // -(command-line arguments)
   static int argc_save;
   static char **argv_save;
   // -(ParticleSystem type)
   static int num_psys;
   static int num_psys_creation;
   static std::vector<Psys_Data> psys_vector;
   // -(DomainInfo type)
   static int num_dinfo;
   static int num_dinfo_creation;
   static std::vector<Dinfo_Data> dinfo_vector;
   // -(TreeForForce type)
   static int num_tree;
   static int num_tree_creation;
   static std::vector<Tree_Data> tree_vector;
   // -(ParticleMesh class)
   static int num_pm;
   static int num_pm_creation;
   static std::vector<PM_Data> pm_vector;

   //=======================
   // Functions
   //=======================
   //-(Initializer)
   void Initialize(int argc, char *argv[]) {
      argc_save = argc;
      argv_save = argv;
      num_psys  = 0;
      num_dinfo = 0;
      num_tree  = 0;
      num_pm    = 0;
      num_psys_creation = 0;
      num_dinfo_creation = 0;
      num_tree_creation = 0;
      num_pm_creation = 0;
   }
   //-(FDPS initilizer/finalizer)
   void PS_Initialize() {
      PS::Initialize(argc_save,argv_save);
      // Display # of MPI processes and threads
      PS::S32 nprocs = PS::Comm::getNumberOfProc();
      PS::S32 nthrds = PS::Comm::getNumberOfThread();
      PS::S32 myrank = PS::Comm::getRank();
      if (myrank == 0) {
      std::cout << "====================================" << std::endl
                << " Paralleization infomation:"          << std::endl
                << "   # of processes is " << nprocs      << std::endl
                << "   # of thread is    " << nthrds      << std::endl
                << "====================================" << std::endl;
      }
   }
   void PS_Finalize() {
      PS::Finalize();
   }
   void PS_Abort(const int err_num) {
      PS::Abort(err_num);
   }
   //--------------------------------
   // ParticleSystem manipulators
   //--------------------------------
   void create_psys(int *psys_num,
                    char *psys_info) {
      std::string psys_info_ = psys_info;
      Psys_Data psys_data;
      //-----------------------------------------------
      // fdps-autogen:create_psys;
      //-----------------------------------------------
      psys_vector.push_back(psys_data);
      *psys_num = psys_data.id;
      num_psys_creation++;
      num_psys++;
   }
   void delete_psys(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            it = psys_vector.erase(it);
            num_psys--;
            break;
         }
      }
   }
   void init_psys(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:init_psys;
            //-----------------------------------------------
         }
      }
   }
   void get_psys_info(const int psys_num,
                      char *psys_info,
                      size_t *charlen) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            std::strcpy(psys_info,it->info.c_str());
            *charlen = std::strlen(psys_info);
            break;
         }
      }
   }
   long long int get_psys_memsize(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_psys_memsize;
            //-----------------------------------------------
         }
      }
   }
   void get_psys_time_prof(const int psys_num,
                           PS::TimeProfile *prof) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_psys_time_prof;
            //-----------------------------------------------
         }
      }
   }
   void clear_psys_time_prof(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:clear_psys_time_prof;
            //-----------------------------------------------
         }
      }
   }
   void set_nptcl_smpl(const int psys_num,
                       const int numPtcl){
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:set_nptcl_smpl;
            //-----------------------------------------------
         }
      }
   }
   void set_nptcl_loc(const int psys_num,
                      const int numPtcl){
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:set_nptcl_loc;
            //-----------------------------------------------
         }
      }
   }
   int get_nptcl_loc(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_nptcl_loc;
            //-----------------------------------------------
         }
      }
   }
   int get_nptcl_glb(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_nptcl_glb;
            //-----------------------------------------------
         }
      }
   }
   void get_psys_cptr(const int psys_num,
                      void **cptr) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_psys_cptr;
            //-----------------------------------------------
         }
      }
   }
   void exchange_particle(const int psys_num,
                          const int dinfo_num) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
            break;
         }
      }
      // Get psys and issue exchangeParticle()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:exchange_particle;
            //-----------------------------------------------
         }
      }
   }
   //-----------------------------------------------------------------------------------
   // fdps-autogen:add_particle;
   //-----------------------------------------------------------------------------------
   void remove_particle(const int psys_num, 
                        const int numPtcl,
                        int *ptcl_indx) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:remove_particle;
            //-----------------------------------------------
         }
      }
   }
   void adjust_pos_into_root_domain(const int psys_num, 
                                    const int dinfo_num) {
      //* Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
            break;
         }
      }
      //* Execute adjustPositionIntoRootDomain()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:adjust_pos_into_root_domain;
            //-----------------------------------------------
         }
      }
   }
   //-----------------------------------------------------------------------------------
   // fdps-autogen:sort_particle;
   //-----------------------------------------------------------------------------------


   //----------------------------
   // DomainInfo manipulators
   //----------------------------
   void create_dinfo(int *dinfo_num) {
      Dinfo_Data dinfo_data;
      PS::DomainInfo *dinfo;
      dinfo = new PS::DomainInfo;
      dinfo_data.id = num_dinfo_creation;
      dinfo_data.ptr = (void *) dinfo;
      dinfo_vector.push_back(dinfo_data);
      *dinfo_num = dinfo_data.id;
      num_dinfo_creation++;
      num_dinfo++;
   }
   void delete_dinfo(const int dinfo_num) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            it = dinfo_vector.erase(it);
            num_dinfo--;
            break;
         }
      }
   }
   void init_dinfo(const int dinfo_num,
                   const float coef_ema) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->initialize(coef_ema);
            break;
         }
      }
   }
   void get_dinfo_time_prof(const int dinfo_num,
                            PS::TimeProfile *prof) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            *prof = dinfo->getTimeProfile();
            break;
         }
      }
   }
   void clear_dinfo_time_prof(const int dinfo_num) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->clearTimeProfile();
            break;
         }
      }
   }
   void set_nums_domain(const int dinfo_num,
                        const int nx,
                        const int ny,
                        const int nz) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->setDomain(nx,ny,nz);
            break;
         }
      }
   }
   void set_boundary_condition(const int dinfo_num,
                               const enum PS_BOUNDARY_CONDITION bc) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            if (bc == BC_OPEN) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
            } else if (bc == BC_PERIODIC_X) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
            } else if (bc == BC_PERIODIC_Y) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Y);
            } else if (bc == BC_PERIODIC_Z) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Z);
            } else if (bc == BC_PERIODIC_XY) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
            } else if (bc == BC_PERIODIC_XZ) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XZ);
            } else if (bc == BC_PERIODIC_YZ) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_YZ);
            } else if (bc == BC_PERIODIC_XYZ) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
            } else if (bc == BC_SHEARING_BOX) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_SHEARING_BOX);
            } else if (bc == BC_USER_DEFINED) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_USER_DEFINED);
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::string errmsg,func_name;
                  errmsg = "Unknown boundary condition is specified.";
                  func_name = "set_boundary_condition";
                  print_errmsg(errmsg,func_name);
               } 
               PS::Abort(-1);
               std::exit(EXIT_FAILURE);
            }
            break;
         }
      }
   }
   void set_pos_root_domain(const int dinfo_num,
                            const PS::F32vec *low,
                            const PS::F32vec *high) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->setPosRootDomain(*low,*high);
            break;
         }
      }
   }
   void collect_sample_particle(const int dinfo_num,
                                const int psys_num,
                                const bool clear,
                                const float weight) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
            break;
         }
      }
      // Get psys and issue collectSampleParticle()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:collect_sample_particle;
            //-----------------------------------------------
         }
      }
   }
   void decompose_domain(const int dinfo_num) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->decomposeDomain();
            break;
         }
      }
   }
   void decompose_domain_all(const int dinfo_num,
                             const int psys_num,
                             const float weight) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
            break;
         }
      }
      // Get psys and issue decomposeDomainAll()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:decompose_domain_all;
            //-----------------------------------------------
         }
      }
   }
   //----------------------------
   // TreeForForce manipulators
   //----------------------------
   void create_tree(int *tree_num,
                    char *tree_info) {
      std::string tree_info_ = tree_info;
      Tree_Data tree_data;
      //------------------------------------------------------
      // fdps-autogen:create_tree;
      //------------------------------------------------------
      tree_vector.push_back(tree_data);
      *tree_num = tree_data.id;
      num_tree_creation++;
      num_tree++;
   }
   void delete_tree(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            it = tree_vector.erase(it);
            num_tree--;
            break;
         }
      }
   }
   void init_tree(const int tree_num,
                  const int numPtcl,
                  const float theta, 
                  const int n_leaf_limit,
                  const int n_group_limit) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:init_tree;
            //------------------------------------------------------
         }
      }
   }
   void get_tree_info(const int tree_num,
                      char *tree_info,
                      size_t *charlen) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            std::strcpy(tree_info,it->info.c_str());
            *charlen  = std::strlen(tree_info);
            break;
         }
      }
   }
   long long int get_tree_memsize(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_tree_memsize;
            //------------------------------------------------------
         }
      }
   }
   void get_tree_time_prof(const int tree_num,
                           PS::TimeProfile *prof) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_tree_time_prof;
            //------------------------------------------------------
         }
      }
   }
   void clear_tree_time_prof(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:clear_tree_time_prof;
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_interact_ep_ep_loc(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_interact_ep_ep_loc;
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_interact_ep_sp_loc(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_interact_ep_sp_loc;
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_interact_ep_ep_glb(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_interact_ep_ep_glb;
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_interact_ep_sp_glb(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_interact_ep_sp_glb;
            //------------------------------------------------------
         }
      }
   }
   void clear_num_interact(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:clear_num_interact;
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_tree_walk_loc(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_tree_walk_loc;
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_tree_walk_glb(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_tree_walk_glb;
            //------------------------------------------------------
         }
      }
   }
   //-----------------------------------------------------------------------------------
   // fdps-autogen:calc_force_all_and_write_back;
   //-----------------------------------------------------------------------------------
   // fdps-autogen:calc_force_all;
   //-----------------------------------------------------------------------------------
   // fdps-autogen:calc_force_making_tree;
   //-----------------------------------------------------------------------------------
   // fdps-autogen:calc_force_and_write_back;
   //-----------------------------------------------------------------------------------
   // fdps-autogen:get_neighbor_list;
   //-----------------------------------------------------------------------------------
   // fdps-autogen:get_epj_from_id;

   //----------------------------
   // Utility functions
   //----------------------------
   void mt_init_genrand(const int s) {
      unsigned long seed;
      seed = (unsigned long) s;
      PS::MT::init_genrand(seed);
   }
   int mt_genrand_int31(void) {
      return (int) PS::MT::genrand_int31();
   }
   double mt_genrand_real1(void){
      return PS::MT::genrand_real1();
   } 
   double mt_genrand_real2(void) {
      return PS::MT::genrand_real2(); 
   } 
   double mt_genrand_real3(void) {
      return PS::MT::genrand_real3();
   }
   double mt_genrand_res53(void) {
      return PS::MT::genrand_res53();
   }

   //----------------------------
   // Particle Mesh
   //----------------------------
#ifdef PARTICLE_SIMULATOR_USE_PM_MODULE
   void create_pm(int *pm_num) {
      PM_Data pm_data;
      PS::PM::ParticleMesh *pm;
      pm = new PS::PM::ParticleMesh;
      pm_data.id = num_pm_creation;
      pm_data.ptr = (void *) pm;
      pm_vector.push_back(pm_data);
      *pm_num = pm_data.id;
      num_pm_creation++;
      num_pm++;
   }
   void delete_pm(const int pm_num) {
      for (std::vector<PM_Data>::iterator it = pm_vector.begin(); it != pm_vector.end(); ++it) {
         if (it->id == pm_num) {
            it = pm_vector.erase(it);
            num_pm--;
            break;
         }
      }
   }
   int get_pm_mesh_num() {
      return SIZE_OF_MESH;
   }
   double get_pm_cutoff_radius() {
      return PS::PM::CUTOFF_RADIUS;
   }
   void set_dinfo_of_pm(const int pm_num,
                        const int dinfo_num) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
            break;
         }
      }
      // Get pm
      PS::PM::ParticleMesh *pm;
      for (std::vector<PM_Data>::iterator it = pm_vector.begin(); it != pm_vector.end(); ++it) {
         if (it->id == pm_num) {
            pm = (PS::PM::ParticleMesh *) it->ptr;
            break;
         }
      }
      // Issue setDomainInfoParticleMesh
      pm->setDomainInfoParticleMesh(*dinfo); 
   }
   void set_psys_of_pm(const int pm_num,
                       const int psys_num,
                       const bool clear) {
      // Get pm
      PS::PM::ParticleMesh *pm;
      for (std::vector<PM_Data>::iterator it = pm_vector.begin(); it != pm_vector.end(); ++it) {
         if (it->id == pm_num) {
            pm = (PS::PM::ParticleMesh *) it->ptr;
            break;
         }
      }
      // Get psys and issue setParticleParticleMesh()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:set_psys_of_pm;
            //-----------------------------------------------
         }
      }
   }
   void get_pm_force(const int pm_num,
                     const PS::F32vec *pos,
                     PS::F32vec *force) {
      // Get pm
      PS::PM::ParticleMesh *pm;
      for (std::vector<PM_Data>::iterator it = pm_vector.begin(); it != pm_vector.end(); ++it) {
         if (it->id == pm_num) {
            pm = (PS::PM::ParticleMesh *) it->ptr;
            break;
         }
      }
      // Issue getForce()
      *force = pm->getForce(*pos);
   }
   void get_pm_potential(const int pm_num,
                         const PS::F32vec *pos,
                         PS::F32 *pot) {
      // Get pm
      PS::PM::ParticleMesh *pm;
      for (std::vector<PM_Data>::iterator it = pm_vector.begin(); it != pm_vector.end(); ++it) {
         if (it->id == pm_num) {
            pm = (PS::PM::ParticleMesh *) it->ptr;
            break;
         }
      }
      // Issue getForce()
      *pot = pm->getPotential(*pos);
   }
   void calc_pm_force_only(const int pm_num) {
      // Get pm
      PS::PM::ParticleMesh *pm;
      for (std::vector<PM_Data>::iterator it = pm_vector.begin(); it != pm_vector.end(); ++it) {
         if (it->id == pm_num) {
            pm = (PS::PM::ParticleMesh *) it->ptr;
            break;
         }
      }
      // Issue calcMeshForceOnly()
      pm->calcMeshForceOnly();
   }
   void calc_pm_force_all_and_write_back(const int pm_num,
                                         const int psys_num,
                                         const int dinfo_num) {
      // Get pm
      PS::PM::ParticleMesh *pm;
      for (std::vector<PM_Data>::iterator it = pm_vector.begin(); it != pm_vector.end(); ++it) {
         if (it->id == pm_num) {
            pm = (PS::PM::ParticleMesh *) it->ptr;
            break;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
            break;
         }
      }
      // Get psys and issue calcForceAllAndWriteBack()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:calc_pm_force_all_and_write_back;
            //-----------------------------------------------
         }
      }
   }
#endif

   //----------------------------
   // Other utility functions
   //----------------------------
   void print_errmsg(const std::string& errmsg,
                     const std::string& func_name) {
    
      std::cout << "*** PS_FTN_IF_ERROR ***" << std::endl; 
      std::cout << "message:  " << errmsg << std::endl;
      std::cout << "function: " << func_name << std::endl;
      std::cout << "file:      FDPS_Manipulators.cpp" << std::endl;

   }

   std::vector<std::string> split(const std::string& s,
                                  const std::string& delim) {
      // This function splits a string by a given delimiter.
      // [Refs.]
      //   1) http://qiita.com/_meki/items/4328c98964ea33b0db0d
      //   2) http://qiita.com/iseki-masaya/items/70b4ee6e0877d12dafa8 and
      std::vector<std::string> results;
      std::string::size_type pos = 0;
      while (pos != std::string::npos) {
         std::string::size_type pos_delim_head = s.find(delim,pos);
         if (pos_delim_head == std::string::npos) { // the end of the target string
            // Register characters between pos and npos as a string
            // return the result to the caller.
            results.push_back(s.substr(pos));
            return results;
         }
         std::string::size_type len = pos_delim_head - pos; 
         results.push_back(s.substr(pos,len));
         pos = pos_delim_head + delim.size();
      }
   }

    void check_split(void) {
       std::string input_string = "interact_type,Force,EPI,EPJ,mode";
       std::string delim = ",";
       std::vector<std::string> results = split(input_string,delim);
       for (std::vector<std::string>::iterator it = results.begin(); it != results.end(); ++it)
          std::cout << *it << delim;
       std::cout << std::endl;
    }


};
