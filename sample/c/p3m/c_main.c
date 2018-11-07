/* Standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
/* FDPS headers */
#include "FDPS_c_if.h"
/* User-defined headers */
#include "user_defined.h"

void setup_NaCl_crystal(int psys_num,
                        int dinfo_num,
                        Crystal_parameters params) {
   // Get # of MPI processes and Rank number
   int nprocs = fdps_get_num_procs();
   int myrank = fdps_get_rank();

   // Define the parameters
   int nptcl = params.nptcl_per_side
             * params.nptcl_per_side
             * params.nptcl_per_side;
   if (params.nptcl_per_side % 2 != 0) {
      if (myrank == 0) {
         fprintf(stderr,"nptcl_per_side is an invalid value: %d\n",
                 params.nptcl_per_side);
      }
      fflush(stderr);
      fdps_finalize();
      exit(EXIT_FAILURE);
   }

   // Compute nptcl_list
   int i,j,k;
   int nptcl_list[nprocs];
   for (i = 0; i < nprocs; i++) nptcl_list[i]=0;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   // MPI-parallel
   int nptcl_loc = nptcl/nprocs;
   int nptcl_rem = nptcl % nprocs;
   if (nptcl_rem != 0) {
      if ((myrank+1) <= nptcl_rem) nptcl_loc++;
   }
   int sendbuf[nprocs];
   for (i = 0; i < nprocs; i++) sendbuf[i] = nptcl_loc;
   MPI_Alltoall(sendbuf,1,MPI_INT,
                nptcl_list,1,MPI_INT,
                MPI_COMM_WORLD);
   int i_start = 0;
   if (myrank > 0) {
      int irank;
      for (irank = 0; irank < myrank; irank++)
         i_start += nptcl_list[irank];
   }
   int i_end = i_start + nptcl_loc - 1;
#else
   // Serial
   int nptcl_loc = nptcl;
   int i_start = 0;
   int i_end = i_start + nptcl_loc - 1;
#endif
   fdps_set_nptcl_loc(psys_num,nptcl_loc);

   // Make a particle distribution
   FP_nbody *ptcl = (FP_nbody *) fdps_get_psys_cptr(psys_num);
   int size_of_mesh = fdps_get_pm_mesh_num();
   double cutoff_radius = fdps_get_pm_cutoff_radius();
   int id = 0;
   for (i = 0; i < params.nptcl_per_side; i++) {
      for (j = 0; j < params.nptcl_per_side; j++) {
         for (k = 0; k < params.nptcl_per_side; k++) {
            // Substitute
            if ((i_start <= id) && (id <= i_end)) {
               int ii = id - i_start;
               ptcl[ii].id = id;
               if ((i+j+k) % 2 == 0) {
                  ptcl[ii].mass =  1.0
                                / (params.nptcl_per_side
                                  *params.nptcl_per_side);
               } else {
                  ptcl[ii].mass = -1.0
                                / (params.nptcl_per_side
                                  *params.nptcl_per_side);
               }
               ptcl[ii].rcut = cutoff_radius/size_of_mesh;
               ptcl[ii].pos.x = (params.pos_vertex.x + i) * (1.0/params.nptcl_per_side);
               ptcl[ii].pos.y = (params.pos_vertex.y + j) * (1.0/params.nptcl_per_side);
               ptcl[ii].pos.z = (params.pos_vertex.z + k) * (1.0/params.nptcl_per_side);
            }
            // Update id
            id++;
         }
      }
   }

   // Domain decomposition & exchange particle
   fdps_set_boundary_condition(dinfo_num,FDPS_BC_PERIODIC_XYZ);
   fdps_f32vec pos_ll, pos_ul;
   pos_ll.x = 0.0; pos_ll.y = 0.0; pos_ll.z = 0.0;
   pos_ul.x = 1.0; pos_ul.y = 1.0; pos_ul.z = 1.0;
   fdps_set_pos_root_domain(dinfo_num,&pos_ll,&pos_ul);
   fdps_decompose_domain_all(dinfo_num,psys_num,-1.0);
   fdps_exchange_particle(psys_num,dinfo_num);

}

void calc_energy_error(int psys_num, double *relerr) {
   const double e0=-1.7475645946332e0;
   // where e0 is the analytical solution for \sum q_{i}*\phi_{i}
   // obtained by the PM^{3} method (computed by K.Nitadori).
   int i, nptcl_loc = fdps_get_nptcl_loc(psys_num);
   FP_nbody *ptcl = (FP_nbody *) fdps_get_psys_cptr(psys_num);
   double ebuf=0.0;
   for (i = 0; i < nptcl_loc; i++) 
      ebuf += ptcl[i].mass * ptcl[i].pot;
   double esum = fdps_get_sum_f64(ebuf);
   *relerr = fabs((esum-e0)/e0);
}

//-----------------------------------------------------------------------
////////////////////////// M A I N   F U N C T I O N ////////////////////
//-----------------------------------------------------------------------
void c_main() {
   // Initialize FDPS 
   fdps_initialize();

   // Make an instance of particle system and initialize it
   int psys_num;
   fdps_create_psys(&psys_num,"fp_nbody");
   fdps_init_psys(psys_num);

   // Make an instance of domain info. and initialize it
   int dinfo_num;
   fdps_create_dinfo(&dinfo_num);
   const float coef_ema = 0.3;
   fdps_init_dinfo(dinfo_num,coef_ema);

   // Make an instance of ParticleMesh object
   int pm_num;
   fdps_create_pm(&pm_num);

   // Initialize Mersenne twister pseudo-random number generator
   int mtts_num;
   fdps_create_mtts(&mtts_num);
   fdps_mtts_init_genrand(mtts_num,0);

   //==================================================================
   // Compute relative energy errors of the Madelung energy
   // due to the P^{3}M method for different # of particles
   // and for different configurations.
   //==================================================================
   //-(Local variables for tree)
   int tree_num;
   bool is_tree_initialized = false;
   const float theta = 0.0;
   const int n_leaf_limit = 8, n_group_limit = 64;
   //-(Local variables for I/O)
   bool is_1st_output = true;
   int nptcl_1d;
   for (nptcl_1d = 4; nptcl_1d <= 32; nptcl_1d += 2) {
      // Information to STDOUT
      if (fdps_get_rank() == 0) 
         printf("Processing %d case...\n",nptcl_1d);

      Crystal_parameters NaCl_params;
      NaCl_params.nptcl_per_side = nptcl_1d;
      double relerr = 0.0;
      const int num_trials = 512;
      int nstep;
      for (nstep = 0; nstep < num_trials; nstep++) {
         // [1] Randomly choose a configuration of the grid
         NaCl_params.pos_vertex.x = fdps_mtts_genrand_res53(mtts_num);
         NaCl_params.pos_vertex.y = fdps_mtts_genrand_res53(mtts_num);
         NaCl_params.pos_vertex.z = fdps_mtts_genrand_res53(mtts_num);

         // [2] Make a NaCl crystal
         setup_NaCl_crystal(psys_num, dinfo_num, NaCl_params);

         // [3] Initialize tree if needed
         if (is_tree_initialized == false) {
            int nptcl_loc = fdps_get_nptcl_loc(psys_num);
            fdps_create_tree(&tree_num,
                             "Long,force_pp,ep_nbody,ep_nbody,MonopoleWithCutoff");
            fdps_init_tree(tree_num,3*nptcl_loc,theta,
                           n_leaf_limit,n_group_limit);
            is_tree_initialized = true;
         }

         // [4] Compute force and potential with P^{3}M method
         // [4-1] Get the pointer to FP and # of local particles
         int nptcl_loc = fdps_get_nptcl_loc(psys_num);
         FP_nbody *ptcl = (FP_nbody *) fdps_get_psys_cptr(psys_num);
         // [4-2] PP part
         fdps_calc_force_all_and_write_back(tree_num,
                                            calc_force_ep_ep,
                                            calc_force_ep_sp,
                                            psys_num,
                                            dinfo_num,
                                            true,
                                            FDPS_MAKE_LIST);
         // [4-3] PM part
         fdps_calc_pm_force_all_and_write_back(pm_num,
                                               psys_num,
                                               dinfo_num);
         int i;
         for (i = 0; i < nptcl_loc; i++) {
            fdps_f32vec pos32;
            pos32.x = ptcl[i].pos.x;
            pos32.y = ptcl[i].pos.y;
            pos32.z = ptcl[i].pos.z;
            fdps_get_pm_potential(pm_num,&pos32,&ptcl[i].pot_pm);
         }
         // [4-4] Compute the total acceleration and potential
         for (i = 0; i < nptcl_loc; i++) {
            ptcl[i].pot -= ptcl[i].pot_pm;
            ptcl[i].acc.x -= ptcl[i].acc_pm.x;
            ptcl[i].acc.y -= ptcl[i].acc_pm.y;
            ptcl[i].acc.z -= ptcl[i].acc_pm.z;
         }

         // [5] Compare with the result of the Ewald summation
         double relerr1;
         calc_energy_error(psys_num,&relerr1);
         relerr += relerr1;

      }
      relerr /= num_trials;

      // Output relative error
      if (fdps_get_rank() == 0) {
         // Output to file
         char file_name[64];
         strcpy(file_name,"EnergyError.dat");
         FILE *fp;
         if (is_1st_output == true) {
            if ((fp = fopen(file_name,"w")) == NULL) {
                fprintf(stderr,"Cannot open file %s.\n",file_name);
                exit(EXIT_FAILURE);
            }
            is_1st_output = false;
         } else {
            if ((fp = fopen(file_name,"a")) == NULL) {
                fprintf(stderr,"Cannot open file %s.\n",file_name);
                exit(EXIT_FAILURE);
            }
         }
         fprintf(fp,"%d %15.7e\n",nptcl_1d,relerr);
         fclose(fp);
         // Output to STDOUT
         printf("********** Result of this experiment **********\n");
         printf("   num_trials     = %d\n",num_trials);
         printf("   nptcl_per_side = %d\n",nptcl_1d);
         printf("   Relative Error = %15.7e\n",relerr);
         printf("***********************************************\n");
      }
   }

   // Finalize FDPS
   fdps_finalize();

}

