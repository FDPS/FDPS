/* C++ standard heasers */
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cfloat>
#include <vector>
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include <omp.h>
#endif
#include <random>
/* FDPS headers */
#include <particle_simulator.hpp>
#include <particle_mesh.hpp>
#include <param.h>
#include <param_fdps.h>

/* Cutoff functions */
PS::F64 S2_pcut(const PS::F64 xi) {
   // This is the potential cutoff function where we used Eq.(8.75)
   // in Hockney & Eastwood (1987).

   if (xi <= 1.0) {
      return 1.0 - xi*(208.0 
                      +(xi*xi)*(-112.0 
                               +(xi*xi)*(56.0 
                                        +xi*(-14.0 
                                            +xi*(-8.0
                                                +3.0*xi)))))/140.0;
   } else if ((1.0 < xi) && (xi < 2.0)) {
      return 1.0 - (12.0 
                   +xi*(128.0
                       +xi*(224.0 
                           +xi*(-448.0 
                               +xi*(280.0 
                                   +xi*(-56.0 
                                       +xi*(-14.0 
                                           +xi*(8.0
                                               -xi))))))))/140.0;
   } else {
      return 0.0;
   }
}

PS::F64 S2_fcut(const PS::F64 xi) {
   // This function returns 1 - R(\xi), where \xi is r/(a/2), a is the
   // scale length of the cutoff function, and R(\xi) is almost the same
   // as the function defined as Eq.(8-72) in Hockney & Eastwood (1987).
   // The only difference is that [1/(r/(a/2))]^2 is factored out 
   // in this function from Eq.(8-72).

   if (xi <= 1.0) {
      return 1.0 - (xi*xi*xi)*(224.0 
                              +(xi*xi)*(-224.0
                                       +xi*(70.0
                                           +xi*(48.0-21.0*xi))))/140.0;
   } else if ((1.0 < xi) && (xi < 2.0)) {
      return 1.0 - (12.0 
                   +(xi*xi)*(-224.0 
                            +xi*(896.0 
                                +xi*(-840.0 
                                    +xi*(224.0 
                                        +xi*(70.0 
                                            +xi*(-48.0+7.0*xi)))))))/140.0;
   } else {
      return 0.0;
   }
}

//#####################
/* Class definitions */
//#####################
//** Force Class
class Nbody_PP_Results
{
   public:
      PS::F64 pot;
      PS::F64vec agrv;
      void clear() {
         pot  = 0.0;
         agrv = 0.0;
      }
};

//** Full Particle Class
class Nbody_FP
{
   public:
      PS::S64 id;
      PS::F64 m;
      PS::F64 rc;
      PS::F64vec x; 
      PS::F64vec v,v_half; 
      PS::F64vec agrv; 
      PS::F64 pot;
      // Member functions required by FDPS
      PS::F64 getCharge() const {
         return m;
      };
      PS::F64 getChargeParticleMesh() const {
         return m;
      };
      PS::F64vec getPos() const {
         return x;
      };
      PS::F64 getRSearch() const {
         return rc;
      };
      void setPos(const PS::F64vec& x) {
         this->x = x;
      };
      void copyFromForce(const Nbody_PP_Results& result) {};
      void copyFromForceParticleMesh(const PS::F64 apm) {};
};

//** Essential Particle Class
class Nbody_EP
{
   public:
      PS::S64 id;
      PS::F64 m;
      PS::F64 rc;
      PS::F64vec x;
      // Member functions required by FDPS
      PS::F64 getCharge() const {
         return m;
      };
      PS::F64vec getPos() const {
         return x;
      };
      PS::F64 getRSearch() const {
         return rc;
      };
      void setPos(const PS::F64vec& x) {
         this->x = x;
      };
      void copyFromFP(const Nbody_FP& FP) {
         id  = FP.id;
         m   = FP.m;
         rc  = FP.rc;
         x   = FP.x;
      };
};

//** Force Functors
class Calc_force_ep_ep{
   public:
      void operator () (const Nbody_EP* const ep_i,
                        const PS::S32 Nip,
                        const Nbody_EP* const ep_j,
                        const PS::S32 Njp,
                        Nbody_PP_Results* const result) {
         for (PS::S32 i=0; i<Nip; i++) {
            for (PS::S32 j=0; j<Njp; j++) {
               PS::F64vec dx = ep_i[i].x - ep_j[j].x;
               PS::F64 rij = std::sqrt(dx * dx);
               if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0)) continue;
               PS::F64 rinv = 1.0/rij;
               PS::F64 rinv3 = rinv*rinv*rinv;
               PS::F64 xi = 2.0*rij/ep_i[i].rc;
               result[i].pot  += ep_j[j].m * S2_pcut(xi) * rinv;
               result[i].agrv += ep_j[j].m * S2_fcut(xi) * rinv3 * dx;
            }
            //* Self-interaction term
            result[i].pot -= ep_i[i].m * (208.0/(70.0*ep_i[i].rc));
         }

      };
};

class Calc_force_ep_sp{
   public:
      void operator () (const Nbody_EP* const ep_i,
                        const PS::S32 Nip,
                        const PS::SPJMonopoleCutoff* const ep_j,
                        const PS::S32 Njp,
                        Nbody_PP_Results* const result) {
         for (PS::S32 i=0; i<Nip; i++) {
            for (PS::S32 j=0; j<Njp; j++) {
               PS::F64vec dx = ep_i[i].x - ep_j[j].pos;
               PS::F64 rij = std::sqrt(dx * dx);
               PS::F64 rinv = 1.0/rij;
               PS::F64 rinv3 = rinv*rinv*rinv;
               PS::F64 xi = 2.0*rij/ep_i[i].rc;
               result[i].pot  += ep_j[j].mass * S2_pcut(xi) * rinv;
               result[i].agrv += ep_j[j].mass * S2_fcut(xi) * rinv3 * dx;
            }
         }

      };
};


//*  Nbody_Objects Class
class Nbody_Objects {
   public:
      PS::ParticleSystem<Nbody_FP> system;
      PS::DomainInfo dinfo;
      PS::TreeForForceLong<Nbody_PP_Results, Nbody_EP, Nbody_EP>::MonopoleWithCutoff pp_tree;
      PS::PM::ParticleMesh pm;
      // Methods for Nbody simulation
      void init_tree() {
         PS::S32 numPtclLocal = system.getNumberOfParticleLocal();
         PS::U64 ntot = 3 * numPtclLocal;
         pp_tree.initialize(ntot,0.0);
      };
      void calc_gravity() {
         //* Local variables
         PS::S32 numPtclLocal = system.getNumberOfParticleLocal();

         //* Reset potential and accelerations
         for (PS::S32 i=0; i<numPtclLocal; i++) {
            system[i].pot  = 0.0;
            system[i].agrv = 0.0;
         }

         //=================
         //* PM part 
         //=================
         pm.setDomainInfoParticleMesh(dinfo);
         pm.setParticleParticleMesh(system);
         pm.calcMeshForceOnly();
         for (PS::S32 i=0; i<numPtclLocal; i++) { 
            PS::F32vec x32 = system[i].x; 
            system[i].pot  -= pm.getPotential(x32);
            system[i].agrv -= pm.getForce(x32);
         }

         //=================
         //* PP part 
         //=================
         pp_tree.calcForceAll(Calc_force_ep_ep(),
                              Calc_force_ep_sp(),
                              system, dinfo);
         for (PS::S32 i=0; i<numPtclLocal; i++) {
            Nbody_PP_Results result = pp_tree.getForce(i);
            system[i].pot  += result.pot;
            system[i].agrv += result.agrv;
         }

      };
};

//** Crystal Parameters Class
class Crystal_Parameters
{
   public:
      PS::S32 numPtcl_per_side;
      PS::F64vec pos_vertex;
};

//* Initial condition
void NaCl_IC(PS::ParticleSystem<Nbody_FP>& system,
             PS::DomainInfo& dinfo,
             Crystal_Parameters& params) {
   //* Local variables
   PS::S32 nprocs,myrank;
   //-(lattice)
   PS::F64 xmin,xmax,ymin,ymax,zmin,zmax;
   PS::F64 m,x,y,z;
   PS::S32 numPtcl,numPtclLocal,numPtclRem;
   PS::S32 i_start,i_end;
   //-(IO)
   std::string filename;
   std::ofstream output_file;
   //-(Time measurement)
   struct timeval start_time,end_time;
   double e_time;

   //* Get # of MPI processes and Rank number
   nprocs = PS::Comm::getNumberOfProc();
   myrank = PS::Comm::getRank();

   //* Define the parameters
   numPtcl = params.numPtcl_per_side
           * params.numPtcl_per_side
           * params.numPtcl_per_side;
   if (params.numPtcl_per_side % 2 != 0) {
      std::cout << "numPtcl_per_side is an invalid value: "
                << params.numPtcl_per_side 
                << std::endl;
      PS::Finalize();
      std::exit(1);
   }

   //* Compute numPtclList
   std::vector<PS::S32> numPtclList(nprocs);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   //* MPI-parallel
   numPtclLocal = numPtcl/nprocs;
   if ((numPtclRem = numPtcl % nprocs) != 0) {
      if ((myrank+1) <= numPtclRem) {
         numPtclLocal++;
      }
   }
   std::vector<PS::S32> sendbuf(nprocs);
   for (PS::S32 i=0; i<nprocs; i++) sendbuf[i]=numPtclLocal;
   MPI::COMM_WORLD.Alltoall(&sendbuf[0],1,MPI::INT, 
                            &numPtclList[0],1,MPI::INT);
   i_start = 0;
   if (myrank > 0) {
      for (PS::S32 irank=0; irank<myrank; irank++)
         i_start += numPtclList[irank];
   }
   i_end = i_start + numPtclLocal - 1;
#else
   numPtclList.push_back(numPtcl);
   numPtclLocal = numPtcl;
   i_start = 0;
   i_end   = i_start + numPtclLocal - 1;
#endif
   system.setNumberOfParticleLocal(numPtclLocal);

   //* Make a particle distribution
   PS::S32 id {0};
   for(int i=0; i<params.numPtcl_per_side; i++) {
      for(int j=0; j<params.numPtcl_per_side; j++) {
         for(int k=0; k<params.numPtcl_per_side; k++){
            //* Substitute
            if ((i_start <= id) && (id <= i_end)) {
               PS::S32 ii = id - i_start;
               system[ii].id  = id;
               system[ii].m   = ((i+j+k)%2 ? 1.0 : -1.0)
                              / (params.numPtcl_per_side
                                *params.numPtcl_per_side);
               system[ii].rc  = PS::ParticleMesh::CUTOFF_RADIUS / SIZE_OF_MESH; // see param.h & param_fdps.h
               system[ii].x.x = (params.pos_vertex.x + i) * (1.0/params.numPtcl_per_side);
               system[ii].x.y = (params.pos_vertex.y + j) * (1.0/params.numPtcl_per_side);
               system[ii].x.z = (params.pos_vertex.z + k) * (1.0/params.numPtcl_per_side);
               system[ii].v   = 0.0;
            } // END of i_start & i_end

            //* Update id
            id++;
         }
      }
   }

   //* Initialize domain info and apply periodic boundary condition
   dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
   dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
                          PS::F64vec(1.0,1.0,1.0));
   dinfo.decomposeDomainAll(system);
   //* Exchange particle
   system.exchangeParticle(dinfo);
   //* Check load balance
   //numPtclLocal = system.getNumberOfParticleLocal();
   //{
   //   PS::S32 Min = PS::Comm::getMinValue(numPtclLocal);
   //   PS::S32 Max = PS::Comm::getMaxValue(numPtclLocal);
   //   PS::S32 Sum = PS::Comm::getSum(numPtclLocal);
   //   if (myrank == 0) {
   //      std::cout << "Max. " << Max << std::endl;
   //      std::cout << "Min. " << Min << std::endl;
   //      std::cout << "Sum. " << Sum << std::endl;
   //   }
   //}
  
   //* Check exchanged particle distribution
//#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
//   std::ostringstream ss;
//   ss << "IC" << std::setfill('0') << std::setw(5) << myrank << ".txt";
//   filename = ss.str();
//#else
//   filename = "IC.txt";
//#endif
//   output_file.open(filename.c_str(),std::ios::trunc);
//   output_file.setf(std::ios_base::scientific,
//                    std::ios_base::floatfield);
//   for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
//      output_file << std::setprecision(15)
//                  << std::showpos
//                  << system[i].x.x << " "
//                  << system[i].x.y << " "
//                  << system[i].x.z << " "
//                  << std::endl;
//   }
//   output_file.close();

   //if (myrank == 0) 
   //   std::cout << "complete to make I.C.!" << std::endl;

}

PS::F64 calc_energy_error(PS::ParticleSystem<Nbody_FP>& system,
                          PS::S32 nstep){

   //* Analytical solution for \sum q_{i}*\phi_{i}
   PS::F64 E0 {-1.7475645946332}; // Ewald method

   //* Get MPI parallel info.
   PS::S32 nprocs = PS::Comm::getNumberOfProc();
   PS::S32 myrank = PS::Comm::getRank();

   //* Particle num.
   PS::S32 numPtclLocal  = system.getNumberOfParticleLocal();

   //* Compute total energy
   PS::F64 Ebuf {0.0}, Esum {0.0};
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      Ebuf += system[i].m * system[i].pot;
   }
   MPI::COMM_WORLD.Allreduce(&Ebuf,&Esum,1,MPI::DOUBLE,MPI::SUM);

   return std::abs((Esum-E0)/E0);
   
}

/*-------------------------------------------------------------------*/
////////////////////// M A I N   F U N C T I O N //////////////////////
/*-------------------------------------------------------------------*/
int main(int argc, char* argv[]) {

   //* Local parameters
   const PS::S32 numTrials=512;

   //* Local variables
   Nbody_Objects Nbody_objs;
   Crystal_Parameters NaCl_params;
   static bool is_tree_initialized=false;

   //* Initialize FDPS 
   PS::Initialize(argc,argv);

   //* Initialize FDPS objects
   Nbody_objs.system.initialize();
   Nbody_objs.dinfo.initialize();

   //* Initialize Mersenne twister pseudo-random number generator
   PS::MT::init_genrand(0);

   //==================================================================
   //* Compute relative energy errors of the Madelung energy
   //  due to the P^{3}M method for different # of particles
   //  and for different configurations.
   //==================================================================
   for (PS::S32 numPtcl1D=4; numPtcl1D <= 32; numPtcl1D += 2) {
      if (PS::Comm::getRank() == 0)
         std::cout << "Processing " << numPtcl1D << " case..." << std::endl;
      NaCl_params.numPtcl_per_side = numPtcl1D;
      PS::F64 relerr {0.0};
      for (PS::S32 nstep=1; nstep <= numTrials; nstep++) {
         //if (PS::Comm::getRank() == 0)
         //   std::cout << "nstep = " << nstep << std::endl;
         //PS::Comm::barrier();

         //* [1] Randomly choose a configuration of the grid
         NaCl_params.pos_vertex.x = PS::MT::genrand_res53();
         NaCl_params.pos_vertex.y = PS::MT::genrand_res53();
         NaCl_params.pos_vertex.z = PS::MT::genrand_res53();

         //* [2] Make a NaCl crystal
         NaCl_IC(Nbody_objs.system,
                 Nbody_objs.dinfo,
                 NaCl_params);

         //* [3] Initialize tree if needed
         if (is_tree_initialized == false) {
            Nbody_objs.init_tree();
            is_tree_initialized = true;
         }

         //* [4] Compute force and potential with P^{3}M method
         Nbody_objs.calc_gravity();

         //* [5] Compare with the result of the Ewald summation
         relerr += calc_energy_error(Nbody_objs.system,nstep);

      }
      relerr /= numTrials;

      //* Output relative error
      if (PS::Comm::getRank() == 0) {
         //* Output to file
         std::string filename {"EnergyError.dat"};
         std::ofstream output_file;
         output_file.open(filename.c_str(),std::ios::out | std::ios::app);
         output_file << std::setprecision(15)
                     << numPtcl1D << " "
                     << relerr << std::endl;
         output_file.close();
         //* Output to STDOUT
         std::cout << "********** Result of this experiment **********" << std::endl;
         std::cout << "   numTrials        = " << numTrials << std::endl;
         std::cout << "   numPtcl_per_side = " << numPtcl1D << std::endl;
         std::cout << "   Relative Error   = " << relerr    << std::endl;
         std::cout << "***********************************************" << std::endl;
      }

   } 

   //* Finalize FDPS
   PS::Finalize();
   return 0;
}
