// Include FDPS header
#include <particle_simulator.hpp>
// Include the standard C++ headers
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <sys/stat.h>

/* Parameters */
const short int Dim = 3;
const PS::F64 SMTH = 1.2;
const PS::U32 OUTPUT_INTERVAL = 10;
const PS::F64 C_CFL = 0.3;

/* Kernel Function */
const PS::F64 pi = atan(1.0) * 4.0;
const PS::F64 kernelSupportRadius = 2.5;

PS::F64 W(const PS::F64vec dr, const PS::F64 h){
   const PS::F64 H = kernelSupportRadius * h;
   const PS::F64 s = sqrt(dr * dr) / H;
   const PS::F64 s1 = (1.0 - s < 0) ? 0 : 1.0 - s;
   const PS::F64 s2 = (0.5 - s < 0) ? 0 : 0.5 - s;
   PS::F64 r_value = pow(s1, 3) - 4.0 * pow(s2, 3);
   //if # of dimension == 3
   r_value *= 16.0 / pi / (H * H * H);
   return r_value;
}

PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h){
   const PS::F64 H = kernelSupportRadius * h;
   const PS::F64 s = sqrt(dr * dr) / H;
   const PS::F64 s1 = (1.0 - s < 0) ? 0 : 1.0 - s;
   const PS::F64 s2 = (0.5 - s < 0) ? 0 : 0.5 - s;
   PS::F64 r_value = - 3.0 * pow(s1, 2) + 12.0 * pow(s2, 2);
   //if # of dimension == 3
   r_value *= 16.0 / pi / (H * H * H);
   return dr * r_value / (sqrt(dr * dr) * H + 1.0e-6 * h);
}

/* Class Definitions */
//** Force Class (Result Class)
class Dens{
   public:
   PS::F64 dens;
   PS::F64 smth;
   void clear(){
      dens = 0;
   }
};
class Hydro{
   public:
   PS::F64vec acc;
   PS::F64 eng_dot;
   PS::F64 dt;
   void clear(){
      acc = 0;
      eng_dot = 0;
   }
};

//** Full Particle Class
struct FP{
   PS::F64 mass;
   PS::F64vec pos;
   PS::F64vec vel;
   PS::F64vec acc;
   PS::F64 dens;
   PS::F64 eng;
   PS::F64 pres;
   PS::F64 smth;
   PS::F64 snds;
   PS::F64 eng_dot;
   PS::F64 dt;
   PS::S64 id;
   PS::F64vec vel_half;
   PS::F64 eng_half;
   void copyFromForce(const Dens& dens){
      this->dens = dens.dens;
   }
   void copyFromForce(const Hydro& force){
      this->acc     = force.acc;
      this->eng_dot = force.eng_dot;
      this->dt      = force.dt;
   }
   PS::F64 getCharge() const{
      return this->mass;
   }
   PS::F64vec getPos() const{
      return this->pos;
   }
   PS::F64 getRSearch() const{
      return kernelSupportRadius * this->smth;
   }
   void setPos(const PS::F64vec& pos){
      this->pos = pos;
   }
   void writeAscii(FILE* fp) const{
      fprintf(fp,
              "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t"
              "%lf\t%lf\t%lf\t%lf\t%lf\n",
              this->id, this->mass,
              this->pos.x, this->pos.y, this->pos.z,
              this->vel.x, this->vel.y, this->vel.z,
              this->dens, this->eng, this->pres);
   }
   void readAscii(FILE* fp){
      fscanf(fp,
             "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t"
             "%lf\t%lf\t%lf\t%lf\t%lf\n",
             &this->id, &this->mass,
             &this->pos.x, &this->pos.y, &this->pos.z,
             &this->vel.x, &this->vel.y, &this->vel.z,
             &this->dens, &this->eng, &this->pres);
   }
   void setPressure(){
      const PS::F64 hcr = 1.4;
      pres = (hcr - 1.0) * dens * eng;
      snds = sqrt(hcr * pres / dens);
   }
};

//** Essential Particle Class
struct EP{
   PS::F64vec pos;
   PS::F64vec vel;
   PS::F64    mass;
   PS::F64    smth;
   PS::F64    dens;
   PS::F64    pres;
   PS::F64    snds;
   void copyFromFP(const FP& rp){
      this->pos  = rp.pos;
      this->vel  = rp.vel;
      this->mass = rp.mass;
      this->smth = rp.smth;
      this->dens = rp.dens;
      this->pres = rp.pres;
      this->snds = rp.snds;
   }
   PS::F64vec getPos() const{
      return this->pos;
   }
   PS::F64 getRSearch() const{
      return kernelSupportRadius * this->smth;
   }
   void setPos(const PS::F64vec& pos){
      this->pos = pos;
   }
};

class FileHeader{
   public:
   PS::S32 Nbody;
   PS::F64 time;
   int readAscii(FILE* fp){
      fscanf(fp, "%lf\n", &time);
      fscanf(fp, "%d\n", &Nbody);
      return Nbody;
   }
   void writeAscii(FILE* fp) const{
      fprintf(fp, "%e\n", time);
      fprintf(fp, "%d\n", Nbody);
   }
};

struct boundary{
   PS::F64 x, y, z;
};


/* Force Functors */
class CalcDensity{
   public:
   void operator () (const EP* const ep_i, const PS::S32 Nip,
                     const EP* const ep_j, const PS::S32 Njp,
                     Dens* const dens){
      for(PS::S32 i = 0 ; i < Nip ; ++i){
         dens[i].clear();
         for(PS::S32 j = 0 ; j < Njp ; ++j){
            const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
            dens[i].dens += ep_j[j].mass * W(dr, ep_i[i].smth);
         }
      }
   }
};

class CalcHydroForce{
   public:
   void operator () (const EP* const ep_i, const PS::S32 Nip,
                     const EP* const ep_j, const PS::S32 Njp,
                     Hydro* const hydro){
      for(PS::S32 i = 0; i < Nip ; ++ i){
         hydro[i].clear();
         PS::F64 v_sig_max = 0.0;
         for(PS::S32 j = 0; j < Njp ; ++j){
            const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
            const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
            const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
            const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;
            v_sig_max = std::max(v_sig_max, v_sig);
            const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens + ep_j[j].dens));
            const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr, ep_j[j].smth));
            hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) * gradW_ij;
            hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;
         }
         hydro[i].dt = C_CFL * 2.0 * ep_i[i].smth / v_sig_max;
      }
   }
};

void makeOutputDirectory(char * dir_name) {
    struct stat st;
    PS::S32 ret;
    if (PS::Comm::getRank() == 0) {
        if (stat(dir_name, &st) != 0) {
            ret = mkdir(dir_name, 0777);
        } else {
            ret = 0; // the directory named dir_name already exists.
        }
    } 
    PS::Comm::broadcast(&ret, 1);
    if (ret == 0) {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
        PS::Abort();
    }
}

void SetupIC(PS::ParticleSystem<FP>& sph_system, PS::F64 *end_time, boundary *box){
   // Place SPH particles
   std::vector<FP> ptcl;
   const PS::F64 dx = 1.0 / 128.0;
   box->x = 1.0;
   box->y = box->z = box->x / 8.0;
   PS::S32 i = 0;
   for(PS::F64 x = 0 ; x < box->x * 0.5 ; x += dx){
      for(PS::F64 y = 0 ; y < box->y ; y += dx){
         for(PS::F64 z = 0 ; z < box->z ; z += dx){
            FP ith;
            ith.pos.x = x;
            ith.pos.y = y;
            ith.pos.z = z;
            ith.dens = 1.0;
            ith.mass = 0.75;
            ith.eng  = 2.5;
            ith.id   = i++;
            ith.smth = 0.012;
            ptcl.push_back(ith);
         }
      }
   }
   for(PS::F64 x = box->x * 0.5 ; x < box->x * 1.0 ; x += dx * 2.0){
      for(PS::F64 y = 0 ; y < box->y ; y += dx){
         for(PS::F64 z = 0 ; z < box->z ; z += dx){
            FP ith;
            ith.pos.x = x;
            ith.pos.y = y;
            ith.pos.z = z;
            ith.dens = 0.5;
            ith.mass = 0.75;
            ith.eng  = 2.5;
            ith.id   = i++;
            ith.smth = 0.012;
            ptcl.push_back(ith);
         }
      }
   }
   for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
      ptcl[i].mass = ptcl[i].mass * box->x * box->y * box->z / (PS::F64)(ptcl.size());
   }
   std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
   // Scatter SPH particles
   assert(ptcl.size() % PS::Comm::getNumberOfProc() == 0);
   const PS::S32 numPtclLocal = ptcl.size() / PS::Comm::getNumberOfProc();
   sph_system.setNumberOfParticleLocal(numPtclLocal);
   const PS::U32 i_head = numPtclLocal * PS::Comm::getRank();
   const PS::U32 i_tail = numPtclLocal * (PS::Comm::getRank() + 1);
   for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
      if(i_head <= i && i < i_tail){
         const PS::U32 ii = i - numPtclLocal * PS::Comm::getRank();
         sph_system[ii] = ptcl[i];
      }
   }
   // Set the end time
   *end_time = 0.12;
   // Fin.
   std::cout << "setup..." << std::endl;
}

void Initialize(PS::ParticleSystem<FP>& sph_system){
   for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
      sph_system[i].setPressure();
   }
}

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<FP>& sph_system){
   PS::F64 dt = 1.0e+30; //set VERY LARGE VALUE
   for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
      dt = std::min(dt, sph_system[i].dt);
   }
   return PS::Comm::getMinValue(dt);
}

void InitialKick(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
   for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
      sph_system[i].vel_half = sph_system[i].vel + 0.5 * dt * sph_system[i].acc;
      sph_system[i].eng_half = sph_system[i].eng + 0.5 * dt * sph_system[i].eng_dot;
   }
}

void FullDrift(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
   // time becomes t + dt;
   for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
      sph_system[i].pos += dt * sph_system[i].vel_half;
   }
}

void Predict(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
   for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
      sph_system[i].vel += dt * sph_system[i].acc;
      sph_system[i].eng += dt * sph_system[i].eng_dot;
   }
}

void FinalKick(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
   for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
      sph_system[i].vel = sph_system[i].vel_half + 0.5 * dt * sph_system[i].acc;
      sph_system[i].eng = sph_system[i].eng_half + 0.5 * dt * sph_system[i].eng_dot;
   }
}

void setPressure(PS::ParticleSystem<FP>& sph_system){
   for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
      sph_system[i].setPressure();
   }
}

void CheckConservativeVariables(const PS::ParticleSystem<FP>& sph_system){
   PS::F64vec Mom=0.0; // total momentum
   PS::F64    Eng=0.0; // total enegry
   for(PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++ i){
      Mom += sph_system[i].vel * sph_system[i].mass;
      Eng += (sph_system[i].eng + 0.5 * sph_system[i].vel * sph_system[i].vel)
            * sph_system[i].mass;
   }
   Eng = PS::Comm::getSum(Eng);
   Mom = PS::Comm::getSum(Mom);
    if(PS::Comm::getRank() == 0){
        printf("%.16e\n", Eng);
        printf("%.16e\n", Mom.x);
        printf("%.16e\n", Mom.y);
        printf("%.16e\n", Mom.z);
    }
}

int main(int argc, char* argv[]){
   // Initialize FDPS
   PS::Initialize(argc, argv);
   // Make a directory
   char dir_name[1024];
   sprintf(dir_name,"./result");
   makeOutputDirectory(dir_name);
   // Display # of MPI processes and threads
   PS::S32 nprocs = PS::Comm::getNumberOfProc();
   PS::S32 nthrds = PS::Comm::getNumberOfThread();
   std::cout << "===========================================" << std::endl
             << " This is a sample program of "               << std::endl
             << " Smoothed Particle Hydrodynamics on FDPS!"   << std::endl
             << " # of processes is " << nprocs               << std::endl
             << " # of thread is    " << nthrds               << std::endl
             << "===========================================" << std::endl;
   // Make an instance of ParticleSystem and initialize it
   PS::ParticleSystem<FP> sph_system;
   sph_system.initialize();
   // Define local variables
   PS::F64 dt, end_time;
   boundary box;
   // Make an initial condition and initialize the particle system
   SetupIC(sph_system, &end_time, &box);
   Initialize(sph_system);
   // Make an instance of DomainInfo and initialize it
   PS::DomainInfo dinfo;
   dinfo.initialize();
   // Set the boundary condition
   dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
   dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0),
                          PS::F64vec(box.x, box.y, box.z));
   // Perform domain decomposition
   dinfo.decomposeDomainAll(sph_system);
   // Exchange the SPH particles between the (MPI) processes
   sph_system.exchangeParticle(dinfo);
   // Make two tree structures
   // (one is for the density calculation and
   //  another is for the force calculation.)
   PS::TreeForForceShort<Dens, EP, EP>::Gather dens_tree;
   dens_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());

   PS::TreeForForceShort<Hydro, EP, EP>::Symmetry hydr_tree;
   hydr_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
   // Compute density, pressure, acceleration due to pressure gradient
   dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
   setPressure(sph_system);
   hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
   // Get timestep
   dt = getTimeStepGlobal(sph_system);
   // Main loop for time integration
   PS::S32 step = 0;
   for(PS::F64 time = 0 ; time < end_time ; time += dt, ++ step){
      // Leap frog: Initial Kick & Full Drift
      InitialKick(sph_system, dt);
      FullDrift(sph_system, dt);
      // Adjust the positions of the SPH particles that run over
      // the computational boundaries.
      sph_system.adjustPositionIntoRootDomain(dinfo);
      // Leap frog: Predict
      Predict(sph_system, dt);
      // Perform domain decomposition again
      dinfo.decomposeDomainAll(sph_system);
      // Exchange the SPH particles between the (MPI) processes
      sph_system.exchangeParticle(dinfo);
      // Compute density, pressure, acceleration due to pressure gradient
      dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
      setPressure(sph_system);
      hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
      // Get a new timestep
      dt = getTimeStepGlobal(sph_system);
      // Leap frog: Final Kick
      FinalKick(sph_system, dt);
      // Output result files
      if(step % OUTPUT_INTERVAL == 0){
         FileHeader header;
         header.time = time;
         header.Nbody = sph_system.getNumberOfParticleGlobal();
         char filename[256];
         sprintf(filename, "result/%04d.txt", step);
         sph_system.writeParticleAscii(filename, header);
         if (PS::Comm::getRank() == 0){
            std::cout << "================================" << std::endl;
            std::cout << "output " << filename << "." << std::endl;
            std::cout << "================================" << std::endl;
         }
      }
      // Output information to STDOUT
      if (PS::Comm::getRank() == 0){
         std::cout << "================================" << std::endl;
         std::cout << "time = " << time << std::endl;
         std::cout << "step = " << step << std::endl;
         std::cout << "================================" << std::endl;
      }
      CheckConservativeVariables(sph_system);
   }
   // Finalize FDPS
   PS::Finalize();
   return 0;
}

