#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<random>

#include <pikg_vector.hpp>
#include "particle.hpp"

#ifndef USE_CUDA_KERNEL
#include "kernel.hpp"
#else
PIKG::S32 DispatchKernel(const PIKG::S32          tag,
			 const PIKG::S32          n_walk,
			 const Particle    *epi[],
			 const PIKG::S32          n_epi[],
			 const Particle    *epj[],
			 const PIKG::S32          n_epj[]);
PIKG::S32 RetrieveKernel(const PIKG::S32 tag,
			 const PIKG::S32 n_walk,
			 const PIKG::S32 ni[],
			 Particle *force[]);
void initialize_Kernel(const double eps2_);

struct Kernel{
  Kernel(const PIKG::F64 eps2){
    initialize_Kernel(eps2);
  }
  void operator()(const Particle* epi,
		  const int nepi,
		  const Particle* epj,
		  const int nepj,
		  Particle* force){
    DispatchKernel(0,1,&epi,&nepi,&epj,&nepj);
    RetrieveKernel(0,1,&nepi,&force);
  }
};
#endif
//#include "kernel_orig.hpp"

void makeColdUniformSphere(const PIKG::F64 mass_glb,
                           const PIKG::S64 n_glb,
                           const PIKG::S64 n_loc,
                           PIKG::F64 *& mass,
                           PIKG::F64vec *& pos,
                           PIKG::F64vec *& vel,
                           const PIKG::F64 eng = -0.25,
                           const PIKG::S32 seed = 0) {
  assert(eng < 0.0);
  {
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::uniform_real_distribution<> dist(-1.0, 1.0);
    for(PIKG::S32 i = 0; i < n_loc; i++){
      mass[i] = mass_glb / n_glb;
      const PIKG::F64 radius = 3.0;
      do {
	pos[i][0] = dist(engine) * radius;
	pos[i][1] = dist(engine) * radius;
	pos[i][2] = dist(engine) * radius;
      }while(pos[i] * pos[i] >= radius * radius);
      vel[i][0] = 0.0;
      vel[i][1] = 0.0;
      vel[i][2] = 0.0;
    }
  }

  PIKG::F64vec cm_pos  = 0.0;
  PIKG::F64vec cm_vel  = 0.0;
  PIKG::F64    cm_mass = 0.0;
  for(PIKG::S32 i = 0; i < n_loc; i++){
    cm_pos  += mass[i] * pos[i];
    cm_vel  += mass[i] * vel[i];
    cm_mass += mass[i];
  }
  cm_pos /= cm_mass;
  cm_vel /= cm_mass;
  for(PIKG::S32 i = 0; i < n_loc; i++){
    pos[i] -= cm_pos;
    vel[i] -= cm_vel;
  }
}

void setParticlesColdUniformSphere(Particle* psys,
                                   const PIKG::S32 n_glb){
  PIKG::F64    * mass = new PIKG::F64[n_glb];
  PIKG::F64vec * pos  = new PIKG::F64vec[n_glb];
  PIKG::F64vec * vel  = new PIKG::F64vec[n_glb];
  const PIKG::F64 m_tot = 1.0;
  const PIKG::F64 eng   = -0.25;
  makeColdUniformSphere(m_tot, n_glb, n_glb, mass, pos, vel, eng);
  for(PIKG::S32 i = 0; i < n_glb; i++){
    psys[i].mass = mass[i];
    psys[i].pos  = pos[i];
    psys[i].vel  = vel[i];
  }
  delete [] mass;
  delete [] pos;
  delete [] vel;
}

void kick(Particle* system,
          const PIKG::F64 dt,
	  const PIKG::S32 n) {
  for(PIKG::S32 i = 0; i < n; i++) {
    system[i].vel  += system[i].acc * dt;
  }
}

void drift(Particle* system,
           const PIKG::F64 dt,
	   const PIKG::S32 n) {
  for(PIKG::S32 i = 0; i < n; i++) {
    system[i].pos  += system[i].vel * dt;
  }
}

void calcEnergy(const Particle* system,
                PIKG::F64 & etot,
                PIKG::F64 & ekin,
                PIKG::F64 & epot,
		const PIKG::S32 n){
  etot = ekin = epot = 0.0;
  for(PIKG::S32 i = 0; i < n; i++){
    ekin += system[i].mass * system[i].vel * system[i].vel;
    epot += system[i].mass * (system[i].pot + system[i].mass / Particle::eps);
  }
  ekin *= 0.5;
  epot *= 0.5;
  etot  = ekin + epot;
}

void printHelp() {
  std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
  std::cerr<<"t: theta (default: 0.5)"<<std::endl;
  std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
  std::cerr<<"s: time_step (default: 1.0 / 128.0)"<<std::endl;
  std::cerr<<"d: dt_diag (default: 1.0 / 8.0)"<<std::endl;
  std::cerr<<"D: dt_snap (default: 1.0)"<<std::endl;
  std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
  std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
  std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
  std::cerr<<"h: help"<<std::endl;
}

void makeOutputDirectory(char * dir_name) {
  struct stat st;
  PIKG::S32 ret;
  if (stat(dir_name, &st) != 0) {
    ret = mkdir(dir_name, 0777);
  } else {
    ret = 0; // the directory named dir_name already exists.
  }
  if (ret == 0) {
    fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
  } else {
    fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
    exit(-1);
  }
}

PIKG::F64 Particle::eps = 1.0/32.0;

int main(int argc, char *argv[]) {
  std::cout<<std::setprecision(15);
  std::cerr<<std::setprecision(15);

  PIKG::F32 theta = 0.5;
  PIKG::S32 n_leaf_limit = 8;
  PIKG::S32 n_group_limit = 64;
  PIKG::F32 time_end = 10.0;
  PIKG::F32 dt = 1.0 / 128.0;
  PIKG::F32 dt_diag = 1.0 / 8.0;
  PIKG::F32 dt_snap = 1.0;
  char dir_name[1024];
  PIKG::S64 nptcl = 1024;
  PIKG::S32 c;
  sprintf(dir_name,"./result");
  opterr = 0;
  while((c=getopt(argc,argv,"i:o:d:D:t:T:l:n:N:hs:")) != -1){
    switch(c){
    case 'o':
      sprintf(dir_name,optarg);
      break;
    case 't':
      theta = atof(optarg);
      std::cerr << "theta =" << theta << std::endl;
      break;
    case 'T':
      time_end = atof(optarg);
      std::cerr << "time_end = " << time_end << std::endl;
      break;
    case 's':
      dt = atof(optarg);
      std::cerr << "time_step = " << dt << std::endl;
      break;
    case 'd':
      dt_diag = atof(optarg);
      std::cerr << "dt_diag = " << dt_diag << std::endl;
      break;
    case 'D':
      dt_snap = atof(optarg);
      std::cerr << "dt_snap = " << dt_snap << std::endl;
      break;
    case 'l':
      n_leaf_limit = atoi(optarg);
      std::cerr << "n_leaf_limit = " << n_leaf_limit << std::endl;
      break;
    case 'n':
      n_group_limit = atoi(optarg);
      std::cerr << "n_group_limit = " << n_group_limit << std::endl;
      break;
    case 'N':
      nptcl = atoi(optarg);
      std::cerr << "nptcl = " << nptcl << std::endl;
      break;
    case 'h':
      printHelp();
      return 0;
    default:
      std::cerr<<"No such option! Available options are here."<<std::endl;
      printHelp();
    }
  }

  makeOutputDirectory(dir_name);

  std::ofstream fout_eng;

  char sout_de[1024];
  sprintf(sout_de, "%s/t-de.dat", dir_name);
  fout_eng.open(sout_de);
  fprintf(stdout, "This is a sample program of N-body simulation using PIKG!\n");

  Particle* system_grav = new Particle[nptcl];
  Kernel calcForce(Particle::eps*Particle::eps);
  setParticlesColdUniformSphere(system_grav, nptcl);
  calcForce(system_grav,nptcl,system_grav,nptcl,system_grav);

  const PIKG::F32 coef_ema = 0.3;
  PIKG::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
  calcEnergy(system_grav, Etot0, Ekin0, Epot0, nptcl);

  PIKG::F32 time_sys = 0.0;
  PIKG::F64 time_diag = 0.0;
  PIKG::S64 n_loop = 0;
  PIKG::S32 id_snap = 0;
  while(time_sys < time_end){
    calcEnergy(system_grav, Etot1, Ekin1, Epot1, nptcl);
    if( (time_sys >= time_diag) || ( (time_sys + dt) - time_diag ) > (time_diag - time_sys) ){
      fout_eng << time_sys << "   " << (Etot1 - Etot0) / Etot0 << std::endl;
      fprintf(stdout, "time: %10.7f energy error: %+e\n",
	      time_sys, (Etot1 - Etot0) / Etot0);
      time_diag += dt_diag;
    }

    kick(system_grav, dt * 0.5, nptcl);
    time_sys += dt;
    drift(system_grav, dt, nptcl);
    for(int i=0;i<nptcl;i++) system_grav[i].clear();
    calcForce(system_grav,nptcl,system_grav,nptcl,system_grav);
    kick(system_grav, dt * 0.5, nptcl);
    n_loop++;
  }

  delete[] system_grav;
  return 0;
}
