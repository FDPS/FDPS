#pragma once

#include<cassert>
#include<fstream>

#include"../ps_defs.hpp"

namespace ParticleSimulator{

    namespace ParticleMesh{

        class ParticleMesh {
           private: 
              PMForce *pm_;
              RunParam this_run_;
              Particle *particle_;
              S32 n_loc_tot_;
              bool is_initialized {false};

              void Initialize() {
                 if (is_initialized == false) {
                    // Initialize private members
                     S32 rank, size;
                     MPI_Comm_rank(MPI_COMM_WORLD,&rank);
                     MPI_Comm_size(MPI_COMM_WORLD,&size);
                    this_run_.inode = rank;
                    this_run_.nnode = size;
                    this_run_.MPI_COMM_INTERNAL = MPI_COMM_WORLD;
                    pm_ = new PMForce(&this_run_);                
                    assert(pm_);
                    pm_->mesh_density_local = NULL;
                    particle_ = NULL;

                    // Update is_initialized
                    is_initialized = true;
                 }
              }
            
           public:
              ParticleMesh() { }

              template<class Tdinfo>
              void setDomainInfoParticleMesh(const Tdinfo & dinfo) {
                 Initialize();
                 S32 rank;
                 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
                 F32ort pos_domain = dinfo.getPosDomain(rank);
                 this_run_.ndiv[0] = dinfo.getNDomain(0);
                 this_run_.ndiv[1] = dinfo.getNDomain(1);
                 this_run_.ndiv[2] = dinfo.getNDomain(2);
                 this_run_.bmin[0] = pos_domain.low_[0];
                 this_run_.bmin[1] = pos_domain.low_[1];
                 this_run_.bmin[2] = pos_domain.low_[2];
                 this_run_.bmax[0] = pos_domain.high_[0];
                 this_run_.bmax[1] = pos_domain.high_[1];
                 this_run_.bmax[2] = pos_domain.high_[2];
              }

              template<class Tpsys>
              void setParticleParticleMesh(const Tpsys & psys,
                                           const bool clear=true) {
                 Initialize();

                 if(clear) {
                     if(particle_ != NULL)
                         delete [] particle_;
                     
                     n_loc_tot_ = psys.getNumberOfParticleLocal();
                     particle_ = new Particle[n_loc_tot_];
                     for(S32 i = 0; i < n_loc_tot_; i++) {
                         particle_[i].id   = 0;
                         particle_[i].xvel = 0.0;
                         particle_[i].mass =  psys[i].getChargeParticleMesh();
                         particle_[i].xpos = (psys[i].getPos())[0];
                         particle_[i].ypos = (psys[i].getPos())[1];
                         particle_[i].zpos = (psys[i].getPos())[2];
                     }
                 } else {
                     Particle *particle_prev = new Particle[n_loc_tot_];
                     for(S32 i = 0; i < n_loc_tot_; i++)
                         particle_prev[i] = particle_[i];
                     delete [] particle_;

                     S32 npsys = psys.getNumberOfParticleLocal();
                     particle_ = new Particle[n_loc_tot_+npsys];
                     for(S32 i = 0; i < n_loc_tot_; i++)
                         particle_[i] = particle_prev[i];
                     for(S32 i = 0; i < npsys; i++) {
                         S32 ii = i + n_loc_tot_;
                         particle_[ii].mass =  psys[i].getChargeParticleMesh();
                         particle_[ii].xpos = (psys[i].getPos())[0];
                         particle_[ii].ypos = (psys[i].getPos())[1];
                         particle_[ii].zpos = (psys[i].getPos())[2];                        
                     }
                     n_loc_tot_ += npsys;

                     delete [] particle_prev;
                 }
              }

              void calcMeshForceOnly() {
                 if (pm_->mesh_density_local != NULL) {
                     pm_->delLocal();
                 }
                 pm_->calcPMMeshForce(particle_, n_loc_tot_);
              }

              F32vec getForce(F32vec pos) {
                 F32 ia[3] = {0., 0., 0.};
                 Particle p;
                 p.xpos = pos[0];
                 p.ypos = pos[1];
                 p.zpos = pos[2];
                 pm_->forceInterpolation(&p, ia);
                 F32vec apm(ia[0], ia[1], ia[2]);
                 return apm;
              }

              F32 getPotential(F32vec pos) {
                 float pot = 0.0;
                 Particle p;
                 p.xpos = pos[0];
                 p.ypos = pos[1];
                 p.zpos = pos[2];
                 pm_->potentialInterpolation(&p, &pot);
                 return (F32) pot;
              }

              template<class Tpsys,
                       class Tdinfo>
              void calcForceAllAndWriteBack(Tpsys & psys,
                                            const Tdinfo & dinfo) {
                  setDomainInfoParticleMesh(dinfo);
                  setParticleParticleMesh(psys);
                  calcMeshForceOnly();
                  for(S32 i = 0; i < n_loc_tot_; i++)
                      psys[i].copyFromForceParticleMesh(getForce(psys[i].getPos()));
              }

        }; // END of class of ParticleMesh
        
    } // END of namespace of ParticleMesh
    
} // END of namespace of ParticleSimulator
