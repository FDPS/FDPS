#pragma once
#include <cstring>
#include <particle_simulator.hpp>

class full_particle {
public:
   PS::S64 id;
   PS::F64 mass;
   PS::F64 eps;
   PS::F64vec pos;
   PS::F64vec vel;
   PS::F64 pot;
   PS::F64vec acc;

   PS::S64 getId() const { 
      return this->id;
   }
   PS::F64vec getPos() const { 
      return this->pos;
   }
   void setPos(const PS::F64vec pos_new) {
      this->pos = pos_new;
   }
   PS::F64 getCharge() const {
      return this->mass;
   }
   PS::F64 getChargeParticleMesh() const {
      return this->mass;
   }
   void copyFromForce(const full_particle & force) {
      this->pot = force.pot;
      this->acc = force.acc;
   }
   void copyFromFP(const full_particle & fp) {
      this->id = fp.id;
      this->mass = fp.mass;
      this->eps = fp.eps;
      this->pos = fp.pos;
   }
   void clear() {
      this->pot = 0.0;
      this->acc = 0.0;
   }
};

