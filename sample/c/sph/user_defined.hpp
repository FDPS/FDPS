#pragma once
#include <cstring>
#include <particle_simulator.hpp>

class force_dens {
public:
   PS::F64 dens;
   PS::F64 smth;

   void clear() {
      this->dens = 0.0;
   }
};

class force_hydro {
public:
   PS::F64vec acc;
   PS::F64 eng_dot;
   PS::F64 dt;

   void clear() {
      this->acc = 0.0;
      this->eng_dot = 0.0;
      this->dt = 0.0;
   }
};

class full_particle {
public:
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
   PS::F64 getRSearch() const {
      return this->smth;
   }
   void copyFromForce(const force_dens & force) {
      this->dens = force.dens;
   }
   void copyFromForce(const force_hydro & force) {
      this->acc = force.acc;
      this->eng_dot = force.eng_dot;
      this->dt = force.dt;
   }
};

class essential_particle {
public:
   PS::S64 id;
   PS::F64vec pos;
   PS::F64vec vel;
   PS::F64 mass;
   PS::F64 smth;
   PS::F64 dens;
   PS::F64 pres;
   PS::F64 snds;

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
   PS::F64 getRSearch() const {
      return this->smth;
   }
   void copyFromFP(const full_particle & fp) {
      this->id = fp.id;
      this->pos = fp.pos;
      this->vel = fp.vel;
      this->mass = fp.mass;
      this->smth = fp.smth;
      this->dens = fp.dens;
      this->pres = fp.pres;
      this->snds = fp.snds;
   }
};

