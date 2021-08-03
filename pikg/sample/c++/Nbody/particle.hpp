#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include <pikg_vector.hpp>
struct Particle{
  PIKG::F64vec pos;
  PIKG::F64vec vel;
  PIKG::F64vec acc;
  PIKG::F64    mass;
  PIKG::F64    pot;
  static PIKG::F64 eps;

  void clear(){
    acc = 0.0;
    pot = 0.0;
  }
};

#endif
