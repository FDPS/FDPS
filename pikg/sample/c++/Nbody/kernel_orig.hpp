#include <pikg_vector.hpp>

struct Kernel{
  const PIKG::F64 eps2;
  Kernel(const PIKG::F64 _eps2) : eps2(_eps2){}
  void operator()(const Particle * ep_i,
		  const PIKG::S32 n_ip,
		  const Particle * ep_j,
		  const PIKG::S32 n_jp,
		  Particle * force) {
    for(PIKG::S32 i = 0; i < n_ip; i++){
      PIKG::F64vec xi = ep_i[i].pos;
      PIKG::F64vec ai = 0.0;
      PIKG::F64 poti = 0.0;
      for(PIKG::S32 j = 0; j < n_jp; j++){
	PIKG::F64vec rij    = xi - ep_j[j].pos;
	PIKG::F64    r3_inv = rij * rij + eps2;
	PIKG::F64    r_inv  = 1.0/sqrt(r3_inv);
	r3_inv  = r_inv * r_inv;
	r_inv  *= ep_j[j].mass;
	r3_inv *= r_inv;
	ai     -= r3_inv * rij;
	poti   -= r_inv;
      }
      force[i].acc += ai;
      force[i].pot += poti;
    }
  }
};
