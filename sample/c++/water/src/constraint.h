#ifndef H_CONSTRAINT
#define H_CONSTRAINT

#include "water_params.h"

class Constraint{
 private:
  int nbond;
  double *length = nullptr;

  const int max_iteration = 10000;
  const double acceptable_error = 1e-12;

  struct Pair{
    int a;
    int b;
  };
  Pair *pair = nullptr;

 public:
  Constraint(){
    nbond = 3;
    length = new double[nbond];
    pair   = new Pair[nbond];

    length[0] = BOND_OH;
    length[1] = BOND_HH;
    length[2] = BOND_OH;

    pair[0].a = 0; pair[0].b = 1;
    pair[1].a = 1; pair[1].b = 2;
    pair[2].a = 2; pair[2].b = 0;
  }

  ~Constraint(){
    if(length == nullptr) delete length;
    if(pair   == nullptr) delete pair;
  }

  template <class FP>
  void Shake(FP *psys,const PS::F64 dt){
    const PS::F64 dti = 1.0 / dt;
    PS::F64 max_error = 0.0;
    for(int iter=1;iter<=max_iteration;iter++){
      bool isEnough = true;
      for(int b=0;b<nbond;b++){
	const int i = pair[b].a;
	const int j = pair[b].b;
	PS::F64vec3 rij = psys[i].pos - psys[j].pos;
	const PS::F64 r2 = rij*rij;
	const PS::F64 b2 = length[b]*length[b];
	const PS::F64 error = fabs((sqrt(r2) - length[b])/length[b]);
	if(error > max_error) max_error = error;
	if(error > acceptable_error){
	  isEnough = false;
	  const PS::F64 mi_i = 1.0 / psys[i].mass;
	  const PS::F64 mi_j = 1.0 / psys[j].mass;
	  const PS::F64vec3 rij_old = psys[i].prev_pos - psys[j].prev_pos;
	  const PS::F64 gij = (r2 - b2) / (2.0*(rij*rij_old) * (mi_i + mi_j));
	  psys[i].pos -= gij * mi_i * rij_old;
	  psys[j].pos += gij * mi_j * rij_old;
	  psys[i].vel -= gij * mi_i * rij_old * dti;
	  psys[j].vel += gij * mi_j * rij_old * dti;
	}
      }
      if(isEnough) break;
      if(iter==max_iteration){
	std::cerr << "warning: SHAKE iteration is too long (>" << max_iteration << "), max error = " << max_error << std::endl;
      }
    }

  }

  template <class FP>
  void Rattle(FP *psys,
	      const PS::F64 dt){
    // iteration parameter
    const int max_iteration = 1000;
    const double max_error = 1e-12;

    for(int iter=1;iter<=max_iteration;iter++){
      bool isEnough = true;
      for(int b=0;b<nbond;b++){
	const int i = pair[b].a;
	const int j = pair[b].b;
	const PS::F64vec3 rij = psys[i].pos - psys[j].pos;
	const PS::F64vec3 vij = psys[i].vel - psys[j].vel;
	const PS::F64 rv = rij*vij;
	const PS::F64 delta = rv*dt/length[b];
	if(fabs(delta) > max_error){
	  isEnough = false;
	  const PS::F64 b2 = length[b]*length[b];
	  const PS::F64 mi_i = 1.0 / psys[i].mass;
	  const PS::F64 mi_j = 1.0 / psys[j].mass;
	  const PS::F64 kij = rv / (b2 * (mi_i + mi_j));
	  psys[i].vel -= kij*mi_i * rij;
	  psys[j].vel += kij*mi_j * rij;
	}
      }
      if(isEnough) break;
      if(iter==max_iteration){
	std::cerr << "warning: RATTLE iteration is too long (>" << max_iteration << ")" << std::endl;
      }
    }
  }
};

#endif
