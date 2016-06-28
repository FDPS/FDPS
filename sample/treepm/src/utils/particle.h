#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <stdint.h>

struct particle {
  int64_t id;
  float   mass;
  float   eps;
  double  xpos, ypos, zpos;
  double  xvel, yvel, zvel;
  double  xacc, yacc, zacc;
};

#endif /* __PARTICLE_H_ */
