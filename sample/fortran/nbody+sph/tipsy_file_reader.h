#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <cstdint>
#include <array>

/* Structure definitions */
// The following two classes are used to read particle data created by MAGI.
typedef struct magi_tipsy_header {
   double time;
   int nbodies;
   int ndim;
   int nsph;
   int ndark;
   int nstar;
} Magi_tipsy_header;

typedef struct magi_tipsy_particle {
   float mass;
   float pos[3];
   float vel[3];
   float eps;
   int idx;
} Magi_tipsy_particle;
