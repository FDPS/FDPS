#pragma once
#include "macro_defs.h"
#include "FDPS_c_if.h"
#include "user_defined.h"

void initial_kick(int psys_num_nbody, 
                  int psys_num_sph,
                  double dt);
void full_drift(int psys_num_nbody,
                int psys_num_sph,
                double dt);
void predict(int psys_num, double dt);
void final_kick(int psys_num_nbody,
                int psys_num_sph,
                double dt);
