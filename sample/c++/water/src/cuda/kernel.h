#ifndef H_FORCE_KERNEL
#define H_FORCE_KERNEL

#include "particle_simulator.hpp"
#include "user_defined_class.h"

#include <cmath>

PS::S32 DispatchKernel(
		       const PS::S32          tag,
		       const PS::S32          n_walk,
		       const EP              *epi[],
		       const PS::S32          n_epi[],
		       const EP              *epj[],
		       const PS::S32          n_epj[]);

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       Force *force[]);

#endif
