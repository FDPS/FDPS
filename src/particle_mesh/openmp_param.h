#ifndef __OPENMP_PARAM__
#define __OPENMP_PARAM__

#include<omp.h>

namespace ParticleSimulator{
    namespace ParticleMesh{


#ifdef KCOMPUTER
const int NUMBER_OF_OMP_THREADS = 8;
#else
const int NUMBER_OF_OMP_THREADS = 4;
#endif

const int NUMBER_OF_LEVEL_QSORT_SERIAL = 10;  // maximum
const int NSORTMAX = 1 << NUMBER_OF_LEVEL_QSORT_SERIAL;
#define CHUNK_QSORT schedule(dynamic, 2) 
#define CHUNK_PM schedule(dynamic) 
#define CHUNK_DECOMPOSITION schedule(dynamic) 
const int NLEVEL_KEY_SORT = 6;
const int NLEVEL_PARTICLE_SORT = 7;
const int NCHUNK_QSORT = 2;

#define CHUNK_FORCELOOP schedule(dynamic, 8)
const int NUMBER_OF_CHUNK_TREE_PARALLEL = 8;
const int NLEAF_PARALLEL_TREE_CONSTRUCTION = 10000;



    } // namespace ParticleMesh
}     // namespace ParticleSimulator

#endif
