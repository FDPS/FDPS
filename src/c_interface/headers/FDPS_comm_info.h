#pragma once

//**** PS::CommInfo
#ifdef __cplusplus
extern "C" {
#endif    
typedef struct {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   MPI_Comm comm;
#endif    
   int rank;
   int n_proc;
}fdps_comm_info;
#ifdef __cplusplus
}
#endif
