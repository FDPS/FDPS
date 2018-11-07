#pragma once

/* 32 bit data types */
typedef int          fdps_s32;
typedef unsigned int fdps_u32;
#ifdef PARTICLE_SIMULATOR_ALL_64BIT_PRECISION
typedef double       fdps_f32;
#else
typedef float        fdps_f32;
#endif

/* 64 bit data types */
typedef long long int          fdps_s64;
typedef unsigned long long int fdps_u64;
typedef double                 fdps_f64;
