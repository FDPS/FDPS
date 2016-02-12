#ifndef __PPPM__
#define __PPPM__

#include "openmp_param.h"
#include "gadget_param.h"

// A. Tanikawa add this.
#include "param_fdps.h"

namespace ParticleSimulator{
    namespace ParticleMesh{

typedef const char *const cchar;


/* Run ---------------------------------------------------------- */
//const int MAX_STEPS = 100000000;        /* max steps for calculation */
const int MAX_STEPS = 0;        /* max steps for calculation */
const double MAX_CPUTIME = 1e30;       /* max cputime (second)*/
//const double Z_FIN = 0.0;
const double Z_FIN = 1e30;
cchar MODEL="Test";
cchar SNAPSHOT="Snapshot";


/* Tree ---------------------------------------------------------- */
const int NJMAX = 130000;           /* j-particle buffer size (g5_get_jmemsize) */
const double Z_SWITCH_THETA = 5.0;  /* redshift to switch opening angle */
#if 1
const double THETA_HIGHZ = 0.5;     /* opening angle at high z*/
const double THETA = 0.5;           /* opening angle */
#else
const double THETA_HIGHZ = 0.0;     /* opening angle at high z*/
const double THETA = 0.0;           /* opening angle */
#endif
const int NCRIT = 300;              /* critical number of particles to share interaction list */
const int NLEAF = 10;               /* minimum number of particle not to divide tree further */
// A. Tanikawa comments out this.
const double CUTOFF_RADIUS = 3.0;


/* Timestep ------------------------------------------------ */
const double Z_SWITCH_ETA = 0.0;
const double ETA_HIGHZ = 0.3;      /* accuracy parameter */
const double ETA_TIMESTEP = 0.3;       /* accuracy parameter */
const double MAX_DLOGA = 0.03;
const double constant_timestep = 5e-3;  


/* Softening ------------------------------------------------ */
const double CONST_Z_VALUE = 0.0;  /* redshift to switch softening shape */


/* load balance  ----------------------------------------------------------------------- */
const int LOADBALANCE_METHOD = 1;   /* 0:interaction 1:pp+pm 2:n 3:1+lowerlimit 4:pp*/
const double SAMPLING_LOWER_LIMIT_FACTOR = 1.2;

const int NRATE_EXCHANGE = 2500;     /* sampling interval to decide area of each domain */
const int SORT_STEP = 10;      /* step interval to sort particle using morton ordering */
const int nstep_smoothing_boundary = 5;
const int nstep_decompose_x = 3;


/* IO --------------------------------------------------------------------------*/

const int NUMBER_OF_COMM_PER_ALLTOALLV_FOR_INPUT_IC = 16;
const int LOADBALANCELOG_ON = 0;
const int BOUNDARYLOG_ON = 0;
cchar STOPFILE = "stop";
const int fmerge = 2;

cchar DUMPDIR = "Dump";
cchar DUMPDIR0 = "Dump0";  // for initial dump;
cchar DUMPFILE = "dump";

const int ONELOGFLAG = 1;            /* 0: each node outputs log  1: root node outputs log */
const int IO_CACHE_SIZE = 2097152;   /* cache number of io (*3) */


/* Buffer ---------------------------------------------------------- */
const int CHARMAX = 256;
const int MAXNNODE = 16384;
const int MAXNNODE3 = MAXNNODE*3;
const int NXMAX = 100;
const int NSAMPMAX = 30000000;
const int bufsize_largemem0 = 10000000;
const int bufsize_largemem = 10000000;


/* PM ---------------------------------------------------------- */
const int _ndiv_fft[3] = {2, 4, 2};
const int _pm_reduce_nx = 1;
const int _pm_reduce_ny = 2;
const int _pm_reduce_nz = 1;

const int _flag_wisdom = 0;
cchar FFTW_WISDOM_FILENAME_F = "wisdom_f";
cchar FFTW_WISDOM_FILENAME_B = "wisdom_b";

const int PMLOG_ON = 0;
cchar PMLOG = "pm.log";

#ifdef FFT3D
#ifndef RMM_PM
#error
#endif
#ifndef FIX_FFTNODE
#error
#endif
#endif


/* Misc ---------------------------------------------------------- */
#define NOACC  /* debug */
// A. Tanikawa comments out this.
//#define UNIFORM  /* uniform mass mode */

/*cannot use multi mass mode in the case of n>2^31 */
#define LONG_ID
#define IDTYPE long long int
#define MPI_IDTYPE MPI_LONG_LONG_INT

cchar TORUSINFOFILE = "torus.dat";



#ifdef TREE2
#ifndef TREE_PARTICLE_CACHE
#error
#endif
#endif



/* on the fly analysis ---------------------------------------------------------- */
cchar IMAGEDIR = "Image";   /* relative path of image directory */
cchar NEWEST_IMAGEFILE = "newest.bmp";
const int IMAGESTEP = 32;   /* step interval to output image (0:no output image) */
const int IMAGEWIDTH = 768;   /* image width (px) */
const int IMAGEHEIGHT= 768;  /* image height (px) */
const double IMAGEFACA = 1.098612;    /* luminosity = s*(density-a) + b */
const double IMAGEFACB = 0.0;
const double IMAGEFACS = 22.77958;
const int COLORMAP = 0;  /* others:gray 1:blue 2:red */




/* NUMBER_OF_PART_ALL -> number of particle of system */
/* NUMBER_OF_PART     -> cache number to allocate particle array */
/* SIZE_OF_MESH       -> number of mesh of pm part (2*NUMBER_OF_PART**(1/3))*/
/* SFT_FOR_PP         -> softening parameter (ref evolve.c:get_eps) (1.0/(10*NUMBER_OF_PART**(1/3))) */

// A. Tanikawa defines this.
#define N32_2H

#ifdef N32_2H
#define NUMBER_OF_PART_ALL  (32768)
#define NUMBER_OF_PART      (32768)
// A. Tanikawa comments out this.
//#define SIZE_OF_MESH        (16)
#define SFT_FOR_PP          (2.5e-4)
#endif

#ifdef N128_2H
#define NUMBER_OF_PART_ALL  (2097152)
#define NUMBER_OF_PART      (2097152)
#define SIZE_OF_MESH        (64)
#define SFT_FOR_PP          (2.5e-4)
#endif

#ifdef N256_H
#define NUMBER_OF_PART_ALL  (16777216)
#define NUMBER_OF_PART      (16777216)
#define SIZE_OF_MESH        (256)
#define SFT_FOR_PP          (2.5e-4) 
#endif

#ifdef N256_2H
#define NUMBER_OF_PART_ALL  (16777216)
#define NUMBER_OF_PART      (16777216)
#define SIZE_OF_MESH        (128)
#define SFT_FOR_PP          (2.5e-4) 
#endif

#ifdef N512_2H
#define NUMBER_OF_PART_ALL  (134217728)
#define NUMBER_OF_PART      (20000000)
#define SIZE_OF_MESH        (256) 
#define SFT_FOR_PP          (6.25e-5) 
#endif

#define SIZE_OF_MESH_P2     (SIZE_OF_MESH+2)
#define SIZE_OF_GREEN       (SIZE_OF_MESH/2+1)
#define SFT_FOR_PM          (CUTOFF_RADIUS/SIZE_OF_MESH)
#define RADIUS_FOR_PP       (SFT_FOR_PM)
#define RADIUS_FOR_PP2      (RADIUS_FOR_PP*RADIUS_FOR_PP)

#define _MIN_(a,b) ( ((a)<(b)) ? (a) : (b) )
#define _MAX_(a,b) ( ((a)>(b)) ? (a) : (b) )

const static long long int nleafmax = 
  (long long int)NUMBER_OF_PART * (long long int)NLEAF / (long long int)NLEAF_PARALLEL_TREE_CONSTRUCTION;


    } // namespace ParticleMesh
}     // namespace ParticleSimulator

#endif /* __PPPM__ */
