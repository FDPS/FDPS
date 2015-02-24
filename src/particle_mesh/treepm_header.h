#ifndef _TREEPM_HEADER_INCLUDED
#define _TREEPM_HEADER_INCLUDED

#include<cstdio>
#include<cmath>
#include<iostream>
#include <fstream>
#include<sys/time.h>
#include<sys/resource.h>

#include "mpi.h"
#include "param.h"
#include "tools.h"

using namespace std;


namespace ParticleSimulator{
    namespace ParticleMesh{

typedef struct Particle{

  /* particle positive ID. 
   * If particle is boundary-superparticle, it is -1.
   */

  double xpos;
  double ypos;
  double zpos;

  float xvel;
  float yvel;
  float zvel;

#ifndef NOACC  // for debug
  float xacc;
  float yacc;
  float zacc;
#endif

  IDTYPE id;

#ifndef UNIFORM
  float mass;
#endif

#ifdef CALCPOT
  float pot;
#endif

}Particle, *pParticle;



typedef struct cosmology{

  float omega0;
  float lambda0;
  float omegab;
  float hubble;
  float tend;

}cosmology, *p_cosmology;



typedef struct TreePMTime{

  double pm;
  double pm_density_assignment;
  double pm_make_comm_density;
  double pm_density_comm;
  double pm_density_comm_reduce;
  double pm_fft_forward;
  double pm_density_phi;
  double pm_fft_backward;
  double pm_make_comm_phi;
  double pm_phi_comm_bcast;
  double pm_phi_comm;
  double pm_mesh_force;
  double pm_interpolation;

  double pp;
  double pp_get_boundary;
  double pp_get_boundary_determine_comm_node;
  double pp_get_boundary_make_key_sort;
  double pp_get_boundary_construct_tree;
  double pp_get_boundary_make_send_list;
  double pp_get_boundary_comm;
  double pp_get_boundary_push;
  double pp_make_key;
  double pp_key_sort;
  double pp_particle_sort;
  double pp_construct_tree;
  double pp_make_intr_list;
  double pp_calc_force;
  double pp_vel_kick;

  double drift;
  double next_dt;

  double exchange;
  double exchange_determine_boundary;
  double exchange_determine_boundary_comm;
  double exchange_make_send_list;
  double exchange_make_send_list_pre;
  double exchange_comm;
  double exchange_push;

  double pp_drift;

  double io_image;
  double io_density;
  double io_snapshot;

  double load_balance_log;

  double cputime;

}TreePMTime, *pTreePMTime;



typedef struct run_param{

  double lunit;
  double munit;
  double tunit;
  
  int   io_ver;
  int   npart;
  int   npart_old;
  int   ngas;
  float tnow;
  float anow;
  float znow;
  float hnow;
  float astart;   /* initial scale factor */
  struct cosmology cosm;

  long long int npart_total;  /* number of particle of system */\

  // mass information
  double uniform_mass;
  double min_mass;  /* for grape(g5_set_range) */


  // snapshot
  int noutput;    /* number of snapshot to output */
  int outflag;   
  double *toutput;  /* output epoch */
  double maxcpu; 
  int restart_flag;

  // log file
  FILE *log_file;
  FILE *profile_file;
  FILE *loadbalance_file;
  FILE *boundary_file;

  // node information
  int inode;      /* rank of own node */
  int nnode;      /* number of node */
  int ndiv[3];    /* division number for 3D */

  double bmin[3];   /* area of this node */
  double bmax[3];
  double global_bmin[3];  /* area of this group */
  double global_bmax[3];
  double other_bmin[3];   /* area of other group */
  double other_bmax[3];
  double bmin_all[MAXNNODE][3];
  double bmax_all[MAXNNODE][3];

  // comm
  int ncomm_tree;
  int nscomm_exchange;  
  int nrcomm_exchange;  
  int nsend_tree;
  int nsend_ex;
  int nrecv_ex;
  int nscomm_pm_l2s;
  int nrcomm_pm_l2s;
  int nsend_pm_l2s;
  int nrecv_pm_l2s;
  int nscomm_pm_s2l;
  int nrcomm_pm_s2l;
  int nsend_pm_s2l;
  int nrecv_pm_s2l;

  int nstep;

  //tree
  int bn;
  int nwalk;
  int nisum;
  int nimin;
  int nimax;
  int njsum;
  int njmin;
  int njmax;
  long long int ninteraction;
  double ninter_sum;

  // stat
  double nintrave;
  double nintrdisp;

  // timestep
  double dtprev;
  double dtime_ini;
  double dtime;
  double a2max;
  double v2max;
  double a2min;
  double v2min;
  double a2ave;
  double v2ave;
  double nrate;

  double dtv;
  double dta;

  TreePMTime t;

  char **argv;
  MPI_Comm MPI_COMM_INTERNAL;

  int inode_yz;
  MPI_Comm MPI_COMM_YZ;

  //PMForce *pm;

}RunParam, *pRunParam;





/**--------------------------------------------------------------------------------
 * libpg5.a
 * --------------------------------------------------------------------------------
 */
#ifndef GRAPE_OFF
extern "C"{
void pg5_set_plummer_table( const double eps, const double eta);
}
#endif



/**--------------------------------------------------------------------------------
 * pp2.c pp
 * --------------------------------------------------------------------------------
 */
void dumpParticle( const Particle *particle, FILE *outstream);


/**--------------------------------------------------------------------------------
 * io.cpp
 * --------------------------------------------------------------------------------
 */
int determineBoxOverlapping( const double cbox1[3],
                             const double box_half_size1[3],
                             const double cbox2[3],
                             const double box_half_size2[3]);
float swapFloat( const float f);
int swapInt( const int i);
double swapDouble( const double d);
long long int swapLLInt( const long long int ll);
void swapParticle( Particle *p);

void inputParam( pRunParam run_param, double *dtime, char const *paramfile);

void inputOneData( pRunParam this_run, pParticle particle, const char *filename, const int *ndiv);
void inputPluralDataContinue( pRunParam this_run, pParticle particle, char *basefile);
void inputPluralData( pRunParam this_run, pParticle particle, const char *filename);



void inputInitWrapper( pRunParam this_run, pParticle particle, 
		       char *init_file, const int nfiles,
		       const char *mass_file);

void inputInitWrapper3( pRunParam this_run, pParticle particle, 
			const int nfiles, const char *filename);

void outputSnapshotDump( const RunParam *this_run, 
			 const Particle *particle, const char *out_file);

void outputSnapshotDumpWrapper0( pParticle particle, pRunParam this_run, const char *filename);

void outputSnapshotWrapper( pParticle particle, pRunParam this_run, const double dtime);

void outputImageWrapper( const Particle *particle, pRunParam this_run);

void dumpAcc( pParticle particle, pRunParam this_run);

void readGadget( RunParam &this_run, pParticle p, const char *infile);
void outputGadget( RunParam &this_run, pParticle p, const char *outfile);




/**--------------------------------------------------------------------------------
 * evolve.c pp
 * --------------------------------------------------------------------------------
 */
double timetoz(double tnow, struct cosmology cosm);
void  funcd(double x, double *f, double *df, double tau);
double rtsafe(double x1, double x2, double xacc, double tau);

void step_pos( pParticle particle, pRunParam this_run, double dtime);
double getEps( const double znow, const double anow);
double getTheta( const double znow);
void getIntegralFactor( const pRunParam this_run, const double dt,
			double *vfac, double *afac);
void kick( pParticle particle, const float *a, 
		   const double vfac, const double afac);
void kick( pParticle particle, const float *a, 
		   const double afac);
void drift( pParticle particle, const double dt);
void drift( pParticle p, const int n, const double dt);
double zToTime( const double znow, cosmology &cosm);
void update_now( pRunParam run_param);
void correctBoundaryCondition( pParticle particle, const int n);
void TreePMTimePrint(const TreePMTime t, pRunParam this_run);
void logPrint( const RunParam *this_run, const double eps, const double global_time);
void initRunParam( pRunParam this_run, const int inode, const int nnode,
		   const int nowstep);
void genCommYZ( pRunParam this_run);
void printLoadBalance( RunParam *this_run, FILE *fout);
void printBoundary( const RunParam *this_run, FILE *fout,
		    const double *bmin, const double *bmax);



/**--------------------------------------------------------------------------------
 * ptobmp.cpp
 * --------------------------------------------------------------------------------
 */
int outputImage( const Particle *particle, const RunParam *this_run, const char *output_image);


/**--------------------------------------------------------------------------------
 * decomposition.cpp
 * --------------------------------------------------------------------------------
 */
void getXYZIndex( const int inode, const int *ndiv, int *x, int *y, int *z);
int getVoxelIndex( const int x, const int y, const int z, const int *ndiv);
void particleSortUsingXpos( pParticle particle, const int n);






inline double getMass( const Particle *particle, const double uniform_mass){

  if( particle->id == -1){
    return particle->xvel;
  }
  else{
#ifndef UNIFORM
    return particle->mass;
#else
    return uniform_mass;
#endif
  }

}



inline void getPos( const Particle *particle, double *pos){

  pos[0] = particle->xpos;
  pos[1] = particle->ypos;
  pos[2] = particle->zpos;

}



inline void getPos2( const Particle *particle, float *pos){

  pos[0] = particle->xpos;
  pos[1] = particle->ypos;
  pos[2] = particle->zpos;

}



/* x-y-z to rank */
inline int getVoxelIndex( const int x, const int y, const int z, const int *ndiv){

  return x*ndiv[1]*ndiv[2] + y*ndiv[2] + z;

}



inline int whichBox( const double *pos,
			    const double bmin[][3],
			    const double bmax[][3], 
			    const int npx,
			    const int npy, 
			    const int npz){

  int p = 0;
  if( pos[0] < bmin[p][0]) return -1;
  for(int ix=0; ix<npx; ix++, p+=npy*npz){
    if(pos[0] < bmax[p][0]) break;
  }
  if(pos[0] > bmax[p][0]) return -1;

  if(pos[1] < bmin[p][1]) return -1;
  for(int iy=0; iy<npy; iy++, p+=npz){
    if(pos[1] < bmax[p][1]) break;
  }
  if(pos[1] > bmax[p][1]) return -1;

  if(pos[2] < bmin[p][2]) return -1;
  for(int iz=0; iz<npz; iz++, p++){
    if(pos[2] < bmax[p][2]) break;
  }
  if(pos[2] > bmax[p][2]) return -1;

  return p;

}




inline double getHubble( const double a, const cosmology &cosm){

  double ainv = 1.0 / a;
  double ainv3 = ainv * ainv * ainv;
  return sqrt( cosm.omega0*ainv3 + cosm.lambda0);

}



class Kick{

 public:

  p_cosmology cosm;

  Kick( p_cosmology _cosm) : cosm(_cosm) {}

  inline double operator()(double anow){

    double hnow = getHubble( anow, *cosm);
    return 1.0 / ( anow * anow * hnow);
  }

};



class Drift{

 public:

  p_cosmology cosm;
  
  Drift( p_cosmology _cosm) : cosm(_cosm) {}

  inline double operator()(double anow){

    double hnow = getHubble( anow, *cosm);

    return 1.0 / ( anow * anow * anow * hnow);
  }

};



inline double getDriftFac( RunParam &this_run, const double anow, const double da){

  Drift d( &(this_run.cosm));
  return integralSimpson<Drift>( anow, anow+da, 1.0e-8, d);


}



inline double getKickFac( RunParam &this_run, const double anow, const double da){

  Kick k( &(this_run.cosm));
  return integralSimpson<Kick>( anow, anow+da, 1.0e-8, k);

}



    } // namespace ParticleMesh
}     // namespace ParticleSimulator

#endif
