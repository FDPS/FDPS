#ifndef _MY_PP_INCLUDED
#define _MY_PP_INCLUDED

#include<cstdio>
#include<unistd.h>
#include<cstdlib>
#include<cassert>
#include<cstring>
#include<cmath>
#include<iostream>

#include "pm_parallel.h"
#include "treepm_header.h"

#ifndef GRAPE_OFF
#include "gp5util.h"
#endif


using namespace std;


namespace ParticleSimulator{
    namespace ParticleMesh{

typedef struct TreeParam{

  double tree_theta;
  double tree_theta2;          /* theta^2 */
  double tree_theta2_quarter;  /* theta^2 * 0.25 */
  int tree_ncrit;
  int tree_nleaf;
  double uniform_mass;

}TreeParam, *pTreeParam;



typedef struct TreeCell{

  int first_index;    
  int n;              // number of particle in this cell
  int next;           // linked-list

  float m;            // total mass
#ifdef KCOMPUTER
  float cm[3];        // gravity center
#else
  double cm[3];        // gravity center
#endif

#ifndef TREE2
  long long int key;  // morton key
  int zeroflag;       // 0x0000__** __(key_level) **(zerflag)
#endif

}TreeCell, *pTreeCell;



/* j-particle list */
typedef struct JList{

#ifdef KCOMPUTER
  double x[NJMAX][4] __attribute ((aligned (16)));
#else
  double m[NJMAX];
  double x[NJMAX][3];
#endif

}JList, *pJList;


typedef struct ICell{

  int id;
  long long int key;
  int key_level;

}ICell, *pICell;


typedef struct TreeTmp{

  int clistmax;         // max of tree cell
  int tree_clistmask;   // mask : key to adr
  pTreeCell tree_cell;  // pointer to tree_cell

  int walklist_size;    // number of i-list

#ifdef TREE2
  pICell walklist;
#else
  int *walklist;        // i-list cell index
#endif

  long long int *key;   // morton key
  int *index;           // local particle index

  double maxx;
  double bmin[3];

}TreeTmp, *pTreeTmp;



/**--------------------------------------------------------------------------------
 * pp2.cpp
 * --------------------------------------------------------------------------------
 */
void dumpAllParticle( const pParticle particle, const int n, FILE *fout);
void dumpTree( const TreeTmp *treetmp, FILE *outstream);

void TreeTmpTreeTmp( pTreeTmp tree_tmp, const pTreeParam treeparam, const int n);
void TreeTmpTreeTmp( pTreeTmp tree_tmp, const int n, const int offset);
void TreeTmpDTreeTmp( pTreeTmp tree_tmp);
void TreeTmpPrintMemsize( pTreeTmp tree_tmp, const int n);

void TreeCellTreeCell(pTreeCell tree_cell, const int n);
void initTreeParam( pTreeParam treeparam, const double theta,
		   const int ncrit, const int nleaf, const double uniform_mass);
int keyToAdr( const long long int key, const pTreeCell tree_cell, 
	      const int tree_clistmask);
void getCellRange( const pTreeCell tree_cell, double *ccell, double *length, 
		   const pParticle particle, const int *index);
void getCellRange2( const pTreeCell tree_cell, double *ccell, double *length);
void getCellRange2( const pTreeCell tree_cell, double *ccell, double *length,
		    const int key_level, const long long int key, const double maxx);
void getCellRange3( const pTreeCell tree_cell, double *ccell, double *length,
		    const pTreeTmp treetmp);
void makeInteractionList(const int ikey, const double cell_length2,
			 const double *cicell, const double *i_half_length,
			 const pTreeCell tree_cell,
			 pJList jlist, int *njlist,
			 const long long int current_key, const int tree_clistmask, 
			 const pParticle particle, const int *index,
			 const pTreeParam treeparam);
int treeConstructionParallel( pTreeTmp tree_tmp,
			      const int n,
			      pICell walklist,
			      const pTreeParam treeparam,
			      double (*p_cache)[4]);
void treeConstruction( pTreeCell tree_cell, const long long int *key, 
		       const int *index, const pParticle particle, 
		       const int first_index, const int n,
		       const int tree_clistmask, int *col_cell_index,
		       long long int current_key, int key_level,
		       int *walklist, int *nwalk, int ipflag,
		       const int parent_cell, const pTreeParam treeparam);
double makeKey( pParticle particle, long long int *key, const int n);

void calculateForce2( pParticle particle, const int *ilist, const int nilist,
		      const pJList jlist, const int njlist, float (*a)[3],
		      double *cicell);
void calculateForcePotSSE( pParticle particle, const int *ilist, const int nilist,
			   const pJList jlist, const int njlist, const double eps);
void calculateForceHost( pParticle particle, const int *ilist, const int nilist,
			 const pJList jlist, const int njlist, float (*a)[3], 
			 const double eps);


void calculateForceUsingTree(const pTreeCell tree_cell, const int ncell,
			     const pParticle particle, const int *index,
			     const int *walklist, const int nwalk,
			     const double maxx, const int tree_clistmask,
			     const double eps, const pTreeParam treeparam,
			     float (*a)[3], const int npart);
void calculateForceUsingMakingTree( pTreeTmp tree_tmp, const pParticle particle, 
				    const double eps,  const pTreeParam treeparam,
				    float (*a)[3], const int n);
void calc_PP_force( pParticle particle, pRunParam this_run,
		    pTreeTmp tree_tmp, TreeParam treeparam);



void calculateForceUsingTreeAndPMForce(const pTreeCell tree_cell, const int ncell,
				       const pParticle particle, const int *index,
				       const int *walklist, const int nwalk,
				       const double maxx, const int tree_clistmask,
				       const double eps, const pTreeParam treeparam,
				       const float *mesh_density, 
				       const int ncalc, const int npart,
				       const double vfac, 
				       const double afac,
				       pRunParam this_run);

void calculateForceUsingMakingTreeAndPMForce( pTreeTmp tree_tmp, 
					      const pParticle particle, 
					      pRunParam this_run,
					      const double eps,  
					      const pTreeParam treeparam,
					      float const *mesh_density,
					      const int ncalc, // n - boundary
					      const int n, //ncalc + boundary
					      const double vfac,
					      const double afac);

void getNextDt( pRunParam this_run, const double dtmid, const float eps);
void calcTreePMForce( pParticle particle, pRunParam this_run,
		      TreeParam treeparam, float const *mesh_density,
		      const float eps);

void calculateForceUsingTreeAndPMForce2(const pTreeCell tree_cell, const int ncell,
					const pParticle particle, const int *index,
					const pTreeTmp tree_tmp, const int nwalk,
					const double maxx, const int tree_clistmask,
					const double eps, const pTreeParam treeparam,
					PMForce &pm,
					const int ncalc, const int npart,
					const double vfac, 
					const double afac,
					pRunParam this_run);

void calculateForceUsingMakingTreeAndPMForce2( pTreeTmp tree_tmp, 
					       const pParticle particle, 
					       pRunParam this_run,
					       const double eps,  
					       const pTreeParam treeparam,
					       PMForce &pm,
					       const int ncalc, // n - boundary
					       const int n, //ncalc + boundary
					       const double vfac,
					       const double afac);
void calcTreePMForce2( pParticle particle, pRunParam this_run,
		       TreeParam treeparam, PMForce &pm,
		       const float eps);


void copyPCache( const pParticle particle, 
		 const int *index, const int n,
		 const double uniform_mass);


#ifndef KCOMPUTER
void makeInteractionList( const pTreeCell tree_cell, const int ipart_cell, 
			  const double *cicell, const double *i_half_length,
			  const int now_cell, const float *peri,
			  const double cell_length2,
			  pJList jlist, int *njlist,
			  const pTreeParam treeparam, 
			  const double (*p_cache)[4]);
#else
void makeInteractionList( const pTreeCell tree_cell, const int ipart_cell, 
			  const double *cicell, const double *i_half_length,
			  const int now_cell, const float *peri,
			  const double cell_length2,
			  pJList jlist, int *njlist,
			  const pTreeParam treeparam, 
			  const double (*p_cache)[4]);
#endif

void treeConstruction( pTreeTmp tree_tmp,
		       const int first_index, const int n,
		       long long int current_key, int key_level,
		       pICell walklist, int *nwalk, int ipflag,
		       const int parent_cell, const pTreeParam treeparam,
		       int *ncell_all, double (*p_cache)[4]);



void calcPMKick( pParticle particle, 
		 pRunParam this_run,
		 PMForce &pm, 
		 const double vfac,
		 const double afac);

void calcPPKick( pParticle particle, 
		 pRunParam this_run,
		 TreeParam treeparam, 
		 PMForce &pm, 
		 const double eps,
		 const double vfac,
		 const double afac);




/**--------------------------------------------------------------------------------
 * decomposition.cpp
 * --------------------------------------------------------------------------------
 */
void createDivision( const int n, int *ndiv);
void getRange( const Particle *particle, const int n,
	       double *bmin, double *bmax);
void setUniformBoundary( pRunParam this_run, 
			 double (*bmin_all_new)[3],
			 double (*bmax_all_new)[3]);
void exchangeParticle( pParticle particle, pRunParam this_run, 
		       double bmin_all[][3], double bmax_all[][3]);
void exchangeParticle( pParticle particle, pRunParam this_run, const int *ndiv);

void determineBoundary( float *xsamp_all, 
			float *ysamp_all,
			float *zsamp_all,
			const int nsamp,
			const pRunParam this_run,
			double *bmin,
			double *bmax);

void dumpSamplingParticle( const pParticle particle, 
			   const pRunParam this_run, 
			   const char *filename);

int getBoundaryParticle( pParticle particle, pRunParam this_run,
			 TreeParam *treeparam
#ifdef BUFFER_FOR_TREE
			 ,TreeTmp &tree_tmp
#endif
);
void createMPIParticle( MPI_Datatype *MPI_PARTICLE);


    } // namespace ParticleMesh
}     // namespace ParticleSimulator

#endif
