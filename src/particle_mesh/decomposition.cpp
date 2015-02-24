#include "pp.h"

namespace ParticleSimulator{
    namespace ParticleMesh{

extern double (*p_cache)[4];

static char cbuf[CHARMAX];


void createDivision( const int n, int *ndiv){


  int nx, ny, nz;
  int n0, n1;
  n0 = (int)pow(n+0.1,0.33333333333333333333);
  while(n%n0)n0--;
  nx = n0;
  n1 = n/nx;
  n0 = (int)sqrt(n1+0.1);
  while(n1%n0)n0++;
  ny = n0; nz = n1/n0;
  int ntmp;
  if (nz > ny){
    ntmp = nz; nz = ny; ny = ntmp;
  }
  if (ny > nx){
    ntmp = nx; nx = ny; ny = ntmp;
  }
  if (nz > ny){
    ntmp = nz; nz = ny; ny = ntmp;
  }

  ndiv[0] = nx;
  ndiv[1] = ny;
  ndiv[2] = nz;

}

void getRange( const Particle *particle, const int n,
	       double *bmin, double *bmax){

  bmin[0] = 1.0e30;
  bmin[1] = 1.0e30;
  bmin[2] = 1.0e30;
  bmax[0] = -1.0e30;
  bmax[1] = -1.0e30;
  bmax[2] = -1.0e30;

  int i;
  for( i=0; i<n; i++){
    if( particle[i].xpos < bmin[0])  bmin[0] = particle[i].xpos;
    if( particle[i].xpos > bmax[0])  bmax[0] = particle[i].xpos;
    if( particle[i].ypos < bmin[1])  bmin[1] = particle[i].ypos;
    if( particle[i].ypos > bmax[1])  bmax[1] = particle[i].ypos;
    if( particle[i].zpos < bmin[2])  bmin[2] = particle[i].zpos;
    if( particle[i].zpos > bmax[2])  bmax[2] = particle[i].zpos;
  }

}



/* rank to x-y-z*/
void getXYZIndex( const int inode, const int *ndiv, int *x, int *y, int *z){

  if( ndiv[1] == 1){
    *x = inode;
    *y = 0;
    *z = 0;
  }
  else{
    if( ndiv[2] == 1){
      *z = 0;
      *x = inode / ndiv[1];
      *y = inode % ndiv[1];
    }
    else{
      int xy = ndiv[1] * ndiv[2];
      *x = inode / xy;
      int inode2 = inode % xy;
      *y = inode2 / ndiv[2];
      *z = inode2 % ndiv[2];
    }
  }

}




void setUniformBoundary( pRunParam this_run, 
			 double (*bmin_all_new)[3],
			 double (*bmax_all_new)[3]){

  double hx = (this_run->global_bmax[0] - this_run->global_bmin[0]) / (double)this_run->ndiv[0];
  double hy = (this_run->global_bmax[1] - this_run->global_bmin[1]) / (double)this_run->ndiv[1];
  double hz = (this_run->global_bmax[2] - this_run->global_bmin[2]) / (double)this_run->ndiv[2];

  double bmin[3], bmax[3];
  int x, y, z;
  getXYZIndex( this_run->inode, this_run->ndiv, &x, &y, &z);
  bmin[0] = (double)x*hx + this_run->global_bmin[0];
  bmin[1] = (double)y*hy + this_run->global_bmin[1];
  bmin[2] = (double)z*hz + this_run->global_bmin[2];

#ifdef SHIFT_INITIAL_BOUNDARY
  double ds = 0.001 / pow( (double)NUMBER_OF_PART_ALL, 0.33333);
  fprintf_verbose( stderr, "shift_initial_boundary %e\n", ds);
  cerr << this_run->npart << "\t" << ds << endl;
  bmin[0] += ds;
  bmin[1] += ds;
  bmin[2] += ds;
#endif

  bmax[0] = bmin[0] + hx;
  bmax[1] = bmin[1] + hy;
  bmax[2] = bmin[2] + hz;
  if( x == 0)  bmin[0] = this_run->global_bmin[0];
  if( y == 0)  bmin[1] = this_run->global_bmin[1];
  if( z == 0)  bmin[2] = this_run->global_bmin[2];
  if( x == this_run->ndiv[0]-1)  bmax[0] = this_run->global_bmax[0];
  if( y == this_run->ndiv[1]-1)  bmax[1] = this_run->global_bmax[1];
  if( z == this_run->ndiv[2]-1)  bmax[2] = this_run->global_bmax[2];
  MPI_Allgather( bmin, 3, MPI_DOUBLE, 
		 bmin_all_new, 3, MPI_DOUBLE, this_run->MPI_COMM_INTERNAL);
  MPI_Allgather( bmax, 3, MPI_DOUBLE, 
		 bmax_all_new, 3, MPI_DOUBLE, this_run->MPI_COMM_INTERNAL);

  this_run->bmin[0] = bmin[0];
  this_run->bmin[1] = bmin[1];
  this_run->bmin[2] = bmin[2];
  this_run->bmax[0] = bmax[0];
  this_run->bmax[1] = bmax[1];
  this_run->bmax[2] = bmax[2];

}



/* sort particle for exchange */
int sortInBoxParticle( Particle *particle, const int n, const double *bmin,
		       const double *bmax){

  /*
   * particle ___________________________________ 
   *                    sort using area
   *          -----------------------------------
   *          \no change node      /\change node/
   */

#ifdef STATIC_ARRAY
  static int inbox_flag[NUMBER_OF_PART];
#else
  int *inbox_flag = (int *) my_malloc( sizeof(int) * n);
#endif

  int oldn = 0;
#pragma omp parallel for reduction(+:oldn)
  for( int i=0; i<n; i++){
    double pos[3];
    getPos( &particle[i], pos);
    inbox_flag[i] = 1;
    if( (bmin[0] <= pos[0]  &&  pos[0] < bmax[0]) &&
	(bmin[1] <= pos[1]  &&  pos[1] < bmax[1]) &&
	(bmin[2] <= pos[2]  &&  pos[2] < bmax[2])){
      inbox_flag[i] = 0;
      oldn ++;
    }
  }

  int i = 0;
  int j = n - 1;
  while(i < j){
    if( inbox_flag[i] == 0){
      i++;
    }
    else if( inbox_flag[j] == 1){
      j--;
    }
    else{
      Particle particle_tmp = particle[i];
      particle[i] = particle[j];
      particle[j] = particle_tmp;
      i++;
      j--;
    }
  }

#ifndef STATIC_ARRAY
  free( inbox_flag);
#endif

  return oldn;

}


void particleSortUsingXpos( pParticle particle, const int n){

  int i = 0;
  int j = n-1;
  double pivot = particle[n/2].xpos;

  while(1){
    while( particle[i].xpos < pivot){
      i++;
    }
    while( particle[j].xpos > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    Particle tmp_p = particle[j];
    particle[j] = particle[i];
    particle[i] = tmp_p;
    i++;
    j--;

  }

  if( 0 < i-1){
    particleSortUsingXpos( particle, i);
  }
  if( j+1 < n-1){
    particleSortUsingXpos( particle+j+1, n-j-1);
  }

}



void particleSortUsingYpos( pParticle particle, const int n){

  int i = 0;
  int j = n-1;
  double pivot = particle[n/2].ypos;

  while(1){
    while( particle[i].ypos < pivot){
      i++;
    }
    while( particle[j].ypos > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    Particle tmp_p = particle[j];
    particle[j] = particle[i];
    particle[i] = tmp_p;
    i++;
    j--;

  }

  if( 0 < i-1){
    particleSortUsingYpos( particle, i);
  }
  if( j+1 < n-1){
    particleSortUsingYpos( particle+j+1, n-j-1);
  }

}



void particleSortUsingZpos( pParticle particle, const int n){

  int i = 0;
  int j = n-1;
  double pivot = particle[n/2].zpos;

  while(1){
    while( particle[i].zpos < pivot){
      i++;
    }
    while( particle[j].zpos > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    Particle tmp_p = particle[j];
    particle[j] = particle[i];
    particle[i] = tmp_p;
    i++;
    j--;

  }

  if( 0 < i-1){
    particleSortUsingZpos( particle, i);
  }
  if( j+1 < n-1){
    particleSortUsingZpos( particle+j+1, n-j-1);
  }

}



void particleSortForExchange( pParticle particle, pRunParam this_run, const int dn,
			      double bmax[][3], int *nsend){

  int n = this_run->npart;
  int oldn = n - dn;
  int nx = this_run->ndiv[0];
  int ny = this_run->ndiv[1];
  int nz = this_run->ndiv[2];

  for( int i=0; i<this_run->nnode; i++){
    nsend[i] = 0;
  }

  if( dn == 0)  return;
  qsortParticleUsingXpos( particle+oldn, dn, NLEVEL_KEY_SORT, NCHUNK_QSORT);

#if 0
  int i, j, k;
  int nxstart = oldn;
  for( i=0; i<nx; i++){
    int bxnode = getVoxelIndex( i, 0, 0, this_run->ndiv);
    double xmax = bmax[bxnode][0];
    int j;
    for( j=nxstart; j<n; j++){
      if( particle[j].xpos >= xmax)  break;
    }
    int nxend = j;
    if( i == nx-1)  nxend = n;
    int dny = nxend - nxstart;
    if( dny == 0)  continue;
    particleSortUsingYpos( &particle[nxstart], dny);

    //fprintf( stderr, "%d %d\n", i, dny);

    int nystart = nxstart;
    for( j=0; j<ny; j++){
      int bynode = getVoxelIndex( i, j, 0, this_run->ndiv);
      double ymax = bmax[bynode][1];
      for( k=nystart; k<nxend; k++){
	if( particle[k].ypos >= ymax)  break;
      }
      int nyend = k;
      if( j == ny-1)  nyend = nxend;
      int dnz = nyend - nystart;
      if( dnz == 0)  continue;
      particleSortUsingZpos( &particle[nystart], dnz);

      int bznode = bynode;
      int lstart = nystart;
      for( k=bznode; k<(bznode+nz); k++){
	//fprintf( stderr, "lstart %d %d %d %d\n", this_run->inode, lstart, nyend, oldn);
	double zmax = bmax[k][2];
	int l;
	for( l=lstart; l<nyend; l++){
	  if( particle[l].zpos < zmax){
	    nsend[k] ++;
	  }
	  else{
	    lstart = l;
	    break;
	  }
	}
	if( l == nyend)  lstart = nyend;
      }

      nystart = nyend;
    }

    nxstart = nxend;

  } 
#else

#ifdef STATIC_ARRAY
  static int nxstart[NXMAX];
  static int nxend[NXMAX];
#else
  int *nxstart = new int[nx];
  int *nxend = new int[nx];
#endif
  nxstart[0] = oldn;
  for( int i=0; i<(nx-1); i++){
    int bxnode = getVoxelIndex( i, 0, 0, this_run->ndiv);
    double xmax = bmax[bxnode][0];
    int j;
    for( j=nxstart[i]; j<n; j++){
      if( particle[j].xpos >= xmax){
	break;
      }
    }
    nxend[i] = j;
    nxstart[i+1] = nxend[i];
  }
  nxend[nx-1] = n;

#pragma omp parallel for CHUNK_DECOMPOSITION
  for( int i=0; i<nx; i++){
    int dny = nxend[i] - nxstart[i];
    if( dny == 0)  continue;
    particleSortUsingYpos( &particle[nxstart[i]], dny);

    int j, k;
    int nystart = nxstart[i];
    for( j=0; j<ny; j++){
      int bynode = getVoxelIndex( i, j, 0, this_run->ndiv);
      double ymax = bmax[bynode][1];
      for( k=nystart; k<nxend[i]; k++){
	if( particle[k].ypos >= ymax)  break;
      }
      int nyend = k;
      if( j == ny-1)  nyend = nxend[i];
      int dnz = nyend - nystart;
      if( dnz == 0)  continue;
      particleSortUsingZpos( &particle[nystart], dnz);

      int bznode = bynode;
      int lstart = nystart;
      for( k=bznode; k<(bznode+nz); k++){
	//fprintf( stderr, "lstart %d %d %d %d\n", this_run->inode, lstart, nyend, oldn);
	double zmax = bmax[k][2];
	int l;
	for( l=lstart; l<nyend; l++){
	  if( particle[l].zpos < zmax){
	    nsend[k] ++;
	  }
	  else{
	    lstart = l;
	    break;
	  }
	}
	if( l == nyend)  lstart = nyend;
      }
      
      nystart = nyend;
    }
  }

#ifndef STATIC_ARRAY
  delete [] nxstart;
  delete [] nxend;
#endif

#endif

}



int samplingParticleX( const pParticle particle, const int n, const int nrate, 
		       float *xsamp_all,
		       const pRunParam this_run){

  int nnode = this_run->nnode;
  int nsamp_local = n / nrate + 1;

  float *xsamp = new float[nsamp_local];
 
  nsamp_local = 0;
#pragma omp parallel for reduction(+:nsamp_local)
  for( int i=0; i<n; i++){
    if( i%nrate != 0)  continue;
    int ip = i/nrate;
    xsamp[ip] = particle[i].xpos;
    nsamp_local ++;
  }

  static int nsamp_all_buf[MAXNNODE];
  static int displs[MAXNNODE];

#ifdef KCOMPUTER
  MPI_Allgather( &nsamp_local, 1, MPI_INT, nsamp_all_buf, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);
#else
  MPI_Gather( &nsamp_local, 1, MPI_INT, nsamp_all_buf, 1, MPI_INT, 0, this_run->MPI_COMM_INTERNAL);
#endif
  displs[0] = 0;
  for( int i=1; i<nnode; i++){
    displs[i] = displs[i-1] + nsamp_all_buf[i-1];
  }
  int nsamp = displs[nnode-1] + nsamp_all_buf[nnode-1];

#ifdef KCOMPUTER
  MPI_Allgatherv( xsamp, nsamp_local, MPI_FLOAT,
		  xsamp_all, nsamp_all_buf, displs, 
		  MPI_FLOAT, this_run->MPI_COMM_INTERNAL);
#else
  MPI_Gatherv( xsamp, nsamp_local, MPI_FLOAT,
	       xsamp_all, nsamp_all_buf, displs, 
	       MPI_FLOAT, 0, this_run->MPI_COMM_INTERNAL);
#endif

  delete [] xsamp;

  return nsamp;

}



int samplingParticleYZ( const pParticle particle, const int n, const int nrate, 
			float *ysamp_all, float *zsamp_all,
			const pRunParam this_run){

  int nnode = this_run->ndiv[1] * this_run->ndiv[2];
  int nsamp_local = n / nrate + 1;

  float *ysamp = new float[nsamp_local];
  float *zsamp = new float[nsamp_local];
 
  nsamp_local = 0;
#pragma omp parallel for reduction(+:nsamp_local)
  for( int i=0; i<n; i++){
    if( i%nrate != 0)  continue;
    int ip = i/nrate;
    ysamp[ip] = particle[i].ypos;
    zsamp[ip] = particle[i].zpos;
    nsamp_local ++;
  }

  static int nsamp_all_buf[MAXNNODE];
  static int displs[MAXNNODE];

  MPI_Allgather( &nsamp_local, 1, MPI_INT, nsamp_all_buf, 1, MPI_INT, 
		 this_run->MPI_COMM_YZ);
  displs[0] = 0;
  for( int i=1; i<nnode; i++){
    displs[i] = displs[i-1] + nsamp_all_buf[i-1];
  }
  int nsamp = displs[nnode-1] + nsamp_all_buf[nnode-1];

  MPI_Allgatherv( ysamp, nsamp_local, MPI_FLOAT,
		  ysamp_all, nsamp_all_buf, displs, 
		  MPI_FLOAT, this_run->MPI_COMM_YZ);
  MPI_Allgatherv( zsamp, nsamp_local, MPI_FLOAT,
		  zsamp_all, nsamp_all_buf, displs, 
		  MPI_FLOAT, this_run->MPI_COMM_YZ);

  delete [] ysamp;
  delete [] zsamp;

  return nsamp;

}




int samplingParticle( const pParticle particle, const int n, const int nrate, 
		      float *xsamp_all, float *ysamp_all, float *zsamp_all,
		      const pRunParam this_run){

  int nnode = this_run->nnode;

  int nsamp_local = n / nrate + 1;

  float *xsamp = new float[nsamp_local];
  float *ysamp = new float[nsamp_local];
  float *zsamp = new float[nsamp_local];
 
  nsamp_local = 0;
#pragma omp parallel for reduction(+:nsamp_local)
  for( int i=0; i<n; i++){
    if( i%nrate != 0)  continue;
    int ip = i/nrate;
    float  pos[3];
    getPos2( &particle[i], pos);
    xsamp[ip] = pos[0];
    ysamp[ip] = pos[1];
    zsamp[ip] = pos[2];
    nsamp_local ++;
  }

  static int nsamp_all_buf[MAXNNODE];
  static int displs[MAXNNODE];

  MPI_Allgather( &nsamp_local, 1, MPI_INT, nsamp_all_buf, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);
  displs[0] = 0;
  for( int i=1; i<nnode; i++){
    displs[i] = displs[i-1] + nsamp_all_buf[i-1];
  }
  int nsamp = displs[nnode-1] + nsamp_all_buf[nnode-1];

  MPI_Allgatherv( xsamp, nsamp_local, MPI_FLOAT,
		  xsamp_all, nsamp_all_buf, displs, 
		  MPI_FLOAT, this_run->MPI_COMM_INTERNAL);
  MPI_Allgatherv( ysamp, nsamp_local, MPI_FLOAT,
		  ysamp_all, nsamp_all_buf, displs, 
		  MPI_FLOAT, this_run->MPI_COMM_INTERNAL);
  MPI_Allgatherv( zsamp, nsamp_local, MPI_FLOAT,
		  zsamp_all, nsamp_all_buf, displs, 
		  MPI_FLOAT, this_run->MPI_COMM_INTERNAL);

  delete [] xsamp;
  delete [] ysamp;
  delete [] zsamp;

  return nsamp;

}



void determineBoundary( float *xsamp_all, 
			float *ysamp_all,
			float *zsamp_all,
			const int nsamp,
			const pRunParam this_run,
			double *bmin,
			double *bmax){

  int inode = this_run->inode;
  int *ndiv = this_run->ndiv;
  int nnode2 = ndiv[1] * ndiv[2];

  qsortPivotWrapper( ysamp_all, zsamp_all, xsamp_all, nsamp, NLEVEL_PARTICLE_SORT, NCHUNK_QSORT);

  int xindex = inode / nnode2;
  if( ndiv[1] == 1)   xindex = inode;
  int xleft = xindex * nsamp/ndiv[0];
  float xmin = xsamp_all[xleft];
  if( xindex == 0)  xmin = this_run->global_bmin[0];
  int xright = (xindex+1) * nsamp/ndiv[0];
  float xmax = xsamp_all[xright];
  if( xindex == (ndiv[0]-1))  xmax = this_run->global_bmax[0];
  
  int nsamp2 = xright - xleft;
  qsortPivotWrapper( zsamp_all+xleft, ysamp_all+xleft, nsamp2, NLEVEL_PARTICLE_SORT, NCHUNK_QSORT);
  int inode2 = inode - xindex*nnode2;
  int yindex = inode2 / ndiv[2];
  if( ndiv[2] == 1)  yindex = inode2;
  int yleft = yindex * nsamp2/ndiv[1] + xleft;
  float ymin = ysamp_all[yleft];
  if( yindex == 0)  ymin = this_run->global_bmin[1];
  int yright = (yindex+1) * nsamp2/ndiv[1] + xleft;
  float ymax = ysamp_all[yright];
  if( yindex == (ndiv[1]-1))  ymax = this_run->global_bmax[1];

  int nsamp3 = yright - yleft;
  qsortBasicWrapper( zsamp_all+yleft, nsamp3, NLEVEL_PARTICLE_SORT, NCHUNK_QSORT);
  int inode3 = inode2 - yindex*ndiv[2];
  int zindex = inode3;
  int zleft = zindex * nsamp3/ndiv[2] + yleft;
  float zmin = zsamp_all[zleft];
  if( zindex == 0)  zmin = this_run->global_bmin[2];
  int zright = (zindex+1) * nsamp3/ndiv[2] + yleft;
  float zmax = zsamp_all[zright];
  if( zindex == (ndiv[2]-1))  zmax = this_run->global_bmax[2];

  bmin[0] = xmin;
  bmin[1] = ymin;
  bmin[2] = zmin;
  bmax[0] = xmax;
  bmax[1] = ymax;
  bmax[2] = zmax;

}



void determineBoundary2( float *xsamp_all, 
			 float *ysamp_all,
			 float *zsamp_all,
			 const int nsamp_x,
			 const int nsamp_yz,
			 const pRunParam this_run,
			 double *bmin,
			 double *bmax,
			 const int sampx_on){

  int inode = this_run->inode;
  int inode_yz = this_run->inode_yz;
  int nnode = this_run->nnode;
  int *ndiv = this_run->ndiv;
  int nnode2 = ndiv[1] * ndiv[2];

  static float bminall[MAXNNODE];
  static float bmaxall[MAXNNODE];
  int nsamp_per_x = nsamp_x / ndiv[0];
  float xmin = this_run->bmin[0];
  float xmax = this_run->bmax[0];
  if( sampx_on == 1){
    if( inode == 0){
      qsortBasicWrapper( xsamp_all, nsamp_x, NLEVEL_PARTICLE_SORT, NCHUNK_QSORT);
#pragma omp parallel for
      for( int i=0; i<nnode; i++){
	int xindex = i / nnode2;
	int xleft = xindex * nsamp_per_x;
	bminall[i] = xsamp_all[xleft];
	if( xindex == 0)  bminall[i] = this_run->global_bmin[0];
	int xright = (xindex+1) * nsamp_per_x;
	bmaxall[i] = xsamp_all[xright];
	if( xindex == (ndiv[0]-1))  bmaxall[i] = this_run->global_bmax[0];
      }
    }
    MPI_Scatter( bminall, 1, MPI_FLOAT, &xmin, 1, MPI_FLOAT, 0, this_run->MPI_COMM_INTERNAL);
    MPI_Scatter( bmaxall, 1, MPI_FLOAT, &xmax, 1, MPI_FLOAT, 0, this_run->MPI_COMM_INTERNAL);
  }

  qsortPivotWrapper( zsamp_all, ysamp_all, nsamp_yz, NLEVEL_PARTICLE_SORT, NCHUNK_QSORT);
  int nsamp_per_y = nsamp_yz / ndiv[1];
  int yindex = inode_yz / ndiv[2];
  int yleft = yindex * nsamp_per_y;
  float ymin = ysamp_all[yleft];
  if( yindex == 0)  ymin = this_run->global_bmin[1];
  int yright = (yindex+1) * nsamp_per_y;
  float ymax = ysamp_all[yright];
  if( yindex == (ndiv[1]-1))  ymax = this_run->global_bmax[1];

  int nsamp_z = yright - yleft;
  qsortBasicWrapper( zsamp_all+yleft, nsamp_z, NLEVEL_PARTICLE_SORT, NCHUNK_QSORT);
  int nsamp_per_z = nsamp_z / ndiv[2];
  int inode_z = inode_yz - yindex*ndiv[2];
  int zleft = inode_z * nsamp_per_z + yleft;
  float zmin = zsamp_all[zleft];
  if( inode_z == 0)  zmin = this_run->global_bmin[2];
  int zright = (inode_z+1) * nsamp_per_z + yleft;
  float zmax = zsamp_all[zright];
  if( inode_z == (ndiv[2]-1))  zmax = this_run->global_bmax[2];

  bmin[0] = xmin;
  bmin[1] = ymin;
  bmin[2] = zmin;
  bmax[0] = xmax;
  bmax[1] = ymax;
  bmax[2] = zmax;

}



void dumpSamplingParticle( const pParticle particle, 
			   const pRunParam this_run, 
			   const char *filename){

  int nrate = NRATE_EXCHANGE;
  int nsamp = this_run->npart_total / nrate + 1;
  int nmem = nsamp + this_run->nnode;
  int n = this_run->npart;

  float *xsamp_all = (float *) my_malloc( sizeof(float) * nmem);
  float *ysamp_all = (float *) my_malloc( sizeof(float) * nmem);
  float *zsamp_all = (float *) my_malloc( sizeof(float) * nmem);

  int nsamp_local = (int)(nsamp * this_run->nrate);
  nrate = n / nsamp_local + 1;
  nsamp = samplingParticle( particle, n, nrate, 
			    xsamp_all, ysamp_all, zsamp_all, this_run);
  assert( nsamp <= nmem);

  if( this_run->inode == 0){
    fprintf( stderr, "writing %s nsamp:%d\n", filename, nsamp);
    FILE *fout = fopen( filename, "wb");
    fwrite( &nsamp, sizeof(int), 1, fout);
    fwrite( xsamp_all, sizeof(float), nsamp, fout);
    fwrite( ysamp_all, sizeof(float), nsamp, fout);
    fwrite( zsamp_all, sizeof(float), nsamp, fout);
    fclose(fout);
  }


  free(xsamp_all);
  free(ysamp_all);
  free(zsamp_all);


}



/* determine own area using sampling particle */
void determineBoundary( pParticle particle, const int n, 
			const long long int npart_all, const int *ndiv, 
			double bmin[][3], double bmax[][3],
			RunParam *this_run){
  double nowtime;
  double time1 = 0.0;
  getTime( &nowtime);

  int nrate = NRATE_EXCHANGE;
  int nsamp = npart_all / nrate + 1;
  int nmem = nsamp + this_run->nnode;

#ifdef STATIC_ARRAY
  static float xsamp_all[NSAMPMAX];
  static float ysamp_all[NSAMPMAX];
  static float zsamp_all[NSAMPMAX];
  nmem = NSAMPMAX;
#else
  float *xsamp_all = new float[nmem];
  float *ysamp_all = new float[nmem];
  float *zsamp_all = new float[nmem];
#endif

  int nsamp_local = (int)(nsamp * this_run->nrate);
  nrate = n / nsamp_local + 1;

  double bminbuf[3], bmaxbuf[3];
#ifdef NEW_DECOMPOSITION
  static int nstep0 = 0;
  if( nstep0 < nstep_decompose_x){
    nsamp = samplingParticle( particle, n, nrate, 
			      xsamp_all, ysamp_all, zsamp_all, this_run);
    My_MPI_Barrier( this_run->MPI_COMM_INTERNAL);
    time1 += getTime( &nowtime);
    assert( nsamp <= nmem);
    determineBoundary( xsamp_all, ysamp_all, zsamp_all, nsamp,
		       this_run, bminbuf, bmaxbuf);
  }
  else{
    int nsamp_x = 0;
    int sampx_on = 0;
    if( nstep0 % nstep_decompose_x == 0){
      nsamp_x = samplingParticleX( particle,n, nrate, xsamp_all, this_run);
      sampx_on = 1;
    }
    int nsamp_yz = samplingParticleYZ( particle,n, nrate, ysamp_all, zsamp_all, this_run);
    determineBoundary2( xsamp_all, ysamp_all, zsamp_all, nsamp_x, nsamp_yz,
			this_run, bminbuf, bmaxbuf, sampx_on);
  }
  nstep0 ++;
#else
  nsamp = samplingParticle( particle, n, nrate, 
			    xsamp_all, ysamp_all, zsamp_all, this_run);
  My_MPI_Barrier( this_run->MPI_COMM_INTERNAL);
  time1 += getTime( &nowtime);
  assert( nsamp <= nmem);

  determineBoundary( xsamp_all, ysamp_all, zsamp_all, nsamp,
		     this_run, bminbuf, bmaxbuf);
#endif

#ifdef BOUNDARY_SMOOTHING
  static int nstep = 0;
  if( this_run->restart_flag == 1){
    if( nstep < (nstep_smoothing_boundary-1)){
      bminbuf[0] = this_run->bmin[0];
      bminbuf[1] = this_run->bmin[1];
      bminbuf[2] = this_run->bmin[2];
      bmaxbuf[0] = this_run->bmax[0];
      bmaxbuf[1] = this_run->bmax[1];
      bmaxbuf[2] = this_run->bmax[2];
    }
  }
#endif

  My_MPI_Barrier( this_run->MPI_COMM_INTERNAL);
  getTime( &nowtime);
  MPI_Allgather( bminbuf, 3, MPI_DOUBLE, bmin, 3, MPI_DOUBLE, this_run->MPI_COMM_INTERNAL);
  MPI_Allgather( bmaxbuf, 3, MPI_DOUBLE, bmax, 3, MPI_DOUBLE, this_run->MPI_COMM_INTERNAL);
  My_MPI_Barrier( this_run->MPI_COMM_INTERNAL);
  time1 += getTime( &nowtime);

  this_run->t.exchange_determine_boundary_comm = time1;


#ifndef STATIC_ARRAY
  delete [] xsamp_all;
  delete [] ysamp_all;
  delete [] zsamp_all;
#endif

#ifdef BOUNDARY_SMOOTHING
  static double t_bmin[nstep_smoothing_boundary][3];
  static int sstep = 0;
  t_bmin[sstep][0] = bminbuf[0];
  t_bmin[sstep][1] = bminbuf[1];
  t_bmin[sstep][2] = bminbuf[2];
  if( nstep >= (nstep_smoothing_boundary-1)){
    double numerator[3] = {0.0, 0.0, 0.0};
    for( int i=0; i<nstep_smoothing_boundary; i++){
      int h = sstep - i;
      if( h < 0)  h += nstep_smoothing_boundary;
      double w = (double)( nstep_smoothing_boundary - i);
      for( int j=0; j<3; j++){
	numerator[j] += t_bmin[h][j] * w;
      }
    }
    double fac = 2.0 / (double)nstep_smoothing_boundary / (double)(nstep_smoothing_boundary+1);
    bminbuf[0] = numerator[0] * fac;
    bminbuf[1] = numerator[1] * fac;
    bminbuf[2] = numerator[2] * fac;

    int inode = this_run->inode;
    int px, py, pz;
    getXYZIndex( inode, ndiv, &px, &py, &pz);
    if( px == 0)  bmin[inode][0] = 0.0;
    if( py == 0)  bmin[inode][1] = 0.0;
    if( pz == 0)  bmin[inode][2] = 0.0;
    MPI_Allgather( bminbuf, 3, MPI_DOUBLE, bmin, 3, MPI_DOUBLE, this_run->MPI_COMM_INTERNAL);
    int px2 = px + 1;
    if( px2 >= ndiv[0])  bmaxbuf[0] = 1.0;
    else  bmaxbuf[0] = bmin[getVoxelIndex(px2, py, pz, ndiv)][0];
    int py2 = py + 1;
    if( py2 >= ndiv[1])  bmaxbuf[1] = 1.0;
    else  bmaxbuf[1] = bmin[getVoxelIndex(px, py2, pz, ndiv)][1];
    int pz2 = pz + 1;
    if( pz2 >= ndiv[2])  bmaxbuf[2] = 1.0;
    else  bmaxbuf[2] = bmin[getVoxelIndex(px, py, pz2, ndiv)][2];
    MPI_Allgather( bmaxbuf, 3, MPI_DOUBLE, bmax, 3, MPI_DOUBLE, this_run->MPI_COMM_INTERNAL);
    //printBoundary( this_run, stderr, bminbuf, bmaxbuf);
  }
  nstep ++;
  sstep ++;
  if( sstep == nstep_smoothing_boundary)  sstep = 0;
#endif

#ifdef BINARY_BOUNDARY
  double scale = (double)(1LL<<16);
  for( int i=0; i<this_run->nnode; i++){
    for( int j=0; j<3; j++){
      bmin[i][j] = (double)((long long int)( bmin[i][j] * scale)) / scale;
      bmax[i][j] = (double)((long long int)( bmax[i][j] * scale)) / scale;
    }
  }
#endif


}



void exchangeParticle( pParticle particle, pRunParam this_run, 
		       double bmin_all[][3], double bmax_all[][3]){

  double nowtime = 0.0;
  double time1 = 0.0;
  double time2 = 0.0;
  double time4 = 0.0;
  getTime(&nowtime);

  int n = this_run->npart;
  long long int npart_all = this_run->npart_total;
  int inode = this_run->inode;
  int nnode = this_run->nnode;
  MPI_Datatype MPI_PARTICLE;
  createMPIParticle( &MPI_PARTICLE);

  fprintf_verbose( stderr, "#exchangeParticle sortInBoxParticle start\n");
  int oldn = sortInBoxParticle( particle, n, bmin_all[inode], bmax_all[inode]);
  int dn = n - oldn;   // number of particle to send
  time4 += getTime( &nowtime);

  if( BOUNDARYLOG_ON == 1){
    sprintf( cbuf, "%e\t%e\t%e\t%e\t%e\t%e\t%d\t%d\n",
	     this_run->bmin[0], this_run->bmax[0], 
	     this_run->bmin[1], this_run->bmax[1], 
	     this_run->bmin[2], this_run->bmax[2],
	     oldn, n);
    mp_print( cbuf, this_run->MPI_COMM_INTERNAL, this_run->boundary_file);
  }

#ifdef USING_MPI_PARTICLE
#ifdef STATIC_ARRAY
  if( dn >= bufsize_largemem0){
    cerr << inode << "\t" << dn << "\t" << bufsize_largemem0 << endl;
  }
  assert( dn < bufsize_largemem0);
  static Particle sendbuf[bufsize_largemem0];
#else  // STATIC_ARRAY
  pParticle sendbuf = (pParticle) my_malloc( sizeof(Particle)*dn + 1);
#endif // STATIC_ARRAY

#else  // USING_MPI_PARTICLE
  double *sendbuf1 = new double[3*dn+1];
  float  *sendbuf2 = new float[3*dn+1];
  IDTYPE *sendbuf3 = new IDTYPE[dn+1];
#endif  // USING_MPI_PARTICLE

  static int nsend[MAXNNODE];
  static int sdispls[MAXNNODE];
  static int nrecv[MAXNNODE];
  static int rdispls[MAXNNODE];
  for( int i=0; i<nnode; i++){
    nsend[i] = sdispls[i] = nrecv[i] = rdispls[i] = 0;
  }

  fprintf_verbose( stderr, "#exchangeParticle particleSortForExchange start\n");
  particleSortForExchange( particle, this_run, dn, bmax_all, nsend);
  for( int i=oldn; i<n; i++){
#ifdef USING_MPI_PARTICLE
    sendbuf[i-oldn] = particle[i];
#else
    sendbuf1[3*(i-oldn)+0] = particle[i].xpos;
    sendbuf1[3*(i-oldn)+1] = particle[i].ypos;
    sendbuf1[3*(i-oldn)+2] = particle[i].zpos;
    sendbuf2[3*(i-oldn)+0] = particle[i].xvel;
    sendbuf2[3*(i-oldn)+1] = particle[i].yvel;
    sendbuf2[3*(i-oldn)+2] = particle[i].zvel;
    sendbuf3[i-oldn]       = particle[i].id;
#endif
  }

  /*
  int offset = 0;
  for( i=0; i<nnode; i++){
    //fprintf( stderr, "%d %d %d\n", inode, nsend[i], dn); 
    int j;
    for( j=0; j<nsend[i]; j++){
      if( sendbuf[offset+j].xpos >= bmax_all[i][0] || 
	  sendbuf[offset+j].ypos >= bmax_all[i][1] || 
	  sendbuf[offset+j].zpos >= bmax_all[i][2] ||
	  sendbuf[offset+j].xpos <  bmin_all[i][0] ||
	  sendbuf[offset+j].ypos <  bmin_all[i][1] ||
	  sendbuf[offset+j].zpos <  bmin_all[i][2]){
	fprintf( stderr, "%d\t%d\t%d\t%d\t%e\t%e\t%e\t%lld\t%e\t%e\n", 
		 inode, i, j, nsend[i],
		 sendbuf[offset+j].xpos, sendbuf[offset+j].ypos, 
		 sendbuf[offset+j].zpos, sendbuf[offset+j].id,
		 bmin_all[i][2], bmax_all[i][2]);
      }
    }
    offset += nsend[i];
  }
  */
  for( int i=1; i<nnode; i++){
    sdispls[i] = sdispls[i-1] + nsend[i-1];
  }
  My_MPI_Barrier(this_run->MPI_COMM_INTERNAL);
  time1 += getTime( &nowtime);

  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);
  int nrecv_total = nrecv[0];
  for( int i=1; i<nnode; i++){
    rdispls[i] = rdispls[i-1] + nrecv[i-1];
    nrecv_total += nrecv[i];
  }

#ifndef USING_MPI_PARTICLE
  double *recvbuf1 = new double[nrecv_total*3];
  float  *recvbuf2 = new float[nrecv_total*3];
  IDTYPE *recvbuf3 = new IDTYPE[nrecv_total];
  static int nsend2[MAXNNODE];
  static int sdispls2[MAXNNODE];
  static int nrecv2[MAXNNODE];
  static int rdispls2[MAXNNODE];
  for( int i=0; i<nnode; i++){
    nsend2[i]   = nsend[i] * 3;
    sdispls2[i] = sdispls[i] * 3;
    nrecv2[i]   = nrecv[i] * 3;
    rdispls2[i] = rdispls[i] * 3;
  }
#endif

  int n_scomm = 0;
  int n_rcomm = 0;
  static int scomm_inode[MAXNNODE];
  static int rcomm_inode[MAXNNODE];
  for( int i=0; i<nnode; i++){
    if( nsend[i] != 0){
      scomm_inode[n_scomm] = i;
      n_scomm ++;
    }
    if( nrecv[i] != 0){
      rcomm_inode[n_rcomm] = i;
      n_rcomm ++;
    }
  }
  this_run->nscomm_exchange = n_scomm;
  this_run->nrcomm_exchange = n_rcomm;
  this_run->nsend_ex = dn;
  this_run->nrecv_ex = nrecv_total;
  fprintf_verbose( stderr, "#exchangeParticle comm start\n");

#ifdef EXCHANGE_COMM_NONBLOCKING
  static MPI_Request req_send[MAXNNODE3];
  static MPI_Request req_recv[MAXNNODE3];
  static MPI_Status  st_send[MAXNNODE3];
  static MPI_Status  st_recv[MAXNNODE3];
  //static int flag_send[MAXNNODE3];
  //static int flag_recv[MAXNNODE3];

#ifdef USING_MPI_PARTICLE
  for( int i=0; i<n_rcomm; i++){
    int src = rcomm_inode[i];
    int e = MPI_Irecv( &particle[oldn+rdispls[src]], nrecv[src], MPI_PARTICLE, src, 0, this_run->MPI_COMM_INTERNAL, req_recv+i);
    assert( e == MPI_SUCCESS);
  }
  for( int i=0; i<n_scomm; i++){
    int dest = scomm_inode[i];
    MPI_Isend( &sendbuf[sdispls[dest]], nsend[dest], MPI_PARTICLE, dest, 0, this_run->MPI_COMM_INTERNAL, req_send+i);
    //    MPI_Testall( n_rcomm, req_recv, flag_recv, st_recv);
    //    MPI_Testall( i+1,  req_send, flag_send, st_send);
  }
  MPI_Waitall( n_scomm, req_send, st_send);
  MPI_Waitall( n_rcomm, req_recv, st_recv);
#else // USING_MPI_PARTICLE
  for( int i=0; i<n_rcomm; i++){
    int src = rcomm_inode[i];
    MPI_Irecv( &recvbuf1[rdispls2[src]], nrecv2[src], MPI_DOUBLE, src, 0, this_run->MPI_COMM_INTERNAL, req_recv+3*i);
    MPI_Irecv( &recvbuf2[rdispls2[src]], nrecv2[src], MPI_FLOAT, src, 1, this_run->MPI_COMM_INTERNAL, req_recv+3*i+1);
    MPI_Irecv( &recvbuf3[rdispls[src]], nrecv[src], MPI_IDTYPE, src, 2, this_run->MPI_COMM_INTERNAL, req_recv+3*i+2);
  }
  for( int i=0; i<n_scomm; i++){
    int dest = scomm_inode[i];
    MPI_Isend( &sendbuf1[sdispls2[dest]], nsend2[dest], MPI_DOUBLE, dest, 0, this_run->MPI_COMM_INTERNAL, req_send+3*i);
    MPI_Isend( &sendbuf2[sdispls2[dest]], nsend2[dest], MPI_FLOAT, dest, 1, this_run->MPI_COMM_INTERNAL, req_send+3*i+1);
    MPI_Isend( &sendbuf3[sdispls[dest]], nsend[dest], MPI_IDTYPE, dest, 2, this_run->MPI_COMM_INTERNAL, req_send+3*i+2);
    //    MPI_Testall( n_rcomm*3, req_recv, flag_recv, st_recv);
    //    MPI_Testall( 3*(i+1),  req_send, flag_send, st_send);
  }
  MPI_Waitall( 3*n_scomm, req_send, st_send);
  MPI_Waitall( 3*n_rcomm, req_recv, st_recv);
#endif // USING_MPI_PARTICLE

#else //EXCHANGE_COMM_NONBLOCKING

#ifndef EXCHANGE_COMM_SENDRECV
#ifdef USING_MPI_PARTICLE
  MPI_Alltoallv( sendbuf, nsend, sdispls, MPI_PARTICLE, 
		 &particle[oldn], nrecv, rdispls, MPI_PARTICLE,
		 this_run->MPI_COMM_INTERNAL);
#else
  MPI_Alltoallv( sendbuf3, nsend, sdispls, MPI_IDTYPE,
		 recvbuf3, nrecv, rdispls, MPI_IDTYPE,
		 this_run->MPI_COMM_INTERNAL);
  MPI_Alltoallv( sendbuf2, nsend2, sdispls2, MPI_FLOAT,
		 recvbuf2, nrecv2, rdispls2, MPI_FLOAT,
		 this_run->MPI_COMM_INTERNAL);
  MPI_Alltoallv( sendbuf1, nsend2, sdispls2, MPI_DOUBLE,
		 recvbuf1, nrecv2, rdispls2, MPI_DOUBLE,
		 this_run->MPI_COMM_INTERNAL);
#endif 
#else //EXCHANGE_COMM_SENDRECV
  MPI_Status stat;
  int inode_src = (inode + 1) % nnode;
  for( int k=(nnode-1); k>0; k--){
    int inode_dest = (k+inode) % nnode;
    if( k != (nnode -1)){
      inode_src = (inode_src+1) % nnode;
      if( inode_src == inode)  inode_src = (inode_src+1) % nnode;
    }
#ifdef USING_MPI_PARTICLE
    MPI_Sendrecv( &sendbuf[sdispls[inode_dest]], nsend[inode_dest],
		  MPI_PARTICLE, inode_dest, inode,
		  &particle[oldn+rdispls[inode_src]], nrecv[inode_src],
		  MPI_PARTICLE, inode_src, inode_src,
		  this_run->MPI_COMM_INTERNAL, &stat);
#else
    MPI_Sendrecv( &sendbuf1[sdispls2[inode_dest]], nsend2[inode_dest],
		  MPI_DOUBLE, inode_dest, inode,
		  &recvbuf1[rdispls2[inode_src]], nrecv2[inode_src],
		  MPI_DOUBLE, inode_src, inode_src,
		  this_run->MPI_COMM_INTERNAL, &stat);
    MPI_Sendrecv( &sendbuf2[sdispls2[inode_dest]], nsend2[inode_dest],
		  MPI_FLOAT, inode_dest, inode,
		  &recvbuf2[rdispls2[inode_src]], nrecv2[inode_src],
		  MPI_FLOAT, inode_src, inode_src,
		  this_run->MPI_COMM_INTERNAL, &stat);
    MPI_Sendrecv( &sendbuf3[sdispls[inode_dest]], nsend[inode_dest],
		  MPI_IDTYPE, inode_dest, inode,
		  &recvbuf3[rdispls[inode_src]], nrecv[inode_src],
		  MPI_IDTYPE, inode_src, inode_src,
		  this_run->MPI_COMM_INTERNAL, &stat);
#endif
  }
#endif //EXCHANGE_COMM_SENDRECV
#endif //EXCHANGE_COMM_NONBLOCKING

  fprintf_verbose( stderr, "#exchangeParticle push start\n");

#ifndef USING_MPI_PARTICLE
  for( int i=0; i<nrecv_total; i++){
    particle[i+oldn].xpos = recvbuf1[i*3+0];
    particle[i+oldn].ypos = recvbuf1[i*3+1];
    particle[i+oldn].zpos = recvbuf1[i*3+2];
    particle[i+oldn].xvel = recvbuf2[i*3+0];
    particle[i+oldn].yvel = recvbuf2[i*3+1];
    particle[i+oldn].zvel = recvbuf2[i*3+2];
    particle[i+oldn].id = recvbuf3[i];
  }
#endif

  this_run->npart = oldn + nrecv_total;
#ifdef USING_MPI_PARTICLE
#ifndef STATIC_ARRAY
  free( sendbuf);
#endif
#else
  delete [] sendbuf1;
  delete [] sendbuf2;
  delete [] sendbuf3;
  delete [] recvbuf1;
  delete [] recvbuf2;
  delete [] recvbuf3;
#endif

  //  fprintf( stderr, "node:%d\toldn:%d\tn:%d\n", inode, oldn, this_run->npart);

  /* for protection */
  long long int lli_npart_all = (long long int)npart_all;
  long long int npart_now = 0;
  long long int this_npart = (long long int)this_run->npart;
  MPI_Allreduce( &this_npart, &npart_now, 1, MPI_LONG_LONG_INT, MPI_SUM, this_run->MPI_COMM_INTERNAL);
  if( npart_now != npart_all){
#ifdef NPART_DIFFERENT_DUMP
    sprintf( cbuf, "npartdiff-%d", inode);
    ofstream fout( cbuf, ios::out);
    cerr << oldn << "\t" << this_run->npart << endl;
    for( int i=oldn; i<this_run->npart; i++){
      fout << particle[i].xpos << "\t"
	   << particle[i].ypos << "\t"
	   << particle[i].zpos << "\t"
	   << particle[i].xvel << "\t"
	   << particle[i].yvel << "\t"
	   << particle[i].zvel << "\t"
	   << particle[i].id << endl;
    }
    fout.close();
    MPI_Barrier( MPI_COMM_WORLD);
#endif
    fprintf( stderr, "npart is different\n");
    fprintf( stderr, "%lld\t%lld\n", lli_npart_all, npart_now);
    exit(0);
  }
  time2 += getTime( &nowtime);

  this_run->t.exchange_make_send_list = time1;
  this_run->t.exchange_make_send_list_pre = time4;
  this_run->t.exchange_comm = time2;


}




/* exchange particle with neighboring node */
void exchangeParticle( pParticle particle, pRunParam this_run, const int *ndiv){

  /*                
   * particle array / old_n              \/ bn        \
   *                ------------------------------------------------
   *                \no change node      /\send buffer/\recv buffer/
   * finished communication
   *                ______________________abcde________edcba
   * move recv buffer particle to send buffer
   * particle array ___________________________
   *                \ new particle array      /
   */

  double nowtime = 0.0;
  getTime(&nowtime);

  int n = this_run->npart;
  long long int npart_all = this_run->npart_total;
  this_run->npart_old = n;
  int inode = this_run->inode;

  fprintf_verbose( stderr, "#exchangeParticle determineBoundary start\n");
#ifndef FIXED_BOUNDARY
  determineBoundary( particle, n, npart_all, ndiv, this_run->bmin_all, 
		     this_run->bmax_all,
		     this_run);
#endif
  this_run->bmin[0] = this_run->bmin_all[inode][0];
  this_run->bmin[1] = this_run->bmin_all[inode][1];
  this_run->bmin[2] = this_run->bmin_all[inode][2];
  this_run->bmax[0] = this_run->bmax_all[inode][0];
  this_run->bmax[1] = this_run->bmax_all[inode][1];
  this_run->bmax[2] = this_run->bmax_all[inode][2];

  this_run->t.exchange_determine_boundary = getTime( &nowtime);

  exchangeParticle( particle, this_run, this_run->bmin_all, this_run->bmax_all);

}



/* make super particle list */
#ifndef TREE2
void makeSendList( const pTreeTmp tree_tmp, const double cell_length2,
		   const double *cicell, const double *i_half_length,
		   float send_buf[][4], int *nsend,
		   const long long int current_key,
		   pParticle particle, pTreeParam treeparam){

  int i, j;
  int adr = keyToAdr( current_key, tree_tmp->tree_cell, tree_tmp->tree_clistmask);
  double rcut2 = RADIUS_FOR_PP2;

  double dr = 0.0;
#pragma loop simd
  for( j=0; j<3; j++){
    double drj = cicell[j] - tree_tmp->tree_cell[adr].cm[j];
    drj = fabs(drj);
    if( drj > 0.5)  drj = 1.0 - drj;
    drj -= i_half_length[j];
    if( drj > 0)  dr += drj * drj;
  }

  double theta2_dr2 = treeparam->tree_theta2 * dr;
  /*fprintf( stdout, "%e\t%e\t%d\t%lld\t%e\t%e\n", 
    dr, sqrt(rcut2), tree_tmp->tree_cell[adr].n, current_key, sqrt(cell_length2),
    sqrt(theta2_dr2));*/
  if( theta2_dr2 > cell_length2){   // theta > d/l
    if( dr <= rcut2){
      send_buf[*nsend][0] = tree_tmp->tree_cell[adr].cm[0];
      send_buf[*nsend][1] = tree_tmp->tree_cell[adr].cm[1];
      send_buf[*nsend][2] = tree_tmp->tree_cell[adr].cm[2];
      send_buf[*nsend][3] = tree_tmp->tree_cell[adr].m;
      (*nsend) ++;
    }
  }
  else{
    if( tree_tmp->tree_cell[adr].n > treeparam->tree_nleaf){  // have child
      for( i=0; i<8; i++){
	int zeroflag = ( tree_tmp->tree_cell[adr].zeroflag>>i) & 1;
	if( zeroflag == 1)  continue;
	long long int tmpkey = ((current_key<<3) | i);
	makeSendList( tree_tmp, cell_length2*0.25, cicell, i_half_length,
		      send_buf, nsend, 
		      tmpkey, particle, treeparam);
      }
    }
    else{  // lowest cell
      int start = tree_tmp->tree_cell[adr].first_index;
      int end = start + tree_tmp->tree_cell[adr].n;
      for( j=start; j<end; j++){
	int ii = tree_tmp->index[j];
#ifdef CLEAN_BOUNDARY_PARTICLE
	double dx = cicell[0] - particle[ii].xpos;
	double dy = cicell[1] - particle[ii].ypos;
	double dz = cicell[2] - particle[ii].zpos;
	dx = fabs( dx);
	dy = fabs( dy);
	dz = fabs( dz);
	if( dx > 0.5)  dx = 1.0 - dx;
	if( dy > 0.5)  dy = 1.0 - dy;
	if( dz > 0.5)  dz = 1.0 - dz;
	dx -= i_half_length[0];
	dy -= i_half_length[1];
	dz -= i_half_length[2];
	dx += fabs( dx);
	dy += fabs( dy);
	dz += fabs( dz);
	double dr2 = 0.25 * (dx*dx + dy*dy + dz*dz);
	if( dr2 <= rcut2){
	  send_buf[*nsend][3] = getMass( &particle[ii], treeparam->uniform_mass);
	  getPos2( &particle[ii], send_buf[*nsend]);
	  (*nsend) ++;
	}
#else
	send_buf[*nsend][3] = getMass( &particle[ii], treeparam->uniform_mass);
	getPos2( &particle[ii], send_buf[*nsend]);
	(*nsend) ++;
#endif
      }
    }
  }

}
#else
/* make super particle list */
void makeSendList( const pTreeTmp tree_tmp, const double cell_length2,
		   const double *cicell, const double *i_half_length,
		   float send_buf[][4], int *nsend,
		   const int now_cell,
		   const pTreeParam treeparam,
		   const double (*p_cache)[4]){

  double rcut2 = RADIUS_FOR_PP2;
  pTreeCell tree_cell = tree_tmp->tree_cell;

  if( tree_cell[now_cell].n > treeparam->tree_nleaf){
    double child_cell_length2 = cell_length2 * 0.25;

    double dr2[8];
    double theta2_dr2[8];
    for( int i=0; i<8; i++){
      int adr = tree_cell[now_cell].next + i;
      double dx = cicell[0] - tree_cell[adr].cm[0];
      double dy = cicell[1] - tree_cell[adr].cm[1];
      double dz = cicell[2] - tree_cell[adr].cm[2];
      dx = fabs( dx);
      dy = fabs( dy);
      dz = fabs( dz);
      if( dx > 0.5)  dx = 1.0 - dx;
      if( dy > 0.5)  dy = 1.0 - dy;
      if( dz > 0.5)  dz = 1.0 - dz;
      dx -= i_half_length[0];
      dy -= i_half_length[1];
      dz -= i_half_length[2];
      dx += fabs( dx);
      dy += fabs( dy);
      dz += fabs( dz);
      dr2[i] = 0.25 * (dx*dx + dy*dy + dz*dz);
      theta2_dr2[i] = treeparam->tree_theta2 * dr2[i];
    }
    for( int i=0; i<8; i++){
      int adr = tree_cell[now_cell].next + i;
      if( tree_cell[adr].n == 0)  continue;
      if( theta2_dr2[i] > child_cell_length2){   // theta > d/l
	if( dr2[i] > rcut2)  continue;
	int ji = *nsend;
	send_buf[ji][0] = tree_tmp->tree_cell[adr].cm[0];
	send_buf[ji][1] = tree_tmp->tree_cell[adr].cm[1];
	send_buf[ji][2] = tree_tmp->tree_cell[adr].cm[2];
	send_buf[ji][3] = tree_tmp->tree_cell[adr].m;
	(*nsend) ++;
      }
      else{
	makeSendList( tree_tmp, child_cell_length2, cicell, i_half_length,
		      send_buf, nsend, adr, treeparam, p_cache);
      }
    }
  }
  else{
    const int start = tree_cell[now_cell].first_index;
    const int n = tree_cell[now_cell].n;
    const int nj = *nsend;
    for( int j=0; j<n; j++){
      int ji = nj + j;
      int pi = start + j;
      send_buf[ji][0] = p_cache[pi][0];
      send_buf[ji][1] = p_cache[pi][1];
      send_buf[ji][2] = p_cache[pi][2];
      send_buf[ji][3] = p_cache[pi][3];
    }
    (*nsend) += n;
  }


}
#endif



void makeSendListWrapper( const pTreeTmp tree_tmp, const double cell_length2,
			  const double *dest_bmin, const double *dest_bmax,
			  float (*send_buf)[4], int *nsend,
			  pParticle particle, pTreeParam treeparam,
			  const double (*p_cache)[4]){

  double dest_center[3], dest_bsize[3];
  for( int j=0; j<3; j++){
    dest_center[j] = 0.5 * (dest_bmax[j] + dest_bmin[j]);
    dest_bsize[j]  = 0.5 * (dest_bmax[j] - dest_bmin[j]);
  }

#ifdef TREE2
  makeSendList( tree_tmp, cell_length2, dest_center, dest_bsize,
		send_buf, nsend, 0, treeparam, p_cache);
#else
  makeSendList( tree_tmp, cell_length2, dest_center, dest_bsize,
		send_buf, nsend, 1, particle, treeparam);
#endif

}




#ifndef TREE2
/*get the number of send list */
int getSendNumber( const pTreeTmp tree_tmp, const double cell_length2,
		   const double *cicell, const double *i_half_length,
		   const long long int current_key,
		   pParticle particle, pTreeParam treeparam){

  int nsend = 0;

  int i, j;
  int adr = keyToAdr( current_key, tree_tmp->tree_cell, tree_tmp->tree_clistmask);
  double rcut2 = RADIUS_FOR_PP2;

  double dr = 0.0;
#pragma loop simd
  for( j=0; j<3; j++){
    double drj = cicell[j] - tree_tmp->tree_cell[adr].cm[j];
    drj = fabs(drj);
    if( drj > 0.5)  drj = 1.0 - drj;
    drj -= i_half_length[j];
    if( drj > 0)  dr += drj * drj;
  }

  double theta2_dr2 = treeparam->tree_theta2 * dr;
  if( theta2_dr2 > cell_length2){   // theta > d/l
    if( dr <= rcut2){
      nsend ++;
    }
  }
  else{
    if( tree_tmp->tree_cell[adr].n > treeparam->tree_nleaf){  // have child
      for( i=0; i<8; i++){
	int zeroflag = ( tree_tmp->tree_cell[adr].zeroflag>>i) & 1;
	if( zeroflag == 1)  continue;
	long long int tmpkey = ((current_key<<3) | i);
	int nsend2 = getSendNumber( tree_tmp, cell_length2*0.25, 
				    cicell, i_half_length, 
				    tmpkey, particle, treeparam);
	nsend += nsend2;
      }
    }
    else{  // lowest cell
#ifdef CLEAN_BOUNDARY_PARTICLE
      for( int i=0; i<tree_tmp->tree_cell[adr].n; i++){
	int ii = tree_tmp->tree_cell[adr].first_index + i;
	int ip = tree_tmp->index[ii];
	double dx = cicell[0] - particle[ip].xpos;
	double dy = cicell[1] - particle[ip].ypos;
	double dz = cicell[2] - particle[ip].zpos;
	dx = fabs( dx);
	dy = fabs( dy);
	dz = fabs( dz);
	if( dx > 0.5)  dx = 1.0 - dx;
	if( dy > 0.5)  dy = 1.0 - dy;
	if( dz > 0.5)  dz = 1.0 - dz;
	dx -= i_half_length[0];
	dy -= i_half_length[1];
	dz -= i_half_length[2];
	dx += fabs( dx);
	dy += fabs( dy);
	dz += fabs( dz);
	double dr2 = 0.25 * (dx*dx + dy*dy + dz*dz);
	if( dr2 <= rcut2) nsend ++;
      }
#else
      nsend += tree_tmp->tree_cell[adr].n;
#endif
    }
  } 

  return nsend;

}
#else
int getSendNumber( const pTreeTmp tree_tmp, const double cell_length2,
		   const double *cicell, const double *i_half_length,
		   const int now_cell,
		   const pTreeParam treeparam){

  int nsend = 0;

  double rcut2 = RADIUS_FOR_PP2;
  pTreeCell tree_cell = tree_tmp->tree_cell;

  if( tree_cell[now_cell].n > treeparam->tree_nleaf){
    double child_cell_length2 = cell_length2 * 0.25;
    double dr2[8];
    double theta2_dr2[8];
    for( int i=0; i<8; i++){
      int adr = tree_cell[now_cell].next + i;
      double dx = cicell[0] - tree_cell[adr].cm[0];
      double dy = cicell[1] - tree_cell[adr].cm[1];
      double dz = cicell[2] - tree_cell[adr].cm[2];
      dx = fabs( dx);
      dy = fabs( dy);
      dz = fabs( dz);
      if( dx > 0.5)  dx = 1.0 - dx;
      if( dy > 0.5)  dy = 1.0 - dy;
      if( dz > 0.5)  dz = 1.0 - dz;
      dx -= i_half_length[0];
      dy -= i_half_length[1];
      dz -= i_half_length[2];
      dx += fabs( dx);
      dy += fabs( dy);
      dz += fabs( dz);
      dr2[i] = 0.25 * (dx*dx + dy*dy + dz*dz);
      theta2_dr2[i] = treeparam->tree_theta2 * dr2[i];
    }
    for( int i=0; i<8; i++){
      int adr = tree_cell[now_cell].next + i;
      if( tree_cell[adr].n == 0)  continue;
      if( theta2_dr2[i] > child_cell_length2){   // theta > d/l
	if( dr2[i] <= rcut2)  nsend ++;
      }
      else{
	nsend += getSendNumber( tree_tmp, child_cell_length2, 
				cicell, i_half_length,
				adr, treeparam);
      }
    }
  }
  else{
    nsend += tree_cell[now_cell].n;
  }


  return nsend;

}
#endif



int getSendNumberWrapper( const pTreeTmp tree_tmp, const double cell_length2,
			  const double *dest_bmin, const double *dest_bmax,
			  pParticle particle, pTreeParam treeparam){

  double dest_center[3], dest_bsize[3];
  for( int j=0; j<3; j++){
    dest_center[j] = 0.5 * (dest_bmax[j] + dest_bmin[j]);
    dest_bsize[j]  = 0.5 * (dest_bmax[j] - dest_bmin[j]);
  }

#ifdef TREE2
  int nsend = getSendNumber( tree_tmp, cell_length2,
			     dest_center, dest_bsize,
			     0, treeparam);
#else
  int nsend = getSendNumber( tree_tmp, cell_length2,
			     dest_center, dest_bsize,
			     1, particle, treeparam);
#endif

    return nsend;


}



/* communicate boundary particle using super-particle */
static float send_buf[bufsize_largemem][4];
static float recv_buf[bufsize_largemem][4];
int sendRecvBoundaryUsingTree( pParticle particle, pTreeParam treeparam,
			       const int n,
			       double bmin_all[][3], double bmax_all[][3],
			       const int *comm_flag,
			       pRunParam this_run, pTreeTmp tree_tmp){

  double nowtime = 0.0;
  double time1 = 0.0;
  double time2 = 0.0;
  double time4 = 0.0;
  getTime(&nowtime);

  int nnode = this_run->nnode;
  double maxx = tree_tmp->maxx;

  // send and receive
  static int nsend[MAXNNODE];
  static int sdispls[MAXNNODE];
  static int nrecv[MAXNNODE];
  static int rdispls[MAXNNODE];
  long long int nsend_total = 0;
  long long int nrecv_total = 0;
  for( int i=0; i<nnode; i++){
    nsend[i] = sdispls[i] = rdispls[i] = 0;
  }

  fprintf_verbose( stderr, "#getSendNumber\n");

  getTime( &nowtime);
#pragma omp parallel for CHUNK_DECOMPOSITION
  for( int i=0; i<nnode; i++){
    if( comm_flag[i] == 0)  continue;
    nsend[i] = getSendNumberWrapper( tree_tmp, maxx*maxx, 
				     this_run->bmin_all[i], this_run->bmax_all[i], 
				     particle, treeparam);
  }
  for( int i=0; i<nnode; i++)  nsend_total += nsend[i];
  assert( nsend_total < bufsize_largemem);
  for( int i=1; i<nnode; i++)  sdispls[i] = sdispls[i-1] + nsend[i-1];

  fprintf_verbose( stderr, "#comm SendNumber\n");

  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);
  nrecv_total = nrecv[0];
  for( int i=1; i<nnode; i++){
    rdispls[i] = rdispls[i-1] + nrecv[i-1];
    nrecv_total += nrecv[i];
  }
  assert( nrecv_total < bufsize_largemem);

#ifdef BOUNDARY_COMM_NONBLOCKING
  fprintf_verbose( stderr, "#start comm\n");
  static int comm_inode[MAXNNODE];
  int ncomm = 0;
  for( int i=0; i<nnode; i++){
    if( comm_flag[i] == 0)  continue;
    comm_inode[ncomm] = i;
    ncomm ++;
  }

  static MPI_Request req_send[MAXNNODE];
  static MPI_Request req_recv[MAXNNODE];
  static MPI_Status st[MAXNNODE];
  for( int i=0; i<ncomm; i++){
    int src = comm_inode[i];
    MPI_Irecv( recv_buf[rdispls[src]], 4*nrecv[src], MPI_FLOAT, src, 0, this_run->MPI_COMM_INTERNAL, req_recv+i);
  }

#ifndef BOUNDARY_COMM_NONBLOCKING_OVERLAP
#pragma omp parallel for CHUNK_DECOMPOSITION
  for( int i=0; i<ncomm; i++){
    int nsend0  = 0;
    int dest = comm_inode[i];
    makeSendListWrapper( tree_tmp, maxx*maxx, this_run->bmin_all[dest],
			 this_run->bmax_all[dest], &send_buf[sdispls[dest]], &nsend0,
			 particle, treeparam, p_cache);
    assert( nsend0 == nsend[dest]);
  }
  fprintf_verbose( stderr, "#end makeSendList\n");
  time1 += getTime(&nowtime);
  for( int i=0; i<ncomm; i++){
    int dest = comm_inode[i];
    MPI_Isend( send_buf[sdispls[dest]], 4*nsend[dest], MPI_FLOAT, dest, 0, this_run->MPI_COMM_INTERNAL, req_send+i);
    //MPI_Testall( ncomm, req_recv, flag_recv, st);
    //MPI_Testall( i+1, req_send, flag_send, st);
  }
#else  //BOUNDARY_COMM_NONBLOCKING_OVERLAP
  for( int i=0; i<ncomm; i++){
    int nsend0  = 0;
    int dest = comm_inode[i];
    makeSendListWrapper( tree_tmp, maxx*maxx, this_run->bmin_all[dest],
			 this_run->bmax_all[dest], &send_buf[sdispls[dest]], &nsend0,
			 particle, treeparam, p_cache);
    assert( nsend0 == nsend[dest]);
    MPI_Isend( send_buf[sdispls[dest]], 4*nsend0, MPI_FLOAT, dest, 0, this_run->MPI_COMM_INTERNAL, req_send+i);
    //MPI_Testall( ncomm, req_recv, flag_recv, st);
    //MPI_Testall( i+1, req_send, flag_send, st);
  }
#endif
  MPI_Waitall( ncomm, req_send, st);
  MPI_Waitall( ncomm, req_recv, st);
  My_MPI_Barrier(this_run->MPI_COMM_INTERNAL);
  fprintf_verbose( stderr, "#end comm\n");
  time2 += getTime(&nowtime);

#else //BOUNDARY_COMM_NONBLOCKING
  MPI_Status stat;
#pragma omp parallel for CHUNK_DECOMPOSITION
  for( int i=0; i<nnode; i++){
    nsend[i] = 0;
    if( comm_flag[i] != 0){
    makeSendListWrapper( tree_tmp, maxx*maxx, this_run->bmin_all[dest],
			 this_run->bmax_all[dest], &send_buf[sdispls[dest]], &nsend0,
			 particle, treeparam, p_cache);
    }
  }
  My_MPI_Barrier(this_run->MPI_COMM_INTERNAL);
  time1 += getTime(&nowtime);

  int inode = this_run->inode;
  int inode_src = (inode + 1) % nnode;
  for( int k=(nnode-1); k>0; k--){
    int inode_dest = (k+inode) % nnode;
    if( k != (nnode -1)){
      inode_src = (inode_src+1) % nnode;
      if( inode_src == inode)  inode_src = (inode_src+1) % nnode;
    }
    MPI_Sendrecv( send_buf[sdispls[inode_dest]], 4*nsend[inode_dest], MPI_FLOAT, inode_dest, inode, 
		  recv_buf[rdispls[inode_src]],  4*nrecv[inode_src] , MPI_FLOAT, inode_src , inode_src,
		  this_run->MPI_COMM_INTERNAL, &stat);
  }
  My_MPI_Barrier(this_run->MPI_COMM_INTERNAL);
  time2 += getTime(&nowtime);
#endif //BOUNDARY_COMM_NONBLOCKING

#pragma omp parallel for
  for( int i=0; i<nrecv_total; i++){
    int ip = n+i;
    particle[ip].xpos = recv_buf[i][0];
    particle[ip].ypos = recv_buf[i][1];
    particle[ip].zpos = recv_buf[i][2];
    particle[ip].xvel = recv_buf[i][3];
    particle[ip].id = -1;
  }
  time4 += getTime(&nowtime);

  this_run->t.pp_get_boundary_make_send_list = time1;
  this_run->t.pp_get_boundary_comm = time2;
  this_run->t.pp_get_boundary_push = time4;
  this_run->nsend_tree = nsend_total;

  return nrecv_total;

}



void determineCommNode( double bmin_all[][3], double bmax_all[][3], 
			int *comm_flag, pRunParam this_run){

  int nnode = this_run->nnode;
  int inode = this_run->inode;

  double icenter[3], ibsize[3], ibsize2[3];;
  for( int j=0; j<3; j++){
    icenter[j] = 0.5 * (bmax_all[inode][j] + bmin_all[inode][j]);
    ibsize[j]  = 0.5 * (bmax_all[inode][j] - bmin_all[inode][j]);
    ibsize2[j] = ibsize[j] + SFT_FOR_PM;
  }

  for( int i=0; i<nnode; i++){
    comm_flag[i] = 1;
    for( int j=0; j<3; j++){
      double center = 0.5 * ( bmax_all[i][j] + bmin_all[i][j] );
      double bsize  = 0.5 * ( bmax_all[i][j] - bmin_all[i][j] );
      double drj = icenter[j] - center;
      drj = fabs(drj);
      if( drj > 0.5)  drj = 1.0 - drj;
      double max_length_j = ibsize2[j] + bsize;
      if( drj > max_length_j)  comm_flag[i] = 0;
    }
  }
  comm_flag[inode] = 0;

  int nnode_comm = 0;
  for( int i=0; i<nnode; i++){
    if( comm_flag[i] == 1)  nnode_comm ++;
  }
  this_run->ncomm_tree = nnode_comm;

}



int getBoundaryParticle( pParticle particle, pRunParam this_run,
			 TreeParam *treeparam
#ifdef BUFFER_FOR_TREE
			 ,TreeTmp &tree_tmp
#endif
			 ){

  double nowtime = 0.0;
  getTime(&nowtime);

  int n = this_run->npart;
  double (*bmin_all)[3] = this_run->bmin_all;
  double (*bmax_all)[3] = this_run->bmax_all;

  static int comm_flag[MAXNNODE];
  determineCommNode( bmin_all, bmax_all, comm_flag, this_run);

  this_run->t.pp_get_boundary_determine_comm_node = getTime(&nowtime);

#ifndef BUFFER_FOR_TREE
  TreeTmp tree_tmp;
  TreeTmpTreeTmp( &tree_tmp, treeparam, n);
#else
  TreeTmpTreeTmp( &tree_tmp, n, 0);
#endif

  fprintf_verbose( stderr, "#PP getBoundaryParticle makeKey\n");
  tree_tmp.maxx = makeKey( particle, tree_tmp.key, n); 

  static int nstep = 0;
  if( nstep % SORT_STEP == 0){
#ifndef CALCPOT
    fprintf_verbose( stderr, "#PP getBoundaryParticle Sort start\n");
    qsortPivotWrapper( particle, tree_tmp.key, n, NLEVEL_PARTICLE_SORT, NCHUNK_QSORT);
    this_run->t.pp_particle_sort = getTime(&nowtime);
#endif
  }
  else{
    fprintf_verbose( stderr, "#PP getBoundaryParticle qsort\n");
    qsortPivotWrapper( tree_tmp.index, tree_tmp.key, n, NLEVEL_KEY_SORT, NCHUNK_QSORT);
    this_run->t.pp_get_boundary_make_key_sort = getTime(&nowtime);
  }
  nstep ++;

#ifdef TREE_PARTICLE_CACHE
  fprintf_verbose( stderr, "#PP getBoundaryParticle copyPCache\n");
  copyPCache( particle, tree_tmp.index, n, treeparam->uniform_mass);
#endif

  fprintf_verbose( stderr, "#PP getBoundaryParticle treeConstruction\n");

  int nwalk = 0;         // for i-particle
  int key_level = 63;
  int ipflag = 0;   //1:i-particle 0:upper 2:lower
#ifdef TREE2
  tree_tmp.tree_cell[0].first_index = 0;
  tree_tmp.tree_cell[0].n = n;
  tree_tmp.tree_cell[0].next = 1;
#ifdef TREECONSTRUCTION_PARALLEL
  nwalk = treeConstructionParallel( &tree_tmp,
				    n,
				    tree_tmp.walklist,
				    treeparam,
				    p_cache);
#else
  int ncell_all = 1;
  treeConstruction( &tree_tmp, 0, n, 1, key_level, 
		    tree_tmp.walklist, &nwalk, ipflag,
		    0, treeparam, &ncell_all, p_cache);
#endif
#else //TREE2
  int col_cell_index = tree_tmp.tree_clistmask + 1;  //for hash collision
  tree_tmp.tree_cell[1].key = 1;
  tree_tmp.tree_cell[1].n = n;
  treeConstruction( tree_tmp.tree_cell, tree_tmp.key, tree_tmp.index, particle, 0, n, 
		    tree_tmp.tree_clistmask, &col_cell_index, 1, key_level,
		    tree_tmp.walklist, &nwalk, ipflag, 1, treeparam);
#endif//TREE2
  this_run->t.pp_get_boundary_construct_tree = getTime(&nowtime);

  fprintf_verbose( stderr, "#PP getBoundaryParticle sendRecvBoundaryUsingTree\n");
  int bn = sendRecvBoundaryUsingTree( particle, treeparam, n,
				      bmin_all, bmax_all, comm_flag,
				      this_run, &tree_tmp);

#ifndef BUFFER_FOR_TREE
  TreeTmpDTreeTmp( &tree_tmp);
#endif

  return bn;

}



void createMPIParticle( MPI_Datatype *MPI_PARTICLE){

  MPI_Type_contiguous( sizeof(Particle), MPI_CHAR, MPI_PARTICLE);
  MPI_Type_commit( MPI_PARTICLE);

}

    } // namespace ParticleMesh
}     // namespace ParticleSimulator


