#include "treepm_header.h"
#include "pp.h"
#include "tools.h"

namespace ParticleSimulator{
    namespace ParticleMesh{


static char cbuf[CHARMAX];
static char cbuf2[CHARMAX];



static void setNTotal( pRunParam this_run){

  this_run->npart_total = 0;
  long long int npart = (long long int)this_run->npart;
  MPI_Allreduce( &npart, &(this_run->npart_total), 
		 1, MPI_LONG_LONG_INT, MPI_SUM, this_run->MPI_COMM_INTERNAL);

  fprintf( stderr, "npart_total= %lld\n", this_run->npart_total);

}



static void posToParticle( pParticle particle, const float *pos){

  particle->xpos = pos[0];
  particle->ypos = pos[1];
  particle->zpos = pos[2];

}



int determineBoxOverlapping( const double cbox1[3],
                             const double box_half_size1[3],
                             const double cbox2[3],
                             const double box_half_size2[3]){

  double rrel[3];
  double length[3];
  double length2[3];
  int flag = 0;
  int flag2 = 0;

  int j;
  for( j=0; j<3; j++){
    rrel[j] = cbox1[j] - cbox2[j];
    if( rrel[j] > 0.50)  rrel[j] -= 1.0;
    if( rrel[j] < -0.50) rrel[j] += 1.0;
    length[j]  = box_half_size1[j] + box_half_size2[j];
    length2[j] = box_half_size1[j] - box_half_size2[j];
    if( fabs(rrel[j]) <= length[j])  flag += 1;
    if( fabs(rrel[j]) <= length2[j])  flag2 += 1;
  }

  if( flag == 3){
    flag = 1;
  }
  else{
    flag = 0;
  }

  return flag;

}



/* convert endian */
float swapFloat( const float f){

  union{
    float fval;
    char cval[4];
  }u;

  u.fval = f;

  char tmp;
  tmp = u.cval[0];
  u.cval[0] = u.cval[3];
  u.cval[3] = tmp;
  tmp = u.cval[1];
  u.cval[1] = u.cval[2];
  u.cval[2] = tmp;

  return u.fval;

}



int swapInt( const int i){

  union{
    int ival;
    char cval[4];
  }u;

  u.ival = i;

  char tmp;
  tmp = u.cval[0];
  u.cval[0] = u.cval[3];
  u.cval[3] = tmp;
  tmp = u.cval[1];
  u.cval[1] = u.cval[2];
  u.cval[2] = tmp;

  return u.ival;

}



double swapDouble( const double d){

  union ud{
    double dval;
    char cval[8];
  }u;

  u.dval = d;

  char tmp;
  tmp = u.cval[0];
  u.cval[0] = u.cval[7];
  u.cval[7] = tmp;
  tmp = u.cval[1];
  u.cval[1] = u.cval[6];
  u.cval[6] = tmp;
  tmp = u.cval[2];
  u.cval[2] = u.cval[5];
  u.cval[5] = tmp;
  tmp = u.cval[3];
  u.cval[3] = u.cval[4];
  u.cval[4] = tmp;

  return u.dval;

}



long long int swapLLInt( const long long int ll){


  union ud{
    char cval[8];
    long long int llval;
  }u;

  u.llval = ll;

  char tmp;
  tmp = u.cval[0];
  u.cval[0] = u.cval[7];
  u.cval[7] = tmp;
  tmp = u.cval[1];
  u.cval[1] = u.cval[6];
  u.cval[6] = tmp;
  tmp = u.cval[2];
  u.cval[2] = u.cval[5];
  u.cval[5] = tmp;
  tmp = u.cval[3];
  u.cval[3] = u.cval[4];
  u.cval[4] = tmp;

  return u.llval;

}



void swapParticle( Particle *p){

  p->xpos = swapDouble(p->xpos);
  p->ypos = swapDouble(p->ypos);
  p->zpos = swapDouble(p->zpos);
  p->xvel = swapFloat(p->xvel);
  p->yvel = swapFloat(p->yvel);
  p->zvel = swapFloat(p->zvel);

#ifdef LONG_ID
  p->id = swapLLInt(p->id);
#else
  p->id = swapInt(p->id);
#endif

#ifndef UNIFORM
  p->mass = swapFloat(p->mass);
#endif


}



void inputParam( pRunParam this_run, double *dtime, char const *paramfile){


  FILE *parameters = my_fopen(paramfile,"r");

  double dum;
  //fscanf(parameters,"%f",&(this_run->maxcpu));
  fscanf(parameters,"%lf",&dum);
  fscanf(parameters,"%lf",dtime);
  fscanf(parameters,"%d",&(this_run->noutput));

  this_run->toutput = (double *) my_malloc( sizeof(double) * (this_run->noutput+1));

  for( int i=0;i<this_run->noutput;i++) {
    fscanf(parameters,"%lf",(this_run->toutput)+i);
  }
  this_run->outflag = this_run->noutput;

  this_run->toutput[this_run->noutput] = -1.0e30;
  for( int i=0;i<this_run->noutput;i++) {
    if (this_run->znow > this_run->toutput[i]) {
      this_run->outflag = i;
      break;
    }
  }

  fclose(parameters);

#ifdef CONSTANT_TIMESTEP
  *dtime = constant_timestep;
#endif


}



/* input snapshot header */
static void inputSnapshotHeader( pRunParam this_run, FILE *inputdata){

  fread( &(this_run->io_ver) , sizeof(int),    1, inputdata);
#ifdef REVERSE_ENDIAN_INPUT
  this_run->io_ver = swapInt( this_run->io_ver);
#endif
  if( this_run->io_ver == 1000){
    fread( &(this_run->npart_total)  , sizeof(long long int), 1, inputdata);
#ifndef REVERSE_ENDIAN_INPUT
    if( this_run->npart_total < 2147483648){
      this_run->npart = (int)(this_run->npart_total);
    }
#endif
    fread( &(this_run->ngas)   , sizeof(long long int),    1, inputdata);
  }
  else{
    fread( &(this_run->npart)  , sizeof(int), 1, inputdata);
    fread( &(this_run->ngas)   , sizeof(int), 1, inputdata);
  }
  fread( &(this_run->cosm.omega0) , sizeof(float),  1, inputdata);
  fread( &(this_run->cosm.omegab) , sizeof(float),  1, inputdata);
  fread( &(this_run->cosm.lambda0), sizeof(float),  1, inputdata);
  fread( &(this_run->cosm.hubble) , sizeof(float),  1, inputdata);
  fread( &(this_run->astart) , sizeof(float),  1, inputdata);
  fread( &(this_run->anow)   , sizeof(float),  1, inputdata);
  fread( &(this_run->tnow)   , sizeof(float),  1, inputdata);
  fread( &(this_run->lunit)  , sizeof(double), 1, inputdata);
  fread( &(this_run->munit)  , sizeof(double), 1, inputdata);
  fread( &(this_run->tunit)  , sizeof(double), 1, inputdata);

#ifdef REVERSE_ENDIAN_INPUT
  if( this_run->io_ver == 1000){
    this_run->npart_total = swapLLInt( this_run->npart_total);
    if( this_run->npart_total < 2147483648){
      this_run->npart = (int)(this_run->npart_total);
    }
    this_run->ngas = swapLLInt( this_run->ngas);
  }
  else{
    this_run->npart = swapInt( this_run->npart);
    this_run->ngas = swapInt( this_run->ngas);
  }
  this_run->cosm.omega0 = swapFloat( this_run->cosm.omega0);
  this_run->cosm.omegab = swapFloat( this_run->cosm.omegab);
  this_run->cosm.lambda0 = swapFloat( this_run->cosm.lambda0);
  this_run->cosm.hubble = swapFloat( this_run->cosm.hubble);
  this_run->astart = swapFloat( this_run->astart);
  this_run->anow = swapFloat( this_run->anow);
  this_run->tnow = swapFloat( this_run->tnow);
  this_run->lunit = swapDouble( this_run->lunit);
  this_run->munit = swapDouble( this_run->munit);
  this_run->tunit = swapDouble( this_run->tunit);
#endif

  this_run->znow = (float)(1.0 / this_run->anow - 1.0);
  if (this_run->ngas > 0) {
    fprintf(stderr,"Warning: Number of gas particles is not zero.\n");
    fprintf(stderr,"Warning: Gas dynamical effect is ignored.");
  }

  fprintf( stderr, "io_ver:%d\n", this_run->io_ver);
  fprintf( stderr, "omega0:%e lambda0:%e hubble:%e\n",
           this_run->cosm.omega0, this_run->cosm.lambda0,
           this_run->cosm.hubble);
  fprintf( stderr, "astart:%e anow:%e tnow:%e znow:%e\n",
           this_run->astart, this_run->anow,
           this_run->tnow, this_run->znow);
  fprintf( stderr, "lunit:%e munit:%e tunit:%e\n",
           this_run->lunit, this_run->munit, this_run->tunit);
  fprintf( stderr, "npart:%d\n",
           this_run->npart);



}



void inputOneData( pRunParam this_run, pParticle particle, 
		   const char *filename, const int *ndiv){

  int i,ii,iii;

  int idum;
  float fdum;
  double ddum;
  long long int llidum;

  int inode = this_run->inode;
  int nnode = this_run->nnode;

  int x, y, z;
  getXYZIndex( inode, ndiv, &x, &y, &z);
  double hx = 1.0/(float)ndiv[0];
  double hy = 1.0/(float)ndiv[1];
  double hz = 1.0/(float)ndiv[2];
  this_run->bmin[0] = (double)x     * hx;
  this_run->bmin[1] = (double)y     * hy;
  this_run->bmin[2] = (double)z     * hz;
  this_run->bmax[0] = (double)(x+1) * hx;
  this_run->bmax[1] = (double)(y+1) * hy;
  this_run->bmax[2] = (double)(z+1) * hz;

  int n_memory = IO_CACHE_SIZE;
  float *xsend  = (float *) my_malloc( sizeof(float) * n_memory * 6);
  int   *ixsend = (int *)   my_malloc( sizeof(int)   * n_memory);	
  float *xrecv  = (float *) my_malloc( sizeof(float) * n_memory * 6);
  int   *ixrecv = (int *)   my_malloc( sizeof(int)   * n_memory);	

  FILE *inputdata = my_fopen(filename,"r");
  inputSnapshotHeader( this_run, inputdata);
  update_now(this_run);

#ifndef UNIFORM
  double uniform_mass = 3.0*(this_run->cosm.omega0)/(8*M_PI*(float)(NUMBER_OF_PART_ALL));
#endif

  float *cache_array = (float *) my_malloc(sizeof(float)*n_memory*3);
  float *cache_array2 = (float *) my_malloc(sizeof(float)*n_memory*3);

  int ipart=0;
  for(i=0;i<(this_run->npart);i+=n_memory){

    int nn = n_memory;
    // i=%last# operation
    if((i+n_memory)>(this_run->npart)) nn = (this_run->npart)-i;
    int ic = i/n_memory;
        
    // rank = 0
    static int ns[MAXNNODE], ns2[MAXNNODE], offset[MAXNNODE], displs1[MAXNNODE], displs2[MAXNNODE], ns3[MAXNNODE];

    if(inode==0){
      if(ic!=0){
	inputdata = my_fopen(filename,"r");
	//for(ii=0;ii<3;ii++)fread(&idum,sizeof(int),1,inputdata);
	fread( &idum,sizeof(int),1,inputdata);
	if( this_run->io_ver == 1000){
	  fread( &llidum,sizeof(long long int),1,inputdata);
	  fread( &llidum,sizeof(long long int),1,inputdata);
	}
	else{
	  fread( &idum,sizeof(int),1,inputdata);
	  fread( &idum,sizeof(int),1,inputdata);
	}
	for(ii=0;ii<7;ii++)fread(&fdum,sizeof(float),1,inputdata);
	for(ii=0;ii<3;ii++)fread(&ddum,sizeof(double),1,inputdata);       
      }
      for(ii=0;ii<ic;ii++){  /* head */
	fread(cache_array,sizeof(float),3*n_memory,inputdata);
      }
      fread(cache_array,sizeof(float),3*nn,inputdata);     /* main */
      for(ii=(i+n_memory);ii<(this_run->npart);ii+=n_memory){ /* tail */
	int nnn = n_memory;
	if((ii+n_memory)>(this_run->npart)) nnn = (this_run->npart)-ii; 
	fread(cache_array2,sizeof(float),3*nnn,inputdata);      
      }
      for(ii=0;ii<ic;ii++){ /* head */
	fread(cache_array2,sizeof(float),3*n_memory,inputdata);
      }
      fread(cache_array2,sizeof(float),3*nn,inputdata); /* main */
      fclose(inputdata);

#ifdef REVERSE_ENDIAN_INPUT
      for( ii=0; ii<nn; ii++){
        cache_array[3*ii+0] = swapFloat(cache_array[3*ii+0]);
        cache_array[3*ii+1] = swapFloat(cache_array[3*ii+1]);
        cache_array[3*ii+2] = swapFloat(cache_array[3*ii+2]);
        cache_array2[3*ii+0] = swapFloat(cache_array2[3*ii+0]);
        cache_array2[3*ii+1] = swapFloat(cache_array2[3*ii+1]);
        cache_array2[3*ii+2] = swapFloat(cache_array2[3*ii+2]);
      }
#endif
      
      for(ii=0;ii<nn;ii++){
	if( cache_array[3*ii+0] >= 1.0)  cache_array[3*ii+0] -= 1.0;
	if( cache_array[3*ii+0] < 0.0)  cache_array[3*ii+0] += 1.0;
	if( cache_array[3*ii+1] >= 1.0)  cache_array[3*ii+1] -= 1.0;
	if( cache_array[3*ii+1] < 0.0)  cache_array[3*ii+1] += 1.0;
	if( cache_array[3*ii+2] >= 1.0)  cache_array[3*ii+2] -= 1.0;
	if( cache_array[3*ii+2] < 0.0)  cache_array[3*ii+2] += 1.0;
	int ix = (int)( cache_array[3*ii+0]/hx);
	int iy = (int)( cache_array[3*ii+1]/hy);
	int iz = (int)( cache_array[3*ii+2]/hz);
	int p_inode = getVoxelIndex( ix, iy, iz, ndiv);
        if(p_inode==0){
	  posToParticle( &particle[ipart], &cache_array[3*ii]);
          particle[ipart].xvel = cache_array2[3*ii+0];
          particle[ipart].yvel = cache_array2[3*ii+1];
          particle[ipart].zvel = cache_array2[3*ii+2];
#ifndef UNIFORM
          particle[ipart].mass = uniform_mass;
#endif
          particle[ipart].id = i+ii;
          ipart++;
	}
      }

      for( iii=0; iii<nnode; iii++){
	ns[iii] = ns2[iii] = offset[iii] = displs1[iii] = displs2[iii] = 0;
	ns3[iii] = 1;
      }
      for( ii=0; ii<nn; ii++){
	int ix = (int)( cache_array[3*ii+0]/hx);
	int iy = (int)( cache_array[3*ii+1]/hy);
	int iz = (int)( cache_array[3*ii+2]/hz);
	int p_inode = getVoxelIndex( ix, iy, iz, ndiv);
	ns[p_inode] ++;
      }
      for( iii=1; iii<nnode; iii++){
	offset[iii] = ns[iii-1] + offset[iii-1];
	displs1[iii] = displs1[iii-1] + 1;
	displs2[iii] = displs2[iii-1] + ns[iii-1];
      }
      for( ii=0; ii<nn; ii++){
	int ix = (int)( cache_array[3*ii+0]/hx);
	int iy = (int)( cache_array[3*ii+1]/hy);
	int iz = (int)( cache_array[3*ii+2]/hz);
	int p_inode = getVoxelIndex( ix, iy, iz, ndiv);
	if( p_inode == 0)  continue;
	int index = ns2[p_inode] + offset[p_inode];
	xsend[6*index+0] = cache_array[3*ii+0];
	xsend[6*index+1] = cache_array[3*ii+1];
	xsend[6*index+2] = cache_array[3*ii+2];
	xsend[6*index+3] = cache_array2[3*ii+0];
	xsend[6*index+4] = cache_array2[3*ii+1];
	xsend[6*index+5] = cache_array2[3*ii+2];
	ixsend[index] = i + ii;
	ns2[p_inode] ++;
      }
    }

    int nrecv = 0;
    MPI_Scatterv( ns2, ns3, displs1, MPI_INT,
		  &nrecv, 1, MPI_INT, 0, this_run->MPI_COMM_INTERNAL);
    MPI_Scatterv( ixsend, ns2, displs2, MPI_INT,
		  ixrecv, nrecv, MPI_INT, 0, this_run->MPI_COMM_INTERNAL);
    for( iii=0; iii<nnode; iii++){
      ns2[iii] *= 6;
      displs2[iii] *= 6;
    }
    MPI_Scatterv( xsend, ns2, displs2, MPI_FLOAT,
		  xrecv, nrecv*6, MPI_FLOAT, 0, this_run->MPI_COMM_INTERNAL);
    if( inode != 0){
      for( ii=0; ii<nrecv; ii++){
	posToParticle( &particle[ipart], &xrecv[6*ii]);
	if( particle[ipart].xpos < 0.0)  particle[ipart].xpos += 1.0;
	if( particle[ipart].ypos < 0.0)  particle[ipart].ypos += 1.0;
	if( particle[ipart].zpos < 0.0)  particle[ipart].zpos += 1.0;
	if( particle[ipart].xpos >= 1.0)  particle[ipart].xpos -= 1.0;
	if( particle[ipart].ypos >= 1.0)  particle[ipart].ypos -= 1.0;
	if( particle[ipart].zpos >= 1.0)  particle[ipart].zpos -= 1.0;
        particle[ipart].xvel = xrecv[6*ii+3];
        particle[ipart].yvel = xrecv[6*ii+4];
        particle[ipart].zvel = xrecv[6*ii+5];
#ifndef UNIFORM
	particle[ipart].mass = uniform_mass;
#endif
        particle[ipart].id = ixrecv[ii];
        ipart++;
      }
    }
  }

  //fprintf(stderr,"node %d: N %d Ntotal %d\n",inode,ipart,(this_run->npart));

  (*this_run).npart_total = (this_run->npart);
  (*this_run).npart = ipart;  

  if (this_run->npart > NUMBER_OF_PART) {
    fprintf(stderr,"Error, Number of particles exceeds the limit.\n");
    exit(EXIT_FAILURE);
  }

  this_run->uniform_mass = 3.0*(this_run->cosm.omega0)/(8*M_PI*(float)(this_run->npart_total));
  //fprintf( stderr, "%.15lf\n", this_run->uniform_mass);

  this_run->min_mass = this_run->uniform_mass;

  free(cache_array); 
  free(cache_array2);
  free(xrecv); 
  free(ixrecv);       
  free(xsend);
  free(ixsend);


}




/* input plural file*/
void inputPluralDataContinue( pRunParam this_run, pParticle particle, char *basefile){

  int inode = this_run->inode;

#ifndef Parallel_IO
  int nnode = this_run->nnode;
  for( int j=0; j<nnode; j++){
    MPI_Barrier(this_run->MPI_COMM_INTERNAL);
    if( j == inode){
#endif
      sprintf( cbuf2, "%s-%d", basefile, inode);
      fprintf( stderr, "node:%d input %s\n", inode, cbuf2);
      fileCheck(cbuf2);
#ifdef GADGET_IO
      readGadget( *this_run, particle, cbuf2);
#else
      inputPluralData( this_run, particle, cbuf2);
#endif
#ifndef Parallel_IO
    }
  }
#endif

  setNTotal( this_run);

}



/* input initial conditions divided into all node in advance*/
void inputPluralData( pRunParam this_run, pParticle particle, const char *filename){


  int i, j;
  this_run->npart_total = NUMBER_OF_PART_ALL;

  FILE *inputdata = my_fopen(filename,"r");
  inputSnapshotHeader( this_run, inputdata);
  update_now(this_run);

  int nmemory = IO_CACHE_SIZE;
  float *cache_array = (float *) my_malloc(sizeof(float)*nmemory*3);
  for( i=0; i<(this_run->npart); i+=nmemory){
    int nread = nmemory;
    if( (i+nmemory) > this_run->npart)  nread = this_run->npart - i;
    fread( cache_array, sizeof(float), nread*3, inputdata);
    for( j=0; j<nread; j++){
      int index = i+j;
#ifdef REVERSE_ENDIAN_INPUT
      particle[index].xpos = swapFloat(cache_array[3*j+0]);
      particle[index].ypos = swapFloat(cache_array[3*j+1]);
      particle[index].zpos = swapFloat(cache_array[3*j+2]);
#else
      particle[index].xpos = cache_array[3*j+0];
      particle[index].ypos = cache_array[3*j+1];
      particle[index].zpos = cache_array[3*j+2];
#endif
    }
  }
  for( i=0; i<(this_run->npart); i+=nmemory){
    int nread = nmemory;
    if( (i+nmemory) > this_run->npart)  nread = this_run->npart - i;
    fread( cache_array, sizeof(float), nread*3, inputdata);
    fflush( inputdata);
    for( j=0; j<nread; j++){
      int index = i+j;
#ifdef REVERSE_ENDIAN_INPUT
      particle[index].xvel = swapFloat(cache_array[3*j+0]);
      particle[index].yvel = swapFloat(cache_array[3*j+1]);
      particle[index].zvel = swapFloat(cache_array[3*j+2]);
#else
      particle[index].xvel = cache_array[3*j+0];
      particle[index].yvel = cache_array[3*j+1];
      particle[index].zvel = cache_array[3*j+2];
#endif
    }
  }
  free( cache_array);

  IDTYPE *cache_array2 = (IDTYPE *) my_malloc(sizeof(IDTYPE)*nmemory);
  for( i=0; i<(this_run->npart); i+=nmemory){
    int nread = nmemory;
    if( (i+nmemory) > this_run->npart)  nread = this_run->npart - i;
    fread( cache_array2, sizeof(IDTYPE), nread, inputdata);
    fflush( inputdata);
    for( j=0; j<nread; j++){
      int index = i+j;
#ifdef REVERSE_ENDIAN_INPUT
#ifdef LONG_ID
      particle[index].id = swapLLInt(cache_array2[j]);
#else
      particle[index].id = swapInt(cache_array2[j]);
#endif
#else
      particle[index].id = cache_array2[j];
#endif
    }
  }
  free( cache_array2);

  if (this_run->npart > NUMBER_OF_PART) {
    fprintf(stderr,"Error, Number of particles exceeds the limit.\n");
    exit(EXIT_FAILURE);
  }

  this_run->uniform_mass = 3.0*(this_run->cosm.omega0)/(8*M_PI*(double)(this_run->npart_total));
  this_run->min_mass = this_run->uniform_mass;


}



static void outputSnapshotEconomyHeader( const RunParam *this_run,
					 const char *filename){

  if( this_run->inode == 0){
    FILE *outstream = my_fopen( filename, "w");
    fprintf( outstream, "step    %d\n", this_run->nstep);
    fprintf( outstream, "znow    %e\n", this_run->znow);
    fprintf( outstream, "anow    %e\n", this_run->anow);
    fprintf( outstream, "tnow    %e\n", this_run->tnow);
    fprintf( outstream, "omega0  %e\n", this_run->cosm.omega0);
    fprintf( outstream, "omegab  %e\n", this_run->cosm.omegab);
    fprintf( outstream, "lambda0 %e\n", this_run->cosm.lambda0);
    fprintf( outstream, "hubble  %e\n", this_run->cosm.hubble);
    fprintf( outstream, "astart  %e\n", this_run->astart);
    fprintf( outstream, "lunit   %e\n", this_run->lunit);
    fprintf( outstream, "munit   %e\n", this_run->munit);
    fprintf( outstream, "tunit   %e\n", this_run->tunit);
    fprintf( outstream, "mass    %e\n", this_run->uniform_mass);
    fclose(outstream);
  }

}






void inputInitWrapper( pRunParam this_run, pParticle particle, 
		       char *init_file, const int nfiles, 
		       const char *mass_file){


  double nowtime;

  getTime( &nowtime);
  if( nfiles == 1){
    inputOneData( this_run, particle, init_file, this_run->ndiv);
  }
  else{
    assert ( nfiles == this_run->nnode);
    inputPluralDataContinue(this_run, particle, init_file);
  }

  MPI_Barrier( MPI_COMM_WORLD);
  getTimeFprint( &nowtime, "Load Initial Condition", this_run->profile_file);


}



static void outputBoundary( const pRunParam this_run, const char *filename,
			    const double bmin_all[][3], const double bmax_all[][3]){

  int nnode = this_run->nnode;
  int idump[4] = { this_run->nnode,
		   this_run->ndiv[0],
		   this_run->ndiv[1],
		   this_run->ndiv[2]};

  if( this_run->inode == 0){    
    fprintf( stderr, "writing %s nnode:%d ndiv:%d %d %d\n",
	     filename, idump[0], idump[1], idump[2], idump[3]);
    FILE *fout = fopen( filename, "wb");
    fwrite( idump, sizeof(int), 4, fout);
    fwrite( bmin_all, sizeof(double), 3*nnode, fout);
    fwrite( bmax_all, sizeof(double), 3*nnode, fout);
    fclose( fout);

  }

}



static void inputBoundary( const pRunParam this_run, const char *filename,
			   const int nfiles,
			   double bmin_all[][3], double bmax_all[][3]){

  int idump[4];
    
  if( this_run->inode == 0){
    FILE *fin = fopen( filename, "rb");
    fread( idump, sizeof(int), 4, fin);
    fprintf( stderr, "reading %s nnode:%d ndiv:%d %d %d\n",
	     filename, idump[0], idump[1], idump[2], idump[3]);
    fread( bmin_all, sizeof(double), 3*nfiles, fin);
    fread( bmax_all, sizeof(double), 3*nfiles, fin);
    fclose( fin);
  }

  MPI_Bcast( bmin_all, 3*nfiles, MPI_DOUBLE, 0, this_run->MPI_COMM_INTERNAL);
  MPI_Bcast( bmax_all, 3*nfiles, MPI_DOUBLE, 0, this_run->MPI_COMM_INTERNAL);


}




typedef struct MPI_Comm3D{

  MPI_Comm comm_x, comm_y, comm_z;

  MPI_Comm3D( int nx, int ny, int nz,
	      int xrank, int yrank, int zrank,
	      MPI_Comm comm_world){

    {
      int color = zrank + nz * yrank;
      int key   = xrank;
      MPI_Comm_split(comm_world, color, key, &comm_x);
      int size, rank;
      MPI_Comm_size(comm_x, &size);
      MPI_Comm_rank(comm_x, &rank);
      assert(size == nx);
      assert(rank == xrank);
    }

    {
      int color = xrank + nx * zrank;
      int key   = yrank;
      MPI_Comm_split(comm_world, color, key, &comm_y);
      int size, rank;
      MPI_Comm_size(comm_y, &size);
      MPI_Comm_rank(comm_y, &rank);
      assert(size == ny);
      assert(rank == yrank);
    }

    {
      int color = yrank + ny * xrank;
      int key   = zrank;
      MPI_Comm_split(comm_world, color, key, &comm_z);
      int size, rank;
      MPI_Comm_size(comm_z, &size);
      MPI_Comm_rank(comm_z, &rank);
      assert(size == nz);
      assert(rank == zrank);
    }
  }

}MPI_Comm3D;



static inline int whichBox3D( int ibox[3],
			      const double *pos,
			      const double bmin[][3],
			      const double bmax[][3],
			      const int npx,
			      const int npy,
			      const int npz){

  int p = 0;
  int ix, iy, iz;
  if( pos[0] < bmin[p][0]) return -1;
  for(ix=0; ix<npx; ix++, p+=npy*npz){
    if(pos[0] < bmax[p][0]) break;
  }
  if(pos[0] > bmax[p][0]) return -1;

  if(pos[1] < bmin[p][1]) return -1;
  for(iy=0; iy<npy; iy++, p+=npz){
    if(pos[1] < bmax[p][1]) break;
  }
  if(pos[1] > bmax[p][1]) return -1;

  if(pos[2] < bmin[p][2]) return -1;
  for(iz=0; iz<npz; iz++, p++){
    if(pos[2] < bmax[p][2]) break;
  }
  if(pos[2] > bmax[p][2]) return -1;

  ibox[0] = ix;
  ibox[1] = iy;
  ibox[2] = iz;

  return p;

}



static int commParticle( const int np, int *nsend, 
			 const pParticle psend, pParticle precv,
			 MPI_Comm comm, MPI_Datatype MPI_PARTICLE){

  static int displs[MAXNNODE];
  static int rdispls[MAXNNODE];
  static int nrecv[MAXNNODE];

  displs[0] = 0;
  for( int i=1; i<np; i++)  displs[i] = displs[i-1] + nsend[i-1];
  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, comm);

  int nrecv_total = nrecv[0];
  rdispls[0] = 0;
  for( int i=1; i<np; i++){
    rdispls[i] = rdispls[i-1] + nrecv[i-1];
    nrecv_total += nrecv[i];
  }

  MPI_Alltoallv( psend, nsend, displs, MPI_PARTICLE,
		 precv, nrecv, rdispls, MPI_PARTICLE, comm);

  return nrecv_total;

}



void readSamplingParticleSetBoundary( pRunParam this_run, 
				      const char *filename){

  int inode = this_run->inode;

  int nsamp = 0;
  FILE *fin = NULL;
  if( inode == 0){
    fin = fopen( filename, "rb");
    fread( &nsamp, sizeof(int), 1, fin);
    fprintf( stderr, "reading %s nsamp:%d\n", filename, nsamp);
  }
  MPI_Bcast( &nsamp, 1, MPI_INT, 0, MPI_COMM_WORLD);

  float *xsamp_all = new float[nsamp];
  float *ysamp_all = new float[nsamp];
  float *zsamp_all = new float[nsamp];

  if( inode == 0){
    fread( xsamp_all, sizeof(float), nsamp, fin);
    fread( ysamp_all, sizeof(float), nsamp, fin);
    fread( zsamp_all, sizeof(float), nsamp, fin);
    fclose( fin);
  }

  MPI_Bcast( xsamp_all, nsamp, MPI_FLOAT, 0, this_run->MPI_COMM_INTERNAL);
  MPI_Bcast( ysamp_all, nsamp, MPI_FLOAT, 0, this_run->MPI_COMM_INTERNAL);
  MPI_Bcast( zsamp_all, nsamp, MPI_FLOAT, 0, this_run->MPI_COMM_INTERNAL);

  determineBoundary( xsamp_all, ysamp_all, zsamp_all, nsamp,
		     this_run, this_run->bmin, this_run->bmax);
  MPI_Allgather( this_run->bmin, 3, MPI_DOUBLE, 
		 this_run->bmin_all, 3, MPI_DOUBLE, this_run->MPI_COMM_INTERNAL);
  MPI_Allgather( this_run->bmax, 3, MPI_DOUBLE, 
		 this_run->bmax_all, 3, MPI_DOUBLE, this_run->MPI_COMM_INTERNAL);

  delete [] xsamp_all;
  delete [] ysamp_all;
  delete [] zsamp_all;



}



void inputInitWrapper3( pRunParam this_run, pParticle particle, 
			const int nfiles, const char *filename){

  int inode = this_run->inode;
  int nnode = this_run->nnode;
  int nx = this_run->ndiv[0];
  int ny = this_run->ndiv[1];
  int nz = this_run->ndiv[2];
  MPI_Datatype MPI_PARTICLE;
  createMPIParticle( &MPI_PARTICLE);

  double nowtime;
  double time1;
  getTime( &time1);

  sprintf( cbuf, "%s.info", filename);
  if( access( cbuf, F_OK) == 0){
    FILE *instream = fopen( cbuf, "r");
    if( instream != NULL)  fscanf( instream, "%d", &this_run->nstep);
    fprintf( stderr, "nstep=%d\n", this_run->nstep);
    fclose(instream);
  }

  sprintf( cbuf, "%s.samp", filename);
  if( access( cbuf, F_OK) == 0){
    readSamplingParticleSetBoundary( this_run, cbuf);
  }
  else{
    setUniformBoundary( this_run, this_run->bmin_all, this_run->bmax_all);
  }

  int ninput = nfiles / nnode;
  int ninput2 = nfiles % nnode;
  if( inode < ninput2)  ninput ++;
  int istart = ninput*inode;
  if( inode >= ninput2)  istart += ninput2;
  int iend = istart + ninput;

  int nparttmp = 0;
  for( int i=istart; i<iend; i++){
    sprintf( cbuf, "%s-%d", filename, i);
    fileCheck( cbuf);
#ifdef GADGET_IO
    readGadget( *this_run, particle, cbuf);
#else
    inputPluralData( this_run, particle+nparttmp, cbuf);
#endif
    fprintf( stderr, "node: %d\t read: %s\t npart: %d\t npart_cum: %d\n", 
	     inode, cbuf, this_run->npart, nparttmp);
    nparttmp += this_run->npart;
  }
  this_run->npart = nparttmp;
  correctBoundaryCondition( particle, this_run->npart);

  if (this_run->npart > NUMBER_OF_PART) {
    fprintf(stderr,"Error: Number of particles in inode=%d exceeds the limit (%d > %d).\n",
	    inode, this_run->npart, NUMBER_OF_PART);
    exit(1);
  }

  MPI_Bcast( &(this_run->cosm.omega0), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->cosm.omegab), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->cosm.lambda0), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->cosm.hubble), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->astart), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->anow), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->znow), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->tnow), 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->lunit), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->munit), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast( &(this_run->tunit), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int npart = this_run->npart;

  int *nsend = new int[nnode];
  pParticle sendbuf = new Particle[NUMBER_OF_PART];
  int *dest = new int[NUMBER_OF_PART];

  if( inode == 0)  cerr << "start comm" << endl;

  int *ndiv = this_run->ndiv;
  int x, y, z;
  getXYZIndex( inode, this_run->ndiv, &x, &y, &z);
  MPI_Comm3D comm3d( ndiv[0], ndiv[1], ndiv[2], x, y, z, MPI_COMM_WORLD);
  int np[3] = { nx, ny, nz};
  MPI_Comm comm_xyz[3] = {comm3d.comm_x, comm3d.comm_y, comm3d.comm_z};

  for( int k=0; k<3; k++){

    for( int i=0; i<npart; i++){
      double pos[3];
      getPos( &particle[i], pos);
      int ibox[3];
      whichBox3D( ibox, pos, this_run->bmin_all, this_run->bmax_all, nx, ny, nz);
      dest[i] = ibox[k];
      assert( dest[i] >= 0);
      assert( dest[i] < np[k]);
      sendbuf[i] = particle[i];
    }

    int nite_comm = NUMBER_OF_COMM_PER_ALLTOALLV_FOR_INPUT_IC;
    int npart00 = npart / nite_comm;
    int offset = 0;
    int nrecv_total = 0;
    for( int i=0; i<nite_comm; i++){
      getTime( &nowtime);
      int npart_this = npart00;
      if( i==0)  npart_this += npart % nite_comm;
      qsortPivotWrapper( sendbuf+offset, dest+offset, npart_this, 
			 NLEVEL_PARTICLE_SORT, NCHUNK_QSORT);
      for( int j=0; j<np[k]; j++)  nsend[j] = 0;
      for( int j=0; j<npart_this; j++)  nsend[dest[j+offset]] ++;
      int nrecv_this = commParticle( np[k], nsend, sendbuf+offset, 
				     particle+nrecv_total, comm_xyz[k], MPI_PARTICLE);
      nrecv_total += nrecv_this;
      offset += npart_this;
      My_MPI_Barrier( MPI_COMM_WORLD);
      if( inode == 0){
	cerr << "alltoallv(load IC) " << k+1 << "/3" << "\t" << i+1 << "/" << nite_comm 
	     << "\t" << getTime(&nowtime) << "sec" << endl;
      }
      My_MPI_Barrier( MPI_COMM_WORLD);
    }
    npart = nrecv_total;
  }
  this_run->npart = npart;

  delete [] nsend;
  delete [] sendbuf;
  delete [] dest;

  setNTotal( this_run);

  this_run->uniform_mass = 3.0*(this_run->cosm.omega0)/(8*M_PI*(double)(this_run->npart_total));
  this_run->min_mass = this_run->uniform_mass;

  MPI_Barrier( MPI_COMM_WORLD);
  getTimeFprint( &time1, "Load Initial Condition", this_run->profile_file);

}




static void mergeParticle( RunParam &this_run, pParticle particle, 
			   const int fmerge){

  MPI_Status stat;
  MPI_Datatype MPI_PARTICLE;
  createMPIParticle( &MPI_PARTICLE);
  int inode = this_run.inode;
  int n = this_run.npart;
  MPI_Comm comm = this_run.MPI_COMM_INTERNAL;

  int ix, iy, iz;
  getXYZIndex( inode, this_run.ndiv, &ix, &iy, &iz);
  int receiver = iz % fmerge;
  receiver = getVoxelIndex( ix, iy, iz-receiver, this_run.ndiv);

  if( receiver == inode){
    for( int i=1; i<fmerge; i++){
      int src = inode + i;
      int nrecv = 0;
      MPI_Recv( &nrecv, 1, MPI_INT, src, 0, comm, &stat);
      cerr << "node:" << inode << " <- " << src << "\t" << nrecv << endl;
      MPI_Recv( particle+n, nrecv, MPI_PARTICLE, src, 0, comm, &stat);
      n += nrecv;
    }
    this_run.npart = n;
  }
  else{
    cerr << "node:" << inode << " -> " << receiver << "\t" << n << endl;
    MPI_Send( &n, 1, MPI_INT, receiver, 0, comm);
    MPI_Send( particle, n, MPI_PARTICLE, receiver, 0, comm);
    this_run.npart = 0;
  }


}



static void outputBoundary2( const pRunParam this_run, const char *filename,
			     const double bmin_all[][3], const double bmax_all[][3],
			     const int fmerge){
  
  int nnode = this_run->nnode;
  int idump[4] = { this_run->nnode,
		   this_run->ndiv[0],
		   this_run->ndiv[1],
		   this_run->ndiv[2]/fmerge};
  double bmin[3], bmax[3];

  if( this_run->inode == 0){    
    fprintf( stderr, "writing %s nnode:%d ndiv:%d %d %d\n",
	     filename, idump[0], idump[1], idump[2], idump[3]);
    FILE *fout = fopen( filename, "wb");
    fwrite( idump, sizeof(int), 4, fout);
    for( int i=0; i<nnode; i++){
      if( i % fmerge != 0)  continue;
      bmin[0] = bmin_all[i][0];
      bmin[1] = bmin_all[i][1];
      bmin[2] = bmin_all[i][2];
      bmax[0] = bmax_all[i][0];
      bmax[1] = bmax_all[i][1];
      bmax[2] = bmax_all[i+fmerge-1][2];
      fwrite( bmin, sizeof(double), 3, fout);
      fwrite( bmax, sizeof(double), 3, fout);
      cout << bmin[0] << "\t" 
	   << bmax[0] << "\t" 
	   << bmin[1] << "\t" 
	   << bmax[1] << "\t" 
	   << bmin[2] << "\t" 
	   << bmax[2] << endl;
    }
    fclose( fout);
  }


}




/* 
   Each node outputs snapshot individually.
 */
void outputSnapshotDump( const RunParam *this_run, 
			 const Particle *particle, const char *out_file){


  int nmemory = IO_CACHE_SIZE;
  int n = this_run->npart;
  double anow2inv = 1.0 / this_run->anow / this_run->anow;

  int io_ver = 0;

  createWritingFile( out_file);
  FILE *outputdata = my_fopen(out_file, "wb");

  //fwrite( &(this_run->io_ver),       sizeof(int),    1, outputdata);
  fwrite( &io_ver,       sizeof(int),    1, outputdata);
  fwrite( &(this_run->npart),        sizeof(int),    1, outputdata);
  fwrite( &(this_run->ngas),         sizeof(int),    1, outputdata);
  fwrite( &(this_run->cosm.omega0),  sizeof(float),  1, outputdata);
  fwrite( &(this_run->cosm.omegab),  sizeof(float),  1, outputdata);
  fwrite( &(this_run->cosm.lambda0), sizeof(float),  1, outputdata);
  fwrite( &(this_run->cosm.hubble),  sizeof(float),  1, outputdata);
  fwrite( &(this_run->astart),       sizeof(float),  1, outputdata);
  fwrite( &(this_run->anow),         sizeof(float),  1, outputdata);
  fwrite( &(this_run->tnow),         sizeof(float),  1, outputdata);
  fwrite( &(this_run->lunit),        sizeof(double), 1, outputdata);
  fwrite( &(this_run->munit),        sizeof(double), 1, outputdata);
  fwrite( &(this_run->tunit),        sizeof(double), 1, outputdata);
  fflush(outputdata);

  float (*cache)[3] = (float (*)[3]) my_malloc( sizeof(float) * nmemory * 3);
  for( int ii=0; ii<n; ii+=nmemory){
    int nwrite = nmemory;
    if( (ii+nmemory) > n)  nwrite = n - ii;
    for( int j=0; j<nwrite; j++){
      int index = ii + j;
      getPos2( &particle[index], cache[j]);
    }
    fwrite( cache, sizeof(float), nwrite*3, outputdata);
  }
  for( int ii=0; ii<n; ii+=nmemory){
    int nwrite = nmemory;
    if( (ii+nmemory) > n)  nwrite = n - ii;
    for( int j=0; j<nwrite; j++){
      int index = ii + j;
      cache[j][0] = particle[index].xvel * anow2inv;
      cache[j][1] = particle[index].yvel * anow2inv;
      cache[j][2] = particle[index].zvel * anow2inv;
    }
    fwrite( cache, sizeof(float), nwrite*3, outputdata);
    fflush(outputdata);
  }
  free( cache);
  
  IDTYPE *i_cache = (IDTYPE *) my_malloc( sizeof(IDTYPE) * nmemory);
  for( int ii=0; ii<n; ii+=nmemory){
    int nwrite = nmemory;
    if( (ii+nmemory) > n)  nwrite = n - ii;
    for( int j=0; j<nwrite; j++){
      int index = ii + j;
      i_cache[j] = particle[index].id;
    }
    fwrite( i_cache, sizeof(IDTYPE), nwrite, outputdata);
    fflush(outputdata);
  }

  free( i_cache);
  fclose( outputdata);
  removeWritingFile(out_file);


}



static int merge_snapshot_on = 0;
void outputSnapshotDumpWrapper0( pParticle particle, pRunParam this_run,
				 const char *filename){

#ifndef MULTI_TIMESTEP
  double dt = -0.5 * this_run->dtprev;
  double driftfac = getDriftFac( *this_run, this_run->anow, dt);
  drift( particle, this_run->npart, driftfac);
  this_run->anow += dt;
  this_run->znow = 1.0/this_run->anow - 1.0;
  this_run->tnow = zToTime( this_run->znow, this_run->cosm);
  exchangeParticle( particle, this_run, this_run->bmin_all, this_run->bmax_all);
#endif

  outputSnapshotEconomyHeader( this_run, filename);

  if( merge_snapshot_on == 1){
    int oldn = this_run->npart;
    mergeParticle( *this_run, particle, fmerge);
    sprintf( cbuf, "%s-%d", filename, this_run->inode/fmerge);
#ifndef Parallel_IO
    for( int j=0; j<this_run->nnode; j++){
      MPI_Barrier(this_run->MPI_COMM_INTERNAL);
      if( j == this_run->inode){
#endif
	if( this_run->npart != 0){
	  outputSnapshotDump( this_run, particle, cbuf);
	}
#ifndef Parallel_IO
      }
    }
#endif
    this_run->npart = oldn;
    sprintf( cbuf, "%s.boundary", filename);
    outputBoundary2( this_run, cbuf, this_run->bmin_all, this_run->bmax_all, fmerge);
  }
  else{
    sprintf( cbuf, "%s-%d", filename, this_run->inode);

#ifndef Parallel_IO
    for( int j=0; j<this_run->nnode; j++){
      MPI_Barrier(this_run->MPI_COMM_INTERNAL);
      if( j == this_run->inode){
#endif
#ifdef GADGET_IO
	outputGadget( *this_run, particle, cbuf);
#else
	outputSnapshotDump( this_run, particle, cbuf);
#endif
#ifndef Parallel_IO
      }
    }
#endif
    sprintf( cbuf, "%s.boundary", filename);
    outputBoundary( this_run, cbuf, this_run->bmin_all, this_run->bmax_all);
  }

  sprintf( cbuf, "%s.samp", filename);
  dumpSamplingParticle( particle, this_run, cbuf);

#ifndef MULTI_TIMESTEP
  drift( particle, this_run->npart, -driftfac);
  this_run->anow -= dt;
  this_run->znow = 1.0/this_run->anow - 1.0;
  this_run->tnow = zToTime( this_run->znow, this_run->cosm);
  exchangeParticle( particle, this_run, this_run->bmin_all, this_run->bmax_all);
#endif

}




void outputSnapshotWrapper( pParticle particle, pRunParam this_run, const double dtime){

  double nowtime = 0.0;
  getTime( &nowtime);

#ifdef MERGE_SNAPSHOT
  merge_snapshot_on = 1;
#endif

  sprintf( cbuf2,"%s/%s.%03d", 
	   SNAPSHOT, SNAPSHOT, this_run->outflag+1);
  outputSnapshotDumpWrapper0(  particle, this_run, cbuf2);

  this_run->outflag++;

  MPI_Barrier( MPI_COMM_WORLD);
  this_run->t.io_snapshot = getTime(&nowtime);

#ifdef MERGE_SNAPSHOT
  merge_snapshot_on = 0;
#endif

}



void outputImageWrapper( const Particle *particle, pRunParam this_run){

  fprintf_verbose( stderr, "#output Image start\n");

  int istep = IMAGESTEP;
  if( istep == 0)  return;
  
  double nowtime = 0.0;
  getTime( &nowtime);
  int nstep = this_run->nstep;
    
  if( nstep % istep == 0){
    int image_number = nstep / istep;
    sprintf( cbuf, "%s/%s_%05d.bmp", IMAGEDIR, MODEL, image_number);
    outputImage( particle, this_run, cbuf);
    MPI_Barrier( MPI_COMM_WORLD);
    this_run->t.io_image = getTime(&nowtime);
  }


}



static void dumpAcc0( pParticle particle, pRunParam this_run, FILE *fout){

#ifndef NOACC
  int i, ii;
  for( ii=0; ii<this_run->nnode; ii++){
    fflush( fout);
    MPI_Barrier(MPI_COMM_WORLD);
    if( ii != this_run->inode)  continue;
    fprintf( stderr, "%d\n", this_run->npart);
    for( i=0; i<this_run->npart; i++){
      fprintf( fout, "%e\t%e\t%e\t%e\t%lld\n",
	       particle[i].xacc, particle[i].yacc, particle[i].zacc,
	       particle[i].xvel, particle[i].id);
      fflush( fout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

}



void dumpAcc( pParticle particle, pRunParam this_run){

  sprintf( cbuf, "acc-%d", this_run->inode);
  FILE *fout = fopen( cbuf, "w");
  dumpAcc0( particle, this_run, fout);
  fclose( fout);

}



void readGadget( RunParam &this_run, pParticle p, const char *infile){

  int nmemory = IO_CACHE_SIZE;
  GadgetHeader gadget_header;

  FILE *fin = fopen( infile, "rb");
  int blksize = 0;
  fread( &blksize, 4, 1, fin);
  fprintf( stderr, "Reading Gadget file blksize= %d\n", blksize);
  fread( &gadget_header, sizeof(gadget_header), 1, fin);
  fread( &blksize, 4, 1, fin);
  fprintf( stderr, "End Reading Gadget file blksize= %d\n", blksize);

  this_run.znow = gadget_header.Redshift;
  this_run.anow = 1.0 / (1.0 + this_run.znow);
  this_run.cosm.omega0 = gadget_header.Omega0;
  this_run.cosm.lambda0 = gadget_header.OmegaLambda;
  this_run.cosm.hubble = gadget_header.HubbleParam;
  this_run.tnow = zToTime( this_run.znow, this_run.cosm);
  int npart = gadget_header.Npart[1];
  this_run.npart = npart;

  this_run.lunit = gadget_header.BoxSize * Gadget_UnitLength_in_Mpc / gadget_header.HubbleParam;
  double tunit = 100.0 * gadget_header.HubbleParam / 3.085678 * 1.0e-19;  ///H0=1.0 in this unit.
  tunit = 1.0 / tunit;
  double lunit = this_run.lunit;
  this_run.munit = lunit * lunit * lunit * 1.0e9 / 4.51703791e-27 / tunit / tunit * 1.0e12;
  double vunit = tunit / (lunit * 1.0e6 * 3.085678e18) * Gadget_UnitVelocity_in_cm_per_s;
  vunit /= sqrt(this_run.anow);
  lunit = 1.0 / gadget_header.BoxSize;

  this_run.tunit = tunit;
  fprintf( stderr, "Boxsize= %lf\n", gadget_header.BoxSize);
  fprintf( stderr, "znow= %lf\n", this_run.znow);
  fprintf( stderr, "anow= %lf\n", this_run.anow);
  fprintf( stderr, "omega0= %f\n", this_run.cosm.omega0);
  fprintf( stderr, "lambda0= %f\n", this_run.cosm.lambda0);
  fprintf( stderr, "hubble0= %f\n", this_run.cosm.hubble);
  fprintf( stderr, "lunit= %e\n", this_run.lunit);
  fprintf( stderr, "munit= %e\n", this_run.munit);
  fprintf( stderr, "tunit= %e\n", this_run.tunit);
  fprintf( stderr, "npart= %d\n", this_run.npart);

  fprintf( stderr, "convert_lunit= %e\n", lunit);
  fprintf( stderr, "convert_vunit= %e\n", vunit);

  float *cache = new float[nmemory*3];

  fread( &blksize, 4, 1, fin);
  fprintf( stderr, "Reading Gadget file blksize= %d\n", blksize);
  for( int i=0; i<npart; i+=nmemory){
    int nread = nmemory;
    if( (i+nmemory) > npart)  nread = npart - i;
    fread( cache, sizeof(float), 3*nread, fin);
    for( int j=0; j<nread; j++){
      p[i+j].xpos = cache[3*j+0] * lunit;
      p[i+j].ypos = cache[3*j+1] * lunit;
      p[i+j].zpos = cache[3*j+2] * lunit;
    }
  }
  fread( &blksize, 4, 1, fin);
  fprintf( stderr, "End Reading Gadget file blksize= %d\n", blksize);

  fread( &blksize, 4, 1, fin);
  fprintf( stderr, "Reading Gadget file blksize= %d\n", blksize);
  for( int i=0; i<npart; i+=nmemory){
    int nread = nmemory;
    if( (i+nmemory) > npart)  nread = npart - i;
    fread( cache, sizeof(float), 3*nread, fin);
    for( int j=0; j<nread; j++){
      p[i+j].xvel = cache[3*j+0] * vunit;
      p[i+j].yvel = cache[3*j+1] * vunit;
      p[i+j].zvel = cache[3*j+2] * vunit;
    }
  }
  fread( &blksize, 4, 1, fin);
  fprintf( stderr, "End Reading Gadget file blksize= %d\n", blksize);

  fread( &blksize, 4, 1, fin);
  fprintf( stderr, "Reading Gadget file blksize= %d\n", blksize);

#ifndef LONGIDS
  int *icache = new int[nmemory];
  for( int i=0; i<npart; i+=nmemory){
    int nread = nmemory;
    if( (i+nmemory) > npart)  nread = npart - i;
    fread( icache, sizeof(int), nread, fin);
    for( int j=0; j<nread; j++){
      p[i+j].id = (long long int)icache[j];
    }
  }
#else
  long long int *icache = new long long int[nmemory];
  for( int i=0; i<npart; i+=nmemory){
    int nread = nmemory;
    if( (i+nmemory) > npart)  nread = npart - i;
    fread( icache, sizeof(long long int), nread, fin);
    for( int j=0; j<nread; j++){
      p[i+j].id = icache[j];
    }
  }
#endif
  fread( &blksize, 4, 1, fin);
  fprintf( stderr, "End Reading Gadget file blksize= %d\n", blksize);
  delete [] icache;

  fclose(fin);
  delete [] cache;

  this_run.npart_total = NUMBER_OF_PART_ALL;
  this_run.uniform_mass = 3.0*(this_run.cosm.omega0)/(8*M_PI*(double)(this_run.npart_total));
  this_run.min_mass = this_run.uniform_mass;

}



void outputGadget( RunParam &this_run, pParticle p, const char *outfile){

  int nmemory = IO_CACHE_SIZE;
  int npart = this_run.npart;

  GadgetHeader gadget_header;
  memset( &gadget_header, 0, sizeof(gadget_header));
  gadget_header.Npart[1] = npart;
  gadget_header.Nall[1] = (unsigned int)this_run.npart_total;
  gadget_header.NallHW[1] = (unsigned int) (this_run.npart_total >> 32);
  gadget_header.Massarr[1] = this_run.uniform_mass * this_run.munit * this_run.cosm.hubble /  Gadget_UnitMass_in_Msun;
  gadget_header.Redshift = this_run.znow;
  gadget_header.Time = this_run.anow;
  gadget_header.Omega0 = this_run.cosm.omega0;
  gadget_header.OmegaLambda = this_run.cosm.lambda0;
  gadget_header.HubbleParam = this_run.cosm.hubble;
  gadget_header.NumFiles = this_run.nnode;
  gadget_header.BoxSize = this_run.lunit * this_run.cosm.hubble / Gadget_UnitLength_in_Mpc;

  double lunit_gadget = gadget_header.BoxSize;
  double vunit = this_run.tunit / (this_run.lunit * 1.0e6 * 3.085678e18) * Gadget_UnitVelocity_in_cm_per_s;
  double vunit_gadget = sqrt(this_run.anow) / vunit;
  fprintf( stderr, "%e\t%e\t%e\n", this_run.lunit, this_run.munit, this_run.tunit);
  fprintf( stderr, "%e\t%e\t%e\n", lunit_gadget, gadget_header.Massarr[1], vunit_gadget);

  double anow2inv = 1.0 / this_run.anow / this_run.anow;
  vunit_gadget *= anow2inv;

  FILE *fout = fopen( outfile, "wb");
  int blksize = sizeof(gadget_header);
  fwrite( &blksize, 4, 1, fout);
  fwrite( &gadget_header, sizeof(gadget_header), 1, fout);
  fwrite( &blksize, 4, 1, fout);

  float *cache = new float[nmemory*3];

  blksize = sizeof(float) * 3 * npart;
  fwrite( &blksize, 4, 1, fout);
  for( int i=0; i<npart; i+=nmemory){
    int nwrite = nmemory;
    if( (i+nmemory) > npart)  nwrite = npart - i;
    for( int j=0; j<nwrite; j++){
      cache[3*j+0] = p[i+j].xpos * lunit_gadget;
      cache[3*j+1] = p[i+j].ypos * lunit_gadget;
      cache[3*j+2] = p[i+j].zpos * lunit_gadget;
    }
    fwrite( cache, sizeof(float), 3*nwrite, fout);
  }
  fwrite( &blksize, 4, 1, fout);

  fwrite( &blksize, 4, 1, fout);
  for( int i=0; i<npart; i+=nmemory){
    int nwrite = nmemory;
    if( (i+nmemory) > npart)  nwrite = npart - i;
    for( int j=0; j<nwrite; j++){
      cache[3*j+0] = p[i+j].xvel * vunit_gadget;
      cache[3*j+1] = p[i+j].yvel * vunit_gadget;
      cache[3*j+2] = p[i+j].zvel * vunit_gadget;
    }
    fwrite( cache, sizeof(float), 3*nwrite, fout);
  }
  fwrite( &blksize, 4, 1, fout);

#ifndef LONGIDS
  blksize = sizeof(int) * npart;
#else
  blksize = sizeof(long long int) * npart;
#endif
  fwrite( &blksize, 4, 1, fout);

#ifndef LONGIDS
  int *icache = new int[nmemory];
  for( int i=0; i<npart; i+=nmemory){
    int nwrite = nmemory;
    if( (i+nmemory) > npart)  nwrite = npart - i;
    for( int j=0; j<nwrite; j++){
      icache[j] = (int)p[i+j].id;
    }
    fwrite( icache, sizeof(int), nwrite, fout);
  }
#else
  long long int *icache = new long long int[nmemory];
  for( int i=0; i<npart; i+=nmemory){
    int nwrite = nmemory;
    if( (i+nmemory) > npart)  nwrite = npart - i;
    for( int j=0; j<nwrite; j++){
      icache[j] = p[i+j].id;
    }
    fwrite( icache, sizeof(long long int), nwrite, fout);
  }
#endif

  fwrite( &blksize, 4, 1, fout);
  delete [] icache;

  fclose(fout);
  delete [] cache;

}

    } // namespace ParticleMesh
}     // namespace ParticleSimulator
