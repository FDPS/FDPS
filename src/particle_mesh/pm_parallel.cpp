#include "pm_parallel.h"

namespace ParticleSimulator{
    namespace ParticleMesh{


static char cbuf[CHARMAX];
static char cbuf2[CHARMAX];


void PMForce::PMForce0_fftw3( MPI_Comm comm){

#ifdef FFTW3_PARALLEL
  ptrdiff_t alloc_local, local_n0, local_0_start;
  alloc_local = fftw_mpi_local_size_3d( SIZE_MESH, 
					SIZE_MESH,
					SIZE_MESH/2+1,
					comm,
					&local_n0,
					&local_0_start);
  alloc_local *= 2;
  local_nx = local_n0;
  local_x_start = local_0_start;
  local_ny_after_transpose = local_nx;
  local_y_start_after_transpose = local_x_start;
  total_local_size = alloc_local;
#endif

}



void PMForce::PMForce1_fftw3( MPI_Comm comm, 
			      const char *fname_wisdom_f, const char *fname_wisdom_b){

#ifdef FFTW3_PARALLEL

  if( inode_fft == 0){
    cerr << "import wisdom: " << fname_wisdom_f << endl;
    fftw_import_wisdom_from_filename( fname_wisdom_f);
  }
  fftw_mpi_broadcast_wisdom( comm);

  fplan = fftw_mpi_plan_dft_r2c_3d( SIZE_MESH,
				    SIZE_MESH,
				    SIZE_MESH,
				    mesh_density_slab,
				    (fftw_complex *)mesh_density_slab,
				    comm,
				    FFTW_WISDOM_ONLY | FFTW_MPI_TRANSPOSED_OUT);

  fftw_forget_wisdom();
  if( inode_fft == 0){
    cerr << "import wisdom: " << fname_wisdom_b << endl;
    fftw_import_wisdom_from_filename( fname_wisdom_b);
  }
  fftw_mpi_broadcast_wisdom( comm);

  bplan = fftw_mpi_plan_dft_c2r_3d( SIZE_MESH,
				    SIZE_MESH,
				    SIZE_MESH,
				    (fftw_complex *)mesh_density_slab,
				    mesh_density_slab,
				    comm,
				    FFTW_WISDOM_ONLY  | FFTW_MPI_TRANSPOSED_IN);
#endif

}



void PMForce::PMForce2_fftw3( MPI_Comm comm){

#ifdef FFTW3_PARALLEL

  if( inode_fft == 0){
    cerr << "wisdom is not available" << endl;
  }

  fplan = fftw_mpi_plan_dft_r2c_3d( SIZE_MESH,
				    SIZE_MESH,
				    SIZE_MESH,
				    mesh_density_slab,
				    (fftw_complex *)mesh_density_slab,
				    comm,
				    FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);

  bplan = fftw_mpi_plan_dft_c2r_3d( SIZE_MESH,
				    SIZE_MESH,
				    SIZE_MESH,
				    (fftw_complex *)mesh_density_slab,
				    mesh_density_slab,
				    comm,
				    FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
#endif

}



void PMForce::initPMReduce( const int x, const int y, const int z){

#ifdef RMM_PM

  int *ndiv = this_run->ndiv;

  reduce_nx = _pm_reduce_nx;
  reduce_ny = _pm_reduce_ny;
  reduce_nz = _pm_reduce_nz;

  ndiv_reduce[0] = ndiv[0] / reduce_nx;
  ndiv_reduce[1] = ndiv[1] / reduce_ny;
  ndiv_reduce[2] = ndiv[2] / reduce_nz;
  int column_x = x / ndiv_reduce[0];
  int column_y = y / ndiv_reduce[1];
  int column_z = z / ndiv_reduce[2];

  int xc = x - column_x*ndiv_reduce[0];
  int yc = y - column_y*ndiv_reduce[1];
  int zc = z - column_z*ndiv_reduce[2];

#if 1
  if( column_x == reduce_nx){
    column_x -= 1;
    xc = x - column_x*ndiv_reduce[0];
  }
  if( column_x == reduce_nx -1){
    ndiv_reduce[0] += ndiv[0] % ndiv_reduce[0];
  }

  if( column_y == reduce_ny){
    column_y -= 1;
    yc = y - column_y*ndiv_reduce[1];
  }
  if( column_y == reduce_ny -1){
    ndiv_reduce[1] += ndiv[1] % ndiv_reduce[1];
  }

  if( column_z == reduce_nz){
    column_z -= 1;
    zc = z - column_z*ndiv_reduce[2];
  }
  if( column_z == reduce_nz -1){
    ndiv_reduce[2] += ndiv[2] % ndiv_reduce[2];
  }
#endif

  int color2 = column_x*reduce_nz*reduce_ny + column_y*reduce_nz + column_z;
  int key2 = getVoxelIndex( xc, yc, zc, ndiv_reduce);
  if( init_gk == 1){
    assert( ndiv_fft[0] <= ndiv_reduce[0]);
    assert( ndiv_fft[1] <= ndiv_reduce[1]);
    assert( ndiv_fft[2] <= ndiv_reduce[2]);
  }

  inode_reduce = key2;
  reduce_node = key2;
  MPI_Comm_split( MPI_COMM_WORLD, color2, key2, &MPI_COMM_COLUMN);
  MPI_Comm_split( MPI_COMM_WORLD, key2, color2, &MPI_COMM_REDUCE);
  nnode_reduce = ndiv_reduce[0] * ndiv_reduce[1] * ndiv_reduce[2];

#if 0
  sprintf( cbuf, "x,y,z,color,key,xc,yc,zc,ndiv_reduce %d %d %d\t%d %d\t%d %d %d %d %d %d\n", x, y, z, color2, key2, 
	   xc, yc, zc, ndiv_reduce[0], ndiv_reduce[1], ndiv_reduce[2]);
  mp_print( cbuf, MPI_COMM_WORLD, stdout);
#endif

#endif

}



void PMForce::PMForce0( const int *ndiv_fft_input){


#ifdef FFTW3_PARALLEL
  fftw_mpi_init();
#endif

  ndiv_fft[0] = ndiv_fft_input[0];
  ndiv_fft[1] = ndiv_fft_input[1];
  ndiv_fft[2] = ndiv_fft_input[2];
#ifdef FIX_FFTNODE
  int nnode_fft0 = ndiv_fft[0] * ndiv_fft[1] * ndiv_fft[2];
#ifndef FFT3D
  assert( SIZE_MESH % nnode_fft0 == 0);
  assert( SIZE_MESH >= nnode_fft0);
#endif
  int x, y, z;
  getXYZIndex( this_run->inode, this_run->ndiv, &x, &y, &z);
  int color = 0;
  int key = 0;
  if( x < ndiv_fft[0] &&  y < ndiv_fft[1]  &&  z < ndiv_fft[2]){
    color = 1;
    key = getVoxelIndex( x, y, z, ndiv_fft);
  }
  else{
    key = this_run->inode;
    local_nx = local_x_start = local_ny_after_transpose = local_y_start_after_transpose = 0;
    total_local_size = 0;
  }
  MPI_Comm_split( MPI_COMM_WORLD, color, key, &MPI_COMM_FFT);
  fft_node = color;

#ifdef RMM_PM
  initPMReduce( x, y, z);
#endif

  MPI_Comm_size(MPI_COMM_FFT,&nnode_fft);
  MPI_Comm_rank(MPI_COMM_FFT,&inode_fft);
  int x2, y2, z2;
  getXYZIndex( inode_fft, ndiv_fft, &x2, &y2, &z2);
  sprintf( cbuf, "%d->%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", this_run->inode, inode_fft, key,
	   x, y, z, x2, y2, z2);
  //mp_print( cbuf, MPI_COMM_FFT, stderr);
#else
  inode_fft = this_run->inode;
  nnode_fft = this_run->nnode;
  fft_node = 1;
  MPI_COMM_FFT = this_run->MPI_COMM_INTERNAL;
#endif

  if( fft_node == 1){
#ifdef FFTW3_PARALLEL
    PMForce0_fftw3( MPI_COMM_FFT);
#else
    fplan = rfftw3d_mpi_create_plan( MPI_COMM_FFT,
				     SIZE_MESH,
				     SIZE_MESH,
				     SIZE_MESH,
				     FFTW_REAL_TO_COMPLEX,
				     FFTW_ESTIMATE);
    // FFTW_MEASURE);
    bplan = rfftw3d_mpi_create_plan( MPI_COMM_FFT,
				     SIZE_MESH,
				     SIZE_MESH,
				     SIZE_MESH,
				     FFTW_COMPLEX_TO_REAL,
				     FFTW_ESTIMATE);
				     //				     FFTW_MEASURE);

    rfftwnd_mpi_local_sizes( fplan, &local_nx, &local_x_start,
			     &local_ny_after_transpose,
			     &local_y_start_after_transpose,
			     &total_local_size);
#endif
  }


#ifdef RMM_PM
  local_nx2 = local_nx;
  local_x_start2 = local_x_start;
  MPI_Bcast( &local_nx2, 1, MPI_INT, 0, MPI_COMM_REDUCE);
  MPI_Bcast( &local_x_start2, 1, MPI_INT, 0, MPI_COMM_REDUCE);
#endif

}




void PMForce::init( RunParam *i_this_run, const int _SIZE_MESH, 
		    const int *_ndiv_fft_input, const int flag_wisdom){

  double nowtime;
  getTime(&nowtime);
  this_run = i_this_run;

  SIZE_MESH = _SIZE_MESH;
  SIZE_MESH_P2 = SIZE_MESH + 2;

#ifdef FIX_FFTNODE
  if( this_run->inode == 0)  cerr << "_ndiv_fft= " 
				  << _ndiv_fft_input[0] << "\t"
				  << _ndiv_fft_input[1] << "\t"
				  << _ndiv_fft_input[2] << endl;
  PMForce0( _ndiv_fft_input);
#else
  if( this_run->inode == 0)  cerr << "_ndiv_fft= " 
				  << this_run->ndiv[0] << "\t"
				  << this_run->ndiv[1] << "\t"
				  << this_run->ndiv[2] << endl;
  PMForce0( this_run->ndiv);
#endif
  getTimePrint( &nowtime, "initFFTW");

#if 0
  sprintf( cbuf, "%d\t%d\t%d\t%d\t%d\n", 
	   local_nx, local_x_start, local_ny_after_transpose,
	   local_y_start_after_transpose, total_local_size);
  mp_print( cbuf, MPI_COMM_WORLD, stderr);
#endif

#ifdef RMM_PM
  l_nx_all2 = new int[nnode_reduce];
  l_x_start_all2 = new int[nnode_reduce];
  MPI_Allgather( &local_nx2, 1, MPI_INT, 
		 l_nx_all2, 1, MPI_INT, MPI_COMM_COLUMN);
  MPI_Allgather( &local_x_start2, 1, MPI_INT, 
		 l_x_start_all2, 1, MPI_INT, MPI_COMM_COLUMN);
#endif

  l_nx_all = new int[this_run->nnode];
  l_x_start_all = new int[this_run->nnode];
  MPI_Allgather( &local_nx, 1, MPI_INT, 
		 l_nx_all, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);
  MPI_Allgather( &local_x_start, 1, MPI_INT, 
		 l_x_start_all, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);

  if( fft3d == 1){
#ifdef FFT3D
    assert( ndiv_fft[0] == ndiv_fft[1]);
    assert( ndiv_fft[0] == ndiv_fft[2]);
    subbox_msize[0] = SIZE_MESH / ndiv_fft[0];
    subbox_msize[1] = subbox_msize[0];
    subbox_msize[2] = subbox_msize[0];
    subbox_total = subbox_msize[0] * subbox_msize[1] * subbox_msize[2];
    msize_slab = subbox_total;
    int x, y, z;
    getXYZIndex( inode_reduce, ndiv_reduce, &x, &y, &z);
    subbox_start[0] = x * subbox_msize[0];
    subbox_start[1] = y * subbox_msize[1];
    subbox_start[2] = z * subbox_msize[2];
    mesh_density_subbox = new TYPE_FFTW[msize_slab][2]; 
    mesh_density_slab0 = new TYPE_FFTW[msize_slab]; 
    fft3d_work = new TYPE_FFTW[msize_slab][2]; 
    MPI_COMM_FFT_F = MPI_Comm_c2f(MPI_COMM_FFT);
#if 0
    sprintf( cbuf, "%d\t%d\t%d\t%d\t%d\t%d\n", 
	     subbox_msize[0], subbox_msize[1], subbox_msize[2],
	     subbox_start[0], subbox_start[1], subbox_start[2]);
    mp_print( cbuf, MPI_COMM_WORLD, stderr);
#endif
#endif
  }
  else{
#ifdef RMM_PM
    msize_slab = local_nx2 * SIZE_MESH * SIZE_MESH_P2;
    mesh_density_slab0 = new TYPE_FFTW[msize_slab]; 
#else  //RMM_PM
    msize_slab = local_nx * SIZE_MESH * SIZE_MESH_P2;
#endif //RMM_PM
  }

  mesh_density_slab = new TYPE_FFTW[msize_slab]; 
  fftwork = new TYPE_FFTW[msize_slab]; 
  getTimePrint( &nowtime, "allocFFT");

#ifdef FFTW3_PARALLEL
  if( flag_wisdom == 1){
    PMForce1_fftw3( MPI_COMM_FFT, FFTW_WISDOM_FILENAME_F, FFTW_WISDOM_FILENAME_B);
  }
  else{
    PMForce2_fftw3( MPI_COMM_FFT);
  }
  getTimePrint( &nowtime, "initFFTW");
#endif

  memset( mesh_density_slab, 0, sizeof(TYPE_FFTW)*msize_slab);

  int sizeg = SIZE_MESH/2 + 1;
  if( fft3d == 1){
    msize_gk = subbox_total;
  }
  else{
    msize_gk = sizeg * sizeg * local_ny_after_transpose;
  }
  gk = new TYPE_FFTW[msize_gk];
  if( init_gk == 1){
    memset( gk, 0, sizeof(TYPE_FFTW)*msize_gk);
    getTime( &nowtime);
    if( fft3d == 1){
      initGK3D();
    }
    else{
      initGK();
    }
    getTimePrint( &nowtime, "initGK");
  }

  output_mesh_flag = 0;

}



PMForce::PMForce( RunParam *i_this_run){

  fft3d = 0;
#ifdef FFT3D
  fft3d = 1;
#endif
  init_gk = 1;
  init( i_this_run, SIZE_OF_MESH, _ndiv_fft, _flag_wisdom);

}






PMForce::~PMForce(){

#ifdef FFTW3_PARALLEL
  fftw_destroy_plan( fplan);
  fftw_destroy_plan( bplan);
#else
  if( fft_node == 1){
    rfftwnd_mpi_destroy_plan( fplan);
    rfftwnd_mpi_destroy_plan( bplan);
  }
#endif

#ifdef RMM_PM
  delete [] mesh_density_slab0;
  delete [] l_nx_all2;
  delete [] l_x_start_all2;
#endif

  if( fft3d == 1){
    delete [] mesh_density_subbox;
    delete [] fft3d_work;
  }

  delete [] mesh_density_slab;
  delete [] fftwork;;
  delete [] l_nx_all;
  delete [] l_x_start_all;
  delete [] gk;

}



double sksq( float wn, float a){
   
   double akb2 = (double)wn*(double)a*5.0e-1;
   double sk   = 12.e0*(pow(2.e0*sin(akb2/2.e0),2.0)-akb2*sin(akb2))/pow(akb2,4.0);
   return sk*sk/(double)(wn*wn);

}



void PMForce::initGK3D(){
#ifdef FFT3D

#ifdef _OPENMP  
  omp_set_num_threads(NUMBER_OF_OMP_THREADS);
#endif

  int sizeg = SIZE_MESH/2 + 1;

  float *sins1 = new float[sizeg];
  float *sins2 = new float[sizeg];
  float *dfi = new float[sizeg];

  double h  = 1.e0/(float)SIZE_MESH;
  double hi = (float)SIZE_MESH;
  double alpha = 4.0/3.0;
  double pi2 = 2.0 * M_PI;
  double pi2hi = 2.0 * M_PI * hi;
  double pi4 = 4.0 * M_PI;

#pragma omp parallel for
  for( int mindx=0; mindx<sizeg; mindx++) {
    float hkb2 = (float)mindx * M_PI * h;
    dfi[mindx] = alpha*sin(hkb2*2.0)*hi + (1.0-alpha)*sin(4.0*hkb2)*hi*0.5;
    sins2[mindx] = sin(hkb2) * hi * 2.0;
    sins1[mindx] = sin(hkb2);
    sins1[mindx] = 1.0 - sins1[mindx]*sins1[mindx] + 0.133333333*pow(sins1[mindx],4);
    sins1[mindx] = sins1[mindx]*sins1[mindx];
  }

   
  for(int i=0; i<subbox_msize[0]; i++){
    int ii = i + subbox_start[0];
    if( ii > SIZE_MESH/2)  ii = SIZE_MESH - ii;
    float xk = (float)(ii) * pi2;
#pragma omp parallel for
    for(int j=0; j<subbox_msize[1]; j++){
      int jj = j + subbox_start[1];
      if( jj > SIZE_MESH/2)  jj = SIZE_MESH - jj;
      float yk = (float)(jj) * pi2;
      for(int k=0; k<subbox_msize[2]; k++){
	int kk = k + subbox_start[2];
	if( kk > SIZE_MESH/2)  kk = SIZE_MESH - kk;
	float zk = (float)(kk) * pi2;
	int li = k + subbox_msize[2]*(j + subbox_msize[1]*i);
	int n1s = -2;
	if(ii==0)  n1s = 1;
	int n2s = -2;
	if(jj==0)  n2s = 1;
	int n3s = -2;
	if(kk==0)  n3s = 1;
	for(int n1=n1s; n1<=1; n1++){
	  for(int n2=n2s; n2<=1; n2++){
	    for(int n3=n3s; n3<=1; n3++){
	      float wnx = xk + (float)(n1) * pi2hi;
	      float wny = yk + (float)(n2) * pi2hi;
	      float wnz = zk + (float)(n3) * pi2hi;
	      if(ii==0)  wnx = 0;
	      if(jj==0)  wny = 0;
	      if(kk==0)  wnz = 0;
	      float wn  = sqrt((wnx*wnx + wny*wny + wnz*wnz));
	      float ukx = sins2[ii] / wnx;
	      float uky = sins2[jj] / wny;
	      float ukz = sins2[kk] / wnz;
	      if(ii==0)  ukx = 1;
	      if(jj==0)  uky = 1;
	      if(kk==0)  ukz = 1;
	      float uknsq = ukx * uky * ukz;
	      uknsq = pow(uknsq,6);
	      float dkkn  = dfi[ii]*wnx + dfi[jj]*wny + dfi[kk]*wnz;
	      gk[li] += dkkn * sksq(wn,SFT_FOR_PM) * uknsq;
	    }
	  }
	} 
	float ukx = sins1[ii];
	float uky = sins1[jj];
	float ukz = sins1[kk];
	if(ii==0)  ukx = 1;
	if(jj==0)  uky = 1;
	if(kk==0)  ukz = 1;
	float usqu = ukx * uky * ukz;

	float dkx = dfi[ii];
	float dky = dfi[jj];
	float dkz = dfi[kk];
	if(ii==0)  dkx = 0.0;
	if(jj==0)  dky = 0.0;
	if(kk==0)  dkz = 0.0;
	float dksqu = dkx*dkx + dky*dky + dkz*dkz;
	gk[li] /= -(dksqu*usqu);
	gk[li] *= pi4;
      }
    }
  }

  if( subbox_start[0] == 0  && subbox_start[1] == 0  && subbox_start[2] ==0){
    gk[0] = 0.0;
  }

  delete [] sins1;
  delete [] sins2;
  delete [] dfi;

#endif

}



void PMForce::initGK(){

#ifdef _OPENMP  
  omp_set_num_threads(NUMBER_OF_OMP_THREADS);
#endif

  int sizeg = SIZE_MESH/2 + 1;

  float *sins1 = new float[sizeg];
  float *sins2 = new float[sizeg];
  float *dfi = new float[sizeg];

  double h  = 1.e0/(float)SIZE_MESH;
  double hi = (float)SIZE_MESH;
  double alpha = 4.0/3.0;
  double pi2 = 2.0 * M_PI;
  double pi2hi = 2.0 * M_PI * hi;
  double pi4 = 4.0 * M_PI;

#pragma omp parallel for
  for( int mindx=0; mindx<sizeg; mindx++) {
    float hkb2 = (float)mindx * M_PI * h;
    dfi[mindx] = alpha*sin(hkb2*2.0)*hi + (1.0-alpha)*sin(4.0*hkb2)*hi*0.5;
    sins2[mindx] = sin(hkb2) * hi * 2.0;
    sins1[mindx] = sin(hkb2);
    sins1[mindx] = 1.0 - sins1[mindx]*sins1[mindx] + 0.133333333*pow(sins1[mindx],4);
    sins1[mindx] = sins1[mindx]*sins1[mindx];
  }

   
  for(int i=0; i<local_ny_after_transpose; i++){
    int ii = i + local_y_start_after_transpose;
    if( ii > SIZE_MESH/2)  ii = SIZE_MESH - ii;
    float xk = (float)(ii) * pi2;
#pragma omp parallel for
    for(int j=0; j<sizeg; j++){
      int jj = j;
      float yk = (float)(jj) * pi2;
      for(int k=0; k<sizeg; k++){
	int kk = k;
	float zk = (float)(kk) * pi2;
	int li = k + sizeg*(j + sizeg*i);
	int n1s = -2;
	if(ii==0)  n1s = 1;
	int n2s = -2;
	if(jj==0)  n2s = 1;
	int n3s = -2;
	if(kk==0)  n3s = 1;
	for(int n1=n1s; n1<=1; n1++){
	  for(int n2=n2s; n2<=1; n2++){
	    for(int n3=n3s; n3<=1; n3++){
	      float wnx = xk + (float)(n1) * pi2hi;
	      float wny = yk + (float)(n2) * pi2hi;
	      float wnz = zk + (float)(n3) * pi2hi;
	      if(ii==0)  wnx = 0;
	      if(jj==0)  wny = 0;
	      if(kk==0)  wnz = 0;
	      float wn  = sqrt((wnx*wnx + wny*wny + wnz*wnz));
	      float ukx = sins2[ii] / wnx;
	      float uky = sins2[jj] / wny;
	      float ukz = sins2[kk] / wnz;
	      if(ii==0)  ukx = 1;
	      if(jj==0)  uky = 1;
	      if(kk==0)  ukz = 1;
	      float uknsq = ukx * uky * ukz;
	      uknsq = pow(uknsq,6);
	      float dkkn  = dfi[ii]*wnx + dfi[jj]*wny + dfi[kk]*wnz;
	      gk[li] += dkkn * sksq(wn,SFT_FOR_PM) * uknsq;
	    }
	  }
	} 
	float ukx = sins1[ii];
	float uky = sins1[jj];
	float ukz = sins1[kk];
	if(ii==0)  ukx = 1;
	if(jj==0)  uky = 1;
	if(kk==0)  ukz = 1;
	float usqu = ukx * uky * ukz;

	float dkx = dfi[ii];
	float dky = dfi[jj];
	float dkz = dfi[kk];
	if(ii==0)  dkx = 0.0;
	if(jj==0)  dky = 0.0;
	if(kk==0)  dkz = 0.0;
	float dksqu = dkx*dkx + dky*dky + dkz*dkz;
	gk[li] /= -(dksqu*usqu);
	gk[li] *= pi4;
      }
    }
  }

  if( local_y_start_after_transpose == 0){
    gk[0] = 0.0;
  }

  delete [] sins1;
  delete [] sins2;
  delete [] dfi;

}



void PMForce::initLocal(){

  local_to_slab_recvbuf_size = slab_to_local_sendbuf_size = 0;

  memset( mesh_density_slab, 0, sizeof(TYPE_FFTW)*msize_slab);

  getLocalMeshParam(3);

  mesh_density_local = new float[ln_total];
  memset( mesh_density_local, 0, sizeof(float)*ln_total);

  mesh_force_local = new float[ln_total][3];

  if( output_mesh_flag == 1){
    mesh_density_local2 = new float[ln_total];
    memset( mesh_density_local2, 0, sizeof(float)*ln_total);
    mesh_density_slab2 = new float[msize_slab]; 
  }

  MPI_Barrier(MPI_COMM_WORLD);
  getLocalMeshParam(1);
  
#ifdef RMM_PM
  memset( mesh_density_slab0, 0, sizeof(TYPE_FFTW)*msize_slab);
#endif


}


void PMForce::delLocal(){

  delete [] mesh_density_local;
  delete [] mesh_force_local;

}



void PMForce::checkMeshSlab(){

  static int msize = SIZE_MESH * SIZE_MESH * SIZE_MESH_P2;
  static int msize_yz = SIZE_MESH * SIZE_MESH_P2;

  float *mesh_density = (float *) malloc( sizeof(float) * msize);
  float *mesh_density_buf = (float *) malloc( sizeof(float) * msize);
  for( int i=0; i<msize; i++){
    mesh_density[i] = mesh_density_buf[i] = 0.0;
  }

  if( fft3d == 1){
    for( int i=0; i<subbox_msize[0]; i++){
      int ii = i + subbox_start[0];
      for( int j=0; j<subbox_msize[1]; j++){
	int jj = j + subbox_start[1];
	for( int k=0; k<subbox_msize[2]; k++){
	  int kk = k + subbox_start[2];
	  int li = kk + SIZE_MESH_P2*( jj + SIZE_MESH*ii);
	  int ri = getSubboxID( i, j, k);
	  mesh_density_buf[li] = mesh_density_slab[ri];
	  //mesh_density_buf[li] = mesh_density_subbox[ri][0];
	}
      }
    }
  }
  else{
    int offset = local_x_start * msize_yz;
    for( int i=0; i<msize_slab; i++){
      mesh_density_buf[i+offset] = mesh_density_slab[i];
    }
  }

#ifdef FIX_FFTNODE
  MPI_Reduce( mesh_density_buf, mesh_density, msize, 
	      MPI_FLOAT, MPI_SUM, 0, MPI_COMM_FFT);
#else
  MPI_Reduce( mesh_density_buf, mesh_density, msize, 
	      MPI_FLOAT, MPI_SUM, 0, this_run->MPI_COMM_INTERNAL);
#endif

  if( this_run->inode == 0){
    for( int i=0; i<msize; i++){
//      fprintf( stdout, "%.8e\n", mesh_density[i]);
      fprintf( stderr, "slab %+.16e\n", mesh_density[i]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  free( mesh_density);
  free( mesh_density_buf);

}



void PMForce::checkMeshGreen(){

  int sizeg = SIZE_MESH/2 + 1;
  int msize = sizeg * sizeg * SIZE_MESH;
  int msize_yz = sizeg * sizeg;
  if( fft3d == 1){
    msize = SIZE_MESH * SIZE_MESH * SIZE_MESH;
    msize_yz = SIZE_MESH * SIZE_MESH;
  }

  float *mesh_density = (float *) malloc( sizeof(float) * msize);
  float *mesh_density_buf = (float *) malloc( sizeof(float) * msize);
  for( int i=0; i<msize; i++){
    mesh_density[i] = mesh_density_buf[i] = 0.0;
  }

  if( fft3d == 1){
    for( int i=0; i<subbox_msize[0]; i++){
      int ii = i + subbox_start[0];
      for( int j=0; j<subbox_msize[1]; j++){
	int jj = j + subbox_start[1];
	for( int k=0; k<subbox_msize[2]; k++){
	  int kk = k + subbox_start[2];
	  int gi = k + subbox_msize[2]*( j + i*subbox_msize[1]);
	  int mi = kk + SIZE_OF_MESH*( jj + ii*SIZE_OF_MESH);
	  mesh_density_buf[mi] = gk[gi];
	}
      }
    }
  }
  else{
    int offset = local_x_start * msize_yz;
    for( int i=0; i<local_nx; i++){
      for( int j=0; j<sizeg; j++){
	for( int k=0; k<sizeg; k++){
	  int gi = k + sizeg*( j + sizeg*i);
	  int mi = offset + k + sizeg*( j + i*sizeg);
	  mesh_density_buf[mi] = gk[gi];
	}
      }
    }
  }

  MPI_Reduce( mesh_density_buf, mesh_density, msize, 
	      MPI_FLOAT, MPI_SUM, 0, this_run->MPI_COMM_INTERNAL);
  
  if( this_run->inode == 0){
    if( fft3d == 1){
      for( int i=0; i<sizeg; i++){
	for( int j=0; j<sizeg; j++){
	  for( int k=0; k<sizeg; k++){
	    int ii = k + SIZE_MESH*( j + SIZE_MESH*i);
	    fprintf( stdout, "%.8e\n", mesh_density[ii]);
	  }
	}
      }
    }
    else{
      /*      for( int i=0; i<msize; i++){
	fprintf( stdout, "%.8e\n", mesh_density[i]);
	}*/
      for( int i=0; i<sizeg; i++){
	for( int j=0; j<sizeg; j++){
	  for( int k=0; k<sizeg; k++){
	    int ii = k + sizeg*( j + sizeg*i);
	    fprintf( stdout, "%.8e\n", mesh_density[ii]);
	  }
	}
      }
    }
  }


  free( mesh_density);
  free( mesh_density_buf);

}



void PMForce::checkMeshLocal(){

  static int msize = SIZE_MESH * SIZE_MESH * SIZE_MESH;

  float *mesh_density = (float *) malloc( sizeof(float) * msize);
  float *mesh_density_buf = (float *) malloc( sizeof(float) * msize);
  for( int i=0; i<msize; i++){
    mesh_density[i] = mesh_density_buf[i] = -1.0e30;
  }

  int ii = 0;
  for( int i=0; i<l_msize[0]; i++){
    for( int j=0; j<l_msize[1]; j++){
      for( int k=0; k<l_msize[2]; k++){
	int i2 = i + g_offset[0];
	if( i2 < SIZE_MESH)  i2 += SIZE_MESH;
	if( i2 >= SIZE_MESH)  i2 -= SIZE_MESH;
	int j2 = j + g_offset[1];
	if( j2 < SIZE_MESH)  j2 += SIZE_MESH;
	if( j2 >= SIZE_MESH)  j2 -= SIZE_MESH;
	int k2 = k + g_offset[2];
	if( k2 < SIZE_MESH)  k2 += SIZE_MESH;
	if( k2 >= SIZE_MESH)  k2 -= SIZE_MESH;
	int iii = k2 + SIZE_MESH*j2 + SIZE_MESH*SIZE_MESH*i2;
	mesh_density_buf[iii] = mesh_density_local[ii];
	ii ++;
      }
    }
  }

  MPI_Reduce( mesh_density_buf, mesh_density, msize, 
	      MPI_FLOAT, MPI_MAX, 0, this_run->MPI_COMM_INTERNAL);

  if( this_run->inode == 0){
    for( int i=0; i<msize; i++){
//      fprintf( stdout, "%.8e\n", mesh_density[i]);
      fprintf( stderr, "hoge %+.16e\n", mesh_density[i]);
    }
  }

  free( mesh_density);
  free( mesh_density_buf);

}



void PMForce::getLocalMeshParam( const double padsize){

  double pad = padsize / SIZE_MESH;

  double gx = (int)(this_run->bmin[0]*(double)SIZE_MESH)/(double)SIZE_MESH-pad;
  double gy = (int)(this_run->bmin[1]*(double)SIZE_MESH)/(double)SIZE_MESH-pad;
  double gz = (int)(this_run->bmin[2]*(double)SIZE_MESH)/(double)SIZE_MESH-pad;

  double gx2 = (int)(this_run->bmax[0]*(double)SIZE_MESH+0.5)/(double)SIZE_MESH+pad;
  double gy2 = (int)(this_run->bmax[1]*(double)SIZE_MESH+0.5)/(double)SIZE_MESH+pad;
  double gz2 = (int)(this_run->bmax[2]*(double)SIZE_MESH+0.5)/(double)SIZE_MESH+pad;

  g_pos[0] = gx;
  g_pos[1] = gy;
  g_pos[2] = gz;
  g_pos_f[0] = gx;
  g_pos_f[1] = gy;
  g_pos_f[2] = gz;

  g_offset[0] = (int)(gx*(double)SIZE_MESH);
  g_offset[1] = (int)(gy*(double)SIZE_MESH);
  g_offset[2] = (int)(gz*(double)SIZE_MESH);

  l_msize[0] = (int)ceil((gx2 - gx)*SIZE_MESH) + 1;
  l_msize[1] = (int)ceil((gy2 - gy)*SIZE_MESH) + 1;
  l_msize[2] = (int)ceil((gz2 - gz)*SIZE_MESH) + 1;

  lnyz = l_msize[1] * l_msize[2];
  ln_total = l_msize[0] * lnyz;

#if 0
  for( int i=0; i<this_run->nnode; i++){
    if( i == this_run->inode){
      cerr << "bmin\t" << this_run->bmin[0] << "\t" << this_run->bmin[1] << "\t" 
	   << this_run->bmin[2] << endl;
      cerr << "gx\t" << g_pos[0] << "\t" << g_pos[1] << "\t" << g_pos[2] << endl;
      cerr << "gix_offset\t" << g_offset[0] 
	   << "\t" << g_offset[1] << "\t" << g_offset[2] << endl;
      cerr << "lnx\t" << l_msize[0] << "\t" << l_msize[1] 
	   << "\t" << l_msize[2] << endl << endl;
    }
    MPI_Barrier( MPI_COMM_WORLD);
  }
#endif

}



static inline int addr(int ix, int iy, int iz,
                       int nx, int ny, int nz){
  return iz + nz * (iy + ny * ix);
}



void PMForce::setLocalMeshDensity( const Particle *particle, const int npart){

#ifdef _OPENMP
  omp_set_num_threads(NUMBER_OF_OMP_THREADS);

  int chunk_size = ln_total * NUMBER_OF_OMP_THREADS;
  float *mesh_density_local_chunk = new float[chunk_size];
  memset( mesh_density_local_chunk, 0, sizeof(float)*chunk_size);
  float *mesh_density_local_thread[NUMBER_OF_OMP_THREADS];
  for( int i=0; i<NUMBER_OF_OMP_THREADS; i++){
    mesh_density_local_thread[i] = &mesh_density_local_chunk[i*ln_total];
  }


#pragma omp parallel 
  {
    int devid = omp_get_thread_num();
#ifdef UNSTABLE
    float g_mass = this_run->uniform_mass;
    const double gx = g_pos[0], gy = g_pos[1], gz = g_pos[2];
    const int nx = l_msize[0];
    const int ny = l_msize[1];
    const int nz = l_msize[2];
    float *rho = mesh_density_local_thread[devid];
#pragma omp for
    for(int p=0; p<npart; p++){
      float x = (float)(particle[p].xpos - gx);
      float y = (float)(particle[p].ypos - gy);
      float z = (float)(particle[p].zpos - gz);

      float xt1 = x * (float)SIZE_MESH;
      float yt1 = y * (float)SIZE_MESH;
      float zt1 = z * (float)SIZE_MESH;

      int ix = (int)(xt1 + 0.5f);
      int iy = (int)(yt1 + 0.5f);
      int iz = (int)(zt1 + 0.5f);

#define LOAD(i,j,k) float rho##i##j##k = rho[ addr(ix+i-1, iy+j-1, iz+k-1, nx, ny, nz) ]
      LOAD(0,0,0);
      LOAD(1,0,0);
      LOAD(2,0,0);
      LOAD(0,1,0);
      LOAD(1,1,0);
      LOAD(2,1,0);
      LOAD(0,2,0);
      LOAD(1,2,0);
      LOAD(2,2,0);
      LOAD(0,0,1);
      LOAD(1,0,1);
      LOAD(2,0,1);
      LOAD(0,1,1);
      LOAD(1,1,1);
      LOAD(2,1,1);
      LOAD(0,2,1);
      LOAD(1,2,1);
      LOAD(2,2,1);
      LOAD(0,0,2);
      LOAD(1,0,2);
      LOAD(2,0,2);
      LOAD(0,1,2);
      LOAD(1,1,2);
      LOAD(2,1,2);
      LOAD(0,2,2);
      LOAD(1,2,2);
      LOAD(2,2,2);
#undef LOAD

      float dx = xt1 - (float)(ix);
      float dy = yt1 - (float)(iy);
      float dz = zt1 - (float)(iz);

      float wx0 = 0.5f  * (0.5f-dx) * (0.5f-dx);
      float wx1 = 0.75f - (     dx) * (     dx);
      float wx2 = 0.5f  * (0.5f+dx) * (0.5f+dx);

      float wy0 = 0.5f  * (0.5f-dy) * (0.5f-dy);
      float wy1 = 0.75f - (     dy) * (     dy);
      float wy2 = 0.5f  * (0.5f+dy) * (0.5f+dy);

      float wz0 = 0.5f  * (0.5f-dz) * (0.5f-dz);
      float wz1 = 0.75f - (     dz) * (     dz);
      float wz2 = 0.5f  * (0.5f+dz) * (0.5f+dz);

      float m = (float(g_mass));

      rho000 += m*wx0*wy0*wz0;
      rho100 += m*wx1*wy0*wz0;
      rho200 += m*wx2*wy0*wz0;
      rho010 += m*wx0*wy1*wz0;
      rho110 += m*wx1*wy1*wz0;
      rho210 += m*wx2*wy1*wz0;
      rho020 += m*wx0*wy2*wz0;
      rho120 += m*wx1*wy2*wz0;
      rho220 += m*wx2*wy2*wz0;
      rho001 += m*wx0*wy0*wz1;
      rho101 += m*wx1*wy0*wz1;
      rho201 += m*wx2*wy0*wz1;
      rho011 += m*wx0*wy1*wz1;
      rho111 += m*wx1*wy1*wz1;
      rho211 += m*wx2*wy1*wz1;
      rho021 += m*wx0*wy2*wz1;
      rho121 += m*wx1*wy2*wz1;
      rho221 += m*wx2*wy2*wz1;
      rho002 += m*wx0*wy0*wz2;
      rho102 += m*wx1*wy0*wz2;
      rho202 += m*wx2*wy0*wz2;
      rho012 += m*wx0*wy1*wz2;
      rho112 += m*wx1*wy1*wz2;
      rho212 += m*wx2*wy1*wz2;
      rho022 += m*wx0*wy2*wz2;
      rho122 += m*wx1*wy2*wz2;
      rho222 += m*wx2*wy2*wz2;
#define STORE(i,j,k) rho[ addr(ix+i-1, iy+j-1, iz+k-1, nx, ny, nz) ] = rho##i##j##k
      STORE(0,0,0);
      STORE(1,0,0);
      STORE(2,0,0);
      STORE(0,1,0);
      STORE(1,1,0);

      STORE(2,1,0);
      STORE(0,2,0);
      STORE(1,2,0);
      STORE(2,2,0);
      STORE(0,0,1);
      STORE(1,0,1);
      STORE(2,0,1);
      STORE(0,1,1);
      STORE(1,1,1);
      STORE(2,1,1);
      STORE(0,2,1);
      STORE(1,2,1);
      STORE(2,2,1);
      STORE(0,0,2);
      STORE(1,0,2);
      STORE(2,0,2);
      STORE(0,1,2);
      STORE(1,1,2);
      STORE(2,1,2);
      STORE(0,2,2);
      STORE(1,2,2);
      STORE(2,2,2);
#undef STORE
    }
#else
#pragma omp for
    for( int p=0; p<npart; p++){
      float pos[3];
      getPos2( &particle[p], pos);
      pos[0] -= g_pos_f[0];
      pos[1] -= g_pos_f[1];
      pos[2] -= g_pos_f[2];
      float mass = getMass( &particle[p], this_run->uniform_mass);
      
      float wi[3][3];
      int iw[3][3];

      for( int j=0; j<3; j++){
	float xt1  = pos[j] * (float)SIZE_MESH;
	iw[1][j] = (int)(xt1 + 0.5);
	float dx1 = xt1 - (float)(iw[1][j]);
	wi[0][j] = 0.5 * (0.5-dx1) * (0.5-dx1);
	wi[1][j] = 0.75 - dx1*dx1;
	wi[2][j] = 0.5 * (0.5 + dx1) * (0.5 + dx1);
	iw[0][j] = iw[1][j] - 1;
	iw[2][j] = iw[1][j] + 1;
      }

     for( int i=0; i<3; i++){
	for( int j=0; j<3; j++){
	  int iwxy = l_msize[2]*( iw[j][1] + l_msize[1]*iw[i][0]);
	  float wixy = mass*wi[i][0]*wi[j][1];
	  int k0 = iw[0][2] + iwxy;
	  int k1 = iw[1][2] + iwxy;
	  int k2 = iw[2][2] + iwxy;
	  mesh_density_local_thread[devid][k0] += wixy*wi[0][2];
	  mesh_density_local_thread[devid][k1] += wixy*wi[1][2];
	  mesh_density_local_thread[devid][k2] += wixy*wi[2][2];
	}
      }
    }
#endif
  }


  for( int i=0; i<NUMBER_OF_OMP_THREADS; i++){
#pragma omp parallel for
    for( int j=0; j<ln_total; j++){
      mesh_density_local[j] += mesh_density_local_thread[i][j];
    }
  }

  delete mesh_density_local_chunk;

#else
  for( int p=0; p<npart; p++){
    float pos[3];
    getPos2( &particle[p], pos);
    pos[0] -= g_pos[0];
    pos[1] -= g_pos[1];
    pos[2] -= g_pos[2];
    float mass = getMass( &particle[p], this_run->uniform_mass);

    float wi[3][3];
    int iw[3][3];

    for( int j=0; j<3; j++){
      float xt1  = pos[j] * (float)SIZE_MESH;
      iw[1][j] = (int)(xt1 + 0.5);
      float dx1 = xt1 - (float)(iw[1][j]);
      wi[0][j] = 0.5 * (0.5-dx1) * (0.5-dx1);
      wi[1][j] = 0.75 - dx1*dx1;
      wi[2][j] = 0.5 * (0.5 + dx1) * (0.5 + dx1);
      iw[0][j] = iw[1][j] - 1;
      iw[2][j] = iw[1][j] + 1;
    }

    for( int i=0; i<3; i++){
      for( int j=0; j<3; j++){
	for( int k=0; k<3; k++){
	  int ii = iw[k][2] + l_msize[2]*( iw[j][1] + l_msize[1]*iw[i][0]);
	  assert( ii >= 0);
	  assert( ii < ln_total);
	  mesh_density_local[ii] += mass*wi[i][0]*wi[j][1]*wi[k][2];
	}
      }
    }
  }
#endif

}



#ifdef RMM_PM
void PMForce::commMeshDensity(){

  double nowtime;
  getTime(&nowtime);

  static int g_offset_all[MAXNNODE][3];
  static int l_msize_all[MAXNNODE][3];

  float *sendbuf = new float[ln_total];
  static int nsend[MAXNNODE];
  static int nrecv[MAXNNODE];
  static int sdispls_sendbuf[MAXNNODE];
  static int rdispls_sendbuf[MAXNNODE];
  static int nsend_slab[MAXNNODE];
  static int nrecv_slab[MAXNNODE];
  static int sdispls_gix[MAXNNODE];
  static int rdispls_gix[MAXNNODE];
  int *gix_send = new int[SIZE_MESH];
  int *dest = new int[SIZE_MESH];
  for( int i=0; i<nnode_reduce; i++){
    nsend[i] = nrecv[i] = sdispls_sendbuf[i] = rdispls_sendbuf[i] = 0;
    nsend_slab[i] = nrecv_slab[i] =  sdispls_gix[i] = rdispls_gix[i] = 0;
  }
  for( int i=0; i<SIZE_MESH; i++){
    gix_send[i] = dest[i] = 0;
  }

  for( int i=0; i<l_msize[0]; i++){
    int gix = g_offset[0] + i;
    if( gix < 0)  gix += SIZE_MESH;
    if( gix >= SIZE_MESH)  gix -= SIZE_MESH;
    for( int j=0; j<nnode_reduce; j++){
      if( l_nx_all2[j] == 0)  continue;
      if( l_x_start_all2[j] <= gix){
	dest[i] = j;
      }
    }
    nsend[dest[i]] += lnyz;
    nsend_slab[dest[i]] ++;
  }

  for( int i=1; i<nnode_reduce; i++){
    sdispls_gix[i] = sdispls_gix[i-1] + nsend_slab[i-1];
    sdispls_sendbuf[i] = sdispls_sendbuf[i-1] + nsend[i-1];
    nsend[i-1] = 0;
    nsend_slab[i-1] = 0;
  }
  nsend[nnode_reduce-1] = 0;
  nsend_slab[nnode_reduce-1] = 0;
  for( int i=0; i<l_msize[0]; i++){
    int l_offset = lnyz * i;
    int buf_offset = sdispls_sendbuf[dest[i]] + nsend[dest[i]];
    for( int j=0; j<lnyz; j++){
      sendbuf[buf_offset+j] = mesh_density_local[l_offset+j];
    }
    nsend[dest[i]] += lnyz;
    int ii = sdispls_gix[dest[i]] + nsend_slab[dest[i]];
    gix_send[ii] = g_offset[0] + i;
    nsend_slab[dest[i]] ++;
  }
  My_MPI_Barrier( MPI_COMM_COLUMN);
  //fprintf_verbose(stderr, "#start comm\n");
  this_run->t.pm_make_comm_density += getTime(&nowtime);

  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_COLUMN);
  MPI_Alltoall( nsend_slab, 1, MPI_INT, nrecv_slab, 1, MPI_INT, MPI_COMM_COLUMN);
  int nrecv_total = nrecv[0];
  int gix_total = nrecv_slab[0];
  for( int i=1; i<nnode_reduce; i++){
    rdispls_sendbuf[i] = rdispls_sendbuf[i-1] + nrecv[i-1];
    rdispls_gix[i] = rdispls_gix[i-1] + nrecv_slab[i-1];
    nrecv_total += nrecv[i];
    gix_total += nrecv_slab[i];
  }
  float *recvbuf = new float[nrecv_total];
  int *gix_recv = new int[gix_total];

  int nscomm_pm_l2s = 0;
  int nrcomm_pm_l2s = 0;
  int nsend_pm_l2s = 0;
  int nrecv_pm_l2s = 0;
  for( int i=0; i<nnode_reduce; i++){
    if( nsend[i] != 0)  nscomm_pm_l2s ++;
    if( nrecv[i] != 0)  nrcomm_pm_l2s ++;
    nsend_pm_l2s += nsend[i];
    nrecv_pm_l2s += nrecv[i];
  }
  this_run->nscomm_pm_l2s = nscomm_pm_l2s;
  this_run->nrcomm_pm_l2s = nrcomm_pm_l2s;
  this_run->nsend_pm_l2s = nsend_pm_l2s;
  this_run->nrecv_pm_l2s = nrecv_pm_l2s;

  MPI_Alltoallv( sendbuf, nsend, sdispls_sendbuf, MPI_FLOAT,
                 recvbuf, nrecv, rdispls_sendbuf, MPI_FLOAT,
                 MPI_COMM_COLUMN);
  MPI_Alltoallv( gix_send, nsend_slab, sdispls_gix, MPI_INT,
                 gix_recv, nrecv_slab, rdispls_gix, MPI_INT,
                 MPI_COMM_COLUMN);
  MPI_Allgather( g_offset, 3, MPI_INT, 
		 g_offset_all, 3, MPI_INT, MPI_COMM_COLUMN);
  MPI_Allgather( l_msize, 3, MPI_INT, 
		 l_msize_all, 3, MPI_INT, MPI_COMM_COLUMN);
  My_MPI_Barrier( MPI_COMM_COLUMN);
  this_run->t.pm_density_comm += getTime(&nowtime);
  //fprintf_verbose( stderr, "#end comm\n");

  for( int p=0; p<nnode_reduce; p++){
    if( nrecv[p] == 0)  continue;
    int offset_this = rdispls_sendbuf[p];
    for( int i=0; i<nrecv_slab[p]; i++){
      int ix0 = gix_recv[rdispls_gix[p]+i];
      if( ix0 < 0)  ix0 += SIZE_MESH;
      if( ix0>= SIZE_MESH)  ix0 -= SIZE_MESH;
      int ix = ix0 - l_x_start_all2[inode_reduce];
#pragma omp parallel for
      for( int j=0; j<l_msize_all[p][1]; j++){
	int iy = g_offset_all[p][1] + j;
	if( iy < 0)  iy += SIZE_MESH;
	if( iy >= SIZE_MESH)  iy -= SIZE_MESH;
	for( int k=0; k<l_msize_all[p][2]; k++){
	  int iz = g_offset_all[p][2] + k;
	  if( iz < 0)  iz += SIZE_MESH;
	  if( iz >= SIZE_MESH)  iz -= SIZE_MESH;
	  long long int ii = k + l_msize_all[p][2]*(j + i*l_msize_all[p][1]);
	  long long int li = offset_this + ii;
	  long long int gi = iz + SIZE_MESH_P2*(iy + SIZE_MESH*ix);
	  assert( gi >= 0);
	  assert( gi < msize_slab);
	  mesh_density_slab0[gi] += recvbuf[li];
	}
      }
    }
  }

  local_to_slab_recvbuf_size = nrecv_total;
  delete [] sendbuf;
  delete [] recvbuf;
  delete [] gix_recv;
  delete [] gix_send;
  delete [] dest;
  this_run->t.pm_make_comm_density += getTime(&nowtime);

  fprintf_verbose( stderr, "#start reduce\n");
  My_MPI_Barrier( MPI_COMM_COLUMN);
  MPI_Reduce( mesh_density_slab0, mesh_density_slab, msize_slab, MPI_TYPE_FFTW,
	      MPI_SUM, 0, MPI_COMM_REDUCE);
  My_MPI_Barrier( MPI_COMM_COLUMN);
  fprintf_verbose( stderr, "#end reduce\n");
  this_run->t.pm_density_comm_reduce += getTime(&nowtime);

}
#else
void PMForce::commMeshDensity(){

  double nowtime;
  getTime(&nowtime);

  static int g_offset_all[MAXNNODE][3];
  static int l_msize_all[MAXNNODE][3];

  float *sendbuf = new float[ln_total];
  static int nsend[MAXNNODE];
  static int nrecv[MAXNNODE];
  static int sdispls_sendbuf[MAXNNODE];
  static int rdispls_sendbuf[MAXNNODE];
  static int nsend_slab[MAXNNODE];
  static int nrecv_slab[MAXNNODE];
  static int sdispls_gix[MAXNNODE];
  static int rdispls_gix[MAXNNODE];

  int *gix_send = new int[SIZE_MESH];
#if 0
  // original 
  int *dest = new int[SIZE_MESH];
  for( int i=0; i<SIZE_MESH; i++){
    gix_send[i] = dest[i] = 0;
  }
#else
  // modified by M.I. 2016 6/7
  const int size_dest = static_cast<int>(l_msize[0]*1.2)+1024;
  int *dest = new int[size_dest];
  for( int i=0; i<SIZE_MESH; i++) gix_send[i] = 0;
  for( int i=0; i<size_dest; i++) dest[i] = 0;
#endif
  for( int i=0; i<this_run->nnode; i++){
    nsend[i] = nrecv[i] = sdispls_sendbuf[i] = rdispls_sendbuf[i] = 0;
    nsend_slab[i] = nrecv_slab[i] =  sdispls_gix[i] = rdispls_gix[i] = 0;
  }  
  for( int i=0; i<l_msize[0]; i++){
    int gix = g_offset[0] + i;
    if( gix < 0)  gix += SIZE_MESH;
    if( gix >= SIZE_MESH)  gix -= SIZE_MESH;
    for( int j=0; j<this_run->nnode; j++){
#ifdef FIX_FFTNODE
      if( l_nx_all[j] == 0)  continue;
#else
      if( l_nx_all[j] == 0)  break;
#endif
      if( l_x_start_all[j] <= gix){
	dest[i] = j;
      }
    }
    nsend[dest[i]] += lnyz;
    nsend_slab[dest[i]] ++;
    /*
    if( this_run->inode == 0){
      cerr << i << "/" << l_msize[0] << "\t" << dest[i] << "\t" 
	   << nsend[dest[i]] << "\t"
	   << gix << "\t" << nsend_slab[dest[i]] << endl;
    }
    */
  }

  for( int i=1; i<this_run->nnode; i++){
    sdispls_gix[i] = sdispls_gix[i-1] + nsend_slab[i-1];
    sdispls_sendbuf[i] = sdispls_sendbuf[i-1] + nsend[i-1];
    nsend[i-1] = 0;
    nsend_slab[i-1] = 0;
  }
  nsend[this_run->nnode-1] = 0;
  nsend_slab[this_run->nnode-1] = 0;
  for( int i=0; i<l_msize[0]; i++){
    int l_offset = lnyz * i;
    int buf_offset = sdispls_sendbuf[dest[i]] + nsend[dest[i]];
    for( int j=0; j<lnyz; j++){
      sendbuf[buf_offset+j] = mesh_density_local[l_offset+j];
    }
    nsend[dest[i]] += lnyz;
    int ii = sdispls_gix[dest[i]] + nsend_slab[dest[i]];
    gix_send[ii] = g_offset[0] + i;
    nsend_slab[dest[i]] ++;
  }
  My_MPI_Barrier( this_run->MPI_COMM_INTERNAL);
  this_run->t.pm_make_comm_density += getTime(&nowtime);

  //fprintf_verbose( stderr, "#start comm\n");
  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);
  MPI_Alltoall( nsend_slab, 1, MPI_INT, nrecv_slab, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);
  int nrecv_total = nrecv[0];
  int gix_total = nrecv_slab[0];
  for( int i=1; i<this_run->nnode; i++){
    rdispls_sendbuf[i] = rdispls_sendbuf[i-1] + nrecv[i-1];
    rdispls_gix[i] = rdispls_gix[i-1] + nrecv_slab[i-1];
    nrecv_total += nrecv[i];
    gix_total += nrecv_slab[i];
  }
  float *recvbuf = new float[nrecv_total];
  int *gix_recv = new int[gix_total];

  int nscomm_pm_l2s = 0;
  int nrcomm_pm_l2s = 0;
  int nsend_pm_l2s = 0;
  int nrecv_pm_l2s = 0;
  for( int i=0; i<this_run->nnode; i++){
    if( nsend[i] != 0)  nscomm_pm_l2s ++;
    if( nrecv[i] != 0)  nrcomm_pm_l2s ++;
    nsend_pm_l2s += nsend[i];
    nrecv_pm_l2s += nrecv[i];
  }
  this_run->nscomm_pm_l2s = nscomm_pm_l2s;
  this_run->nrcomm_pm_l2s = nrcomm_pm_l2s;
  this_run->nsend_pm_l2s = nsend_pm_l2s;
  this_run->nrecv_pm_l2s = nrecv_pm_l2s;

  MPI_Alltoallv( sendbuf, nsend, sdispls_sendbuf, MPI_FLOAT,
                 recvbuf, nrecv, rdispls_sendbuf, MPI_FLOAT,
                 this_run->MPI_COMM_INTERNAL);
  MPI_Alltoallv( gix_send, nsend_slab, sdispls_gix, MPI_INT,
                 gix_recv, nrecv_slab, rdispls_gix, MPI_INT,
                 this_run->MPI_COMM_INTERNAL);
  MPI_Allgather( g_offset, 3, MPI_INT, 
		 g_offset_all, 3, MPI_INT, this_run->MPI_COMM_INTERNAL);
  MPI_Allgather( l_msize, 3, MPI_INT, 
		 l_msize_all, 3, MPI_INT, this_run->MPI_COMM_INTERNAL);
  My_MPI_Barrier( this_run->MPI_COMM_INTERNAL);
  //fprintf_verbose( stderr, "#end comm\n");
  this_run->t.pm_density_comm += getTime(&nowtime);

  for( int p=0; p<this_run->nnode; p++){
    if( nrecv[p] == 0)  continue;
    int offset_this = rdispls_sendbuf[p];
    for( int i=0; i<nrecv_slab[p]; i++){
      int ix0 = gix_recv[rdispls_gix[p]+i];
      if( ix0 < 0)  ix0 += SIZE_MESH;
      if( ix0>= SIZE_MESH)  ix0 -= SIZE_MESH;
      int ix = ix0 - l_x_start_all[this_run->inode];
#pragma omp parallel for
      for( int j=0; j<l_msize_all[p][1]; j++){
	int iy = g_offset_all[p][1] + j;
	if( iy < 0)  iy += SIZE_MESH;
	if( iy >= SIZE_MESH)  iy -= SIZE_MESH;
	for( int k=0; k<l_msize_all[p][2]; k++){
	  int iz = g_offset_all[p][2] + k;
	  if( iz < 0)  iz += SIZE_MESH;
	  if( iz >= SIZE_MESH)  iz -= SIZE_MESH;
	  long long int ii = k + l_msize_all[p][2]*(j + i*l_msize_all[p][1]);
	  long long int li = offset_this + ii;
	  long long int gi = iz + SIZE_MESH_P2*(iy + SIZE_MESH*ix);
	  assert( gi >= 0);
	  assert( gi < msize_slab);
	  mesh_density_slab[gi] += recvbuf[li];

	}
      }
    }
  }

  local_to_slab_recvbuf_size = nrecv_total;
  delete [] sendbuf;
  delete [] recvbuf;
  delete [] gix_recv;
  delete [] gix_send;
  delete [] dest;

  this_run->t.pm_make_comm_density += getTime(&nowtime);

}
#endif



void PMForce::commMeshDensity3D(){
#ifdef FFT3D

  double nowtime;
  getTime(&nowtime);

  static int nsend[MAXNNODE];
  static int nrecv[MAXNNODE];
  static int sdispls_sendbuf[MAXNNODE];
  static int rdispls_sendbuf[MAXNNODE];
  int *sendbuf_subbox_dest = new int[ln_total];
  int *dest = new int[ln_total];
  for( int i=0; i<nnode_reduce; i++)  nsend[i] = 0;

#pragma omp parallel for
  for( int i=0; i<l_msize[0]; i++){
    int gix = g_offset[0] + i;
    if( gix < 0)  gix += SIZE_MESH;
    if( gix >= SIZE_MESH)  gix -= SIZE_MESH;
    int xnode = gix / subbox_msize[0];
    int lix = gix - xnode*subbox_msize[0];

    for( int j=0; j<l_msize[1]; j++){
      int giy = g_offset[1] + j;
      if( giy < 0)  giy += SIZE_MESH;
      if( giy >= SIZE_MESH)  giy -= SIZE_MESH;
      int ynode = giy / subbox_msize[1];
      int liy = giy - ynode*subbox_msize[1];

      for( int k=0; k<l_msize[2]; k++){
	int giz = g_offset[2] + k;
	if( giz < 0)  giz += SIZE_MESH;
	if( giz >= SIZE_MESH)  giz -= SIZE_MESH;
	int znode = giz / subbox_msize[2];
	int liz = giz - znode*subbox_msize[2];

	int im= getLocalID( i, j, k);
	sendbuf_subbox_dest[im] = getSubboxID( lix, liy, liz);
	dest[im] = getVoxelIndex( xnode, ynode, znode, ndiv_reduce);
      }
    }
  }
  qsortPivotWrapper( mesh_density_local, sendbuf_subbox_dest, dest, 
		     ln_total, NLEVEL_KEY_SORT, NCHUNK_QSORT);
  for( int i=0; i<ln_total; i++)  nsend[dest[i]] ++;

  My_MPI_Barrier( MPI_COMM_COLUMN);
  //fprintf_verbose( stderr, "#start comm\n");
  this_run->t.pm_make_comm_density += getTime(&nowtime);

  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_COLUMN);
  int nrecv_total = nrecv[0];
  sdispls_sendbuf[0] = rdispls_sendbuf[0] = 0;
  for( int i=1; i<nnode_reduce; i++){
    sdispls_sendbuf[i] = sdispls_sendbuf[i-1] + nsend[i-1];
    rdispls_sendbuf[i] = rdispls_sendbuf[i-1] + nrecv[i-1];
    nrecv_total += nrecv[i];
  }
  local_to_slab_recvbuf_size = nrecv_total;

  float *recvbuf = new float[nrecv_total];
  int *recvbuf_subbox_dest = new int[nrecv_total];

  int nscomm_pm_l2s = 0;
  int nrcomm_pm_l2s = 0;
  int nsend_pm_l2s = 0;
  int nrecv_pm_l2s = 0;
  for( int i=0; i<nnode_reduce; i++){
    if( nsend[i] != 0)  nscomm_pm_l2s ++;
    if( nrecv[i] != 0)  nrcomm_pm_l2s ++;
    nsend_pm_l2s += nsend[i];
    nrecv_pm_l2s += nrecv[i];
  }
  this_run->nscomm_pm_l2s = nscomm_pm_l2s;
  this_run->nrcomm_pm_l2s = nrcomm_pm_l2s;
  this_run->nsend_pm_l2s = nsend_pm_l2s;
  this_run->nrecv_pm_l2s = nrecv_pm_l2s;

  MPI_Alltoallv( mesh_density_local, nsend, sdispls_sendbuf, MPI_FLOAT,
                 recvbuf, nrecv, rdispls_sendbuf, MPI_FLOAT,
                 MPI_COMM_COLUMN);
  MPI_Alltoallv( sendbuf_subbox_dest, nsend, sdispls_sendbuf, MPI_INT,
                 recvbuf_subbox_dest, nrecv, rdispls_sendbuf, MPI_INT,
                 MPI_COMM_COLUMN);
  My_MPI_Barrier( MPI_COMM_COLUMN);
  this_run->t.pm_density_comm += getTime(&nowtime);
  //fprintf_verbose( stderr, "#end comm\n");

  for( int i=0; i<nrecv_total; i++){
    int im = recvbuf_subbox_dest[i];
    mesh_density_slab0[im] += recvbuf[i];
  }

  delete [] sendbuf_subbox_dest;
  delete [] dest;
  delete [] recvbuf;
  delete [] recvbuf_subbox_dest;
  this_run->t.pm_make_comm_density += getTime(&nowtime);

  fprintf_verbose( stderr, "#start reduce\n");
  My_MPI_Barrier( MPI_COMM_COLUMN);
  MPI_Reduce( mesh_density_slab0, mesh_density_slab, msize_slab, MPI_TYPE_FFTW,
	      MPI_SUM, 0, MPI_COMM_REDUCE);
  My_MPI_Barrier( MPI_COMM_COLUMN);
  fprintf_verbose( stderr, "#end reduce\n");
  this_run->t.pm_density_comm_reduce += getTime(&nowtime);

#endif
}



void PMForce::fftForward(){

  if( fft_node != 1)  return;

  if( fft3d == 1){
#pragma omp parallel for
    for( int i=0; i<subbox_total; i++){
      mesh_density_subbox[i][0] = mesh_density_slab[i];
      mesh_density_subbox[i][1] = 0.0;
    }

#ifdef KCOMPUTER
    int nm = SIZE_MESH;
    int isn = 1;
    int icon = 0;
    int nnode0 = 0;
    MPI_Comm_size(MPI_COMM_FFT, &nnode0);
    static int first_call = 0;
    if( first_call == 0){
      sprintf( cbuf, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	       nnode0, subbox_msize[0], subbox_msize[1], subbox_msize[2],
	       nm, ndiv_fft[0], ndiv_fft[1], ndiv_fft[2], subbox_total);
      //mp_print( cbuf, MPI_COMM_FFT, stderr);
      first_call = 1;
    }

    ds_v3dcft3_( (double *)mesh_density_subbox,
		 subbox_msize, subbox_msize+1,
		 subbox_msize, subbox_msize+1, subbox_msize+2,
		 &nm, &nm, &nm,
		 &ndiv_fft[0], &ndiv_fft[1], &ndiv_fft[2],
		 (double *)fft3d_work, &subbox_total, &isn, &MPI_COMM_FFT_F, &icon);

    if( inode_fft == 0){
      if( icon != 0){
	if( this_run->inode == 0)  cerr << "fftForward icon=" << icon << endl;
      }
    }
#else
    if( this_run->inode == 0)  cerr << "fft3d is not available" << endl;
#endif
  }

  else{

#ifdef FFTW3_PARALLEL
    fftw_execute(fplan);
#else
    rfftwnd_mpi( fplan, 1, mesh_density_slab, fftwork, FFTW_TRANSPOSED_ORDER);
#endif
  }

}



void PMForce::calcMeshPhi(){


  if( fft3d == 1){
#pragma omp parallel for
    for( int i=0; i<subbox_total; i++){
      mesh_density_subbox[i][0] *= gk[i];
      mesh_density_subbox[i][1] *= gk[i];
    }
  }

  else{

    static int nz = SIZE_MESH/2 + 1;
    static int nz2 = SIZE_MESH/2;

    fftw_complex *mesh = (fftw_complex*) mesh_density_slab;

#pragma omp parallel for
    for( int i=0; i<local_ny_after_transpose; i++){
      for( int j=0; j<nz2; j++){
	for( int k=0; k<nz; k++){
	  int mi = k + nz*( j + SIZE_MESH*i);
	  int gi = k + nz*( j + nz*i);
#ifdef FFTW3_PARALLEL
	  mesh[mi][0] *= gk[gi];
	  mesh[mi][1] *= gk[gi];
#else
	  mesh[mi].re *= gk[gi];
	  mesh[mi].im *= gk[gi];
#endif
	  mi = k + nz*( nz2+j + SIZE_MESH*i);
	  gi = k + nz*( nz2-j + nz*i);
#ifdef FFTW3_PARALLEL
	  mesh[mi][0] *= gk[gi];
	  mesh[mi][1] *= gk[gi];
#else
	  mesh[mi].re *= gk[gi];
	  mesh[mi].im *= gk[gi];
#endif
	}
      }
    }
  }

}



void PMForce::fftBackward(){

  if( fft_node != 1)  return;

  if( fft3d == 1){
#ifdef KCOMPUTER
    int nm = SIZE_MESH;
    int isn = -1;
    int icon = 0;
    ds_v3dcft3_( (double *)mesh_density_subbox,
		 subbox_msize, subbox_msize+1,
		 subbox_msize, subbox_msize+1, subbox_msize+2,
		 &nm, &nm, &nm,
		 ndiv_fft, ndiv_fft+1, ndiv_fft+2,
		 (double *)fft3d_work, &subbox_total, &isn, &MPI_COMM_FFT_F, &icon);
    if( inode_fft == 0){
      if( icon != 0){
	cerr << "fftBackward icon=" << icon << endl;
      }
    }
#endif

#pragma omp parallel for
    for( int i=0; i<subbox_total; i++){
      mesh_density_slab[i] = mesh_density_subbox[i][0];
    }

  }
  else{

#ifdef FFTW3_PARALLEL
    fftw_execute(bplan);
#else
    rfftwnd_mpi( bplan, 1, mesh_density_slab, fftwork, FFTW_TRANSPOSED_ORDER);
#endif
  }


}



#ifdef RMM_PM
void PMForce::commMeshPhi(){

  double nowtime;
  fprintf_verbose( stderr, "#start bcast\n");
  My_MPI_Barrier( MPI_COMM_REDUCE);
  getTime(&nowtime);
  MPI_Bcast( mesh_density_slab, msize_slab, MPI_TYPE_FFTW,
	     0, MPI_COMM_REDUCE);
  My_MPI_Barrier( MPI_COMM_REDUCE);
  this_run->t.pm_phi_comm_bcast += getTime(&nowtime);
  fprintf_verbose( stderr, "#end bcast\n");

  getLocalMeshParam(3);
  static int g_offset_all[MAXNNODE][3];
  static int l_msize_all[MAXNNODE][3];
  MPI_Allgather( g_offset, 3, MPI_INT, 
		 g_offset_all, 3, MPI_INT, MPI_COMM_COLUMN);
  MPI_Allgather( l_msize, 3, MPI_INT, 
		 l_msize_all, 3, MPI_INT, MPI_COMM_COLUMN);

  float *recvbuf = new float[ln_total];
  int *gix_recv = new int[SIZE_MESH];
  static int nsend[MAXNNODE];
  static int nrecv[MAXNNODE];
  static int sdispls_sendbuf[MAXNNODE];
  static int rdispls_sendbuf[MAXNNODE];
  static int nsend_slab[MAXNNODE];
  static int nrecv_slab[MAXNNODE];
  static int sdispls_gix[MAXNNODE];
  static int rdispls_gix[MAXNNODE];
  for( int i=0; i<nnode_reduce; i++){
    nsend[i] = nrecv[i] = sdispls_sendbuf[i] = rdispls_sendbuf[i] = 0;
    nsend_slab[i] = nrecv_slab[i] =  sdispls_gix[i] = rdispls_gix[i] = 0;
  }
  for( int i=0; i<SIZE_MESH; i++){
    gix_recv[i] = 0;
  }

  int nsend_total = 0;
  int nsend_slab_total = 0;
  for( int p=0; p<nnode_reduce; p++){
    int p_xstart = g_offset_all[p][0];
    for( int i=0; i<l_msize_all[p][0]; i++){
      int p_ix = p_xstart + i;
      if( p_ix < 0)  p_ix += SIZE_MESH;
      if( p_ix >= SIZE_MESH)  p_ix -= SIZE_MESH;
      if( p_ix < local_x_start2)  continue;
      if( p_ix >= local_x_start2+local_nx2)  continue;
      nsend_slab[p] ++;
      nsend[p] += l_msize_all[p][1]*l_msize_all[p][2];
    }
    nsend_total += nsend[p];
    nsend_slab_total += nsend_slab[p];
  }

  for( int i=1; i<this_run->nnode; i++){
    sdispls_sendbuf[i] = sdispls_sendbuf[i-1] + nsend[i-1];
    sdispls_gix[i] = sdispls_gix[i-1] + nsend_slab[i-1];
    nsend[i-1] = 0;
    nsend_slab[i-1] = 0;
  }
  nsend[nnode_reduce-1] = 0;
  nsend_slab[nnode_reduce-1] = 0;

  float *sendbuf = new float[nsend_total];
  int *gix_send = new int[nsend_slab_total];

#pragma omp parallel for CHUNK_PM
  for( int p=0; p<nnode_reduce; p++){
    int l_offset = sdispls_gix[p];
    int buf_offset = sdispls_sendbuf[p];
    int p_xstart = g_offset_all[p][0];

    int ii = 0;
    for( int i=0; i<l_msize_all[p][0]; i++){
      int p_ix0 = p_xstart + i;
      int p_ix = p_ix0;
      if( p_ix < 0)  p_ix += SIZE_MESH;
      if( p_ix >= SIZE_MESH)  p_ix -= SIZE_MESH;
      if( p_ix < local_x_start2)  continue;
      if( p_ix >= local_x_start2+local_nx2)  continue;
      p_ix -= local_x_start2;
      gix_send[l_offset + nsend_slab[p]] = p_ix0;
      nsend_slab[p] ++;
      nsend[p] += l_msize_all[p][1]*l_msize_all[p][2];

      int p_ystart = g_offset_all[p][1];
      int p_zstart = g_offset_all[p][2];
      for( int j=0; j<l_msize_all[p][1]; j++){
	int p_iy = p_ystart + j;
	if( p_iy < 0)  p_iy += SIZE_MESH;
	if( p_iy >= SIZE_MESH)  p_iy -= SIZE_MESH;
	for( int k=0; k<l_msize_all[p][2]; k++){
	  int p_iz = p_zstart + k;
	  if( p_iz < 0)  p_iz += SIZE_MESH;
	  if( p_iz >= SIZE_MESH)  p_iz -= SIZE_MESH;
	  int li = p_iz + SIZE_MESH_P2*( p_iy + SIZE_MESH*p_ix);
	  sendbuf[buf_offset+ii] = mesh_density_slab[li];
	  ii ++;
	}
      }
    }
  }
  My_MPI_Barrier( MPI_COMM_COLUMN);
  this_run->t.pm_make_comm_phi += getTime(&nowtime);

  //fprintf_verbose( stderr, "#start comm0\n");
  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_COLUMN);
  MPI_Alltoall( nsend_slab, 1, MPI_INT, nrecv_slab, 1, MPI_INT, MPI_COMM_COLUMN);
  //fprintf_verbose( stderr, "#end comm0\n");

  int nrecv_total = nrecv[0];
  int nrecv_slab_total = nrecv_slab[0];
  for( int i=1; i<nnode_reduce; i++){
    rdispls_sendbuf[i] = rdispls_sendbuf[i-1] + nrecv[i-1];
    rdispls_gix[i] = rdispls_gix[i-1] + nrecv_slab[i-1];
    nrecv_total += nrecv[i];
    nrecv_slab_total += nrecv_slab[i];
  }

  int nscomm_pm_s2l = 0;
  int nrcomm_pm_s2l = 0;
  int nsend_pm_s2l = 0;
  int nrecv_pm_s2l = 0;
  for( int i=0; i<nnode_reduce; i++){
    if( nsend[i] != 0)  nscomm_pm_s2l ++;
    if( nrecv[i] != 0)  nrcomm_pm_s2l ++;
    nsend_pm_s2l += nsend[i];
    nrecv_pm_s2l += nrecv[i];
  }
  this_run->nscomm_pm_s2l = nscomm_pm_s2l;
  this_run->nrcomm_pm_s2l = nrcomm_pm_s2l;
  this_run->nsend_pm_s2l = nsend_pm_s2l;
  this_run->nrecv_pm_s2l = nrecv_pm_s2l;

  //fprintf_verbose( stderr, "#start comm\n");
  MPI_Alltoallv( sendbuf, nsend, sdispls_sendbuf, MPI_FLOAT,
                 recvbuf, nrecv, rdispls_sendbuf, MPI_FLOAT,
                 MPI_COMM_COLUMN);
  MPI_Alltoallv( gix_send, nsend_slab, sdispls_gix, MPI_INT,
                 gix_recv, nrecv_slab, rdispls_gix, MPI_INT,
                 MPI_COMM_COLUMN);
  My_MPI_Barrier( MPI_COMM_COLUMN);
  this_run->t.pm_phi_comm += getTime(&nowtime);
  //fprintf_verbose( stderr, "#end comm\n");

#pragma omp parallel for CHUNK_PM
  for( int p=0; p<nnode_reduce; p++){
    for( int i=0; i<nrecv_slab[p]; i++){
      int p_ix = gix_recv[rdispls_gix[p]+i] - g_offset[0];
      int buf_offset = rdispls_sendbuf[p] + i*lnyz;
      for( int j=0; j<lnyz; j++){
	int gi = p_ix*lnyz + j;
	int li = buf_offset + j;
	mesh_density_local[gi] = recvbuf[li];
      }
    }
  }

  slab_to_local_sendbuf_size = nsend_total;
  delete [] recvbuf;
  delete [] sendbuf;
  delete [] gix_send;
  delete [] gix_recv;

  this_run->t.pm_make_comm_phi += getTime(&nowtime);

}

#else
void PMForce::commMeshPhi(){

  double nowtime;
  getTime(&nowtime);

  getLocalMeshParam(3);
  static int g_offset_all[MAXNNODE][3];
  static int l_msize_all[MAXNNODE][3];
  MPI_Allgather( g_offset, 3, MPI_INT, 
		 g_offset_all, 3, MPI_INT, this_run->MPI_COMM_INTERNAL);
  MPI_Allgather( l_msize, 3, MPI_INT, 
		 l_msize_all, 3, MPI_INT, this_run->MPI_COMM_INTERNAL);

  float *recvbuf = new float[ln_total];
  //int *gix_recv = new int[SIZE_MESH]; // comment out M.I. 2016/02/15
  static int nsend[MAXNNODE];
  static int nrecv[MAXNNODE];
  static int sdispls_sendbuf[MAXNNODE];
  static int rdispls_sendbuf[MAXNNODE];
  static int nsend_slab[MAXNNODE];
  static int nrecv_slab[MAXNNODE];
  static int sdispls_gix[MAXNNODE];
  static int rdispls_gix[MAXNNODE];
  for( int i=0; i<this_run->nnode; i++){
    nsend[i] = nrecv[i] = sdispls_sendbuf[i] = rdispls_sendbuf[i] = 0;
    nsend_slab[i] = nrecv_slab[i] =  sdispls_gix[i] = rdispls_gix[i] = 0;
  }
  /*
    comment out M.I. 2016/02/15
  for( int i=0; i<SIZE_MESH; i++){
    gix_recv[i] = 0;
  }
  */
  int nsend_total = 0;
  int nsend_slab_total = 0;
  for( int p=0; p<this_run->nnode; p++){
    int p_xstart = g_offset_all[p][0];
    for( int i=0; i<l_msize_all[p][0]; i++){
      int p_ix = p_xstart + i;
      if( p_ix < 0)  p_ix += SIZE_MESH;
      if( p_ix >= SIZE_MESH)  p_ix -= SIZE_MESH;
      if( p_ix < local_x_start)  continue;
      if( p_ix >= local_x_start+local_nx)  continue;
      nsend_slab[p] ++;
      nsend[p] += l_msize_all[p][1]*l_msize_all[p][2];
    }
    nsend_total += nsend[p];
    nsend_slab_total += nsend_slab[p];
  }

  for( int i=1; i<this_run->nnode; i++){
    sdispls_sendbuf[i] = sdispls_sendbuf[i-1] + nsend[i-1];
    sdispls_gix[i] = sdispls_gix[i-1] + nsend_slab[i-1];
    nsend[i-1] = 0;
    nsend_slab[i-1] = 0;
  }
  nsend[this_run->nnode-1] = 0;
  nsend_slab[this_run->nnode-1] = 0;

  float *sendbuf = new float[nsend_total];
  int *gix_send = new int[nsend_slab_total];

#pragma omp parallel for CHUNK_PM
  for( int p=0; p<this_run->nnode; p++){
    int l_offset = sdispls_gix[p];
    int buf_offset = sdispls_sendbuf[p];
    int p_xstart = g_offset_all[p][0];

    int ii = 0;
    for( int i=0; i<l_msize_all[p][0]; i++){
      int p_ix0 = p_xstart + i;
      int p_ix = p_ix0;
      if( p_ix < 0)  p_ix += SIZE_MESH;
      if( p_ix >= SIZE_MESH)  p_ix -= SIZE_MESH;
      if( p_ix < local_x_start)  continue;
      if( p_ix >= local_x_start+local_nx)  continue;
      p_ix -= local_x_start;
      gix_send[l_offset + nsend_slab[p]] = p_ix0;
      nsend_slab[p] ++;
      nsend[p] += l_msize_all[p][1]*l_msize_all[p][2];

      int p_ystart = g_offset_all[p][1];
      int p_zstart = g_offset_all[p][2];
      for( int j=0; j<l_msize_all[p][1]; j++){
	int p_iy = p_ystart + j;
	if( p_iy < 0)  p_iy += SIZE_MESH;
	if( p_iy >= SIZE_MESH)  p_iy -= SIZE_MESH;
	for( int k=0; k<l_msize_all[p][2]; k++){
	  int p_iz = p_zstart + k;
	  if( p_iz < 0)  p_iz += SIZE_MESH;
	  if( p_iz >= SIZE_MESH)  p_iz -= SIZE_MESH;
	  int li = p_iz + SIZE_MESH_P2*( p_iy + SIZE_MESH*p_ix);
	  sendbuf[buf_offset+ii] = mesh_density_slab[li];
	  ii ++;
	}
      }
    }
  }
  My_MPI_Barrier( this_run->MPI_COMM_INTERNAL);
  this_run->t.pm_make_comm_phi += getTime(&nowtime);
  //fprintf_verbose( stderr, "#start comm\n");

  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);
  MPI_Alltoall( nsend_slab, 1, MPI_INT, nrecv_slab, 1, MPI_INT, this_run->MPI_COMM_INTERNAL);

  int nrecv_total = nrecv[0];
  int nrecv_slab_total = nrecv_slab[0];
  for( int i=1; i<this_run->nnode; i++){
    rdispls_sendbuf[i] = rdispls_sendbuf[i-1] + nrecv[i-1];
    rdispls_gix[i] = rdispls_gix[i-1] + nrecv_slab[i-1];
    nrecv_total += nrecv[i];
    nrecv_slab_total += nrecv_slab[i];
  }

  // add by M.I. 2016/02/15
  // bug fix by M.I. 2016/03/01
  int *gix_recv = new int[static_cast<int>(nrecv_slab_total*1.2)+1024];
  for(int i=0; i<nrecv_slab_total; i++){
      gix_recv[i] = 0;
  }  

  int nscomm_pm_s2l = 0;
  int nrcomm_pm_s2l = 0;
  int nsend_pm_s2l = 0;
  int nrecv_pm_s2l = 0;
  for( int i=0; i<this_run->nnode; i++){
    if( nsend[i] != 0)  nscomm_pm_s2l ++;
    if( nrecv[i] != 0)  nrcomm_pm_s2l ++;
    nsend_pm_s2l += nsend[i];
    nrecv_pm_s2l += nrecv[i];
  }
  this_run->nscomm_pm_s2l = nscomm_pm_s2l;
  this_run->nrcomm_pm_s2l = nrcomm_pm_s2l;
  this_run->nsend_pm_s2l = nsend_pm_s2l;
  this_run->nrecv_pm_s2l = nrecv_pm_s2l;

  MPI_Alltoallv( sendbuf, nsend, sdispls_sendbuf, MPI_FLOAT,
                 recvbuf, nrecv, rdispls_sendbuf, MPI_FLOAT,
                 this_run->MPI_COMM_INTERNAL);

  MPI_Alltoallv( gix_send, nsend_slab, sdispls_gix, MPI_INT,
                 gix_recv, nrecv_slab, rdispls_gix, MPI_INT,
                 this_run->MPI_COMM_INTERNAL);
  My_MPI_Barrier( this_run->MPI_COMM_INTERNAL);
  this_run->t.pm_phi_comm += getTime(&nowtime);
  //fprintf_verbose( stderr, "#end comm\n");

#pragma omp parallel for CHUNK_PM
  for( int p=0; p<this_run->nnode; p++){
    for( int i=0; i<nrecv_slab[p]; i++){
      int p_ix = gix_recv[rdispls_gix[p]+i] - g_offset[0];
      int buf_offset = rdispls_sendbuf[p] + i*lnyz;

      /*
      if( this_run->inode == 0){
	cerr << p_ix << "\t" << l_msize[0] << "\t" 
	     << buf_offset << "\t" << nrecv_total << "\t" << lnyz 
	     << "\t" << rdispls_gix[p] << "\t" << nrecv_slab[p]
	     << "\t" << p_ix*lnyz << "\t" << ln_total << endl;
      }
      */

      for( int j=0; j<lnyz; j++){
	int gi = p_ix*lnyz + j;
	int li = buf_offset + j;
	mesh_density_local[gi] = recvbuf[li];
      }
    }
  }

  slab_to_local_sendbuf_size = nsend_total;
  delete [] recvbuf;
  delete [] sendbuf;
  delete [] gix_send;
  delete [] gix_recv;

  this_run->t.pm_make_comm_phi += getTime(&nowtime);

}
#endif



int determineBoxOverlapping2( const int *subbox_start,
			      const int *subbox_msize, 
			      const int (*g_offset_all)[3],
			      const int (*l_msize_all)[3],
			      const int SIZE_MESH,
			      const int nnode_reduce,
			      int *inode_list){

  double dminv = 1.0 / (double)SIZE_MESH;

  int nnode_comm = 0;
  double cbox1[3], box_half_size1[3];
  for( int j=0; j<3; j++){
    cbox1[j] = dminv * ( (double)subbox_start[j] + 0.5*(double)subbox_msize[j]);
    box_half_size1[j] = 0.5 * dminv * (double)subbox_msize[j];
  }

#pragma omp parallel for
  for( int p=0; p<nnode_reduce; p++){
    double cbox2[3], box_half_size2[3];
    for( int j=0; j<3; j++){
      cbox2[j] = dminv * ( (double)g_offset_all[p][j] + 0.5*(double)l_msize_all[p][j]);
      box_half_size2[j] = 0.5 * dminv * (double)l_msize_all[p][j];
    }
    inode_list[p] = 0;
    if( determineBoxOverlapping( cbox1, box_half_size1, 
				 cbox2, box_half_size2) == 1){

      inode_list[p] = 1;
      nnode_comm ++;
    }
  }

  return nnode_comm;
  
}



void getOverlap( const int x0_s, const int x0_e,
		 const int x1_s, const int x1_e,
		 int &xs, int &xe){

  xs = -10;
  xe = -10;

  if( x1_s <= x0_s){
    xs = x0_s;
  }
  else if( x1_s <= x1_e){
    xs = x1_s;
  }
  if( x1_e >= x0_e){
    xe = x0_e;
  }
  else if( x1_e > x0_s){
    xe = x1_e;
  }

}



int getOverlapBox( const int *g_offset_all,
		   const int *l_msize_all,
		   const int *subbox_start,
		   const int *subbox_msize,
		   const int SIZE_MESH,
		   int (*ms)[2], int (*me)[2],
		   int (*m_corr)[2]){

  for( int j=0; j<3; j++){
    m_corr[j][0] = m_corr[j][1] = 0;

    int gs = g_offset_all[j];
    int ge = gs + l_msize_all[j];
    int s1 = 0;  
    int e1= 0;
    int s2 = gs;  
    int e2 = ge;
    if( gs < 0){
      s1 = gs + SIZE_MESH;
      m_corr[j][0] = -SIZE_MESH;
      e1 = SIZE_MESH;
      s2 = 0;
    }
    if( ge >= SIZE_MESH){
      s1 = 0;
      e1 = ge - SIZE_MESH;
      m_corr[j][0] = SIZE_MESH;
      e2 = SIZE_MESH;
    }
    getOverlap( s1, e1, subbox_start[j], subbox_start[j]+subbox_msize[j], ms[j][0], me[j][0]);
    getOverlap( s2, e2, subbox_start[j], subbox_start[j]+subbox_msize[j], ms[j][1], me[j][1]);
  }

  int nsend = 0;
  int dmx0 = me[0][0] - ms[0][0];
  if( dmx0 < 0)  dmx0 = 0;
  int dmx1 = me[0][1] - ms[0][1];
  if( dmx1 < 0)  dmx1 = 0;
  int dmy0 = me[1][0] - ms[1][0];
  if( dmy0 < 0)  dmy0 = 0;
  int dmy1 = me[1][1] - ms[1][1];
  if( dmy1 < 0)  dmy1 = 0;
  int dmz0 = me[2][0] - ms[2][0];
  if( dmz0 < 0)  dmz0 = 0;
  int dmz1 = me[2][1] - ms[2][1];
  if( dmz1 < 0)  dmz1 = 0;
  nsend += dmx0 * dmy0 * dmz0;
  nsend += dmx0 * dmy0 * dmz1;
  nsend += dmx0 * dmy1 * dmz0;
  nsend += dmx0 * dmy1 * dmz1;
  nsend += dmx1 * dmy0 * dmz0;
  nsend += dmx1 * dmy0 * dmz1;
  nsend += dmx1 * dmy1 * dmz0;
  nsend += dmx1 * dmy1 * dmz1;

  return nsend;

}



void PMForce::commMeshPhi3D(){
#ifdef FFT3D

  double nowtime;
  fprintf_verbose( stderr, "#start bcast\n");
  My_MPI_Barrier( MPI_COMM_REDUCE);
  getTime(&nowtime);
  MPI_Bcast( mesh_density_slab, msize_slab, MPI_TYPE_FFTW,
	     0, MPI_COMM_REDUCE);
  My_MPI_Barrier( MPI_COMM_REDUCE);
  this_run->t.pm_phi_comm_bcast += getTime(&nowtime);
  fprintf_verbose( stderr, "#end bcast\n");

  getLocalMeshParam(3);
  static int g_offset_all[MAXNNODE][3];
  static int l_msize_all[MAXNNODE][3];
  MPI_Allgather( g_offset, 3, MPI_INT, 
		 g_offset_all, 3, MPI_INT, MPI_COMM_COLUMN);
  MPI_Allgather( l_msize, 3, MPI_INT, 
		 l_msize_all, 3, MPI_INT, MPI_COMM_COLUMN);

  static int nsend[MAXNNODE];
  static int nrecv[MAXNNODE];
  static int sdispls_sendbuf[MAXNNODE];
  static int rdispls_sendbuf[MAXNNODE];
  float *recvbuf = new float[ln_total];
  int *recvbuf_local_dest = new int[ln_total];
  for( int i=0; i<nnode_reduce; i++)  nsend[i] = 0;

  int *inode_list = new int[nnode_reduce];
  int nnode_comm = determineBoxOverlapping2( subbox_start, subbox_msize,
					     g_offset_all, l_msize_all,
					     SIZE_MESH, nnode_reduce,
					     inode_list);

  int p_nsend_all = 0;
  for( int p=0; p<nnode_reduce; p++){
    if( inode_list[p] == 0)  continue;
    int ms[3][2], me[3][2], m_corr[3][2];
    nsend[p] = getOverlapBox( g_offset_all[p], l_msize_all[p],
			      subbox_start, subbox_msize,
			      SIZE_MESH, ms, me, m_corr);
    p_nsend_all += nsend[p];
  }

  float *sendbuf = new float[p_nsend_all];
  int *sendbuf_local_dest = new int[p_nsend_all];
  int *dest = new int [p_nsend_all];
  sdispls_sendbuf[0] = 0;
  for( int i=1; i<nnode_reduce; i++){
    sdispls_sendbuf[i] = sdispls_sendbuf[i-1] + nsend[i-1];
  }

#pragma omp parallel for CHUNK_PM
  for( int p=0; p<nnode_reduce; p++){
    if( inode_list[p] == 0)  continue;
    int p_ln_total = l_msize_all[p][0] * l_msize_all[p][1] * l_msize_all[p][2];
    int ms[3][2], me[3][2], m_corr[3][2];
    getOverlapBox( g_offset_all[p], l_msize_all[p],
		   subbox_start, subbox_msize,
		   SIZE_MESH, ms, me, m_corr);

    int ib = sdispls_sendbuf[p];
    for( int ni=0; ni<2; ni++){
      for( int i=ms[0][ni]; i<me[0][ni]; i++){
	int li = i + m_corr[0][ni] - g_offset_all[p][0];
	for( int nj=0; nj<2; nj++){
	  for( int j=ms[1][nj]; j<me[1][nj]; j++){
	    int lj = j + m_corr[1][nj] - g_offset_all[p][1];
	    for( int nk=0; nk<2; nk++){
	      for( int k=ms[2][nk]; k<me[2][nk]; k++){
		int lk = k + m_corr[2][nk] - g_offset_all[p][2];
		int im = getSubboxID( i-subbox_start[0], j-subbox_start[1], k-subbox_start[2]);
		assert( im < subbox_total);
		sendbuf[ib] = mesh_density_slab[im];
		sendbuf_local_dest[ib] = lk + l_msize_all[p][2]*( lj + l_msize_all[p][1]*li);
		assert( sendbuf_local_dest[ib] < p_ln_total);
		dest[ib] = p;
		ib ++;
	      }
	    }
	  }
	}	    
      }
    }
  }

  delete [] inode_list;

  My_MPI_Barrier( MPI_COMM_COLUMN);
  this_run->t.pm_make_comm_phi += getTime(&nowtime);

  My_MPI_Barrier( MPI_COMM_COLUMN);
  //fprintf_verbose( stderr, "#start comm\n");
  My_MPI_Barrier( MPI_COMM_COLUMN);
  MPI_Alltoall( nsend, 1, MPI_INT, nrecv, 1, MPI_INT, MPI_COMM_COLUMN);
  My_MPI_Barrier( MPI_COMM_COLUMN);
  //fprintf_verbose( stderr, "#start comm0\n");
  My_MPI_Barrier( MPI_COMM_COLUMN);
  int nrecv_total = nrecv[0];
  rdispls_sendbuf[0] = 0;
  for( int i=1; i<nnode_reduce; i++){
    rdispls_sendbuf[i] = rdispls_sendbuf[i-1] + nrecv[i-1];
    nrecv_total += nrecv[i];
  }
  //cerr << nrecv_total << "\t" << ln_total << endl;
  assert( nrecv_total == ln_total);

  int nscomm_pm_s2l = 0;
  int nrcomm_pm_s2l = 0;
  int nsend_pm_s2l = 0;
  int nrecv_pm_s2l = 0;
  for( int i=0; i<nnode_reduce; i++){
    if( nsend[i] != 0)  nscomm_pm_s2l ++;
    if( nrecv[i] != 0)  nrcomm_pm_s2l ++;
    nsend_pm_s2l += nsend[i];
    nrecv_pm_s2l += nrecv[i];
  }
  this_run->nscomm_pm_s2l = nscomm_pm_s2l;
  this_run->nrcomm_pm_s2l = nrcomm_pm_s2l;
  this_run->nsend_pm_s2l = nsend_pm_s2l;
  this_run->nrecv_pm_s2l = nrecv_pm_s2l;

  My_MPI_Barrier( MPI_COMM_COLUMN);
  //fprintf_verbose( stderr, "#end comm0\n");
  My_MPI_Barrier( MPI_COMM_COLUMN);
  MPI_Alltoallv( sendbuf, nsend, sdispls_sendbuf, MPI_FLOAT,
                 recvbuf, nrecv, rdispls_sendbuf, MPI_FLOAT,
                 MPI_COMM_COLUMN);
  MPI_Alltoallv( sendbuf_local_dest, nsend, sdispls_sendbuf, MPI_INT,
                 recvbuf_local_dest, nrecv, rdispls_sendbuf, MPI_INT,
                 MPI_COMM_COLUMN);
  My_MPI_Barrier( MPI_COMM_COLUMN);
  this_run->t.pm_phi_comm += getTime(&nowtime);
  //fprintf_verbose( stderr, "#end comm\n");

#pragma omp parallel for
  for( int i=0; i<nrecv_total; i++){
    int im = recvbuf_local_dest[i];
    mesh_density_local[im] = recvbuf[i];
  }

  slab_to_local_sendbuf_size = p_nsend_all;
  delete [] sendbuf;
  delete [] sendbuf_local_dest;
  delete [] recvbuf;
  delete [] recvbuf_local_dest;
  delete [] dest;

  this_run->t.pm_make_comm_phi += getTime(&nowtime);
#endif

}





#define DRHO(a,b) (mesh_density_local[a]-mesh_density_local[b])
void PMForce::calcMeshForce(){

  static double fac1 = -2.0 * SIZE_MESH / 3.0;
  static double fac2 = SIZE_MESH / 12.0;

#pragma omp parallel for
  for( int i=2; i<(l_msize[0]-2); i++){
    for( int j=2; j<(l_msize[1]-2); j++){
      for( int k=2; k<(l_msize[2]-2); k++){
	int li = k + l_msize[2]*( j + l_msize[1]*i);
	assert( li >= 0);
	assert( li < ln_total);
	int lix2p = li + lnyz*2;
	int lix1p = li + lnyz;
	int lix2m = li - lnyz*2;
	int lix1m = li - lnyz;
	mesh_force_local[li][0] = fac1*DRHO(lix1p,lix1m) + fac2*DRHO(lix2p,lix2m);
	int liy2p = li + l_msize[2]*2;
	int liy1p = li + l_msize[2];
	int liy2m = li - l_msize[2]*2;
	int liy1m = li - l_msize[2];
	mesh_force_local[li][1] = fac1*DRHO(liy1p,liy1m) + fac2*DRHO(liy2p,liy2m);
	int liz2p = li + 2;
	int liz1p = li + 1;
	int liz2m = li - 2;
	int liz1m = li - 1;
	mesh_force_local[li][2] = fac1*DRHO(liz1p,liz1m) + fac2*DRHO(liz2p,liz2m);
      }
    }
  }

}



void PMForce::forceInterpolation(const Particle *particle, float *a){

  float pos[3];
  getPos2( particle, pos);
  pos[0] -= g_pos_f[0];
  pos[1] -= g_pos_f[1];
  pos[2] -= g_pos_f[2];

  float wi[3][3];
  int iw[3][3];

  for( int j=0; j<3; j++){
    float xt1  = pos[j] * (float)SIZE_MESH;
    iw[1][j] = (int)(xt1 + 0.5);
    float dx1 = xt1 - (float)(iw[1][j]);
    wi[0][j] = 0.5 * (0.5-dx1) * (0.5-dx1);
    wi[1][j] = 0.75 - dx1*dx1;
    wi[2][j] = 0.5 * (0.5 + dx1) * (0.5 + dx1);
    iw[0][j] = iw[1][j] - 1;
    iw[2][j] = iw[1][j] + 1;
  }

  float ax = 0.0;
  float ay = 0.0;
  float az = 0.0;
  for( int i=0; i<3; i++){
    for( int j=0; j<3; j++){
      for( int k=0; k<3; k++){
         int ii = iw[k][2] + l_msize[2]*( iw[j][1] + l_msize[1]*iw[i][0]);
         float w = wi[i][0]*wi[j][1]*wi[k][2];
         ax += mesh_force_local[ii][0]*w;
         ay += mesh_force_local[ii][1]*w;
         az += mesh_force_local[ii][2]*w;
      }
    }
  }

  a[0] = ax;
  a[1] = ay;
  a[2] = az;


}

void PMForce::potentialInterpolation(const Particle *particle, float *pot){

   float pos[3];
   getPos2( particle, pos);
   pos[0] -= g_pos_f[0];
   pos[1] -= g_pos_f[1];
   pos[2] -= g_pos_f[2];

   float wi[3][3];
   int iw[3][3];
   for(int j=0; j<3; j++){
     float xt1 = pos[j] * (float)SIZE_MESH;
     iw[1][j] = (int)(xt1 + 0.5);
     float dx1 = xt1 - (float)(iw[1][j]);
     wi[0][j] = 0.5 * (0.5-dx1) * (0.5-dx1);
     wi[1][j] = 0.75 - dx1*dx1;
     wi[2][j] = 0.5 * (0.5 + dx1) * (0.5 + dx1);
     iw[0][j] = iw[1][j] - 1;
     iw[2][j] = iw[1][j] + 1;
   }

   float ptnl = 0.0;
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         for(int k=0; k<3; k++){
            int ii = iw[k][2] + l_msize[2]*( iw[j][1] + l_msize[1]*iw[i][0]);
            float w = wi[i][0]*wi[j][1]*wi[k][2];
            ptnl += mesh_density_local[ii]*w;
         }
      }
   }

   *pot = ptnl;

}


void PMForce::outputLog( FILE *outstream){

  double memsize = 1.0e-6 * sizeof(TYPE_FFTW) * ln_total * 4;
  double memsize2 = 1.0e-6 * sizeof(TYPE_FFTW) * msize_slab * 2;
  double memsize3 = 1.0e-6 * sizeof(TYPE_FFTW) * msize_gk;
  double recvbufsize = 1.0e-6 * sizeof(float) * local_to_slab_recvbuf_size;
  double sendbufsize = 1.0e-6 * sizeof(float) * slab_to_local_sendbuf_size;

  sprintf( cbuf, "%e\t%e\t%e\t%d\t%d\t%d\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n", 
	   g_pos[0], g_pos[1], g_pos[2], 
	   l_msize[0], l_msize[1], l_msize[2],
	   memsize, memsize2+memsize3,
	   recvbufsize, sendbufsize);
  mp_print( cbuf, this_run->MPI_COMM_INTERNAL, outstream);


}






void PMForce::calcPMMeshForce(const Particle *particle, const int npart){

   static int first_flag = 0;
   double nowtime = 0.0;
   getTime(&nowtime);
   initLocal();

   //fprintf_verbose( stderr, "#PM setLocalMeshDensity start\n");
   setLocalMeshDensity( particle, npart);
   this_run->t.pm_density_assignment = getTime(&nowtime);

   //fprintf_verbose( stderr, "#PM commMeshDensity start\n");
#ifdef FFT3D
   commMeshDensity3D();
#else
   commMeshDensity();
#endif

   getTime(&nowtime);

   //fprintf_verbose( stderr, "#PM fftForward start\n");
   fftForward();
   this_run->t.pm_fft_forward = getTime(&nowtime);


   //fprintf_verbose( stderr, "#PM calcMeshPhi start\n");
   calcMeshPhi();
   this_run->t.pm_density_phi = getTime(&nowtime);

   //fprintf_verbose( stderr, "#PM fftBackward start\n");
   fftBackward();
   this_run->t.pm_fft_backward = getTime(&nowtime);

   //fprintf_verbose( stderr, "#PM commMeshPhi start\n");
   if (this_run->nnode != 1){
#ifdef FFT3D
      commMeshPhi3D();
#else
      commMeshPhi();
#endif
   }

   //fprintf_verbose( stderr, "#PM calcMeshForce start\n");
   getTime(&nowtime);
   calcMeshForce();
   this_run->t.pm_mesh_force = getTime(&nowtime);

   if (PMLOG_ON == 1){
      FILE *fout = my_fopen( PMLOG, "a+");
      if (first_flag == 0){
         if (this_run->inode == 0){
            fprintf(fout, "node, g_pos[0,1,2], msize[0,1,2], memsize_local[Mbyte], memsize_slab, recvmemsize, sendmemsize\n");
         }
         sprintf(cbuf,"%d\t%d\t%d\t%d\t%d\n",
                 local_nx, local_x_start,
                 local_ny_after_transpose,
                 local_y_start_after_transpose,
                 total_local_size);
         mp_print(cbuf, this_run->MPI_COMM_INTERNAL, fout);
         first_flag = 1;
      }
      outputLog(fout);
      fclose(fout);
   }

}


    } // END of namespace ParticleMesh
} // END of namespace ParticleSimulator
