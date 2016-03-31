#ifndef _PM_PARALLEL_INCLUDED
#define _PM_PARALLEL_INCLUDED

#include <iostream>
#include <cassert>

using namespace std;

#include "treepm_header.h"
#include "qsort_omp.h"



#ifdef FFTW3_PARALLEL

#include <fftw3-mpi.h>
#ifdef FFTW_DOUBLE
#define TYPE_FFTW double
#define MPI_TYPE_FFTW MPI_DOUBLE
#else
#define TYPE_FFTW float
#define MPI_TYPE_FFTW MPI_FLOAT
#define fftw_complex fftwf_complex
#define fftw_plan fftwf_plan
#define fftw_mpi_init fftwf_mpi_init
#define fftw_mpi_local_size_3d fftwf_mpi_local_size_3d
#define fftw_plan_dft_r2c_3d fftwf_plan_dft_r2c_3d
#define fftw_plan_dft_c2r_3d fftwf_plan_dft_c2r_3d
#define fftw_mpi_plan_dft_r2c_3d fftwf_mpi_plan_dft_r2c_3d
#define fftw_mpi_plan_dft_c2r_3d fftwf_mpi_plan_dft_c2r_3d
#define fftw_destroy_plan fftwf_destroy_plan
#define fftw_execute fftwf_execute
#define fftw_mpi_broadcast_wisdom fftwf_mpi_broadcast_wisdom
#define fftw_forget_wisdom fftwf_forget_wisdom
#define fftw_import_wisdom_from_filename fftwf_import_wisdom_from_filename
#endif

#else // if FFTW3_PARALLEL is not defined

#ifdef FFTW_DOUBLE
#define TYPE_FFTW double
#define MPI_TYPE_FFTW MPI_DOUBLE
#ifdef KCOMPUTER
#include <rfftw_mpi.h>
#else
#include <drfftw_mpi.h>
#endif
#else
#define TYPE_FFTW float
#define MPI_TYPE_FFTW MPI_FLOAT
#include <srfftw_mpi.h>
#endif

#endif //FFTW3_PARALLEL

#include<complex>
namespace ParticleSimulator{
    namespace ParticleMesh{

#ifdef KCOMPUTER
extern "C"{
  void ds_v3dcft3_( double *x,
                    int *kx1, int *kx2,
		    int *kx1p, int *kx2p, int *kx3p,
                    int *n1, int *n2, int *n3,
                    int *nd1, int *nd2, int *nd3,
                    double *w, int *nw,
                    int *isn, MPI_Fint *comm, int *icon);
}
#endif



class PMForce{

 public:

  /* initialized once  */
  RunParam *this_run;

  int SIZE_MESH;
  int SIZE_MESH_P2;

  int local_nx;            // size of mesh for X direction
  int local_x_start;       // offset of mesh for X direction
  int local_ny_after_transpose;
  int local_y_start_after_transpose;
  int total_local_size;
  int *l_nx_all;        // size of mesh for all processes
  int *l_x_start_all;   // offset of mesh for all processes
  int msize_slab;       // slab-mesh size i.e. local_nx*SIZE_OF_MESH*SIZE_OF_MESH_P2
  int msize_gk;         // gk table size i.e. local_nx*sizeg*sizeg

#ifdef FFTW3_PARALLEL
  fftw_plan fplan;
  fftw_plan bplan;
#else
  rfftwnd_mpi_plan fplan;
  rfftwnd_mpi_plan bplan;
#endif

  TYPE_FFTW *mesh_density_slab;
  TYPE_FFTW *fftwork;
  TYPE_FFTW *gk;            // green function table

  int init_gk;
  int fft_node;
  int ndiv_fft[3];
  int inode_fft;
  int nnode_fft;
  MPI_Comm MPI_COMM_FFT;

#ifdef RMM_PM
  int reduce_nx;
  int reduce_ny;
  int reduce_nz;
  int ndiv_reduce[3];
  int reduce_node;
  int inode_reduce;
  int nnode_reduce;
  MPI_Comm MPI_COMM_COLUMN;
  MPI_Comm MPI_COMM_REDUCE;
  int local_nx2;
  int local_x_start2;
  int *l_nx_all2;
  int *l_x_start_all2;
#endif
  TYPE_FFTW *mesh_density_slab0;  // rmm

  /* FFT3D */
  int fft3d;
  int subbox_total;
  int subbox_msize[3];  // common to all fft nodes
  int subbox_start[3];
  TYPE_FFTW (*mesh_density_subbox)[2];
  TYPE_FFTW (*fft3d_work)[2];
  MPI_Fint MPI_COMM_FFT_F;

  /* updated every step */
  double g_pos[3];      // position of top-left mesh
  float  g_pos_f[3];      // position of top-left mesh
  int g_offset[3];      // offset of local mesh from origin of entire mesh (0,0,0)
  int l_msize[3];       // size of local mesh
  int lnyz;             // l_msize[1] * l_msize[2]
  int ln_total;         // l_msize[0] * l_msize[1] * l_msize[2]
  float *mesh_density_local;
  float (*mesh_force_local)[3];

  /* log */
  int local_to_slab_recvbuf_size;  // the size of recv buffer(float) 
  int slab_to_local_sendbuf_size;  // the size of send buffer(float) 

  /* for visualization */
  int output_mesh_flag;
  float *mesh_density_local2;
  float *mesh_density_slab2;


  PMForce( RunParam *i_this_run);
  virtual ~PMForce();

  void PMForce0_fftw3( MPI_Comm comm);
  void PMForce1_fftw3( MPI_Comm comm,
		       const char *fname_wisdom_f, const char *fname_wisdom_b);
  void PMForce2_fftw3( MPI_Comm comm);
  void PMForce0( const int *ndiv_fft_input);
  void initPMReduce( const int x, const int y, const int z);
  void init( RunParam *i_this_run, const int _SIZE_MESH, 
	     const int *ndiv_fft_input, const int flag_wisdom);
  virtual void initGK();
  virtual void initGK3D();
  virtual void initLocal();
  virtual void delLocal();
  virtual void getLocalMeshParam( const double padsize);
  virtual void checkMeshSlab();
  virtual void checkMeshLocal();
  virtual void checkMeshGreen();
  virtual void setLocalMeshDensity( const Particle *particle, const int npart);
  virtual void commMeshDensity();
  virtual void commMeshDensity3D();
  virtual void fftForward();
  virtual void calcMeshPhi();
  virtual void fftBackward();
  virtual void commMeshPhi();
  virtual void commMeshPhi3D();
  virtual void calcMeshForce();
  virtual void forceInterpolation(const Particle *particle, float *a);
  virtual void forceInterpolation( const double *_pos, float *a);
  virtual void potentialInterpolation(const Particle *particle, float *pot);
  virtual void outputLog( FILE *outstream);
  virtual void calcPMMeshForce(const Particle *particle, const int npart);

  int getLocalID( const int x, const int y, const int z) const;
  int getSubboxID( const int x, const int y, const int z) const;

};



inline void PMForce::forceInterpolation( const double *_pos, float *a){

  float pos[3];
  pos[0] = _pos[0] - g_pos_f[0];
  pos[1] = _pos[1] - g_pos_f[1];
  pos[2] = _pos[2] - g_pos_f[2];

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



inline int PMForce::getLocalID( const int x, const int y, const int z) const{

  return z + l_msize[2]*( y + l_msize[1]*x);

}



inline int PMForce::getSubboxID( const int x, const int y, const int z) const{

  int r = 0;

  r = z + subbox_msize[0]*( y + subbox_msize[1]*x);

  return r;

}


    } // namespace ParticleMesh
}     // namespace ParticleSimulator

#endif

