#ifndef _TOOLS_INCLUDED
#define _TOOLS_INCLUDED

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cstdarg>
#include<cstring>
#include<iostream>
#include<dirent.h>
#include<unistd.h>
#include<sys/time.h>
#include<sys/resource.h>
#include<sys/stat.h>
#include<sys/types.h>

#include "mpi.h"
#include "param.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace std;

namespace ParticleSimulator{
    namespace ParticleMesh{

/* get time function (return nowtime-*now_time) */
double getTimeDiff( double now_time);
double getTimePrint( double *now_time, const char *command_name);
double getTimeFprint( double *now_time, const char *command_name, FILE *fp);

/* file function */
void fileCheck( char const *fname);
void dirCheck( char *const dname);
FILE *my_fopen(const char *path, const char *mode);
void createWritingFile( const char *filename);
void removeWritingFile( const char *filename);
void fprintf_verbose( FILE *fout, const char *format, ...);
void fprintf_verbose2( MPI_Comm comm, FILE *fout, const char *format, ...);


/* memory allocation */
void *my_malloc( const size_t size);

/* parallel */
void mpiCheck( const int error_class, const int inode);
int mp_print(char *msg, MPI_Comm comm, FILE *outstream);
void exitMPI();

/* OpenMP */
int setNthreadsOpenMP();
int getDevidOpenMP();



inline void My_MPI_Barrier( MPI_Comm comm){

#ifdef MY_MPI_BARRIER
  MPI_Barrier( comm);
#endif

}



inline double getTime( double *now_time){

  struct timeval tv;
  gettimeofday( &tv, NULL);
  double diff_time;

  diff_time = (tv.tv_sec + (double)tv.tv_usec*1e-6) - *now_time;
  *now_time = (tv.tv_sec + (double)tv.tv_usec*1e-6);

  return diff_time;

}




template<typename type> type retMax(const type *array, const int n){

  type max_array = array[0];

  for( int i=1; i<n; i++){
    if( max_array < array[i])  max_array = array[i];
  }

  return max_array;

}




template<typename type> type retMin(const type *array, const int n){

  type min_array = array[0];

  for( int i=1; i<n; i++){
    if( min_array > array[i])  min_array = array[i];
  }

  return min_array;

}



template<class Func>
inline double integralTrapezoid( const double xmin, const double xmax, const int nbin,
				 Func &func){

  double dx = (xmax - xmin)/(double)nbin;

  double y = 0.5 * func(xmin);
  double x = xmin + dx;
  for( int i=1; i<nbin; i++){
    y += func(x);
    x += dx;
  }
  y += 0.5 * func(x);

  return dx * y;

}



template<class Func>
inline double integralSimpson( const double xmin, const double xmax, const double err_threshold,
			       Func &func){

  double err = 1.0e30;
  int nbin = 1000;
  double t1 = integralTrapezoid(xmin, xmax, nbin, func);
  double t2 = integralTrapezoid(xmin, xmax, nbin*2, func);
  double yprev = 4.0*t2/3.0 - t1/3.0;
  t1 = t2;

  while( fabs(err) > err_threshold){
    nbin *= 2;
    double t2 = integralTrapezoid(xmin, xmax, nbin, func);
    double y = 4.0*t2/3.0 - t1/3.0;
    err = yprev - y;
    //    fprintf( stderr, "%d\t%e\t%e\t%e\t%e\t%e\n", nbin, t1, t2, y, yprev, err);
    t1 = t2;
    yprev = y;
  }

  return yprev;

}





    } // namespace ParticleMesh
}     // namespace ParticleSimulator

#endif
