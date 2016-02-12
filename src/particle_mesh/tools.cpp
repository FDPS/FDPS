#include "tools.h"

namespace ParticleSimulator{
    namespace ParticleMesh{


static char cbuf[256];



double getTimeDiff( double now_time){

  struct timeval tv;
  gettimeofday( &tv, NULL);
  double diff_time;

  diff_time = (tv.tv_sec + (double)tv.tv_usec*1e-6) - now_time;

  return diff_time;

}



double getTimePrint( double *now_time, const char *command_name){

  struct timeval tv;
  gettimeofday( &tv, NULL);
  double diff_time;

  diff_time = (tv.tv_sec + (double)tv.tv_usec*1e-6) - *now_time;
  *now_time = (tv.tv_sec + (double)tv.tv_usec*1e-6);

  fprintf( stderr, "-%s: %lf sec\n", command_name, diff_time);

  return diff_time;

}



double getTimeFprint( double *now_time, const char *command_name, FILE *fp){

  struct timeval tv;
  gettimeofday( &tv, NULL);
  double diff_time;

  diff_time = (tv.tv_sec + (double)tv.tv_usec*1e-6) - *now_time;
  *now_time = (tv.tv_sec + (double)tv.tv_usec*1e-6);

  fprintf( fp, "%s: %lf sec\n", command_name, diff_time);

  return diff_time;

}



void fileCheck( char const *fname){

  FILE *fp;
  if( (fp = fopen( fname, "r")) == NULL){
    perror( fname);
    exit(1);
  }
  fclose(fp);

}



void dirCheck( char *const dname){

  DIR *dp;
  if( (dp = opendir( dname)) == NULL){
    sprintf( cbuf, "mkdir %s", dname);
    fprintf( stderr, "*%s\n", cbuf);
    //system( cbuf);
    mode_t mode = S_IRWXU;
    mkdir( dname, mode);
  }
  closedir(dp);

}



FILE *my_fopen(const char *path, const char *mode){

  FILE *fp;
  fp = fopen( path, mode);

  if( fp == NULL){
    perror(path);
    exit(EXIT_FAILURE);
  }

  return fp;

}



void createWritingFile( const char *filename){

  sprintf( cbuf, "%s.writing", filename);

  FILE *outstream = fopen( cbuf, "w");
  fclose(outstream);

}



void removeWritingFile( const char *filename){

  sprintf( cbuf, "%s.writing", filename);
  remove(cbuf);

}



void fprintf_verbose( FILE *fout, const char *format, ...){

#ifdef VERBOSE_MODE2
  MPI_Barrier( MPI_COMM_WORLD);
#endif

#ifdef VERBOSE_MODE
  va_list ap;
  va_start(ap, format);
  vfprintf( fout, format, ap);
  va_end(ap);
#endif

}



void fprintf_verbose2( MPI_Comm comm, FILE *fout, const char *format, ...){

#ifdef VERBOSE_MODE2
  MPI_Barrier( comm);
  va_list ap;
  va_start(ap, format);
  vfprintf( fout, format, ap);
  va_end(ap);
#endif

}




void *my_malloc( const size_t size){

  void *p;
  p = malloc( size);
  if( p == NULL){
    fprintf( stderr, "malloc error size %ld\n", sizeof(size));
    exit(1);
  }

  return p;

}




void mpiCheck( const int error_class, const int inode){

  if( error_class == MPI_ERR_BUFFER){
    fprintf( stderr, "node:%d MPI_ERR_BUFFER", inode);
    exit(1);
  }

  if( error_class == MPI_ERR_COUNT){
    fprintf( stderr, "node:%d MPI_ERR_COUNT", inode);
    exit(1);
  }

  if( error_class == MPI_ERR_TYPE){
    fprintf( stderr, "node:%d MPI_ERR_TYPE", inode);
    exit(1);
  }

  if( error_class == MPI_ERR_TAG){
    fprintf( stderr, "node:%d MPI_ERR_TAG", inode);
    exit(1);
  }

  if( error_class == MPI_ERR_RANK){
    fprintf( stderr, "node:%d MPI_ERR_RANK", inode);
    exit(1);
  }

  if( error_class == MPI_ERR_COMM){
    fprintf( stderr, "node:%d MPI_ERR_COMM", inode);
    exit(1);
  }

}



#if 1
int mp_print(char *msg, MPI_Comm comm, FILE *outstream){

  static char buf[MAXNNODE][CHARMAX];

  if( (int)strlen(msg) >= CHARMAX){
    cerr << "message size is larger than CHARMAX in mp_print" << endl;
    cerr << msg << endl;
  }

  int inode, nnode;
  MPI_Comm_rank(comm, &inode);
  MPI_Comm_size(comm, &nnode);
  int err = MPI_Gather( msg, CHARMAX, MPI_CHAR, buf, CHARMAX, MPI_CHAR, 0, comm);

  if( inode == 0){
    for( int i=0; i<nnode; i++){
      char *buf2 = buf[i];
      buf2[CHARMAX-1] = '\0';
      fprintf( outstream, "%d\t%s", i, buf2);
    }
  }

  return err;

}
#else
int mp_print(char *msg, MPI_Comm comm, FILE *outstream){

  int err = 0;
  int i;
  int len;
  int nprocs, myrank;
  static char buf[256];
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myrank);

  if(myrank == 0){
    fprintf(outstream, "%d\t%s", 0, msg);
    for(i=1; i<nprocs; i++){
      MPI_Status stat;
      err |= MPI_Recv(&len, 1, MPI_INT, i, 1, comm, &stat);
      err |= MPI_Recv(buf, len+1, MPI_CHAR, i, 2, comm, &stat);
      buf[255] = '\0';
      fprintf(outstream, "%d\t%s", i, buf);
    }
  }else{
    len = (int)strlen(msg);
    len = len<256 ? len : 255;
    err |= MPI_Send(&len, 1, MPI_INT, 0, 1, comm);
    err |= MPI_Send(msg, len+1, MPI_CHAR, 0, 2, comm);
  }
  return err;
}
#endif



void exitMPI(){

  MPI_Finalize();
  exit(1);

}



int setNthreadsOpenMP(){

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = NUMBER_OF_OMP_THREADS;
  omp_set_num_threads(NUMBER_OF_OMP_THREADS);
#endif

  return nthreads;

}



int getDevidOpenMP(){

  int devid = 0;
#ifdef _OPENMP
  devid = omp_get_thread_num();
#endif

  return devid;

}


    } // namespace ParticleMesh
}     // namespace ParticleSimulator
