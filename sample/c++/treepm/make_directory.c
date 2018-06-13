#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

#define MAXPATHLEN (1024)

void make_directory(char *directory_name)
{
  // Set the absolute PATH of the directory
  static char cwd_path[MAXPATHLEN];
  if (getcwd(cwd_path, sizeof(cwd_path)) == NULL) {
      fprintf(stderr,"Error in make_directory");
      MPI_Finalize();
      exit(EXIT_FAILURE);
  } 
  strcat(cwd_path,"/");
  strcat(cwd_path, directory_name);

  // Make the directory at rank 0
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  if (my_rank == 0) {
     mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
     mkdir(cwd_path, mode);
  }

  // Sleep just for the safety
  struct timespec ts;
  ts.tv_sec = 0;
  ts.tv_nsec = 1000; // 1[\mu s]
  clock_nanosleep(CLOCK_MONOTONIC, 0, &ts, NULL);
  MPI_Barrier(MPI_COMM_WORLD);

}
