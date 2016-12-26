#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>

#define MAXPATHLEN (1024)

void make_directory(char *directory_name)
{
  mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;

  //  char *cwd_path;
  //  cwd_path = get_current_dir_name();

  static char cwd_path[MAXPATHLEN];
  if(getcwd(cwd_path, sizeof(cwd_path)) == NULL) {
    fprintf(stderr,"Error in make_directory");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } 
  

  strcat(cwd_path,"/");
  strcat(cwd_path, directory_name);

  mkdir(cwd_path, mode);

  //  free(cwd_path);
}
