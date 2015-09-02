#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "timeprof.h"

double get_dtime(void){
  struct timeval tv;

  gettimeofday(&tv, NULL);
  return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 0.001 * 0.001);
}
