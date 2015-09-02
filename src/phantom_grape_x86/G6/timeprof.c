#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "banana01.h"

REAL get_dtime(void){
  struct timeval tv;

  gettimeofday(&tv, NULL);
  return ((REAL)(tv.tv_sec) + (REAL)(tv.tv_usec) * 0.001 * 0.001);
}
