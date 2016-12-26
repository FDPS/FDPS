#include <time.h>
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>

#ifndef CLK_TCK
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif

float timing(struct tms from, struct tms to) 
{
  return ((float)(to.tms_utime+to.tms_stime)-(float)(from.tms_utime+from.tms_stime))/(float)CLK_TCK;
}

float wallclock_timing(struct timeval from, struct timeval to)
{
  return ((to.tv_sec+to.tv_usec*1.e-6)-(from.tv_sec+from.tv_usec*1.e-6));
}
