#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>

#if 0

/* CPU time */
void
get_cputime(double *laptime, double *sprittime)
{
    struct rusage x;
    double sec,microsec;

    getrusage(RUSAGE_SELF,&x);
    sec = x.ru_utime.tv_sec + x.ru_stime.tv_sec ;
    microsec = x.ru_utime.tv_usec + x.ru_stime.tv_usec ;

    *laptime = sec + microsec / 1000000.0 - *sprittime;
    *sprittime = sec + microsec / 1000000.0;
}

#else

/* elapsed time in real world */
void
get_cputime(double *laptime, double *sprittime)
{
    struct timeval x;
    double sec,microsec;

    gettimeofday(&x, NULL);
    sec = x.tv_sec;
    microsec = x.tv_usec;

    *laptime = sec + microsec / 1000000.0 - *sprittime;
    *sprittime = sec + microsec / 1000000.0;
}

#endif
