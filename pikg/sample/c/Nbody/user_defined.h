#ifndef H_USER_DEFINED
#define H_USER_DEFINED

#include <pikg_vector.h>

typedef struct {
    long long id;
    double  mass;
    double  eps;
    pikg_f64vec pos;
    pikg_f64vec vel;
    double  pot;
    pikg_f64vec acc;
}particle;

#endif
