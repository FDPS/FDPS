#pragma once
/* Standard headers */
#include <math.h>
/* FDPS headers */
#include "FDPS_c_if.h"

typedef struct full_particle { //$fdps FP,EPI,EPJ,Force
    //$fdps copyFromForce full_particle (pot,pot) (acc,acc)
    //$fdps copyFromFP full_particle (id,id) (mass,mass) (eps,eps) (pos,pos) 
    //$fdps clear id=keep, mass=keep, eps=keep, pos=keep, vel=keep
    long long id; // $fdps id
    double  mass; //$fdps charge
    double  eps;
    fdps_f64vec pos; //$fdps position
    fdps_f64vec vel; //$fdps velocity
    double  pot;
    fdps_f64vec acc;
} Full_particle;

void  calc_gravity_ep_ep(Full_particle *ep_i,
                         int n_ip,
                         Full_particle *ep_j,
                         int n_jp,
                         Full_particle *f);

void  calc_gravity_ep_sp(Full_particle *ep_i,
                         int n_ip,
                         fdps_spj_monopole *ep_j,
                         int n_jp,
                         Full_particle *f);
