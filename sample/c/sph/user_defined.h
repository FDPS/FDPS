#pragma once
/* Standard headers */
#include <math.h>
/* FDPS headers */
#include "FDPS_c_if.h"
/* User-defined headers */
#include "mathematical_constants.h"

/* Force types */
typedef struct force_dens { //$fdps Force
    //$fdps clear smth=keep
    double dens;
    double smth;
} Force_dens;

typedef struct force_hydro { //$fdps Force
    //$fdps clear 
    fdps_f64vec acc;
    double eng_dot;
    double dt;
} Force_hydro;

/* Full particle type */
typedef struct full_particle { //$fdps FP
    //$fdps copyFromForce force_dens (dens,dens)
    //$fdps copyFromForce force_hydro (acc,acc) (eng_dot,eng_dot) (dt,dt)
    double mass; //$fdps charge
    fdps_f64vec pos; //$fdps position
    fdps_f64vec vel; 
    fdps_f64vec acc;
    double dens;
    double eng;
    double pres;
    double smth; //$fdps rsearch
    double snds;
    double eng_dot;
    double dt;
    long long int id;
    fdps_f64vec vel_half;
    double eng_half;
} Full_particle;

/* Essential particle type */
typedef struct essential_particle { //$fdps EPI,EPJ
    //$fdps copyFromFP full_particle (id,id) (pos,pos) (vel,vel) (mass,mass) (smth,smth) (dens,dens) (pres,pres) (snds,snds)
    long long int id;
    fdps_f64vec pos; //$fdps position
    fdps_f64vec vel;
    double mass; //$fdps charge
    double smth; //$fdps rsearch
    double dens;
    double pres;
    double snds;
} Essential_particle;

/* Prototype declarations */
double W(fdps_f64vec dr, double h);
fdps_f64vec gradW(fdps_f64vec dr, double h);

void calc_density(Essential_particle *ep_i,
                  int n_ip,
                  Essential_particle *ep_j,
                  int n_jp,
                  Force_dens *f);
void calc_hydro_force(Essential_particle *ep_i,
                      int n_ip,
                      Essential_particle *ep_j,
                      int n_jp,
                      Force_hydro *f);

/* Gloabl variable */
extern const double kernel_support_radius;
