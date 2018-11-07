#pragma once
/* Standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
/* FDPS headers */
#include "FDPS_c_if.h"
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"

/* Force types */
typedef struct force_grav { //$fdps Force
   //$fdps clear
   fdps_f64vec acc;
   double pot;
} Force_grav;

typedef struct force_dens { //$fdps Force
   //$fdps clear smth=keep
   int flag;
   double dens;
   double smth;
   double gradh;
   double divv;
   fdps_f64vec rotv;
} Force_dens;

typedef struct force_hydro { //$fdps Force
   //$fdps clear 
   fdps_f64vec acc;
   double eng_dot;
   double ent_dot;
   double dt;
} Force_hydro;

/* Full particle types */
typedef struct fp_nbody { //$fdps FP
   //$fdps copyFromForce force_grav (acc,acc) (pot,pot)
   long long int id; //$fdps id
   double mass; //$fdps charge
   fdps_f64vec pos; //$fdps position
   fdps_f64vec vel;
   fdps_f64vec acc;
   double pot;
} FP_nbody;

typedef struct fp_sph { //$fdps FP
   //$fdps copyFromForce force_grav (acc,acc_grav) (pot,pot_grav)
   //$fdps copyFromForce force_dens (flag,flag) (dens,dens) (smth,smth) (gradh,gradh) (divv,divv) (rotv,rotv)
   //$fdps copyFromForce force_hydro (acc,acc_hydro) (eng_dot,eng_dot) (ent_dot,ent_dot) (dt,dt)
   long long id; //$fdps id
   double mass; //$fdps charge
   fdps_f64vec pos; //$fdps position
   fdps_f64vec vel; 
   fdps_f64vec acc_grav;
   double pot_grav;
   fdps_f64vec acc_hydro;
   int flag;
   double dens;
   double eng;
   double ent;
   double pres;
   double smth;
   double gradh;
   double divv;
   fdps_f64vec rotv;
   double balsw;
   double snds;
   double eng_dot;
   double ent_dot;
   double dt;
   fdps_f64vec vel_half;
   double eng_half;
   double ent_half;
} FP_sph;

/* Essential particle types */
typedef struct ep_grav { //$fdps EPI,EPJ
   //$fdps copyFromFP fp_nbody (id,id) (mass,mass) (pos,pos)
   //$fdps copyFromFP fp_sph (id,id) (mass,mass) (pos,pos)
   long long id; //$fdps id
   double mass; //$fdps charge
   fdps_f64vec pos; //$fdps position
} EP_grav;

typedef struct ep_hydro { //$fdps EPI,EPJ
   //$fdps copyFromFP fp_sph (id,id) (pos,pos) (vel,vel) (mass,mass) (smth,smth) (dens,dens) (pres,pres) (gradh,gradh) (snds,snds) (balsw,balsw)
   long long id; //$fdps id
   fdps_f64vec pos; //$fdps position
   fdps_f64vec vel;
   double mass; //$fdps charge
   double smth; //$fdps rsearch
   double dens;
   double pres;
   double gradh;
   double snds;
   double balsw;
} EP_hydro;

/* Prototype declarations */
double W(double r, double h);
fdps_f64vec gradW(fdps_f64vec dr, double h);
double dWdh(double r, double h);

void calc_gravity_ep_ep(struct ep_grav *ep_i,
                        int n_ip,
                        struct ep_grav *ep_j,
                        int n_jp,
                        struct force_grav *f);

void calc_gravity_ep_sp(struct ep_grav *ep_i,
                        int n_ip,
                        fdps_spj_monopole *ep_j,
                        int n_jp,
                        struct force_grav *f);

void calc_density(struct ep_hydro *ep_i,
                  int n_ip,
                  struct ep_hydro *ep_j,
                  int n_jp,
                  struct force_dens *f);

void calc_hydro_force(struct ep_hydro *ep_i,
                      int n_ip,
                      struct ep_hydro *ep_j,
                      int n_jp,
                      struct force_hydro *f);

/* Global variables */
extern double specific_heat_ratio;
extern double CFL_dyn;
extern double CFL_hydro;
extern double SCF_smth;
extern int N_neighbor;
extern double eps_grav;
extern double mass_avg;
extern double dt_max;

