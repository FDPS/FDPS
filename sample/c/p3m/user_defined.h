#pragma once
/* Standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
/* FDPS headers */
#include "FDPS_c_if.h"

/* Force class */
typedef struct force_pp { //$fdps Force
   //$fdps clear
   double pot;
   fdps_f64vec acc;
} Force_pp;

/* Full Particle Class */
typedef struct fp_nbody { //$fdps FP
   //$fdps copyFromForce force_pp (pot,pot) (acc,acc)
   //$fdps copyFromForcePM acc_pm
   long long int id;
   double mass; //$fdps charge
   double rcut; //$fdps rsearch
   fdps_f64vec pos; //$fdps position
   fdps_f64vec acc;
   double pot;
   fdps_f32vec acc_pm;
   float pot_pm;
} FP_nbody;

/* Essential Particle Class */
typedef struct ep_nbody { //$fdps EPI,EPJ
   //$fdps copyFromFP fp_nbody (id,id) (mass,mass) (rcut,rcut) (pos,pos)
   long long int id;
   double mass; //$fdps charge
   double rcut; //$fdps rsearch
   fdps_f64vec pos; //$fdps position
} EP_nbody;

/* Crystal Parameters class */
typedef struct crystal_parameters {
   int nptcl_per_side;
   fdps_f64vec pos_vertex;
} Crystal_parameters;

/* Prototype declarations */
double S2_pcut(double xi);
double S2_fcut(double xi);
void calc_force_ep_ep(EP_nbody *ep_i,
                      int n_ip,
                      EP_nbody *ep_j,
                      int n_jp,
                      Force_pp *f);
void calc_force_ep_sp(EP_nbody *ep_i,
                      int n_ip,
                      fdps_spj_monopole_cutoff *ep_j,
                      int n_jp,
                      Force_pp *f);
