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
/* User-defined headers */
#include "macro_defs.h"
#include "user_defined.h"
#include "mathematical_constants.h"
#include "physical_constants.h"

/* Structure definitions */
// The following two classes are used in function readTipsyFile,
// which is used to read particle data created by MAGI.
typedef struct magi_tipsy_header {
   double time;
   int nbodies;
   int ndim;
   int nsph;
   int ndark;
   int nstar;
} Magi_tipsy_header;

typedef struct magi_tipsy_particle {
   float mass;
   float pos[3];
   float vel[3];
   float eps;
   int idx;
} Magi_tipsy_particle;

/* Prototype declarations */
void read_tipsy_file(char file_name[], int psys_num);
void galaxy_IC(int psys_num_nbody,
               int psys_num_sph,
               int *bc,
               fdps_f32vec *pos_root_domain_low,
               fdps_f32vec *pos_root_domain_high,
               double *time_dum,
               double *dt_dump,
               double *time_end);
void cold_collapse_test_IC(int psys_num_nbody,
                           int psys_num_sph,
                           int *bc,
                           fdps_f32vec *pos_root_domain_low,
                           fdps_f32vec *pos_root_domain_high,
                           double *time_dum,
                           double *dt_dump,
                           double *time_end);
void Evrard_test_IC(int psys_num_nbody,
                    int psys_num_sph,
                    int *bc,
                    fdps_f32vec *pos_root_domain_low,
                    fdps_f32vec *pos_root_domain_high,
                    double *time_dum,
                    double *dt_dump,
                    double *time_end,
                    int gen_mode);
double get_pos_cell_center(double pos_left_bound,
                           double pos_right_bound,
                           int number_of_cells,
                           int i);
void make_glass_IC(int psys_num_nbody,
                   int psys_num_sph,
                   int *bc,
                   fdps_f32vec *pos_root_domain_low,
                   fdps_f32vec *pos_root_domain_high,
                   double *time_dum,
                   double *dt_dump,
                   double *time_end);

