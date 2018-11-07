/* Standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
/* FDPS headers */
#include "FDPS_c_if.h"
/* user-defined headers */
#include "mathematical_constants.h"
#include "user_defined.h"

void setup_IC(int psys_num,
              double *end_time,
              fdps_f32vec *pos_ll,
              fdps_f32vec *pos_ul) {
    // Get # of MPI processes and rank number
    int nprocs = fdps_get_num_procs();
    int myrank = fdps_get_rank();

    // Set the box size
    pos_ll->x = 0.0;
    pos_ll->y = 0.0;
    pos_ll->z = 0.0;
    pos_ul->x = 1.0;
    pos_ul->y = pos_ul->x / 8.0;
    pos_ul->z = pos_ul->x / 8.0;

    // Make an initial condition at RANK 0
    if (myrank == 0) {
        // Set the left and right states
        const double dens_L = 1.0;
        const double eng_L  = 2.5;
        const double dens_R = 0.5;
        const double eng_R  = 2.5;
        // Set the separation of particle of the left state
        const double dx = 1.0 / 128.0;
        const double dy = dx;
        const double dz = dx;
        // Set the number of local particles
        int nptcl_glb = 0;
        // (1) Left-half
        const int nx_L = 0.5*pos_ul->x/dx;
        const int ny_L = pos_ul->y/dy;
        const int nz_L = pos_ul->z/dz;
        nptcl_glb += nx_L * ny_L * nz_L;
        printf("nptcl_glb(L)   = %d\n",nptcl_glb);
        // (2) Right-half
        const int nx_R = 0.5*pos_ul->x/((dens_L/dens_R)*dx);
        const int ny_R = ny_L;
        const int nz_R = nz_L;
        nptcl_glb += nx_R * ny_R * nz_R;
        printf("nptcl_glb(L+R) = %d\n",nptcl_glb);
        // Place SPH particles
        fdps_set_nptcl_loc(psys_num,nptcl_glb);
        Full_particle *ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
        int id = -1;
        // (1) Left-half
        int i,j,k;
        for (i = 0; i < nx_L; i++) {
            for (j = 0; j < ny_L; j++) {
                for (k = 0; k < nz_L; k++) {
                    id++;
                    ptcl[id].id    = id;
                    ptcl[id].pos.x = dx * i;
                    ptcl[id].pos.y = dy * j;
                    ptcl[id].pos.z = dz * k;
                    ptcl[id].dens  = dens_L;
                    ptcl[id].eng   = eng_L;
                }
            }
        }
        // (2) Right-half
        for (i = 0; i < nx_R; i++) {
            for (j = 0; j < ny_R; j++) {
                for (k = 0; k < nz_R; k++) {
                    id++;
                    ptcl[id].id    = id;
                    ptcl[id].pos.x = 0.5*pos_ul->x + ((dens_L/dens_R)*dx)*i;
                    ptcl[id].pos.y = dy * j;
                    ptcl[id].pos.z = dz * k;
                    ptcl[id].dens  = dens_R;
                    ptcl[id].eng   = eng_R;
                }
            }
        }
        printf("nptcl(L+R) = %d\n",id+1);
        // Set particle mass and smoothing length
        for (i = 0; i < nptcl_glb; i++) {
            ptcl[i].mass = 0.5*(dens_L+dens_R)
                         * (pos_ul->x*pos_ul->y*pos_ul->z)
                         / nptcl_glb;
            ptcl[i].smth = kernel_support_radius * 0.012;
        }
    } else {
        fdps_set_nptcl_loc(psys_num,0);
    }

    // Set the end time
    *end_time = 0.12;

    // Inform to STDOUT
    if (fdps_get_rank() == 0) printf("setup...completed!\n");
    //fdps_finalize();
    //exit(0);

}

double get_timestep(int psys_num) {
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    double dt_loc = 1.0e30;
    int i;
    for (i = 0; i < nptcl_loc; i++) 
        if (ptcl[i].dt < dt_loc)
            dt_loc = ptcl[i].dt;
    return fdps_get_min_value_f64(dt_loc);
}

void initial_kick(int psys_num, double dt) {
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        ptcl[i].vel_half.x = ptcl[i].vel.x + 0.5 * dt * ptcl[i].acc.x;
        ptcl[i].vel_half.y = ptcl[i].vel.y + 0.5 * dt * ptcl[i].acc.y;
        ptcl[i].vel_half.z = ptcl[i].vel.z + 0.5 * dt * ptcl[i].acc.z;
        ptcl[i].eng_half = ptcl[i].eng + 0.5 * dt * ptcl[i].eng_dot;
    }
}

void full_drift(int psys_num, double dt) {
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        ptcl[i].pos.x += dt * ptcl[i].vel_half.x;
        ptcl[i].pos.y += dt * ptcl[i].vel_half.y;
        ptcl[i].pos.z += dt * ptcl[i].vel_half.z;
    }
}

void predict(int psys_num, double dt) {
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        ptcl[i].vel.x += dt * ptcl[i].acc.x;
        ptcl[i].vel.y += dt * ptcl[i].acc.y;
        ptcl[i].vel.z += dt * ptcl[i].acc.z;
        ptcl[i].eng += dt * ptcl[i].eng_dot;
    }
}

void final_kick(int psys_num, double dt) {
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        ptcl[i].vel.x = ptcl[i].vel_half.x + 0.5 * dt * ptcl[i].acc.x;
        ptcl[i].vel.y = ptcl[i].vel_half.y + 0.5 * dt * ptcl[i].acc.y;
        ptcl[i].vel.z = ptcl[i].vel_half.z + 0.5 * dt * ptcl[i].acc.z;
        ptcl[i].eng = ptcl[i].eng_half + 0.5 * dt * ptcl[i].eng_dot;
    }
}

void set_pressure(int psys_num) {
    const double hcr=1.4;
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        ptcl[i].pres = (hcr - 1.0) * ptcl[i].dens * ptcl[i].eng;
        ptcl[i].snds = sqrt(hcr * ptcl[i].pres / ptcl[i].dens);
    }
}

void output(int psys_num, int nstep) {
    int myrank = fdps_get_rank();
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    char filename[64] = {'\0'};
    sprintf(filename,"./result/snap%05d-proc%05d.txt",nstep,myrank);
    FILE *fp;
    if ((fp = fopen(filename,"w")) == NULL) {
        fprintf(stderr,"Cannot open file %s\n",filename);
        exit(EXIT_FAILURE);
    }
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        fprintf(fp,"%ld %15.7e %15.7e %15.7e %15.7e",
                ptcl[i].id,ptcl[i].mass, 
                ptcl[i].pos.x,ptcl[i].pos.y,ptcl[i].pos.z);
        fprintf(fp,"%15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n",
                ptcl[i].vel.x,ptcl[i].vel.y,ptcl[i].vel.z,
                ptcl[i].dens,ptcl[i].eng,ptcl[i].pres);
    }
    fclose(fp);
}

void check_cnsrvd_vars(int psys_num){
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl = fdps_get_psys_cptr(psys_num);
    fdps_f64vec mom_loc;
    mom_loc.x = 0.0; mom_loc.y = 0.0; mom_loc.z = 0.0;
    double eng_loc = 0.0;
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        mom_loc.x += ptcl[i].vel.x * ptcl[i].mass;
        mom_loc.y += ptcl[i].vel.y * ptcl[i].mass;
        mom_loc.z += ptcl[i].vel.z * ptcl[i].mass;
        eng_loc += ptcl[i].mass *(ptcl[i].eng
                                 + 0.5*(ptcl[i].vel.x * ptcl[i].vel.x
                                       +ptcl[i].vel.y * ptcl[i].vel.y
                                       +ptcl[i].vel.z * ptcl[i].vel.z));
    }
    double eng = fdps_get_sum_f64(eng_loc);
    fdps_f64vec mom;
    mom.x = fdps_get_sum_f64(mom_loc.x);
    mom.y = fdps_get_sum_f64(mom_loc.y);
    mom.z = fdps_get_sum_f64(mom_loc.z);
    if (fdps_get_rank() == 0) {
        printf("eng   = %15.7e\n",eng);
        printf("mom.x = %15.7e\n",mom.x);
        printf("mom.y = %15.7e\n",mom.y);
        printf("mom.z = %15.7e\n",mom.z);
    }
}

void c_main() {
    // Initialize some global variables
    setup_math_const();

    // Initialize FDPS
    fdps_initialize();

    // Make an instance of ParticleSystem and initialize it
    int psys_num;
    fdps_create_psys(&psys_num,"full_particle");
    fdps_init_psys(psys_num);

    // Make an initial condition and initialize the particle system
    double end_time;
    fdps_f32vec pos_ll, pos_ul;
    setup_IC(psys_num,&end_time,&pos_ll,&pos_ul);

    // Make an instance of DomainInfo and initialize it
    int dinfo_num;
    fdps_create_dinfo(&dinfo_num);
    float coef_ema = 0.3;
    fdps_init_dinfo(dinfo_num,coef_ema);
    fdps_set_boundary_condition(dinfo_num,FDPS_BC_PERIODIC_XYZ);
    fdps_set_pos_root_domain(dinfo_num,&pos_ll,&pos_ul);

    // Perform domain decomposition and exchange particles
    fdps_decompose_domain_all(dinfo_num,psys_num,-1.0);
    fdps_exchange_particle(psys_num,dinfo_num);

    // Make two tree structures
    int ntot = fdps_get_nptcl_glb(psys_num);
    // tree_dens (used for the density calculation)
    int tree_num_dens;
    fdps_create_tree(&tree_num_dens,
                     "Short,force_dens,essential_particle,essential_particle,Gather");
    float theta = 0.5;
    int n_leaf_limit = 8;
    int n_group_limit = 64;
    fdps_init_tree(tree_num_dens,3*ntot,theta,n_leaf_limit,n_group_limit);

    // tree_hydro (used for the force calculation)
    int tree_num_hydro;
    fdps_create_tree(&tree_num_hydro,
                     "Short,force_hydro,essential_particle,essential_particle,Symmetry");
    fdps_init_tree(tree_num_hydro,3*ntot,theta,n_leaf_limit,n_group_limit);

    // Compute density, pressure, acceleration due to pressure gradient
    fdps_calc_force_all_and_write_back(tree_num_dens, 
                                       calc_density,
                                       NULL,
                                       psys_num,
                                       dinfo_num,
                                       true,
                                       FDPS_MAKE_LIST);
    set_pressure(psys_num);
    fdps_calc_force_all_and_write_back(tree_num_hydro,
                                       calc_hydro_force,
                                       NULL,
                                       psys_num,
                                       dinfo_num,
                                       true,
                                       FDPS_MAKE_LIST);
    // Get timestep
    double dt = get_timestep(psys_num);

    // Main loop for time integration
    int nstep = 0; double time = 0.0;
    for (;;) {
        // Leap frog: Initial Kick & Full Drift
        initial_kick(psys_num,dt);
        full_drift(psys_num,dt);

        // Adjust the positions of the SPH particles that run over
        // the computational boundaries.
        fdps_adjust_pos_into_root_domain(psys_num,dinfo_num);

        // Leap frog: Predict
        predict(psys_num,dt);

        // Perform domain decomposition and exchange particles again
        fdps_decompose_domain_all(dinfo_num,psys_num,-1.0);
        fdps_exchange_particle(psys_num,dinfo_num);

        // Compute density, pressure, acceleration due to pressure gradient
        fdps_calc_force_all_and_write_back(tree_num_dens,
                                           calc_density,
                                           NULL,
                                           psys_num,
                                           dinfo_num,
                                           true,
                                           FDPS_MAKE_LIST);
        set_pressure(psys_num);
        fdps_calc_force_all_and_write_back(tree_num_hydro,
                                           calc_hydro_force,
                                           NULL,
                                           psys_num,
                                           dinfo_num,
                                           true,
                                           FDPS_MAKE_LIST);

        // Get a new timestep
        dt = get_timestep(psys_num);

        // Leap frog: Final Kick
        final_kick(psys_num,dt);

        // Output result files
        int output_interval = 10;
        if (nstep % output_interval == 0) {
           output(psys_num,nstep);
           check_cnsrvd_vars(psys_num);
        }

        // Output information to STDOUT
        if (fdps_get_rank() == 0) {
           printf("================================\n");
           printf("time  = %15.7e\n",time);
           printf("nstep = %d\n",nstep);
           printf("================================\n");
        }

        // Termination condition
        if (time >= end_time) break;

        // Update time & step
        time  += dt;
        nstep++;
    }
    fdps_finalize();
}
