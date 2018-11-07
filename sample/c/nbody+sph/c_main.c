/* Standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
/* FDPS headers */
#include "FDPS_c_if.h"
/* User-defined headers */
#include "macro_defs.h"
#include "user_defined.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "leapfrog.h"
#include "ic.h"

void setup_IC(int psys_num_nbody,
              int psys_num_sph,
              int dinfo_num,
              double *time_dump,
              double *dt_dump,
              double *time_end) {
    int i, bc;
    fdps_f32vec pos_root_domain_low, pos_root_domain_high;
    if (fdps_get_rank() == 0) {
#if (INITIAL_CONDITION == 0)
        galaxy_IC(psys_num_nbody,psys_num_sph,
                  &bc,
                  &pos_root_domain_low,
                  &pos_root_domain_high,
                  time_dump,dt_dump,time_end);
#elif (INITIAL_CONDITION == 1)
        cold_collapse_test_IC(psys_num_nbody,psys_num_sph,
                              &bc,
                              &pos_root_domain_low,
                              &pos_root_domain_high,
                              time_dump,dt_dump,time_end);
#elif (INITIAL_CONDITION == 2)
        const int gen_mode = 1; // 0 or 1
        Evrard_test_IC(psys_num_nbody,psys_num_sph,
                       &bc,
                       &pos_root_domain_low,
                       &pos_root_domain_high,
                       time_dump,dt_dump,time_end,
                       gen_mode);
#elif (INITIAL_CONDITION == 3)
        make_glass_IC(psys_num_nbody,psys_num_sph,
                      &bc,
                      &pos_root_domain_low,
                      &pos_root_domain_high,
                      time_dump,dt_dump,time_end);
#else
#error Invalid IC number is specified.
#endif
        // Check the initial distribution
        char file_name[64] = {'\0'};
        strcpy(file_name,"IC.txt");
        FILE *fp;
        if ((fp = fopen(file_name,"w")) == NULL) {
            fprintf(stderr,"Cannot open file %s.\n",file_name);
            exit(EXIT_FAILURE);
        }
        int nptcl_loc = fdps_get_nptcl_loc(psys_num_nbody);
        FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
        for (i = 0; i < nptcl_loc; i++) {
            fprintf(fp,"%15.7e %15.7e %15.7e\n",
                    ptcl_nbody[i].pos.x,ptcl_nbody[i].pos.y,ptcl_nbody[i].pos.z);
        }
        fprintf(fp,"\n\n");
        nptcl_loc = fdps_get_nptcl_loc(psys_num_sph);
        FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
        for (i = 0; i < nptcl_loc; i++) {
            fprintf(fp,"%15.7e %15.7e %15.7e\n",
                    ptcl_sph[i].pos.x,ptcl_sph[i].pos.y,ptcl_sph[i].pos.z);
        }
        fclose(fp);
    } else  {
        fdps_set_nptcl_loc(psys_num_nbody,0);
        fdps_set_nptcl_loc(psys_num_sph,0);
    }

    // Broadcast from RANK 0 
    fdps_broadcast_s32(&bc, 1, 0);
    fdps_broadcast_f32(&pos_root_domain_low.x, 1, 0);
    fdps_broadcast_f32(&pos_root_domain_low.y, 1, 0);
    fdps_broadcast_f32(&pos_root_domain_low.z, 1, 0);
    fdps_broadcast_f32(&pos_root_domain_high.x, 1, 0);
    fdps_broadcast_f32(&pos_root_domain_high.y, 1, 0);
    fdps_broadcast_f32(&pos_root_domain_high.z, 1, 0);
    fdps_broadcast_f64(&eps_grav, 1, 0);
    fdps_broadcast_f64(time_dump, 1, 0);
    fdps_broadcast_f64(dt_dump, 1, 0);
    fdps_broadcast_f64(time_end, 1, 0);
    fdps_broadcast_f64(&dt_max, 1, 0);
   
    // Set the boundary condition and the size of the computational domain if needed.
    fdps_set_boundary_condition(dinfo_num,bc);
    if (bc != FDPS_BC_OPEN) {
        fdps_set_pos_root_domain(dinfo_num,
                                 &pos_root_domain_low,
                                 &pos_root_domain_high);
    }

    // Compute the average mass of SPH particles
    int nptcl_glb = fdps_get_nptcl_glb(psys_num_sph);
    int nptcl_loc = fdps_get_nptcl_loc(psys_num_sph);
    FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
    double m_sum_loc = 0.0, m_sum = 0.0;
    for (i = 0; i < nptcl_loc; i++) m_sum_loc += ptcl_sph[i].mass;
    m_sum = fdps_get_sum_f64(m_sum_loc);
    mass_avg = m_sum / nptcl_glb;

    // Inform to STDOUT
    if (fdps_get_rank() == 0)
        printf("setup_IC() is completed.\n");
    //fdps_finalize();
    //exit(EXIT_SUCCESS);

}

double get_timestep(int psys_num_nbody,
                    int psys_num_sph) {
    double dt_loc = DBL_MAX;
    if (dt_max > 0.0) dt_loc = dt_max;

    int i,nptcl_loc;
    // Timescale for N-body system
#if defined(ENABLE_GRAVITY_INTERACT)
    nptcl_loc = fdps_get_nptcl_loc(psys_num_nbody);
    FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
    for (i = 0; i < nptcl_loc; i++) {
        double acc = sqrt(ptcl_nbody[i].acc.x * ptcl_nbody[i].acc.x
                         +ptcl_nbody[i].acc.y * ptcl_nbody[i].acc.y
                         +ptcl_nbody[i].acc.z * ptcl_nbody[i].acc.z);
        if (acc > 0.0) { 
            double tmp = CFL_dyn * sqrt(eps_grav / acc);
            if (tmp < dt_loc) dt_loc = tmp;
        }
    }
#endif

    // Timescale for SPH system
    nptcl_loc = fdps_get_nptcl_loc(psys_num_sph);
    FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
    for (i = 0; i < nptcl_loc; i++) {
#if defined(ENABLE_GRAVITY_INTERACT)
        fdps_f64vec acc_tot;
        acc_tot.x = ptcl_sph[i].acc_grav.x + ptcl_sph[i].acc_hydro.x;
        acc_tot.y = ptcl_sph[i].acc_grav.y + ptcl_sph[i].acc_hydro.y;
        acc_tot.z = ptcl_sph[i].acc_grav.z + ptcl_sph[i].acc_hydro.z;
        double acc = sqrt(acc_tot.x * acc_tot.x
                         +acc_tot.y * acc_tot.y
                         +acc_tot.z * acc_tot.z);
        if (acc > 0.0) {
            double tmp = CFL_dyn * sqrt(eps_grav / acc);
            if (tmp < dt_loc) dt_loc = tmp;
        }
#endif
#if defined(ENABLE_HYDRO_INTERACT)
        if (ptcl_sph[i].dt < dt_loc) dt_loc = ptcl_sph[i].dt;
#endif
    }
    // Reduction
    return fdps_get_min_value_f64(dt_loc);

}

void calc_density_wrapper(int psys_num,
                          int dinfo_num,
                          int tree_num) {
#if defined(ENABLE_VARIABLE_SMOOTHING_LENGTH)
   int nptcl_loc = fdps_get_nptcl_loc(psys_num);
   int nptcl_glb = fdps_get_nptcl_glb(psys_num);
   FP_sph *ptcl = (FP_sph *) fdps_get_psys_cptr(psys_num);
   // Determine the density and the smoothing length
   // so that Eq.(6) in Springel (2005) holds within a specified accuracy.
   for (;;) { 
       // Increase smoothing length 
       int i;
       for (i = 0; i < nptcl_loc; i++) ptcl[i].smth *= SCF_smth;
       // Compute density, etc.
       fdps_calc_force_all_and_write_back(tree_num,
                                          calc_density,
                                          NULL,
                                          psys_num,
                                          dinfo_num,
                                          true,
                                          FDPS_MAKE_LIST);
       // Check convergence
       int n_compl_loc = 0;
       for (i = 0; i < nptcl_loc; i++)
           if (ptcl[i].flag == 1) n_compl_loc++;
       int n_compl = fdps_get_sum_s32(n_compl_loc);
       if (n_compl == nptcl_glb) break;
   }
#else
   fdps_calc_force_all_and_write_back(tree_num,
                                      calc_density,
                                      NULL,
                                      psys_num,
                                      dinfo_num,
                                      true,
                                      FDPS_MAKE_LIST);
#endif
}

void set_pressure(int psys_num) {
    int i, nptcl_loc = fdps_get_nptcl_loc(psys_num);
    FP_sph *ptcl = (FP_sph *) fdps_get_psys_cptr(psys_num);
    for (i = 0; i < nptcl_loc; i++) {
#if defined(ISOTHERMAL_EOS)
        // In this case, eng = const.
        ptcl[i].pres = (specific_heat_ratio - 1.0) * ptcl[i].dens * ptcl[i].eng;
        ptcl[i].ent  = ptcl[i].pres / pow(ptcl[i].dens,specific_heat_ratio);
#else
#if defined(USE_ENTROPY)
        ptcl[i].pres = ptcl[i].ent * pow(ptcl[i].dens,specific_heat_ratio);
        ptcl[i].eng  = ptcl[i].pres / ((specific_heat_ratio - 1.0) * ptcl[i].dens);
#else
        ptcl[i].pres = (specific_heat_ratio - 1.0) * ptcl[i].dens * ptcl[i].eng;
        ptcl[i].ent  = ptcl[i].pres / pow(ptcl[i].dens,specific_heat_ratio);
#endif
#endif
        ptcl[i].snds = sqrt(specific_heat_ratio * ptcl[i].pres / ptcl[i].dens);
#if defined(USE_BALSARA_SWITCH)
        double rotv_abs = sqrt(ptcl[i].rotv.x * ptcl[i].rotv.x
                              +ptcl[i].rotv.y * ptcl[i].rotv.y
                              +ptcl[i].rotv.z * ptcl[i].rotv.z);
        ptcl[i].balsw = fabs(ptcl[i].divv) / (fabs(ptcl[i].divv) + rotv_abs
                                              +1.0e-4 * ptcl[i].snds / ptcl[i].smth);
#else
        ptcl[i].balsw = 1.0;
#endif
    }
}

void set_entropy(int psys_num) {
   int i, nptcl_loc = fdps_get_nptcl_loc(psys_num);
   FP_sph *ptcl = (FP_sph *) fdps_get_psys_cptr(psys_num);
   for (i = 0; i < nptcl_loc; i++)
       ptcl[i].ent = (specific_heat_ratio - 1.0) * ptcl[i].eng
                   / pow(ptcl[i].dens,specific_heat_ratio - 1.0);

}

void set_circular_velocity(int psys_num) {
    int i, nptcl_loc = fdps_get_nptcl_loc(psys_num);
    FP_sph *ptcl = (FP_sph *) fdps_get_psys_cptr(psys_num);
    for (i = 0; i < nptcl_loc; i++) {
        fdps_f64vec acc;
        acc.x = ptcl[i].acc_grav.x + ptcl[i].acc_hydro.x;
        acc.y = ptcl[i].acc_grav.y + ptcl[i].acc_hydro.y;
        acc.z = ptcl[i].acc_grav.z + ptcl[i].acc_hydro.z;
        double r = sqrt(ptcl[i].pos.x * ptcl[i].pos.x
                       +ptcl[i].pos.y * ptcl[i].pos.y
                       +ptcl[i].pos.z * ptcl[i].pos.z);
        double rcyl = sqrt(ptcl[i].pos.x * ptcl[i].pos.x
                          +ptcl[i].pos.y * ptcl[i].pos.y);
        double phi = atan2(ptcl[i].pos.y, ptcl[i].pos.x);
        double theta = atan2(rcyl, ptcl[i].pos.z);
        fdps_f64vec base_vect_r;
        base_vect_r.x = sin(theta) * cos(phi);
        base_vect_r.y = sin(theta) * sin(phi);
        base_vect_r.z = cos(theta);
        double vel_circ = sqrt(r * fabs(acc.x * base_vect_r.x
                                       +acc.y * base_vect_r.y
                                       +acc.z * base_vect_r.z));
        ptcl[i].vel.x = - vel_circ * sin(phi);
        ptcl[i].vel.y =   vel_circ * cos(phi);
        ptcl[i].vel.z = 0.0;
    }
#if 0
   // [for debug] 
   // Output the velocity field on the x-y plane
   char file_name[64];
   sprintf(file_name,"velc_fld%05d.txt",fdps_get_rank());
   FILE *fp;
   if ((fp = fopen(file_name,"w")) == NULL) {
       fprintf(stderr,"Cannot open file %s.\n",file_name);
       exit(EXIT_FAILURE);
   }
   for (i = 0; i < nptcl_loc; i++) {
       fprintf(fp,"%15.7e %15.7e %15.7e %15.7e\n",
               ptcl[i].pos.x,ptcl[i].pos.y,
               ptcl[i].vel.x,ptcl[i].pos.y);
   }
   fclose(fp);
   // Output the rotation curve
   memset(file_name, '\0', strlen(file_name));
   sprintf(file_name,"rot_curve%05d.txt",fdps_get_rank());
   if ((fp = fopen(file_name,"w")) == NULL) {
       fprintf(stderr,"Cannot open file %s.\n",file_name);
       exit(EXIT_FAILURE);
   }
   for (i = 0; i < nptcl_loc; i++) {
       fdps_f64vec pos = ptcl[i].pos;
       double r = sqrt(pos.x * pos.x 
                      +pos.y * pos.y
                      +pos.z * pos.z);
       double rcyl = sqrt(pos.x * pos.x 
                         +pos.y * pos.y);
       fdps_f64vec vel = ptcl[i].vel;
       double vel_circ = sqrt(vel.x * vel.x 
                             +vel.y * vel.y
                             +vel.z * vel.z);
       double T_rot = 2.0 * pi * r / vel_circ;
       fprintf(fp,"%15.7e %15.7e %15.7e\n",
               rcyl,vel_circ,T_rot);
   }
   fclose(fp);
   //call fdps_finalize();
   //exit(EXIT_SUCCESS);
#endif

}

void output(int psys_num_nbody,int psys_num_sph) {
    static int ndump = 1;
    int myrank = fdps_get_rank();

    // Output N-body data
    int i,nptcl_loc = fdps_get_nptcl_loc(psys_num_nbody);
    FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
    char file_name[64] = {'\0'};
    sprintf(file_name,"result/nbody%05d-proc%05d.dat",ndump,myrank);
    FILE *fp;
    if ((fp = fopen(file_name,"w")) == NULL) {
        fprintf(fp,"Cannot open file %s.\n",file_name);
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < nptcl_loc; i++) {
        fprintf(fp,"%ld %15.7e %15.7e %15.7e %15.7e ",
                ptcl_nbody[i].id,ptcl_nbody[i].mass,
                ptcl_nbody[i].pos.x,ptcl_nbody[i].pos.y,ptcl_nbody[i].pos.z);
        fprintf(fp,"%15.7e %15.7e %15.7e\n",
                ptcl_nbody[i].vel.x,ptcl_nbody[i].vel.y,ptcl_nbody[i].vel.z);
    }
    fclose(fp);

    // Output SPH data
    nptcl_loc = fdps_get_nptcl_loc(psys_num_sph);
    FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
    memset(file_name,'\0',strlen(file_name));
    sprintf(file_name,"result/sph%05d-proc%05d.dat",ndump,myrank);
    if ((fp = fopen(file_name,"w")) == NULL) {
        fprintf(fp,"Cannot open file %s.\n",file_name);
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < nptcl_loc; i++) {
        fprintf(fp,"%ld %15.7e %15.7e %15.7e %15.7e ",
                ptcl_sph[i].id,ptcl_sph[i].mass,
                ptcl_sph[i].pos.x,ptcl_sph[i].pos.y,ptcl_sph[i].pos.z);
        fprintf(fp,"%15.7e %15.7e %15.7e ",
                ptcl_sph[i].vel.x,ptcl_sph[i].vel.y,ptcl_sph[i].vel.z);
        fprintf(fp,"%15.7e %15.7e %15.7e\n",
                ptcl_sph[i].dens,ptcl_sph[i].eng,ptcl_sph[i].pres);
    }
    fclose(fp);

    // Update ndump
    ndump++;
}

void check_cnsrvd_vars(int psys_num_nbody,
                       int psys_num_sph,
                       double time) {
   static bool is_initialized = false;
   static double emech_ini,etot_ini;

   double ekin_loc = 0.0;
   double epot_loc = 0.0;
   double eth_loc  = 0.0;
   fdps_f64vec mom_loc;
   mom_loc.x = 0.0; mom_loc.y = 0.0; mom_loc.z = 0.0;
   // Sum over N-body system
   int i, nptcl_loc = fdps_get_nptcl_loc(psys_num_nbody);
   FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
   for (i = 0; i < nptcl_loc; i++) {
      ekin_loc += 0.5 * ptcl_nbody[i].mass * (ptcl_nbody[i].vel.x * ptcl_nbody[i].vel.x
                                             +ptcl_nbody[i].vel.y * ptcl_nbody[i].vel.y
                                             +ptcl_nbody[i].vel.z * ptcl_nbody[i].vel.z);
      epot_loc += 0.5 * ptcl_nbody[i].mass * (ptcl_nbody[i].pot
                                             +ptcl_nbody[i].mass / eps_grav);
      mom_loc.x += ptcl_nbody[i].mass * ptcl_nbody[i].vel.x;
      mom_loc.y += ptcl_nbody[i].mass * ptcl_nbody[i].vel.y;
      mom_loc.z += ptcl_nbody[i].mass * ptcl_nbody[i].vel.z;
   }
   // Sum over SPH system
   nptcl_loc = fdps_get_nptcl_loc(psys_num_sph);
   FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
   for (i = 0; i < nptcl_loc; i++) {
      ekin_loc += 0.5 * ptcl_sph[i].mass * (ptcl_sph[i].vel.x * ptcl_sph[i].vel.x
                                           +ptcl_sph[i].vel.y * ptcl_sph[i].vel.y
                                           +ptcl_sph[i].vel.z * ptcl_sph[i].vel.z);
      epot_loc += 0.5 * ptcl_sph[i].mass * (ptcl_sph[i].pot_grav
                                           +ptcl_sph[i].mass / eps_grav);
      eth_loc  += ptcl_sph[i].mass * ptcl_sph[i].eng;
      mom_loc.x += ptcl_sph[i].mass * ptcl_sph[i].vel.x;
      mom_loc.y += ptcl_sph[i].mass * ptcl_sph[i].vel.y;
      mom_loc.z += ptcl_sph[i].mass * ptcl_sph[i].vel.z;
   }

   // Reduction 
   double ekin  = fdps_get_sum_f64(ekin_loc);
   double epot  = fdps_get_sum_f64(epot_loc);
   double eth   = fdps_get_sum_f64(eth_loc);
   fdps_f64vec mom;
   mom.x = fdps_get_sum_f64(mom_loc.x);
   mom.y = fdps_get_sum_f64(mom_loc.y);
   mom.z = fdps_get_sum_f64(mom_loc.z);

   if (is_initialized == false){
       emech_ini = ekin + epot;
       etot_ini  = ekin + epot + eth;
       is_initialized = true;
   }

   // Output
   if (fdps_get_rank() == 0) {
       double emech = ekin + epot;
       double etot  = ekin + epot + eth;
       double relerr_mech = fabs((emech - emech_ini)/emech_ini);
       double relerr_tot  = fabs((etot  - etot_ini)/etot_ini);
       printf("-------------------------\n");
       printf("E_kin   = %15.7e\n",ekin);
       printf("E_pot   = %15.7e\n",epot);
       printf("E_th    = %15.7e\n",eth);
       printf("E_mech  = %15.7e (time = %15.7e, rel.err = %15.7e)\n",emech,time,relerr_mech);
       printf("E_tot   = %15.7e (time = %15.7e, rel.err = %15.7e)\n",etot,time,relerr_tot);
       printf("Mom (x) = %15.7e\n",mom.x);
       printf("Mom (y) = %15.7e\n",mom.y);
       printf("Mom (z) = %15.7e\n",mom.z);
       printf("-------------------------\n");
   }

}

void check_dens_fluc(int psys_num) {
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    int nptcl_glb = fdps_get_nptcl_glb(psys_num);
    FP_sph *ptcl = (FP_sph *) fdps_get_psys_cptr(psys_num);
    // Compute the average of density
    double tmp = 0.0;
    int i;
    for (i = 0; i < nptcl_loc; i++) tmp += ptcl[i].dens;
    double dens_avg = fdps_get_sum_f64(tmp) / nptcl_glb;
    // Compute the dispersion of density
    tmp = 0.0;
    for (i = 0; i < nptcl_loc; i++)
        tmp += pow(ptcl[i].dens - dens_avg,2.0);
    double dens_disp = sqrt(fdps_get_sum_f64(tmp)/nptcl_glb);
    // Output the status 
    double fluc_str = dens_disp / dens_avg;
    if (fdps_get_rank() == 0) {
        printf("###########################\n");
        printf("avg.       = %15.7e\n",dens_avg);
        printf("disp.      = %15.7e\n",dens_disp);
        printf("disp./avg. = %15.7e\n",fluc_str);
        printf("###########################\n");
    }
    // Output data and end the simulation if the fluctuation is small
    const double eps=3.05e-3;
    if (fluc_str < eps) {
        if (fdps_get_rank() == 0) {
           char file_name[64] = {'\0'};
           sprintf(file_name,"result/glass_data_header.dat");
           FILE *fp;
           if ((fp = fopen(file_name,"wb")) == NULL) {
               fprintf(stderr,"Cannot open file %s.\n",file_name);
               exit(EXIT_FAILURE);
           }
           fwrite(&nptcl_glb,sizeof(int),1,fp);
           fclose(fp);
        }
        char file_name[64] = {'\0'};
        sprintf(file_name,"result/glass_data%05d.dat",fdps_get_rank());
        FILE *fp;
        if ((fp = fopen(file_name,"wb")) == NULL) {
            fprintf(stderr,"Cannot open file %s.\n",file_name);
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < nptcl_loc; i++)
           fwrite(&ptcl[i].pos,sizeof(fdps_f64vec),1,fp);
        fclose(fp);
        if (fdps_get_rank() == 0) {
            printf("A glass-like distribution is obtained.\n");
            printf("The particle position data is output as files %s, etc.\n",
                   file_name);
        }
        fdps_finalize();
        exit(EXIT_SUCCESS);
    }
}

//-----------------------------------------------------------------------
////////////////////// < M A I N   F U N C T I O N > ////////////////////
//-----------------------------------------------------------------------
void c_main() {
    int i;

    // Initialize some global variables
    setup_math_const();
    setup_phys_const();

    // Initialize FDPS
    fdps_initialize();

    // Make instances of ParticleSystems and initialize them
    int psys_num_nbody;
    fdps_create_psys(&psys_num_nbody,"fp_nbody");
    fdps_init_psys(psys_num_nbody);
    int psys_num_sph;
    fdps_create_psys(&psys_num_sph,"fp_sph");
    fdps_init_psys(psys_num_sph);

    // Make an instance of DomainInfo and initialize it
    int dinfo_num;
    fdps_create_dinfo(&dinfo_num);
    const float coef_ema=0.3;
    fdps_init_dinfo(dinfo_num,coef_ema);

    // Make an initial condition and initialize the particle system
    double time_dump, dt_dump, time_end;
    setup_IC(psys_num_nbody, psys_num_sph, dinfo_num,
             &time_dump, &dt_dump, &time_end);

    // Perform domain decomposition
    fdps_collect_sample_particle(dinfo_num, psys_num_nbody, true, -1.0);
    fdps_collect_sample_particle(dinfo_num, psys_num_sph, false, -1.0);
    fdps_decompose_domain(dinfo_num);

    // Perform particle exchange
    fdps_exchange_particle(psys_num_nbody,dinfo_num);
    fdps_exchange_particle(psys_num_sph,dinfo_num);

    // Make three tree structures
    int nptcl_loc_sph = 1;
    if (fdps_get_nptcl_loc(psys_num_sph) > 1)
        nptcl_loc_sph = fdps_get_nptcl_loc(psys_num_sph);
    int nptcl_loc_nbody = fdps_get_nptcl_loc(psys_num_nbody);
    int nptcl_loc_all   = nptcl_loc_nbody + nptcl_loc_sph;
    // tree for gravity calculation
    int tree_num_grav;
    fdps_create_tree(&tree_num_grav,
                     "Long,force_grav,ep_grav,ep_grav,Monopole");
    const float theta=0.5;
    const int n_leaf_limit=8, n_group_limit=64;
    fdps_init_tree(tree_num_grav, 3*nptcl_loc_all, theta,
                   n_leaf_limit, n_group_limit);
    // tree for the density calculation
    int tree_num_dens;
    fdps_create_tree(&tree_num_dens,
                     "Short,force_dens,ep_hydro,ep_hydro,Gather");
    fdps_init_tree(tree_num_dens, 3*nptcl_loc_sph, theta,
                   n_leaf_limit, n_group_limit);
    // tree for the hydrodynamic force calculation
    int tree_num_hydro;
    fdps_create_tree(&tree_num_hydro,
                     "Short,force_hydro,ep_hydro,ep_hydro,Symmetry");
    fdps_init_tree(tree_num_hydro, 3*nptcl_loc_sph, theta,
                   n_leaf_limit, n_group_limit);

    // Perform force calculations
#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_open();
    g5_set_eps_to_all(eps_grav);
#endif
    // Gravity calculation
    double t_start = fdps_get_wtime();
#if defined(ENABLE_GRAVITY_INTERACT)
    fdps_set_particle_local_tree(tree_num_grav, psys_num_nbody, true);
    fdps_set_particle_local_tree(tree_num_grav, psys_num_sph, false);
    fdps_calc_force_making_tree(tree_num_grav,
                                calc_gravity_ep_ep,
                                calc_gravity_ep_sp,
                                dinfo_num,
                                true);
    nptcl_loc_nbody = fdps_get_nptcl_loc(psys_num_nbody);
    FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
    for (i = 0; i < nptcl_loc_nbody; i++) {
        Force_grav f_grav;
        void *pforce = (void *) &f_grav;
        fdps_get_force(tree_num_grav, i, pforce);
        ptcl_nbody[i].acc.x = f_grav.acc.x;
        ptcl_nbody[i].acc.y = f_grav.acc.y;
        ptcl_nbody[i].acc.z = f_grav.acc.z;
        ptcl_nbody[i].pot   = f_grav.pot;
    }
    int offset = nptcl_loc_nbody;
    nptcl_loc_sph = fdps_get_nptcl_loc(psys_num_sph);
    FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
    for (i = 0; i < nptcl_loc_sph; i++) {
        Force_grav f_grav;
        fdps_get_force(tree_num_grav, i + offset, (void *)&f_grav);
        ptcl_sph[i].acc_grav.x = f_grav.acc.x;
        ptcl_sph[i].acc_grav.y = f_grav.acc.y;
        ptcl_sph[i].acc_grav.z = f_grav.acc.z;
        ptcl_sph[i].pot_grav   = f_grav.pot;
    }
#endif
    double t_grav = fdps_get_wtime() - t_start;
    // SPH calculations
    t_start = fdps_get_wtime();
#if defined(ENABLE_HYDRO_INTERACT)
    calc_density_wrapper(psys_num_sph, dinfo_num, tree_num_dens);
    set_entropy(psys_num_sph);
    set_pressure(psys_num_sph);
    fdps_calc_force_all_and_write_back(tree_num_hydro,
                                       calc_hydro_force,
                                       NULL,
                                       psys_num_sph,
                                       dinfo_num,
                                       true,
                                       FDPS_MAKE_LIST);
#endif
    double t_hydro = fdps_get_wtime() - t_start;

    // Set the initial velocity of gas particle
#if defined(SET_CIRCULAR_VELOCITY) 
    set_circular_velocity(psys_num_sph);
#endif

    // Get timestep
    double dt = get_timestep(psys_num_nbody, psys_num_sph);

    // Check the conserved variables
    check_cnsrvd_vars(psys_num_nbody,psys_num_sph,0.0);

    // Main loop for time integration
    int nstep = 1; double time = 0.0;
    for (;;) { 
        if (fdps_get_rank() == 0) {
            printf("nstep = %d, dt = %15.7e, time = %15.7e, time_end = %15.7e\n",
                   nstep,dt,time,time_end);
        }
        double t_start_1step = fdps_get_wtime();

        // Leap frog: Initial Kick & Full Drift
        initial_kick(psys_num_nbody,psys_num_sph,dt);
        full_drift(psys_num_nbody,psys_num_sph,dt);
        if (fdps_get_boundary_condition(dinfo_num) != FDPS_BC_OPEN) {
            fdps_adjust_pos_into_root_domain(psys_num_nbody,dinfo_num);
            fdps_adjust_pos_into_root_domain(psys_num_sph,dinfo_num);
        }

        // Leap frog: Predict
        predict(psys_num_sph,dt);

        // Perform domain decomposition again
        fdps_collect_sample_particle(dinfo_num, psys_num_nbody, true, -1.0);
        fdps_collect_sample_particle(dinfo_num, psys_num_sph, false, -1.0);
        fdps_decompose_domain(dinfo_num);

        // Exchange the particles between the (MPI) processes
        fdps_exchange_particle(psys_num_nbody,dinfo_num);
        fdps_exchange_particle(psys_num_sph,dinfo_num);

        // Perform force calculations
        // Gravity calculation
        t_start = fdps_get_wtime();
#if defined(ENABLE_GRAVITY_INTERACT)
        fdps_set_particle_local_tree(tree_num_grav, psys_num_nbody, true);
        fdps_set_particle_local_tree(tree_num_grav, psys_num_sph, false);
        fdps_calc_force_making_tree(tree_num_grav,
                                    calc_gravity_ep_ep,
                                    calc_gravity_ep_sp,
                                    dinfo_num,
                                    true);
        nptcl_loc_nbody = fdps_get_nptcl_loc(psys_num_nbody);
        FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
        for (i = 0; i < nptcl_loc_nbody; i++) {
            Force_grav f_grav;
            fdps_get_force(tree_num_grav, i, (void *)&f_grav);
            ptcl_nbody[i].acc.x = f_grav.acc.x;
            ptcl_nbody[i].acc.y = f_grav.acc.y;
            ptcl_nbody[i].acc.z = f_grav.acc.z;
            ptcl_nbody[i].pot   = f_grav.pot;
        }
        offset = nptcl_loc_nbody;
        nptcl_loc_sph = fdps_get_nptcl_loc(psys_num_sph);
        FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
        for (i = 0; i < nptcl_loc_sph; i++) {
            Force_grav f_grav;
            fdps_get_force(tree_num_grav, i + offset, (void *)&f_grav);
            ptcl_sph[i].acc_grav.x = f_grav.acc.x;
            ptcl_sph[i].acc_grav.y = f_grav.acc.y;
            ptcl_sph[i].acc_grav.z = f_grav.acc.z;
            ptcl_sph[i].pot_grav   = f_grav.pot;
        }
#endif
        t_grav = fdps_get_wtime() - t_start;
        // SPH calculations
        t_start = fdps_get_wtime();
#if defined(ENABLE_HYDRO_INTERACT)
        calc_density_wrapper(psys_num_sph, dinfo_num, tree_num_dens);
        set_pressure(psys_num_sph);
        fdps_calc_force_all_and_write_back(tree_num_hydro,
                                           calc_hydro_force,
                                           NULL,
                                           psys_num_sph,
                                           dinfo_num,
                                           true,
                                           FDPS_MAKE_LIST);
#endif
        t_hydro = fdps_get_wtime() - t_start;

        // Get a new timestep
        dt = get_timestep(psys_num_nbody,psys_num_sph);

        // Leap frog: Final Kick
        final_kick(psys_num_nbody,psys_num_sph,dt);

        // Get the elapsed time for this step
        double t_1step = fdps_get_wtime() - t_start_1step;
        double t_1step_max = fdps_get_max_value_f64(t_1step);
        if (fdps_get_rank() == 0)
            printf("t_1step_max = %15.7e\n",t_1step_max);

        // Output result files
        if (time > time_dump) {
            output(psys_num_nbody,psys_num_sph);
            if (fdps_get_rank() == 0) {
                printf("============================================\n");
                printf("output a file at time = %15.7e\n",time);
                printf("============================================\n");
            }
            time_dump += dt_dump;
        }

        // Check the conserved variables
        check_cnsrvd_vars(psys_num_nbody,psys_num_sph,time);

        // Check the amplitude of density fluctuation
#if defined(CHECK_DENSITY_FLUCTUATION)
        if (nstep % 100 == 0) check_dens_fluc(psys_num_sph);
#endif

        // Termination condition
        if (time >= time_end) break;

        // Update time & step
        time += dt;
        nstep++;

    }

#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_close();
#endif

    // Finalize FDPS
    fdps_finalize();

}

