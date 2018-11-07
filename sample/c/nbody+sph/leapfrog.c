#include "leapfrog.h"

void initial_kick(int psys_num_nbody,
                  int psys_num_sph,
                  double dt) {
    // Update N-body system
    int nptcl_loc = fdps_get_nptcl_loc(psys_num_nbody);
    FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        ptcl_nbody[i].vel.x += 0.5 * dt * ptcl_nbody[i].acc.x;
        ptcl_nbody[i].vel.y += 0.5 * dt * ptcl_nbody[i].acc.y;
        ptcl_nbody[i].vel.z += 0.5 * dt * ptcl_nbody[i].acc.z;
    }

    // Update SPH system
    nptcl_loc = fdps_get_nptcl_loc(psys_num_sph);
    FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
    for (i = 0; i < nptcl_loc; i++) {
        fdps_f64vec acc_tot;
        acc_tot.x = ptcl_sph[i].acc_grav.x + ptcl_sph[i].acc_hydro.x;
        acc_tot.y = ptcl_sph[i].acc_grav.y + ptcl_sph[i].acc_hydro.y;
        acc_tot.z = ptcl_sph[i].acc_grav.z + ptcl_sph[i].acc_hydro.z;
        ptcl_sph[i].vel_half.x = ptcl_sph[i].vel.x + 0.5 * dt * acc_tot.x;
        ptcl_sph[i].vel_half.y = ptcl_sph[i].vel.y + 0.5 * dt * acc_tot.y;
        ptcl_sph[i].vel_half.z = ptcl_sph[i].vel.z + 0.5 * dt * acc_tot.z;
#if !defined(ISOTHERMAL_EOS)
        ptcl_sph[i].eng_half = ptcl_sph[i].eng + 0.5 * dt * ptcl_sph[i].eng_dot;
        ptcl_sph[i].ent_half = ptcl_sph[i].ent + 0.5 * dt * ptcl_sph[i].ent_dot;
#endif
    }
}

void full_drift(int psys_num_nbody, 
                int psys_num_sph,
                double dt) {
    // Update N-body system
    int nptcl_loc = fdps_get_nptcl_loc(psys_num_nbody);
    FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        ptcl_nbody[i].pos.x += dt * ptcl_nbody[i].vel.x;
        ptcl_nbody[i].pos.y += dt * ptcl_nbody[i].vel.y;
        ptcl_nbody[i].pos.z += dt * ptcl_nbody[i].vel.z;
    }

    // Update SPH system
    nptcl_loc = fdps_get_nptcl_loc(psys_num_sph);
    FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
    for (i = 0; i < nptcl_loc; i++) {
        ptcl_sph[i].pos.x += dt * ptcl_sph[i].vel_half.x;
        ptcl_sph[i].pos.y += dt * ptcl_sph[i].vel_half.y;
        ptcl_sph[i].pos.z += dt * ptcl_sph[i].vel_half.z;
    }
}

void predict(int psys_num, double dt) {
    // Update SPH system
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    FP_sph *ptcl = (FP_sph *) fdps_get_psys_cptr(psys_num);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        fdps_f64vec acc_tot;
        acc_tot.x = ptcl[i].acc_grav.x + ptcl[i].acc_hydro.x;
        acc_tot.y = ptcl[i].acc_grav.y + ptcl[i].acc_hydro.y;
        acc_tot.z = ptcl[i].acc_grav.z + ptcl[i].acc_hydro.z;
        ptcl[i].vel.x += dt * acc_tot.x;
        ptcl[i].vel.y += dt * acc_tot.y;
        ptcl[i].vel.z += dt * acc_tot.z;
#if !defined(ISOTHERMAL_EOS)
        ptcl[i].eng += dt * ptcl[i].eng_dot;
        ptcl[i].ent += dt * ptcl[i].ent_dot;
#endif
    }
}

void final_kick(int psys_num_nbody,
                int psys_num_sph,
                double dt) {
    // Update N-body system
    int nptcl_loc = fdps_get_nptcl_loc(psys_num_nbody);
    FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
    int i;
    for (i = 0; i < nptcl_loc; i++) {
        ptcl_nbody[i].vel.x += 0.5 * dt * ptcl_nbody[i].acc.x;
        ptcl_nbody[i].vel.y += 0.5 * dt * ptcl_nbody[i].acc.y;
        ptcl_nbody[i].vel.z += 0.5 * dt * ptcl_nbody[i].acc.z;
    }
    // Update SPH system
    nptcl_loc = fdps_get_nptcl_loc(psys_num_sph);
    FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
    for (i = 0; i < nptcl_loc; i++) {
        fdps_f64vec acc_tot;
        acc_tot.x = ptcl_sph[i].acc_grav.x + ptcl_sph[i].acc_hydro.x;
        acc_tot.y = ptcl_sph[i].acc_grav.y + ptcl_sph[i].acc_hydro.y;
        acc_tot.z = ptcl_sph[i].acc_grav.z + ptcl_sph[i].acc_hydro.z;
        ptcl_sph[i].vel.x = ptcl_sph[i].vel_half.x + 0.5 * dt * acc_tot.x;
        ptcl_sph[i].vel.y = ptcl_sph[i].vel_half.y + 0.5 * dt * acc_tot.y;
        ptcl_sph[i].vel.z = ptcl_sph[i].vel_half.z + 0.5 * dt * acc_tot.z;
#if !defined(ISOTHERMAL_EOS)
        ptcl_sph[i].eng = ptcl_sph[i].eng_half + 0.5 * dt * ptcl_sph[i].eng_dot;
        ptcl_sph[i].ent = ptcl_sph[i].ent_half + 0.5 * dt * ptcl_sph[i].ent_dot;
#endif
    }
#if defined(DUMP_VELOCITY_OF_SPH_PARTICLE)
    const double coeff = 0.1;
    for (i = 0; i < nptcl_loc; i++) {
        double frac_dump = exp(- coeff * (CFL_hydro/0.1) 
                              * ptcl_sph[i].snds * dt / ptcl_sph[i].smth);
        ptcl_sph[i].vel.x *= frac_dump;
        ptcl_sph[i].vel.y *= frac_dump;
        ptcl_sph[i].vel.z *= frac_dump;
    }
#endif
}
