#include "user_defined.h"

/* Global variables */
double specific_heat_ratio = 5.0e0/3.0e0;
double CFL_dyn = 0.3e0;
double CFL_hydro = 0.3e0;
double SCF_smth = 1.25e0;
int N_neighbor = 50;
double eps_grav;
double mass_avg;
double dt_max;

/* Kernel functions */
double W(double r, double h) {
    // M4 Cubic spline kernel
    // (see Eq. (4) in Springel (2005)[MNRAS,364,1105])
    double u = r/h;
    double cc=8.0/(pi*h*h*h);
    if (u <= 0.5) {
        double u2 = u*u;
        return cc*(1.0+u2*6.0*(-1.0+u));
    } else if ((0.5 < u) && (u <= 1.0)) {
        double s = 1.0-u;
        return cc*2.0*s*s*s;
    } else {
        return 0.0;
    }
}
fdps_f64vec gradW(fdps_f64vec dr, double h) {
    // This subroutine gives \nabla W(r,h), i.e.,
    // \dfrac{\partial W(r,h)}{\partial r}\dfrac{dr}{r}.
    double r = sqrt(dr.x * dr.x
                   +dr.y * dr.y
                   +dr.z * dr.z);
    double u=r/h;
    double cc = -48.0/(pi*h*h*h*h);
#if defined(USE_PRESCR_OF_THOMAS_COUCHMAN_1992)
    if (u <= 1.0/3.0) {
        cc = cc*(1.0/3.0)/(r);
    } else if ((1.0/3.0 < u) && (u <= 0.5)) {
        cc = cc*u*(2.0-3.0*u)/(r);
    } else if ((0.5 < u) && (u < 1.0)) {
        cc = cc*(1.0-u)*(1.0-u)/(r);
    } else {
        cc = 0.0;
    }
#else
    if ((0.0 < u) && (u <= 0.5)) {
        cc = cc*u*(2.0-3.0*u)/(r);
    } else if ((0.5 < u) && (u < 1.0)) {
        cc = cc*(1.0-u)*(1.0-u)/(r);
    } else {
        // r=0 case is included in this branch
        cc = 0.0;
    }
#endif
    fdps_f64vec v;
    v.x = dr.x * cc;
    v.y = dr.y * cc;
    v.z = dr.z * cc;
    return v;
}

double dWdh(double r, double h) {
    // This subroutine gives dW(r,h)/dh, i.e.,
    // \dfrac{\partial W(r,h)}{\partial h}.
    double u=r/h;
    double cc=-24.0/(pi*h*h*h*h);
    if (u <= 0.5) {
        double u2 = u*u;
        return cc*(1.0
                  +u2*(-10.0
                       +12.0*u));
    } else if ((0.5 < u) && (u < 1.0)) {
        double s = 1.0-u;
        return cc*2.0*s*s*(1.0-2.0*u);
    } else {
        return 0.0;
    }
}

/* Interaction functions */
#if defined(ENABLE_PHANTOM_GRAPE_X86)
void calc_gravity_ep_ep(struct ep_grav *ep_i,
                        int n_ip,
                        struct ep_grav *ep_j,
                        int n_jp,
                        struct force_grav *f) {
    int i,j;
    int nipipe = n_ip;
    int njpipe = n_jp;
    double (*xi)[3] = (double (*)[3])malloc(sizeof(double) * nipipe * 3);
    double (*ai)[3] = (double (*)[3])malloc(sizeof(double) * nipipe * 3);
    double  *pi     = (double  *    )malloc(sizeof(double) * nipipe);
    double (*xj)[3] = (double (*)[3])malloc(sizeof(double) * njpipe * 3);
    double  *mj     = (double  *    )malloc(sizeof(double) * njpipe);
    for (i = 0; i < n_ip; i++) {
        xi[i][0] = ep_i[i].pos.x;
        xi[i][1] = ep_i[i].pos.y;
        xi[i][2] = ep_i[i].pos.z;
        ai[i][0] = 0.0;
        ai[i][1] = 0.0;
        ai[i][2] = 0.0;
        pi[i]    = 0.0;
    }
    for (j = 0; j < n_jp; j++) {
        xj[j][0] = ep_j[j].pos.x;
        xj[j][1] = ep_j[j].pos.y;
        xj[j][2] = ep_j[j].pos.z;
        mj[j]    = ep_j[j].mass;
    }
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
    int devid = omp_get_thread_num();
    // [IMPORTANT NOTE]
    //   The function calc_gravity_ep_ep is called by a OpenMP thread
    //   in the FDPS. This means that here is already in the parallel region.
    //   So, you can use omp_get_thread_num() without !$OMP parallel directives.
    //   If you use them, a nested parallel resions is made and the gravity
    //   calculation will not be performed correctly.
#else
    int devid = 0;
#endif
    g5_set_xmjMC(devid, 0, n_jp, xj, mj);
    g5_set_nMC(devid, n_jp);
    g5_calculate_force_on_xMC(devid, xi, ai, pi, n_ip);
    for (i = 0; i < n_ip; i++) {
       f[i].acc.x += ai[i][0];
       f[i].acc.y += ai[i][1];
       f[i].acc.z += ai[i][2];
       f[i].pot   -= pi[i];
    }
    free(xi);
    free(ai);
    free(pi);
    free(xj);
    free(mj);
}

void calc_gravity_ep_sp(struct ep_grav *ep_i,
                        int n_ip,
                        fdps_spj_monopole *ep_j,
                        int n_jp,
                        struct force_grav *f) {
    int i,j;
    int nipipe = n_ip;
    int njpipe = n_jp;
    double (*xi)[3] = (double (*)[3])malloc(sizeof(double) * nipipe * 3);
    double (*ai)[3] = (double (*)[3])malloc(sizeof(double) * nipipe * 3);
    double  *pi     = (double  *    )malloc(sizeof(double) * nipipe);
    double (*xj)[3] = (double (*)[3])malloc(sizeof(double) * njpipe * 3);
    double  *mj     = (double  *    )malloc(sizeof(double) * njpipe);
    for (i = 0; i < n_ip; i++) {
        xi[i][0] = ep_i[i].pos.x;
        xi[i][1] = ep_i[i].pos.y;
        xi[i][2] = ep_i[i].pos.z;
        ai[i][0] = 0.0;
        ai[i][1] = 0.0;
        ai[i][2] = 0.0;
        pi[i]    = 0.0;
    }
    for (j = 0; j < n_jp; j++) {
        xj[j][0] = ep_j[j].pos.x;
        xj[j][1] = ep_j[j].pos.y;
        xj[j][2] = ep_j[j].pos.z;
        mj[j]    = ep_j[j].mass;
    }
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
    int devid = omp_get_thread_num();
    // [IMPORTANT NOTE]
    //   The function calc_gravity_ep_ep is called by a OpenMP thread
    //   in the FDPS. This means that here is already in the parallel region.
    //   So, you can use omp_get_thread_num() without !$OMP parallel directives.
    //   If you use them, a nested parallel resions is made and the gravity
    //   calculation will not be performed correctly.
#else
    int devid = 0;
#endif
    g5_set_xmjMC(devid, 0, n_jp, xj, mj);
    g5_set_nMC(devid, n_jp);
    g5_calculate_force_on_xMC(devid, xi, ai, pi, n_ip);
    for (i = 0; i < n_ip; i++) {
       f[i].acc.x += ai[i][0];
       f[i].acc.y += ai[i][1];
       f[i].acc.z += ai[i][2];
       f[i].pot   -= pi[i];
    }
    free(xi);
    free(ai);
    free(pi);
    free(xj);
    free(mj);
}
#else
void calc_gravity_ep_ep(struct ep_grav *ep_i,
                        int n_ip,
                        struct ep_grav *ep_j,
                        int n_jp,
                        struct force_grav *f) {
    int i,j;
    double eps2 = eps_grav * eps_grav;
    for (i = 0; i < n_ip; i++) {
        fdps_f64vec xi,ai;
        double poti;
        xi.x = ep_i[i].pos.x;
        xi.y = ep_i[i].pos.y;
        xi.z = ep_i[i].pos.z;
        ai.x = 0.0;
        ai.y = 0.0;
        ai.z = 0.0;
        poti = 0.0;
        for (j = 0; j < n_jp; j++) {
            fdps_f64vec rij;
            rij.x  = xi.x - ep_j[j].pos.x;
            rij.y  = xi.y - ep_j[j].pos.y;
            rij.z  = xi.z - ep_j[j].pos.z;
            double r3_inv = rij.x * rij.x 
                          + rij.y * rij.y 
                          + rij.z * rij.z 
                          + eps2;
            double r_inv  = 1.0/sqrt(r3_inv);
            r3_inv = r_inv * r_inv;
            r_inv  = r_inv * ep_j[j].mass;
            r3_inv = r3_inv * r_inv;
            ai.x -= r3_inv * rij.x;
            ai.y -= r3_inv * rij.y;
            ai.z -= r3_inv * rij.z;
            poti -= r_inv;
        }
        f[i].acc.x += ai.x;
        f[i].acc.y += ai.y;
        f[i].acc.z += ai.z;
        f[i].pot   += poti;
    }
}

void calc_gravity_ep_sp(struct ep_grav *ep_i,
                        int n_ip,
                        fdps_spj_monopole *ep_j,
                        int n_jp,
                        struct force_grav *f) {
    int i,j;
    double eps2 = eps_grav * eps_grav;
    for (i = 0; i < n_ip; i++) {
        fdps_f64vec xi,ai;
        double poti;
        xi.x = ep_i[i].pos.x;
        xi.y = ep_i[i].pos.y;
        xi.z = ep_i[i].pos.z;
        ai.x = 0.0;
        ai.y = 0.0;
        ai.z = 0.0;
        poti = 0.0;
        for (j = 0; j < n_jp; j++) {
            fdps_f64vec rij;
            rij.x  = xi.x - ep_j[j].pos.x;
            rij.y  = xi.y - ep_j[j].pos.y;
            rij.z  = xi.z - ep_j[j].pos.z;
            double r3_inv = rij.x * rij.x 
                          + rij.y * rij.y 
                          + rij.z * rij.z 
                          + eps2;
            double r_inv  = 1.0/sqrt(r3_inv);
            r3_inv = r_inv * r_inv;
            r_inv  = r_inv * ep_j[j].mass;
            r3_inv = r3_inv * r_inv;
            ai.x -= r3_inv * rij.x;
            ai.y -= r3_inv * rij.y;
            ai.z -= r3_inv * rij.z;
            poti -= r_inv;
        }
        f[i].acc.x += ai.x;
        f[i].acc.y += ai.y;
        f[i].acc.z += ai.z;
        f[i].pot   += poti;
    }
}
#endif

void calc_density(struct ep_hydro *ep_i,
                  int n_ip,
                  struct ep_hydro *ep_j,
                  int n_jp,
                  struct force_dens *f) {
#if defined(ENABLE_VARIABLE_SMOOTHING_LENGTH)
    // Local parameters
    const double eps=1.0e-6;
    // Local variables
    int i,j;
    int n_unchanged;
    double M,M_trgt;
    double dens,drho_dh;
    double h,h_max_alw,h_L,h_U,dh,dh_prev;
    fdps_f64vec dr,dv,gradW_i;
    double *mj  = (double *)malloc(sizeof(double) * n_jp);
    double *rij = (double *)malloc(sizeof(double) * n_jp);
    M_trgt = mass_avg * N_neighbor;
    for (i = 0; i < n_ip; i++) {
        dens = 0.0;
        h_max_alw = ep_i[i].smth; // maximum allowance
        h = h_max_alw / SCF_smth;
        // Note that we increase smth by a factor of scf_smth
        // before calling calc_density().
        h_L = 0.0;
        h_U = h_max_alw;
        dh_prev = 0.0;
        n_unchanged = 0;
        // Software cache
        for (j = 0; j < n_jp; j++) {
            mj[j] = ep_j[j].mass;
            dr.x = ep_i[i].pos.x - ep_j[j].pos.x;
            dr.y = ep_i[i].pos.y - ep_j[j].pos.y;
            dr.z = ep_i[i].pos.z - ep_j[j].pos.z;
            rij[j] = sqrt(dr.x * dr.x
                         +dr.y * dr.y
                         +dr.z * dr.z);
        }
        for(;;) {
            // Calculate density
            dens = 0.0;
            for (j = 0; j < n_jp; j++)
               dens += mj[j] * W(rij[j], h);
            // Check if the current value of the smoohting length satisfies 
            // Eq.(5) in Springel (2005).
            M = 4.0 * pi * h * h * h * dens / 3.0;
            if ((h < h_max_alw) && (fabs(M/M_trgt - 1.0) < eps)) {
                // In this case, Eq.(5) holds within a specified accuracy.
                f[i].flag = 1;
                f[i].dens = dens;
                f[i].smth = h;
                break;
            }
            if (((h == h_max_alw) && (M < M_trgt)) || (n_unchanged == 4)) {
                // In this case, we skip this particle forcibly.
                // In order to determine consistently the density
                // and the smoohting length for this particle,
                // we must re-perform calcForceAllAndWriteBack().
                f[i].flag = 0;
                f[i].dens = dens;
                f[i].smth = h_max_alw;
                break;
            }
            // Update h_L & h_U
            if (M < M_trgt) {
                if (h_L < h) h_L = h;
            } else if (M_trgt < M) {
                if (h < h_U) h_U = h;
            }
            dh = h_U - h_L;
            if (dh == dh_prev) {
               n_unchanged++;
            } else {
               dh_prev = dh;
               n_unchanged = 0;
            }
            // Update smoothing length
            h = pow((3.0 * M_trgt)/(4.0 * pi * dens),1.0/3.0);
            if ((h <= h_L) || (h == h_U)) {
               // In this case, we switch to the bisection search.
               // The inclusion of '=' in the if statement is very
               // important to escape a limit cycle.
               h = 0.5 * (h_L + h_U);
            } else if (h_U < h) {
               h = h_U;
            }
        }
        // Calculate grad-h term
        if (f[i].flag == 1) {
            drho_dh = 0.0;
            for (j = 0; j < n_jp; j++)
               drho_dh += mj[j] * dWdh(rij[j], h);
            f[i].gradh = 1.0 / (1.0 + (h * drho_dh) / (3.0 * dens));
        } else {
            f[i].gradh = 1.0; // dummy value
        }
        // Compute \div v & \rot v for Balsara switch
#if defined(USE_BALSARA_SWITCH)
        for (j = 0; j < n_jp; j++) {
           dr.x = ep_i[i].pos.x - ep_j[j].pos.x;
           dr.y = ep_i[i].pos.y - ep_j[j].pos.y;
           dr.z = ep_i[i].pos.z - ep_j[j].pos.z;
           dv.x = ep_i[i].vel.x - ep_j[j].vel.x;
           dv.y = ep_i[i].vel.y - ep_j[j].vel.y;
           dv.z = ep_i[i].vel.z - ep_j[j].vel.z;
           gradW_i = gradW(dr, f[i].smth);
           f[i].divv -= mj[j] * (dv.x * gradW_i.x
                                +dv.y * gradW_i.y
                                +dv.z * gradW_i.z);
           f[i].rotv.x -= mj[j] * (dv.y * gradW_i.z - dv.z * gradW_i.y);
           f[i].rotv.y -= mj[j] * (dv.z * gradW_i.x - dv.x * gradW_i.z);
           f[i].rotv.z -= mj[j] * (dv.x * gradW_i.y - dv.y * gradW_i.x);
        }
        f[i].divv   /= f[i].dens;
        f[i].rotv.x /= f[i].dens;
        f[i].rotv.y /= f[i].dens;
        f[i].rotv.z /= f[i].dens;
#endif
    }
    free(mj);
    free(rij);
#else
    int i,j;
    for (i = 0; i < n_ip; i++) {
        f[i].dens = 0.0;
        for (j = 0; j < n_jp; j++) {
            fdps_f64vec dr;
            dr.x = ep_j[j].pos.x - ep_i[i].pos.x;
            dr.y = ep_j[j].pos.y - ep_i[i].pos.y;
            dr.z = ep_j[j].pos.z - ep_i[i].pos.z;
            double rij = sqrt(dr.x * dr.x
                             +dr.y * dr.y
                             +dr.z * dr.z);
            f[i].dens += ep_j[j].mass * W(rij,ep_i[i].smth);
        }
        f[i].smth = ep_i[i].smth;
        f[i].gradh = 1.0;
        // Compute \div v & \rot v for Balsara switch
#if defined(USE_BALSARA_SWITCH)
        for (j = 0; j < n_jp; j++) {
            double mj = ep_j[j].mass;
            fdps_f64vec dr,dv,gradW_i;
            dr.x = ep_i[i].pos.x - ep_j[j].pos.x;
            dr.y = ep_i[i].pos.y - ep_j[j].pos.y;
            dr.z = ep_i[i].pos.z - ep_j[j].pos.z;
            dv.x = ep_i[i].vel.x - ep_j[j].vel.x;
            dv.y = ep_i[i].vel.y - ep_j[j].vel.y;
            dv.z = ep_i[i].vel.z - ep_j[j].vel.z;
            gradW_i = gradW(dr, f[i].smth);
            f[i].divv -= mj * (dv.x * gradW_i.x
                              +dv.y * gradW_i.y
                              +dv.z * gradW_i.z);
            f[i].rotv.x -= mj * (dv.y * gradW_i.z - dv.z * gradW_i.y);
            f[i].rotv.y -= mj * (dv.z * gradW_i.x - dv.x * gradW_i.z);
            f[i].rotv.z -= mj * (dv.x * gradW_i.y - dv.y * gradW_i.x);
        }
        f[i].divv   /= f[i].dens;
        f[i].rotv.x /= f[i].dens;
        f[i].rotv.y /= f[i].dens;
        f[i].rotv.z /= f[i].dens;
#endif
    }
#endif
} 

void calc_hydro_force(struct ep_hydro *ep_i,
                      int n_ip,
                      struct ep_hydro *ep_j,
                      int n_jp,
                      struct force_hydro *f) {
    // Local variables
    int i,j;
    double mass_i,mass_j,smth_i,smth_j,
           dens_i,dens_j,pres_i,pres_j,
           gradh_i,gradh_j,balsw_i,balsw_j,
           snds_i,snds_j;
    double povrho2_i,povrho2_j,
           v_sig_max,dr_dv,w_ij,v_sig,AV;
    fdps_f64vec pos_i,pos_j,vel_i,vel_j,
                dr,dv,gradW_i,gradW_j,gradW_ij;
    for (i = 0; i < n_ip; i++) {
        // Zero-clear
        v_sig_max = 0.0;
        // Extract i-particle info.
        pos_i = ep_i[i].pos;
        vel_i = ep_i[i].vel;
        mass_i  = ep_i[i].mass;
        smth_i  = ep_i[i].smth;
        dens_i  = ep_i[i].dens;
        pres_i  = ep_i[i].pres;
        gradh_i = ep_i[i].gradh;
        balsw_i = ep_i[i].balsw;
        snds_i  = ep_i[i].snds;
        povrho2_i = pres_i/(dens_i*dens_i);
        for (j = 0; j < n_jp; j++) {
            // Extract j-particle info.
            pos_j.x = ep_j[j].pos.x;
            pos_j.y = ep_j[j].pos.y;
            pos_j.z = ep_j[j].pos.z;
            vel_j.x = ep_j[j].vel.x;
            vel_j.y = ep_j[j].vel.y;
            vel_j.z = ep_j[j].vel.z;
            mass_j  = ep_j[j].mass;
            smth_j  = ep_j[j].smth;
            dens_j  = ep_j[j].dens;
            pres_j  = ep_j[j].pres;
            gradh_j = ep_j[j].gradh;
            balsw_j = ep_j[j].balsw;
            snds_j  = ep_j[j].snds;
            povrho2_j = pres_j/(dens_j*dens_j);
            // Compute dr & dv
            dr.x = pos_i.x - pos_j.x;
            dr.y = pos_i.y - pos_j.y;
            dr.z = pos_i.z - pos_j.z;
            dv.x = vel_i.x - vel_j.x;
            dv.y = vel_i.y - vel_j.y;
            dv.z = vel_i.z - vel_j.z;
            // Compute the signal velocity
            dr_dv = dr.x * dv.x + dr.y * dv.y + dr.z * dv.z;
            if (dr_dv < 0.0) {
               w_ij = dr_dv / sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
            } else {
               w_ij = 0.0;
            }
            v_sig = snds_i + snds_j - 3.0 * w_ij;
            if (v_sig > v_sig_max) v_sig_max = v_sig;
            // Compute the artificial viscosity
            AV = - 0.5*v_sig*w_ij / (0.5*(dens_i+dens_j)) * 0.5*(balsw_i+balsw_j);
            // Compute the average of the gradients of kernel
            gradW_i  = gradW(dr,smth_i);
            gradW_j  = gradW(dr,smth_j);
            gradW_ij.x = 0.5 * (gradW_i.x + gradW_j.x);
            gradW_ij.y = 0.5 * (gradW_i.y + gradW_j.y);
            gradW_ij.z = 0.5 * (gradW_i.z + gradW_j.z);
            // Compute the acceleration and the heating rate
            f[i].acc.x -= mass_j*(gradh_i * povrho2_i * gradW_i.x
                                 +gradh_j * povrho2_j * gradW_j.x
                                 +AV * gradW_ij.x);
            f[i].acc.y -= mass_j*(gradh_i * povrho2_i * gradW_i.y
                                 +gradh_j * povrho2_j * gradW_j.y
                                 +AV * gradW_ij.y);
            f[i].acc.z -= mass_j*(gradh_i * povrho2_i * gradW_i.z
                                 +gradh_j * povrho2_j * gradW_j.z
                                 +AV * gradW_ij.z);
            f[i].eng_dot += mass_j * gradh_i * povrho2_i * (dv.x * gradW_i.x
                                                           +dv.y * gradW_i.y
                                                           +dv.z * gradW_i.z);
                         + mass_j * 0.5 * AV * (dv.x * gradW_ij.x
                                               +dv.y * gradW_ij.y
                                               +dv.z * gradW_ij.z);
            f[i].ent_dot += 0.5 * mass_j * AV * (dv.x * gradW_ij.x
                                                +dv.y * gradW_ij.y
                                                +dv.z * gradW_ij.z);
        }
        f[i].ent_dot *= ((specific_heat_ratio - 1.0) 
                        /pow(dens_i,specific_heat_ratio - 1.0));
        f[i].dt = CFL_hydro*2.0*smth_i/v_sig_max;
    }
}
