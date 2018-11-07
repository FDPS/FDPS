#include "user_defined.h"

/* Global variable */ 
const double kernel_support_radius=2.5;

/* Kernel functions */
double W(fdps_f64vec dr, double h) {
    double s,s1,s2,ret;
    s = sqrt(dr.x * dr.x 
            +dr.y * dr.y 
            +dr.z * dr.z)/h;
    s1 = 1.0 - s;
    if (s1 < 0.0) s1 = 0.0;
    s2 = 0.5 - s;
    if (s2 < 0.0) s2 = 0.0;
    ret = (s1*s1*s1) - 4.0*(s2*s2*s2);
    ret = ret * 16.0e0/(pi*h*h*h);
    return ret;
}

fdps_f64vec gradW(fdps_f64vec dr, double h) {
    double dr_abs,s,s1,s2,coef;
    fdps_f64vec ret;
    dr_abs = sqrt(dr.x * dr.x 
                 +dr.y * dr.y 
                 +dr.z * dr.z);
    s = dr_abs/h;
    s1 = 1.0 - s;
    if (s1 < 0.0) s1 = 0.0;
    s2 = 0.5 - s;
    if (s2 < 0.0) s2 = 0.0;
    coef = - 3.0*(s1*s1) + 12.0*(s2*s2);
    coef = coef * 16.0/(pi*h*h*h);
    coef = coef / (dr_abs*h + 1.0e-6*h);
    ret.x = dr.x * coef;
    ret.y = dr.y * coef;
    ret.z = dr.z * coef;
    return ret;
}

/* Interaction functions */
void calc_density(Essential_particle *ep_i,
                  int n_ip,
                  Essential_particle *ep_j,
                  int n_jp,
                  Force_dens *f) {
    int i,j;
    fdps_f64vec dr;
    for (i = 0; i< n_ip; i++) {
        for (j = 0; j < n_jp; j++) {
            dr.x = ep_j[j].pos.x - ep_i[i].pos.x;
            dr.y = ep_j[j].pos.y - ep_i[i].pos.y;
            dr.z = ep_j[j].pos.z - ep_i[i].pos.z;
            f[i].dens += ep_j[j].mass * W(dr,ep_i[i].smth);
        }
    }

}

void calc_hydro_force(Essential_particle *ep_i,
                      int n_ip, 
                      Essential_particle *ep_j,
                      int n_jp,
                      Force_hydro *f) {
    // Local parameters
    const double C_CFL=0.3;
    // Local variables
    int i,j;
    double mass_i,mass_j,smth_i,smth_j,
           dens_i,dens_j,pres_i,pres_j,
           snds_i,snds_j;
    double povrho2_i,povrho2_j, 
           v_sig_max,dr_dv,w_ij,v_sig,AV;
    fdps_f64vec pos_i,pos_j,vel_i,vel_j,
                dr,dv,gradW_i,gradW_j,gradW_ij;

    for (i = 0; i < n_ip; i++) {
        // Zero-clear
        v_sig_max = 0.0;
        // Extract i-particle info.
        pos_i.x = ep_i[i].pos.x;
        pos_i.y = ep_i[i].pos.y;
        pos_i.z = ep_i[i].pos.z;
        vel_i.x = ep_i[i].vel.x;
        vel_i.y = ep_i[i].vel.y;
        vel_i.z = ep_i[i].vel.z;
        mass_i  = ep_i[i].mass;
        smth_i  = ep_i[i].smth;
        dens_i  = ep_i[i].dens;
        pres_i  = ep_i[i].pres;
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
            if (v_sig > v_sig_max) v_sig_max=v_sig;
            // Compute the artificial viscosity
            AV = - 0.5*v_sig*w_ij / (0.5*(dens_i+dens_j));
            // Compute the average of the gradients of kernel
            gradW_i  = gradW(dr,smth_i);
            gradW_j  = gradW(dr,smth_j);
            gradW_ij.x = 0.5 * (gradW_i.x + gradW_j.x);
            gradW_ij.y = 0.5 * (gradW_i.y + gradW_j.y);
            gradW_ij.z = 0.5 * (gradW_i.z + gradW_j.z);
            // Compute the acceleration and the heating rate
            f[i].acc.x -= mass_j*(povrho2_i+povrho2_j+AV)*gradW_ij.x;
            f[i].acc.y -= mass_j*(povrho2_i+povrho2_j+AV)*gradW_ij.y;
            f[i].acc.z -= mass_j*(povrho2_i+povrho2_j+AV)*gradW_ij.z;
            f[i].eng_dot += mass_j * (povrho2_i + 0.5*AV)
                            *(dv.x * gradW_ij.x
                             +dv.y * gradW_ij.y
                             +dv.z * gradW_ij.z);
        }
        f[i].dt = C_CFL*2.0*smth_i/(v_sig_max*kernel_support_radius);
    }
}
