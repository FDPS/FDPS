#include "user_defined.h"

#ifdef USE_PIKG_KERNEL
#include <stdlib.h>
#include "kernel_epep.h"
void calc_gravity_ep_ep(Full_particle *ep_i,
                        int n_ip,
                        Full_particle *ep_j,
                        int n_jp,
                        Full_particle* f)
{
    if (n_ip > 0) {
        Epi_grav * ep_i_tmp = (Epi_grav *) malloc(sizeof(Epi_grav) * n_ip);
        Force_grav * f_tmp  = (Force_grav *) malloc(sizeof(Force_grav) * n_ip);
        int i;
        for (i = 0; i < n_ip; i++) {
          ep_i_tmp[i].pos.x = ep_i[i].pos.x - ep_i[0].pos.x;
          ep_i_tmp[i].pos.y = ep_i[i].pos.y - ep_i[0].pos.y;
          ep_i_tmp[i].pos.z = ep_i[i].pos.z - ep_i[0].pos.z;
          f_tmp[i].acc.x = 0.0;
          f_tmp[i].acc.y = 0.0;
          f_tmp[i].acc.z = 0.0;
          f_tmp[i].pot = 0.0;
        }
        Epj_grav * ep_j_tmp = (Epj_grav *) malloc(sizeof(Epj_grav) * n_jp);
        for (i = 0; i < n_jp; i++) {
          ep_j_tmp[i].pos.x = ep_j[i].pos.x - ep_i[0].pos.x;
          ep_j_tmp[i].pos.y = ep_j[i].pos.y - ep_i[0].pos.y;
          ep_j_tmp[i].pos.z = ep_j[i].pos.z - ep_i[0].pos.z;
          ep_j_tmp[i].mass = ep_j[i].mass;
        }
        pikg_calc_grav_ep_ep(ep_i_tmp, n_ip, ep_j_tmp, n_jp, f_tmp);
        for (i = 0; i < n_ip; i++) {
          f[i].acc.x += f_tmp[i].acc.x;
          f[i].acc.y += f_tmp[i].acc.y;
          f[i].acc.z += f_tmp[i].acc.z;
          f[i].pot += f_tmp[i].pot;
        }
        free(ep_i_tmp);
        free(f_tmp);
        free(ep_j_tmp);
    }
}

void calc_gravity_ep_sp(Full_particle *ep_i,
                        int n_ip,
                        fdps_spj_monopole *ep_j,
                        int n_jp,
                        Full_particle *f)
{
    if (n_ip > 0) {
        Epi_grav * ep_i_tmp = (Epi_grav *) malloc(sizeof(Epi_grav) * n_ip);
        Force_grav * f_tmp  = (Force_grav *) malloc(sizeof(Force_grav) * n_ip);
        int i;
        for (i = 0; i < n_ip; i++) {
          ep_i_tmp[i].pos.x = ep_i[i].pos.x - ep_i[0].pos.x;
          ep_i_tmp[i].pos.y = ep_i[i].pos.y - ep_i[0].pos.y;
          ep_i_tmp[i].pos.z = ep_i[i].pos.z - ep_i[0].pos.z;
          f_tmp[i].acc.x = 0.0;
          f_tmp[i].acc.y = 0.0;
          f_tmp[i].acc.z = 0.0;
          f_tmp[i].pot = 0.0;
        }
        Epj_grav * ep_j_tmp = (Epj_grav *) malloc(sizeof(Epj_grav) * n_jp);
        for (i = 0; i < n_jp; i++) {
          ep_j_tmp[i].pos.x = ep_j[i].pos.x - ep_i[0].pos.x;
          ep_j_tmp[i].pos.y = ep_j[i].pos.y - ep_i[0].pos.y;
          ep_j_tmp[i].pos.z = ep_j[i].pos.z - ep_i[0].pos.z;
          ep_j_tmp[i].mass = ep_j[i].mass;
        }
        pikg_calc_grav_ep_ep(ep_i_tmp, n_ip, ep_j_tmp, n_jp, f_tmp);
        for (i = 0; i < n_ip; i++) {
          f[i].acc.x += f_tmp[i].acc.x;
          f[i].acc.y += f_tmp[i].acc.y;
          f[i].acc.z += f_tmp[i].acc.z;
          f[i].pot += f_tmp[i].pot;
        }
        free(ep_i_tmp);
        free(f_tmp);
        free(ep_j_tmp);
    }
}

#else // USE_PIKG_KERNEL
void calc_gravity_ep_ep(Full_particle *ep_i,
                        int n_ip,
                        Full_particle *ep_j,
                        int n_jp,
                        Full_particle* f)
{
    int i, j;
    for (i=0; i<n_ip ;i++) {
        Full_particle *pi = ep_i + i;
        double eps2 = pi->eps * pi->eps;
        double xi = pi->pos.x;
        double yi = pi->pos.y;
        double zi = pi->pos.z;
        double ax, ay, az, pot;
        ax = ay = az = pot = 0;
        for (j=0; j<n_jp; j++) {
            Full_particle *pj = ep_j + j;
            double dx = xi - pj->pos.x;
            double dy = yi - pj->pos.y;
            double dz = zi - pj->pos.z;
            double r2 = dx*dx+dy*dy+dz*dz+eps2;
            double rinv = 1.0/sqrt(r2);
            double mrinv = pj->mass* rinv;
            double mr3inv = mrinv*rinv*rinv;
            ax -= dx*mr3inv;
            ay -= dy*mr3inv;
            az -= dz*mr3inv;
            pot = pot - mrinv;
        }
        Full_particle *pfi = f+i;
        pfi->pot += pot;
        pfi->acc.x += ax;
        pfi->acc.y += ay;
        pfi->acc.z += az;
    }
}

void calc_gravity_ep_sp(Full_particle *ep_i,
                        int n_ip,
                        fdps_spj_monopole *ep_j,
                        int n_jp,
                        Full_particle *f)
{
    int i, j;
    for (i=0; i<n_ip; i++) {
        Full_particle *pi = ep_i + i;
        double eps2 = pi->eps*pi->eps;
        double xi = pi->pos.x;
        double yi = pi->pos.y;
        double zi = pi->pos.z;
        double ax, ay, az, pot;
        ax = ay = az = pot = 0;
        for (j=0; j<n_jp; j++) {
            fdps_spj_monopole *pj = ep_j + j;
            double dx = xi - pj->pos.x;
            double dy = yi - pj->pos.y;
            double dz = zi - pj->pos.z;
            double r2 = dx*dx+dy*dy+dz*dz+eps2;
            double rinv = 1.0/sqrt(r2);
            double mrinv = pj->mass* rinv;
            double mr3inv = mrinv*rinv*rinv;
            ax -= dx*mr3inv;
            ay -= dy*mr3inv;
            az -= dz*mr3inv;
            pot = pot - mrinv;
        }
        Full_particle *pfi = f+i;
        pfi->pot += pot;
        pfi->acc.x += ax;
        pfi->acc.y += ay;
        pfi->acc.z += az;
    }
}
#endif // USE_PIKG_KERNEL
