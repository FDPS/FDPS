#include "user_defined.h"

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
