#include "user_defined.h"

double S2_pcut(double xi) {
    // This is the potential cutoff function where we used Eq.(8.75)
    // in Hockney & Eastwood (1987).
    if (xi <= 1.0) {
       return 1.0
              - xi*(208.0
              + (xi*xi)*(-112.0
              + (xi*xi)*(56.0
              + xi*(-14.0
              + xi*(-8.0
              + 3.0*xi)))))/140.0;
    } else if ((1.0 < xi) && (xi < 2.0)) {
       return 1.0
              - (12.0
              + xi*(128.0
              + xi*(224.0
              + xi*(-448.0
              + xi*(280.0
              + xi*(-56.0
              + xi*(-14.0
              + xi*(8.0
              - xi))))))))/140.0;
    } else {
       return 0.0;
    }
}

double S2_fcut(double xi) {
    // This function returns 1 - R(\xi), where \xi is r/(a/2), a is the
    // scale length of the cutoff function, and R(\xi) is almost the same
    // as the function defined as Eq.(8-72) in Hockney & Eastwood (1987).
    // The only difference is that [1/(r/(a/2))]^2 is factored out 
    // in this function from Eq.(8-72).
    if (xi <= 1.0) {
       return 1.0
              - (xi*xi*xi)*(224.0
              + (xi*xi)*(-224.0
              + xi*(70.0
              + xi*(48.0
              - 21.0*xi))))/140.0;
    } else if ((1.0 < xi) && (xi < 2.0)) {
       return 1.0
              - (12.0
              + (xi*xi)*(-224.0
              + xi*(896.0
              + xi*(-840.0
              + xi*(224.0
              + xi*(70.0
              + xi*(-48.0
              + 7.0*xi)))))))/140.0;
    } else {
       return 0.0;
    }
}
  
void calc_force_ep_ep(EP_nbody *ep_i,
                      int n_ip,
                      EP_nbody *ep_j,
                      int n_jp,
                      Force_pp *f) {
    int i,j;
    for (i = 0; i < n_ip; i++) {
        for (j = 0; j < n_jp; j++) {
            fdps_f64vec dr;
            dr.x = ep_i[i].pos.x - ep_j[j].pos.x;
            dr.y = ep_i[i].pos.y - ep_j[j].pos.y;
            dr.z = ep_i[i].pos.z - ep_j[j].pos.z;
            double rij  = sqrt(dr.x * dr.x
                              +dr.y * dr.y
                              +dr.z * dr.z);
            if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0)) continue;
            double rinv = 1.0/rij;
            double rinv3 = rinv*rinv*rinv;
            double xi = 2.0*rij/ep_i[i].rcut;
            f[i].pot   += ep_j[j].mass * S2_pcut(xi) * rinv;
            f[i].acc.x += ep_j[j].mass * S2_fcut(xi) * rinv3 * dr.x;
            f[i].acc.y += ep_j[j].mass * S2_fcut(xi) * rinv3 * dr.y;
            f[i].acc.z += ep_j[j].mass * S2_fcut(xi) * rinv3 * dr.z;
        }
        // Self-interaction term
        f[i].pot -= ep_i[i].mass * (208.0/(70.0*ep_i[i].rcut));
    }
}

void calc_force_ep_sp(EP_nbody *ep_i,
                      int n_ip,
                      fdps_spj_monopole_cutoff *ep_j,
                      int n_jp,
                      Force_pp *f) {
   int i,j;
   for (i = 0; i < n_ip; i++) {
      for (j = 0; j < n_jp; j++) {
         fdps_f64vec dr;
         dr.x = ep_i[i].pos.x - ep_j[j].pos.x;
         dr.y = ep_i[i].pos.y - ep_j[j].pos.y;
         dr.z = ep_i[i].pos.z - ep_j[j].pos.z;
         double rij = sqrt(dr.x * dr.x
                          +dr.y * dr.y
                          +dr.z * dr.z);
         double rinv = 1.0/rij;
         double rinv3 = rinv*rinv*rinv;
         double xi = 2.0*rij/ep_i[i].rcut;
         f[i].pot   += ep_j[j].mass * S2_pcut(xi) * rinv;
         f[i].acc.x += ep_j[j].mass * S2_fcut(xi) * rinv3 * dr.x;
         f[i].acc.y += ep_j[j].mass * S2_fcut(xi) * rinv3 * dr.y;
         f[i].acc.z += ep_j[j].mass * S2_fcut(xi) * rinv3 * dr.z;
      }
   }
}
