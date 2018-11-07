/* Standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
/* FDPS headers */
#include "user_defined.h"
#include "FDPS_c_if.h"

void dump_fullp(Full_particle p)
{
    printf("%lld %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e",
           p.id, p.mass, p.pos.x, p.pos.y, p.pos.z,
           p.vel.x, p.vel.y, p.vel.z);
    printf("%15.7e %15.7e %15.7e %15.7e\n",
           p.acc.x, p.acc.y, p.acc.z, p.pot);
}
void dump_fullpsys(Full_particle *p, int n)
{
    int i;
    for (i=0;i<n;i++)dump_fullp(p[i]);
}

void dump_particles(int psys_num)
{
    Full_particle *ptcl;
    ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int n = fdps_get_nptcl_loc(psys_num);
    dump_fullpsys(ptcl, n);
}

void setup_IC(int psys_num,
              int nptcl_glb)
{

    double m_tot=1.0;
    double rmax=3.0;
    double r2max=rmax*rmax;
    //   Get # of MPI processes and rank number
    int nprocs = fdps_get_num_procs();
    int myrank = fdps_get_rank();
    // Make an initial condition at RANK 0
    if (myrank == 0 ){
        //Set # of local particles
        fdps_set_nptcl_loc(psys_num,nptcl_glb);
        Full_particle *ptcl;
        ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
        //** initialize Mersenne twister
        int mtts_num;
        fdps_create_mtts(&mtts_num);
        fdps_mtts_init_genrand(mtts_num,0);
        int i;
        for (i=0; i < nptcl_glb; i++){
            Full_particle *q = ptcl+i;
            q->id = i;
            q->mass = m_tot/nptcl_glb;
            double r2 = r2max*2;
            fdps_f64vec pos;
            while (r2 >= r2max){
                pos.x= (2*fdps_mtts_genrand_res53(mtts_num)-1.0) * rmax;
                pos.y= (2*fdps_mtts_genrand_res53(mtts_num)-1.0) * rmax;
                pos.z= (2*fdps_mtts_genrand_res53(mtts_num)-1.0) * rmax;
                r2 = pos.x*pos.x + pos.y*pos.y + pos.z*pos.z;
            }
            q->pos = pos;
            q->vel.x = 0.0;
            q->vel.y = 0.0;
            q->vel.z = 0.0;
            q->eps = 1.0/32.0;
        }
        fdps_f64vec cm_pos;
        fdps_f64vec cm_vel;
        cm_pos.x = 0.0; cm_pos.y = 0.0; cm_pos.z = 0.0;
        cm_vel.x = 0.0; cm_vel.y = 0.0; cm_vel.z = 0.0;
        double cm_mass = 0;
        for (i=0; i < nptcl_glb; i++){
            Full_particle *pi = ptcl+i;
            cm_pos.x +=  pi->pos.x* pi->mass;
            cm_pos.y +=  pi->pos.y* pi->mass;
            cm_pos.z +=  pi->pos.z* pi->mass;
            cm_vel.x +=  pi->vel.x* pi->mass;
            cm_vel.y +=  pi->vel.y* pi->mass;
            cm_vel.z +=  pi->vel.z* pi->mass;
            cm_mass += pi->mass;
        }
        cm_pos.x /= cm_mass;
        cm_pos.y /= cm_mass;
        cm_pos.z /= cm_mass;
        cm_vel.x /= cm_mass;
        cm_vel.y /= cm_mass;
        cm_vel.z /= cm_mass;
        for (i=0; i < nptcl_glb; i++){
            Full_particle* q = ptcl+i;
            q->pos.x -= cm_pos.x;
            q->pos.y -= cm_pos.y;
            q->pos.z -= cm_pos.z;
            q->vel.x -= cm_vel.x;
            q->vel.y -= cm_vel.y;
            q->vel.z -= cm_vel.z;
        }
        //dump_fullpsys(ptcl, nptcl_glb);
    } else{
        fdps_set_nptcl_loc(psys_num,0);
    }
}

void calc_energy(int psys_num,
                 double *etot,
                 double *ekin,
                 double *epot)
{
    *etot = *ekin = *epot = 0;
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    Full_particle *ptcl;
    ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    
    double  ekin_loc = 0.0;
    double  epot_loc = 0.0;
    int i;
    for (i=0;i < nptcl_loc; i++){
        Full_particle *pi = ptcl+i;
        fdps_f64vec v = pi->vel;
        ekin_loc += pi->mass * (v.x*v.x+v.y*v.y+v.z*v.z);
        epot_loc += pi->mass * (pi->pot + pi->mass/pi->eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    double etot_loc = ekin_loc + epot_loc;
    *ekin = fdps_get_sum_f64(ekin_loc);
    *epot = fdps_get_sum_f64(epot_loc);
    *etot = fdps_get_sum_f64(etot_loc);
}

void  kick(int psys_num, double dt)
{
    Full_particle *ptcl;
    ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int n = fdps_get_nptcl_loc(psys_num);
    int i;
    for (i=0;i < n; i++){
        Full_particle *pi = ptcl+i;
        fdps_f64vec *pv, *pa;
        pv = &(pi->vel);
        pa = &(pi->acc);
        pv->x += pa->x * dt;
        pv->y += pa->y * dt;
        pv->z += pa->z * dt;
    }
}

void  drift(int psys_num, double dt)
{
    Full_particle *ptcl;
    ptcl = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int n = fdps_get_nptcl_loc(psys_num);
    int i;
    for (i=0;i < n; i++){
        Full_particle *pi = ptcl+i;
        fdps_f64vec *px, *pv; 
        pv = &(pi->vel);
        px = &(pi->pos);
        px->x += pv->x * dt;
        px->y += pv->y * dt;
        px->z += pv->z * dt;
    }
}

int c_main()
{
    fprintf(stderr, "FDPS on C test code\n");
    fdps_initialize();
    // Create and initialize dinfo object
    int dinfo_num;
    float coef_ema=0.3;
    fdps_create_dinfo(&dinfo_num);
    fdps_init_dinfo(dinfo_num,coef_ema);
    // Create and initialize psys object
    int psys_num;
    fdps_create_psys(&psys_num,"full_particle");
    fdps_init_psys(psys_num);
    // Create and initialize tree object
    int tree_num;
    fdps_create_tree(&tree_num, 
                     "Long,full_particle,full_particle,full_particle,Monopole");
    int ntot=1024;
    double theta = 0.5;
    int n_leaf_limit = 8;
    int n_group_limit = 64;
    fdps_init_tree(tree_num, ntot, theta, n_leaf_limit, n_group_limit);
    // Make an initial condition
    setup_IC(psys_num,ntot);
    // Domain decomposition and exchange particle
    fdps_decompose_domain_all(dinfo_num,psys_num,-1.0);
    fdps_exchange_particle(psys_num,dinfo_num);
    // Compute force at the initial time
    fdps_calc_force_all_and_write_back(tree_num, 
                                       calc_gravity_ep_ep,
                                       calc_gravity_ep_sp,
                                       psys_num,
                                       dinfo_num,
                                       true,
                                       FDPS_MAKE_LIST);
    //dump_particles(psys_num);
    // Compute energies at the initial time
    double etot0,ekin0,epot0;
    calc_energy(psys_num, &etot0, &ekin0,&epot0);
    printf("Energies = %21.14e  %21.14e  %21.14e\n",etot0,ekin0,epot0);
    // Time integration
    double time_diag = 0;
    double time_snap = 0;
    double time_sys  = 0;
    double time_end = 10.0;
    double dt = 1.0/128.0;
    double dt_diag = 1.0;
    double dt_snap = 1.0;
    int    num_loop = 0;
    while (time_sys <= time_end){
        if (time_sys + dt/2 >= time_snap){
            // output(psys_num)
            time_snap += dt_snap;
        }
        double etot1, ekin1, epot1;
        calc_energy(psys_num, &etot1,&ekin1,&epot1);
        //printf( "Energies = %21.14e  %21.14e  %21.14e\n",etot1,ekin1,epot1);
        //dump_particles(psys_num);
        if (fdps_get_rank() == 0){
            if (time_sys + dt/2 >= time_diag){
                printf ("time: %10.3f, energy error: %15.7e\n",
                        time_sys, (etot1-etot0)/etot0);
                time_diag = time_diag + dt_diag;
            }
        }
        kick(psys_num,0.5*dt);
        time_sys +=  dt;
        drift(psys_num,dt);
        // Domain decomposition & exchange particle
        if (num_loop %4  == 0) {
            fdps_decompose_domain_all(dinfo_num,psys_num,-1.0);
        }
        fdps_exchange_particle(psys_num,dinfo_num);
        // Force calculation
        fdps_calc_force_all_and_write_back(tree_num,
                                           calc_gravity_ep_ep,
                                           calc_gravity_ep_sp,
                                           psys_num,
                                           dinfo_num,
                                           true,
                                           FDPS_MAKE_LIST);
        kick(psys_num,0.5*dt);
        num_loop += 1;
    }
    fdps_finalize();
    return 0;
}
