#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "user_defined.h"
#include "kernel.h"

void dump_fullp(particle p)
{
    printf("%lld %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e",
           p.id, p.mass, p.pos.x, p.pos.y, p.pos.z,
           p.vel.x, p.vel.y, p.vel.z);
    printf("%15.7e %15.7e %15.7e %15.7e\n",
           p.acc.x, p.acc.y, p.acc.z, p.pot);
}
void dump_particles(particle *p, int n)
{
    int i;
    for (i=0;i<n;i++)dump_fullp(p[i]);
}

void setup_IC(particle *ptcl,
              int nptcl)
{
  double m_tot=1.0;
  double rmax=3.0;
  double r2max=rmax*rmax;
  int i;
  for (i=0; i < nptcl; i++){
    particle *q = ptcl+i;
    q->id = i;
    q->mass = m_tot/nptcl;
    double r2 = r2max*2;
    pikg_f64vec pos;
    while (r2 >= r2max){
      pos.x= (2*((double)rand()/(double)RAND_MAX)) * rmax;
      pos.y= (2*((double)rand()/(double)RAND_MAX)) * rmax;
      pos.z= (2*((double)rand()/(double)RAND_MAX)) * rmax;
      r2 = pos.x*pos.x + pos.y*pos.y + pos.z*pos.z;
    }
    q->pos = pos;
    q->vel.x = 0.0;
    q->vel.y = 0.0;
    q->vel.z = 0.0;
    q->eps = 1.0/32.0;
  }
  pikg_f64vec cm_pos;
  pikg_f64vec cm_vel;
  cm_pos.x = 0.0; cm_pos.y = 0.0; cm_pos.z = 0.0;
  cm_vel.x = 0.0; cm_vel.y = 0.0; cm_vel.z = 0.0;
  double cm_mass = 0;
  for (i=0; i < nptcl; i++){
    particle *pi = ptcl+i;
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
  for (i=0; i < nptcl; i++){
    particle* q = ptcl+i;
    q->pos.x -= cm_pos.x;
    q->pos.y -= cm_pos.y;
    q->pos.z -= cm_pos.z;
    q->vel.x -= cm_vel.x;
    q->vel.y -= cm_vel.y;
    q->vel.z -= cm_vel.z;
  }
  //dump_fullpsys(ptcl, nptcl);
}

void calc_energy(particle* ptcl,
		 int nptcl,
                 double *etot,
                 double *ekin,
                 double *epot)
{
    *etot = *ekin = *epot = 0.0;
    int i;
    for (i=0;i < nptcl; i++){
      particle *pi = ptcl+i;
      pikg_f64vec v = pi->vel;
      *ekin += pi->mass * (v.x*v.x+v.y*v.y+v.z*v.z);
      *epot += pi->mass * (pi->pot + pi->mass/pi->eps);
    }
    *ekin *= 0.5;
    *epot *= 0.5;
    *etot = *ekin + *epot;
}

void  kick(particle* ptcl,int n, double dt)
{
  int i;
  for (i=0;i < n; i++){
    particle *pi = ptcl+i;
    pikg_f64vec *pv, *pa;
    pv = &(pi->vel);
    pa = &(pi->acc);
    pv->x += pa->x * dt;
    pv->y += pa->y * dt;
    pv->z += pa->z * dt;
  }
}

void  drift(particle* ptcl,int n, double dt)
{
  int i;
  for (i=0;i < n; i++){
    particle *pi = ptcl+i;
    pikg_f64vec *px, *pv;
    pv = &(pi->vel);
    px = &(pi->pos);
    px->x += pv->x * dt;
    px->y += pv->y * dt;
    px->z += pv->z * dt;
  }
}

void clear_ptcl(particle* ptcl,int n)
{
  int i;
  for(i=0;i<n;i++){
    particle *pi = ptcl + i;
    pi->acc.x = 0.0;
    pi->acc.y = 0.0;
    pi->acc.z = 0.0;
    pi->pot   = 0.0;
  }
}

int main()
{
  fprintf(stderr,"PIKG sapmple of N-body\n");
  int ntot=1024;
  particle* ptcl = (particle*)malloc(sizeof(particle)*ntot);
  if(ptcl == NULL){
    fprintf(stderr,"error: malloc ptcl failed\n");
    exit(EXIT_FAILURE);
  }
  setup_IC(ptcl,ntot);
  clear_ptcl(ptcl,ntot);
  calc_gravity(ptcl,ntot,ptcl,ntot,ptcl);

  //dump_particles(psys_num);
  // Compute energies at the initial time
  double etot0,ekin0,epot0;
  calc_energy(ptcl,ntot, &etot0, &ekin0,&epot0);
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
      // output(ptcl)
      time_snap += dt_snap;
    }
    double etot1, ekin1, epot1;
    calc_energy(ptcl,ntot, &etot1,&ekin1,&epot1);
    //printf( "Energies = %21.14e  %21.14e  %21.14e\n",etot1,ekin1,epot1);
    //dump_particles(ptcl);
    if (time_sys + dt/2 >= time_diag){
      printf ("time: %10.3f, energy error: %15.7e\n",
	      time_sys, (etot1-etot0)/etot0);
      time_diag = time_diag + dt_diag;
    }
    kick(ptcl,ntot,0.5*dt);
    time_sys +=  dt;
    drift(ptcl,ntot,dt);
    clear_ptcl(ptcl,ntot);
    calc_gravity(ptcl,ntot,ptcl,ntot,ptcl);
    kick(ptcl,ntot,0.5*dt);
    num_loop += 1;
  }
  if(ptcl != NULL) free(ptcl);
  return 0;
}
