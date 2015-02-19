#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "treepm_header.h"
#include "pp.h"

namespace ParticleSimulator{
    namespace ParticleMesh{


static char cbuf[CHARMAX];

#define MAXIT 100
double rtsafe(double x1, double x2, double xacc, double tau){

  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;

  funcd(x1,&fl,&df,tau);
  funcd(x2,&fh,&df,tau);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    fprintf(stderr,"Root must be bracketed in rtsafe\n");
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  funcd(rts,&f,&df,tau);
  for (j=1;j<=MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
	|| (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    funcd(rts,&f,&df,tau);
    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
  }
  fprintf(stderr,"Maximum number of iterations exceeded in rtsafee\n");
  return 0.0;
}
#undef MAXIT



double timetoz( double tnow, struct cosmology cosm){

  double t1 = 1.e-2;
  double t2 = 5.e0;
  double xacc = 1.e-6;

  if(cosm.omega0 > 0.99 && cosm.lambda0 < 0.01){
    double zred = (double)(pow(tnow*3.0/2.0,-2.0/3.0) - 1.0);
    return zred;
  }

  if(cosm.omega0 < 1.e0 && cosm.lambda0 < 0.01){
    double tau = tnow*2.e0*pow(1.e0-cosm.omega0,3.e0/2.e0)/cosm.omega0;
    double theta = rtsafe(t1,t2,xacc,tau);
    double zred = (double)(2.0*(1.0/cosm.omega0-1.0)/(5.e-1*(exp(theta)+exp(-theta))-1.0)-1.0);
    return zred;
  }

  if(cosm.omega0 < 1.e0 && cosm.lambda0 > 0.01){
    double tau = tnow*3.0*sqrt(1.0-cosm.omega0);
    double theta = (exp(tau)-1.0)/(exp(tau)+1.0);
    double zred = (1.0/cosm.omega0-1.0)*(1.0/(theta*theta)-1.0);
    zred = pow(zred,1.e0/3.e0)-1.e0;
    return zred;
  }
   
   fprintf(stderr,"Errors in timetoz(..)\n");
   exit(EXIT_FAILURE);
}



void funcd(double x, double *f, double *df, double tau){

  *f  = 0.5*(exp(x)-exp(-x)) - x - tau;
  *df = 0.5*(exp(x)+exp(-x)) - 1.e0;

  return;

}



double zToTime( const double znow, cosmology &cosm){

  double om = cosm.omega0;
  double ov = cosm.lambda0;

  double ov2 = sqrt(ov);
  double znow3 = 1.0 + znow;
  znow3 = znow3 * znow3 * znow3;

  double x = sqrt(om*znow3 + ov);

  return -0.333333 * log( (x-ov2)/(x+ov2)) / ov2;


}



/* update redshift, scale factor, hubble */
void update_now( pRunParam this_run){

  double om = this_run->cosm.omega0;
  double ov = this_run->cosm.lambda0;

  double zred = timetoz(this_run->tnow, this_run->cosm);

  double ascale = 1.e0/(1.e0+zred);

  this_run->znow = zred;
  this_run->anow = ascale;
  this_run->hnow = sqrt(1.e0+om*(1.e0/ascale-1.e0)
			+ov*(ascale*ascale-1.e0))/ascale;

}



void getIntegralFactor( const pRunParam this_run, const double dt,
			double *vfac, double *afac){

  double om = this_run->cosm.omega0;
  double ov = this_run->cosm.lambda0;
  double ascale = this_run->anow;

  double at = sqrt(1.0+om*(1.0/ascale-1.0)+ov*(ascale*ascale-1.0))/(ascale);
  double bt = 1.0/(ascale*ascale*ascale);
  double atdt1 = 1.0+at*dt;
  (*vfac)  = (2.0-atdt1)/atdt1;
  (*afac)  = bt*dt/atdt1;

}



void kick( pParticle particle, const float *a, 
		   const double vfac, const double afac){

  particle->xvel = vfac*particle->xvel + afac*a[0];
  particle->yvel = vfac*particle->yvel + afac*a[1];
  particle->zvel = vfac*particle->zvel + afac*a[2];

}



void kick( pParticle particle, const float *a, 
		   const double afac){

  particle->xvel += afac*a[0];
  particle->yvel += afac*a[1];
  particle->zvel += afac*a[2];

}




void drift( pParticle particle, const double dt){

  particle->xpos += particle->xvel * dt;
  particle->ypos += particle->yvel * dt;
  particle->zpos += particle->zvel * dt;

  if(particle->xpos < 0.0) particle->xpos += 1.0;
  if(particle->ypos < 0.0) particle->ypos += 1.0;
  if(particle->zpos < 0.0) particle->zpos += 1.0;
  if(particle->xpos >= 1.0) particle->xpos -= 1.0;
  if(particle->ypos >= 1.0) particle->ypos -= 1.0;
  if(particle->zpos >= 1.0) particle->zpos -= 1.0;

}



void drift( pParticle p, const int n, const double dt){

#pragma omp parallel for
  for( int i=0; i<n; i++){
    p[i].xpos += p[i].xvel * dt;
    p[i].ypos += p[i].yvel * dt;
    p[i].zpos += p[i].zvel * dt;

    if( p[i].xpos < 0.0) p[i].xpos += 1.0;
    if( p[i].ypos < 0.0) p[i].ypos += 1.0;
    if( p[i].zpos < 0.0) p[i].zpos += 1.0;
    if( p[i].xpos >= 1.0) p[i].xpos -= 1.0;
    if( p[i].ypos >= 1.0) p[i].ypos -= 1.0;
    if( p[i].zpos >= 1.0) p[i].zpos -= 1.0;
  }

}



void step_pos( pParticle particle, pRunParam this_run, double dtime){

#pragma omp parallel for
  for( int i=0; i<this_run->npart; i++) {

    particle[i].xpos += particle[i].xvel*dtime;
    particle[i].ypos += particle[i].yvel*dtime;
    particle[i].zpos += particle[i].zvel*dtime;

    if(particle[i].xpos >= 1.0) particle[i].xpos -= 1.0;
    if(particle[i].ypos >= 1.0) particle[i].ypos -= 1.0;
    if(particle[i].zpos >= 1.0) particle[i].zpos -= 1.0;
    if(particle[i].xpos < 0.0) particle[i].xpos += 1.0;
    if(particle[i].ypos < 0.0) particle[i].ypos += 1.0;
    if(particle[i].zpos < 0.0) particle[i].zpos += 1.0;
  }

  this_run->tnow += dtime;
  update_now(this_run);

}



double getEps( const double znow, const double anow){

  double epstmp;
  if( znow < CONST_Z_VALUE){
    epstmp = SFT_FOR_PP/anow;
  }
  else{
    epstmp = (1+CONST_Z_VALUE) * SFT_FOR_PP;
  }


  return epstmp;

}



double getTheta(const double znow){

  double theta;

  if( znow < Z_SWITCH_THETA){
    theta = THETA;
  }
  else{
    theta = THETA_HIGHZ;
  }

  return theta;


}



void correctBoundaryCondition( pParticle particle, const int n){

#pragma omp parallel for
  for( int ii=0; ii<n; ii++){
    if( particle[ii].xpos <  0.0)  particle[ii].xpos += 1.0;
    if( particle[ii].xpos >= 1.0)  particle[ii].xpos -= 1.0;
    if( particle[ii].ypos <  0.0)  particle[ii].ypos += 1.0;
    if( particle[ii].ypos >= 1.0)  particle[ii].ypos -= 1.0;
    if( particle[ii].zpos <  0.0)  particle[ii].zpos += 1.0;
    if( particle[ii].zpos >= 1.0)  particle[ii].zpos -= 1.0;
  }

}



void TreePMTimePrint(const TreePMTime t, pRunParam this_run){

  fprintf( this_run->profile_file,"------ step %d ------\n",this_run->nstep-1);

  fprintf( this_run->profile_file, 
	   "node: %d\n",
	   this_run->inode);

  fprintf( this_run->profile_file, "global_boundary: %e %e %e %e %e %e\n",
	   this_run->global_bmin[0], this_run->global_bmax[0],
	   this_run->global_bmin[1], this_run->global_bmax[1],
	   this_run->global_bmin[2], this_run->global_bmax[2]);

  fprintf( this_run->profile_file, "boundary: %e %e %e %e %e %e\n",
	   this_run->bmin[0], this_run->bmax[0],
	   this_run->bmin[1], this_run->bmax[1],
	   this_run->bmin[2], this_run->bmax[2]);

  fprintf( this_run->profile_file, 
	   "n_total: %lld\n",
	   this_run->npart_total);

  fprintf( this_run->profile_file, 
	   "n: %d\n",
	   this_run->npart_old);

  fprintf( this_run->profile_file, 
	   "bn: %d\n",
	   this_run->bn);

  fprintf( this_run->profile_file, 
	   "nwalk: %d\n",
	   this_run->nwalk);

  fprintf( this_run->profile_file, 
	   "nisum: %d\n",
	   this_run->nisum);

  fprintf( this_run->profile_file, 
	   "nimin: %d\n",
	   this_run->nimin);

  fprintf( this_run->profile_file, 
	   "nimax: %d\n",
	   this_run->nimax);

  fprintf( this_run->profile_file, 
	   "njsum: %d\n",
	   this_run->njsum);

  fprintf( this_run->profile_file, 
	   "njmin: %d\n",
	   this_run->njmin);

  fprintf( this_run->profile_file, 
	   "njmax: %d\n",
	   this_run->njmax);

  fprintf( this_run->profile_file, 
	   "ninteraction: %lld\n",
	   this_run->ninteraction);

  fprintf( this_run->profile_file, 
	   "ninteraction_sum: %.15e\n",
	   this_run->ninter_sum);

  fprintf( this_run->profile_file, 
	   "load_balance: %f / %f\n\n",
	   this_run->nrate, 1.0/(double)this_run->nnode);

  fprintf( this_run->profile_file, 
	   "CPUTIME_for_this_step: %lf\n",
	   t.cputime);

  fprintf( this_run->profile_file, 
	   "  PM: %lf\n",
	   t.pm);

  fprintf( this_run->profile_file, 
	   "    PM_Density_Assignment: %lf\n",
	   t.pm_density_assignment);

  fprintf( this_run->profile_file, 
	   "    PM_Make_Comm_Density: %lf\n",
	   t.pm_make_comm_density);

  fprintf( this_run->profile_file, 
	   "    PM_Density_Communication: %lf\n",
	   t.pm_density_comm);

#ifdef RMM_PM
  fprintf( this_run->profile_file, 
	   "    PM_Density_Communication_Reduce: %lf\n",
	   t.pm_density_comm_reduce);
#endif

  fprintf( this_run->profile_file, 
	   "    PM_FFT_F: %lf\n",
	   t.pm_fft_forward);

  fprintf( this_run->profile_file, 
	   "    PM_FFT_B: %lf\n",
	   t.pm_fft_backward);

  fprintf( this_run->profile_file, 
	   "    PM_Density->Phi: %lf\n",
	   t.pm_density_phi);

  fprintf( this_run->profile_file, 
	   "    PM_Make_Comm_Phi: %lf\n",
	   t.pm_make_comm_phi);

  fprintf( this_run->profile_file, 
	   "    PM_Phi_Communication: %lf\n",
	   t.pm_phi_comm);

#ifdef RMM_PM
  fprintf( this_run->profile_file, 
	   "    PM_Phi_Communication_Bcast: %lf\n",
	   t.pm_phi_comm_bcast);
#endif

  fprintf( this_run->profile_file, 
	   "    PM_Mesh_Force: %lf\n",
	   t.pm_mesh_force);

  fprintf( this_run->profile_file, 
	   "    PM_Interpolation+Kick: %lf\n",
	   t.pm_interpolation);

  fprintf( this_run->profile_file, 
	   "  PP_and_Drift: %lf\n",
	   t.pp_drift);

  fprintf( this_run->profile_file, 
	   "    PP_on_this_node: %lf\n",
	   t.pp);

  fprintf( this_run->profile_file, 
	   "      PP_Get_Boundary_Particle: %lf\n",
	   t.pp_get_boundary);

  fprintf( this_run->profile_file, 
	   "        Get_Boundary_Determine_Comm_Node: %lf\n",
	   t.pp_get_boundary_determine_comm_node);

  fprintf( this_run->profile_file, 
	   "        Get_Boundary_Make_Key_Sort: %lf\n",
	   t.pp_get_boundary_make_key_sort);

  fprintf( this_run->profile_file, 
	   "        Get_Boundary_Make_Tree: %lf\n",
	   t.pp_get_boundary_construct_tree);

  fprintf( this_run->profile_file, 
	   "        Get_Boundary_Make_Send_List: %lf\n",
	   t.pp_get_boundary_make_send_list);

  fprintf( this_run->profile_file, 
	   "        Get_Boundary_Comm: %lf\n",
	   t.pp_get_boundary_comm);

  fprintf( this_run->profile_file, 
	   "        Get_Boundary_Push: %lf\n",
	   t.pp_get_boundary_push);

  fprintf( this_run->profile_file, 
	   "      PP_Make_Key: %lf\n",
	   t.pp_make_key);

  fprintf( this_run->profile_file, 
	   "      PP_Particle_Sort: %lf\n",
	   t.pp_particle_sort);

  fprintf( this_run->profile_file, 
	   "      PP_Key_Sort: %lf\n",
	   t.pp_key_sort);

  fprintf( this_run->profile_file, 
	   "      PP_Construct_Tree: %lf\n",
	   t.pp_construct_tree);

  fprintf( this_run->profile_file, 
	   "      PP_Make_Interaction_List: %lf\n",
	   t.pp_make_intr_list);

  fprintf( this_run->profile_file, 
	   "      PP_Calc_Force: %lf\n",
	   t.pp_calc_force);

  fprintf( this_run->profile_file, 
	   "      PP_Vel_Kick: %lf\n",
	   t.pp_vel_kick);

  fprintf( this_run->profile_file, 
	   "    Drift: %lf\n",
	   t.drift);

#ifdef MULTI_TIMESTEP
  fprintf( this_run->profile_file, 
	   "  Next_Dt: %lf\n",
	   t.next_dt);
#else
  fprintf( this_run->profile_file, 
	   "    Next_Dt: %lf\n",
	   t.next_dt);
#endif

  fprintf( this_run->profile_file, 
	   "  Exchange: %lf\n",
	   t.exchange);

  fprintf( this_run->profile_file, 
	   "    Exchange_Determine_Boundary: %lf\n",
	   t.exchange_determine_boundary);

  fprintf( this_run->profile_file, 
	   "      Exchange_Determine_Boundary_Comm: %lf\n",
	   t.exchange_determine_boundary_comm);

  fprintf( this_run->profile_file, 
	   "    Exchange_Make_Send_List_Pre: %lf\n",
	   t.exchange_make_send_list_pre);

  fprintf( this_run->profile_file, 
	   "    Exchange_Make_Send_List: %lf\n",
	   t.exchange_make_send_list);

  fprintf( this_run->profile_file, 
	   "    Exchange_Comm: %lf\n",
	   t.exchange_comm);

  fprintf( this_run->profile_file, 
	   "  IO_Image: %lf\n",
	   t.io_image);

  fprintf( this_run->profile_file, 
	   "  IO_Density: %lf\n",
	   t.io_density);

  fprintf( this_run->profile_file, 
	   "  IO_Snapshot: %lf\n",
	   t.io_snapshot);

  fprintf( this_run->profile_file, 
	   "  LOAD_BALANCE_LOG: %lf\n\n",
	   t.load_balance_log);

  memset( &(this_run->t), 0, sizeof(TreePMTime));

}



void logPrint( const RunParam *this_run, const double eps, const double global_time){


#ifdef KCOMPUTER
  double perf = this_run->nintrave * 51.0 / 128e9 / this_run->t.cputime;
  double perf2 = this_run->nintrdisp * 51.0 / 128e9 / this_run->t.cputime;
  double perf3 = this_run->ninter_sum * 51.0 / 128e9 / global_time / this_run->nnode;
  fprintf( this_run->log_file,"  %5d\t%8.4e\t%8.4e\t%5.4f\t%5.4e\t%e\t%5.4f\t%8.4e\t%e\t%e\t%e\n", 
	   this_run->nstep, this_run->dtime, this_run->tnow, this_run->znow, 
	   this_run->anow, eps, 
	   this_run->t.cputime, global_time, perf, perf2, perf3);
#else
  fprintf( this_run->log_file,"  %5d\t%8.4e\t%8.4e\t%5.4f\t%5.4e\t%e\t%5.4f\t%8.4e\n", 
	   this_run->nstep, this_run->dtime, this_run->tnow, this_run->znow, 
	   this_run->anow, eps, this_run->t.cputime, global_time);
#endif

}



void initRunParam( pRunParam this_run, const int inode, const int nnode,
		   const int nowstep){

  memset( this_run, 0, sizeof(RunParam));

  /* init run_param */
  this_run->inode = inode;
  this_run->nnode = nnode;
  this_run->nrate = 1.0 / (double)nnode;
  this_run->global_bmin[0] = 0.0;
  this_run->global_bmin[1] = 0.0;
  this_run->global_bmin[2] = 0.0;
  this_run->global_bmax[0] = 1.0;
  this_run->global_bmax[1] = 1.0;
  this_run->global_bmax[2] = 1.0;
  this_run->other_bmin[0] = 0.0;
  this_run->other_bmin[1] = 0.0;
  this_run->other_bmin[2] = 0.0;
  this_run->other_bmax[0] = 1.0;
  this_run->other_bmax[1] = 1.0;
  this_run->other_bmax[2] = 1.0;
  this_run->MPI_COMM_INTERNAL = MPI_COMM_WORLD;

  //log, profile setting
  char log_name[CHARMAX], profile_name[CHARMAX];
  sprintf( log_name,"%s.log%d", MODEL, inode);
  sprintf( profile_name,"%s.out%d", MODEL, inode);
  if( ONELOGFLAG == 1){
    if( inode != 0){
      sprintf( log_name,"/dev/null");
      sprintf( profile_name,"/dev/null");
    }
  }
  this_run->log_file     = fopen((const char*)log_name, "a");
  this_run->profile_file = fopen((const char*)profile_name, "a");

  if( this_run->inode == 0){
    if( LOADBALANCELOG_ON == 1){
      sprintf( cbuf, "%s.loadbalance", MODEL);
      this_run->loadbalance_file = fopen( cbuf, "a");
    }
    if( BOUNDARYLOG_ON == 1){
      sprintf( cbuf, "%s.bound", MODEL);
      this_run->boundary_file = fopen( cbuf, "a");
    }
  }

  this_run->dtime = 0.0;
  this_run->dtime_ini = 1.0e30;
  this_run->nstep = 0;

  this_run->ninter_sum = 0.0;

  /* g5_setup */
#ifndef GRAPE_OFF
  g5_open();
  g5_set_range( 0.0, 1.0, this_run->min_mass);
  //g5_set_cutoff_table(force_cutoff_func,1.9,1.0,pot_cutoff_func,1.8,1.0);  // for g6a
  double eta = SFT_FOR_PM * 0.5;
  g5_set_eta(eta);
#endif

}



void genCommYZ( pRunParam this_run){

  int x, y, z;
  int *ndiv = this_run->ndiv;
  getXYZIndex( this_run->inode, ndiv, &x, &y, &z);
  int ndiv_ic[3] = { 1, ndiv[1], ndiv[2]};
  int color = x;
  int key = getVoxelIndex( 0, y, z, ndiv_ic);
  MPI_Comm_split( MPI_COMM_WORLD, color, key, &this_run->MPI_COMM_YZ);
  this_run->inode_yz = key;

}



void printLoadBalance( RunParam *this_run, FILE *fout){

#define NDATA_STAT (24)
  double dnpart = this_run->npart;
  double dbn = this_run->bn;
  double dnint = this_run->ninteraction;
  static char cdata[NDATA_STAT][CHARMAX] = {
    "npart", "bnpart", "nintr", 
    "ncomm_tree", "nsend_tree",
    "nscomm_ex", "nsend_ex",
    "nrcomm_ex", "nrecv_ex",
    "nscomm_pm_l2s", "nsend_pm_l2s",
    "nrcomm_pm_l2s", "nrrecv_pm_l2s",
    "nscomm_pm_s2l", "nsend_pm_s2l",
    "nrcomm_pm_s2l", "nrrecv_pm_s2l",
    "t_PM", "t_PP", "t_boundary", 
    "t_makeInt", "t_force", "t_exchange", "t_PPdrift"
  };

  double data[NDATA_STAT] = { (double)dnpart, (double)dbn, (double)dnint, 
			      (double)this_run->ncomm_tree, (double)this_run->nsend_tree,
			      (double)this_run->nscomm_exchange, (double)this_run->nsend_ex,
			      (double)this_run->nrcomm_exchange, (double)this_run->nrecv_ex,
			      (double)this_run->nscomm_pm_l2s, (double)this_run->nsend_pm_l2s,
			      (double)this_run->nrcomm_pm_l2s, (double)this_run->nrecv_pm_l2s,
			      (double)this_run->nscomm_pm_s2l, (double)this_run->nsend_pm_s2l,
			      (double)this_run->nrcomm_pm_s2l, (double)this_run->nrecv_pm_s2l,
			      this_run->t.pm, this_run->t.pp, this_run->t.pp_get_boundary,
			      this_run->t.pp_make_intr_list, this_run->t.pp_calc_force,
			      this_run->t.exchange, this_run->t.pp_drift};
  double ninv = 1.0 / (double)this_run->nnode;

  if( this_run->inode == 0){
    fprintf( fout, "------ step %d  (ave, disp, min, max)------\n", this_run->nstep);
  }

  for( int i=0; i<NDATA_STAT; i++){
    double min = 0.0;
    double max = 0.0;
    double sum = 0.0;
    double sum2 = 0.0;
    double data2 = data[i] * data[i];
    MPI_Reduce( data+i,  &min, 1, MPI_DOUBLE, MPI_MIN, 0, this_run->MPI_COMM_INTERNAL);
    MPI_Reduce( data+i,  &max, 1, MPI_DOUBLE, MPI_MAX, 0, this_run->MPI_COMM_INTERNAL);
    MPI_Reduce( data+i,  &sum, 1, MPI_DOUBLE, MPI_SUM, 0, this_run->MPI_COMM_INTERNAL);
    MPI_Reduce( &data2, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, this_run->MPI_COMM_INTERNAL);
    double ave = sum * ninv;
    double disp = sum2*ninv - ave*ave;
    if( disp < 1.0e-20)  disp = 0.0;
    else disp = sqrt(disp);
    if( this_run->inode != 0)  continue;
    fprintf( fout, "%-15s: %e\t%e\t%e\t%e\n", 
	     cdata[i], ave, disp, min, max);
    if( i==2){
      this_run->nintrave = ave;
      this_run->nintrdisp = disp;
      this_run->ninter_sum += sum;
    }
  }

}



void printBoundary( const RunParam *this_run, FILE *fout,
		    const double *bmin, const double *bmax){

  sprintf( cbuf, "%e\t%e\t%e\t%e\t%e\t%e\t%d\n",
	   bmin[0], bmax[0], 
	   bmin[1], bmax[1], 
	   bmin[2], bmax[2], 
	   this_run->npart);
  mp_print( cbuf, this_run->MPI_COMM_INTERNAL, fout);

}

    } // namespace ParticleMesh
}     // namespace ParticleSimulator
