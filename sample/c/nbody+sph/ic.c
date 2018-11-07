#include "ic.h"

void read_tipsy_file(char file_name[], int psys_num) {
   // This function is used to read particle data created by
   // MAGI (https://bitbucket.org/ymiki/magi). The particle
   // data must be in the TIPSY format.
   FILE *fp;
   if ((fp = fopen(file_name,"rb")) == NULL) {
      fprintf(stderr,"Cannot open file %s\n.",file_name);
   }
   Magi_tipsy_header header;
   fread(&header,sizeof(Magi_tipsy_header),1,fp);
   fprintf(stdout,"nbodies = %d\n",header.nbodies);
   fdps_set_nptcl_loc(psys_num,header.nbodies);
   FP_nbody *ptcl = (FP_nbody *) fdps_get_psys_cptr(psys_num);
   int i;
   for (i = 0; i < header.nbodies; i++) {
      Magi_tipsy_particle ptcl_tipsy;
      fread(&ptcl_tipsy,sizeof(Magi_tipsy_particle),1,fp);
      ptcl[i].mass  = ptcl_tipsy.mass;
      ptcl[i].pos.x = ptcl_tipsy.pos[0];
      ptcl[i].pos.y = ptcl_tipsy.pos[1];
      ptcl[i].pos.z = ptcl_tipsy.pos[2];
      ptcl[i].vel.x = ptcl_tipsy.vel[0];
      ptcl[i].vel.y = ptcl_tipsy.vel[1];
      ptcl[i].vel.z = ptcl_tipsy.vel[2];
   }
   fclose(fp);
}

void galaxy_IC(int psys_num_nbody,
               int psys_num_sph,
               int *bc,
               fdps_f32vec *pos_root_domain_low,
               fdps_f32vec *pos_root_domain_high,
               double *time_dump,
               double *dt_dump,
               double *time_end) {
   // Local parameters
   //- Definitions of the code units of MAGI
   //    [Important]
   //    (1) The values MUST BE consistent with "the computational units"
   //        written in the file ./magi_data/doc/unit.txt, which is output
   //        by MAGI when we create a particle data with MAGI.
   //    (2) The MAGI's code units are DIFFERENT for unit systems
   //        a user choose in the file ./magi_data/cfg/Galaxy.tipsy.
   //        For detail, read Section "Unit systems in inc/constants.[c h]"
   //        in https://bitbucket.org/ymiki/magi.
   //    (3) In this sample code, "Galactic scale" unit is adopted.
   //        It is consistent with ./magi_data/cfg/Galaxy.cfg.
   double magi_unit_mass = 1.0e8 * Msolar;
   double magi_unit_leng = kpc;
   double magi_unit_time = 1.0e2 * Myr;
   double magi_unit_velc = magi_unit_leng/magi_unit_time;
   //- Definitions of the model parameters for a gaseous exponential disk
   int nptcl_sph = pow(2,18);
   double m_gas = 1.0e10 * Msolar;
   double Rs = 7.0e0 * kpc; // scale radius
   double Rt = 12.5e0 * kpc; // truncation radius
   double zd = 4.0e2 * pc; // scale height
   double zt = 1.0e0 * kpc; // truncation height
   double temp = 1.0e4; // gas temperature 
   double mu = 0.5e0; // mean molecular weight relative to the mass of hydrogen
   double eps = 1.0e-6; // iteraction accuracy 
   //- Definitions of parameters for output
   int nbin=64;
   double safety=1.001e0;

   // Initialize pseudorandom number generator
   int mtts_num;
   fdps_create_mtts(&mtts_num);
   fdps_mtts_init_genrand(mtts_num,0);
   // Place Nbody particles
   char file_name[64] = {'\0'};
   strcpy(file_name,"./magi_data/dat/Galaxy.tipsy");
   read_tipsy_file(file_name, psys_num_nbody);
   FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
   int nptcl_nbody = fdps_get_nptcl_loc(psys_num_nbody);
   int i;
   for (i = 0; i < nptcl_nbody; i++) {
      ptcl_nbody[i].id    = i;
      ptcl_nbody[i].acc.x = 0.0;
      ptcl_nbody[i].acc.y = 0.0;
      ptcl_nbody[i].acc.z = 0.0;
      ptcl_nbody[i].pot   = 0.0;
   }
   // Place SPH particles
   fdps_set_nptcl_loc(psys_num_sph,nptcl_sph);
   FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
   for (i = 0; i < nptcl_sph; i++) {
      // First make a uniform disk with a finite thickness
      double x,y,z,r,r2;
      for (;;) { 
         x = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * Rt;
         y = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * Rt;
         z = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * zt;
         r2 = (x * x + y * y);
         if ((r2 < Rt * Rt) && (fabs(z) < zt)) break;
      }
      r = sqrt(r2);
      // Then re-scale particle position to generate an exponential disk
      double r_low = 0.0;
      double r_high = Rt;
      double val, val_trgt = (r/Rt) * (r/Rt);
      double r_new;
      for (;;) {
         r_new = 0.5 * (r_low + r_high);
         val = (1.0 - (r_new/Rs + 1.0)*exp(-r_new/Rs))
             / (1.0 - (   Rt/Rs + 1.0)*exp(-   Rt/Rs));
         if (val < val_trgt) r_low  = r_new;
         if (val > val_trgt) r_high = r_new;
         double reldiff = 2.0 * fabs(r_low - r_high)/(r_low + r_high);
         if (reldiff < eps) {
            r_new = 0.5 * (r_low + r_high);
            break;
         }
      }
      double z_new;
      if (z >= 0.0) {
         z_new = - zd * log(1.0 - (z/zt) * (1.0 - exp(-zt/zd)));
      } else {
         z_new = zd * log(1.0 + (z/zt) * (1.0 - exp(-zt/zd)));
      }
      // Set  
      ptcl_sph[i].id   = i + nptcl_nbody;
      ptcl_sph[i].mass = m_gas / nptcl_sph;
      ptcl_sph[i].pos.x = (r_new / r) * x;
      ptcl_sph[i].pos.y = (r_new / r) * y;
      ptcl_sph[i].pos.z = z_new;
      //printf("%15.7e %15.7e\n",ptcl_sph[i].pos.x,ptcl_sph[i].pos.y);
      ptcl_sph[i].vel.x = 0.0;
      ptcl_sph[i].vel.y = 0.0;
      ptcl_sph[i].vel.z = 0.0;
      ptcl_sph[i].acc_grav.x = 0.0;
      ptcl_sph[i].acc_grav.y = 0.0;
      ptcl_sph[i].acc_grav.z = 0.0;
      ptcl_sph[i].pot_grav = 0.0;
      ptcl_sph[i].acc_hydro.x = 0.0;
      ptcl_sph[i].acc_hydro.y = 0.0;
      ptcl_sph[i].acc_hydro.z = 0.0;
      ptcl_sph[i].eng = (kBoltz * temp)/((specific_heat_ratio - 1.0) * mu * Mhydrogen);
      ptcl_sph[i].smth = pow(Rt*Rt*zt,1.0/3.0) * pow(((double)N_neighbor)/((double)nptcl_sph),1.0/3.0);
   }
   // Unit convertion (MAGI unit, CGS unit -> G=M=R=1 system)
   double unit_mass = 0.0;
   double unit_leng = 0.0;
   for (i = 0; i < nptcl_nbody; i++) {
      ptcl_nbody[i].mass *= magi_unit_mass;
      ptcl_nbody[i].pos.x *= magi_unit_leng;
      ptcl_nbody[i].pos.y *= magi_unit_leng;
      ptcl_nbody[i].pos.z *= magi_unit_leng;
      ptcl_nbody[i].vel.x *= magi_unit_velc;
      ptcl_nbody[i].vel.y *= magi_unit_velc;
      ptcl_nbody[i].vel.z *= magi_unit_velc;
      unit_mass += ptcl_nbody[i].mass;
      fdps_f64vec pos = ptcl_nbody[i].pos;
      double r = sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
      if (r > unit_leng) unit_leng = r;
   }
   printf("Total mass in N-body particles = %15.7e [Msolar]\n",unit_mass/Msolar);
   for (i = 0; i < nptcl_sph; i++) {
      unit_mass += ptcl_sph[i].mass;
      fdps_f64vec *pos = &ptcl_nbody[i].pos;
      double r = sqrt(pos->x * pos->x + pos->y * pos->y + pos->z * pos->z);
      if (r > unit_leng) unit_leng = r;
   }
   double unit_time = sqrt(pow(unit_leng,3.0)/(Ggrav * unit_mass));
   double unit_velc = unit_leng / unit_time;
   double unit_eng  = unit_mass * unit_velc * unit_velc;
   printf("unit_mass = %15.7e [g]    = %15.7e [Msolar]\n",unit_mass,unit_mass/Msolar);
   printf("unit_leng = %15.7e [cm]   = %15.7e [kpc]\n",unit_leng,unit_leng/kpc);
   printf("unit_time = %15.7e [s]    = %15.7e [Gyr]\n",unit_time,unit_time/Gyr);
   printf("unit_velc = %15.7e [cm/s] = %15.7e [km/s]\n",unit_velc,unit_velc/km);
   for (i = 0; i < nptcl_nbody; i++) {
      ptcl_nbody[i].mass /= unit_mass;
      ptcl_nbody[i].pos.x /= unit_leng;
      ptcl_nbody[i].pos.y /= unit_leng;
      ptcl_nbody[i].pos.z /= unit_leng;
      ptcl_nbody[i].vel.x /= unit_velc;
      ptcl_nbody[i].vel.y /= unit_velc;
      ptcl_nbody[i].vel.z /= unit_velc;
   }
   for (i = 0; i < nptcl_sph; i++) {
      ptcl_sph[i].mass /= unit_mass;
      ptcl_sph[i].pos.x /= unit_leng;
      ptcl_sph[i].pos.y /= unit_leng;
      ptcl_sph[i].pos.z /= unit_leng;
      ptcl_sph[i].vel.x /= unit_velc;
      ptcl_sph[i].vel.y /= unit_velc;
      ptcl_sph[i].vel.z /= unit_velc;
      ptcl_sph[i].smth /= unit_leng;
      ptcl_sph[i].eng  /= (unit_eng/unit_mass);
   }
   // Set boundary condition
   *bc = FDPS_BC_OPEN;
   // Set gravitational softening
   eps_grav = 1.0e-3;
   // Set I/O intervals
   *dt_dump = 0.01;
   *time_dump = *dt_dump;
   *time_end = 1.0;
   // Set maximum timestep
   dt_max = 1.0e-4;
   // Output
   printf("An initial condition for isolated galaxy simulation is made.\n");
   printf("N_nbody = %d, N_sph = %d\n",nptcl_nbody,nptcl_sph);

   // In the following, we compute the surface gas density and
   // the distribution function of particle's z coordinates.
   // Compute the maxium clyndrical radius of gas particles
   double rmax = 0.0;
   for (i = 0; i < nptcl_sph; i++) {
      fdps_f64vec *pos = &ptcl_sph[i].pos;
      double r = sqrt(pos->x * pos->x + pos->y * pos->y);
      if (r > rmax) rmax = r;
   }
   // Compute the surface density
   double dr = (safety * rmax)/nbin;
   double sigma[nbin];
   for (i = 0; i < nbin; i++) sigma[i]=0.0;
   for (i = 0; i < nptcl_sph; i++) {
      fdps_f64vec *pos = &ptcl_sph[i].pos;
      double r = sqrt(pos->x * pos->x + pos->y * pos->y);
      int indx = r/dr;
      assert(indx >= 0 && indx < nbin);
      sigma[indx] += ptcl_sph[i].mass;
   }
   for (i = 0; i < nbin; i++)
       sigma[i] /= (2.0 * pi * (2*i + 1) * dr*dr);
   // Output the surface density
   char file_num[64] = {'\0'};
   strcpy(file_name,"Sigma.txt");
   FILE *fp;
   if ((fp = fopen(file_name,"w")) == NULL) {
       fprintf(stderr,"Cannot open file %s\n",file_name);
       exit(EXIT_FAILURE);
   }
   for (i = 0; i < nbin; i++)
       fprintf(fp,"%15.7e %15.7e\n",(i + 0.5) * dr, sigma[i]);
   fclose(fp);
   // Compute the distribution function 
   double zmax = zt/unit_leng;
   double dz = (safety * 2.0 * zmax)/nbin;
   int dist_func[nbin];
   for (i = 0; i < nbin; i++) dist_func[i] = 0;
   for (i = 0; i < nptcl_sph; i++) {
      int indx = (ptcl_sph[i].pos.z + safety * zmax)/dz;
      assert(indx >= 0 && indx < nbin);
      dist_func[indx]++;
   }
   // Output the distribution function
   memset(file_num,'\0',strlen(file_name));
   strcpy(file_name,"dist_func.txt");
   if ((fp = fopen(file_name,"w")) == NULL) {
      fprintf(stderr,"Cannot open file %s\n",file_name);
      exit(EXIT_FAILURE);
   }
   for (i = 0; i < nbin; i++) {
       double z = (i + 0.5) * dz - safety * zmax;
       fprintf(fp,"%15.7e %d\n",z,dist_func[i]);
   }
   fclose(fp);
}

void cold_collapse_test_IC(int psys_num_nbody,
                           int psys_num_sph,
                           int *bc,
                           fdps_f32vec *pos_root_domain_low,
                           fdps_f32vec *pos_root_domain_high,
                           double *time_dump, 
                           double *dt_dump,
                           double *time_end) {
   // Local parameters
   const double M_nbody = 0.5;
   const double M_sph   = 0.5;
   const int N_nbody = 512;
   const int N_sph   = 512;
   const double R_nbody = 3.0;
   const double R_sph   = 3.0;
   const double E_bind_nbody = - 3.0 * M_nbody * M_nbody / 5.0 / R_nbody;
   const double E_bind_sph   = - 3.0 * M_sph * M_sph / 5.0 / R_sph;
   const double virial_ratio = 0.5;

   // Initialize pseudorandom number generator
   int mtts_num;
   fdps_create_mtts(&mtts_num);
   fdps_mtts_init_genrand(mtts_num,0);
   // Place Nbody particles
   fdps_set_nptcl_loc(psys_num_nbody,N_nbody);
   FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
   int i;
   for (i = 0; i < N_nbody; i++) {
       ptcl_nbody[i].id = i;
       ptcl_nbody[i].mass = M_nbody / N_nbody;
       fdps_f64vec pos;
       for (;;) {
           pos.x = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * R_nbody;
           pos.y = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * R_nbody;
           pos.z = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * R_nbody;
           double r2 = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
           if (r2 <= R_nbody * R_nbody) break;
       }
       ptcl_nbody[i].pos.x = pos.x;
       ptcl_nbody[i].pos.y = pos.y;
       ptcl_nbody[i].pos.z = pos.z;
       ptcl_nbody[i].vel.x = 0.0;
       ptcl_nbody[i].vel.y = 0.0;
       ptcl_nbody[i].vel.z = 0.0;
       ptcl_nbody[i].acc.x = 0.0;
       ptcl_nbody[i].acc.y = 0.0;
       ptcl_nbody[i].acc.z = 0.0;
       ptcl_nbody[i].pot   = 0.0;
   }
   // Place SPH particles
   fdps_set_nptcl_loc(psys_num_sph,N_sph);
   FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
   for (i = 0; i < N_sph; i++) {
       ptcl_sph[i].id = i + N_nbody;
       ptcl_sph[i].mass = M_sph / N_sph;
       fdps_f64vec pos;
       for (;;) {
           pos.x = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * R_sph;
           pos.y = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * R_sph;
           pos.z = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0) * R_sph;
           double r2 = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z;
           if (r2 <= R_sph * R_sph) break;
       }
       ptcl_sph[i].pos.x = pos.x;
       ptcl_sph[i].pos.y = pos.y;
       ptcl_sph[i].pos.z = pos.z;
       ptcl_sph[i].vel.x = 0.0;
       ptcl_sph[i].vel.y = 0.0;
       ptcl_sph[i].vel.z = 0.0;
       ptcl_sph[i].acc_grav.x = 0.0;
       ptcl_sph[i].acc_grav.y = 0.0;
       ptcl_sph[i].acc_grav.z = 0.0;
       ptcl_sph[i].pot_grav = 0.0;
       ptcl_sph[i].acc_hydro.x = 0.0;
       ptcl_sph[i].acc_hydro.y = 0.0;
       ptcl_sph[i].acc_hydro.z = 0.0;
       ptcl_sph[i].eng = virial_ratio * fabs(E_bind_sph) / M_sph;  // specific thermal energy
       double dens = M_sph / (4.0 * pi * R_sph * R_sph * R_sph / 3.0); 
       double p = specific_heat_ratio - 1.0;
       ptcl_sph[i].ent = p * ptcl_sph[i].eng / pow(dens,p);
       ptcl_sph[i].smth = R_sph * pow(((double)N_neighbor)/((double)N_sph),1.0/3.0); // smoothing length
   }
   // Set boundary condition
   *bc = FDPS_BC_OPEN;
   // Set other parameters
   eps_grav = 0.01 * R_nbody;
   // Set I/O intervals
   double dens_nbody = 3.0 * M_nbody / (4.0 * pi * R_nbody * R_nbody * R_nbody);
   double dens_sph   = 3.0 * M_sph / (4.0 * pi * R_sph * R_sph * R_sph);
   double unit_dens = dens_nbody + dens_sph;
   double unit_time = sqrt(3.0 * pi/ (32.0 * unit_dens));
   printf("unit_dens = %15.7e\n",unit_dens);
   printf("unit_time = %15.7e\n",unit_time);
   *dt_dump   = 0.1 * unit_time;
   *time_dump = *dt_dump;
   *time_end  = 2.0 * unit_time;
   // Set the maximum timestep
   dt_max = 1.0e-3;

}

void Evrard_test_IC(int psys_num_nbody,
                    int psys_num_sph,
                    int *bc,
                    fdps_f32vec *pos_root_domain_low,
                    fdps_f32vec *pos_root_domain_high,
                    double *time_dump,
                    double *dt_dump,
                    double *time_end,
                    int gen_mode) {
   // Local parameters
   const int N_nbody=1;
   const int N_ptcl_per_side=64;
   const double radius_of_sphere = 1.0;
   const double M_sph = 1.0;
   int N_sph;

   // Place Nbody particles
   // here, we place only one dummy, mass-less particle.
   fdps_set_nptcl_loc(psys_num_nbody,N_nbody);
   FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
   int i,j,k;
   for (i = 0; i < N_nbody; i++) {
       ptcl_nbody[i].id = i;
       ptcl_nbody[i].mass = 0.0;
       ptcl_nbody[i].pos.x = 0.0;
       ptcl_nbody[i].pos.y = 0.0;
       ptcl_nbody[i].pos.z = 0.0;
       ptcl_nbody[i].vel.x = 0.0;
       ptcl_nbody[i].vel.y = 0.0;
       ptcl_nbody[i].vel.z = 0.0;
       ptcl_nbody[i].acc.x = 0.0;
       ptcl_nbody[i].acc.y = 0.0;
       ptcl_nbody[i].acc.z = 0.0;
       ptcl_nbody[i].pot = 0.0;
   }
   // Place SPH particles
   if (gen_mode == 0) {
      // In this mode, we create an initial distribution of particles
      // by rescaling the positions of particles which are placed in a grid.
      // (1) Count # of particles in a sphere of radius 1
      N_sph = 0;
      printf("OK1\n");
      // x-loop
      for (i = 0; i < N_ptcl_per_side; i++) {
          double x = get_pos_cell_center(-radius_of_sphere, radius_of_sphere,
                                         N_ptcl_per_side, i);
          // y-loop
          for (j = 0; j < N_ptcl_per_side; j++) {
              double y = get_pos_cell_center(-radius_of_sphere, radius_of_sphere,
                                             N_ptcl_per_side, j);
              // z-loop
              for (k = 0; k < N_ptcl_per_side; k++) {
                  double z = get_pos_cell_center(-radius_of_sphere, radius_of_sphere,
                                                 N_ptcl_per_side, k);
                  double r = sqrt(x*x + y*y + z*z);
                  if (r <= radius_of_sphere) N_sph++;
              }
          }
      }
      printf("N_sph = %d\n",N_sph);
      printf("OK2\n");
      if (N_sph == 0) {
         printf("There are no SPH particles in the specified sphere!\n");
         fdps_abort(-1);
         exit(EXIT_FAILURE);
      }
      fdps_set_nptcl_loc(psys_num_sph,N_sph);
      FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
      // (2) Actually place particles
      int id = 0;
      // x-loop
      for (i = 0; i < N_ptcl_per_side; i++) {
          double x = get_pos_cell_center(-radius_of_sphere, radius_of_sphere,
                                         N_ptcl_per_side, i);
          // y-loop
          for (j = 0; j < N_ptcl_per_side; j++) {
              double y = get_pos_cell_center(-radius_of_sphere, radius_of_sphere,
                                             N_ptcl_per_side, j);
              // z-loop
              for (k = 0; k < N_ptcl_per_side; k++) {
                  double z = get_pos_cell_center(-radius_of_sphere, radius_of_sphere,
                                                 N_ptcl_per_side, k);
                  double r = sqrt(x*x + y*y + z*z);
                  if (r <= radius_of_sphere) {
                      double r_new = radius_of_sphere * pow(r/radius_of_sphere,1.5);
                      ptcl_sph[id].id    = id + N_nbody;
                      ptcl_sph[id].mass  = M_sph / N_sph;
                      ptcl_sph[id].pos.x = (r_new / r) * x;
                      ptcl_sph[id].pos.y = (r_new / r) * y;
                      ptcl_sph[id].pos.z = (r_new / r) * z;
                      ptcl_sph[id].vel.x = 0.0;
                      ptcl_sph[id].vel.y = 0.0;
                      ptcl_sph[id].vel.z = 0.0;
                      ptcl_sph[id].acc_grav.x = 0.0;
                      ptcl_sph[id].acc_grav.y = 0.0;
                      ptcl_sph[id].acc_grav.z = 0.0;
                      ptcl_sph[id].pot_grav = 0.0;
                      ptcl_sph[id].acc_hydro.x = 0.0;
                      ptcl_sph[id].acc_hydro.y = 0.0;
                      ptcl_sph[id].acc_hydro.z = 0.0;
                      ptcl_sph[id].eng = 0.05 * M_sph / radius_of_sphere;
                      double dens = M_sph / (2.0 * pi * radius_of_sphere * radius_of_sphere * r_new);
                      double p = specific_heat_ratio - 1.0;
                      ptcl_sph[id].ent = p * ptcl_sph[id].eng / pow(dens,p);
                      ptcl_sph[id].smth = radius_of_sphere * pow(((double)N_neighbor)/((double)N_sph),1.0/3.0);
                      id++;
                  }
              }
          }
      }
   } else {
      // In this mode, we set an initial distribution of particles by reading a file.
      // (1) Read particle data
      char file_name[64] = {'\0'};
      strcpy(file_name,"result/glass_data_header.dat");
      FILE *fp;
      if ((fp = fopen(file_name,"rb")) == NULL) {
          fprintf(stderr,"Cannot open file %s.",file_name);
          exit(EXIT_FAILURE);
      }
      int nptcl_in;
      fread(&nptcl_in,sizeof(int),1,fp);
      fclose(fp);
      fdps_f64vec *pos = (fdps_f64vec *)malloc(sizeof(fdps_f64vec) * nptcl_in);
      for (i = 0; i < strlen(file_name); i++) file_name[i] = '\0'; // clear
      strcpy(file_name,"result/glass_data_body.dat");
      if ((fp = fopen(file_name,"rb")) == NULL) {
         fprintf(stderr,"Cannot open file %s.\n",file_name);
         exit(EXIT_FAILURE);
      }
      for (i = 0; i < nptcl_in; i++)
         fread(pos+i,sizeof(fdps_f64vec),1,fp);
      fclose(fp);
      printf("%d particles are read.\n",nptcl_in);
#if 0
      // For debug
      memset(file_name,'\0',strlen(file_name));
      strcpy(file_name,"glass_data.txt");
      if ((fp = fopen(file_name,"w")) == NULL) {
         fprintf(stderr,"Cannot open file %s.\n",file_name);
         exit(EXIT_FAILURE);
      }
      for (i = 0; i< nptcl_in; i++)
          fprintf(fp,"%15.7e %15.7e %15.7e\n",
                  pos[i].x,pos[i].y,pos[i].z);
      fclose(fp);
#endif
      // (2) Count # of particles in a sphere of radius 1
      N_sph = 0;
      for (i = 0; i < nptcl_in; i++) {
         double r2 = pos[i].x * pos[i].x 
                   + pos[i].y * pos[i].y 
                   + pos[i].z * pos[i].z;
         if (r2 < 1.0) N_sph++;
      }
      printf("%d particles of them will be used to make a Evrard sphere.\n",N_sph);
      // (3) Place SPH particles
      fdps_set_nptcl_loc(psys_num_sph,N_sph);
      FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
      int j = -1;
      for (i = 0; i < N_sph; i++) {
         ptcl_sph[i].id = i;
         ptcl_sph[i].mass = M_sph / N_sph;
         double r2;
         for (;;) {
            j++;
            r2 = pos[j].x * pos[j].x 
               + pos[j].y * pos[j].y 
               + pos[j].z * pos[j].z;
            if (r2 < 1.0) break;
         }
         double r = sqrt(r2);
         if (r > 0.0) {
             double r_new = radius_of_sphere * pow(r,1.5);
             ptcl_sph[i].pos.x = (r_new / r) * pos[j].x;
             ptcl_sph[i].pos.y = (r_new / r) * pos[j].y;
             ptcl_sph[i].pos.z = (r_new / r) * pos[j].z;
         } else {
             ptcl_sph[i].pos.x = pos[j].x;
             ptcl_sph[i].pos.y = pos[j].y;
             ptcl_sph[i].pos.z = pos[j].z;
         }
         ptcl_sph[i].vel.x = 0.0;
         ptcl_sph[i].vel.y = 0.0;
         ptcl_sph[i].vel.z = 0.0;
         ptcl_sph[i].acc_grav.x = 0.0;
         ptcl_sph[i].acc_grav.y = 0.0;
         ptcl_sph[i].acc_grav.z = 0.0;
         ptcl_sph[i].pot_grav = 0.0;
         ptcl_sph[i].acc_hydro.x = 0.0;
         ptcl_sph[i].acc_hydro.y = 0.0;
         ptcl_sph[i].acc_hydro.z = 0.0;
         ptcl_sph[i].eng = 0.05 * M_sph / radius_of_sphere;
         ptcl_sph[i].smth = radius_of_sphere * pow(((double)N_neighbor)/((double)N_sph),1.0/3.0);
         // Note that entropy is determined in setEntropy().
      }
   }
   // Set boundary condition
   *bc = FDPS_BC_OPEN;
   // Set the other parameters
   eps_grav = 1.0e-4 * radius_of_sphere;
   printf("An initial condition for the Evrard test is made.\n");
   printf("(N_sph = %d\n)",N_sph);
   // Set I/O intervals
   double unit_dens = 3.0 * M_sph / (4.0 * pi * pow(radius_of_sphere,3.0));
   double unit_time = sqrt(pi * pi / 8.0) * pow(radius_of_sphere,1.5) / sqrt(M_sph);
   printf("unit_dens = %15.7e\n",unit_dens);
   printf("unit_time = %15.7e\n",unit_time);
   *dt_dump   = 0.1 * unit_time;
   *time_dump = *dt_dump;
   *time_end  = unit_time;
   // Set the maximum timestep
   dt_max = 1.0e-3;

}

double get_pos_cell_center(double pos_left_bound,
                           double pos_right_bound,
                           int number_of_cells,
                           int i) {
    // Error check
    if ((i < 0) || (number_of_cells <= i) ||
        (pos_left_bound > pos_right_bound)) {
       fprintf(stderr,"Given arguments are out of ranges!\n");
       fprintf(stderr,"pos_left_bound  = %15.7e\n",pos_left_bound);
       fprintf(stderr,"pos_right_bound = %15.7e\n",pos_right_bound);
       fprintf(stderr,"number_of_cells = %d\n",number_of_cells);
       fprintf(stderr,"i               = %d\n",i);
       fflush(stderr);
       exit(EXIT_FAILURE);
    }

    double dx = (pos_right_bound - pos_left_bound) / number_of_cells;
    if ( number_of_cells % 2 == 0) {
        // # of cells is even.
        if (i < number_of_cells/2) {
            return pos_left_bound + dx * (i + 0.5);
        } else {
            return pos_right_bound - dx * (number_of_cells - i - 0.5);
        }
    } else {
        // # of cells is odd.
        double center = 0.5 * (pos_left_bound + pos_right_bound);
        if (i < number_of_cells/2) {
            return center - dx * (number_of_cells/2 - i);
        } else {
            return center + dx * (i - number_of_cells/2);
        }
    }
}

void make_glass_IC(int psys_num_nbody,
                   int psys_num_sph,
                   int *bc,
                   fdps_f32vec *pos_root_domain_low,
                   fdps_f32vec *pos_root_domain_high,
                   double *time_dump,
                   double *dt_dump,
                   double *time_end) {
   // Local parameters
   const int N_nbody = 1; // dummy value
   const int N_sph   = pow(2,18);

   // Initialize pseudorandom number generator
   int mtts_num;
   fdps_create_mtts(&mtts_num);
   fdps_mtts_init_genrand(mtts_num,0);
   // Place Nbody particles
   fdps_set_nptcl_loc(psys_num_nbody,N_nbody);
   FP_nbody *ptcl_nbody = (FP_nbody *) fdps_get_psys_cptr(psys_num_nbody);
   int i;
   for (i = 0; i < N_nbody; i++) {
       ptcl_nbody[i].id   = i;
       ptcl_nbody[i].mass = 0.0;
       ptcl_nbody[i].pos.x = 0.0;
       ptcl_nbody[i].pos.y = 0.0;
       ptcl_nbody[i].pos.z = 0.0;
       ptcl_nbody[i].vel.x = 0.0;
       ptcl_nbody[i].vel.y = 0.0;
       ptcl_nbody[i].vel.z = 0.0;
       ptcl_nbody[i].acc.x = 0.0;
       ptcl_nbody[i].acc.y = 0.0;
       ptcl_nbody[i].acc.z = 0.0;
       ptcl_nbody[i].pot = 0.0;
   }
   // Place SPH particles
   fdps_set_nptcl_loc(psys_num_sph,N_sph);
   FP_sph *ptcl_sph = (FP_sph *) fdps_get_psys_cptr(psys_num_sph);
   for (i = 0; i < N_sph; i++) {
       ptcl_sph[i].id = i + N_nbody;
       ptcl_sph[i].mass = 8.0 / N_sph;
       double x,y,z;
       for (;;) {
           x = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0); 
           y = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0);
           z = (2.0 * fdps_mtts_genrand_res53(mtts_num) - 1.0);
           if ((-1.0 <= x) && (x < 1.0) &&
               (-1.0 <= y) && (y < 1.0) &&
               (-1.0 <= z) && (z < 1.0)) break;
       }
       ptcl_sph[i].pos.x = x;
       ptcl_sph[i].pos.y = y;
       ptcl_sph[i].pos.z = z;
       ptcl_sph[i].vel.x = 0.0;
       ptcl_sph[i].vel.y = 0.0;
       ptcl_sph[i].vel.z = 0.0;
       ptcl_sph[i].acc_grav.x = 0.0;
       ptcl_sph[i].acc_grav.y = 0.0;
       ptcl_sph[i].acc_grav.z = 0.0;
       ptcl_sph[i].pot_grav = 0.0;
       ptcl_sph[i].acc_hydro.x = 0.0;
       ptcl_sph[i].acc_hydro.y = 0.0;
       ptcl_sph[i].acc_hydro.z = 0.0;
       ptcl_sph[i].eng = 1.0;
       ptcl_sph[i].smth = 2.0 * (((double)N_neighbor)/((double)N_sph),1.0/3.0);
       // [Notes]
       //   (1) The value of the specific thermal energy is chosen 
       //       so that the sound speed is nearly equal to 1.
       //   (2) The value of the entropy is determined in set_entropy().
   }
   // Set boundary condition
   *bc = FDPS_BC_PERIODIC_XYZ;
   pos_root_domain_low->x  = - 1.0;
   pos_root_domain_low->y  = - 1.0;
   pos_root_domain_low->z  = - 1.0;
   pos_root_domain_high->x = 1.0;
   pos_root_domain_high->y = 1.0;
   pos_root_domain_high->z = 1.0;
   // Set gravitational softening
   eps_grav = 0.01;
   // Set I/O intervals
   double tcross = 2.0 * sqrt(3.0);
   printf("The sound crossing time = %15.7e\n",tcross);
   *dt_dump   = 4.0 * tcross;
   *time_dump = *dt_dump;
   *time_end  = 64.0 * tcross;
   // Set maximum timestep
   dt_max = 1.0e-2;

}
