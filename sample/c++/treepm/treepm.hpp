#ifndef __TREEPM_HPP__
#define __TREEPM_HPP__

#include <sys/times.h>
#include <sys/time.h>
#include <param_fdps.h>
#include <mpi.h>

#include "run_param.hpp"
#include "cosmology.hpp"

#ifdef ENABLE_PHANTOM_GRAPE_X86
#include "gp5util.h"
#include "pg5_table.h"
#endif

#define TINY (1.0e-30)

class Result_treepm {
public:
  PS::F32vec acc;
  PS::F32    pot;

  void clear() {
    acc = 0.0;
    pot = 0.0;
  }
};


class FPtreepm {
private:
    template<class T>
    T reverseEndian(T value){
        char * first = reinterpret_cast<char*>(&value);
        char * last = first + sizeof(T);
        std::reverse(first, last);
        return value;
    }
public:
    PS::S64    id;
    PS::F32    mass;
    PS::F32    eps;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64vec acc_pm;
    
    //static PS::F64vec low_boundary;
    //static PS::F64vec high_boundary;
    //static PS::F64 unit_l;
    //static PS::F64 unit_m;
    static PS::F64 H0;
    static PS::F64 Lbnd;
    
    PS::F64vec getPos() const {
        return pos;
    }

    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromForce(const Result_treepm & force) {
        this->acc = force.acc;
    }

    PS::F64 getRSearch() const {
        PS::F64 rcut = 3.0/SIZE_OF_MESH;
        return rcut;
    }

    void setPos(const PS::F64vec pos_new) {
        this->pos = pos_new;
    }

    PS::F64 getChargeParticleMesh() const {
        return this->mass;
    }

    void copyFromForceParticleMesh(const PS::F64vec & acc_pm) {
        this->acc_pm = acc_pm;
    }

    /*    
    void writeParticleBinary(FILE *fp) {
        int count;
        count = 0;

        count += fwrite(&mass,   sizeof(PS::F32),1,fp);
        count += fwrite(&eps,    sizeof(PS::F32),1,fp);
        count += fwrite(&pos[0], sizeof(PS::F64),1,fp);
        count += fwrite(&pos[1], sizeof(PS::F64),1,fp);
        count += fwrite(&pos[2], sizeof(PS::F64),1,fp);
        count += fwrite(&vel[0], sizeof(PS::F64),1,fp);
        count += fwrite(&vel[1], sizeof(PS::F64),1,fp);
        count += fwrite(&vel[2], sizeof(PS::F64),1,fp);
    }
    */
    /*
    int readParticleBinary(FILE *fp) {
        int count;
        count = 0;
        
        count += fread(&mass,   sizeof(PS::F32),1,fp);
        count += fread(&eps,    sizeof(PS::F32),1,fp);
        count += fread(&pos[0], sizeof(PS::F64),1,fp);
        count += fread(&pos[1], sizeof(PS::F64),1,fp);
        count += fread(&pos[2], sizeof(PS::F64),1,fp);
        count += fread(&vel[0], sizeof(PS::F64),1,fp);
        count += fread(&vel[1], sizeof(PS::F64),1,fp);
        count += fread(&vel[2], sizeof(PS::F64),1,fp);
        
        return count;
    }
    */

    void writeParticleBinary(FILE *fp) {
        PS::F32 x = pos[0];
        PS::F32 y = pos[1];
        PS::F32 z = pos[2];
        PS::F32 vx = vel[0];
        PS::F32 vy = vel[1];
        PS::F32 vz = vel[2];    
        PS::S32 i = id;
        PS::S32 m = mass;
        fwrite(&x, sizeof(PS::F32),1,fp);
        fwrite(&vx, sizeof(PS::F32),1,fp);
        fwrite(&y, sizeof(PS::F32),1,fp);
        fwrite(&vy, sizeof(PS::F32),1,fp);
        fwrite(&z, sizeof(PS::F32),1,fp);
        fwrite(&vz, sizeof(PS::F32),1,fp);    
        //fwrite(&mass,   sizeof(PS::F32),1,fp);
	fwrite(&m,   sizeof(PS::F32),1,fp);
	fwrite(&i,   sizeof(PS::F32),1,fp);
        //fwrite(&id,   sizeof(PS::F32),1,fp);
    }


    // for API of FDPS 
    // in snapshot, L unit is Mpc/h, M unit is Msun, v unit is km/s
    void readBinary(FILE *fp){
        static PS::S32 ONE = 1;
        static bool is_little_endian = *reinterpret_cast<char*>(&ONE) == ONE;
        static const PS::F64 Mpc_m = 3.08567e22; // unit is m
        static const PS::F64 Mpc_km = 3.08567e19; // unit is km    
        static const PS::F64 Msun_kg = 1.9884e30; // unit is kg
        static const PS::F64 G = 6.67428e-11; // m^3*kg^-1*s^-2
        static const PS::F64 Cl = 1.0 / FPtreepm::Lbnd;
        static const PS::F64 Cv = 1.0 / (FPtreepm::Lbnd * FPtreepm::H0);
        static const PS::F64 Cm = 1.0 / (pow(Mpc_m*FPtreepm::Lbnd, 3.0) / pow(Mpc_km/FPtreepm::H0, 2.0) / G / Msun_kg);
        PS::F32 x, y, z, vx, vy, vz, m;
        PS::S32 i;
        fread(&x, 4, 1, fp);
        fread(&vx, 4, 1, fp);
        fread(&y, 4, 1, fp);
        fread(&vy, 4, 1, fp);
        fread(&z, 4, 1, fp);
        fread(&vz, 4, 1, fp);
        fread(&m,   4, 1, fp);
        fread(&i,   4, 1, fp);
        if( is_little_endian){
            pos.x = x * Cl;
            pos.y = y * Cl;
            pos.z = z * Cl;
            vel.x = vx * Cv;
            vel.y = vy * Cv;
            vel.z = vz * Cv;
            mass = m * Cm;
            //mass = m / 1.524e17;
            id = i;
        }
        else{
            pos.x = reverseEndian(x) * Cl;
            pos.y = reverseEndian(y) * Cl;
            pos.z = reverseEndian(z) * Cl;
            vel.x = reverseEndian(vx) * Cv;
            vel.y = reverseEndian(vy) * Cv;
            vel.z = reverseEndian(vz) * Cv;
            mass = reverseEndian(m) * Cm;
            //mass = reverseEndian(m) / 1.524e17;
            id = reverseEndian(i);
        }
    }

    // for API of FDPS 
    void writeBinary(FILE *fp){
        static const PS::F64 Mpc_m = 3.08567e22; // unit is m
        static const PS::F64 Mpc_km = 3.08567e19; // unit is km    
        static const PS::F64 Msun_kg = 1.9884e30; // unit is kg
        static const PS::F64 G = 6.67428e-11; // m^3*kg^-1*s^-2
        static const PS::F64 Cl = FPtreepm::Lbnd;
        static const PS::F64 Cv = (FPtreepm::Lbnd * FPtreepm::H0);
        static const PS::F64 Cm = (pow(Mpc_m*FPtreepm::Lbnd, 3.0) / pow(Mpc_km/FPtreepm::H0, 2.0) / G / Msun_kg);    
        PS::F32vec x = pos * Cl;
        PS::F32vec v = vel * Cv;
        PS::F32 m = mass * Cm;
        PS::S32 i = id;
        fwrite(&x.x,   sizeof(PS::F32), 1, fp);
        fwrite(&v.x,   sizeof(PS::F32), 1, fp);
        fwrite(&x.y,   sizeof(PS::F32), 1, fp);
        fwrite(&v.y,   sizeof(PS::F32), 1, fp);
        fwrite(&x.z,   sizeof(PS::F32), 1, fp);
        fwrite(&v.z,   sizeof(PS::F32), 1, fp);
        fwrite(&m,   sizeof(PS::F32), 1, fp);
        fwrite(&i,   sizeof(PS::S32), 1, fp);
    }

    PS::F64 calcDtime(run_param &this_run) {
    PS::F64 dtime_v, dtime_a, dtime;
    PS::F64 vnorm, anorm;
    vnorm = sqrt(SQR(this->vel))+TINY;
    anorm = sqrt(SQR(this->acc+this->acc_pm))+TINY;

    dtime_v = this->eps/vnorm;
    dtime_a = sqrt(this->eps/anorm)*CUBE(this_run.anow);

    dtime = fmin(0.5*dtime_v, dtime_a);

    return dtime;
  }
};

//PS::F64vec FPtreepm::low_boundary;
//PS::F64vec FPtreepm::high_boundary;
//PS::F64 FPtreepm::unit_l;
//PS::F64 FPtreepm::unit_m;
PS::F64 FPtreepm::H0;
PS::F64 FPtreepm::Lbnd;

class EPItreepm {
public:
  PS::S64    id;
  PS::F32    eps;
  PS::F64vec pos;

  PS::F64vec getPos() const {
    return this->pos;
  }

  void copyFromFP(const FPtreepm & fp) {
    this->id = fp.id;
    this->eps = fp.eps;
    this->pos = fp.pos;
  }
  
};

class EPJtreepm {
public:
  PS::S64    id;
  PS::F64vec pos;
  PS::F64    mass;
  //  PS::F64    rcut;

  PS::F64vec getPos() const {
    return this->pos;
  }

  PS::F64 getCharge() const {
    return this->mass;
  }

  void copyFromFP(const FPtreepm & fp) {
    this->id = fp.id;
    this->mass = fp.mass;
    this->pos = fp.pos;
  }

  PS::F64 getRSearch() const {
    PS::F64 rcut = 3.0/SIZE_OF_MESH;
    return rcut;
  }

  void setPos(const PS::F64vec pos_new) {
    this->pos = pos_new;
  }
};

inline PS::F64 gfactor_S2(const PS::F64 rad, const PS::F64 eps_pm) 
{
  PS::F64 R;
  PS::F64 g;
  PS::F64 S;

  R = 2.0*rad/eps_pm;
  R = (R > 2.0) ? 2.0 : R;
  S = R-1.0;
  S = (S > 0.0) ? S : 0.0;

  g = 1.0 + CUBE(R)*(-1.6+SQR(R)*(1.6+R*(-0.5+R*(0.15*R-12.0/35.0))))
    -CUBE(S)*CUBE(S)*(3.0/35.0+R*(18.0/35.0+0.2*R));

  return g;
}

#ifdef ENABLE_PHANTOM_GRAPE_X86
template <class TPJ>
class calc_pp_force {
public:
    void operator()(EPItreepm *iptcl,
                    const PS::S32 ni,
                    TPJ *jptcl,
                    const PS::S32 nj,
                    Result_treepm *ppforce){
        static const PS::S32 IPTCL_MAX = 4096;
        static const PS::S32 JPTCL_MAX = 16384;
        static __thread PS::F64vec xi[IPTCL_MAX];
        static __thread PS::F64vec ai[IPTCL_MAX];
        static __thread PS::F64 pi[IPTCL_MAX];
        static __thread PS::F64vec xj[JPTCL_MAX];
        static __thread PS::F64 mj[JPTCL_MAX];
        const int j_loop_max = ((nj-1) / JPTCL_MAX) + 1;
        assert(ni <= IPTCL_MAX);
        for (int i=0; i<ni; i++){
            xi[i] = iptcl[i].getPos();
            ai[i] = 0.0;
            pi[i] = 0.0;
        }
        for (int j_loop=0; j_loop<j_loop_max; j_loop++){
            const PS::F32 eps = iptcl[0].eps;
            const int j_head = j_loop*JPTCL_MAX;
            const int nj_tmp = std::min( (nj-j_head), JPTCL_MAX );
            const int j_tail = j_head + nj_tmp;
            int j_cnt = 0;
            for (int j=j_head; j<j_tail; j++){
                mj[j_cnt] = jptcl[j].getCharge();
                xj[j_cnt] = jptcl[j].getPos();
                j_cnt++;
            }

            assert(j_cnt == nj_tmp);
            g5_set_xmj(0, nj_tmp, (double (*)[3])xj, mj);
            g5_set_eps_to_all(eps);
            g5_set_n(nj);
            g5_calculate_force_on_x((double(*)[3])xi, (double(*)[3])ai, pi, ni);
        }
        for (PS::S32 i = 0; i < ni; i++) {
            ppforce[i].acc = ai[i];
        }
    }
};
#else
template <class TPJ>
class calc_pp_force {
public:
    void operator () (EPItreepm *iptcl,
                      const PS::S32 ni,
                      TPJ *jptcl,
                      const PS::S32 nj,
                      Result_treepm *ppforce) {
        for (PS::S32 i=0;i < ni;i++) {
            PS::F64 eps2 = SQR(iptcl[i].eps);
            for (PS::S32 j=0; j < nj;j++) {
                PS::F64vec dr = iptcl[i].pos - jptcl[j].pos;
                PS::F64 rsq = dr*dr;
                PS::F64 rad = sqrt(rsq+eps2);
                PS::F64 gfact = gfactor_S2(rad, 3.0/SIZE_OF_MESH);
                PS::F64 rinv  = 1.0/rad;
                PS::F64 mrinv3 = jptcl[j].mass*CUBE(rinv);
                ppforce[i].acc -= dr*gfact*mrinv3;
            }
        }
    }
};
#endif

class Map2D{
private:
    PS::F64 * dens_loc;
    PS::F64 * dens_glb;
    PS::S32 n_mesh;
    PS::F64 l_mesh;
    PS::S32 snp_id;    
    std::string file_name_prefix;
    std::string dir_name_prefix;
public:
    template<class Tptcl>    
    void initialize(const Tptcl & ptcl,
            const run_param & this_run,
            const PS::F64 res = 1.0){
	//PS::S32 n_loc = ptcl.getNumberOfParticleLocal();
    PS::S32 n_tot = ptcl.getNumberOfParticleGlobal();
    n_mesh = cbrt((PS::F64)n_tot*1.0000001) * res;
    l_mesh = 1.0 / n_mesh;
    if(PS::Comm::getRank() == 0)std::cerr<<"n_mesh="<<n_mesh<<" l_mesh="<<l_mesh<<std::endl;
    dens_loc = new PS::F64[n_mesh*n_mesh];
    dens_glb = new PS::F64[n_mesh*n_mesh];
    snp_id = 0;
    dir_name_prefix = this_run.model_name;
    dir_name_prefix += "_";
    file_name_prefix = "map_2d";
    file_name_prefix += "_";    
    if(PS::Comm::getRank() == 0)std::cerr<<"dir_name_prefix="<<dir_name_prefix<<" file_name_prefix="<<file_name_prefix<<std::endl;
    }

    template<class Tptcl>    
    void write(const Tptcl & ptcl,
           const run_param & this_run){
    for(PS::S32 i=0; i<n_mesh*n_mesh; i++){
        dens_loc[i] = dens_glb[i] = 0.0;        
    }
    PS::S32 n_loc = ptcl.getNumberOfParticleLocal();
    PS::F64 m_tot_loc = 0.0;    
    for(PS::S32 i=0; i<n_loc; i++){    
        PS::S32 idx = (PS::S32)(ptcl[i].pos.x / l_mesh);
        PS::S32 idy = (PS::S32)(ptcl[i].pos.y / l_mesh);
        dens_loc[idx*n_mesh+idy] += ptcl[i].mass;
        m_tot_loc += ptcl[i].mass;
    }
    MPI_Allreduce((PS::F64*)&dens_loc[0], (PS::F64*)&dens_glb[0], n_mesh*n_mesh, PS::GetDataType<PS::F64>(), MPI_SUM, MPI_COMM_WORLD);
    PS::F64 m_tot_glb = 0.0;
    MPI_Allreduce(&m_tot_loc, &m_tot_glb, 1, PS::GetDataType<PS::F64>(), MPI_SUM, MPI_COMM_WORLD);
    for(PS::S32 i=0; i<n_mesh*n_mesh; i++){
        dens_glb[i] /= (m_tot_glb*l_mesh*l_mesh);
    }
    if(PS::Comm::getRank() == 0){
        std::ostringstream snp_id_str;
        //snp_id_str.setf(std::ios::right);
        //snp_id_str.width(5);
        //snp_id_str.fill('0');
        snp_id_str<<snp_id;
        std::ostringstream input;
        input<<dir_name_prefix<<snp_id_str.str()<<"/"<<file_name_prefix<<snp_id_str.str();
        std::cerr<<"input="<<input.str()<<std::endl;
        std::ofstream fout;
        fout.open(input.str());
        for(PS::S32 ix=0; ix<n_mesh; ix++){
        fout<<ix*l_mesh<<"  "<<"0.0"<<"  "<<dens_glb[ix*n_mesh]<<std::endl;
        for(PS::S32 iy=1; iy<n_mesh; iy++){
            fout<<ix*l_mesh<<"  "<<iy*l_mesh<<"  "<<dens_glb[ix*n_mesh+iy]<<std::endl;
        }
        fout<<std::endl;
        }
        fout.close();
    }
    snp_id++;
    }

};

void output_data(PS::ParticleSystem<FPtreepm> &ptcl, run_param &this_run, 
         char *filename)
{
  FILE *output_fp = fopen(filename,"w");

  if(output_fp == NULL) {
    fprintf(stderr, "File %s cannot be written.",filename);
    exit(EXIT_FAILURE);
  }

  this_run.mpi_rank = PS::Comm::getRank();
  this_run.mpi_nproc = PS::Comm::getNumberOfProc();
  this_run.npart_local = ptcl.getNumberOfParticleLocal();
  this_run.npart_total = ptcl.getNumberOfParticleGlobal();
  this_run.write_header(output_fp);

  for(PS::S64 i=0;i<this_run.npart_local;i++) {
    ptcl[i].writeParticleBinary(output_fp);
  }

  fclose(output_fp);
}

/*
void input_data(PS::ParticleSystem<FPtreepm> &ptcl, run_param &this_run, 
        char *filename)
{
  FILE *input_fp = fopen(filename,"r");

  if(input_fp == NULL) {
    fprintf(stderr, "File %s not found.\n", filename);
    exit(EXIT_FAILURE);
  }

  this_run.read_header(input_fp);

  ptcl.setNumberOfParticleLocal(this_run.npart_local);
  for(PS::S64 i=0;i<ptcl.getNumberOfParticleLocal();i++) {
    ptcl[i].readParticleBinary(input_fp);
  }

  fclose(input_fp);
}
*/

void output_data_in_run(PS::ParticleSystem<FPtreepm> &ptcl, run_param &this_run)
{
  static char filename[256], directory_name[256];

  sprintf(directory_name,"%s_%d",
      this_run.model_name, this_run.output_indx);
  sprintf(filename, "%s/%s_%d-%d", 
      directory_name, this_run.model_name, 
      this_run.output_indx, this_run.mpi_rank);

  make_directory(directory_name);
  
  if(this_run.znow < this_run.output_timing[this_run.output_indx]+0.001) {
    output_data(ptcl, this_run, filename);
    this_run.output_indx++;
  }
}

void output_data_in_run(PS::ParticleSystem<FPtreepm> &ptcl,
            run_param &this_run,
            Map2D & map)
{
  static char filename[256], directory_name[256];

  sprintf(directory_name,"%s_%d",
      this_run.model_name, this_run.output_indx);
  sprintf(filename, "%s/%s_%d-%d", 
      directory_name, this_run.model_name, 
      this_run.output_indx, this_run.mpi_rank);

  make_directory(directory_name);
  
  if(this_run.znow < this_run.output_timing[this_run.output_indx]+0.001) {
    output_data(ptcl, this_run, filename);
    map.write(ptcl, this_run);
    this_run.output_indx++;
  }
}

void drift_ptcl(PS::ParticleSystem<FPtreepm> &ptcl, PS::DomainInfo &dinfo, 
        const PS::F64 dtime)
{
  PS::S32 npart_local = ptcl.getNumberOfParticleLocal();
  for(PS::S64 i=0;i<npart_local;i++) 
    ptcl[i].pos += ptcl[i].vel*dtime;

  ptcl.adjustPositionIntoRootDomain(dinfo);
}

void reverse_ptcl_acc(PS::ParticleSystem<FPtreepm> &ptcl)
{
  PS::S64 npart_local = ptcl.getNumberOfParticleLocal();

  for(PS::S64 i=0;i<npart_local;i++) {
    ptcl[i].acc *= -1.0;
    ptcl[i].acc_pm *= -1.0;
  }

}

void kick_ptcl(PS::ParticleSystem<FPtreepm> &ptcl, const PS::F64 dtime, 
           run_param &this_run)
{
  PS::S64 npart_local = ptcl.getNumberOfParticleLocal();

  PS::F64 om = this_run.cosm.omegam;
  PS::F64 ov = this_run.cosm.omegav;

  PS::F64 anow = this_run.cosm.timetoa(this_run.tnow);
  
  PS::F64 at = sqrt(1.e0+om*(1.e0/anow-1.e0)+ov*(SQR(anow)-1.e0))/anow;
  PS::F64 bt = 1.0/CUBE(anow);

  PS::F64 atdt1 = 1.0+at*dtime;
  PS::F64 vfact = (2.0-atdt1)/atdt1;
  PS::F64 afact = bt*dtime/atdt1;

  for(PS::S64 i=0;i<npart_local;i++) {
    ptcl[i].vel = vfact*ptcl[i].vel + afact*(ptcl[i].acc+ptcl[i].acc_pm);
    //    ptcl[i].vel = vfact*ptcl[i].vel - afact*(ptcl[i].acc+ptcl[i].acc_pm);
  }
}


PS::F64 calc_dtime(PS::ParticleSystem<FPtreepm> &ptcl, run_param &this_run)
{
  PS::F64 dtime;

  dtime = DBL_MAX;
  for(PS::S64 i=0;i<ptcl.getNumberOfParticleLocal();i++) {
    dtime = fmin(dtime, ptcl[i].calcDtime(this_run));
  }

  if(this_run.noutput >= 1) {
    COSM::REAL zred_next;
    COSM::REAL zred_next_output;
    COSM::REAL time_next_output;

    if(this_run.znow < this_run.output_timing[this_run.output_indx] && 
               this_run.output_indx+1 < this_run.noutput) {
     zred_next_output = this_run.output_timing[this_run.output_indx+1];
    }else{
      zred_next_output = this_run.output_timing[this_run.output_indx];
    }
    zred_next = this_run.cosm.timetoz(this_run.tnow+dtime);

    if(zred_next < zred_next_output){
      time_next_output = this_run.cosm.ztotime(zred_next_output/1.0001);
      dtime = time_next_output - this_run.tnow;
    }
    
  }

  MPI_Allreduce(MPI_IN_PLACE, &dtime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  
  return dtime;
}


#endif /* __TREEPM_HPP__ */
