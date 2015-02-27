#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>

#ifdef USEPHANTOMGRAPE
#include"phantomgrape.hpp"
#endif

class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp){
		fscanf(fp, "%e\n", &time);
		fscanf(fp, "%d\n", &n_body);
		return n_body;
    }
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%e\n", time);
		fprintf(fp, "%d\n", n_body);
	}
};


class ForceGrav{
public:
    PS::F64vec acc;
    PS::F64 pot;
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
};

class FPGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const ForceGrav & force){
        acc = force.acc;
        pot = force.pot;
    }
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z);
	}

	void readAscii(FILE* fp){
		fscanf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z);
	}

};

class EPIGrav{
public:
    PS::S64 id;
    PS::F64vec pos;
    static PS::F64 eps;
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FPGrav & fp){ 
        pos = fp.pos;
        id = fp.id;
    }
};

PS::F64 EPIGrav::eps = 1.0/32.0;

class EPJGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    void copyFromFP(const FPGrav & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & pos_new){ pos = pos_new;}
    PS::F64 getCharge() const { return mass; }
};

#ifdef USEPHANTOMGRAPE

struct CalcForceEpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const EPJGrav * ep_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        double xi[PhantomGRAPE::PG_NIMAX][3];
        double mxj[PhantomGRAPE::PG_NJMAX][4];
        double ai[PhantomGRAPE::PG_NIMAX][3];
        double pi[PhantomGRAPE::PG_NIMAX];
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        const PS::F64 pot_corr = 1.0 / sqrt(eps2) * ep_j[0].getCharge(); // potential correction for equal mass
        for(PS::S32 i=0; i<n_ip; i++){
            xi[i][0] = ep_i[i].getPos().x;
            xi[i][1] = ep_i[i].getPos().y;
            xi[i][2] = ep_i[i].getPos().z;

        }
        for(PS::S32 j=0; j<n_jp; j++){
            mxj[j][0] = ep_j[j].getPos().x;
            mxj[j][1] = ep_j[j].getPos().y;
            mxj[j][2] = ep_j[j].getPos().z;
            mxj[j][3] = ep_j[j].getCharge();
        }
        PhantomGRAPE pg;
        pg.set_eps2(eps2);
        pg.set_xj(n_jp, mxj);
        pg.set_xi(n_ip, xi);
        pg.run(n_ip, n_jp);
        pg.get_ai(n_ip, ai, pi);

        for(PS::S32 i=0; i<n_ip; i++){
            force[i].acc.x += ai[i][0];
            force[i].acc.y += ai[i][1];
            force[i].acc.z += ai[i][2];
            force[i].pot += pi[i] + pot_corr;
        }
    }
};

#else
struct CalcForceEpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const EPJGrav * ep_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){

        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S64 idi = ep_i[i].id;
            for(PS::S32 j=0; j<n_jp; j++){
                if( idi == ep_j[j].id ) continue;
                PS::F64vec rij = xi - ep_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#endif // USEPHANTOMGRAPE

#ifdef MONOPOLE

#ifdef USEPHANTOMGRAPE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        double xi[PhantomGRAPE::PG_NIMAX][3];
        double mxj[PhantomGRAPE::PG_NJMAX][4];
        double ai[PhantomGRAPE::PG_NIMAX][3];
        double pi[PhantomGRAPE::PG_NIMAX];
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            xi[i][0] = ep_i[i].pos.x;
            xi[i][1] = ep_i[i].pos.y;
            xi[i][2] = ep_i[i].pos.z;
        }
        for(PS::S32 j=0; j<n_jp; j++){
            mxj[j][0] = sp_j[j].getPos().x;
            mxj[j][1] = sp_j[j].getPos().y;
            mxj[j][2] = sp_j[j].getPos().z;
            mxj[j][3] = sp_j[j].getCharge();
        }
        PhantomGRAPE pg;
        pg.set_eps2(eps2);
        pg.set_xj(n_jp, mxj);
        pg.set_xi(n_ip, xi);
        pg.run(n_ip, n_jp);
        pg.get_ai(n_ip, ai, pi);

        for(PS::S32 i=0; i<n_ip; i++){
            force[i].acc.x += ai[i][0];
            force[i].acc.y += ai[i][1];
            force[i].acc.z += ai[i][2];
            force[i].pot += pi[i];
        }

    }
};
#else
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#endif // #ifdef USEPHANTOMGRAPE

#elif QUADRUPOLE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJQuadrupole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F64vec xi = ep_i[ip].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64 mj = sp_j[jp].mass;
                PS::F64vec xj= sp_j[jp].pos;
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};
#elif MONOPOLEGEOMETRICCENTER
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopoleGeometricCenter * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].charge;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#elif DIPOLEGEOMETRICCENTER
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJDipoleGeometricCenter * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec di = sp_j[j].dipole;
                const PS::F64vec rij = xi - sp_j[j].pos;
                const PS::F64 r2 = rij * rij + eps2;
                const PS::F64 r_inv = 1.0/sqrt(r2);
                const PS::F64 r2_inv = r_inv * r_inv;
                const PS::F64 r3_inv = r_inv * r2_inv;
                const PS::F64vec hij = rij * r_inv;
                const PS::F64 dihij = di * hij;
                const PS::F64 mj = sp_j[j].charge;
                poti -= mj * r_inv + dihij* r2_inv;
                ai -= (mj*r2_inv + 3.0*dihij*r3_inv) * hij - di;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#elif QUADRUPOLEGEOMETRICCENTER
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJQuadrupoleGeometricCenter * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            const PS::F64vec xi = ep_i[ip].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64 mj = sp_j[jp].charge;
                PS::F64vec xj= sp_j[jp].pos;
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quadrupole;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
                const PS::F64vec di = sp_j[jp].dipole;
                const PS::F64 dirij = di * rij;
                poti -= dirij* r3_inv;
                ai -= 3.0*dirij*r5_inv*rij - di;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};
#endif

template<class Tpsys>
void ReadNemoAscii(Tpsys & psys,
                   PS::S32 & n_glb,
                   PS::S32 & n_loc,  
                   PS::F32 & t_sys,
                   const char * ifile){
    std::ifstream finput;
    finput.open(ifile);
    assert(finput);
    PS::S32 dim;
    finput>>n_glb>>dim>>t_sys;
    std::cerr<<"ifile:"<<ifile<<std::endl;
    std::cerr<<"n_glb="<<n_glb<<std::endl;
    std::cerr<<"dim="<<dim<<std::endl;
    std::cerr<<"t_sys="<<t_sys<<std::endl;

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb/n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;

    psys.setNumberOfParticleLocal(n_loc);

    PS::F32vec pos_shift = 0.0;

    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    const PS::S32 i_t = i_h+n_loc;
    PS::F32 xf32;
    PS::F32vec vf32;

    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++)psys[n].id = i;

    for(PS::S32 i=0; i<i_h; i++) finput>>xf32;
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++){
        finput>>psys[n].mass;
        //psys[n].mass *= (PS::MT::genrand_real1()+0.5);
    }
    for(PS::S32 i=i_t; i<n_glb; i++) finput>>xf32;

    for(PS::S32 i=0; i<i_h; i++) finput>>vf32;
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++){
        finput>>psys[n].pos;
        psys[n].pos += pos_shift;
    }
    for(PS::S32 i=i_t; i<n_glb; i++) finput>>vf32;

    for(PS::S32 i=0; i<i_h; i++) finput>>vf32;
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++) finput>>psys[n].vel;
    for(PS::S32 i=i_t; i<n_glb; i++) finput>>vf32;
    finput.close();
}

template<class Tpsys>
void WriteNemoAscii(const Tpsys & psys,
                    const PS::F32 time_sys,
                    const PS::S32 snp_id,
                    const char * dir_name){
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::S32 n_glb = 0;
    FPGrav * fp;
    PS::AllGatherParticle(fp, n_glb, &psys[0], n_loc);
    if(PS::Comm::getRank () == 0){
        const PS::S32 STRINGSIZE = 1024;
        char sout[STRINGSIZE];
        sprintf(sout,"%s/snap%04d.dat", dir_name, snp_id);
        for(int i=0;i<STRINGSIZE;i++)if(sout[i]==' ')sout[i]='0';
        std::ofstream foutput;
        foutput.open(sout);
        foutput<<std::setprecision(15);
        foutput<<n_glb<<std::endl;
        foutput<<"3"<<std::endl;
        foutput<<time_sys<<std::endl;
        for(PS::S32 i=0; i<n_glb; i++) foutput<<fp[i].mass<<std::endl;
        for(PS::S32 i=0; i<n_glb; i++) foutput<<fp[i].pos<<std::endl;
        for(PS::S32 i=0; i<n_glb; i++) foutput<<fp[i].vel<<std::endl;
        foutput.close();
    }
    delete [] fp;
}

void MakeColdUniformSphere(const double mass_glb,
                           const long long int n_glb,
                           const long long int n_loc,
                           double *& mass,
                           PS::F64vec *& pos,
                           PS::F64vec *& vel,
                           const double eng = -0.25,
                           const int seed = 0){
    
    assert(eng < 0.0);
    static const double PI = atan(1.0) * 4.0;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];
#pragma omp parallel
    {
        PS::MTTS mt;
        mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());
#pragma omp for
        for(int i=0; i<n_loc; i++){
            mass[i] = mass_glb / n_glb;
            //double m_tmp = PS::MT::genrand_real3();
            const double radius = 3.0;
            do {
                pos[i][0] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][1] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][2] = (2. * mt.genrand_res53() - 1.) * radius;
            }while(pos[i] * pos[i] >= radius * radius);
            while(1){
                const double v_max = 0.1;
                const double v_try = mt.genrand_res53();
                //const double v_crit = v_max * PS::MT::genrand_real3();
                const double v_crit = v_max * mt.genrand_res53();
                if(v_crit < v_try * v_try * pow( (1.0 - v_try * v_try), 3.5) ){
                    vel[i][0] = 0.0;
                    vel[i][1] = 0.0;
                    vel[i][2] = 0.0;
                    break;
                }
            }
        }
    }

    PS::F64vec cm_pos = 0.0;
    PS::F64vec cm_vel = 0.0;
    double  cm_mass = 0.0;
    for(size_t i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(size_t i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    /*
    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(size_t i=0; i<n_loc; i++){
        pos[i] *= r_scale;
        vel[i] *= coef;
    }
    */

}


template<class Tpsys>
void SetParticlesPlummer(Tpsys & psys,
                         const PS::S32 n_glb,
                         PS::S32 & n_loc,  
                         PS::F32 & t_sys){

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;

    const PS::F64 m_tot = 1.0;
    const PS::F64 eng = -0.25;
//    MakePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);
//    MakePlummerModelLike(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    MakeColdUniformSphere(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

template<class Tpsys>
void Kick(Tpsys & system,
          const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
        system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys>
void Drift(Tpsys & system,
           const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
        system[i].pos  += system[i].vel * dt;
    }
}

template<class Tpsys>
void CalcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
    if(clear){
        etot = ekin = epot = 0.0;
    }
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * system[i].pot;
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc = ekin_loc + epot_loc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI::COMM_WORLD.Allreduce(&etot_loc, &etot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&epot_loc, &epot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&ekin_loc, &ekin, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
    etot = etot_loc;
    epot = epot_loc;
    ekin = ekin_loc;
#endif
}

int main(int argc, char *argv[]){
    PS::F64 Tbegin = PS::GetWtime();
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);
    PS::F32 theta = 0.5;
    const PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F32 time_end = 10.0;
    char sinput[1024];
    char dir_name[1024];
    long long int n_tot = 1024;
    int c;
    sprintf(dir_name,"./result");
    while((c=getopt(argc,argv,"i:o:t:T:n:N:h")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 't':
            theta = atof(optarg);
            std::cerr<<"theta="<<theta<<std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr<<"time_end="<<time_end<<std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'N':
            n_tot = atoi(optarg);
            std::cerr<<"n_tot="<<n_tot<<std::endl;
            break;
        case 'h':
            std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
            std::cerr<<"t: theta (default: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
            std::cerr<<"n: n_group_limit (default: 64.0)"<<std::endl;
            std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
            return 0;
        }
    }

    std::ofstream fout_eng;
    std::ofstream fout_tcal;

#if 0
    if(PS::Comm::getRank() == 0){
#else
    if(1){
#endif
        char sout_de[1024];
        char sout_tcal[1024];
        sprintf(sout_de, "%s/t-de.dat", dir_name);
        sprintf(sout_tcal, "%s/t-tcal.dat", dir_name);
        std::cerr<<sout_de<<std::endl;
        std::cerr<<sout_tcal<<std::endl;
        fout_eng.open(sout_de);
        fout_tcal.open(sout_tcal);
    }


    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    PS::S32 n_grav_loc;
    PS::F32 time_sys;
#if 0
    PS::S32 n_grav_glb;
    ReadNemoAscii(system_grav, n_grav_glb, n_grav_loc, time_sys, sinput);
#else
    PS::S32 n_grav_glb = n_tot;
    SetParticlesPlummer(system_grav, n_grav_glb, n_grav_loc, time_sys);
#endif

    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);
    n_grav_loc = system_grav.getNumberOfParticleLocal();
#ifdef MONOPOLE
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::Monopole tree_grav;
#elif QUADRUPOLE
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::Quadrupole tree_grav;
#elif MONOPOLEGEOMETRICCENTER
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::MonopoleGeometricCenter tree_grav;
#elif DIPOLEGEOMETRICCENTER
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::DipoleGeometricCenter tree_grav;
#elif QUADRUPOLEGEOMETRICCENTER
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::QuadrupoleGeometricCenter tree_grav;
#endif

    tree_grav.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    
    tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo);

    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
    CalcEnergy(system_grav, Etot0, Ekin0, Epot0);

    const PS::F32 dt = 1.0/128.0;


    Kick(system_grav, dt*0.5);

    PS::F64 Tloop = 0.0;

    PS::S32 snp_id = 0;
    while(time_sys < time_end){

        PS::Timer timer;
        timer.reset();
        timer.start();
        if( fmod(time_sys, 1.0) == 0.0){
            //WriteNemoAscii(system_grav, time_sys, snp_id++, dir_name);
            FileHeader header;
            header.time = time_sys;
            header.n_body = system_grav.getNumberOfParticleLocal();
            struct stat st;
            if(stat(dir_name, &st) != 0) {
                PS::S32 rank = PS::Comm::getRank();
                PS::S32 ret_loc, ret;
                if(rank == 0)
                    ret_loc = mkdir(dir_name, 0777);
                PS::Comm::broadcast(&ret_loc, ret);
                if(ret == 0) {
                    if(rank == 0)
                        fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
                } else {
                    fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
                    PS::Finalize();
                    exit(0);
                }
            }
			char filename[256];
			sprintf(filename, "%s/%04d.dat", dir_name, snp_id++);
			system_grav.writeParticleAscii(filename, header);
        }
        timer.restart("WriteNemoAscii");
	
        time_sys += dt;
        Drift(system_grav, dt);

        timer.restart("Drift");
	
        if( fmod(time_sys, 1.0/32.0) == 0.0){
            dinfo.collectSampleParticle(system_grav, Tloop);
            timer.restart("collect");
            dinfo.decomposeDomain();
        }
        else{
            timer.restart("collect");
        }

        timer.restart("decompose");
	
        system_grav.exchangeParticle(dinfo);

        timer.restart("exchangeParticle");

        Tloop = PS::GetWtime();
	    
        tree_grav.calcForceAllAndWriteBackWithTimer
            (CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, timer, true);

        Tloop = PS::GetWtime() - Tloop;
	
	
        Kick(system_grav, dt*0.5);

        timer.stop("Kick");
	
        fout_tcal<<"time_sys= "<<time_sys<<std::endl;
        fout_tcal<<"tree_grav.getMemSizeUsed()= "<<tree_grav.getMemSizeUsed()*1e-9<<" [Gbyte]";
        fout_tcal<<" system_grav.getMemSizeUsed()= "<<system_grav.getMemSizeUsed()*1e-9<<" [Gbyte]"<<std::endl;
        tree_grav.dump_calc_cost(PS::Comm::getMaxValue(Tloop), fout_tcal);
        fout_tcal<<"Tloop= "<<Tloop<<" Ttot="<<PS::GetWtime()-Tbegin<<std::endl;
        timer.dump(fout_tcal);
        fout_tcal<<std::endl;

        CalcEnergy(system_grav, Etot1, Ekin1, Epot1);
        if(PS::Comm::getRank() == 0){
            fout_eng<<time_sys<<"   "<<(Etot1-Etot0)/Etot0<<std::endl;
        }
        Kick(system_grav, dt*0.5);
    }

    PS::Finalize();
    return 0;
}
