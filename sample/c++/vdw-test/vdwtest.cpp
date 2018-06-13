//
// vdwtest.cpp
//
// simple short-range MD test program
//
// Jun Makino March 9 2015
//
// Known problem as of Mar 13 2015
//   -- acc and phi are consistent only for the case of m=1
// This has been fixed as of Mar 15 2015


#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>


class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp){
        fscanf(fp, "%lf\n", &time);
        fscanf(fp, "%lld\n", &n_body);
        return n_body;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%e\n", time);
        fprintf(fp, "%lld\n", n_body);
    }
};


class FPLJ{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64 search_radius;
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
    PS::F64 getRSearch() const{
        return this->search_radius;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const FPLJ & force){
        acc = force.acc;
        pot = force.pot;
    }
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z);
	}

	void readAscii(FILE* fp){
		fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z);
	}
    void copyFromFP(const FPLJ & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
        search_radius = fp.search_radius;
    }
    
};


void CalcForceFpFp(const FPLJ * ep_i,
                   const PS::S32 n_ip,
                   const FPLJ * ep_j,
                   const PS::S32 n_jp,
                   FPLJ * force){
    const PS::F64 r0 = 3;
    const PS::F64 r0sq = r0*r0;
    const PS::F64 r0inv = 1/r0;
    const PS::F64 r0invp6 = 1/(r0sq*r0sq*r0sq);
    const PS::F64 r0invp7 = r0invp6*r0inv;
    const PS::F64 foffset = -12.0*r0invp6*r0invp7+6*r0invp7;
    const PS::F64 poffset = -13.0*r0invp6*r0invp6+7*r0invp6;
    for(PS::S32 i=0; i<n_ip; i++){
        PS::F64vec xi = ep_i[i].pos;
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        PS::S64 idi = ep_i[i].id;
        for(PS::S32 j=0; j<n_jp; j++){
            if( idi == ep_j[j].id ) continue;
            PS::F64vec rij = xi - ep_j[j].pos;
            PS::F64 r2 = rij * rij;
            if (r2 < r0sq){
                PS::F64 r_inv = 1.0/sqrt(r2);
                PS::F64 r = r2*r_inv;
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r6_inv = r2_inv * r2_inv * r2_inv;
                PS::F64 r12_inv = r6_inv * r6_inv;
                poti += r12_inv - r6_inv-foffset*r+poffset;
                ai += (12*r12_inv*r2_inv - 6*r6_inv*r2_inv+ foffset*r_inv)
                    * rij;
            }
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}


template<class Tpsys>
void ReadNemoAscii(Tpsys & psys,
                   PS::S32 & n_glb,
                   PS::S32 & n_loc,  
                   PS::F64 & t_sys,
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

    PS::F64vec pos_shift = 0.0;

    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    const PS::S32 i_t = i_h+n_loc;
    PS::F64 xf32;
    PS::F64vec vf32;

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
                    const PS::F64 time_sys,
                    const PS::S32 snp_id,
                    const char * dir_name){
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::S32 n_glb = 0;
    FPLJ * fp;
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

void MakeUniformCube(const double mass_glb,
                           const long long int n_glb,
                           const long long int n_loc,
                           double *& mass,
                           PS::F64vec *& pos,
                           PS::F64vec *& vel,
                           const double eng = -0.25,
                           const int seed = 0){
    
    assert(eng < 0.0);
    //static const double PI = atan(1.0) * 4.0;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];

    PS::MTTS mt;

    if (n_loc != n_glb){
        std::cerr << "MakeUniformCube: n_loc and n_glb must be the same" << n_loc << "!= " << n_glb <<std::endl;
        PS::Abort();
    }
    int n_1d = pow(n_loc+0.0,0.33333334)+0.1;
    std::cerr << "MakeUniformCube: n_loc="<<n_loc
              << " and n_1d="<<n_1d<<std::endl;
    if (n_loc != n_1d*n_1d*n_1d){
        std::cerr << "MakeUniformCube: n_loc and n_1d^3 must be the same"
                  << n_loc << "!= " << n_1d*n_1d*n_1d <<std::endl;
        PS::Abort();
    }
    
    mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());

    //double offset = (n_1d -1.0)*0.5;
    int ip=0;
    for(int i=0; i<n_1d; i++){
        for(int j=0; j<n_1d; j++){
            for(int k=0; k<n_1d; k++){
                pos[ip][0] = i;
                pos[ip][1] = j;
                pos[ip][2] = k;
                ip++;
            }
        }
    }
    
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        const double v_max = 0.1;
        do {
            vel[i][0] = (2. * mt.genrand_res53() - 1.) * v_max;
            vel[i][1] = (2. * mt.genrand_res53() - 1.) * v_max;
            vel[i][2] = (2. * mt.genrand_res53() - 1.) * v_max;
        } while(vel[i] * vel[i] >= v_max * v_max);
    }

    PS::F64vec cm_pos = 0.0;
    PS::F64vec cm_vel = 0.0;
    double  cm_mass = 0.0;
    for(int i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(int i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }
}


template<class Tpsys>
void SetParticlesPlummer(Tpsys & psys,
                         const PS::S32 n_glb,
                         PS::S32 & n_loc,  
                         PS::F64 & t_sys){

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;

    const PS::F64 m_tot = n_glb;
    const PS::F64 eng = -0.25;
    MakeUniformCube(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(int i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
        psys[i].search_radius = 4.0;
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
           const PS::F64 dt,
           const PS::F64 boxdh){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
        system[i].pos  += system[i].vel * dt;
	for(int k=0;k<3;k++){
	    if (system[i].pos[k] < -boxdh)system[i].pos[k] += boxdh*2;
	    if (system[i].pos[k] > boxdh)system[i].pos[k] -= boxdh*2;
	}
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
    PS::F64 theta = 0.5;
    const PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F64 time_end = 10.0;
    PS::F64 dt = 1.0/128.0;
    PS::F64 dt_diag = 1.0/16.0;
    PS::F64 boxdh=5;
    //char sinput[1024];
    char dir_name[1024];
    long long int n_tot = 512;
    int c;
    sprintf(dir_name,"./result");
    while((c=getopt(argc,argv,"b:d:D:i:o:t:T:n:N:h")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 'b':
            boxdh = atof(optarg);
            std::cerr<<"boxdh="<<boxdh<<std::endl;
            break;
        case 'd':
            dt = atof(optarg);
            std::cerr<<"dt="<<dt<<std::endl;
            break;
        case 'D':
            dt_diag = atof(optarg);
            std::cerr<<"dt_diag="<<dt_diag<<std::endl;
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
            std::cerr<<"b: half of the box size (default: 5.0)"<<std::endl;
            std::cerr<<"d: time step (default: 1/128)"<<std::endl;
            std::cerr<<"D: time step for diag(default: 1/16)"<<std::endl;
            std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
            std::cerr<<"t: theta (default: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
            std::cerr<<"n: n_group_limit (default: 64.0)"<<std::endl;
            std::cerr<<"N: n_tot (default: 512)"<<std::endl;
            return 0;
        }
    }


    struct stat st;
    if (stat(dir_name, &st) != 0) {
        PS::S32 rank = PS::Comm::getRank();
        PS::S32 ret_loc=0, ret;
        if (rank == 0)
            ret_loc = mkdir(dir_name, 0777);
        ret=ret_loc;
        PS::Comm::broadcast(&ret, 1);
        if (ret == 0) {
            if(rank == 0)
                fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
        } else {
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
            PS::Abort();
            exit(0);
        }
    }
    
    std::ofstream fout_eng;
    std::ofstream fout_tcal;
    {
        char sout_de[1024];
        char sout_tcal[1024];
        sprintf(sout_de, "%s/t-de.dat", dir_name);
        sprintf(sout_tcal, "%s/t-tcal.dat", dir_name);
        std::cerr<<sout_de<<std::endl;
        std::cerr<<sout_tcal<<std::endl;
        fout_eng.open(sout_de);
        fout_tcal.open(sout_tcal);
    }


    PS::ParticleSystem<FPLJ> system_grav;
    system_grav.initialize();
    PS::S32 n_grav_loc;
    PS::F64 time_sys;
    PS::F64 time_diag;
#if 0
    PS::S32 n_grav_glb;
    ReadNemoAscii(system_grav, n_grav_glb, n_grav_loc, time_sys, sinput);
#else
    PS::S32 n_grav_glb = n_tot;
    SetParticlesPlummer(system_grav, n_grav_glb, n_grav_loc, time_sys);
#endif

    const PS::F64 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(PS::F64vec(-boxdh,-boxdh,-boxdh),
                           PS::F64vec(boxdh,boxdh,boxdh));
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);
    n_grav_loc = system_grav.getNumberOfParticleLocal();
    PS::TreeForForceShort<FPLJ, FPLJ, FPLJ>::Scatter tree_grav;

    tree_grav.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    
    tree_grav.calcForceAllAndWriteBack(CalcForceFpFp,  system_grav, dinfo);

    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
    CalcEnergy(system_grav, Etot0, Ekin0, Epot0);



    Kick(system_grav, dt*0.5);

    PS::F64 Tloop = 0.0;

    PS::S32 snp_id = 0;
    time_diag = time_sys+dt_diag;
    while(time_sys < time_end){

        PS::Timer timer;
        timer.reset();
        timer.start();
        if( fmod(time_sys, 1.0) == 0.0){
            //WriteNemoAscii(system_grav, time_sys, snp_id++, dir_name);
            FileHeader header;
            header.time = time_sys;
            header.n_body = system_grav.getNumberOfParticleLocal();
            char filename[256];
            sprintf(filename, "%s/%04d.dat", dir_name, snp_id++);
            system_grav.writeParticleAscii(filename, header);
        }
        timer.restart("WriteNemoAscii");

        time_sys += dt;
        Drift(system_grav, dt,boxdh);

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
        //tree_grav.calcForceAllAndWriteBackWithTimer
        //            (CalcForceEpEp(),  system_grav, dinfo, timer, true);
        tree_grav.calcForceAllAndWriteBack
                    (CalcForceFpFp,  system_grav, dinfo);
        tree_grav.calcForceAllAndWriteBack
                    (CalcForceFpFp,  system_grav, dinfo, true);
        Tloop = PS::GetWtime() - Tloop;


        Kick(system_grav, dt*0.5);

        timer.stop("Kick");

        fout_tcal<<"time_sys= "<<time_sys<<std::endl;
        fout_tcal<<"tree_grav.getMemSizeUsed()= "<<tree_grav.getMemSizeUsed()*1e-9<<" [Gbyte]";
        fout_tcal<<" system_grav.getMemSizeUsed()= "<<system_grav.getMemSizeUsed()*1e-9<<" [Gbyte]"<<std::endl;
//      tree_grav.dump_calc_cost(PS::Comm::getMaxValue(Tloop), fout_tcal);
        fout_tcal<<"Tloop= "<<Tloop<<" Ttot="<<PS::GetWtime()-Tbegin<<std::endl;
        timer.dump(fout_tcal);
        fout_tcal<<std::endl;

        CalcEnergy(system_grav, Etot1, Ekin1, Epot1);
        if (PS::Comm::getRank() == 0){
            fout_eng<<time_sys<<"   "<<(Etot1-Etot0)/Etot0<<std::endl;
            if (time_sys >= time_diag) {
                fprintf(stderr, "time: %10.7f energy error: %+e\n",
                        time_sys, (Etot1 - Etot0) / Etot0);
                time_diag += dt_diag;
            }            
        }
        Kick(system_grav, dt*0.5);
    }

    PS::Finalize();
    return 0;
}
