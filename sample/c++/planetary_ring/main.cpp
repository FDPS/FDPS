#include<iostream>
#include<cstdio>
#include<unistd.h>
#include<random>
#include<particle_simulator.hpp>
#include"user_defined.hpp"
#include"./fdps-util/my_lib.hpp"
#include"./fdps-util/kepler.hpp"


PS::F64 FP_t::kappa;
PS::F64 FP_t::eta;
PS::F64 FP_t::eps;

using namespace MY_LIB::LITERAL;
void DEBUG_PRINT_RING(){
#if defined(DEBUG_PRINT_RING)
    PS::Comm::barrier();
    std::cerr<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
}

constexpr PS::F64 MY_PI = MY_LIB::CONSTANT::pi;
FP_t PLANET;

template<typename Tpsys, typename Tsat>
void RemoveParticlesAroundSatellite(Tpsys & psys, Tsat & sat){
    const auto n_loc = psys.getNumberOfParticleLocal();
    const auto n_sat = sat.getNumberOfSat();
    std::vector<PS::S32> id_remove;
    for(auto i=0; i<n_loc; i++){
        const auto pos_ptcl = psys[i].pos_car;
	const auto r_coll_ptcl = psys[i].r_coll;
	for(auto j=0; j<n_sat; j++){
	    const auto rij = pos_ptcl -sat[j].pos_car;
	    const auto r_coll = sat[j].r_coll + r_coll_ptcl;
	    if( (rij*rij) < (r_coll*r_coll) ){
	        id_remove.push_back(i);
	    }
	}
    }
    if(id_remove.size() > 0)
        psys.removeParticle(&id_remove[0], id_remove.size());

    if(id_remove.size() > 0){
        std::cerr<<"rank= "<<PS::Comm::getRank()<<" id_remove.size()= "<<id_remove.size()<<std::endl;
    }
}

template<typename Tpsys>
void UpdateCyl(Tpsys & psys){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
      psys[i].pos_cyl = ConvertCar2Cyl(psys[i].pos_car);
    }
}


template<typename Tpsys>
void Rotate(Tpsys & psys, const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        const auto pos = psys[i].pos_car;
        const auto vel = psys[i].vel;
        const auto r = sqrt(pos*pos);
        const auto mm = sqrt(PLANET.mass/(r*r*r)); // mean mortion
        const auto theta = mm*dt;
        const auto cth = std::cos(theta);
        const auto sth = std::sin(theta);
        const auto l = pos ^ vel;
        const auto n = l / sqrt(l*l);

        const auto ex_pos = pos / sqrt(pos*pos);
        const auto ey_pos = n ^ ex_pos;
        const auto pos_new = sqrt(pos*pos)*(cth*ex_pos + sth*ey_pos);

        const auto ex_vel = vel / sqrt(vel*vel);
        const auto ey_vel = n ^ ex_vel;
        const auto vel_new = sqrt(vel*vel)*(cth*ex_vel + sth*ey_vel);
        //const auto l_new = pos_new ^ vel_new;

        psys[i].pos_car = pos_new;
	psys[i].pos_cyl = ConvertCar2Cyl(pos_new);
        psys[i].vel     = vel_new;
    }
}

template<typename Tsys>
void RigidRotation(Tsys & sys, const PS::F64 dt){
    const auto n = sys.getNumberOfParticleLocal();
    const auto ax = 1.0;
    const auto mm = sqrt(PLANET.mass/(ax*ax*ax)); // mean mortion
    const auto theta = mm*dt;
    const auto cth = std::cos(theta);
    const auto sth = std::sin(theta);    
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        const auto x_new  = cth*sys[i].pos_car.x - sth*sys[i].pos_car.y;
        const auto y_new  = sth*sys[i].pos_car.x + cth*sys[i].pos_car.y;
        const auto vx_new = cth*sys[i].vel.x - sth*sys[i].vel.y;
        const auto vy_new = sth*sys[i].vel.x + cth*sys[i].vel.y;
        sys[i].pos_car = PS::F64vec(x_new, y_new, sys[i].pos_car.z);
	sys[i].pos_cyl = ConvertCar2Cyl(sys[i].pos_car);
        sys[i].vel     = PS::F64vec(vx_new, vy_new, sys[i].vel.z);
    }
}



void DivideNProc(PS::S32 & nx,
                 PS::S32 & ny,
                 PS::S32 & nz,
                 const PS::S32 n_proc,
                 const PS::F64 delta_ax){
    nz = 1;
    ny = 1;
    nx = n_proc / ny;
    if(n_proc == 1) return;
    const PS::F64 dx = 2.0*MY_PI   / nx;
    double dy = delta_ax / ny;
    double ratio = (dx < dy) ? dx / dy : dy / dx;
    PS::S32 ny_tmp = ny;
    PS::S32 nx_tmp = nx;
    double dx_tmp = dx;
    double dy_tmp = dy;
    double ratio_tmp = ratio;
    do{
        ny = ny_tmp;
        nx = nx_tmp;
        ratio = ratio_tmp;
        ny_tmp += 1;
        while( n_proc % ny_tmp != 0) ny_tmp++;
        nx_tmp = n_proc / ny_tmp;
        dx_tmp = 2.0*MY_PI   / nx_tmp;
        dy_tmp = delta_ax / ny_tmp;
        ratio_tmp = (dx_tmp < dy_tmp) ? dx_tmp / dy_tmp : dy_tmp / dx_tmp;
    }while( fabs(ratio_tmp-1.0) < fabs(ratio-1.0));
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
    }
    assert(n_proc == nx*ny*nz);
}

PS::F64ort GetPosDomainCyl(const PS::F64 delta_ax,
			   const PS::S32 nx,
			   const PS::S32 ny){
    constexpr PS::F64 len_x = 2.0 * MY_PI;
    PS::F64ort pos_domain;
    const auto my_rank = PS::Comm::getRank();
    const auto rank_x = my_rank / ny;
    const auto rank_y = my_rank % ny;
    const auto dx = len_x / nx;
    const auto dy = delta_ax / ny;
    const auto dy_offset = 1.0-delta_ax*0.5;
    const auto x_offset = -MY_PI;
    pos_domain.low_.x = dx*rank_x + x_offset;
    pos_domain.low_.y = dy*rank_y + dy_offset;
    pos_domain.low_.z = -MY_PI;
    pos_domain.high_.x = dx*(rank_x+1) + x_offset;
    pos_domain.high_.y = dy*(rank_y+1) + dy_offset;
    pos_domain.high_.z = MY_PI;
    return pos_domain;
}

PS::F64vec GetVel(const PS::F64vec & pos){
    const PS::F64 dr = sqrt(pos*pos);
    const PS::F64 v_kep = sqrt(1.0/(dr));
    const PS::F64 theta = atan2(pos.y, pos.x) + MY_PI*0.5;
    const PS::F64 vx = v_kep * cos(theta);
    const PS::F64 vy = v_kep * sin(theta);
    const PS::F64 vz = 0.0;
    return PS::F64vec(vx, vy, vz);
}

#if defined(QUAD)
using SPJ_t    = MySPJQuadrupole;
using Moment_t = MyMomentQuadrupole;
using CalcForceSp = CalcForceSpQuad<EPI_t, SPJ_t, Force_t>;
#else
using SPJ_t    = MySPJMonopole;
using Moment_t = MyMomentMonopole;
using CalcForceSp = CalcForceSpMono<EPI_t, SPJ_t, Force_t>;
#endif

using MY_SEARCH_MODE = PS::SEARCH_MODE_LONG_SCATTER;
using Tree_t = PS::TreeForForce<MY_SEARCH_MODE, Force_t, EPI_t, EPJ_t, Moment_t, Moment_t, SPJ_t, PS::CALC_DISTANCE_TYPE_NEAREST_X>;

template<typename Tpsys>
void CalcForceFromPlanet(Tpsys & psys, const FP_t & pla){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        const auto rij = psys[i].pos_car - pla.pos_car;
        const auto r_sq = rij*rij;
        const auto r_inv = 1.0 / sqrt(r_sq);
        const auto pot   = pla.mass * r_inv;
        psys[i].acc -= pot * r_inv * r_inv * rij;
        psys[i].pot -= pot;
    }
}

template<typename Tpsys>
void Kick(Tpsys & psys, const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        psys[i].vel += psys[i].acc*dt;
        psys[i].vel_full = psys[i].vel + psys[i].acc*dt;
    }
}
template<typename Tpsys>
void Drift(Tpsys & psys,
           const PS::F64 dt){
    const auto n = psys.getNumberOfParticleLocal();
#pragma omp parallel for
    for(auto i=0; i<n; i++){
        psys[i].pos_car += psys[i].vel*dt;
    }
}



template<typename Tpsys>
void SetID(Tpsys & psys, const PS::S32 start_id=0){
    PS::S32 n = psys.getNumberOfParticleLocal();
    PS::S32 n_proc  = PS::Comm::getNumberOfProc();
    PS::S32 my_rank = PS::Comm::getRank();
    PS::ReallocatableArray<PS::S32> n_ar;
    n_ar.resizeNoInitialize(n_proc);
    PS::Comm::allGather(&n, 1, n_ar.getPointer());
    PS::S32 offset = start_id;
    for(auto i=0; i<my_rank; i++){
        offset += n_ar[i];
    }
    for(auto i=0; i<n; i++){
        psys[i].id = i + offset;
    }
}


class DiskInfo{
    PS::F64 ax_;
    PS::F64 delta_ax_;
    PS::F64 r_hill_;
    PS::F64 r_phy_;
    PS::F64 dens_;
    PS::F64 m_ptcl_;
    PS::F64 kappa_;
    PS::F64 eta_;
    PS::F64 e_refl_;
    PS::F64 t_dur_;
    PS::F64 tau_;
    PS::S64 n_glb_;
    void setParticlesOnLayer(std::vector<PS::F64vec> & pos, const PS::S64 id_x_head, const PS::S64 id_x_tail, const PS::S64 id_y_head, const PS::S64 id_y_tail, const PS::F64 dtheta, const PS::F64 dr, const PS::F64ort box,
			     const PS::F64 offset_theta = 0.0, const PS::F64 offset_r = 0.0, const PS::F64 offset_z = 0.0, const PS::F64 eps = 0.0){
	static std::mt19937 mt(PS::Comm::getRank());
	static std::uniform_real_distribution<double> dist(0.0,1.0);
	PS::F64 pos_z = offset_z;
	for(auto i=id_x_head-3; i<=id_x_tail+3; i++){
	    PS::F64 pos_x = ((PS::F64)i)*dtheta + offset_theta;
	    if(pos_x < box.low_.x || pos_x >= box.high_.x) continue;
	    for(auto j=id_y_head-3; j<=id_y_tail+3; j++){
		PS::F64 pos_y = j*dr + offset_r;
		if(pos_y < box.low_.y || pos_y >= box.high_.y) continue;
		PS::F64 eps_x = eps*(dist(mt)-0.5)*2.0;
		PS::F64 eps_y = eps*(dist(mt)-0.5)*2.0;
		PS::F64 eps_z = eps*(dist(mt)-0.5)*2.0;
		pos.push_back(PS::F64vec(pos_x+eps_x, pos_y+eps_y, pos_z+eps_z));
	    }
	}
    }

    void calcKappa() {
	kappa_ = std::pow((2.0*MY_PI/t_dur_), 2.0); // k^{'};
    }
    void calcEta() {
	const auto ln_e_refl = std::log(e_refl_);
	eta_ = 4.0*MY_PI/(t_dur_*std::sqrt(1.0+std::pow((MY_PI/ln_e_refl),2.0))); // \eta^{'}
	//return 4.0*MY_PI*ln_e_refl/(t_dur_*std::sqrt(MY_PI*MY_PI+ln_e_refl)); // \eta^{'}
    }
    
public:
    void setParams(const PS::F64 delta_ax, const PS::F64 e_refl, const PS::F64 t_dur, const PS::F64 tau, const PS::F64 rphy_over_rhill, const PS::S64 n_glb){
	ax_ = 1.0;
	delta_ax_ = delta_ax;
	e_refl_ = e_refl;
	t_dur_  = t_dur;
	tau_    = tau;
	n_glb_ = n_glb;
	PS::F64 ax_in  = ax_ - 0.5*delta_ax_;
	PS::F64 ax_out = ax_ + 0.5*delta_ax_;
	r_phy_  = sqrt(tau_*(ax_out*ax_out - ax_in*ax_in) / n_glb_);
	r_hill_ = r_phy_ / rphy_over_rhill;
	m_ptcl_ = (r_hill_/((ax_out+ax_in)*0.5))*(r_hill_/((ax_out+ax_in)*0.5))*(r_hill_/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5;
	dens_ = m_ptcl_ * n_glb / (MY_PI*(ax_out*ax_out - ax_in*ax_in));
	calcKappa();
	calcEta();
    }
    DiskInfo(){}    
    DiskInfo(const PS::F64 delta_ax, const PS::F64 e_refl, const PS::F64 t_dur, const PS::F64 tau, const PS::F64 rphy_over_rhill, const PS::S64 n_glb){
	setParams(delta_ax, e_refl, t_dur, tau, rphy_over_rhill, n_glb);
    }
    void writeAscii(FILE * fp){
        fprintf(fp, "%12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e %lld \n",
		ax_, delta_ax_, r_hill_, r_phy_, dens_, m_ptcl_, kappa_, eta_, e_refl_, t_dur_, tau_, n_glb_);
    }
    void readAscii(FILE * fp){
        auto tmp = fscanf(fp, "%lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf %lld",
			  &ax_, &delta_ax_, &r_hill_, &r_phy_, &dens_, &m_ptcl_, &kappa_, &eta_, &e_refl_, &t_dur_, &tau_, &n_glb_);
	assert(tmp == 12);
	FP_t::kappa = kappa_;
	FP_t::eta   = eta_;
    }


    void writeBinary(FILE * fp){
	fwrite(this, sizeof(*this), 1, fp);
    }
    void readBinary(FILE * fp){
	auto tmp = fread(this, sizeof(*this), 1, fp);
	assert(tmp == 1);
	FP_t::kappa = kappa_;
	FP_t::eta   = eta_;
    }
    
    template<class Tpsys>
    PS::S64 setParticles(Tpsys & psys, const PS::F64ort box, const bool layer=true, const bool random_shift=true){
	const auto ax_in  = ax_ - 0.5*delta_ax_;
	const auto ax_out = ax_ + 0.5*delta_ax_;
	const auto area = MY_PI*(ax_out*ax_out-ax_in*ax_in);
	PS::F64	dS = area / n_glb_;
	if(layer){
	    dS *= 3.0;
	}
	const auto dl = sqrt(dS);
	assert(delta_ax_ > dl);
	const auto dz = dl;
	const PS::S64 n_theta = (2.0*MY_PI*ax_ / dl);
	const auto dtheta = 2.0*MY_PI / n_theta;
	const PS::S64 n_r = (delta_ax_ / dl);
	const auto dr = delta_ax_ / n_r;

	const PS::S64 id_x_head = box.low_.x / dtheta;
	const PS::S64 id_x_tail = box.high_.x / dtheta;
	const PS::S64 id_y_head = box.low_.y / dr;
	const PS::S64 id_y_tail = box.high_.y / dr;

	PS::F64 eps = 0.0;
	if(random_shift){
	    eps = (dl > 2.0*r_phy_) ? (0.5*(dl-2.0*r_phy_))*0.9 : 0.0;
	}
    
	PS::Comm::barrier();
	if(PS::Comm::getRank()==0){
	    std::cerr<<"delta_ax_= "<<delta_ax_
		     <<" dS = "<<dS
		     <<" dl = "<<dl
		     <<" n_r= "<<n_r
		     <<" dr= "<<dr
		     <<" eps= "<<eps
		     <<std::endl;
	    std::cerr<<"n_theta= "<<n_theta
		     <<" dtheta= "<<dtheta
		     <<std::endl;
	}
	PS::Comm::barrier();
	std::vector<PS::F64vec> pos;
	PS::F64 offset_theta = dtheta*(sqrt(2.0)-1.0)*0.5;
	PS::F64 offset_r = 0.0;
	PS::F64 offset_z = 0.0;
	setParticlesOnLayer(pos, id_x_head, id_x_tail, id_y_head, id_y_tail, dtheta, dr, box, offset_theta, offset_r, offset_z, eps);
	if(layer == true){
	    offset_theta = 0.5*dtheta + dtheta*(sqrt(2.0)-1.0)*0.5;
	    offset_r = 0.5*dr;
	    offset_z = -0.5*dz;
	    setParticlesOnLayer(pos, id_x_head, id_x_tail, id_y_head, id_y_tail, dtheta, dr, box, offset_theta, offset_r, offset_z, eps);
	    offset_z = 0.5*dz;
	    setParticlesOnLayer(pos, id_x_head, id_x_tail, id_y_head, id_y_tail, dtheta, dr, box, offset_theta, offset_r, offset_z, eps);

	}
	
	PS::S64 n_loc = pos.size();
	psys.setNumberOfParticleLocal(n_loc);
	for(PS::S64 i=0; i<n_loc; i++){
	    psys[i].pos_cyl = pos[i];
	    psys[i].pos_car = ConvertCyl2Car(pos[i]);
	    psys[i].vel     = GetVel(psys[i].pos_car);
	    psys[i].mass    = m_ptcl_;
	    psys[i].r_coll  = r_phy_;
	    psys[i].r_search = 6.0*r_hill_;
	    psys[i].eps = 0.0;
	    psys[i].kappa = kappa_;
	    psys[i].eta = eta_;
	}
	return n_loc;
    }
};


int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    
    PS::Initialize(argc, argv);
    const auto my_rank = PS::Comm::getRank();
    const auto n_proc  = PS::Comm::getNumberOfProc();

    PLANET.mass = 1.0;
    PLANET.pos_car = PLANET.vel = 0.0;

    MY_LIB::UnitManager unit_man;
    unit_man.initializeByLMG(1e5_km2cm, MY_LIB::CONSTANT::mass_saturn_cgs, 1.0);
    unit_man.dump(std::cerr);

    FDPS_UTIL::ComandLineOptionManager cmd(argc, argv);
    cmd.append("delta_ax",  "a", "the width of semi^major axis", "1.0e-3");
    cmd.append("inv_dt",    "d", "inverse of dt", "256");
    cmd.append("theta",     "t", "iopening criterion", "0.5");
    cmd.append("time_end",  "T", "ending time", "10.0");
    cmd.append("n_leaf_limit",    "l", "n_leaf_limit", "8");
    cmd.append("n_group_limit",     "n", "maximum # of particles in a leaf cell", "64");
    cmd.append("n_step_share",  "L", "# of steps sharing the same interaction list", "1");
    cmd.append("n_loc",     "N", "# of particles per process", "32768");
    cmd.append("ex_let_mode",   "e", "exchange LET mode \n  0: PS::EXCHANGE_LET_A2A \n  1: PS::EXCHANGE_LET_P2P_EXACT \n  2: PS::EXCHANGE_LET_P2P_FAST \n", "0");
    cmd.append("e_refl",    "E", "restitution coefficient", "0.5");
    cmd.append("t_dur",     "D", "inverse of duration time", "4");
    cmd.append("tau",       "O", "optical depth", "1.0");
    cmd.append("rphy_over_rhill", "R", "rphy / rhill", "1.0");
    cmd.appendFlag("flag_mtm", "m", "use mid point tree method");
    cmd.appendFlag("flag_rot", "r", "use rotation method");
    cmd.appendFlag("flag_para_out", "output files in parallel");
    cmd.appendFlag("flag_bin_out", "output binary format");
    cmd.appendNoDefault("read_file",  "i", "write file name base", false);
    cmd.appendNoDefault("write_file", "o", "write file name base", false);
    cmd.append("sat_mode",  "satellite mode\n  0: no satellite\n  1: PAN\n  2: few satellites\n", "0");
    cmd.append("sat_mass_ratio",  "satellite-mass / particle-mass (this option is available ONLY IF sat_mode=2)", "8.0");
    cmd.append("sat_num",  "satellite number per process (this option is available ONLY IF sat_mode=2)", "10");
    cmd.append("n_smp",  "# of sample particles per process", "100");
    cmd.append("dt_snp",  "the interval time of snapshot", "1.0");
    cmd.appendNoDefault("log_file",  "logfile name", false);
    cmd.read();
    
    Param param;
    param.setTimeSys(0.0);
    PS::S64 id_snp = 0;
    PS::F64 delta_ax = cmd.get("delta_ax");
    PS::S32 inv_dt = cmd.get("inv_dt");
    param.setDt( 1.0 / inv_dt );
    PS::F64 theta    = cmd.get("theta");
    PS::F64 time_end = cmd.get("time_end");
    PS::S32 n_leaf_limit = cmd.get("n_leaf_limit");
    PS::S32 n_group_limit = cmd.get("n_group_limit");
    PS::S32 n_step_share = cmd.get("n_step_share");
    PS::S64 n_loc        = cmd.get("n_loc");
    PS::S64 n_glb        = n_loc * PS::Comm::getNumberOfProc();
    PS::S32 ex_let_mode_n = cmd.get("ex_let_mode");
    auto ex_let_mode = PS::EXCHANGE_LET_A2A;
    if(ex_let_mode_n == 0)
	ex_let_mode = PS::EXCHANGE_LET_A2A;
    else if(ex_let_mode_n == 1)
	ex_let_mode = PS::EXCHANGE_LET_P2P_EXACT;
    else if(ex_let_mode_n == 2)
	ex_let_mode = PS::EXCHANGE_LET_P2P_FAST;
    PS::F64 e_refl = cmd.get("e_refl");
    PS::F64 t_dur  = cmd.get("t_dur");
    PS::F64 tau  = cmd.get("tau");
    PS::F64 rphy_over_rhill  = cmd.get("rphy_over_rhill");
    bool    flag_mtm = cmd.get("flag_mtm");
    bool    flag_rot = cmd.get("flag_rot");
    bool    flag_para_out = cmd.get("flag_para_out");
    bool    flag_bin_out = cmd.get("flag_bin_out");
    PS::S32 sat_mode = cmd.get("sat_mode");
    PS::F64 sat_mass_ratio = cmd.get("sat_mass_ratio");
    PS::S32 sat_num = cmd.get("sat_num");
    PS::S32 n_smp    = cmd.get("n_smp");
    PS::F64 dt_snp = cmd.get("dt_snp");


    std::string write_file_name_base;
    bool flag_snapshot   = false; // output snapshot    
    if(cmd.hasOption("write_file") || cmd.hasOption("o")){
	write_file_name_base = cmd.get<std::string>("write_file");
	flag_snapshot   = true;
    }
    
    std::string read_file_name_base;
    bool flag_start_file = false; // resume simulations using snapshot
    if(cmd.hasOption("read_file") || cmd.hasOption("i")){
	read_file_name_base = cmd.get<std::string>("read_file");
	flag_start_file = true;
    }

    std::ofstream fout_log;
    if(cmd.hasOption("log_file")){
	std::string log_file_name = cmd.get<std::string>("log_file");
	//std::ofstream tmp;
	//tmp.open(log_file_name);
	//fout_log.rdbuf(tmp.rdbuf);
	fout_log.open(log_file_name);
    }
    
    
    PS::S32 nx, ny, nz;
    DivideNProc(nx, ny, nz, n_proc, delta_ax);
    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.setDomain(nx, ny, nz);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setPosRootDomainX(-MY_PI, MY_PI);
    PS::ParticleSystem<FP_t> system;
    system.initialize();
    system.setAverageTargetNumberOfSampleParticlePerProcess(n_smp);
    const auto pos_domain = GetPosDomainCyl(delta_ax, nx, ny);
    if(my_rank==0){
	fout_log<<"my_rank= "<<my_rank<<" pos_domain= "<<pos_domain<<std::endl;
    }
    DiskInfo disk_info;
    SatelliteSystem sat_system_glb;

    FDPS_UTIL::SnapshotManager snapshot_manager(write_file_name_base, dt_snp, flag_para_out, flag_bin_out); // single ascii
    Energy eng_init;
    PS::F64 eng_disp_glb = 0;
    if(flag_start_file){
	PS::F64 t;
	snapshot_manager.read< PS::ParticleSystem<FP_t>, DiskInfo, Energy, SatelliteSystem>(read_file_name_base, system, t, disk_info, eng_init, sat_system_glb);
	param.setTimeSys(t);
	eng_disp_glb = eng_init.getDisp();
    } else {
	disk_info.setParams(delta_ax, e_refl, t_dur, tau, rphy_over_rhill, n_glb);
	n_loc = disk_info.setParticles(system, pos_domain);
	n_glb = system.getNumberOfParticleGlobal();
	for(auto i=0; i<n_loc; i++){
	    if(!dinfo.getPosRootDomain().contained(system[i].getPos())){
		std::cerr<<"n_loc= "<<n_loc<<std::endl;
		std::cerr<<"dinfo.getPosRootDomain()= "<<dinfo.getPosRootDomain()<<std::endl;
		std::cerr<<"i= "<<i<<" system[i].getPos()= "<<system[i].getPos()<<std::endl;
		std::cerr<<"i= "<<i<<" system[i].pos_car= "<<system[i].pos_car<<std::endl;
		std::cerr<<"i= "<<i<<" system[i].pos_cyl= "<<system[i].pos_cyl<<std::endl;
	    }
	    assert(dinfo.getPosRootDomain().contained(system[i].getPos()));
	}
	snapshot_manager.setIdSnp(id_snp);
	snapshot_manager.setTimeSnp(param.getTimeSys());

	SatelliteSystem sat_system_loc;
	sat_system_loc.clear();
	if (sat_mode == 1){
	  if(my_rank==0){
		FP_t sat_tmp = system[0];
		using namespace MY_LIB::LITERAL;
		PS::F64 ax  = 1.0;
		// Daphnis
		/*
		PS::F64 sat_mass = 1e-13;
		PS::F64 ecc = 3e-5;
		PS::F64 inc = 4e-3_deg;
		*/
		// PAN
		PS::F64 sat_mass = 1e-11;
		PS::F64 ecc = 1e-5;
		PS::F64 inc = 1e-4_deg;
                /*
		PS::F64 sat_mass = system[0].mass * 64.0;
		PS::F64 ecc = 0.0;
		PS::F64 inc = 0.0;
                */
		PS::F64 OMG = M_PI; // pericenter is in the direction of -x
		PS::F64 omg = 0.0;
		PS::F64 u = M_PI; // satellite is at the apocenter
		MY_LIB::OrbParam2PosVel(PLANET.pos_car, sat_tmp.pos_car, PLANET.vel, sat_tmp.vel, PLANET.mass, 0.0,
					ax, ecc, inc, OMG, omg, u);
		sat_tmp.mass = system[0].mass;
		sat_tmp.changeMass(sat_mass);
		std::cerr<<"sat_tmp.pos_car= "<<sat_tmp.pos_car
			 <<" vel= "<<sat_tmp.vel
			 <<" mass= "<<sat_tmp.mass
			 <<std::endl;
		std::cerr<<"system[0].mass= "<<system[0].mass<<std::endl;
		sat_tmp.pos_cyl = ConvertCar2Cyl(sat_tmp.pos_car);
		sat_system_loc.setSatellite(sat_tmp);	    
	    }
	} else if (sat_mode == 2){
            auto n_loc_tmp = system.getNumberOfParticleLocal();
            assert(sat_num < n_loc_tmp);
	    const PS::S32 n_remove = sat_num;
	    PS::S32 id_remove[n_remove];
            MY_LIB::GetUniqueID(id_remove, n_remove, 0, n_loc_tmp, 0);
            /*
            if(my_rank==0){
                std::cerr<<"n_remove= "<<n_remove
                         <<" n_loc_tmp= "<<n_loc_tmp
                         <<std::endl;
                for(auto i=0; i<n_remove; i++){
                    std::cerr<<"i= "<<i<<" id_remove[i]= "<<id_remove[i]
                             <<" system[id_remove[i]].pos_car= "<<system[id_remove[i]].pos_car
                             <<" system[n_loc_tmp-1-i].pos_car= "<<system[n_loc_tmp-1-i].pos_car
                             <<std::endl;
                }
            }
            */
	    for(auto i=0; i<n_remove; i++){
                const auto id_tmp = id_remove[i];
		system[id_tmp].changeMass(system[id_tmp].mass * sat_mass_ratio);
		sat_system_loc.setSatellite(system[id_tmp]);
		id_remove[i] = id_tmp;
	    }
	    system.removeParticle(id_remove, n_remove);

            /*
            if(my_rank==0){
                std::cerr<<"n_remove= "<<n_remove
                         <<" n_loc_tmp= "<<n_loc_tmp
                         <<std::endl;
                for(auto i=0; i<n_remove; i++){
                    std::cerr<<"i= "<<i<<" id_remove[i]= "<<id_remove[i]
                             <<" sat_system_loc[i].pos= "<<sat_system_loc[i].pos_car
                             <<" system[id_remove[i]].pos= "<<system[id_remove[i]].pos_car
                             <<std::endl;
                }
            }
            exit(1);
            */
            
	}
	sat_system_glb.allgatherSystem(sat_system_loc);
	sat_system_glb.setId();
    }
    if(sat_system_glb.getNumberOfSat() > 0){
	std::cerr<<"sat_system_glb[0].mass= "<<sat_system_glb[0].mass
		 <<" pos_car= "<<sat_system_glb[0].pos_car
		 <<" pos_cyl= "<<sat_system_glb[0].pos_cyl
		 <<" vel= "<<sat_system_glb[0].vel
		 <<std::endl;
    }
    RemoveParticlesAroundSatellite(system, sat_system_glb); // remove overlapped particles
    //exit(1);
    if(my_rank==1){
	std::cerr<<"sat_system_glb.getNumberOfSat()= "<<sat_system_glb.getNumberOfSat()<<std::endl;
	for(auto i=0; i<sat_system_glb.getNumberOfSat(); i++){
	    std::cerr<<"sat_system_glb[i].id= "<<sat_system_glb[i].id<<" pos= "<<sat_system_glb[i].pos_car<<std::endl;
	}
    }
    n_loc = system.getNumberOfParticleLocal();
    n_glb = system.getNumberOfParticleGlobal();    
    if(my_rank==0){
	std::cerr<<"my_rank= "<<my_rank<<" n_loc= "<<n_loc<<" n_glb= "<<n_glb<<std::endl;
	system[0].dump(unit_man, fout_log);
	unit_man.dump(fout_log);
    }
    dinfo.decomposeDomainAll(system);
    if(my_rank==0){
        std::cerr<<"dinfo.getPosRootDomain()= "<<dinfo.getPosRootDomain()<<std::endl;
        for(auto i=0; i<n_proc; i++){
            std::cerr<<"i= "<<i<<" dinfo.getPosDomain(i)= "<<dinfo.getPosDomain(i)<<std::endl;
        }
    }
    system.adjustPositionIntoRootDomain(dinfo);
    n_loc = system.getNumberOfParticleLocal();
    system.exchangeParticle(dinfo);
    n_loc = system.getNumberOfParticleLocal();
    SetID(system, sat_system_glb.getNumberOfSat());

    PS::F64 mass_loc = 0.0;
    for(auto i=0; i<system.getNumberOfParticleLocal(); i++){
      mass_loc += system[i].mass;
    }
    PS::F64 mass_glb = PS::Comm::getSum(mass_loc);
    if(my_rank==0){
      std::cerr<<"mass_glb= "<<mass_glb<<std::endl;
    }
    
    Tree_t tree;
    tree.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    tree.setExchangeLETMode(ex_let_mode);
    tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
    //tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::MAKE_LIST);

    sat_system_glb.CalcForce(system);
    AddSatToSys(system, sat_system_glb);

    CalcForceFromPlanet(system, PLANET);

    if(my_rank != 0){
	RemoveSatFromSys(system, sat_system_glb);
    }
    eng_init.calc(system);
    if(my_rank != 0){
	AddSatToSys(system, sat_system_glb);
    }
    //exit(1);

    if(flag_snapshot){
	RemoveSatFromSys(system, sat_system_glb);
	snapshot_manager.write< PS::ParticleSystem<FP_t>, DiskInfo, Energy, SatelliteSystem>(system, param.getTimeSys(), disk_info, eng_init, sat_system_glb);
	AddSatToSys(system, sat_system_glb);
    }

    
    PS::S64 n_loop = 0;
    PS::S64 n_loop_prev = 0;
    PS::F64 eng_disp_loc = 0.0;
    //while(time_sys <= time_end){
    while(param.getTimeSys() <= time_end){

        PS::Comm::barrier();
	
        n_loc = system.getNumberOfParticleLocal();
	
        PS::F64 va0_loc = 0.0;
        for(int i=0; i<n_loc; i++){
            va0_loc += system[i].mass*system[i].acc_dash*system[i].vel;
        }
	
        Kick(system, 0.5*param.getDt());
        Drift(system, param.getDt());
	RemoveSatFromSys(system, sat_system_glb);

	param.integrateTime();
        n_loop++;
        if(n_loop % n_step_share == 0){
	    DEBUG_PRINT_RING();
            UpdateCyl(system);
	    if(flag_rot){
		RigidRotation(system, -param.getDt()*(n_loop-n_loop_prev) );
		n_loop_prev = n_loop;
	    }
	    DEBUG_PRINT_RING();
	    if(flag_mtm){
		Rotate(system, param.getDt()*n_step_share*0.5);
	    }
	    DEBUG_PRINT_RING();
            dinfo.decomposeDomainAll(system);

	    DEBUG_PRINT_RING();
            system.adjustPositionIntoRootDomain(dinfo);

	    DEBUG_PRINT_RING();
            system.exchangeParticle(dinfo);

	    DEBUG_PRINT_RING();

            n_loc = system.getNumberOfParticleLocal();
            n_glb = system.getNumberOfParticleGlobal();
        }
        tree.clearCounterAll();
	
	if(flag_mtm){	
	    if(n_loop % n_step_share == 0){
	        tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);

		DEBUG_PRINT_RING();

		Rotate(system, -param.getDt()*n_step_share*0.5);

		DEBUG_PRINT_RING();
		
	    }
	    tree.clearCounterAll();
	    tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, PS::REUSE_LIST);
	} else {
	  const auto int_mode = (n_loop % n_step_share == 0) ? PS::MAKE_LIST_FOR_REUSE : PS::REUSE_LIST;
	  //const auto int_mode = PS::MAKE_LIST_FOR_REUSE;
	    tree.calcForceAllAndWriteBack(CalcForceEp<EPI_t, EPJ_t, Force_t>(), CalcForceSp(), system, dinfo, true, int_mode);
	}
	
	DEBUG_PRINT_RING();

        PS::MemoryPool::checkEmpty();
        sat_system_glb.CalcForce(system);
	AddSatToSys(system, sat_system_glb);	
        CalcForceFromPlanet(system, PLANET);
	
        Kick(system, 0.5*param.getDt());

        n_loc = system.getNumberOfParticleLocal();
        PS::F64 va1_loc = 0.0;
        for(int i=0; i<n_loc; i++){
            va1_loc += system[i].mass*system[i].acc_dash*system[i].vel;
        }
        eng_disp_loc = (va0_loc+va1_loc)*param.getDt()*0.5;
        eng_disp_glb += PS::Comm::getSum(eng_disp_loc);
        Energy eng_now;
	if(my_rank != 0){
	    RemoveSatFromSys(system, sat_system_glb);
	}
        eng_now.calc(system);
	if(my_rank != 0){
	    AddSatToSys(system, sat_system_glb);
	}
	eng_now.setEngDispGlb(eng_disp_glb);
	
        auto n_walk = tree.getNumberOfWalkGlobal();
        auto n_int_ep_ep = tree.getNumberOfInteractionEPEPGlobal();
        auto n_int_ep_sp = tree.getNumberOfInteractionEPSPGlobal();
        if(PS::Comm::getRank()==0){
            eng_init.dump(std::cout);
            eng_now.dump(std::cout);
            std::cout<<"eng_now.disp= "<<eng_now.disp<<std::endl;
            std::cout<<"param.getTimeSys()= "<<param.getTimeSys()<<" n_loop= "<<n_loop<<" (eng_now.tot-eng_init.tot)/eng_init.tot= "<<(eng_now.tot-eng_init.tot)/eng_init.tot
                     <<" (eng_now.tot-eng_init.tot-eng_now.disp)/eng_init.tot= "<<(eng_now.tot-eng_init.tot-eng_now.disp)/eng_init.tot
		//<<" (eng_now.tot-eng_init.tot+eng_now.disp)/eng_init.tot= "<<(eng_now.tot-eng_init.tot+eng_now.disp)/eng_init.tot
                     <<std::endl;
            std::cout<<"n_int_ep_ep= "<<n_int_ep_ep
                     <<" n_int_ep_sp= "<<n_int_ep_sp
                     <<" <n_epi>= "<<((PS::F64)n_glb) / n_walk
                     <<" <n_epj>= "<<((PS::F64)n_int_ep_ep) / n_glb
                     <<" <n_spj>= "<<((PS::F64)n_int_ep_sp) / n_glb
                     <<std::endl;


            eng_init.dump(fout_log);
            eng_now.dump(fout_log);
            fout_log<<"eng_now.disp= "<<eng_now.disp<<std::endl;
            fout_log<<"param.getTimeSys()= "<<param.getTimeSys()<<" n_loop= "<<n_loop<<" (eng_now.tot-eng_init.tot)/eng_init.tot= "<<(eng_now.tot-eng_init.tot)/eng_init.tot
		    <<" (eng_now.tot-eng_init.tot-eng_now.disp)/eng_init.tot= "<<(eng_now.tot-eng_init.tot-eng_now.disp)/eng_init.tot
		//<<" (eng_now.tot-eng_init.tot+eng_now.disp)/eng_init.tot= "<<(eng_now.tot-eng_init.tot+eng_now.disp)/eng_init.tot
		    <<std::endl;
            fout_log<<"n_int_ep_ep= "<<n_int_ep_ep
		    <<" n_int_ep_sp= "<<n_int_ep_sp
		    <<" <n_epi>= "<<((PS::F64)n_glb) / n_walk
		    <<" <n_epj>= "<<((PS::F64)n_int_ep_ep) / n_glb
		    <<" <n_spj>= "<<((PS::F64)n_int_ep_sp) / n_glb
		    <<std::endl;
	    
        }
        tree.clearCounterAll();
	/*
        PS::S64 size_used_loc = tree.getUsedMemorySize();
        PS::S64 size_used_glb = PS::Comm::getMaxValue(size_used_loc);
        if(size_used_loc == size_used_glb){
            std::cerr<<"my_rank= "<<my_rank<<" tree.getUsedMemorySize()= "<<tree.getUsedMemorySize()<<std::endl;
            std::cerr<<" tree.dumpMemSizeUsed():"<<std::endl;
            tree.dumpMemSizeUsed(std::cerr);
        }
	*/
	if(flag_snapshot){
	    RemoveSatFromSys(system, sat_system_glb);
	    snapshot_manager.write< PS::ParticleSystem<FP_t>, DiskInfo, Energy, SatelliteSystem>(system, param.getTimeSys(), disk_info, eng_now, sat_system_glb);
	    AddSatToSys(system, sat_system_glb);
	}
    }
    
    PS::Finalize();
    return 0;
}
