#define SANITY_CHECK_REALLOCATABLE_ARRAY
#define VERBOSE

#include "header.h"

template <typename real> real plus(const real a){
	return (a > 0) ? a : 0;
}

template <typename real> real cube(const real a){
	return a * a * a;
}

const PS::F64 PI = atan(1.0) * 4.0;

struct kernel_t{
	const PS::F64 width;
	kernel_t(): width(2.0){
	}
	//W
	PS::F64 W(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = width * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
		r_value = plus(cube(1.0 - s)) - 4.0 * plus(cube(0.5 - s));
		//#if N_DIM == 1
		//r_value *= 8.0 / 3.0 / H;
		//#elif N_DIM == 2
		//r_value *= 80.0 / (7.0 * PI) / (H * H);
		//#elif N_DIM == 3
		r_value *= 16.0 / PI / (H * H * H);
		//#endif
		return r_value;
	}
	//gradW
	PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = width * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
		r_value = - 3.0 * plus(1.0 - s) * plus(1.0 - s) + 12.0 * plus(0.5 - s) * plus(0.5 - s);
		//#if N_DIM == 1
		//r_value *= 8.0 / 3.0 / H;
		//#elif N_DIM == 2
		//r_value *= 80.0 / (7.0 * PI) / (H * H);
		//#elif N_DIM == 3
		r_value *= 16.0 / PI / (H * H * H);
		//#endif
		return dr * r_value / (sqrt(dr * dr) * H  + 1.0e-6 * h);
	}
};

class CalcDensity{
	kernel_t kernel;
	public:
	void operator () (const EPI::Dens* const ep_i, const PS::S32 Nip, const EPJ::Dens* const ep_j, const PS::S32 Njp, RESULT::Dens* const dens){
		for(PS::S32 i = 0 ; i < Nip ; ++ i){
			dens[i].clear();
			for(PS::S32 j = 0 ; j < Njp ; ++ j){
				const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
				dens[i].dens += ep_j[j].mass * kernel.W(dr, ep_i[i].smth);
			}
			#warning TEMPORARY
			dens[i].smth = 1.2 * pow(ep_i[i].mass / dens[i].dens, 1.0/3.0);
		}
	}
};

class CalcDerivative{
	kernel_t kernel;
	public:
	void operator () (const EPI::Drvt* const ep_i, const PS::S32 Nip, const EPJ::Drvt* const ep_j, const PS::S32 Njp, RESULT::Drvt* const drvt){
		for(PS::S32 i = 0; i < Nip ; ++ i){
			drvt[i].clear();
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
				const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
				drvt[i].div_v += - ep_j[j].mass * dv * kernel.gradW(dr, ep_i[i].smth);
				drvt[i].rot_v += - ep_j[j].mass * dv ^ kernel.gradW(dr, ep_i[i].smth);
			}
			drvt[i].div_v /= ep_i[i].dens;
			drvt[i].rot_v /= ep_i[i].dens;
		}
	}
};

class CalcHydroForce{
	const kernel_t kernel;
	public:
	void operator () (const EPI::Hydro* const ep_i, const PS::S32 Nip, const EPJ::Hydro* const ep_j, const PS::S32 Njp, RESULT::Hydro* const hydro){
		const PS::F64 C_CFL = 0.3;
		for(PS::S32 i = 0; i < Nip ; ++ i){
			hydro[i].clear();
			PS::F64 v_sig_max = 0.0;
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
				const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
				const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
				const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;
				v_sig_max = std::max(v_sig_max, v_sig);
				const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens + ep_j[j].dens));
				const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ep_i[i].smth) + kernel.gradW(dr, ep_j[j].smth));
				hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) * gradW;
				hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + AV) * dv * gradW;
			}
			hydro[i].dt = C_CFL * 2.0 * ep_i[i].smth / v_sig_max;
		}
	}
};

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system){
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		#warning TEMPORARY
		sph_system[i].smth = 1.2 * pow(sph_system[i].mass / sph_system[i].dens, 1.0/3.0);
		sph_system[i].setPressure();
	}
}

template <class FDPSTree, class Functor> void PlantTreeByFDPS(Functor Force, PS::ParticleSystem<RealPtcl>& sph_system, PS::DomainInfo& dinfo){
	FDPSTree tree;
	std::cout << PS::Comm::getRank() << " = " << sph_system.getNumberOfParticleGlobal() << std::endl;
	tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
	//tree.initializeLocalTree(sph_system.getHalfLength());
	printf("setParticleLT\n");
	tree.setParticleLocalTree(sph_system);
	
	printf("setRC\n");
	tree.setRootCell(dinfo);
	printf("mSLTO\n");
	tree.mortonSortLocalTreeOnly();
	#ifdef VERBOSE
	tree.checkMortonSortLocalTreeOnly();
	#endif

	printf("linkCellLocalTreeOnly\n");
	tree.linkCellLocalTreeOnly();
	#ifdef VERBOSE
	tree.checkMakeLocalTree();
	#endif

	printf("calcMomentLocalTreeOnly\n");
	tree.calcMomentLocalTreeOnly();
	#ifdef VERBOSE
	tree.checkCalcMomentLocalTree();
	#endif
	printf("exchangeLocalEssentialTree\n");
	tree.exchangeLocalEssentialTree(dinfo);
	#ifdef VERBOSE
	tree.checkExchangeLocalEssentialTree(dinfo);
	#endif

	printf("exchangeLocalEssentialGlobalTree\n");
	tree.setLocalEssentialTreeToGlobalTree();

	printf("mortonSortGlobalTreeOnly\n");
	tree.mortonSortGlobalTreeOnly();
	#ifdef VERBOSE
	tree.checkMortonSortGlobalTreeOnly();
	#endif

	printf("linkCellGlobalTreeOnly\n");
	tree.linkCellGlobalTreeOnly();
	#ifdef VERBOSE
	tree.checkMakeGlobalTree();
	#endif

	printf("calcMomentGlobalTreeOnly\n");
	tree.calcMomentGlobalTreeOnly();
	#ifdef VERBOSE
	tree.checkCalcMomentGlobalTree();
	#endif

	printf("makeIPGroup\n");
	tree.makeIPGroup();
	#ifdef VERBOSE
	tree.checkMakeIPGroup();
	#endif

	printf("calcForceAndWriteBack\n");
	tree.calcForceAndWriteBack(Force, sph_system);
}

void TimeUpdate(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt){
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].pos += dt * sph_system[i].vel;
		sph_system[i].vel += dt * sph_system[i].hydro.acc;
		sph_system[i].eng += dt * sph_system[i].eng_dot;
	}
}

void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt){
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].vel_half = sph_system[i].vel + 0.5 * dt * sph_system[i].hydro.acc;
		sph_system[i].eng_half = sph_system[i].eng + 0.5 * dt * sph_system[i].eng_dot;
	}
}

void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt){
	//time becomes t + dt;
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].pos += dt * sph_system[i].vel_half;
	}
}

void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt){
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].vel += dt * sph_system[i].hydro.acc;
		sph_system[i].eng += dt * sph_system[i].eng_dot;
	}
}

void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt){
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].vel = sph_system[i].vel_half + 0.5 * dt * sph_system[i].hydro.acc;
		sph_system[i].eng = sph_system[i].eng_half + 0.5 * dt * sph_system[i].eng_dot;
	}
}

int main(int argc, char* argv[]){
	PS::Initialize(argc, argv);
	DisplayInfo();
	//
	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	PS::F64 dt, end_time;
	boundary box;
	//
	//Setup Initial
	if(1){//WRITING
		SetupIC(sph_system, &end_time, &box);
		Initialize(sph_system);
		#if 1
		if(0){//ASCII WRITING
			{
				FileHeader header;
				header.Nbody = sph_system.getNumberOfParticleGlobal();
				sph_system.writeParticleAscii("initial_with_header.txt", header);
			}
			/*
			{
				FileHeader header;
				header.Nbody = sph_system.getNumberOfParticleLocal();
				sph_system.writeParticleAscii("initial_with_header", "%s_%04d_%04d.txt", header);
			}
			{
				sph_system.writeParticleAscii("initial.txt");
			}
			{
				sph_system.writeParticleAscii("initial", "%s_%04d_%04d.txt");
			}
			*/
		}else{//BINARY WRITING
		}
		#endif
	}else{//READING
		if(1){//ASCII READING
			/*
			{//multiple file w/ header
				FileHeader header;
				sph_system.readParticleAscii("initial_with_header", "%s_%04d_%04d.txt", header);
				sph_system.writeParticleAscii("initial2_with_header", "%s_%04d_%04d.txt", header);
			}
			*/
			{//single file w/ header
				FileHeader header;
				sph_system.readParticleAscii("initial_with_header.txt", header);
				//sph_system.writeParticleAscii("initial2_with_header.txt", header);
				end_time = 0.2;
				box.x = 1.0;
				box.y = box.z = box.x / 1.0;
			}
			/*
			{//multiple file w/o header
				sph_system.readParticleAscii("initial", "%s_%04d_%04d.txt");
				sph_system.writeParticleAscii("initial2", "%s_%04d_%04d.txt");
			}
			{//single file w/o header
				sph_system.readParticleAscii("initial.txt");
				sph_system.writeParticleAscii("initial2.txt");
			}
			*/
		}else{//BINARY READING
		}
		Initialize(sph_system);
	}
	PS::DomainInfo dinfo;
	dinfo.initialize();
	dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
	dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), PS::F64vec(box.x, box.y, box.z));
	dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
	std::cerr << "# of ptcls in this process is..." << sph_system.getNumberOfParticleLocal() << std::flush << std::endl;
	PS::S32 cnt = 0;
	sph_system.adjustPositionIntoRootDomain(dinfo);
	for(PS::F64 time = 0 ; time < end_time ; time += dt, ++ cnt){
		std::cerr << "collectSampleParticle" << std::flush << std::endl;
		dinfo.collectSampleParticle(sph_system);
		std::cout << "decomposeDomain" << std::endl;
		dinfo.decomposeDomain();
		std::cout << "exchangeParticle" << std::endl;
		sph_system.exchangeParticle(dinfo);
		sph_system.adjustPositionIntoRootDomain(dinfo);
		//
		std::cout << "Density tree" << std::endl;
		#if 0
		PlantTreeByFDPS<PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather>(CalcDensity()   , sph_system, dinfo);
		#else
		PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather dens_tree;
		dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
		#endif

		std::cout << "Derivative tree" << std::endl;
		PlantTreeByFDPS<PS::TreeForForceShort<RESULT::Drvt, EPI::Drvt, EPJ::Drvt>::Gather>(CalcDerivative(), sph_system, dinfo);
		//set pres.
		{
			for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
				sph_system[i].setPressure();
			}
		}
		std::cout << "Force tree" << std::endl;
		PlantTreeByFDPS<PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry>(CalcHydroForce(), sph_system, dinfo);
		{
			dt = 1.0e+30;
			for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
				dt = std::min(dt, sph_system[i].dt);
			}
			dt = PS::Comm::getMinValue(dt);
		}
		TimeUpdate(sph_system, dt);
		sph_system.adjustPositionIntoRootDomain(dinfo);
		{
			FileHeader header;
			header.time = time;
			header.Nbody = sph_system.getNumberOfParticleGlobal();
			char filename[256];
			sprintf(filename, "result/%04d.txt", cnt);
			sph_system.writeParticleAscii(filename, header);
		}
		std::cout << "time = " << time << std::endl;
	}

	PS::Finalize();
	return 0;
}

