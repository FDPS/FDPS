#define SANITY_CHECK_REALLOCATABLE_ARRAY
//#define VERBOSE
#include "header.h"

int main(int argc, char* argv[]){
	PS::Initialize(argc, argv);
	DisplayInfo();
	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	PS::F64 dt, end_time;
	boundary box;
	//
	//Setup Initial
	SetupIC(sph_system, &end_time, &box);
	Initialize(sph_system);
	//Dom. info
	PS::DomainInfo dinfo;
	dinfo.initialize();
	dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
	dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), PS::F64vec(box.x, box.y, box.z));
	dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
	dinfo.decomposeDomain();
	sph_system.exchangeParticle(dinfo);
	{
		PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather dens_tree;
		dens_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
		dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
	}
	{
		PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
		hydr_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
		hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
	}
	dt = getTimeStepGlobal(sph_system);

	PS::S32 step = 0;
	for(PS::F64 time = 0 ; time < end_time ; time += dt, ++ step){
		InitialKick(sph_system, dt);
		FullDrift(sph_system, dt);
		sph_system.adjustPositionIntoRootDomain(dinfo);
		Predict(sph_system, dt);
		dinfo.decomposeDomain();
		sph_system.exchangeParticle(dinfo);
		{
			PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather dens_tree;
			dens_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
			dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
		}
		for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure();
		}
		{
			PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
			hydr_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
			hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
		}
		dt = getTimeStepGlobal(sph_system);

		FinalKick(sph_system, dt);

		if(step % PARAM::OUTPUT_INTERVAL == 0){
			FileHeader header;
			header.time = time;
			header.Nbody = sph_system.getNumberOfParticleGlobal();
			char filename[256];
			sprintf(filename, "result/%04d.txt", step);
			sph_system.writeParticleAscii(filename, header);
			if(PS::Comm::getRank() == 0){
				std::cout << "//================================" << std::endl;
				std::cout << "output " << filename << "." << std::endl;
				std::cout << "//================================" << std::endl;
			}
		}
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << "time = " << time << std::endl;
			std::cout << "step = " << step << std::endl;
			std::cout << "//================================" << std::endl;
		}
	}

	PS::Finalize();
	return 0;
}

