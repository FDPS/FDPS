//#define SANITY_CHECK_REALLOCATABLE_ARRAY

#include <sys/stat.h>
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
	printf("Dom info\n");
	PS::DomainInfo dinfo;
	printf("init.\n");
	dinfo.initialize();
	printf("setBC\n");
	dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
	printf("setPRD\n");
	dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), PS::F64vec(box.x, box.y, box.z));
	printf("setD\n");
	dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
	printf("decomposeDomain\n");
	dinfo.decomposeDomain();
	printf("ex.P\n");
	sph_system.exchangeParticle(dinfo);
	printf("TreeFFS\n");
	{
		PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather dens_tree;
		printf("Tree initialize\n");
		dens_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
		printf("Tree calcForceAllAndWriteBack\n");
		dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
	}
	printf("TreeFFS\n");
	{
		PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
		hydr_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
		hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
	}
	printf("getDT\n");
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
            struct stat st;
            if(stat("result", &st) != 0) {
                PS::S32 rank = PS::Comm::getRank();
                PS::S32 ret_loc, ret;
                if(rank == 0)
                    ret_loc = mkdir("result", 0777);
                PS::Comm::broadcast(&ret_loc, ret);
                if(ret == 0) {
                    if(rank == 0)
                        fprintf(stderr, "Directory \"result\" is successfully made.\n");
                } else {
                    fprintf(stderr, "Directory \"result\" fails to be made.\n");
                    PS::Finalize();
                    exit(0);
                }
            }
			sprintf(filename, "result/%04d.dat", step);
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
		}
			CheckConservativeVariables(sph_system);
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
        }
	}

	PS::Finalize();
	return 0;
}

