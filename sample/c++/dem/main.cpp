#include <particle_simulator.hpp>
#include <random>

#include "force.h"
#include "class.h"
#include "surface.h"
#include "init/InfiniteCylinder2.h"
#include "integral.h"
#include "io.h"

int main(int argc, char *argv[]){
	PS::Initialize(argc, argv);
	PS::ParticleSystem<FP> ptcl;
	ptcl.initialize();
	PS::TreeForForceShort<Force, FP, FP>::Symmetry tree_coll;

	PS::DomainInfo dinfo;
	dinfo.initialize();
	system_t sysinfo;
	FileIO<FP> io;

	if(argc == 1){
		Problem<FP>::SetIC(ptcl, dinfo, sysinfo);
		tree_coll.initialize(ptcl.getNumberOfParticleLocal(), 0.5, 1, 1);
	}else{
		io.Restore(ptcl, sysinfo);
		tree_coll.initialize(ptcl.getNumberOfParticleLocal(), 0.5, 1, 1);
		goto savepoint;
	}

	dinfo.decomposeDomainAll(ptcl);
	ptcl.exchangeParticle(dinfo);
	tree_coll.calcForceAllAndWriteBack(Force(), ptcl, dinfo);
	Problem<FP>::externalForce(ptcl, sysinfo);

	for(sysinfo.time = 0, sysinfo.step = 0 ; sysinfo.time < sysinfo.end_time ; sysinfo.time += sysinfo.dt, ++ sysinfo.step){
		dinfo.decomposeDomainAll(ptcl);
		ptcl.exchangeParticle(dinfo);
		sysinfo.dt = GetGlobalTimestep<FP>(ptcl);
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			ptcl[i].kick(sysinfo.dt);
			ptcl[i].drift(sysinfo.dt);
			ptcl[i].clear();
		}
		tree_coll.calcForceAllAndWriteBack(Force(), ptcl, dinfo);
		Problem<FP>::externalForce(ptcl, sysinfo);
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			ptcl[i].kick2(sysinfo.dt);
		}
		Problem<FP>::postTimestep(ptcl, sysinfo);
		ptcl.adjustPositionIntoRootDomain(dinfo);
		if(PS::Comm::getRank() == 0 && sysinfo.step % 100 == 0) std::cerr << "time = " << sysinfo.time << " (dt = " << sysinfo.dt << ")" << std::endl;
		io.OutputFileWithTimeInterval(ptcl, sysinfo);
		if(sysinfo.step % 10000 == 0) io.Create(ptcl, sysinfo);
		savepoint:;
	}

	
	PS::Finalize();
	return 0;
}

