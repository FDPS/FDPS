#include "header.h"

void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::U32* Nptcl){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	const PS::F64 dx = 1.0 / 16.0;
	const PS::F64 rad = 1.0;
	PS::S32 i = 0;
	for(PS::F64 x = 0 ; x < rad ; x += dx){
		for(PS::F64 y = 0 ; y < rad ; y += dx){
			for(PS::F64 z = 0 ; z < rad ; z += dx){
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.dens = 1.0;
				ith.eng  = 1.0;
				ith.id   = i++;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		ptcl[i].mass = 1.0 / (PS::F64)(ptcl.size());
	}
	*Nptcl = ptcl.size();
	std::cout << "# of ptcls is... " << *Nptcl << std::endl; 
	/////////
	//scatter ptcls
	/////////
	assert(ptcl.size() % PS::Comm::getNumberOfProc() == 0);
	const PS::U32 numPtclLocal = ptcl.size() / PS::Comm::getNumberOfProc();
	sph_system.createParticle(numPtclLocal*8+10000);
	const PS::U32 i_head = numPtclLocal * PS::Comm::getRank();
	const PS::U32 i_tail = numPtclLocal * (PS::Comm::getRank() + 1);
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		if(i_head <= i && i < i_tail){
			const PS::U32 ii = i - numPtclLocal * PS::Comm::getRank();
			sph_system[ii] = ptcl[i];
		}
	}
	sph_system.setNumberOfParticleLocal(numPtclLocal);
	//Fin.
	std::cout << "setup..." << std::endl;
}
