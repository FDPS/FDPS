#include "header.h"

void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time, boundary *box){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	const PS::F64 dx = 1.0 / 128.0;
	box->x = 1.0;
	box->y = box->z = box->x / 8.0;
	PS::S32 i = 0;
	for(PS::F64 x = 0 ; x < box->x * 0.5 ; x += dx){
		for(PS::F64 y = 0 ; y < box->y ; y += dx){
			for(PS::F64 z = 0 ; z < box->z ; z += dx){
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.dens = 1.0;
				ith.mass = 0.75;
				ith.eng  = 2.5;
				ith.id   = i++;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::F64 x = box->x * 0.5 ; x < box->x * 1.0 ; x += dx * 2.0){
		for(PS::F64 y = 0 ; y < box->y ; y += dx){
			for(PS::F64 z = 0 ; z < box->z ; z += dx){
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.dens = 0.5;
				ith.mass = 0.75;
				ith.eng  = 2.5;
				ith.id   = i++;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		ptcl[i].mass = ptcl[i].mass * box->x * box->y * box->z / (PS::F64)(ptcl.size());
	}
	std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
	/////////
	//scatter ptcls
	/////////↲
	assert(ptcl.size() % PS::Comm::getNumberOfProc() == 0);
	const PS::S32 numPtclLocal = ptcl.size() / PS::Comm::getNumberOfProc();
	sph_system.setNumberOfParticleLocal(numPtclLocal);
	const PS::U32 i_head = numPtclLocal * PS::Comm::getRank();
	const PS::U32 i_tail = numPtclLocal * (PS::Comm::getRank() + 1);
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		if(i_head <= i && i < i_tail){
			const PS::U32 ii = i - numPtclLocal * PS::Comm::getRank();
			sph_system[ii] = ptcl[i];
		}
	}
	/////////
	*end_time = 0.1;
	//Fin.
	std::cout << "setup..." << std::endl;
}

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system){
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0/3.0);
		sph_system[i].setPressure();
	}
}
