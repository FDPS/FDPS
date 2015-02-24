#include "header.h"

void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time, boundary *box){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	const PS::F64 dx = 1.0 / 32.0;
	box->x = 1.0;
	box->y = box->z = box->x / 1.0;
	PS::S32 i = 0;
	for(PS::F64 x = 0 ; x < box->x * 0.5 ; x += dx){
		for(PS::F64 y = 0 ; y < box->y ; y += dx){
			for(PS::F64 z = 0 ; z < box->z ; z += dx){
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.dens = 1.0;
				ith.mass = ith.dens;
				ith.eng  = 1.0;
				ith.id   = i++;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::F64 x = box->x * 0.5 ; x < box->x * 1.0 ; x += dx){
		for(PS::F64 y = 0 ; y < box->y ; y += dx){
			for(PS::F64 z = 0 ; z < box->z ; z += dx){
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.dens = 0.5;
				ith.mass = ith.dens;
				ith.eng  = 1.0;
				ith.id   = i++;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		ptcl[i].mass /= (PS::F64)(ptcl.size());
	}
	std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
	/////////
	//scatter ptcls↲
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
	*end_time = 0.2;
	//Fin.
	std::cout << "setup..." << std::endl;
}

/*
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time, boundary *box){
	/////////
	//place ptcls
	/////////
	if(PS::Comm::getRank() == 0){
		std::vector<RealPtcl> ptcl;
		const PS::F64 dx = 1.0 / 16.0;
		box->x = 1.0;
		box->y = box->z = box->x / 1.0;
		PS::S32 i = 0;
		for(PS::F64 x = 0 ; x < box->x ; x += dx){
			for(PS::F64 y = 0 ; y < box->y ; y += dx){
				for(PS::F64 z = 0 ; z < box->z ; z += dx){
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
			ptcl[i].mass = 1.0 * (box->x * box->y * box->z) / (PS::F64)(ptcl.size());
		}
		std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
		const PS::U32 numPtclLocal = ptcl.size();
		sph_system.setNumberOfParticleLocal(numPtclLocal);
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			sph_system[i] = ptcl[i];
		}
	}else{
		sph_system.setNumberOfParticleLocal(0);
	}
	*end_time = 0.2;
	//Fin.
	std::cout << "setup..." << std::endl;
}
*/