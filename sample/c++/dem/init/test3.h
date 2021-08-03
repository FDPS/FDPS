#pragma once

template <class ThisPtcl> class Problem{
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		const PS::S64 N = 100;
		ptcl.setNumberOfParticleLocal(N);
		//dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
		//dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), PS::F64vec(1.0, 1.0, 1.0));
		sysinfo.end_time = 5.0;
		if(PS::Comm::getRank() != 0) return;
		for(int i = 0 ; i < N ; ++ i){
			ptcl[i].id = i;
			ptcl[i].pos.x = 0.9 * rand() / (double)(RAND_MAX) + 0.05;
			ptcl[i].pos.y = 0.9 * rand() / (double)(RAND_MAX) + 0.05;
			ptcl[i].pos.z = 3.0 * rand() / (double)(RAND_MAX) + 1.0;
			ptcl[i].vel.x = 0.0;
			ptcl[i].vel.y = 0.0;
			ptcl[i].vel.z = 0.0;
			if(ptcl[i].id % 5 == 0){
				ptcl[i].rad = 0.05;
				ptcl[i].mat = Basalt;
			}else{
				ptcl[i].rad = 0.02;
				ptcl[i].mat = Ice;
			}
			ptcl[i].mass = 4.0 * M_PI / 3.0 * pow(ptcl[i].rad, 3) * ptcl[i].mat.getDensity();
		}
	}
	static void externalForce(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		Plane floor(0, 0, 1, 0.0, Basalt);
		Plane wall1(1, 0, 0,  0.0, Basalt);
		Plane wall2(1, 0, 0, -1.0, Basalt);
		Plane wall3(0, 1, 0,  0.0, Basalt);
		Plane wall4(0, 1, 0, -1.0, Basalt);
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			//Earth's gravity
			ptcl[i].acc.z -= 9.8;
			ptcl[i].acc += floor.Force(ptcl[i]);
			ptcl[i].acc += wall1.Force(ptcl[i]);
			ptcl[i].acc += wall2.Force(ptcl[i]);
			ptcl[i].acc += wall3.Force(ptcl[i]);
			ptcl[i].acc += wall4.Force(ptcl[i]);
		}
	}
};

