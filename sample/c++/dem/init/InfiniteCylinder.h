#pragma once

template <class ThisPtcl> class Problem{
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		const PS::S64 N = 2;
		ptcl.setNumberOfParticleLocal(N);
		sysinfo.end_time = 10.0;
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
		dinfo.setPosRootDomain(PS::F64vec(-3.0, -3.0, -3.0), PS::F64vec(3.0, 3.0, 3.0));
		if(PS::Comm::getRank() != 0) return;
		int i = 0;
		ptcl[i].id = i;
		ptcl[i].pos.x = 0.0;
		ptcl[i].pos.y = 1.5;
		ptcl[i].pos.z = 0.0;
		ptcl[i].vel.x = -1.0;
		ptcl[i].vel.y = 0.0;
		ptcl[i].vel.z = 0.0;
		ptcl[i].avel.x = 0.0;
		ptcl[i].avel.y = 0.0;
		ptcl[i].avel.z = 0.0;
		ptcl[i].rad  = 0.5;
		ptcl[i].mat  = Test;
		ptcl[i].mass = 4.0 * M_PI / 3.0 * pow(ptcl[i].rad, 3) * ptcl[i].mat.getDensity();
		ptcl[i].iner = 0.4 * ptcl[i].mass * ptcl[i].rad * ptcl[i].rad;

		i = 1;
		ptcl[i].id = i;
		ptcl[i].pos.x = 1.5;
		ptcl[i].pos.y = 0.0;
		ptcl[i].pos.z = 0.0;
		ptcl[i].vel.x = -1.0;
		ptcl[i].vel.y = 0.0;
		ptcl[i].vel.z = 0.0;
		ptcl[i].avel.x = 0.0;
		ptcl[i].avel.y = 0.0;
		ptcl[i].avel.z = 0.0;
		ptcl[i].rad  = 0.5;
		ptcl[i].mat  = Test;
		ptcl[i].mass = 4.0 * M_PI / 3.0 * pow(ptcl[i].rad, 3) * ptcl[i].mat.getDensity();
		ptcl[i].iner = 0.4 * ptcl[i].mass * ptcl[i].rad * ptcl[i].rad;


	}
	static void externalForce(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		InfiniteCylinder<ThisPtcl> cylinder(Test, PS::F64vec(0.0, 1.0, 0.0), PS::F64vec(0.0, 0.0, 0.0), 1.0);
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			//Earth's gravity
			cylinder.setForce(ptcl[i]);
		}
	}
	static void postTimestep(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
	}
};

