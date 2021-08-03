#pragma once

template <class ThisPtcl> class Problem{
	static constexpr double rad = 0.125;
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		sysinfo.end_time = 1.0;
		ptcl.setNumberOfParticleLocal(0);
		if(PS::Comm::getRank() != 0) return;
		ptcl.setNumberOfParticleLocal(3);
		
		ThisPtcl ith;

		ith.pos.x = 0.5;
		ith.pos.y = 0.0;
		ith.pos.z = 2.0 * rad;
		ith.rad   = rad;
		ith.mat   = Test;
		ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
		ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
		ith.id    = 0;
		ptcl[0]   = ith;

		ith.pos.x = 1.0 + rad / 2.0;
		ith.pos.y = 0.0;
		ith.pos.z = 2.0 * rad;
		ith.rad   = rad;
		ith.mat   = Test;
		ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
		ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
		ith.id    = 1;
		ptcl[1]   = ith;

		ith.pos.x = 0.0 + 0.09;
		ith.pos.y = 0.0;
		ith.pos.z = 1.1 * rad;
		ith.rad   = rad * 0.25;
		ith.mat   = Test;
		ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
		ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
		ith.id    = 1;
		ptcl[2]   = ith;
		
	}
	static void externalForce(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		Disc<ThisPtcl> disc(Test, PS::F64vec(0.0, 0.0, 1.0), PS::F64vec(0.0, 0.0, 0.0), 0.1, 1.0);
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			//Earth's gravity
			ptcl[i].acc.z -= 9.8;
			disc.setForce(ptcl[i]);
		}
	}
	static void postTimestep(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
	}
};

