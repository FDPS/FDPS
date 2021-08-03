#pragma once

template <class ThisPtcl> class Problem{
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		const PS::S64 N = 1;
		ptcl.setNumberOfParticleLocal(N);
		sysinfo.end_time = 60.0;
		//dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
		//dinfo.setPosRootDomain(PS::F64vec(-1.0, -1.0, -1.0), PS::F64vec(1.0, 1.0, 1.0));
		if(PS::Comm::getRank() != 0) return;
		const int i = 0;
		ptcl[i].id = i;
		ptcl[i].rad  = 1.0/16.0;
		ptcl[i].pos.x = 0.0;
		ptcl[i].pos.y = 0.0;
		ptcl[i].pos.z = 0.0 + 1.5 *  ptcl[i].rad;
		ptcl[i].vel.x = 0.0;
		ptcl[i].vel.y = 0.0;
		ptcl[i].vel.z = 0.0;
		ptcl[i].avel.x = 0.0;
		ptcl[i].avel.y = 0.0;
		ptcl[i].avel.z = 0.0;
		ptcl[i].mat  = Test;
		ptcl[i].mass = 4.0 * M_PI / 3.0 * pow(ptcl[i].rad, 3) * ptcl[i].mat.getDensity();
		ptcl[i].iner = 0.4 * ptcl[i].mass * ptcl[i].rad * ptcl[i].rad;
	}
	static void externalForce(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		InfiniteCylinder<ThisPtcl> cylinder(Test, PS::F64vec(0.0, 0.0, 1.0), PS::F64vec(0.0, 0.0, 0.0), 1.0);
		Disc<ThisPtcl> disc1(Test, PS::F64vec(0.0, 0.0, 1.0), PS::F64vec(0.0, 0.0, 0.0), 0.0, 1.0);
		Disc<ThisPtcl> disc2(Test, PS::F64vec(0.0, 0.0, 1.0), PS::F64vec(0.0, 0.0, 0.0), 0.25, 1.0);
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			//Earth's gravity
			ptcl[i].acc.z -= 9.8;
			cylinder.setForce(ptcl[i]);
			if(sysinfo.time < 30){
				disc1.setForce(ptcl[i]);
			}else{
				disc2.setForce(ptcl[i]);
			}
		}
	}
	static void postTimestep(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			if(ptcl[i].pos.z < - 1.0){
				ptcl[i] = ptcl[ptcl.getNumberOfParticleLocal() - 1];
				ptcl.setNumberOfParticleLocal(ptcl.getNumberOfParticleLocal() - 1);
			}
		}
		if(PS::Comm::getRank() != 0) return;
		//radius
		const PS::F64 rad = 1.0/16.0;
		//inflow duration
		const PS::F64 duration = 30;
		//inflow speed
		const PS::F64 velocity = 2.0;
		//inflow span [sec]
		const PS::F64 inflow_span = (-velocity + sqrt(velocity * velocity + 2.0 * 9.8 * rad)) / 9.8 * 2.0;

		if(sysinfo.time >= duration) return;
		static PS::F64 time = inflow_span;
		
		if(sysinfo.time >= time){
			const std::size_t inflow_number = 9;
			const std::size_t tail = ptcl.getNumberOfParticleLocal() - 1;
			ptcl.setNumberOfParticleLocal(ptcl.getNumberOfParticleLocal() + inflow_number);
			//
			ThisPtcl ith;
			{
				ith.pos.x = 0.0;
				ith.pos.y = 0.0;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 1;
				ptcl[tail + 1] = ith;
			}
			{
				ith.pos.x = 0.5;
				ith.pos.y = 0.5;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 2;
				ptcl[tail + 2] = ith;
			}
			{
				ith.pos.x = -0.5;
				ith.pos.y = 0.5;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 3;
				ptcl[tail + 3] = ith;
			}
			{
				ith.pos.x = 0.5;
				ith.pos.y = -0.5;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 4;
				ptcl[tail + 4] = ith;
			}
			{
				ith.pos.x = -0.5;
				ith.pos.y = -0.5;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 5;
				ptcl[tail + 5] = ith;
			}
			{
				ith.pos.x = 0.25;
				ith.pos.y = 0.25;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 6;
				ptcl[tail + 6] = ith;
			}
			{
				ith.pos.x = -0.25;
				ith.pos.y = 0.25;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 7;
				ptcl[tail + 7] = ith;
			}
			{
				ith.pos.x = 0.25;
				ith.pos.y = -0.25;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 8;
				ptcl[tail + 8] = ith;
			}
			{
				ith.pos.x = -0.25;
				ith.pos.y = -0.25;
				ith.pos.z = 4.0;
				ith.vel.x = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.y = 0.1 * (2.0 * (double)rand() / RAND_MAX - 1.0);
				ith.vel.z = - velocity;
				ith.rad   = rad;
				ith.mat   = Test;
				ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
				ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
				ith.id    = tail + 9;
				ptcl[tail + 9] = ith;
			}

			time += inflow_span;
		}
	}
};

