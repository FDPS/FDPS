#pragma once

template <class ThisPtcl> class Problem{
	static constexpr int N = 16;
	static constexpr double rad = 0.25 / N;
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		sysinfo.end_time = 30.0;
		ptcl.setNumberOfParticleLocal(0);
		if(PS::Comm::getRank() != 0) return;
		ptcl.setNumberOfParticleLocal(1);
		
		ThisPtcl ith;
		ith.pos.x = (rand() / double(RAND_MAX)) * (1.0 - 2.0 * rad) + rad;
		ith.pos.y = (rand() / double(RAND_MAX)) * (1.0 - 2.0 * rad) + rad;
		ith.pos.z = 1.0;
		ith.rad   = rad;
		ith.mat   = Pyrex;
		ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
		ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
		ith.id    = 0;
		ptcl[0]   = ith;
		
	}
	static void externalForce(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		const PS::F64 duration = 10.0;
		if(sysinfo.time >= duration){
			InfinitePlane<ThisPtcl> floor2(Pyrex, PS::F64vec(0.0, 0.0, 1.0), PS::F64vec(0.0, 0.0, 0.0));
			InfinitePlane<ThisPtcl> wall1(Pyrex, PS::F64vec(1.0, 0.0, 0.0), PS::F64vec(0.0, 0.0, 0.0));
			InfinitePlane<ThisPtcl> wall2(Pyrex, PS::F64vec(0.0, 1.0, 0.0), PS::F64vec(0.0, 0.0, 0.0));
			InfinitePlane<ThisPtcl> wall3(Pyrex, PS::F64vec(-1.0, 0.0, 0.0), PS::F64vec(2.0, 1.0, 0.0));
			InfinitePlane<ThisPtcl> wall4(Pyrex, PS::F64vec(0.0, -1.0, 0.0), PS::F64vec(2.0, 1.0, 0.0));
			for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
				//Earth's gravity
				ptcl[i].acc.z -= 9.8;
				floor2.setForce(ptcl[i]);
				wall1.setForce(ptcl[i]);
				wall2.setForce(ptcl[i]);
				wall3.setForce(ptcl[i]);
				wall4.setForce(ptcl[i]);
			}
		}else{
			InfinitePlane<ThisPtcl> floor2(Pyrex, PS::F64vec(0.0, 0.0, 1.0), PS::F64vec(0.0, 0.0, 0.0));
			InfinitePlane<ThisPtcl> wall1(Pyrex, PS::F64vec(1.0, 0.0, 0.0), PS::F64vec(0.0, 0.0, 0.0));
			InfinitePlane<ThisPtcl> wall2(Pyrex, PS::F64vec(0.0, 1.0, 0.0), PS::F64vec(0.0, 0.0, 0.0));
			InfinitePlane<ThisPtcl> wall3(Pyrex, PS::F64vec(-1.0, 0.0, 0.0), PS::F64vec(1.0, 1.0, 0.0));
			InfinitePlane<ThisPtcl> wall4(Pyrex, PS::F64vec(0.0, -1.0, 0.0), PS::F64vec(1.0, 1.0, 0.0));
			for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
				//Earth's gravity
				ptcl[i].acc.z -= 9.8;
				floor2.setForce(ptcl[i]);
				wall1.setForce(ptcl[i]);
				wall2.setForce(ptcl[i]);
				wall3.setForce(ptcl[i]);
				wall4.setForce(ptcl[i]);
			}
		}
	}
	static void postTimestep(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		if(PS::Comm::getRank() != 0) return;
		//radnom number
		static std::mt19937 mt(0);
		//inflow span [sec]
		const PS::F64 inflow_span = 2.0 * sqrt(rad / 9.8);
		//inflow duration
		const PS::F64 duration = 10.0;

		if(sysinfo.time >= duration) return;
		static PS::F64 time = inflow_span;
		
		if(sysinfo.time >= time){
			const std::size_t inflow_number = N * N;
			const std::size_t tail = ptcl.getNumberOfParticleLocal() - 1;
			ptcl.setNumberOfParticleLocal(ptcl.getNumberOfParticleLocal() + inflow_number);
			//
			ThisPtcl ith;
			int cnt = 0;
			for(int Nx = 0 ; Nx < N ; ++ Nx){
				for(int Ny = 0 ; Ny < N ; ++ Ny){
					++ cnt;
					const PS::F64 min_x = double(Nx) / N + rad;
					const PS::F64 max_x = double(Nx + 1) / N - rad;
					const PS::F64 min_y = double(Ny) / N + rad;
					const PS::F64 max_y = double(Ny + 1) / N - rad;
					std::uniform_real_distribution<double> x(min_x, max_x);
					std::uniform_real_distribution<double> y(min_y, max_y);
					ith.pos.x = x(mt);
					ith.pos.y = y(mt);
					ith.pos.z = 1.0;
					ith.rad   = rad;
					ith.mat   = Pyrex;
					ith.mass  = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
					ith.iner  = 0.4 * ith.mass * ith.rad * ith.rad;
					ith.id    = tail + cnt;
					ptcl[tail + cnt] = ith;
				}
			}
			time += inflow_span;
		}
	}
};

