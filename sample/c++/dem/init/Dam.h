#pragma once

template <class ThisPtcl> class Problem{
	public:
	static void SetIC(PS::ParticleSystem<ThisPtcl>& ptcl, PS::DomainInfo& dinfo, system_t& sysinfo){
		std::vector<ThisPtcl> tmp;
		const double dx = 1.0/50.0;
		sysinfo.end_time = 5.0;
		if(PS::Comm::getRank() != 0) return;
		int cnt = 0;
		for(double x = 0.5 * dx ; x < 0.2 ; x += dx){
			for(double y = 0.8 + 0.5 * dx ; y < 1.0 ; y += dx){
				for(double z = 0.5 * dx ; z < 0.2 ; z += dx){
					ThisPtcl ith;
					ith.id = cnt;
					ith.pos.x = x;
					ith.pos.y = y;
					ith.pos.z = z;
					ith.vel.x = 0.01 * sqrt(2 * 9.8) * (1.0 - 2.0 * rand()/(double)(RAND_MAX));
					ith.vel.y = 0.01 * sqrt(2 * 9.8) * (1.0 - 2.0 * rand()/(double)(RAND_MAX));
					ith.vel.z = 0.01 * sqrt(2 * 9.8) * (1.0 - 2.0 * rand()/(double)(RAND_MAX));
					ith.rad = 0.9 * dx / 2.0;
					ith.mat = Test;
					ith.mass = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
					ith.iner = 0.4 * ith.mass * ith.rad * ith.rad;
					tmp.push_back(ith);
					++ cnt;
				}
			}
		}
		#if 0
		for(double x = - dx ; x <= 1.0 + dx ; x += dx){
			for(double y = - dx ; y <= 1.0 + dx ; y += dx){
				for(double z = - dx ; z <= 1.0 + dx ; z += dx){
					if(0 <= x || x <= 1.0){
						continue;
					}
					ThisPtcl ith;
					ith.id = cnt;
					ith.pos.x = x;
					ith.pos.y = y;
					ith.pos.z = z;
					ith.vel.x = 0.01 * sqrt(2 * 9.8) * (1.0 - 2.0 * rand()/(double)(RAND_MAX));
					ith.vel.y = 0.01 * sqrt(2 * 9.8) * (1.0 - 2.0 * rand()/(double)(RAND_MAX));
					ith.vel.z = 0.01 * sqrt(2 * 9.8) * (1.0 - 2.0 * rand()/(double)(RAND_MAX));
					ith.rad = 0.9 * dx / 2.0;
					ith.mat = Test;
					ith.mass = 4.0 * M_PI / 3.0 * pow(ith.rad, 3) * ith.mat.getDensity();
					ith.iner = 0.4 * ith.mass * ith.rad * ith.rad;
					tmp.push_back(ith);
					++ cnt;
				}
			}
		}
		#endif

		ptcl.setNumberOfParticleLocal(tmp.size());
		for(int i = 0 ; i < tmp.size() ; ++ i){
			ptcl[i] = tmp[i];
		}
	}
	static void externalForce(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		InfinitePlane<ThisPtcl> floor1(Test, PS::F64vec(0.0, 0.0, 1.0), PS::F64vec(0.0, 0.0, 0.0));
		InfinitePlane<ThisPtcl> wall1(Test, PS::F64vec(1.0, 0.0, 0.0), PS::F64vec(0.0, 0.0, 0.0));
		InfinitePlane<ThisPtcl> wall2(Test, PS::F64vec(0.0, 1.0, 0.0), PS::F64vec(0.0, 0.0, 0.0));
		InfinitePlane<ThisPtcl> wall3(Test, PS::F64vec(1.0, 0.0, 0.0), PS::F64vec(1.0, 1.0, 0.0));
		InfinitePlane<ThisPtcl> wall4(Test, PS::F64vec(0.0, 1.0, 0.0), PS::F64vec(1.0, 1.0, 0.0));
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			//Earth's gravity
			ptcl[i].acc.z -= 9.8;
			floor1.setForce(ptcl[i]);
			wall1.setForce(ptcl[i]);
			wall2.setForce(ptcl[i]);
			wall3.setForce(ptcl[i]);
			wall4.setForce(ptcl[i]);
		}
	}
	static void postTimestep(PS::ParticleSystem<ThisPtcl>& ptcl, const system_t& sysinfo){
		static double time = 0;
		if(sysinfo.time < time) return;
		PS::F64 vel = 0;
		PS::F64 avel = 0;
		for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
			vel += ptcl[i].vel * ptcl[i].vel;
			avel += ptcl[i].avel * ptcl[i].avel;
		}
		std::cout << "RMS of vel and avel at t = " << time << " : " << vel << " " << avel << std::endl;
		time += sysinfo.end_time / 50.0;
	}
};

