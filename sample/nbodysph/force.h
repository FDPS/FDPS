#pragma once

class CalcDensity{
	kernel_t kernel;
	public:
	void operator () (const EPI::Dens* const ep_i, const PS::S32 Nip, const EPJ::Dens* const ep_j, const PS::S32 Njp, RESULT::Dens* const dens){
		for(PS::S32 i = 0 ; i < Nip ; ++ i){
			dens[i].clear();
			const EPI::Dens& ith = ep_i[i];
			for(PS::S32 j = 0 ; j < Njp ; ++ j){
				const EPJ::Dens& jth = ep_j[j];
				const PS::F64vec dr = jth.pos - ith.pos;
				dens[i].dens += jth.mass * kernel.W(dr, ith.smth);
			}
			dens[i].smth = PARAM::SMTH * pow(ith.mass / dens[i].dens, 1.0/(PS::F64)(PARAM::Dim));
		}
	}
};

class CalcDerivative{
	kernel_t kernel;
	public:
	void operator () (const EPI::Drvt* ep_i, const PS::S32 Nip, const EPJ::Drvt* ep_j, const PS::S32 Njp, RESULT::Drvt* const drvt){
		for(PS::S32 i = 0; i < Nip ; ++ i){
			drvt[i].clear();
			const EPI::Drvt& ith = ep_i[i];
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const EPJ::Drvt& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64vec dv = ith.vel - jth.vel;
				drvt[i].div_v += - jth.mass * dv * kernel.gradW(dr, ith.smth);
				drvt[i].rot_v += - jth.mass * dv ^ kernel.gradW(dr, ith.smth);
			}
			drvt[i].div_v /= ith.dens;
			drvt[i].rot_v /= ith.dens;
		}
	}
};

class CalcHydroForce{
	const kernel_t kernel;
	public:
	void operator () (const EPI::Hydro* const ep_i, const PS::S32 Nip, const EPJ::Hydro* const ep_j, const PS::S32 Njp, RESULT::Hydro* const hydro){
		for(PS::S32 i = 0; i < Nip ; ++ i){
			hydro[i].clear();
			PS::F64 v_sig_max = 0.0;
			const EPI::Hydro& ith = ep_i[i];
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const EPJ::Hydro& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64vec dv = ith.vel - jth.vel;
				const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
				const PS::F64 v_sig = ith.snds + jth.snds - 3.0 * w_ij;
				v_sig_max = std::max(v_sig_max, v_sig);
				#if 1 //With Balsara Switch
				const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ith.dens + jth.dens)) * 0.5 * (ith.Bal + jth.Bal);
				//if(Bal_i > 0.1 && ep_i[i].pos.x < 0.6)std::cout << ep_i[i].pos.x << ", " << Bal_i << std::endl;
				#else
				const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ith.dens + jth.dens));
				#endif
				const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ith.smth) + kernel.gradW(dr, jth.smth));
				hydro[i].acc     -= jth.mass * (ith.pres / (ith.dens * ith.dens) + jth.pres / (jth.dens * jth.dens) + AV) * gradW;
				hydro[i].eng_dot += jth.mass * (ith.pres / (ith.dens * ith.dens) + 0.5 * AV) * dv * gradW;
			}
			hydro[i].dt = PARAM::C_CFL * 2.0 * ith.smth / v_sig_max;
		}
	}
};

template <class TParticleJ> class CalcGravityForce{
	public:
	void operator () (const EPI::Grav* const ep_i, const PS::S32 Nip, const TParticleJ* const ep_j, const PS::S32 Njp, RESULT::Grav* const grav){
		for(PS::S32 i = 0; i < Nip ; ++ i){
			const EPI::Grav& ith = ep_i[i];
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const TParticleJ& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64 dr2 = dr * dr;
				const PS::F64 dr_inv = 1.0 / sqrt(dr2 + ith.getEps2());
				const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
				grav[i].acc += - m_dr3_inv * dr;
				grav[i].pot += - jth.mass * dr_inv;
			}
		}
	}
};

/*
void CalcGravEpEp(const ThisRun::EPI::Grav* ep_i, const PS::S32 n_ip, const ThisRun::EPJ::Grav* ep_j, const PS::S32 n_jp, ThisRun::Result::Grav* grav){
	for(PS::S32 i = 0; i < n_ip ; ++ i){
		grav[i].clear();
		for(PS::S32 j = 0; j < n_jp ; ++ j){
			if(ep_i[i].id == ep_j[j].id) continue;
			const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
			const PS::F64 dr2 = dr * dr;
			const PS::F64 dr_inv = 1.0 / sqrt(dr2 + 0.5 * (ep_i[i].getEps2() + ep_j[j].getEps2()));
			const PS::F64 m_dr3_inv = ep_j[j].mass * math::pow3(dr_inv);
			grav[i].acc += - m_dr3_inv * dr;
			grav[i].pot += - ep_j[j].mass * dr_inv;
		}
	}
}
*/
/*
void CalcGravEpSp(const ThisRun::EPI::Grav* ep_i, const PS::S32 n_ip, const ThisRun::SPJ::Grav* ep_j, const PS::S32 n_jp, ThisRun::Result::Grav* grav){
	for(PS::S32 i = 0; i < n_ip ; ++ i){
		//grav[i].clear();
		for(PS::S32 j = 0; j < n_jp ; ++ j){
			const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
			const PS::F64 dr2 = dr * dr;
			const PS::F64 dr_inv = 1.0 / sqrt(dr2 + Rtmp * Rtmp * ep_i[i].getEps2());
			const PS::F64 m_dr3_inv = ep_j[j].mass * math::pow3(dr_inv);
			grav[i].acc += - Gtmp * m_dr3_inv * dr;
			grav[i].pot += - Gtmp * ep_j[j].mass * dr_inv;
		}
	}
}
*/
