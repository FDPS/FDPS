#pragma once

class CalcDensity{
	kernel_t kernel;
	public:
	void operator () (const EPI::Dens* const ep_i, const PS::S32 Nip, const EPJ::Dens* const ep_j, const PS::S32 Njp, RESULT::Dens* const dens){
		for(PS::S32 i = 0 ; i < Nip ; ++ i){
			dens[i].clear();
			for(PS::S32 j = 0 ; j < Njp ; ++ j){
				const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
				dens[i].dens += ep_j[j].mass * kernel.W(dr, ep_i[i].smth);
			}
			dens[i].smth = PARAM::SMTH * pow(ep_i[i].mass / dens[i].dens, 1.0/(PS::F64)(PARAM::Dim));
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
			for(PS::S32 j = 0; j < Njp ; ++ j){
				const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
				const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
				const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
				const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;
				v_sig_max = std::max(v_sig_max, v_sig);
				const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens + ep_j[j].dens));
				const PS::F64vec gradW = 0.5 * (kernel.gradW(dr, ep_i[i].smth) + kernel.gradW(dr, ep_j[j].smth));
				hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) * gradW;
				hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens * ep_i[i].dens) + 0.5 * AV) * dv * gradW;
			}
			hydro[i].dt = PARAM::C_CFL * 2.0 * ep_i[i].smth / v_sig_max;
		}
	}
};
