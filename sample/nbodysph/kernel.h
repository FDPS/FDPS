#pragma once

struct kernel_t{
	kernel_t(){
	}
	//W
	PS::F64 W(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = supportRadius() * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
		r_value = math::plus(math::pow3(1.0 - s)) - 4.0 * math::plus(math::pow3(0.5 - s));
		//#if N_DIM == 1
		//r_value *= 8.0 / 3.0 / H;
		//#elif N_DIM == 2
		//r_value *= 80.0 / (7.0 * PI) / (H * H);
		//#elif N_DIM == 3
		r_value *= 16.0 / math::pi / (H * H * H);
		//#endif
		return r_value;
	}
	//gradW
	PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h) const{
		const PS::F64 H = supportRadius() * h;
		const PS::F64 s = sqrt(dr * dr) / H;
		PS::F64 r_value;
		r_value = - 3.0 * math::pow2(math::plus(1.0 - s)) + 12.0 * math::pow2(math::plus(0.5 - s));
		//#if N_DIM == 1
		//r_value *= 8.0 / 3.0 / H;
		//#elif N_DIM == 2
		//r_value *= 80.0 / (7.0 * PI) / (H * H);
		//#elif N_DIM == 3
		r_value *= 16.0 / math::pi / (H * H * H);
		//#endif
		return dr * r_value / (sqrt(dr * dr) * H  + 1.0e-6 * h);
	}
	static PS::F64 supportRadius(){
		return 2.5;
	}
};
