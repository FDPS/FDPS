#pragma once
template <int MaxSize, typename real = double> class Hash{
	public:
	real value[MaxSize];
	std::size_t key[MaxSize];
	std::size_t tail;
	Hash(): tail(0){
	}
	void clear(){
		tail = 0;
	}
	real& operator [](const std::size_t key_){
		if(tail == MaxSize){
			std::cout << "range over 0:" << tail << std::endl;
			dump();
			exit(0);
		}
		for(std::size_t idx = 0 ; idx < tail ; ++ idx){
			if(key_ == key[idx]) return value[idx];
		}
		key[tail] = key_;
		return value[tail++];
	}
	real operator [](const std::size_t key_) const{
		for(std::size_t idx = 0 ; idx < tail ; ++ idx){
			if(key_ == key[idx]) return value[idx];
		}
		return +0.0;
	}
	real getValue(const std::size_t key_) const{
		for(std::size_t idx = 0 ; idx < tail ; ++ idx){
			if(key_ == key[idx]) return value[idx];
		}
		return +0.0/0.0;
	}
	std::size_t size() const{
		return tail;
	}
	void dump(){
		std::cout << "tail: " << tail  << std::endl;
		for(int i = 0 ; i < tail ; ++ i){
			std::cout << "   idx: " << i << ", key: " << key[i] << ", value: " << value[i] << std::endl;
		}
	}
};


class Force{
	static constexpr PS::F64 cn = 0.4;
	static constexpr PS::F64 ct = 0.4;
	static constexpr PS::F64 mu_s = 0.6;//coefficient for sliding friction
	static constexpr PS::F64 mu_r = 0.5;//coefficient for rolling friction in metre
	public:
	PS::F64vec acc;
	PS::F64vec tor;
	Hash<16, PS::F64vec> t_vel;//tangential velocity
	PS::F64 dt;
	void clear(){
		acc = tor = 0;
		t_vel.clear();
	}
	template <class ThisPtcl> void operator () (const ThisPtcl* const ep_i, const PS::S32 Nip, const ThisPtcl* const ep_j, const PS::S32 Njp, Force* const force){
		for(PS::S32 i = 0 ; i < Nip ; ++ i){
			const ThisPtcl& ith = ep_i[i];
			force[i].dt = M_PI / 50.0 * sqrt(ith.mass / (2.0 / 3.0 * ith.mat.getYoungModulus() * sqrt(ith.rad * 0.5)));
			for(PS::S32 j = 0 ; j < Njp ; ++ j){
				const ThisPtcl& jth = ep_j[j];
				const PS::F64vec dr = jth.pos - ith.pos;
				const PS::F64    pd = ith.rad + jth.rad - sqrt(dr * dr);
				if(pd < 0 || dr * dr <= 0.0) continue;
				//unit normal vec.
				const PS::F64vec n  = dr / sqrt(dr * dr);
				const PS::F64vec dv = jth.vel - ((jth.rad * n) ^ jth.avel) - ith.vel + ((ith.rad * n) ^ ith.avel);
				//unit tangential vec.
				const PS::F64vec dv_n = (dv * n) * n;
				const PS::F64vec dv_t = dv - dv_n;
				//max tangential displacement
				const PS::F64 pr  = 2.0 * (1.0 / ith.mat.getPoissonRatio() + 1.0 / jth.mat.getPoissonRatio());//Poisson ratio
				const PS::F64 Eij = 1.0 / ((1.0 - pow(ith.mat.getPoissonRatio(), 2)) / ith.mat.getYoungModulus() + (1.0 - pow(jth.mat.getPoissonRatio(), 2)) / jth.mat.getYoungModulus());
				const PS::F64 rij = 1.0 / (1.0 / ith.rad + 1.0 / jth.rad);
				const PS::F64 mij = 1.0 / (1.0 / ith.mass + 1.0 / jth.mass);
				//repulsive normal force
				const PS::F64 knr = 4.0 / 3.0 * Eij * sqrt(rij);
				const PS::F64vec Fnr = - knr * powf(pd, 1.5) * n;
				//damp normal force
				const PS::F64 Cn = cn * sqrt(6.0 * mij * Eij * sqrt(rij * pd));
				const PS::F64vec Fnd = Cn * dv_n;
				//repulsive tangential force
				const PS::F64 t_displ_max = mu_s * (2.0 - pr) / (1.0 - pr) * 0.5 * pd;
				const PS::F64 ktr = mu_s * sqrt(Fnr * Fnr) / sqrt(ith.t_displ[jth.id] * ith.t_displ[jth.id] + 1.0e-4) * (1.0 - powf(1.0 - std::min(sqrt(ith.t_displ[jth.id] * ith.t_displ[jth.id]), t_displ_max) / t_displ_max, 1.5));
				const PS::F64vec Ftr = ktr * ith.t_displ[jth.id];
				//damp tangential force
				const PS::F64 Ct = ct * sqrt(6.0 * mij * mu_s * sqrt(Fnr * Fnr) * sqrt(1.0 - std::min(sqrt(ith.t_displ[jth.id] * ith.t_displ[jth.id]), t_displ_max) / t_displ_max) / t_displ_max);
				const PS::F64vec Ftd = Ct * dv_t;
				//tangential vel.
				force[i].t_vel[jth.id] = dv_t;
				//resulting force
				force[i].acc += Fnr + Fnd + Ftr + Ftd;
				force[i].tor += ith.rad * (n ^ (Ftr + Ftd)) - mu_r * ith.rad * sqrt(Fnr * Fnr) * ith.avel / sqrt(ith.avel * ith.avel + 1.0e-16);
				force[i].dt  = std::min(M_PI / 50.0 * sqrt(mij / (knr * (1.0 - Cn * Cn / (knr * knr)))), force[i].dt);
			}
			force[i].acc /= ith.mass;
			force[i].tor /= ith.iner;
		}
	}
};


