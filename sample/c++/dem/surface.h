#pragma once
template <class ThisPtcl> class Wall{
	public:
	virtual void setForce(ThisPtcl&, const PS::F64) const = 0;
};

template <class ThisPtcl> class InfinitePlane : public Wall<ThisPtcl>{
	const PS::F64vec n;//normal vector
	const PS::F64vec O;//origin vector
	const PS::F64vec V;//velocity vector of the plane
	static constexpr PS::F64 cn = 0.4;
	static constexpr PS::F64 mu_s = 0.6;//coefficient for sliding friction
	static constexpr PS::F64 mu_r = 0.5;//coefficient for rolling friction
	const Material mat;
	public:
	InfinitePlane(const Material mat_, const PS::F64vec n_, const PS::F64vec O_, const PS::F64vec V_ = 0): mat(mat_), n(n_ / sqrt(n_ * n_)), O(O_), V(V_){
	}
	void dispParameter(){
		std::cout << "See Schwartz et al. (2012)" << std::endl;
		std::cout << "normal vector  : " << n << std::endl;
		std::cout << "origin vector  : " << O << std::endl;
		std::cout << "velocity vector: " << V << " (currently untested)" << std::endl;
		return ;
	}
	void setForce(ThisPtcl& ith, const PS::F64 time = 0) const{
		const PS::F64vec dr = ith.pos - (O + V * time);
		const PS::F64vec dv = ith.vel + ((ith.rad * n) ^ ith.avel) - V;
		const PS::F64vec dr_n = (dr * n) * n;
		const PS::F64vec dv_n = (dv * n) * n;
		const PS::F64vec dr_t = dr - dr_n;
		const PS::F64vec dv_t = dv - dv_n;
		const PS::F64 pd = ith.rad - sqrt(dr_n * dr_n);
		if(pd < 0) return;
		const PS::F64    Eij = 1.0 / ((1.0 - pow(ith.mat.getPoissonRatio(), 2)) / ith.mat.getYoungModulus() + (1.0 - pow(mat.getPoissonRatio(), 2)) / mat.getYoungModulus());
		const PS::F64    rij = ith.rad;
		const PS::F64    mij = ith.mass;
		const PS::F64    knr = 4.0 / 3.0 * Eij * sqrt(rij);
		const PS::F64vec Fnr = (dr * n < 0 ? -1.0 : 1.0) * knr * powf(pd, 1.5) * n;
		const PS::F64vec Fnd = - cn * sqrt(6.0 * mij * Eij * sqrt(rij * pd)) * dv_n;
		//temporary
		const PS::F64vec Ftd = - mu_s * sqrt(Fnr * Fnr) * dv_t / sqrt(dv_t * dv_t + 1.0e-16);
		ith.acc += (Fnr + Fnd + Ftd) / ith.mass;
		ith.tor += (- ith.rad * (n ^ Ftd) - mu_r * ith.rad * sqrt(Fnr * Fnr) * ith.avel / sqrt(ith.avel * ith.avel + 1.0e-4)) / ith.iner;
		return;
	}
};

template <class ThisPtcl> class InfiniteCylinder : public Wall<ThisPtcl>{
	/**************************************************
	* force to Z-axis and torque haven't imposed yet. *
	**************************************************/

	const PS::F64vec Z;//Axis-directional vector
	const PS::F64vec O;//origin vector
	const PS::F64    R;//radius
	const PS::F64vec V;//velocity vector
	/*const*/ PS::F64    Vang;//rotation velocity
	static constexpr PS::F64 cn = 0.3;
	static constexpr PS::F64 mu_s = 0.5;//coefficient for sliding friction
	static constexpr PS::F64 mu_r = 0.5;//coefficient for rolling friction
	const Material mat;
	public:
	InfiniteCylinder(const Material mat_, const PS::F64vec Z_, const PS::F64vec O_, const PS::F64 R_, const PS::F64vec V_ = 0, const PS::F64 Vang_ = 0): mat(mat_), Z(Z_ / sqrt(Z_ * Z_)), O(O_), R(R_), V(V_), Vang(Vang_){
		//Vang = 2.0 * M_PI / 50.0;
	}
	void dispParameter(){
		std::cout << "See Schwartz et al. (2012)" << std::endl;
		std::cout << "Axis vector    : " << Z << std::endl;
		std::cout << "origin vector  : " << O << std::endl;
		std::cout << "radius         : " << R << std::endl;
		std::cout << "velocity vector: " << V << " (currently not tested)" << std::endl;
		return ;
	}
	void setForce(ThisPtcl& ith, const PS::F64 time = 0) const{
		const PS::F64vec dr   = ith.pos - O;
		const PS::F64vec dr_Z = (dr * Z) * Z;//along with Z-axis
		const PS::F64vec dr_r = dr - dr_Z;   //along with R-axis
		const PS::F64vec dr_t = dr_r ^ Z;    //along with theta-axis
		//const PS::F64vec dv   = ith.vel + (ith.rad * ith.avel) ^ dr_r / sqrt(dr_r * dr_r) - (V + (R * Vang * Z) ^ (- dr_r) / sqrt(dr_r * dr_r));
		const PS::F64vec dv   = ith.vel + ith.rad * (ith.avel ^ dr_r) / sqrt(dr_r * dr_r) - (V + R * Vang * (Z ^ (-dr_r)) / sqrt(dr_r * dr_r));
		const PS::F64vec dv_Z = (dv * Z) * Z;
		const PS::F64vec dv_r = (dv * dr_r / sqrt(dr_r * dr_r)) * dr_r / sqrt(dr_r * dr_r);
		const PS::F64vec dv_t = dv - dv_r - dv_Z;//theta
		//inside the cylinder & outside it
		const PS::F64 pd      = std::min(ith.rad + sqrt(dr_r * dr_r) - R, R + ith.rad - sqrt(dr_r * dr_r));
		const PS::F64 sign    = (ith.rad + sqrt(dr_r * dr_r) - R < R + ith.rad - sqrt(dr_r * dr_r)) ? -1.0 : 1.0;
		if(!(std::max(R - ith.rad, 0.0) < sqrt(dr_r * dr_r) && sqrt(dr_r * dr_r) < R + ith.rad)) return;
		const PS::F64    Eij = 1.0 / ((1.0 - pow(ith.mat.getPoissonRatio(), 2)) / ith.mat.getYoungModulus() + (1.0 - pow(mat.getPoissonRatio(), 2)) / mat.getYoungModulus());
		const PS::F64    rij = ith.rad;
		const PS::F64    mij = ith.mass;
		const PS::F64    knr = 4.0 / 3.0 * Eij * sqrt(rij);
		const PS::F64vec Fnr = sign * knr * powf(pd, 1.5) * dr_r;
		const PS::F64vec Fnd = - cn * sqrt(6.0 * mij * Eij * sqrt(rij * pd)) * (dv_r);
		const PS::F64vec Ftd = - mu_s * sqrt(Fnr * Fnr) * dv_t / sqrt(dv_t * dv_t + 1.0e-4);
		#if 0
		std::cout << ith.pos << std::endl;
		std::cout << "dr  : " << dr   << std::endl;
		std::cout << "dr_r: " << dr_r << std::endl;
		std::cout << "dr_Z: " << dr_Z << std::endl;
		std::cout << "dr_t: " << dr_t << std::endl;
		std::cout << "dv  : " << dv   << std::endl;
		std::cout << "dv_r: " << dv_r << std::endl;
		std::cout << "dv_t: " << dv_t << std::endl;
		std::cout << "dv_Z: " << dv_Z << std::endl;
		std::cout << "Fnr :" << Fnr << std::endl;
		std::cout << "Fnd :" << Fnd << std::endl;
		std::cout << "Ftd :" << Ftd << std::endl;
		getchar();
		#endif
		ith.acc += (Fnr + Fnd + Ftd) / ith.mass;
		//ith.tor += (/*- ith.rad * (dr_r / sqrt(dr_r * dr_r) ^ Ftd)*/ - mu_r * ith.rad * sqrt(Fnr * Fnr) * ith.avel / sqrt(ith.avel * ith.avel + 1.0e-4)) / ith.iner;
		return;
	}
};

template <class ThisPtcl> class Disc : public Wall<ThisPtcl>{
	const PS::F64vec n;//normal vector
	const PS::F64vec O;//origin vector
	const PS::F64vec V;//velocity vector of the disc
	const PS::F64    R_in;//inner edge
	const PS::F64    R_out;//outer edge
	static constexpr PS::F64 cn = 0.3;
	static constexpr PS::F64 mu_s = 0.1;//coefficient for sliding friction
	static constexpr PS::F64 mu_r = 0.1;//coefficient for rolling friction
	const Material mat;
	public:
	Disc(const Material mat_, const PS::F64vec n_, const PS::F64vec O_, const PS::F64 R_in_, const PS::F64 R_out_, const PS::F64vec V_ = 0): mat(mat_), n(n_ / sqrt(n_ * n_)), O(O_), R_in(R_in_), R_out(R_out_), V(V_){
		assert(R_in < R_out);
		assert(R_in >= 0);
		assert(R_out > 0);
	}
	void dispParameter(){
		std::cout << "See Schwartz et al. (2012)" << std::endl;
		std::cout << "normal vector  : " << n << std::endl;
		std::cout << "origin vector  : " << O << std::endl;
		std::cout << "inner edge     : " << R_in << std::endl;
		std::cout << "outer edge     : " << R_out << std::endl;
		std::cout << "velocity vector: " << V << " (currently untested)" << std::endl;
		return ;
	}
	void setForce(ThisPtcl& ith, const PS::F64 time = 0) const{
		const PS::F64vec dr = ith.pos - (O + V * time);
		const PS::F64vec dv = ith.vel + ((ith.rad * n) ^ ith.avel) - V;
		const PS::F64vec dr_n = (dr * n) * n;
		const PS::F64vec dv_n = (dv * n) * n;
		const PS::F64vec dr_t = dr - dr_n;
		const PS::F64vec dv_t = dv - dv_n;
		const PS::F64 pd = ith.rad - sqrt(dr_n * dr_n);
		//first cut
		if(pd < 0) return;
		const PS::F64 dr_t_norm = sqrt(dr_t * dr_t);
		if(R_in <= dr_t_norm && dr_t_norm <= R_out){
			const PS::F64    Eij = 1.0 / ((1.0 - pow(ith.mat.getPoissonRatio(), 2)) / ith.mat.getYoungModulus() + (1.0 - pow(mat.getPoissonRatio(), 2)) / mat.getYoungModulus());
			const PS::F64    rij = ith.rad;
			const PS::F64    mij = ith.mass;
			const PS::F64    knr = 4.0 / 3.0 * Eij * sqrt(rij);
			const PS::F64vec Fnr = (dr * n < 0 ? -1.0 : 1.0) * knr * powf(pd, 1.5) * n;
			const PS::F64vec Fnd = - cn * sqrt(6.0 * mij * Eij * sqrt(rij * pd)) * dv_n;
			//temporary
			const PS::F64vec Ftd = - mu_s * sqrt(Fnr * Fnr) * dv_t / sqrt(dv_t * dv_t + 1.0e-4);
			//the particle touches the flat portion of the disc
			ith.acc += (Fnr + Fnd + Ftd) / ith.mass;
			ith.tor += (- ith.rad * (n ^ Ftd) - mu_r * ith.rad * sqrt(Fnr * Fnr) * ith.avel / sqrt(ith.avel * ith.avel + 1.0e-4)) / ith.iner;
		}else if(R_out < dr_t_norm){
			if(sqrt((dr - R_out * dr_t / dr_t_norm) * (dr - R_out * dr_t / dr_t_norm)) > ith.rad) return;
			const PS::F64vec contact = O + R_out * dr_t / dr_t_norm;
			const PS::F64vec force_dir = (ith.pos - contact) / sqrt((ith.pos - contact) * (ith.pos - contact));
			const PS::F64    Eij = 1.0 / ((1.0 - pow(ith.mat.getPoissonRatio(), 2)) / ith.mat.getYoungModulus() + (1.0 - pow(mat.getPoissonRatio(), 2)) / mat.getYoungModulus());
			const PS::F64    rij = ith.rad;
			const PS::F64    mij = ith.mass;
			const PS::F64    knr = 4.0 / 3.0 * Eij * sqrt(rij);
			const PS::F64vec Fnr = (dr * n < 0 ? -1.0 : 1.0) * knr * powf(pd, 1.5) * force_dir;
			const PS::F64vec Fnd = - cn * sqrt(6.0 * mij * Eij * sqrt(rij * pd)) * dv_n;
			//temporary
			const PS::F64vec Ftd = - mu_s * sqrt(Fnr * Fnr) * dv_t / sqrt(dv_t * dv_t + 1.0e-4);
			//the particle touches the flat portion of the disc
			ith.acc += (Fnr + Fnd + Ftd) / ith.mass;
			ith.tor += (- ith.rad * (n ^ Ftd) - mu_r * ith.rad * sqrt(Fnr * Fnr) * ith.avel / sqrt(ith.avel * ith.avel + 1.0e-4)) / ith.iner;
		}else if(dr_t_norm < R_in){
			if(dr_t_norm == 0.0){
				if(dr * dr + R_in * R_in > ith.rad * ith.rad) return;
				//phantom contact
				const PS::F64vec contact = - (ith.rad - sqrt(ith.rad * ith.rad - R_in * R_in)) * dr_n / sqrt(dr_n * dr_n);
				const PS::F64    pd_phantom = sqrt(ith.rad * ith.rad - R_in * R_in) - sqrt(dr_n * dr_n);
				const PS::F64vec force_dir = n;
				const PS::F64    Eij = 1.0 / ((1.0 - pow(ith.mat.getPoissonRatio(), 2)) / ith.mat.getYoungModulus() + (1.0 - pow(mat.getPoissonRatio(), 2)) / mat.getYoungModulus());
				const PS::F64    rij = ith.rad;
				const PS::F64    mij = ith.mass;
				const PS::F64    knr = 4.0 / 3.0 * Eij * sqrt(rij);
				const PS::F64vec Fnr = (dr * n < 0 ? -1.0 : 1.0) * knr * powf(pd_phantom, 1.5) * force_dir;
				const PS::F64vec Fnd = - cn * sqrt(6.0 * mij * Eij * sqrt(rij * pd_phantom)) * dv_n;
				//temporary
				const PS::F64vec Ftd = - mu_s * sqrt(Fnr * Fnr) * dv_t / sqrt(dv_t * dv_t + 1.0e-4);
				//the particle touches the flat portion of the disc
				ith.acc += (Fnr + Fnd + Ftd) / ith.mass;
				ith.tor += (- ith.rad * (n ^ Ftd) - mu_r * ith.rad * sqrt(Fnr * Fnr) * ith.avel / sqrt(ith.avel * ith.avel + 1.0e-4)) / ith.iner;
			}else{
				if(sqrt((dr - R_in * dr_t / dr_t_norm) * (dr - R_in * dr_t / dr_t_norm)) > ith.rad) return;
				const PS::F64vec contact = O + R_in * dr_t / dr_t_norm;
				const PS::F64vec force_dir = (ith.pos - contact) / sqrt((ith.pos - contact) * (ith.pos - contact));
				const PS::F64    Eij = 1.0 / ((1.0 - pow(ith.mat.getPoissonRatio(), 2)) / ith.mat.getYoungModulus() + (1.0 - pow(mat.getPoissonRatio(), 2)) / mat.getYoungModulus());
				const PS::F64    rij = ith.rad;
				const PS::F64    mij = ith.mass;
				const PS::F64    knr = 4.0 / 3.0 * Eij * sqrt(rij);
				const PS::F64vec Fnr = (dr * n < 0 ? -1.0 : 1.0) * knr * powf(pd, 1.5) * force_dir;
				const PS::F64vec Fnd = - cn * sqrt(6.0 * mij * Eij * sqrt(rij * pd)) * dv_n;
				//temporary
				const PS::F64vec Ftd = - mu_s * sqrt(Fnr * Fnr) * dv_t / sqrt(dv_t * dv_t + 1.0e-4);
				//the particle touches the flat portion of the disc
				ith.acc += (Fnr + Fnd + Ftd) / ith.mass;
				ith.tor += (- ith.rad * (n ^ Ftd) - mu_r * ith.rad * sqrt(Fnr * Fnr) * ith.avel / sqrt(ith.avel * ith.avel + 1.0e-4)) / ith.iner;
			}
			
		}
		return;
	}
};

