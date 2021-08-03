#pragma once

class FileHeader{
	public:
	int Nbody;
	double time;
	int readAscii(FILE* fp){
		fscanf(fp, "%lf\n", &time);
		fscanf(fp, "%d\n", &Nbody);
		return Nbody;
	}
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%e\n", time);
		fprintf(fp, "%d\n", Nbody);
	}
};

class Material{
	//bulk modulus, shear modulus and reference density
	PS::F64 bulk, shear, dens;
	public:
	Material(): bulk(0.0), shear(0.0), dens(0.0){
	}
	Material(const PS::F64 bulk_, const PS::F64 shear_, const PS::F64 dens_): bulk(bulk_), shear(shear_), dens(dens_){
	}
	PS::F64 getBulkModulus() const{
		return bulk;
	}
	PS::F64 getShearModulus() const{
		return shear;
	}
	PS::F64 getDensity() const{
		return dens;
	}
	PS::F64 getYoungModulus() const{
		return 9.0 * shear * bulk / (3.0 * bulk + shear);
	}
	PS::F64 getPoissonRatio() const{
		return (3.0 * bulk - 2.0 * shear) / (6.0 * bulk + 2.0 * shear);
	}
};

const Material Basalt(26.70e+9, 22.7e+9, 2700.0);
const Material Ice   ( 9.47e+9,  2.8e+9,  917.0);
const Material Pyrex (21.20e+9, 28.2e+9, 2235.0);
const Material Test ( 5.56e+6, 4.17e+6, 2500.0);
//const Material Test  ( 1.67e+6, 0.77e+6, 2500.0);

class FP{
	public:
	//id
	PS::S64    id;
	//tag
	PS::S64    tag;
	//mass
	PS::F64    mass;
	//inertia
	PS::F64    iner;
	//position & angle
	PS::F64vec pos;
	PS::F64vec ang;
	//velocity
	PS::F64vec vel, vel_half;
	PS::F64vec avel, avel_half;
	//acceleration
	PS::F64vec acc;
	//torque
	PS::F64vec tor;
	//radius
	PS::F64    rad;
	//timestep
	PS::F64    dt;
	//material parameter
	Material   mat;
	//tangential displacement between surrounding ptcls.
	Hash<16, PS::F64vec> t_displ;
	Hash<16, PS::F64vec> t_vel;

	FP(){
		id = 0;
		avel = 0;
	}
	PS::F64vec getPos() const{
		return pos;
	}
	PS::F64 getCharge() const{
		return mass;
	}
	PS::F64 getRSearch() const{
		return 2.0 * rad;
	}
	void setPos(const PS::F64vec& pos){
		this->pos = pos;
	}
	void copyFromForce(const Force& force){
		this->acc += force.acc;
		this->tor += force.tor;
		this->dt    = force.dt;
		this->t_vel = force.t_vel;
	}
	void copyFromFP(const FP& fp){
		*this = fp;
	}
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", id, mass, rad, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, avel.x, avel.y ,avel.z, ang.x, ang.y, ang.z);
	}
	//KDK
	void kick(const PS::F64 dt){
		vel_half = vel + 0.5 * acc * dt;
		avel_half = avel + 0.5 * tor * dt;
	}
	void drift(const PS::F64 dt){
		pos += vel_half * dt;
		ang += avel_half * dt;
		vel += acc * dt;
		avel += tor * dt;
		//update tangential displacement
		Hash<16, PS::F64vec> t_displ_tmp;
		for(int i = 0 ; i < t_vel.size() ; ++ i){
			std::size_t key = t_vel.key[i];
			t_displ_tmp[key] = t_displ[key] + t_vel[key] * dt;
		}
		t_displ = t_displ_tmp;
	}
	void kick2(const PS::F64 dt){
		vel = vel_half + 0.5 * acc * dt;
		avel = avel_half + 0.5 * tor * dt;
	}
	//
	void clear(){
		acc = 0.0;
		tor = 0.0;
	}
};

struct system_t{
	PS::F64 time, dt, end_time;
	unsigned int step;
	system_t() : time(0.0), step(0), end_time(0.0), dt(1.0e+30){
	}
};


