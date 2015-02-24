#pragma once

struct boundary{
	PS::F64 x, y, z;
};

class FileHeader{
public:
	int Nbody;
	double time;
	int readAscii(FILE* fp){
		fscanf(fp, "%e\n", &time);
		fscanf(fp, "%d\n", &Nbody);
		return Nbody;
	}
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%e\n", time);
		fprintf(fp, "%d\n", Nbody);
	}
	/*
	int readBinary(FILE* fp){
		fread(this, sizeof(*this), 1, fp);
		return Nbody;
	}
	void writeBinary(FILE* fp) const{
		//std::cout << Nbody << std::endl;
		fwrite(this, sizeof(*this), 1, fp);
	}
	*/
};

namespace RESULT{
	//Density summation
	class Dens{
		public:
		PS::F64 dens;
		PS::F64 smth;
		bool isNotConverged;
		void clear(){
			dens = smth = 0;
			isNotConverged = true;
		}
	};
	class Drvt{
		public:
		PS::F64    div_v;
		PS::F64vec rot_v;
		void clear(){
			div_v = 0;
			rot_v = 0;
		}
	};
	class Hydro{
		public:
		PS::F64vec acc;
		PS::F64 eng_dot;
		PS::F64 dt;
		void clear(){
			acc = 0;
			eng_dot = 0;
			dt = 1.0e+30;
		}
	};
}

class RealPtcl{
	public:
	PS::F64 mass;
	PS::F64vec pos, vel;
	struct{
		PS::F64vec acc;
	}hydro, grav, ext;
	PS::F64 dens;//DENSity
	PS::F64 eng; //ENerGy
	PS::F64 pres;//PRESsure
	PS::F64 smth;//SMooTHing length

	PS::F64 div_v;
	PS::F64vec rot_v;//rot_v.z
	PS::F64 snds; //SouND Speed
	PS::F64 pot;

	PS::F64 eng_dot;
	PS::F64vec vel_half;
	PS::F64 eng_half;
	PS::F64 dt;
	PS::S64 id;
	bool isNotConverged;
	//Copy functions
	void copyFromForce(const RESULT::Dens& dens){
		this->dens = dens.dens;
		this->smth = dens.smth;
		this->isNotConverged = dens.isNotConverged;
	}
	void copyFromForce(const RESULT::Drvt& drvt){
		this->div_v = drvt.div_v;
		this->rot_v = drvt.rot_v;
	}
	void copyFromForce(const RESULT::Hydro& force){
		this->hydro.acc = force.acc;
		this->eng_dot   = force.eng_dot;
		this->dt        = force.dt;
	}
	//Give necessary values to FDPS
	PS::F64 getCharge() const{
		return this->mass;
	}
	PS::F64vec getPos() const{
		return this->pos;
	}
	PS::F64 getRSearch() const{
		#warning TEMPORARY
		const PS::F64 KernelSearchRadius = 2.0;
		return KernelSearchRadius * this->smth;
	}
	void setPos(const PS::F64vec& pos){
		this->pos = pos;
	}
	void writeAscii(FILE* fp) const{
		fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z, this->dens, this->eng, this->smth);
	}
	void readAscii(FILE* fp){
		fscanf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z, &this->dens, &this->eng, &this->smth);
	}
	/*
	void writeBinary(FILE* fp) const{
		fwrite(this, sizeof(this), 1, fp);
	}
	void readBinary(FILE* fp){
		fread(this, sizeof(this), 1, fp);
	}
	*/
	void setPressure(){
		PS::F64 hcr = 1.4;//heat capacity ratio
		pres = (hcr - 1.0) * dens * eng;
		snds = sqrt(hcr * pres / dens);
	}
};

namespace EPI{
	class Dens{
	public:
		PS::F64vec pos;
		PS::F64    mass;
		PS::F64    smth;
		bool       isNotConverged;
		void copyFromFP(const RealPtcl& rp){
			this->pos  = rp.pos;
			this->mass = rp.mass;
			this->smth = rp.smth;
			this->isNotConverged = rp.isNotConverged;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		PS::F64 getRSearch() const{
			#warning TEMPORARY
			return 2.0 * this->smth;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
	};
	class Drvt{
	public:
		PS::F64vec pos;
		PS::F64vec vel;
		PS::F64    smth;
		PS::F64    dens;
		void copyFromFP(const RealPtcl& rp){
			pos  = rp.pos;
			vel  = rp.vel;
			dens = rp.dens;
			smth = rp.smth;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		PS::F64 getRSearch() const{
			#warning TEMPORARY
			return 2.0 * this->smth;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
	};
	class Hydro{
	public:
		PS::F64vec pos;
		PS::F64vec vel;
		PS::F64    div_v;
		PS::F64vec rot_v;
		PS::F64    smth;
		PS::F64    dens;
		PS::F64    pres;
		PS::F64    snds;
		PS::S64    id;///DEBUG
		void copyFromFP(const RealPtcl& rp){
			this->pos   = rp.pos;
			this->vel   = rp.vel;
			this->div_v = rp.div_v;
			this->rot_v = rp.rot_v;
			this->smth  = rp.smth;
			this->dens  = rp.dens;
			this->pres  = rp.pres;
			this->snds = rp.snds;
			this->id   = rp.id;///DEBUG
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		PS::F64 getRSearch() const{
			#warning TEMPORARY
			return 2.0 * this->smth;
		}
	};
}

namespace EPJ{
	class Dens{
	public:
		PS::F64    mass;
		PS::F64vec pos;
		void copyFromFP(const RealPtcl& rp){
			this->mass = rp.mass;
			this->pos  = rp.pos;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
	};
	class Drvt{
	public:
		PS::F64    mass;
		PS::F64vec pos;
		PS::F64vec vel;
		void copyFromFP(const RealPtcl& rp){
			this->mass = rp.mass;
			this->pos  = rp.pos;
			this->vel  = rp.vel;
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
	};
	class Hydro{
	public:
		PS::F64vec pos;
		PS::F64vec vel;
		PS::F64    div_v;
		PS::F64vec rot_v;
		PS::F64    dens;
		PS::F64    pres;
		PS::F64    smth;
		PS::F64    mass;
		PS::F64    snds;
		PS::S64    id;///DEBUG
		void copyFromFP(const RealPtcl& rp){
			this->pos   = rp.pos;
			this->vel   = rp.vel;
			this->div_v = rp.div_v;
			this->rot_v = rp.rot_v;
			this->dens  = rp.dens;
			this->pres  = rp.pres;
			this->smth  = rp.smth;
			this->mass  = rp.mass;
			this->snds  = rp.snds;
			this->id    = rp.id;///DEBUG
		}
		PS::F64vec getPos() const{
			return this->pos;
		}
		PS::F64 getRSearch() const{
			#warning TEMPORARY
			return 2.0 * this->smth;
		}
		void setPos(const PS::F64vec& pos){
			this->pos = pos;
		}
	};
}


