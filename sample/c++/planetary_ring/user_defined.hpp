#include<sstream>
#include"./fdps-util/my_lib.hpp"
#include"./fdps-util/fdps-util.hpp"

inline PS::F64vec ConvertCar2Cyl(const PS::F64vec & pos){
    const auto pos_r   = sqrt(pos.x*pos.x + pos.y*pos.y);
    const auto pos_phi = atan2(pos.y, pos.x);
    return PS::F64vec(pos_phi, pos_r, pos.z);
}

inline PS::F64vec ConvertCyl2Car(const PS::F64vec & pos){
    const auto cth = cos(pos.x);
    const auto sth = sin(pos.x);
    const auto r = pos.y;
    const auto pos_x = r*cth;
    const auto pos_y = r*sth;
    return PS::F64vec(pos_x, pos_y, pos.z);
}

class Force_t{
public:
    PS::F64vec acc; // total
    PS::F64    pot;
    PS::F64vec acc_dash; // dashpot
    void clear(){
        acc = 0.0;
        pot = 0.0;
        acc_dash = 0.0;
    }
};

void SumForce(Force_t in[], Force_t out[], const int n){
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
    MPI_Allreduce(in, out, (sizeof(Force_t)/sizeof(PS::F64))*n, MPI_DOUBLE, MPI_SUM, PS::Comm::getCommunicator());
#else
    for(auto i=0; i<n; i++){
	out[i] = in[i];
    }
#endif
}

class FP_t{
public:
    PS::F64vec pos_car; // cartesian
    PS::F64vec pos_cyl; // cyl
    PS::F64    mass;
    PS::F64vec vel; // cartesian
    PS::F64vec vel_full; // for calculating disipation energy
    PS::S64    id;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64 r_coll;
    PS::F64vec acc_dash;
    PS::F64 r_search;
    static PS::F64 kappa;
    static PS::F64 eta;
    static PS::F64 eps;
    //static inline PS::F64 kappa;
    //static inline PS::F64 eta;
    //static inline PS::F64 eps;
    
    PS::F64vec getPos() const {
        return pos_cyl;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void setPos(const PS::F64vec & pos_new){
        pos_cyl = pos_new;
    }

    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64 getRSearch() const {
        return r_search;
    }
    void copyFromForce(const Force_t & force) {
        acc = force.acc;
        pot = force.pot;
        acc_dash = force.acc_dash;
    }
    void accumFromForce(const Force_t & force) {
        acc += force.acc;
        pot += force.pot;
        acc_dash += force.acc_dash;
    }
    void clearForce(){
	acc = acc_dash = 0.0;
	pot = 0.0;
    }
    
    void writeAscii(FILE* fp) const {
        fprintf(fp, "%lld   %12.11e  %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e   %12.11e\n", 
                this->id, this->mass,
                this->pos_car.x, this->pos_car.y, this->pos_car.z,
                this->vel.x, this->vel.y, this->vel.z,
		this->pot, this->r_coll, this->r_search);
    }
    void readAscii(FILE* fp) {
        auto tmp = fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
			  &this->id, &this->mass,
			  &this->pos_car.x, &this->pos_car.y, &this->pos_car.z,
			  &this->vel.x, &this->vel.y, &this->vel.z,
			  &this->pot, &this->r_coll, &this->r_search);
	assert(tmp==11);
	this->pos_cyl = ConvertCar2Cyl(this->pos_car);
    }

    void writeBinary(FILE* fp) const {
	fwrite(&id, sizeof(id), 1, fp);
	fwrite(&mass, sizeof(mass), 1, fp);
	fwrite(&pos_car.x, sizeof(pos_car.x), 1, fp);
	fwrite(&pos_car.y, sizeof(pos_car.y), 1, fp);
	fwrite(&pos_car.z, sizeof(pos_car.z), 1, fp);
	fwrite(&vel.x, sizeof(vel.x), 1, fp);
	fwrite(&vel.y, sizeof(vel.y), 1, fp);
	fwrite(&vel.z, sizeof(vel.z), 1, fp);
	fwrite(&pot, sizeof(pot), 1, fp);
	fwrite(&r_coll, sizeof(r_coll), 1, fp);
	fwrite(&r_search, sizeof(r_search), 1, fp);
    }
    void readBinary(FILE* fp) {
	size_t tmp = 0;
	tmp += fread(&id, sizeof(id), 1, fp);
	tmp += fread(&mass, sizeof(mass), 1, fp);
	tmp += fread(&pos_car.x, sizeof(pos_car.x), 1, fp);
	tmp += fread(&pos_car.y, sizeof(pos_car.y), 1, fp);
	tmp += fread(&pos_car.z, sizeof(pos_car.z), 1, fp);
	tmp += fread(&vel.x, sizeof(vel.x), 1, fp);
	tmp += fread(&vel.y, sizeof(vel.y), 1, fp);
	tmp += fread(&vel.z, sizeof(vel.z), 1, fp);
	tmp += fread(&pot, sizeof(pot), 1, fp);
	tmp += fread(&r_coll, sizeof(r_coll), 1, fp);
	tmp += fread(&r_search, sizeof(r_search), 1, fp);
	assert(tmp == 11);
	this->pos_cyl = ConvertCar2Cyl(this->pos_car);
    }
    void changeMass(const PS::F64 m_new){
        const auto m_ratio = m_new / mass;
	mass = m_new;
	r_coll *= std::cbrt(m_ratio);
	r_search *= std::cbrt(m_ratio);
    }
    template<typename Tunit>
    void dump(const Tunit & unit, std::ostream & fout=std::cerr){
	const auto dens = mass / (4.0/3.0*MY_LIB::CONSTANT::pi*r_coll*r_coll*r_coll);
	fout<<"r_coll= "<<r_coll<<" ( "<<unit.getLen(r_coll)
	    <<" ),  mass= "<<mass<<" ( "<<unit.getMass(mass)
	    <<" ),  vel= "<<vel<<" ( "<<unit.getVel(vel)
	    <<" ),  dens= "<<dens<<" ( "<<unit.getDens(dens)
	    <<std::endl;
    }
};

class EPI_t{
public:
    PS::F64vec pos_car; // cartesian
    PS::F64vec pos_cyl; // cyl
    PS::F64    mass;
    PS::F64vec vel_full; // for calculating disipation energy
    PS::S64    id;
    PS::F64 r_coll;
    PS::F64 r_search;
    PS::F64vec getPos() const {
        return pos_cyl;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void setPos(const PS::F64vec & pos_new){
        pos_cyl = pos_new;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64 getRSearch() const {
        return r_search;
    }
    void copyFromFP(const FP_t & fp){ 
        pos_car  = fp.pos_car;
	pos_cyl  = fp.pos_cyl;
        mass = fp.mass;
	r_search = fp.r_search;
        vel_full  = fp.vel_full;
        id   = fp.id;
        r_coll = fp.r_coll;
    }
};

using EPJ_t = EPI_t;

class MyMomentMonopole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64ort boundary;
    MyMomentMonopole() : mass(0.0), pos(PS::F64vec(0.0)), pos_car(PS::F64vec(0.0)), boundary(PS::F64ort(PS::F64vec(-100.0), PS::F64vec(100.0))){}
    MyMomentMonopole(const PS::F64 m, const PS::F64vec & p, const PS::F64vec & p_car, const PS::F64ort & b) : mass(m), pos(p), pos_car(p_car), boundary(b){}
    void init(){
        mass = 0.0;
        pos = 0.0;
        pos_car = 0.0;
        boundary.init();
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        mass += epj.getCharge();
        pos += epj.getCharge() * epj.getPos();
        pos_car += epj.getCharge() * epj.getPosCar();
        boundary.merge(epj.getPos());
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){}
    void set(){
        pos = pos / mass;
        pos_car = pos_car / mass;
    }
    void accumulate(const MyMomentMonopole & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        pos_car += mom.mass * mom.pos_car;
        boundary.merge(mom.boundary);
    }
    void accumulate2(const MyMomentMonopole & mom){}
    // for DEBUG 
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"pos_car="<<pos_car<<std::endl;
    }
};

class MySPJMonopole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64ort boundary;
    template<class Tmom>
    void copyFromMoment(const Tmom & mom){
        this->mass     = mom.mass;
        this->pos      = mom.pos;
        this->pos_car  = mom.pos_car;
        this->boundary = mom.boundary;
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        pos_car = 0.0;
        boundary.init();
    }
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    MyMomentMonopole convertToMoment() const {
        return MyMomentMonopole(mass, pos, pos_car, boundary);
    }
};

class MyMomentQuadrupole{
public:
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64mat quad;
    PS::F64vec pos_car;
    void init(){
        pos = 0.0;
        mass = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
    MyMomentQuadrupole(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
    MyMomentQuadrupole(const PS::F64 m, const PS::F64vec & p, const PS::F64mat & q, const PS::F64vec & p_car){
        mass = m;
        pos = p;
        quad = q;
        pos_car = p_car;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        mass += epj.getCharge();
        pos  += epj.getCharge() * epj.getPos();
        pos_car += epj.getCharge() * epj.getPosCar();
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){
        PS::F64 ctmp = epj.getCharge();
        PS::F64vec ptmp = epj.getPosCar() - this->pos_car;
        PS::F64 cx = ctmp * ptmp.x;
        PS::F64 cy = ctmp * ptmp.y;
        PS::F64 cz = ctmp * ptmp.z;
        this->quad.xx += cx * ptmp.x;
        this->quad.yy += cy * ptmp.y;
        this->quad.zz += cz * ptmp.z;
        this->quad.xy += cx * ptmp.y;
        this->quad.xz += cx * ptmp.z;
        this->quad.yz += cy * ptmp.z;
    }
    void set(){
        pos = pos / mass;
        pos_car = pos_car / mass;
    }
    void accumulate(const MyMomentQuadrupole & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        pos_car += mom.mass * mom.pos_car;
    }
    void accumulate2(const MyMomentQuadrupole & mom){
        PS::F64 mtmp = mom.mass;
        PS::F64vec ptmp = mom.pos_car - this->pos_car;
        PS::F64 cx = mtmp * ptmp.x;
        PS::F64 cy = mtmp * ptmp.y;
        PS::F64 cz = mtmp * ptmp.z;
        this->quad.xx += cx * ptmp.x + mom.quad.xx;
        this->quad.yy += cy * ptmp.y + mom.quad.yy;
        this->quad.zz += cz * ptmp.z + mom.quad.zz;
        this->quad.xy += cx * ptmp.y + mom.quad.xy;
        this->quad.xz += cx * ptmp.z + mom.quad.xz;
        this->quad.yz += cy * ptmp.z + mom.quad.yz;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass= "<<mass<<std::endl;
        fout<<"pos= "<<pos<<std::endl;
        fout<<"quad= "<<quad<<std::endl;
    }
};

class MySPJQuadrupole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64mat quad;
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void copyFromMoment(const MyMomentQuadrupole & mom){
        mass = mom.mass;
        pos = mom.pos;
        quad = mom.quad;
        pos_car = mom.pos_car;
    }
    MyMomentQuadrupole convertToMoment() const {
        return MyMomentQuadrupole(mass, pos, quad, pos_car);
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
};


template<typename Tpi, typename Tpj, typename Tforce>
void CalcForceFromArray(const Tpi & pi, const Tpj * pj, const PS::S32 nj, Tforce & force){
    const auto eps2  = FP_t::eps*FP_t::eps;
    const auto kappa = FP_t::kappa;
    const auto eta   = FP_t::eta;
    PS::F64vec xj[nj];
    for(auto j=0; j<nj; j++){
	xj[j] = pj[j].getPosCar();
    }
    PS::F64vec ai_grav = 0.0;
    PS::F64vec ai_sprg = 0.0;
    PS::F64vec ai_dash = 0.0;
    PS::F64 poti_grav = 0.0;
    PS::F64 poti_sprg = 0.0;    
    for(auto j=0; j<nj; j++){
	if(pi.id == pj[j].id) continue;
	const auto r_coll    = (pi.r_coll + pj[j].r_coll);
	const auto r_coll_sq = r_coll * r_coll;
	const auto rij       = pi.getPosCar() - xj[j];
	const auto r_real_sq = rij * rij + eps2;
	const auto over_r_real = 1.0 / sqrt(r_real_sq);
	const auto over_r_real_sq = over_r_real * over_r_real;
	if(r_coll_sq > r_real_sq){
	    const auto r_coll_cu  = r_coll_sq * r_coll;
	    const auto over_r_coll_cu = 1.0 / r_coll_cu;
	    ai_grav -= pj[j].getCharge() * over_r_coll_cu * rij;
	    const auto pot_offset = -1.5 / r_coll;
	    poti_grav += 0.25 * pj[j].mass * (r_real_sq * over_r_coll_cu + pot_offset);

	    const auto m_red = pj[j].mass / (pi.mass + pj[j].mass);
	    const auto r_real = r_real_sq * over_r_real;
	    const auto dr = r_coll - r_real;
	    ai_sprg += kappa * m_red * dr * over_r_real * rij;
	    poti_sprg += 0.25 * kappa * m_red * dr * dr;

	    const auto vij = pi.vel_full - pj[j].vel_full;
	    const auto rv  = rij * vij;
	    ai_dash -= eta * m_red * rv * over_r_real_sq * rij;
	}
	else{
	    const auto m_over_r_real = pj[j].mass * over_r_real;
	    ai_grav   -= m_over_r_real * over_r_real_sq * rij;
	    poti_grav -= 0.5 * m_over_r_real;
	}
    }
    force.acc      += ai_grav + ai_sprg + ai_dash;
    force.acc_dash += ai_dash;
    force.pot      += poti_grav + poti_sprg;    
}

#if defined(USE_PIKG_KERNEL)

#include"user_defined_kernel.hpp"

#else


template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceEp{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
#if 1
	for(auto i=0; i<ni; i++){
	    CalcForceFromArray(pi[i], pj, nj, force[i]);
	}
#else
        const auto eps2  = FP_t::eps * FP_t::eps;
        const auto kappa = FP_t::kappa;
        const auto eta   = FP_t::eta;
	PS::F64vec xj[nj];
	for(auto j=0; j<nj; j++){
	  xj[j] = pj[j].getPosCar();
	}
        for(auto i=0; i<ni; i++){
            const PS::F64vec xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64vec ai_dash = 0.0;
            PS::F64 poti = 0.0;
            for(auto j=0; j<nj; j++){
                const auto r_coll    = (pi[i].r_coll + pj[j].r_coll);
                const auto r_coll_sq = r_coll * r_coll;
		PS::F64vec rij       = xi - xj[j];
                if(pi[i].id == pj[j].id) continue;
                PS::F64 r2_real = rij * rij + eps2;
                PS::F64 r2      = std::max(r2_real, r_coll_sq);
                PS::F64 r_inv   = 1.0/sqrt(r2);
                PS::F64 r2_inv  = r_inv * r_inv;
                PS::F64 pot = r_inv * pj[j].getCharge();
                if(r_coll_sq > r2_real){
                    ai     -= pj[j].getCharge() / (r_coll_sq*r_coll) * rij;
                    PS::F64 pot_offset = -1.5/r_coll;
                    poti   += 0.5*pj[j].getCharge()*(0.5*r2_real/(r_coll_sq*r_coll) + pot_offset);
                    PS::F64 m_r = pj[j].mass / (pi[i].mass+pj[j].mass);
                    PS::F64 r   = sqrt(r2_real);
                    PS::F64 dr  = r_coll-r ;
                    ai += kappa * m_r * dr/r * rij;
                    poti += 0.5*kappa*m_r*dr*dr * 0.5;
                    PS::F64vec vij = pi[i].vel_full - pj[j].vel_full;
                    PS::F64 rv = rij*vij;
                    PS::F64vec a_eta = eta * m_r * rv / r2_real * rij;
                    ai_dash += a_eta;
                    ai += a_eta;
                }
                else{
                    ai     -= pot * r2_inv * rij;
                    poti   -= 0.5 * pot;
                }
            }
            force[i].acc += ai;
            force[i].acc_dash += ai_dash;
            force[i].pot += poti;
        }
#endif
    }
};

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpMono{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2 = FP_t::eps*FP_t::eps;
	PS::F64vec xj[nj];
	for(auto j=0; j<nj; j++){
	  xj[j] = pj[j].getPosCar();
	}
        for(auto i=0; i<ni; i++){
            const auto xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64 poti  = 0.0;
            for(auto j=0; j<nj; j++){
	      //PS::F64vec rij    = xi - pj[j].getPosCar();
	        const auto rij    = xi - xj[j];
                auto r3_inv = rij * rij + eps2;
                auto r_inv  = 1.0/sqrt(r3_inv);
                r3_inv  = r_inv * r_inv;
                r_inv  *= pj[j].getCharge();
                r3_inv *= r_inv;
                ai     -= r3_inv * rij;
                poti   -= 0.5*r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpQuad{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2 = FP_t::eps*FP_t::eps;
	PS::F64vec xj[nj];
	for(auto j=0; j<nj; j++){
	  xj[j] = pj[j].getPosCar();
	}
        for(auto i=0; i<ni; i++){
            PS::F64vec xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(auto j=0; j<nj; j++){
                PS::F64 mj = pj[j].getCharge();
                //PS::F64vec xj= pj[j].getPosCar();
                //PS::F64vec rij= xi - xj;
		PS::F64vec rij = xi - xj[j];
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = pj[j].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5.0*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[i].acc += ai;
            force[i].pot += 0.5*poti;
        }
    }
};

#endif


class Energy{
public:
    PS::F64 tot;
    PS::F64 pot;
    PS::F64 kin;
    PS::F64 disp; // glb
    Energy(){
        pot = kin = tot = disp = 0.0;
    }
    void setEngDispGlb(const PS::F64 e){
	disp = e;
    }
    template<typename Tpsys>
    void calc(const Tpsys & psys){
        const auto n = psys.getNumberOfParticleLocal();
        PS::F64 pot_loc = 0.0;
        PS::F64 kin_loc = 0.0;
        for(auto i=0; i<n; i++){
            kin_loc += 0.5*psys[i].mass*psys[i].vel*psys[i].vel;
            pot_loc += psys[i].mass*psys[i].pot;
        }
        pot = PS::Comm::getSum(pot_loc);
        kin = PS::Comm::getSum(kin_loc);
        tot = pot + kin;
    }
    void writeAscii(FILE * fp){
	fprintf(fp, "%12.11e   %12.11e   %12.11e   %12.11e \n", tot, pot, kin, disp);
    }
    void readAscii(FILE * fp){
	auto tmp = fscanf(fp, "%lf   %lf   %lf   %lf", &tot, &pot, &kin, &disp);
	assert(tmp == 4);
    }

    void writeBinary(FILE * fp){
	fwrite(this, sizeof(*this), 1, fp);
    }
    void readBinary(FILE * fp){
	auto tmp = fread(this, sizeof(*this), 1, fp);
	assert(tmp == 1);
    }
    
    PS::F64 getKin() const {
	return kin;
    }
    PS::F64 getPot() const {
	return pot;
    }
    PS::F64 getDisp() const {
	return disp;
    }
    void dump(std::ostream & fout = std::cerr){
        fout<<"tot= "<<tot<<" pot= "<<pot<<" kin= "<<kin<<std::endl;
    }
};

class SatelliteSystem{
    std::vector<FP_t> sat;
public:
    SatelliteSystem(){
	sat.clear();
    }
    PS::S32 getNumberOfSat() const { return sat.size();}
    
    const SatelliteSystem & operator = (const SatelliteSystem & s){
	const PS::S32 n_sat = s.getNumberOfSat();
	this->clear();
	for(auto i=0; i<n_sat; i++){
	    this->setSatellite(s[i]);
	}
	return (*this);
    }
    
    void setSatellite(const FP_t & s){
	sat.push_back(s);
    }
    void clear(){
	sat.clear();
    }
    void copy(const SatelliteSystem & s){
	const PS::S32 n_sat = s.getNumberOfSat();
	this->clear();
	for(auto i=0; i<n_sat; i++){
	    this->setSatellite(s[i]);
	}
    }
    
    template<typename Tpsys>
    void CalcForce(Tpsys & psys){
	const auto n = psys.getNumberOfParticleLocal();
	const auto n_sat = sat.size();
	Force_t f_sat_loc[n_sat];
	Force_t f_sat_glb[n_sat];
#pragma omp parallel for
	for(size_t i=0; i<n_sat; i++){
	    f_sat_loc[i].clear();
	    f_sat_glb[i].clear();
	    sat[i].clearForce();
	    CalcForceFromArray(sat[i], &psys[0], n, f_sat_loc[i]); // ptcl  -> satellite 
	}
	assert(sizeof(Force_t) % sizeof(PS::F64) == 0);
	//MPI_Allreduce(f_sat_loc, f_sat_glb, (sizeof(Force_t)/sizeof(PS::F64))*n_sat, MPI_DOUBLE, MPI_SUM, PS::Comm::getCommunicator());
	SumForce(f_sat_loc, f_sat_glb, n_sat);
#pragma omp parallel for
	for(auto i=0; i<n; i++){
	    CalcForceFromArray(psys[i], &sat[0], n_sat, psys[i]); // satellite -> ptcl
	}

#pragma omp parallel for
	for(size_t i=0; i<n_sat; i++){
	    CalcForceFromArray(sat[i], &sat[0], n_sat, sat[i]); // satellite -> satellite
	}

#pragma omp parallel for
	//for(auto i=0; i<n_sat; i++){
	for(size_t i=0; i<n_sat; i++){
	    sat[i].accumFromForce(f_sat_glb[i]);
	}
    }

    void writeAscii(FILE* fp){
	PS::S64 n_sat = (PS::S64)sat.size();
        fprintf(fp, "%lld\n", n_sat);
	for(auto i=0; i<n_sat; i++){
	    sat[i].writeAscii(fp);
	}
    }
    void readAscii(FILE * fp){
	PS::S64 n_sat = (PS::S64)sat.size();
        auto tmp = fscanf(fp, "%lld", &n_sat);
	assert(tmp == 1);
	sat.resize(n_sat);
	for(auto i=0; i<n_sat; i++){
	    sat[i].readAscii(fp);
	}
    }
    void writeBinary(FILE * fp){
	PS::S64 n_sat = (PS::S64)sat.size();
        fwrite(&n_sat, sizeof(n_sat), 1, fp);
	for(auto i=0; i<n_sat; i++){
	    sat[i].writeBinary(fp);
	}
    }    
    void readBinary(FILE * fp){
	PS::S64 n_sat = (PS::S64)sat.size();
	auto tmp = fread(&n_sat, sizeof(n_sat), 1, fp);
	assert(tmp == 1);
	sat.resize(n_sat);
	for(auto i=0; i<n_sat; i++){
	    sat[i].readBinary(fp);
	}
    }


    FP_t & operator [] (const PS::S32 id) {return sat[id];}
    const FP_t & operator [] (const PS::S32 id) const {return sat[id];}
    
    void allgatherSystem(SatelliteSystem & sys){
	const auto n_proc = PS::Comm::getNumberOfProc();
	PS::S32 n_sats[n_proc];
	PS::S32 n_sat  = sys.getNumberOfSat();
	PS::Comm::allGather(&n_sat, 1, n_sats);
	/*
	if(PS::Comm::getRank()==1){
	    for(auto i=0; i<n_proc; i++){
		std::cerr<<"n_sats[i]= "<<n_sats[i]<<std::endl;
	    }
	}
	*/
	PS::S32 n_sat_glb = 0;
	PS::S32 n_disp_sats[n_proc+1];
	n_disp_sats[0] = 0;
	for(auto i=0; i<n_proc; i++){
	    n_sat_glb += n_sats[i];
	    n_disp_sats[i+1] = n_disp_sats[i] + n_sats[i];
	}
	/*
	if(PS::Comm::getRank()==1){
	    std::cerr<<"n_sat_glb= "<<n_sat_glb<<std::endl;
	}
	*/
	FP_t sat_send[n_sat];
	for(auto i=0; i<n_sat; i++){
	    sat_send[i] = sys.sat[i];
	}
	sat.resize(n_sat_glb);
	if(n_sat_glb > 0){
	    PS::Comm::allGatherV(sat_send, n_sat, &sat[0], n_sats, n_disp_sats);
	}
    }

    void broadcast(){
	const PS::S32 n_sat = sat.size();
	if(n_sat > 0){
	    PS::Comm::broadcast(&sat[0], n_sat);
	}
    }
    
    void setId(){
	const PS::S64 n_sat  = sat.size();
	for(auto i=0; i<n_sat; i++){
	    sat[i].id = i;
	}
    }
};


template<typename Tpsys>
void AddSatToSys(Tpsys & psys, const SatelliteSystem & sat){
    const PS::S32 n_sat = sat.getNumberOfSat();
    for(auto i=0; i<n_sat; i++){
	psys.addOneParticle(sat[i]);
    }    
    /*
    if(PS::Comm::getRank()==0){
	const PS::S32 n_sat = sat.getNumberOfSat();
	for(auto i=0; i<n_sat; i++){
	    psys.addOneParticle(sat[i]);
	}
    }
    */
}

template<typename Tpsys>
void RemoveSatFromSys(Tpsys & psys, SatelliteSystem & sat){
    PS::S32 n_remove = sat.getNumberOfSat();
    PS::S32 n_ptcl = psys.getNumberOfParticleLocal() - n_remove;
    PS::S32 id_remove[n_remove];
    for(auto i=0; i<n_remove; i++){
	sat[i] = psys[n_ptcl+i];
	id_remove[i] = n_ptcl+i;
    }
    psys.removeParticle(id_remove, n_remove);    
    /*
    if(PS::Comm::getRank()==0){
	PS::S32 n_remove = sat.getNumberOfSat();
	PS::S32 n_ptcl = psys.getNumberOfParticleLocal() - n_remove;
	PS::S32 id_remove[n_remove];
	for(auto i=0; i<n_remove; i++){
	    sat[i] = psys[n_ptcl+i];
	    id_remove[i] = n_ptcl+i;
	}
	psys.removeParticle(id_remove, n_remove);
    }
    if(PS::Comm::getRank()==1){
	std::cerr<<"A) sat[0].acc= "<<sat[0].acc<<" pot= "<<sat[0].pot<<std::endl;
    }    
    sat.broadcast();
    if(PS::Comm::getRank()==1){
	std::cerr<<"B) sat[0].acc= "<<sat[0].acc<<" pot= "<<sat[0].pot<<std::endl;
    }
    */
}


class Param{
    PS::F64 time_sys;
    PS::F64 dt;
public:
    Param(){}
    
    void integrateTime(){
	time_sys += dt;
    }
    PS::F64 getDt() const {
	return dt;
    }
    PS::F64 getTimeSys() const {
	return time_sys;
    }
    
    void setDt(const PS::F64 t){
	dt = t;
    }
    void setTimeSys(const PS::F64 t){
	time_sys = t;
    }    
};

