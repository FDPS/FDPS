#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>

///////////////
/// GRAVITY ///
class ForceGrav{
public:
    PS::F64vec acc;
    PS::F64 pot;
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
};

class FPGrav{
public:
    PS::F64vec pos;
    PS::F64vec getPos() const {
        return this->pos;
    }
    void copyFromForce(const ForceGrav & force){
        this->acc = force.acc;
        this->pot = force.pot;
    }
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
        fout<<"acc="<<acc<<std::endl;
        fout<<"pot="<<pot<<std::endl;
    }
};

class EPIGrav{
public:
    PS::S64 id;
    PS::F64vec pos;
    static PS::F64 eps;
    PS::F64vec getPos() const {
        return this->pos;
    }
    void copyFromFP(const FPGrav & fp){ 
        pos = fp.pos;
        id = fp.id;
    }
};

PS::F64 EPIGrav::eps = 1.0/32.0;

class EPJGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    void copyFromFP(const FPGrav & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

class MomentGrav{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64mat quad;

    MomentGrav(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
    }

    MomentGrav(const PS::F64 m, const PS::F64vec & p, const PS::F64mat & q){
        mass = m;
        pos = p;
        quad = q;
    }

    void init(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    void accumulateAtLeaf(const EPJGrav & epj){
        mass += epj.mass;
        pos += epj.mass * epj.pos;
    }
    void accumulateAtLeaf2(const EPJGrav & epj){
        PS::F64 ctmp = epj.mass;
        PS::F64vec ptmp = epj.pos - this->pos;
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
    }
    void accumulate(const MomentGrav & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
    }
    void accumulate2(const MomentGrav & mom){
        PS::F64 mtmp = mom.mass;
        PS::F64vec ptmp = mom.pos - this->pos;
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
    // for DEBUG 
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"quad="<<quad<<std::endl;
    }
};

class SPJ{
public:
    void copyFromMoment(const MomentGrav & mom){
        mass = mom.mass;
        pos = mom.pos;
        quad = mom.quad;
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    MomentGrav convertToMoment() const {
        return MomentGrav(mass, pos, quad);
    }

    PS::F64 mass;
    PS::F64vec pos;
    PS::F64mat quad;
};

#if 0
void CalcForceEpEp(const EPIGrav * ep_i,
                   const PS::S32 n_ip,
                   const EPJGrav * ep_j,
                   const PS::S32 n_jp,
                   ForceGrav * force){
    PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
    for(PS::S32 i=0; i<n_ip; i++){
        PS::F64vec xi = ep_i[i].pos;
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        PS::S64 idi = ep_i[i].id;
        for(PS::S32 j=0; j<n_jp; j++){
            if( idi == ep_j[j].id ) continue;
            PS::F64vec rij = xi - ep_j[j].pos;
            PS::F64 r3_inv = rij * rij + eps2;
            PS::F64 r_inv = 1.0/sqrt(r3_inv);
            r3_inv = r_inv * r_inv;
            r_inv *= ep_j[j].mass;
            r3_inv *= r_inv;
            ai -= r3_inv * rij;
            poti -= r_inv;
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}
#else
struct CalcForceEpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const EPJGrav * ep_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S64 idi = ep_i[i].id;
            for(PS::S32 j=0; j<n_jp; j++){
                if( idi == ep_j[j].id ) continue;
                PS::F64vec rij = xi - ep_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#endif

#if 0
void CalcForceSpEp(const EPIGrav * ep_i,
                   const PS::S32 n_ip,
                   const SPJ * sp_j,
                   const PS::S32 n_jp,
                   ForceGrav * force){
    PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
    for(PS::S32 ip=0; ip<n_ip; ip++){
        PS::F64vec xi = ep_i[ip].pos;
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 jp=0; jp<n_jp; jp++){
            PS::F64 mj = sp_j[jp].mass;
            PS::F64vec xj= sp_j[jp].pos;
            PS::F64vec rij= xi - xj;
            PS::F64 r2 = rij * rij + eps2;
            PS::F64mat qj = sp_j[jp].quad;
            PS::F64 tr = qj.getTrace();
            PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                           (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                           (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
            PS::F64 qrr = qr * rij;
            PS::F64 r_inv = 1.0f/sqrt(r2);
            PS::F64 r2_inv = r_inv * r_inv;
            PS::F64 r3_inv = r2_inv * r_inv;
            PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
            PS::F64 qrr_r5 = r5_inv * qrr;
            PS::F64 qrr_r7 = r2_inv * qrr_r5;
            PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
            PS::F64 B = -2.0*r5_inv;
            ai -= A*rij + B*qr;
            poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
        }
        force[ip].acc += ai;
        force[ip].pot += poti;
    }
}
#else
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const SPJ * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
        /*
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F64vec xi = ep_i[ip].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64 mj = sp_j[jp].mass;
                PS::F64vec xj= sp_j[jp].pos;
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
        */
    }
};
#endif

struct CompareGrav{
    void operator () (ForceGrav grav0, ForceGrav grav1){
        PS::F64 dpot = std::abs( (grav0.pot - grav1.pot) / grav0.pot );
        PS::F64vec dacc_vec = grav0.acc - grav1.acc;
        PS::F64 dacc = sqrt( (dacc_vec*dacc_vec) / (grav0.acc*grav0.acc) );
        if( dpot > 1e-3 || dacc > 1e-1){
            std::cout<<"Compare Grav: FAIL"<<std::endl;
            std::cout<<"grav0.pot="<<grav0.pot<<" grav1.pot="<<grav1.pot<<std::endl;
            std::cout<<"grav0.acc="<<grav0.acc<<" grav1.acc="<<grav1.acc<<std::endl;
        }
    }
};


///////////
/// SPH ///

class ResultDens{
public:
    void clear(){
        dens = r_search_next = kernel_length = 0.0;
        n_neighbour = 0;
    }
    PS::F64 dens;
    PS::F64 r_search_next;
    PS::S32 n_neighbour;
    PS::F64 kernel_length;
};

class FPSPH{
public:
    void copyFromForce(const ResultDens & dens){
        this->dens = dens.dens;
        this->kernel_length = dens.kernel_length;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getRSearch() const {
        return this->r_search;
    }

    PS::F64 mass;
    PS::F64vec pos;
    PS::F64 r_search;
    PS::F64 kernel_length;
    PS::F64 dens;
    PS::F64 divv;
    PS::F64 pres;
    PS::F64 vel_sound;
    PS::F64vec acc;
    PS::F64 engdot;
    PS::S64 id;
    PS::F64vec vel;
    PS::F64 eng;
    PS::F64vec vel_half;
    PS::F64 eng_half;
    PS::S64 n_neighbour;
    void dump(std::ostream & fout=std::cout) const {
        fout<<"pos="<<pos<<std::endl;
        fout<<"dens="<<dens<<std::endl;
        fout<<"kernel_length="<<kernel_length<<std::endl;
    }
};

class EPIGather{
public:
    enum{
        n_neighbour_crit = 7,
    };
    PS::F64vec pos;
    PS::F64 r_search;
    PS::S64 id;
    PS::S64 getIndex() const { return this->id;}
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getRSearch() const {
        return this->r_search;
    }
    void copyFromFP(const FPSPH & fp){
        this->pos = fp.getPos();
        this->r_search = fp.getRSearch();
        id = fp.id;
    }
};

class EPJGather{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::S64 id;
    void copyFromFP(const FPSPH & fp){
        this->mass = fp.getCharge();
        this->pos = fp.getPos();
        this->id = fp.id;
    }
};

class EPIScatter{
public:
    PS::S64 id;
    PS::F64vec pos;
    void copyFromFP(const FPSPH & fp){
        this->pos = fp.getPos();
        id = fp.id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
};

class EPJScatter{
public:
    enum{
        n_neighbour_crit = 7,
    };
    PS::F64 mass;
    PS::F64 r_search;
    PS::F64vec pos;
    PS::S64 id;
    PS::S64 getIndex() const { return this->id; }
    PS::F64 getCharge() const { return this->mass; }
    PS::F64vec getPos() const { return this->pos; }
    PS::F64 getRSearch() const { return this->r_search; }
    void copyFromFP(const FPSPH & fp){
        this->mass = fp.getCharge();
        this->r_search = fp.getRSearch();
        this->pos = fp.getPos();
        this->id = fp.id;
    }
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"r_search="<<r_search<<std::endl;
    }
};

PS::F64 CubicSpline(const PS::F64 r_sq,
                    const PS::F64 h_inv){
    PS::F64 xi = sqrt(r_sq) * h_inv;
    PS::F64 xi10 = (1.0-xi > 0.0) ? 1.0-xi : 0.0;
    PS::F64 xi05 = (0.5-xi > 0.0) ? 0.5-xi : 0.0;
    return xi10*xi10*xi10 - 4.0*xi05*xi05*xi05;
}

PS::F64vec CubicSpline_ri(const PS::F64vec & rij,
                          const PS::F64 h_inv){
    PS::F64 r_sq = rij * rij;
    PS::F64 r = sqrt(r_sq);
    PS::F64 xi = r * h_inv;
    PS::F64 xi10 = (1.0-xi > 0.0) ? 1.0-xi : 0.0;
    PS::F64 xi05 = (0.5-xi > 0.0) ? 0.5-xi : 0.0;
    PS::F64 C = (3.0*xi10*xi10 - 12.0*xi05*xi05) * h_inv / r;
    return -C * rij;
}


struct CalcDensityEpEp{
    void operator () (const EPIScatter * ep_i,
                      const PS::S32 n_ip,
                      const EPJScatter * ep_j,
                      const PS::S32 n_jp,
                      ResultDens * dens){
        static PS::F64 Cnorm = 8.0/3.0; // 1D
        //static PS::F64 Cnorm = 80.0/(7.0*M_PI); // 2D
        //static PS::F64 Cnorm = 16.0/M_PI; // 3D
        for(PS::S32 i=0; i<n_ip; i++){
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64 r_crit_sq = ep_j[j].getRSearch() * ep_j[j].getRSearch();
                const PS::F64vec dr = ep_j[j].getPos() - ep_i[i].getPos();
                const PS::F64 dr_sq = dr * dr;
                if( r_crit_sq > dr_sq ){
                    const PS::F64 h_inv = 1.0 / ep_j[j].getRSearch();
                    dens[i].dens += ep_j[j].getCharge() * CubicSpline(dr_sq, h_inv);
                }
            }
            dens[i].dens *= Cnorm;
        }
    }
};

struct CompareDensity{
    void operator () (ResultDens dens0, ResultDens dens1){
        if( std::abs(dens0.dens - dens1.dens) > 1e-5){
            std::cout<<"FAIL"<<std::endl;
            std::cout<<"desn0.dens="<<dens0.dens<<" dens1.dens="<<dens1.dens<<std::endl;
        }
    }
};

template<class Tpsys>
void ReadNemoAscii(Tpsys & psys,
                   PS::S32 & n_glb,
                   PS::S32 & n_loc,  
                   PS::F32 & t_sys,
                   const char * ifile){
    std::ifstream finput;
    finput.open(ifile);
    assert(finput);
    PS::S32 dim;
    finput>>n_glb>>dim>>t_sys;
    std::cerr<<"ifile:"<<ifile<<std::endl;
    std::cerr<<"n_glb="<<n_glb<<std::endl;
    std::cerr<<"dim="<<dim<<std::endl;
    std::cerr<<"t_sys="<<t_sys<<std::endl;

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb/n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;

    psys.createParticle((n_glb/n_proc)*4+1000);
    psys.setNumberOfParticleLocal(n_loc);

    PS::F32vec pos_shift(1.0, 2.0, 3.0);

    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    const PS::S32 i_t = i_h+n_loc;
    PS::F32 xf32;
    PS::F32vec vf32;

    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++)psys[n].id = i;

    for(PS::S32 i=0; i<i_h; i++) finput>>xf32;
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++)finput>>psys[n].mass;
    for(PS::S32 i=i_t; i<n_glb; i++) finput>>xf32;

    for(PS::S32 i=0; i<i_h; i++) finput>>vf32;
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++){
        finput>>psys[n].pos;
        psys[n].pos += pos_shift;
    }
    for(PS::S32 i=i_t; i<n_glb; i++) finput>>vf32;

    for(PS::S32 i=0; i<i_h; i++) finput>>vf32;
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++) finput>>psys[n].vel;
    for(PS::S32 i=i_t; i<n_glb; i++) finput>>vf32;

    finput.close();

}

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);

    char sinput[1024];
    int c;
    while((c=getopt(argc,argv,"i:h")) != -1){
        switch(c){
        case 'i':
            sprintf(sinput,optarg);
            break;
        case 'h':
            std::cerr<<"i: input file name (nemo ascii)"<<std::endl;
            return 0;
        }
    }


    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    PS::S32 n_grav_glb, n_grav_loc;
    PS::F32 time_sys;
    ReadNemoAscii(system_grav, n_grav_glb, n_grav_loc, time_sys, sinput);

    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomain();

    system_grav.exchangeParticle(dinfo);
    n_grav_loc = system_grav.getNumberOfParticleLocal();

    PS::F64 half_len_grav_glb = system_grav.getHalfLength();
    std::cout<<"half_len_grav_glb="<<half_len_grav_glb<<std::endl;

#if 1
    PS::TreeForForce<PS::SEARCH_MODE_LONG, ForceGrav, EPIGrav, EPJGrav, MomentGrav, SPJ> tree_grav;
    PS::F32 theta = 0.4;
    tree_grav.initialize(n_grav_glb, theta);
    tree_grav.initializeLocalTree(half_len_grav_glb);
    tree_grav.setParticleLocalTree(system_grav);

    tree_grav.mortonSortLocalTreeOnly();
    //tree_grav.checkMortonSortLocalTreeOnly();

    tree_grav.linkCellLocalTreeOnly();
    //tree_grav.checkMakeLocalTree();

    tree_grav.calcMomentLocalTreeOnly();
    //tree_grav.checkCalcMomentLocalTree();

    tree_grav.exchangeLocalEssentialTree(dinfo);
    //tree_grav.checkExchangeLocalEssentialTree(dinfo);

    tree_grav.setLocalEssentialTreeToGlobalTree();

    tree_grav.mortonSortGlobalTreeOnly();
    tree_grav.checkMortonSortGlobalTreeOnly();

    tree_grav.linkCellGlobalTreeOnly();
    tree_grav.checkMakeGlobalTree();

    tree_grav.calcMomentGlobalTreeOnly();
    tree_grav.checkCalcMomentGlobalTree();

    tree_grav.makeIPGroup();
    tree_grav.checkMakeIPGroup();

    PS::S32 n_ipg_grav = tree_grav.getNumberOfIPG();
    bool err_grav = false;
    for(PS::S32 i=0; i<n_ipg_grav; i++){
        tree_grav.makeInteractionList(i, err_grav);
        //tree_grav.checkMakeInteractionList();
        tree_grav.calcForceOnly(CalcForceEpEp(), CalcForceSpEp(), i);
    }
    tree_grav.copyForceOriginalOrder();
    ForceGrav * force_direct = new ForceGrav[n_grav_loc];
    tree_grav.checkForce( CalcForceEpEp(), force_direct, CompareGrav());

    /*
    //tree_grav.calcForceDirect(CalcForceEpEp(), force_direct);
    tree_grav.checkForce(CalcForceEpEp(), force_direct);
    PS::F64 Epot_loc_direct = 0.0;
    for(PS::S32 i=0; i<n_grav_loc; i++){
        Epot_loc_direct += 0.5 * system_grav[i].mass * force_direct[i].pot;
    }
    PS::F64 Epot_glb_direct = 0.0;
    MPI::COMM_WORLD.Allreduce(&Epot_loc_direct, &Epot_glb_direct, 1, PS::GetDataType<PS::F64>(), MPI::SUM);

    PS::F64 Epot_loc = 0.0;
    for(PS::S32 i=0; i<n_grav_loc; i++){
        Epot_loc += 0.5 * system_grav[i].mass * tree_grav.getForce(i).pot;
    }
    PS::F64 Epot_glb = 0.0;
    MPI::COMM_WORLD.Allreduce(&Epot_loc, &Epot_glb, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
    std::cout<<"Epot_loc="<<Epot_loc<<" Epot_loc_direct="<<Epot_loc_direct<<std::endl;
    std::cout<<"Epot_glb="<<Epot_glb<<" Epot_glb_direct="<<Epot_glb_direct<<std::endl;
    */
#endif

#if 1
    ///////////////
    //// SHORT ////
    PS::ParticleSystem<FPSPH> system_sph;
    system_sph.initialize();
    PS::S32 n_sph_glb, n_sph_loc;
    ReadNemoAscii(system_sph, n_sph_glb, n_sph_loc, time_sys, sinput);
    PS::F64 half_len_sph_glb = system_sph.getHalfLength();
    std::cout<<"half_len_sph_glb="<<half_len_sph_glb<<std::endl;

    for(PS::S32 i=0; i<n_sph_loc; i++){
        system_sph[i].r_search = pow( (n_sph_glb/(half_len_sph_glb*half_len_sph_glb*half_len_sph_glb)), -0.333333) * (system_sph[i].id%10)*0.2;
    }

    system_sph.exchangeParticle(dinfo);
    n_sph_loc = system_sph.getNumberOfParticleLocal();

    //////////////////////
    //// SCATTER MODE ////
    PS::TreeForForce<PS::SEARCH_MODE_SCATTER, ResultDens, EPIScatter, EPJScatter, PS::MomentSearch, PS::SuperParticleBase> tree_scatter;
    tree_scatter.initialize(n_sph_glb);
    tree_scatter.initializeLocalTree(half_len_sph_glb);
    tree_scatter.setParticleLocalTree(system_sph);

    tree_scatter.mortonSortLocalTreeOnly();
    tree_scatter.checkMortonSortLocalTreeOnly();

    tree_scatter.linkCellLocalTreeOnly();
    tree_scatter.checkMakeLocalTree();

    tree_scatter.calcMomentLocalTreeOnly();
    tree_scatter.checkCalcMomentLocalTree();

    tree_scatter.exchangeLocalEssentialTree(dinfo);
    tree_scatter.checkExchangeLocalEssentialTree(dinfo);


    tree_scatter.setLocalEssentialTreeToGlobalTree();

    tree_scatter.mortonSortGlobalTreeOnly();
    tree_scatter.checkMortonSortGlobalTreeOnly();

    tree_scatter.linkCellGlobalTreeOnly();
    tree_scatter.checkMakeGlobalTree();

    tree_scatter.calcMomentGlobalTreeOnly();
    tree_scatter.checkCalcMomentGlobalTree();

    tree_scatter.makeIPGroup();
    tree_scatter.checkMakeIPGroup();


    PS::S32 n_ipg_sph = tree_scatter.getNumberOfIPG();
    bool err_sph = false;
    for(PS::S32 i=0; i<n_ipg_sph; i++){
        tree_scatter.makeInteractionList(i, err_sph);
        //tree_scatter.checkMakeInteractionList(i);
        tree_scatter.calcForceOnly( CalcDensityEpEp(), i);
    }
    tree_scatter.copyForceOriginalOrder();

    ResultDens * dens_direct = new ResultDens[n_sph_loc];
    tree_scatter.checkForce( CalcDensityEpEp(), dens_direct, CompareDensity());
    /*
    tree_scatter.calcForceDirect( CalcDensityEpEp(), dens_direct );
    for(PS::S32 i=0; i<n_sph_loc; i++){
        if( std::abs(tree_scatter.getForce(i).dens - dens_direct[i].dens) > 1e-5){
            std::cout<<"FAIL"<<std::endl;
            std::cout<<"tree_scatter.getForce(i).dens="<<tree_scatter.getForce(i).dens<<
                " dens_direct[i].dens="<<dens_direct[i].dens<<std::endl;
        }
    }
    */

#endif
#if 0
    /////////////////////
    //// GATHER MODE ////
    PS::TreeForForce<PS::SEARCH_MODE_GATHER, ResultDens, EPIGather, EPJGather, PS::MomentSearch, PS::SuperParticleBase> tree_gather;
    tree_gather.initialize(n_sph_glb);
    tree_gather.initializeLocalTree(half_len_sph_glb);
    tree_gather.setParticleLocalTree(system_sph);

    tree_gather.mortonSortLocalTreeOnly();
    tree_scatter.checkMortonSortLocalTreeOnly();

    tree_gather.linkCellLocalTreeOnly();
    tree_gather.checkMakeLocalTree();

    tree_gather.calcMomentLocalTreeOnly();
    tree_gather.checkCalcMomentLocalTree();
#endif

    PS::Finalize();
    return 0;

}
