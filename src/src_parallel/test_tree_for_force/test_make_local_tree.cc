#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>

class Force{
public:
    PS::F64vec acc;
    PS::F64 pot;
};

class FP{
public:
    PS::F64vec pos;
    PS::F64vec getPos() const {
        return this->pos;
    }
    void copyFromForce(const Force & force){
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


class EPI{
public:
    PS::F64vec pos;
    PS::S64 id;
    static PS::F64 eps;
    void copyFromFP(const FP & fp){ 
        pos = fp.pos;
        id = fp.id;
    }
};

PS::F64 EPI::eps = 1.0/32.0;

class EPJ{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    void copyFromFP(const FP & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

class Moment{
public:
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64mat quad;
    void init(){
        pos = 0.0;
        mass = 0.0;
        quad = 0.0;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    void accumulateAtLeaf(const EPJ & epj){
        mass += epj.mass;
        pos += epj.mass * epj.pos;
    }
    void accumulateAtLeaf2(const EPJ & epj){
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
    void accumulate(const Moment & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
    }
    void accumulate2(const Moment & mom){
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
    void copyFromMoment(const Moment & mom){
        mass = mom.mass;
        pos = mom.pos;
        quad = mom.quad;
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
    }
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64mat quad;
};

int main(int argc, char *argv[]){
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
    std::ifstream finput;
    finput.open(sinput);
    PS::S32 n_ptcl, dim;
    PS::F64 time_sys;
    finput>>n_ptcl>>dim>>time_sys;

    const PS::S32 n_ptcl_copy = 20;
    //const PS::S32 n_ptcl_copy = 0;
    PS::ParticleSystem<FP> system;
    system.createParticle(n_ptcl*4+1000);
    system.setNumberOfParticleLocal(n_ptcl+n_ptcl_copy);

    PS::F32vec pos_shift(100.0, 200.0, 300.0);
    //PS::F32vec pos_shift(0.0, 0.0, 0.0);

    for(int i=0; i<n_ptcl; i++){
        system[i].id = i;
        finput>>system[i].mass;
    }
    for(int i=0; i<n_ptcl; i++){
        finput>>system[i].pos;
        system[i].pos += pos_shift;
    }
    for(int i=0; i<n_ptcl; i++){
        finput>>system[i].vel;
    }
    for(int i=n_ptcl; i<n_ptcl+n_ptcl_copy; i++){
        system[i].id = i;
        system[i] = system[0];
    }

    PS::F64 half_len = -1.0;
    for(PS::S32 i=0; i<n_ptcl; i++){
        PS::F64 tmp = system[i].pos.applyEach( PS::Abs<PS::F64>() ).getMax();
        half_len = tmp > half_len ? tmp : half_len;
    }
    std::cerr<<"half_len="<<half_len<<std::endl;

    PS::TreeForForce<PS::SEARCH_MODE_LONG, Force, EPI, EPJ, Moment, SPJ> tree;
    tree.initialize(n_ptcl+n_ptcl_copy);

    tree.initializeLocalTree(half_len);

    tree.setParticleLocalTree(system);

    tree.mortonSortLocalTreeOnly();

    tree.linkCellLocalTreeOnly();

    tree.checkMakeLocalTree();

    tree.calcMomentLocalTreeOnly();

    tree.checkCalcMomentLocalTree();




    PS::Finalize();
    return 0;

}
