#include "particle_simulator.hpp"
//#include "particle_mesh/particle_mesh.hpp"
#include "particle_mesh.hpp"

class Force{
public:
    PS::F64vec acc;
    PS::F64 pot;
    PS::S64 nj;
    PS::S64 njreal;
    void clear(){
        acc = 0.0;
        pot = 0.0;
        nj  = 0;
        njreal = 0;
    }
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
        this->nj  = force.nj;
        this->njreal = force.njreal;
    }
    // This setter is required for Particle Mesh...
    void copyFromForceParticleMesh(const PS::F32vec & force) {
        this->apm = force;
    }
    PS::S64 id;
    PS::F64 mass;
    // This getter is required for Particle Mesh...
    PS::F64 getChargeParticleMesh() const {
        return this->mass;
    }
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F32vec apm;
    PS::F64 pot;
    PS::S64 nj;
    PS::S64 njreal;
    static PS::F64 cutoff;
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
        fout<<"acc="<<acc<<std::endl;
        fout<<"pot="<<pot<<std::endl;
    }
};

PS::F64 FP::cutoff = 0.1875;

class EPI{
public:
    PS::F64vec pos;
    PS::S64 id;
    PS::S64 nj;
    PS::S64 njreal;
    static PS::F64 eps;
    void copyFromFP(const FP & fp){ 
        pos = fp.pos;
        id  = fp.id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
};

PS::F64 EPI::eps = 2.5e-4;

class EPJ{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    static PS::F64 cutoff;
    void copyFromFP(const FP & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
        id   = fp.id;
        cutoff = fp.cutoff;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64vec setPos(const PS::F64vec pos_new) {
        pos = pos_new;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64 getRSearch() {
        return this->cutoff;
    }
    PS::F64 getRSearch() const {
        return this->cutoff;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

PS::F64 EPJ::cutoff = 0.;

inline PS::F64 getP3M(PS::F64 xi) {
    PS::F64 p3m;

    xi = (xi > 2.) ? 2. : xi;
    PS::F64 zt  = (xi <= 1.) ? 0. : xi - 1.;
    
    PS::F64 xi2 = xi * xi;
    PS::F64 xi3 = xi * xi2;
    PS::F64 zt6 = zt * zt * zt * zt * zt * zt;

    p3m = 1. + xi3 * (-1.6 + xi2 * (1.6 + xi * (-0.5 + xi * (-0.34285714285 + xi * 0.15))))
        - zt6 * (0.08571428571 + xi * (0.51428571428 + xi * 0.2));

    return p3m;
}

#include <limits>

struct calcForceEpEp {
    void operator () (const EPI *epi,
                      const PS::S32 nip,
                      const EPJ *epj,
                      const PS::S32 njp,
                      Force *force) {

        PS::F64 eps2 = epi[0].eps * epi[0].eps;
        PS::F64 cinv = 1. / FP::cutoff;
        for(PS::S32 i = 0; i < nip; i++) {
            PS::S64 idi     = epi[i].id;
            PS::F64vec posi = epi[i].pos;
            PS::F64vec acci = 0.;
            PS::F64 poti    = 0.;
            PS::S64 njreal  = 0;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::S64 idj     = epj[j].id;
                PS::F64 mj      = epj[j].mass;
                PS::F64vec posj = epj[j].pos;
                PS::F64vec dx   = posj - posi;

                PS::F64 r2    = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + eps2;
                PS::F64 r1    = sqrt(r2);
                PS::F64 r3inv = 1. / (r1 * r1 * r1);

                PS::F64 xi    = 2. * cinv * r1;
                PS::F64 p3m   = getP3M(xi);
                
                acci += mj * r3inv * p3m * dx;
            }
            force[i].acc += acci;
            force[i].pot += poti;
            force[i].nj  += njp;
            force[i].njreal += njreal;
        }

    }
};

struct calcForceEpSp {
    void operator () (const EPI *epi,
                      const PS::S32 nip,
                      const PS::SPJMonopolePeriodic *spj,
                      const PS::S32 njp,
                      Force *force) {

        PS::F64 eps2 = epi[0].eps * epi[0].eps;
        PS::F64 cinv = 1. / FP::cutoff;
        for(PS::S32 i = 0; i < nip; i++) {
            PS::F64vec posi = epi[i].pos;
            PS::F64vec acci = 0.;
            PS::F64 poti    = 0.;
            PS::S64 njreal  = 0;
            for(PS::S32 j = 0; j < njp; j++) {
                PS::F64 mj      = spj[j].mass;
                PS::F64vec posj = spj[j].pos;
                PS::F64vec dx   = posj - posi;

                PS::F64 r2    = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + eps2;
                PS::F64 r1    = sqrt(r2);
                PS::F64 r3inv = 1. / (r1 * r1 * r1);

                PS::F64 xi    = 2. * cinv * r1;
                PS::F64 p3m   = getP3M(xi);

                acci += mj * r3inv * p3m * dx;
            }
            force[i].acc += acci;
            force[i].pot += poti;
            force[i].nj  += njp;
            force[i].njreal += njreal;
        }
    }
};

