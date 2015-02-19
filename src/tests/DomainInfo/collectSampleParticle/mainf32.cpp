#include <iostream>
#include <vector>
#include <cassert>

#include <particle_simulator.hpp>
#include "../../basic_particle.hpp"
#include "../check_domain_info.hpp"

template <class Tptcl>
void generateSphere(PS::U32 seed,
                    PS::S32 ntot,
                    Tptcl & system,
                    PS::F64 radius,
                    PS::F64vec center)
{
    PS::S32 rank = PS::Comm::getRank();
    PS::S32 size = PS::Comm::getNumberOfProc();

    PS::S32 ibgn = (ntot * rank) / size;
    PS::S32 iend = (ntot * (rank + 1)) / size;
    PS::S32 nloc = iend - ibgn;    
    PS::MT::init_genrand(seed);
    for(PS::S32 i = 0; i < nloc; i++)
        system[i].generateSphere(i+ibgn, radius, center);
    system.setNumberOfParticleLocal(nloc);
    
    return;
}                    

template <class Tptcl>
void writeAscii(Tptcl & system)
{
    PS::S32 rank = PS::Comm::getRank();
    char filename[1024];
    sprintf(filename, "snap.init.%04d", rank);

    PS::S32 nloc = system.getNumberOfParticleLocal();
    FILE *fp = fopen(filename, "w");
    for(PS::S32 i = 0; i < nloc; i++)
        system[i].writeAscii(fp);
    fclose(fp);

    return;
}

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32    nmem = 131072;
    PS::S32    ntot = 83927;
    PS::F64    prad = 3.1d;
    PS::F64vec pcen(1.d, -10.d, 199.8d);
    PS::U32    seed = PS::Comm::getRank();
    
    PS::S32    code  = 0;
    bool success_loc = false;
    bool success_glb = false;
    
    PS::DomainInfo dinfo;
    PS::ParticleSystem<BasicParticle32> bp;

    dinfo.initialize();
    bp.initialize();
    bp.createParticle(nmem);    
    generateSphere(seed, ntot, bp, prad, pcen);

    dinfo.collectSampleParticle(bp);

    success_loc = dinfo.checkCollectSampleParticleSubset(bp);
    success_glb = PS::Comm::synchronizeConditionalBranchAND(success_loc);
    code = (success_glb) ? code : (code | 1);

    success_loc = dinfo.checkCollectSampleParticleAverage(bp);
    success_glb = PS::Comm::synchronizeConditionalBranchAND(success_loc);
    code = (success_glb) ? code : (code | (1 << 1));

    PS::Finalize();       

    return code;
}
