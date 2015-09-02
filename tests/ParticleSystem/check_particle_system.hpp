namespace ParticleSimulator {

    template<class Tptcl>
    bool ParticleSystem<Tptcl>::checkExchangeParticleAllParticleInside(DomainInfo & dinfo)
    {
        bool success_loc = true;
        bool success_glb = false;

        PS::S32    nloc   = this->getNumberOfParticleLocal();
        PS::S32    rank   = PS::Comm::getRank();
        PS::F64ort domain = dinfo.getPosDomain(rank);

        for(S32 i = 0; i < nloc; i++) {
            for(S32 k = 0; k < DIMENSION; k++) {
                success_loc *= ((this->ptcl_[i].getPos())[k] >= domain.low_[k]) * ((this->ptcl_[i].getPos())[k] < domain.high_[k]);
            }
        }

        success_glb = PS::Comm::synchronizeConditionalBranchAND(success_loc);
        
        return success_glb;
    }

    template<class Tptcl>
    bool ParticleSystem<Tptcl>::checkExchangeParticleSumOfNumberOfParticle(DomainInfo & dinfo,
                                                                           S32 ntot_init)
    {
        PS::S32 rank = Comm::getRank();
        PS::S32 nloc = this->getNumberOfParticleLocal();
        PS::S32 ntot = Comm::getSum(nloc);

        bool success_glb = (ntot == ntot_init);

        return success_glb;
    }
}
