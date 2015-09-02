namespace ParticleSimulator {

    namespace CheckDomainInfo {
        bool compareF64vecx(const F64vec & left,
                            const F64vec & right)
        {
            return left[0] < right[0];
        }

        F64vec calcPosArithmeticMean(const std::vector<F64vec> & array)
        {
            S32    n    = array.size();
            F64vec pave = 0.;
            for(S32 i = 0; i < n; i++)
                pave += array[i];
            pave /= (F64)n;
            return pave;
        }

        F64vec calcPosStandardDeviation(const std::vector<F64vec> & array,
                                        const F64vec pave)
        {
            S32    n    = array.size();
            F64vec pdis = 0.;
            for(S32 i = 0; i < n; i++)
                for(S32 k = 0; k < DIMENSION; k++)
                    pdis[k] += (array[i][k] - pave[k]) * (array[i][k] - pave[k]);
            for(S32 k = 0; k < DIMENSION; k++)
                pdis[k] = sqrt(pdis[k] / (F64)n);
            return pdis;
        }

    }

    template<class Tpsys>
    bool DomainInfo::checkCollectSampleParticleSubset(Tpsys & psys)
    {
        S32 rank = PS::Comm::getRank();
        std::vector<F64vec> cp_pos_psys;
        std::vector<F64vec> cp_pos_sample_loc;

        S32 nloc = psys.getNumberOfParticleLocal();
        for(S32 i = 0; i < nloc; i++)
            cp_pos_psys.push_back(psys[i].getPos());        
        std::sort(cp_pos_psys.begin(),
                  cp_pos_psys.end(),
                  CheckDomainInfo::compareF64vecx);

        for(S32 i = 0; i < number_of_sample_particle_loc_; i++)
            cp_pos_sample_loc.push_back(pos_sample_loc_[i]);
        std::sort(cp_pos_sample_loc.begin(),
                  cp_pos_sample_loc.end(),
                  CheckDomainInfo::compareF64vecx);

        S32 nsuc = 0;
        for(S32 i = 0; i < number_of_sample_particle_loc_; i++) {
            S32 ism  = 0;
            S32 ilg  = nloc - 1;
            S32 imd  = (ism + ilg) / 2;
            S32 ntry = 0;
            while(ntry < nloc) {
                if(cp_pos_sample_loc[i][0] == cp_pos_psys[imd][0]) {
                    if(cp_pos_sample_loc[i][1] == cp_pos_psys[imd][1]
                        && cp_pos_sample_loc[i][2] == cp_pos_psys[imd][2]) {
                        nsuc++;
                        break;
                    }
                }
                if(cp_pos_sample_loc[i][0] < cp_pos_psys[imd][0]) {
                    ilg = imd;
                    imd = (ism + ilg) / 2;
                } else {
                    ism = imd;
                    imd = (ism + ilg) / 2;
                }
                ntry++;
            }
        }
        
        bool success_loc = false;        
        if(nsuc == number_of_sample_particle_loc_)
            success_loc = true;
        bool success_glb = Comm::synchronizeConditionalBranchAND(success_loc);

        return success_glb;
    }

    template<class Tpsys>
    bool DomainInfo::checkCollectSampleParticleAverage(Tpsys & psys)
    {
        std::vector<F64vec> cp_pos_psys;
        std::vector<F64vec> cp_pos_sample_loc;
        F64vec pos_psys_ave;
        F64vec pos_psys_dis;
        F64vec pos_sample_ave;

        S32 nloc = psys.getNumberOfParticleLocal();
        for(S32 i = 0; i < nloc; i++)
            cp_pos_psys.push_back(psys[i].getPos());

        for(S32 i = 0; i < number_of_sample_particle_loc_; i++)
            cp_pos_sample_loc.push_back(pos_sample_loc_[i]);

        pos_psys_ave   = CheckDomainInfo::calcPosArithmeticMean(cp_pos_psys);
        pos_psys_dis   = CheckDomainInfo::calcPosStandardDeviation(cp_pos_psys, pos_psys_ave);
        pos_sample_ave = CheckDomainInfo::calcPosArithmeticMean(cp_pos_sample_loc);        

        F64vec low  = pos_psys_ave - 3.d * pos_psys_dis;
        F64vec high = pos_psys_ave + 3.d * pos_psys_dis;

        bool success_loc = false;        
        if(low[0] <= pos_sample_ave[0] && pos_sample_ave[0] <= high[0]
           && low[1] <= pos_sample_ave[1] && pos_sample_ave[1] <= high[1]
           && low[2] <= pos_sample_ave[2] && pos_sample_ave[2] <= high[2])
            success_loc = true;
        bool success_glb = Comm::synchronizeConditionalBranchAND(success_loc);

        return success_glb;
    }

    template<class Tpsys>
    bool DomainInfo::checkDecomposeDomain(Tpsys & psys)
    {
        S32 size = Comm::getNumberOfProc();
        S32 nloc = psys.getNumberOfParticleLocal();
        S32 ntot = Comm::getSum(nloc);

        bool success_loc = false;
        if(nloc == ntot / size)
            success_loc = true;
        bool success_glb = Comm::synchronizeConditionalBranchAND(success_loc);

        return success_glb;
    }

}
