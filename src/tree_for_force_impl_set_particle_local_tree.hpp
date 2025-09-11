#pragma once
namespace ParticleSimulator {

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tpsys>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setParticleLocalTree(const Tpsys &psys, const bool clear) {
    const F64 t0 = GetWtime();
    const auto n_loc = psys.getNumberOfParticleLocal();
    epi_sorted_.resizeNoInitialize(n_loc);
    setParticleLocalTreeImpl(psys, clear);
    //time_profile_.calc_force_all__set_particle_local_tree += GetWtime() - t0;
    //std::cout<<"BBBBBBBBBBBBBBB"<<std::endl;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tpsys>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setParticleLocalTreeImpl(const Tpsys &psys,
                                                                                                                 const bool clear) {
    const S32 nloc = psys.getNumberOfParticleLocal();
    if (clear) {
        n_loc_tot_ = 0;
        inner_boundary_of_local_tree_.init();
        outer_boundary_of_local_tree_.init();
    }
    const S32 offset = n_loc_tot_;
    n_loc_tot_ += nloc;
    epj_org_.resizeNoInitialize(n_loc_tot_);
    epi_org_.resizeNoInitialize(n_loc_tot_);

#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
    const S32 n_thread = Comm::getNumberOfThread();
#else
    const S32 n_thread = 1;
#endif
    F64ort inner[n_thread];
    F64ort outer[n_thread];
    PS_OMP_PARALLEL {
        const S32 ith = Comm::getThreadNum();
        inner[ith].init();
        outer[ith].init();
        PS_OMP_FOR
        for (S32 i = 0; i < nloc; i++) {
            epi_org_[i + offset].copyFromFP(psys[i]);
            epj_org_[i + offset].copyFromFP(psys[i]);
            inner[ith].merge(psys[i].getPos());
            if constexpr (HasgetRSearchMethod<Tepj>::value) {
                outer[ith].merge(psys[i].getPos(), epj_org_[i + offset].getRSearch());
            } else if constexpr (HasgetRSearchMethod<Tepi>::value) {
                outer[ith].merge(psys[i].getPos(), epi_org_[i + offset].getRSearch());
            } else {
                outer[ith].merge(psys[i].getPos());
            }
        }
    }
    for (S32 i = 0; i < n_thread; i++) {
        inner_boundary_of_local_tree_.merge(inner[i]);
        outer_boundary_of_local_tree_.merge(outer[i]);
    }
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT)
    PARTICLE_SIMULATOR_PRINT_LINE_INFO();
    std::cout << "nloc=" << nloc << std::endl;
    std::cout << "n_loc_tot_=" << n_loc_tot_ << std::endl;
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename... Args>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setParticleLocalTreeVA(Args &&...args) {
    auto t0 = GetWtime();
    n_loc_tot_ = 0;
    countNumberOfParticleSystems(args...);
    epi_sorted_.resizeNoInitialize(n_loc_tot_);
    epj_sorted_.resizeNoInitialize(n_loc_tot_);
    epi_org_.resizeNoInitialize(n_loc_tot_);
    epj_org_.resizeNoInitialize(n_loc_tot_);
    inner_boundary_of_local_tree_.init();
    outer_boundary_of_local_tree_.init();
    n_loc_tot_ = 0;
    setParticleLocalTreeVAImpl(args...);
    time_profile_.set_particle_lt += GetWtime() - t0;
    //time_profile_.calc_force_all__set_particle_local_tree += GetWtime() - t0;
    //std::cout<<"time_profile_.calc_force_all__set_particle_local_tree= "<<time_profile_.calc_force_all__set_particle_local_tree<<std::endl;
    //std::cout<<"AAAAAAAAAAAAAAA"<<std::endl;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Head, typename... Tail>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setParticleLocalTreeVAImpl(Head &head, Tail &...tail) {
    if constexpr (IsParticleSystem<Head>::value) {
        const auto n_thread = Comm::getNumberOfThread();
        F64ort inner[n_thread];
        F64ort outer[n_thread];
        const auto n0 = n_loc_tot_;
        const auto n_loc = head.getNumberOfParticleLocal();
        n_loc_tot_ += n_loc;
        PS_OMP_PARALLEL {
            auto ith = Comm::getThreadNum();
            inner[ith].init();
            outer[ith].init();
            S32 idx0, idx1;
            CalcAdrToSplitData(idx0, idx1, ith, n_thread, n_loc);
            for (int i = idx0; i < idx1; i++) {
                epi_org_[i + n0].copyFromFP(head[i]);
                epj_org_[i + n0].copyFromFP(head[i]);
                inner[ith].merge(head[i].getPos());
                if constexpr (HasgetRSearchMethod<Tepj>::value) {
                    outer[ith].merge(epj_org_[i + n0].getPos(), epj_org_[i + n0].getRSearch() * 1.000001);
                } else if constexpr (HasgetRSearchMethod<Tepi>::value) {
                    outer[ith].merge(epi_org_[i + n0].getPos(), epi_org_[i + n0].getRSearch() * 1.000001);
                } else {
                    outer[ith].merge(epi_org_[i + n0].getPos());
                }
            }
        }  // PS_OMP_PARALLEL
        for (S32 i = 0; i < n_thread; i++) {
            inner_boundary_of_local_tree_.merge(inner[i]);
            outer_boundary_of_local_tree_.merge(outer[i]);
        }
    }
    if constexpr (sizeof...(tail) > 0) {
        setParticleLocalTreeVAImpl(tail...);
    }
}

// helper functions
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Head, typename... Tail>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::countNumberOfParticleSystems(Head &head, Tail &&...tail) {
    if constexpr (IsParticleSystem<Head>::value == true) {
        n_loc_tot_ += head.getNumberOfParticleLocal();
    }
    if constexpr (sizeof...(tail) > 0) {
        countNumberOfParticleSystems(tail...);
    }
}
}  // namespace ParticleSimulator
