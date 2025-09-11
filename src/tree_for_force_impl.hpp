////////////////////////////////////////////////
/// implementaion of methods of TreeForForce ///

#include "tree_walk.hpp"
#include "tree_for_force_impl_exlet.hpp"
#if defined(PARTICLE_SIMULATOR_USE_SAMPLE_SORT)
#include "samplesortlib.hpp"
#endif

namespace ParticleSimulator {
/////////////
// WRITE BACK
#if 1

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Head, typename... Tail>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::writeBackImpl(S32 &offset, Head &head, Tail &&...tail) {
    if constexpr (IsParticleSystem<Head>::value) {
        const F64 time_offset = GetWtime();
        const auto n_loc = head.getNumberOfParticleLocal();
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_loc; i++) head[i].copyFromForce(force_org_[i + offset]);
        offset += n_loc;
        time_profile_.write_back += GetWtime() - time_offset;
    }
    if constexpr (sizeof...(tail) > 0) {
        writeBackImpl(offset, std::forward<Tail>(tail)...);
    }
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Head, typename... Tail>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::writeBack(Head &head, Tail &&...tail) {
    S32 offset = 0;
    auto t0 = GetWtime();
    writeBackImpl(offset, head, std::forward<Tail>(tail)...);
    time_profile_.write_back += GetWtime() - t0;
}

#else

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tpsys>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::writeBack(Tpsys &psys) {
    const F64 time_offset = GetWtime();
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
    time_profile_.write_back += GetWtime() - time_offset;
}

#endif

// WRITE BACK
/////////////

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setCommInfo(const CommInfo &c) {
    comm_info_ = c;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setDomainInfo(const DomainInfo &dinfo) {
    p_domain_info_ = const_cast<DomainInfo *>(&dinfo);
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setDomainInfo(DomainInfo &dinfo) {
    p_domain_info_ = &dinfo;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
Tepj *TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getEpjFromId(const S64 id, const Tepj *epj_tmp) {
    PS_OMP_CRITICAL {
        if (map_id_to_epj_.empty()) {
            S64 n_epj = epj_sorted_.size();
            for (S32 i = 0; i < n_epj; i++) {
                if (GetMSB(tp_glb_[i].adr_ptcl_) == 1) continue;
                Tepj *epj_tmp = epj_sorted_.getPointer(i);
                S64 id_tmp = epj_tmp->getId();
                map_id_to_epj_.insert(std::pair<S64, Tepj *>(id_tmp, epj_tmp));
            }
        }
    }
    Tepj *epj = NULL;
    typename MyMap::iterator it = map_id_to_epj_.find(id);
    if (it != map_id_to_epj_.end()) epj = it->second;
    return epj;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getMemSizeUsed() const {
    size_t tmp = 0;
    for (int i = 0; i < Comm::getNumberOfThread(); i++) {
        //tmp = epj_for_force_[i].getMemSize() + spj_for_force_[i].getMemSize() + epjr_send_buf_[i].getMemSize() +
        //      epjr_send_buf_for_scatter_[i].getMemSize() + epjr_recv_1st_sorted_[i].getMemSize();
        tmp = epjr_send_buf_[i].getMemSize() +
              epjr_send_buf_for_scatter_[i].getMemSize() + epjr_recv_1st_sorted_[i].getMemSize();        
    }
    return tmp + tp_glb_.getMemSize() + tc_loc_.getMemSize() + tc_glb_.getMemSize() + epi_sorted_.getMemSize() + epi_org_.getMemSize() +
           epj_sorted_.getMemSize() + epj_org_.getMemSize() + spj_sorted_.getMemSize() + spj_org_.getMemSize() + ipg_.getMemSize() +
           epj_send_.getMemSize() + spj_send_.getMemSize() + force_sorted_.getMemSize() + force_org_.getMemSize() + epjr_sorted_.getMemSize() +
           epjr_send_.getMemSize() + epjr_recv_.getMemSize() + epjr_recv_1st_buf_.getMemSize() + epjr_recv_2nd_buf_.getMemSize();
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getUsedMemorySize() const {
    return getMemSizeUsed();
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::initialize(const U64 n_glb_tot, const F64 theta,
                                                                                                   const U32 n_leaf_limit, const U32 n_group_limit) {
    if (is_initialized_ == true) {
        PARTICLE_SIMULATOR_PRINT_ERROR("Do not initialize the tree twice");
        std::cerr << "SEARCH_MODE: " << typeid(TSM).name() << std::endl;
        std::cerr << "Force: " << typeid(Tforce).name() << std::endl;
        std::cerr << "EPI: " << typeid(Tepi).name() << std::endl;
        std::cerr << "EPJ: " << typeid(Tepj).name() << std::endl;
        std::cerr << "SPJ: " << typeid(Tspj).name() << std::endl;
        Abort(-1);
    }
    lev_group_limit_ = lev_leaf_limit_ = 0;
    coef_ipg_size_limit_ = 1.0;
    comm_info_.setCommunicator();
    map_id_to_epj_.clear();
    is_initialized_ = true;
    n_glb_tot_ = n_glb_tot;
    theta_ = theta;
    n_leaf_limit_ = n_leaf_limit;
    n_group_limit_ = n_group_limit;
    lev_max_loc_ = lev_max_glb_ = 0;
    const S32 n_thread = Comm::getNumberOfThread();

    n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
    n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
    ni_ave_ = nj_ave_ = 0;
    // Comm::barrier();
    comm_info_.barrier();
    bool err = false;
    if (n_leaf_limit_ <= 0) {
        err = true;
        PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of the particles in the leaf cell must be > 0");
        std::cout << "n_leaf_limit_= " << n_leaf_limit_ << std::endl;
    }
    if (n_group_limit_ < n_leaf_limit_) {
        err = true;
        PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of particles in ip graoups msut be >= that in leaf cells");
        std::cout << "n_group_limit_= " << n_group_limit_ << std::endl;
        std::cout << "n_leaf_limit_= " << n_leaf_limit_ << std::endl;
    }
    if (typeid(TSM) == typeid(SEARCH_MODE_LONG) || typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF)) {
        if (theta_ < 0.0) {
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The opening criterion of the tree must be >= 0.0");
            std::cout << "theta_= " << theta_ << std::endl;
        }
    }
    if (err) {
        std::cout << "SEARCH_MODE: " << typeid(TSM).name() << std::endl;
        std::cout << "Force: " << typeid(Tforce).name() << std::endl;
        std::cout << "EPI: " << typeid(Tepi).name() << std::endl;
        std::cout << "SPJ: " << typeid(Tspj).name() << std::endl;
        Abort(-1);
        //ParticleSimulator::Abort(-1);
    }

#if defined(PARTICLE_SIMULATOR_USE_MEMORY_POOL)
    epi_sorted_.setAllocMode(1);  // msortLT --- final
    spj_sorted_.setAllocMode(1);  //  --- final
    epi_org_.setAllocMode(1);     // setPtclLT ---
    spj_org_.setAllocMode(1);     // insted of it, use spj_recv
    epj_org_.setAllocMode(1);
    force_sorted_.setAllocMode(1);  // -- final
    epj_send_.setAllocMode(1);
    spj_send_.setAllocMode(1);
    epj_for_force_ = new ReallocatableArray<Tepj>[n_thread];
    spj_for_force_ = new ReallocatableArray<Tspj>[n_thread];
    for (S32 i = 0; i < n_thread; i++) {
        epj_for_force_[i].setAllocMode(1);
        spj_for_force_[i].setAllocMode(1);
    }
#else
    /*
    epj_org_.setAllocMode(MemoryAllocMode::Pool);
    epi_org_.setAllocMode(MemoryAllocMode::Pool);
    epj_send_.setAllocMode(MemoryAllocMode::Pool);
    spj_send_.setAllocMode(MemoryAllocMode::Pool);
    epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
    spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
    for(S32 i=0; i<n_thread; i++){
        epj_for_force_[i].setAllocMode(MemoryAllocMode::Pool);
        spj_for_force_[i].setAllocMode(MemoryAllocMode::Pool);
    }
    */
    epj_org_.setAllocMode(MemoryAllocMode::Default);
    epi_org_.setAllocMode(MemoryAllocMode::Default);
    epj_send_.setAllocMode(MemoryAllocMode::Default);
    spj_send_.setAllocMode(MemoryAllocMode::Default);
    //epj_for_force_ = new ReallocatableArray<Tepj>[n_thread];
    //spj_for_force_ = new ReallocatableArray<Tspj>[n_thread];
    for (S32 i = 0; i < n_thread; i++) {
        //epj_for_force_[i].setAllocMode(MemoryAllocMode::Default);
        //spj_for_force_[i].setAllocMode(MemoryAllocMode::Default);
    }
#endif

    adr_epj_for_force_ = new ReallocatableArray<S32>[n_thread];
    adr_spj_for_force_ = new ReallocatableArray<S32>[n_thread];
    adr_ipg_for_force_ = new ReallocatableArray<S32>[n_thread];

    epjr_send_buf_ = new ReallocatableArray<EPJWithR>[n_thread];
    epjr_send_buf_for_scatter_ = new ReallocatableArray<EPJWithR>[n_thread];
    epjr_recv_1st_sorted_ = new ReallocatableArray<EPJWithR>[n_thread];
    epj_neighbor_ = new ReallocatableArray<Tepj>[n_thread];
    n_cell_open_ = new CountT[n_thread];
    comm_info_.barrier();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    const S64 n_proc = comm_info_.getNumberOfProc();
    req_send_ = new MPI_Request[n_proc];
    req_recv_ = new MPI_Request[n_proc];
    status_ = new MPI_Status[n_proc];
#endif
    SET_VAR_NAME(epi_sorted_);
    SET_VAR_NAME(epj_sorted_);
    SET_VAR_NAME(force_sorted_);
    SET_VAR_NAME(epi_org_);
    SET_VAR_NAME(epj_org_);
    SET_VAR_NAME(spj_org_);
    SET_VAR_NAME(epj_send_);
    SET_VAR_NAME(spj_send_);
    for (S32 i = 0; i < n_thread; i++) {
        //SET_VAR_NAME(epj_for_force_[i]);
        //SET_VAR_NAME(spj_for_force_[i]);
    }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
    MemoryPool::dump();
    if (comm_info_.getRank() == 0) {
        std::cerr << "used mem size for tree= " << this->getMemSizeUsed() * 1e-9 << "[GB]" << std::endl;
    }
#endif
}

/*
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMomentLocalTreeOnlyImpl(TagCalcMomLongEpjLt) {
    CalcMomentLongLocalTree(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_, morton_key_);
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMomentLocalTreeOnlyImpl(TagCalcMomShortEpjLt) {
    CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMomentLocalTreeOnlyImpl(TagCalcMomShortEpiLt) {
    CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epi_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
}
*/

/////////////////////////////
/// MAKE INTERACTION LIST ///
// pj is epj or spj
template <class Tpj>
void CopyPjForForceST(const ReallocatableArray<S32> &adr_pj, const ReallocatableArray<Tpj> &pj_sorted, ReallocatableArray<Tpj> &pj_for_force) {
    const S32 n_pj = adr_pj.size();
    pj_for_force.resizeNoInitialize(n_pj);
    for (S32 i = 0; i < n_pj; i++) {
        const S32 adr_pj_src = adr_pj[i];
        pj_for_force[i] = pj_sorted[adr_pj_src];
    }
}
template <class Tpj>
void CopyPjForForceST(const ReallocatableArray<S32> &adr_pj, const ReallocatableArray<Tpj> &pj_sorted, const S32 n_head, const S32 n_tail,
                      ReallocatableArray<Tpj> &pj_for_force) {
    pj_for_force.resizeNoInitialize(n_tail);
    for (S32 i = n_head; i < n_tail; i++) {
        const S32 adr_pj_src = adr_pj[i];
        pj_for_force[i] = pj_sorted[adr_pj_src];
    }
}
template <class Tpj, class Ttc>
void CopyPjForForceST(const ReallocatableArray<S32> &adr_pj, const ReallocatableArray<Tpj> &pj_sorted, const ReallocatableArray<Ttc> &tc,
                      const S32 n_let_sp, const S32 n_head, const S32 n_tail, ReallocatableArray<Tpj> &pj_for_force) {
    pj_for_force.resizeNoInitialize(n_tail);
    for (S32 i = n_head; i < n_tail; i++) {
        const auto adr_pj_src = adr_pj[i];
        if (adr_pj_src < n_let_sp) {
            pj_for_force[i] = pj_sorted[adr_pj_src];
        } else {
            const auto adr_tc_src = adr_pj_src - n_let_sp;
            pj_for_force[i].copyFromMoment(tc[adr_tc_src].mom_);
        }
    }
}

template <class T>
struct TraitsForCutoff {
    typedef TagWithoutCutoff type_cutoff;
};
template <>
struct TraitsForCutoff<TagSearchLongCutoff> {
    typedef TagWithCutoff type_cutoff;
};
/*
// commented out 2025/02/22
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeInteractionListLongForZeroTheta(TagWithoutCutoff,
                                                                                                                            const S32 adr_ipg) {
    const S32 ith = Comm::getThreadNum();
    const S32 n_tmp = tc_glb_[0].n_ptcl_;
    S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
    epj_for_force_[ith].reserveEmptyAreaAtLeast(n_tmp);
    for (S32 ip = 0; ip < n_tmp; ip++) {
        if (GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0) {
            epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
        } else {
            adr_ptcl_tmp++;
        }
    }
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeInteractionListLongForZeroTheta(TagWithCutoff,
                                                                                                                            const S32 adr_ipg) {
    const S32 ith = Comm::getThreadNum();
    if (!tc_glb_[0].isLeaf(n_leaf_limit_)) {
        const F64 r_cut_sq = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_in_;
        const F64ort cell_box = pos_root_cell_;
        MakeInteractionListLongCutoffEPForZeroTheta(tc_glb_, tc_glb_[0].adr_tc_, tp_glb_, epj_sorted_, epj_for_force_[ith], cell_box, pos_target_box,
                                                    r_cut_sq, n_leaf_limit_);
    } else {
        const S32 n_tmp = tc_glb_[0].n_ptcl_;
        S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
        epj_for_force_[ith].reserveEmptyAreaAtLeast(n_tmp);
        for (S32 ip = 0; ip < n_tmp; ip++) {
            if (GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0) {
                epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
            } else {
                adr_ptcl_tmp++;
            }
        }
    }
}
// commented out 2025/02/22
*/

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::copyForceOriginalOrder() {
    epi_org_.freeMem(1);
    force_org_.resizeNoInitialize(n_loc_tot_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
    for (S32 i = 0; i < n_loc_tot_; i++) {
        const S32 adr = adr_org_from_adr_sorted_loc_[i];
        force_org_[adr] = force_sorted_[i];
    }
}

// return forces in original order
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tsys>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceDirectParallel(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                                Tforce force[], const Tsys &system,
                                                                                                                const DomainInfo &dinfo,
                                                                                                                const bool clear) {
    if (clear) {
        for (S32 i = 0; i < n_loc_tot_; i++) force[i].clear();
    }
    // S32 n_epj_loc_max = Comm::getMaxValue(n_loc_tot_);
    ReallocatableArray<Tepi> my_epi(n_loc_tot_, n_loc_tot_, MemoryAllocMode::Pool);
    ReallocatableArray<Tepj> my_epj(n_loc_tot_, n_loc_tot_, MemoryAllocMode::Pool);
    for (S32 i = 0; i < n_loc_tot_; i++) {
        const S32 adr = adr_org_from_adr_sorted_loc_[i];
        my_epi[adr] = epi_sorted_[i];
        my_epj[i].copyFromFP(system[i]);
    }
    // const S32 n_proc  = Comm::getNumberOfProc();
    const S32 n_proc = comm_info_.getNumberOfProc();
    // const S32 my_rank = Comm::getRank();
    const S32 my_rank = comm_info_.getRank();
    ReallocatableArray<F64vec> shift;
    F64ort pos_root_domain = dinfo.getPosRootDomain();
    F64vec domain_size = pos_root_domain.getFullLength();
    bool pa[DIMENSION];
    dinfo.getPeriodicAxis(pa);
    CalcNumberAndShiftOfImageDomain(shift, domain_size, pos_root_cell_, pos_root_domain, pa);
    const S32 n_image = shift.size();
    for (S32 i = 0; i < n_proc; i++) {
        auto rank_send = (my_rank + n_proc + i) % n_proc;
        auto rank_recv = (my_rank + n_proc - i) % n_proc;
        S32 n_recv;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Status stat;
        MPI_Sendrecv(&n_loc_tot_, 1, GetDataType<S32>(), rank_send, 0, &n_recv, 1, GetDataType<S32>(), rank_recv, 0, MPI_COMM_WORLD, &stat);
        ReallocatableArray<Tepj> epj(n_recv * n_image, n_recv * n_image, MemoryAllocMode::Pool);
        MPI_Sendrecv(my_epj.getPointer(), n_loc_tot_, GetDataType<Tepj>(), rank_send, 1, epj.getPointer(), n_recv, GetDataType<Tepj>(), rank_recv, 1,
                     MPI_COMM_WORLD, &stat);
#else
        n_recv = n_loc_tot_;
        ReallocatableArray<Tepj> epj(n_recv * n_image, n_recv * n_image, MemoryAllocMode::Pool);
        for (S32 j = 0; j < n_recv; j++) {
            epj[j] = my_epj[j];
        }
#endif
        for (S32 j = 0; j < n_image; j++) {
            for (S32 k = 0; k < n_recv; k++) {
                epj[j * n_image + k] = epj[k];
                const F64vec pos_new = epj[k].getPos() + shift[j];
                epj[j * n_image + k].pos = pos_new;
            }
        }
        pfunc_ep_ep(my_epi.getPointer(), n_loc_tot_, epj.getPointer(), n_recv * n_image, force);
    }
}

// return forces in original order
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceDirect(Tfunc_ep_ep pfunc_ep_ep, Tforce force[],
                                                                                                        const DomainInfo &dinfo, const bool clear) {
    if (clear) {
        for (S32 i = 0; i < n_loc_tot_; i++) force[i].clear();
    }
    Tepj *epj_tmp;
    S32 n_epj_tmp;
    bool pa[DIMENSION];
    dinfo.getPeriodicAxis(pa);
    AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain(), pos_root_cell_, pa, comm_info_);
    pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force);
    delete[] epj_tmp;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceDirectAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                                    const DomainInfo &dinfo,
                                                                                                                    const bool clear) {
    force_org_.resizeNoInitialize(n_loc_tot_);
    if (clear) {
        for (S32 i = 0; i < n_loc_tot_; i++) force_org_[i].clear();
    }
    Tepj *epj_tmp;
    S32 n_epj_tmp;
    bool pa[DIMENSION];
    dinfo.getPeriodicAxis(pa);
    AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa,
                      comm_info_);
    pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force_org_.getPointer());
    delete[] epj_tmp;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tptcl>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getNeighborListOneParticle(const Tptcl &ptcl, Tepj *&epj) {
    return getNeighborListOneParticleImpl(typename TSM::neighbor_search_type(), ptcl, epj);
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tptcl>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getNeighborListOneParticleImpl(TagNeighborSearchScatter,
                                                                                                                      const Tptcl &ptcl, Tepj *&epj) {
    const S32 id_thread = Comm::getThreadNum();
    epj_neighbor_[id_thread].clearSize();
    const F64vec pos_target = ptcl.getPos();
    const S32 adr = 0;
    SearchNeighborListOneParticleScatter(pos_target, tc_glb_.getPointer(), tp_glb_.getPointer(), adr, epj_sorted_, epj_neighbor_[id_thread],
                                         n_leaf_limit_);
    S32 nnp = epj_neighbor_[id_thread].size();
    epj = epj_neighbor_[id_thread].getPointer();
    return nnp;
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tptcl>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getNeighborListOneParticleImpl(TagNeighborSearchGather,
                                                                                                                      const Tptcl &ptcl, Tepj *&epj) {
    const S32 id_thread = Comm::getThreadNum();
    epj_neighbor_[id_thread].clearSize();
    const F64vec pos_target = ptcl.getPos();
    F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
    const S32 adr = 0;
    SearchNeighborListOneParticleGather(pos_target, r_search_sq, tc_glb_.getPointer(), tp_glb_.getPointer(), adr, epj_sorted_,
                                        epj_neighbor_[id_thread], n_leaf_limit_);
    S32 nnp = epj_neighbor_[id_thread].size();
    epj = epj_neighbor_[id_thread].getPointer();
    return nnp;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tptcl>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getNeighborListOneParticleImpl(TagNeighborSearchSymmetry,
                                                                                                                      const Tptcl &ptcl, Tepj *&epj) {
    const S32 id_thread = Comm::getThreadNum();
    epj_neighbor_[id_thread].clearSize();
    const F64vec pos_target = ptcl.getPos();
    F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
    const S32 adr = 0;
    SearchNeighborListOneParticleSymmetry(pos_target, r_search_sq, tc_glb_.getPointer(), tp_glb_.getPointer(), adr, epj_sorted_,
                                          epj_neighbor_[id_thread], n_leaf_limit_);
    S32 nnp = epj_neighbor_[id_thread].size();
    epj = epj_neighbor_[id_thread].getPointer();
    return nnp;
}
// New neighbor search mode (this does not require getRSearch)
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tptcl>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getNeighborListOneParticle(const Tptcl &ptcl, const S32 n_ngb,
                                                                                                                  Tepj *&epj) {
    const S32 id_thread = Comm::getThreadNum();
    epj_neighbor_[id_thread].clearSize();
    const F64vec pos_target = ptcl.getPos();
    std::vector<std::pair<Tepj, F64>> epj_md_list;  // md := metadata
    const S32 adr = 0;
    F64ort vertex;
    SearchNeighborListOneParticleNumber(pos_target, n_ngb, tc_glb_.getPointer(), tp_glb_.getPointer(), adr, epj_sorted_, epj_md_list, vertex,
                                        n_leaf_limit_);
    for (const auto &e : epj_md_list) epj_neighbor_[id_thread].push_back(e.first);
    S32 nnp = epj_neighbor_[id_thread].size();
    epj = epj_neighbor_[id_thread].getPointer();
    return nnp;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::clearSizeOfArray() {
    tp_glb_.clearSize();
    tc_loc_.clearSize();
    tc_glb_.clearSize();
    epi_sorted_.clearSize();
    epi_org_.clearSize();
    epj_sorted_.clearSize();
    epj_org_.clearSize();
    spj_sorted_.clearSize();
    spj_org_.clearSize();
    ipg_.clearSize();
    epj_send_.clearSize();
    spj_send_.clearSize();
    force_sorted_.clearSize();
    force_org_.clearSize();
    epjr_sorted_.clearSize();
    epjr_send_.clearSize();
    epjr_recv_.clearSize();
    epjr_recv_1st_buf_.clearSize();
    epjr_recv_2nd_buf_.clearSize();
    const S32 n_thread = Comm::getNumberOfThread();
    for (S32 i = 0; i < n_thread; i++) {
        //epj_for_force_[i].clearSize();
        //spj_for_force_[i].clearSize();
        epjr_send_buf_[i].clearSize();
        epjr_send_buf_for_scatter_[i].clearSize();
        epjr_recv_1st_sorted_[i].clearSize();
    }
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getOuterBoundaryOfLocalTreeImpl(TagSearchLongSymmetry) {
    return tc_loc_[0].geo_.getVertexOut();
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::getOuterBoundaryOfLocalTreeImpl(TagSearchShortSymmetry) {
    return tc_loc_[0].geo_.getVertexOut();
}

template <class TSM>
F64ort GetIpgBoxForInteractionList(const IPGroup<TSM> &ipg) {
    return ipg.vertex_;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::freeMem() {
    tp_glb_.freeMem();
    tc_loc_.freeMem();
    tc_glb_.freeMem();
    epi_sorted_.freeMem();
    epi_org_.freeMem();
    epj_sorted_.freeMem();
    epj_org_.freeMem();
    spj_sorted_.freeMem();
    spj_org_.freeMem();
    ipg_.freeMem();
    epj_send_.freeMem();
    spj_send_.freeMem();
    force_sorted_.freeMem();
    force_org_.freeMem();
    epjr_sorted_.freeMem();
    epjr_send_.freeMem();
    epjr_recv_.freeMem();
    epjr_recv_1st_buf_.freeMem();
    epjr_recv_2nd_buf_.freeMem();
    const S32 n_thread = Comm::getNumberOfThread();
    for (S32 i = 0; i < n_thread; i++) {
        //epj_for_force_[i].freeMem();
        //spj_for_force_[i].freeMem();
        epjr_send_buf_[i].freeMem();
        epjr_send_buf_for_scatter_[i].freeMem();
        epjr_recv_1st_sorted_[i].freeMem();
    }
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::dump(std::ostream &fout) {
    fout << "n_loc_tot_=" << n_loc_tot_ << std::endl;
    fout << "n_glb_tot_=" << n_glb_tot_ << std::endl;
    fout << "length_=" << length_ << std::endl;
    fout << "center_=" << center_ << std::endl;
    fout << "pos_root_cell_=" << pos_root_cell_ << std::endl;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::dumpMemSizeUsed(std::ostream &fout) {
    S32 n_thread = Comm::getNumberOfThread();
    if (comm_info_.getRank() == 0) {
        fout << "tp_glb_.getMemSize()= " << tp_glb_.getMemSize() << std::endl;
        fout << "tc_loc_.getMemSize()= " << tc_loc_.getMemSize() << std::endl;
        fout << "tc_glb_.getMemSize()= " << tc_glb_.getMemSize() << std::endl;
        fout << "epi_sorted_.getMemSize()= " << epi_sorted_.getMemSize() << std::endl;
        fout << "epi_org_.getMemSize()= " << epi_org_.getMemSize() << std::endl;
        fout << "epj_sorted_.getMemSize()= " << epj_sorted_.getMemSize() << std::endl;
        fout << "epj_org_.getMemSize()= " << epj_org_.getMemSize() << std::endl;
        fout << "spj_sorted_.getMemSize()= " << spj_sorted_.getMemSize() << std::endl;
        fout << "spj_org_.getMemSize()= " << spj_org_.getMemSize() << std::endl;
        fout << "ipg_.getMemSize()= " << ipg_.getMemSize() << std::endl;
        fout << "epj_send_.getMemSize()= " << epj_send_.getMemSize() << std::endl;
        fout << "spj_send_.getMemSize()= " << spj_send_.getMemSize() << std::endl;
        fout << "force_org_.getMemSize()= " << force_org_.getMemSize() << std::endl;
        fout << "force_sorted_.getMemSize()= " << force_sorted_.getMemSize() << std::endl;
        //for (S32 i = 0; i < n_thread; i++) fout << "epj_for_force_[" << i << "].getMemSize()= " << epj_for_force_[i].getMemSize() << std::endl;
        //for (S32 i = 0; i < n_thread; i++) fout << "spj_for_force_[" << i << "].getMemSize()= " << spj_for_force_[i].getMemSize() << std::endl;
        fout << "epjr_sorted_.getMemSie()= " << epjr_sorted_.getMemSize() << std::endl;
        fout << "epjr_send_.getMemSie()= " << epjr_send_.getMemSize() << std::endl;
        fout << "epjr_recv_.getMemSie()= " << epjr_recv_.getMemSize() << std::endl;
        fout << "epjr_recv_1st_buf_.getMemSie()= " << epjr_recv_1st_buf_.getMemSize() << std::endl;
        fout << "epjr_recv_2nd_buf_.getMemSie()= " << epjr_recv_2nd_buf_.getMemSize() << std::endl;
        for (S32 i = 0; i < n_thread; i++) fout << "epjr_send_buf_[" << i << "].getMemSize()= " << epjr_send_buf_[i].getMemSize() << std::endl;
        for (S32 i = 0; i < n_thread; i++)
            fout << "epjr_send_buf_for_scatter_[" << i << "].getMemSize()= " << epjr_send_buf_for_scatter_[i].getMemSize() << std::endl;
        for (S32 i = 0; i < n_thread; i++)
            fout << "epjr_recv_1st_sorted_[" << i << "].getMemSize()= " << epjr_recv_1st_sorted_[i].getMemSize() << std::endl;
        size_t size_epj_for_force_tot = 0;
        size_t size_spj_for_force_tot = 0;
        size_t size_epjr_send_buf_tot = 0;
        size_t size_epjr_send_buf_for_scatter_tot = 0;
        size_t size_epjr_recv_1st_sorted_tot = 0;
        for (S32 i = 0; i < n_thread; i++) {
            //size_epj_for_force_tot += epj_for_force_[i].getMemSize();
            //size_spj_for_force_tot += spj_for_force_[i].getMemSize();
            size_epjr_send_buf_tot += epjr_send_buf_[i].getMemSize();
            size_epjr_send_buf_for_scatter_tot += epjr_send_buf_for_scatter_[i].getMemSize();
            size_epjr_recv_1st_sorted_tot += epjr_recv_1st_sorted_[i].getMemSize();
        }

        fout << "sum= "
             << (double)(tp_glb_.getMemSize() + tc_loc_.getMemSize() + tc_glb_.getMemSize() + epi_sorted_.getMemSize() + epi_org_.getMemSize() +
                         epj_sorted_.getMemSize() + epj_org_.getMemSize() + spj_sorted_.getMemSize() + spj_org_.getMemSize() + ipg_.getMemSize() +
                         epj_send_.getMemSize() + spj_send_.getMemSize() + force_org_.getMemSize() + force_sorted_.getMemSize() +
                         size_epj_for_force_tot + size_spj_for_force_tot + epjr_sorted_.getMemSize() + epjr_send_.getMemSize() +
                         epjr_recv_.getMemSize() + epjr_recv_1st_buf_.getMemSize() + epjr_recv_2nd_buf_.getMemSize() + size_epjr_send_buf_tot +
                         size_epjr_send_buf_for_scatter_tot + size_epjr_recv_1st_sorted_tot) /
                    1e9
             << " [GB]" << std::endl;
    }
}
}  // namespace ParticleSimulator
#include "tree_for_force_impl_set_particle_local_tree.hpp"
#include "tree_for_force_impl_make_ipg.hpp"
#include "tree_for_force_impl_force.hpp"
