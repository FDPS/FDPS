#pragma once

#include <tree_for_force_utils_exlet.hpp>

namespace ParticleSimulator {

// These classes are used only for exchanging LET with P2P and Allgather
template <typename Tspj, typename Ttc, typename Tmode>
struct TopMoment {
    Tspj spj_;
    F64ort vertex_in_;
    void set(const Tspj& spj, const F64ort vertex_in, const Ttc& tc) {
        spj_ = spj;
        vertex_in_ = vertex_in;
    }
    template <typename Ttm>
    bool isOverlapped(const Ttm& tm) {
        PARTICLE_SIMULATOR_PRINT_ERROR("FDPS error. This search mode is not correct.");
        Abort();
        return false;
    }
    F64ort getVertexOut() const {
        PARTICLE_SIMULATOR_PRINT_ERROR("FDPS error. This search mode is not correct.");
        Abort();
        return F64ort(-1234.5, -9876.5);
    }
};

template <typename Tspj, typename Ttc>
struct TopMoment<Tspj, Ttc, SEARCH_MODE_LONG_SYMMETRY> {
    Tspj spj_;
    F64ort vertex_in_;
    F64ort vertex_out_;
    void set(const Tspj& spj, const F64ort vertex_in, const Ttc& tc) {
        spj_ = spj;
        vertex_in_ = vertex_in;
        vertex_out_ = tc.geo_.getVertexOut();
    }
    bool isOverlapped(const TopMoment<Tspj, Ttc, SEARCH_MODE_LONG_SYMMETRY>& tm) {
        return vertex_in_.overlapped(tm.vertex_out_) || vertex_out_.overlapped(tm.vertex_in_);
    }
    F64ort getVertexOut() const { return vertex_out_; }
};

template <typename Tspj, typename Ttc>
struct TopMoment<Tspj, Ttc, SEARCH_MODE_LONG_SCATTER> {
    Tspj spj_;
    F64ort vertex_in_;
    F64ort vertex_out_;  // To check if p2p is ok or not
    void set(const Tspj& spj, const F64ort vertex_in, const Ttc& tc) {
        spj_ = spj;
        vertex_in_ = vertex_in;
        vertex_out_ = tc.geo_.getVertexOut();
    }
    bool isOverlapped(const TopMoment<Tspj, Ttc, SEARCH_MODE_LONG_SCATTER>& tm) { return vertex_out_.overlapped(tm.vertex_in_); }
    F64ort getVertexOut() const { return vertex_out_; }
};

template <typename Ttc, typename Tmode>
struct TopGeometry {
    F64ort vertex_in_;
    void set(const F64ort vertex_in, const Ttc& tc) { vertex_in_ = vertex_in; }
};
template <typename Ttc>
struct TopGeometry<Ttc, SEARCH_MODE_LONG_SYMMETRY> {
    F64ort vertex_in_;
    F64ort vertex_out_;
    void set(const F64ort vertex_in, const Ttc& tc) {
        vertex_in_ = vertex_in;
        vertex_out_ = tc.geo_.getVertexOut();
    }
    F64ort getVertexOut() const { return vertex_out_; }
};

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::exchangeLocalEssentialTree(const DomainInfo& dinfo,
                                                                                                                   const bool flag_reuse) {
    PS_COMPILE_TIME_MESSAGE("in exchangeLocalEssentialTree");
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTree)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START exchangeLocalEssentialTree ****", 0);
#endif
    if (typeid(TSM) == typeid(SEARCH_MODE_LONG) && dinfo.getBoundaryCondition() != BOUNDARY_CONDITION_OPEN) {
        PARTICLE_SIMULATOR_PRINT_ERROR("The forces w/o cutoff can be evaluated only under the open boundary condition");
        Abort(-1);
    }
    if (!flag_reuse) {
        comm_table_.clearSize();
    }
    auto t0 = GetWtime();
    if constexpr (std::is_same_v<typename TSM::search_type, TagSearchLong> || std::is_same_v<typename TSM::search_type, TagSearchLongScatter> ||
                  std::is_same_v<typename TSM::search_type, TagSearchLongSymmetry> ||
                  std::is_same_v<typename TSM::search_type, TagSearchLongCutoff>) {
        PS_COMPILE_TIME_MESSAGE("TagSearchLong@exchangeLocalEssentialTree@exchangeLocalEssentialTree");
        if (exchange_let_mode_ == EXCHANGE_LET_A2A) {
            // exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
            exchangeLocalEssentialTreeLongNew(dinfo, flag_reuse);
        } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) {
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
        } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_FAST) {
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
            comm_table_.setPosDomainAllgather(dinfo, pos_root_cell_);
        }
    } else if constexpr (std::is_same_v<typename TSM::search_type, TagSearchLongParticleMeshMultipole>) {
        PS_COMPILE_TIME_MESSAGE("TagSearchLongParticleMeshMultipole@exchangeLocalEssentialTree");
        if (exchange_let_mode_ == EXCHANGE_LET_A2A) {
            // exchangeLocalEssentialTreeLongPMM(dinfo, flag_reuse);
            exchangeLocalEssentialTreeLongNew(dinfo, flag_reuse);
        } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) {
            // under construction
            // exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
        } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_FAST) {
            // under construction
            // exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
            // comm_table_.setPosDomainAllgather(dinfo, pos_root_cell_);
        }
    } else if constexpr (std::is_same_v<typename TSM::search_type, TagSearchShortScatter>) {
        PS_COMPILE_TIME_MESSAGE("TagSearchShortScatter@exchangeLocalEssentialTree");
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTree)
        PARTICLE_SIMULATOR_PRINT_MSG("- TagSearchShortScatter", 0);
#endif
        if (exchange_let_mode_ == EXCHANGE_LET_A2A) {
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTree)
            PARTICLE_SIMULATOR_PRINT_MSG("- before exchangeLocalEssentialTreeShortScatter", 0);
#endif
            exchangeLocalEssentialTreeShortScatter(dinfo, flag_reuse);

        } else {
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTree)
            PARTICLE_SIMULATOR_PRINT_MSG("- before exchangeLocalEssentialTreeShortScatterP2P", 0);
#endif
            exchangeLocalEssentialTreeShortScatterP2P(dinfo, flag_reuse);
        }
    } else if constexpr (std::is_same_v<typename TSM::search_type, TagSearchShortGather> ||
                         std::is_same_v<typename TSM::search_type, TagSearchShortSymmetry>) {
        if (exchange_let_mode_ == EXCHANGE_LET_A2A) {
            exchangeLocalEssentialTreeShort2StepA2A(dinfo, flag_reuse);
        } else {
            exchangeLocalEssentialTreeShort2StepP2P(dinfo, flag_reuse);
        }
    } else {
        static_assert([] { return false; }());
    }
    time_profile_.construct_exchange_let = GetWtime() - t0;
    //time_profile_.exchange_LET_tot =
    //    time_profile_.make_LET_1st + time_profile_.exchange_LET_1st + time_profile_.make_LET_2nd + time_profile_.exchange_LET_2nd;
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTree)
    PARTICLE_SIMULATOR_PRINT_MSG("END exchangeLocalEssentialTree", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::exchangeLocalEssentialTreeLongNew(const DomainInfo& dinfo,
                                                                                                                          const bool flag_reuse) {
    F64 t0 = 0.0;
    //const auto time_offset = GetWtime();
    PS_COMPILE_TIME_MESSAGE("exchangeLocalEssentialTreeLongNew");
    const auto n_proc = comm_info_.getNumberOfProc();
    if (!flag_reuse) {
        ReallocatableArray<TopGeometry<TreeCellLoc, TSM>> top_geo(n_proc, n_proc, MemoryAllocMode::Pool);
        TopGeometry<TreeCellLoc, TSM> my_top_geo;
        my_top_geo.set(inner_boundary_of_local_tree_, tc_loc_[0]);
        BarrierForMeasure(comm_info_);
        t0 = GetWtime();
        comm_info_.allGather(&my_top_geo, 1, top_geo.getPointer());
        time_profile_.let_communication_1 += GetWtime() - t0;
        BarrierForMeasure(comm_info_);
        t0 = GetWtime();
        FindScatterParticleLong<TSM, TreeCellLoc, TreeParticle, Tepj, Tspj, TopGeometry<TreeCellLoc, TSM>, MortonKey<typename TSM::mkey_type>,
                                CALC_DISTANCE_TYPE>(tc_loc_, tp_glb_, epj_sorted_, comm_table_.n_ep_send_, comm_table_.adr_ep_send_, dinfo,
                                                    n_leaf_limit_, lev_leaf_limit_, comm_table_.n_sp_send_, comm_table_.adr_sp_send_,
                                                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_, comm_table_.n_ep_per_image_,
                                                    comm_table_.n_sp_per_image_, top_geo, theta_, i_cut_, lev_pm_cell_, morton_key_);
        time_profile_.construct_let_1 += GetWtime() - t0;
        t0 = GetWtime();
        ExchangeNumberLong(comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_sp_send_, comm_table_.n_sp_recv_, comm_info_);
        time_profile_.let_communication_1 += GetWtime() - t0;
        comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
        comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
        comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
        S32 n_ep_recv_tot_tmp = 0;
        S32 n_sp_recv_tot_tmp = 0;
        PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp))
        for (S32 i = 0; i < n_proc; i++) {
            n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
            n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[i];
        }
        comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        comm_table_.n_sp_recv_tot_ = n_sp_recv_tot_tmp;
    }
    t0 = GetWtime();
#if defined(SERIALIZE_EX_LET)
    PS_COMPILE_TIME_MESSAGE("EXCHANGE LET in SERIALIZE MODE");
    ExchangeLetSerialize<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                          comm_table_.adr_ep_send_, epj_org_, n_loc_tot_, tc_loc_, comm_table_.n_sp_send_, comm_table_.n_sp_recv_,
                                          comm_table_.n_sp_per_image_, comm_table_.adr_sp_send_, spj_org_, comm_table_.shift_per_image_,
                                          comm_table_.n_image_per_proc_, comm_info_, comm_a2a_epj_spj_2d_, n_thread, param);
#else
    PS_COMPILE_TIME_MESSAGE("EXCHANGE LET in NON-SERIALIZE MODE");
    ExchangeLet<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_, comm_table_.adr_ep_send_,
                                 epj_org_, n_loc_tot_, tc_loc_, comm_table_.n_sp_send_, comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
                                 comm_table_.adr_sp_send_, spj_org_, comm_table_.shift_per_image_, comm_table_.n_image_per_proc_, comm_info_);
#endif
    time_profile_.let_communication_1 += GetWtime() - t0;
    //time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    // comm_table_.dump();

#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeLongPMM)
    S32 n_ep_send_loc = 0;
    S32 n_sp_send_loc = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_ep_send_loc += comm_table_.n_ep_send_[i];
        n_sp_send_loc += comm_table_.n_sp_send_[i];
    }
    PS_PRINT_VARIABLE2(n_ep_send_loc, n_sp_send_loc);
    PARTICLE_SIMULATOR_PRINT_MSG("***** END exchangeLocalEssentialTreeLongPMM *****", 0);
#endif

    /*
    if (exchange_let_mode_ == EXCHANGE_LET_A2A) {
        const auto n_proc = comm_info_.getNumberOfProc();
    } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) {
        exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
    } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_FAST) {
        exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
        comm_table_.setPosDomainAllgather(dinfo, pos_root_cell_);
    }
    */
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::exchangeLocalEssentialTreeLongP2P(const DomainInfo& dinfo,
                                                                                                                          const bool flag_reuse) {
    F64 time_offset = GetWtime();
    // const auto n_proc  = Comm::getNumberOfProc();
    const auto n_proc = comm_info_.getNumberOfProc();
    F64ort pos_root_cell = getPosRootCell();
    ReallocatableArray<Tspj> top_spj;
    top_spj.setAllocMode(MemoryAllocMode::Pool);
    ReallocatableArray<TopMoment<Tspj, TreeCellLoc, TSM>> top_moment;
    top_moment.setAllocMode(MemoryAllocMode::Pool);
    Tspj my_top_spj;
    my_top_spj.copyFromMoment(tc_loc_[0].mom_);
    if (flag_reuse) {
        top_spj.resizeNoInitialize(n_proc);
        const F64 time_offset_exchange_top_moment = GetWtime();
        comm_info_.allGather(&my_top_spj, 1, top_spj.getPointer());
        time_profile_.let_communication_1 += GetWtime() - time_offset_exchange_top_moment;
        //time_profile_.make_LET_1st__exchange_top_moment += GetWtime() - time_offset_exchange_top_moment;
    }
    if (!flag_reuse) {
        top_moment.resizeNoInitialize(n_proc);
        TopMoment<Tspj, TreeCellLoc, TSM> my_top_moment;
        my_top_moment.set(my_top_spj, inner_boundary_of_local_tree_, tc_loc_[0]);
        const F64 time_offset_exchange_top_moment = GetWtime();
        // Comm::allGather(&my_top_moment, 1, top_moment.getPointer());
        comm_info_.allGather(&my_top_moment, 1, top_moment.getPointer());
        time_profile_.let_communication_1 += GetWtime() - time_offset_exchange_top_moment;
        //time_profile_.make_LET_1st__exchange_top_moment += GetWtime() - time_offset_exchange_top_moment;
        const F64 time_offset_find_particle = GetWtime();
        FindScatterParticleP2P<TSM, TreeCellLoc, TreeParticle, Tepj, Tspj, TopMoment<Tspj, TreeCellLoc, TSM>, MortonKey<typename TSM::mkey_type>,
                               CALC_DISTANCE_TYPE>(tc_loc_, tp_glb_, epj_sorted_, comm_table_.n_ep_send_, comm_table_.adr_ep_send_, dinfo,
                                                   n_leaf_limit_, comm_table_.n_sp_send_, comm_table_.adr_sp_send_, comm_table_.shift_per_image_,
                                                   comm_table_.n_image_per_proc_, comm_table_.n_ep_per_image_, comm_table_.n_sp_per_image_,
                                                   comm_table_.rank_send_, comm_table_.rank_recv_, top_moment, pos_root_cell, theta_,
                                                   comm_table_.rank_recv_allgather_, morton_key_, exchange_let_mode_);
        time_profile_.construct_let_1 += GetWtime() - time_offset_find_particle;
        //time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
        const F64 time_offset_exchange_n = GetWtime();
        ExchangeNumberLong(comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_sp_send_, comm_table_.n_sp_recv_, comm_table_.rank_send_,
                           comm_table_.rank_recv_, comm_info_);
        time_profile_.let_communication_1 += GetWtime() - time_offset_exchange_n;
        //time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;

        comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
        comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
        comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
        S32 n_ep_recv_tot_tmp = 0;
        S32 n_sp_recv_tot_tmp = 0;
        PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp))
        for (S32 i = 0; i < comm_table_.rank_recv_.size(); i++) {
            S32 rank = comm_table_.rank_recv_[i];
            n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[rank];
            n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[rank];
        }
        comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        comm_table_.n_sp_recv_tot_ = n_sp_recv_tot_tmp;
    }
    //time_profile_.make_LET_1st += GetWtime() - time_offset;
    //F64 wtime_ex_let_comm = 0.0;
    time_offset = GetWtime();
    ExchangeLetP2P<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                    comm_table_.adr_ep_send_, epj_org_, n_loc_tot_, tc_loc_, comm_table_.n_sp_send_, comm_table_.n_sp_recv_,
                                    comm_table_.n_sp_per_image_, comm_table_.adr_sp_send_, spj_org_, comm_table_.shift_per_image_,
                                    comm_table_.n_image_per_proc_, comm_table_.rank_send_, comm_table_.rank_recv_, comm_info_);
    time_profile_.let_communication_1 += GetWtime() - time_offset;
    //time_profile_.exchange_LET_1st__icomm_ptcl += wtime_ex_let_comm;
    const auto spj_offset = spj_org_.size();
    const auto n_write = comm_table_.rank_recv_allgather_.size();
    spj_org_.increaseSize(n_write);
    PS_OMP_PARALLEL {
        const S32 ith = Comm::getThreadNum();
        const S32 n_thread = Comm::getNumberOfThread();
        S32 head, end;
        CalcAdrToSplitData(head, end, ith, n_thread, n_write);
        if (flag_reuse) {
            for (S32 i = head; i < end; i++) {
                spj_org_[i + spj_offset] = top_spj[comm_table_.rank_recv_allgather_[i]];
            }
        } else {
            for (S32 i = head; i < end; i++) {
                spj_org_[i + spj_offset] = top_moment[comm_table_.rank_recv_allgather_[i]].spj_;
            }
        }
    }
    //time_profile_.exchange_LET_1st += GetWtime() - time_offset;

#ifdef PS_DEBUG_EXCHANGE_LET
    // #if 1
    F64 mass = 0.0;
    for (int i = 0; i < epj_org_.size(); i++) {
        mass += epj_org_[i].getCharge();
    }
    for (int i = 0; i < spj_org_.size(); i++) {
        mass += spj_org_[i].getCharge();
    }
    // Comm::barrier();
    comm_info_.barrier();
    if (mass != 1.0) {
        // std::cerr<<"my_rank= "<<Comm::getRank()
        std::cerr << "my_rank= " << comm_info_.getRank() << " mass= " << mass << std::endl;
    }
    // Comm::barrier();
    comm_info_.barrier();
    assert(mass == 1.0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::exchangeLocalEssentialTreeLongCutoffP2P(
    const DomainInfo& dinfo, const bool flag_reuse) {
    F64 time_offset = GetWtime();
    //time_profile_.make_LET_1st__exchange_top_moment = 0.0;
    if (!flag_reuse) {
        const F64 time_offset_find_particle = GetWtime();
        FindScatterParticleLongCutoffP2P<TSM, TreeCellLoc, TreeParticle, Tepj, Tspj, TopMoment<Tspj, TreeCellLoc, TSM>, CALC_DISTANCE_TYPE>(
            tc_loc_, tp_glb_, epj_sorted_, comm_table_.n_ep_send_, comm_table_.adr_ep_send_, dinfo, n_leaf_limit_, comm_table_.n_sp_send_,
            comm_table_.adr_sp_send_, comm_table_.shift_per_image_, comm_table_.n_image_per_proc_, comm_table_.n_ep_per_image_,
            comm_table_.n_sp_per_image_, comm_table_.rank_send_, comm_table_.rank_recv_, theta_, exchange_let_mode_);
        time_profile_.construct_let_1 += GetWtime() - time_offset_find_particle;
        //time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
        const F64 time_offset_exchange_n = GetWtime();
        ExchangeNumberLong(comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_sp_send_, comm_table_.n_sp_recv_, comm_table_.rank_send_,
                           comm_table_.rank_recv_, comm_info_);
        time_profile_.let_communication_1 += GetWtime() - time_offset_exchange_n;
        //time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
        comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
        comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
        comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
        S32 n_ep_recv_tot_tmp = 0;
        S32 n_sp_recv_tot_tmp = 0;
        PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp))
        for (S32 i = 0; i < comm_table_.rank_recv_.size(); i++) {
            S32 rank = comm_table_.rank_recv_[i];
            n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[rank];
            n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[rank];
        }
        comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        comm_table_.n_sp_recv_tot_ = n_sp_recv_tot_tmp;
    }
    //time_profile_.make_LET_1st += GetWtime() - time_offset;
    F64 wtime_ex_let_comm = 0.0;
    time_offset = GetWtime();
    ExchangeLetP2P<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                    comm_table_.adr_ep_send_, epj_org_, n_loc_tot_, tc_loc_, comm_table_.n_sp_send_, comm_table_.n_sp_recv_,
                                    comm_table_.n_sp_per_image_, comm_table_.adr_sp_send_, spj_org_, comm_table_.shift_per_image_,
                                    comm_table_.n_image_per_proc_, comm_table_.rank_send_, comm_table_.rank_recv_, wtime_ex_let_comm);
    time_profile_.let_communication_1 += GetWtime() - time_offset;
    //time_profile_.exchange_LET_1st__icomm_ptcl += wtime_ex_let_comm;
    //time_profile_.exchange_LET_1st += GetWtime() - time_offset;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::exchangeLocalEssentialTreeShortScatter(
    const DomainInfo& dinfo, const bool flag_reuse) {
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShortScatter)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START exchangeLocalEssentialTreeShortScatter ***", 0);
#endif
    auto time_offset = GetWtime();
    if (!flag_reuse) {
        const F64 time_offset_find_particle = GetWtime();
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShortScatter)
        PARTICLE_SIMULATOR_PRINT_MSG("*** beforce FindScatterParticleShort ***", 0);
#endif
        FindScatterParticleShort<TreeCellLoc, TreeParticle, Tepj, Tspj, CALC_DISTANCE_TYPE>(
            tc_loc_, tp_glb_, epj_sorted_, comm_table_.n_ep_send_, comm_table_.adr_ep_send_, dinfo, n_leaf_limit_, comm_table_.shift_per_image_,
            comm_table_.n_image_per_proc_, comm_table_.n_ep_per_image_);
        time_profile_.construct_let_1 += GetWtime() - time_offset_find_particle;
        //time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
        const F64 time_offset_exchange_n = GetWtime();

#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShortScatter)
        PARTICLE_SIMULATOR_PRINT_MSG("*** beforce ExchangeNumberShort ***", 0);
#endif
        ExchangeNumberShort(comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_info_);
        time_profile_.let_communication_1 += GetWtime() - time_offset_exchange_n;
        //time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
        // const auto n_proc = Comm::getNumberOfProc();
        const auto n_proc = comm_info_.getNumberOfProc();
        comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
        comm_table_.n_ep_recv_tot_ = 0;
        S32 n_ep_recv_tot_tmp = 0;
        PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp))
        for (S32 i = 0; i < n_proc; i++) {
            n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
        }
        comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
    }
    //time_profile_.make_LET_1st += GetWtime() - time_offset;

    time_offset = GetWtime();
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShortScatter)
    /*
        if (Comm::getRank() == 0) {
            comm_table_.dump();
            for(S32 i=0; i<comm_table_.adr_ep_send_.size(); i++){
                std::cerr<<"i= "<<i<<" adr_ep_send= "<<comm_table_.adr_ep_send_[i]<<std::endl;
            }
            for(S32 i=0; i<comm_table_.n_ep_send_.size(); i++){
                std::cerr<<"i= "<<i<<" n_ep_send= "<<comm_table_.n_ep_send_[i]<<std::endl;
            }
            for(S32 i=0; i<comm_table_.n_ep_per_image_.size(); i++){
                std::cerr<<"i= "<<i<<" n_ep_per_image= "<<comm_table_.n_ep_per_image_[i]<<std::endl;
            }
            for(S32 i=0; i<comm_table_.shift_per_image_.size(); i++){
                std::cerr<<"i= "<<i<<" shift_per_image= "<<comm_table_.shift_per_image_[i]<<std::endl;
            }
            for(S32 i=0; i<comm_table_.n_image_per_proc_.size(); i++){
                std::cerr<<"i= "<<i<<" n_image_per_proc= "<<comm_table_.n_image_per_proc_[i]<<std::endl;
            }
        }
    */
    PARTICLE_SIMULATOR_PRINT_MSG("*** beforce ExchangeLetShortA2A ***", 0);
#endif
    ExchangeLetShortA2A<TSM, Tepj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_, comm_table_.adr_ep_send_,
                                   epj_org_, n_loc_tot_, comm_table_.shift_per_image_, comm_table_.n_image_per_proc_, comm_info_);
    time_profile_.let_communication_1 += GetWtime() - time_offset;
    //time_profile_.exchange_LET_1st += GetWtime() - time_offset;
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShortScatter)
    PARTICLE_SIMULATOR_PRINT_MSG("*** END exchangeLocalEssentialTreeShortScatter ***", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::exchangeLocalEssentialTreeShortScatterP2P(
    const DomainInfo& dinfo, const bool flag_reuse) {
    auto time_offset = GetWtime();
    if (!flag_reuse) {
        const F64 time_offset_find_particle = GetWtime();
        FindScatterParticleP2P<TreeCellLoc, TreeParticle, Tepj, Tspj, CALC_DISTANCE_TYPE>(
            tc_loc_, tp_glb_, epj_sorted_, comm_table_.n_ep_send_, comm_table_.adr_ep_send_, dinfo, n_leaf_limit_, comm_table_.shift_per_image_,
            comm_table_.n_image_per_proc_, comm_table_.n_ep_per_image_, comm_table_.rank_send_, comm_table_.rank_recv_);
        time_profile_.construct_let_1 += GetWtime() - time_offset_find_particle;
        //time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
        const F64 time_offset_exchange_n = GetWtime();
        ExchangeNumberShort(comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.rank_send_, comm_table_.rank_recv_, comm_info_);
        time_profile_.let_communication_1 += GetWtime() - time_offset_exchange_n;
        //time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
        // const auto n_proc = Comm::getNumberOfProc();
        const auto n_proc = comm_info_.getNumberOfProc();
        comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
        comm_table_.n_ep_recv_tot_ = 0;
        S32 n_ep_recv_tot_tmp = 0;
        PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp))
        for (S32 i = 0; i < n_proc; i++) {
            n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
        }
        comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
    }
    //time_profile_.make_LET_1st += GetWtime() - time_offset;

    time_offset = GetWtime();
    ExchangeLetShortP2P<TSM, Tepj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_, comm_table_.adr_ep_send_,
                                   epj_org_, n_loc_tot_, comm_table_.shift_per_image_, comm_table_.n_image_per_proc_, comm_table_.rank_send_,
                                   comm_table_.rank_recv_, comm_info_);
    time_profile_.let_communication_1 += GetWtime() - time_offset;
    //time_profile_.exchange_LET_1st += GetWtime() - time_offset;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::exchangeLocalEssentialTreeShort2StepA2A(
    const DomainInfo& dinfo, const bool flag_reuse) {
    // PS_COMPILE_TIME_MESSAGE("exchangeLocalEssentialTreeShort2StepA2A");
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShort2StepA2A)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START exchangeLocalEssentialTreeShort2StepA2A *****", 0);
#endif
    auto time_offset = GetWtime();
    ReallocatableArray<S32> n_ep_send_per_proc_1st;
    ReallocatableArray<S32> n_ep_send_per_proc_2nd;
    ReallocatableArray<S32> n_ep_recv_per_proc_1st;
    ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
    ReallocatableArray<S32> adr_ep_send_1st;
    ReallocatableArray<S32> adr_ep_send_2nd;
    ReallocatableArray<F64vec> shift_per_image_1st;
    ReallocatableArray<F64vec> shift_per_image_2nd;
    ReallocatableArray<S32> n_image_per_proc_1st;
    ReallocatableArray<S32> n_image_per_proc_2nd;
    ReallocatableArray<S32> n_ep_send_per_image_1st;
    ReallocatableArray<S32> n_ep_send_per_image_2nd;
    ReallocatableArray<Tepj> epj_recv_1st;
    ReallocatableArray<EssentialParticleBase> epi_base_sorted;    // for gather mode
    ReallocatableArray<EssentialParticleBase> epi_base_recv_1st;  // for gather mode
    if constexpr (std::is_same_v<typename TSM::neighbor_search_type, TagNeighborSearchGather>) {
        // ReallocatableArray<EssentialParticleBase> ep_recv_1st;
        const S32 n_epi_sorted = epi_sorted_.size();
        epi_base_sorted.resizeNoInitialize(n_epi_sorted);
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_epi_sorted; i++) {
            epi_base_sorted[i].pos = epi_sorted_[i].getPos();
            epi_base_sorted[i].r_search = epi_sorted_[i].getRSearch();
        }
    } else {
        // ReallocatableArray<Tepj> ep_recv_1st;
    }
    if (!flag_reuse) {
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShort2StepA2A)
        PARTICLE_SIMULATOR_PRINT_MSG("- before FindScatterParticle", 0);
#endif
        ////////////
        // 1st STEP (send epi_base)
        auto time_offset_find_particle = GetWtime();
        if constexpr (std::is_same_v<typename TSM::neighbor_search_type, TagNeighborSearchGather>) {
            FindScatterParticleShort<TreeCellLoc, TreeParticle, EssentialParticleBase, Tspj, CALC_DISTANCE_TYPE>(
                tc_loc_, tp_glb_, epi_base_sorted, n_ep_send_per_proc_1st, adr_ep_send_1st, dinfo, n_leaf_limit_, shift_per_image_1st,
                n_image_per_proc_1st, n_ep_send_per_image_1st);
        } else {
            FindScatterParticleShort<TreeCellLoc, TreeParticle, Tepj, Tspj, CALC_DISTANCE_TYPE>(
                tc_loc_, tp_glb_, epj_sorted_, n_ep_send_per_proc_1st, adr_ep_send_1st, dinfo, n_leaf_limit_, shift_per_image_1st,
                n_image_per_proc_1st, n_ep_send_per_image_1st);
        }
        time_profile_.construct_let_1 += GetWtime() - time_offset_find_particle;
        //time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;

#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShort2StepA2A)
        PARTICLE_SIMULATOR_PRINT_MSG("- before ExchangeNumberShort", 0);
#endif
        auto time_offset_exchange_n = GetWtime();
        ExchangeNumberShort(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st, comm_info_);
        time_profile_.let_communication_1 += GetWtime() - time_offset_exchange_n;
        //time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
        auto t0 = GetWtime();
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShort2StepA2A)
        PARTICLE_SIMULATOR_PRINT_MSG("- before ExchangeLET", 0);
#endif
        if constexpr (std::is_same_v<typename TSM::neighbor_search_type, TagNeighborSearchGather>) {
            ExchangeLetShortA2A<TSM, EssentialParticleBase>(epi_base_sorted, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st, n_ep_send_per_image_1st,
                                                            adr_ep_send_1st, epi_base_recv_1st, 0, shift_per_image_1st, n_image_per_proc_1st,
                                                            comm_info_);
        } else {
            ExchangeLetShortA2A<TSM, Tepj>(epj_sorted_, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st, n_ep_send_per_image_1st, adr_ep_send_1st,
                                           epj_recv_1st, 0, shift_per_image_1st, n_image_per_proc_1st, comm_info_);
        }
        time_profile_.let_communication_1 += GetWtime() - t0;
        ////////////
        // 2nd STEP (find j particle)
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShort2StepA2A)
        PARTICLE_SIMULATOR_PRINT_MSG("- before FindExchangeParticleDoubleWalk", 0);
#endif
        t0 = GetWtime();
        if constexpr (std::is_same_v<typename TSM::neighbor_search_type, TagNeighborSearchGather>) {
            FindExchangeParticleDoubleWalk<TreeCellLoc, TreeParticle, EssentialParticleBase>(
                epi_base_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo, n_leaf_limit_, n_ep_send_per_proc_2nd,
                n_ep_send_per_image_2nd, n_image_per_proc_2nd, adr_ep_send_2nd, shift_per_image_2nd, epi_base_sorted, center_, length_, morton_key_);
        } else {
            FindExchangeParticleDoubleWalk<TreeCellLoc, TreeParticle, Tepj>(
                epj_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo, n_leaf_limit_, n_ep_send_per_proc_2nd,
                n_ep_send_per_image_2nd, n_image_per_proc_2nd, adr_ep_send_2nd, shift_per_image_2nd, epj_sorted_, center_, length_, morton_key_);
        }
        time_profile_.construct_let_2 += GetWtime() - t0;
        /////////////////////
        // 3rd STEP (exchange # of particles again)
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShort2StepA2A)
        PARTICLE_SIMULATOR_PRINT_MSG("- before ExchangeNumberShort", 0);
#endif
        t0 = GetWtime();
        ExchangeNumberShort(n_ep_send_per_proc_2nd, n_ep_recv_per_proc_2nd, comm_info_);
        time_profile_.let_communication_2 += GetWtime() - t0;
        /////////////////////
        // 4th STEP (make communication table)
        const auto n_proc = comm_info_.getNumberOfProc();
        t0 = GetWtime();
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShort2StepA2A)
        PARTICLE_SIMULATOR_PRINT_MSG("- before MakeCommTableFor2StepCommuniction", 0);
#endif
        MakeCommTableFor2StepCommuniction(comm_table_, n_ep_send_per_proc_1st, n_image_per_proc_1st, shift_per_image_1st, n_ep_send_per_image_1st,
                                          adr_ep_send_1st, n_ep_recv_per_proc_1st, n_ep_send_per_proc_2nd, n_image_per_proc_2nd, shift_per_image_2nd,
                                          n_ep_send_per_image_2nd, adr_ep_send_2nd, n_ep_recv_per_proc_2nd, n_proc);
        time_profile_.let_communication_2 += GetWtime() - t0;
    }  // end of reuse flag
    //time_profile_.make_LET_1st += GetWtime() - time_offset;

#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShort2StepA2A)
    PARTICLE_SIMULATOR_PRINT_MSG("- before ExchangeLet", 0);
#endif
    time_offset = GetWtime();
    ExchangeLetShortA2A<TSM, Tepj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_, comm_table_.adr_ep_send_,
                                   epj_org_, n_loc_tot_, comm_table_.shift_per_image_, comm_table_.n_image_per_proc_, comm_info_);
    time_profile_.let_communication_2 += GetWtime() - time_offset;
    //time_profile_.exchange_LET_1st += GetWtime() - time_offset;
#if defined(PS_DEBUG_PRINT_exchangeLocalEssentialTreeShortGather)
    PARTICLE_SIMULATOR_PRINT_MSG("***** END exchangeLocalEssentialTreeShortGather *****", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::exchangeLocalEssentialTreeShort2StepP2P(
    const DomainInfo& dinfo, const bool flag_reuse) {
    // PS_COMPILE_TIME_MESSAGE("exchangeLocalEssentialTreeShort2StepP2P");
    F64 time_offset = GetWtime();
    ReallocatableArray<S32> n_ep_send_per_proc_1st;
    ReallocatableArray<S32> n_ep_send_per_proc_2nd;
    ReallocatableArray<S32> n_ep_recv_per_proc_1st;
    ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
    ReallocatableArray<S32> adr_ep_send_1st;
    ReallocatableArray<S32> adr_ep_send_2nd;
    ReallocatableArray<F64vec> shift_per_image_1st;
    ReallocatableArray<F64vec> shift_per_image_2nd;
    ReallocatableArray<S32> n_image_per_proc_1st;
    ReallocatableArray<S32> n_image_per_proc_2nd;
    ReallocatableArray<S32> n_ep_send_per_image_1st;
    ReallocatableArray<S32> n_ep_send_per_image_2nd;
    ReallocatableArray<Tepj> epj_recv_1st;                        // symmetry
    ReallocatableArray<EssentialParticleBase> epi_base_recv_1st;  // gather
    ReallocatableArray<EssentialParticleBase> epi_base_sorted;    // gather
    if constexpr (std::is_same_v<typename TSM::neighbor_search_type, TagNeighborSearchGather>) {
        const S32 n_epi_sorted = epi_sorted_.size();
        epi_base_sorted.resizeNoInitialize(n_epi_sorted);
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_epi_sorted; i++) {
            epi_base_sorted[i].pos = epi_sorted_[i].getPos();
            epi_base_sorted[i].r_search = epi_sorted_[i].getRSearch();
        }
    }
    if (!flag_reuse) {
        ////////////
        // 1st STEP (send epi_base)
        const F64 time_offset_find_particle = GetWtime();
        FindScatterParticleP2P<TreeCellLoc, TreeParticle, EssentialParticleBase, Tspj, CALC_DISTANCE_TYPE>(
            tc_loc_, tp_glb_, epi_base_sorted, n_ep_send_per_proc_1st, adr_ep_send_1st, dinfo, n_leaf_limit_, shift_per_image_1st,
            n_image_per_proc_1st, n_ep_send_per_image_1st, comm_table_.rank_send_, comm_table_.rank_recv_);
        time_profile_.construct_let_1 += GetWtime() - time_offset_find_particle;
        //time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
        const F64 time_offset_exchange_n = GetWtime();
        ExchangeNumberShort(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st, comm_table_.rank_send_, comm_table_.rank_recv_, comm_info_);
        time_profile_.let_communication_1 += GetWtime() - time_offset_exchange_n;
        //time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
        auto t0 = GetWtime();
        if constexpr (std::is_same_v<typename TSM::neighbor_search_type, TagNeighborSearchGather>) {
            ExchangeLetShortP2P<TSM, EssentialParticleBase>(epi_base_sorted, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st, n_ep_send_per_image_1st,
                                                            adr_ep_send_1st, epi_base_recv_1st, 0, shift_per_image_1st, n_image_per_proc_1st,
                                                            comm_table_.rank_send_, comm_table_.rank_recv_, comm_info_);
        } else {
            ExchangeLetShortP2P<TSM, Tepj>(epj_sorted_, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st, n_ep_send_per_image_1st, adr_ep_send_1st,
                                           epj_recv_1st, 0, shift_per_image_1st, n_image_per_proc_1st, comm_table_.rank_send_, comm_table_.rank_recv_,
                                           comm_info_);
        }
        time_profile_.let_communication_1 += GetWtime() - t0;
        ////////////
        // 2nd STEP (find j particle)
        t0 = GetWtime();
        if constexpr (std::is_same_v<typename TSM::neighbor_search_type, TagNeighborSearchGather>) {
            FindExchangeParticleDoubleWalk<TreeCellLoc, TreeParticle, EssentialParticleBase>(
                epi_base_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo, n_leaf_limit_, n_ep_send_per_proc_2nd,
                n_ep_send_per_image_2nd, n_image_per_proc_2nd, adr_ep_send_2nd, shift_per_image_2nd, epi_base_sorted, center_, length_, morton_key_);
        } else {
            FindExchangeParticleDoubleWalk<TreeCellLoc, TreeParticle, Tepj>(
                epj_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo, n_leaf_limit_, n_ep_send_per_proc_2nd,
                n_ep_send_per_image_2nd, n_image_per_proc_2nd, adr_ep_send_2nd, shift_per_image_2nd, epj_sorted_, center_, length_, morton_key_);
        }
        time_profile_.construct_let_2 += GetWtime() - t0;
        /////////////////////
        // 3rd STEP (exchange # of particles again)
        t0 = GetWtime();
        ExchangeNumberShort(n_ep_send_per_proc_2nd, n_ep_recv_per_proc_2nd, comm_table_.rank_send_, comm_table_.rank_recv_, comm_info_);
        time_profile_.construct_let_2 += GetWtime() - t0;
        /////////////////////
        // 4th STEP (make communication table)
        const auto n_proc = comm_info_.getNumberOfProc();
        t0 = GetWtime();
        MakeCommTableFor2StepCommuniction(comm_table_, n_ep_send_per_proc_1st, n_image_per_proc_1st, shift_per_image_1st, n_ep_send_per_image_1st,
                                          adr_ep_send_1st, n_ep_recv_per_proc_1st, n_ep_send_per_proc_2nd, n_image_per_proc_2nd, shift_per_image_2nd,
                                          n_ep_send_per_image_2nd, adr_ep_send_2nd, n_ep_recv_per_proc_2nd, n_proc);
        time_profile_.let_communication_2 += GetWtime() - time_offset;                                          
    }  // end of reuse flag
    //time_profile_.make_LET_1st += GetWtime() - time_offset;
    time_offset = GetWtime();
    ExchangeLetShortP2P<TSM, Tepj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_, comm_table_.adr_ep_send_,
                                   epj_org_, n_loc_tot_, comm_table_.shift_per_image_, comm_table_.n_image_per_proc_, comm_table_.rank_send_,
                                   comm_table_.rank_recv_, comm_info_);
    time_profile_.let_communication_2 += GetWtime() - time_offset;
    //time_profile_.exchange_LET_1st += GetWtime() - time_offset;
}

}  // namespace ParticleSimulator