#pragma once

extern int N_STEP;

namespace ParticleSimulator {

template <typename Tmorton_key>
inline F64ort GetCorrespondingTreeCell(const F64ort &box, const Tmorton_key &morton_key) {
    return morton_key.getCorrespondingTreeCell(box);
}

///////////////////////////////
// small functions for dispatch
template <class Ttc>
inline F64ort GetOuterBoundaryOfMyTree(TagGeometrySize, const ReallocatableArray<Ttc> &tc_first) {
    PS_COMPILE_TIME_MESSAGE("TagGeometrySize@GetOuterBoundaryOfMyTree")
    return F64ort(-1234.5, -9876.5);
}
template <class Ttc>
inline F64ort GetOuterBoundaryOfMyTree(TagGeometrySizeOut, const ReallocatableArray<Ttc> &tc_first) {
    PS_COMPILE_TIME_MESSAGE("TagGeometrySizeOut@GetOuterBoundaryOfMyTree")
    return tc_first[0].geo_.getVertexOut();
}
template <class Ttc>
inline F64ort GetOuterBoundaryOfMyTree(TagGeometrySizeInOut, const ReallocatableArray<Ttc> &tc_first) {
    PS_COMPILE_TIME_MESSAGE("TagGeometrySizeInOut@GetOuterBoundaryOfMyTree")
    return tc_first[0].geo_.getVertexOut();
}
template <class Ttc>
inline F64ort GetOuterBoundaryOfMyTree(TagGeometryInAndOut, const ReallocatableArray<Ttc> &tc_first) {
    PS_COMPILE_TIME_MESSAGE("TagGeometrySizeInAndOut@GetOuterBoundaryOfMyTree")
    return tc_first[0].geo_.getVertexOut();
}
template <class Ttc>
inline F64ort GetOuterBoundaryOfMyTree(TagGeometrySizeIn,  // for PMMM
                                       const ReallocatableArray<Ttc> &tc_first) {
    PS_COMPILE_TIME_MESSAGE("TagGeometrySizeIn@GetOuterBoundaryOfMyTree")
    return tc_first[0].geo_.getVertexIn();
}

////////////////////
// for exchange LET
// only for open boundary
template <class Tptcl, class Ttc>
inline void CopyPtclFromTreeToSendBuf(TagSearchBoundaryConditionOpenOnly, ReallocatableArray<Tptcl> &ptcl_send, const ReallocatableArray<Ttc> &tc,
                                      const ReallocatableArray<S32> &adr_ptcl_send, const S32 n_ptcl, const S32 n_ptcl_offset, const F64vec &shift) {
    for (S32 j = 0; j < n_ptcl; j++) {
        S32 adr = adr_ptcl_send[n_ptcl_offset + j];
        ptcl_send[n_ptcl_offset + j].copyFromMoment(tc[adr].mom_);
    }
}

// for periodic boundary (not only for periodic but also open boudnary)
template <class Tptcl, class Ttc>
inline void CopyPtclFromTreeToSendBuf(TagSearchBoundaryConditionOpenPeriodic, ReallocatableArray<Tptcl> &ptcl_send, const ReallocatableArray<Ttc> &tc,
                                      const ReallocatableArray<S32> &adr_ptcl_send, const S32 n_ptcl, const S32 n_ptcl_offset, const F64vec &shift) {
    for (S32 j = 0; j < n_ptcl; j++) {
        S32 adr = adr_ptcl_send[n_ptcl_offset + j];
        ptcl_send[n_ptcl_offset + j].copyFromMoment(tc[adr].mom_);
        const F64vec pos_new = ptcl_send[n_ptcl_offset + j].getPos() - shift;
        ptcl_send[n_ptcl_offset + j].setPos(pos_new);
    }
}
template <class Tptcl>
inline void CopyPtclToSendBuf(TagSearchBoundaryConditionOpenOnly, ReallocatableArray<Tptcl> &ptcl_send, const ReallocatableArray<Tptcl> &ptcl,
                              const ReallocatableArray<S32> &adr_ptcl_send, const S32 n_ptcl, const S32 n_ptcl_offset, const F64vec &shift) {
    for (S32 j = 0; j < n_ptcl; j++) {
        S32 adr = adr_ptcl_send[n_ptcl_offset + j];
        ptcl_send[n_ptcl_offset + j] = ptcl[adr];
    }
}
template <class Tptcl>
inline void CopyPtclToSendBuf(TagSearchBoundaryConditionOpenPeriodic, ReallocatableArray<Tptcl> &ptcl_send, const ReallocatableArray<Tptcl> &ptcl,
                              const ReallocatableArray<S32> &adr_ptcl_send, const S32 n_ptcl, const S32 n_ptcl_offset, const F64vec &shift) {
    for (S32 j = 0; j < n_ptcl; j++) {
        S32 adr = adr_ptcl_send[n_ptcl_offset + j];
        ptcl_send[n_ptcl_offset + j] = ptcl[adr];
        ptcl_send[n_ptcl_offset + j].setPos(ptcl[adr].getPos() - shift);
    }
}
// for exchange LET
////////////////////

// small functions for dispatch
///////////////////////////////

////////////////
// for long mode
template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, typename Tmoment, typename Tmorton_key, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
inline void FindScatterParticleP2P(const ReallocatableArray<Ttc> &tc_first, const ReallocatableArray<Ttp> &tp_first,
                                   const ReallocatableArray<Tep> &ep_first, ReallocatableArray<S32> &n_ep_send, ReallocatableArray<S32> &adr_ep_send,
                                   const DomainInfo &dinfo, const S32 n_leaf_limit, ReallocatableArray<S32> &n_sp_send,
                                   ReallocatableArray<S32> &adr_sp_send, ReallocatableArray<F64vec> &shift_per_image,
                                   ReallocatableArray<S32> &n_image_per_proc, ReallocatableArray<S32> &n_ep_per_image,
                                   ReallocatableArray<S32> &n_sp_per_image, ReallocatableArray<S32> &rank_send, ReallocatableArray<S32> &rank_recv,
                                   ReallocatableArray<Tmoment> &top_moment, const F64ort &pos_root_cell, const F64 theta,
                                   ReallocatableArray<S32> &rank_recv_a2a, const Tmorton_key &morton_key,
                                   const enum EXCHANGE_LET_MODE exchange_let_mode_) {
    PS_COMPILE_TIME_MESSAGE("FindScatterParticleP2P(long)");
    const auto n_thread = Comm::getNumberOfThread();
    const auto n_proc = dinfo.getCommInfo().getNumberOfProc();
    const auto my_rank = dinfo.getCommInfo().getRank();
    std::vector<ReallocatableArray<S32>> rank_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<F64vec>> shift_per_image_tmp(n_thread, ReallocatableArray<F64vec>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> adr_ep_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> n_ep_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> adr_sp_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> n_sp_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> rank_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> rank_recv_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> rank_recv_allgather_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    ReallocatableArray<F64ort> minimum_tree_boundary_of_domain(n_proc, n_proc, MemoryAllocMode::Pool);
    n_ep_send.resizeNoInitialize(n_proc);
    n_sp_send.resizeNoInitialize(n_proc);
    n_image_per_proc.resizeNoInitialize(n_proc);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        n_ep_send[i] = n_sp_send[i] = n_image_per_proc[i] = 0;
    }
    bool pa[DIMENSION];
    dinfo.getPeriodicAxis(pa);
    F64ort pos_root_domain = dinfo.getPosRootDomain();
    F64vec len_peri = pos_root_domain.getFullLength();
    F64ort outer_boundary_of_my_tree = GetOuterBoundaryOfMyTree(typename TSM::tree_cell_loc_geometry_type(), tc_first);
    // exchange outer boundary
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        minimum_tree_boundary_of_domain[i] = GetCorrespondingTreeCell(top_moment[i].vertex_in_, morton_key);
    }
    F64 r_crit_send = 0.0;
    if (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) {
        r_crit_send = minimum_tree_boundary_of_domain[my_rank].getFullLength().getMax();
    } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_FAST) {
        r_crit_send = top_moment[my_rank].vertex_in_.getFullLength().getMax();
    } else {
        assert(0);
    }
    if (theta > 0.0) {
        r_crit_send = r_crit_send * r_crit_send / (theta * theta);
    } else {
        r_crit_send = -1.0;
    }
    PS_OMP_PARALLEL {
        const S32 ith = Comm::getThreadNum();
        const S32 n_thread = Comm::getNumberOfThread();
        const S32 head = (n_proc / n_thread) * ith + std::min(n_proc % n_thread, ith);
        const S32 end = (n_proc / n_thread) * (ith + 1) + std::min(n_proc % n_thread, (ith + 1));
        rank_send_tmp[ith].clearSize();
        rank_recv_tmp[ith].clearSize();
        rank_recv_allgather_tmp[ith].clearSize();
        for (S32 ib = head; ib < end; ib++) {
            const S32 id0 = std::min(ib, my_rank);
            const S32 id1 = std::max(ib, my_rank);
            bool flag_recv_allgather = true;
            if (((typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER) || typeid(TSM) == typeid(SEARCH_MODE_LONG_SYMMETRY)) &&
                 (top_moment[id0].isOverlapped(top_moment[id1]))) ||
                (theta <= 0.0)) {
                rank_send_tmp[ith].push_back(ib);
                rank_recv_tmp[ith].push_back(ib);
                flag_recv_allgather = false;
            } else {
                if ((GetDistanceMinSq<CALC_DISTANCE_TYPE>(top_moment[ib].vertex_in_, top_moment[my_rank].spj_.getPos(), len_peri) <= r_crit_send) ||
                    ((GetDistanceMinSq<CALC_DISTANCE_TYPE>(top_moment[ib].vertex_in_, minimum_tree_boundary_of_domain[my_rank], len_peri) <= 0.0) &&
                     (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT))) {
                    rank_send_tmp[ith].push_back(ib);
                }
                F64 r_crit_recv = 0.0;
                if (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) {
                    r_crit_recv = minimum_tree_boundary_of_domain[ib].getFullLength().getMax();
                } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_FAST) {
                    r_crit_recv = top_moment[ib].vertex_in_.getFullLength().getMax();
                } else {
                    assert(0);
                }
                if (theta > 0.0) {
                    r_crit_recv = r_crit_recv * r_crit_recv / (theta * theta);
                } else {
                    r_crit_recv = -1.0;
                }
                if ((GetDistanceMinSq<CALC_DISTANCE_TYPE>(top_moment[my_rank].vertex_in_, top_moment[ib].spj_.getPos(), len_peri) <= r_crit_recv) ||
                    ((GetDistanceMinSq<CALC_DISTANCE_TYPE>(top_moment[my_rank].vertex_in_, minimum_tree_boundary_of_domain[ib], len_peri) <= 0.0) &&
                     (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT))) {
                    rank_recv_tmp[ith].push_back(ib);
                    flag_recv_allgather = false;
                }
            }
            if (flag_recv_allgather) {
                rank_recv_allgather_tmp[ith].push_back(ib);
            }
        }
    }  // end of OMP scope
    // Comm::barrier();
    // exit(1);
    PackData(&rank_send_tmp[0], Comm::getNumberOfThread(), rank_send);
    PackData(&rank_recv_tmp[0], Comm::getNumberOfThread(), rank_recv);
    PackData(&rank_recv_allgather_tmp[0], Comm::getNumberOfThread(), rank_recv_a2a);

    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        n_image_per_proc[i] = 0;
    }
    S32 adr_tc = 0;
    S32 adr_tree_sp_first = 0;
    PS_OMP_PARALLEL {
        S32 ith = Comm::getThreadNum();
        rank_tmp[ith].clearSize();
        shift_per_image_tmp[ith].clearSize();
        adr_ep_send_tmp[ith].clearSize();
        n_ep_per_image_tmp[ith].clearSize();
        adr_sp_send_tmp[ith].clearSize();
        n_sp_per_image_tmp[ith].clearSize();
        ReallocatableArray<Tsp> sp_first;
PS_OMP(omp for schedule(dynamic, 4))
for (S32 ib = 0; ib < rank_send.size(); ib++) {
    const S32 rank = rank_send[ib];
    const S32 n_image_tmp_prev = shift_per_image_tmp[ith].size();
    rank_tmp[ith].push_back(rank);
    CalcNumberAndShiftOfImageDomain(shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(), outer_boundary_of_my_tree,
                                    dinfo.getPosDomain(rank), pa, false);
    const S32 n_image_tmp = shift_per_image_tmp[ith].size();
    n_image_per_proc[rank] = n_image_tmp - n_image_tmp_prev;
    S32 n_ep_prev = adr_ep_send_tmp[ith].size();
    S32 n_sp_prev = adr_sp_send_tmp[ith].size();
    for (S32 j = n_image_tmp_prev; j < n_image_tmp; j++) {
        S32 n_ep_prev_2 = adr_ep_send_tmp[ith].size();
        S32 n_sp_prev_2 = adr_sp_send_tmp[ith].size();
        if (my_rank == rank && j == n_image_tmp_prev) {
            // self image
            n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
            n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
            continue;
        }
        TargetBox<TSM> target_box;
        target_box.setForExLet(top_moment[rank], shift_per_image_tmp[ith][j]);
        MakeListUsingTreeRecursiveTop<TSM, Ttc, TreeParticle, Tep, Tsp, TagChopLeafFalse, TagCopyInfoCloseNoSp, CALC_DISTANCE_TYPE>(
            tc_first, adr_tc, tp_first, ep_first, adr_ep_send_tmp[ith], sp_first, adr_sp_send_tmp[ith], target_box, n_leaf_limit, adr_tree_sp_first,
            len_peri, theta);
        n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
        n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
    }
    n_ep_send[rank] = adr_ep_send_tmp[ith].size() - n_ep_prev;
    n_sp_send[rank] = adr_sp_send_tmp[ith].size() - n_sp_prev;
}
    }  // end of OMP scope
    ReallocatableArray<S32> n_disp_image_per_proc(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    n_disp_image_per_proc[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_image_per_proc[i + 1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
    }
    const S32 n_image_tot = n_disp_image_per_proc[n_proc];
    shift_per_image.resizeNoInitialize(n_image_tot);
    n_ep_per_image.resizeNoInitialize(n_image_tot);
    n_sp_per_image.resizeNoInitialize(n_image_tot);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_thread; i++) {
        S32 n_cnt = 0;
        for (S32 j = 0; j < rank_tmp[i].size(); j++) {
            S32 rank = rank_tmp[i][j];
            S32 offset = n_disp_image_per_proc[rank];
            for (S32 k = 0; k < n_image_per_proc[rank]; k++) {
                shift_per_image[offset + k] = shift_per_image_tmp[i][n_cnt];
                n_ep_per_image[offset + k] = n_ep_per_image_tmp[i][n_cnt];
                n_sp_per_image[offset + k] = n_sp_per_image_tmp[i][n_cnt];
                n_cnt++;
            }
        }
    }
    ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    n_disp_ep_per_image[0] = 0;
    n_disp_sp_per_image[0] = 0;
    for (S32 i = 0; i < n_image_tot; i++) {
        n_disp_ep_per_image[i + 1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
        n_disp_sp_per_image[i + 1] = n_disp_sp_per_image[i] + n_sp_per_image[i];
    }
    const S32 n_ep_send_tot = n_disp_ep_per_image[n_image_tot];
    const S32 n_sp_send_tot = n_disp_sp_per_image[n_image_tot];
    adr_ep_send.resizeNoInitialize(n_ep_send_tot);
    adr_sp_send.resizeNoInitialize(n_sp_send_tot);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_thread; i++) {
        S32 n_cnt_ep = 0;
        S32 n_cnt_sp = 0;
        for (S32 j = 0; j < rank_tmp[i].size(); j++) {
            S32 rank = rank_tmp[i][j];
            const S32 adr_image_head = n_disp_image_per_proc[rank];
            const S32 adr_image_end = n_disp_image_per_proc[rank + 1];
            for (S32 k = adr_image_head; k < adr_image_end; k++) {
                const S32 adr_ep_head = n_disp_ep_per_image[k];
                const S32 adr_ep_end = n_disp_ep_per_image[k + 1];
                for (S32 l = adr_ep_head; l < adr_ep_end; l++) {
                    adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                }
                const S32 adr_sp_head = n_disp_sp_per_image[k];
                const S32 adr_sp_end = n_disp_sp_per_image[k + 1];
                for (S32 l = adr_sp_head; l < adr_sp_end; l++) {
                    adr_sp_send[l] = adr_sp_send_tmp[i][n_cnt_sp++];
                }
            }
        }
    }
}

template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, typename Tmoment, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
inline void FindScatterParticleLongCutoffP2P(const ReallocatableArray<Ttc> &tc_first, const ReallocatableArray<Ttp> &tp_first,
                                             const ReallocatableArray<Tep> &ep_first, ReallocatableArray<S32> &n_ep_send,
                                             ReallocatableArray<S32> &adr_ep_send, const DomainInfo &dinfo, const S32 n_leaf_limit,
                                             ReallocatableArray<S32> &n_sp_send, ReallocatableArray<S32> &adr_sp_send,
                                             ReallocatableArray<F64vec> &shift_per_image, ReallocatableArray<S32> &n_image_per_proc,
                                             ReallocatableArray<S32> &n_ep_per_image, ReallocatableArray<S32> &n_sp_per_image,
                                             ReallocatableArray<S32> &rank_send, ReallocatableArray<S32> &rank_recv, const F64 theta,
                                             const enum EXCHANGE_LET_MODE exchange_let_mode_) {
    PS_COMPILE_TIME_MESSAGE("FindScatterParticleLongCutoffP2P");
    const auto n_thread = Comm::getNumberOfThread();
    const auto n_proc = dinfo.getCommInfo().getNumberOfProc();
    const auto my_rank = dinfo.getCommInfo().getRank();
    std::vector<ReallocatableArray<S32>> rank_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<F64vec>> shift_per_image_tmp(n_thread, ReallocatableArray<F64vec>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> adr_ep_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> n_ep_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> adr_sp_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> n_sp_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> rank_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> rank_recv_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    const auto pos_root_domain = dinfo.getPosRootDomain();
    const auto len_peri = pos_root_domain.getFullLength();
    PS_OMP_PARALLEL {
        const auto ith = Comm::getThreadNum();
        const auto n_thread = Comm::getNumberOfThread();
        S32 head, end;
        CalcAdrToSplitData(head, end, ith, n_thread, n_proc);
        rank_send_tmp[ith].clearSize();
        rank_recv_tmp[ith].clearSize();
        for (S32 ib = head; ib < end; ib++) {
            const auto id0 = std::min(ib, my_rank);
            const auto id1 = std::max(ib, my_rank);
            const auto vertex0 = dinfo.getPosDomain(id0);
            const auto vertex1 = dinfo.getPosDomain(id1);
            const auto dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex0, vertex1, len_peri);
        }
    }  // end of OMP scope
    PackData(&rank_send_tmp[0], Comm::getNumberOfThread(), rank_send);
    PackData(&rank_recv_tmp[0], Comm::getNumberOfThread(), rank_recv);

#if 0
        std::vector<ReallocatableArray<S32>> rank_recv_allgather_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        ReallocatableArray<F64ort> minimum_tree_boundary_of_domain(n_proc, n_proc, MemoryAllocMode::Pool);
        n_ep_send.resizeNoInitialize(n_proc);
        n_sp_send.resizeNoInitialize(n_proc);
        n_image_per_proc.resizeNoInitialize(n_proc);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_ep_send[i] = n_sp_send[i] = n_image_per_proc[i] = 0;
        }
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        F64ort pos_root_domain = dinfo.getPosRootDomain();
        F64vec len_peri = pos_root_domain.getFullLength();
        F64ort outer_boundary_of_my_tree = GetOuterBoundaryOfMyTree(typename TSM::tree_cell_loc_geometry_type(), tc_first);
        // exchange outer boundary
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            minimum_tree_boundary_of_domain[i] = GetCorrespondingTreeCell(top_moment[i].vertex_in_);
        }
        F64 r_crit_send = 0.0;
        if(exchange_let_mode_==EXCHANGE_LET_P2P_EXACT){
            r_crit_send = minimum_tree_boundary_of_domain[my_rank].getFullLength().getMax();
        }else if(exchange_let_mode_==EXCHANGE_LET_P2P_FAST){
            r_crit_send = top_moment[my_rank].vertex_in_.getFullLength().getMax();
        }else{assert(0);}
        if(theta > 0.0){
            r_crit_send = r_crit_send*r_crit_send / (theta*theta);
        }else{
            r_crit_send = -1.0;
        }
PS_OMP_PARALLEL
        {
            const S32 ith = Comm::getThreadNum();
            const S32 n_thread = Comm::getNumberOfThread();
            const S32 head = (n_proc/n_thread)*ith + std::min(n_proc%n_thread, ith);
            const S32 end  = (n_proc/n_thread)*(ith+1) + std::min(n_proc%n_thread, (ith+1));
            rank_send_tmp[ith].clearSize();
            rank_recv_tmp[ith].clearSize();
            rank_recv_allgather_tmp[ith].clearSize();
            for(S32 ib=head; ib<end; ib++){
                const S32 id0 = std::min(ib, my_rank);
                const S32 id1 = std::max(ib, my_rank);
                bool flag_recv_allgather = true;
                if( ((typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER) || typeid(TSM) == typeid(SEARCH_MODE_LONG_SYMMETRY))
                         && (top_moment[id0].isOverlapped(top_moment[id1])))
                    || (theta <= 0.0) ){
                    rank_send_tmp[ith].push_back(ib);
                    rank_recv_tmp[ith].push_back(ib);
                    flag_recv_allgather = false;
                }
                else{
                    if( (GetDistanceMinSq<CALC_DISTANCE_TYPE>(top_moment[ib].vertex_in_, top_moment[my_rank].spj_.getPos(), len_peri) <= r_crit_send)
                        || ( (GetDistanceMinSq<CALC_DISTANCE_TYPE>(top_moment[ib].vertex_in_, minimum_tree_boundary_of_domain[my_rank], len_peri) <= 0.0)
                            && (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) ) ){
                        rank_send_tmp[ith].push_back(ib);
                    }
                    F64 r_crit_recv = 0.0;
                    if(exchange_let_mode_==EXCHANGE_LET_P2P_EXACT){
                        r_crit_recv = minimum_tree_boundary_of_domain[ib].getFullLength().getMax();
                    } else if(exchange_let_mode_==EXCHANGE_LET_P2P_FAST){
                        r_crit_recv = top_moment[ib].vertex_in_.getFullLength().getMax();
                    } else{assert(0);}
                    if(theta > 0.0){
                        r_crit_recv = r_crit_recv*r_crit_recv / (theta*theta);
                    }else{
                        r_crit_recv = -1.0;
                    }
                    if( (GetDistanceMinSq<CALC_DISTANCE_TYPE>(top_moment[my_rank].vertex_in_, top_moment[ib].spj_.getPos(), len_peri) <= r_crit_recv)
                       || ( (GetDistanceMinSq<CALC_DISTANCE_TYPE>(top_moment[my_rank].vertex_in_, minimum_tree_boundary_of_domain[ib], len_peri) <= 0.0)
                            && (exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) ) ){
                        rank_recv_tmp[ith].push_back(ib);
                        flag_recv_allgather = false;
                    }
                }
                if(flag_recv_allgather){
                    rank_recv_allgather_tmp[ith].push_back(ib);
                }
            }
        } // end of OMP scope
        //Comm::barrier();
        //exit(1);
        PackData(&rank_send_tmp[0], Comm::getNumberOfThread(), rank_send);
        PackData(&rank_recv_tmp[0], Comm::getNumberOfThread(), rank_recv);
        PackData(&rank_recv_allgather_tmp[0], Comm::getNumberOfThread(), rank_recv_a2a);
        
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_image_per_proc[i] = 0;
        }
        S32 adr_tc = 0;
        S32 adr_tree_sp_first = 0;
PS_OMP_PARALLEL
        {
            S32 ith = Comm::getThreadNum();
            rank_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ep_send_tmp[ith].clearSize();
            n_ep_per_image_tmp[ith].clearSize();
            adr_sp_send_tmp[ith].clearSize();
            n_sp_per_image_tmp[ith].clearSize();
            ReallocatableArray<Tsp> sp_first;
PS_OMP(omp for schedule(dynamic, 4))
            for(S32 ib=0; ib<rank_send.size(); ib++){
                const S32 rank = rank_send[ib];
                const S32 n_image_tmp_prev = shift_per_image_tmp[ith].size();
                rank_tmp[ith].push_back(rank);
                CalcNumberAndShiftOfImageDomain
                    (shift_per_image_tmp[ith],  dinfo.getPosRootDomain().getFullLength(),
                     outer_boundary_of_my_tree, dinfo.getPosDomain(rank), pa, false);
                const S32 n_image_tmp = shift_per_image_tmp[ith].size();
                n_image_per_proc[rank] = n_image_tmp - n_image_tmp_prev;
                S32 n_ep_prev = adr_ep_send_tmp[ith].size();
                S32 n_sp_prev = adr_sp_send_tmp[ith].size();
                for(S32 j=n_image_tmp_prev; j<n_image_tmp; j++){
                    S32 n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                    S32 n_sp_prev_2 = adr_sp_send_tmp[ith].size();
                    if(my_rank==rank && j==n_image_tmp_prev){
                        // self image
                        n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                        n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
                        continue;
                    }
                    //F64ort pos_target_domain = dinfo.getPosDomain(rank).shift(shift_per_image_tmp[ith][j]);
                    //TargetBox<TSM> target_box;
                    //SetTargetBoxExLet(target_box, pos_target_domain, outer_boundary_of_tree[rank].shift(shift_per_image_tmp[ith][j]));
                    TargetBox<TSM> target_box;
                    target_box.setForExLet(top_moment[rank], shift_per_image_tmp[ith][j]);
                    MakeListUsingTreeRecursiveTop
                        <TSM, Ttc, TreeParticle, Tep, Tsp, TagChopLeafFalse, TagCopyInfoCloseNoSp, CALC_DISTANCE_TYPE>
                        (tc_first, adr_tc, tp_first,
                         ep_first, adr_ep_send_tmp[ith],
                         sp_first, adr_sp_send_tmp[ith],
                         target_box,
                         n_leaf_limit,
                         adr_tree_sp_first, len_peri, theta);
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                    n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
                }
                n_ep_send[rank] = adr_ep_send_tmp[ith].size() - n_ep_prev;
                n_sp_send[rank] = adr_sp_send_tmp[ith].size() - n_sp_prev;
            }
        } // end of OMP scope
        ReallocatableArray<S32> n_disp_image_per_proc(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
        n_disp_image_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_per_proc[i+1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
        }
        const S32 n_image_tot = n_disp_image_per_proc[n_proc];
        shift_per_image.resizeNoInitialize( n_image_tot );
        n_ep_per_image.resizeNoInitialize( n_image_tot);
        n_sp_per_image.resizeNoInitialize( n_image_tot);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                S32 offset = n_disp_image_per_proc[rank];
                for(S32 k=0; k<n_image_per_proc[rank]; k++){
                    shift_per_image[offset+k] = shift_per_image_tmp[i][n_cnt];
                    n_ep_per_image[offset+k]  = n_ep_per_image_tmp[i][n_cnt];
                    n_sp_per_image[offset+k]  = n_sp_per_image_tmp[i][n_cnt];
                    n_cnt++;
                }
            }
        }
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1, n_image_tot+1, MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot+1, n_image_tot+1, MemoryAllocMode::Pool);
        n_disp_ep_per_image[0] = 0;
        n_disp_sp_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
            n_disp_sp_per_image[i+1] = n_disp_sp_per_image[i] + n_sp_per_image[i];
        }
        const S32 n_ep_send_tot = n_disp_ep_per_image[ n_image_tot ];
        const S32 n_sp_send_tot = n_disp_sp_per_image[ n_image_tot ];
        adr_ep_send.resizeNoInitialize( n_ep_send_tot );
        adr_sp_send.resizeNoInitialize( n_sp_send_tot );
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_ep = 0;
            S32 n_cnt_sp = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                const S32 adr_image_head = n_disp_image_per_proc[rank];
                const S32 adr_image_end = n_disp_image_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++){
                    const S32 adr_ep_head = n_disp_ep_per_image[k];
                    const S32 adr_ep_end = n_disp_ep_per_image[k+1];
                    for(S32 l=adr_ep_head; l<adr_ep_end; l++){
                        adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                    }
                    const S32 adr_sp_head = n_disp_sp_per_image[k];
                    const S32 adr_sp_end  = n_disp_sp_per_image[k+1];
                    for(S32 l=adr_sp_head; l<adr_sp_end; l++){
                        adr_sp_send[l] = adr_sp_send_tmp[i][n_cnt_sp++];
                    }
                }
            }
        }
#endif
}

template <typename TSM, typename Ttc, typename Ttp, typename Tep, typename Tsp, typename Tgeometry, typename Tmorton_key,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
inline void FindScatterParticleLong(const ReallocatableArray<Ttc> &tc_first, const ReallocatableArray<Ttp> &tp_first,
                                    const ReallocatableArray<Tep> &ep_first, ReallocatableArray<S32> &n_ep_send, ReallocatableArray<S32> &adr_ep_send,
                                    const DomainInfo &dinfo, const S32 n_leaf_limit, const S32 lev_leaf_limit, ReallocatableArray<S32> &n_sp_send,
                                    ReallocatableArray<S32> &adr_sp_send, ReallocatableArray<F64vec> &shift_per_image,
                                    ReallocatableArray<S32> &n_image_per_proc, ReallocatableArray<S32> &n_ep_per_image,
                                    ReallocatableArray<S32> &n_sp_per_image, const ReallocatableArray<Tgeometry> &top_geometry, const F64 theta,
                                    const S32 icut, const S32 lev_pm_cell, const Tmorton_key &morton_key) {
    PS_COMPILE_TIME_MESSAGE("FindScatterParticleLong");
#if defined(PS_DEBUG_PRINT_FindScatterParticleLong)
    PARTICLE_SIMULATOR_PRINT_MSG("****** START FindScatterParticleLong ******", 0);
#endif
    const auto n_thread = Comm::getNumberOfThread();
    const auto n_proc = dinfo.getCommInfo().getNumberOfProc();
    const auto my_rank = dinfo.getCommInfo().getRank();
    ReallocatableArray<S32> adr_ep_send_tmp[n_thread];
    ReallocatableArray<S32> adr_sp_send_tmp[n_thread];
    ReallocatableArray<S32> rank_tmp[n_thread];
    ReallocatableArray<F64vec> shift_per_image_tmp[n_thread];
    ReallocatableArray<S32> n_ep_per_image_tmp[n_thread];
    ReallocatableArray<S32> n_sp_per_image_tmp[n_thread];
    for (int i = 0; i < n_thread; i++) {
        rank_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        shift_per_image_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        n_ep_per_image_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        n_sp_per_image_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        if (Comm::getRank() == 14) {  // not OK
            adr_ep_send_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        }
        adr_sp_send_tmp[i].setAllocMode(MemoryAllocMode::Pool);
    }

    n_ep_send.resizeNoInitialize(n_proc);
    n_sp_send.resizeNoInitialize(n_proc);
    n_image_per_proc.resizeNoInitialize(n_proc);

    for (int i = 0; i < n_thread; i++) {
        // rank_tmp[i].reserve(n_proc*2);
        // shift_per_image_tmp[i].reserve(n_proc*2);
        // adr_ep_send_tmp[i].reserve(100000000);
        // n_ep_per_image_tmp[i].reserve(n_proc*2);
        // adr_sp_send_tmp[i].reserve(1000000);
        // n_sp_per_image_tmp[i].reserve(n_proc*2);
    }
    /*
#if defined(NO_PUSH_BACK)
    for(int i=0; i<n_thread; i++){
        rank_tmp[i].reserve(n_proc*2);
        shift_per_image_tmp[i].reserve(n_proc*2);
        adr_ep_send_tmp[i].reserve(10000000);
        n_ep_per_image_tmp[i].reserve(n_proc*2);
        adr_sp_send_tmp[i].reserve(1000000);
        n_sp_per_image_tmp[i].reserve(n_proc*2);
    }
#endif
    */

    bool pa[DIMENSION];
    dinfo.getPeriodicAxis(pa);
    F64ort pos_root_domain = dinfo.getPosRootDomain();
    F64vec len_peri = pos_root_domain.getFullLength();
    F64ort outer_boundary_of_my_tree = GetOuterBoundaryOfMyTree(typename TSM::tree_cell_loc_geometry_type(), tc_first);

#if PS_DEBUG_PRINT_FindScatterParticlePMM >= 2
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE(outer_boundary_of_my_tree);
    }
#endif

    S32 adr_tc = 0;
    S32 adr_tree_sp_first = 0;
    PS_OMP_PARALLEL {
        S32 ith = Comm::getThreadNum();
        rank_tmp[ith].clearSize();
        shift_per_image_tmp[ith].clearSize();
        adr_ep_send_tmp[ith].clearSize();
        n_ep_per_image_tmp[ith].clearSize();
        adr_sp_send_tmp[ith].clearSize();
        n_sp_per_image_tmp[ith].clearSize();
        ReallocatableArray<Tsp> sp_first;
        PS_OMP_FOR_D4
        for (S32 i = 0; i < n_proc; i++) {
            const S32 n_image_tmp_prev = shift_per_image_tmp[ith].size();
            rank_tmp[ith].push_back(i);
            // if (typeid(TSM) != typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE)) {
            if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
                CalcNumberAndShiftOfImageDomain(shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(), top_geometry[my_rank].vertex_in_,
                                                top_geometry[i].vertex_in_, pa, icut, morton_key, false);
            } else {
                CalcNumberAndShiftOfImageDomain(shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(), outer_boundary_of_my_tree,
                                                dinfo.getPosDomain(i), pa, false);
            }
            const S32 n_image_tmp = shift_per_image_tmp[ith].size();
            n_image_per_proc[i] = n_image_tmp - n_image_tmp_prev;
            /*
            if(Comm::getRank()==0){
                std::cerr<<"i= "<<" n_image_per_proc[i]= "<<n_image_per_proc[i]
                         <<" shift_per_image_tmp[0].size()= "<<shift_per_image_tmp[0].size()
                         <<std::endl;
                for(S32 j=0; j<shift_per_image_tmp[0].size(); j++){
                    std::cerr<<"j= "<<j<<" shift_per_image_tmp[0][j]= "<<shift_per_image_tmp[0][j]<<std::endl;
                }
            }
            */
            S32 n_ep_prev = adr_ep_send_tmp[ith].size();
            S32 n_sp_prev = adr_sp_send_tmp[ith].size();
            for (S32 j = n_image_tmp_prev; j < n_image_tmp; j++) {
                S32 n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                S32 n_sp_prev_2 = adr_sp_send_tmp[ith].size();
                if (my_rank == i && j == n_image_tmp_prev) {
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);  // is 0
                    n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);  // is 0
                    continue;
                }

                TargetBox<TSM> target_box;
                if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
                    target_box.setForExLet(top_geometry[i], shift_per_image_tmp[ith][j], icut, morton_key);
                } else {
                    target_box.setForExLet(top_geometry[i], shift_per_image_tmp[ith][j]);
                }

                // if(Comm::getRank()==0){
                //     std::cerr<<"i= "<<i<<" j= "<<j
                //              <<" shift_per_image_tmp[ith][j]= "<<shift_per_image_tmp[ith][j]
                //              <<" target_box.vertex_= "<<target_box.vertex_<<std::endl;
                // }
                /*
                if(N_STEP == 1431){
                    Comm::barrier();
                    std::cout<<"check B "<<Comm::getRank()<<std::endl;
                    Comm::barrier();
                }
                */
#if 1
                // if (typeid(TSM) != typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE)) {
                if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
                    MakeListUsingTreeRecursiveTopPMM<TSM, Ttc, TreeParticle, Tep, Tsp, Tmorton_key, TagChopLeafFalse, TagChopNonleafFalse,
                                                     TagCopyInfoCloseNoSp, CALC_DISTANCE_TYPE>(
                        tc_first, adr_tc, tp_first, ep_first, adr_ep_send_tmp[ith], sp_first, adr_sp_send_tmp[ith], target_box, n_leaf_limit,
                        lev_leaf_limit, adr_tree_sp_first, len_peri, theta, lev_pm_cell, morton_key);
                } else {
                    MakeListUsingTreeRecursiveTop<TSM, Ttc, TreeParticle, Tep, Tsp, TagChopLeafFalse, TagCopyInfoCloseNoSp, CALC_DISTANCE_TYPE>(
                        tc_first, adr_tc, tp_first, ep_first, adr_ep_send_tmp[ith], sp_first, adr_sp_send_tmp[ith], target_box, n_leaf_limit,
                        adr_tree_sp_first, len_peri, theta);
                }
#else
                MakeListUsingTreeRecursiveTopPMM<TSM, Ttc, TreeParticle, Tep, Tsp, Tmorton_key, TagChopLeafFalse, TagChopNonleafFalse,
                                                 TagCopyInfoCloseNoSp, CALC_DISTANCE_TYPE>(
                    tc_first, adr_tc, tp_first, ep_first, adr_ep_send_tmp[ith], sp_first, adr_sp_send_tmp[ith], target_box, n_leaf_limit,
                    lev_leaf_limit, adr_tree_sp_first, len_peri, theta, lev_pm_cell, morton_key);
#endif
                /*
                if(N_STEP == 1431){
                    std::cout<<"check C "<<Comm::getRank()<<std::endl;
                }
                */
                n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
            }
            n_ep_send[i] = adr_ep_send_tmp[ith].size() - n_ep_prev;
            n_sp_send[i] = adr_sp_send_tmp[ith].size() - n_sp_prev;
        }
    }  // end of OMP scope
#if PS_DEBUG_PRINT_FindScatterParticlePMM >= 2
    if (Comm::getRank() == 0) {
        for (int i = 0; i < n_proc; i++) {
            PS_PRINT_VARIABLE4(i, n_ep_send[i], n_sp_send[i], n_image_per_proc[i]);
        }
    }
#endif

    // ReallocatableArray<S32> n_disp_image_per_proc(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
    S32 n_disp_image_per_proc[n_proc + 1];
    // n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
    n_disp_image_per_proc[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_image_per_proc[i + 1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
    }
    const S32 n_image_tot = n_disp_image_per_proc[n_proc];
    shift_per_image.resizeNoInitialize(n_image_tot);
    n_ep_per_image.resizeNoInitialize(n_image_tot);
    n_sp_per_image.resizeNoInitialize(n_image_tot);
    // PS_OMP(omp parallel for num_threads(n_thread))
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_thread; i++) {
        S32 n_cnt = 0;
        for (S32 j = 0; j < rank_tmp[i].size(); j++) {
            S32 rank = rank_tmp[i][j];
            S32 offset = n_disp_image_per_proc[rank];
            for (S32 k = 0; k < n_image_per_proc[rank]; k++) {
                shift_per_image[offset + k] = shift_per_image_tmp[i][n_cnt];
                n_ep_per_image[offset + k] = n_ep_per_image_tmp[i][n_cnt];
                n_sp_per_image[offset + k] = n_sp_per_image_tmp[i][n_cnt];
                n_cnt++;
            }
        }
    }

#if PS_DEBUG_PRINT_FindScatterParticlePMM >= 2
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE(n_image_tot);
    }
#endif

    ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    n_disp_ep_per_image[0] = 0;
    n_disp_sp_per_image[0] = 0;
    for (S32 i = 0; i < n_image_tot; i++) {
        n_disp_ep_per_image[i + 1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
        n_disp_sp_per_image[i + 1] = n_disp_sp_per_image[i] + n_sp_per_image[i];
    }
    const S32 n_ep_send_tot = n_disp_ep_per_image[n_image_tot];
    const S32 n_sp_send_tot = n_disp_sp_per_image[n_image_tot];
    adr_ep_send.resizeNoInitialize(n_ep_send_tot);
    adr_sp_send.resizeNoInitialize(n_sp_send_tot);
    // PS_OMP(omp parallel for num_threads(n_thread))
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_thread; i++) {
        S32 n_cnt_ep = 0;
        S32 n_cnt_sp = 0;
        for (S32 j = 0; j < rank_tmp[i].size(); j++) {
            S32 rank = rank_tmp[i][j];
            const S32 adr_image_head = n_disp_image_per_proc[rank];
            const S32 adr_image_end = n_disp_image_per_proc[rank + 1];
            for (S32 k = adr_image_head; k < adr_image_end; k++) {
                const S32 adr_ep_head = n_disp_ep_per_image[k];
                const S32 adr_ep_end = n_disp_ep_per_image[k + 1];
                for (S32 l = adr_ep_head; l < adr_ep_end; l++) {
                    adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                }
                const S32 adr_sp_head = n_disp_sp_per_image[k];
                const S32 adr_sp_end = n_disp_sp_per_image[k + 1];
                for (S32 l = adr_sp_head; l < adr_sp_end; l++) {
                    adr_sp_send[l] = adr_sp_send_tmp[i][n_cnt_sp++];
                }
            }
        }
    }

    for (S32 i = 0; i < n_thread; i++) {
        rank_tmp[i].freeMem();
        shift_per_image_tmp[i].freeMem();
        n_ep_per_image_tmp[i].freeMem();
        n_sp_per_image_tmp[i].freeMem();
        adr_ep_send_tmp[i].freeMem();
        adr_sp_send_tmp[i].freeMem();
    }

#if defined(PS_DEBUG_PRINT_FindScatterParticlePMM)
#if PS_DEBUG_PRINT_FindScatterParticlePMM >= 2
    if (Comm::getRank() == 0) {
        //
    }
#endif
    PARTICLE_SIMULATOR_PRINT_MSG("****** END FindScatterParticlePMM ******", 0);
#endif
}

/////////////////
// for short mode
template <class Ttc, class Ttp, class Tep, class Tsp, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
inline void FindScatterParticleP2P(const ReallocatableArray<Ttc> &tc_first, const ReallocatableArray<Ttp> &tp_first,
                                   const ReallocatableArray<Tep> &ep_first, ReallocatableArray<S32> &n_ep_send, ReallocatableArray<S32> &adr_ep_send,
                                   const DomainInfo &dinfo, const S32 n_leaf_limit, ReallocatableArray<F64vec> &shift_per_image,
                                   ReallocatableArray<S32> &n_image_per_proc, ReallocatableArray<S32> &n_ep_per_image,
                                   ReallocatableArray<S32> &rank_send, ReallocatableArray<S32> &rank_recv) {
    PS_COMPILE_TIME_MESSAGE("FindScatterParticleP2P(short)");
    const auto n_thread = Comm::getNumberOfThread();
    const auto n_proc = dinfo.getCommInfo().getNumberOfProc();
    const auto my_rank = dinfo.getCommInfo().getRank();
    rank_send.resizeNoInitialize(n_proc);
    rank_recv.resizeNoInitialize(n_proc);
    n_ep_send.resizeNoInitialize(n_proc);
    n_image_per_proc.resizeNoInitialize(n_proc);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        rank_send[i] = rank_recv[i] = n_ep_send[i] = n_image_per_proc[i] = 0;
    }
    std::vector<ReallocatableArray<S32>> rank_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<F64vec>> shift_per_image_tmp(n_thread, ReallocatableArray<F64vec>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> adr_ep_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> n_ep_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> rank_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> rank_recv_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));

    ReallocatableArray<F64ort> outer_boundary_of_tree(n_proc, n_proc, MemoryAllocMode::Pool);
    bool pa[DIMENSION];
    dinfo.getPeriodicAxis(pa);
    const auto pos_root_domain = dinfo.getPosRootDomain();
    auto len_peri = pos_root_domain.getFullLength();
    // for(S32 i=0; i<DIMENSION; i++){
    //     if(pa[i]==false) len_peri[i] = 0.0;
    // }
    // const auto outer_boundary_of_my_tree = tc_first[0].mom_.vertex_out_;
    const auto outer_boundary_of_my_tree = tc_first[0].geo_.getVertexOut();
    // Comm::allGather(&outer_boundary_of_my_tree, 1, outer_boundary_of_tree.getPointer());
    dinfo.getCommInfo().allGather(&outer_boundary_of_my_tree, 1, outer_boundary_of_tree.getPointer());
    PS_OMP_PARALLEL {
        const auto ith = Comm::getThreadNum();
        const auto n_thread = Comm::getNumberOfThread();
        const auto head = (n_proc / n_thread) * ith + std::min(n_proc % n_thread, ith);
        const auto end = (n_proc / n_thread) * (ith + 1) + std::min(n_proc % n_thread, (ith + 1));
        rank_send_tmp[ith].clearSize();
        rank_recv_tmp[ith].clearSize();
        for (S32 ib = head; ib < end; ib++) {
            const S32 id0 = std::min(ib, my_rank);
            const S32 id1 = std::max(ib, my_rank);
            if (GetDistanceMinSq<CALC_DISTANCE_TYPE>(outer_boundary_of_tree[id0], outer_boundary_of_tree[id1], len_peri) <= 0.0) {
                rank_send_tmp[ith].push_back(ib);
                rank_recv_tmp[ith].push_back(ib);
            }
        }
    }  // end of OMP scope
    PackData(&rank_send_tmp[0], Comm::getNumberOfThread(), rank_send);
    PackData(&rank_recv_tmp[0], Comm::getNumberOfThread(), rank_recv);

    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        n_image_per_proc[i] = 0;
    }
    S32 adr_tc = 0;
    S32 adr_tree_sp_first = 0;
    PS_OMP_PARALLEL {
        auto ith = Comm::getThreadNum();
        rank_tmp[ith].clearSize();
        shift_per_image_tmp[ith].clearSize();
        adr_ep_send_tmp[ith].clearSize();
        n_ep_per_image_tmp[ith].clearSize();
        ReallocatableArray<S32> adr_sp_send_tmp;
        ReallocatableArray<Tsp> sp_first;
        // F64 r_crit_sq = 1.0; // dummy
        PS_OMP(omp for schedule(dynamic, 4))
        for (S32 ib = 0; ib < rank_send.size(); ib++) {
            const auto rank = rank_send[ib];
            const auto n_image_tmp_prev = shift_per_image_tmp[ith].size();
            rank_tmp[ith].push_back(rank);
            CalcNumberAndShiftOfImageDomain(shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(), outer_boundary_of_my_tree,
                                            dinfo.getPosDomain(rank), pa, false);
            const auto n_image_tmp = shift_per_image_tmp[ith].size();
            n_image_per_proc[rank] = n_image_tmp - n_image_tmp_prev;
            const auto n_ep_prev = adr_ep_send_tmp[ith].size();
            for (S32 j = n_image_tmp_prev; j < n_image_tmp; j++) {
                const auto n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                if (my_rank == rank && j == n_image_tmp_prev) {
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                    continue;
                }
                const auto pos_target_domain = dinfo.getPosDomain(rank).shift(shift_per_image_tmp[ith][j]);
                TargetBox<SEARCH_MODE_SCATTER> target_box;
                target_box.vertex_in_ = pos_target_domain;
                const F64 theta_tmp = 1.0;
                MakeListUsingTreeRecursiveTop<SEARCH_MODE_SCATTER, Ttc, TreeParticle, Tep, Tsp, TagChopLeafFalse, TagCopyInfoCloseNoSp,
                                              CALC_DISTANCE_TYPE>(tc_first, adr_tc, tp_first, ep_first, adr_ep_send_tmp[ith], sp_first,
                                                                  adr_sp_send_tmp, target_box, n_leaf_limit, adr_tree_sp_first, len_peri, theta_tmp);
                n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
            }
            n_ep_send[rank] = adr_ep_send_tmp[ith].size() - n_ep_prev;
        }
    }  // end of OMP scope
    ReallocatableArray<S32> n_disp_image_per_proc(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    n_disp_image_per_proc[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_image_per_proc[i + 1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
    }
    const auto n_image_tot = n_disp_image_per_proc[n_proc];
    shift_per_image.resizeNoInitialize(n_image_tot);
    n_ep_per_image.resizeNoInitialize(n_image_tot);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_thread; i++) {
        auto n_cnt = 0;
        for (S32 j = 0; j < rank_tmp[i].size(); j++) {
            auto rank = rank_tmp[i][j];
            auto offset = n_disp_image_per_proc[rank];
            for (S32 k = 0; k < n_image_per_proc[rank]; k++) {
                shift_per_image[offset + k] = shift_per_image_tmp[i][n_cnt];
                n_ep_per_image[offset + k] = n_ep_per_image_tmp[i][n_cnt];
                n_cnt++;
            }
        }
    }
    ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    n_disp_ep_per_image[0] = 0;
    for (S32 i = 0; i < n_image_tot; i++) {
        n_disp_ep_per_image[i + 1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
    }
    const auto n_ep_send_tot = n_disp_ep_per_image[n_image_tot];
    adr_ep_send.resizeNoInitialize(n_ep_send_tot);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_thread; i++) {
        auto n_cnt_ep = 0;
        for (S32 j = 0; j < rank_tmp[i].size(); j++) {
            auto rank = rank_tmp[i][j];
            const auto adr_image_head = n_disp_image_per_proc[rank];
            const auto adr_image_end = n_disp_image_per_proc[rank + 1];
            for (S32 k = adr_image_head; k < adr_image_end; k++) {
                const auto adr_ep_head = n_disp_ep_per_image[k];
                const auto adr_ep_end = n_disp_ep_per_image[k + 1];
                for (S32 l = adr_ep_head; l < adr_ep_end; l++) {
                    adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                }
            }
        }
    }
}

template <class Ttc, class Ttp, class Tep, class Tsp, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
inline void FindScatterParticleShort(const ReallocatableArray<Ttc> &tc_first, const ReallocatableArray<Ttp> &tp_first,
                                     const ReallocatableArray<Tep> &ep_first, ReallocatableArray<S32> &n_ep_send,
                                     ReallocatableArray<S32> &adr_ep_send, const DomainInfo &dinfo, const S32 n_leaf_limit,
                                     ReallocatableArray<F64vec> &shift_per_image, ReallocatableArray<S32> &n_image_per_proc,
                                     ReallocatableArray<S32> &n_ep_per_image) {
    PS_COMPILE_TIME_MESSAGE("FindScatterParticleShort");
    const S32 n_thread = Comm::getNumberOfThread();
    std::vector<ReallocatableArray<S32>> rank_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<F64vec>> shift_per_image_tmp(n_thread, ReallocatableArray<F64vec>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> adr_ep_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    std::vector<ReallocatableArray<S32>> n_ep_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));

    const auto n_proc = dinfo.getCommInfo().getNumberOfProc();
    const auto my_rank = dinfo.getCommInfo().getRank();
    n_ep_send.resizeNoInitialize(n_proc);
    n_image_per_proc.resizeNoInitialize(n_proc);
    bool pa[DIMENSION];
    dinfo.getPeriodicAxis(pa);
    auto pos_root_domain = dinfo.getPosRootDomain();
    auto len_peri = pos_root_domain.getFullLength();
    const auto outer_boundary_of_my_tree = tc_first[0].geo_.getVertexOut();
    auto adr_tc = 0;
    auto adr_tree_sp_first = 0;
    PS_OMP_PARALLEL {
        auto ith = Comm::getThreadNum();
        rank_tmp[ith].clearSize();
        shift_per_image_tmp[ith].clearSize();
        adr_ep_send_tmp[ith].clearSize();
        n_ep_per_image_tmp[ith].clearSize();
        ReallocatableArray<S32> adr_sp_send_tmp;
        ReallocatableArray<Tsp> sp_first;
        // F64 r_crit_sq = 1.0; // dummy
        F64 theta = 1.0;  // dummy
        PS_OMP(omp for schedule(dynamic, 4))
        for (S32 i = 0; i < n_proc; i++) {
            const auto n_image_tmp_prev = shift_per_image_tmp[ith].size();
            rank_tmp[ith].push_back(i);
            CalcNumberAndShiftOfImageDomain(shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(), outer_boundary_of_my_tree,
                                            dinfo.getPosDomain(i), pa, false);

            // if(my_rank == 0){
            //     std::cerr<<"outer_boundary_of_my_tree= "<<outer_boundary_of_my_tree
            //<<" dinfo.getPosDomain(i)= "<<dinfo.getPosDomain(i)<<std::endl;
            //}

            const auto n_image_tmp = shift_per_image_tmp[ith].size();
            n_image_per_proc[i] = n_image_tmp - n_image_tmp_prev;
            auto n_ep_prev = adr_ep_send_tmp[ith].size();
            for (S32 j = n_image_tmp_prev; j < n_image_tmp; j++) {
                auto n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                if (my_rank == i && j == n_image_tmp_prev) {
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);  // is 0
                    continue;
                }
                auto pos_target_domain = dinfo.getPosDomain(i).shift(shift_per_image_tmp[ith][j]);
                TargetBox<SEARCH_MODE_SCATTER> target_box;
                target_box.vertex_in_ = pos_target_domain;
                MakeListUsingTreeRecursiveTop<SEARCH_MODE_SCATTER, Ttc, TreeParticle, Tep, Tsp, TagChopLeafFalse, TagCopyInfoCloseNoSp,
                                              CALC_DISTANCE_TYPE>(tc_first, adr_tc, tp_first, ep_first, adr_ep_send_tmp[ith], sp_first,
                                                                  adr_sp_send_tmp, target_box, n_leaf_limit, adr_tree_sp_first, len_peri, theta);
                n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
            }
            n_ep_send[i] = adr_ep_send_tmp[ith].size() - n_ep_prev;
        }
    }  // end of OMP scope
    ReallocatableArray<S32> n_disp_image_per_proc(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    n_disp_image_per_proc[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_image_per_proc[i + 1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
    }
    const auto n_image_tot = n_disp_image_per_proc[n_proc];
    shift_per_image.resizeNoInitialize(n_image_tot);
    n_ep_per_image.resizeNoInitialize(n_image_tot);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_thread; i++) {
        S32 n_cnt = 0;
        for (S32 j = 0; j < rank_tmp[i].size(); j++) {
            auto rank = rank_tmp[i][j];
            auto offset = n_disp_image_per_proc[rank];
            for (S32 k = 0; k < n_image_per_proc[rank]; k++) {
                shift_per_image[offset + k] = shift_per_image_tmp[i][n_cnt];
                n_ep_per_image[offset + k] = n_ep_per_image_tmp[i][n_cnt];
                n_cnt++;
            }
        }
    }
    ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    n_disp_ep_per_image[0] = 0;
    for (S32 i = 0; i < n_image_tot; i++) {
        n_disp_ep_per_image[i + 1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
    }
    const auto n_ep_send_tot = n_disp_ep_per_image[n_image_tot];
    adr_ep_send.resizeNoInitialize(n_ep_send_tot);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_thread; i++) {
        auto n_cnt_ep = 0;
        for (S32 j = 0; j < rank_tmp[i].size(); j++) {
            auto rank = rank_tmp[i][j];
            const auto adr_image_head = n_disp_image_per_proc[rank];
            const auto adr_image_end = n_disp_image_per_proc[rank + 1];
            for (S32 k = adr_image_head; k < adr_image_end; k++) {
                const auto adr_ep_head = n_disp_ep_per_image[k];
                const auto adr_ep_end = n_disp_ep_per_image[k + 1];
                for (S32 l = adr_ep_head; l < adr_ep_end; l++) {
                    adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                }
            }
        }
    }
}

///////////////////
// exchange # of LET
////// FOR LONG SEARCH
inline void ExchangeNumberLong(const ReallocatableArray<S32> &n_ep_send, ReallocatableArray<S32> &n_ep_recv, const ReallocatableArray<S32> &n_sp_send,
                               ReallocatableArray<S32> &n_sp_recv, const ReallocatableArray<S32> &rank_send, const ReallocatableArray<S32> &rank_recv,
                               const CommInfo &comm_info) {
    ReallocatableArray<S32> n_ep_sp_send(rank_send.size() * 2, rank_send.size() * 2, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_ep_sp_recv(rank_recv.size() * 2, rank_recv.size() * 2, MemoryAllocMode::Pool);
    PS_OMP_PARALLEL_FOR
    for (int i = 0; i < rank_send.size(); i++) {
        const S32 rank = rank_send[i];
        n_ep_sp_send[i * 2] = n_ep_send[rank];
        n_ep_sp_send[i * 2 + 1] = n_sp_send[rank];
    }
    // Comm::sendIrecv(n_ep_sp_send.getPointer(), rank_send.getPointer(), 2, rank_send.size(), n_ep_sp_recv.getPointer(), rank_recv.getPointer(), 2,
    // rank_recv.size());
    comm_info.sendIrecv(n_ep_sp_send.getPointer(), rank_send.getPointer(), 2, rank_send.size(), n_ep_sp_recv.getPointer(), rank_recv.getPointer(), 2,
                        rank_recv.size());
    const S32 n_proc = comm_info.getNumberOfProc();
    n_ep_recv.resizeNoInitialize(n_proc);
    n_sp_recv.resizeNoInitialize(n_proc);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        n_ep_recv[i] = n_sp_recv[i] = 0;
    }
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < rank_recv.size(); i++) {
        const S32 rank = rank_recv[i];
        n_ep_recv[rank] = n_ep_sp_recv[i * 2];
        n_sp_recv[rank] = n_ep_sp_recv[i * 2 + 1];
    }
}
inline void ExchangeNumberLong(const ReallocatableArray<S32> &n_ep_send, ReallocatableArray<S32> &n_ep_recv, const ReallocatableArray<S32> &n_sp_send,
                               ReallocatableArray<S32> &n_sp_recv, const CommInfo &comm_info) {
    const S32 n_proc = comm_info.getNumberOfProc();
    ReallocatableArray<S32> n_ep_sp_send(n_proc * 2, n_proc * 2, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_ep_sp_recv(n_proc * 2, n_proc * 2, MemoryAllocMode::Pool);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        n_ep_sp_send[i * 2] = n_ep_send[i];
        n_ep_sp_send[i * 2 + 1] = n_sp_send[i];
    }
#ifdef FAST_ALL_TO_ALL_FOR_K
    CommForAllToAll<S32, 2> comm_a2a_2d;
    comm_a2a_2d.execute(n_ep_sp_send, 2, n_ep_sp_recv);
#else
    // Comm::allToAll(n_ep_sp_send.getPointer(), 2, n_ep_sp_recv.getPointer());
    comm_info.allToAll(n_ep_sp_send.getPointer(), 2, n_ep_sp_recv.getPointer());
#endif
    n_ep_recv.resizeNoInitialize(n_proc);
    n_sp_recv.resizeNoInitialize(n_proc);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        n_ep_recv[i] = n_ep_sp_recv[i * 2];
        n_sp_recv[i] = n_ep_sp_recv[i * 2 + 1];
    }
}

///////////////////////
////// FOR SHORT SEARCH
inline void ExchangeNumberShort(const ReallocatableArray<S32> &n_ep_send, ReallocatableArray<S32> &n_ep_recv,
                                const ReallocatableArray<S32> &rank_send, const ReallocatableArray<S32> &rank_recv, const CommInfo &comm_info) {
    ReallocatableArray<S32> n_ep_send_tmp(rank_send.size(), rank_send.size(), MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_ep_recv_tmp(rank_recv.size(), rank_recv.size(), MemoryAllocMode::Pool);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < rank_send.size(); i++) {
        const auto rank = rank_send[i];
        n_ep_send_tmp[i] = n_ep_send[rank];
    }
    // Comm::sendIrecv(n_ep_send_tmp.getPointer(), rank_send.getPointer(), 1, rank_send.size(), n_ep_recv_tmp.getPointer(), rank_recv.getPointer(), 1,
    // rank_recv.size());
    comm_info.sendIrecv(n_ep_send_tmp.getPointer(), rank_send.getPointer(), 1, rank_send.size(), n_ep_recv_tmp.getPointer(), rank_recv.getPointer(),
                        1, rank_recv.size());
    const auto n_proc = comm_info.getNumberOfProc();
    n_ep_recv.resizeNoInitialize(n_proc);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        n_ep_recv[i] = 0;
    }
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < rank_recv.size(); i++) {
        const auto rank = rank_recv[i];
        n_ep_recv[rank] = n_ep_recv_tmp[i];
    }
}

inline void ExchangeNumberShort(const ReallocatableArray<S32> &n_ep_send, ReallocatableArray<S32> &n_ep_recv, const CommInfo &comm_info) {
    const auto n_proc = comm_info.getNumberOfProc();
    n_ep_recv.resizeNoInitialize(n_proc);
#ifdef FAST_ALL_TO_ALL_FOR_K
    CommForAllToAll<S32, 2> comm_a2a_2d;
    comm_a2a_2d.execute(n_ep_send, 1, n_ep_recv);
#else
    // Comm::allToAll(n_ep_send.getPointer(), 1, n_ep_recv.getPointer());
    comm_info.allToAll(n_ep_send.getPointer(), 1, n_ep_recv.getPointer());
#endif
}
// exchange # of LET
//////////////

//////////////////
// EXCHANGE LET //
//////////////////
// FOR LONG SEARCH
template <class TSM, class Tep, class Tsp, class Ttc>
inline void ExchangeLetP2P(const ReallocatableArray<Tep> &ep, const ReallocatableArray<S32> &n_ep_send, const ReallocatableArray<S32> &n_ep_recv,
                           const ReallocatableArray<S32> &n_ep_per_image, const ReallocatableArray<S32> &adr_ep_send, ReallocatableArray<Tep> &ep_org,
                           const S32 n_ep_offset,  // ep_org[n_ep_offset] = n_ep_recv[0]
                           const ReallocatableArray<Ttc> &tc, const ReallocatableArray<S32> &n_sp_send, const ReallocatableArray<S32> &n_sp_recv,
                           const ReallocatableArray<S32> &n_sp_per_image, const ReallocatableArray<S32> &adr_sp_send, ReallocatableArray<Tsp> &sp_org,
                           const ReallocatableArray<F64vec> &shift_image_domain, const ReallocatableArray<S32> &n_image_per_proc,
                           const ReallocatableArray<S32> &rank_send, const ReallocatableArray<S32> &rank_recv, const CommInfo &comm_info) {
    ReallocatableArray<S32> n_disp_ep_send(rank_send.size() + 1, rank_send.size() + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_sp_send(rank_send.size() + 1, rank_send.size() + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_ep_recv(rank_recv.size() + 1, rank_recv.size() + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_sp_recv(rank_recv.size() + 1, rank_recv.size() + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_ep_send_tmp(rank_send.size(), rank_send.size(), MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_sp_send_tmp(rank_send.size(), rank_send.size(), MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_ep_recv_tmp(rank_recv.size(), rank_recv.size(), MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_sp_recv_tmp(rank_recv.size(), rank_recv.size(), MemoryAllocMode::Pool);

    n_disp_ep_send[0] = n_disp_sp_send[0] = n_disp_ep_recv[0] = n_disp_sp_recv[0] = 0;
    for (int i = 0; i < rank_send.size(); i++) {
        const S32 rank = rank_send[i];
        n_disp_ep_send[i + 1] = n_ep_send[rank] + n_disp_ep_send[i];
        n_disp_sp_send[i + 1] = n_sp_send[rank] + n_disp_sp_send[i];
        n_ep_send_tmp[i] = n_ep_send[rank];
        n_sp_send_tmp[i] = n_sp_send[rank];
    }
    for (int i = 0; i < rank_recv.size(); i++) {
        const S32 rank = rank_recv[i];
        n_disp_ep_recv[i + 1] = n_ep_recv[rank] + n_disp_ep_recv[i];
        n_disp_sp_recv[i + 1] = n_sp_recv[rank] + n_disp_sp_recv[i];
        n_ep_recv_tmp[i] = n_ep_recv[rank];
        n_sp_recv_tmp[i] = n_sp_recv[rank];
    }
    const S32 n_image_tot = n_ep_per_image.size();
    ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    n_disp_ep_per_image[0] = n_disp_sp_per_image[0] = 0;
    for (S32 i = 0; i < n_image_tot; i++) {
        n_disp_ep_per_image[i + 1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
        n_disp_sp_per_image[i + 1] = n_disp_sp_per_image[i] + n_sp_per_image[i];
    }
    ReallocatableArray<Tep> ep_send(n_disp_ep_send[rank_send.size()], n_disp_ep_send[rank_send.size()], MemoryAllocMode::Pool);
    ReallocatableArray<Tsp> sp_send(n_disp_sp_send[rank_send.size()], n_disp_sp_send[rank_send.size()], MemoryAllocMode::Pool);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif
    for (S32 i = 0; i < n_image_tot; i++) {
        F64vec shift = shift_image_domain[i];
        S32 n_ep = n_ep_per_image[i];
        S32 n_ep_offset = n_disp_ep_per_image[i];
        CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep, n_ep_offset, shift);
        S32 n_sp = n_sp_per_image[i];
        S32 n_sp_offset = n_disp_sp_per_image[i];
        CopyPtclFromTreeToSendBuf(typename TSM::search_boundary_type(), sp_send, tc, adr_sp_send, n_sp, n_sp_offset, shift);
    }
    ep_org.resizeNoInitialize(n_disp_ep_recv[rank_recv.size()] + n_ep_offset);
    sp_org.resizeNoInitialize(n_disp_sp_recv[rank_recv.size()]);
#ifdef PS_MEASURE_BARRIER
    comm_info_.barrier();
#endif
    // const auto wtime_offset = GetWtime();
    comm_info.sendIrecvV(ep_send.getPointer(), rank_send.getPointer(), n_ep_send_tmp.getPointer(), n_disp_ep_send.getPointer(), rank_send.size(),
                         sp_send.getPointer(), rank_send.getPointer(), n_sp_send_tmp.getPointer(), n_disp_sp_send.getPointer(), rank_send.size(),
                         ep_org.getPointer(n_ep_offset), rank_recv.getPointer(), n_ep_recv_tmp.getPointer(), n_disp_ep_recv.getPointer(),
                         rank_recv.size(), sp_org.getPointer(), rank_recv.getPointer(), n_sp_recv_tmp.getPointer(), n_disp_sp_recv.getPointer(),
                         rank_recv.size());
#ifdef PS_MEASURE_BARRIER
    comm_info.barrier();
#endif
    // wtime_comm = GetWtime() - wtime_offset;
}

template <class TSM, class Tep, class Tsp, class Ttc>
inline void ExchangeLet(const ReallocatableArray<Tep> &ep, const ReallocatableArray<S32> &n_ep_send, const ReallocatableArray<S32> &n_ep_recv,
                        const ReallocatableArray<S32> &n_ep_per_image, const ReallocatableArray<S32> &adr_ep_send, ReallocatableArray<Tep> &ep_org,
                        const S32 n_ep_offset, const ReallocatableArray<Ttc> &tc, const ReallocatableArray<S32> &n_sp_send,
                        const ReallocatableArray<S32> &n_sp_recv, const ReallocatableArray<S32> &n_sp_per_image,
                        const ReallocatableArray<S32> &adr_sp_send, ReallocatableArray<Tsp> &sp_org,
                        const ReallocatableArray<F64vec> &shift_image_domain, const ReallocatableArray<S32> &n_image_per_proc,
                        const CommInfo &comm_info) {
    // std::cerr<<"ExchangeLet"<<std::endl;
    const S32 n_proc = comm_info.getNumberOfProc();
    ReallocatableArray<S32> n_disp_ep_send(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_sp_send(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_ep_recv(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_sp_recv(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    n_disp_ep_send[0] = n_disp_sp_send[0] = n_disp_ep_recv[0] = n_disp_sp_recv[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_ep_send[i + 1] = n_ep_send[i] + n_disp_ep_send[i];
        n_disp_sp_send[i + 1] = n_sp_send[i] + n_disp_sp_send[i];
        n_disp_ep_recv[i + 1] = n_ep_recv[i] + n_disp_ep_recv[i];
        n_disp_sp_recv[i + 1] = n_sp_recv[i] + n_disp_sp_recv[i];
    }
    const S32 n_image_tot = n_ep_per_image.size();
    // std::cerr<<"n_image_tot= "<<n_image_tot<<std::endl;
    ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    n_disp_ep_per_image[0] = n_disp_sp_per_image[0] = 0;
    for (S32 i = 0; i < n_image_tot; i++) {
        n_disp_ep_per_image[i + 1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
        n_disp_sp_per_image[i + 1] = n_disp_sp_per_image[i] + n_sp_per_image[i];
    }
    ReallocatableArray<Tep> ep_send(n_disp_ep_send[n_proc], n_disp_ep_send[n_proc], MemoryAllocMode::Pool);
    ReallocatableArray<Tsp> sp_send(n_disp_sp_send[n_proc], n_disp_sp_send[n_proc], MemoryAllocMode::Pool);

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif
    for (S32 i = 0; i < n_image_tot; i++) {
#if 0
            F64 cm_mass = 0.0;
            F64vec cm_pos = 0.0;
#endif
        F64vec shift = shift_image_domain[i];

        S32 n_ep = n_ep_per_image[i];
        S32 n_ep_offset = n_disp_ep_per_image[i];
        CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep, n_ep_offset, shift);
#if 0
            for(S32 i2=n_ep_offset; i2<n_ep_offset+n_ep; i2++){
                cm_mass += ep_send[i2].mass;
                cm_pos += ep_send[i2].mass * ep_send[i2].pos;
            }
#endif
        S32 n_sp = n_sp_per_image[i];
        S32 n_sp_offset = n_disp_sp_per_image[i];
        CopyPtclFromTreeToSendBuf(typename TSM::search_boundary_type(), sp_send, tc, adr_sp_send, n_sp, n_sp_offset, shift);
#if 0
            for(S32 i2=n_sp_offset; i2<n_sp_offset+n_sp; i2++){
                cm_mass += sp_send[i2].mass;
                cm_pos += sp_send[i2].mass * sp_send[i2].pos;
            }
            if(Comm::getRank()==0){
                std::cerr<<"i= "<<i
                         <<" cm_mass= "<<cm_mass
                         <<" cm_pos= "<<cm_pos / cm_mass
                         <<std::endl;
            }
#endif
    }
    // exit(1);
    ep_org.resizeNoInitialize(n_disp_ep_recv[n_proc] + n_ep_offset);
    sp_org.resizeNoInitialize(n_disp_sp_recv[n_proc]);

#if 0
        if(Comm::getRank()==0){
            std::cout<<"DEBUG EXCHANGE_LET in tree_for_force_utils2.hpp"<<std::endl;
            if(Comm::getNumberOfProc() > 1){
                std::cout<<"n_sp_send[1]= "<<n_sp_send[1]<<std::endl;
                for(int i=0; i<n_sp_send[1]; i++){
                    std::cout<<"i= "<<i<<" sp_send[i].pos= "<<sp_send[i].pos<<" sp_send[i].mass= "<<sp_send[i].mass<<std::endl;
                }
            }
        }
#endif

#ifdef FAST_ALL_TO_ALL_FOR_K
    CommForAllToAll<Tep, 2> comm_a2a_epj_2d;
    comm_a2a_epj_2d.executeV(ep_send, ep_org, n_ep_send.getPointer(), n_ep_recv.getPointer(), 0, n_ep_offset);
    CommForAllToAll<Tsp, 2> comm_a2a_spj_2d;
    comm_a2a_spj_2d.executeV(sp_send, sp_org, n_sp_send.getPointer(), n_sp_recv.getPointer(), 0, 0);
#else
    comm_info.allToAllV(ep_send.getPointer(), n_ep_send.getPointer(), n_disp_ep_send.getPointer(), ep_org.getPointer(n_ep_offset),
                        n_ep_recv.getPointer(), n_disp_ep_recv.getPointer());

    /*
    if(Comm::getRank()==0){
        std::cerr<<"sp_send.size()= "<<sp_send.size()<<std::endl;
        for(auto i=0; i<sp_send.size(); i++){
            std::cerr<<"i= "<<i<<std::endl;
            sp_send[i].dump(std::cerr);
        }
    }
    */
    comm_info.allToAllV(sp_send.getPointer(), n_sp_send.getPointer(), n_disp_sp_send.getPointer(), sp_org.getPointer(), n_sp_recv.getPointer(),
                        n_disp_sp_recv.getPointer());

#endif

#if 0
        F64 cm_mass_glb = 0.0;
        F64vec cm_pos_glb = 0.0;        
        for(int i=0; i<ep_org.size(); i++){
            cm_mass_glb += ep_org[i].getCharge();
            cm_pos_glb  += ep_org[i].getCharge()*ep_org[i].getPos();
        }
        cm_pos_glb /= cm_mass_glb;
        std::cerr<<"cm_mass_glb= "<<cm_mass_glb
                 <<" cm_pos_glb= "<<cm_pos_glb
                 <<std::endl;
#endif
}

//////////////////
// FOR SHORT SEARCH
template <typename TSM, typename Tep>
inline void ExchangeLetShortP2P(const ReallocatableArray<Tep> &ep, const ReallocatableArray<S32> &n_ep_send, const ReallocatableArray<S32> &n_ep_recv,
                                const ReallocatableArray<S32> &n_ep_per_image, const ReallocatableArray<S32> &adr_ep_send,
                                ReallocatableArray<Tep> &ep_org, const S32 n_ep_offset, const ReallocatableArray<F64vec> &shift_image_domain,
                                const ReallocatableArray<S32> &n_image_per_proc, const ReallocatableArray<S32> &rank_send,
                                const ReallocatableArray<S32> &rank_recv, const CommInfo &comm_info) {
    PS_COMPILE_TIME_MESSAGE("ExchangeLetShortP2P");
    ReallocatableArray<S32> n_ep_send_tmp(rank_send.size(), rank_send.size(), MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_ep_recv_tmp(rank_recv.size(), rank_recv.size(), MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_ep_send(rank_send.size() + 1, rank_send.size() + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_ep_recv(rank_recv.size() + 1, rank_recv.size() + 1, MemoryAllocMode::Pool);
    n_disp_ep_send[0] = n_disp_ep_recv[0] = 0;
    for (int i = 0; i < rank_send.size(); i++) {
        const S32 rank = rank_send[i];
        n_disp_ep_send[i + 1] = n_ep_send[rank] + n_disp_ep_send[i];
        n_ep_send_tmp[i] = n_ep_send[rank];
    }
    for (int i = 0; i < rank_recv.size(); i++) {
        const S32 rank = rank_recv[i];
        n_disp_ep_recv[i + 1] = n_ep_recv[rank] + n_disp_ep_recv[i];
        n_ep_recv_tmp[i] = n_ep_recv[rank];
    }
    const S32 n_image_tot = n_ep_per_image.size();
    ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    n_disp_ep_per_image[0] = 0;
    for (S32 i = 0; i < n_image_tot; i++) {
        n_disp_ep_per_image[i + 1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
    }
    ReallocatableArray<Tep> ep_send(n_disp_ep_send[rank_send.size()], n_disp_ep_send[rank_send.size()], MemoryAllocMode::Pool);
    PS_OMP(omp parallel for schedule(dynamic, 4))
    for (S32 i = 0; i < n_image_tot; i++) {
        F64vec shift = shift_image_domain[i];
        S32 n_ep = n_ep_per_image[i];
        S32 n_ep_offset = n_disp_ep_per_image[i];
        CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep, n_ep_offset, shift);
    }
    ep_org.resizeNoInitialize(n_disp_ep_recv[rank_recv.size()] + n_ep_offset);
    comm_info.sendIrecvV(ep_send.getPointer(), rank_send.getPointer(), n_ep_send_tmp.getPointer(), n_disp_ep_send.getPointer(), rank_send.size(),
                         ep_org.getPointer(n_ep_offset), rank_recv.getPointer(), n_ep_recv_tmp.getPointer(), n_disp_ep_recv.getPointer(),
                         rank_recv.size());
}

template <typename TSM, typename Tep>
inline void ExchangeLetShortA2A(const ReallocatableArray<Tep> &ep, const ReallocatableArray<S32> &n_ep_send, const ReallocatableArray<S32> &n_ep_recv,
                                const ReallocatableArray<S32> &n_ep_per_image, const ReallocatableArray<S32> &adr_ep_send,
                                ReallocatableArray<Tep> &ep_org, const S32 n_ep_offset, const ReallocatableArray<F64vec> &shift_image_domain,
                                const ReallocatableArray<S32> &n_image_per_proc, const CommInfo &comm_info) {
#if defined(PS_DEBUG_PRINT_ExchangeLetShortA2A)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START ExchangeLetShortA2A ***", 0);    
#endif
    PS_COMPILE_TIME_MESSAGE("ExchangeLetShortA2A");
    const S32 n_proc = comm_info.getNumberOfProc();
    ReallocatableArray<S32> n_disp_ep_send(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_ep_recv(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    n_disp_ep_send[0] = n_disp_ep_recv[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_ep_send[i + 1] = n_ep_send[i] + n_disp_ep_send[i];
        n_disp_ep_recv[i + 1] = n_ep_recv[i] + n_disp_ep_recv[i];
    }
    const S32 n_image_tot = n_ep_per_image.size();
    ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot + 1, n_image_tot + 1, MemoryAllocMode::Pool);
    n_disp_ep_per_image[0] = 0;
    for (S32 i = 0; i < n_image_tot; i++) {
        n_disp_ep_per_image[i + 1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
    }
    ReallocatableArray<Tep> ep_send(n_disp_ep_send[n_proc], n_disp_ep_send[n_proc], MemoryAllocMode::Pool);

#if defined(PS_DEBUG_PRINT_ExchangeLetShortA2A)
    PARTICLE_SIMULATOR_PRINT_MSG("*** before CopyPtclToSendBuf ***", 0);
#endif
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
//#pragma omp parallel for schedule(dynamic, 4)
#endif
    for (S32 i = 0; i < n_image_tot; i++) {
        F64vec shift = shift_image_domain[i];
        S32 n_ep = n_ep_per_image[i];
        S32 n_ep_offset = n_disp_ep_per_image[i];
#if defined(PS_DEBUG_PRINT_ExchangeLetShortA2A)
        std::cerr<<"i= "<<i<<" n_ep= "<<n_ep<<" n_ep_offset= "<<n_ep_offset<<std::endl;
#endif
        CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep, n_ep_offset, shift);
    }
    ep_org.resizeNoInitialize(n_disp_ep_recv[n_proc] + n_ep_offset);
#if defined(PS_DEBUG_PRINT_ExchangeLetShortA2A)
    PARTICLE_SIMULATOR_PRINT_MSG("*** before exchange let with a2a ***", 0);
#endif    
#if defined(FAST_ALL_TO_ALL_FOR_K)
    CommForAllToAll<Tep, 2> comm_a2a_epj_2d;
    comm_a2a_epj_2d.executeV(ep_send, ep_org, n_ep_send.getPointer(), n_ep_recv.getPointer(), 0, n_ep_offset);
#else
    comm_info.allToAllV(ep_send.getPointer(), n_ep_send.getPointer(), n_disp_ep_send.getPointer(), ep_org.getPointer(n_ep_offset),
                        n_ep_recv.getPointer(), n_disp_ep_recv.getPointer());
#endif
#if defined(PS_DEBUG_PRINT_ExchangeLetShortA2A)
    PARTICLE_SIMULATOR_PRINT_MSG("*** END ExchangeLetShortA2A ***", 0);    
#endif
}

template <typename TSM, typename Tep>
inline void ExchangeLetShort(const ReallocatableArray<Tep> &ep, const ReallocatableArray<S32> &n_ep_send, const ReallocatableArray<S32> &n_ep_recv,
                             const ReallocatableArray<S32> &n_ep_per_image, const ReallocatableArray<S32> &adr_ep_send,
                             ReallocatableArray<Tep> &ep_org, const S32 n_ep_offset, const ReallocatableArray<F64vec> &shift_image_domain,
                             const ReallocatableArray<S32> &n_image_per_proc, const CommInfo &comm_info) {
    ExchangeLetShortA2A(ep, n_ep_send, n_ep_recv, n_ep_per_image, adr_ep_send, ep_org, n_ep_offset, shift_image_domain, n_image_per_proc, comm_info);
}
template <typename TSM, typename Tep>
inline void ExchangeLet(const ReallocatableArray<Tep> &ep, const ReallocatableArray<S32> &n_ep_send, const ReallocatableArray<S32> &n_ep_recv,
                        const ReallocatableArray<S32> &n_ep_per_image, const ReallocatableArray<S32> &adr_ep_send,
                        ReallocatableArray<Tep> &ep_org, const S32 n_ep_offset, const ReallocatableArray<F64vec> &shift_image_domain,
                        const ReallocatableArray<S32> &n_image_per_proc, const CommInfo &comm_info) {
    ExchangeLetShort(ep, n_ep_send, n_ep_recv, n_ep_per_image, adr_ep_send, ep_org, n_ep_offset, shift_image_domain, n_image_per_proc, comm_info);
}
// exchange LET
//////////////

//////////////////
// add moment as sp
template <class Ttreecell, class Tspj>
inline void AddMomentAsSp(TagForceLong, const ReallocatableArray<Ttreecell> &_tc, const S32 offset, ReallocatableArray<Tspj> &_spj) {
    _spj.resizeNoInitialize(offset + _tc.size());
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < _tc.size(); i++) {
        _spj[offset + i].copyFromMoment(_tc[i].mom_);
    }
}
template <class Ttreecell, class Tspj>
inline void AddMomentAsSp(TagForceShort, const ReallocatableArray<Ttreecell> &_tc, const S32 offset, ReallocatableArray<Tspj> &_spj) {
    // do nothing
}

// for long force
template <class Ttp, class Tepj, class Tspj, typename Tmorton_key>
inline void SetLocalEssentialTreeToGlobalTreeLong(const ReallocatableArray<Tepj> &epj_org, const ReallocatableArray<Tspj> &spj_org, const S32 n_loc,
                                                  ReallocatableArray<Ttp> &tp_glb, const Tmorton_key &morton_key, const bool flag_reuse = false) {
    PS_COMPILE_TIME_MESSAGE("SetLocalEssentialTreeToGlobalTreeLong");
#if defined(PS_DEBUG_PRINT_SetLocalEssentialTreeToGlobalTreeLong)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START SetLocalEssentialTreeToGlobalTreeLong *****", 0);
#endif
    const S32 n_loc_ep = epj_org.size();
    const S32 n_loc_ep_sp = n_loc_ep + spj_org.size();
#if defined(PS_DEBUG_PRINT_SetLocalEssentialTreeToGlobalTreeLong)
PS_PRINT_VARIABLE3(n_loc, n_loc_ep, n_loc_ep_sp);
#endif
    tp_glb.resizeNoInitialize(n_loc_ep_sp);
    if (!flag_reuse) {
        PS_OMP_PARALLEL {
            PS_OMP_FOR
            for (S32 i = n_loc; i < n_loc_ep; i++) {
                tp_glb[i].setFromEP(epj_org[i], i, morton_key);
            }
            PS_OMP_FOR
            for (S32 i = n_loc_ep; i < n_loc_ep_sp; i++) {
                const S32 i_src = i - n_loc_ep;
                tp_glb[i].setFromSP(spj_org[i_src], i_src, morton_key);
            }
        }
    }
#if defined(PS_DEBUG_PRINT_SetLocalEssentialTreeToGlobalTreeLong)
    PS_PRINT_VARIABLE(tp_glb.size());
    PARTICLE_SIMULATOR_PRINT_MSG("***** END SetLocalEssentialTreeToGlobalTreeLong *****", 0);
#endif
}

template <class Ttp, class Tepj, typename Tmorton_key>
inline void SetLocalEssentialTreeToGlobalTreeShort(const ReallocatableArray<Tepj> &epj_org, const S32 n_loc, ReallocatableArray<Ttp> &tp_glb,
                                                   const Tmorton_key &morton_key, const bool flag_reuse = false) {
    PS_COMPILE_TIME_MESSAGE("SetLocalEssentialTreeToGlobalTreeShort");
    const S32 n_loc_ep = epj_org.size();
    tp_glb.resizeNoInitialize(n_loc_ep);
    if (!flag_reuse) {
        PS_OMP_PARALLEL_FOR
        for (S32 i = n_loc; i < n_loc_ep; i++) {
            tp_glb[i].setFromEP(epj_org[i], i, morton_key);
        }
    }
}

// for long force
template <class Ttp, class Tepj, class Tspj, class Tmorton_key>
inline void SetLocalEssentialTreeToLETTreeLong(const ReallocatableArray<Tepj> &epj_org, const ReallocatableArray<Tspj> &spj_org, const S32 n_loc,
                                               ReallocatableArray<Ttp> &tp_let, const Tmorton_key &morton_key, const bool flag_reuse = false) {
    const S32 n_let_ep = epj_org.size() - n_loc;
    assert(n_let_ep >= 0);
    const S32 n_let_ep_sp = n_let_ep + spj_org.size();
    tp_let.resizeNoInitialize(n_let_ep_sp);
    if (!flag_reuse) {
        PS_OMP_PARALLEL {
            PS_OMP_FOR
            for (S32 i = 0; i < n_let_ep; i++) {
                const S32 i_src = i + n_loc;
                tp_let[i].setFromEP(epj_org[i_src], i_src, morton_key);
            }
            PS_OMP_FOR
            for (S32 i = n_let_ep; i < n_let_ep_sp; i++) {
                const S32 i_src = i - n_let_ep;
                tp_let[i].setFromSP(spj_org[i_src], i_src, morton_key);
            }
        }
    }
}

template <class Ttp, class Tepj, class Tmorton_key>
inline void SetLocalEssentialTreeToLETTreeShort(const ReallocatableArray<Tepj> &epj_org, const S32 n_loc, ReallocatableArray<Ttp> &tp_let,
                                                const Tmorton_key &morton_key, const bool flag_reuse = false) {
    const S32 n_let_ep = epj_org.size() - n_loc;
    assert(n_let_ep >= 0);
    tp_let.resizeNoInitialize(n_let_ep);
    if (!flag_reuse) {
        PS_OMP_PARALLEL {
            PS_OMP_FOR
            for (S32 i = 0; i < n_let_ep; i++) {
                const S32 i_src = i + n_loc;
                tp_let[i].setFromEP(epj_org[i_src], i_src, morton_key);
            }
        }
    }
}

template <class Ttc>
inline void SetOuterBoxGlobalTreeForLongCutoffRecursive(Ttc tc[], const S32 n_leaf_limit, const F64 r_cut, const S32 adr_tc, const F64 tc_hlen,
                                                        const F64vec tc_cen) {
    F64 child_hlen = tc_hlen * 0.5;
    for (S32 i = 0; i < N_CHILDREN; i++) {
        F64vec child_cen = tc_cen + SHIFT_CENTER[i] * tc_hlen;
        // tc[adr_tc+i].mom_.vertex_out_.high_ = child_cen + (child_hlen+r_cut);
        // tc[adr_tc+i].mom_.vertex_out_.low_  = child_cen - (child_hlen+r_cut);
        tc[adr_tc + i].geo_.vertex_out_.high_ = child_cen + (child_hlen + r_cut);
        tc[adr_tc + i].geo_.vertex_out_.low_ = child_cen - (child_hlen + r_cut);
        S32 child_adr_tc = tc[adr_tc + i].adr_tc_;
        if (tc[adr_tc + i].n_ptcl_ <= 0)
            continue;
        else if (tc[adr_tc + i].isLeaf(n_leaf_limit))
            continue;
        else {
            SetOuterBoxGlobalTreeForLongCutoffRecursive(tc, n_leaf_limit, r_cut, child_adr_tc, child_hlen, child_cen);
        }
    }
}

template <class Ttc, class Tepj>
inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongCutoff, Ttc tc[], Tepj epj[], const S32 n_leaf_limit, const F64 tc_hlen,
                                                  const F64vec tc_cen) {
    F64 r_cut = epj[0].getRSearch();
    // tc[0].mom_.vertex_out_.high_ = tc_cen + (tc_hlen+r_cut);
    // tc[0].mom_.vertex_out_.low_  = tc_cen - (tc_hlen+r_cut);
    tc[0].geo_.vertex_out_.high_ = tc_cen + (tc_hlen + r_cut);
    tc[0].geo_.vertex_out_.low_ = tc_cen - (tc_hlen + r_cut);
    // for(S32 i=1; i<N_CHILDREN; i++) tc[i].mom_.vertex_out_.init();
    for (S32 i = 1; i < N_CHILDREN; i++) tc[i].geo_.vertex_out_.init();
    if (tc[0].n_ptcl_ < 0 || tc[0].isLeaf(n_leaf_limit)) return;
    S32 adr_tc = N_CHILDREN;
    SetOuterBoxGlobalTreeForLongCutoffRecursive(tc, n_leaf_limit, r_cut, adr_tc, tc_hlen, tc_cen);
}
template <class Ttc, class Tepj>
inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLong, Ttc tc[], Tepj epj[], const S32 n_leaf_limit, const F64 tc_hlen,
                                                  const F64vec tc_cen) {
    // do nothing
}
template <class Ttc, class Tepj>
inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongScatter, Ttc tc[], Tepj epj[], const S32 n_leaf_limit, const F64 tc_hlen,
                                                  const F64vec tc_cen) {
    // do nothing
}
template <class Ttc, class Tepj>
inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongSymmetry, Ttc tc[], Tepj epj[], const S32 n_leaf_limit, const F64 tc_hlen,
                                                  const F64vec tc_cen) {
    // do nothing
}
template <class Ttc, class Tepj>
inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortGather, Ttc tc[], Tepj epj[], const S32 n_leaf_limit, const F64 tc_hlen,
                                                  const F64vec tc_cen) {
    // do nothing
}
template <class Ttc, class Tepj>
inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortScatter, Ttc tc[], Tepj epj[], const S32 n_leaf_limit, const F64 tc_hlen,
                                                  const F64vec tc_cen) {
    // do nothing
}
template <class Ttc, class Tepj>
inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortSymmetry, Ttc tc[], Tepj epj[], const S32 n_leaf_limit, const F64 tc_hlen,
                                                  const F64vec tc_cen) {
    // do nothing
}

// it may works only if tp has one component.
template <class Tsys, class Ttp, class Tepi, class Tepj>
inline void CopyFpToEpSortedLocalTree(const Tsys &sys, const ReallocatableArray<Ttp> &tp, ReallocatableArray<Tepi> &epi_sorted,
                                      ReallocatableArray<Tepj> &epj_sorted) {
    const S32 n_loc = sys.getNumberOfParticleLocal();
    for (S32 i = 0; i < n_loc; i++) {
        const S32 adr = tp[i].adr_ptcl_;
        epi_sorted[i].copyFromFP(sys[adr]);
        epj_sorted[i].copyFromFP(sys[adr]);
    }
}

template <class Tcomm>
void MakeCommTableFor2StepCommuniction(Tcomm &comm_table, const ReallocatableArray<S32> &n_ep_send_per_proc_1st,
                                       const ReallocatableArray<S32> &n_image_per_proc_1st, const ReallocatableArray<F64vec> &shift_per_image_1st,
                                       const ReallocatableArray<S32> &n_ep_send_per_image_1st, const ReallocatableArray<S32> &adr_ep_send_1st,
                                       const ReallocatableArray<S32> &n_ep_recv_per_proc_1st, const ReallocatableArray<S32> &n_ep_send_per_proc_2nd,
                                       const ReallocatableArray<S32> &n_image_per_proc_2nd, const ReallocatableArray<F64vec> &shift_per_image_2nd,
                                       const ReallocatableArray<S32> &n_ep_send_per_image_2nd, const ReallocatableArray<S32> &adr_ep_send_2nd,
                                       const ReallocatableArray<S32> &n_ep_recv_per_proc_2nd, const S32 n_proc) {
    ReallocatableArray<S32> n_disp_ep_send_per_proc(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    n_disp_ep_send_per_proc[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_ep_send_per_proc[i + 1] = n_disp_ep_send_per_proc[i] + n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
    }
    ReallocatableArray<S32> n_disp_image_per_proc(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_image_per_proc_1st(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_image_per_proc_2nd(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    n_disp_image_per_proc[0] = n_disp_image_per_proc_1st[0] = n_disp_image_per_proc_2nd[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_image_per_proc[i + 1] = n_disp_image_per_proc[i] + n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
        n_disp_image_per_proc_1st[i + 1] = n_disp_image_per_proc_1st[i] + n_image_per_proc_1st[i];
        n_disp_image_per_proc_2nd[i + 1] = n_disp_image_per_proc_2nd[i] + n_image_per_proc_2nd[i];
    }
    const S32 n_image_tot_1st = shift_per_image_1st.size();
    const S32 n_image_tot_2nd = shift_per_image_2nd.size();
    ReallocatableArray<S32> n_disp_ep_send_per_image_1st(n_image_tot_1st + 1, n_image_tot_1st + 1, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_disp_ep_send_per_image_2nd(n_image_tot_2nd + 1, n_image_tot_2nd + 1, MemoryAllocMode::Pool);
    n_disp_ep_send_per_image_1st[0] = n_disp_ep_send_per_image_2nd[0] = 0;
    for (S32 i = 0; i < n_image_tot_1st; i++) {
        n_disp_ep_send_per_image_1st[i + 1] = n_disp_ep_send_per_image_1st[i] + n_ep_send_per_image_1st[i];
    }
    for (S32 i = 0; i < n_image_tot_2nd; i++) {
        n_disp_ep_send_per_image_2nd[i + 1] = n_disp_ep_send_per_image_2nd[i] + n_ep_send_per_image_2nd[i];
    }
    const S32 n_image_tot = n_disp_image_per_proc_1st[n_proc] + n_disp_image_per_proc_2nd[n_proc];
    comm_table.shift_per_image_.resizeNoInitialize(n_image_tot);
    comm_table.n_ep_per_image_.resizeNoInitialize(n_image_tot);
    const S32 n_send_tot = n_disp_ep_send_per_proc[n_proc];
    comm_table.adr_ep_send_.resizeNoInitialize(n_send_tot);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        S32 n_ep_cnt = 0;
        S32 n_image_cnt = 0;
        const S32 ep_head = n_disp_ep_send_per_proc[i];
        const S32 image_head = n_disp_image_per_proc[i];
        const S32 image_head_1st = n_disp_image_per_proc_1st[i];
        const S32 image_end_1st = n_disp_image_per_proc_1st[i + 1];
        for (S32 j = image_head_1st; j < image_end_1st; j++, n_image_cnt++) {
            const S32 ep_head_1st = n_disp_ep_send_per_image_1st[j];
            const S32 ep_end_1st = n_disp_ep_send_per_image_1st[j + 1];
            comm_table.shift_per_image_[image_head + n_image_cnt] = shift_per_image_1st[j];
            comm_table.n_ep_per_image_[image_head + n_image_cnt] = n_ep_send_per_image_1st[j];
            for (S32 k = ep_head_1st; k < ep_end_1st; k++, n_ep_cnt++) {
                comm_table.adr_ep_send_[ep_head + n_ep_cnt] = adr_ep_send_1st[k];
            }
        }
        const S32 image_head_2nd = n_disp_image_per_proc_2nd[i];
        const S32 image_end_2nd = n_disp_image_per_proc_2nd[i + 1];
        for (S32 j = image_head_2nd; j < image_end_2nd; j++, n_image_cnt++) {
            const S32 ep_head_2nd = n_disp_ep_send_per_image_2nd[j];
            const S32 ep_end_2nd = n_disp_ep_send_per_image_2nd[j + 1];
            comm_table.shift_per_image_[image_head + n_image_cnt] = shift_per_image_2nd[j];
            comm_table.n_ep_per_image_[image_head + n_image_cnt] = n_ep_send_per_image_2nd[j];
            for (S32 k = ep_head_2nd; k < ep_end_2nd; k++, n_ep_cnt++) {
                comm_table.adr_ep_send_[ep_head + n_ep_cnt] = adr_ep_send_2nd[k];
            }
        }
    }
    comm_table.n_ep_send_.resizeNoInitialize(n_proc);
    comm_table.n_ep_recv_.resizeNoInitialize(n_proc);
    comm_table.n_image_per_proc_.resizeNoInitialize(n_proc);
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        comm_table.n_ep_send_[i] = n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
        comm_table.n_ep_recv_[i] = n_ep_recv_per_proc_1st[i] + n_ep_recv_per_proc_2nd[i];
        comm_table.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
    }
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        comm_table.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
    }
    comm_table.n_ep_send_tot_ = comm_table.adr_ep_send_.size();
    S32 n_ep_recv_tot_tmp = 0;
PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp))
for (S32 i = 0; i < n_proc; i++) {
    n_ep_recv_tot_tmp += comm_table.n_ep_recv_[i];
}
comm_table.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
}
}  // namespace ParticleSimulator
