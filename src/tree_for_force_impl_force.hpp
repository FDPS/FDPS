#pragma once
#include "tree_walk.hpp"

namespace ParticleSimulator {
struct KERNEL_TYPE_INDEX {};
struct KERNEL_TYPE_PTCL {};
struct KERNEL_TYPE_NORMAL {};
template <typename Tforce>
void CalcForceMultiWalkInitialize(ReallocatableArray<Tforce>& force_sorted, const S32 n_loc_tot, const bool tag_max, const bool clear) {
    if (tag_max <= 0) {
        PARTICLE_SIMULATOR_PRINT_ERROR("tag_max is illegal. In currente version, tag_max must be 1");
        Abort(-1);
    }
    force_sorted.resizeNoInitialize(n_loc_tot);
    if (clear) {
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_loc_tot; i++) {
            force_sorted[i].clear();
        }
    }
}

template <typename Tforce>
void CalcForceInitialize(ReallocatableArray<Tforce>& force_sorted, const S32 n_loc_tot, const bool tag_max, const bool clear) {
    CalcForceMultiWalkInitialize(force_sorted, n_loc_tot, tag_max, clear);
}

///////////////////////
// calcForceMakingGroup
template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tkernel_type, typename Tfunc_dispatch, typename Tfunc_retrieve>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingGroupMultiWalk(
    Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit, const INTERACTION_LIST_MODE list_mode,
    const bool copy_force, const bool clear) {
    F64 t0 = GetWtime();
    makeIPGroupOnly(list_mode);
    time_profile_.construct_ipg += GetWtime() - t0;
    S32 ret = 0;
    ret = calcForceMultiWalk<Tkernel_type>(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, clear);
    return ret;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tfunc_ep_sp>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingGroupLong(
    Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear_force, const bool copy_force, const INTERACTION_LIST_MODE list_mode,
    const TREE_TYPE tree_type) {
#if defined(PS_DEBUG_PRINT_calcForceMakingGroup)
    PARTICLE_SIMULATOR_PRINT_MSG_B("**** START calcForceMakingGroupLong ****", 0);
#endif
    F64 t0 = GetWtime();
    makeIPGroupOnly(list_mode);
    time_profile_.construct_ipg += GetWtime() - t0;
#if 1
    auto n_ipg_this_rnd_max = 100;
    auto n_loop = (ipg_.size() + n_ipg_this_rnd_max - 1) / n_ipg_this_rnd_max;
    if (list_mode == MAKE_LIST_FOR_REUSE || list_mode == REUSE_LIST) {
        n_ipg_this_rnd_max = ipg_.size();
        n_loop = 1;
    }
    // std::cout<<"n_loop= "<<n_loop<<std::endl;
    for (int i = 0; i < n_loop; i++) {
        auto adr_ipg_first = i * n_ipg_this_rnd_max;
        auto n_ipg_this_rnd = std::min(n_ipg_this_rnd_max, ipg_.size() - adr_ipg_first);
        F64 t1 = GetWtimeNoBarrier();
        makeJPListOnly(adr_ipg_first, n_ipg_this_rnd, list_mode, tree_type);
        F64 t2 = GetWtimeNoBarrier();
        calcForceNoWalkPMM(pfunc_ep_ep, pfunc_ep_sp, clear_force, copy_force, tree_type, adr_ipg_first, n_ipg_this_rnd);
        F64 t3 = GetWtimeNoBarrier();
        // time_profile_.interaction_list_construction += t2 - t1;
        time_profile_.construct_interaction_list += t2 - t1;
        time_profile_.calculate_force += t3 - t2;
        // time_profile_.hoge += t3 - t1;
    }
#else
    makeJPListOnly(0, ipg_.size(), list_mode, tree_type);
    calcForceNoWalkPMM(pfunc_ep_ep, pfunc_ep_sp, clear_force, copy_force, tree_type, 0, ipg_.size());
    const S32 n_thread = Comm::getNumberOfThread();
    for (int i = 0; i < n_thread; i++) {
        // epj_for_force_[i].freeMem();
        // spj_for_force_[i].freeMem();
    }
#endif
#if defined(PS_DEBUG_PRINT_calcForceMakingGroup)
    PARTICLE_SIMULATOR_PRINT_MSG("**** END calcForceMakingGroupLong ****", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_ep_ep>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingGroupShort(
    Tfunc_ep_ep pfunc_ep_ep, const bool clear_force, const bool copy_force, const INTERACTION_LIST_MODE list_mode, const TREE_TYPE tree_type) {
#if defined(PS_DEBUG_PRINT_calcForceMakingGroup)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START calcForceMakingGroupShort", 0);
#endif

    F64 t0 = GetWtime();
    makeIPGroupOnly(list_mode);
    time_profile_.construct_ipg += GetWtime() - t0;

#if 1
    auto n_ipg_this_rnd_max = 10;
    auto n_loop = (ipg_.size() + n_ipg_this_rnd_max - 1) / n_ipg_this_rnd_max;
    if (list_mode == MAKE_LIST_FOR_REUSE || list_mode == REUSE_LIST) {
        n_ipg_this_rnd_max = ipg_.size();
        n_loop = 1;
    }
    // std::cout<<"n_loop= "<<n_loop<<std::endl;
    for (int i = 0; i < n_loop; i++) {
        auto adr_ipg_first = i * n_ipg_this_rnd_max;
        auto n_ipg_this_rnd = std::min(n_ipg_this_rnd_max, ipg_.size() - adr_ipg_first);
        // std::cout<<"n_ipg_this_rnd= "<<n_ipg_this_rnd<<std::endl;
        F64 t1 = GetWtimeNoBarrier();
        makeJPListOnly(adr_ipg_first, n_ipg_this_rnd, list_mode, tree_type);
        F64 t2 = GetWtimeNoBarrier();
        calcForceNoWalk(pfunc_ep_ep, clear_force, adr_ipg_first, n_ipg_this_rnd);
        F64 t3 = GetWtimeNoBarrier();
        time_profile_.construct_interaction_list += t2 - t1;
        time_profile_.calculate_force += t3 - t2;
    }
#else

    t0 = GetWtime();
    makeJPListOnly(0, ipg_.size(), list_mode, tree_type);
    time_profile_.construct_interaction_list += GetWtime() - t0;
    t0 = GetWtime();
    calcForceNoWalkOrg(pfunc_ep_ep, clear_force);
    time_profile_.calculate_force += GetWtime() - t0;
    // const S32 n_thread = Comm::getNumberOfThread();

#endif
#if defined(PS_DEBUG_PRINT_calcForceMakingGroup)
    PARTICLE_SIMULATOR_PRINT_MSG("**** END calcForceMakingGroupShort", 0);
#endif
}
// calcForceMakingGroup
///////////////////////

#if 1
template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tkernel_type, typename Tfunc_dispatch, typename Tfunc_retrieve>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                                                                                          Tfunc_retrieve pfunc_retrieve,
                                                                                                          const S32 tag_max, const S32 n_walk_limit,
                                                                                                          const INTERACTION_LIST_MODE list_mode,
                                                                                                          const bool clear) {
    S32 ret = 0;
    CalcForceMultiWalkInitialize(force_sorted_, n_loc_tot_, tag_max, clear);
    // const auto wtime_offset = GetWtime();
    ret = calcForceMultiWalkImpl<Tkernel_type>(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, clear);
    // time_profile_.calculate_force += GetWtime() - wtime_offset;
    return ret;
}
#else
//////////// Walk+Force, Kernel:Index ////////////
template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_dispatch, typename Tfunc_retrieve>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMultiWalkIndex(
    Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit, const INTERACTION_LIST_MODE list_mode,
    const bool clear) {
    S32 ret = 0;
    CalcForceMultiWalkInitialize(force_sorted_, n_loc_tot_, tag_max, clear);
    const auto wtime_offset = GetWtime();
    ret = calcForceMultiWalkImpl<KERNEL_TYPE_INDEX>(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, clear);
    // time_profile_.calc_force += GetWtime() - wtime_offset;
    return ret;
}

//////////// Walk+Force, Kernel:Ptcl //////////////
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_dispatch, class Tfunc_retrieve>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMultiWalkPtcl(
    Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit, const INTERACTION_LIST_MODE list_mode,
    const bool clear) {
    S32 ret = 0;
    CalcForceMultiWalkInitialize(force_sorted_, n_loc_tot_, tag_max, clear);
    const F64 time_offset = GetWtime();
    ret = calcForceMultiWalkImpl<KERNEL_TYPE_PTCL>(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, clear);
    // time_profile_.calc_force += GetWtime() - time_offset;
    return ret;
}
#endif

template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_ep_ep, typename Tfunc_ep_sp>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceNoWalkPMM(
    Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear_force, const bool copy_force, const TREE_TYPE tree_type,
    const PS::S32 adr_ipg_first, const PS::S32 n_ipg_this_rnd) {
    //F64 time_offset = GetWtime();
    if (tree_type == GLOBAL_TREE) {
        // std::cerr << "CHECK AA" << std::endl;
        CalcForceNoWalkPMM<TreeCellGlb, IPGroup<typename TSM::ipg_type>, InteractionList, Tepi, Tepj, Tspj, Tforce, Tfunc_ep_ep, Tfunc_ep_sp>(
            tc_glb_, ipg_, interaction_list_glb_, n_loc_tot_, n_let_sp_, epi_sorted_, epj_sorted_, spj_sorted_, force_sorted_, n_group_limit_,
            pfunc_ep_ep, pfunc_ep_sp, adr_ipg_first, n_ipg_this_rnd, clear_force, false, n_interaction_ep_ep_local_, n_interaction_ep_sp_local_,
            time_profile_);
        // std::cerr << "CHECK BB" << std::endl;
        /*
        CalcForceNoWalkPMM<TreeCellGlb, IPGroup<typename TSM::ipg_type>, InteractionList, Tepi, Tepj, Tspj, Tforce, Tfunc_ep_ep, Tfunc_ep_sp>(
            tc_glb_, ipg_, interaction_list_glb_, n_loc_tot_, n_let_sp_, epi_sorted_, epj_sorted_, spj_sorted_, force_sorted_, n_group_limit_,
            pfunc_ep_ep, pfunc_ep_sp, 0, ipg_.size(), clear_force, false, n_interaction_ep_ep_local_, n_interaction_ep_sp_local_, time_profile_);
            */
    } else if (tree_type == LOCAL_TREE) {
        /*
          CalcForceNoWalk
          <TreeCellLoc,
          IPGroup< typename TSM::ipg_type >,
          InteractionList,
          Tepi, Tepj, Tspj, Tforce,
          Tfunc_ep_ep, Tfunc_ep_sp, N_POP>
          (tc_loc_, ipg_[N_POP],
          interaction_list_loc_[N_POP],
          n_loc_tot_, 0,
          epi_sorted_, epj_sorted_,
          spj_sorted_, force_sorted_,
          n_group_limit_,
          pfunc_ep_ep, pfunc_ep_sp,
          N_POP, 0, ipg_[N_POP].size(),
          clear_force, false,
          n_epi_processed_local_,
          n_epj_processed_local_,
          n_spj_processed_local_,
          n_interaction_ep_ep_local_,
          n_interaction_ep_sp_local_,
          time_profile_,
          n_thread);
        */
    } else if (tree_type == LET_TREE) {
        /*
          CalcForceNoWalk
          <TreeCellGlb,
          IPGroup< typename TSM::ipg_type >,
          InteractionList,
          Tepi, Tepj, Tspj, Tforce,
          Tfunc_ep_ep, Tfunc_ep_sp, N_POP>
          (tc_let_, ipg_[N_POP],
          interaction_list_let_[N_POP],
          n_loc_tot_, n_let_sp_,
          epi_sorted_, epj_sorted_,
          spj_sorted_, force_sorted_,
          n_group_limit_,
          pfunc_ep_ep, pfunc_ep_sp,
          N_POP, 0, ipg_[N_POP].size(),
          clear_force, false,
          n_epi_processed_local_,
          n_epj_processed_local_,
          n_spj_processed_local_,
          n_interaction_ep_ep_local_,
          n_interaction_ep_sp_local_,
          time_profile_,
          n_thread);
        */
    }
    if (copy_force) {
        //const F64 offset_copy_original_order = GetWtime();
        // std::cerr << "CHECK CC" << std::endl;
        copyForceOriginalOrder();
        // std::cerr << "CHECK DD" << std::endl;
        //time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
    }
    // time_profile_.calc_force += GetWtime() - time_offset;
}

///////////////////////////////////////////////////
//
// FUNCTIONS OF WALK+FORCE WITH DOUBLE BUFFERING
//
///////////////////////////////////////////////////
// RELEASE V8 functions are
// * calcForceMultiWalkIndex(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit, const bool
// flag_keep_list, const bool clear);
// * calcForceMultiWalkPtcl(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit, const bool
// flag_keep_list, const bool clear);
// * calcForceNoWalkForMultiWalkIndexImpl(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 n_walk_limit, const bool clear);
// * calcForceNoWalkForMultiWalkPtclImpl(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 n_walk_limit, const bool clear);
// *************************** //
// *************************** //
// below functions are not implemented
// * calcForce(Tfunc_ep_ep pfunc_ep_ep, const bool clear);
// * calcForceOnly(Tfunc_ep_ep pfunc_ep_ep, const S32 adr_ipg, const bool clear);
// * makeInteractionListImpl(TagForceShort, const S32 adr_ipg, const bool clear);
// * calcForce(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear);
// * calcForceOnly(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const S32 adr_ipg, const bool clear);
// * makeInteractionListImpl(TagForceLong, const S32 adr_ipg, const bool clear);
// * calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep, const bool clear);
// * calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear);
// * makeInteractionList(const S32 adr_ipg, const bool clear); // this just call both modes of makeInteractionListImpl

//////////// Walk+Force, Kernel:Ptcl and Index, List:Index //////////////
template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename KERNEL_TYPE, typename Tfunc_dispatch, typename Tfunc_retrieve>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMultiWalkImpl(
    Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit, const INTERACTION_LIST_MODE list_mode,
    const bool clear) {
    //const F64 offset_core = GetWtime();
    if (clear) {
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_loc_tot_; i++) {
            force_sorted_[i].clear();
        }
    }
    // send all epj and spj
    Tepi** epi_dummy = nullptr;
    S32* n_epi_dummy = nullptr;
    S32** id_epj_dummy = nullptr;
    S32* n_epj_dummy = nullptr;
    if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_INDEX>) {
        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
            S32** id_spj_dummy = nullptr;
            S32* n_spj_dummy = nullptr;
            pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy, (const S32**)id_epj_dummy, n_epj_dummy, (const S32**)id_spj_dummy, n_spj_dummy,
                           epj_sorted_.getPointer(), epj_sorted_.size(), spj_sorted_.getPointer(), spj_sorted_.size(), true);
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
            pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy, (const S32**)id_epj_dummy, n_epj_dummy, epj_sorted_.getPointer(),
                           epj_sorted_.size(), true);
        } else {
            static_assert([] { return false; }());
        }
    }
    S32 ret = 0;
    S32 tag = 0;
    const S32 n_thread = Comm::getNumberOfThread();
    force_sorted_.resizeNoInitialize(n_loc_tot_);
    const S32 n_ipg = ipg_.size();
    if (n_ipg <= 0) return 0;
    // n_walk_local_ += n_ipg;
    if (list_mode == MAKE_LIST_FOR_REUSE) {
        interaction_list_glb_.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list_glb_.n_disp_ep_.resizeNoInitialize(n_ipg + 1);
        interaction_list_glb_.adr_ep_.clearSize();
        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
            interaction_list_glb_.n_sp_.resizeNoInitialize(n_ipg);
            interaction_list_glb_.n_disp_sp_.resizeNoInitialize(n_ipg + 1);
            interaction_list_glb_.adr_sp_.clearSize();
        }
    }
    const S32 n_loop_max = n_ipg / n_walk_limit + ((n_ipg % n_walk_limit) == 0 ? 0 : 1);
    ReallocatableArray<Tforce*> ptr_force_per_walk[2];
    ReallocatableArray<S32> n_epi_per_walk[2];
    ReallocatableArray<Tepi*> ptr_epi_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_epj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);
    ReallocatableArray<S32*> ptr_adr_epj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_epj_disp_per_thread[n_thread];
    ReallocatableArray<S32> n_spj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);
    ReallocatableArray<S32*> ptr_adr_spj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);
    ReallocatableArray<S32> n_spj_disp_per_thread[n_thread];                                        // for long search case only
    ReallocatableArray<Tepj*> ptr_epj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);  // for ptcl kernel
    ReallocatableArray<Tspj*> ptr_spj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);  // for ptcl kernel
    ptr_force_per_walk[0].initialize(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);            // array of pointer *[n_walk]
    ptr_force_per_walk[1].initialize(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);            // array of pointer *[n_walk]
    n_epi_per_walk[0].initialize(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);
    n_epi_per_walk[1].initialize(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);

    ReallocatableArray<Tepj> epj_for_force[n_thread];
    ReallocatableArray<Tspj> spj_for_force[n_thread];

    for (int i = 0; i < n_thread; i++) {
        n_epj_disp_per_thread[i].initialize(n_walk_limit + 1, n_walk_limit + 1, MemoryAllocMode::Pool);
    }
    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
        for (int i = 0; i < n_thread; i++) {
            n_spj_disp_per_thread[i].initialize(n_walk_limit + 1, n_walk_limit + 1, MemoryAllocMode::Pool);
        }
    }
    const auto adr_tree_sp_first = n_let_sp_;
    bool first_loop = true;
    S32 n_walk_prev = 0;
    S64 n_interaction_ep_ep_tmp = 0;
    S64 n_interaction_ep_sp_tmp = 0;
    const auto len_peri = p_domain_info_->getLenRootDomain();
    if (n_ipg > 0) {
        for (int wg = 0; wg < n_loop_max; wg++) {
            const S32 n_walk = n_ipg / n_loop_max + (((n_ipg % n_loop_max) > wg) ? 1 : 0);
            const S32 walk_grp_head = (n_ipg / n_loop_max) * wg + std::min((n_ipg % n_loop_max), wg);
            //const auto offset_calc_force__core__walk_tree = GetWtime();
            const S32 lane_now = wg % 2;
            const S32 lane_old = (wg + 1) % 2;
            if (list_mode == REUSE_LIST) {
#if 1
                if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL>) {
                    const S64 n_ep_head = interaction_list_glb_.n_disp_ep_[walk_grp_head];
                    const S64 n_ep_tail = interaction_list_glb_.n_disp_ep_[walk_grp_head + n_walk];
                    // epj_for_force_[0].resizeNoInitialize((n_ep_tail - n_ep_head));
                    epj_for_force[0].resizeNoInitialize((n_ep_tail - n_ep_head));
                    PS_OMP_PARALLEL_FOR
                    for (S32 jp = 0; jp < (n_ep_tail - n_ep_head); jp++) {
                        // epj_for_force_[0][jp] = epj_sorted_[interaction_list_glb_.adr_ep_[jp + n_ep_head]];
                        epj_for_force[0][jp] = epj_sorted_[interaction_list_glb_.adr_ep_[jp + n_ep_head]];
                    }
                    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                        const S64 n_sp_head = interaction_list_glb_.n_disp_sp_[walk_grp_head];
                        const S64 n_sp_tail = interaction_list_glb_.n_disp_sp_[walk_grp_head + n_walk];
                        // spj_for_force_[0].resizeNoInitialize((n_sp_tail - n_sp_head));
                        spj_for_force[0].resizeNoInitialize((n_sp_tail - n_sp_head));
                        PS_OMP_PARALLEL_FOR
                        for (S32 jp = 0; jp < (n_sp_tail - n_sp_head); jp++) {
                            // spj_for_force_[0][jp] = spj_sorted_[interaction_list_glb_.adr_sp_[jp + n_sp_head]];
                            spj_for_force[0][jp] = spj_sorted_[interaction_list_glb_.adr_sp_[jp + n_sp_head]];
                        }
                    }
                }
                PS_OMP_PARALLEL_FOR
                for (int iw = 0; iw < n_walk; iw++) {
                    const S32 iw_glb = walk_grp_head + iw;
                    const S32 first_adr_ip = ipg_[iw_glb].adr_ptcl_;
                    n_epi_per_walk[lane_now][iw] = ipg_[iw_glb].n_ptcl_;
                    ptr_epi_per_walk[iw] = epi_sorted_.getPointer(first_adr_ip);
                    ptr_force_per_walk[lane_now][iw] = force_sorted_.getPointer(first_adr_ip);
                    n_epj_per_walk[iw] = interaction_list_glb_.n_ep_[iw_glb];
                    if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL>) {
                        const S32 offset_ep = interaction_list_glb_.n_disp_ep_[iw_glb] - interaction_list_glb_.n_disp_ep_[walk_grp_head];
                        // ptr_epj_per_walk[iw] = epj_for_force_[0].getPointer(offset_ep);
                        ptr_epj_per_walk[iw] = epj_for_force[0].getPointer(offset_ep);
                    } else if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_INDEX>) {
                        ptr_adr_epj_per_walk[iw] = interaction_list_glb_.adr_ep_.getPointer(interaction_list_glb_.n_disp_ep_[iw_glb]);
                    }
                    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                        n_spj_per_walk[iw] = interaction_list_glb_.n_sp_[iw_glb];
                        if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL>) {
                            const S32 offset_sp = interaction_list_glb_.n_disp_sp_[iw_glb] - interaction_list_glb_.n_disp_sp_[walk_grp_head];
                            // ptr_spj_per_walk[iw] = spj_for_force_[0].getPointer(offset_sp);
                            ptr_spj_per_walk[iw] = spj_for_force[0].getPointer(offset_sp);
                        } else if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_INDEX>) {
                            ptr_adr_spj_per_walk[iw] = interaction_list_glb_.adr_sp_.getPointer(interaction_list_glb_.n_disp_sp_[iw_glb]);
                        }
                    }
                }
            }
#else
                if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL>) {
                    const S64 n_ep_head = interaction_list_glb_.n_disp_ep_[walk_grp_head];
                    const S64 n_ep_tail = interaction_list_glb_.n_disp_ep_[walk_grp_head + n_walk];
                    epj_for_force_[0].resizeNoInitialize((n_ep_tail - n_ep_head));
                    PS_OMP_PARALLEL_FOR
                    for (S32 jp = 0; jp < (n_ep_tail - n_ep_head); jp++) {
                        epj_for_force_[0][jp] = epj_sorted_[interaction_list_glb_.adr_ep_[jp + n_ep_head]];
                    }
                    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                        const S64 n_sp_head = interaction_list_glb_.n_disp_sp_[walk_grp_head];
                        const S64 n_sp_tail = interaction_list_glb_.n_disp_sp_[walk_grp_head + n_walk];
                        spj_for_force_[0].resizeNoInitialize((n_sp_tail - n_sp_head));
                        PS_OMP_PARALLEL_FOR
                        for (S32 jp = 0; jp < (n_sp_tail - n_sp_head); jp++) {
                            spj_for_force_[0][jp] = spj_sorted_[interaction_list_glb_.adr_sp_[jp + n_sp_head]];
                        }
                    }
                    PS_OMP_PARALLEL_FOR
                    for (int iw = 0; iw < n_walk; iw++) {
                        const S32 iw_glb = walk_grp_head + iw;
                        const S32 first_adr_ip = ipg_[iw_glb].adr_ptcl_;
                        n_epi_per_walk[lane_now][iw] = ipg_[iw_glb].n_ptcl_;
                        ptr_epi_per_walk[iw] = epi_sorted_.getPointer(first_adr_ip);
                        ptr_force_per_walk[lane_now][iw] = force_sorted_.getPointer(first_adr_ip);
                        n_epj_per_walk[iw] = interaction_list_glb_.n_ep_[iw_glb];
                        const S32 offset_ep = interaction_list_glb_.n_disp_ep_[iw_glb] - interaction_list_glb_.n_disp_ep_[walk_grp_head];
                        ptr_epj_per_walk[iw] = epj_for_force_[0].getPointer(offset_ep);
                        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                            n_spj_per_walk[iw] = interaction_list_glb_.n_sp_[iw_glb];
                            const S32 offset_sp = interaction_list_glb_.n_disp_sp_[iw_glb] - interaction_list_glb_.n_disp_sp_[walk_grp_head];
                            ptr_spj_per_walk[iw] = spj_for_force_[0].getPointer(offset_sp);
                        }
                    }
                } else if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_INDEX>) {
                    PS_OMP_PARALLEL_FOR
                    for (int iw = 0; iw < n_walk; iw++) {
                        const S32 iw_glb = walk_grp_head + iw;
                        const S32 first_adr_ip = ipg_[iw_glb].adr_ptcl_;
                        n_epi_per_walk[lane_now][iw] = ipg_[iw_glb].n_ptcl_;
                        ptr_epi_per_walk[iw] = epi_sorted_.getPointer(first_adr_ip);
                        ptr_force_per_walk[lane_now][iw] = force_sorted_.getPointer(first_adr_ip);
                        n_epj_per_walk[iw] = interaction_list_glb_.n_ep_[iw_glb];
                        ptr_adr_epj_per_walk[iw] = interaction_list_glb_.adr_ep_.getPointer(interaction_list_glb_.n_disp_ep_[iw_glb]);
                        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                            n_spj_per_walk[iw] = interaction_list_glb_.n_sp_[iw_glb];
                            ptr_adr_spj_per_walk[iw] = interaction_list_glb_.adr_sp_.getPointer(interaction_list_glb_.n_disp_sp_[iw_glb]);
                        }
                    }
                }  // end of KERNEL TYPE scope
            }
#endif
            else if (list_mode != REUSE_LIST) {
                PS_OMP_PARALLEL {
                    const S32 ith = Comm::getThreadNum();
                    ReallocatableArray<S32> iwloc2iw(n_walk_limit, n_walk_limit, MemoryAllocMode::Pool);
                    S32 n_walk_loc = 0;
                    S32 n_ep_cum = 0;
                    S32 n_sp_cum = 0;
                    n_epj_disp_per_thread[ith][0] = 0;
                    adr_epj_for_force_[ith].clearSize();
                    adr_ipg_for_force_[ith].clearSize();
                    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                        n_spj_disp_per_thread[ith][0] = 0;
                        adr_spj_for_force_[ith].clearSize();
                    }
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
#pragma omp for schedule(dynamic, 4) reduction(+ : n_interaction_ep_ep_tmp) reduction(+ : n_interaction_ep_sp_tmp)
#endif
                    for (int iw = 0; iw < n_walk; iw++) {
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_for_force_[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_;
                        const S32 ith = Comm::getThreadNum();
                        n_epi_per_walk[lane_now][iw] = ipg_[id_ipg].n_ptcl_;
                        ptr_epi_per_walk[iw] = epi_sorted_.getPointer(first_adr_ip);
                        ptr_force_per_walk[lane_now][iw] = force_sorted_.getPointer(first_adr_ip);
                        TargetBox<TSM> target_box;
                        if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
                            target_box.set(ipg_[id_ipg], i_cut_, morton_key_);
                        } else {
                            target_box.set(ipg_[id_ipg]);
                        }
                        // target_box.set(ipg_[id_ipg]);
                        S32 adr_tc = 0;
                        auto t0 = GetWtime();
                        if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
                            MakeListUsingTreeRecursiveTopPMM<TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, MortonKey<typename TSM::mkey_type>,
                                                             TagChopLeafFalse, TagChopNonleafFalse, TagCopyInfoCloseWithTpAdrptcl,
                                                             CALC_DISTANCE_TYPE>(
                                tc_glb_, adr_tc, tp_glb_, epj_sorted_, adr_epj_for_force_[ith], spj_sorted_, adr_spj_for_force_[ith], target_box,
                                n_leaf_limit_, lev_leaf_limit_, adr_tree_sp_first, len_peri, theta_, lev_pm_cell_, morton_key_);
                        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                            MakeListUsingTreeRecursiveTop<TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl,
                                                          CALC_DISTANCE_TYPE>(tc_glb_, adr_tc, tp_glb_, epj_sorted_, adr_epj_for_force_[ith],
                                                                              spj_sorted_, adr_spj_for_force_[ith], target_box, n_leaf_limit_,
                                                                              adr_tree_sp_first, len_peri, theta_);
                        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
                            MakeListUsingTreeShortRecursiveTop<TSM, TreeCellGlb, TreeParticle, Tepj, TagChopLeafFalse, TagCopyInfoCloseNoSp,
                                                               CALC_DISTANCE_TYPE>(tc_glb_, adr_tc, tp_glb_, epj_sorted_, adr_epj_for_force_[ith],
                                                                                   target_box, n_leaf_limit_, len_peri);
                        }
                        // time_profile_.interaction_list_construction += GetWtime() - t0;
                        time_profile_.construct_interaction_list += GetWtime() - t0;
                        n_epj_per_walk[iw] = adr_epj_for_force_[ith].size() - n_ep_cum;
                        n_ep_cum = adr_epj_for_force_[ith].size();
                        n_epj_disp_per_thread[ith][n_walk_loc + 1] = n_ep_cum;
                        n_interaction_ep_ep_tmp += ((S64)n_epj_per_walk[iw] * (S64)n_epi_per_walk[lane_now][iw]);
                        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                            n_spj_per_walk[iw] = adr_spj_for_force_[ith].size() - n_sp_cum;
                            n_sp_cum = adr_spj_for_force_[ith].size();
                            n_spj_disp_per_thread[ith][n_walk_loc + 1] = n_sp_cum;
                            n_interaction_ep_sp_tmp += ((S64)n_spj_per_walk[iw] * (S64)n_epi_per_walk[lane_now][iw]);
                        }
                        iwloc2iw[n_walk_loc] = iw;
                        n_walk_loc++;
                    }  // end of OMP for
                    if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL>) {
                        // epj_for_force_[ith].resizeNoInitialize(adr_epj_for_force_[ith].size());
                        epj_for_force[ith].resizeNoInitialize(adr_epj_for_force_[ith].size());
                        for (S32 iwloc = 0; iwloc < n_walk_loc; iwloc++) {
                            S32 iw = iwloc2iw[iwloc];
                            for (int i = 0; i < n_epj_per_walk[iw]; i++) {
                                S32 offset = n_epj_disp_per_thread[ith][iwloc];
                                S32 adr = adr_epj_for_force_[ith][offset + i];
                                // epj_for_force_[ith][offset + i] = epj_sorted_[adr];
                                epj_for_force[ith][offset + i] = epj_sorted_[adr];
                            }
                            // ptr_epj_per_walk[iw] = epj_for_force_[ith].getPointer(n_epj_disp_per_thread[ith][iwloc]);
                            ptr_epj_per_walk[iw] = epj_for_force[ith].getPointer(n_epj_disp_per_thread[ith][iwloc]);
                            if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                                // spj_for_force_[ith].resizeNoInitialize(adr_spj_for_force_[ith].size());
                                spj_for_force[ith].resizeNoInitialize(adr_spj_for_force_[ith].size());
                                for (int i = 0; i < n_spj_per_walk[iw]; i++) {
                                    S32 offset = n_spj_disp_per_thread[ith][iwloc];
                                    S32 adr = adr_spj_for_force_[ith][offset + i];
                                    // spj_for_force_[ith][offset + i] = spj_sorted_[adr];
                                    spj_for_force[ith][offset + i] = spj_sorted_[adr];
                                }
                                // ptr_spj_per_walk[iw] = spj_for_force_[ith].getPointer(n_spj_disp_per_thread[ith][iwloc]);
                                ptr_spj_per_walk[iw] = spj_for_force[ith].getPointer(n_spj_disp_per_thread[ith][iwloc]);
                            }
                        }
                    } else if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_INDEX>) {
                        for (S32 iwloc = 0; iwloc < n_walk_loc; iwloc++) {
                            S32 iw = iwloc2iw[iwloc];
                            ptr_adr_epj_per_walk[iw] = adr_epj_for_force_[ith].getPointer(n_epj_disp_per_thread[ith][iwloc]);
                            if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                                ptr_adr_spj_per_walk[iw] = adr_spj_for_force_[ith].getPointer(n_spj_disp_per_thread[ith][iwloc]);
                            }
                        }
                    }
                    if (list_mode == MAKE_LIST_FOR_REUSE) {
                        PS_OMP_FOR
                        for (int iw = 0; iw < n_walk; iw++) {
                            const S32 id_ipg = walk_grp_head + iw;
                            interaction_list_glb_.n_ep_[id_ipg] = n_epj_per_walk[iw];
                            if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                                interaction_list_glb_.n_sp_[id_ipg] = n_spj_per_walk[iw];
                            }
                        }
                        PS_OMP_SINGLE {
                            interaction_list_glb_.n_disp_ep_[0] = interaction_list_glb_.n_disp_sp_[0] = 0;
                            for (int iw = 0; iw < n_walk; iw++) {
                                const S32 id_ipg = walk_grp_head + iw;
                                interaction_list_glb_.n_disp_ep_[id_ipg + 1] =
                                    interaction_list_glb_.n_disp_ep_[id_ipg] + interaction_list_glb_.n_ep_[id_ipg];
                                if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                                    interaction_list_glb_.n_disp_sp_[id_ipg + 1] =
                                        interaction_list_glb_.n_disp_sp_[id_ipg] + interaction_list_glb_.n_sp_[id_ipg];
                                }
                            }
                            interaction_list_glb_.adr_ep_.resizeNoInitialize(interaction_list_glb_.n_disp_ep_[walk_grp_head + n_walk]);
                            if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                                interaction_list_glb_.adr_sp_.resizeNoInitialize(interaction_list_glb_.n_disp_sp_[walk_grp_head + n_walk]);
                            }
                        }
                        for (S32 j = 0; j < adr_ipg_for_force_[ith].size(); j++) {
                            const S32 adr_ipg = adr_ipg_for_force_[ith][j];
                            S32 adr_ep = interaction_list_glb_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_per_thread[ith][j];
                            const S32 k_ep_e = n_epj_disp_per_thread[ith][j + 1];
                            for (S32 k = k_ep_h; k < k_ep_e; k++, adr_ep++) {
                                interaction_list_glb_.adr_ep_[adr_ep] = adr_epj_for_force_[ith][k];
                            }
                            if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                                S32 adr_sp = interaction_list_glb_.n_disp_sp_[adr_ipg];
                                const S32 k_sp_h = n_spj_disp_per_thread[ith][j];
                                const S32 k_sp_e = n_spj_disp_per_thread[ith][j + 1];
                                for (S32 k = k_sp_h; k < k_sp_e; k++, adr_sp++) {
                                    interaction_list_glb_.adr_sp_[adr_sp] = adr_spj_for_force_[ith][k];
                                }
                            }
                        }
                    }
                }  // end of OMP parallel scope
            }
            //time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
            auto t0 = GetWtime();
            if (!first_loop) {
                ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[lane_old].getPointer(), ptr_force_per_walk[lane_old].getPointer());
            }  // retrieve
            if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_INDEX>) {
                if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                    pfunc_dispatch(0, n_walk, (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(),
                                   (const S32**)ptr_adr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(),
                                   (const S32**)ptr_adr_spj_per_walk.getPointer(), n_spj_per_walk.getPointer(), epj_sorted_.getPointer(),
                                   epj_sorted_.size(), spj_sorted_.getPointer(), spj_sorted_.size(), false);
                } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
                    pfunc_dispatch(0, n_walk, (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(),
                                   (const S32**)ptr_adr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(), epj_sorted_.getPointer(),
                                   epj_sorted_.size(), false);
                }
            } else if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL>) {
                if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                    ret += pfunc_dispatch(tag, n_walk, (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(),
                                          (const Tepj**)ptr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(),
                                          (const Tspj**)ptr_spj_per_walk.getPointer(), n_spj_per_walk.getPointer());
                } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
                    ret += pfunc_dispatch(tag, n_walk, (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(),
                                          (const Tepj**)ptr_epj_per_walk.getPointer(),
                                          n_epj_per_walk.getPointer());  // new
                }
            }
            time_profile_.calculate_force += GetWtime() - t0;
            first_loop = false;
            n_walk_prev = n_walk;
        }  // end of walk group loop
        auto t0 = GetWtime();
        ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[(n_loop_max + 1) % 2].getPointer(),
                              ptr_force_per_walk[(n_loop_max + 1) % 2].getPointer());
        time_profile_.calculate_force += GetWtime() - t0;
        n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
    }  // if(n_ipg > 0)
    else {
        ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
        n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
    }
    //time_profile_.calc_force__core += GetWtime() - offset_core;
    //const F64 offset_copy_original_order = GetWtime();
    copyForceOriginalOrder();
    //time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
    ptr_force_per_walk[0].freeMem();
    ptr_force_per_walk[1].freeMem();
    n_epi_per_walk[0].freeMem();
    n_epi_per_walk[1].freeMem();
    for (int i = 0; i < n_thread; i++) {
        n_epj_disp_per_thread[i].freeMem();
        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
            n_spj_disp_per_thread[i].freeMem();
        }
    }
    return ret;
}

#if defined(UNDER_CONSTRUCTION)
template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_dispatch, typename Tfunc_retrieve>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::sendJptcl(Tfunc_dispatch pfunc_dispatch) {
    // send all epj and spj
    Tepi** epi_dummy = nullptr;
    S32* n_epi_dummy = nullptr;
    S32** id_epj_dummy = nullptr;
    S32* n_epj_dummy = nullptr;
    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
        S32** id_spj_dummy = nullptr;
        S32* n_spj_dummy = nullptr;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy, (const S32**)id_epj_dummy, n_epj_dummy, (const S32**)id_spj_dummy, n_spj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(), spj_sorted_.getPointer(), spj_sorted_.size(), true);
    } else if (std::is_same_v<typename TSM::force_type, TagForceShort>) {
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy, (const S32**)id_epj_dummy, n_epj_dummy, epj_sorted_.getPointer(),
                       epj_sorted_.size(), true);
    } else {
        static_assert([] { return false; }());
    }
}

#endif  // UNDER_CONSTRUCTION)

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                        const bool clear_force,
                                                                                                        const S32 adr_ipg_first,
                                                                                                        const S32 n_ipg_this_rnd) {
    bool update_epi = false;
    CalcForceNoWalkShort<TreeCellGlb, IPGroup<typename TSM::ipg_type>, InteractionList, Tepi, Tepj, Tforce, Tfunc_ep_ep>(
        tc_glb_, ipg_, interaction_list_glb_, n_loc_tot_, epi_sorted_, epj_sorted_, force_sorted_, n_group_limit_, pfunc_ep_ep, adr_ipg_first,
        n_ipg_this_rnd, clear_force, update_epi, n_interaction_ep_ep_local_, time_profile_);
    bool copy_force = true;
    if (copy_force) {
        //const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        //time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
    }
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceNoWalkOrg(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                           const bool clear) {
#if defined(PS_DEBUG_PRINT_calcForceOnlyNoWalk)
    PARTICLE_SIMULATOR_PRINT_MSG_B("***** START calcForceOnlyNoWalk *****", 0);
#endif
    F64 time_offset = GetWtime();
    const S32 n_thread = Comm::getNumberOfThread();
    ReallocatableArray<Tepj> epj_for_force[n_thread];
    force_sorted_.resizeNoInitialize(n_loc_tot_);
    if (clear) {
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_loc_tot_; i++) {
            force_sorted_[i].clear();
        }
    }
    S64 n_interaction_ep_ep_tmp = 0;
    const S64 n_ipg = ipg_.size();
    if (n_ipg > 0) {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_interaction_ep_ep_tmp)
#endif
        for (S32 i = 0; i < n_ipg; i++) {
            const S32 ith = Comm::getThreadNum();
            const S32 n_epi = ipg_[i].n_ptcl_;
            const S32 adr_epi_head = ipg_[i].adr_ptcl_;
            const S32 n_epj = interaction_list_glb_.n_ep_[i];
            // if (Comm::getRank() == 0) {
            //     std::cerr << "i= " << i << " n_epj=" << n_epj << std::endl;
            // }
            const S32 adr_epj_head = interaction_list_glb_.n_disp_ep_[i];
            const S32 adr_epj_end = interaction_list_glb_.n_disp_ep_[i + 1];
            n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * n_epj;
            // epj_for_force_[ith].resizeNoInitialize(n_epj);
            epj_for_force[ith].resizeNoInitialize(n_epj);
            S32 n_cnt = 0;
            for (S32 j = adr_epj_head; j < adr_epj_end; j++, n_cnt++) {
                const S32 adr_epj = interaction_list_glb_.adr_ep_[j];
                // epj_for_force_[ith][n_cnt] = epj_sorted_[adr_epj];
                epj_for_force[ith][n_cnt] = epj_sorted_[adr_epj];
            }
            // pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head), n_epi, epj_for_force_[ith].getPointer(), n_epj,
            // force_sorted_.getPointer(adr_epi_head));
            pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head), n_epi, epj_for_force[ith].getPointer(), n_epj, force_sorted_.getPointer(adr_epi_head));
        }
    }
    n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
    copyForceOriginalOrder();
    // time_profile_.calc_force += GetWtime() - time_offset;
#if defined(PS_DEBUG_PRINT_calcForceOnlyNoWalk)
    PARTICLE_SIMULATOR_PRINT_MSG_B("***** END calcForceOnlyNoWalk *****", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tfunc_ep_sp>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                        Tfunc_ep_sp pfunc_ep_sp, const bool clear) {
#if defined(PS_DEBUG_PRINT_calcForceOnlyNoWalk)
    PARTICLE_SIMULATOR_PRINT_MSG_B("***** START calcForceOnlyNoWalk *****", 0);
    PS_PRINT_VARIABLE2(interaction_list_glb_.n_ep_.size(), interaction_list_glb_.n_sp_.size());
#endif
    F64 time_offset = GetWtime();
    const S32 n_thread = Comm::getNumberOfThread();
    ReallocatableArray<Tepj> epj_for_force[n_thread];
    ReallocatableArray<Tspj> spj_for_force[n_thread];
    force_sorted_.resizeNoInitialize(n_loc_tot_);
    if (clear) {
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_loc_tot_; i++) {
            force_sorted_[i].clear();
        }
    }
    S64 n_interaction_ep_ep_tmp = 0;
    S64 n_interaction_ep_sp_tmp = 0;
    const S64 n_ipg = ipg_.size();
    if (n_ipg > 0) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_interaction_ep_ep_tmp, n_interaction_ep_sp_tmp)
#endif
        for (S32 i = 0; i < n_ipg; i++) {
            const S32 ith = Comm::getThreadNum();
            const S32 n_epi = ipg_[i].n_ptcl_;
            const S32 adr_epi_head = ipg_[i].adr_ptcl_;
            const S32 n_epj = interaction_list_glb_.n_ep_[i];
            const S32 adr_epj_head = interaction_list_glb_.n_disp_ep_[i];
            const S32 adr_epj_end = interaction_list_glb_.n_disp_ep_[i + 1];
            const S32 n_spj = interaction_list_glb_.n_sp_[i];
            const S32 adr_spj_head = interaction_list_glb_.n_disp_sp_[i];
            const S32 adr_spj_end = interaction_list_glb_.n_disp_sp_[i + 1];
            n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * n_epj;
            n_interaction_ep_sp_tmp += ipg_[i].n_ptcl_ * n_spj;
            // epj_for_force_[ith].resizeNoInitialize(n_epj);
            epj_for_force[ith].resizeNoInitialize(n_epj);
            // spj_for_force_[ith].resizeNoInitialize(n_spj);
            spj_for_force[ith].resizeNoInitialize(n_spj);
            S32 n_ep_cnt = 0;
            for (S32 j = adr_epj_head; j < adr_epj_end; j++, n_ep_cnt++) {
                const S32 adr_epj = interaction_list_glb_.adr_ep_[j];
                // epj_for_force_[ith][n_ep_cnt] = epj_sorted_[adr_epj];
                epj_for_force[ith][n_ep_cnt] = epj_sorted_[adr_epj];
            }
            // pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head), n_epi, epj_for_force_[ith].getPointer(), n_epj,
            // force_sorted_.getPointer(adr_epi_head));
            pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head), n_epi, epj_for_force[ith].getPointer(), n_epj, force_sorted_.getPointer(adr_epi_head));
            S32 n_sp_cnt = 0;
#if 0
                // must use this combinig with AddMomentAsSp()
                for(S32 j=adr_spj_head; j<adr_spj_end; j++, n_sp_cnt++){
                    const S32 adr_spj = interaction_list_glb_.adr_sp_[j];
                    spj_for_force_[ith][n_sp_cnt] = spj_sorted_[adr_spj];
                }
#else
            for (S32 j = adr_spj_head; j < adr_spj_end; j++, n_sp_cnt++) {
                const S32 adr_spj = interaction_list_glb_.adr_sp_[j];
                if (adr_spj < n_let_sp_) {
                    // spj_for_force_[ith][n_sp_cnt] = spj_sorted_[adr_spj];
                    spj_for_force[ith][n_sp_cnt] = spj_sorted_[adr_spj];
                } else {
                    // spj_for_force_[ith][n_sp_cnt].copyFromMoment(tc_glb_[adr_spj - n_let_sp_].mom_);
                    spj_for_force[ith][n_sp_cnt].copyFromMoment(tc_glb_[adr_spj - n_let_sp_].mom_);
                }
            }
#endif
            // pfunc_ep_sp(epi_sorted_.getPointer(adr_epi_head), n_epi, spj_for_force_[ith].getPointer(), n_spj,
            // force_sorted_.getPointer(adr_epi_head));
            pfunc_ep_sp(epi_sorted_.getPointer(adr_epi_head), n_epi, spj_for_force[ith].getPointer(), n_spj, force_sorted_.getPointer(adr_epi_head));
        }
    }
    n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
    n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
    copyForceOriginalOrder();
    // time_profile_.calc_force += GetWtime() - time_offset;
#if defined(PS_DEBUG_PRINT_calcForceOnlyNoWalk)
    PARTICLE_SIMULATOR_PRINT_MSG_B("***** END calcForceOnlyNoWalk *****", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeInteractionList(const S32 adr_ipg, const bool clear) {
    makeInteractionListImpl(typename TSM::force_type(), adr_ipg, clear);
}

}  // namespace ParticleSimulator
