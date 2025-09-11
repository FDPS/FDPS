namespace ParticleSimulator {
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeLocalTree(DomainInfo& dinfo) {
#if 1
    setDomainInfo(dinfo);
    setRootCell(dinfo);
    mortonSortLocalTreeOnly(false);
    epi_org_.freeMem();
    linkCellLocalTreeOnly();
    // calcMomentLocalTreeOnly();
#else
    setRootCell(dinfo);
    mortonSortLocalTreeOnly();
    linkCellLocalTreeOnly();
#endif
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeLocalTree(const F64 l, F64vec& c) {
#if 1
    setRootCell(l, c);
    mortonSortLocalTreeOnly(false);
    epi_org_.freeMem();
    linkCellLocalTreeOnly();
#else
    setRootCell(l, c);
    mortonSortLocalTreeOnly();
    linkCellLocalTreeOnly();
#endif
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeGlobalTree(DomainInfo& dinfo) {
#if 1
    calcMomentLocalTreeOnly();
    setLocalEssentialTreeToGlobalTree(false);
    epj_send_.freeMem();
    spj_send_.freeMem();
    mortonSortGlobalTreeOnly(false);
    linkCellGlobalTreeOnly();
    epj_org_.freeMem();
    spj_org_.freeMem();
    // calcMomentGlobalTreeOnly();
    // calcAdrTCFromParticleMeshCellIndex(GLOBAL_TREE);
#else
    calcMomentLocalTreeOnly();
    exchangeLocalEssentialTree(dinfo);
    setLocalEssentialTreeToGlobalTree();
    mortonSortGlobalTreeOnly();
    linkCellGlobalTreeOnly();
#endif
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMomentGlobalTree() {
#if 1
    calcMomentGlobalTreeOnly();
    calcAdrTCFromParticleMeshCellIndex(GLOBAL_TREE);
    AddMomentAsSp(typename TSM::force_type(), tc_glb_, spj_sorted_.size(), spj_sorted_);
#else
    calcMomentGlobalTreeOnly();
    makeIPGroup();
#endif
}

/*
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tfunc_ep_sp>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceNormal(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                        Tfunc_ep_sp pfunc_ep_sp, const bool clear,
                                                                                                        const INTERACTION_LIST_MODE list_mode) {
    S32 tag_max_dummy = 1;
    S32 n_walk_limit = 10000;
    CalcForceInitialize(force_sorted_, n_loc_tot_, tag_max_dummy, clear);
    makeListAndCalcForceImpl<KERNEL_TYPE_NORMAL>(pfunc_ep_ep, pfunc_ep_sp, tag_max_dummy, n_walk_limit, list_mode, clear);

}
*/
/*
// commented out 2025/02/22
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tfunc_ep_sp>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForce(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp,
                                                                                                  const bool clear) {
#if 1
    calcForceMakingGroupLong(pfunc_ep_ep, pfunc_ep_sp, clear, true, MAKE_LIST, GLOBAL_TREE);
#else
    // calcForceNormal(pfunc_ep_ep, pfunc_ep_sp, clear, MAKE_LIST);
    makeJPListOnly(0, ipg_.size(), MAKE_LIST, GLOBAL_TREE);
    calcForceNoWalkPMM(pfunc_ep_ep, pfunc_ep_sp, clear, false, GLOBAL_TREE);
#endif
    const S32 n_thread = Comm::getNumberOfThread();
    for (int i = 0; i < n_thread; i++) {
        epj_for_force_[i].freeMem();
        spj_for_force_[i].freeMem();
    }

}
    
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForce(Tfunc_ep_ep pfunc_ep_ep, const bool clear) {
#if 1
    calcForceMakingGroupShort(pfunc_ep_ep, clear, true, MAKE_LIST, GLOBAL_TREE);
#else    
    // class DummyClass {
    // } func_dummy;
    // calcForceNormal(pfunc_ep_ep, func_dummy, clear, MAKE_LIST);
    makeJPListOnly(0, ipg_.size(), MAKE_LIST, GLOBAL_TREE);
    calcForceNoWalk(pfunc_ep_ep, clear);

#endif
    const S32 n_thread = Comm::getNumberOfThread();
    for (int i = 0; i < n_thread; i++) {
        epj_for_force_[i].freeMem();
    }
}
// commented out 2025/02/22    
    */

/*
// commented out 2025/02/22    
//////////// Walk+Force, Kernel:Ptcl and Index, List:Index //////////////
template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename KERNEL_TYPE, typename Tfunc0, typename Tfunc1>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeListAndCalcForceImpl(
    Tfunc0 func0, Tfunc1 func1, const S32 tag_max, const S32 n_walk_limit, const INTERACTION_LIST_MODE list_mode, const bool clear) {
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
    PARTICLE_SIMULATOR_PRINT_MSG_B("***** START makeListAndCalcForceImpl", 0);
#endif
    const F64 offset_core = GetWtime();
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
            func0(0, 0, (const Tepi**)epi_dummy, n_epi_dummy, (const S32**)id_epj_dummy, n_epj_dummy, (const S32**)id_spj_dummy, n_spj_dummy,
                  epj_sorted_.getPointer(), epj_sorted_.size(), spj_sorted_.getPointer(), spj_sorted_.size(), true);
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
            func0(0, 0, (const Tepi**)epi_dummy, n_epi_dummy, (const S32**)id_epj_dummy, n_epj_dummy, epj_sorted_.getPointer(), epj_sorted_.size(),
                  true);
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
    ReallocatableArray<Tepi*> ptr_epi_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);
    ReallocatableArray<S32> n_epj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);
    ReallocatableArray<S32*> ptr_adr_epj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);
    // ReallocatableArray<S32>* n_epj_disp_per_thread;
    ReallocatableArray<S32> n_epj_disp_per_thread[n_thread];
    ReallocatableArray<S32> n_spj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);
    ReallocatableArray<S32*> ptr_adr_spj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);
    // ReallocatableArray<S32>* n_spj_disp_per_thread = nullptr;
    ReallocatableArray<S32> n_spj_disp_per_thread[n_thread];                                           // for long search case only
    ReallocatableArray<Tepj*> ptr_epj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);  // for ptcl kernel
    ReallocatableArray<Tspj*> ptr_spj_per_walk(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);  // for ptcl kernel
    ptr_force_per_walk[0].initialize(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);            // array of pointer *[n_walk]
    ptr_force_per_walk[1].initialize(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);            // array of pointer *[n_walk]
    n_epi_per_walk[0].initialize(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);
    n_epi_per_walk[1].initialize(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);
    // n_epj_disp_per_thread = new ReallocatableArray<S32>[n_thread];
    for (int i = 0; i < n_thread; i++) {
        n_epj_disp_per_thread[i].initialize(n_walk_limit + 1, n_walk_limit + 1, MemoryAllocMode::Default);
    }
    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
        // n_spj_disp_per_thread = new ReallocatableArray<S32>[n_thread];
        for (int i = 0; i < n_thread; i++) {
            n_spj_disp_per_thread[i].initialize(n_walk_limit + 1, n_walk_limit + 1, MemoryAllocMode::Default);
        }
    }
    const auto adr_tree_sp_first = n_let_sp_;
    bool first_loop = true;
    S32 n_walk_prev = 0;
    S64 n_interaction_ep_ep_tmp = 0;
    S64 n_interaction_ep_sp_tmp = 0;
    const auto len_peri = p_domain_info_->getLenRootDomain();
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
    PARTICLE_SIMULATOR_PRINT_MSG_B("- berofr walk loop", 0);
#endif
    if (n_ipg > 0) {
        for (int wg = 0; wg < n_loop_max; wg++) {
            const S32 n_walk = n_ipg / n_loop_max + (((n_ipg % n_loop_max) > wg) ? 1 : 0);
            const S32 walk_grp_head = (n_ipg / n_loop_max) * wg + std::min((n_ipg % n_loop_max), wg);
            const F64 offset_calc_force__core__walk_tree = GetWtime();
            const S32 lane_now = wg % 2;
            const S32 lane_old = (wg + 1) % 2;
            if (list_mode == REUSE_LIST) {
                if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL> || std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_NORMAL>) {
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
                        const S32 offset_ep = interaction_list_glb_.n_disp_ep_[iw_glb] - interaction_list_glb_.n_disp_ep_[walk_grp_head];
                        ptr_epj_per_walk[iw] = epj_for_force_[0].getPointer(offset_ep);
                        n_epj_per_walk[iw] = interaction_list_glb_.n_ep_[iw_glb];
                        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                            const S32 offset_sp = interaction_list_glb_.n_disp_sp_[iw_glb] - interaction_list_glb_.n_disp_sp_[walk_grp_head];
                            if (Comm::getRank() == 0) {
                                std::cout << "C) offset_sp= " << offset_sp << std::endl;
                            }
                            ptr_spj_per_walk[iw] = spj_for_force_[0].getPointer(offset_sp);
                            n_spj_per_walk[iw] = interaction_list_glb_.n_sp_[iw_glb];
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
            } else if (list_mode != REUSE_LIST) {
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                PARTICLE_SIMULATOR_PRINT_MSG_B("- list_mode != REUSE_LIST", 0);
#endif
                PS_OMP_PARALLEL {
                    const S32 ith = Comm::getThreadNum();
                    ReallocatableArray<S32> iwloc2iw(n_walk_limit, n_walk_limit, MemoryAllocMode::Default);
                    S32 n_walk_loc = 0;
                    S32 n_ep_cum = 0;
                    S32 n_sp_cum = 0;
                    n_epj_disp_per_thread[ith][0] = 0;
                    adr_epj_for_force_[ith].clearSize();
                    adr_ipg_for_force_[ith].clearSize();
                    adr_epj_for_force_[ith].reserve(n_walk);
                    adr_ipg_for_force_[ith].reserve(n_walk);
                    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                        n_spj_disp_per_thread[ith][0] = 0;
                        adr_spj_for_force_[ith].clearSize();
                        adr_spj_for_force_[ith].reserve(n_walk);
                    }
                    // #if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                    // #pragma omp for schedule(dynamic, 4) reduction(+: n_interaction_ep_ep_tmp) reduction(+: n_interaction_ep_sp_tmp)
                    // #endif
                    PS_OMP(omp for reduction(+: n_interaction_ep_ep_tmp), reduction(+: n_interaction_ep_sp_tmp))
                    for (int iw = 0; iw < n_walk; iw++) {
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                        PS_PRINT_VARIABLE3(ith, iw, n_walk);
#endif
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
                        S32 adr_tc = 0;
                        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                            MakeListUsingTreeRecursiveTop<TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl,
                                                          CALC_DISTANCE_TYPE>(tc_glb_, adr_tc, tp_glb_, epj_sorted_, adr_epj_for_force_[ith],
                                                                              spj_sorted_, adr_spj_for_force_[ith], target_box, n_leaf_limit_,
                                                                              adr_tree_sp_first, len_peri, theta_);
                        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
                            MakeListUsingTreeShortRecursiveTop<TSM, TreeCellGlb, TreeParticle, Tepj, TagChopLeafFalse, TagCopyInfoCloseNoSp,
                                                               CALC_DISTANCE_TYPE>(tc_glb_, adr_tc, tp_glb_, epj_sorted_, adr_epj_for_force_[ith],
                                                                                   target_box, n_leaf_limit_, len_peri);
                        }
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                        // PS_PRINT_VARIABLE4(iw, first_adr_ip, ptr_epi_per_walk[iw][0].id, ptr_epi_per_walk[iw][1].id);
                        PS_PRINT_VARIABLE4(iw, first_adr_ip, ptr_epi_per_walk[iw][0].pos.x, ptr_epi_per_walk[iw][1].pos.x);
#endif
                        n_epj_per_walk[iw] = adr_epj_for_force_[ith].size() - n_ep_cum;
                        n_ep_cum = adr_epj_for_force_[ith].size();
                        n_epj_disp_per_thread[ith][n_walk_loc + 1] = n_ep_cum;
                        n_interaction_ep_ep_tmp += ((S64)n_epj_per_walk[iw] * (S64)n_epi_per_walk[lane_now][iw]);
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                        PS_PRINT_VARIABLE2(iw, n_interaction_ep_ep_tmp);
#endif
                        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                            n_spj_per_walk[iw] = adr_spj_for_force_[ith].size() - n_sp_cum;
                            n_sp_cum = adr_spj_for_force_[ith].size();
                            n_spj_disp_per_thread[ith][n_walk_loc + 1] = n_sp_cum;
                            n_interaction_ep_sp_tmp += ((S64)n_spj_per_walk[iw] * (S64)n_epi_per_walk[lane_now][iw]);
                        }
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                        PS_PRINT_VARIABLE2(iw, n_interaction_ep_sp_tmp);
#endif
                        iwloc2iw[n_walk_loc] = iw;
                        n_walk_loc++;
                    }  // end of OMP for
                    if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL> || std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_NORMAL>) {
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                        PARTICLE_SIMULATOR_PRINT_MSG_B("- KERNEL_TYPE is KERNEL_TYPE_PTCL or KERNEL_TYPE_NORMAL", 0);
#endif
                        // PS_OMP_CRITICAL
                        epj_for_force_[ith].resizeNoInitialize(adr_epj_for_force_[ith].size());
                        for (S32 iwloc = 0; iwloc < n_walk_loc; iwloc++) {
                            S32 iw = iwloc2iw[iwloc];
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                            PS_PRINT_VARIABLE3(iw, iwloc, n_epj_per_walk[iw]);
#endif
                            for (int i = 0; i < n_epj_per_walk[iw]; i++) {
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                                PS_PRINT_VARIABLE(i);
#endif
                                S32 offset = n_epj_disp_per_thread[ith][iwloc];
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                                PS_PRINT_VARIABLE(offset);
#endif
                                S32 adr = adr_epj_for_force_[ith][offset + i];
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                                PS_PRINT_VARIABLE(adr);
                                PS_PRINT_VARIABLE(adr);
#endif
                                epj_for_force_[ith][offset + i] = epj_sorted_[adr];
                            }
                            ptr_epj_per_walk[iw] = epj_for_force_[ith].getPointer(n_epj_disp_per_thread[ith][iwloc]);
#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
                            PS_PRINT_VARIABLE(ptr_epj_per_walk[iw]);
#endif
                            if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                                // PS_OMP_CRITICAL
                                spj_for_force_[ith].resizeNoInitialize(adr_spj_for_force_[ith].size());
                                for (int i = 0; i < n_spj_per_walk[iw]; i++) {
                                    S32 offset = n_spj_disp_per_thread[ith][iwloc];
                                    S32 adr = adr_spj_for_force_[ith][offset + i];
                                    spj_for_force_[ith][offset + i] = spj_sorted_[adr];
                                }
                                ptr_spj_per_walk[iw] = spj_for_force_[ith].getPointer(n_spj_disp_per_thread[ith][iwloc]);
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

            time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
            if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_NORMAL>) {
                PS_OMP_PARALLEL_FOR
                for (int iw = 0; iw < n_walk; iw++) {
                    const S32 iw_glb = walk_grp_head + iw;
                    const S32 first_adr_ip = ipg_[iw_glb].adr_ptcl_;
                    // PS_PRINT_VARIABLE3(iw, iw_glb, first_adr_ip);
                    // PS_PRINT_VARIABLE2(n_epi_per_walk[lane_now][iw], n_epj_per_walk[iw]);
                    // PS_PRINT_VARIABLE(ptr_epi_per_walk[iw][0].id);
                    func0(ptr_epi_per_walk[iw], n_epi_per_walk[lane_now][iw], ptr_epj_per_walk[iw], n_epj_per_walk[iw],
                          force_sorted_.getPointer(first_adr_ip));
                    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                        // PS_PRINT_VARIABLE2(n_epi_per_walk[lane_now][iw], n_spj_per_walk[iw]);
                        func1(ptr_epi_per_walk[iw], n_epi_per_walk[lane_now][iw], ptr_spj_per_walk[iw], n_spj_per_walk[iw],
                              force_sorted_.getPointer(first_adr_ip));
                    }
                }
            } else {
                if (!first_loop) {
                    ret += func1(tag, n_walk_prev, n_epi_per_walk[lane_old].getPointer(), ptr_force_per_walk[lane_old].getPointer());
                }
                if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_INDEX>) {
                    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                        func0(0, n_walk, (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(),
                              (const S32**)ptr_adr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(),
                              (const S32**)ptr_adr_spj_per_walk.getPointer(), n_spj_per_walk.getPointer(), epj_sorted_.getPointer(),
                              epj_sorted_.size(), spj_sorted_.getPointer(), spj_sorted_.size(), false);
                    } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
                        func0(0, n_walk, (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(),
                              (const S32**)ptr_adr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(), epj_sorted_.getPointer(),
                              epj_sorted_.size(), false);
                    }
                } else if constexpr (std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_PTCL>) {
                    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                        ret += func0(tag, n_walk, (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(),
                                     (const Tepj**)ptr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(),
                                     (const Tspj**)ptr_spj_per_walk.getPointer(), n_spj_per_walk.getPointer());
                    } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
                        ret += func0(tag, n_walk, (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(),
                                     (const Tepj**)ptr_epj_per_walk.getPointer(),
                                     n_epj_per_walk.getPointer());  // new
                    }
                }
            }
            first_loop = false;
            n_walk_prev = n_walk;
        }  // end of walk group loop

#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
        PARTICLE_SIMULATOR_PRINT_MSG_B("- end of walk group loop", 0);
#endif

        if constexpr (!std::is_same_v<KERNEL_TYPE, KERNEL_TYPE_NORMAL>) {
            ret += func1(tag, n_walk_prev, n_epi_per_walk[(n_loop_max + 1) % 2].getPointer(), ptr_force_per_walk[(n_loop_max + 1) % 2].getPointer());
        }
        n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
    }  // if(n_ipg > 0)
    else {
        ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
        n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
    }
    time_profile_.calc_force__core += GetWtime() - offset_core;
    const F64 offset_copy_original_order = GetWtime();
    copyForceOriginalOrder();
    time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
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
    // delete[] n_epj_disp_per_thread;
    // delete[] n_spj_disp_per_thread;

#if defined(PS_DEBUG_PRINT_makeListAndCalcForceImpl)
    PARTICLE_SIMULATOR_PRINT_MSG_B("***** END makeListAndCalcForceImpl", 0);
#endif

    return ret;
}
// commented out 2025/02/22    
*/


}  // namespace ParticleSimulator