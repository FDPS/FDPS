#pragma once
namespace ParticleSimulator {
////////////////////
/// MAKE IPGROUP ///
// #if defined(RELEASE_V8)
#if 1
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeIPGroup() {
#if defined(PS_DEBUG_PRINT_makeIPGroup)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START makeIPGroup ****", 0);
#endif
    ipg_.clearSize();
    const F64 ipg_size_limit = coef_ipg_size_limit_ * inner_boundary_of_local_tree_.getFullLength().getMax();
#if defined(PS_DEBUG_PRINT_makeIPGroup)
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE5(ipg_size_limit, n_group_limit_, lev_group_limit_, n_leaf_limit_, lev_leaf_limit_);
    }
#endif
    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
        MakeIPGroupUseGLBTreeLong(ipg_, tc_loc_, tc_glb_, epi_sorted_, 0, 0, n_group_limit_, lev_group_limit_, n_leaf_limit_, lev_leaf_limit_,
                                  ipg_size_limit);
        const S32 n_ipg = ipg_.size();
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_ipg; i++) {
            const S32 n = ipg_[i].n_ptcl_;
            const S32 adr = ipg_[i].adr_ptcl_;
            ipg_[i].n_epi_ = n;
            ipg_[i].vertex_in_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
        }
    } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
        MakeIPGroupShort(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
    } else {
        static_assert([] { return false; }());
    }
    n_walk_local_ += ipg_.size();
#if defined(PS_DEBUG_PRINT_makeIPGroup)
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE(ipg_.size());
    }
#endif
#if defined(PS_DEBUG_PRINT_makeIPGroup)
    PARTICLE_SIMULATOR_PRINT_MSG("**** END makeIPGroup ****", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeIPGroupOnly(const INTERACTION_LIST_MODE list_mode) {
#if defined(PS_DEBUG_PRINT_makeIPGroupOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START makeIPGroupOnly ***", 0);
#endif
    // if (tree_type == GLOBAL_TREE) assert(is_global_tree_constructed_);
    // else if (tree_type == LOCAL_TREE) assert(is_local_tree_constructed_);
    // else if (tree_type == LET_TREE) assert(is_LET_tree_constructed_);
    if (list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE) {
        makeIPGroup();
        // if (!is_ipg_constructed_[N_POP]) {
        //     makeIPGroup(n_thread);
        //     is_ipg_constructed_[N_POP] = true;
        //     }
    }
#if defined(PS_DEBUG_PRINT_makeIPGroupOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** END makeIPGroupOnly ***", 0);
#endif
}

#elif 1
// PMMM
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeIPGroup() {
    const F64 time_offset = GetWtime();
    ipg_.clearSize();
    makeIPGroupImpl(typename TSM::force_type());
    n_walk_local_ += ipg_.size();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
    PARTICLE_SIMULATOR_PRINT_LINE_INFO();
    std::cout << "ipg_.size()=" << ipg_.size() << std::endl;
#endif
    time_profile_.calc_force += GetWtime() - time_offset;
    time_profile_.calc_force__make_ipgroup += GetWtime() - time_offset;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeIPGroupImpl(TagForceLong) {
    const F64 ipg_size_limit = coef_ipg_size_limit_ * inner_boundary_of_local_tree_.getFullLength().getMax();
    MakeIPGroupUseGLBTreeLong(ipg_, tc_loc_, tc_glb_, epi_sorted_, 0, 0, n_group_limit_, lev_group_limit_, n_leaf_limit_, lev_leaf_limit_,
                              ipg_size_limit);
    const S32 n_ipg = ipg_.size();
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_ipg; i++) {
        const S32 n = ipg_[i].n_ptcl_;
        const S32 adr = ipg_[i].adr_ptcl_;
        ipg_[i].n_epi_ = n;
        ipg_[i].vertex_in_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
    }
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeIPGroupImpl(TagForceShort) {
    const F64 ipg_size_limit = coef_ipg_size_limit_ * inner_boundary_of_local_tree_.getFullLength().getMax();
    MakeIPGroupUseGLBTreeShort(ipg_, tc_loc_, tc_glb_, epi_sorted_, 0, 0, n_group_limit_, lev_group_limit_, n_leaf_limit_, lev_leaf_limit_,
                               ipg_size_limit);
}

#else
// original
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeIPGroup() {
    const F64 time_offset = GetWtime();
    ipg_.clearSize();
    makeIPGroupImpl(typename TSM::force_type());
    n_walk_local_ += ipg_.size();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
    PARTICLE_SIMULATOR_PRINT_LINE_INFO();
    std::cout << "ipg_.size()=" << ipg_.size() << std::endl;
#endif
    time_profile_.calc_force += GetWtime() - time_offset;
    time_profile_.calc_force__make_ipgroup += GetWtime() - time_offset;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeIPGroupImpl(TagForceLong) {
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
    PARTICLE_SIMULATOR_PRINT_LINE_INFO();
    std::cout << "n_group_limit_=" << n_group_limit_ << std::endl;
#endif

#ifdef PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
    MakeIPGroupLongGLBTreeCellAsIPGBox(ipg_, tc_loc_, tc_glb_, epj_sorted_, 0, 0, n_group_limit_, n_leaf_limit_);
#else  // PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
#if 1
    MakeIPGroupUseGLBTreeLong(ipg_, tc_loc_, tc_glb_, epi_sorted_, 0, 0, n_group_limit_, n_leaf_limit_);  // NEW
#else
    MakeIPGroupLong(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
#endif

    const S32 n_ipg = ipg_.size();
    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_ipg; i++) {
        const S32 n = ipg_[i].n_ptcl_;
        const S32 adr = ipg_[i].adr_ptcl_;
        ipg_[i].vertex_in_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
    }
#endif  // PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeIPGroupImpl(TagForceShort) {
    MakeIPGroupShort(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
}
#endif

}  // namespace ParticleSimulator
