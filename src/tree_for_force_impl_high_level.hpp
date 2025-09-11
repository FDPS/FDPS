#pragma once

#include <unordered_map>
#include <bitset>

#include <sort.hpp>
#include <tree.hpp>
#include <comm_table.hpp>
#include <interaction_list.hpp>
#include <tree_walk.hpp>
#include <tree_for_force_utils.hpp>
#include <tree_for_force_utils_force.hpp>

namespace ParticleSimulator {
//////////////////
// VARIADIC TEMPLATE VERSION

// NO MULTI WALK
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_ep_ep, typename Head, typename... Tails>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Head &&head,
                                                                                                                 Tails &&...tails) {
    calcForceAll(pfunc_ep_ep, std::forward<Head>(head), tails...);
    writeBack(head, tails...);
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_ep_ep, typename Head, typename... Tails>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAll(Tfunc_ep_ep pfunc_ep_ep, Head &&head,
                                                                                                     Tails &&...tails) {
    auto t0 = GetWtime();
    auto clear_force = true;
    auto list_mode = FDPS_DFLT_VAL_LIST_MODE;
    auto flag_serialize = false;
    setCalcForceParams(clear_force, list_mode, flag_serialize, tails...);

#if defined(PS_DEBUG_PRINT_calcForceAll)
    std::cout << "calcForceAll: clear_force=" << clear_force << ", list_mode=" << list_mode << ", flag_serialize=" << flag_serialize << std::endl;
#endif

    checkModeConsistency(*p_domain_info_);
#if defined(PS_DEBUG_PRINT_calcForceAll)
    std::cout << "after checkModeConsistency" << std::endl;
#endif
    setParticleLocalTreeVA(head, tails...);
    
#if defined(PS_DEBUG_PRINT_calcForceAll)
    std::cout << "after setParticleLocalTree" << std::endl;
#endif
    if constexpr (IsParticleSystem<Head>::value == true) {  // SHORT SEARCH MODE
        calcForceMakingTreeShort(pfunc_ep_ep, *p_domain_info_, clear_force, list_mode, flag_serialize);
    } else {  // LONG SEARCH MODE
        calcForceMakingTreeLong(pfunc_ep_ep, head, *p_domain_info_, clear_force, list_mode, flag_serialize);
    }
#if defined(PS_DEBUG_PRINT_calcForceAll)
    std::cout << "calcForceAll: finished calcForceMakingTree" << std::endl;
#endif
    //std::cout<<"time_profile_.calc_force_all__set_particle_local_tree: "<<time_profile_.calc_force_all__set_particle_local_tree<<std::endl;
    time_profile_.calc_force_all += GetWtime() - t0;
}

// MULTI WALK PTCL
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAllAndWriteBackMultiWalk(
    Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args) {
    S32 ret = calcForceAllMultiWalk(pfunc_dispatch, pfunc_retrieve, tag_max, args...);
    writeBack(args...);
    return ret;
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAllMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                                                                                             Tfunc_retrieve pfunc_retrieve,
                                                                                                             const S32 tag_max, Args &&...args) {
    auto t0 = GetWtime();
    S32 n_walk_limit = 10;
    bool clear_force = true;
    INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE;
    bool flag_serialize = false;
    setCalcForceParamsMultiWalk(n_walk_limit, clear_force, list_mode, flag_serialize, args...);
    checkModeConsistency(*p_domain_info_);
    setParticleLocalTreeVA(args...);
    S32 ret =
        calcForceMakingTreeMultiWalk(pfunc_dispatch, pfunc_retrieve, tag_max, *p_domain_info_, n_walk_limit, clear_force, list_mode, flag_serialize);
    time_profile_.calc_force_all += GetWtime() - t0;
    return ret;
}

// MULTI WALK INDEX
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAllAndWriteBackMultiWalkIndex(
    Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args) {
    S32 ret = calcForceAllMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, tag_max, args...);
    writeBack(args...);
    return ret;
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAllMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                                                                                                  Tfunc_retrieve pfunc_retrieve,
                                                                                                                  const S32 tag_max, Args &&...args) {
    auto t0 = GetWtime();
    S32 n_walk_limit = 10;
    bool clear_force = true;
    INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE;
    bool flag_serialize = false;
    setCalcForceParamsMultiWalk(n_walk_limit, clear_force, list_mode, flag_serialize, args...);
    checkModeConsistency(*p_domain_info_);
    setParticleLocalTreeVA(args...);
    S32 ret = calcForceMakingTreeMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, tag_max, *p_domain_info_, n_walk_limit, clear_force, list_mode,
                                                flag_serialize);
    time_profile_.calc_force_all += GetWtime() - t0;
    return ret;
}
// VARIADIC TEMPLATE VERSION
//////////////////

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tfunc_ep_sp>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                            Tfunc_ep_sp pfunc_ep_sp,
                                                                                                            DomainInfo &dinfo, const bool clear_force,
                                                                                                            const INTERACTION_LIST_MODE list_mode,
                                                                                                            const bool flag_serialize) {
    calcForceMakingTreeLong(pfunc_ep_ep, pfunc_ep_sp, dinfo, clear_force, list_mode, flag_serialize);
}
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                            DomainInfo &dinfo, const bool clear_force,
                                                                                                            const INTERACTION_LIST_MODE list_mode,
                                                                                                            const bool flag_serialize) {
    calcForceMakingTreeShort(pfunc_ep_ep, dinfo, clear_force, list_mode, flag_serialize);
}

//////////////////
// FOR LONG FORCE
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tfunc_ep_sp>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingTreeLong(
    Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, DomainInfo &dinfo, const bool clear_force, const INTERACTION_LIST_MODE list_mode,
    const bool flag_serialize) {
    PS_COMPILE_TIME_MESSAGE("calcForceMakingTreeLong");
#if defined(PS_DEBUG_PRINT_calcForceMakingTree)
    PARTICLE_SIMULATOR_PRINT_MSG_B("*** START calcForceMakingTree (Long) ***", 0);
#endif
    makeTreeAll(dinfo, list_mode, flag_serialize);
#if defined(PS_DEBUG_PRINT_calcForceMakingTree)
    PARTICLE_SIMULATOR_PRINT_MSG_B("*** after makeTreeAll (Long) ***", 0);
#endif
    AddMomentAsSp(typename TSM::force_type(), tc_glb_, spj_sorted_.size(), spj_sorted_);
#if defined(PS_DEBUG_PRINT_calcForceMakingTree)
    PARTICLE_SIMULATOR_PRINT_MSG_B("*** after AddMomentAsSp (Long) ***", 0);
#endif
    calcForceMakingGroupLong(pfunc_ep_ep, pfunc_ep_sp, clear_force, true, list_mode, GLOBAL_TREE);
#if defined(PS_DEBUG_CHECK_MEMORY_POOL)
    MemoryPool::checkEmpty();
#endif
#if defined(PS_DEBUG_PRINT_calcForceMakingTree)
    PARTICLE_SIMULATOR_PRINT_MSG_B("*** END calcForceMakingTree (Long) ***", 0);
#endif
}

//////////////////
// FOR SHORT FORCE
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingTreeShort(
    Tfunc_ep_ep pfunc_ep_ep, DomainInfo &dinfo, const bool clear_force, const INTERACTION_LIST_MODE list_mode, const bool flag_serialize) {
    PS_COMPILE_TIME_MESSAGE("calcForceMakingTreeShort");
    if (flag_serialize == true) {
        PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
        Abort(-1);
    }
    makeTreeAll(dinfo, list_mode, flag_serialize);
    calcForceMakingGroupShort(pfunc_ep_ep, clear_force, true, list_mode, GLOBAL_TREE);
#if defined(PS_DEBUG_CHECK_MEMORY_POOL)
    MemoryPool::checkEmpty();
#endif
}

///////////
// for multiwalk index
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_dispatch, class Tfunc_retrieve>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingTreeMultiWalkIndex(
    Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, DomainInfo &dinfo, const S32 n_walk_limit, const bool clear,
    const INTERACTION_LIST_MODE list_mode, const bool flag_serialize) {
#if defined(PS_DEBUG_PRINT_calcForceMakingTreeMultiWalkIndex)
    PARTICLE_SIMULATOR_PRINT_MSG_B("** START calcForceMakingTreeMultiWalkIndex **", 0);
#endif
    makeTreeAll(dinfo, list_mode, flag_serialize);
    AddMomentAsSp(typename TSM::force_type(), tc_glb_, spj_sorted_.size(), spj_sorted_);
#if defined(PS_DEBUG_PRINT_calcForceMakingTreeMultiWalk)
    PARTICLE_SIMULATOR_PRINT_MSG_B("** before makeIPGroup", 0);
#endif

#if 1
    S32 ret = 0;
    bool copy_force = true;
    // ret = calcForceMakingGroupMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, copy_force, clear);
    ret = calcForceMakingGroupMultiWalk<KERNEL_TYPE_INDEX>(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, copy_force, clear);
#else
    makeIPGroupOnly(list_mode);
    S32 ret = calcForceMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, clear);
#endif

#if defined(PS_DEBUG_CHECK_MEMORY_POOL)
    MemoryPool::checkEmpty();
#endif
#if defined(PS_DEBUG_PRINT_calcForceMakingTreeMultiWalkIndex)
    PARTICLE_SIMULATOR_PRINT_MSG_B("** END calcForceMakingTreeMultiWalkIndex", 0);
#endif
    return ret;
}

/////////////////////////////
// for multiwalk
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_dispatch, class Tfunc_retrieve>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceMakingTreeMultiWalk(
    Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, DomainInfo &dinfo, const S32 n_walk_limit, const bool clear,
    const INTERACTION_LIST_MODE list_mode, const bool flag_serialize) {
#if defined(PS_DEBUG_PRINT_calcForceMakingTreeMultiWalk)
    PARTICLE_SIMULATOR_PRINT_MSG_B("** START calcForceMakingTreeMultiWalk **", 0);
#endif
    makeTreeAll(dinfo, list_mode, flag_serialize);
    AddMomentAsSp(typename TSM::force_type(), tc_glb_, spj_sorted_.size(), spj_sorted_);
#if defined(PS_DEBUG_PRINT_calcForceMakingTreeMultiWalk)
    PARTICLE_SIMULATOR_PRINT_MSG_B("** before makeIPGroup", 0);
#endif

#if 1
    S32 ret = 0;
    bool copy_force = true;
    // ret = calcForceMakingGroupMultiWalkPtcl(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, copy_force, clear);
    ret = calcForceMakingGroupMultiWalk<KERNEL_TYPE_PTCL>(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, copy_force, clear);
#else
    makeIPGroupOnly(list_mode);
    S32 ret = 0;
    ret = calcForceMultiWalkPtcl(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, list_mode, clear);
#endif

#if defined(PS_DEBUG_CHECK_MEMORY_POOL)
    MemoryPool::checkEmpty();
#endif
#if defined(PS_DEBUG_PRINT_calcForceMakingTreeMultiWalk)
    PARTICLE_SIMULATOR_PRINT_MSG_B("** END calcForceMakingTreeMultiWalk", 0);
#endif
    return ret;
}

#if 1
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tfunc_ep_ep, typename Head, typename... Tails>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Head &&head, Tails &&...tails){
    //assert(0);
    // under construction
    auto clear_force = true;
    auto list_mode = FDPS_DFLT_VAL_LIST_MODE;
    auto flag_serialize = true;
    setCalcForceParams2(clear_force, list_mode, flag_serialize, std::forward<Tails>(tails)...);

    const auto copy_force = true; // copy force_org
    if constexpr (IsParticleSystem<Head>::value == true) {  // SHORT SEARCH MODE
        calcForceMakingGroupShort(pfunc_ep_ep, clear_force, copy_force, list_mode, GLOBAL_TREE);
    } else {  // LONG SEARCH MODE
        calcForceMakingGroupLong(pfunc_ep_ep, head, clear_force, copy_force, list_mode, GLOBAL_TREE);
    }
    writeBack(head, tails...);
}
#else
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                                                                                              Tfunc_ep_sp pfunc_ep_sp, Tpsys &psys,
                                                                                                              const bool clear,
                                                                                                              const bool flag_serialize) {
    if (flag_serialize == true) {
        PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
        Abort(-1);
    }
    // calcForce(pfunc_ep_ep, pfunc_ep_sp, clear);
    calcForceMakingGroupLong(pfunc_ep_ep, pfunc_ep_sp, clear, true, MAKE_LIST, GLOBAL_TREE);
    writeBack(psys);
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <class Tfunc_ep_ep, class Tpsys>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Tpsys &psys,
                                                                                                              const bool clear,
                                                                                                              const bool flag_serialize) {
    if (flag_serialize == true) {
        PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
        Abort(-1);
    }
    // calcForce(pfunc_ep_ep, clear);
    calcForceMakingGroupShort(pfunc_ep_ep, clear, true, MAKE_LIST, GLOBAL_TREE);
    writeBack(psys);
}
#endif

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeTreeAll(DomainInfo &dinfo,
                                                                                                    const INTERACTION_LIST_MODE list_mode,
                                                                                                    const bool flag_serialize) {
#if defined(PS_DEBUG_PRINT_makeTreeAll)
    PARTICLE_SIMULATOR_PRINT_MSG_B("**** START makeTreeAll ****", 0);
    PARTICLE_SIMULATOR_PRINT_MSG_B("- before makeLocalTreeOnly", 0);
#endif

    makeLocalTreeOnly(dinfo, list_mode, flag_serialize);

#if defined(PS_DEBUG_PRINT_makeTreeAll)
    PARTICLE_SIMULATOR_PRINT_MSG_B("- before exchangeLETOnly **", 0);
#endif
    makeAndExchangeLETOnly(dinfo, list_mode, flag_serialize);

#if defined(PS_DEBUG_PRINT_makeTreeAll)
    PARTICLE_SIMULATOR_PRINT_MSG_B("- before makeGlobalTreeOnly", 0);
#endif
    makeGlobalTreeOnly(dinfo, list_mode, flag_serialize);
#if defined(PS_DEBUG_PRINT_makeTreeAll)
    PARTICLE_SIMULATOR_PRINT_MSG_B("**** END makeTreeAll ****", 0);
#endif
}
}  // namespace ParticleSimulator
