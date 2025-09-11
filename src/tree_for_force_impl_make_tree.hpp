#pragma once
#include <tree_for_force_utils_make_tree.hpp>
namespace ParticleSimulator {

////////////////////////
// initialize morton key
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::initializeMortonKey(TagMortonKeyNormal) {
    PS_COMPILE_TIME_MESSAGE("initializeMortonKey(TagMortonKeyNormal)");
    morton_key_.initialize(length_ * 0.5, center_);
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::initializeMortonKey(TagMortonKeyMeshBased) {
    PS_COMPILE_TIME_MESSAGE("initializeMortonKey(TagMortonKeyMeshBased)");
    // morton_key_.initialize(pos_root_cell_.getHalfLength(), center_, pos_unit_cell_.low_, idx_unit_cell_, lev_pm_cell_, width_pm_cell_);
    morton_key_.initialize(pos_root_cell_.getHalfLength(), center_, pos_unit_cell_.low_, idx_pm3_unit_cell_, lev_pm_cell_, width_pm_cell_);

    //std::cout << "**************" << std::endl;
    //PS_PRINT_VARIABLE(pos_root_cell_);
    //PS_PRINT_VARIABLE4(center_, pos_unit_cell_, idx_pm3_unit_cell_, lev_pm_cell_);
    // PS_PRINT_VARIABLE4(center_, pos_unit_cell_, idx_unit_cell_, lev_pm_cell_);
    // PS_PRINT_VARIABLE2(((width_pm_cell_.x * 11) + pos_root_cell_.low_.x), ((width_pm_cell_.x * 12) + pos_root_cell_.low_.x));
    // auto idx_offset = (0.0 - pos_root_cell_.low_.x) / width_pm_cell_.x;
    // PS_PRINT_VARIABLE(idx_offset);
    //std::cout << "**************" << std::endl;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::initializeMortonKey(TagMortonKeyMeshBased2) {
    PS_COMPILE_TIME_MESSAGE("initializeMortonKey(TagMortonKeyMeshBased2)");

    //std::cout << "**************" << std::endl;
    //PS_PRINT_VARIABLE(pos_root_cell_);
    // PS_PRINT_VARIABLE4(center_, pos_unit_cell_, idx_unit_cell_, lev_pm_cell_);
    //PS_PRINT_VARIABLE4(center_, pos_unit_cell_, idx_pm3_unit_cell_, lev_pm_cell_);
    // PS_PRINT_VARIABLE2(((width_pm_cell_.x * 11) + pos_root_cell_.low_.x), ((width_pm_cell_.x * 12) + pos_root_cell_.low_.x));
    // auto idx_offset = (0.0 - pos_root_cell_.low_.x) / width_pm_cell_.x;
    // PS_PRINT_VARIABLE(idx_offset);
    //std::cout << "**************" << std::endl;

    // Comm::barrier();
    // exit(1);

    // morton_key_.initialize(pos_root_cell_.getHalfLength(), center_, pos_unit_cell_.low_, idx_unit_cell_, lev_pm_cell_, width_pm_cell_);
    morton_key_.initialize(pos_root_cell_.getHalfLength(), center_, pos_unit_cell_.low_, idx_pm3_unit_cell_, lev_pm_cell_, width_pm_cell_);
}
// initialize morton key
////////////////////////

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeLocalTreeOnly(DomainInfo& dinfo,
                                                                                                          const INTERACTION_LIST_MODE list_mode,
                                                                                                          const bool flag_serialize) {
    auto wtime0 = GetWtime();
#if defined(PS_DEBUG_PRINT_makeLocalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START makeLocalTreeOnly ***", 0);
    PARTICLE_SIMULATOR_PRINT_MSG("*** before setDomainInfo ***", 0);
#endif

    if (flag_serialize == true) {
        PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
        Abort(-1);
    }
    Check_INTERACTION_LIST_MODE(list_mode);
    setDomainInfo(dinfo);
    /*
    if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>){
        if( p_domain_info_->getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
            p_domain_info_->setPosRootDomain( inner_boundary_of_local_tree_.getExpandedBoxSize() );
            std::cout<<"inner_boundary_of_local_tree_= "<<inner_boundary_of_local_tree_<<std::endl;
            //exit(1);
            std::cout<<"p_domain_info_->getPosRootDomain()= "<<p_domain_info_->getPosRootDomain()<<std::endl;
        } else if (p_domain_info_->getBoundaryCondition() != BOUNDARY_CONDITION_PERIODIC_XYZ){
            assert(0);
        }
    }
    */

    if (list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE) {
#if defined(PS_DEBUG_PRINT_makeLocalTreeOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("*** before setRootCell ***", 0);
#endif
        auto t0 = GetWtime();
        setRootCell(dinfo);
        time_profile_.set_root_cell += GetWtime() - t0;
        /*
          if (param.set_root_cell_mode == SET_ROOT_CELL_MODE::NORMAL) {
          setRootCell(dinfo);
          } else if (param.set_root_cell_mode == SET_ROOT_CELL_MODE::USE_PRECALCULATED_OUTER_BOUNDARY) {
          //setRootCell(dinfo, param.outer_boundary);
          } else if (param.set_root_cell_mode == SET_ROOT_CELL_MODE::USE_PRECALCULATED_CENTER_AND_LENGTH) {
          //setRootCell(dinfo, param.length, param.center);
          }
        */

#if defined(PS_DEBUG_PRINT_makeLocalTreeOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("*** before mortonSortLocalTreeOnly ***", 0);
#endif
        t0 = GetWtime();
        mortonSortLocalTreeOnly(false);
        time_profile_.sort_particle_lt += GetWtime() - t0;
        /*
#if defined(DEBUG_PRINT_makeLocalTreeOnly)
        if(Comm::getRank()==0){
            std::cerr<<"A) epj_sorted_.size= "<<epj_sorted_.size()<<std::endl;
            for(auto i=0; i<epj_sorted_.size(); i++){
                std::cerr<<"i= "<<i<<" mass= "<<epj_sorted_[i].mass<<" pos= "<<epj_sorted_[i].pos<<std::endl;
            }
        }
#endif
        */
        epi_org_.freeMem();
#if defined(PS_DEBUG_PRINT_makeLocalTreeOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("*** before linkCellLocalTreeOnly ***", 0);
#endif
        t0 = GetWtime();
        linkCellLocalTreeOnly();
        time_profile_.link_cell_lt += GetWtime() - t0;

#if defined(PS_DEBUG_PRINT_makeLocalTreeOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("*** before calcMomentLocalTreeOnly ***", 0);
#endif
        t0 = GetWtime();
        calcMomentLocalTreeOnly();
        time_profile_.set_moment_lt += GetWtime() - t0;
    } else if (list_mode == REUSE_LIST) {
        auto t0 = GetWtime();
        mortonSortLocalTreeOnly(true);
        time_profile_.sort_particle_lt += GetWtime() - t0;
        epi_org_.freeMem();
        // DEBUG_PRINT_MAKING_TREE(comm_info_);
        t0 = GetWtime();
        if (typeid(typename TSM::force_type) == typeid(TagForceLong)) {
            // std::cerr<<"XXXXXXXX"<<std::endl;
            calcMomentLocalTreeOnly();
        } else if (typeid(typename TSM::force_type) == typeid(TagForceShort)) {
            // std::cerr<<"YYYYYYYYYY"<<std::endl;
        } else {
            assert(0);
        }
        time_profile_.set_moment_lt += GetWtime() - t0;
    } else {
        std::cout << "INTERACTION_LIST_MODE: " << list_mode << std::endl;
        Abort(-1);
    }
#if defined(PS_DEBUG_PRINT_makeLocalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** before decomposeParticleMeshCell ***", 0);
#endif
    // for PMMM
    auto t0 = GetWtime();
    correspondPM3Cells2TreeCellAndLocalDomain(LOCAL_TREE);
    time_profile_.map_pm_cell_lt += GetWtime() - t0;
    time_profile_.construct_lt += GetWtime() - wtime0;
    
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeAndExchangeLETOnly(DomainInfo& dinfo,
                                                                                                               const INTERACTION_LIST_MODE list_mode,
                                                                                                               const bool flag_serialize) {
#if defined(PS_DEBUG_PRINT_exchangeLETOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START exchangeLETOnly ***", 0);
#endif
    if (flag_serialize == true) {
        PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
        Abort(-1);
    }

    Check_INTERACTION_LIST_MODE(list_mode);
    auto flag_reuse = false;
    if (list_mode == REUSE_LIST) {
        flag_reuse = true;
    }
#if defined(PS_DEBUG_PRINT_exchangeLETOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** before exchangeLocalEssentialTree ***", 0);
#endif
    exchangeLocalEssentialTree(dinfo, flag_reuse);

#if defined(PS_DEBUG_PRINT_exchangeLETOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** END exchangeLETOnly ***", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeGlobalTreeOnly(DomainInfo& dinfo,
                                                                                                           const INTERACTION_LIST_MODE list_mode,
                                                                                                           const bool flag_serialize) {
    const auto wtime0 = GetWtime();
#if defined(PS_DEBUG_PRINT_makeGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START makeGlobalTreeOnly ***", 0);
#endif
    if (flag_serialize == true) {
        PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
        Abort(-1);
    }
    Check_INTERACTION_LIST_MODE(list_mode);
    if (list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE) {
#if defined(PS_DEBUG_PRINT_makeGlobalTreeOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("*** before setLocalEssentialTreeToGlobalTreeOnly ***", 0);
#endif
        auto t0 = GetWtime();
        setLocalEssentialTreeToGlobalTree(false);
        time_profile_.set_particle_gt += GetWtime() - t0;
        epj_send_.freeMem();
        spj_send_.freeMem();
#if defined(PS_DEBUG_PRINT_makeGlobalTreeOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("*** before mortonSortGlobalTreeOnly ***", 0);
#endif
        t0 = GetWtime();
        mortonSortGlobalTreeOnly(false);
        time_profile_.sort_particle_gt += GetWtime() - t0;
#if defined(PS_DEBUG_PRINT_makeGlobalTreeOnly)
        PS_PRINT_VARIABLE(tp_glb_.size());
#if PS_DEBUG_PRINT_makeGlobalTreeOnly >= 2
        if (comm_info_.getRank() == 0) {
            for (int i = 0; i < tp_glb_.size(); i++) {
                tp_glb_[i].dump(std::cout);
                auto adr = tp_glb_[i].adr_ptcl_;
                if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                    if (GetMSB(adr) == 0) {
                        PS_PRINT_VARIABLE4(i, adr, epj_sorted_[adr].mass, epj_sorted_[adr].pos);
                    } else {
                        adr = ClearMSB(adr);
                        PS_PRINT_VARIABLE4(i, adr, spj_sorted_[adr].mass, spj_sorted_[adr].pos);
                        // PS_PRINT_VARIABLE(spj_sorted_[adr].quad);
                    }
                } else {
                    PS_PRINT_VARIABLE4(i, adr, epj_sorted_[adr].mass, epj_sorted_[adr].pos);
                }
            }
        }
#endif
        for (int i = 1; i < tp_glb_.size(); i++) {
            assert(tp_glb_[i - 1].key_ <= tp_glb_[i].key_);
        }
        PARTICLE_SIMULATOR_PRINT_MSG("*** before linkCellGlobalTreeOnly ***", 0);
#endif
        t0 = GetWtime();
        linkCellGlobalTreeOnly();
        time_profile_.link_cell_gt += GetWtime() - t0;
        epj_org_.freeMem();
        spj_org_.freeMem();
#if defined(PS_DEBUG_PRINT_makeGlobalTreeOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("*** before calcMomentGlobalTreeOnly ***", 0);
        // exit(1);
#endif
        t0 = GetWtime();
        calcMomentGlobalTreeOnly();
        time_profile_.set_moment_gt += GetWtime() - t0;
#if defined(PS_DEBUG_PRINT_makeGlobalTreeOnly)
        if (Comm::getRank() == 0) {
            tc_glb_[0].dump(std::cout);
            if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
                PS_PRINT_VARIABLE(spj_sorted_.size());
#if PS_DEBUG_PRINT_makeGlobalTreeOnly >= 2
                for (auto i = 0; i < spj_sorted_.size(); i++) {
                    PS_PRINT_VARIABLE3(i, spj_sorted_[i].mass, spj_sorted_[i].pos);
                }
#endif
            }
        }
#endif
    } else if (list_mode == REUSE_LIST) {
        auto t0 = GetWtime();
        setLocalEssentialTreeToGlobalTree(true);
        time_profile_.set_particle_gt += GetWtime() - t0;
        epj_send_.freeMem();
        t0 = GetWtime();
        mortonSortGlobalTreeOnly(true);
        time_profile_.sort_particle_gt += GetWtime() - t0;
        epj_org_.freeMem();
        t0 = GetWtime();
        if (typeid(typename TSM::force_type) == typeid(TagForceLong)) {
            calcMomentGlobalTreeOnly();
        }
        time_profile_.set_moment_gt += GetWtime() - t0;
    } else {
        PARTICLE_SIMULATOR_PRINT_ERROR("INVALID INTERACTION_LIST_MODE.");
        std::cerr << "INTERACTION_LIST_MODE: " << list_mode << std::endl;
        Abort(-1);
    }

    auto t0 = GetWtime();

#if defined(PS_DEBUG_PRINT_makeGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** before calcAdrTCFromParticleMeshCellIndex(GLOBAL_TREE) ***", 0);
#endif
    calcAdrTCFromParticleMeshCellIndex(GLOBAL_TREE);

#if defined(PS_DEBUG_PRINT_makeGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** END makeGlobalTreeOnly ***", 0);
#endif
    time_profile_.map_pm_cell_gt += GetWtime() - t0;
    time_profile_.construct_gt += GetWtime() - wtime0;
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::mortonSortLocalTreeOnly(const bool reuse) {
#if defined(PS_DEBUG_PRINT_mortonSortLocalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("START mortonSortLocalTreeOnly", 0);
#endif
    //const F64 wtime_offset = GetWtime();
    //F64 wtime_make_key = 0.0;
    //F64 wtime_sort = 0.0;
    // F64 wtime_offset2 = 0.0;
    epi_sorted_.resizeNoInitialize(n_loc_tot_);
    epj_sorted_.resizeNoInitialize(n_loc_tot_);
    adr_org_from_adr_sorted_loc_.resizeNoInitialize(n_loc_tot_);
    tp_glb_.resizeNoInitialize(n_loc_tot_);
    if (!reuse) {
        //F64 wtime_offset2 = 0.0;
        initializeMortonKey(typename TSM::mkey_type());
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_loc_tot_; i++) {
            tp_glb_[i].setFromEP(epj_org_[i], i, morton_key_);
        }
        //wtime_make_key = GetWtime() - wtime_offset2;
#if defined(PS_DEBUG_PRINT_mortonSortLocalTreeOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("before mortons sort in Local Tree", 0);
#endif
        //wtime_offset2 = GetWtime();
#if defined(PARTICLE_SIMULATOR_USE_RADIX_SORT)
        ReallocatableArray<TreeParticle> tp_buf(n_loc_tot_, n_loc_tot_, MemoryAllocMode::Pool);
        rs_.lsdSort(tp_glb_.getPointer(), tp_buf.getPointer(), 0, n_loc_tot_ - 1);
        tp_buf.freeMem();
#elif defined(PARTICLE_SIMULATOR_USE_STD_SORT)
        std::sort(tp_glb_.getPointer(), tp_glb_.getPointer() + n_loc_tot_,
                  [](const TreeParticle& l, const TreeParticle& r) -> bool { return l.getKey() < r.getKey(); });
#elif defined(PARTICLE_SIMULATOR_USE_SAMPLE_SORT)
        SampleSortLib::samplesort_bodies(tp_glb_.getPointer(), n_loc_tot_, [](const TreeParticle& l) -> auto { return l.getKey(); });
#else
        MergeSortOmp(tp_glb_, 0, n_loc_tot_, [](const TreeParticle& l, const TreeParticle& r) -> bool { return l.getKey() < r.getKey(); });
#endif

        //wtime_sort = GetWtime() - wtime_offset2;

        ///////
        // PMMM
        tp_loc_.resizeNoInitialize(n_loc_tot_);
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_loc_tot_; i++) {
            tp_loc_[i] = tp_glb_[i];
        }
        // PMMM
        ///////

        PS_OMP_PARALLEL_FOR_S
        for (S32 i = 0; i < n_loc_tot_; i++) {
            const S32 adr = tp_glb_[i].adr_ptcl_;
            adr_org_from_adr_sorted_loc_[i] = adr;
            epi_sorted_[i] = epi_org_[adr];
            epj_sorted_[i] = epj_org_[adr];
        }
    } else {
        PS_OMP_PARALLEL_FOR_S
        for (S32 i = 0; i < n_loc_tot_; i++) {
            const S32 adr = adr_org_from_adr_sorted_loc_[i];
            epi_sorted_[i] = epi_org_[adr];
            epj_sorted_[i] = epj_org_[adr];
        }
    }  // end of if() no reuse
#ifdef SORTLIB_MEASURE_TIME
    std::cout << "After sort: " << GetWtime() - wtime_offset2 << "\n";
    auto time3 = GetWtime();
#endif

#ifdef SORTLIB_MEASURE_TIME
    std::cout << "After copy: " << GetWtime() - time3 << "\n";
#endif
#if defined(PS_DEBUG_PRINT_mortonSortLocalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("END mortonSortLocalTreeOnly", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::linkCellLocalTreeOnly() {
#if defined(PS_DEBUG_PRINT_linkCellLocalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START linkCellLocalTreeOnly ***", 0);
#endif
    LinkCell(tc_loc_, adr_tc_level_partition_loc_, tp_glb_.getPointer(), lev_max_loc_, n_loc_tot_, n_leaf_limit_, lev_leaf_limit_, morton_key_);
    /*
    S32 err = 0;
    tc_loc_[0].checkTree(epj_sorted_.getPointer(), tc_loc_.getPointer(), center_, length_*0.5,
                         n_leaf_limit_, lev_leaf_limit_, 1e-8, err);
    */

#if defined(PS_DEBUG_PRINT_linkCellLocalTreeOnly)
#if PS_DEBUG_PRINT_linkCellLocalTreeOnly >= 2
    if (Comm::getRank() == 0) PS_PRINT_VARIABLE2(tc_loc_.size(), lev_max_loc);
#endif
    PARTICLE_SIMULATOR_PRINT_MSG("*** END linkCellLocalTreeOnly ***", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMomentLocalTreeOnly() {
    if constexpr (std::is_same_v<typename TSM::calc_moment_local_tree_type, TagCalcMomLongEpjLt>) {
        PS_COMPILE_TIME_MESSAGE("TagCalcMomLongEpjLt@calcMomentLocalTreeOnly");
        CalcMomentLongLocalTree(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_,
                                lev_leaf_limit_, morton_key_);
    } else if constexpr (std::is_same_v<typename TSM::calc_moment_local_tree_type, TagCalcMomShortEpjLt>) {
        PS_COMPILE_TIME_MESSAGE("TagCalcMomShortEpjLt@calcMomentLocalTreeOnly");
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    } else if constexpr (std::is_same_v<typename TSM::calc_moment_local_tree_type, TagCalcMomShortEpiLt>) {
        PS_COMPILE_TIME_MESSAGE("TagCalcMomShortEpiLt@calcMomentLocalTreeOnly");
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epi_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
    // calcMomentLocalTreeOnlyImpl(typename TSM::calc_moment_local_tree_type());
#if PS_DEBUG_PRINT_calcMomentLocalTreeOnly >= 2
    if (Comm::getRank() == 0) {
        for (int i = 0; i < tc_loc_.size(); i++) {
            std::cout << "***************" << std::endl;
            std::cout << "i= " << i << std::endl;
            tc_loc_[i].dump(std::cout);
        }
    }
#endif
}

///////////////////////////////
/// morton sort global tree ///
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::mortonSortGlobalTreeOnly(const bool reuse) {
    //F64 time_offset = GetWtime();
    if (!map_id_to_epj_.empty()) {
        map_id_to_epj_.clear();
    }
    assert(map_id_to_epj_.empty());
    tp_glb_.resizeNoInitialize(n_glb_tot_);
    const S32 n_ep_tot = epj_org_.size();
    epj_sorted_.resizeNoInitialize(n_ep_tot);
    adr_org_from_adr_sorted_glb_.resizeNoInitialize(n_glb_tot_);
    //std::cout<<"n_glb_tot_= "<<n_glb_tot_<<std::endl;
    if (!reuse) {
        ReallocatableArray<TreeParticle> tp_buf(n_glb_tot_, n_glb_tot_, MemoryAllocMode::Default);
#if defined(PARTICLE_SIMULATOR_USE_RADIX_SORT)
        rs_.lsdSort(tp_glb_.getPointer(), tp_buf.getPointer(), 0, n_glb_tot_ - 1);
#elif defined(PARTICLE_SIMULATOR_USE_STD_SORT)
        std::sort(tp_glb_.getPointer(), tp_glb_.getPointer() + n_glb_tot_,
                  [](const TreeParticle& l, const TreeParticle& r) -> bool { return l.getKey() < r.getKey(); });
#elif defined(PARTICLE_SIMULATOR_USE_SAMPLE_SORT)
        SampleSortLib::samplesort_bodies(tp_glb_.getPointer(), n_glb_tot_, [](const TreeParticle& l) -> auto { return l.getKey(); });
#else
        const S32 n_add = n_glb_tot_ - n_loc_tot_;
        ReallocatableArray<TreeParticle> tp_add_buf(n_add, n_add, MemoryAllocMode::Default);
        ReallocatableArray<TreeParticle> tp_loc_buf(n_loc_tot_, n_loc_tot_, MemoryAllocMode::Default);
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_add; i++) {
            tp_add_buf[i] = tp_glb_[n_loc_tot_ + i];
        }
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_loc_tot_; i++) {
            tp_loc_buf[i] = tp_glb_[i];
        }
        // for(S32 i=1; i<n_loc_tot_; i++){ assert(tp_loc_buf[i-1].getKey() <= tp_loc_buf[i].getKey()); }
        MergeSortOmp(tp_add_buf, 0, n_add, [](const TreeParticle& l, const TreeParticle& r) -> bool { return l.getKey() < r.getKey(); });
        // for(S32 i=1; i<n_add; i++){ assert(tp_add_buf[i-1].getKey() <= tp_add_buf[i].getKey()); }
        MergeSortOmpImpl(tp_loc_buf.getPointer(0), tp_loc_buf.getPointer(n_loc_tot_), tp_add_buf.getPointer(0), tp_add_buf.getPointer(n_add),
                         tp_glb_.getPointer(0), [](const TreeParticle& l, const TreeParticle& r) -> bool { return l.getKey() < r.getKey(); });
        // for(S32 i=1; i<n_add; i++){ assert(tp_glb_[i-1].getKey() <= tp_glb_[i].getKey()); }
#endif
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_glb_tot_; i++) {
            adr_org_from_adr_sorted_glb_[i] = tp_glb_[i].adr_ptcl_;
        }
        tp_buf.freeMem();
    }

    // F64 timeaftersort = GetWtime();

    if (typeid(typename TSM::force_type) == typeid(TagForceLong)) {  // long mode
        if (!reuse) {
            const S32 n_sp_tot = spj_org_.size();
            assert(n_ep_tot + n_sp_tot == n_glb_tot_);
            spj_sorted_.resizeNoInitialize(n_sp_tot);
            adr_ep_org_from_adr_ep_sorted_glb_.resizeNoInitialize(n_ep_tot);
            adr_sp_org_from_adr_sp_sorted_glb_.resizeNoInitialize(n_sp_tot);
            //std::cout<<"n_ep_tot= "<<n_ep_tot<<std::endl;
            //std::cout<<"n_sp_tot= "<<n_sp_tot<<std::endl;
            const S32 n_thread = Comm::getNumberOfThread();
            U32 n_cnt_ep = 0;
            U32 n_cnt_sp = 0;
            if (n_thread < 3) {
                for (S32 i = 0; i < n_glb_tot_; i++) {
                    const U32 adr = adr_org_from_adr_sorted_glb_[i];
                    if (GetMSB(adr) == 0) {
                        epj_sorted_[n_cnt_ep] = epj_org_[adr];
                        tp_glb_[i].adr_ptcl_ = n_cnt_ep;
                        adr_ep_org_from_adr_ep_sorted_glb_[n_cnt_ep] = adr;
                        n_cnt_ep++;
                    } else {
                        spj_sorted_[n_cnt_sp] = spj_org_[ClearMSB(adr)];
                        tp_glb_[i].adr_ptcl_ = SetMSB(n_cnt_sp);
                        adr_sp_org_from_adr_sp_sorted_glb_[n_cnt_sp] = ClearMSB(adr);
                        n_cnt_sp++;
                    }
                }
            } else {
                U32 nwork[n_thread];
                U32 istart[n_thread];
                U32 iend[n_thread];
                U32 nep[n_thread];
                U32 nsp[n_thread];
                U32 epstart[n_thread];
                U32 spstart[n_thread];
                /*
                auto nwork0 = (n_glb_tot_ + n_thread - 1) / n_thread;
                for (auto i = 0; i < n_thread; i++) {
                    istart[i] = nwork0 * i;
                    nwork[i] = nwork0;
                    if (istart[i] + nwork[i] > n_glb_tot_) {
                        nwork[i] = n_glb_tot_ - istart[i];
                    }
                    nep[i] = 0;
                    nsp[i] = 0;
                }
                */

                PS_OMP_PARALLEL {
                    auto it = Comm::getThreadNum();
                    auto [is, ie] = CalcAdrToSplitData(it, n_thread, n_glb_tot_);
                    istart[it] = is;
                    iend[it] = ie;
                    nwork[it] = iend[it] - istart[it];
                    nep[it] = nsp[it] = 0;
                    for (U32 i = 0; i < nwork[it]; i++) {
                        auto adr = adr_org_from_adr_sorted_glb_[i + istart[it]];
                        if (GetMSB(adr) == 0) {
                            nep[it]++;
                        } else {
                            nsp[it]++;
                        }
                    }
                    PS_OMP_BARRIER
                    if (it == 0) {
                        epstart[0] = 0;
                        spstart[0] = 0;
                        for (auto it = 1; it < n_thread; it++) {
                            epstart[it] = epstart[it - 1] + nep[it - 1];
                            spstart[it] = spstart[it - 1] + nsp[it - 1];
                        }
                        auto itlast = n_thread - 1;
                        n_cnt_sp = spstart[itlast] + nsp[itlast];
                        n_cnt_ep = epstart[itlast] + nep[itlast];
                    }
                    PS_OMP_BARRIER
                    auto isp = 0;
                    auto iep = 0;
                    for (U32 i = 0; i < nwork[it]; i++) {
                        const U32 adr = adr_org_from_adr_sorted_glb_[i + istart[it]];
                        if (GetMSB(adr) == 0) {
                            epj_sorted_[epstart[it] + iep] = epj_org_[adr];
                            tp_glb_[i + istart[it]].adr_ptcl_ = epstart[it] + iep;
                            adr_ep_org_from_adr_ep_sorted_glb_[epstart[it] + iep] = adr;
                            iep++;
                        } else {
                            spj_sorted_[spstart[it] + isp] = spj_org_[ClearMSB(adr)];
                            tp_glb_[i + istart[it]].adr_ptcl_ = SetMSB(spstart[it] + isp);
                            adr_sp_org_from_adr_sp_sorted_glb_[spstart[it] + isp] = ClearMSB(adr);
                            isp++;
                        }
                    }
                }
            }
        } else {  // long & reuse
            const auto n_ep_tot = epj_org_.size();
            const auto n_sp_tot = spj_org_.size();
            PS_OMP_PARALLEL_FOR
            for (S32 i = 0; i < n_ep_tot; i++) {
                epj_sorted_[i] = epj_org_[adr_ep_org_from_adr_ep_sorted_glb_[i]];
            }
            PS_OMP_PARALLEL_FOR
            for (S32 i = 0; i < n_sp_tot; i++) {
                spj_sorted_[i] = spj_org_[adr_sp_org_from_adr_sp_sorted_glb_[i]];
            }
        }
    } else {  // short mode
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n_glb_tot_; i++) {
            const U32 adr = adr_org_from_adr_sorted_glb_[i];
            epj_sorted_[i] = epj_org_[adr];
            tp_glb_[i].adr_ptcl_ = i;
        }
    }

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT)
    PARTICLE_SIMULATOR_PRINT_LINE_INFO();
    std::cout << "tp_glb_.size()=" << tp_glb_.size() << std::endl;
    std::cout << "epj_sorted_.size()=" << epj_sorted_.size() << " spj_sorted_.size()=" << spj_sorted_.size() << std::endl;
#endif
}

/////////////////////////////
/// link cell global tree ///
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::linkCellGlobalTreeOnly() {
#if defined(PS_DEBUG_PRINT_linkCellGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START linkCellGlobalTreeOnly ****", 0);
    PARTICLE_SIMULATOR_PRINT_MSG("**** before LinkCell ****", 0);
#endif
#if 1
    // PMMM
    LinkCell(tc_glb_, adr_tc_level_partition_glb_, tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_, n_leaf_limit_, lev_leaf_limit_, morton_key_);
    // LinkCell(tc_glb_, adr_tc_level_partition_glb_, tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_, n_leaf_limit_, lev_leaf_limit_, morton_key_,
    // is_being_debugged_, comm_info_, n_thread);
#else
    LinkCell(tc_glb_, adr_tc_level_partition_glb_, tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_, n_leaf_limit_, morton_key_);
#endif
#if defined(PS_DEBUG_PRINT_linkCellGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("**** after LinkCell ****", 0);
    PS_PRINT_VARIABLE2(tc_glb_.size(), lev_max_glb_);
#endif
#if defined(PS_DEBUG_PRINT_linkCellGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("**** END linkCellGlobalTreeOnly ****", 0);
#endif
}

//////////////////////////
// CALC MOMENT GLOBAL TREE
#if 1
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMomentGlobalTreeOnly() {
#if defined(PS_DEBUG_PRINT_calcMomentGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START calcMomentGlobalTreeOnly ****", 0);
#endif
    //const auto time_offset = GetWtime();
    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
        if (exchange_let_mode_ == EXCHANGE_LET_A2A || exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) {
            CalcMomentLongGlobalTree(adr_tc_level_partition_glb_, tc_glb_.getPointer(), tp_glb_.getPointer(), epj_sorted_.getPointer(),
                                     spj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_, lev_leaf_limit_, morton_key_);
        } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_FAST) {
            CalcMomentLongGlobalTreeP2P(adr_tc_level_partition_glb_, tc_glb_.getPointer(), tp_glb_.getPointer(), epj_sorted_.getPointer(),
                                        spj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_, comm_table_, pos_root_cell_,
                                        adr_org_from_adr_sorted_glb_.getPointer(), morton_key_);
        }
    } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(), epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    } else {
        static_assert([] { return false; }());
    }
    if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_CUTOFF>) {
        SetOuterBoxGlobalTreeForLongCutoffTop(typename TSM::search_type(), tc_glb_.getPointer(),
                                              epj_sorted_.getPointer(),  // to get rcut
                                              n_leaf_limit_, length_ * 0.5, center_);
    }
#if defined(PS_DEBUG_PRINT_calcMomentGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("**** END calcMomentGlobalTreeOnly ****", 0);
#endif
}
#else
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMomentGlobalTreeOnly() {
#if defined(PS_DEBUG_PRINT_calcMomentGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START calcMomentGlobalTreeOnly ****", 0);
#endif
    const F64 time_offset = GetWtime();
    // calcMomentGlobalTreeOnlyImpl(typename TSM::force_type());
    calcMomentGlobalTreeOnlyImpl();
    if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_CUTOFF>) {
        SetOuterBoxGlobalTreeForLongCutoffTop(typename TSM::search_type(), tc_glb_.getPointer(),
                                              epj_sorted_.getPointer(),  // to get rcut
                                              n_leaf_limit_, length_ * 0.5, center_);
    }
    time_profile_.calc_moment_global_tree += GetWtime() - time_offset;
    time_profile_.make_global_tree_tot = time_profile_.calc_moment_global_tree + time_profile_.make_global_tree;
#if defined(PS_DEBUG_PRINT_calcMomentGlobalTreeOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("**** END calcMomentGlobalTreeOnly ****", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMomentGlobalTreeOnlyImpl() {
#if defined(PS_DEBUG_PRINT_calcMomentGlobalTreeOnlyImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START calcMomentGlobalTreeOnlyImpl *****", 0);
#endif
    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
        if (exchange_let_mode_ == EXCHANGE_LET_A2A || exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT) {
            /*
            if(Comm::getRank()==0){
                //std::cerr<<"adr_tc_level_partition_glb_.size()= "<<adr_tc_level_partition_glb_.size()<<std::endl;
                std::cerr<<" tc_glb_.size()= "<<tc_glb_.size()<<std::endl;
                std::cerr<<" tp_glb_.size()= "<<tp_glb_.size()<<std::endl;
                std::cerr<<" epj_sorted_.size()= "<<epj_sorted_.size()<<std::endl;
                std::cerr<<" spj_sorted_.size()= "<<spj_sorted_.size()<<std::endl;
                std::cerr<<" lev_max_glb_= "<<lev_max_glb_<<std::endl;
                std::cerr<<" n_leaf_limit_= "<<n_leaf_limit_<<std::endl;
                std::cerr<<" lev_leaf_limit_= "<<lev_leaf_limit_<<std::endl;
            }
            */
#if 1
            // PMM
            CalcMomentLongGlobalTree(adr_tc_level_partition_glb_, tc_glb_.getPointer(), tp_glb_.getPointer(), epj_sorted_.getPointer(),
                                     spj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_, lev_leaf_limit_, morton_key_);
            /*
            for(int i=0; i<4; i++){
                if(Comm::getRank()==i)   {
                    std::cerr<<"rank="<<i<<std::endl;
            CalcMomentLongGlobalTree
                (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
                 tp_glb_.getPointer(),     epj_sorted_.getPointer(),
                 spj_sorted_.getPointer(), lev_max_glb_,
                 n_leaf_limit_, lev_leaf_limit_, morton_key_);
            std::cerr<<"rank="<<i<<std::endl;
            tc_glb_[0].dump(std::cerr);
            }
                Comm::barrier();
                Comm::barrier();
            }
            */
#else
            // original
            CalcMomentLongGlobalTree(adr_tc_level_partition_glb_, tc_glb_.getPointer(), tp_glb_.getPointer(), epj_sorted_.getPointer(),
                                     spj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_, morton_key_);
#endif
        } else if (exchange_let_mode_ == EXCHANGE_LET_P2P_FAST) {
            CalcMomentLongGlobalTreeP2P(adr_tc_level_partition_glb_, tc_glb_.getPointer(), tp_glb_.getPointer(), epj_sorted_.getPointer(),
                                        spj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_, comm_table_, pos_root_cell_,
                                        adr_org_from_adr_sorted_glb_.getPointer(), morton_key_);
        }
    } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(), epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    } else {
        static_assert([] { return false; }());
    }
}
#endif

/*
template<class TSM, class Tforce, class Tepi, class Tepj,
         class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
calcMomentGlobalTreeOnlyImpl(TagForceLong){
    if(exchange_let_mode_ == EXCHANGE_LET_A2A
       || exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT){
#if 1
        // PMM
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_, lev_leaf_limit_, morton_key_);
#else
        // original
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_, morton_key_);
#endif
    }
    else if(exchange_let_mode_ == EXCHANGE_LET_P2P_FAST){
        CalcMomentLongGlobalTreeP2P
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_, comm_table_, pos_root_cell_,
             adr_org_from_adr_sorted_glb_.getPointer(), morton_key_);
    }
}
template<class TSM, class Tforce, class Tepi, class Tepj,
         class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
calcMomentGlobalTreeOnlyImpl(TagForceShort){
    CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
               epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
}
*/

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setLocalEssentialTreeToGlobalTree(const bool flag_reuse) {
    if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
        PS_COMPILE_TIME_MESSAGE("TagForceLong@setLocalEssentialTreeToGlobalTree");
        SetLocalEssentialTreeToGlobalTreeLong(epj_org_, spj_org_, n_loc_tot_, tp_glb_, morton_key_, flag_reuse);
        n_let_sp_ = spj_org_.size();
    } else {
        PS_COMPILE_TIME_MESSAGE("TagForceShort@setLocalEssentialTreeToGlobalTree");
        SetLocalEssentialTreeToGlobalTreeShort(epj_org_, n_loc_tot_, tp_glb_, morton_key_, flag_reuse);
    }
    this->n_glb_tot_ = tp_glb_.size();
}

}  // namespace ParticleSimulator
