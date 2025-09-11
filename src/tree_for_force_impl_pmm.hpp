#pragma once
namespace ParticleSimulator {

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::correspondPM3Cells2TreeCellAndLocalDomain(const TREE_TYPE tree_type){
    if constexpr (std::is_same_v<typename TSM::search_type, TagSearchLongParticleMeshMultipole>) {
        //decomposeParticleMeshCellImpl(dinfo);
        //decomposeParticleMeshCellImpl(*p_domain_info_);
        decomposeParticleMeshCell();
        calcAdrTCFromParticleMeshCellIndex(tree_type);
        calcIPInfoFromParticleMeshCellIndex();
    }
}

// 
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::decomposeParticleMeshCell(){
#if defined(PS_DEBUG_PRINT_decomposeParticleMeshCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START decomposeParticleMeshCellImpl *****", 0);
#endif
    const S32 n_proc = comm_info_.getNumberOfProc();
    const S32 my_rank = comm_info_.getRank();
    // const S32 n_cell_tot = n_cell_pm3_.x * n_cell_pm3_.y * n_cell_pm3_.z;
    const S32 n_thread = Comm::getNumberOfThread();
    // Calculate a list of S32ort describing the size of PM cells
    // that overlap with the domain of each rank
    pos_pm_domain_.resize(n_proc);
#if defined(PS_DEBUG_PRINT_decomposeParticleMeshCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** IN decomposeParticleMeshCellImpl *****", 0);
#endif

    PS_OMP_PARALLEL_FOR
    for (S32 i = 0; i < n_proc; i++) {
        const F64ort pos_domain = p_domain_info_->getPosDomain(i);
        const S32vec low = morton_key_.getPMCellID(pos_domain.low_);
        // const S32vec low = morton_key_.getPMCellID(pos_domain.low_);
        S32vec high = morton_key_.getPMCellID(pos_domain.high_);
#if PS_DEBUG_PRINT_decomposeParticleMeshCellImpl >= 2
        if (Comm::getRank() == 0) {
            PS_PRINT_VARIABLE4(i, pos_domain, low, high);
        }
#endif
        //PS_PRINT_VARIABLE4(i, pos_domain, low, high);
        assert((0 <= low.x) && (0 <= low.y) && (0 <= low.z));
        assert((high.x <= n_cell_pm3_.x) && (high.y <= n_cell_pm3_.y) && (high.z <= n_cell_pm3_.z));
        if (high.x == n_cell_pm3_.x) high.x = n_cell_pm3_.x - 1;
        if (high.y == n_cell_pm3_.y) high.y = n_cell_pm3_.y - 1;
        if (high.z == n_cell_pm3_.z) high.z = n_cell_pm3_.z - 1;
        pos_pm_domain_[i].low_ = low;
        pos_pm_domain_[i].high_ = high;
        // In the implementation above, it could happen that `high` is exactly
        // located on the boundary of the PM cell, even though there is no local
        // particle that is exactly located on the boundary of the PM cell.
        // If this happens, `pos_pm_domain_` includes unnecessary PM cells.
        // You can get around this by using `inner_boundary_of_local_tree` of
        // all the MPI processes, but which requires MPI_Allgather.
    }

#if defined(PS_DEBUG_PRINT_decomposeParticleMeshCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** IN decomposeParticleMeshCellImpl *****", 0);
#if PS_DEBUG_PRINT_decomposeParticleMeshCellImpl >= 2
    if (Comm::getRank() == 0) {
        for (S32 i = 0; i < n_proc; i++) {
            PS_PRINT_VARIABLE2(i, pos_pm_domain_[i]);
        }
    }
#endif
#endif

    // Calculate a list of indices of PM cells that overlap with the local domain
    // Note that inner_boundary_of_trees
    pm_cell_idx_loc_.clear();
    const S32ort pos_my_pm_domain = pos_pm_domain_[my_rank];
    const S32vec low = pos_my_pm_domain.low_;
    const S32vec high = pos_my_pm_domain.high_;
    for (S32 z = low.z; z <= high.z; z++) {
        for (S32 y = low.y; y <= high.y; y++) {
            for (S32 x = low.x; x <= high.x; x++) {
                const S32 idx_1d = x + n_cell_pm3_.x * y + (n_cell_pm3_.x * n_cell_pm3_.y) * z;
                pm_cell_idx_loc_.push_back(idx_1d);
            }
        }
    }
#if defined(PS_DEBUG_PRINT_decomposeParticleMeshCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** IN decomposeParticleMeshCellImpl *****", 0);
    const auto n_pm_cell_idx_loc = pm_cell_idx_loc_.size();
    const auto n_pm_cell_idx_glb = comm_info_.getSum(n_pm_cell_idx_loc);
    const auto n_cell_tot = n_cell_pm3_.x * n_cell_pm3_.y * n_cell_pm3_.z;
    if (comm_info_.getRank() == 0) {
        std::cout << "n_pm_cell_idx_glb= " << n_pm_cell_idx_glb << " n_cell_tot= " << n_cell_tot << std::endl;
    }
#endif
    // assert(0);
    //  Make a list of indices of PM cells for which this process
    //  is responsible for the calculation of multipole moment.
    std::vector<std::vector<S32>> idx_tmp(n_thread);
    std::vector<std::unordered_map<S32, std::vector<S32>>> rank_list_tmp(n_thread);
#if defined(PS_DEBUG_PRINT_decomposeParticleMeshCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** IN decomposeParticleMeshCellImpl *****", 0);
#endif

    PS_OMP_PARALLEL_FOR
    for (decltype(pm_cell_idx_loc_.size()) i = 0; i < pm_cell_idx_loc_.size(); i++) {
        const S32 ith = Comm::getThreadNum();
        const S32 idx_1d = pm_cell_idx_loc_[i];
        // First, calculate three-dimensional cell index of the PM cell
        S32vec idx_3d;
        idx_3d.z = idx_1d / (n_cell_pm3_.x * n_cell_pm3_.y);
        idx_3d.y = (idx_1d - (n_cell_pm3_.x * n_cell_pm3_.y) * idx_3d.z) / n_cell_pm3_.x;
        idx_3d.x = idx_1d - (n_cell_pm3_.x * n_cell_pm3_.y) * idx_3d.z - n_cell_pm3_.x * idx_3d.y;
        // Next, calculate the box size of the PM cell
        const F64ort pos_cell = GetBoxOfParticleMeshCell(idx_3d, pos_unit_cell_.low_, width_pm_cell_);
        // Then, find MPI ranks whose domain intersects with `box`
        using intxn_info_t = std::pair<F64, S32>;  // volume & rank
        std::vector<intxn_info_t> info_list;
        for (S32 j = 0; j < n_proc; j++) {
            const F64ort& pos_domain = p_domain_info_->getPosDomain(j);
            if (pos_cell.overlapped(pos_domain)) {
                const F64ort pos_intxn = pos_cell.getIntersectionNoCheck(pos_domain);
                const F64 vol_intxn = pos_intxn.getVolume();
                intxn_info_t tmp;
                tmp.first = vol_intxn;
                tmp.second = j;
                info_list.push_back(tmp);
            }
        }
        // Finaly, find an MPI rank that has the largest volume of
        // the intersection between its domain and the cell.
        if (info_list.size() > 0) {
            std::sort(info_list.begin(), info_list.end(),
                      [](const intxn_info_t lhs, const intxn_info_t rhs) -> bool { return lhs.first > rhs.first; });
            const S32 rank_largest = info_list.front().second;
            if (rank_largest == my_rank) {
                idx_tmp[ith].push_back(idx_1d);
            }
            for (const auto& info : info_list) {
                rank_list_tmp[ith][idx_1d].push_back(info.second);
            }
        }
    }
#if defined(PS_DEBUG_PRINT_decomposeParticleMeshCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** IN decomposeParticleMeshCellImpl *****", 0);
#endif
    pm_cell_idx_mm_.clear();
    rank_list_from_pm_cell_idx_loc_.clear();
    for (S32 ith = 0; ith < n_thread; ith++) {
        for (decltype(idx_tmp[ith].size()) i = 0; i < idx_tmp[ith].size(); i++) {
            pm_cell_idx_mm_.push_back(idx_tmp[ith][i]);
        }
        for (const auto& kv : rank_list_tmp[ith]) {
            rank_list_from_pm_cell_idx_loc_[kv.first] = kv.second;
        }
    }

#if defined(PS_DEBUG_PRINT_decomposeParticleMeshCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** IN decomposeParticleMeshCellImpl *****", 0);
  #if PS_DEBUG_PRINT_decomposeParticleMeshCellImpl >= 3
    // Check
    {
        std::ostringstream ss;
        ss << "pm_cell_idx_loc_" << std::setfill('0') << std::setw(5) << my_rank << ".txt";
        const std::string file_name = ss.str();
        std::ofstream ofs;
        ofs.open(file_name.c_str(), std::ios::trunc);
        for (S32 i = 0; i < pm_cell_idx_loc_.size(); i++) {
            const S32 idx = pm_cell_idx_loc_[i];
            ofs << idx << std::endl;
        }
        ofs.close();
    }
    {
        std::ostringstream ss;
        ss << "pm_cell_index_mm_" << std::setfill('0') << std::setw(5) << my_rank << ".txt";
        const std::string file_name = ss.str();
        std::ofstream ofs;
        ofs.open(file_name.c_str(), std::ios::trunc);
        for (S32 i = 0; i < pm_cell_idx_mm_.size(); i++) {
            const S32 idx = pm_cell_idx_mm_[i];
            ofs << idx << std::endl;
        }
        ofs.close();
    }
  #endif
    const auto n_pm_cell_idx_mm_loc = pm_cell_idx_mm_.size();
    const auto n_rank_list_from_pm_cell_idx_loc = rank_list_from_pm_cell_idx_loc_.size();
    const auto n_pm_cell_idx_mm_glb = comm_info_.getSum(n_pm_cell_idx_mm_loc);
    const auto n_rank_list_from_pm_cell_idx_glb = comm_info_.getSum(n_rank_list_from_pm_cell_idx_loc);
    if (comm_info_.getRank() == 0) {
        PS_PRINT_VARIABLE3(n_pm_cell_idx_mm_glb, n_rank_list_from_pm_cell_idx_glb, n_cell_tot);
    }
    // exit(1);
    assert((S32)n_pm_cell_idx_mm_glb == n_cell_tot);
    PARTICLE_SIMULATOR_PRINT_MSG("END decomposeParticleMeshCellImpl", 0);
#endif
}


template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcAdrTCFromParticleMeshCellIndex(
    const TREE_TYPE tree_type) {

    if constexpr (std::is_same_v<typename TSM::search_type, TagSearchLongParticleMeshMultipole>) {  // OK
#if defined(PS_DEBUG_PRINT_calcAdrTCFromParticleMeshCellIndex)
        PARTICLE_SIMULATOR_PRINT_MSG("**** START calcAdrTCFromParticleMeshCellIndex ****", 0);
#endif
        if (tree_type == GLOBAL_TREE) {
#if defined(PS_DEBUG_PRINT_calcAdrTCFromParticleMeshCellIndex)
            PARTICLE_SIMULATOR_PRINT_MSG("**** tree_type == GLOBAL_TREE in calcAdrTCFromParticleMeshCellIndex ****", 0);
#endif
            // assert(is_global_tree_constructed_);
            adr_tc_glb_from_pm_cell_idx_.clear();
            CalcAdrTCFromParticleMeshCellIndex<TreeCellGlb, TreeParticle, Tepj, Tspj, MortonKey<typename TSM::mkey_type>>(
                adr_tc_level_partition_glb_, tc_glb_, tp_glb_, epj_sorted_,
                spj_sorted_,  // for the convenience of debugging
                n_cell_pm3_, lev_pm_cell_, morton_key_, adr_tc_glb_from_pm_cell_idx_);
        } else if (tree_type == LOCAL_TREE) {
#if defined(PS_DEBUG_PRINT_calcAdrTCFromParticleMeshCellIndex)
            PARTICLE_SIMULATOR_PRINT_MSG("**** tree_type == LOCAL_TREE in calcAdrTCFromParticleMeshCellIndex ****", 0);
#endif
            // assert(is_local_tree_constructed_);
            adr_tc_loc_from_pm_cell_idx_.clear();
            CalcAdrTCFromParticleMeshCellIndex<TreeCellLoc, TreeParticle, Tepj, Tspj, MortonKey<typename TSM::mkey_type>>(
                adr_tc_level_partition_loc_, tc_loc_, tp_loc_, epj_sorted_,
                spj_sorted_,  // just for compatibility; DO NOT ACCESS!
                n_cell_pm3_, lev_pm_cell_, morton_key_, adr_tc_loc_from_pm_cell_idx_);
        } else if (tree_type == LET_TREE) {
            // assert(is_LET_tree_constructed_);
            /*
            // under construction
            adr_tc_let_from_pm_cell_idx_.clear();
            CalcAdrTCFromParticleMeshCellIndex
            <TreeCellGlb, TreeParticle, Tepj, Tspj,
            MortonKey<typename TSM::mkey_type>>
            (adr_tc_level_partition_let_,
            tc_let_,
            tp_let_,
            epj_sorted_,
            spj_sorted_, // for the convenience of debugging
            n_cell_pm3_,
            lev_pm_cell_,
            morton_key_,
            adr_tc_let_from_pm_cell_idx_);
            */
        } else {
            PARTICLE_SIMULATOR_PRINT_ERROR("Necessary tree is not constructed!");
            Abort(-1);
        }
#if defined(PS_DEBUG_PRINT_calcAdrTCFromParticleMeshCellIndex)
        PARTICLE_SIMULATOR_PRINT_MSG("**** END calcAdrTCFromParticleMeshCellIndex ****", 0);
#endif
    }
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcIPInfoFromParticleMeshCellIndex() {
#if defined(PS_DEBUG_PRINT_calcIPInfoFromParticleMeshCellIndex)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START calcIPInfoFromParticleMeshCellIndex ****", 0);
#endif
    if constexpr (std::is_same_v<typename TSM::search_type, TagSearchLongParticleMeshMultipole>) {  // OK
        PS_COMPILE_TIME_MESSAGE("TagSearchLongParticleMeshMultipole@calcIPInfoFromParticleMeshCellIndex");
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        ip_info_from_pm_cell_idx_.clear();
        for (const S32 idx : pm_cell_idx_loc_) {
            auto itr = adr_tc_loc_from_pm_cell_idx_.find(idx);
            if (itr != adr_tc_loc_from_pm_cell_idx_.end()) {  // found
                // This PM cell contains local particles
                const S32 adr_tc = itr->second;
                const S32 n_ptcl = tc_loc_[adr_tc].n_ptcl_;
                const U32 adr_ptcl = tc_loc_[adr_tc].adr_ptcl_;
                assert(n_ptcl > 0);
                assert(GetMSB(adr_ptcl) == 0);  // just in case
                IParticleInfo<Tepi, Tforce> tmp;
                tmp.n_epi = n_ptcl;
                tmp.adr_epi_first = adr_ptcl;
                tmp.epi_first = epi_sorted_.getPointer(adr_ptcl);
                tmp.force_first = force_sorted_.getPointer(adr_ptcl);
                ip_info_from_pm_cell_idx_[idx] = tmp;
            } else {  // not found
                // This PM cell does not contain local particles
                IParticleInfo<Tepi, Tforce> tmp;
                tmp.n_epi = 0;
                tmp.adr_epi_first = -1;
                tmp.epi_first = nullptr;
                tmp.force_first = nullptr;
                ip_info_from_pm_cell_idx_[idx] = tmp;
            }
        }
    }

#if defined(PS_DEBUG_PRINT_calcIPInfoFromParticleMeshCellIndex)
#if PS_DEBUG_PRINT_calcIPInfoFromParticleMeshCellIndex >= 2
    PARTICLE_SIMULATOR_PRINT_MSG("IN calcIPInfoFromParticleMeshCellIndex, dump all elements of ip_info_from_pm_cell_idx_", 0);
    for (const auto& kv : ip_info_from_pm_cell_idx_) {
        std::cout << "pm cell index = " << kv.first << " ip information in the coresponding cell= " << kv.second << std::endl;
    }
#endif
    PARTICLE_SIMULATOR_PRINT_MSG("**** END calcIPInfoFromParticleMeshCellIndex ****", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename real_t, typename cplx_t, template<typename, typename> typename mm_t>
S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcMultipoleMomentOfParticleMeshCell(
    const S32 idx, const F64vec& center, const F64ort& pos_my_domain, mm_t<real_t, cplx_t>& mm, const S32 p_spj2mm) const {
#if defined(PS_DEBUG_PRINT_calcMultipoleMomentOfParticleMeshCell)
    PARTICLE_SIMULATOR_PRINT_MSG("START calcMultipoleMomentOfParticleMeshCell", 0);
#endif
    // This function calculates multipole moments for PM cell whose one-
    // dimensional cell index is `idx` assuming that `mm`, variable to
    // store multipole moments, is cleared by the caller.
    // The definitions of the other arguments are as follows.
    // `center` : the center of the PM cell, which is not calculated in
    //            TreeForForce class.
    // `pos_my_domain` : the part of the computational domain that this
    //                   process is responsible for.
    // `p_spj2mm` : the maximum order of the moment of SPJ that are
    //              converted to `mm`.
    // assert(is_global_tree_constructed_ || is_LET_tree_constructed_);
    S32 ret = -1;
    //const S32 my_rank = comm_info_.getRank();
    auto itr_outer = std::find(pm_cell_idx_mm_.begin(), pm_cell_idx_mm_.end(), idx);
    if (itr_outer != pm_cell_idx_mm_.end()) {
        F64 r_crit_sq;
        if (theta_ > 0.0) {
            r_crit_sq = (length_ * length_) / (theta_ * theta_);
        } else {
            r_crit_sq = -1.0;
        }
        S32vec idx_3d;
        idx_3d.z = idx / (n_cell_pm3_.x * n_cell_pm3_.y);
        idx_3d.y = (idx - (n_cell_pm3_.x * n_cell_pm3_.y) * idx_3d.z) / n_cell_pm3_.x;
        idx_3d.x = idx - (n_cell_pm3_.x * n_cell_pm3_.y) * idx_3d.z - n_cell_pm3_.x * idx_3d.y;
#if 1
        auto itr_inner = adr_tc_glb_from_pm_cell_idx_.find(idx);
        if (itr_inner != adr_tc_glb_from_pm_cell_idx_.end()) {
            const S32 adr_tc = itr_inner->second;
            const S32 adr_tree_sp_first = n_let_sp_;
            CalcMultipoleMomentOfParticleMeshCell<TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, real_t, cplx_t>(
                idx_3d, adr_tc, tc_glb_, tp_glb_, epj_sorted_, spj_sorted_, n_leaf_limit_, lev_leaf_limit_, adr_tree_sp_first, r_crit_sq,
                pos_unit_cell_.low_, i_cut_, lev_pm_cell_, width_pm_cell_, center, pos_my_domain, mm, p_spj2mm);
        } else {
            // In this case, no epj & spj exist in this PM cell.
            // mm = 0.
        }
        ret = 0;
#else
        // original
        if (is_global_tree_constructed_) {
            auto itr_inner = adr_tc_glb_from_pm_cell_idx_.find(idx);
            if (itr_inner != adr_tc_glb_from_pm_cell_idx_.end()) {
                const S32 adr_tc = itr_inner->second;
                const S32 adr_tree_sp_first = n_let_sp_;
                CalcMultipoleMomentOfParticleMeshCell<TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, real_t, cplx_t>(
                    idx_3d, adr_tc, tc_glb_, tp_glb_, epj_sorted_, spj_sorted_, n_leaf_limit_, lev_leaf_limit_, adr_tree_sp_first, r_crit_sq,
                    pos_unit_cell_.low_, icut_, lev_pm_cell_, width_pm_cell_, center, pos_my_domain, mm, p_spj2mm);
            } else {
                // In this case, no epj & spj exist in this PM cell.
                // mm = 0.
            }
            ret = 0;
        } else if (is_LET_tree_constructed_) {
            // Contribution from local tree
            {
                auto itr_inner = adr_tc_loc_from_pm_cell_idx_.find(idx);
                if (itr_inner != adr_tc_loc_from_pm_cell_idx_.end()) {
                    const S32 adr_tc = itr_inner->second;
                    const S32 adr_tree_sp_first = 0;  // No spj
                    CalcMultipoleMomentOfParticleMeshCell<TSM, TreeCellLoc, TreeParticle, Tepj, Tspj, real_t, cplx_t>(
                        idx_3d, adr_tc, tc_loc_, tp_loc_, epj_sorted_, spj_sorted_, n_leaf_limit_, lev_leaf_limit_, adr_tree_sp_first, r_crit_sq,
                        pos_unit_cell_.low_, icut_, lev_pm_cell_, width_pm_cell_, center, pos_my_domain, mm, p_spj2mm);
                } else {
                    // In this case, no local particle in this PM cell.
                }
            }
            // Contribution from LET tree
            {
                auto itr_inner = adr_tc_let_from_pm_cell_idx_.find(idx);
                if (itr_inner != adr_tc_let_from_pm_cell_idx_.end()) {
                    const S32 adr_tc = itr_inner->second;
                    const S32 adr_tree_sp_first = n_let_sp_;
                    CalcMultipoleMomentOfParticleMeshCell<TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, real_t, cplx_t>(
                        idx_3d, adr_tc, tc_let_, tp_let_, epj_sorted_, spj_sorted_, n_leaf_limit_, lev_leaf_limit_, adr_tree_sp_first, r_crit_sq,
                        pos_unit_cell_.low_, icut_, lev_pm_cell_, width_pm_cell_, center, pos_my_domain, mm, p_spj2mm);
                } else {
                    // In this case, no epj(LET) & spj(LET) exist in this PM cell.
                }
            }
            ret = 0;
        }
#endif
    }
#if defined(PS_DEBUG_PRINT_calcMultipoleMomentOfParticleMeshCell)
    PARTICLE_SIMULATOR_PRINT_MSG("END calcMultipoleMomentOfParticleMeshCell", 0);
#endif
    return ret;
}

}  // namespace ParticleSimulator
