namespace ParticleSimulator {

template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeInteractionListIndexLong(const TREE_TYPE tree_type,
                                                                                                                     const S32 adr_ipg_first,
                                                                                                                     const S32 n_ipg_this_rnd) {
#if defined(PS_DEBUG_PRINT_makeInteractionListIndexLong2)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START makeInteractionListIndexLong2 ****", 0);
#endif
    //const F64 time_offset = GetWtime();
    F64vec len_peri = p_domain_info_->getPosRootDomain().getFullLength();
    // std::cerr<<"XXXXXX"<<std::endl;
    //  [Note by D.Namekata]
    //     To make the code more compact, all of the definitions of structures used as
    //     "tag" such as TagChop* should be moved to ps_defs.hpp and we should declare
    //     aliases in structures SEARCH_MODE_* necessary to control how we walk tree
    //     when making interaction list and exchanging LETs. If this is realized, we
    //     can eliminate the if branch below.
    // if (typeid(TSM) != typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE)) {
    if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
#if defined(PS_DEBUG_PRINT_makeInteractionListIndexLong2)
        PARTICLE_SIMULATOR_PRINT_MSG("**** TSM == SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, in makeInteractionListIndexLong2 ****", 0);
#endif
        if (tree_type == GLOBAL_TREE) {
            // std::cerr << "CHECK A" << std::endl;
            MakeInteractionListIndexLong<TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, IPGroup<typename TSM::ipg_type>,
                                         MortonKey<typename TSM::mkey_type>, TagChopLeafFalse, TagChopNonleafFalse, TagCopyInfoCloseWithTpAdrptcl,
                                         InteractionList, CALC_DISTANCE_TYPE>(tc_glb_, 0, tp_glb_, n_let_sp_, epj_sorted_, spj_sorted_, ipg_,
                                                                              adr_ipg_first, n_ipg_this_rnd, n_leaf_limit_, lev_leaf_limit_, theta_,
                                                                              i_cut_, 0, morton_key_, interaction_list_glb_, len_peri);
            // std::cerr << "CHECK B" << std::endl;
        } else if (tree_type == LOCAL_TREE) {
            /*
              MakeInteractionListIndexLong
              <TSM, TreeCellLoc, TreeParticle, Tepj, Tspj,
              IPGroup<typename TSM::ipg_type>,
              MortonKey<typename TSM::mkey_type>,
              TagWalkModeNormal, TagChopLeafFalse, TagCopyInfoCloseWithTpAdrptcl, TagChopNonleafFalse,
              InteractionList>
              (tc_loc_,
              0,
              tp_loc_,
              0,
              epj_sorted_,
              spj_sorted_,
              ipg_[pop_id],
              adr_ipg_first,
              n_ipg_this_rnd,
              n_leaf_limit_,
              lev_leaf_limit_,
              theta_,
              icut_,
              0,
              morton_key_,
              interaction_list_loc_[pop_id],
              n_thread);
            */
        } else if (tree_type == LET_TREE) {
            /*
              MakeInteractionListIndexLong
              <TSM, TreeCellGlb, TreeParticle, Tepj, Tspj,
              IPGroup<typename TSM::ipg_type>,
              MortonKey<typename TSM::mkey_type>,
              TagWalkModeNormal, TagChopLeafFalse, TagCopyInfoCloseWithTpAdrptcl, TagChopNonleafFalse,
              InteractionList>
              (tc_let_,
              0,
              tp_let_,
              n_let_sp_,
              epj_sorted_,
              spj_sorted_,
              ipg_[pop_id],
              adr_ipg_first,
              n_ipg_this_rnd,
              n_leaf_limit_,
              lev_leaf_limit_,
              theta_,
              icut_,
              0,
              morton_key_,
              interaction_list_let_[pop_id],
              n_thread);
            */
        }
    } else {
#if defined(PS_DEBUG_PRINT_makeInteractionListIndexLong2)
        PARTICLE_SIMULATOR_PRINT_MSG("**** TSM != SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, in makeInteractionListIndexLong2 ****", 0);
#endif
        if (tree_type == GLOBAL_TREE) {
            // std::cerr << "CHECK C" << std::endl;
            MakeInteractionListIndexLong<TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, IPGroup<typename TSM::ipg_type>,
                                         MortonKey<typename TSM::mkey_type>, TagChopLeafTrue, TagChopNonleafTrue, TagCopyInfoCloseWithTpAdrptcl,
                                         InteractionList, CALC_DISTANCE_TYPE>(tc_glb_, 0, tp_glb_, n_let_sp_, epj_sorted_, spj_sorted_, ipg_,
                                                                              adr_ipg_first, n_ipg_this_rnd, n_leaf_limit_, lev_leaf_limit_, theta_,
                                                                              i_cut_, lev_pm_cell_, morton_key_, interaction_list_glb_, len_peri);
            // std::cerr << "CHECK D" << std::endl;
        } else if (tree_type == LOCAL_TREE) {
            /*
              MakeInteractionListIndexLong
              <TSM, TreeCellLoc, TreeParticle, Tepj, Tspj,
              IPGroup<typename TSM::ipg_type>,
              MortonKey<typename TSM::mkey_type>,
              TagWalkModeNormal, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl, TagChopNonleafTrue,
              InteractionList>
              (tc_loc_,
              0,
              tp_loc_,
              0,
              epj_sorted_,
              spj_sorted_,
              ipg_[pop_id],
              adr_ipg_first,
              n_ipg_this_rnd,
              n_leaf_limit_,
              lev_leaf_limit_,
              theta_,
              icut_,
              lev_pm_cell_,
              morton_key_,
              interaction_list_loc_[pop_id],
              n_thread);
            */
        } else if (tree_type == LET_TREE) {
            /*
              MakeInteractionListIndexLong
              <TSM, TreeCellGlb, TreeParticle, Tepj, Tspj,
              IPGroup<typename TSM::ipg_type>,
              MortonKey<typename TSM::mkey_type>,
              TagWalkModeNormal, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl, TagChopNonleafTrue,
              InteractionList>
              (tc_let_,
              0,
              tp_let_,
              n_let_sp_,
              epj_sorted_,
              spj_sorted_,
              ipg_[pop_id],
              adr_ipg_first,
              n_ipg_this_rnd,
              n_leaf_limit_,
              lev_leaf_limit_,
              theta_,
              icut_,
              lev_pm_cell_,
              morton_key_,
              interaction_list_let_[pop_id],
              n_thread);
            */
        }
    }
    //const F64 time_elapsed = GetWtime() - time_offset;
    // time_profile_.calc_force += time_elapsed;
    //  time_profile_.calc_force__make_list += time_elapsed;
    //  n_walk_local_ += ipg_[pop_id].size();
    n_walk_local_ += ipg_.size();
#if defined(PS_DEBUG_PRINT_makeInteractionListIndexLong2)
    PARTICLE_SIMULATOR_PRINT_MSG("**** END makeInteractionListIndexLong2 ****", 0);
#endif
}

template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeInteractionListIndexShort(const TREE_TYPE tree_type,
                                                                                                                      const S32 adr_ipg_first,
                                                                                                                      const S32 n_ipg_this_rnd) {
#if defined(PS_DEBUG_PRINT_makeInteractionListIndexShort2)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START makeInteractionListIndexShort2 ****", 0);
#endif
    //const F64 time_offset = GetWtime();
    if (tree_type == GLOBAL_TREE) {
        MakeInteractionListIndexShort<TSM, TreeCellGlb, TreeParticle, Tepj, IPGroup<typename TSM::ipg_type>, MortonKey<typename TSM::mkey_type>,
                                      TagChopLeafFalse, TagChopNonleafFalse, TagCopyInfoCloseWithTpAdrptcl, InteractionList, CALC_DISTANCE_TYPE>(
            tc_glb_, 0, tp_glb_, epj_sorted_, ipg_, adr_ipg_first, n_ipg_this_rnd, n_leaf_limit_, lev_leaf_limit_, theta_, morton_key_,
            interaction_list_glb_);
    } else if (tree_type == LOCAL_TREE) {
        /*
          MakeInteractionListIndexLong
          <TSM, TreeCellLoc, TreeParticle, Tepj, Tspj,
          IPGroup<typename TSM::ipg_type>,
          MortonKey<typename TSM::mkey_type>,
          TagWalkModeNormal, TagChopLeafFalse, TagCopyInfoCloseWithTpAdrptcl, TagChopNonleafFalse,
          InteractionList>
          (tc_loc_,
          0,
          tp_loc_,
          0,
          epj_sorted_,
          spj_sorted_,
          ipg_[pop_id],
          adr_ipg_first,
          n_ipg_this_rnd,
          n_leaf_limit_,
          lev_leaf_limit_,
          theta_,
          icut_,
          0,
          morton_key_,
          interaction_list_loc_[pop_id],
          n_thread);
        */
    } else if (tree_type == LET_TREE) {
        /*
          MakeInteractionListIndexLong
          <TSM, TreeCellGlb, TreeParticle, Tepj, Tspj,
          IPGroup<typename TSM::ipg_type>,
          MortonKey<typename TSM::mkey_type>,
          TagWalkModeNormal, TagChopLeafFalse, TagCopyInfoCloseWithTpAdrptcl, TagChopNonleafFalse,
          InteractionList>
          (tc_let_,
          0,
          tp_let_,
          n_let_sp_,
          epj_sorted_,
          spj_sorted_,
          ipg_[pop_id],
          adr_ipg_first,
          n_ipg_this_rnd,
          n_leaf_limit_,
          lev_leaf_limit_,
          theta_,
          icut_,
          0,
          morton_key_,
          interaction_list_let_[pop_id],
          n_thread);
        */
    }
    //const F64 time_elapsed = GetWtime() - time_offset;
    // time_profile_.calc_force += time_elapsed;
    n_walk_local_ += ipg_.size();
#if defined(PS_DEBUG_PRINT_makeInteractionListIndexShort2)
    PARTICLE_SIMULATOR_PRINT_MSG("**** END makeInteractionListIndexShort2 ****", 0);
#endif
}

#if 0
template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeJPListOnly20241019(const INTERACTION_LIST_MODE list_mode,
                                                                                                               const TREE_TYPE tree_type) {
#if defined(PS_DEBUG_PRINT_makeJPListOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START makeJPListOnly ***", 0);
#endif
    if (list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE) {
#if defined(PS_DEBUG_PRINT_makeJPListOnly)
        PARTICLE_SIMULATOR_PRINT_MSG("***  list mode is MAKE_LIST or MAKE_LIST_FOR_REUSE in makeJPListOnly ***", 0);
#endif
        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
#if defined(PS_DEBUG_PRINT_makeJPListOnly)
            PARTICLE_SIMULATOR_PRINT_MSG("***  force_type is TagForceLong in makeJPListOnly ***", 0);
#endif
            makeInteractionListIndexLong(tree_type, 0, ipg_.size());
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
#if defined(PS_DEBUG_PRINT_makeJPListOnly)
            PARTICLE_SIMULATOR_PRINT_MSG("***  force_type is TagForceShort in makeJPListOnly ***", 0);
#endif
            makeInteractionListIndexShort(tree_type, 0, ipg_.size());
        } else {
            static_assert([] { return false; }());
        }
    }
#if defined(PS_DEBUG_PRINT_makeJPListOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** END makeJPListOnly ***", 0);
#endif
}
#endif

template <typename TSM, typename Tforce, typename Tepi, typename Tepj, typename Tmomloc, typename Tmomglb, typename Tspj,
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::makeJPListOnly(const S32 adr_ipg_first,
                                                                                                       const S32 n_ipg_this_rnd,
                                                                                                       const INTERACTION_LIST_MODE list_mode,
                                                                                                       const TREE_TYPE tree_type) {
#if defined(PS_DEBUG_PRINT_makeJPListOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** START makeJPListOnly ***", 0);
#endif
    if (list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE) {
        if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
            makeInteractionListIndexLong(tree_type, adr_ipg_first, n_ipg_this_rnd);
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
            makeInteractionListIndexShort(tree_type, adr_ipg_first, n_ipg_this_rnd);
        } else {
            static_assert([] { return false; }());
        }
    }
#if defined(PS_DEBUG_PRINT_makeJPListOnly)
    PARTICLE_SIMULATOR_PRINT_MSG("*** END makeJPListOnly ***", 0);
#endif
}

}  // namespace ParticleSimulator
