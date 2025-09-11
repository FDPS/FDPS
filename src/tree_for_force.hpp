#pragma once

// #if __cplusplus <= 199711L
// #include<map>
// #else
// #include<unordered_map>
// #endif

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

///////////////////////////////////////
/// TREE FOR FORCE CLASS DEFINITION ///
template <class TSM,      // search mode
          class Tforce,   // USER def
          class Tepi,     // USER def
          class Tepj,     // USER def
          class Tmomloc,  // PS or USER def
          class Tmomglb,  // PS or USER def
          class Tspj,     // PS or USER def
          enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE = CALC_DISTANCE_TYPE_NORMAL>
class TreeForForce {
#if defined(DEV_CODE)
   public:
#else
   private:
#endif

    DomainInfo *p_domain_info_;
    CommInfo comm_info_;
    // MortonKey morton_key_;
    MortonKey<typename TSM::mkey_type> morton_key_;
    F64ort inner_boundary_of_local_tree_;
    F64ort outer_boundary_of_local_tree_;
    F64 length_;     // length of a side of the root cell
    F64vec center_;  // new member (not used)
    F64ort pos_root_cell_;
    TimeProfile time_profile_;
    using MyMap = std::unordered_map<S64, Tepj *>;

    MyMap map_id_to_epj_;

    enum EXCHANGE_LET_MODE exchange_let_mode_ = EXCHANGE_LET_A2A;
    enum SET_ROOT_CELL_MODE set_root_cell_mode_ = SET_ROOT_CELL_AUTO;

    CountT n_interaction_ep_ep_local_, n_interaction_ep_sp_local_, n_walk_local_;
    CountT *n_cell_open_;

    bool is_initialized_;

    S64 n_interaction_ep_ep_;
    S64 n_interaction_ep_sp_;
    S32 ni_ave_;
    S32 nj_ave_;

    // RadixSort<KeyT, 8> rs_;
    S32 n_loc_tot_;  // # of all kinds of assigned particles in local proc
    S64 n_let_sp_;
    S64 n_glb_tot_;  // n_loc_tot_ + LETs
    S32 n_leaf_limit_;
    S32 n_group_limit_;
    F64 theta_;

    using TreeCellLoc = TreeCell<Tmomloc, Geometry<typename TSM::tree_cell_loc_geometry_type>>;
    using TreeCellGlb = TreeCell<Tmomglb, Geometry<typename TSM::tree_cell_glb_geometry_type>>;

    ReallocatableArray<TreeParticle> tp_glb_;  // not removed for neighbour search
    ReallocatableArray<TreeCellLoc> tc_loc_;   // not removed for reusing method (calc moment LT)
    ReallocatableArray<TreeCellGlb> tc_glb_;   // not removed for neighbour search

    ReallocatableArray<Tepj> epj_sorted_;  // not removed for neighbour search
    ReallocatableArray<Tepj> epj_send_;    // not removed (sometimes needed (no API))
    ReallocatableArray<Tspj> spj_send_;    // not removed (sometimes needed (no API))

    ReallocatableArray<IPGroup<typename TSM::ipg_type>> ipg_;  // not removed for reusing method
    ReallocatableArray<Tforce> force_org_;
    ReallocatableArray<S32> adr_org_from_adr_sorted_loc_;
    ReallocatableArray<U32> adr_org_from_adr_sorted_glb_;

    // added in 2025/02/22
    ReallocatableArray<S32> adr_ep_org_from_adr_ep_sorted_glb_;  // for reuse to avoid conditional branch to distinguish epj and spj
    ReallocatableArray<S32> adr_sp_org_from_adr_sp_sorted_glb_;  // for reuse to avoid conditional branch to distinguish epj and spj
    // added in 2025/02/22

    ReallocatableArray<Tepi> epi_sorted_;  // msortLT --- final
    ReallocatableArray<Tspj> spj_sorted_;  //  --- final
    ReallocatableArray<Tepi> epi_org_;     // setPtclLT ---
    ReallocatableArray<Tspj> spj_org_;     // insted of it, use spj_recv
    ReallocatableArray<Tepj> epj_org_;
    ReallocatableArray<Tforce> force_sorted_;  // -- final

    // ReallocatableArray<Tepj> *epj_for_force_;
    // ReallocatableArray<Tspj> *spj_for_force_;
    ReallocatableArray<S32> *adr_epj_for_force_;
    ReallocatableArray<S32> *adr_spj_for_force_;
    ReallocatableArray<S32> *adr_ipg_for_force_;

    // for gather mode
    class EPJWithR {
        Tepj epj_;
        F64 r_search_;

       public:
        Tepj getEPJ() const { return epj_; }
        F64vec getPos() const { return epj_.getPos(); }
        F64 getCharge() const { return epj_.getCharge(); }
        F64 getRSearch() const { return r_search_; }
        void copyFromEPJ(const Tepj &epj) { epj_ = epj; }
        void setRSearch(const F64 r_search) { r_search_ = r_search; }
        void setPos(const F64vec &pos) { epj_.setPos(pos); }
    };
    ReallocatableArray<EPJWithR> epjr_sorted_;  // cirectly copied from EPJ + RSearch
    ReallocatableArray<EPJWithR> epjr_send_;
    ReallocatableArray<EPJWithR> epjr_recv_;
    ReallocatableArray<EPJWithR> epjr_recv_1st_buf_;
    ReallocatableArray<EPJWithR> epjr_recv_2nd_buf_;
    ReallocatableArray<EPJWithR> *epjr_send_buf_;  // for 1st communication
    ReallocatableArray<EPJWithR> *epjr_send_buf_for_scatter_;
    ReallocatableArray<EPJWithR> *epjr_recv_1st_sorted_;

#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
    MPI_Request *req_send_;
    MPI_Request *req_recv_;
    MPI_Status *status_;
#endif

    void checkModeConsistency(const DomainInfo &dinfo) {
        bool fail_flag = false;
        if (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE> || std::is_same_v<TSM, SEARCH_MODE_LONG_CUTOFF>) {
            if (dinfo.getBoundaryCondition() != BOUNDARY_CONDITION_PERIODIC_XYZ) {
                PARTICLE_SIMULATOR_PRINT_ERROR(
                    "SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE or SEARCH_MODE_LONG_CUTOFF is only supported for the periodic boundary condition.");
                fail_flag = true;
            }
            for (S32 k = 0; k < DIMENSION; k++) {
                if (dinfo.getFlagSetPosRootDomain(k) == false) {
                    PARTICLE_SIMULATOR_PRINT_ERROR("The root domain is not set.");
                    fail_flag = true;
                }
            }
        }
        if (fail_flag) {
            PARTICLE_SIMULATOR_PRINT_ERROR("The mode is not consistent with the domain info.");
            Abort(-1);
        }
    }

    ////////////////////////////
    // for use calcForceAllAndWriteBack
    template <typename Head, typename... Tails>
    void setCalcForceParams(bool &clear_force, INTERACTION_LIST_MODE &list_mode, bool &flag_serialize, Head &head, Tails &&...tails) {
        if constexpr (IsDomainInfo<Head>::value) {
            setCalcForceParamsTail(clear_force, list_mode, flag_serialize, head, tails...);
        } else {
            setCalcForceParams(clear_force, list_mode, flag_serialize, tails...);
        }
    }
    void setCalcForceParamsTail(bool &clear_force, INTERACTION_LIST_MODE &list_mode, bool &flag_serialize, DomainInfo &dinfo,
                                const bool clear_force_ = true, const INTERACTION_LIST_MODE list_mode_ = FDPS_DFLT_VAL_LIST_MODE,
                                const bool flag_serialize_ = false) {
        checkModeConsistency(dinfo);
        p_domain_info_ = const_cast<DomainInfo *>(&dinfo);
        clear_force = clear_force_;
        list_mode = list_mode_;
        flag_serialize = flag_serialize_;
        if (flag_serialize == true) {
            PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
            Abort(-1);
        }
    }

    template <typename Head, typename... Tails>
    void setCalcForceParamsMultiWalk(S32 &n_walk_limit, bool &clear_force, INTERACTION_LIST_MODE &list_mode, bool &flag_serialize, Head &head,
                                     Tails &&...tails) {
        if constexpr (IsDomainInfo<Head>::value) {
            setCalcForceParamsTailMultiWalk(n_walk_limit, clear_force, list_mode, flag_serialize, head, tails...);
        } else {
            setCalcForceParamsMultiWalk(n_walk_limit, clear_force, list_mode, flag_serialize, tails...);
        }
    }
    void setCalcForceParamsTailMultiWalk(S32 &n_walk_limit, bool &clear_force, INTERACTION_LIST_MODE &list_mode, bool &flag_serialize,
                                         DomainInfo &dinfo, const S32 &n_walk_limit_, const bool clear_force_ = true,
                                         const INTERACTION_LIST_MODE list_mode_ = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize_ = false) {
        checkModeConsistency(dinfo);
        p_domain_info_ = const_cast<DomainInfo *>(&dinfo);
        n_walk_limit = n_walk_limit_;
        clear_force = clear_force_;
        list_mode = list_mode_;
        flag_serialize = flag_serialize_;
        if (flag_serialize == true) {
            PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
            Abort(-1);
        }
    }

    // for use calcForceAndWriteBack
    template <typename Head, typename... Tails>
    void setCalcForceParams2(bool &clear_force, INTERACTION_LIST_MODE &list_mode, bool &flag_serialize, Head &&head, Tails &&...tails) {
        if constexpr (IsParticleSystem<Head>::value) {
            setCalcForceParams2(clear_force, list_mode, flag_serialize, std::forward<Tails>(tails)...);
        } else {
            setCalcForceParamsTail2(clear_force, list_mode, flag_serialize, std::forward<Head>(head), std::forward<Tails>(tails)...);
        }
    }
    void setCalcForceParamsTail2(bool &clear_force, INTERACTION_LIST_MODE &list_mode, bool &flag_serialize, const bool clear_force_ = true,
                                 const INTERACTION_LIST_MODE list_mode_ = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize_ = false) {
        clear_force = clear_force_;
        list_mode = list_mode_;
        flag_serialize = flag_serialize_;
        if (flag_serialize == true) {
            PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
            Abort(-1);
        }
    }
    // for use calcForceAllAndWriteBackVA
    ////////////////////////////

   public:
    // added 2025/02/22
    LET_USAGE_MODE let_usage_mode_ = LET_USAGE_MODE::CONSTRUCT_GLOBAL_TREE;
    // LET_USAGE_MODE let_usage_mode_ = LET_USAGE_MODE::CONSTRUCT_LET_TREE;
    //  added 2025/02/22

    ////////////////////////////
    /// HIGH LEVEL FUNCTIONS
    /// in tree_for_force_high_level.hpp
    // Theses VA functions are not open to the user.
    /*
    template <typename Tfunc_ep_ep, typename Head, typename... Tails>
    void calcForceAllAndWriteBackVA(Tfunc_ep_ep pfunc_ep_ep, Head && head, Tails &&...tails);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
    S32 calcForceAllAndWriteBackMultiWalkVA(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
    S32 calcForceAllAndWriteBackMultiWalkIndexVA(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args);
    */
// #if defined(PARTICLE_SIMULATOR_USE_VARIADIC_TEMPLATE)
#if 1
    template <typename Tfunc_ep_ep, typename Head, typename... Tails>
    void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Head &&head, Tails &&...tails);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
    S32 calcForceAllAndWriteBackMultiWalk(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
    S32 calcForceAllAndWriteBackMultiWalkIndex(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args);
#else
    template <class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
    void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, Tpsys &psys, DomainInfo &dinfo, const bool clear_force = true,
                                  const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_ep_ep, class Tpsys>
    void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Tpsys &psys, DomainInfo &dinfo, const bool clear_force = true,
                                  const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
    S32 calcForceAllAndWriteBackMultiWalk(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Tpsys &psys,
                                          DomainInfo &dinfo, const S32 n_walk_limit, const bool clear = true,
                                          const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
    S32 calcForceAllAndWriteBackMultiWalkIndex(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Tpsys &psys,
                                               DomainInfo &dinfo, const S32 n_walk_limit, const bool clear = true,
                                               const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
#endif

/*
    template <typename Tfunc_ep_ep, typename Head, typename... Tails>
    void calcForceAllVA(Tfunc_ep_ep pfunc_ep_ep, Head && head, Tails &&...tails);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
    S32 calcForceAllMultiWalkVA(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
    S32 calcForceAllMultiWalkIndexVA(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args);
*/
// #if defined(PARTICLE_SIMULATOR_USE_VARIADIC_TEMPLATE)
#if 1
    // calcForceAll* (set FP & make tree & make ipg & make list & calc force) (wihtout write back)
    template <typename Tfunc_ep_ep, typename Head, typename... Tails>
    void calcForceAll(Tfunc_ep_ep pfunc_ep_ep, Head &&head, Tails &&...tails);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
    S32 calcForceAllMultiWalk(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve, typename... Args>
    S32 calcForceAllMultiWalkIndex(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Args &&...args);
#else
    template <class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
    void calcForceAll(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, Tpsys &psys, DomainInfo &dinfo, const bool clear_force = true,
                      const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_ep_ep, class Tpsys>
    void calcForceAll(Tfunc_ep_ep pfunc_ep_ep, Tpsys &psys, DomainInfo &dinfo, const bool clear_force = true,
                      const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
    S32 calcForceAllMultiWalk(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Tpsys &psys, DomainInfo &dinfo,
                              const S32 n_walk_limit, const bool clear = true, const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                              const bool flag_serialize = false);
    template <class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
    S32 calcForceAllMultiWalkIndex(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, Tpsys &psys, DomainInfo &dinfo,
                                   const S32 n_walk_limit, const bool clear = true, const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                                   const bool flag_serialize = false);
    // calcForceAll* (set FP & make tree & make ipg & make list & calc force)
#endif

    // calcFoceMakingTree* (make tree & make ipg & make list & calc force)
    template <class Tfunc_ep_ep, class Tfunc_ep_sp>
    void calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, DomainInfo &dinfo, const bool clear_force = true,
                             const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_ep_ep>
    void calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep, DomainInfo &dinfo, const bool clear_force = true,
                             const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_ep_ep, class Tfunc_ep_sp>
    void calcForceMakingTreeLong(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, DomainInfo &dinfo, const bool clear_force = true,
                                 const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_ep_ep>
    void calcForceMakingTreeShort(Tfunc_ep_ep pfunc_ep_ep, DomainInfo &dinfo, const bool clear_force = true,
                                  const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Tfunc_dispatch, class Tfunc_retrieve>
    S32 calcForceMakingTreeMultiWalk(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, DomainInfo &dinfo,
                                     const S32 n_walk_limit, const bool clear = true, const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                                     const bool flag_serialize = false);
    template <class Tfunc_dispatch, class Tfunc_retrieve>
    S32 calcForceMakingTreeMultiWalkIndex(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, DomainInfo &dinfo,
                                          const S32 n_walk_limit, const bool clear = true,
                                          const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    // calcFoceMakingTree* (make tree & make ipg & make list & calc force)

#if 1
    // calcForceAndWriteBack* (calc force & write back)
    template <typename Tfunc_ep_ep, typename Head, typename... Tails>
    void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Head &&head, Tails &&...tails);
    // calcForceAndWriteBack* (calc force & write back)
#else
    // calcForceAndWriteBack* (calc force & write back)
    template <class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
    void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, Tpsys &psys, const bool clear = true,
                               const bool flag_serialize = false);
    template <class Tfunc_ep_ep, class Tpsys>
    void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, Tpsys &psys, const bool clear = true, const bool flag_serialize = false);
    // calcForceAndWriteBack* (calc force & write back)
#endif
    // PRIVATE FUNCTIONS
    // make all tree. This function is invoked in the functions calcForceMakingTreeLong and calcForceMakingTreeShort.
    void makeTreeAll(DomainInfo &dinfo, const INTERACTION_LIST_MODE list_mode, const bool flag_serialize);
    /// HIGH LEVEL FUNCTIONS
    ////////////////////////////

    ////////////////////////////
    /// LOW LEVEL FUNCTIONS
    Tforce getForce(const S32 i) const { return force_org_[i]; }
    Tepi getEPISorted(const S32 i) const { return epi_sorted_[i]; }
    void copyLocalTreeStructure();  // not yet
    void repeatLocalCalcForce();    // not yet
    /// LOW LEVEL FUNCTIONS
    ////////////////////////////

    //////////////////////////////
    /// TREE WALK & CALC FORCE
    /// in tree_for_force_impl_force.hpp
    // make ipg & make list & calc force
    template <class Tfunc_ep_ep, class Tfunc_ep_sp>
    void calcForceMakingGroupLong(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear_force = true, const bool copy_force = true,
                                  const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const TREE_TYPE tree_type = GLOBAL_TREE);
    template <typename Tfunc_ep_ep>
    void calcForceMakingGroupShort(Tfunc_ep_ep pfunc_ep_ep, const bool clear_force = true, const bool copy_force = true,
                                   const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const TREE_TYPE tree_type = GLOBAL_TREE);
    template <typename Tkernel_type, typename Tfunc_dispatch, typename Tfunc_retrieve>
    S32 calcForceMakingGroupMultiWalk(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit,
                                      const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool copy_force = true,
                                      const bool clear = true);
    template <typename Tkernel_type, typename Tfunc_dispatch, typename Tfunc_retrieve>
    S32 calcForceMultiWalk(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit,
                           const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool clear = true);
    // walk + force. This function is invoked in the functions calcForceMultiWalkIndex and calcForceMultiWalkPtcl.
    template <typename KERNEL_TYPE, class Tfunc_dispatch, class Tfunc_retrieve>
    S32 calcForceMultiWalkImpl(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 tag_max, const S32 n_walk_limit,
                               const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool clear = true);

    /// CALC FORCE & TREE WALK
    ///////////////////////////

    ////////////////////////////
    // SET ROOT CELL
    void setRootCell(const DomainInfo &dinfo);
    void setRootCell(const F64 l, const F64vec &c = F64vec(0.0));
    void setRootCellMode(enum SET_ROOT_CELL_MODE rcm) { set_root_cell_mode_ = rcm; }
    // SET ROOT CELL
    ////////////////////////////

    ////////////////////////////
    // SET PARTICLE LOCAL TREE (not public)
    template <class Tpsys>
    void setParticleLocalTreeImpl(const Tpsys &psys, const bool clear);
    template <typename Head, typename... Tail>
    void setParticleLocalTreeVAImpl(Head &head, Tail &...tail);
    template <typename... Args>
    void setParticleLocalTreeVA(Args &&...args);
    // SET PARTICLE LOCAL TREE
    ////////////////////////////

    ////////////////////////////
    // helper functons
    template <typename Head, typename... Tail>
    void countNumberOfParticleSystems(Head &head, Tail &&...tail);
    // SET PARTICLE LOCAL TREE
    ////////////////////////////

    /////////////
    // WRITE BACK (not public)
    template <typename Head, typename... Tail>
    void writeBackImpl(S32 &offset, Head &head, Tail &&...tails);
    // WRITE BACK
    /////////////

    ///////////////
    // MAKE LOCAL TREE
    void makeLocalTreeOnly(DomainInfo &dinfo, const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);
    template <class Ttree>
    void copyRootCell(const Ttree &tree);
    void mortonSortLocalTreeOnly(const bool reuse = false);
    void linkCellLocalTreeOnly();
    void calcMomentLocalTreeOnly();
    void calcMomentGlobalTreeOnlyImpl();
    // MAKE LOCAL TREE
    ///////////////

    void makeInteractionListLongForZeroTheta(TagWithoutCutoff, const S32 adr_ipg);
    void makeInteractionListLongForZeroTheta(TagWithCutoff, const S32 adr_ipg);
    void makeInteractionListImpl(TagForceLong, const S32 adr_ipg, const bool clear);
    void makeInteractionListImpl(TagForceShort, const S32 adr_ipg, const bool clear);
    void makeInteractionListIndexImpl(TagSearchLong, const S32 adr_ipg, const bool clear);

    void calcCenterAndLengthOfRootCell(const DomainInfo &dinfo);
    // void calcCenterAndLengthOfRootCell(const DomainInfo & dinfo, const F64ort & outer_boundary);

    // for compile
    template <typename Tep2>
    void calcCenterAndLengthOfRootCellImpl(const Tep2 ep[], const DomainInfo &dinfo);
    template <typename Tep2>
    void calcCenterAndLengthOfRootCellImpl(const Tep2 ep[], const DomainInfo &dinfo, const F64ort &outer_boundary);

    void checkMortonSortGlobalTreeOnlyImpl(TagForceLong, std::ostream &fout);
    void checkMortonSortGlobalTreeOnlyImpl(TagForceShort, std::ostream &fout);

    void checkMakeInteractionListImpl(TagSearchLong, const DomainInfo &dinfo, const S32 adr_ipg, const S32 ith, const F64 tolerance,
                                      std::ostream &fout);
    void checkMakeInteractionListImpl(TagSearchLongCutoff, const DomainInfo &dinfo, const S32 adr_ipg, const S32 ith, const F64 tolerance,
                                      std::ostream &fout);
    void checkMakeInteractionListImpl(TagSearchShortScatter, const DomainInfo &dinfo, const S32 adr_ipg, const S32 ith, const F64 tolerance,
                                      std::ostream &fout);
    void checkMakeInteractionListImpl(TagSearchShortGather, const DomainInfo &dinfo, const S32 adr_ipg, const S32 ith, const F64 tolerance,
                                      std::ostream &fout);
    void checkMakeInteractionListImpl(TagSearchShortSymmetry, const DomainInfo &dinfo, const S32 adr_ipg, const S32 ith, const F64 tolerance,
                                      std::ostream &fout);

    void freeObjectFromMemoryPool() {
#ifdef PARTICLE_SIMULATOR_USE_MEMORY_POOL
        epi_sorted_.freeMem();
        spj_sorted_.freeMem();
        force_sorted_.freeMem();
        for (int i = 0; i < Comm::getNumberOfThread(); i++) {
            epj_for_force_[i].freeMem();
            spj_for_force_[i].freeMem();
        }
        epi_org_.freeMem();
        spj_org_.freeMem();
        epj_org_.freeMem();
        epj_send_.freeMem();
        spj_send_.freeMem();
        assert(MemoryPool::getSize() == 0);
#else
        for (int i = 0; i < Comm::getNumberOfThread(); i++) {
            // spj_for_force_[i].freeMem();
            // epj_for_force_[i].freeMem();
        }
        FreeMemVA(spj_send_, epj_send_, epi_org_);
        // spj_send_.freeMem();
        // epj_send_.freeMem();
        // epi_org_.freeMem();
#endif
    }

   public:
    using TypeEpi = Tepi;
    using TypeEpj = Tepj;

    void setCoefIpgSizeLimit(const F64 c) { coef_ipg_size_limit_ = c; }
    ///////
    // PMMM
    F64 coef_ipg_size_limit_;
    ReallocatableArray<TreeParticle> tp_loc_;  // [NEW] not removed for the reuse of tp_loc_ in other TreeForForce
    S32vec n_cell_pm3_;
    S32 i_cut_;
    std::unordered_map<S32, S32> adr_tc_loc_from_pm_cell_idx_;  // Data to obtain an array index of tc_*_ from
    std::unordered_map<S32, S32> adr_tc_let_from_pm_cell_idx_;  // a given 1D PM cell index. These data are used
    std::unordered_map<S32, S32> adr_tc_glb_from_pm_cell_idx_;  // in calcMultipoleMomentOfParticleMeshCell() &
    InteractionList interaction_list_glb_;
    void correspondPM3Cells2TreeCellAndLocalDomain(const TREE_TYPE tree_type);
    // PUBLIC API
    void setParamPMMM(const S32vec n_cell_pm3, const S32 i_cut) {
        n_cell_pm3_ = n_cell_pm3;
        i_cut_ = i_cut;
        bool err{false};
        if (i_cut_ < 1) {
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The cell separation must be >= 1.");
            std::cout << "i_cut_ = " << i_cut_ << std::endl;
        }
        if ((n_cell_pm3_.x <= 0) || (n_cell_pm3_.y <= 0) || (n_cell_pm3_.z <= 0)) {
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The number of the cells in each dimension must be > 0.");
            std::cout << "n_cell_pm3_ = " << n_cell_pm3_ << std::endl;
        }
        if (err) {
            Abort(-1);
        }
    }
    void setParamPMMM(const S32 n_cell_pm3_1d, const S32 i_cut) { setParamPMMM(S32vec(n_cell_pm3_1d), i_cut); }
    S32 getICut() const { return i_cut_; }
    S32vec getNCell() const { return n_cell_pm3_; }
    void repositionPosRootCellToEmbedParticleMesh();

    F64ort pos_unit_cell_;      // the domain size of the part of the particle mesh used for PMMM
                                // calculation and it is also used as the reference position for
                                // Morton key calculation.
    S32vec idx_pm3_unit_cell_;  // If we divide pos_root_cell into PM cells and PM cells are
                                // indexed such that (0,0,0) corresponds to the PM cell at
                                // the lower bound of pos_root_cell, idx_unit_cell_ describes
                                // the position of a PM cell at the lower bound of the particle
                                // mesh used for the PMMM calculation.
    F64vec width_pm_cell_;      // the width of a PM cell
    S32 lev_pm_cell_;           // the tree level at which a PM cell matches with a tree cell
    S32 lev_leaf_limit_;        // the level of leaf cells must be equal to or greater than this value.
    S32 lev_group_limit_;       // the level of tree cells corresponding to i-particle group must be
                                // equal to or greather than this value.

    std::vector<S32ort> pos_pm_domain_;
    std::vector<S32> pm_cell_idx_loc_;  // a list of 1-d indices of PM cells that overlap with the local domain
    std::vector<S32> pm_cell_idx_mm_;
    std::unordered_map<S32, std::vector<S32>> rank_list_from_pm_cell_idx_loc_;
    std::unordered_map<S32, IParticleInfo<Tepi, Tforce>> ip_info_from_pm_cell_idx_;  // i particle information of PM cell
    void repositionPosRootCellToEmbedParticleMeshImpl();
    void initializeMortonKey(TagMortonKeyNormal);
    void initializeMortonKey(TagMortonKeyMeshBased);
    void initializeMortonKey(TagMortonKeyMeshBased2);
    void decomposeParticleMeshCell();
    void calcIPInfoFromParticleMeshCellIndex();
    // PMMM
    ///////

    void makeAndExchangeLETOnly(DomainInfo &dinfo, const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                                const bool flag_serialize = false);

    void makeGlobalTreeOnly(DomainInfo &dinfo, const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const bool flag_serialize = false);

    void calcAdrTCFromParticleMeshCellIndex(const TREE_TYPE tree_type);

    void makeIPGroupOnly(const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE);

    void makeInteractionListIndexLong(const TREE_TYPE tree_type, const S32 adr_ipg_first, const S32 n_ipg_this_rnd);

    void makeInteractionListIndexShort(const TREE_TYPE tree_type, const S32 adr_ipg_first, const S32 n_ipg_this_rnd);

    void makeJPListOnly20241019(const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE, const TREE_TYPE tree_type = GLOBAL_TREE);
    void makeJPListOnly(const S32 adr_ipg_first, const S32 n_ipg_this_rnd, const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                        const TREE_TYPE tree_type = GLOBAL_TREE);

    /*
template <typename Tfunc_ep_ep, typename Tfunc_ep_sp>
void calcForceNoWalkPMM(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear_force, const bool copy_force,
        const TREE_TYPE tree_type);
*/
    template <typename Tfunc_ep_ep, typename Tfunc_ep_sp>
    void calcForceNoWalkPMM(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear_force, const bool copy_force,
                            const TREE_TYPE tree_type, const PS::S32 first_adr_ipg, const PS::S32 n_ipg_this_rnd);

    F64ort getPosUnitCell() const { return pos_unit_cell_; }

    std::vector<S32> getParticleMeshCellIndexLocal() const { return pm_cell_idx_loc_; }

    S32 getIParticleInfoOfParticleMeshCell(const S32 idx, IParticleInfo<Tepi, Tforce> &info) {
        auto itr = ip_info_from_pm_cell_idx_.find(idx);
        if (itr != ip_info_from_pm_cell_idx_.end()) {
            info = ip_info_from_pm_cell_idx_[idx];
            // [Note (tag: #fe161b03)]
            //     In order to add a const qulifier to this function,
            //     we have to use std::unordered_map::at instead of operator [].
            //     However, Fujitsu C++ compiler does not support `at`.
            //     Hence, we gave up on making this function a const function.
            //     (see also Note #cee527db in particle_mesh_multipole/M2L_engine.hpp)
            return 0;
        } else {
            return -1;
        }
    }

    void calcTotalChargeAndDispersionOfParticleMeshCellFull(const S32 idx, const F64vec &center, const F64ort &pos_my_domain, F64 &msum,
                                                            F64 &quad0) const {
        // This function calculates the sum of charges and dispersion for
        // PM cell whose one-dimensional cell index is `idx` (for details,
        // see the quantities in the brackets of the 2nd and 3rd terms in
        // the right hand side of Eq.(B.7) in Nitadori 2014).
        // The calculation involves local particles only. Hence, essentially,
        // we need the information of the local tree only. In the case that
        // the combination of the local tree and the LET tree is used for the
        // interaction calculation, we use the local tree only.
        // But, if the global tree is used for the interaction calculation,
        // epj_sorted_ is sorted for the global tree. So, we use the global
        // tree with `pos_my_domain`, which is used to identify local particles.

        // assert(is_global_tree_constructed_ || is_LET_tree_constructed_);
        //  This assertion is to force the users to use this function after
        //  the PP part is completed.
        msum = quad0 = 0;
#if 1
        auto itr = adr_tc_glb_from_pm_cell_idx_.find(idx);
        if (itr != adr_tc_glb_from_pm_cell_idx_.end()) {
            const S32 adr_tc = itr->second;
            CalcTotalChargeAndDispersionOfParticleMeshCellFull<TreeCellGlb, TreeParticle, Tepj>(adr_tc, tc_glb_, tp_glb_, epj_sorted_, center,
                                                                                                pos_my_domain, msum, quad0);
        }
#else
        if (is_global_tree_constructed_) {
            auto itr = adr_tc_glb_from_pm_cell_idx_.find(idx);
            if (itr != adr_tc_glb_from_pm_cell_idx_.end()) {
                const S32 adr_tc = itr->second;
                CalcTotalChargeAndDispersionOfParticleMeshCellFull<TreeCellGlb, TreeParticle, Tepj>(adr_tc, tc_glb_, tp_glb_, epj_sorted_, center,
                                                                                                    pos_my_domain, msum, quad0);
            }
        } else if (is_LET_tree_constructed_) {
            auto itr = adr_tc_loc_from_pm_cell_idx_.find(idx);
            if (itr != adr_tc_loc_from_pm_cell_idx_.end()) {
                const S32 adr_tc = itr->second;
                CalcTotalChargeAndDispersionOfParticleMeshCellFull<TreeCellLoc, TreeParticle, Tepj>(adr_tc, tc_loc_, tp_loc_, epj_sorted_, center,
                                                                                                    pos_my_domain, msum, quad0);
            }
        }
#endif
    }

    void calcTotalChargeAndDispersionOfParticleMeshCellStripe(const S32 idx, const F64vec &center, const F64ort &pos_my_domain, F64 &msum, F64 &quad0,
                                                              const S32 n_thread, const S32 thread_id) const {
        // Multi-thread version
        // assert(is_global_tree_constructed_ || is_LET_tree_constructed_);
        // This assertion is to force the users to use this function after
        // the PP part is completed.
        msum = quad0 = 0;
#if 1
        auto itr = adr_tc_glb_from_pm_cell_idx_.find(idx);
        if (itr != adr_tc_glb_from_pm_cell_idx_.end()) {
            const S32 adr_tc = itr->second;
            CalcTotalChargeAndDispersionOfParticleMeshCellStripe<TreeCellGlb, TreeParticle, Tepj>(adr_tc, tc_glb_, tp_glb_, epj_sorted_, center,
                                                                                                  pos_my_domain, msum, quad0, n_thread, thread_id);
        }
#else
        // original
        if (is_global_tree_constructed_) {
            auto itr = adr_tc_glb_from_pm_cell_idx_.find(idx);
            if (itr != adr_tc_glb_from_pm_cell_idx_.end()) {
                const S32 adr_tc = itr->second;
                CalcTotalChargeAndDispersionOfParticleMeshCellStripe<TreeCellGlb, TreeParticle, Tepj>(
                    adr_tc, tc_glb_, tp_glb_, epj_sorted_, center, pos_my_domain, msum, quad0, n_thread, thread_id);
            }
        } else if (is_LET_tree_constructed_) {
            auto itr = adr_tc_loc_from_pm_cell_idx_.find(idx);
            if (itr != adr_tc_loc_from_pm_cell_idx_.end()) {
                const S32 adr_tc = itr->second;
                CalcTotalChargeAndDispersionOfParticleMeshCellStripe<TreeCellLoc, TreeParticle, Tepj>(
                    adr_tc, tc_loc_, tp_loc_, epj_sorted_, center, pos_my_domain, msum, quad0, n_thread, thread_id);
            }
        }
#endif
    }

    //template <class real_t, class cplx_t>
    //S32 calcMultipoleMomentOfParticleMeshCell(const S32 idx, const F64vec &center, const F64ort &pos_my_domain, MultipoleMoment<real_t, cplx_t> &mm,
    //                                          const S32 p_spj2mm) const;

    template <typename real_t, typename cplx_t, template<typename, typename> typename mm_t>
    S32 calcMultipoleMomentOfParticleMeshCell(const S32 idx, const F64vec &center, const F64ort &pos_my_domain, mm_t<real_t, cplx_t> &mm,
                                              const S32 p_spj2mm) const;    
    std::unordered_map<S32, std::vector<S32>> getRankListFromParticleMeshCellIndexLocal() const { return rank_list_from_pm_cell_idx_loc_; }

    S32ort getPosParticleMeshDomain(const S32 i) const { return pos_pm_domain_[i]; }

    std::vector<S32ort> getPosParticleMeshDomain() const { return pos_pm_domain_; }

    // PMMM
    ///////

    void setCommInfo(const CommInfo &c);
    F64ort getPosRootCell() const { return pos_root_cell_; }

    void setPrefixOfProfile(const char *str) {
        // profile.setPrefix(str);
    }

    void setExchangeLETMode(enum EXCHANGE_LET_MODE elm) { exchange_let_mode_ = elm; }

    // for neighbour search
    ReallocatableArray<Tepj> *epj_neighbor_;

    // new API from V8
    TimeProfile getTimeProfile() const { return time_profile_; }
    void dumpTimeProfile(std::ostream &fout = std::cout, const S32 level = 0) const { time_profile_.dump_tree_for_force(fout, level); }
    void clearTimeProfile() { time_profile_.clear(); }

    CountT getNumberOfWalkLocal() const { return n_walk_local_; }
    CountT getNumberOfInteractionEPEPLocal() const { return n_interaction_ep_ep_local_; }
    CountT getNumberOfInteractionEPSPLocal() const { return n_interaction_ep_sp_local_; }
    // CountT getNumberOfWalkGlobal() const { return Comm::getSum(n_walk_local_); }
    CountT getNumberOfWalkGlobal() const { return comm_info_.getSum(n_walk_local_); }
    // CountT getNumberOfInteractionEPEPGlobal() const { return Comm::getSum(n_interaction_ep_ep_local_); }
    CountT getNumberOfInteractionEPEPGlobal() const { return comm_info_.getSum(n_interaction_ep_ep_local_); }
    // CountT getNumberOfInteractionEPSPGlobal() const { return Comm::getSum(n_interaction_ep_sp_local_); }
    CountT getNumberOfInteractionEPSPGlobal() const { return comm_info_.getSum(n_interaction_ep_sp_local_); }

    CountT getNumberOfCellOpenLocal() const { return n_cell_open_[0]; }
    // CountT getNumberOfCellOpenGlobal() const {return Comm::getSum(n_cell_open_[0]); }
    CountT getNumberOfCellOpenGlobal() const { return comm_info_.getSum(n_cell_open_[0]); }
    CountT getNumberOfCellGlobal() const { return tc_glb_.size(); }

    CountT getNumberOfExLetAllgatherLocal() const { return comm_table_.rank_recv_allgather_.size(); }
    // CountT getNumberOfExLetAllgatherGlobal() const {return Comm::getSum(comm_table_.rank_recv_allgather_.size());}
    CountT getNumberOfExLetAllgatherGlobal() const { return comm_info_.getSum(comm_table_.rank_recv_allgather_.size()); }

    CountT getNumberOfExLetEpSendLocal() const { return comm_table_.n_ep_send_tot_; }
    CountT getNumberOfExLetSpSendLocal() const { return comm_table_.n_sp_send_tot_; }
    CountT getNumberOfExLetEpRecvLocal() const { return comm_table_.n_ep_recv_tot_; }
    CountT getNumberOfExLetSpRecvLocal() const { return comm_table_.n_sp_recv_tot_; }
    // CountT getNumberOfExLetEpGlobal() const {return Comm::getSum(comm_table_.n_ep_send_tot_);}
    CountT getNumberOfExLetEpGlobal() const { return comm_info_.getSum(comm_table_.n_ep_send_tot_); }
    // CountT getNumberOfExLetSpGlobal() const {return Comm::getSum(comm_table_.n_sp_send_tot_);}
    CountT getNumberOfExLetSpGlobal() const { return comm_info_.getSum(comm_table_.n_sp_send_tot_); }

    void clearNumberOfInteraction() { n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0; }
    void clearCounterAll() {
        // clearCommTable();
        clearNumberOfInteraction();
        clearTimeProfile();
    }

    TreeForForce() : is_initialized_(false) {}
    ~TreeForForce() {
        delete[] n_cell_open_;
        delete[] epjr_send_buf_;
        delete[] epjr_send_buf_for_scatter_;
        delete[] epjr_recv_1st_sorted_;
        delete[] epj_neighbor_;
        // delete[] spj_for_force_;
        // delete[] epj_for_force_;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        delete[] req_send_;
        delete[] req_recv_;
#endif
    }

    size_t getMemSizeUsed() const;
    size_t getUsedMemorySize() const;  // same function as getMemSizeUsed()

    void setNInteractionEPEP(const S64 n_ep_ep) { n_interaction_ep_ep_ = n_ep_ep; }
    void setNInteractionEPSP(const S64 n_ep_sp) { n_interaction_ep_sp_ = n_ep_sp; }

    S64 getNInteractionEPEP() const { return n_interaction_ep_ep_; }
    S64 getNInteractionEPSP() const { return n_interaction_ep_sp_; }

    void initialize(const U64 n_glb_tot, const F64 theta = FDPS_DFLT_VAL_THETA, const U32 n_leaf_limit = FDPS_DFLT_VAL_N_LEAF_LIMIT,
                    const U32 n_group_limit = FDPS_DFLT_VAL_N_GROUP_LIMIT);

    void reallocMem();
    void freeMem();
    void clearSizeOfArray();

    void setDomainInfo(const DomainInfo &dinfo);
    void setDomainInfo(DomainInfo &dinfo);
    void linkCellGlobalTreeOnly();
    void calcMomentGlobalTreeOnly();
    void makeIPGroup();
    void mortonSortGlobalTreeOnly(const bool reuse = false);
    CountT getNumberOfIPG() const { return (CountT)ipg_.size(); }
    CountT getNumberOfIPGLocal() const { return getNumberOfIPG(); }
    CountT getNumberOfIPGGlobal() const {
        CountT n_loc = getNumberOfIPGLocal();
        CountT n_glb = comm_info_.getSum(n_loc);
        return n_glb;
    }
    void makeInteractionList(const S32 adr_ipg, const bool clear = true);
    template <class Tfunc_ep_ep>
    void calcForceOnly(Tfunc_ep_ep pfunc_ep_ep, const S32 adr_ipg, const bool clear = true);
    template <class Tfunc_ep_ep, class Tfunc_ep_sp>
    void calcForceOnly(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const S32 adr_ipg, const bool clear = true);

    template <class Tfunc_ep_ep>
    void calcForceWalkOnly(Tfunc_ep_ep pfunc_ep_ep, const bool clear = true);

    ////////////
    // for neighbour search APIs
    template <class Tptcl>
    S32 getNeighborListOneParticle(const Tptcl &ptcl, Tepj *&epj);
    template <class Tptcl>
    S32 getNeighborListOneParticleImpl(TagNeighborSearchScatter, const Tptcl &ptcl, Tepj *&epj);
    template <class Tptcl>
    S32 getNeighborListOneParticleImpl(TagNeighborSearchGather, const Tptcl &ptcl, Tepj *&epj);
    template <class Tptcl>
    S32 getNeighborListOneParticleImpl(TagNeighborSearchSymmetry, const Tptcl &ptcl, Tepj *&epj);
#ifdef PARTICLE_SIMULATOR_CHECK_SEARCH_MODE
    // for test_all_search_mode
    template <class Tptcl>
    S32 getNeighborListOneParticleImpl(TagNeighborSearchNo, const Tptcl &ptcl, Tepj *&epj);
#endif
    template <class Tptcl>
    S32 getNeighborListOneParticle(const Tptcl &ptcl, const S32 n_ngb, Tepj *&epj);

    F64ort getOuterBoundaryOfLocalTree() { return getOuterBoundaryOfLocalTreeImpl(typename TSM::search_type()); }
    F64ort getInnerBoundaryOfLocalTree() { return inner_boundary_of_local_tree_; }
    F64ort getOuterBoundaryOfLocalTreeImpl(TagSearchLongSymmetry);
    F64ort getOuterBoundaryOfLocalTreeImpl(TagSearchShortSymmetry);

    ///////////////
    // UTILS
    Tepj *getEpjFromId(const S64 id, const Tepj *epj_tmp = NULL);

    ///////////////
    // DEBUG
    S32 getNumberOfEpiSorted() const { return epi_sorted_.size(); }
    S32 getNumberOfEpjSorted() const { return epj_sorted_.size(); }
    S32 getNumberOfSpjSorted() const { return spj_sorted_.size(); }

    /////////////////////////////////
    // FOR REUSING INTERACTION LIST
    void addMomentAsSp() {
        AddMomentAsSp(typename TSM::force_type(), tc_glb_, n_let_sp_, spj_sorted_);
    }

    /////////////////////////////////
    // EXCHANGE LET
    void exchangeLocalEssentialTree(const DomainInfo &dinfo, const bool flag_reuse = false);
    // LONG
    void exchangeLocalEssentialTreeImpl(TagSearchLong, const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeImpl(TagSearchLongCutoff, const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeImpl(TagSearchLongScatter, const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeImpl(TagSearchLongSymmetry, const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeImpl(TagSearchShortScatter, const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry, const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeImpl(TagSearchShortGather, const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeImpl(TagSearchLongParticleMeshMultipole, const DomainInfo &dinfo, const bool flag_reuse = false);

    void exchangeLocalEssentialTreeLongNew(const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeLongP2P(const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeLongCutoffP2P(const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeLongP2P2(const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeShortScatter(const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeShortScatterP2P(const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeShort2StepA2A(const DomainInfo &dinfo, const bool flag_reuse = false);
    void exchangeLocalEssentialTreeShort2StepP2P(const DomainInfo &dinfo, const bool flag_reuse = false);

    CommTable comm_table_;
    InteractionList interaction_list_loc_;

    template <class Tfunc_ep_ep>
    void calcForceNoWalkOrg(Tfunc_ep_ep pfunc_ep_ep, const bool clear = true);

    template <class Tfunc_ep_ep>
    void calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep, const bool clear, const PS::S32 first_adr_ipg, const PS::S32 n_ipg_this_rnd);
    template <class Tfunc_ep_ep, class Tfunc_ep_sp>
    void calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear = true);

    template <class Tfunc_dispatch, class Tfunc_retrieve>
    void calcForceNoWalkForMultiWalkIndex(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 n_walk_limit,
                                          const bool clear = true) {
        //F64 wtime_offset = GetWtime();
        calcForceNoWalkForMultiWalkIndexImpl(typename TSM::force_type(), pfunc_dispatch, pfunc_retrieve, n_walk_limit, clear);
        //time_profile_.calc_force += GetWtime() - wtime_offset;
    }

    template <typename Tfunc_dispatch, typename Tfunc_retrieve>
    void calcForceNoWalkForMultiWalkIndexImpl(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 n_walk_limit, const bool clear);
    template <typename Tfunc_dispatch, typename Tfunc_retrieve>
    void calcForceNoWalkForMultiWalkPtclImpl(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 n_walk_limit, const bool clear);

    template <class Tfunc_dispatch, class Tfunc_retrieve>
    void calcForceNoWalkForMultiWalkPtcl(Tfunc_dispatch pfunc_dispatch, Tfunc_retrieve pfunc_retrieve, const S32 n_walk_limit,
                                         const bool clear = true) {
        //F64 wtime_offset = GetWtime();
        calcForceNoWalkForMultiWalkPtclImpl(typename TSM::force_type(), pfunc_dispatch, pfunc_retrieve, n_walk_limit, clear);
        //time_profile_.calc_force += GetWtime() - wtime_offset;
    }

    void copyForceOriginalOrder();

    S32 adr_tc_level_partition_loc_[TREE_LEVEL_LIMIT + 2];
    S32 lev_max_loc_;
    S32 adr_tc_level_partition_glb_[TREE_LEVEL_LIMIT + 2];
    S32 lev_max_glb_;

    ///////
    // SET LET
    void setLocalEssentialTreeToGlobalTree(const bool flag_reuse = false);
    void setLocalEssentialTreeToGlobalTreeImpl(TagForceLong, const bool flag_reuse = false);
    void setLocalEssentialTreeToGlobalTreeImpl(TagForceShort, const bool flag_reuse = false);

    template <class Tfunc_ep_ep, class Tsys>
    void calcForceDirectParallel(Tfunc_ep_ep pfunc_ep_ep, Tforce force[], const Tsys &system, const DomainInfo &dinfo, const bool clear = true);

    template <class Tfunc_ep_ep>
    void calcForceDirect(Tfunc_ep_ep pfunc_ep_ep, Tforce force[], const DomainInfo &dinfo, const bool clear = true);

    template <class Tfunc_ep_ep>
    void calcForceDirectAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, const DomainInfo &dinfo, const bool clear = true);

    void exchangeLocalEssentialTreeUsingCommTable() {}

    void dump(std::ostream &fout = std::cout);
    void dumpMemSizeUsed(std::ostream &fout = std::cout);

    //////////////////////////////
    // DEPRECATED FUNCTIONS
    ///// public API /////
    template <typename Tpsys>
    void setParticleLocalTree(const Tpsys &psys, const bool clear = true);  // deprecated
    void makeLocalTree(DomainInfo &dinfo);                     // deprecated
    void makeLocalTree(const F64 l, F64vec &c = F64vec(0.0));  // deprecated
    void calcMomentGlobalTree();                               // deprecated
    void makeGlobalTree(DomainInfo &dinfo);                    // deprecated
    template <class Tfunc_ep_ep>
    void calcForce(Tfunc_ep_ep pfunc_ep_ep, const bool clear = true);  // deprecated
    template <class Tfunc_ep_ep, class Tfunc_ep_sp>
    void calcForce(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear = true);  // deprecated
    ///// private API /////
    template <class Tfunc_ep_ep, class Tfunc_ep_sp>
    void calcForceNormal(Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const bool clear = true,
                         const INTERACTION_LIST_MODE list_mode = MAKE_LIST);  // deprecated
    template <typename KERNEL_TYPE, typename Tfunc0, typename Tfunc1>
    S32 makeListAndCalcForceImpl(Tfunc0 func0, Tfunc1 func1, const S32 tag_max, const S32 n_walk_limit, const INTERACTION_LIST_MODE list_mode,
                                 const bool clear);
    template <typename Head, typename... Tail>
    void writeBack(Head &head, Tail &&...tails);
    // DEPRECATED FUNCTIONS
    //////////////////////////////
};

}  // namespace ParticleSimulator

#include "tree_for_force_type.hpp"
#include "tree_for_force_impl_set_root_cell.hpp"
#include "tree_for_force_impl.hpp"
#include "tree_for_force_impl_make_tree.hpp"
#include "tree_for_force_impl_make_jplist.hpp"
#include "tree_for_force_impl_pmm.hpp"
#include "tree_for_force_impl_high_level.hpp"
// #include"tree_for_force_check_impl.hpp"
#include "tree_for_force_deprecated.hpp"