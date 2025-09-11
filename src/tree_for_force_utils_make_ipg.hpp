#pragma once
namespace ParticleSimulator {
template <class Tipg, class Ttc, class Tepi>
inline void MakeIPGroupLong(ReallocatableArray<Tipg> &ipg_first, const ReallocatableArray<Ttc> &tc_first, const ReallocatableArray<Tepi> &epi_first,
                            const S32 adr_tc, const S32 n_grp_limit) {
    const Ttc *tc_tmp = tc_first.getPointer(adr_tc);
    const S32 n_tmp = tc_tmp->n_ptcl_;
    if (n_tmp == 0)
        return;
    else if (tc_tmp->isLeaf(n_grp_limit)) {
        ipg_first.increaseSize();
        ipg_first.back().copyFromTC(*tc_tmp);
        // ipg_first.back().vertex_ = GetMinBoxSingleThread(epi_first.data()+(tc_tmp->adr_ptcl_), n_tmp);
        return;
    } else {
        S32 adr_tc_tmp = tc_tmp->adr_tc_;
        for (S32 i = 0; i < N_CHILDREN; i++) {
            MakeIPGroupLong(ipg_first, tc_first, epi_first, adr_tc_tmp + i, n_grp_limit);
            // MakeIPGroupLong<Tipg, Ttc, Tepi, Ti0, Ti1, Ti2>
            //     (ipg_first, tc_first, epi_first, adr_tc_tmp+i, n_grp_limit);
        }
    }
}

// NOTE: This FUnction mimic sirial code
template <class Tipg, class Ttcloc, class Ttcglb, class Tepj>
inline void MakeIPGroupLongGLBTreeCellAsIPGBox(ReallocatableArray<Tipg> &ipg_first, const ReallocatableArray<Ttcloc> &tc_loc_first,
                                               const ReallocatableArray<Ttcglb> &tc_glb_first, const ReallocatableArray<Tepj> &epj_first,
                                               const S32 adr_tc_loc, const S32 adr_tc_glb, const S32 n_grp_limit, const S32 n_leaf_limit) {
    const Ttcloc *tc_loc_tmp = tc_loc_first.getPointer(adr_tc_loc);
    const Ttcglb *tc_glb_tmp = tc_glb_first.getPointer(adr_tc_glb);
    const S32 n_loc_tmp = tc_loc_tmp->n_ptcl_;
    const S32 n_glb_tmp = tc_glb_tmp->n_ptcl_;
    if (n_loc_tmp == 0 || n_glb_tmp == 0)
        return;
    else if (tc_glb_tmp->isLeaf(n_grp_limit) || tc_loc_tmp->isLeaf(n_leaf_limit)) {
        ipg_first.increaseSize();
        ipg_first.back().copyFromTC(*tc_loc_tmp);
        const S32 adr_tmp = tc_glb_tmp->adr_ptcl_;
        ipg_first.back().vertex_ = GetMinBoxSingleThread(epj_first.getPointer(adr_tmp), n_glb_tmp);
        return;
    } else {
        S32 adr_tc_loc_tmp = tc_loc_tmp->adr_tc_;
        S32 adr_tc_glb_tmp = tc_glb_tmp->adr_tc_;
        for (S32 i = 0; i < N_CHILDREN; i++) {
            MakeIPGroupLongGLBTreeCellAsIPGBox(ipg_first, tc_loc_first, tc_glb_first, epj_first, adr_tc_loc_tmp + i, adr_tc_glb_tmp + i, n_grp_limit,
                                               n_leaf_limit);
        }
    }
}

// PMMM
template <class Tipg, class Ttcloc, class Ttcglb, class Tepi>
inline void MakeIPGroupUseGLBTreeLong(ReallocatableArray<Tipg> &ipg_first, const ReallocatableArray<Ttcloc> &tc_loc_first,
                                      const ReallocatableArray<Ttcglb> &tc_glb_first, const ReallocatableArray<Tepi> &epi_first, const S32 adr_tc_loc,
                                      const S32 adr_tc_glb, const S32 n_grp_limit, const S32 lev_grp_limit, const S32 n_leaf_limit,
                                      const S32 lev_leaf_limit, const F64 ipg_size_limit) {
    const Ttcloc *tc_loc_tmp = tc_loc_first.getPointer(adr_tc_loc);
    const Ttcglb *tc_glb_tmp = tc_glb_first.getPointer(adr_tc_glb);
    const S32 n_loc_tmp = tc_loc_tmp->n_ptcl_;
    const S32 n_glb_tmp = tc_glb_tmp->n_ptcl_;
    if (n_loc_tmp == 0 || n_glb_tmp == 0) {
        return;
    } else if (tc_loc_tmp->isLeafPMM(n_leaf_limit, lev_leaf_limit) || tc_glb_tmp->isGroup(n_grp_limit, lev_grp_limit, ipg_size_limit)) {
        ipg_first.increaseSize();
        ipg_first.back().copyFromTC(*tc_loc_tmp, adr_tc_loc);
        // const S32 adr_tmp = tc_loc_tmp->adr_ptcl_;
        // ipg_first.back().vertex_ = GetMinBoxSingleThread(epi_first.getPointer(adr_tmp), n_loc_tmp);
        return;
    } else {
        S32 adr_tc_loc_tmp = tc_loc_tmp->adr_tc_;
        S32 adr_tc_glb_tmp = tc_glb_tmp->adr_tc_;
        for (S32 i = 0; i < N_CHILDREN; i++) {
            MakeIPGroupUseGLBTreeLong(ipg_first, tc_loc_first, tc_glb_first, epi_first, adr_tc_loc_tmp + i, adr_tc_glb_tmp + i, n_grp_limit,
                                      lev_grp_limit, n_leaf_limit, lev_leaf_limit, ipg_size_limit);
        }
    }
}

template <class Tipg, class Ttcloc, class Ttcglb, class Tepi>
inline void MakeIPGroupUseGLBTreeShort(ReallocatableArray<Tipg> &ipg_first, const ReallocatableArray<Ttcloc> &tc_loc_first,
                                       const ReallocatableArray<Ttcglb> &tc_glb_first, const ReallocatableArray<Tepi> &epi_first,
                                       const S32 adr_tc_loc, const S32 adr_tc_glb, const S32 n_grp_limit, const S32 lev_grp_limit,
                                       const S32 n_leaf_limit, const S32 lev_leaf_limit, const F64 ipg_size_limit) {
    MakeIPGroupUseGLBTreeLong<Tipg, Ttcloc, Ttcglb, Tepi>(ipg_first, tc_loc_first, tc_glb_first, epi_first, adr_tc_loc, adr_tc_glb, n_grp_limit,
                                                          lev_grp_limit, n_leaf_limit, lev_leaf_limit, ipg_size_limit);
}

#if 1

#else
template <class Tipg, class Ttcloc, class Ttcglb, class Tepi>
inline void MakeIPGroupUseGLBTreeLong(ReallocatableArray<Tipg> &ipg_first, const ReallocatableArray<Ttcloc> &tc_loc_first,
                                      const ReallocatableArray<Ttcglb> &tc_glb_first, const ReallocatableArray<Tepi> &epi_first, const S32 adr_tc_loc,
                                      const S32 adr_tc_glb, const S32 n_grp_limit, const S32 n_leaf_limit) {
    const Ttcloc *tc_loc_tmp = tc_loc_first.getPointer(adr_tc_loc);
    const Ttcglb *tc_glb_tmp = tc_glb_first.getPointer(adr_tc_glb);
    const S32 n_loc_tmp = tc_loc_tmp->n_ptcl_;
    const S32 n_glb_tmp = tc_glb_tmp->n_ptcl_;
    if (n_loc_tmp == 0 || n_glb_tmp == 0)
        return;
    else if (tc_glb_tmp->isLeaf(n_grp_limit) || tc_loc_tmp->isLeaf(n_leaf_limit)) {
        ipg_first.increaseSize();
        ipg_first.back().copyFromTC(*tc_loc_tmp);
        return;
    } else {
        S32 adr_tc_loc_tmp = tc_loc_tmp->adr_tc_;
        S32 adr_tc_glb_tmp = tc_glb_tmp->adr_tc_;
        for (S32 i = 0; i < N_CHILDREN; i++) {
            MakeIPGroupUseGLBTreeLong(ipg_first, tc_loc_first, tc_glb_first, epi_first, adr_tc_loc_tmp + i, adr_tc_glb_tmp + i, n_grp_limit,
                                      n_leaf_limit);
        }
    }
}
#endif

template <class Tipg, class Ttc, class Tepi>
inline void MakeIPGroupShort(ReallocatableArray<Tipg> &ipg_first, const ReallocatableArray<Ttc> &tc_first, const ReallocatableArray<Tepi> &epi_first,
                             const S32 adr_tc, const S32 n_grp_limit) {
    const Ttc *tc_tmp = tc_first.getPointer() + adr_tc;
    const S32 n_tmp = tc_tmp->n_ptcl_;
    if (n_tmp == 0)
        return;
    else if (tc_tmp->isLeaf(n_grp_limit)) {
        ipg_first.increaseSize();
        // ipg_first.back().copyFromTC(*tc_tmp);
        ipg_first.back().copyFromTC(*tc_tmp, adr_tc);
        return;
    } else {
        S32 adr_tc_tmp = tc_tmp->adr_tc_;
        for (S32 i = 0; i < N_CHILDREN; i++) {
            MakeIPGroupShort<Tipg, Ttc, Tepi>(ipg_first, tc_first, epi_first, adr_tc_tmp + i, n_grp_limit);
        }
    }
}
}  // namespace ParticleSimulator