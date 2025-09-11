#pragma once
namespace ParticleSimulator {
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
#if defined(PS_DEBUG_PRINT_FindScatterParticlePMM)
    PARTICLE_SIMULATOR_PRINT_MSG("****** START FindScatterParticlePMM ******", 0);
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
            if (typeid(TSM) != typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE)) {
                CalcNumberAndShiftOfImageDomain(shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(), outer_boundary_of_my_tree,
                                                dinfo.getPosDomain(i), pa, false);
            } else {
                CalcNumberAndShiftOfImageDomain(shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(), top_geometry[my_rank].vertex_in_,
                                                top_geometry[i].vertex_in_, pa, icut, morton_key, false);
            }
            const S32 n_image_tmp = shift_per_image_tmp[ith].size();
            n_image_per_proc[i] = n_image_tmp - n_image_tmp_prev;
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
                target_box.setForExLet(top_geometry[i], shift_per_image_tmp[ith][j], icut, morton_key);

                if (typeid(TSM) != typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE)) {
                    MakeListUsingTreeRecursiveTopPMM<TSM, Ttc, TreeParticle, Tep, Tsp, Tmorton_key, TagChopLeafFalse, TagChopNonleafFalse,
                                                     TagCopyInfoCloseNoSp, CALC_DISTANCE_TYPE>(
                        tc_first, adr_tc, tp_first, ep_first, adr_ep_send_tmp[ith], sp_first, adr_sp_send_tmp[ith], target_box, n_leaf_limit,
                        lev_leaf_limit, adr_tree_sp_first, len_peri, theta, lev_pm_cell, morton_key);
                } else {
                    MakeListUsingTreeRecursiveTop<TSM, Ttc, TreeParticle, Tep, Tsp, TagChopLeafFalse, TagCopyInfoCloseNoSp, CALC_DISTANCE_TYPE>(
                        tc_first, adr_tc, tp_first, ep_first, adr_ep_send_tmp[ith], sp_first, adr_sp_send_tmp[ith], target_box, n_leaf_limit,
                        adr_tree_sp_first, len_peri, theta);
                }
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

    S32 n_disp_image_per_proc[n_proc + 1];
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
    PARTICLE_SIMULATOR_PRINT_MSG("****** END FindScatterParticlePMM ******", 0);
#endif
}
}  // namespace ParticleSimulator