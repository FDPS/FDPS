namespace ParticleSimulator {
// PMMM and long mode
#if 1
template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tipg, class Tmorton_key, class Tchopleaf, class Tchopnonleaf,
          class Tcopyinfoclose, class Tinteraction_list, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void MakeInteractionListIndexLong(const ReallocatableArray<Ttc> &tc, const S32 adr_tc, const ReallocatableArray<Ttp> &tp, const S32 adr_tree_sp_first,
                                  const ReallocatableArray<Tep> &ep, const ReallocatableArray<Tsp> &sp, const ReallocatableArray<Tipg> &ipg,
                                  const S32 adr_ipg_first, const S32 n_ipg_this_rnd, const S32 n_leaf_limit, const S32 lev_leaf_limit,
                                  const F64 theta, const S32 icut, const S32 lev_pm_cell, const Tmorton_key &morton_key,
                                  Tinteraction_list &interaction_list, const F64vec len_peri) {
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START MakeInteractionListIndexLong *****", 0);
#endif
    const S32 n_ipg = ipg.size();
    const auto n_thread = Comm::getNumberOfThread();
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
    PS_PRINT_VARIABLE2(n_ipg, n_thread);
#endif

    if (n_ipg_this_rnd > 0) {
        assert(0 <= adr_ipg_first && adr_ipg_first < n_ipg);
        assert(n_ipg_this_rnd <= n_ipg);
        interaction_list.n_ep_.resizeNoInitialize(n_ipg_this_rnd);
        interaction_list.n_disp_ep_.resizeNoInitialize(n_ipg_this_rnd + 1);
        interaction_list.n_sp_.resizeNoInitialize(n_ipg_this_rnd);
        interaction_list.n_disp_sp_.resizeNoInitialize(n_ipg_this_rnd + 1);

        ReallocatableArray<S32> adr_epj_tmp[n_thread];
        ReallocatableArray<S32> adr_spj_tmp[n_thread];
        ReallocatableArray<S32> adr_ipg_tmp[n_thread];
        ReallocatableArray<S32> n_disp_epj_tmp[n_thread];
        ReallocatableArray<S32> n_disp_spj_tmp[n_thread];
        // std::cerr<<"adr_epj_tmp[0].data_= "<<adr_epj_tmp[0].data()<<std::endl;
        // std::cerr<<"adr_spj_tmp[0].data_= "<<adr_spj_tmp[0].data()<<std::endl;
        // std::cerr<<"adr_ipg_tmp[0].data_= "<<adr_ipg_tmp[0].data()<<std::endl;
        // std::cerr<<"n_disp_epj_tmp[0].data_= "<<n_disp_epj_tmp[0].data()<<std::endl;
        // std::cerr<<"n_disp_spj_tmp[0].data_= "<<n_disp_spj_tmp[0].data()<<std::endl;
        // const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
        // std::cerr<<"adr_ipg_first= "<<adr_ipg_first<<" adr_ipg_last= "<<adr_ipg_last<<std::endl;
        PS_OMP(omp parallel num_threads(n_thread)) {
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
            PARTICLE_SIMULATOR_PRINT_MSG("- start OMP scope", 0);
#endif
            const S32 ith = Comm::getThreadNum();
            adr_epj_tmp[ith].clearSize();
            adr_spj_tmp[ith].clearSize();
            adr_ipg_tmp[ith].clearSize();
            n_disp_epj_tmp[ith].clearSize();
            n_disp_spj_tmp[ith].clearSize();
            // std::cerr<<"adr_epj_tmp[ith].data_= "<<adr_epj_tmp[ith].data()<<std::endl;
            // std::cerr<<"adr_spj_tmp[ith].data_= "<<adr_spj_tmp[ith].data()<<std::endl;
            // std::cerr<<"adr_ipg_tmp[ith].data_= "<<adr_ipg_tmp[ith].data()<<std::endl;
            // std::cerr<<"n_disp_epj_tmp[ith].data_= "<<n_disp_epj_tmp[ith].data()<<std::endl;
            // std::cerr<<"n_disp_spj_tmp[ith].data_= "<<n_disp_spj_tmp[ith].data()<<std::endl;

            S32 n_ep_cum_prev = 0;
            S32 n_sp_cum_prev = 0;
            PS_OMP(omp for schedule(dynamic, 4))
            for (S32 i = 0; i < n_ipg_this_rnd; i++) {
                // std::cerr<<"A) n_disp_epj_tmp[ith].size()= "<<n_disp_epj_tmp[ith].size()<<std::endl;
                // std::cerr<<"i= "<<i<<" adr_ipg_first= "<<adr_ipg_first<<" adr_ipg_last= "<<adr_ipg_last<<std::endl;
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
                // std::cerr<<"B) n_disp_epj_tmp[ith].size()= "<<n_disp_epj_tmp[ith].size()<<std::endl;
                TargetBox<TSM> target_box;
                target_box.set(ipg[i + adr_ipg_first], icut, morton_key);
                if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
                    PS_COMPILE_TIME_MESSAGE("SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE @ MakeInteractionListIndexLong");
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                    PARTICLE_SIMULATOR_PRINT_MSG("- before MakeListUsingTreeRecursiveTopPMM", 0);
#endif
                    MakeListUsingTreeRecursiveTopPMM<TSM, Ttc, Ttp, Tep, Tsp, Tmorton_key, Tchopleaf, Tchopnonleaf, Tcopyinfoclose,
                                                     CALC_DISTANCE_TYPE>(tc, adr_tc, tp, ep, adr_epj_tmp[ith], sp, adr_spj_tmp[ith], target_box,
                                                                         n_leaf_limit, lev_leaf_limit, adr_tree_sp_first, len_peri, theta,
                                                                         lev_pm_cell, morton_key);
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                    PARTICLE_SIMULATOR_PRINT_MSG("- after MakeListUsingTreeRecursiveTopPMM", 0);
#endif
                } else {
                    PS_COMPILE_TIME_MESSAGE("other than SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE @ MakeInteractionListIndexLong");
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                    PARTICLE_SIMULATOR_PRINT_MSG("- before MakeListUsingTreeRecursiveTop", 0);
                    target_box.dump(std::cout);
#endif
                    MakeListUsingTreeRecursiveTop<TSM, Ttc, Ttp, Tep, Tsp, Tchopleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>(
                        tc, adr_tc, tp, ep, adr_epj_tmp[ith], sp, adr_spj_tmp[ith], target_box, n_leaf_limit, adr_tree_sp_first, len_peri, theta);
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                    PARTICLE_SIMULATOR_PRINT_MSG("- after MakeListUsingTreeRecursiveTop", 0);
#endif
                }
                // std::cerr<<"C) n_disp_epj_tmp[ith].size()= "<<n_disp_epj_tmp[ith].size()<<std::endl;
                interaction_list.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                interaction_list.n_sp_[i] = adr_spj_tmp[ith].size() - n_sp_cum_prev;
                n_ep_cum_prev = adr_epj_tmp[ith].size();
                n_sp_cum_prev = adr_spj_tmp[ith].size();
            }  // end of OMP for
            // std::cerr<<"ith= "<<ith<<" n_ep_cum_prev= "<<n_ep_cum_prev<<" n_sp_cum_prev= "<<n_sp_cum_prev<<std::endl;

            // std::cerr<<"D) n_disp_epj_tmp[ith].size()= "<<n_disp_epj_tmp[ith].size()<<std::endl;
            n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
            n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
            // std::cerr<<"E) n_disp_epj_tmp[ith].size()= "<<n_disp_epj_tmp[ith].size()<<std::endl;
            // std::cerr<<"E2) n_disp_epj_tmp[0].size()= "<<n_disp_epj_tmp[0].size()<<std::endl;
        }  // end of OMP
        // std::cerr<<"F) n_disp_epj_tmp[0].size()= "<<n_disp_epj_tmp[0].size()<<std::endl;
        interaction_list.n_disp_ep_[0] = 0;
        interaction_list.n_disp_sp_[0] = 0;
        // for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
        for (S32 i = 0; i < n_ipg_this_rnd; i++) {
            interaction_list.n_disp_ep_[i + 1] = interaction_list.n_disp_ep_[i] + interaction_list.n_ep_[i];
            interaction_list.n_disp_sp_[i + 1] = interaction_list.n_disp_sp_[i] + interaction_list.n_sp_[i];
        }
        // std::cerr<<"G) n_disp_epj_tmp[0].size()= "<<n_disp_epj_tmp[0].size()<<std::endl;
        interaction_list.adr_ep_.resizeNoInitialize(interaction_list.n_disp_ep_[n_ipg_this_rnd] + interaction_list.n_ep_[n_ipg_this_rnd - 1]);
        interaction_list.adr_sp_.resizeNoInitialize(interaction_list.n_disp_sp_[n_ipg_this_rnd] + interaction_list.n_sp_[n_ipg_this_rnd - 1]);
        // std::cerr<<"H) n_disp_epj_tmp[0].size()= "<<n_disp_epj_tmp[0].size()<<std::endl;
        PS_OMP(omp parallel for num_threads(n_thread))
        for (S32 i = 0; i < n_thread; i++) {
            // std::cerr<<"ith= "<<i<<" adr_ipg_tmp[i].size()= "<<adr_ipg_tmp[i].size()<<std::endl;
            for (S32 j = 0; j < adr_ipg_tmp[i].size(); j++) {
                // std::cerr<<"i= "<<i<<" j= "<<j<<" adr_ipg_tmp[i][j]= "<<adr_ipg_tmp[i][j]<<std::endl;
                const S32 adr_ipg = adr_ipg_tmp[i][j];
                S32 adr_ep = interaction_list.n_disp_ep_[adr_ipg];
                // std::cerr<<"adr_ep= "<<adr_ep<<std::endl;
                const S32 k_ep_h = n_disp_epj_tmp[i][j];
                const S32 k_ep_e = n_disp_epj_tmp[i][j + 1];
                // std::cerr<<"k_ep_h= "<<k_ep_h<<" k_ep_e= "<<k_ep_e<<std::endl;
                for (S32 k = k_ep_h; k < k_ep_e; k++, adr_ep++) {
                    interaction_list.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                }
                // std::cerr<<"check 1"<<std::endl;
                S32 adr_sp = interaction_list.n_disp_sp_[adr_ipg];
                // std::cerr<<"adr_sp= "<<adr_sp<<std::endl;
                const S32 k_sp_h = n_disp_spj_tmp[i][j];
                const S32 k_sp_e = n_disp_spj_tmp[i][j + 1];
                // std::cerr<<"k_sp_h= "<<k_sp_h<<" k_sp_e= "<<k_sp_e<<std::endl;
                for (S32 k = k_sp_h; k < k_sp_e; k++, adr_sp++) {
                    interaction_list.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
                }
                // std::cerr<<"check 2"<<std::endl;
            }
        }
        // std::cerr<<"I) n_disp_epj_tmp[0].size()= "<<n_disp_epj_tmp[0].size()<<std::endl;

        for (int i = 0; i < n_thread; i++) {
            /*
            std::cerr << "adr_epj_tmp[" << i << "].size()=" << adr_epj_tmp[i].size() << std::endl;
            std::cerr << "adr_spj_tmp[" << i << "].size()=" << adr_spj_tmp[i].size() << std::endl;
            std::cerr << "adr_ipg_tmp[" << i << "].size()=" << adr_ipg_tmp[i].size() << std::endl;
            std::cerr << "n_disp_epj_tmp[" << i << "].size()=" << n_disp_epj_tmp[i].size() << std::endl;
            std::cerr << "n_disp_spj_tmp[" << i << "].size()=" << n_disp_spj_tmp[i].size() << std::endl;
            std::cerr << " adr_spj_tmp[" << i << "].capacity()=" << adr_spj_tmp[i].capacity() << std::endl;
            std::cerr << " adr_epj_tmp[" << i << "].capacity()=" << adr_epj_tmp[i].capacity() << std::endl;
            std::cerr << " adr_ipg_tmp[" << i << "].capacity()=" << adr_ipg_tmp[i].capacity() << std::endl;
            std::cerr << " n_disp_epj_tmp[" << i << "].capacity()=" << n_disp_epj_tmp[i].capacity() << std::endl;
            std::cerr << " n_disp_spj_tmp[" << i << "].capacity()=" << n_disp_spj_tmp[i].capacity() << std::endl;
            std::cerr << "n_disp_epj_tmp[" << i << "].getPointer()=" << n_disp_epj_tmp[i].getPointer() << std::endl;
            std::cerr << "n_disp_spj_tmp[" << i << "].getPointer()=" << n_disp_spj_tmp[i].getPointer() << std::endl;
            std::cerr << "adr_epj_tmp[" << i << "].getPointer()=" << adr_epj_tmp[i].getPointer() << std::endl;
            std::cerr << "adr_spj_tmp[" << i << "].getPointer()=" << adr_spj_tmp[i].getPointer() << std::endl;
            std::cerr << "adr_ipg_tmp[" << i << "].getPointer()=" << adr_ipg_tmp[i].getPointer() << std::endl;
            std::cerr<<" adr_epj_tmp["<<i<<"].data_= "<<adr_epj_tmp[i].data_<<std::endl;
            std::cerr<<" adr_spj_tmp["<<i<<"].data_= "<<adr_spj_tmp[i].data_<<std::endl;
            std::cerr<<" adr_ipg_tmp["<<i<<"].data_= "<<adr_ipg_tmp[i].data_<<std::endl;
            std::cerr<<" n_disp_epj_tmp["<<i<<"].data_= "<<n_disp_epj_tmp[i].data_<<std::endl;
            std::cerr<<" n_disp_spj_tmp["<<i<<"].data_= "<<n_disp_spj_tmp[i].data_<<std::endl;
            */
            n_disp_epj_tmp[i].freeMem();
            n_disp_spj_tmp[i].freeMem();
            adr_epj_tmp[i].freeMem();
            adr_spj_tmp[i].freeMem();
            adr_ipg_tmp[i].freeMem();
            // n_disp_spj_tmp[i].freeMem();
            // delete [] n_disp_spj_tmp[i].data_;
            // std::cerr << "check 5" << std::endl;

            // adr_epj_tmp[i].freeMem();
            // delete [] n_disp_epj_tmp[i].data_;
            // std::cerr << "check 1" << std::endl;

            // adr_spj_tmp[i].freeMem();
            // delete [] adr_spj_tmp[i].data_;
            // std::cerr << "check 2" << std::endl;

            // adr_ipg_tmp[i].freeMem();
            // delete [] adr_ipg_tmp[i].data_;
            // std::cerr << "check 3" << std::endl;

            // n_disp_epj_tmp[i].freeMem();
            // delete [] n_disp_epj_tmp[i].data_;
            // std::cerr << "check 4" << std::endl;
        }
    }
}
#else
template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tipg, class Tmorton_key, class Tchopleaf, class Tchopnonleaf,
          class Tcopyinfoclose, class Tinteraction_list, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void MakeInteractionListIndexLong(const ReallocatableArray<Ttc> &tc, const S32 adr_tc, const ReallocatableArray<Ttp> &tp, const S32 adr_tree_sp_first,
                                  const ReallocatableArray<Tep> &ep, const ReallocatableArray<Tsp> &sp, const ReallocatableArray<Tipg> &ipg,
                                  const S32 adr_ipg_first, const S32 n_ipg_this_rnd, const S32 n_leaf_limit, const S32 lev_leaf_limit,
                                  const F64 theta, const S32 icut, const S32 lev_pm_cell, const Tmorton_key &morton_key,
                                  Tinteraction_list &interaction_list, const F64vec len_peri) {
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START MakeInteractionListIndexLong *****", 0);
#endif
    const S32 n_ipg = ipg.size();
    const auto n_thread = Comm::getNumberOfThread();
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
    PS_PRINT_VARIABLE2(n_ipg, n_thread);
#endif
    if (n_ipg > 0) {
        interaction_list.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list.n_disp_ep_.resizeNoInitialize(n_ipg + 1);
        interaction_list.n_sp_.resizeNoInitialize(n_ipg);
        interaction_list.n_disp_sp_.resizeNoInitialize(n_ipg + 1);
        if (n_ipg_this_rnd > 0) {
            assert(0 <= adr_ipg_first && adr_ipg_first < n_ipg);
            assert(n_ipg_this_rnd <= n_ipg);
            const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
            ReallocatableArray<S32> adr_epj_tmp[n_thread];
            ReallocatableArray<S32> adr_spj_tmp[n_thread];
            ReallocatableArray<S32> adr_ipg_tmp[n_thread];
            ReallocatableArray<S32> n_disp_epj_tmp[n_thread];
            ReallocatableArray<S32> n_disp_spj_tmp[n_thread];
            PS_OMP(omp parallel num_threads(n_thread)) {
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                PARTICLE_SIMULATOR_PRINT_MSG("- start OMP scope", 0);
#endif
                const S32 ith = Comm::getThreadNum();
                adr_epj_tmp[ith].clearSize();
                adr_spj_tmp[ith].clearSize();
                adr_ipg_tmp[ith].clearSize();
                n_disp_epj_tmp[ith].clearSize();
                n_disp_spj_tmp[ith].clearSize();
                S32 n_ep_cum_prev = 0;
                S32 n_sp_cum_prev = 0;
                PS_OMP(omp for schedule(dynamic, 4))
                for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
                    adr_ipg_tmp[ith].push_back(i);
                    n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                    n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
                    TargetBox<TSM> target_box;
                    target_box.set(ipg[i], icut, morton_key);
                    if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
                        PS_COMPILE_TIME_MESSAGE("SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE @ MakeInteractionListIndexLong");
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                        PARTICLE_SIMULATOR_PRINT_MSG("- before MakeListUsingTreeRecursiveTopPMM", 0);
#endif
                        MakeListUsingTreeRecursiveTopPMM<TSM, Ttc, Ttp, Tep, Tsp, Tmorton_key, Tchopleaf, Tchopnonleaf, Tcopyinfoclose,
                                                         CALC_DISTANCE_TYPE>(tc, adr_tc, tp, ep, adr_epj_tmp[ith], sp, adr_spj_tmp[ith], target_box,
                                                                             n_leaf_limit, lev_leaf_limit, adr_tree_sp_first, len_peri, theta,
                                                                             lev_pm_cell, morton_key);
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                        PARTICLE_SIMULATOR_PRINT_MSG("- after MakeListUsingTreeRecursiveTopPMM", 0);
#endif
                    } else {
                        PS_COMPILE_TIME_MESSAGE("other than SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE @ MakeInteractionListIndexLong");
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                        PARTICLE_SIMULATOR_PRINT_MSG("- before MakeListUsingTreeRecursiveTop", 0);
                        target_box.dump(std::cout);
#endif
                        MakeListUsingTreeRecursiveTop<TSM, Ttc, Ttp, Tep, Tsp, Tchopleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>(
                            tc, adr_tc, tp, ep, adr_epj_tmp[ith], sp, adr_spj_tmp[ith], target_box, n_leaf_limit, adr_tree_sp_first, len_peri, theta);
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                        PARTICLE_SIMULATOR_PRINT_MSG("- after MakeListUsingTreeRecursiveTop", 0);
#endif
                    }
                    interaction_list.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                    interaction_list.n_sp_[i] = adr_spj_tmp[ith].size() - n_sp_cum_prev;
                    n_ep_cum_prev = adr_epj_tmp[ith].size();
                    n_sp_cum_prev = adr_spj_tmp[ith].size();
                }  // end of OMP for
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
            }  // end of OMP

            interaction_list.n_disp_ep_[0] = 0;
            interaction_list.n_disp_sp_[0] = 0;
            for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
                interaction_list.n_disp_ep_[i + 1] = interaction_list.n_disp_ep_[i] + interaction_list.n_ep_[i];
                interaction_list.n_disp_sp_[i + 1] = interaction_list.n_disp_sp_[i] + interaction_list.n_sp_[i];
            }
            interaction_list.adr_ep_.resizeNoInitialize(interaction_list.n_disp_ep_[adr_ipg_last + 1]);
            interaction_list.adr_sp_.resizeNoInitialize(interaction_list.n_disp_sp_[adr_ipg_last + 1]);
            PS_OMP(omp parallel for num_threads(n_thread))
            for (S32 i = 0; i < n_thread; i++) {
                for (S32 j = 0; j < adr_ipg_tmp[i].size(); j++) {
                    const S32 adr_ipg = adr_ipg_tmp[i][j];
                    S32 adr_ep = interaction_list.n_disp_ep_[adr_ipg];
                    const S32 k_ep_h = n_disp_epj_tmp[i][j];
                    const S32 k_ep_e = n_disp_epj_tmp[i][j + 1];
                    for (S32 k = k_ep_h; k < k_ep_e; k++, adr_ep++) {
                        interaction_list.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                    }
                    S32 adr_sp = interaction_list.n_disp_sp_[adr_ipg];
                    const S32 k_sp_h = n_disp_spj_tmp[i][j];
                    const S32 k_sp_e = n_disp_spj_tmp[i][j + 1];
                    for (S32 k = k_sp_h; k < k_sp_e; k++, adr_sp++) {
                        interaction_list.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
                    }
                }
            }
            for (int i = 0; i < n_thread; i++) {
                adr_epj_tmp[i].freeMem();
                adr_spj_tmp[i].freeMem();
                adr_ipg_tmp[i].freeMem();
                n_disp_epj_tmp[i].freeMem();
                n_disp_spj_tmp[i].freeMem();
            }
        }
    }

#if defined(PS_DEBUG_PRINT_makeInteractionListIndexLong)
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE(interaction_list.adr_ep_.size());
#if PS_DEBUG_PRINT_makeInteractionListIndexLong >= 2
        for (auto i = 0; i < interaction_list.adr_ep_.size(); i++) {
            PS_PRINT_VARIABLE2(i, interaction_list.adr_ep_[i]);
        }
#endif
        PS_PRINT_VARIABLE(interaction_list.adr_sp_.size());
#if PS_DEBUG_PRINT_makeInteractionListIndexLong >= 2
        for (auto i = 0; i < interaction_list.adr_sp_.size(); i++) {
            PS_PRINT_VARIABLE2(i, interaction_list.adr_sp_[i]);
        }
#endif
    }
    PARTICLE_SIMULATOR_PRINT_MSG("***** END makeInteractionListIndexLong *****", 0);
#endif
}
#endif

#if 1
/*
template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tipg, class Tmorton_key, class Tchopleaf, class Tchopnonleaf,
          class Tcopyinfoclose, class Tinteraction_list, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void MakeInteractionListIndexShort(const ReallocatableArray<Ttc> &tc, const S32 adr_tc, const ReallocatableArray<Ttp> &tp,
                                   const S32 adr_tree_sp_first, const ReallocatableArray<Tep> &ep, const ReallocatableArray<Tsp> &sp,
                                   const ReallocatableArray<Tipg> &ipg, const S32 adr_ipg_first, const S32 n_ipg_this_rnd, const S32 n_leaf_limit,
                                   const S32 lev_leaf_limit, const F64 theta, const S32 icut, const S32 lev_pm_cell, const Tmorton_key &morton_key,
                                   Tinteraction_list &interaction_list, const F64vec len_peri) {
*/
template <class TSM, class Ttc, class Ttp, class Tep, class Tipg, class Tmorton_key, class Tchopleaf, class Tchopnonleaf, class Tcopyinfoclose,
          class Tinteraction_list, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void MakeInteractionListIndexShort(const ReallocatableArray<Ttc> &tc, const S32 adr_tc, const ReallocatableArray<Ttp> &tp,
                                   const ReallocatableArray<Tep> &ep, const ReallocatableArray<Tipg> &ipg, const S32 adr_ipg_first,
                                   const S32 n_ipg_this_rnd, const S32 n_leaf_limit, const S32 lev_leaf_limit, const F64 theta,
                                   const Tmorton_key &morton_key, Tinteraction_list &interaction_list) {
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START MakeInteractionListIndexLong *****", 0);
#endif
    constexpr S32 icut = 0;
    const S32 n_ipg = ipg.size();
    const auto n_thread = Comm::getNumberOfThread();
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
    PS_PRINT_VARIABLE2(n_ipg, n_thread);
#endif
    if (n_ipg_this_rnd > 0) {
        assert(0 <= adr_ipg_first && adr_ipg_first < n_ipg);
        assert(n_ipg_this_rnd <= n_ipg);
        interaction_list.n_ep_.resizeNoInitialize(n_ipg_this_rnd);
        interaction_list.n_disp_ep_.resizeNoInitialize(n_ipg_this_rnd + 1);
        //interaction_list.n_sp_.resizeNoInitialize(n_ipg_this_rnd);
        //interaction_list.n_disp_sp_.resizeNoInitialize(n_ipg_this_rnd + 1);

        ReallocatableArray<S32> adr_epj_tmp[n_thread];
        //ReallocatableArray<S32> adr_spj_tmp[n_thread];
        ReallocatableArray<S32> adr_ipg_tmp[n_thread];
        ReallocatableArray<S32> n_disp_epj_tmp[n_thread];
        //ReallocatableArray<S32> n_disp_spj_tmp[n_thread];
        PS_OMP(omp parallel num_threads(n_thread)) {
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
            PARTICLE_SIMULATOR_PRINT_MSG("- start OMP scope", 0);
#endif
            const S32 ith = Comm::getThreadNum();
            adr_epj_tmp[ith].clearSize();
            //adr_spj_tmp[ith].clearSize();
            adr_ipg_tmp[ith].clearSize();
            n_disp_epj_tmp[ith].clearSize();
            //n_disp_spj_tmp[ith].clearSize();

            S32 n_ep_cum_prev = 0;
            //S32 n_sp_cum_prev = 0;
            PS_OMP(omp for schedule(dynamic, 4))
            for (S32 i = 0; i < n_ipg_this_rnd; i++) {
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                //n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
                TargetBox<TSM> target_box;
                target_box.set(ipg[i + adr_ipg_first], icut, morton_key);

                MakeListUsingTreeShortRecursiveTop<TSM, Ttc, Ttp, Tep, Tchopleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>(
                    tc, adr_tc, tp, ep, adr_epj_tmp[ith], target_box, n_leaf_limit, F64vec(0.0));

                /*
                if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
                    PS_COMPILE_TIME_MESSAGE("SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE @ MakeInteractionListIndexLong");
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                    PARTICLE_SIMULATOR_PRINT_MSG("- before MakeListUsingTreeRecursiveTopPMM", 0);
#endif
                    MakeListUsingTreeRecursiveTopPMM<TSM, Ttc, Ttp, Tep, Tsp, Tmorton_key, Tchopleaf, Tchopnonleaf, Tcopyinfoclose,
                                                     CALC_DISTANCE_TYPE>(tc, adr_tc, tp, ep, adr_epj_tmp[ith], sp, adr_spj_tmp[ith], target_box,
                                                                         n_leaf_limit, lev_leaf_limit, adr_tree_sp_first, len_peri, theta,
                                                                         lev_pm_cell, morton_key);
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                    PARTICLE_SIMULATOR_PRINT_MSG("- after MakeListUsingTreeRecursiveTopPMM", 0);
#endif
                } else {
                    PS_COMPILE_TIME_MESSAGE("other than SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE @ MakeInteractionListIndexLong");
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                    PARTICLE_SIMULATOR_PRINT_MSG("- before MakeListUsingTreeRecursiveTop", 0);
                    target_box.dump(std::cout);
#endif
                    MakeListUsingTreeRecursiveTop<TSM, Ttc, Ttp, Tep, Tsp, Tchopleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>(
                        tc, adr_tc, tp, ep, adr_epj_tmp[ith], sp, adr_spj_tmp[ith], target_box, n_leaf_limit, adr_tree_sp_first, len_peri, theta);
#if defined(PS_DEBUG_PRINT_MakeInteractionListIndexLong)
                    PARTICLE_SIMULATOR_PRINT_MSG("- after MakeListUsingTreeRecursiveTop", 0);
#endif
                }
                */

                interaction_list.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                //interaction_list.n_sp_[i] = adr_spj_tmp[ith].size() - n_sp_cum_prev;
                n_ep_cum_prev = adr_epj_tmp[ith].size();
                //n_sp_cum_prev = adr_spj_tmp[ith].size();
            }  // end of OMP for
            n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
            //n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
        }  // end of OMP
        interaction_list.n_disp_ep_[0] = 0;
        //interaction_list.n_disp_sp_[0] = 0;
        for (S32 i = 0; i < n_ipg_this_rnd; i++) {
            interaction_list.n_disp_ep_[i + 1] = interaction_list.n_disp_ep_[i] + interaction_list.n_ep_[i];
            //interaction_list.n_disp_sp_[i + 1] = interaction_list.n_disp_sp_[i] + interaction_list.n_sp_[i];
        }
        interaction_list.adr_ep_.resizeNoInitialize(interaction_list.n_disp_ep_[n_ipg_this_rnd] + interaction_list.n_ep_[n_ipg_this_rnd - 1]);
        //interaction_list.adr_sp_.resizeNoInitialize(interaction_list.n_disp_sp_[n_ipg_this_rnd] + interaction_list.n_sp_[n_ipg_this_rnd - 1]);
        PS_OMP(omp parallel for num_threads(n_thread))
        for (S32 i = 0; i < n_thread; i++) {
            for (S32 j = 0; j < adr_ipg_tmp[i].size(); j++) {
                const S32 adr_ipg = adr_ipg_tmp[i][j];
                S32 adr_ep = interaction_list.n_disp_ep_[adr_ipg];
                const S32 k_ep_h = n_disp_epj_tmp[i][j];
                const S32 k_ep_e = n_disp_epj_tmp[i][j + 1];
                for (S32 k = k_ep_h; k < k_ep_e; k++, adr_ep++) {
                    interaction_list.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                }
                /*
                S32 adr_sp = interaction_list.n_disp_sp_[adr_ipg];
                const S32 k_sp_h = n_disp_spj_tmp[i][j];
                const S32 k_sp_e = n_disp_spj_tmp[i][j + 1];
                for (S32 k = k_sp_h; k < k_sp_e; k++, adr_sp++) {
                    interaction_list.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
                }
                    */
            }
        }

        for (int i = 0; i < n_thread; i++) {
            n_disp_epj_tmp[i].freeMem();
            //n_disp_spj_tmp[i].freeMem();
            adr_epj_tmp[i].freeMem();
            //adr_spj_tmp[i].freeMem();
            adr_ipg_tmp[i].freeMem();
        }
    }
}
#else
// short
template <class TSM, class Ttc, class Ttp, class Tep, class Tipg, class Tmorton_key, class Tchopleaf, class Tchopnonleaf, class Tcopyinfoclose,
          class Tinteraction_list, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void MakeInteractionListIndexShort(const ReallocatableArray<Ttc> &tc, const S32 adr_tc, const ReallocatableArray<Ttp> &tp,
                                   const ReallocatableArray<Tep> &ep, const ReallocatableArray<Tipg> &ipg, const S32 adr_ipg_first,
                                   const S32 n_ipg_this_rnd, const S32 n_leaf_limit, const S32 lev_leaf_limit, const F64 theta,
                                   const Tmorton_key &morton_key, Tinteraction_list &interaction_list) {
#if defined(PS_DEBUG_PRINT_makeInteractionListIndexShort)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START makeInteractionListIndexShort *****", 0);
#endif
    const S32 n_ipg = ipg.size();
    const auto n_thread = Comm::getNumberOfThread();
    ReallocatableArray<SuperParticleBase> spj_dummy;
    ReallocatableArray<S32> adr_spj;
    constexpr S32 icut = 0;
    if (n_ipg > 0) {
        interaction_list.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list.n_disp_ep_.resizeNoInitialize(n_ipg + 1);
        if (n_ipg_this_rnd > 0) {
            assert(0 <= adr_ipg_first && adr_ipg_first < n_ipg);
            assert(n_ipg_this_rnd <= n_ipg);
            const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
            // std::vector<ReallocatableArray<S32>> adr_epj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
            // std::vector<ReallocatableArray<S32>> adr_ipg_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
            // std::vector<ReallocatableArray<S32>> n_disp_epj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
            ReallocatableArray<S32> adr_epj_tmp[n_thread];
            ReallocatableArray<S32> adr_ipg_tmp[n_thread];
            ReallocatableArray<S32> n_disp_epj_tmp[n_thread];
            PS_OMP_PARALLEL {
                const S32 ith = Comm::getThreadNum();
                adr_epj_tmp[ith].clearSize();
                adr_ipg_tmp[ith].clearSize();
                n_disp_epj_tmp[ith].clearSize();
                S32 n_ep_cum_prev = 0;
                PS_OMP_PARALLEL_FOR_D4
                for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
                    adr_ipg_tmp[ith].push_back(i);
                    n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                    TargetBox<TSM> target_box;
                    target_box.set(ipg[i], icut, morton_key);

                    MakeListUsingTreeShortRecursiveTop<TSM, Ttc, Ttp, Tep, Tchopleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>(
                        tc, adr_tc, tp, ep, adr_epj_tmp[ith], target_box, n_leaf_limit, F64vec(0.0));

                    interaction_list.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                    n_ep_cum_prev = adr_epj_tmp[ith].size();
                }
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
            }  // end of OMP
            interaction_list.n_disp_ep_[0] = 0;
            for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
                interaction_list.n_disp_ep_[i + 1] = interaction_list.n_disp_ep_[i] + interaction_list.n_ep_[i];
            }
            interaction_list.adr_ep_.resizeNoInitialize(interaction_list.n_disp_ep_[adr_ipg_last + 1]);
            PS_OMP_PARALLEL
            for (S32 i = 0; i < n_thread; i++) {
                for (S32 j = 0; j < adr_ipg_tmp[i].size(); j++) {
                    const S32 adr_ipg = adr_ipg_tmp[i][j];
                    S32 adr_ep = interaction_list.n_disp_ep_[adr_ipg];
                    const S32 k_ep_h = n_disp_epj_tmp[i][j];
                    const S32 k_ep_e = n_disp_epj_tmp[i][j + 1];
                    for (S32 k = k_ep_h; k < k_ep_e; k++, adr_ep++) {
                        interaction_list.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                    }
                }
            }
            for (int i = 0; i < n_thread; i++) {
                adr_epj_tmp[i].freeMem();
                adr_ipg_tmp[i].freeMem();
                n_disp_epj_tmp[i].freeMem();
            }
        }
    }
}  // end of MakeInteractionListIndexShort
#endif

}  // namespace ParticleSimulator