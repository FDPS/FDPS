#pragma once
// In this file, function templates used in TreeForForce::calcForceNoWalk are implemeted.

namespace ParticleSimulator {

// for short
#if 1
// PMMM
template <class Tipg, class Tinteraction_list, class Tepi, class Tepj, class Tforce, class Tfunc_ep_ep>
void CalcForceNoWalkPMM(const ReallocatableArray<Tipg>& ipg, const Tinteraction_list& interaction_list, const S32 n_loc_tot,
                        ReallocatableArray<Tepi>& epi_sorted, const ReallocatableArray<Tepj>& epj_sorted, ReallocatableArray<Tforce>& force_sorted,
                        const S32 n_group_limit, Tfunc_ep_ep pfunc_ep_ep, const U32 pop_id, const S32 adr_ipg_first, const S32 n_ipg_this_rnd,
                        const bool clear_force, const bool update_epi, std::vector<CountT>& n_interaction_ep_ep_local, TimeProfile& time_profile) {
    /*
    // Check the arguments
    assert(0 <= pop_id && pop_id <= N_POP);
    assert(0 <= adr_ipg_first);
    assert(0 <= n_ipg_this_rnd);
    const S32 n_thread_max = Comm::getMaxThreads();
    // Allocate local buffers
    static bool is_first_call {true};
    static ReallocatableArray<S32> * adr_epi_for_force {nullptr};
    static ReallocatableArray<Tepi> * epi_for_force {nullptr};
    static ReallocatableArray<Tepj> * epj_for_force {nullptr};
    static ReallocatableArray<Tforce> * force_for_force {nullptr};
    if (is_first_call) {
        adr_epi_for_force = new ReallocatableArray<S32>[n_thread_max];
        epi_for_force = new ReallocatableArray<Tepi>[n_thread_max];
        epj_for_force = new ReallocatableArray<Tepj>[n_thread_max];
        force_for_force = new ReallocatableArray<Tforce>[n_thread_max];
        is_first_call = false;
    }
    force_sorted.resizeNoInitialize(n_loc_tot);
    const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
    const S64 n_ipg = ipg.size();
    assert(adr_ipg_last < n_ipg);
    if (n_ipg_this_rnd > 0) {
        if (!(0 <= pop_id && pop_id < N_POP)) {
            if (clear_force) {
PS_OMP(omp parallel for num_threads(n_thread))
                for(S32 i=0; i<n_loc_tot; i++){
                    force_sorted[i].clear();
                }
            }
//PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 4))
PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 1))
            for(S32 i=adr_ipg_first; i<=adr_ipg_last; i++){
                const S32 ith = Comm::getThreadNum();
                const S32 n_epi = ipg[i].n_ptcl_;
                const S32 adr_epi_head = ipg[i].adr_ptcl_;
                const S32 n_epj = interaction_list.n_ep_[i];
                const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
                const S32 adr_epj_end  = interaction_list.n_disp_ep_[i+1];
                epj_for_force[ith].resizeNoInitialize(n_epj);
                S32 n_cnt = 0;
                for(S32 j=adr_epj_head; j<adr_epj_end; j++, n_cnt++){
                    const S32 adr_epj = interaction_list.adr_ep_[j];
                    epj_for_force[ith][n_cnt] = epj_sorted[adr_epj];
                }
                pfunc_ep_ep(epi_sorted.getPointer(adr_epi_head),     n_epi,
                            epj_for_force[ith].getPointer(),   n_epj,
                            force_sorted.getPointer(adr_epi_head));
                n_epi_processed_local[ith] += n_epi;
                n_epj_processed_local[ith] += n_epj;
                n_interaction_ep_ep_local[ith] += n_epi * n_epj;
            }
        } else {
//PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 4))
PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 1))
            for(S32 i=adr_ipg_first; i<=adr_ipg_last; i++){
                const S32 ith = Comm::getThreadNum();
                const S32 offset = ipg[i].adr_ptcl_;
                const S32 n_ptcl = ipg[i].n_ptcl_;
                adr_epi_for_force[ith].clearSize();
                adr_epi_for_force[ith].reserveEmptyAreaAtLeast(n_group_limit);
                const S32 n_tail = offset + n_ptcl;
                for (S32 k=offset; k<n_tail; k++) {
                    if (GetMyPopId(epi_sorted[k]) == pop_id) {
                        adr_epi_for_force[ith].push_back(k);
                    }
                }
                const S32 n_epi = adr_epi_for_force[ith].size();
                epi_for_force[ith].resizeNoInitialize(n_epi);
                for (S32 k=0; k<n_epi; k++) {
                    const S32 adr_epi = adr_epi_for_force[ith][k];
                    epi_for_force[ith][k] = epi_sorted[adr_epi];
                }
                if(clear_force){
                    for (S32 k=0; k<n_epi; k++) {
                        const S32 adr_epi = adr_epi_for_force[ith][k];
                        force_sorted[adr_epi].clear();
                    }
                }
                force_for_force[ith].resizeNoInitialize(n_epi);
                for (S32 k=0; k<n_epi; k++) {
                    const S32 adr_epi = adr_epi_for_force[ith][k];
                    force_for_force[ith][k] = force_sorted[adr_epi];
                }
                const S32 n_epj = interaction_list.n_ep_[i];
                const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
                const S32 adr_epj_end  = interaction_list.n_disp_ep_[i+1];
                epj_for_force[ith].resizeNoInitialize(n_epj);
                S32 n_cnt = 0;
                for(S32 k=adr_epj_head; k<adr_epj_end; k++, n_cnt++){
                    const S32 adr_epj = interaction_list.adr_ep_[k];
                    epj_for_force[ith][n_cnt] = epj_sorted[adr_epj];
                }
                pfunc_ep_ep(epi_for_force[ith].getPointer(),   n_epi,
                            epj_for_force[ith].getPointer(),   n_epj,
                            force_for_force[ith].getPointer());
                for (S32 k=0; k<n_epi; k++) {
                    const S32 adr_epi = adr_epi_for_force[ith][k];
                    force_sorted[adr_epi] = force_for_force[ith][k];
                }
                if (update_epi) {
                    for (S32 k=0; k<n_epi; k++) {
                        const S32 adr_epi = adr_epi_for_force[ith][k];
                        force_sorted[adr_epi] = force_for_force[ith][k];
                    }
                }
                n_epi_processed_local[ith] += n_epi;
                n_epj_processed_local[ith] += n_epj;
                n_interaction_ep_ep_local[ith] += n_epi * n_epj;
            }
        }
    }
    */
}
#else
// original
template <class Tipg, class Tinteraction_list, class Tepi, class Tepj, class Tforce, class Tfunc_ep_ep, U32 N_POP>
void CalcForceNoWalk(const ReallocatableArray<Tipg>& ipg, const Tinteraction_list& interaction_list, const S32 n_loc_tot,
                     ReallocatableArray<Tepi>& epi_sorted, const ReallocatableArray<Tepj>& epj_sorted, ReallocatableArray<Tforce>& force_sorted,
                     const S32 n_group_limit, Tfunc_ep_ep pfunc_ep_ep, const U32 pop_id, const S32 adr_ipg_first, const S32 n_ipg_this_rnd,
                     const bool clear_force, const bool update_epi, std::vector<CountT>& n_epi_processed_local,
                     std::vector<CountT>& n_epj_processed_local, std::vector<CountT>& n_interaction_ep_ep_local, TimeProfile& time_profile,
                     const S32 n_thread) {
    // Check the arguments
    assert(0 <= pop_id && pop_id <= N_POP);
    assert(0 <= adr_ipg_first);
    assert(0 <= n_ipg_this_rnd);
    const S32 n_thread_max = Comm::getMaxThreads();
    // Allocate local buffers
    static bool is_first_call{true};
    static ReallocatableArray<S32>* adr_epi_for_force{nullptr};
    static ReallocatableArray<Tepi>* epi_for_force{nullptr};
    static ReallocatableArray<Tepj>* epj_for_force{nullptr};
    static ReallocatableArray<Tforce>* force_for_force{nullptr};
    if (is_first_call) {
        adr_epi_for_force = new ReallocatableArray<S32>[n_thread_max];
        epi_for_force = new ReallocatableArray<Tepi>[n_thread_max];
        epj_for_force = new ReallocatableArray<Tepj>[n_thread_max];
        force_for_force = new ReallocatableArray<Tforce>[n_thread_max];
        is_first_call = false;
    }
    force_sorted.resizeNoInitialize(n_loc_tot);
    const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
    const S64 n_ipg = ipg.size();
    assert(adr_ipg_last < n_ipg);
    if (n_ipg_this_rnd > 0) {
        if (!(0 <= pop_id && pop_id < N_POP)) {
            if (clear_force) {
PS_OMP(omp parallel for num_threads(n_thread))
for (S32 i = 0; i < n_loc_tot; i++) {
    force_sorted[i].clear();
}
            }
            // PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 4))
PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 1))
for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
    const S32 ith = Comm::getThreadNum();
    const S32 n_epi = ipg[i].n_ptcl_;
    const S32 adr_epi_head = ipg[i].adr_ptcl_;
    const S32 n_epj = interaction_list.n_ep_[i];
    const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
    const S32 adr_epj_end = interaction_list.n_disp_ep_[i + 1];
    epj_for_force[ith].resizeNoInitialize(n_epj);
    S32 n_cnt = 0;
    for (S32 j = adr_epj_head; j < adr_epj_end; j++, n_cnt++) {
        const S32 adr_epj = interaction_list.adr_ep_[j];
        epj_for_force[ith][n_cnt] = epj_sorted[adr_epj];
    }
    pfunc_ep_ep(epi_sorted.getPointer(adr_epi_head), n_epi, epj_for_force[ith].getPointer(), n_epj, force_sorted.getPointer(adr_epi_head));
    n_epi_processed_local[ith] += n_epi;
    n_epj_processed_local[ith] += n_epj;
    n_interaction_ep_ep_local[ith] += n_epi * n_epj;
}
        } else {
            // PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 4))
PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 1))
for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
    const S32 ith = Comm::getThreadNum();
    const S32 offset = ipg[i].adr_ptcl_;
    const S32 n_ptcl = ipg[i].n_ptcl_;
    adr_epi_for_force[ith].clearSize();
    adr_epi_for_force[ith].reserveEmptyAreaAtLeast(n_group_limit);
    const S32 n_tail = offset + n_ptcl;
    for (S32 k = offset; k < n_tail; k++) {
        if (GetMyPopId(epi_sorted[k]) == pop_id) {
            adr_epi_for_force[ith].push_back(k);
        }
    }
    const S32 n_epi = adr_epi_for_force[ith].size();
    epi_for_force[ith].resizeNoInitialize(n_epi);
    for (S32 k = 0; k < n_epi; k++) {
        const S32 adr_epi = adr_epi_for_force[ith][k];
        epi_for_force[ith][k] = epi_sorted[adr_epi];
    }
    if (clear_force) {
        for (S32 k = 0; k < n_epi; k++) {
            const S32 adr_epi = adr_epi_for_force[ith][k];
            force_sorted[adr_epi].clear();
        }
    }
    force_for_force[ith].resizeNoInitialize(n_epi);
    for (S32 k = 0; k < n_epi; k++) {
        const S32 adr_epi = adr_epi_for_force[ith][k];
        force_for_force[ith][k] = force_sorted[adr_epi];
    }
    const S32 n_epj = interaction_list.n_ep_[i];
    const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
    const S32 adr_epj_end = interaction_list.n_disp_ep_[i + 1];
    epj_for_force[ith].resizeNoInitialize(n_epj);
    S32 n_cnt = 0;
    for (S32 k = adr_epj_head; k < adr_epj_end; k++, n_cnt++) {
        const S32 adr_epj = interaction_list.adr_ep_[k];
        epj_for_force[ith][n_cnt] = epj_sorted[adr_epj];
    }
    pfunc_ep_ep(epi_for_force[ith].getPointer(), n_epi, epj_for_force[ith].getPointer(), n_epj, force_for_force[ith].getPointer());
    for (S32 k = 0; k < n_epi; k++) {
        const S32 adr_epi = adr_epi_for_force[ith][k];
        force_sorted[adr_epi] = force_for_force[ith][k];
    }
    if (update_epi) {
        for (S32 k = 0; k < n_epi; k++) {
            const S32 adr_epi = adr_epi_for_force[ith][k];
            force_sorted[adr_epi] = force_for_force[ith][k];
        }
    }
    n_epi_processed_local[ith] += n_epi;
    n_epj_processed_local[ith] += n_epj;
    n_interaction_ep_ep_local[ith] += n_epi * n_epj;
}
        }
    }
}
#endif

// for long
#if 1
template <class Ttc, class Tipg, class Tinteraction_list, class Tepi, class Tepj, class Tspj, class Tforce, class Tfunc_ep_ep, class Tfunc_ep_sp>
void CalcForceNoWalkPMM(const ReallocatableArray<Ttc>& tc, const ReallocatableArray<Tipg>& ipg, const Tinteraction_list& interaction_list,
                        const S32 n_loc_tot, const S32 n_let_sp, ReallocatableArray<Tepi>& epi_sorted, const ReallocatableArray<Tepj>& epj_sorted,
                        const ReallocatableArray<Tspj>& spj_sorted, ReallocatableArray<Tforce>& force_sorted, const S32 n_group_limit,
                        Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const S32 adr_ipg_first, const S32 n_ipg_this_rnd, const bool clear_force,
                        const bool update_epi, CountT& n_interaction_ep_ep_local, CountT& n_interaction_ep_sp_local, TimeProfile& time_profile) {
    // Check the arguments
    // assert(0 <= pop_id && pop_id <= N_POP);
    // assert(0 <= adr_ipg_first);
    // assert(0 <= n_ipg_this_rnd);
    // const S32 n_thread_max = Comm::getMaxThreads();
    const S32 n_thread = Comm::getNumberOfThread();
    // Allocate local buffers
    ReallocatableArray<Tepj> epj_for_force[n_thread];
    ReallocatableArray<Tspj> spj_for_force[n_thread];
    force_sorted.resizeNoInitialize(n_loc_tot);
    const S64 n_ipg = ipg.size();
    //std::cerr << "n_ipg= " << n_ipg << std::endl;
    const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
    //std::cerr << "n_ipg_this_rnd=" << n_ipg_this_rnd << std::endl;
    //std::cerr << "adr_ipg_first=" << adr_ipg_first << std::endl;
    //std::cerr << "adr_ipg_last=" << adr_ipg_last << std::endl;
    assert(adr_ipg_last < n_ipg);
    if (n_ipg_this_rnd > 0) {
        if (clear_force) {
            auto adr_ipg_tail = adr_ipg_first + n_ipg_this_rnd - 1;
            auto adr_epi_head = ipg[adr_ipg_first].adr_ptcl_;
            auto adr_epi_end  = ipg[adr_ipg_tail].adr_ptcl_ + ipg[adr_ipg_tail].n_ptcl_;
            //std::cerr<<"adr_epi_head="<<adr_epi_head<<" adr_epi_end="<<adr_epi_end<<std::endl;
            PS_OMP_PARALLEL_FOR
            for (S32 i = adr_epi_head; i < adr_epi_end; i++) {
                force_sorted[i].clear();
            }
        }
        PS_OMP_PARALLEL_FOR_D
        // for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
        for (S32 i = 0; i < n_ipg_this_rnd; i++) {
            //std::cerr << "i= " << i << " adr_ipg_first=" << adr_ipg_first << " adr_ipg_last=" << adr_ipg_last << std::endl;
            const S32 ith = Comm::getThreadNum();
            const S32 adr_ipg = i + adr_ipg_first;
            const S32 n_epi = ipg[adr_ipg].n_ptcl_;
            const S32 adr_epi_head = ipg[adr_ipg].adr_ptcl_;
            const S32 n_epj = interaction_list.n_ep_[i];
            const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
            const S32 adr_epj_end = interaction_list.n_disp_ep_[i + 1];
            const S32 n_spj = interaction_list.n_sp_[i];
            const S32 adr_spj_head = interaction_list.n_disp_sp_[i];
            const S32 adr_spj_end = interaction_list.n_disp_sp_[i + 1];
            epj_for_force[ith].resizeNoInitialize(n_epj);
            spj_for_force[ith].resizeNoInitialize(n_spj);
            S32 n_ep_cnt = 0;
            for (S32 j = adr_epj_head; j < adr_epj_end; j++, n_ep_cnt++) {
                const S32 adr_epj = interaction_list.adr_ep_[j];
                epj_for_force[ith][n_ep_cnt] = epj_sorted[adr_epj];
            }
            //std::cerr << "n_epi=" << n_epi << " n_epj=" << n_epj << " n_spj=" << n_spj << std::endl;
            //std::cerr << "adr_epi_head=" << adr_epi_head << std::endl;
            //std::cerr << "n_epj=" << n_epj << std::endl;
            pfunc_ep_ep(epi_sorted.getPointer(adr_epi_head), n_epi, epj_for_force[ith].getPointer(), n_epj, force_sorted.getPointer(adr_epi_head));
            //std::cerr << "check 1" << std::endl;
            S32 n_sp_cnt = 0;
            for (S32 j = adr_spj_head; j < adr_spj_end; j++, n_sp_cnt++) {
                const S32 adr_spj = interaction_list.adr_sp_[j];
                if (adr_spj < n_let_sp) {
                    spj_for_force[ith][n_sp_cnt] = spj_sorted[adr_spj];
                } else {
                    spj_for_force[ith][n_sp_cnt].copyFromMoment(tc[adr_spj - n_let_sp].mom_);
                }
            }
            //std::cerr << "check 2" << std::endl;
            pfunc_ep_sp(epi_sorted.getPointer(adr_epi_head), n_epi, spj_for_force[ith].getPointer(), n_spj, force_sorted.getPointer(adr_epi_head));
            //std::cerr << "check 3" << std::endl;
            n_interaction_ep_ep_local += n_epi * n_epj;
            n_interaction_ep_sp_local += n_epi * n_spj;
            //std::cerr << "n_interaction_ep_ep_local=" << n_interaction_ep_ep_local << std::endl;
            //std::cerr << "n_interaction_ep_sp_local=" << n_interaction_ep_sp_local << std::endl;
        }
    }
}


template <typename Ttc, typename Tipg, typename Tinteraction_list, typename Tepi, typename Tepj, typename Tforce, typename Tfunc_ep_ep>
void CalcForceNoWalkShort(const ReallocatableArray<Ttc>& tc, const ReallocatableArray<Tipg>& ipg, const Tinteraction_list& interaction_list,
                            const S32 n_loc_tot, ReallocatableArray<Tepi>& epi_sorted, const ReallocatableArray<Tepj>& epj_sorted,
                            ReallocatableArray<Tforce>& force_sorted, const S32 n_group_limit,
                            Tfunc_ep_ep pfunc_ep_ep, const S32 adr_ipg_first, const S32 n_ipg_this_rnd, const bool clear_force,
                            const bool update_epi, CountT& n_interaction_ep_ep_local, TimeProfile& time_profile) {
    // Check the arguments
    // assert(0 <= pop_id && pop_id <= N_POP);
    // assert(0 <= adr_ipg_first);
    // assert(0 <= n_ipg_this_rnd);
    // const S32 n_thread_max = Comm::getMaxThreads();
    const S32 n_thread = Comm::getNumberOfThread();
    // Allocate local buffers
    ReallocatableArray<Tepj> epj_for_force[n_thread];
    force_sorted.resizeNoInitialize(n_loc_tot);
    const S64 n_ipg = ipg.size();
    const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
    assert(adr_ipg_last < n_ipg);
    if (n_ipg_this_rnd > 0) {
        if (clear_force) {
            auto adr_ipg_tail = adr_ipg_first + n_ipg_this_rnd - 1;
            auto adr_epi_head = ipg[adr_ipg_first].adr_ptcl_;
            auto adr_epi_end  = ipg[adr_ipg_tail].adr_ptcl_ + ipg[adr_ipg_tail].n_ptcl_;
            PS_OMP_PARALLEL_FOR
            for (S32 i = adr_epi_head; i < adr_epi_end; i++) {
                force_sorted[i].clear();
            }
        }
        PS_OMP_PARALLEL_FOR_D
        for (S32 i = 0; i < n_ipg_this_rnd; i++) {
            const S32 ith = Comm::getThreadNum();
            const S32 adr_ipg = i + adr_ipg_first;
            const S32 n_epi = ipg[adr_ipg].n_ptcl_;
            const S32 adr_epi_head = ipg[adr_ipg].adr_ptcl_;
            const S32 n_epj = interaction_list.n_ep_[i];
            const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
            const S32 adr_epj_end = interaction_list.n_disp_ep_[i + 1];
            epj_for_force[ith].resizeNoInitialize(n_epj);
            S32 n_ep_cnt = 0;
            for (S32 j = adr_epj_head; j < adr_epj_end; j++, n_ep_cnt++) {
                const S32 adr_epj = interaction_list.adr_ep_[j];
                epj_for_force[ith][n_ep_cnt] = epj_sorted[adr_epj];
            }
            pfunc_ep_ep(epi_sorted.getPointer(adr_epi_head), n_epi, epj_for_force[ith].getPointer(), n_epj, force_sorted.getPointer(adr_epi_head));
            //S32 n_sp_cnt = 0;
            n_interaction_ep_ep_local += n_epi * n_epj;
        }
    }
}

#elif 0
// PMMM
template <class Ttc, class Tipg, class Tinteraction_list, class Tepi, class Tepj, class Tspj, class Tforce, class Tfunc_ep_ep, class Tfunc_ep_sp>
void CalcForceNoWalkPMM(const ReallocatableArray<Ttc>& tc, const ReallocatableArray<Tipg>& ipg, const Tinteraction_list& interaction_list,
                        const S32 n_loc_tot, const S32 n_let_sp, ReallocatableArray<Tepi>& epi_sorted, const ReallocatableArray<Tepj>& epj_sorted,
                        const ReallocatableArray<Tspj>& spj_sorted, ReallocatableArray<Tforce>& force_sorted, const S32 n_group_limit,
                        Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const S32 adr_ipg_first, const S32 n_ipg_this_rnd, const bool clear_force,
                        const bool update_epi, CountT& n_interaction_ep_ep_local, CountT& n_interaction_ep_sp_local, TimeProfile& time_profile) {
    // Check the arguments
    // assert(0 <= pop_id && pop_id <= N_POP);
    // assert(0 <= adr_ipg_first);
    // assert(0 <= n_ipg_this_rnd);
    // const S32 n_thread_max = Comm::getMaxThreads();
    const S32 n_thread = Comm::getNumberOfThread();
    // Allocate local buffers
    static bool is_first_call{true};
    // static ReallocatableArray<S32> * adr_epi_for_force {nullptr};
    // static ReallocatableArray<Tepi> * epi_for_force {nullptr};
    static ReallocatableArray<Tepj>* epj_for_force{nullptr};
    static ReallocatableArray<Tspj>* spj_for_force{nullptr};
    // static ReallocatableArray<Tforce> * force_for_force {nullptr};
    if (is_first_call) {
        // adr_epi_for_force = new ReallocatableArray<S32>[n_thread];
        // epi_for_force = new ReallocatableArray<Tepi>[n_thread];
        epj_for_force = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force = new ReallocatableArray<Tspj>[n_thread];
        // force_for_force = new ReallocatableArray<Tforce>[n_thread];
        is_first_call = false;
    }
    force_sorted.resizeNoInitialize(n_loc_tot);
    const S64 n_ipg = ipg.size();
    std::cerr << "n_ipg= " << n_ipg << std::endl;
    const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
    std::cerr << "n_ipg_this_rnd=" << n_ipg_this_rnd << std::endl;
    std::cerr << "adr_ipg_first=" << adr_ipg_first << std::endl;
    std::cerr << "adr_ipg_last=" << adr_ipg_last << std::endl;
    assert(adr_ipg_last < n_ipg);
    if (n_ipg_this_rnd > 0) {
        if (clear_force) {
            PS_OMP_PARALLEL_FOR
            for (S32 i = 0; i < n_loc_tot; i++) {
                force_sorted[i].clear();
            }
        }
        PS_OMP_PARALLEL_FOR_D
        for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
            std::cerr << "i= " << i << " adr_ipg_first=" << adr_ipg_first << " adr_ipg_last=" << adr_ipg_last << std::endl;
            const S32 ith = Comm::getThreadNum();
            const S32 n_epi = ipg[i].n_ptcl_;
            const S32 adr_epi_head = ipg[i].adr_ptcl_;
            const S32 n_epj = interaction_list.n_ep_[i];
            const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
            const S32 adr_epj_end = interaction_list.n_disp_ep_[i + 1];
            const S32 n_spj = interaction_list.n_sp_[i];
            const S32 adr_spj_head = interaction_list.n_disp_sp_[i];
            const S32 adr_spj_end = interaction_list.n_disp_sp_[i + 1];
            epj_for_force[ith].resizeNoInitialize(n_epj);
            spj_for_force[ith].resizeNoInitialize(n_spj);
            S32 n_ep_cnt = 0;
            for (S32 j = adr_epj_head; j < adr_epj_end; j++, n_ep_cnt++) {
                const S32 adr_epj = interaction_list.adr_ep_[j];
                epj_for_force[ith][n_ep_cnt] = epj_sorted[adr_epj];
            }
            std::cerr << "n_epi=" << n_epi << " n_epj=" << n_epj << " n_spj=" << n_spj << std::endl;
            std::cerr << "adr_epi_head=" << adr_epi_head << std::endl;
            std::cerr << "n_epj=" << n_epj << std::endl;
            pfunc_ep_ep(epi_sorted.getPointer(adr_epi_head), n_epi, epj_for_force[ith].getPointer(), n_epj, force_sorted.getPointer(adr_epi_head));
            std::cerr << "check 1" << std::endl;
            S32 n_sp_cnt = 0;
            for (S32 j = adr_spj_head; j < adr_spj_end; j++, n_sp_cnt++) {
                const S32 adr_spj = interaction_list.adr_sp_[j];
                if (adr_spj < n_let_sp) {
                    spj_for_force[ith][n_sp_cnt] = spj_sorted[adr_spj];
                } else {
                    spj_for_force[ith][n_sp_cnt].copyFromMoment(tc[adr_spj - n_let_sp].mom_);
                }
            }
            std::cerr << "check 2" << std::endl;
            pfunc_ep_sp(epi_sorted.getPointer(adr_epi_head), n_epi, spj_for_force[ith].getPointer(), n_spj, force_sorted.getPointer(adr_epi_head));
            std::cerr << "check 3" << std::endl;
            n_interaction_ep_ep_local += n_epi * n_epj;
            n_interaction_ep_sp_local += n_epi * n_spj;
            std::cerr << "n_interaction_ep_ep_local=" << n_interaction_ep_ep_local << std::endl;
            std::cerr << "n_interaction_ep_sp_local=" << n_interaction_ep_sp_local << std::endl;
        }
    } else {
        /*
//PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 4))
//PS_OMP(omp parallel for schedule(dynamic, 1))
PS_OMP_PARALLEL_FOR_D
            for(S32 i=adr_ipg_first; i<=adr_ipg_last; i++){
                const S32 ith = Comm::getThreadNum();
                const S32 offset = ipg[i].adr_ptcl_;
                const S32 n_ptcl = ipg[i].n_ptcl_;
                adr_epi_for_force[ith].clearSize();
                adr_epi_for_force[ith].reserveEmptyAreaAtLeast(n_group_limit);
                const S32 n_tail = offset + n_ptcl;
                for (S32 k=offset; k<n_tail; k++) {
                    if (GetMyPopId(epi_sorted[k]) == pop_id) {
                        adr_epi_for_force[ith].push_back(k);
                    }
                }
                const S32 n_epi = adr_epi_for_force[ith].size();
                epi_for_force[ith].resizeNoInitialize(n_epi);
                for (S32 k=0; k<n_epi; k++) {
                    const S32 adr_epi = adr_epi_for_force[ith][k];
                    epi_for_force[ith][k] = epi_sorted[adr_epi];
                }
                if(clear_force){
                    for (S32 k=0; k<n_epi; k++) {
                        const S32 adr_epi = adr_epi_for_force[ith][k];
                        force_sorted[adr_epi].clear();
                    }
                }
                force_for_force[ith].resizeNoInitialize(n_epi);
                for (S32 k=0; k<n_epi; k++) {
                    const S32 adr_epi = adr_epi_for_force[ith][k];
                    force_for_force[ith][k] = force_sorted[adr_epi];
                }
                const S32 n_epj = interaction_list.n_ep_[i];
                const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
                const S32 adr_epj_end  = interaction_list.n_disp_ep_[i+1];
                const S32 n_spj = interaction_list.n_sp_[i];
                const S32 adr_spj_head = interaction_list.n_disp_sp_[i];
                const S32 adr_spj_end  = interaction_list.n_disp_sp_[i+1];
                epj_for_force[ith].resizeNoInitialize(n_epj);
                spj_for_force[ith].resizeNoInitialize(n_spj);
                S32 n_ep_cnt = 0;
                for(S32 k=adr_epj_head; k<adr_epj_end; k++, n_ep_cnt++){
                    const S32 adr_epj = interaction_list.adr_ep_[k];
                    epj_for_force[ith][n_ep_cnt] = epj_sorted[adr_epj];
                }
                pfunc_ep_ep(epi_for_force[ith].getPointer(),   n_epi,
                            epj_for_force[ith].getPointer(),   n_epj,
                            force_for_force[ith].getPointer());
                S32 n_sp_cnt = 0;
                for(S32 k=adr_spj_head; k<adr_spj_end; k++, n_sp_cnt++){
                    const S32 adr_spj = interaction_list.adr_sp_[k];
                    if(adr_spj < n_let_sp){
                        spj_for_force[ith][n_sp_cnt] = spj_sorted[adr_spj];
                    }
                    else{
                        spj_for_force[ith][n_sp_cnt].copyFromMoment(tc[adr_spj-n_let_sp].mom_);
                    }
                }
                pfunc_ep_sp(epi_for_force[ith].getPointer(),   n_epi,
                            spj_for_force[ith].getPointer(),   n_spj,
                            force_for_force[ith].getPointer());
                for (S32 k=0; k<n_epi; k++) {
                    const S32 adr_epi = adr_epi_for_force[ith][k];
                    force_sorted[adr_epi] = force_for_force[ith][k];
                }
                if (update_epi) {
                    for (S32 k=0; k<n_epi; k++) {
                        const S32 adr_epi = adr_epi_for_force[ith][k];
                        epi_sorted[adr_epi] = epi_for_force[ith][k];
                    }
                }
                n_interaction_ep_ep_local[ith] += n_epi * n_epj;
                n_interaction_ep_sp_local[ith] += n_epi * n_spj;
            }
        }
        */
    }
}
#else
// original
template <class Ttc, class Tipg, class Tinteraction_list, class Tepi, class Tepj, class Tspj, class Tforce, class Tfunc_ep_ep, class Tfunc_ep_sp,
          U32 N_POP>
void CalcForceNoWalk(const ReallocatableArray<Ttc>& tc, const ReallocatableArray<Tipg>& ipg, const Tinteraction_list& interaction_list,
                     const S32 n_loc_tot, const S32 n_let_sp, ReallocatableArray<Tepi>& epi_sorted, const ReallocatableArray<Tepj>& epj_sorted,
                     const ReallocatableArray<Tspj>& spj_sorted, ReallocatableArray<Tforce>& force_sorted, const S32 n_group_limit,
                     Tfunc_ep_ep pfunc_ep_ep, Tfunc_ep_sp pfunc_ep_sp, const U32 pop_id, const S32 adr_ipg_first, const S32 n_ipg_this_rnd,
                     const bool clear_force, const bool update_epi, std::vector<CountT>& n_epi_processed_local,
                     std::vector<CountT>& n_epj_processed_local, std::vector<CountT>& n_spj_processed_local,
                     std::vector<CountT>& n_interaction_ep_ep_local, std::vector<CountT>& n_interaction_ep_sp_local, TimeProfile& time_profile,
                     const S32 n_thread) {
    // Check the arguments
    assert(0 <= pop_id && pop_id <= N_POP);
    assert(0 <= adr_ipg_first);
    assert(0 <= n_ipg_this_rnd);
    const S32 n_thread_max = Comm::getMaxThreads();
    // Allocate local buffers
    static bool is_first_call{true};
    static ReallocatableArray<S32>* adr_epi_for_force{nullptr};
    static ReallocatableArray<Tepi>* epi_for_force{nullptr};
    static ReallocatableArray<Tepj>* epj_for_force{nullptr};
    static ReallocatableArray<Tspj>* spj_for_force{nullptr};
    static ReallocatableArray<Tforce>* force_for_force{nullptr};
    if (is_first_call) {
        adr_epi_for_force = new ReallocatableArray<S32>[n_thread_max];
        epi_for_force = new ReallocatableArray<Tepi>[n_thread_max];
        epj_for_force = new ReallocatableArray<Tepj>[n_thread_max];
        spj_for_force = new ReallocatableArray<Tspj>[n_thread_max];
        force_for_force = new ReallocatableArray<Tforce>[n_thread_max];
        is_first_call = false;
    }
    force_sorted.resizeNoInitialize(n_loc_tot);
    const S64 n_ipg = ipg.size();
    const S32 adr_ipg_last = adr_ipg_first + n_ipg_this_rnd - 1;
    assert(adr_ipg_last < n_ipg);
    if (n_ipg_this_rnd > 0) {
        if (!(0 <= pop_id && pop_id < N_POP)) {
            if (clear_force) {
PS_OMP(omp parallel for num_threads(n_thread))
for (S32 i = 0; i < n_loc_tot; i++) {
    force_sorted[i].clear();
}
            }
            // PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 4))
PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 1))
for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
    const S32 ith = Comm::getThreadNum();
    const S32 n_epi = ipg[i].n_ptcl_;
    const S32 adr_epi_head = ipg[i].adr_ptcl_;
    const S32 n_epj = interaction_list.n_ep_[i];
    const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
    const S32 adr_epj_end = interaction_list.n_disp_ep_[i + 1];
    const S32 n_spj = interaction_list.n_sp_[i];
    const S32 adr_spj_head = interaction_list.n_disp_sp_[i];
    const S32 adr_spj_end = interaction_list.n_disp_sp_[i + 1];
    epj_for_force[ith].resizeNoInitialize(n_epj);
    spj_for_force[ith].resizeNoInitialize(n_spj);
    S32 n_ep_cnt = 0;
    for (S32 j = adr_epj_head; j < adr_epj_end; j++, n_ep_cnt++) {
        const S32 adr_epj = interaction_list.adr_ep_[j];
        epj_for_force[ith][n_ep_cnt] = epj_sorted[adr_epj];
    }
    pfunc_ep_ep(epi_sorted.getPointer(adr_epi_head), n_epi, epj_for_force[ith].getPointer(), n_epj, force_sorted.getPointer(adr_epi_head));
    S32 n_sp_cnt = 0;
    for (S32 j = adr_spj_head; j < adr_spj_end; j++, n_sp_cnt++) {
        const S32 adr_spj = interaction_list.adr_sp_[j];
        if (adr_spj < n_let_sp) {
            spj_for_force[ith][n_sp_cnt] = spj_sorted[adr_spj];
        } else {
            spj_for_force[ith][n_sp_cnt].copyFromMoment(tc[adr_spj - n_let_sp].mom_);
        }
    }
    pfunc_ep_sp(epi_sorted.getPointer(adr_epi_head), n_epi, spj_for_force[ith].getPointer(), n_spj, force_sorted.getPointer(adr_epi_head));
    n_epi_processed_local[ith] += n_epi;
    n_epj_processed_local[ith] += n_epj;
    n_spj_processed_local[ith] += n_spj;
    n_interaction_ep_ep_local[ith] += n_epi * n_epj;
    n_interaction_ep_sp_local[ith] += n_epi * n_spj;
}
        } else {
            // PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 4))
PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 1))
for (S32 i = adr_ipg_first; i <= adr_ipg_last; i++) {
    const S32 ith = Comm::getThreadNum();
    const S32 offset = ipg[i].adr_ptcl_;
    const S32 n_ptcl = ipg[i].n_ptcl_;
    adr_epi_for_force[ith].clearSize();
    adr_epi_for_force[ith].reserveEmptyAreaAtLeast(n_group_limit);
    const S32 n_tail = offset + n_ptcl;
    for (S32 k = offset; k < n_tail; k++) {
        if (GetMyPopId(epi_sorted[k]) == pop_id) {
            adr_epi_for_force[ith].push_back(k);
        }
    }
    const S32 n_epi = adr_epi_for_force[ith].size();
    epi_for_force[ith].resizeNoInitialize(n_epi);
    for (S32 k = 0; k < n_epi; k++) {
        const S32 adr_epi = adr_epi_for_force[ith][k];
        epi_for_force[ith][k] = epi_sorted[adr_epi];
    }
    if (clear_force) {
        for (S32 k = 0; k < n_epi; k++) {
            const S32 adr_epi = adr_epi_for_force[ith][k];
            force_sorted[adr_epi].clear();
        }
    }
    force_for_force[ith].resizeNoInitialize(n_epi);
    for (S32 k = 0; k < n_epi; k++) {
        const S32 adr_epi = adr_epi_for_force[ith][k];
        force_for_force[ith][k] = force_sorted[adr_epi];
    }
    const S32 n_epj = interaction_list.n_ep_[i];
    const S32 adr_epj_head = interaction_list.n_disp_ep_[i];
    const S32 adr_epj_end = interaction_list.n_disp_ep_[i + 1];
    const S32 n_spj = interaction_list.n_sp_[i];
    const S32 adr_spj_head = interaction_list.n_disp_sp_[i];
    const S32 adr_spj_end = interaction_list.n_disp_sp_[i + 1];
    epj_for_force[ith].resizeNoInitialize(n_epj);
    spj_for_force[ith].resizeNoInitialize(n_spj);
    S32 n_ep_cnt = 0;
    for (S32 k = adr_epj_head; k < adr_epj_end; k++, n_ep_cnt++) {
        const S32 adr_epj = interaction_list.adr_ep_[k];
        epj_for_force[ith][n_ep_cnt] = epj_sorted[adr_epj];
    }
    pfunc_ep_ep(epi_for_force[ith].getPointer(), n_epi, epj_for_force[ith].getPointer(), n_epj, force_for_force[ith].getPointer());
    S32 n_sp_cnt = 0;
    for (S32 k = adr_spj_head; k < adr_spj_end; k++, n_sp_cnt++) {
        const S32 adr_spj = interaction_list.adr_sp_[k];
        if (adr_spj < n_let_sp) {
            spj_for_force[ith][n_sp_cnt] = spj_sorted[adr_spj];
        } else {
            spj_for_force[ith][n_sp_cnt].copyFromMoment(tc[adr_spj - n_let_sp].mom_);
        }
    }
    pfunc_ep_sp(epi_for_force[ith].getPointer(), n_epi, spj_for_force[ith].getPointer(), n_spj, force_for_force[ith].getPointer());
    for (S32 k = 0; k < n_epi; k++) {
        const S32 adr_epi = adr_epi_for_force[ith][k];
        force_sorted[adr_epi] = force_for_force[ith][k];
    }
    if (update_epi) {
        for (S32 k = 0; k < n_epi; k++) {
            const S32 adr_epi = adr_epi_for_force[ith][k];
            epi_sorted[adr_epi] = epi_for_force[ith][k];
        }
    }
    n_epi_processed_local[ith] += n_epi;
    n_epj_processed_local[ith] += n_epj;
    n_spj_processed_local[ith] += n_spj;
    n_interaction_ep_ep_local[ith] += n_epi * n_epj;
    n_interaction_ep_sp_local[ith] += n_epi * n_spj;
}
        }
    }
}
#endif

}  // End of namespace ParticleSimulator
