#pragma once

#include"tree_walk.hpp"

namespace ParticleSimulator{
    ///////////////////////////////////////////////////
    //
    // FUNCTIONS OF WALK+FORCE WITH DOUBLE BUFFERING 
    //
    ///////////////////////////////////////////////////

    template<class Tforce>
    void CalcForceMultiWalkInitialize(ReallocatableArray<Tforce> & force_sorted,
                                      const S32 n_loc_tot,
                                      const bool tag_max,
                                      const bool clear){
        if(tag_max <= 0){
            PARTICLE_SIMULATOR_PRINT_ERROR("tag_max is illegal. In currente version, tag_max must be 1");
            Abort(-1);
        }
        force_sorted.resizeNoInitialize(n_loc_tot);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot; i++){
                force_sorted[i].clear();
            }
        }
    }
    
    //////////////////////////////////////////////////////////////
    //////////// Walk+Force, Kernel:Index, List:Index ////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                            Tfunc_retrieve pfunc_retrieve,
                            const S32 tag_max,
                            const S32 n_walk_limit,
                            const bool flag_keep_list,
                            const bool clear){
        CalcForceMultiWalkInitialize(force_sorted_, n_loc_tot_, tag_max, clear);
        S32 ret = 0;
        const F64 wtime_offset = GetWtime();
        ret = calcForceMultiWalkIndexImpl(typename TSM::force_type(),
                                          pfunc_dispatch,
                                          pfunc_retrieve,
                                          tag_max,
                                          n_walk_limit,
                                          flag_keep_list,
                                          clear);
        time_profile_.calc_force += GetWtime() - wtime_offset;
        return ret;
    }

    //////////// Walk+Force, Kernel:Index, List:Index, Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexImpl(TagForceLong,
                                Tfunc_dispatch pfunc_dispatch,
                                Tfunc_retrieve pfunc_retrieve,
                                const S32 tag_max,
                                const S32 n_walk_limit,
                                const bool flag_keep_list,
                                const bool clear){
        const F64 offset_core = GetWtime();
        // send all epj and spj
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;
        S32 ** id_spj_dummy = NULL;
        S32 *  n_spj_dummy = NULL;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       (const S32**)id_spj_dummy, n_spj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       spj_sorted_.getPointer(), spj_sorted_.size(),
                       true);        
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        n_walk_local_ += n_ipg;
        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
            interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_sp_.clearSize();
        }
        const S32 n_loop_max = n_ipg/n_walk_limit + ((n_ipg%n_walk_limit)==0 ? 0 : 1);
        ReallocatableArray<Tforce*> ptr_force_per_walk[2];
        ReallocatableArray<S32> n_epi_per_walk[2];
        ReallocatableArray<Tepi*> ptr_epi_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32> n_epj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32> n_spj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32*> ptr_adr_epj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32*> ptr_adr_spj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32> * n_epj_disp_per_thread;
        ReallocatableArray<S32> * n_spj_disp_per_thread;
        ptr_force_per_walk[0].initialize(n_walk_limit, n_walk_limit, 1); // array of pointer *[n_walk]
        ptr_force_per_walk[1].initialize(n_walk_limit, n_walk_limit, 1); // array of pointer *[n_walk]
        n_epi_per_walk[0].initialize(n_walk_limit, n_walk_limit, 1);
        n_epi_per_walk[1].initialize(n_walk_limit, n_walk_limit, 1);        
        n_epj_disp_per_thread = new ReallocatableArray<S32>[n_thread];
        n_spj_disp_per_thread = new ReallocatableArray<S32>[n_thread];
        for(int i=0; i<n_thread; i++){
            n_epj_disp_per_thread[i].initialize(n_walk_limit+1, n_walk_limit+1, 1);
            n_spj_disp_per_thread[i].initialize(n_walk_limit+1, n_walk_limit+1, 1);
        }
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        const S32 adr_tree_sp_first = comm_table_.n_sp_recv_tot_;
        bool first_loop = true;
        S32 n_walk_prev = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_ipg/n_loop_max + (((n_ipg%n_loop_max) > wg) ? 1 : 0);
                const S32 walk_grp_head = (n_ipg/n_loop_max)*wg + std::min((n_ipg%n_loop_max), wg);
                const F64 offset_calc_force__core__walk_tree = GetWtime();
                const S32 lane_now = wg%2;
                const S32 lane_old = (wg+1)%2;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    ReallocatableArray<S32> iwloc2iw(n_walk_limit, n_walk_limit, 1);
                    S32 n_walk_loc = 0;
                    S32 n_ep_cum = 0;
                    S32 n_sp_cum = 0;
                    n_epj_disp_per_thread[ith][0] = 0;
                    n_spj_disp_per_thread[ith][0] = 0;
                    adr_epj_for_force_[ith].clearSize();
                    adr_spj_for_force_[ith].clearSize();
                    adr_ipg_for_force_[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) reduction(+: n_interaction_ep_ep_tmp) reduction(+: n_interaction_ep_sp_tmp)
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_for_force_[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_per_walk[lane_now][iw] = ipg_[id_ipg].n_ptcl_;
                        ptr_epi_per_walk[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        ptr_force_per_walk[lane_now][iw] = force_sorted_.getPointer(first_adr_ip);
                        TargetBox<TSM> target_box;
                        target_box.set(ipg_[id_ipg]);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, 
                             Tspj, WALK_MODE_NORMAL, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_for_force_[ith],
                             spj_sorted_, adr_spj_for_force_[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             adr_tree_sp_first, F64vec(0.0));
                        n_epj_per_walk[iw] = adr_epj_for_force_[ith].size() - n_ep_cum;
                        n_spj_per_walk[iw] = adr_spj_for_force_[ith].size() - n_sp_cum;
                        n_ep_cum = adr_epj_for_force_[ith].size();
                        n_sp_cum = adr_spj_for_force_[ith].size();
                        n_epj_disp_per_thread[ith][n_walk_loc+1] = n_ep_cum;
                        n_spj_disp_per_thread[ith][n_walk_loc+1] = n_sp_cum;
                        n_interaction_ep_ep_tmp += ((S64)n_epj_per_walk[iw]*(S64)n_epi_per_walk[lane_now][iw]);
                        n_interaction_ep_sp_tmp += ((S64)n_spj_per_walk[iw]*(S64)n_epi_per_walk[lane_now][iw]);
                        iwloc2iw[n_walk_loc] = iw;
                        n_walk_loc++;
                    } // end of OMP for
                    for(S32 iwloc=0; iwloc<n_walk_loc; iwloc++){
                        S32 iw = iwloc2iw[iwloc];
                        ptr_adr_epj_per_walk[iw] = adr_epj_for_force_[ith].getPointer(n_epj_disp_per_thread[ith][iwloc]);
                        ptr_adr_spj_per_walk[iw] = adr_spj_for_force_[ith].getPointer(n_spj_disp_per_thread[ith][iwloc]);
                    }
                } // end of OMP parallel scope
                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = interaction_list_.n_disp_sp_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_per_walk[iw];
                        interaction_list_.n_sp_[id_ipg] = n_spj_per_walk[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                        interaction_list_.n_disp_sp_[id_ipg+1] = interaction_list_.n_disp_sp_[id_ipg] + interaction_list_.n_sp_[id_ipg];
                    }
                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
                    interaction_list_.adr_sp_.resizeNoInitialize( interaction_list_.n_disp_sp_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_for_force_[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_for_force_[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_per_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_per_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_for_force_[i][k];
                            }
                            S32 adr_sp = interaction_list_.n_disp_sp_[adr_ipg];
                            const S32 k_sp_h = n_spj_disp_per_thread[i][j];
                            const S32 k_sp_e = n_spj_disp_per_thread[i][j+1];
                            for(S32 k=k_sp_h; k<k_sp_e; k++, adr_sp++){
                                interaction_list_.adr_sp_[adr_sp] = adr_spj_for_force_[i][k];
                            }
                            
                        }
                    }
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[lane_old].getPointer(), ptr_force_per_walk[lane_old].getPointer());
                } // retrieve
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)ptr_epi_per_walk.getPointer(),   n_epi_per_walk[lane_now].getPointer(),
                               (const S32**)ptr_adr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(),
                               (const S32**)ptr_adr_spj_per_walk.getPointer(), n_spj_per_walk.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               spj_sorted_.getPointer(), spj_sorted_.size(),
                               false);
                first_loop = false;
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[(n_loop_max+1)%2].getPointer(), ptr_force_per_walk[(n_loop_max+1)%2].getPointer());
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
        
        ptr_force_per_walk[0].freeMem();
        ptr_force_per_walk[1].freeMem();
        n_epi_per_walk[0].freeMem();
        n_epi_per_walk[1].freeMem();
        for(int i=0; i<n_thread; i++){
            n_epj_disp_per_thread[i].freeMem();
            n_spj_disp_per_thread[i].freeMem();
        }
        delete [] n_epj_disp_per_thread;
        delete [] n_spj_disp_per_thread;
        return ret;
    }
    
    //////////// Walk+Force, Kernel:Index, List:Index, Force:Short //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexImpl(TagForceShort,
                                Tfunc_dispatch pfunc_dispatch,
                                Tfunc_retrieve pfunc_retrieve,
                                const S32 tag_max,
                                const S32 n_walk_limit,
                                const bool flag_keep_list,
                                const bool clear){
        const F64 offset_core = GetWtime();
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       true);

        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        //const S32 n_ipg_amari = (n_ipg > 0) ? n_ipg%n_walk_limit : 0;
        n_walk_local_ += n_ipg;
        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
        }
        const S32 n_loop_max = n_ipg/n_walk_limit + ((n_ipg%n_walk_limit)==0 ? 0 : 1);
        ReallocatableArray<Tforce*> ptr_force_per_walk[2];
        ReallocatableArray<S32> n_epi_per_walk[2];
        ReallocatableArray<Tepi*> epi_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32> n_epj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32*> adr_epj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32> * n_epj_disp_per_thread;

        ptr_force_per_walk[0].initialize(n_walk_limit, n_walk_limit, 1); // array of pointer *[n_walk]
        ptr_force_per_walk[1].initialize(n_walk_limit, n_walk_limit, 1); // array of pointer *[n_walk]
        n_epi_per_walk[0].initialize(n_walk_limit, n_walk_limit, 1);
        n_epi_per_walk[1].initialize(n_walk_limit, n_walk_limit, 1);        
        n_epj_disp_per_thread = new ReallocatableArray<S32>[n_thread];
        for(int i=0; i<n_thread; i++){
            n_epj_disp_per_thread[i].initialize(n_walk_limit+1, n_walk_limit+1, 1);
        }
        
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        bool first_loop = true;
        S32 n_walk_prev = 0;
        S64 n_interaction_ep_ep_tmp = 0;

        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_ipg/n_loop_max + (((n_ipg%n_loop_max) > wg) ? 1 : 0);
                const S32 walk_grp_head = (n_ipg/n_loop_max)*wg + std::min((n_ipg%n_loop_max), wg);
                const F64 offset_calc_force__core__walk_tree = GetWtime();
                const S32 lane_now = wg%2;
                const S32 lane_old = (wg+1)%2;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    ReallocatableArray<S32> iwloc2iw(n_walk_limit, n_walk_limit, 1);
                    S32 n_walk_loc = 0;
                    S32 n_ep_cum = 0;
                    n_epj_disp_per_thread[ith][0] = 0;
                    adr_epj_for_force_[ith].clearSize();
                    adr_ipg_for_force_[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) reduction(+: n_interaction_ep_ep_tmp)
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_for_force_[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_per_walk[lane_now][iw] = ipg_[id_ipg].n_ptcl_;
                        epi_per_walk[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        ptr_force_per_walk[lane_now][iw] = force_sorted_.getPointer(first_adr_ip);

                        TargetBox<TSM> target_box;
                        target_box.set(ipg_[id_ipg]);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj,
                             WALK_MODE_NORMAL, TagChopLeafTrue>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_for_force_[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             F64vec(0.0));
                        n_epj_per_walk[iw] = adr_epj_for_force_[ith].size() - n_ep_cum;
                        n_ep_cum = adr_epj_for_force_[ith].size();
                        n_epj_disp_per_thread[ith][n_walk_loc+1] = n_ep_cum;
                        n_interaction_ep_ep_tmp += ((S64)n_epj_per_walk[iw]*(S64)n_epi_per_walk[lane_now][iw]);
                        iwloc2iw[n_walk_loc] = iw;
                        n_walk_loc++;
                    } // end of OMP for
                    for(S32 iwloc=0; iwloc<n_walk_loc; iwloc++){
                        S32 iw = iwloc2iw[iwloc];
                        adr_epj_per_walk[iw] = adr_epj_for_force_[ith].getPointer(n_epj_disp_per_thread[ith][iwloc]);
                    }
                } // end of OMP parallel scope
                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_per_walk[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                    }
                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_for_force_[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_for_force_[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_per_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_per_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_for_force_[i][k];
                            }
                        }
                    }
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[lane_old].getPointer(), ptr_force_per_walk[lane_old].getPointer());
                } // retrieve
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_per_walk.getPointer(),   n_epi_per_walk[lane_now].getPointer(),
                               (const S32**)adr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               false);
                first_loop = false;
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[(n_loop_max+1)%2].getPointer(), ptr_force_per_walk[(n_loop_max+1)%2].getPointer());
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = 0;
            n_interaction_ep_ep_local_ = 0;
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
        ptr_force_per_walk[0].freeMem();
        ptr_force_per_walk[1].freeMem();
        n_epi_per_walk[0].freeMem();
        n_epi_per_walk[1].freeMem();
        for(int i=0; i<n_thread; i++){
            n_epj_disp_per_thread[i].freeMem();
        }
        delete [] n_epj_disp_per_thread;
        return ret;
    }
    
    //////////////////////////////////////////////////
    //////////// Walk+Force, Kernel:Ptcl, List:Index //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkPtcl(Tfunc_dispatch pfunc_dispatch,
                           Tfunc_retrieve pfunc_retrieve,
                           const S32 tag_max,
                           const S32 n_walk_limit,
                           const bool flag_keep_list,
                           const bool clear){
        CalcForceMultiWalkInitialize(force_sorted_, n_loc_tot_, tag_max, clear);
        S32 ret = 0;
        const F64 time_offset = GetWtime();
        ret = calcForceMultiWalkPtclImpl(typename TSM::force_type(),
                                             pfunc_dispatch,
                                             pfunc_retrieve,
                                             tag_max,
                                             n_walk_limit,
                                             flag_keep_list,
                                             clear);
        time_profile_.calc_force += GetWtime() - time_offset;
        return ret;
    }

    //////////// Walk+Force, Kernel:Ptcl, List:Index Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkPtclImpl(TagForceLong,
                               Tfunc_dispatch pfunc_dispatch,
                               Tfunc_retrieve pfunc_retrieve,
                               const S32 tag_max,
                               const S32 n_walk_limit,
                               const bool flag_keep_list,
                               const bool clear){
        const F64 offset_core = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        n_walk_local_ += n_ipg;
        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
            interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_sp_.clearSize();
        }
        const S32 n_loop_max = n_ipg/n_walk_limit + ((n_ipg%n_walk_limit)==0 ? 0 : 1);
        ReallocatableArray<Tforce*> ptr_force_per_walk[2];
        ReallocatableArray<S32> n_epi_per_walk[2];
        ReallocatableArray<Tepi*> ptr_epi_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32> n_epj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32> n_spj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<Tepj*> ptr_epj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<Tspj*> ptr_spj_per_walk(n_walk_limit, n_walk_limit, 1);        
        ReallocatableArray<S32> * n_epj_disp_per_thread;
        ReallocatableArray<S32> * n_spj_disp_per_thread;
        ptr_force_per_walk[0].initialize(n_walk_limit, n_walk_limit, 1); // array of pointer *[n_walk]
        ptr_force_per_walk[1].initialize(n_walk_limit, n_walk_limit, 1); // array of pointer *[n_walk]
        n_epi_per_walk[0].initialize(n_walk_limit, n_walk_limit, 1);
        n_epi_per_walk[1].initialize(n_walk_limit, n_walk_limit, 1);        
        n_epj_disp_per_thread = new ReallocatableArray<S32>[n_thread];
        n_spj_disp_per_thread = new ReallocatableArray<S32>[n_thread];
        for(int i=0; i<n_thread; i++){
            n_epj_disp_per_thread[i].initialize(n_walk_limit+1, n_walk_limit+1, 1);
            n_spj_disp_per_thread[i].initialize(n_walk_limit+1, n_walk_limit+1, 1);
        }
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        const S32 adr_tree_sp_first = comm_table_.n_sp_recv_tot_;
        bool first_loop = true;
        S32 n_walk_prev = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_ipg/n_loop_max + (((n_ipg%n_loop_max) > wg) ? 1 : 0);
                const S32 walk_grp_head = (n_ipg/n_loop_max)*wg + std::min((n_ipg%n_loop_max), wg);
                const F64 offset_calc_force__core__walk_tree = GetWtime();
                const S32 lane_now = wg%2;
                const S32 lane_old = (wg+1)%2;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    ReallocatableArray<S32> iwloc2iw(n_walk_limit, n_walk_limit, 1);
                    S32 n_walk_loc = 0;
                    S32 n_ep_cum = 0;
                    S32 n_sp_cum = 0;
                    n_epj_disp_per_thread[ith][0] = 0;
                    n_spj_disp_per_thread[ith][0] = 0;
                    adr_epj_for_force_[ith].clearSize();
                    adr_spj_for_force_[ith].clearSize();
                    adr_ipg_for_force_[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) reduction(+: n_interaction_ep_ep_tmp) reduction(+: n_interaction_ep_sp_tmp)
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_for_force_[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_per_walk[lane_now][iw] = ipg_[id_ipg].n_ptcl_;
                        ptr_epi_per_walk[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        ptr_force_per_walk[lane_now][iw] = force_sorted_.getPointer(first_adr_ip);
                        TargetBox<TSM> target_box;
                        target_box.set(ipg_[id_ipg]);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, 
                             Tspj, WALK_MODE_NORMAL, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_for_force_[ith],
                             spj_sorted_, adr_spj_for_force_[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             adr_tree_sp_first, F64vec(0.0));
                        n_epj_per_walk[iw] = adr_epj_for_force_[ith].size() - n_ep_cum;
                        n_spj_per_walk[iw] = adr_spj_for_force_[ith].size() - n_sp_cum;
                        n_ep_cum = adr_epj_for_force_[ith].size();
                        n_sp_cum = adr_spj_for_force_[ith].size();
                        n_epj_disp_per_thread[ith][n_walk_loc+1] = n_ep_cum;
                        n_spj_disp_per_thread[ith][n_walk_loc+1] = n_sp_cum;
                        n_interaction_ep_ep_tmp += ((S64)n_epj_per_walk[iw]*(S64)n_epi_per_walk[lane_now][iw]);
                        n_interaction_ep_sp_tmp += ((S64)n_spj_per_walk[iw]*(S64)n_epi_per_walk[lane_now][iw]);
                        iwloc2iw[n_walk_loc] = iw;
                        n_walk_loc++;
                    } // end of OMP for
                    epj_for_force_[ith].resizeNoInitialize(adr_epj_for_force_[ith].size());
                    spj_for_force_[ith].resizeNoInitialize(adr_spj_for_force_[ith].size());
                    for(S32 jp=0; jp<adr_epj_for_force_[ith].size(); jp++){
                        S32 adr = adr_epj_for_force_[ith][jp];
                        epj_for_force_[ith][jp] = epj_sorted_[adr];
                    }
                    for(S32 jp=0; jp<adr_spj_for_force_[ith].size(); jp++){
                        S32 adr = adr_spj_for_force_[ith][jp];
                        spj_for_force_[ith][jp] = spj_sorted_[adr];
                    }
                    for(S32 iw_loc=0; iw_loc<n_walk_loc; iw_loc++){
                        S32 iw = iwloc2iw[iw_loc];
                        ptr_epj_per_walk[iw] = epj_for_force_[ith].getPointer(n_epj_disp_per_thread[ith][iw_loc]);
                        ptr_spj_per_walk[iw] = spj_for_force_[ith].getPointer(n_spj_disp_per_thread[ith][iw_loc]);
                    }
                } // end of OMP parallel scope
                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = interaction_list_.n_disp_sp_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_per_walk[iw];
                        interaction_list_.n_sp_[id_ipg] = n_spj_per_walk[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                        interaction_list_.n_disp_sp_[id_ipg+1] = interaction_list_.n_disp_sp_[id_ipg] + interaction_list_.n_sp_[id_ipg];
                    }
                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
                    interaction_list_.adr_sp_.resizeNoInitialize( interaction_list_.n_disp_sp_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_for_force_[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_for_force_[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_per_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_per_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_for_force_[i][k];
                            }
                            S32 adr_sp = interaction_list_.n_disp_sp_[adr_ipg];
                            const S32 k_sp_h = n_spj_disp_per_thread[i][j];
                            const S32 k_sp_e = n_spj_disp_per_thread[i][j+1];
                            for(S32 k=k_sp_h; k<k_sp_e; k++, adr_sp++){
                                interaction_list_.adr_sp_[adr_sp] = adr_spj_for_force_[i][k];
                            }
                            
                        }
                    }
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[lane_old].getPointer(), ptr_force_per_walk[lane_old].getPointer());
                } // retrieve
                //pfunc_dispatch(0, n_walk, 
                //               (const Tepi**)ptr_epi_per_walk.getPointer(),   n_epi_per_walk[lane_now].getPointer(),
                //               (const S32**)ptr_adr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(),
                //               (const S32**)ptr_adr_spj_per_walk.getPointer(), n_spj_per_walk.getPointer(),
                //               epj_sorted_.getPointer(), epj_sorted_.size(),
                //               spj_sorted_.getPointer(), spj_sorted_.size(),
                //               false);
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(), 
                                      (const Tepj**)ptr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer(), 
                                      (const Tspj**)ptr_spj_per_walk.getPointer(), n_spj_per_walk.getPointer());
                first_loop = false;
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[(n_loop_max+1)%2].getPointer(), ptr_force_per_walk[(n_loop_max+1)%2].getPointer());
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
        ptr_force_per_walk[0].freeMem();
        ptr_force_per_walk[1].freeMem();
        n_epi_per_walk[0].freeMem();
        n_epi_per_walk[1].freeMem();
        for(int i=0; i<n_thread; i++){
            n_epj_disp_per_thread[i].freeMem();
            n_spj_disp_per_thread[i].freeMem();
        }
        delete [] n_epj_disp_per_thread;
        delete [] n_spj_disp_per_thread;
        return ret;
    }
    //////////// Walk+Force, Kernel:Ptcl, List:Index Force:Short //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkPtclImpl(TagForceShort,
                               Tfunc_dispatch pfunc_dispatch,
                               Tfunc_retrieve pfunc_retrieve,
                               const S32 tag_max,
                               const S32 n_walk_limit,
                               const bool flag_keep_list,
                               const bool clear){
        const F64 offset_core = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        const S32 n_ipg_amari = (n_ipg > 0) ? n_ipg%n_walk_limit : 0;
        n_walk_local_ += n_ipg;
        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
        }
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg_amari==0 ? 0 : 1);
        ReallocatableArray<Tforce*> ptr_force_per_walk[2];
        ReallocatableArray<S32> n_epi_per_walk[2];
        ReallocatableArray<Tepi*> ptr_epi_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<S32> n_epj_per_walk(n_walk_limit, n_walk_limit, 1);
        ReallocatableArray<Tepj*> ptr_epj_per_walk(n_walk_limit, n_walk_limit, 1); // new
        ReallocatableArray<S32> * n_epj_disp_per_thread;
        ptr_force_per_walk[0].initialize(n_walk_limit, n_walk_limit, 1); // array of pointer *[n_walk]
        ptr_force_per_walk[1].initialize(n_walk_limit, n_walk_limit, 1); // array of pointer *[n_walk]
        n_epi_per_walk[0].initialize(n_walk_limit, n_walk_limit, 1);
        n_epi_per_walk[1].initialize(n_walk_limit, n_walk_limit, 1);        
        n_epj_disp_per_thread = new ReallocatableArray<S32>[n_thread];
        for(int i=0; i<n_thread; i++){
            n_epj_disp_per_thread[i].initialize(n_walk_limit+1, n_walk_limit+1, 1);
        }
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        bool first_loop = true;
        S32 n_walk_prev = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_ipg/n_loop_max + (((n_ipg%n_loop_max) > wg) ? 1 : 0);
                const S32 walk_grp_head = (n_ipg/n_loop_max)*wg + std::min((n_ipg%n_loop_max), wg);
                const F64 offset_calc_force__core__walk_tree = GetWtime();
                const S32 lane_now = wg%2;
                const S32 lane_old = (wg+1)%2;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    ReallocatableArray<S32> iwloc2iw(n_walk_limit, n_walk_limit, 1);
                    S32 n_walk_loc = 0;
                    S32 n_ep_cum = 0;
                    n_epj_disp_per_thread[ith][0] = 0;
                    adr_epj_for_force_[ith].clearSize();
                    adr_ipg_for_force_[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) reduction(+: n_interaction_ep_ep_tmp)
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_for_force_[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_per_walk[lane_now][iw] = ipg_[id_ipg].n_ptcl_;
                        ptr_epi_per_walk[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        ptr_force_per_walk[lane_now][iw] = force_sorted_.getPointer(first_adr_ip);

                        TargetBox<TSM> target_box;
                        target_box.set(ipg_[id_ipg]);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj,
                             WALK_MODE_NORMAL, TagChopLeafTrue>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_for_force_[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             F64vec(0.0));
                        n_epj_per_walk[iw] = adr_epj_for_force_[ith].size() - n_ep_cum;
                        n_ep_cum = adr_epj_for_force_[ith].size();
                        n_epj_disp_per_thread[ith][n_walk_loc+1] = n_ep_cum;
                        n_interaction_ep_ep_tmp += ((S64)n_epj_per_walk[iw]*(S64)n_epi_per_walk[lane_now][iw]);
                        iwloc2iw[n_walk_loc] = iw;
                        n_walk_loc++;
                    } // end of OMP for
                    epj_for_force_[ith].resizeNoInitialize(adr_epj_for_force_[ith].size());
                    for(S32 jp=0; jp<adr_epj_for_force_[ith].size(); jp++){
                        S32 adr = adr_epj_for_force_[ith][jp];
                        epj_for_force_[ith][jp] = epj_sorted_[adr];
                    }
                    for(S32 iw_loc=0; iw_loc<n_walk_loc; iw_loc++){
                        S32 iw = iwloc2iw[iw_loc];
                        ptr_epj_per_walk[iw] = epj_for_force_[ith].getPointer(n_epj_disp_per_thread[ith][iw_loc]);
                    }
                } // end of OMP parallel scope
                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_per_walk[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                    }
                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_for_force_[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_for_force_[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_per_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_per_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_for_force_[i][k];
                            }
                        }
                    }
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[lane_old].getPointer(), ptr_force_per_walk[lane_old].getPointer());
                } // retrieve
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)ptr_epi_per_walk.getPointer(), n_epi_per_walk[lane_now].getPointer(), 
                                      (const Tepj**)ptr_epj_per_walk.getPointer(), n_epj_per_walk.getPointer()); // new
                
                first_loop = false;
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_per_walk[(n_loop_max+1)%2].getPointer(), ptr_force_per_walk[(n_loop_max+1)%2].getPointer());
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = 0;
            n_interaction_ep_ep_local_ = 0;
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
        ptr_force_per_walk[0].freeMem();
        ptr_force_per_walk[1].freeMem();
        n_epi_per_walk[0].freeMem();
        n_epi_per_walk[1].freeMem();
        for(int i=0; i<n_thread; i++){
            n_epj_disp_per_thread[i].freeMem();
        }
        delete [] n_epj_disp_per_thread;
        return ret;
    }
    //////////// Kernel:Ptcl, List:Index //////////////
    //////////////////////////////////////////////////


    ///////////////////////////////////////////////////
    //
    // FUNCTIONS OF FORCE WITHOUT WALK
    // (MUST BE USED AFTER WALK)
    // 
    ///////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    //////////// Force Only, Kernel:Index, List:Index //////////////

    //////////// Force Only, Kernel:Index, List:Index, Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkIndexImpl(TagForceLong,
                                         Tfunc_dispatch pfunc_dispatch,
                                         Tfunc_retrieve pfunc_retrieve,
                                         const S32 n_walk_limit,
                                         const bool clear){
        F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_epj_ar;
        id_epj_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_spj_ar;
        n_spj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_spj_ar;
        id_spj_ar.resizeNoInitialize(n_walk_limit);

        force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_sorted_[i].clear();
                //force_org_[i].clear();
            }
        }
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;
        S32 ** id_spj_dummy = NULL;
        S32 *  n_spj_dummy = NULL;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       (const S32**)id_spj_dummy, n_spj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       spj_sorted_.getPointer(), spj_sorted_.size(),
                       true);
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        //F64 time_offset = GetWtime();
        if(n_ipg > 0){
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_; 
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    id_epj_ar[iw] = interaction_list_.adr_ep_.getPointer(interaction_list_.n_disp_ep_[id_walk]);

                    n_spj_ar[iw]  = interaction_list_.n_sp_[id_walk];
                    id_spj_ar[iw] = interaction_list_.adr_sp_.getPointer(interaction_list_.n_disp_sp_[id_walk]);

                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                    n_interaction_ep_sp_local_ += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar[lane_0].getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               (const S32**)id_spj_ar.getPointer(), n_spj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               spj_sorted_.getPointer(), spj_sorted_.size(),
                               false);
                n_walk_prev = n_walk;
                first_loop = false;
            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
        }
	time_profile_.calc_force__core += GetWtime() - time_offset;
	const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
	time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
    }

    //////////// Force Only, Kernel:Index, List:Index, Force:Short //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkIndexImpl(TagForceShort,
                                         Tfunc_dispatch pfunc_dispatch,
                                         Tfunc_retrieve pfunc_retrieve,
                                         const S32 n_walk_limit,
                                         const bool clear){
        F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
        // overlape version
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_epj_ar;
        id_epj_ar.resizeNoInitialize(n_walk_limit);

        force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_sorted_[i].clear();
                //force_org_[i].clear();
            }
        }
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;

        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       true);
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        //F64 time_offset = GetWtime();
        if(n_ipg > 0){
            // overlape version
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    id_epj_ar[iw] = interaction_list_.adr_ep_.getPointer(interaction_list_.n_disp_ep_[id_walk]);

                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar[lane_0].getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               false);
                n_walk_prev = n_walk;
                first_loop = false;
            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
        }
        time_profile_.calc_force__core += GetWtime() - time_offset;
	const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
	time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
    }

    ///////////////////////////////////////////////////
    //////////// Force Only, Kernel:Ptcl, List:Index //////////////

    //////////// Force Only, Kernel:Ptcl, List:Index, Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkPtclImpl(TagForceLong,
                                        Tfunc_dispatch pfunc_dispatch,
                                        Tfunc_retrieve pfunc_retrieve,
                                        const S32 n_walk_limit,
                                        const bool clear){

        //F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
        
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepj*> epj_ar;
        epj_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_spj_ar;
        n_spj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tspj*> spj_ar;
        spj_ar.resizeNoInitialize(n_walk_limit);

        force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_sorted_[i].clear();
                //force_org_[i].clear();
            }
        }
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        if(n_ipg > 0){
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                n_walk_prev = (wg>0) ? n_walk_ar[wg-1] : 0;
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);
                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    n_spj_ar[iw]  = interaction_list_.n_sp_[id_walk];
                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                    n_interaction_ep_sp_local_ += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }

                const S64 n_ep_head = interaction_list_.n_disp_ep_[n_walk_head];
                const S64 n_ep_end  = interaction_list_.n_disp_ep_[n_walk_head+n_walk];
                const S64 n_epj_tot = n_ep_end - n_ep_head;
                epj_for_force_[0].resizeNoInitialize(n_epj_tot);
                for(S32 jp=0; jp<n_epj_tot; jp++){
                    epj_for_force_[0][jp] = epj_sorted_[ interaction_list_.adr_ep_[jp+n_ep_head] ];
                }
                const S64 n_sp_head = interaction_list_.n_disp_sp_[n_walk_head];
                const S64 n_sp_end  = interaction_list_.n_disp_sp_[n_walk_head+n_walk];
                const S64 n_spj_tot = n_sp_end - n_sp_head;
                spj_for_force_[0].resizeNoInitialize(n_spj_tot);
                for(S32 jp=0; jp<n_spj_tot; jp++){
                    spj_for_force_[0][jp] = spj_sorted_[ interaction_list_.adr_sp_[jp+n_sp_head] ];
                }
                S64 n_epj_cnt = 0;
                S64 n_spj_cnt = 0;
                epj_ar.resizeNoInitialize(n_walk);
                spj_ar.resizeNoInitialize(n_walk);
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    spj_ar[iw] = spj_for_force_[0].getPointer(n_spj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                    n_spj_cnt += n_spj_ar[iw];
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)epi_ar.getPointer(), n_epi_ar[lane_0].getPointer(), 
                                      (const Tepj**)epj_ar.getPointer(), n_epj_ar.getPointer(), 
                                      (const Tspj**)spj_ar.getPointer(), n_spj_ar.getPointer());
                first_loop = false;
            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
        }
        copyForceOriginalOrder();
        //time_profile_.calc_force += GetWtime() - time_offset;
    }

    //////////// Force Only, Kernel:Ptcl, List:Index, Force:Short //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkPtclImpl(TagForceShort,
                                        Tfunc_dispatch pfunc_dispatch,
                                        Tfunc_retrieve pfunc_retrieve,
                                        const S32 n_walk_limit,
                                        const bool clear){
        F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepj*> epj_ar;
        epj_ar.resizeNoInitialize(n_walk_limit);

        force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_sorted_[i].clear();
                //force_org_[i].clear();
            }
        }
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        if(n_ipg > 0){
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                n_walk_prev = (wg>0) ? n_walk_ar[wg-1] : 0;
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);
                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }
                const S64 n_ep_head = interaction_list_.n_disp_ep_[n_walk_head];
                const S64 n_ep_end  = interaction_list_.n_disp_ep_[n_walk_head+n_walk];
                const S64 n_epj_tot = n_ep_end - n_ep_head;
                epj_for_force_[0].resizeNoInitialize(n_epj_tot);
                for(S32 jp=0; jp<n_epj_tot; jp++){
                    epj_for_force_[0][jp] = epj_sorted_[ interaction_list_.adr_ep_[jp+n_ep_head] ];
                }
                S64 n_epj_cnt = 0;
                epj_ar.resizeNoInitialize(n_walk);
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)epi_ar.getPointer(), n_epi_ar[lane_0].getPointer(), 
                                      (const Tepj**)epj_ar.getPointer(), n_epj_ar.getPointer());
                first_loop = false;

            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
        }
	time_profile_.calc_force__core += GetWtime() - time_offset;
	const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
	time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
    }


    //////////////////////////////////////////////////////////////
    //////////// Walk+Force, Kernel:Ptcl, List:  ////////////    
    // SHORT
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForce(Tfunc_ep_ep pfunc_ep_ep,
              const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
        F64 offset_walk_tree,offset_dispatch;

        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                offset_walk_tree = GetWtime();
                const S32 ith = Comm::getThreadNum();
                epj_for_force_[ith].clearSize();
                adr_epj_for_force_[ith].clearSize();
                makeInteractionList(i);
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_walk_tree;
                ni_tmp += ipg_[i].n_ptcl_;
                nj_tmp += epj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * epj_for_force_[Comm::getThreadNum()].size();
                offset_dispatch = GetWtime();
                calcForceOnly( pfunc_ep_ep, i, clear);
                time_profile_.calc_force__core__dispatch += GetWtime() - offset_dispatch;
            }
            ni_ave_ = ni_tmp / n_ipg;
            nj_ave_ = nj_tmp / n_ipg;
            n_interaction_ep_ep_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
        }
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = 0;
            n_walk_local_ = n_interaction_ep_ep_local_ = 0;
        }
        copyForceOriginalOrder();
        
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
        std::cout<<"force_sorted_.size()="<<force_sorted_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceOnly(Tfunc_ep_ep pfunc_ep_ep,
                  const S32 adr_ipg,
                  const bool clear){
        const S32 offset = ipg_[adr_ipg].adr_ptcl_;
        const S32 n_epi = ipg_[adr_ipg].n_ptcl_;
        const S32 ith = Comm::getThreadNum();
        const S32 n_epj = epj_for_force_[ith].size();
        const S32 n_tail = offset + n_epi;
        if(clear){
            for(S32 i=offset; i<n_tail; i++) force_sorted_[i].clear();
        }
        pfunc_ep_ep(epi_sorted_.getPointer(offset),     n_epi,
                    epj_for_force_[ith].getPointer(),   n_epj,
                    force_sorted_.getPointer(offset));
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagForceShort, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        if (clear){
            epj_for_force_[ith].clearSize();
        }
        const F64 r_crit_sq = 9999.9;
        TargetBox<TSM> target_box;
        target_box.set(ipg_[adr_ipg]);
        S32 adr_tc = 0;
        S32 n_head = adr_epj_for_force_[ith].size();
        MakeListUsingTreeRecursiveTop
            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, WALK_MODE_NORMAL, TagChopLeafTrue>
            (tc_glb_,  adr_tc, tp_glb_,
             epj_sorted_, adr_epj_for_force_[ith],
             target_box,
             r_crit_sq, n_leaf_limit_,
             F64vec(0.0));
        S32 n_tail = adr_epj_for_force_[ith].size();
        CopyPjForForceST(adr_epj_for_force_[ith], epj_sorted_, n_head, n_tail, epj_for_force_[ith]);
    }


    ///// LONG

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForce(Tfunc_ep_ep pfunc_ep_ep,
              Tfunc_ep_sp pfunc_ep_sp,
              const bool clear){
        const F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp, n_interaction_ep_sp_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                const S32 ith = Comm::getThreadNum();
                epj_for_force_[ith].clearSize();
                spj_for_force_[ith].clearSize();                
                adr_epj_for_force_[ith].clearSize();
                adr_spj_for_force_[ith].clearSize();
                makeInteractionList(i);
                ni_tmp += ipg_[i].n_ptcl_;
                nj_tmp += epj_for_force_[Comm::getThreadNum()].size();
                nj_tmp += spj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * epj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_sp_tmp += ipg_[i].n_ptcl_ * spj_for_force_[Comm::getThreadNum()].size();
                calcForceOnly( pfunc_ep_ep, pfunc_ep_sp, i, clear);
            }
            ni_ave_ = ni_tmp / n_ipg;
            nj_ave_ = nj_tmp / n_ipg;
            n_interaction_ep_ep_ = n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_ = n_interaction_ep_sp_tmp;
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
        }
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
        //PROFILE::Stop(profile.calc_force);
        copyForceOriginalOrder();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
        std::cout<<"force_sorted_.size()="<<force_sorted_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceOnly(Tfunc_ep_ep pfunc_ep_ep,
                  Tfunc_ep_sp pfunc_ep_sp,
                  const S32 adr_ipg,
                  const bool clear){
        const S32 offset = ipg_[adr_ipg].adr_ptcl_;
        const S32 n_epi = ipg_[adr_ipg].n_ptcl_;
        const S32 ith = Comm::getThreadNum();
        const S32 n_epj = epj_for_force_[ith].size();
        const S32 n_spj = spj_for_force_[ith].size();
        const S32 n_tail = offset + n_epi;
        if(clear){
            for(S32 i=offset; i<n_tail; i++) force_sorted_[i].clear();
        }
        pfunc_ep_ep(epi_sorted_.getPointer(offset),     n_epi,
                    epj_for_force_[ith].getPointer(),   n_epj,
                    force_sorted_.getPointer(offset));
        pfunc_ep_sp(epi_sorted_.getPointer(offset),     n_epi,
                    spj_for_force_[ith].getPointer(),   n_spj,
                    force_sorted_.getPointer(offset));
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagForceLong,
                            const S32 adr_ipg,
                            const bool clear){
        const S32 ith = Comm::getThreadNum();
        if (clear) {
            epj_for_force_[ith].clearSize();
            spj_for_force_[ith].clearSize();
            adr_epj_for_force_[ith].clearSize();
            adr_spj_for_force_[ith].clearSize();
        }
        if (theta_ > 0.0){
            const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
            S32 adr_tree_sp_first = spj_sorted_.size() - tc_glb_.size();
            S32 adr_tc = 0;
            TargetBox<TSM> target_box;
            target_box.set(ipg_[adr_ipg]);
            S32 n_epj_head = adr_epj_for_force_[ith].size();
            S32 n_spj_head = adr_spj_for_force_[ith].size();
            MakeListUsingTreeRecursiveTop
                <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                 WALK_MODE_NORMAL, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl>
                (tc_glb_,  adr_tc, tp_glb_,
                 epj_sorted_, adr_epj_for_force_[ith],
                 spj_sorted_, adr_spj_for_force_[ith],
                 target_box,
                 r_crit_sq, n_leaf_limit_,
                 adr_tree_sp_first, F64vec(0.0));
            S32 n_epj_tail = adr_epj_for_force_[ith].size();
            S32 n_spj_tail = adr_spj_for_force_[ith].size();
            CopyPjForForceST(adr_epj_for_force_[ith], epj_sorted_, n_epj_head, n_epj_tail, epj_for_force_[ith]);
            CopyPjForForceST(adr_spj_for_force_[ith], spj_sorted_, n_spj_head, n_spj_tail, spj_for_force_[ith]);
        } else {
            // theta_ = 0 case
            makeInteractionListLongForZeroTheta(typename TraitsForCutoff<typename TSM::force_type>::type_cutoff(), adr_ipg);
        }
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                    const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_sorted_[i].clear();
            }
        }
        S64 n_interaction_ep_ep_tmp = 0;
        const S64 n_ipg = ipg_.size();
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_interaction_ep_ep_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                const S32 ith = Comm::getThreadNum();
                const S32 n_epi = ipg_[i].n_ptcl_;
                const S32 adr_epi_head = ipg_[i].adr_ptcl_;
                const S32 n_epj = interaction_list_.n_ep_[i];
                const S32 adr_epj_head = interaction_list_.n_disp_ep_[i];
                const S32 adr_epj_end  = interaction_list_.n_disp_ep_[i+1];
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * n_epj;
                epj_for_force_[ith].resizeNoInitialize(n_epj);
                S32 n_cnt = 0;
                for(S32 j=adr_epj_head; j<adr_epj_end; j++, n_cnt++){
                    const S32 adr_epj = interaction_list_.adr_ep_[j];
                    epj_for_force_[ith][n_cnt] = epj_sorted_[adr_epj];
                }
                pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            epj_for_force_[ith].getPointer(),   n_epj,
                            force_sorted_.getPointer(adr_epi_head));
            }
        }
        n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        copyForceOriginalOrder();
        time_profile_.calc_force += GetWtime() - time_offset;
    }

    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                    Tfunc_ep_sp pfunc_ep_sp,
                    const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_sorted_[i].clear();
            }
        }
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        const S64 n_ipg = ipg_.size();
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_interaction_ep_ep_tmp, n_interaction_ep_sp_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                const S32 ith = Comm::getThreadNum();
                const S32 n_epi = ipg_[i].n_ptcl_;
                const S32 adr_epi_head = ipg_[i].adr_ptcl_;
                const S32 n_epj = interaction_list_.n_ep_[i];
                const S32 adr_epj_head = interaction_list_.n_disp_ep_[i];
                const S32 adr_epj_end  = interaction_list_.n_disp_ep_[i+1];
                const S32 n_spj = interaction_list_.n_sp_[i];
                const S32 adr_spj_head = interaction_list_.n_disp_sp_[i];
                const S32 adr_spj_end  = interaction_list_.n_disp_sp_[i+1];
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * n_epj;
                n_interaction_ep_sp_tmp += ipg_[i].n_ptcl_ * n_spj;
                epj_for_force_[ith].resizeNoInitialize(n_epj);
                spj_for_force_[ith].resizeNoInitialize(n_spj);
                S32 n_ep_cnt = 0;
                for(S32 j=adr_epj_head; j<adr_epj_end; j++, n_ep_cnt++){
                    const S32 adr_epj = interaction_list_.adr_ep_[j];
                    epj_for_force_[ith][n_ep_cnt] = epj_sorted_[adr_epj];
                }
                pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            epj_for_force_[ith].getPointer(),   n_epj,
                            force_sorted_.getPointer(adr_epi_head));
                S32 n_sp_cnt = 0;
                for(S32 j=adr_spj_head; j<adr_spj_end; j++, n_sp_cnt++){
                    const S32 adr_spj = interaction_list_.adr_sp_[j];
                    spj_for_force_[ith][n_sp_cnt] = spj_sorted_[adr_spj];
                }
                pfunc_ep_sp(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            spj_for_force_[ith].getPointer(),   n_spj,
                            force_sorted_.getPointer(adr_epi_head));
            }
        }
        n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
        copyForceOriginalOrder();
        time_profile_.calc_force += GetWtime() - time_offset;
    }
    
    
}



