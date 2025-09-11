namespace ParticleSimulator{
    using MomentMonopoleInAndOut = MomentMonopole;
    using MomentMonopoleSymmetry = MomentMonopole;
    using MomentMonopoleScatter  = MomentMonopole;
    using MomentMonopoleCutoff   = MomentMonopole;
    using MomentQuadrupoleScatter = MomentQuadrupole;
    using MomentSearchInAndOut    = MomentShort;
    using MomentSearchInOnly      = MomentShort;
    using SPJMonopoleInAndOut    = SPJMonopole;
    using SPJQuadrupoleInAndOut  = SPJQuadrupole;
    using SPJMonopoleScatter     = SPJMonopole;
    using SPJQuadrupoleScatter   = SPJQuadrupole;
    using SPJMonopoleCutoff      = SPJMonopole;


    //////////////////////////////////////////////////////////////
    //////////// Walk+Force, Kernel:Ptcl, List:  ////////////    
    // SHORT
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcForce_before_v8(Tfunc_ep_ep pfunc_ep_ep, const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        //n_walk_local_ += n_ipg;
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


    //////////
    ///// LONG
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcForce_before_v8(Tfunc_ep_ep pfunc_ep_ep,
                        Tfunc_ep_sp pfunc_ep_sp,
                        const bool clear){
        const F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
        if(n_ipg > 0){
#ifdef LOOP_TREE
            const auto n_thread = Comm::getNumberOfThread();
            const F64 inv_theta_sq = 1.0 / (theta_*theta_);
            constexpr S32 size_vec_limit = SIMD_VEC_LEN*2;
            ReallocatableArray<S32> adr_epj[N_THREAD_LIMIT][size_vec_limit];
            ReallocatableArray<S32> adr_spj[N_THREAD_LIMIT][size_vec_limit];
            for(S32 i=0; i<n_thread; i++){
                for(S32 j=0; j<size_vec_limit; j++){
                    adr_epj[i][j].reserve(100000);
                    adr_spj[i][j].reserve(100000);
                }
            }
            TargetBox<TSM> target_box[N_THREAD_LIMIT][size_vec_limit];
            S32 n_epj[N_THREAD_LIMIT][size_vec_limit];
            S32 n_spj[N_THREAD_LIMIT][size_vec_limit];
            S32 n_disp_epj[N_THREAD_LIMIT][size_vec_limit+1];
            S32 n_disp_spj[N_THREAD_LIMIT][size_vec_limit+1];
            const S32 adr_tree_sp_first = n_let_sp_;
            const S32 n_n_ipg = ((n_ipg-1) / size_vec_limit)+1;
            //std::cerr<<"n_ipg= "<<n_ipg<<" n_n_ipg= "<<n_n_ipg<<std::endl;
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_interaction_ep_ep_tmp, n_interaction_ep_sp_tmp)
            for(S32 i=0; i<n_n_ipg; i++){
                const S32 ith = Comm::getThreadNum();
                const S32 ipg_head = i*size_vec_limit;
                const S32 size_vec = std::min( (n_ipg-ipg_head), size_vec_limit);
                const S32 ipg_end  = ipg_head + size_vec;
                //std::cerr<<"ith= "<<ith<<" size_vec= "<<size_vec<<" ipg_head= "<<ipg_head<<std::endl;
                for(S32 j=0; j<size_vec; j++){
                    const S32 adr_ipg = ipg_head + j;
                    target_box[ith][j].set(ipg_[adr_ipg]);
                    adr_epj[ith][j].clearSize();
                    adr_spj[ith][j].clearSize();
                    //if(ith==0){
                    //std::cerr<<"check A"<<std::endl;
                    //std::cerr<<"adr_epj[ith][j].size()= "<<adr_epj[ith][j].size()<<std::endl;
                    //std::cerr<<"adr_spj[ith][j].size()= "<<adr_spj[ith][j].size()<<std::endl;
                        //}
                }
                MakeListUsingTreeLoop(tc_glb_, tp_glb_, target_box[ith], adr_epj[ith], adr_spj[ith], size_vec, n_leaf_limit_, adr_tree_sp_first, inv_theta_sq);
                for(S32 j=0; j<size_vec; j++){
                    //if(ith==0){
                    //std::cerr<<"check B"<<std::endl;
                    //std::cerr<<"adr_epj[ith][j].size()= "<<adr_epj[ith][j].size()<<std::endl;
                    //std::cerr<<"adr_spj[ith][j].size()= "<<adr_spj[ith][j].size()<<std::endl;
                        //}
                    const S32 ith = Comm::getThreadNum();
                    const S32 n_epj_head = 0;
                    const S32 n_epj_tail = adr_epj[ith][j].size();
                    const S32 n_epj = n_epj_tail - n_epj_head;
                    const S32 n_spj_head = 0;
                    const S32 n_spj_tail = adr_spj[ith][j].size();
                    const S32 n_spj = n_spj_tail - n_spj_head;
                    const S32 adr_ipg = ipg_head + j;
                    n_interaction_ep_ep_tmp += ipg_[adr_ipg].n_ptcl_ * n_epj;
                    n_interaction_ep_sp_tmp += ipg_[adr_ipg].n_ptcl_ * n_spj;
                    epj_for_force_[ith].resizeNoInitialize(n_epj);
                    spj_for_force_[ith].resizeNoInitialize(n_spj);
                    CopyPjForForceST(adr_epj[ith][j], epj_sorted_, n_epj_head, n_epj_tail, epj_for_force_[ith]);
                    CopyPjForForceST(adr_spj[ith][j], spj_sorted_, tc_glb_, adr_tree_sp_first, n_spj_head, n_spj_tail, spj_for_force_[ith]);
                    F64 cm_mass = 0.0;
                    for(auto j=0; j<n_epj; j++){
                        cm_mass += epj_for_force_[ith][j].getCharge();
                    }
                    for(auto j=0; j<n_spj; j++){
                        cm_mass += spj_for_force_[ith][j].getCharge();
                    }
                    assert(cm_mass == 1.0);
                    //std::cerr<<"cm_mass= "<<cm_mass<<std::endl;
                    const S32 offset = ipg_[adr_ipg].adr_ptcl_;
                    const S32 n_epi = ipg_[adr_ipg].n_ptcl_;
                    if(clear){
                        for(S32 i=offset; i<offset+n_epi; i++) force_sorted_[i].clear();
                    }
                    pfunc_ep_ep(epi_sorted_.getPointer(offset),     n_epi,
                                epj_for_force_[ith].getPointer(),   n_epj,
                                force_sorted_.getPointer(offset));
                    pfunc_ep_sp(epi_sorted_.getPointer(offset),     n_epi,
                                spj_for_force_[ith].getPointer(),   n_spj,
                                force_sorted_.getPointer(offset));
                }
            }
            //sleep(1);
            //Abort();
            n_interaction_ep_ep_ = n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_ = n_interaction_ep_sp_tmp;
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
#else
            auto n_thread = Comm::getNumberOfThread();
            F64 maketreetime[n_thread];
            F64 calcforcetime[n_thread];
            for(auto i=0; i< n_thread; i++){
                maketreetime[i]=0;
                calcforcetime[i]=0;
            }
            auto wtime_start=GetWtime();
            //	    std::cerr << "Start parallel force calc\n";
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic,4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp, n_interaction_ep_sp_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
#ifdef PARTICLE_SIMULATOR_TREEWALK_DETAILED_TIMING		
                F64 wtime_offset = GetWtime();
#endif		
                const S32 ith = Comm::getThreadNum();
                epj_for_force_[ith].clearSize();
                spj_for_force_[ith].clearSize();                
                adr_epj_for_force_[ith].clearSize();
                adr_spj_for_force_[ith].clearSize();
                makeInteractionList(i);
#ifdef PARTICLE_SIMULATOR_TREEWALK_DETAILED_TIMING		
                F64 wtime_offset2 = GetWtime();
                maketreetime[ith] += wtime_offset2 - wtime_offset;
#endif		
                ni_tmp += ipg_[i].n_ptcl_;
                nj_tmp += epj_for_force_[Comm::getThreadNum()].size();
                nj_tmp += spj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * epj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_sp_tmp += ipg_[i].n_ptcl_ * spj_for_force_[Comm::getThreadNum()].size();
                calcForceOnly( pfunc_ep_ep, pfunc_ep_sp, i, clear);
#ifdef PARTICLE_SIMULATOR_TREEWALK_DETAILED_TIMING		
                calcforcetime[ith] += GetWtime() - wtime_offset2;
#endif
            }
            F64 makemax=0;
            F64 calcmax=0;
#ifdef PARTICLE_SIMULATOR_TREEWALK_DETAILED_TIMING		
            for(S32 i=0; i<n_thread; i++){
                if (makemax < maketreetime[i])makemax = maketreetime[i];
                if (calcmax < calcforcetime[i])calcmax = calcforcetime[i];
            }
#else
            calcmax = GetWtime()-wtime_start;
#endif	    
            time_profile_.calc_force__core__walk_tree  += makemax;
            time_profile_.calc_force__core += calcmax;
            ni_ave_ = ni_tmp / n_ipg;
            nj_ave_ = nj_tmp / n_ipg;
            n_interaction_ep_ep_ = n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_ = n_interaction_ep_sp_tmp;
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
#endif
        }
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcForceOnly_before_v8(Tfunc_ep_ep pfunc_ep_ep,
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    makeInteractionListImpl_before_v8(TagForceShort, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        if (clear){
            epj_for_force_[ith].clearSize();
        }
        //const F64 r_crit_sq = 9999.9;
        TargetBox<TSM> target_box;
        target_box.set(ipg_[adr_ipg]);
        S32 adr_tc = 0;
        S32 n_head = adr_epj_for_force_[ith].size();
        const auto len_peri = p_domain_info_->getLenRootDomain();
        MakeListUsingTreeRecursiveTop
            <TSM, TreeCellGlb, TreeParticle, Tepj, TagChopLeafTrue, CALC_DISTANCE_TYPE>
            (tc_glb_,  adr_tc, tp_glb_,
             epj_sorted_, adr_epj_for_force_[ith],
             target_box,
             n_leaf_limit_,
             len_peri);
        S32 n_tail = adr_epj_for_force_[ith].size();
        CopyPjForForceST(adr_epj_for_force_[ith], epj_sorted_, n_head, n_tail, epj_for_force_[ith]);
    }



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcForceOnly_before_v8(Tfunc_ep_ep pfunc_ep_ep,
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    makeInteractionListImpl_before_v8(TagForceLong,
                            const S32 adr_ipg,
                            const bool clear){
        const S32 ith = Comm::getThreadNum();
	//	std::cerr<<"makeInteractionListImpl(TagForceLong, ...\n";
        if (clear) {
            epj_for_force_[ith].clearSize();
            spj_for_force_[ith].clearSize();
            adr_epj_for_force_[ith].clearSize();
            adr_spj_for_force_[ith].clearSize();
        }
        if (theta_ > 0.0){
            //S32 adr_tree_sp_first = spj_sorted_.size() - tc_glb_.size();
            const auto adr_tree_sp_first = n_let_sp_;
            S32 adr_tc = 0;
            TargetBox<TSM> target_box;
            target_box.set(ipg_[adr_ipg]);
            S32 n_epj_head = adr_epj_for_force_[ith].size();
            S32 n_spj_head = adr_spj_for_force_[ith].size();
            const auto len_peri = p_domain_info_->getLenRootDomain();
            MakeListUsingTreeRecursiveTop
                <TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, TagChopLeafTrue,
                 TagCopyInfoCloseWithTpAdrptcl, CALC_DISTANCE_TYPE>
                (tc_glb_,  adr_tc, tp_glb_,
                 epj_sorted_, adr_epj_for_force_[ith],
                 spj_sorted_, adr_spj_for_force_[ith],
                 target_box,
                 n_leaf_limit_,
                 adr_tree_sp_first, len_peri, theta_);
            S32 n_epj_tail = adr_epj_for_force_[ith].size();
            S32 n_spj_tail = adr_spj_for_force_[ith].size();
            CopyPjForForceST(adr_epj_for_force_[ith], epj_sorted_, n_epj_head, n_epj_tail, epj_for_force_[ith]);
            CopyPjForForceST(adr_spj_for_force_[ith], spj_sorted_, tc_glb_, n_let_sp_, n_spj_head, n_spj_tail, spj_for_force_[ith]);
        } else {
            makeInteractionListLongForZeroTheta(typename TraitsForCutoff<typename TSM::search_type>::type_cutoff(), adr_ipg);
        }
    }
    
}
