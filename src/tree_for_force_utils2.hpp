namespace ParticleSimulator{

    ///////////////////////////////
    // small functions for dispatch
    
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchLong,
                                           const ReallocatableArray<Ttc> & tc_first){
        return F64ort(-1234.5, -9876.5);
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchLongCutoff,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchLongScatter,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchLongSymmetry,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchShortGather,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchShortSymmetry,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchShortScatter,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }

    inline void ExchangeOuterBoundary(TagSearchLong,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){/* do nothing */ }
    inline void ExchangeOuterBoundary(TagSearchLongCutoff,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){/* do nothing */ }
    inline void ExchangeOuterBoundary(TagSearchLongScatter,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){/* do nothing */ }
    inline void ExchangeOuterBoundary(TagSearchLongSymmetry,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
        Comm::allGather(my_outer, 1, outer.getPointer());
    }
    inline void ExchangeOuterBoundary(TagSearchShortScatter,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){/* do nothing */ }
    inline void ExchangeOuterBoundary(TagSearchShortGather,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){/* do nothing */ }
    inline void ExchangeOuterBoundary(TagSearchShortSymmetry,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){/* do nothing */ }

    ////////////////////
    // for exchange LET
    // only for open boundary
    template<class Tptcl, class Ttc>
    inline void CopyPtclFromTreeToSendBuf(TagSearchBoundaryConditionOpenOnly,
                                          ReallocatableArray<Tptcl> & ptcl_send,
                                          const ReallocatableArray<Ttc> & tc,
                                          const ReallocatableArray<S32> & adr_ptcl_send,
                                          const S32 n_ptcl,
                                          const S32 n_ptcl_offset,
                                          const F64vec & shift){
        for(S32 j=0; j<n_ptcl; j++){
            S32 adr = adr_ptcl_send[n_ptcl_offset+j];
            ptcl_send[n_ptcl_offset+j].copyFromMoment(tc[adr].mom_);
        }
    }

    // for periodic boundary (not only for periodic but also open boudnary)
    template<class Tptcl, class Ttc>
    inline void CopyPtclFromTreeToSendBuf(TagSearchBoundaryConditionOpenPeriodic,
                                          ReallocatableArray<Tptcl> & ptcl_send,
                                          const ReallocatableArray<Ttc> & tc,
                                          const ReallocatableArray<S32> & adr_ptcl_send,
                                          const S32 n_ptcl,
                                          const S32 n_ptcl_offset,
                                          const F64vec & shift){
        for(S32 j=0; j<n_ptcl; j++){
            S32 adr = adr_ptcl_send[n_ptcl_offset+j];
            ptcl_send[n_ptcl_offset+j].copyFromMoment(tc[adr].mom_);
            const F64vec pos_new = ptcl_send[n_ptcl_offset+j].getPos() - shift;
            ptcl_send[n_ptcl_offset+j].setPos(pos_new);
        }
    }
    template<class Tptcl>
    inline void CopyPtclToSendBuf(TagSearchBoundaryConditionOpenOnly,
                                  ReallocatableArray<Tptcl> & ptcl_send,
                                  const ReallocatableArray<Tptcl> & ptcl,
                                  const ReallocatableArray<S32> & adr_ptcl_send,
                                  const S32 n_ptcl,
                                  const S32 n_ptcl_offset,
                                  const F64vec & shift){
        for(S32 j=0; j<n_ptcl; j++){
            S32 adr = adr_ptcl_send[n_ptcl_offset+j];
            ptcl_send[n_ptcl_offset+j] = ptcl[adr];
        }
    }
    template<class Tptcl>
    inline void CopyPtclToSendBuf(TagSearchBoundaryConditionOpenPeriodic,
                                  ReallocatableArray<Tptcl> & ptcl_send,
                                  const ReallocatableArray<Tptcl> & ptcl,
                                  const ReallocatableArray<S32> & adr_ptcl_send,
                                  const S32 n_ptcl,
                                  const S32 n_ptcl_offset,
                                  const F64vec & shift){
        for(S32 j=0; j<n_ptcl; j++){
            S32 adr = adr_ptcl_send[n_ptcl_offset+j];
            ptcl_send[n_ptcl_offset+j] = ptcl[adr];
            ptcl_send[n_ptcl_offset+j].setPos(ptcl[adr].getPos() - shift);
        }
    }
    // for exchange LET
    ////////////////////
    
    // small functions for dispatch
    ///////////////////////////////


    // for long mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Twalkmode>
    inline void FindScatterParticle(const ReallocatableArray<Ttc> & tc_first,
                                    const ReallocatableArray<Ttp> & tp_first,
                                    const ReallocatableArray<Tep> & ep_first,
                                    ReallocatableArray<S32> & n_ep_send,
                                    ReallocatableArray<S32> & adr_ep_send,
                                    const DomainInfo & dinfo,
                                    const S32 n_leaf_limit,
                                    ReallocatableArray<S32> & n_sp_send,
                                    ReallocatableArray<S32> & adr_sp_send,
                                    ReallocatableArray<F64vec> & shift_per_image,
                                    ReallocatableArray<S32> & n_image_per_proc,
                                    ReallocatableArray<S32> & n_ep_per_image,
                                    ReallocatableArray<S32> & n_sp_per_image,
                                    const F64 r_crit_sq){
        const S32 n_thread = Comm::getNumberOfThread();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();        
        static bool first = true;
        static ReallocatableArray<S32> * rank_tmp;
        static ReallocatableArray<F64vec> * shift_per_image_tmp;
        static ReallocatableArray<S32> * adr_ep_send_tmp;
        static ReallocatableArray<S32> * n_ep_per_image_tmp;
        static ReallocatableArray<S32> * adr_sp_send_tmp;
        static ReallocatableArray<S32> * n_sp_per_image_tmp;
        static ReallocatableArray<F64ort> outer_boundary_of_tree;
        if(first){
            rank_tmp = new ReallocatableArray<S32>[n_thread];
            shift_per_image_tmp = new ReallocatableArray<F64vec>[n_thread];
            adr_ep_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_ep_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            adr_sp_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_sp_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        outer_boundary_of_tree.resizeNoInitialize(n_proc);
        n_ep_send.resizeNoInitialize(n_proc);
        n_sp_send.resizeNoInitialize(n_proc);
        n_image_per_proc.resizeNoInitialize(n_proc);
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        F64ort pos_root_domain = dinfo.getPosRootDomain();
        F64vec len_peri = pos_root_domain.getFullLength();
        for(S32 i=0; i<DIMENSION; i++){
            if(pa[i]==false) len_peri[i] = 0.0;
        }
        F64ort outer_boundary_of_my_tree = GetOuterBoundaryOfMyTree(typename TSM::search_type(), tc_first);
        S32 adr_tc = 0;
        S32 adr_tree_sp_first = 0;
        ExchangeOuterBoundary(typename TSM::search_type(), &outer_boundary_of_my_tree, outer_boundary_of_tree);
        // for long symmetry
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            S32 ith = Comm::getThreadNum();
            rank_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ep_send_tmp[ith].clearSize();
            n_ep_per_image_tmp[ith].clearSize();
            adr_sp_send_tmp[ith].clearSize();
            n_sp_per_image_tmp[ith].clearSize();
            ReallocatableArray<Tsp> sp_first;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 i=0; i<n_proc; i++){
                const S32 n_image_tmp_prev = shift_per_image_tmp[ith].size();
                rank_tmp[ith].push_back(i);
                CalcNumberAndShiftOfImageDomain
                    (shift_per_image_tmp[ith],  dinfo.getPosRootDomain().getFullLength(),
                     outer_boundary_of_my_tree, dinfo.getPosDomain(i), pa, false);
                const S32 n_image_tmp = shift_per_image_tmp[ith].size();
                n_image_per_proc[i] = n_image_tmp - n_image_tmp_prev;
                S32 n_ep_prev = adr_ep_send_tmp[ith].size();
                S32 n_sp_prev = adr_sp_send_tmp[ith].size();
                for(S32 j=n_image_tmp_prev; j<n_image_tmp; j++){
                    S32 n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                    S32 n_sp_prev_2 = adr_sp_send_tmp[ith].size();
                    if(my_rank==i && j==n_image_tmp_prev){
                        n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2); // is 0
                        n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2); // is 0
                        continue;
                    }
                    F64ort pos_target_domain = dinfo.getPosDomain(i).shift(shift_per_image_tmp[ith][j]);
                    TargetBox<TSM> target_box;
                    SetTargetBoxExLet(target_box, pos_target_domain, outer_boundary_of_tree[i].shift(shift_per_image_tmp[ith][j]));
                    MakeListUsingTreeRecursiveTop
                        <TSM, Ttc, TreeParticle, Tep, Tsp, Twalkmode, TagChopLeafFalse, TagCopyInfoCloseNoSp>
                        (tc_first, adr_tc, tp_first,
                         ep_first, adr_ep_send_tmp[ith],
                         sp_first, adr_sp_send_tmp[ith],
                         target_box,
                         r_crit_sq, n_leaf_limit,
                         adr_tree_sp_first, F64vec(0.0));
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                    n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
                }
                n_ep_send[i] = adr_ep_send_tmp[ith].size() - n_ep_prev;
                n_sp_send[i] = adr_sp_send_tmp[ith].size() - n_sp_prev;
            }
        } // end of OMP scope

        ReallocatableArray<S32> n_disp_image_per_proc;
        n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
        n_disp_image_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_per_proc[i+1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
        }
        const S32 n_image_tot = n_disp_image_per_proc[n_proc];
        shift_per_image.resizeNoInitialize( n_image_tot );
        n_ep_per_image.resizeNoInitialize( n_image_tot);
        n_sp_per_image.resizeNoInitialize( n_image_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt = 0;
            //S32 n_cnt_ep = 0;
            //S32 n_cnt_sp = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                S32 offset = n_disp_image_per_proc[rank];
                for(S32 k=0; k<n_image_per_proc[rank]; k++){
                    shift_per_image[offset+k] = shift_per_image_tmp[i][n_cnt];
                    n_ep_per_image[offset+k] = n_ep_per_image_tmp[i][n_cnt];
                    n_sp_per_image[offset+k] = n_sp_per_image_tmp[i][n_cnt];
                    n_cnt++;
                }
            }
        }
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1);
        ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot+1);
        n_disp_ep_per_image[0] = 0;
        n_disp_sp_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
            n_disp_sp_per_image[i+1] = n_disp_sp_per_image[i] + n_sp_per_image[i];
        }
        const S32 n_ep_send_tot = n_disp_ep_per_image[ n_image_tot ];
        const S32 n_sp_send_tot = n_disp_sp_per_image[ n_image_tot ];
        adr_ep_send.resizeNoInitialize( n_ep_send_tot );
        adr_sp_send.resizeNoInitialize( n_sp_send_tot );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_ep = 0;
            S32 n_cnt_sp = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                const S32 adr_image_head = n_disp_image_per_proc[rank];
                const S32 adr_image_end = n_disp_image_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++){
                    const S32 adr_ep_head = n_disp_ep_per_image[k];
                    const S32 adr_ep_end = n_disp_ep_per_image[k+1];
                    for(S32 l=adr_ep_head; l<adr_ep_end; l++){
                        adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                    }
                    const S32 adr_sp_head = n_disp_sp_per_image[k];
                    const S32 adr_sp_end  = n_disp_sp_per_image[k+1];
                    for(S32 l=adr_sp_head; l<adr_sp_end; l++){
                        adr_sp_send[l] = adr_sp_send_tmp[i][n_cnt_sp++];
                    }
                }
            }
        }
    }

    // for short mode
    template<class Ttc, class Ttp, class Tep, class Tsp, class Twalkmode>
    inline void FindScatterParticle(const ReallocatableArray<Ttc> & tc_first,
                                    const ReallocatableArray<Ttp> & tp_first,
                                    const ReallocatableArray<Tep> & ep_first,
                                    ReallocatableArray<S32> & n_ep_send,
                                    ReallocatableArray<S32> & adr_ep_send,
                                    const DomainInfo & dinfo,
                                    const S32 n_leaf_limit,
                                    ReallocatableArray<F64vec> & shift_per_image,
                                    ReallocatableArray<S32> & n_image_per_proc,
                                    ReallocatableArray<S32> & n_ep_per_image){
        static bool first = true;
        static const S32 n_thread = Comm::getNumberOfThread();
        static ReallocatableArray<S32> * rank_tmp;
        static ReallocatableArray<F64vec> * shift_per_image_tmp;
        static ReallocatableArray<S32> * adr_ep_send_tmp;
        static ReallocatableArray<S32> * n_ep_per_image_tmp;
        if(first){
            rank_tmp = new ReallocatableArray<S32>[n_thread];
            shift_per_image_tmp = new ReallocatableArray<F64vec>[n_thread];
            adr_ep_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_ep_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        n_ep_send.resizeNoInitialize(n_proc);
        n_image_per_proc.resizeNoInitialize(n_proc);
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        F64ort pos_root_domain = dinfo.getPosRootDomain();
        F64vec len_peri = pos_root_domain.getFullLength();
        for(S32 i=0; i<DIMENSION; i++){
            if(pa[i]==false) len_peri[i] = 0.0;
        }
        const F64ort outer_boundary_of_my_tree = tc_first[0].mom_.vertex_out_;
        S32 adr_tc = 0;
        S32 adr_tree_sp_first = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            S32 ith = Comm::getThreadNum();
            rank_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ep_send_tmp[ith].clearSize();
            n_ep_per_image_tmp[ith].clearSize();
            ReallocatableArray<S32> adr_sp_send_tmp;
            ReallocatableArray<Tsp> sp_first;
            F64 r_crit_sq = 1.0; // dummy
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 i=0; i<n_proc; i++){
                const S32 n_image_tmp_prev = shift_per_image_tmp[ith].size();
                rank_tmp[ith].push_back(i);
                CalcNumberAndShiftOfImageDomain
                    (shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(),
                     outer_boundary_of_my_tree, dinfo.getPosDomain(i), pa, false);
                const S32 n_image_tmp = shift_per_image_tmp[ith].size();
                n_image_per_proc[i] = n_image_tmp - n_image_tmp_prev;
                S32 n_ep_prev = adr_ep_send_tmp[ith].size();
                for(S32 j=n_image_tmp_prev; j<n_image_tmp; j++){
                    S32 n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                    if(my_rank==i && j==n_image_tmp_prev){
                        n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2); // is 0
                        continue;
                    }
                    F64ort pos_target_domain = dinfo.getPosDomain(i).shift(shift_per_image_tmp[ith][j]);
                    TargetBox<SEARCH_MODE_SCATTER> target_box;
                    target_box.vertex_in_ = pos_target_domain;
                    MakeListUsingTreeRecursiveTop
                        <SEARCH_MODE_SCATTER, Ttc, TreeParticle, Tep, Tsp, Twalkmode, TagChopLeafFalse, TagCopyInfoCloseNoSp>
                        (tc_first, adr_tc, tp_first,
                         ep_first, adr_ep_send_tmp[ith],
                         sp_first, adr_sp_send_tmp,
                         target_box,
                         r_crit_sq, n_leaf_limit,
                         adr_tree_sp_first, F64vec(0.0));

                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                }
                n_ep_send[i] = adr_ep_send_tmp[ith].size() - n_ep_prev;
            }
        } // end of OMP scope
        ReallocatableArray<S32> n_disp_image_per_proc;
        n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
        n_disp_image_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_per_proc[i+1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
        }
        const S32 n_image_tot = n_disp_image_per_proc[n_proc];
        shift_per_image.resizeNoInitialize( n_image_tot );
        n_ep_per_image.resizeNoInitialize( n_image_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt = 0;
            //S32 n_cnt_ep = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                S32 offset = n_disp_image_per_proc[rank];
                for(S32 k=0; k<n_image_per_proc[rank]; k++){
                    shift_per_image[offset+k] = shift_per_image_tmp[i][n_cnt];
                    n_ep_per_image[offset+k] = n_ep_per_image_tmp[i][n_cnt];
                    n_cnt++;
                }
            }
        }
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1);
        n_disp_ep_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
        }
        const S32 n_ep_send_tot = n_disp_ep_per_image[ n_image_tot ];
        adr_ep_send.resizeNoInitialize( n_ep_send_tot );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_ep = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                const S32 adr_image_head = n_disp_image_per_proc[rank];
                const S32 adr_image_end = n_disp_image_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++){
                    const S32 adr_ep_head = n_disp_ep_per_image[k];
                    const S32 adr_ep_end = n_disp_ep_per_image[k+1];
                    for(S32 l=adr_ep_head; l<adr_ep_end; l++){
                        adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                    }
                }
            }
        }
    }    
    
    //////////////
    // exchange # of LET
    inline void ExchangeNumber(const ReallocatableArray<S32> & n_ep_send,
                               ReallocatableArray<S32> & n_ep_recv,
                               const ReallocatableArray<S32> & n_sp_send,
                               ReallocatableArray<S32> & n_sp_recv){
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<S32> n_ep_sp_send(n_proc*2);
        static ReallocatableArray<S32> n_ep_sp_recv(n_proc*2);
        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send[i*2]   = n_ep_send[i];
            n_ep_sp_send[i*2+1] = n_sp_send[i];
        }
        Comm::allToAll(n_ep_sp_send.getPointer(), 2, n_ep_sp_recv.getPointer());
        n_ep_recv.resizeNoInitialize(n_proc);
        n_sp_recv.resizeNoInitialize(n_proc);
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv[i] = n_ep_sp_recv[i*2];
            n_sp_recv[i] = n_ep_sp_recv[i*2+1];
        }
    }

    inline void ExchangeNumber(const ReallocatableArray<S32> & n_ep_send,
                               ReallocatableArray<S32> & n_ep_recv){
        const S32 n_proc = Comm::getNumberOfProc();
        n_ep_recv.resizeNoInitialize(n_proc);
        Comm::allToAll(n_ep_send.getPointer(), 1, n_ep_recv.getPointer());
    }
    // exchange # of LET
    //////////////

    ///////////////
    // exchange LET
    template<class TSM, class Tep, class Tsp, class Ttc>
    inline void ExchangeLet(const ReallocatableArray<Tep> & ep,
                            const ReallocatableArray<S32> & n_ep_send,
                            const ReallocatableArray<S32> & n_ep_recv,
                            const ReallocatableArray<S32> & n_ep_per_image,
                            const ReallocatableArray<S32> & adr_ep_send,
                            ReallocatableArray<Tep> & ep_org,
                            const S32 n_ep_offset, // ep_org[n_ep_offset] = n_ep_recv[0]
                            const ReallocatableArray<Ttc> & tc,
                            const ReallocatableArray<S32> & n_sp_send,
                            const ReallocatableArray<S32> & n_sp_recv,
                            const ReallocatableArray<S32> & n_sp_per_image,
                            const ReallocatableArray<S32> & adr_sp_send,
                            ReallocatableArray<Tsp> & sp_org,
                            const ReallocatableArray<F64vec> & shift_image_domain,
                            const ReallocatableArray<S32> & n_image_per_proc){
        const S32 n_proc = Comm::getNumberOfProc();
        //const S32 my_rank = Comm::getRank();
        ReallocatableArray<S32> n_disp_ep_send(n_proc+1);
        ReallocatableArray<S32> n_disp_sp_send(n_proc+1);
        ReallocatableArray<S32> n_disp_ep_recv(n_proc+1);
        ReallocatableArray<S32> n_disp_sp_recv(n_proc+1);
        n_disp_ep_send[0] = n_disp_sp_send[0] = n_disp_ep_recv[0] = n_disp_sp_recv[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_ep_send[i+1] = n_ep_send[i] + n_disp_ep_send[i];
            n_disp_sp_send[i+1] = n_sp_send[i] + n_disp_sp_send[i];
            n_disp_ep_recv[i+1] = n_ep_recv[i] + n_disp_ep_recv[i];
            n_disp_sp_recv[i+1] = n_sp_recv[i] + n_disp_sp_recv[i];
        }
        const S32 n_image_tot = n_ep_per_image.size();
        //if(my_rank==0) std::cerr<<"n_image_tot= "<<n_image_tot<<std::endl;
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1);
        ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot+1);
        n_disp_ep_per_image[0] = n_disp_sp_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] +  n_ep_per_image[i];
            n_disp_sp_per_image[i+1] = n_disp_sp_per_image[i] +  n_sp_per_image[i];
        }
        ReallocatableArray<Tep> ep_send( n_disp_ep_send[n_proc] );
        ReallocatableArray<Tsp> sp_send( n_disp_sp_send[n_proc] );

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for(S32 i=0; i<n_image_tot; i++){
            //F64 cm_mass = 0.0;
            //F64vec cm_pos = 0.0;
            F64vec shift = shift_image_domain[i];

            S32 n_ep = n_ep_per_image[i];
            S32 n_ep_offset = n_disp_ep_per_image[i];
            CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep,
                              n_ep_offset, shift);
            //for(S32 i2=n_ep_offset; i2<n_ep_offset+n_ep; i2++){
            //    cm_mass += ep_send[i2].mass;
            //    cm_pos += ep_send[i2].mass * ep_send[i2].pos;
            //}
            
            S32 n_sp = n_sp_per_image[i];
            S32 n_sp_offset = n_disp_sp_per_image[i];
            CopyPtclFromTreeToSendBuf(typename TSM::search_boundary_type(), sp_send, tc, adr_sp_send, n_sp, n_sp_offset, shift);
            /*
            for(S32 i2=n_sp_offset; i2<n_sp_offset+n_sp; i2++){
                cm_mass += sp_send[i2].mass;
                cm_pos += sp_send[i2].mass * sp_send[i2].pos;
            }
            if(Comm::getRank()==0){
                std::cerr<<"i= "<<i
                         <<" cm_mass= "<<cm_mass
                         <<" cm_pos= "<<cm_pos / cm_mass
                         <<std::endl;
            }
            */

        }
        //exit(1);
        ep_org.resizeNoInitialize( n_disp_ep_recv[n_proc] + n_ep_offset);
        sp_org.resizeNoInitialize( n_disp_sp_recv[n_proc] );
        Comm::allToAllV(ep_send.getPointer(), n_ep_send.getPointer(), n_disp_ep_send.getPointer(),
                        ep_org.getPointer(n_ep_offset), n_ep_recv.getPointer(), n_disp_ep_recv.getPointer());
        Comm::allToAllV(sp_send.getPointer(), n_sp_send.getPointer(), n_disp_sp_send.getPointer(),
                        sp_org.getPointer(), n_sp_recv.getPointer(), n_disp_sp_recv.getPointer());
    }
    
    template<class Tep>
    inline void ExchangeLet(const ReallocatableArray<Tep> & ep,
                            const ReallocatableArray<S32> & n_ep_send,
                            const ReallocatableArray<S32> & n_ep_recv,
                            const ReallocatableArray<S32> & n_ep_per_image,
                            const ReallocatableArray<S32> & adr_ep_send,
                            ReallocatableArray<Tep> & ep_org,
                            const S32 n_ep_offset, 
                            const ReallocatableArray<F64vec> & shift_image_domain,
                            const ReallocatableArray<S32> & n_image_per_proc){
        const S32 n_proc = Comm::getNumberOfProc();
        ReallocatableArray<S32> n_disp_ep_send(n_proc+1);
        ReallocatableArray<S32> n_disp_ep_recv(n_proc+1);
        n_disp_ep_send[0] = n_disp_ep_recv[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_ep_send[i+1] = n_ep_send[i] + n_disp_ep_send[i];
            n_disp_ep_recv[i+1] = n_ep_recv[i] + n_disp_ep_recv[i];
        }
        const S32 n_image_tot = n_ep_per_image.size();
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1);
        n_disp_ep_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] +  n_ep_per_image[i];
        }
        ReallocatableArray<Tep> ep_send( n_disp_ep_send[n_proc] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for(S32 i=0; i<n_image_tot; i++){
            F64vec shift = shift_image_domain[i];
            S32 n_ep = n_ep_per_image[i];
            S32 n_ep_offset = n_disp_ep_per_image[i];
            for(S32 j=0; j<n_ep; j++){
                S32 adr = adr_ep_send[n_ep_offset+j];
                ep_send[n_ep_offset+j] = ep[adr];
                ep_send[n_ep_offset+j].setPos(ep[adr].getPos() - shift);
            }
        }
        ep_org.resizeNoInitialize( n_disp_ep_recv[n_proc] + n_ep_offset);
        Comm::allToAllV(ep_send.getPointer(), n_ep_send.getPointer(), n_disp_ep_send.getPointer(),
                        ep_org.getPointer(n_ep_offset), n_ep_recv.getPointer(), n_disp_ep_recv.getPointer());
    }
    // exchange LET
    //////////////

    template<class Tep>
    inline void ExchangeParticle(const ReallocatableArray<Tep> & ep,
                                 const ReallocatableArray<S32> & n_ep_send,
                                 const ReallocatableArray<S32> & n_ep_recv,
                                 const ReallocatableArray<S32> & n_ep_per_image,
                                 const ReallocatableArray<S32> & adr_ep_send,
                                 ReallocatableArray<Tep> & ep_recv,
                                 const ReallocatableArray<F64vec> & shift_image_domain,
                                 const ReallocatableArray<S32> & n_image_per_proc){
        const S32 n_proc = Comm::getNumberOfProc();
        ReallocatableArray<S32> n_disp_ep_send(n_proc+1);
        ReallocatableArray<S32> n_disp_ep_recv(n_proc+1);
        n_disp_ep_send[0] = n_disp_ep_recv[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_ep_send[i+1] = n_ep_send[i] + n_disp_ep_send[i];
            n_disp_ep_recv[i+1] = n_ep_recv[i] + n_disp_ep_recv[i];
        }
        const S32 n_image_tot = n_ep_per_image.size();
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1);
        n_disp_ep_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] +  n_ep_per_image[i];
        }
        ReallocatableArray<Tep> ep_send( n_disp_ep_send[n_proc] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for(S32 i=0; i<n_image_tot; i++){
            F64vec shift = shift_image_domain[i];
            S32 n_ep = n_ep_per_image[i];
            S32 n_ep_offset = n_disp_ep_per_image[i];
            for(S32 j=0; j<n_ep; j++){
                S32 adr = adr_ep_send[n_ep_offset+j];
                ep_send[n_ep_offset+j] = ep[adr];
                ep_send[n_ep_offset+j].setPos(ep[adr].getPos() - shift);
            }
        }
        ep_recv.resizeNoInitialize( n_disp_ep_recv[n_proc] );
        Comm::allToAllV(ep_send.getPointer(), n_ep_send.getPointer(), n_disp_ep_send.getPointer(),
                        ep_recv.getPointer(), n_ep_recv.getPointer(), n_disp_ep_recv.getPointer());
    }

    template<class Ttc, class Ttp, class Tepj>
    inline void FindExchangeParticleDoubleWalk(const ReallocatableArray<Tepj> & epj_A, // received particles
                                               const ReallocatableArray<Ttc> & tc_first_B,
                                               const ReallocatableArray<S32> & n_epj_src_per_proc,
                                               const ReallocatableArray<S32> & n_image_send_per_proc_irnai, // not needed
                                               const DomainInfo & dinfo,
                                               const S32 n_leaf_limit_B,
                                               ReallocatableArray<S32> & n_epj_send_per_proc,
                                               ReallocatableArray<S32> & n_epj_send_per_image,
                                               ReallocatableArray<S32> & n_image_send_per_proc,
                                               ReallocatableArray<S32> & adr_ep_send,
                                               ReallocatableArray<F64vec> & shift_per_image,
                                               const ReallocatableArray<Tepj> & epj_B, // assigned
                                               const F64vec & center_tree,
                                               const F64    & full_len_tree
                                               ){
        const S32 n_proc = Comm::getNumberOfProc();
        //const S32 my_rank = Comm::getRank();
        /*
        if(my_rank==0){
            std::cerr<<"epj_A.size()= "<<epj_A.size()
                     <<" tc_first_B.size()= "<<tc_first_B.size()
                     <<std::endl;
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"n_epj_src_per_proc[i]= "<<n_epj_src_per_proc[i]<<std::endl;
            }
        }
        */
        const S32 n_thread = Comm::getNumberOfThread();
        const F64ort pos_root_domain = dinfo.getPosRootDomain();
        const F64vec len_root_domain = pos_root_domain.getFullLength();
        n_image_send_per_proc.resizeNoInitialize(n_proc);
        n_epj_send_per_proc.resizeNoInitialize(n_proc);
        ReallocatableArray<S32> n_disp_epj_per_proc;
        n_disp_epj_per_proc.resizeNoInitialize(n_proc+1);
        n_disp_epj_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_epj_per_proc[i+1] = n_disp_epj_per_proc[i] + n_epj_src_per_proc[i];
        }
        for(S32 i=0; i<n_proc; i++){
            n_epj_send_per_proc[i] = 0;
        }
        
        static ReallocatableArray<Tepj> * epj_sorted_tmp;
        static ReallocatableArray<F64vec> * shift_per_image_tmp;
        static ReallocatableArray<Ttp> * tp_tmp;
        static ReallocatableArray<Ttc> * tc_tmp;
        static ReallocatableArray<S32> * adr_tc_level_partition_tmp;
        static ReallocatableArray<S32> * adr_ptcl_send_tmp;
        static ReallocatableArray<S32> * rank_dst_tmp;
        static ReallocatableArray<S32> * adr_epj_src_per_image_tmp;
        static ReallocatableArray<S32> * n_epj_src_per_image_tmp;
        static ReallocatableArray<S32> * n_epj_send_per_image_tmp;
        static bool first = true;
        if(first){
            epj_sorted_tmp = new ReallocatableArray<Tepj>[n_thread];
            shift_per_image_tmp = new ReallocatableArray<F64vec>[n_thread];
            tp_tmp = new ReallocatableArray<Ttp>[n_thread];
            tc_tmp = new ReallocatableArray<Ttc>[n_thread];
            adr_tc_level_partition_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ptcl_send_tmp = new ReallocatableArray<S32>[n_thread];
            for(S32 i=0; i<n_thread; i++){
                adr_tc_level_partition_tmp[i].resizeNoInitialize(TREE_LEVEL_LIMIT+2);
            }
            rank_dst_tmp = new ReallocatableArray<S32>[n_thread];
            adr_epj_src_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            n_epj_src_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            n_epj_send_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        ReallocatableArray<S32> rank_src(n_proc);
        for(S32 i=0; i<n_proc; i++){
            n_image_send_per_proc[i] = 0;
            if(n_epj_src_per_proc[i] > 0) rank_src.push_back(i);
        }
        ReallocatableArray<F64> len_peri(DIMENSION);
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        for(S32 i=0; i<DIMENSION; i++){
            if(pa[i]==false) len_peri[i] = 0.0;
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            epj_sorted_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ptcl_send_tmp[ith].clearSize();
            rank_dst_tmp[ith].clearSize();
            tc_tmp[ith].clearSize();
            tp_tmp[ith].clearSize();
            adr_epj_src_per_image_tmp[ith].clearSize();
            n_epj_src_per_image_tmp[ith].clearSize();
            n_epj_send_per_image_tmp[ith].clearSize();
            //S32 n_ep_send_cum_old = 0;
            //bool first_loop = true;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 1)
            //#pragma omp for schedule(static)
#endif
            for(S32 ib=0; ib<rank_src.size(); ib++){
                const S32 rank_tmp = rank_src[ib];
                rank_dst_tmp[ith].push_back(rank_tmp);
                /*
                if(Comm::getRank() == 0){
                    std::cerr<<"rank_tmp= "<<rank_tmp
                             <<" n_epj_per_proc[rank_tmp]= "
                             <<n_epj_per_proc[rank_tmp]
                             <<std::endl;
                }
                */
                if( n_epj_src_per_proc[rank_tmp] <= 0) continue;
                const S32 adr_ptcl_head = n_disp_epj_per_proc[rank_tmp];
                const S32 adr_ptcl_end = n_disp_epj_per_proc[rank_tmp+1];
                S32vec id_image_old = -9999;
                S32 n_image = 0;
                for(S32 ip=adr_ptcl_head; ip<adr_ptcl_end; ip++){
                    const F64vec pos_target = epj_A[ip].getPos();
                    const S32vec id_image_new = CalcIDOfImageDomain(pos_root_domain, pos_target, pa);
                    if(id_image_old != id_image_new){
                        adr_epj_src_per_image_tmp[ith].push_back(ip);
                        n_epj_src_per_image_tmp[ith].push_back(0);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
                        shift_per_image_tmp[ith].push_back(F64vec(id_image_new.x*len_root_domain.x, id_image_new.y*len_root_domain.y));
#else
                        shift_per_image_tmp[ith].push_back(F64vec(id_image_new.x*len_root_domain.x, id_image_new.y*len_root_domain.y, id_image_new.z*len_root_domain.z));
#endif
                        id_image_old = id_image_new;
                        n_image++;
                    }
                    n_epj_src_per_image_tmp[ith].back()++;
                }
                /*
                if(Comm::getRank()==0 && rank_tmp==1){
                    std::cerr<<"adr_ptcl_head= "<<adr_ptcl_head
                             <<" adr_ptcl_end= "<<adr_ptcl_end
                             <<std::endl;
                }
                */
                /*
                if(Comm::getRank()==0){
                    std::cerr<<"rank_tmp= "<<rank_tmp
                             <<" n_image= "<<n_image
                             <<" n_disp_epj_per_image_tmp[ith].size()= "<<n_disp_epj_per_image_tmp[ith].size()
                             <<std::endl;
                    for(S32 ip=adr_ptcl_head; ip<adr_ptcl_end; ip++){
                        std::cerr<<"ip, pos= "<<ip<<", "<<epj_A[ip].getPos()<<std::endl;
                    }
                }
                */
                n_image_send_per_proc[rank_tmp] = n_image;
                const S32 adr_image_end  = n_epj_src_per_image_tmp[ith].size();
                const S32 adr_image_head = adr_image_end - n_image;
                const S32 n_epj_send_cum_prev = adr_ptcl_send_tmp[ith].size();
                for(S32 ii=adr_image_head; ii<adr_image_end; ii++){
                    const F64vec shift = shift_per_image_tmp[ith][ii];
                    /*
                    if(Comm::getRank()==0){
                        std::cerr<<"rank_tmp= "<<rank_tmp
                                 <<" ii= "<<ii
                                 <<" shift= "<<shift
                                 <<" n_epj_src_per_image_tmp[ith][ii]= "<<n_epj_src_per_image_tmp[ith][ii]
                                 <<std::endl;
                    }
                    */
                    const F64ort pos_domain = dinfo.getPosDomain(rank_tmp).shift(shift);
                    //const F64ort pos_domain = dinfo.getPosDomain(rank_tmp).shift(-shift);
                    /*
                    if(Comm::getRank() == 0 && rank_tmp==1){
                        std::cerr<<"rank_tmp= "<<rank_tmp
                                 <<" shift= "<<shift
                                 <<" pos_domain= "<<pos_domain
                                 <<" n_disp_epj_per_image_tmp[ith][ii]= "
                                 <<n_disp_epj_per_image_tmp[ith][ii]
                                 <<" n_disp_epj_per_image_tmp[ith][ii+1]= "
                                 <<n_disp_epj_per_image_tmp[ith][ii+1]
                                 <<std::endl;
                    }
                    */
                    ///////////////
                    // MAKE TREE A
                    S32 n_cnt = 0;
                    tp_tmp[ith].resizeNoInitialize(n_epj_src_per_image_tmp[ith][ii]);
                    for(S32 ip=adr_epj_src_per_image_tmp[ith][ii]; ip<adr_epj_src_per_image_tmp[ith][ii]+n_epj_src_per_image_tmp[ith][ii]; ip++, n_cnt++){
                        tp_tmp[ith][n_cnt].setFromEP(epj_A[ip], ip);
                        //if(Comm::getRank()==0) std::cerr<<"epj_A[ip].pos= "<<epj_A[ip].pos<<std::endl;
                    }
                    std::sort(tp_tmp[ith].getPointer(), tp_tmp[ith].getPointer()+n_cnt, LessOPKEY());
                    epj_sorted_tmp[ith].resizeNoInitialize(n_cnt);
                    for(S32 ip=0; ip<n_cnt; ip++){
                        const S32 adr = tp_tmp[ith][ip].adr_ptcl_;
                        epj_sorted_tmp[ith][ip] = epj_A[adr];
                        tp_tmp[ith][ip].adr_ptcl_ = ip;
                        /*
                        if(Comm::getRank()==0 && rank_tmp == 1){
                            std::cout<<"ip= "<<ip<<" epj_A[adr].pos= "<<epj_A[adr].pos<<std::endl;
                        }
                        */
                    }
                    /*
                    if(Comm::getRank() == 0){
                        for(S32 ip=0; ip<n_cnt; ip++){
                            std::cerr<<"epj_sorted_tmp[ith][ip].pos= "<<epj_sorted_tmp[ith][ip].pos<<std::endl;
                        }
                    }
                    */
                    const S32 n_leaf_limit_A = 1;
                    S32 lev_max_A = 0;
                    LinkCellST(tc_tmp[ith], adr_tc_level_partition_tmp[ith].getPointer(),
                               tp_tmp[ith].getPointer(), lev_max_A, n_cnt, n_leaf_limit_A);
                    CalcMomentST(adr_tc_level_partition_tmp[ith].getPointer(),
                                 tc_tmp[ith].getPointer(), 
                                 epj_sorted_tmp[ith].getPointer(), lev_max_A, n_leaf_limit_A);
                    /*
                    if(Comm::getRank() == 0 && rank_tmp == 1){
                        std::cerr<<"tc_tmp[ith][0].mom_.getVertexOut()= "<<tc_tmp[ith][0].mom_.getVertexOut()
                                 <<"tc_tmp[ith][0].mom_.getVertexIn()= "<<tc_tmp[ith][0].mom_.getVertexIn()
                                 <<" epj_sorted_tmp[ith][0].pos= "<<epj_sorted_tmp[ith][0].pos
                                 <<std::endl;
                    }
                    */
                    /*
                    if(Comm::getRank() == 0){
                        S32 err = 0;
                        tc_tmp[ith].getPointer()->checkTree(epj_sorted_tmp[ith].getPointer(),
                                                            tc_tmp[ith].getPointer(),
                                                            center_tree, full_len_tree*0.5,
                                                            n_leaf_limit_A, 1e-4,
                                                            err);
                        tc_tmp[ith].getPointer()->dumpTree(epj_sorted_tmp[ith].getPointer(),
                                                           tc_tmp[ith].getPointer(),
                                                           center_tree, full_len_tree*0.5,
                                                           n_leaf_limit_A);
                    }
                    */
                    const S32 n_epj_send_per_image_prev = adr_ptcl_send_tmp[ith].size();
                    MakeListDoubleWalkTop(tc_tmp[ith], tc_first_B,
                                          epj_A,    pos_domain,
                                          n_leaf_limit_A, n_leaf_limit_B,
                                          adr_ptcl_send_tmp[ith]);
                    n_epj_send_per_image_tmp[ith].push_back(adr_ptcl_send_tmp[ith].size() - n_epj_send_per_image_prev);

                    if(Comm::getRank() == 0 && rank_tmp == 1){
                        /*
                        std::cout<<"rank_tmp= "<<rank_tmp
                                 <<" shift= "<<shift
                                 <<" pos_domain= "<<pos_domain
                                 <<" tc_tmp[ith].size()= "<<tc_tmp[ith].size()
                                 <<" tc_first_B.size()= "<<tc_first_B.size()
                                 <<" n_epj_send_per_image_tmp[ith].back()= "<<n_epj_send_per_image_tmp[ith].back()
                                 <<" adr_ptcl_send_tmp[ith].size()= "<<adr_ptcl_send_tmp[ith].size()
                                 <<std::endl;
                        */
                        /*
                        for(S32 iii=n_epj_send_per_image_prev; iii<adr_ptcl_send_tmp[ith].size(); iii++){
                            std::cerr<<"adr_ptcl_send_tmp[ith][iii]= "<<adr_ptcl_send_tmp[ith][iii]
                                     <<" epj_B[adr_ptcl_send_tmp[ith][iii]].id= "<<epj_B[adr_ptcl_send_tmp[ith][iii]].id
                                     <<" pos= "<<epj_B[adr_ptcl_send_tmp[ith][iii]].pos
                                     <<std::endl;
                        }
                        */
                    }

                    /*
                    const S32 n_ep_per_image = adr_ptcl_send_tmp[ith].size();
                    n_ep_per_image_tmp[ith].push_back(n_ep_per_image);
                    for(S32 jp=0; jp<n_ep_per_image; jp++){
                        const S32 adr = adr_ptcl_send_tmp[ith][jp];
                        //const F64vec pos_j = epj[adr].getPos();
                        //if( pos_root_cell_.notOverlapped(pos_j-shift) ) continue;
                    }
                    */
                } // end of for ii=0 to n_image
                n_epj_send_per_proc[rank_tmp] = adr_ptcl_send_tmp[ith].size() - n_epj_send_cum_prev;
            } // end of OMP for
        } // end of OMP scope
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"n_image_send_per_proc[i]= "<<n_image_send_per_proc[i]
                         <<" n_epj_send_per_proc[i]= "<<n_epj_send_per_proc[i]
                         <<std::endl;
            }
        }
        */
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_thread; i++){
                std::cerr<<"rank_send_tmp[i].size()= "<<rank_send_tmp[i].size()<<std::endl;
                for(S32 j=0; j<rank_send_tmp[i].size(); j++){
                    std::cerr<<"rank_send_tmp[i][j]= "<<rank_send_tmp[i][j]<<std::endl;
                }
            }
        }
        */

        ReallocatableArray<S32> n_disp_image_send_per_proc;
        n_disp_image_send_per_proc.resizeNoInitialize(n_proc+1);
        n_disp_image_send_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_send_per_proc[i+1] = n_disp_image_send_per_proc[i] + n_image_send_per_proc[i];
        }
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"n_image_send_per_proc[i]= "<<n_image_send_per_proc[i]
                         <<" n_disp_image_send_per_proc[i]= "<<n_disp_image_send_per_proc[i]
                         <<std::endl;
            }
        }
        */
        
        /*
        S32 n_image_send_tot = 0;
        for(S32 i=0; i<n_thread; i++){
            n_image_send_tot += shift_per_image_tmp[i].size();
        }
        if(Comm::getRank()==0){
            std::cerr<<"n_image_send_tot= "<<n_image_send_tot<<std::endl;
        }
        */

        const S32 n_image_send_tot = n_disp_image_send_per_proc[n_proc];
        n_epj_send_per_image.resizeNoInitialize(n_image_send_tot);
        shift_per_image.resizeNoInitialize(n_image_send_tot);
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_image = 0;
            //S32 n_cnt_ep = 0;
            for(S32 j=0; j<rank_dst_tmp[i].size(); j++){
                const S32 rank = rank_dst_tmp[i][j];
                const S32 adr_image_head = n_disp_image_send_per_proc[rank];
                const S32 adr_image_end  = n_disp_image_send_per_proc[rank+1];
                /*
                if(Comm::getRank()==0){
                    std::cerr<<"rank= "<<rank
                             <<" adr_image_head= "<<adr_image_head
                             <<" adr_image_end= "<<adr_image_end
                             <<std::endl;
                }
                */
                for(S32 k=adr_image_head; k<adr_image_end; k++, n_cnt_image++){
                    n_epj_send_per_image[k] = n_epj_send_per_image_tmp[i][n_cnt_image];
                    shift_per_image[k] = shift_per_image_tmp[i][n_cnt_image];
                    /*
                    if(Comm::getRank()==0){
                        std::cerr<<"k= "<<k
                                 <<" n_epj_send_per_image[k]= "<<n_epj_send_per_image[k]
                                 <<std::endl;
                    }
                    */
                }
            }
        }
        ReallocatableArray<S32> n_disp_epj_send_per_image;
        n_disp_epj_send_per_image.resizeNoInitialize(n_image_send_tot+1);
        n_disp_epj_send_per_image[0] = 0;
        for(S32 i=0; i<n_image_send_tot; i++){
            n_disp_epj_send_per_image[i+1] = n_disp_epj_send_per_image[i] + n_epj_send_per_image[i];
        }
        S32 n_epj_send_tot = 0;
        for(S32 i=0; i<n_proc; i++) n_epj_send_tot += n_epj_send_per_proc[i];
        adr_ep_send.resizeNoInitialize(n_epj_send_tot);
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_image = 0;
            S32 n_cnt_ep = 0;
            for(S32 j=0; j<rank_dst_tmp[i].size(); j++){
                const S32 rank = rank_dst_tmp[i][j];
                const S32 adr_image_head = n_disp_image_send_per_proc[rank];
                const S32 adr_image_end  = n_disp_image_send_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++, n_cnt_image++){
                    n_epj_send_per_image[k] = n_epj_send_per_image_tmp[i][n_cnt_image];
                    shift_per_image[k] = shift_per_image_tmp[i][n_cnt_image];
                    const S32 adr_epj_head = n_disp_epj_send_per_image[k];
                    const S32 adr_epj_end  = n_disp_epj_send_per_image[k+1];
                    for(S32 l=adr_epj_head; l<adr_epj_end; l++, n_cnt_ep++){
                        adr_ep_send[l] = adr_ptcl_send_tmp[i][n_cnt_ep];
                        /*
                        if(Comm::getRank()==0){
                            std::cerr<<"l= "<<l
                                     <<" adr_ep_send[l]= "<<adr_ep_send[l]
                                     <<std::endl;
                        }
                        */
                    }
                }
            }
        }
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_proc; i++){
                std::cout<<"rank= "<<i<<std::endl;
                const S32 adr_image_head = n_disp_image_send_per_proc[i];
                const S32 adr_image_end  = n_disp_image_send_per_proc[i+1];
                for(S32 j=adr_image_head; j<adr_image_end; j++){
                    std::cout<<"image= "<<j
                             <<" shift_per_image[j]= "<<shift_per_image[j]
                             <<std::endl;
                    const S32 adr_epj_head = n_disp_epj_send_per_image[j];
                    const S32 adr_epj_end  = n_disp_epj_send_per_image[j+1];
                    for(S32 k=adr_epj_head; k<adr_epj_end; k++){
                        std::cout<<"k= "<<k
                                 <<" adr_ep_send[k]= "<<adr_ep_send[k]
                                 <<std::endl;
                    }
                }
            }
        }
        */
        //Comm::barrier();
        //exit(1);        
    }

    //////////////////
    // add moment as sp
    template<class Ttreecell, class Tspj>
    inline void AddMomentAsSpImpl(TagForceLong,
                                  const ReallocatableArray<Ttreecell> & _tc,
                                  const S32 offset,
                                  ReallocatableArray<Tspj> & _spj){
        _spj.resizeNoInitialize(offset+_tc.size());
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif //PARTICLE_SIMULATOR_THREAD_PARALLEL
        for(S32 i=0; i<_tc.size(); i++){
            _spj[offset+i].copyFromMoment(_tc[i].mom_);
        }
    }    
    template<class Ttreecell, class Tspj>
    inline void AddMomentAsSpImpl(TagForceShort,
                                  const ReallocatableArray<Ttreecell> & _tc,
                                  const S32 offset,
                                  ReallocatableArray<Tspj> & _spj){
        // do nothing
    }


    
    // for long force
    template<class Ttp, class Tepj, class Tspj>
    inline void SetLocalEssentialTreeToGlobalTreeImpl
    (const ReallocatableArray<Tepj> & epj_org,
     const ReallocatableArray<Tspj> & spj_org,
     const S32 n_loc,
     ReallocatableArray<Ttp> & tp_glb,
     const bool flag_reuse = false){
        const S32 n_loc_ep = epj_org.size();
        const S32 n_loc_ep_sp = n_loc_ep + spj_org.size(); 
        tp_glb.resizeNoInitialize( n_loc_ep_sp );
        if(!flag_reuse){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
            {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
                for(S32 i=n_loc; i<n_loc_ep; i++){
                    tp_glb[i].setFromEP(epj_org[i], i);
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
                for(S32 i=n_loc_ep; i<n_loc_ep_sp; i++){
                    const S32 i_src = i-n_loc_ep;
                    tp_glb[i].setFromSP(spj_org[i_src], i_src);
                }
            }
        }
    }

    template<class Ttp, class Tepj>
    inline void SetLocalEssentialTreeToGlobalTreeImpl
    (
     const ReallocatableArray<Tepj> & epj_org,
     const S32 n_loc,
     ReallocatableArray<Ttp> & tp_glb,
     const bool flag_reuse = false){
        const S32 n_loc_ep = epj_org.size();
        tp_glb.resizeNoInitialize( n_loc_ep );
        if(!flag_reuse){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
            {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
                for(S32 i=n_loc; i<n_loc_ep; i++){
                    tp_glb[i].setFromEP(epj_org[i], i);
                }
            }
        }
    }

    template<class Ttc>
    inline void SetOuterBoxGlobalTreeForLongCutoffRecursive(Ttc tc[],
                                                            const S32 n_leaf_limit,
                                                            const F64 r_cut,
                                                            const S32 adr_tc,
                                                            const F64 tc_hlen,
                                                            const F64vec tc_cen){
        F64 child_hlen = tc_hlen*0.5;
        for(S32 i=0; i<N_CHILDREN; i++){
            F64vec child_cen = tc_cen + SHIFT_CENTER[i]*tc_hlen;
            tc[adr_tc+i].mom_.vertex_out_.high_ = child_cen + (child_hlen+r_cut);
            tc[adr_tc+i].mom_.vertex_out_.low_  = child_cen - (child_hlen+r_cut);
            S32 child_adr_tc = tc[adr_tc+i].adr_tc_;
            if(tc[adr_tc+i].n_ptcl_ <= 0) continue;
            else if(tc[adr_tc+i].isLeaf(n_leaf_limit)) continue;
            else{
                SetOuterBoxGlobalTreeForLongCutoffRecursive(tc, n_leaf_limit, r_cut, child_adr_tc,
                                                            child_hlen, child_cen);
            }
        }
    }
    
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongCutoff,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        F64 r_cut = epj[0].getRSearch();
        tc[0].mom_.vertex_out_.high_ = tc_cen + (tc_hlen+r_cut);
        tc[0].mom_.vertex_out_.low_  = tc_cen - (tc_hlen+r_cut);
        for(S32 i=1; i<N_CHILDREN; i++) tc[i].mom_.vertex_out_.init();
        if( tc[0].n_ptcl_ < 0 || tc[0].isLeaf(n_leaf_limit) ) return;
        S32 adr_tc = N_CHILDREN;
        SetOuterBoxGlobalTreeForLongCutoffRecursive(tc, n_leaf_limit, r_cut, adr_tc,
                                                    tc_hlen, tc_cen);
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLong,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongScatter,
                                               Ttc tc[],
                                               Tepj epj[],
                                               const S32 n_leaf_limit,
                                               const F64 tc_hlen,
                                               const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongSymmetry,
                                               Ttc tc[],
                                               Tepj epj[],
                                               const S32 n_leaf_limit,
                                               const F64 tc_hlen,
                                               const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortGather,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortScatter,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortSymmetry,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        // do nothing
    }

    // it may works only if tp has one component.
    template<class Tsys, class Ttp, class Tepi, class Tepj>
    inline void CopyFpToEpSortedLocalTree(const Tsys & sys,
                                   const ReallocatableArray<Ttp> & tp,
                                   ReallocatableArray<Tepi> &epi_sorted,
                                   ReallocatableArray<Tepj> &epj_sorted){
        const S32 n_loc = sys.getNumberOfParticleLocal();
        for(S32 i=0; i<n_loc; i++){
            const S32 adr = tp[i].adr_ptcl_;
            epi_sorted[i].copyFromFP(sys[adr]);
            epj_sorted[i].copyFromFP(sys[adr]);
        }
    }
}
