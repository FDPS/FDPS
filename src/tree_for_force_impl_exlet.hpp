#pragma once

namespace ParticleSimulator{

    // These classes are used only for exchanging LET with P2P and Allgather
    template<typename Tspj, typename Ttc, typename Tmode>
    struct TopMoment{
        Tspj spj_;
        F64ort vertex_in_;
        void set(const Tspj & spj,
                 const F64ort vertex_in,
                 const Ttc & tc){
            spj_ = spj;
            vertex_in_ = vertex_in;
        }
        template<typename Ttm>
        bool isOverlapped(const Ttm & tm){
            PARTICLE_SIMULATOR_PRINT_ERROR("FDPS error. This search mode is not correct.");
            Abort();
            return false;
        }
        F64ort getVertexOut() const {
            PARTICLE_SIMULATOR_PRINT_ERROR("FDPS error. This search mode is not correct.");
            Abort();
            return F64ort(-1234.5, -9876.5);
        }
    };

    template<typename Tspj, typename Ttc>
    struct TopMoment<Tspj, Ttc, SEARCH_MODE_LONG_SYMMETRY>{
        Tspj spj_;
        F64ort vertex_in_;
        F64ort vertex_out_;
        void set(const Tspj & spj,
                 const F64ort vertex_in,
                 const Ttc & tc){
            spj_ = spj;
            vertex_in_ = vertex_in;
            vertex_out_ = tc.geo_.getVertexOut();
        }
        bool isOverlapped(const TopMoment<Tspj, Ttc, SEARCH_MODE_LONG_SYMMETRY> & tm){
            return vertex_in_.overlapped(tm.vertex_out_) || vertex_out_.overlapped(tm.vertex_in_);
        }
        F64ort getVertexOut() const {
            return vertex_out_;
        }
    };

    template<typename Tspj, typename Ttc>
    struct TopMoment<Tspj, Ttc, SEARCH_MODE_LONG_SCATTER>{
        Tspj spj_;
        F64ort vertex_in_;
        F64ort vertex_out_; // To check if p2p is ok or not
        void set(const Tspj & spj,
                 const F64ort vertex_in,
                 const Ttc & tc){
            spj_ = spj;
            vertex_in_ = vertex_in;
            vertex_out_ = tc.geo_.getVertexOut();
        }
        bool isOverlapped(const TopMoment<Tspj, Ttc, SEARCH_MODE_LONG_SCATTER> & tm){
            return vertex_out_.overlapped(tm.vertex_in_);
        }
        F64ort getVertexOut() const {
            return vertex_out_;
        }
    };

    template<typename Ttc, typename Tmode>
    struct TopGeometry{
        F64ort vertex_in_;
        void set(const F64ort vertex_in,
                 const Ttc & tc){
            vertex_in_ = vertex_in;
        }
    };
    template<typename Ttc>
    struct TopGeometry<Ttc, SEARCH_MODE_LONG_SYMMETRY>{
        F64ort vertex_in_;
        F64ort vertex_out_;
        void set(const F64ort vertex_in,
                 const Ttc & tc){
            vertex_in_  = vertex_in;
            vertex_out_ = tc.geo_.getVertexOut();
        }
        F64ort getVertexOut() const {
            return vertex_out_;
        }
    };
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTree(const DomainInfo & dinfo,
                               const bool flag_reuse){
        if(typeid(TSM) == typeid(SEARCH_MODE_LONG)
           && dinfo.getBoundaryCondition() != BOUNDARY_CONDITION_OPEN){
            PARTICLE_SIMULATOR_PRINT_ERROR("The forces w/o cutoff can be evaluated only under the open boundary condition");
            Abort(-1);
        }
        if(!flag_reuse){ comm_table_.clearSize(); }
        exchangeLocalEssentialTreeImpl(typename TSM::search_type(), dinfo, flag_reuse);
        time_profile_.exchange_LET_tot = time_profile_.make_LET_1st
            + time_profile_.exchange_LET_1st
            + time_profile_.make_LET_2nd
            + time_profile_.exchange_LET_2nd;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeImpl(TagSearchLong,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        if(exchange_let_mode_ == EXCHANGE_LET_A2A){
            exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
        } else if(exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT){
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
        } else if(exchange_let_mode_ == EXCHANGE_LET_P2P_FAST){
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
            comm_table_.setPosDomainAllgather(dinfo, pos_root_cell_);
        }
    }
  
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeImpl(TagSearchLongCutoff,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        if(exchange_let_mode_ == EXCHANGE_LET_A2A){
            exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
        } else if(exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT){
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
        } else if(exchange_let_mode_ == EXCHANGE_LET_P2P_FAST){
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
            comm_table_.setPosDomainAllgather(dinfo, pos_root_cell_);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeImpl(TagSearchLongScatter,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        if(exchange_let_mode_ == EXCHANGE_LET_A2A){
            exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
        } else if(exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT){
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
        } else if(exchange_let_mode_ == EXCHANGE_LET_P2P_FAST){
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
            comm_table_.setPosDomainAllgather(dinfo, pos_root_cell_);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeImpl(TagSearchLongSymmetry,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        if(exchange_let_mode_ == EXCHANGE_LET_A2A){
            exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
        } else if(exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT){
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
        } else if(exchange_let_mode_ == EXCHANGE_LET_P2P_FAST){
            exchangeLocalEssentialTreeLongP2P(dinfo, flag_reuse);
            comm_table_.setPosDomainAllgather(dinfo, pos_root_cell_);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeImpl(TagSearchShortScatter,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        if(exchange_let_mode_ == EXCHANGE_LET_A2A){
            exchangeLocalEssentialTreeShortScatter(dinfo, flag_reuse);
        }else{
            exchangeLocalEssentialTreeShortScatterP2P(dinfo, flag_reuse);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        if(exchange_let_mode_ == EXCHANGE_LET_A2A){
            exchangeLocalEssentialTreeShortSymmetry(dinfo, flag_reuse);
        }else{
            exchangeLocalEssentialTreeShortSymmetryP2P(dinfo, flag_reuse);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeImpl(TagSearchShortGather,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        if(exchange_let_mode_ == EXCHANGE_LET_A2A){
            exchangeLocalEssentialTreeShortGather(dinfo, flag_reuse);
        } else{
            exchangeLocalEssentialTreeShortGatherP2P(dinfo, flag_reuse);
        }
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeLongP2P(const DomainInfo & dinfo,
                                      const bool flag_reuse){
        F64 time_offset = GetWtime();
        //const auto n_proc  = Comm::getNumberOfProc();
        const auto n_proc  = comm_info_.getNumberOfProc();
        F64ort pos_root_cell = getPosRootCell();
        ReallocatableArray<Tspj> top_spj;
        top_spj.setAllocMode(MemoryAllocMode::Pool);
        ReallocatableArray<TopMoment<Tspj, TreeCellLoc, TSM>> top_moment;
        top_moment.setAllocMode(MemoryAllocMode::Pool);
        Tspj my_top_spj;
        my_top_spj.copyFromMoment(tc_loc_[0].mom_);
        if(flag_reuse){
            top_spj.resizeNoInitialize(n_proc);
            const F64 time_offset_exchange_top_moment = GetWtime();
            comm_info_.allGather(&my_top_spj, 1, top_spj.getPointer());
            time_profile_.make_LET_1st__exchange_top_moment += GetWtime() - time_offset_exchange_top_moment;
        }
        if(!flag_reuse){
            top_moment.resizeNoInitialize(n_proc);
            TopMoment<Tspj, TreeCellLoc, TSM> my_top_moment;
            my_top_moment.set(my_top_spj, inner_boundary_of_local_tree_, tc_loc_[0]);
            const F64 time_offset_exchange_top_moment = GetWtime();
            //Comm::allGather(&my_top_moment, 1, top_moment.getPointer());
            comm_info_.allGather(&my_top_moment, 1, top_moment.getPointer());
            time_profile_.make_LET_1st__exchange_top_moment += GetWtime() - time_offset_exchange_top_moment;
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticleP2P<TSM, TreeCellLoc, TreeParticle, Tepj, Tspj, TopMoment<Tspj, TreeCellLoc, TSM>, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_,
                 epj_sorted_,
                 comm_table_.n_ep_send_,   comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.n_sp_send_,   comm_table_.adr_sp_send_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.n_sp_per_image_,
                 comm_table_.rank_send_,
                 comm_table_.rank_recv_,
                 top_moment,
                 pos_root_cell,
                 theta_,
                 comm_table_.rank_recv_allgather_,
                 morton_key_,
                 exchange_let_mode_);
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberLong(comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                               comm_table_.n_sp_send_, comm_table_.n_sp_recv_,
                               comm_table_.rank_send_, comm_table_.rank_recv_,
                               comm_info_);
            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
            comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
            S32 n_sp_recv_tot_tmp = 0;
PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp))
            for(S32 i=0; i<comm_table_.rank_recv_.size(); i++){
                S32 rank = comm_table_.rank_recv_[i];
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[rank];
                n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[rank];                
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
            comm_table_.n_sp_recv_tot_ = n_sp_recv_tot_tmp;
        }
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        F64 wtime_ex_let_comm = 0.0;
        time_offset = GetWtime();
        ExchangeLetP2P<TSM, Tepj, Tspj>
            (epj_sorted_, comm_table_.n_ep_send_,
             comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
             comm_table_.adr_ep_send_,
             epj_org_, n_loc_tot_,
             tc_loc_, comm_table_.n_sp_send_,
             comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
             comm_table_.adr_sp_send_,
             spj_org_,
             comm_table_.shift_per_image_,
             comm_table_.n_image_per_proc_,
             comm_table_.rank_send_,
             comm_table_.rank_recv_,
             comm_info_);
        time_profile_.exchange_LET_1st__icomm_ptcl += wtime_ex_let_comm;
        const auto spj_offset = spj_org_.size();
        const auto n_write = comm_table_.rank_recv_allgather_.size();
        spj_org_.increaseSize(n_write);
PS_OMP_PARALLEL
        {
            const S32 ith = Comm::getThreadNum();
            const S32 n_thread = Comm::getNumberOfThread();
            S32 head, end;
            CalcAdrToSplitData(head, end, ith, n_thread, n_write);
            if(flag_reuse){            
                for(S32 i=head; i<end; i++){
                    spj_org_[i+spj_offset] = top_spj[comm_table_.rank_recv_allgather_[i]];
                }
            }
            else{
                for(S32 i=head; i<end; i++){
                    spj_org_[i+spj_offset] = top_moment[comm_table_.rank_recv_allgather_[i]].spj_;
                }
            }
        }
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
        
#ifdef PS_DEBUG_EXCHANGE_LET
        //#if 1
        F64 mass = 0.0;
        for(int i=0; i<epj_org_.size(); i++){
            mass += epj_org_[i].getCharge();
        }
        for(int i=0; i<spj_org_.size(); i++){
            mass += spj_org_[i].getCharge();
        }
        //Comm::barrier();
        comm_info_.barrier();
        if(mass != 1.0){
            //std::cerr<<"my_rank= "<<Comm::getRank()
            std::cerr<<"my_rank= "<<comm_info_.getRank()
                     <<" mass= "<<mass
                     <<std::endl;
        }
        //Comm::barrier();
        comm_info_.barrier();
        assert(mass == 1.0);
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeLong(const DomainInfo & dinfo,
                                   const bool flag_reuse){
        F64 time_offset = GetWtime();
        if(!flag_reuse){
            const auto n_proc  = comm_info_.getNumberOfProc();
            ReallocatableArray<TopGeometry<TreeCellLoc, TSM>> top_geo(n_proc, n_proc, MemoryAllocMode::Pool);
            TopGeometry<TreeCellLoc, TSM> my_top_geo;
            my_top_geo.set(inner_boundary_of_local_tree_, tc_loc_[0]);
            const F64 time_offset_exchange_top_moment = GetWtime();
            comm_info_.allGather(&my_top_geo, 1, top_geo.getPointer());
            time_profile_.make_LET_1st__exchange_top_moment += GetWtime() - time_offset_exchange_top_moment;
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticle<TSM, TreeCellLoc, TreeParticle, Tepj, Tspj, TopGeometry<TreeCellLoc, TSM>, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_,
                 epj_sorted_,
                 comm_table_.n_ep_send_,   comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.n_sp_send_,   comm_table_.adr_sp_send_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.n_sp_per_image_,
                 top_geo,
                 theta_);
            
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberLong(comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                               comm_table_.n_sp_send_, comm_table_.n_sp_recv_,
                               comm_info_);

            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
            comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
            S32 n_sp_recv_tot_tmp = 0;
PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp))
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
                n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
            comm_table_.n_sp_recv_tot_ = n_sp_recv_tot_tmp;
        }
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();
        
        ExchangeLet<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_,
                                     comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                     comm_table_.adr_ep_send_,
                                     epj_org_, n_loc_tot_,
                                     tc_loc_, comm_table_.n_sp_send_,
                                     comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
                                     comm_table_.adr_sp_send_,
                                     spj_org_,
                                     comm_table_.shift_per_image_,
                                     comm_table_.n_image_per_proc_,
                                     comm_info_);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }


#if 0

    /*
    void  ExchangeNumberLong2(const ReallocatableArray<S32> & n_ep_send,
			      ReallocatableArray<S32> * n_ep_recv,
			      const ReallocatableArray<S32> * n_sp_send,
			      ReallocatableArray<S32> & n_sp_recv,
			      const ReallocatableArray<S32> & rank_send,
			      const ReallocatableArray<S32> & rank_recv,
			      const CommInfo & comm_info,
			      const int dim){
	
    }
    */
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeLong2(const DomainInfo & dinfo,
				    const bool flag_reuse){
        F64 time_offset = GetWtime();
        if(!flag_reuse){
            const auto n_proc  = comm_info_.getNumberOfProc();
            ReallocatableArray<TopGeometry<TreeCellLoc, TSM>> top_geo(n_proc, n_proc, MemoryAllocMode::Pool);
            TopGeometry<TreeCellLoc, TSM> my_top_geo;
            my_top_geo.set(inner_boundary_of_local_tree_, tc_loc_[0]);
            const F64 time_offset_exchange_top_moment = GetWtime();
            comm_info_.allGather(&my_top_geo, 1, top_geo.getPointer());
            time_profile_.make_LET_1st__exchange_top_moment += GetWtime() - time_offset_exchange_top_moment;
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticle<TSM, TreeCellLoc, TreeParticle, Tepj, Tspj, TopGeometry<TreeCellLoc, TSM>, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_,
                 epj_sorted_,
                 comm_table_.n_ep_send_,   comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.n_sp_send_,   comm_table_.adr_sp_send_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.n_sp_per_image_,
                 top_geo,
                 theta_);
            
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberLong2(comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                               comm_table_.n_sp_send_, comm_table_.n_sp_recv_,
                               comm_info_);

            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
            comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
            S32 n_sp_recv_tot_tmp = 0;
PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp))
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
                n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
            comm_table_.n_sp_recv_tot_ = n_sp_recv_tot_tmp;
        }
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();
        
        ExchangeLet2<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_,
                                     comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                     comm_table_.adr_ep_send_,
                                     epj_org_, n_loc_tot_,
                                     tc_loc_, comm_table_.n_sp_send_,
                                     comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
                                     comm_table_.adr_sp_send_,
                                     spj_org_,
                                     comm_table_.shift_per_image_,
                                     comm_table_.n_image_per_proc_,
                                     comm_info_);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }    
#endif
  
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeLongCutoffP2P(const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 time_offset = GetWtime();
        time_profile_.make_LET_1st__exchange_top_moment = 0.0;
        if(!flag_reuse){
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticleLongCutoffP2P<TSM, TreeCellLoc, TreeParticle, Tepj, Tspj, TopMoment<Tspj, TreeCellLoc, TSM>, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_,
                 epj_sorted_,
                 comm_table_.n_ep_send_,   comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.n_sp_send_,   comm_table_.adr_sp_send_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.n_sp_per_image_,
                 comm_table_.rank_send_,
                 comm_table_.rank_recv_,
                 theta_,
                 exchange_let_mode_);
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberLong(comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                               comm_table_.n_sp_send_, comm_table_.n_sp_recv_,
                               comm_table_.rank_send_, comm_table_.rank_recv_,
                               comm_info_);
            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
            comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
            S32 n_sp_recv_tot_tmp = 0;
PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp))
            for(S32 i=0; i<comm_table_.rank_recv_.size(); i++){
                S32 rank = comm_table_.rank_recv_[i];
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[rank];
                n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[rank];                
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
            comm_table_.n_sp_recv_tot_ = n_sp_recv_tot_tmp;
        }
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        F64 wtime_ex_let_comm = 0.0;
        time_offset = GetWtime();
        ExchangeLetP2P<TSM, Tepj, Tspj>
            (epj_sorted_, comm_table_.n_ep_send_,
             comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
             comm_table_.adr_ep_send_,
             epj_org_, n_loc_tot_,
             tc_loc_, comm_table_.n_sp_send_,
             comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
             comm_table_.adr_sp_send_,
             spj_org_,
             comm_table_.shift_per_image_,
             comm_table_.n_image_per_proc_,
             comm_table_.rank_send_,
             comm_table_.rank_recv_,
             wtime_ex_let_comm);
        time_profile_.exchange_LET_1st__icomm_ptcl += wtime_ex_let_comm;
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }

    
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeShortScatterP2P(const DomainInfo & dinfo,
                                              const bool flag_reuse){
        auto time_offset = GetWtime();
        if(!flag_reuse){
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticleP2P<TreeCellLoc, TreeParticle, Tepj, Tspj, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_, epj_sorted_,
                 comm_table_.n_ep_send_,  comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.rank_send_,
                 comm_table_.rank_recv_);
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberShort(comm_table_.n_ep_send_,
                                comm_table_.n_ep_recv_,
                                comm_table_.rank_send_,
                                comm_table_.rank_recv_,
                                comm_info_);
            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            //const auto n_proc = Comm::getNumberOfProc();
            const auto n_proc = comm_info_.getNumberOfProc();
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_ep_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp))
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        }
        time_profile_.make_LET_1st += GetWtime() - time_offset;

        time_offset = GetWtime();
        ExchangeLetAllgather<TSM, Tepj>(epj_sorted_,
                                        comm_table_.n_ep_send_,
                                        comm_table_.n_ep_recv_,
                                        comm_table_.n_ep_per_image_,
                                        comm_table_.adr_ep_send_,
                                        epj_org_, n_loc_tot_, 
                                        comm_table_.shift_per_image_,
                                        comm_table_.n_image_per_proc_,
                                        comm_table_.rank_send_,
                                        comm_table_.rank_recv_,
                                        comm_info_);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeShortScatter(const DomainInfo & dinfo,
                                           const bool flag_reuse){
        auto time_offset = GetWtime();
        if(!flag_reuse){
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticle<TreeCellLoc, TreeParticle, Tepj, Tspj, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_, epj_sorted_,
                 comm_table_.n_ep_send_,  comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_);
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberShort(comm_table_.n_ep_send_,
                                comm_table_.n_ep_recv_,
                                comm_info_);
            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            //const auto n_proc = Comm::getNumberOfProc();
            const auto n_proc = comm_info_.getNumberOfProc();
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_ep_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:n_ep_recv_tot_tmp)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        }
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        
        time_offset = GetWtime();
        ExchangeLet<TSM, Tepj>(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                               comm_table_.n_ep_per_image_,
                               comm_table_.adr_ep_send_,
                               epj_org_, n_loc_tot_, 
                               comm_table_.shift_per_image_,
                               comm_table_.n_image_per_proc_,
                               comm_info_);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }

    
    ////////////////
    // SHORT SYMMETRY
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeShortSymmetryP2P
    (const DomainInfo & dinfo,
     const bool flag_reuse){
        F64 time_offset = GetWtime();
        const auto n_proc = comm_info_.getNumberOfProc();
        ReallocatableArray<S32> n_ep_send_per_proc_1st(MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_ep_send_per_proc_2nd(MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_ep_recv_per_proc_1st(MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_ep_recv_per_proc_2nd(MemoryAllocMode::Pool);
        ReallocatableArray<S32> adr_ep_send_1st(MemoryAllocMode::Pool);
        ReallocatableArray<S32> adr_ep_send_2nd(MemoryAllocMode::Pool);
        ReallocatableArray<F64vec> shift_per_image_1st(MemoryAllocMode::Pool);
        ReallocatableArray<F64vec> shift_per_image_2nd(MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_image_per_proc_1st(MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_image_per_proc_2nd(MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_ep_send_per_image_1st(MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_ep_send_per_image_2nd(MemoryAllocMode::Pool);
        ReallocatableArray<Tepj> ep_recv_1st(MemoryAllocMode::Pool);
        
        if(!flag_reuse){
            comm_table_.rank_send_.resizeNoInitialize(n_proc);
            comm_table_.rank_recv_.resizeNoInitialize(n_proc);
            ////////////
            // 1st STEP (send j particles)
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticleP2P<TreeCellLoc, TreeParticle, Tepj, Tspj, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_, epj_sorted_,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st,
                 comm_table_.rank_send_,
                 comm_table_.rank_recv_);
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberShort(n_ep_send_per_proc_1st,
                                n_ep_recv_per_proc_1st,
                                comm_table_.rank_send_,
                                comm_table_.rank_recv_,
                                comm_info_);
            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            ExchangeLetAllgather<TSM, Tepj>
                (epj_sorted_,
                 n_ep_send_per_proc_1st,
                 n_ep_recv_per_proc_1st,
                 n_ep_send_per_image_1st,
                 adr_ep_send_1st,
                 ep_recv_1st,
                 0, 
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 comm_table_.rank_send_,
                 comm_table_.rank_recv_,
                 comm_info_);

            ////////////
            // 2nd STEP
            FindExchangeParticleDoubleWalk<TreeCellLoc, TreeParticle, Tepj>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epj_sorted_,
                 center_, length_, morton_key_);
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumberShort(n_ep_send_per_proc_2nd,
                                n_ep_recv_per_proc_2nd,
                                comm_table_.rank_send_,
                                comm_table_.rank_recv_,
                                comm_info_);
            
            /////////////////////
            // 4th STEP (make communication table)
            MakeCommTableFor2StepCommuniction(comm_table_,
                                              n_ep_send_per_proc_1st,
                                              n_image_per_proc_1st,
                                              shift_per_image_1st,
                                              n_ep_send_per_image_1st,
                                              adr_ep_send_1st,
                                              n_ep_recv_per_proc_1st,
                                              n_ep_send_per_proc_2nd,
                                              n_image_per_proc_2nd,
                                              shift_per_image_2nd,
                                              n_ep_send_per_image_2nd,
                                              adr_ep_send_2nd,
                                              n_ep_recv_per_proc_2nd,
                                              n_proc);
        } // end of reuse
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        
        time_offset = GetWtime();
        ExchangeLetAllgather<TSM, Tepj>(epj_sorted_,
                               comm_table_.n_ep_send_,
                               comm_table_.n_ep_recv_,
                               comm_table_.n_ep_per_image_,
                               comm_table_.adr_ep_send_,
                               epj_org_, n_loc_tot_, 
                               comm_table_.shift_per_image_,
                               comm_table_.n_image_per_proc_,
                               comm_table_.rank_send_,
                                        comm_table_.rank_recv_,
                                        comm_info_);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeShortSymmetry(const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 time_offset = GetWtime();
        ReallocatableArray<S32> n_ep_send_per_proc_1st;
        ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        ReallocatableArray<S32> adr_ep_send_1st;
        ReallocatableArray<S32> adr_ep_send_2nd;
        ReallocatableArray<F64vec> shift_per_image_1st;
        ReallocatableArray<F64vec> shift_per_image_2nd;
        ReallocatableArray<S32> n_image_per_proc_1st;
        ReallocatableArray<S32> n_image_per_proc_2nd;
        ReallocatableArray<S32> n_ep_send_per_image_1st;
        ReallocatableArray<S32> n_ep_send_per_image_2nd;
        ReallocatableArray<Tepj> ep_recv_1st;
        
        if(!flag_reuse){
            ////////////
            // 1st STEP (send j particles)
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticle<TreeCellLoc, TreeParticle, Tepj, Tspj, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_, epj_sorted_,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberShort(n_ep_send_per_proc_1st,
                                n_ep_recv_per_proc_1st,
                                comm_info_);
            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            ExchangeLet<TSM, Tepj>
                (epj_sorted_,
                 n_ep_send_per_proc_1st,
                 n_ep_recv_per_proc_1st,
                 n_ep_send_per_image_1st,
                 adr_ep_send_1st,
                 ep_recv_1st,
                 0,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 comm_info_);
            
            ////////////
            // 2nd STEP (send j particles)
            FindExchangeParticleDoubleWalk<TreeCellLoc, TreeParticle, Tepj>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epj_sorted_,
                 center_, length_, morton_key_);
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumberShort(n_ep_send_per_proc_2nd,
                                n_ep_recv_per_proc_2nd,
                                comm_info_);

            /////////////////////
            // 4th STEP (make communication table)
            const auto n_proc = comm_info_.getNumberOfProc();
            MakeCommTableFor2StepCommuniction(comm_table_,
                                              n_ep_send_per_proc_1st,
                                              n_image_per_proc_1st,
                                              shift_per_image_1st,
                                              n_ep_send_per_image_1st,
                                              adr_ep_send_1st,
                                              n_ep_recv_per_proc_1st,
                                              n_ep_send_per_proc_2nd,
                                              n_image_per_proc_2nd,
                                              shift_per_image_2nd,
                                              n_ep_send_per_image_2nd,
                                              adr_ep_send_2nd,
                                              n_ep_recv_per_proc_2nd,
                                              n_proc);
        } // end of reuse
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        
        time_offset = GetWtime();
        ExchangeLet<TSM, Tepj>
            (epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
             comm_table_.n_ep_per_image_,
             comm_table_.adr_ep_send_,
             epj_org_, n_loc_tot_,
             comm_table_.shift_per_image_, comm_table_.n_image_per_proc_, comm_info_);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }

    ////////////////
    // SHORT GATHER MODE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeShortGather(const DomainInfo & dinfo,
                                          const bool flag_reuse){
        F64 time_offset = GetWtime();
        ReallocatableArray<S32> n_ep_send_per_proc_1st;
        ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        ReallocatableArray<S32> adr_ep_send_1st;
        ReallocatableArray<S32> adr_ep_send_2nd;
        ReallocatableArray<F64vec> shift_per_image_1st;
        ReallocatableArray<F64vec> shift_per_image_2nd;
        ReallocatableArray<S32> n_image_per_proc_1st;
        ReallocatableArray<S32> n_image_per_proc_2nd;
        ReallocatableArray<S32> n_ep_send_per_image_1st;
        ReallocatableArray<S32> n_ep_send_per_image_2nd;
        ReallocatableArray<EssentialParticleBase> ep_recv_1st;
        ReallocatableArray<EssentialParticleBase> epi_base_sorted;
        
        const S32 n_epi_sorted = epi_sorted_.size();
        epi_base_sorted.resizeNoInitialize(n_epi_sorted);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_epi_sorted; i++){
            epi_base_sorted[i].pos      = epi_sorted_[i].getPos();
            epi_base_sorted[i].r_search = epi_sorted_[i].getRSearch();
        }
        if(!flag_reuse){
            ////////////
            // 1st STEP (send epi_base)
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticle<TreeCellLoc, TreeParticle, EssentialParticleBase, Tspj, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_, epi_base_sorted,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberShort(n_ep_send_per_proc_1st,
                                n_ep_recv_per_proc_1st,
                                comm_info_);
            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            ExchangeLet<TSM, EssentialParticleBase>
                (epi_base_sorted,
                 n_ep_send_per_proc_1st,
                 n_ep_recv_per_proc_1st,
                 n_ep_send_per_image_1st,
                 adr_ep_send_1st,
                 ep_recv_1st,
                 0,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 comm_info_);
            ////////////
            // 2nd STEP (find j particle)
            FindExchangeParticleDoubleWalk<TreeCellLoc, TreeParticle, EssentialParticleBase>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epi_base_sorted,
                 center_, length_, morton_key_);
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumberShort(n_ep_send_per_proc_2nd,
                                n_ep_recv_per_proc_2nd,
                                comm_info_);
            /////////////////////
            // 4th STEP (make communication table)
            const auto n_proc = comm_info_.getNumberOfProc();
            MakeCommTableFor2StepCommuniction(comm_table_,
                                              n_ep_send_per_proc_1st,
                                              n_image_per_proc_1st,
                                              shift_per_image_1st,
                                              n_ep_send_per_image_1st,
                                              adr_ep_send_1st,
                                              n_ep_recv_per_proc_1st,
                                              n_ep_send_per_proc_2nd,
                                              n_image_per_proc_2nd,
                                              shift_per_image_2nd,
                                              n_ep_send_per_image_2nd,
                                              adr_ep_send_2nd,
                                              n_ep_recv_per_proc_2nd,
                                              n_proc);
        } // end of reuse flag
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();
        ExchangeLet<TSM, Tepj>
            (epj_sorted_,
             comm_table_.n_ep_send_,
             comm_table_.n_ep_recv_,
             comm_table_.n_ep_per_image_,
             comm_table_.adr_ep_send_,
             epj_org_,
             n_loc_tot_,
             comm_table_.shift_per_image_,
             comm_table_.n_image_per_proc_,
             comm_info_);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    exchangeLocalEssentialTreeShortGatherP2P(const DomainInfo & dinfo,
                                             const bool flag_reuse){
        F64 time_offset = GetWtime();
        ReallocatableArray<S32> n_ep_send_per_proc_1st;
        ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        ReallocatableArray<S32> adr_ep_send_1st;
        ReallocatableArray<S32> adr_ep_send_2nd;
        ReallocatableArray<F64vec> shift_per_image_1st;
        ReallocatableArray<F64vec> shift_per_image_2nd;
        ReallocatableArray<S32> n_image_per_proc_1st;
        ReallocatableArray<S32> n_image_per_proc_2nd;
        ReallocatableArray<S32> n_ep_send_per_image_1st;
        ReallocatableArray<S32> n_ep_send_per_image_2nd;
        ReallocatableArray<EssentialParticleBase> ep_recv_1st;
        ReallocatableArray<EssentialParticleBase> epi_base_sorted;
        
        const S32 n_epi_sorted = epi_sorted_.size();
        epi_base_sorted.resizeNoInitialize(n_epi_sorted);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_epi_sorted; i++){
            epi_base_sorted[i].pos      = epi_sorted_[i].getPos();
            epi_base_sorted[i].r_search = epi_sorted_[i].getRSearch();
        }
        if(!flag_reuse){
            ////////////
            // 1st STEP (send epi_base)
            const F64 time_offset_find_particle = GetWtime();
            FindScatterParticleP2P<TreeCellLoc, TreeParticle, EssentialParticleBase, Tspj, CALC_DISTANCE_TYPE>
                (tc_loc_, tp_glb_, epi_base_sorted,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st,
                 comm_table_.rank_send_,
                 comm_table_.rank_recv_);
            time_profile_.make_LET_1st__find_particle += GetWtime() - time_offset_find_particle;
            const F64 time_offset_exchange_n = GetWtime();
            ExchangeNumberShort(n_ep_send_per_proc_1st,
                                n_ep_recv_per_proc_1st,
                                comm_table_.rank_send_,
                                comm_table_.rank_recv_,
                                comm_info_);
            time_profile_.make_LET_1st__exchange_n += GetWtime() - time_offset_exchange_n;
            ExchangeLetAllgather<TSM, EssentialParticleBase>
                (epi_base_sorted,
                 n_ep_send_per_proc_1st,
                 n_ep_recv_per_proc_1st,
                 n_ep_send_per_image_1st,
                 adr_ep_send_1st,
                 ep_recv_1st,
                 0,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 comm_table_.rank_send_,
                 comm_table_.rank_recv_,
                 comm_info_);
            ////////////
            // 2nd STEP (find j particle)
            FindExchangeParticleDoubleWalk<TreeCellLoc, TreeParticle, EssentialParticleBase>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epi_base_sorted,
                 center_, length_, morton_key_);
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumberShort(n_ep_send_per_proc_2nd,
                                n_ep_recv_per_proc_2nd,
                                comm_table_.rank_send_,
                                comm_table_.rank_recv_,
                                comm_info_);
            /////////////////////
            // 4th STEP (make communication table)
            const auto n_proc = comm_info_.getNumberOfProc();
            MakeCommTableFor2StepCommuniction(comm_table_,
                                              n_ep_send_per_proc_1st,
                                              n_image_per_proc_1st,
                                              shift_per_image_1st,
                                              n_ep_send_per_image_1st,
                                              adr_ep_send_1st,
                                              n_ep_recv_per_proc_1st,
                                              n_ep_send_per_proc_2nd,
                                              n_image_per_proc_2nd,
                                              shift_per_image_2nd,
                                              n_ep_send_per_image_2nd,
                                              adr_ep_send_2nd,
                                              n_ep_recv_per_proc_2nd,
                                              n_proc);
        } // end of reuse flag
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();
        ExchangeLetAllgather<TSM, Tepj>
            (epj_sorted_,
             comm_table_.n_ep_send_,
             comm_table_.n_ep_recv_,
             comm_table_.n_ep_per_image_,
             comm_table_.adr_ep_send_,
             epj_org_,
             n_loc_tot_,
             comm_table_.shift_per_image_,
             comm_table_.n_image_per_proc_,
             comm_table_.rank_send_,
             comm_table_.rank_recv_,
             comm_info_);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }
    
}
