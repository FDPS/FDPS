////////////////////////////////////////////////
/// implementaion of methods of TreeForForce ///

#include"tree_walk.hpp"

namespace ParticleSimulator{
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    Tepj * TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getEpjFromId(const S64 id, const Tepj * epj_tmp){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
        {
            if(map_id_to_epj_.empty()){
                S64 n_epj = epj_sorted_.size();
                for(S32 i=0; i<n_epj; i++){
                    if(GetMSB(tp_glb_[i].adr_ptcl_) == 1) continue;
                    Tepj * epj_tmp = epj_sorted_.getPointer(i);
                    S64 id_tmp = epj_tmp->getId();
                    map_id_to_epj_.insert( std::pair<S64, Tepj*>(id_tmp, epj_tmp) );
                }
            }
        }
        Tepj * epj = NULL;
        typename MyMap::iterator it = map_id_to_epj_.find(id);
        if( it != map_id_to_epj_.end() ) epj = it->second;
        return epj;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::getMemSizeUsed() const {
        size_t tmp = 0;
        for(int i=0; i<Comm::getNumberOfThread(); i++){
            tmp =
                epj_for_force_[i].getMemSize() + spj_for_force_[i].getMemSize()
                + epjr_send_buf_[i].getMemSize() + epjr_send_buf_for_scatter_[i].getMemSize() 
                + epjr_recv_1st_sorted_[i].getMemSize();
        }
        return tmp 
            + tp_glb_.getMemSize()
            + tc_loc_.getMemSize()
            + tc_glb_.getMemSize()
            + epi_sorted_.getMemSize() + epi_org_.getMemSize()
            + epj_sorted_.getMemSize() + epj_org_.getMemSize()
            + spj_sorted_.getMemSize() + spj_org_.getMemSize()
            + ipg_.getMemSize()
            + epj_send_.getMemSize()
            + spj_send_.getMemSize()
            + force_sorted_.getMemSize() + force_org_.getMemSize()
            + epjr_sorted_.getMemSize() + epjr_send_.getMemSize()
            + epjr_recv_.getMemSize() + epjr_recv_1st_buf_.getMemSize()
            + epjr_recv_2nd_buf_.getMemSize();
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::getUsedMemorySize() const {
        return getMemSizeUsed();
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    initialize(const U64 n_glb_tot,
               const F64 theta,
               const U32 n_leaf_limit,
               const U32 n_group_limit) {
        if(is_initialized_ == true){
            PARTICLE_SIMULATOR_PRINT_ERROR("Do not initialize the tree twice");
            std::cerr<<"SEARCH_MODE: "<<typeid(TSM).name()<<std::endl;
            std::cerr<<"Force: "<<typeid(Tforce).name()<<std::endl;
            std::cerr<<"EPI: "<<typeid(Tepi).name()<<std::endl;
            std::cerr<<"EPJ: "<<typeid(Tepj).name()<<std::endl;
            std::cerr<<"SPJ: "<<typeid(Tspj).name()<<std::endl;
            Abort(-1);
        }
        map_id_to_epj_.clear();
        is_initialized_ = true;
        n_glb_tot_ = n_glb_tot;
        theta_ = theta;
        n_leaf_limit_ = n_leaf_limit;
        n_group_limit_ = n_group_limit;
        lev_max_loc_ = lev_max_glb_ = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        const S64 n_proc = Comm::getNumberOfProc();
        n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
        wtime_exlet_comm_ = wtime_exlet_a2a_ = wtime_exlet_a2av_ = 0.0;
        wtime_walk_LET_1st_ = wtime_walk_LET_2nd_ = 0.0;
        ni_ave_ = nj_ave_ = 0;
        Comm::barrier();
        bool err = false;
        if(n_leaf_limit_ <= 0){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of the particles in the leaf cell must be > 0");
            std::cout<<"n_leaf_limit_= "<<n_leaf_limit_<<std::endl;
        }
        if(n_group_limit_ < n_leaf_limit_){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of particles in ip graoups msut be >= that in leaf cells");
            std::cout<<"n_group_limit_= "<<n_group_limit_<<std::endl;
            std::cout<<"n_leaf_limit_= "<<n_leaf_limit_<<std::endl;
        }
        if( typeid(TSM) == typeid(SEARCH_MODE_LONG) || 
            typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) ){
            if(theta_ < 0.0){
                err = true;
                PARTICLE_SIMULATOR_PRINT_ERROR("The opening criterion of the tree must be >= 0.0");
                std::cout<<"theta_= "<<theta_<<std::endl;
            }
        }
        if(err){
            std::cout<<"SEARCH_MODE: "<<typeid(TSM).name()<<std::endl;
            std::cout<<"Force: "<<typeid(Tforce).name()<<std::endl;
            std::cout<<"EPI: "<<typeid(Tepi).name()<<std::endl;
            std::cout<<"SPJ: "<<typeid(Tspj).name()<<std::endl;
            ParticleSimulator::Abort(-1);
        }

        /*
        epi_sorted_.setAllocMode(1); // msortLT --- final
        spj_sorted_.setAllocMode(1); //  --- final
        epi_org_.setAllocMode(1); // setPtclLT ---
        spj_org_.setAllocMode(1); // insted of it, use spj_recv
        epj_org_.setAllocMode(1);
        force_sorted_.setAllocMode(1); // -- final
        epj_send_.setAllocMode(1);
        spj_send_.setAllocMode(1);
        
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].setAllocMode(1);
            spj_for_force_[i].setAllocMode(1);
        }
        */

#ifdef PARTICLE_SIMULATOR_USE_MEMORY_POOL
        epi_sorted_.setAllocMode(1); // msortLT --- final
        spj_sorted_.setAllocMode(1); //  --- final
        epi_org_.setAllocMode(1); // setPtclLT ---
        spj_org_.setAllocMode(1); // insted of it, use spj_recv
        epj_org_.setAllocMode(1);
        force_sorted_.setAllocMode(1); // -- final
        epj_send_.setAllocMode(1);
        spj_send_.setAllocMode(1);
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].setAllocMode(1);
            spj_for_force_[i].setAllocMode(1);
        }
#else
        epi_org_.setAllocMode(1);
        epj_send_.setAllocMode(1);
        spj_send_.setAllocMode(1);
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].setAllocMode(1);
            spj_for_force_[i].setAllocMode(1);
        }
#endif

        
        adr_epj_for_force_ = new ReallocatableArray<S32>[n_thread];
        adr_spj_for_force_ = new ReallocatableArray<S32>[n_thread];
        adr_ipg_for_force_ = new ReallocatableArray<S32>[n_thread];
        
        epjr_send_buf_ = new ReallocatableArray<EPJWithR>[n_thread];
        epjr_send_buf_for_scatter_ = new ReallocatableArray<EPJWithR>[n_thread];
        epjr_recv_1st_sorted_ = new ReallocatableArray<EPJWithR>[n_thread];
        epj_neighbor_ = new ReallocatableArray<Tepj>[n_thread];
        n_cell_open_ = new CountT[n_thread];
        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        req_send_ = new MPI_Request[n_proc];
        req_recv_ = new MPI_Request[n_proc];
        status_   = new MPI_Status[n_proc];
#endif
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree= "<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tpsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setParticleLocalTreeImpl(const Tpsys & psys,
                             const bool clear){
        const F64 time_offset = GetWtime();
        const S32 nloc = psys.getNumberOfParticleLocal();
        if(clear){ n_loc_tot_ = 0;}
        //        const S32 offset = 0;
        const S32 offset = n_loc_tot_;
        n_loc_tot_ += nloc;
        epj_org_.resizeNoInitialize(n_loc_tot_);
        epi_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<nloc; i++){
                epi_org_[i].copyFromFP( psys[i] );
                epj_org_[i].copyFromFP( psys[i] );
            }
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<nloc; i++){
                epi_org_[i+offset].copyFromFP( psys[i] );
                epj_org_[i+offset].copyFromFP( psys[i] );
            }
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"nloc="<<nloc<<std::endl;
        std::cout<<"n_loc_tot_="<<n_loc_tot_<<std::endl;
#endif
        time_profile_.set_particle_local_tree += GetWtime() - time_offset;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setRootCell(const DomainInfo & dinfo){
        const F64 time_offset = GetWtime();
        if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
            calcCenterAndLengthOfRootCellOpenImpl(typename TSM::search_type());
        }
        else{
            calcCenterAndLengthOfRootCellPeriodicImpl(typename TSM::search_type());
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"length_="<<length_<<" center_="<<center_<<std::endl;
            std::cout<<"pos_root_cell_="<<pos_root_cell_<<std::endl;
#endif
        time_profile_.set_root_cell += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setRootCell(const F64 l, const F64vec & c){
        const F64 time_offset = GetWtime();
        center_ = c;
        length_ = l;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
        time_profile_.set_root_cell += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortLocalTreeOnly(const bool reuse){
        const F64 wtime_offset = GetWtime();
        epi_sorted_.resizeNoInitialize(n_loc_tot_);
        epj_sorted_.resizeNoInitialize(n_loc_tot_);
        adr_org_from_adr_sorted_loc_.resizeNoInitialize(n_loc_tot_);
        tp_glb_.resizeNoInitialize(n_loc_tot_);
        if(!reuse){
            //tp_loc_.resizeNoInitialize(n_loc_tot_);
            //tp_glb_.resizeNoInitialize(n_loc_tot_);
            ReallocatableArray<TreeParticle> tp_buf(n_loc_tot_, n_loc_tot_, 1);
            MortonKey::initialize( length_ * 0.5, center_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //tp_loc_[i].setFromEP(epj_org_[i], i);
                tp_glb_[i].setFromEP(epj_org_[i], i);
            }
#ifdef USE_STD_SORT
            //std::sort(tp_loc_.getPointer(), tp_loc_.getPointer()+n_loc_tot_, 
            //          [](const TreeParticle & l, const TreeParticle & r )
            //          ->bool{return l.getKey() < r.getKey();} );
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_loc_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );            

#else
            //rs_.lsdSort(tp_loc_.getPointer(), tp_buf.getPointer(), 0, n_loc_tot_-1);
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf.getPointer(), 0, n_loc_tot_-1);
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //const S32 adr = tp_loc_[i].adr_ptcl_;
                const S32 adr = tp_glb_[i].adr_ptcl_;
                adr_org_from_adr_sorted_loc_[i] = adr;
            }
            tp_buf.freeMem(1);
        } // end of if() no reuse
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = adr_org_from_adr_sorted_loc_[i];
            epi_sorted_[i] = epi_org_[adr];
            epj_sorted_[i] = epj_org_[adr];
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epi_sorted_.size()="<<epi_sorted_.size()<<" epj_sorted_.size()="<<epj_sorted_.size()<<std::endl;
#endif
        time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset;
        time_profile_.make_local_tree += GetWtime() - wtime_offset;

        /*
        if(Comm::getRank()==0){
            std::cerr<<"epi_sorted_.size()= "<<epi_sorted_.size()
                     <<" epi_org_.size()= "<<epi_org_.size()
                     <<" epj_sorted_.size()= "<<epj_sorted_.size()
                     <<" epj_org_.size()= "<<epj_org_.size()
                     <<std::endl;
            for(S32 i=0; i<epi_sorted_.size(); i += epi_sorted_.size()/10){
                std::cerr<<"i= "<<i
                         <<" epi_sorted_[i].id= "<<epi_sorted_[i].id
                         <<" epi_sorted_[i].pos= "<<epi_sorted_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<epj_sorted_.size(); i += epj_sorted_.size()/10){
                std::cerr<<"i= "<<i
                         <<"epj_sorted_[i].id= "<<epj_sorted_[i].id
                         <<" epj_sorted_[i].pos= "<<epj_sorted_[i].pos
                         <<std::endl;
            }
        }
        */

    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellLocalTreeOnly(){
        const F64 time_offset = GetWtime();
        LinkCell(tc_loc_,  adr_tc_level_partition_loc_, tp_glb_.getPointer(), lev_max_loc_, n_loc_tot_, n_leaf_limit_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tc_loc_.size()="<<tc_loc_.size()<<std::endl;
        std::cout<<"lev_max_loc_="<<lev_max_loc_<<std::endl;
#endif
        time_profile_.link_cell_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnly(){
        F64 time_offset = GetWtime();
        calcMomentLocalTreeOnlyImpl(typename TSM::search_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        time_profile_.calc_moment_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree_tot = time_profile_.calc_moment_local_tree + time_profile_.make_local_tree;
    }


    // FOR P^3T
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongSymmetry){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
    
    // FOR P^3T
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongCutoffScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLong){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongCutoff){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortGather){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epi_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortSymmetry){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epi_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    ///////////////////////////////
    /// morton sort global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortGlobalTreeOnly(const bool reuse){
        F64 time_offset = GetWtime();
        //ReallocatableArray< TreeParticle > tp_buf;
        if(!map_id_to_epj_.empty()){
            map_id_to_epj_.clear();
        }
        assert(map_id_to_epj_.empty());
        tp_glb_.resizeNoInitialize(n_glb_tot_);
        const S32 n_ep_tot = epj_org_.size();
        epj_sorted_.resizeNoInitialize(n_ep_tot);
        adr_org_from_adr_sorted_glb_.resizeNoInitialize(n_glb_tot_);
        if(!reuse){
            ReallocatableArray<TreeParticle> tp_buf(n_glb_tot_, n_glb_tot_, 1);
            //tp_buf.resizeNoInitialize(n_glb_tot_);
#ifdef USE_STD_SORT
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );
#else
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf.getPointer(), 0, n_glb_tot_-1);
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                adr_org_from_adr_sorted_glb_[i] = tp_glb_[i].adr_ptcl_;
            }
            tp_buf.freeMem(1);
        }
        if( typeid(typename TSM::force_type) == typeid(TagForceLong)){
            //const S32 n_sp_tot = spj_recv_.size();
            const S32 n_sp_tot = spj_org_.size();
            assert(n_ep_tot+n_sp_tot == n_glb_tot_);
            
            spj_sorted_.resizeNoInitialize(n_sp_tot);

#if 0
            // thread parallelized, but not fast.
            const S32 n_thread = Comm::getNumberOfThread();
            static ReallocatableArray<U32> n_cnt_ep;
            static ReallocatableArray<U32> n_cnt_sp;
            n_cnt_ep.resizeNoInitialize(n_thread);
            n_cnt_sp.resizeNoInitialize(n_thread);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
#endif
            {
                S32 id_thread = Comm::getThreadNum();
                n_cnt_ep[id_thread] = n_cnt_sp[id_thread] = 0;
                U32 id_tp_head = (n_glb_tot_/n_thread) * id_thread + std::min(id_thread, (S32)n_glb_tot_%n_thread);
                U32 id_tp_tail = (n_glb_tot_/n_thread) * (id_thread+1) + std::min( (id_thread+1), (S32)n_glb_tot_%n_thread);
                for(U32 i=id_tp_head; i<id_tp_tail; i++){
                    const U32 adr = adr_org_from_adr_sorted_glb_[i];
                    //const U32 adr = tp_glb_[i].adr_ptcl_;
                    if( GetMSB(adr) == 0){
                        n_cnt_ep[id_thread]++;
                    }
                    else{
                        n_cnt_sp[id_thread]++;
                    }
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp barrier
#endif
                U32 id_ep_head = 0;
                U32 id_sp_head = 0;
                for(U32 i=0; i<id_thread; i++){
                    id_ep_head += n_cnt_ep[i];
                    id_sp_head += n_cnt_sp[i];
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp barrier
#endif
                n_cnt_ep[id_thread] = n_cnt_sp[id_thread] = 0;
                for(U32 i=id_tp_head; i<id_tp_tail; i++){
                    const U32 adr = adr_org_from_adr_sorted_glb_[i];
                    //const U32 adr = tp_glb_[i].adr_ptcl_;
                    if( GetMSB(adr) == 0){
                        const U32 id_dst = id_ep_head + n_cnt_ep[id_thread];
                        epj_sorted_[id_dst] = epj_org_[adr];
                        tp_glb_[i].adr_ptcl_ = id_dst;
                        n_cnt_ep[id_thread]++;
                    }
                    else{
                        const U32 id_dst = id_sp_head + n_cnt_sp[id_thread];
                        spj_sorted_[id_dst] = spj_org_[ClearMSB(adr)];
                        tp_glb_[i].adr_ptcl_ = SetMSB(id_dst);
                        n_cnt_sp[id_thread]++;
                    }
                }
            }
#else
            U32 n_cnt_ep = 0;
            U32 n_cnt_sp = 0;
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = adr_org_from_adr_sorted_glb_[i];
                if( GetMSB(adr) == 0){
                    epj_sorted_[n_cnt_ep] = epj_org_[adr];
                    tp_glb_[i].adr_ptcl_ = n_cnt_ep;
                    n_cnt_ep++;
                }
                else{
                    spj_sorted_[n_cnt_sp] = spj_org_[ClearMSB(adr)];
                    tp_glb_[i].adr_ptcl_ = SetMSB(n_cnt_sp);
                    n_cnt_sp++;
                }
            }
#endif
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = adr_org_from_adr_sorted_glb_[i];
                epj_sorted_[i] = epj_org_[adr];
                tp_glb_[i].adr_ptcl_ = i;
            }
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tp_glb_.size()="<<tp_glb_.size()<<std::endl;
        std::cout<<"epj_sorted_.size()="<<epj_sorted_.size()<<" spj_sorted_.size()="<<spj_sorted_.size()<<std::endl;
#endif
        time_profile_.morton_sort_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree += GetWtime() - time_offset;
    }

    /////////////////////////////
    /// link cell global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellGlobalTreeOnly(){
        const F64 time_offset = GetWtime();
        LinkCell(tc_glb_, adr_tc_level_partition_glb_,
                 tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_, n_leaf_limit_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tc_glb_.size()="<<tc_glb_.size()<<std::endl;
        std::cout<<"lev_max_glb_="<<lev_max_glb_<<std::endl;
#endif
        time_profile_.link_cell_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree += GetWtime() - time_offset;
    }

    
    //////////////////////////
    // CALC MOMENT GLOBAL TREE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnly(){
        const F64 time_offset = GetWtime();
        calcMomentGlobalTreeOnlyImpl(typename TSM::search_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        time_profile_.calc_moment_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree_tot = time_profile_.calc_moment_global_tree + time_profile_.make_global_tree;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLong){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongScatter){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongSymmetry){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }    

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongCutoff){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortScatter){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortGather){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortSymmetry){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    }
    
    ////////////////////
    /// MAKE IPGROUP ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroup(){
        const F64 time_offset = GetWtime();
        ipg_.clearSize();
        makeIPGroupImpl(typename TSM::force_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
        time_profile_.calc_force__make_ipgroup += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroupImpl(TagForceLong){
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_group_limit_="<<n_group_limit_<<std::endl;
#endif

#ifdef PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
        MakeIPGroupLongGLBTreeCellAsIPGBox(ipg_, tc_loc_, tc_glb_, epj_sorted_, 0, 0, n_group_limit_, n_leaf_limit_);
#else //PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
#if 1
        MakeIPGroupUseGLBTreeLong(ipg_, tc_loc_, tc_glb_, epi_sorted_, 0, 0, n_group_limit_, n_leaf_limit_); // NEW
#else
        MakeIPGroupLong(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
#endif

        const S32 n_ipg = ipg_.size();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif //PARTICLE_SIMULATOR_THREAD_PARALLEL
        for(S32 i=0; i<n_ipg; i++){
            const S32 n = ipg_[i].n_ptcl_;
            const S32 adr = ipg_[i].adr_ptcl_;
            ipg_[i].vertex_in_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
        }
#endif //PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroupImpl(TagForceShort){
        MakeIPGroupShort(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
    }
    

    /////////////////////////////
    /// MAKE INTERACTION LIST ///
    // pj is epj or spj
    template<class Tpj>
    void CopyPjForForceST(const ReallocatableArray<S32> & adr_pj,
                          const ReallocatableArray<Tpj> & pj_sorted,
                          ReallocatableArray<Tpj> & pj_for_force){
        const S32 n_pj = adr_pj.size();
        pj_for_force.resizeNoInitialize(n_pj);
        for(S32 i=0; i<n_pj; i++){
            const S32 adr_pj_src = adr_pj[i];
            pj_for_force[i] = pj_sorted[adr_pj_src];
        }
    }
    template<class Tpj>
    void CopyPjForForceST(const ReallocatableArray<S32> & adr_pj,
                          const ReallocatableArray<Tpj> & pj_sorted,
                          const S32 n_head,
                          const S32 n_tail,
                          ReallocatableArray<Tpj> & pj_for_force){
        pj_for_force.resizeNoInitialize(n_tail);
        for(S32 i=n_head; i<n_tail; i++){
            const S32 adr_pj_src = adr_pj[i];
            pj_for_force[i] = pj_sorted[adr_pj_src];
        }
    }
    
    template<class T>
    struct TraitsForCutoff{
        typedef TagWithoutCutoff type_cutoff;
    };
    template<>
    struct TraitsForCutoff<TagSearchLongCutoff>{
        typedef TagWithCutoff type_cutoff;
    };
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListLongForZeroTheta(TagWithoutCutoff, const S32 adr_ipg){
        const S32 ith = Comm::getThreadNum();
        const S32 n_tmp = tc_glb_[0].n_ptcl_;
        S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
        epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
        for(S32 ip=0; ip<n_tmp; ip++){
            if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
            }
            else {
                adr_ptcl_tmp++;
            }
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListLongForZeroTheta(TagWithCutoff, const S32 adr_ipg){
        const S32 ith = Comm::getThreadNum();
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
            const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_in_;
            const F64ort cell_box = pos_root_cell_;
            MakeInteractionListLongCutoffEPForZeroTheta
                (tc_glb_, tc_glb_[0].adr_tc_, tp_glb_, 
                 epj_sorted_, epj_for_force_[ith],
                 cell_box,
                 pos_target_box, r_cut_sq, n_leaf_limit_); 
        }
        else{
            const S32 n_tmp = tc_glb_[0].n_ptcl_;
            S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
            epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
            for(S32 ip=0; ip<n_tmp; ip++){
                if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                    epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                }
                else {
                    adr_ptcl_tmp++;
                }
            }
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionList(const S32 adr_ipg, const bool clear){
        makeInteractionListImpl(typename TSM::force_type(), adr_ipg, clear);
    }




    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcCenterAndLengthOfRootCellOpenNoMargenImpl(const Tep2 ep[]){

        const F64ort min_box  = GetMinBox(ep, n_loc_tot_);
        center_ = min_box.getCenter();
        const F64 tmp0 = (min_box.high_ - center_).getMax();
        const F64 tmp1 = (center_ - min_box.low_).getMax();
        length_ = std::max(tmp0, tmp1) * 2.0 * 1.001;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcCenterAndLengthOfRootCellOpenWithMargenImpl(const Tep2 ep[]){
        const F64ort min_box  = GetMinBoxWithMargen(ep, n_loc_tot_);
        center_ = min_box.getCenter();
        const F64 tmp0 = (min_box.high_ - center_).getMax();
        const F64 tmp1 = (center_ - min_box.low_).getMax();
        length_ = std::max(tmp0, tmp1) * 2.0 * 1.001;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcCenterAndLengthOfRootCellPeriodicImpl2(const Tep2 ep[]){
        //F64 rsearch_max_loc = std::numeric_limits<F64>::max() * -0.25;
        F64 rsearch_max_loc = -LARGE_FLOAT;
        F64ort box_loc;
        box_loc.initNegativeVolume();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            //F64 rsearch_max_loc_tmp = std::numeric_limits<F64>::max() * -0.25;
            F64 rsearch_max_loc_tmp = -LARGE_FLOAT;
            F64ort box_loc_tmp;
            box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
            for(S32 ip=0; ip<n_loc_tot_; ip++){
                rsearch_max_loc_tmp = (rsearch_max_loc_tmp > ep[ip].getRSearch()) ? rsearch_max_loc_tmp : ep[ip].getRSearch();
                box_loc_tmp.merge(ep[ip].getPos());
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
            {
                rsearch_max_loc = rsearch_max_loc > rsearch_max_loc_tmp ? rsearch_max_loc : rsearch_max_loc_tmp;
                box_loc.merge(box_loc_tmp);
            }
        }
        F64 rsearch_max_glb = 1.001 * Comm::getMaxValue(rsearch_max_loc);
        // The factor 1.001 is for the safety.
        F64vec xlow_loc = box_loc.low_;
        F64vec xhigh_loc = box_loc.high_;
        F64vec xlow_glb = Comm::getMinValue(xlow_loc);
        F64vec xhigh_glb = Comm::getMaxValue(xhigh_loc);

        xlow_glb -= F64vec(rsearch_max_glb);
        xhigh_glb += F64vec(rsearch_max_glb);
        center_ = (xhigh_glb + xlow_glb) * 0.5;
        F64 tmp0 = (xlow_glb - center_).applyEach(Abs<F64>()).getMax();
        F64 tmp1 = (xhigh_glb - center_).applyEach(Abs<F64>()).getMax();
        length_ = std::max(tmp0, tmp1) * 2.0 * 2.0;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    copyForceOriginalOrder(){
        epi_org_.freeMem(1);
        force_org_.resizeNoInitialize(n_loc_tot_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = adr_org_from_adr_sorted_loc_[i];
            force_org_[adr] = force_sorted_[i];
        }
    }
    
    // return forces in original order
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceDirect(Tfunc_ep_ep pfunc_ep_ep,
                    Tforce force[],
                    const DomainInfo & dinfo,
                    const bool clear){
        if(clear){
            for(S32 i=0; i<n_loc_tot_; i++) force[i].clear();
        }
        Tepj * epj_tmp;
        S32 n_epj_tmp;
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain(), pos_root_cell_, pa);
        pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force);
        delete [] epj_tmp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceDirectAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                const DomainInfo & dinfo,
                                const bool clear){
        force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
            for(S32 i=0; i<n_loc_tot_; i++)force_org_[i].clear();
        }
        Tepj * epj_tmp;
        S32 n_epj_tmp;
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa);
        pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force_org_.getPointer());
        delete [] epj_tmp;
    }

#if 1
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        return getNeighborListOneParticleImpl(typename TSM::neighbor_search_type(), ptcl, epj);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagNeighborSearchScatter, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagNeighborSearchGather, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        const S32 adr = 0;
        SearchNeighborListOneParticleGather(pos_target,
                                            r_search_sq,
                                            tc_glb_.getPointer(),
                                            tp_glb_.getPointer(), adr, 
                                            epj_sorted_, epj_neighbor_[id_thread], n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagNeighborSearchSymmetry, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        const S32 adr = 0;
        SearchNeighborListOneParticleSymmetry(pos_target,
                                              r_search_sq,
                                              tc_glb_.getPointer(),
                                              tp_glb_.getPointer(), adr, 
                                              epj_sorted_, epj_neighbor_[id_thread],
                                              n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagNeighborSearchNo, const Tptcl & ptcl, Tepj * & epj){
        return -1;
        // std::cerr<<"not implemented"<<std::endl;
    }

#else
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        return getNeighborListOneParticleImpl(typename TSM::search_type(), ptcl, epj);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchShortScatter, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchLongScatter, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleScatter(pos_target, tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_, epj_neighbor_[id_thread], n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchShortSymmetry, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        const S32 adr = 0;
        SearchNeighborListOneParticleSymmetry(pos_target,
                                              r_search_sq,
                                              tc_glb_.getPointer(),
                                              tp_glb_.getPointer(), adr, 
                                              epj_sorted_, epj_neighbor_[id_thread],
                                              n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchLongSymmetry, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        const F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        SearchNeighborListOneParticleSymmetry(pos_target,
                                              r_search_sq,
                                              tc_glb_.getPointer(),
                                              tp_glb_.getPointer(), adr,
                                              epj_sorted_,   epj_neighbor_[id_thread],
                                              n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchShortGather, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        const S32 adr = 0;
        SearchNeighborListOneParticleGather(pos_target,
                                            r_search_sq,
                                            tc_glb_.getPointer(),
                                            tp_glb_.getPointer(), adr, 
                                            epj_sorted_, epj_neighbor_[id_thread], n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
#endif

    

    
    // 2016 02/05
    /*
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleScatter(pos_target,  tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_, epj_neighbor_[id_thread], n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    */
    /*
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){

        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleGather(pos_target,  tc_glb_.getPointer(),
                                            tp_glb_.getPointer(), adr, 
                                            epj_sorted_, epj_neighbor_[id_thread], n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;

    }
    */
    /*
    template<class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<SEARCH_MODE_SYMMETRY, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleSymmetry(pos_target,  tc_glb_.getPointer(),
                                              tp_glb_.getPointer(), adr, 
                                              epj_sorted_, epj_neighbor_[id_thread], n_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }    
    */


    
#if 0
    // under construction
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneIPGroupImpl(TagSearchLongScatter, const Tptcl & ptcl, S32 & nnp){
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        const S32 size_old = epj_neighbor_.size();
        SearchNeighborListOneIPGroupScatter(pos_target,    tc_glb_.getPointer(),       adr, 
                                            epj_sorted_,   epj_neighbor_, n_leaf_limit_);
        nnp = epj_neighbor_.size() - size_old;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneIPGroup(const S32 iipg, S32 & nip, 
                              const Tepi * & epi, S32 & nnp, Tepj * & epj){
        nip = ipg_[iipg].n_ptcl_;
        const S32 adr = ipg_[iipg].adr_ptcl_;
        epi = epi_sorted_[adr].getPointer();
        const S32 head = epj_neighbor_.size();
        getNeighborListOneIPGroupImpl(typename TSM::search_type(), ptcl, nnp);
        epj = epj_neighbor_.getPointer(head);
    }
#endif

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    clearSizeOfArray(){
        tp_glb_.clearSize();
        tc_loc_.clearSize();
        tc_glb_.clearSize();
        epi_sorted_.clearSize();
        epi_org_.clearSize();
        epj_sorted_.clearSize();
        epj_org_.clearSize();
        spj_sorted_.clearSize();
        spj_org_.clearSize();
        ipg_.clearSize();
        epj_send_.clearSize();
        spj_send_.clearSize();
        force_sorted_.clearSize();
        force_org_.clearSize();
        epjr_sorted_.clearSize();
        epjr_send_.clearSize();
        epjr_recv_.clearSize();
        epjr_recv_1st_buf_.clearSize();
        epjr_recv_2nd_buf_.clearSize();
        const S32 n_thread = Comm::getNumberOfThread();
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].clearSize();
            spj_for_force_[i].clearSize();
            epjr_send_buf_[i].clearSize();
            epjr_send_buf_for_scatter_[i].clearSize();
            epjr_recv_1st_sorted_[i].clearSize();
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getOuterBoundaryOfLocalTreeImpl(TagSearchLongSymmetry){
        return tc_loc_[0].mom_.vertex_out_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getInnerBoundaryOfLocalTreeImpl(TagSearchLongSymmetry){
        return tc_loc_[0].mom_.vertex_in_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getOuterBoundaryOfLocalTreeImpl(TagSearchShortSymmetry){
        return tc_loc_[0].mom_.vertex_out_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getInnerBoundaryOfLocalTreeImpl(TagSearchShortSymmetry){
        return tc_loc_[0].mom_.vertex_in_;
    }
    
    ////////////////////
    // for reuse
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseList(const DomainInfo & dinfo,const bool flag_reuse){
        if(typeid(TSM) == typeid(SEARCH_MODE_LONG)
           && dinfo.getBoundaryCondition() != BOUNDARY_CONDITION_OPEN){
            PARTICLE_SIMULATOR_PRINT_ERROR("The forces w/o cutoff can be evaluated only under the open boundary condition");
            Abort(-1);
        }
        if(!flag_reuse){ comm_table_.clearSize(); }
        //comm_table_.clear();
        exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo, flag_reuse);
        time_profile_.exchange_LET_tot = time_profile_.make_LET_1st
            + time_profile_.exchange_LET_1st
            + time_profile_.make_LET_2nd
            + time_profile_.exchange_LET_2nd;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListLong(const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 r_crit_sq = LARGE_FLOAT;
        if (theta_ > 0) r_crit_sq = (length_ * length_) / (theta_ * theta_);
        else r_crit_sq = - 1.0; // negative value is used to represent theta_ = 0
        if(!flag_reuse){
            FindScatterParticle<TSM, TreeCell<Tmomloc>, TreeParticle,
                                Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_glb_,
                 epj_sorted_,
                 comm_table_.n_ep_send_,   comm_table_.adr_ep_send_, 
                 dinfo,          n_leaf_limit_,
                 comm_table_.n_sp_send_,   comm_table_.adr_sp_send_, 
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.n_sp_per_image_,
                 r_crit_sq);            
            ExchangeNumber(comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                           comm_table_.n_sp_send_, comm_table_.n_sp_recv_);
            const S32 n_proc = Comm::getNumberOfProc();
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
            comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
            S32 n_sp_recv_tot_tmp = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
                n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ += n_ep_recv_tot_tmp;
            comm_table_.n_sp_recv_tot_ += n_sp_recv_tot_tmp;
        }
        ExchangeLet<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_,
                                     comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                     comm_table_.adr_ep_send_,
                                     //epj_recv_,
                                     epj_org_, n_loc_tot_,
                                     tc_loc_, comm_table_.n_sp_send_,
                                     comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
                                     comm_table_.adr_sp_send_,
                                     //spj_recv_,
                                     spj_org_,
                                     comm_table_.shift_per_image_,
                                     comm_table_.n_image_per_proc_);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLong,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 time_offset = GetWtime();
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }    
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLongCutoff,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLongScatter,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLongSymmetry,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchShortScatter,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 time_offset = GetWtime();

        if(!flag_reuse){
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_glb_, epj_sorted_,
                 comm_table_.n_ep_send_,  comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_);

            ExchangeNumber(comm_table_.n_ep_send_, comm_table_.n_ep_recv_);

            const S32 n_proc = Comm::getNumberOfProc();
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
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_,
                    //epj_recv_,
                    epj_org_, n_loc_tot_, 
                    comm_table_.shift_per_image_,
                    comm_table_.n_image_per_proc_);

        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }

    ////////////////
    // SYMMETRY
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchShortSymmetry,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<S32> n_ep_send_per_proc_1st;
        static ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        static ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        static ReallocatableArray<S32> adr_ep_send_1st;
        static ReallocatableArray<S32> adr_ep_send_2nd;
        static ReallocatableArray<F64vec> shift_per_image_1st;
        static ReallocatableArray<F64vec> shift_per_image_2nd;
        static ReallocatableArray<S32> n_image_per_proc_1st;
        static ReallocatableArray<S32> n_image_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_send_per_image_1st;
        static ReallocatableArray<S32> n_ep_send_per_image_2nd;
        static ReallocatableArray<Tepj> ep_recv_1st;
        
        if(!flag_reuse){
            ////////////
            // 1st STEP (send j particles)
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_glb_, epj_sorted_,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            ExchangeNumber(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st);
            ExchangeParticle(epj_sorted_, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            ////////////
            // 2nd STEP (send j particles)
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, Tepj>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epj_sorted_,
                 center_, length_);
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumber(n_ep_send_per_proc_2nd, n_ep_recv_per_proc_2nd);

            ReallocatableArray<S32> n_disp_ep_send_per_proc;
            n_disp_ep_send_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_ep_send_per_proc[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_per_proc[i+1] = n_disp_ep_send_per_proc[i]
                    + n_ep_send_per_proc_1st[i]
                    + n_ep_send_per_proc_2nd[i];
            }

            ReallocatableArray<S32> n_disp_image_per_proc;
            ReallocatableArray<S32> n_disp_image_per_proc_1st;
            ReallocatableArray<S32> n_disp_image_per_proc_2nd;
            n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_1st.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_2nd.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc[0] = n_disp_image_per_proc_1st[0] = n_disp_image_per_proc_2nd[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_image_per_proc[i+1]     = n_disp_image_per_proc[i]     + n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
                n_disp_image_per_proc_1st[i+1] = n_disp_image_per_proc_1st[i] + n_image_per_proc_1st[i];
                n_disp_image_per_proc_2nd[i+1] = n_disp_image_per_proc_2nd[i] + n_image_per_proc_2nd[i];
            }

            const S32 n_image_tot_1st = shift_per_image_1st.size();
            const S32 n_image_tot_2nd = shift_per_image_2nd.size();
            ReallocatableArray<S32> n_disp_ep_send_per_image_1st;
            ReallocatableArray<S32> n_disp_ep_send_per_image_2nd;
            n_disp_ep_send_per_image_1st.resizeNoInitialize(n_image_tot_1st+1);
            n_disp_ep_send_per_image_2nd.resizeNoInitialize(n_image_tot_2nd+1);
            n_disp_ep_send_per_image_1st[0] = n_disp_ep_send_per_image_2nd[0] = 0;
            for(S32 i=0; i<n_image_tot_1st; i++){
                n_disp_ep_send_per_image_1st[i+1] = n_disp_ep_send_per_image_1st[i] + n_ep_send_per_image_1st[i];
            }
            for(S32 i=0; i<n_image_tot_2nd; i++){
                n_disp_ep_send_per_image_2nd[i+1] = n_disp_ep_send_per_image_2nd[i] + n_ep_send_per_image_2nd[i];
            }
            
            const S32 n_image_tot = n_disp_image_per_proc_1st[n_proc]  + n_disp_image_per_proc_2nd[n_proc];
            comm_table_.shift_per_image_.resizeNoInitialize(n_image_tot);
            comm_table_.n_ep_per_image_.resizeNoInitialize(n_image_tot);
            const S32 n_send_tot = n_disp_ep_send_per_proc[n_proc];
            comm_table_.adr_ep_send_.resizeNoInitialize(n_send_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                S32 n_ep_cnt = 0;
                S32 n_image_cnt = 0;
                const S32 ep_head    = n_disp_ep_send_per_proc[i];
                const S32 image_head = n_disp_image_per_proc[i];
                const S32 image_head_1st = n_disp_image_per_proc_1st[i];
                const S32 image_end_1st  = n_disp_image_per_proc_1st[i+1];
                for(S32 j=image_head_1st; j<image_end_1st; j++, n_image_cnt++){
                    const S32 ep_head_1st = n_disp_ep_send_per_image_1st[j];
                    const S32 ep_end_1st  = n_disp_ep_send_per_image_1st[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_1st[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_1st[j];
                    for(S32 k=ep_head_1st; k<ep_end_1st; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_1st[k];
                    }
                }
                const S32 image_head_2nd = n_disp_image_per_proc_2nd[i];
                const S32 image_end_2nd  = n_disp_image_per_proc_2nd[i+1];
                for(S32 j=image_head_2nd; j<image_end_2nd; j++, n_image_cnt++){
                    const S32 ep_head_2nd = n_disp_ep_send_per_image_2nd[j];
                    const S32 ep_end_2nd  = n_disp_ep_send_per_image_2nd[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_2nd[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_2nd[j];
                    for(S32 k=ep_head_2nd; k<ep_end_2nd; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_2nd[k];
                    }
                }
            }
            comm_table_.n_ep_send_.resizeNoInitialize(n_proc);
            comm_table_.n_ep_recv_.resizeNoInitialize(n_proc);
            comm_table_.n_image_per_proc_.resizeNoInitialize(n_proc);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                comm_table_.n_ep_send_[i] = n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
                comm_table_.n_ep_recv_[i] = n_ep_recv_per_proc_1st[i] + n_ep_recv_per_proc_2nd[i];
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_image_tot; i++){
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            S32 n_ep_recv_tot_tmp = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:n_ep_recv_tot_tmp)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        } // end of reuse
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_,
                    //epj_recv_,
                    epj_org_, n_loc_tot_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);
    }


    ////////////////
    // GATHER MODE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchShortGather,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<S32> n_ep_send_per_proc_1st;
        static ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        static ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        static ReallocatableArray<S32> adr_ep_send_1st;
        static ReallocatableArray<S32> adr_ep_send_2nd;
        static ReallocatableArray<F64vec> shift_per_image_1st;
        static ReallocatableArray<F64vec> shift_per_image_2nd;
        static ReallocatableArray<S32> n_image_per_proc_1st;
        static ReallocatableArray<S32> n_image_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_send_per_image_1st;
        static ReallocatableArray<S32> n_ep_send_per_image_2nd;
        //static ReallocatableArray<Tepi> ep_recv_1st;
        static ReallocatableArray<EssentialParticleBase> ep_recv_1st;
        static ReallocatableArray<EssentialParticleBase> epi_base_sorted;
        const S32 n_epi_sorted = epi_sorted_.size();
        epi_base_sorted.resizeNoInitialize(n_epi_sorted);
        for(S32 i=0; i<n_epi_sorted; i++){
            epi_base_sorted[i].pos      = epi_sorted_[i].getPos();
            epi_base_sorted[i].r_search = epi_sorted_[i].getRSearch();
        }
        
        if(!flag_reuse){
            ////////////
            // 1st STEP (send epi_base)
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, EssentialParticleBase, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_glb_, epi_base_sorted,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            ExchangeNumber(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st);
            // exchange epi_base particles
            ExchangeParticle(epi_base_sorted, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            
            ////////////
            // 2nd STEP (find j particle)
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, EssentialParticleBase>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epi_base_sorted,
                 center_, length_);
            
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumber(n_ep_send_per_proc_2nd, n_ep_recv_per_proc_2nd);

            ReallocatableArray<S32> n_disp_ep_send_per_proc;
            n_disp_ep_send_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_ep_send_per_proc[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_per_proc[i+1] = n_disp_ep_send_per_proc[i]
                    + n_ep_send_per_proc_1st[i]
                    + n_ep_send_per_proc_2nd[i];
            }

            ReallocatableArray<S32> n_disp_image_per_proc;
            ReallocatableArray<S32> n_disp_image_per_proc_1st;
            ReallocatableArray<S32> n_disp_image_per_proc_2nd;
            n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_1st.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_2nd.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc[0] = n_disp_image_per_proc_1st[0] = n_disp_image_per_proc_2nd[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_image_per_proc[i+1]     = n_disp_image_per_proc[i]     + n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
                n_disp_image_per_proc_1st[i+1] = n_disp_image_per_proc_1st[i] + n_image_per_proc_1st[i];
                n_disp_image_per_proc_2nd[i+1] = n_disp_image_per_proc_2nd[i] + n_image_per_proc_2nd[i];
            }

            const S32 n_image_tot_1st = shift_per_image_1st.size();
            const S32 n_image_tot_2nd = shift_per_image_2nd.size();
            ReallocatableArray<S32> n_disp_ep_send_per_image_1st;
            ReallocatableArray<S32> n_disp_ep_send_per_image_2nd;
            n_disp_ep_send_per_image_1st.resizeNoInitialize(n_image_tot_1st+1);
            n_disp_ep_send_per_image_2nd.resizeNoInitialize(n_image_tot_2nd+1);
            n_disp_ep_send_per_image_1st[0] = n_disp_ep_send_per_image_2nd[0] = 0;
            for(S32 i=0; i<n_image_tot_1st; i++){
                n_disp_ep_send_per_image_1st[i+1] = n_disp_ep_send_per_image_1st[i] + n_ep_send_per_image_1st[i];
            }
            for(S32 i=0; i<n_image_tot_2nd; i++){
                n_disp_ep_send_per_image_2nd[i+1] = n_disp_ep_send_per_image_2nd[i] + n_ep_send_per_image_2nd[i];
            }
            
            const S32 n_image_tot = n_disp_image_per_proc_1st[n_proc]  + n_disp_image_per_proc_2nd[n_proc];
            comm_table_.shift_per_image_.resizeNoInitialize(n_image_tot);
            comm_table_.n_ep_per_image_.resizeNoInitialize(n_image_tot);
            const S32 n_send_tot = n_disp_ep_send_per_proc[n_proc];
            comm_table_.adr_ep_send_.resizeNoInitialize(n_send_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                S32 n_ep_cnt = 0;
                S32 n_image_cnt = 0;
                const S32 ep_head    = n_disp_ep_send_per_proc[i];
                const S32 image_head = n_disp_image_per_proc[i];
                const S32 image_head_1st = n_disp_image_per_proc_1st[i];
                const S32 image_end_1st  = n_disp_image_per_proc_1st[i+1];
                for(S32 j=image_head_1st; j<image_end_1st; j++, n_image_cnt++){
                    const S32 ep_head_1st = n_disp_ep_send_per_image_1st[j];
                    const S32 ep_end_1st  = n_disp_ep_send_per_image_1st[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_1st[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_1st[j];
                    for(S32 k=ep_head_1st; k<ep_end_1st; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_1st[k];
                    }
                }
                const S32 image_head_2nd = n_disp_image_per_proc_2nd[i];
                const S32 image_end_2nd  = n_disp_image_per_proc_2nd[i+1];
                for(S32 j=image_head_2nd; j<image_end_2nd; j++, n_image_cnt++){
                    const S32 ep_head_2nd = n_disp_ep_send_per_image_2nd[j];
                    const S32 ep_end_2nd  = n_disp_ep_send_per_image_2nd[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_2nd[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_2nd[j];
                    for(S32 k=ep_head_2nd; k<ep_end_2nd; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_2nd[k];
                    }
                }
            }
            comm_table_.n_ep_send_.resizeNoInitialize(n_proc);
            comm_table_.n_ep_recv_.resizeNoInitialize(n_proc);
            comm_table_.n_image_per_proc_.resizeNoInitialize(n_proc);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                comm_table_.n_ep_send_[i] = n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
                comm_table_.n_ep_recv_[i] = n_ep_recv_per_proc_1st[i] + n_ep_recv_per_proc_2nd[i];
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_image_tot; i++){
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            S32 n_ep_recv_tot_tmp = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:n_ep_recv_tot_tmp)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        } // end of reuse flag
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_,
                    //epj_recv_,
                    epj_org_, n_loc_tot_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);
    }
    
    template<class TSM>
    F64ort GetIpgBoxForInteractionList(const IPGroup<TSM> & ipg){
        return ipg.vertex_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListIndexShort(){
        static bool first = true;
        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;
        static ReallocatableArray<S32> * n_disp_epj_tmp;
        const S32 n_thread = Comm::getNumberOfThread();
        if(first){
            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];
            n_disp_epj_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        const S32 n_ipg = ipg_.size();
        const S32 adr_tc = 0;
        ReallocatableArray<Tspj> spj_dummy;
        ReallocatableArray<S32> adr_spj;
        interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_ep_cum_prev = 0;
            adr_epj_tmp[ith].clearSize();
            adr_ipg_tmp[ith].clearSize();
            n_disp_epj_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 i=0; i<n_ipg; i++){
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                const F64 r_crit_sq = 999.9;
                TargetBox<TSM> target_box;
                target_box.set(ipg_[i]);
                MakeListUsingTreeRecursiveTop
                    <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj,
                     WALK_MODE_NORMAL, TagChopLeafTrue>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     target_box,
                     r_crit_sq,
                     n_leaf_limit_,
                     F64vec(0.0));

                interaction_list_.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                n_ep_cum_prev = adr_epj_tmp[ith].size();
            }
            n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
        }
        interaction_list_.n_disp_ep_[0] = 0;
        for(S32 i=0; i<n_ipg; i++){
            interaction_list_.n_disp_ep_[i+1] = interaction_list_.n_disp_ep_[i] + interaction_list_.n_ep_[i];
        }
        interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[n_ipg] );

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                const S32 adr_ipg = adr_ipg_tmp[i][j];
                S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                const S32 k_h = n_disp_epj_tmp[i][j];
                const S32 k_e = n_disp_epj_tmp[i][j+1];
                for(S32 k=k_h; k<k_e; k++, adr_ep++){
                    interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                }
            }
        }
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListIndexLong(){
        static bool first = true;
        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_spj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;
        static ReallocatableArray<S32> * n_disp_epj_tmp;
        static ReallocatableArray<S32> * n_disp_spj_tmp;
        const S32 n_thread = Comm::getNumberOfThread();
        if(first){
            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_spj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];
            n_disp_epj_tmp = new ReallocatableArray<S32>[n_thread];
            n_disp_spj_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        const S32 n_ipg = ipg_.size();
        const S32 adr_tc = 0;
        const S32 adr_tree_sp_first = spj_sorted_.size() - tc_glb_.size();
        F64 r_crit_sq = LARGE_FLOAT;
        if (theta_ > 0.0) r_crit_sq = (length_ * length_) / (theta_ * theta_);
        else r_crit_sq = - 1.0; // negative value is used to represent theta_ = 0
        interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
        interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_ep_cum_prev = 0;
            S32 n_sp_cum_prev = 0;
            adr_epj_tmp[ith].clearSize();
            adr_spj_tmp[ith].clearSize();
            adr_ipg_tmp[ith].clearSize();
            n_disp_epj_tmp[ith].clearSize();
            n_disp_spj_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 i=0; i<n_ipg; i++){
                //F64 m_tmp = 0.0;
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
                TargetBox<TSM> target_box;
                target_box.set(ipg_[i]);
                MakeListUsingTreeRecursiveTop
                    <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                     WALK_MODE_NORMAL, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     spj_sorted_, adr_spj_tmp[ith],
                     target_box,
                     r_crit_sq, n_leaf_limit_,
                     adr_tree_sp_first, F64vec(0.0));
                interaction_list_.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                interaction_list_.n_sp_[i] = adr_spj_tmp[ith].size() - n_sp_cum_prev;
                n_ep_cum_prev = adr_epj_tmp[ith].size();
                n_sp_cum_prev = adr_spj_tmp[ith].size();
            }
            n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
            n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
        } // end of OMP

        interaction_list_.n_disp_ep_[0] = 0;
        interaction_list_.n_disp_sp_[0] = 0;
        for(S32 i=0; i<n_ipg; i++){
            interaction_list_.n_disp_ep_[i+1] = interaction_list_.n_disp_ep_[i] + interaction_list_.n_ep_[i];
            interaction_list_.n_disp_sp_[i+1] = interaction_list_.n_disp_sp_[i] + interaction_list_.n_sp_[i];
        }
        interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[n_ipg] );
        interaction_list_.adr_sp_.resizeNoInitialize( interaction_list_.n_disp_sp_[n_ipg] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                const S32 adr_ipg = adr_ipg_tmp[i][j];
                S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                const S32 k_ep_h = n_disp_epj_tmp[i][j];
                const S32 k_ep_e = n_disp_epj_tmp[i][j+1];
                for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                    interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                }
                S32 adr_sp = interaction_list_.n_disp_sp_[adr_ipg];
                const S32 k_sp_h = n_disp_spj_tmp[i][j];
                const S32 k_sp_e = n_disp_spj_tmp[i][j+1];
                for(S32 k=k_sp_h; k<k_sp_e; k++, adr_sp++){
                    interaction_list_.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
                }
            }
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::    
    freeMem(){
        tp_glb_.freeMem();
        tc_loc_.freeMem();
        tc_glb_.freeMem();
        epi_sorted_.freeMem();
        epi_org_.freeMem();
        epj_sorted_.freeMem();
        epj_org_.freeMem();
        spj_sorted_.freeMem();
        spj_org_.freeMem();
        ipg_.freeMem();
        epj_send_.freeMem();
        spj_send_.freeMem();
        force_sorted_.freeMem();
        force_org_.freeMem();
        epjr_sorted_.freeMem();
        epjr_send_.freeMem();
        epjr_recv_.freeMem();
        epjr_recv_1st_buf_.freeMem();
        epjr_recv_2nd_buf_.freeMem();
        const S32 n_thread = Comm::getNumberOfThread();
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].freeMem();
            spj_for_force_[i].freeMem();
            epjr_send_buf_[i].freeMem();
            epjr_send_buf_for_scatter_[i].freeMem();
            epjr_recv_1st_sorted_[i].freeMem();
        }
    }
    
}

#include"tree_for_force_impl_force.hpp"
