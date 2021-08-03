////////////////////////////////////////////////
/// implementaion of methods of TreeForForce ///

#include"tree_walk.hpp"
#include"tree_for_force_impl_exlet.hpp"

namespace ParticleSimulator{

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    setCommInfo(const CommInfo & c){
        comm_info_ = c;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    Tepj * TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    getEpjFromId(const S64 id, const Tepj * epj_tmp){
PS_OMP_CRITICAL      
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>
    ::getUsedMemorySize() const {
        return getMemSizeUsed();
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
        comm_info_.setCommunicator();
        map_id_to_epj_.clear();
        is_initialized_ = true;
        n_glb_tot_ = n_glb_tot;
        theta_ = theta;
        n_leaf_limit_ = n_leaf_limit;
        n_group_limit_ = n_group_limit;
        lev_max_loc_ = lev_max_glb_ = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        const S64 n_proc = comm_info_.getNumberOfProc();
        n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
        ni_ave_ = nj_ave_ = 0;
        //Comm::barrier();
        comm_info_.barrier();
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
        epi_org_.setAllocMode(MemoryAllocMode::Pool);
        epj_send_.setAllocMode(MemoryAllocMode::Pool);
        spj_send_.setAllocMode(MemoryAllocMode::Pool);
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].setAllocMode(MemoryAllocMode::Pool);
            spj_for_force_[i].setAllocMode(MemoryAllocMode::Pool);
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
        comm_info_.barrier();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        req_send_ = new MPI_Request[n_proc];
        req_recv_ = new MPI_Request[n_proc];
        status_   = new MPI_Status[n_proc];
#endif
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(comm_info_.getRank() == 0){
            std::cerr<<"used mem size for tree= "<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tpsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    setParticleLocalTreeImpl(const Tpsys & psys,
                             const bool clear){
        const F64 time_offset = GetWtime();
        const S32 nloc = psys.getNumberOfParticleLocal();
        if(clear){
            n_loc_tot_ = 0;
            inner_boundary_of_local_tree_.init();
        }
        const S32 offset = n_loc_tot_;
        n_loc_tot_ += nloc;
        epj_org_.resizeNoInitialize(n_loc_tot_);
        epi_org_.resizeNoInitialize(n_loc_tot_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        const S32 n_thread = Comm::getNumberOfThread();
        ReallocatableArray<F64ort> inner(n_thread, n_thread, MemoryAllocMode::Pool);
#pragma omp parallel
        {
            const S32 ith = Comm::getThreadNum();
            inner[ith].init();
            if(clear){
#pragma omp for
                for(S32 i=0; i<nloc; i++){
                    epi_org_[i].copyFromFP( psys[i] );
                    epj_org_[i].copyFromFP( psys[i] );
                    inner[ith].merge(psys[i].getPos());
                }
            }
            else{
#pragma omp for
                for(S32 i=0; i<nloc; i++){
                    epi_org_[i+offset].copyFromFP( psys[i] );
                    epj_org_[i+offset].copyFromFP( psys[i] );
                    inner[ith].merge(psys[i].getPos());
                }
            }
        } // end of OMP scope
        for(S32 i=0; i<n_thread; i++){
            inner_boundary_of_local_tree_.merge(inner[i]);
        }
#else //PARTICLE_SIMULATOR_THREAD_PARALLEL
        if(clear){
            for(S32 i=0; i<nloc; i++){
                epi_org_[i].copyFromFP( psys[i] );
                epj_org_[i].copyFromFP( psys[i] );
                inner_boundary_of_local_tree_.merge(psys[i].getPos());
            }
        }
        else{
            for(S32 i=0; i<nloc; i++){
                epi_org_[i+offset].copyFromFP( psys[i] );
                epj_org_[i+offset].copyFromFP( psys[i] );
                inner_boundary_of_local_tree_.merge(psys[i].getPos());
            }
        }
#endif //PARTICLE_SIMULATOR_THREAD_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"nloc="<<nloc<<std::endl;
        std::cout<<"n_loc_tot_="<<n_loc_tot_<<std::endl;
#endif
        time_profile_.set_particle_local_tree += GetWtime() - time_offset;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    setRootCell(const DomainInfo & dinfo){
        const F64 time_offset = GetWtime();
        calcCenterAndLengthOfRootCell(typename TSM::search_type(), dinfo);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        //if(Comm::getRank()==0){
        if(comm_info_.getRank()==0){
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"length_="<<length_<<" center_="<<center_<<std::endl;
            std::cout<<"pos_root_cell="<<getPosRootCell<<std::endl;
        }
#endif
        time_profile_.set_root_cell += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    setRootCell(const F64 l, const F64vec & c){
        const F64 time_offset = GetWtime();
        center_ = c;
        length_ = l;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
        time_profile_.set_root_cell += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    mortonSortLocalTreeOnly(const bool reuse){
        const F64 wtime_offset = GetWtime();
        epi_sorted_.resizeNoInitialize(n_loc_tot_);
        epj_sorted_.resizeNoInitialize(n_loc_tot_);
        adr_org_from_adr_sorted_loc_.resizeNoInitialize(n_loc_tot_);
        tp_glb_.resizeNoInitialize(n_loc_tot_);
        if(!reuse){
            ReallocatableArray<TreeParticle> tp_buf(n_loc_tot_, n_loc_tot_, MemoryAllocMode::Pool);
            morton_key_.initialize( length_ * 0.5, center_);
PS_OMP_PARALLEL_FOR
            for(S32 i=0; i<n_loc_tot_; i++){
                tp_glb_[i].setFromEP(epj_org_[i], i, morton_key_);
            }
#if defined(PARTICLE_SIMULATOR_USE_RADIX_SORT)
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf.getPointer(), 0, n_loc_tot_-1);
#elif defined(PARTICLE_SIMULATOR_USE_STD_SORT)
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_loc_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );
#else
	    MergeSortOmp(tp_glb_, 0, n_loc_tot_,
			 [](const TreeParticle & l, const TreeParticle & r )
			 ->bool{return l.getKey() < r.getKey();} );
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //const S32 adr = tp_loc_[i].adr_ptcl_;
                const S32 adr = tp_glb_[i].adr_ptcl_;
                adr_org_from_adr_sorted_loc_[i] = adr;
            }
            tp_buf.freeMem(1);
        } // end of if() no reuse
PS_OMP_PARALLEL_FOR
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
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    linkCellLocalTreeOnly(){
        const F64 time_offset = GetWtime();
        LinkCell(tc_loc_,  adr_tc_level_partition_loc_, tp_glb_.getPointer(), lev_max_loc_, n_loc_tot_, n_leaf_limit_, morton_key_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tc_loc_.size()="<<tc_loc_.size()<<std::endl;
        std::cout<<"lev_max_loc_="<<lev_max_loc_<<std::endl;
#endif
        time_profile_.link_cell_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcMomentLocalTreeOnly(){
        F64 time_offset = GetWtime();
        calcMomentLocalTreeOnlyImpl(typename TSM::calc_moment_local_tree_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        time_profile_.calc_moment_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree_tot = time_profile_.calc_moment_local_tree + time_profile_.make_local_tree;
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcMomentLocalTreeOnlyImpl(TagCalcMomLongEpjLt){
        CalcMomentLongLocalTree(adr_tc_level_partition_loc_,
                                tc_loc_.getPointer(),
                                epj_sorted_.getPointer(),
                                lev_max_loc_,
                                n_leaf_limit_,
                                morton_key_);
    }
  
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcMomentLocalTreeOnlyImpl(TagCalcMomShortEpjLt){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
  
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcMomentLocalTreeOnlyImpl(TagCalcMomShortEpiLt){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epi_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
    
    ///////////////////////////////
    /// morton sort global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    mortonSortGlobalTreeOnly(const bool reuse){
        F64 time_offset = GetWtime();
        if(!map_id_to_epj_.empty()){
            map_id_to_epj_.clear();
        }
        assert(map_id_to_epj_.empty());
        tp_glb_.resizeNoInitialize(n_glb_tot_);
        const S32 n_ep_tot = epj_org_.size();
        epj_sorted_.resizeNoInitialize(n_ep_tot);
        adr_org_from_adr_sorted_glb_.resizeNoInitialize(n_glb_tot_);
        if(!reuse){
            ReallocatableArray<TreeParticle> tp_buf(n_glb_tot_, n_glb_tot_, MemoryAllocMode::Pool);
#if defined(PARTICLE_SIMULATOR_USE_RADIX_SORT)
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf.getPointer(), 0, n_glb_tot_-1);
#elif defined(PARTICLE_SIMULATOR_USE_STD_SORT)
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );
#else
            const S32 n_add = n_glb_tot_ - n_loc_tot_;
	    ReallocatableArray<TreeParticle> tp_add_buf(n_add, n_add, MemoryAllocMode::Pool);
	    ReallocatableArray<TreeParticle> tp_loc_buf(n_loc_tot_, n_loc_tot_, MemoryAllocMode::Pool);
PS_OMP_PARALLEL_FOR
	    for(S32 i=0; i<n_add; i++){
	      tp_add_buf[i] = tp_glb_[n_loc_tot_+i];
	    }
PS_OMP_PARALLEL_FOR
	    for(S32 i=0; i<n_loc_tot_; i++){
	      tp_loc_buf[i] = tp_glb_[i];
	    }
	    //for(S32 i=1; i<n_loc_tot_; i++){ assert(tp_loc_buf[i-1].getKey() <= tp_loc_buf[i].getKey()); }
	    MergeSortOmp(tp_add_buf, 0, n_add,
			 [](const TreeParticle & l, const TreeParticle & r )
			 ->bool{return l.getKey() < r.getKey();} );
	    //for(S32 i=1; i<n_add; i++){ assert(tp_add_buf[i-1].getKey() <= tp_add_buf[i].getKey()); }

	    MergeSortOmpImpl(tp_loc_buf.getPointer(0), tp_loc_buf.getPointer(n_loc_tot_),
			     tp_add_buf.getPointer(0), tp_add_buf.getPointer(n_add),
			     tp_glb_.getPointer(0),
			     [](const TreeParticle & l, const TreeParticle & r )
			     ->bool{return l.getKey() < r.getKey();} );
	    //for(S32 i=1; i<n_add; i++){ assert(tp_glb_[i-1].getKey() <= tp_glb_[i].getKey()); }
#endif
PS_OMP_PARALLEL_FOR
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
            ReallocatableArray<U32> n_cnt_ep(n_thread, n_thread, MemoryAllocMode::Pool);
            ReallocatableArray<U32> n_cnt_sp(n_thread, n_thread, MemoryAllocMode::Pool);
PS_OMP_PARALLEL
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
PS_OMP_BARRIER
                U32 id_ep_head = 0;
                U32 id_sp_head = 0;
                for(U32 i=0; i<id_thread; i++){
                    id_ep_head += n_cnt_ep[i];
                    id_sp_head += n_cnt_sp[i];
                }
PS_OMP_BARRIER
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
PS_OMP_PARALLEL_FOR
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    linkCellGlobalTreeOnly(){
        const F64 time_offset = GetWtime();
        LinkCell(tc_glb_, adr_tc_level_partition_glb_, tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_, n_leaf_limit_, morton_key_);
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcMomentGlobalTreeOnly(){
        const F64 time_offset = GetWtime();
        calcMomentGlobalTreeOnlyImpl(typename TSM::force_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        time_profile_.calc_moment_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree_tot = time_profile_.calc_moment_global_tree + time_profile_.make_global_tree;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcMomentGlobalTreeOnlyImpl(TagForceLong){
        if(exchange_let_mode_ == EXCHANGE_LET_A2A
           || exchange_let_mode_ == EXCHANGE_LET_P2P_EXACT){
            CalcMomentLongGlobalTree
                (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
                 tp_glb_.getPointer(),     epj_sorted_.getPointer(),
                 spj_sorted_.getPointer(), lev_max_glb_,
                 n_leaf_limit_, morton_key_);
        }
        else if(exchange_let_mode_ == EXCHANGE_LET_P2P_FAST){
            CalcMomentLongGlobalTreeP2P
                (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
                 tp_glb_.getPointer(),     epj_sorted_.getPointer(),
                 spj_sorted_.getPointer(), lev_max_glb_,
                 n_leaf_limit_, comm_table_, pos_root_cell_,
                 adr_org_from_adr_sorted_glb_.getPointer(), morton_key_);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcMomentGlobalTreeOnlyImpl(TagForceShort){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    }
    
    ////////////////////
    /// MAKE IPGROUP ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    makeIPGroup(){
        const F64 time_offset = GetWtime();
        ipg_.clearSize();
        makeIPGroupImpl(typename TSM::force_type());
        n_walk_local_ += ipg_.size();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
        time_profile_.calc_force__make_ipgroup += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_ipg; i++){
            const S32 n = ipg_[i].n_ptcl_;
            const S32 adr = ipg_[i].adr_ptcl_;
            ipg_[i].vertex_in_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
        }
#endif //PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
    template<class Tpj, class Ttc>
    void CopyPjForForceST(const ReallocatableArray<S32> & adr_pj,
                          const ReallocatableArray<Tpj> & pj_sorted,
                          const ReallocatableArray<Ttc> & tc,
                          const S32 n_let_sp,
                          const S32 n_head,
                          const S32 n_tail,
                          ReallocatableArray<Tpj> & pj_for_force){
        pj_for_force.resizeNoInitialize(n_tail);
        for(S32 i=n_head; i<n_tail; i++){
            const auto adr_pj_src = adr_pj[i];
            if(adr_pj_src < n_let_sp){
                pj_for_force[i] = pj_sorted[adr_pj_src];
            }
            else{
                const auto adr_tc_src = adr_pj_src - n_let_sp;
                pj_for_force[i].copyFromMoment(tc[adr_tc_src].mom_);
            }
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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



    inline void CalcNumShift(const F64ort root_domain,
                             const F64ort my_outer_boundary,
                             const F64vec shift,
                             S32 & n_shift){
        F64ort root_domain_tmp = root_domain;
        root_domain_tmp.high_ += shift;
        root_domain_tmp.low_ += shift;
        while(my_outer_boundary.overlapped(root_domain_tmp)){
            root_domain_tmp.high_ += shift;
            root_domain_tmp.low_ += shift;
            n_shift++;
        }        
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcCenterAndLengthOfRootCellImpl(const Tep2 ep[],
                                      const DomainInfo & dinfo){
        F64ort box_loc;
        box_loc.initNegativeVolume();
        if (typeid(TSM) == typeid(SEARCH_MODE_LONG)) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
            {
                F64ort box_loc_tmp;
                box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
                for(S32 ip=0; ip<n_loc_tot_; ip++){
                    box_loc_tmp.merge(ep[ip].getPos());
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                {
                    box_loc.merge(box_loc_tmp);
                }
            }
        } else {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
            {
                F64ort box_loc_tmp;
                box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
                for(S32 ip=0; ip<n_loc_tot_; ip++){
                    box_loc_tmp.merge(ep[ip].getPos(), GetMyRSearch(ep[ip])*1.000001);
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                {
                    box_loc.merge(box_loc_tmp);
                }
            }
        }
        F64vec xlow_loc  = box_loc.low_;
        F64vec xhigh_loc = box_loc.high_;
        F64vec xlow_glb  = comm_info_.getMinValue(xlow_loc);
        F64vec xhigh_glb = comm_info_.getMaxValue(xhigh_loc);
        const F64ort my_outer_boundary(xlow_glb, xhigh_glb);
        const F64ort root_domain = dinfo.getPosRootDomain();
        const F64vec shift = root_domain.high_ - root_domain.low_;
        bool pa[DIMENSION_LIMIT];
        dinfo.getPeriodicAxis(pa);
        S32 num_shift_p[DIMENSION_LIMIT];
        S32 num_shift_n[DIMENSION_LIMIT];
        S32 num_shift_max[DIMENSION_LIMIT];
        for(S32 cid=0; cid<DIMENSION; cid++){
            if (pa[cid]) {
                num_shift_p[cid] = num_shift_n[cid] = 0;
                F64vec shift_tmp(0.0);
                shift_tmp[cid] = shift[cid];
                CalcNumShift(root_domain, my_outer_boundary, shift_tmp, num_shift_p[cid]);
                CalcNumShift(root_domain, my_outer_boundary, -shift_tmp, num_shift_n[cid]);
                num_shift_max[cid] = std::max(num_shift_p[cid], num_shift_n[cid]);
            } else {
                num_shift_max[cid] = 0;
            }
        }
        length_ = 0.0;
        for(S32 cid=0; cid<DIMENSION; cid++){
            if (pa[cid]) {
                F64 length_tmp = (2*num_shift_max[cid]+1)*shift[cid];
                if(length_tmp > length_) length_ = length_tmp;
                center_[cid] = root_domain.getCenter()[cid];
            } else {
                F64 length_tmp = my_outer_boundary.getFullLength()[cid];
                if (length_tmp > length_) length_ = length_tmp;
                center_[cid] = my_outer_boundary.getCenter()[cid];
            }
        }
        length_ *= 1.000001;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tfunc_ep_ep, class Tsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    calcForceDirectParallel(Tfunc_ep_ep pfunc_ep_ep,
                            Tforce force[],
                            const Tsys & system,
                            const DomainInfo & dinfo,
                            const bool clear){
        if(clear){
            for(S32 i=0; i<n_loc_tot_; i++) force[i].clear();
        }
        //S32 n_epj_loc_max = Comm::getMaxValue(n_loc_tot_);
        ReallocatableArray<Tepi> my_epi(n_loc_tot_, n_loc_tot_, MemoryAllocMode::Pool);
        ReallocatableArray<Tepj> my_epj(n_loc_tot_, n_loc_tot_, MemoryAllocMode::Pool);
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = adr_org_from_adr_sorted_loc_[i];
            my_epi[adr] = epi_sorted_[i];
            my_epj[i].copyFromFP(system[i]);
        }
        //const S32 n_proc  = Comm::getNumberOfProc();
        const S32 n_proc  = comm_info_.getNumberOfProc();
        //const S32 my_rank = Comm::getRank();
        const S32 my_rank = comm_info_.getRank();
        ReallocatableArray<F64vec> shift;
        F64ort pos_root_domain = dinfo.getPosRootDomain();
        F64vec domain_size = pos_root_domain.getFullLength();
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        CalcNumberAndShiftOfImageDomain(shift, domain_size, 
                                        pos_root_cell_, pos_root_domain, pa);
        const S32 n_image = shift.size();
        for(S32 i=0; i<n_proc; i++){
            auto rank_send = (my_rank + n_proc + i) % n_proc;
            auto rank_recv = (my_rank + n_proc - i) % n_proc;
            S32 n_recv;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Status stat;
            MPI_Sendrecv(&n_loc_tot_, 1, GetDataType<S32>(), rank_send, 0,
                         &n_recv, 1, GetDataType<S32>(), rank_recv, 0,
                         MPI_COMM_WORLD, &stat);
            ReallocatableArray<Tepj> epj(n_recv*n_image, n_recv*n_image, MemoryAllocMode::Pool);
            MPI_Sendrecv(my_epj.getPointer(), n_loc_tot_, GetDataType<Tepj>(), rank_send, 1,
                         epj.getPointer(), n_recv, GetDataType<Tepj>(), rank_recv, 1,
                         MPI_COMM_WORLD, &stat);
#else
            n_recv = n_loc_tot_;
            ReallocatableArray<Tepj> epj(n_recv*n_image, n_recv*n_image, MemoryAllocMode::Pool);
            for(S32 j=0; j<n_recv; j++){
                epj[j] = my_epj[j];
            }
#endif
            for(S32 j=0; j<n_image; j++){
                for(S32 k=0; k<n_recv; k++){
                    epj[j*n_image+k] = epj[k];
                    const F64vec pos_new = epj[k].getPos() + shift[j];
                    epj[j*n_image+k].pos = pos_new;
                }
            }
            pfunc_ep_ep(my_epi.getPointer(), n_loc_tot_, epj.getPointer(), n_recv*n_image, force);
        }
    }

    
    // return forces in original order
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
        AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain(), pos_root_cell_, pa, comm_info_);
        pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force);
        delete [] epj_tmp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
        AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa, comm_info_);
        pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force_org_.getPointer());
        delete [] epj_tmp;
    }

#if 1
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        return getNeighborListOneParticleImpl(typename TSM::neighbor_search_type(), ptcl, epj);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
    #ifdef PARTICLE_SIMULATOR_CHECK_SEARCH_MODE
    // for test_all_force_mode
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    getNeighborListOneParticleImpl(TagNeighborSearchNo, const Tptcl & ptcl, Tepj * & epj){
        return -1;
        // std::cerr<<"not implemented"<<std::endl;
    }
    #endif
    // New neighbor search mode (this does not require getRSearch)
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    getNeighborListOneParticle(const Tptcl & ptcl, const S32 n_ngb, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        std::vector< std::pair<Tepj, F64> > epj_md_list; // md := metadata
        const S32 adr = 0;
        F64ort vertex;
        SearchNeighborListOneParticleNumber(pos_target,
                                            n_ngb,
                                            tc_glb_.getPointer(),
                                            tp_glb_.getPointer(), adr, 
                                            epj_sorted_, epj_md_list, vertex,
                                            n_leaf_limit_);
        for (const auto & e : epj_md_list) epj_neighbor_[id_thread].push_back( e.first );
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
#else
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        return getNeighborListOneParticleImpl(typename TSM::search_type(), ptcl, epj);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    getOuterBoundaryOfLocalTreeImpl(TagSearchLongSymmetry){
        return tc_loc_[0].geo_.getVertexOut();
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    getOuterBoundaryOfLocalTreeImpl(TagSearchShortSymmetry){
        return tc_loc_[0].geo_.getVertexOut();
    }
    
    template<class TSM>
    F64ort GetIpgBoxForInteractionList(const IPGroup<TSM> & ipg){
        return ipg.vertex_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    makeInteractionListIndexShort(DomainInfo & dinfo){
        const S32 n_thread = Comm::getNumberOfThread();
        std::vector<ReallocatableArray<S32>> adr_epj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        std::vector<ReallocatableArray<S32>> adr_ipg_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        std::vector<ReallocatableArray<S32>> n_disp_epj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        
        const S32 n_ipg = ipg_.size();
        const S32 adr_tc = 0;
        ReallocatableArray<Tspj> spj_dummy;
        ReallocatableArray<S32> adr_spj;
        interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
	//const auto len_peri = dinfo.getPosRootDomain().getFullLength();
	const auto len_peri = dinfo.getLenRootDomain();
PS_OMP_PARALLEL
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_ep_cum_prev = 0;
            adr_epj_tmp[ith].clearSize();
            adr_ipg_tmp[ith].clearSize();
            n_disp_epj_tmp[ith].clearSize();
PS_OMP(omp for schedule(dynamic, 4))
            for(S32 i=0; i<n_ipg; i++){
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                //const F64 r_crit_sq = 999.9;
                TargetBox<TSM> target_box;
                target_box.set(ipg_[i]);
                MakeListUsingTreeRecursiveTop
		    <TSM, TreeCellGlb, TreeParticle, Tepj, TagChopLeafTrue, CALC_DISTANCE_TYPE>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     target_box,
                     n_leaf_limit_,
                     len_peri);
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
    makeInteractionListIndexLong(DomainInfo & dinfo){
        const S32 n_thread = Comm::getNumberOfThread();
        std::vector<ReallocatableArray<S32>> adr_epj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        std::vector<ReallocatableArray<S32>> adr_spj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        std::vector<ReallocatableArray<S32>> adr_ipg_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        std::vector<ReallocatableArray<S32>> n_disp_epj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        std::vector<ReallocatableArray<S32>> n_disp_spj_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
        const S32 n_ipg = ipg_.size();
        const S32 adr_tc = 0;
        const S32 adr_tree_sp_first = spj_sorted_.size() - tc_glb_.size();
        interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
        interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
        const auto len_peri = dinfo.getPosRootDomain().getFullLength();
PS_OMP_PARALLEL
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_ep_cum_prev = 0;
            S32 n_sp_cum_prev = 0;
            adr_epj_tmp[ith].clearSize();
            adr_spj_tmp[ith].clearSize();
            adr_ipg_tmp[ith].clearSize();
            n_disp_epj_tmp[ith].clearSize();
            n_disp_spj_tmp[ith].clearSize();
PS_OMP(omp for schedule(dynamic, 4))
       
            for(S32 i=0; i<n_ipg; i++){
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
                TargetBox<TSM> target_box;
                target_box.set(ipg_[i]);
                MakeListUsingTreeRecursiveTop
		    <TSM, TreeCellGlb, TreeParticle, Tepj, Tspj, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl, CALC_DISTANCE_TYPE>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     spj_sorted_, adr_spj_tmp[ith],
                     target_box,
                     n_leaf_limit_,
                     adr_tree_sp_first, len_peri, theta_);

                interaction_list_.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                interaction_list_.n_sp_[i] = adr_spj_tmp[ith].size() - n_sp_cum_prev;

		/*
		if(Comm::getRank()==0){
		  std::cerr<<"i= "<<i<<" interaction_list_.n_ep_[i]= "<<interaction_list_.n_ep_[i]<<" interaction_list_.n_sp_[i]= "<<interaction_list_.n_sp_[i]<<std::endl;
		  target_box.dump(std::cerr);
		}
		*/
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
PS_OMP_PARALLEL_FOR
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
             class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::
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
