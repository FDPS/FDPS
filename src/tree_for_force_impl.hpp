////////////////////////////////////////////////
/// implementaion of methods of TreeForForce ///

#include"tree_walk.hpp"

namespace ParticleSimulator{
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    Tepj * TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getEpjFromId(const S64 id, const Tepj * epj_tmp){
#if 1
#pragma omp critical
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
#else
        // original
        if(map_id_to_epj_.empty()){
            S64 n_epj = epj_sorted_.size();
            for(S32 i=0; i<n_epj; i++){
                if(GetMSB(tp_glb_[i].adr_ptcl_) == 1) continue;
                Tepj * epj_tmp = epj_sorted_.getPointer(i);
                S64 id_tmp = epj_tmp->getId();
                map_id_to_epj_.insert( std::pair<S64, Tepj*>(id_tmp, epj_tmp) );
            }
        }
#endif
        Tepj * epj = NULL;
        typename MyMap::iterator it = map_id_to_epj_.find(id);
        if( it != map_id_to_epj_.end() ) epj = it->second;
        return epj;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2, class Tep3>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::scatterEP(S32 n_send[],
                S32 n_send_disp[],
                S32 n_recv[],
                S32 n_recv_disp[],
                ReallocatableArray<Tep2> & ep_send,  // send buffer
                ReallocatableArray<Tep2> & ep_recv,  // recv buffer
                ReallocatableArray<Tep2> * ep_send_buf,  // send buffer
                const ReallocatableArray<Tep3> & ep_org, // original
                const DomainInfo & dinfo){
        F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        const F64ort outer_boundary_of_my_tree = tc_loc_[0].mom_.vertex_out_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"outer_boundary_of_my_tree="<<outer_boundary_of_my_tree<<std::endl;
#endif
        wtime_walk_LET_1st_ = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            bool pa[DIMENSION];
            dinfo.getPeriodicAxis(pa);
            ep_send_buf[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc; ib++){
                n_send[ib] = n_recv[ib] = 0;
                shift_image_domain_[ith].clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_domain_[ith].push_back( (F64vec)(0.0));
                }
                else{
                    CalcNumberAndShiftOfImageDomain
                        (shift_image_domain_[ith], 
                         dinfo.getPosRootDomain().getFullLength(),
                         outer_boundary_of_my_tree,
                         dinfo.getPosDomain(ib),
                         pa);
                }
                S32 n_ep_offset = ep_send_buf[ith].size();
                const S32 n_image = shift_image_domain_[ith].size();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
                PARTICLE_SIMULATOR_PRINT_LINE_INFO();
                std::cout<<"n_image="<<n_image<<std::endl;
                std::cout<<"dinfo.getPosDomain(ib)="<<dinfo.getPosDomain(ib)<<std::endl;
#endif

                for(S32 ii=0; ii<n_image; ii++){
                    if(my_rank == ib && ii == 0) continue; // skip self image
                    F64ort pos_domain = dinfo.getPosDomain(ib).shift(shift_image_domain_[ith][ii]);
                    const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#if 0
                    const S32 adr_root_cell = 0;
                    MakeListUsingOuterBoundaryIteration
                        (tc_loc_.getPointer(),  adr_root_cell,
                         ep_org.getPointer(),   ep_send_buf[ith],
                         pos_domain,            n_leaf_limit_,
                         -shift_image_domain_[ith][ii]);
#else
                    if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
                        /*
                          MakeListUsingOuterBoundary
                          (tc_loc_.getPointer(),  adr_tc_tmp,
                          ep_org.getPointer(),   ep_send_buf[ith],
                          pos_domain,            n_leaf_limit_,
                          -shift_image_domain_[ith][ii]);
                        */
                        // add M.I. 2016/03/23
                        MakeListUsingOuterBoundaryWithRootCellCheck
                            (tc_loc_.getPointer(),  adr_tc_tmp,
                             ep_org.getPointer(),   ep_send_buf[ith],
                             pos_domain,            n_leaf_limit_,
                             pos_root_cell_,
                             -shift_image_domain_[ith][ii]);
                    }
                    else{
                        const S32 n_tmp = tc_loc_[0].n_ptcl_;
                        S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                        for(S32 iii=0; iii<n_tmp; iii++, adr_tmp++){
                            const F64vec pos_tmp = ep_org[adr_tmp].getPos();
                            const F64 size_tmp = ep_org[adr_tmp].getRSearch();
                            const F64 dis_sq_tmp = pos_domain.getDistanceMinSQ(pos_tmp);
#ifdef ORIGINAL_SCATTER_MODE
                            if(dis_sq_tmp > size_tmp*size_tmp) continue;
#endif
                            //if( pos_root_cell_.notOverlapped(pos_tmp-shift_image_domain_[ith][ii]) ) continue;  // add M.I. 2016/03/12
                            if( pos_root_cell_.notContained(pos_tmp-shift_image_domain_[ith][ii]) ) continue;  // add M.I. 2016/03/12
                            ep_send_buf[ith].increaseSize();
                            ep_send_buf[ith].back() = ep_org[adr_tmp];
                            const F64vec pos_new = ep_send_buf[ith].back().getPos() - shift_image_domain_[ith][ii];
                            ep_send_buf[ith].back().setPos(pos_new);
                        }
                    }
#endif
                }
                n_send[ib] = ep_send_buf[ith].size() - n_ep_offset;
                id_proc_send_[ith][n_proc_cum++] = ib;  // new
            } // omp for
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                n_send_disp[0] = 0;
                for(S32 i=0; i<n_proc; i++){
                    n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                }
                ep_send.resizeNoInitialize(n_send_disp[n_proc]);
            }
            S32 n_ep_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                const S32 adr_ep_tmp = n_send_disp[id];
                const S32 n_ep_tmp = n_send[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    ep_send[adr_ep_tmp+ip] = ep_send_buf[ith][n_ep_cnt++];
                }
            }
        } // omp parallel scope
        wtime_walk_LET_1st_ = GetWtime() - wtime_walk_LET_1st_;
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        //time_profile_.exchange_LET_tot += GetWtime() - time_offset;
        time_offset = GetWtime();

        Tcomm_scatterEP_tmp_ = GetWtime();
        Comm::allToAll(n_send, 1, n_recv); // TEST
        n_recv_disp[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
        }

        ep_recv.resizeNoInitialize( n_recv_disp[n_proc] );

        Comm::allToAllV(ep_send.getPointer(), n_send, n_send_disp,
                        ep_recv.getPointer(), n_recv, n_recv_disp);

        Tcomm_scatterEP_tmp_ = GetWtime() - Tcomm_scatterEP_tmp_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ep_send.size()="<<ep_send.size()<<std::endl;
        std::cout<<"ep_recv.size()="<<ep_recv.size()<<std::endl;
#endif
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
        //time_profile_.exchange_LET_tot += GetWtime() - time_offset;
    }





    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::getMemSizeUsed() const {
        size_t tmp = 0;
        for(int i=0; i<Comm::getNumberOfThread(); i++){
            tmp = id_ep_send_buf_[i].getMemSize() + id_sp_send_buf_[i].getMemSize()
                + epj_for_force_[i].getMemSize() + spj_for_force_[i].getMemSize()
                + ep_send_buf_for_scatter_[i].getMemSize() + shift_image_domain_[i].getMemSize()
                + epjr_send_buf_[i].getMemSize() + epjr_send_buf_for_scatter_[i].getMemSize() 
                + epjr_recv_1st_sorted_[i].getMemSize() + epj_send_buf_[i].getMemSize() 
                + id_ptcl_send_[i].getMemSize() + shift_image_box_[i].getMemSize()
                + ip_disp_[i].getMemSize() + tp_scatter_[i].getMemSize()
                + tc_recv_1st_[i].getMemSize() + epj_recv_1st_sorted_[i].getMemSize();
        }
        return tmp 
            + tp_buf_.getMemSize() + tp_loc_.getMemSize() + tp_glb_.getMemSize()
            + tc_loc_.getMemSize() + tc_glb_.getMemSize()
            + epi_sorted_.getMemSize() + epi_org_.getMemSize()
            + epj_sorted_.getMemSize() + epj_org_.getMemSize()
            + spj_sorted_.getMemSize() + spj_org_.getMemSize()
            + ipg_.getMemSize()
            + epj_send_.getMemSize() + epj_recv_.getMemSize()
            + spj_send_.getMemSize() + spj_recv_.getMemSize()
            + force_sorted_.getMemSize() + force_org_.getMemSize()
            + epjr_sorted_.getMemSize() + epjr_send_.getMemSize()
            + epjr_recv_.getMemSize() + epjr_recv_1st_buf_.getMemSize()
            + epjr_recv_2nd_buf_.getMemSize()
            + epj_recv_1st_buf_.getMemSize() + epj_recv_2nd_buf_.getMemSize();
    }
    

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    initialize(const U64 n_glb_tot,
               const F64 theta,
               const U32 n_leaf_limit,
               const U32 n_group_limit){
        if(is_initialized_ == true){
            PARTICLE_SIMULATOR_PRINT_ERROR("Do not initialize the tree twice");
            std::cerr<<"SEARCH_MODE: "<<typeid(TSM).name()<<std::endl;
            std::cerr<<"Force: "<<typeid(Tforce).name()<<std::endl;
            std::cerr<<"EPI: "<<typeid(Tepi).name()<<std::endl;
            std::cerr<<"EPJ: "<<typeid(Tepj).name()<<std::endl;
            std::cerr<<"SPJ: "<<typeid(Tspj).name()<<std::endl;
            Abort(-1);
        }
        //comm_table_.initialize();

        map_id_to_epj_.clear();
        is_initialized_ = true;
        n_glb_tot_ = n_glb_tot;
        theta_ = theta;
        n_leaf_limit_ = n_leaf_limit;
        n_group_limit_ = n_group_limit;
        lev_max_loc_ = lev_max_glb_ = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        const S64 n_proc = Comm::getNumberOfProc();
        const S64 np_ave = (n_glb_tot_ / n_proc);
        n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
        wtime_exlet_comm_ = wtime_exlet_a2a_ = wtime_exlet_a2av_ = 0.0;
        wtime_walk_LET_1st_ = wtime_walk_LET_2nd_ = 0.0;
        ni_ave_ = nj_ave_ = 0;
        Comm::barrier();
        //if(Comm::getRank() == 0) std::cerr<<"np_ave="<<np_ave<<std::endl;

        const F64 np_one_dim = pow( ((F64)np_ave)*1.0001, 1.0/DIMENSION) + 4;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        n_surface_for_comm_ = (4*np_one_dim)*6;
#else
        n_surface_for_comm_ = (6*np_one_dim*np_one_dim+8*np_one_dim)*6;
        //n_surface_for_comm_ = (6*np_one_dim*np_one_dim+8*np_one_dim)*2;
#endif
        //epi_org_.reserve( np_ave*3 + 100 );
        epi_org_.reserve( np_ave+(np_ave/3) + 100 );
        epi_sorted_.reserve( epi_org_.capacity() );



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
                //PARTICLE_SIMULATOR_PRINT_ERROR("theta must be >= 0.0");
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

        id_ep_send_buf_   = new ReallocatableArray<S32>[n_thread];
        id_sp_send_buf_   = new ReallocatableArray<S32>[n_thread];
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        id_epj_for_force_ = new ReallocatableArray<S32>[n_thread];
        id_spj_for_force_ = new ReallocatableArray<S32>[n_thread];

        id_proc_send_ = new S32*[n_thread];
        shift_image_domain_ = new ReallocatableArray<F64vec>[n_thread];
        id_ptcl_send_ = new ReallocatableArray<S32>[n_thread];
        ip_disp_ = new ReallocatableArray<S32>[n_thread];
        tp_scatter_ = new ReallocatableArray< TreeParticle >[n_thread];
        tc_recv_1st_ = new ReallocatableArray< TreeCell< Tmomloc > >[n_thread];
        shift_image_box_ = new ReallocatableArray<F64vec>[n_thread];// note
        ep_send_buf_for_scatter_ = new ReallocatableArray<Tepj>[n_thread];
        epj_send_buf_ = new ReallocatableArray<Tepj>[n_thread]; // note
        epj_recv_1st_sorted_ = new ReallocatableArray< Tepj >[n_thread];
        epjr_send_buf_ = new ReallocatableArray<EPJWithR>[n_thread];
        epjr_send_buf_for_scatter_ = new ReallocatableArray<EPJWithR>[n_thread];
        epjr_recv_1st_sorted_ = new ReallocatableArray<EPJWithR>[n_thread];
        epj_neighbor_ = new ReallocatableArray<Tepj>[n_thread];
        n_cell_open_ = new CountT[n_thread];

        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(0)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
        if( typeid(TSM) == typeid(SEARCH_MODE_LONG) || 
            typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) || 
            typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER) ||
            typeid(TSM) == typeid(SEARCH_MODE_LONG_SYMMETRY)){
            if(theta_ > 0.0){
                const S64 n_tmp = epi_org_.capacity() + (S64)2000 * pow((0.5 / theta_), DIMENSION);
                const S64 n_new = std::min( std::min(n_tmp, n_glb_tot_+(S64)(100)), (S64)(20000) );
                //epj_org_.reserve( n_new );
                //epj_sorted_.reserve( n_new );
                epj_org_.reserve( n_new + 10000);
                epj_sorted_.reserve( n_new + 10000);
                spj_org_.reserve( epj_org_.capacity() );
                spj_sorted_.reserve( epj_sorted_.capacity() );
            }
            else if(theta == 0.0) {
                epj_org_.reserve(n_glb_tot_ + 100);
                epj_sorted_.reserve(n_glb_tot_ + 100);
                spj_org_.reserve(n_glb_tot_ + 100);
                spj_sorted_.reserve(n_glb_tot_ + 100);
            }
            for(S32 i=0; i<n_thread; i++){
                id_ep_send_buf_[i].reserve(100);
                id_sp_send_buf_[i].reserve(100);
            }
            spj_send_.reserve(n_proc+1000);
            spj_recv_.reserve(n_proc+1000);
            n_sp_send_ = new S32[n_proc];
            n_sp_recv_ = new S32[n_proc];
            n_sp_send_disp_ = new S32[n_proc+1];
            n_sp_recv_disp_ = new S32[n_proc+1];
            for(S32 i=0; i<n_thread; i++){
                epj_for_force_[i].reserve(10000);
                spj_for_force_[i].reserve(10000);
            }
            n_ep_sp_send_ = new S32[n_proc * 2];
            n_ep_sp_recv_ = new S32[n_proc * 2];
        }
        else{
            // FOR SHORT MODE
            epj_org_.reserve( epi_org_.capacity() + n_surface_for_comm_ );
            epj_sorted_.reserve( epj_org_.capacity() );
            for(S32 i=0; i<n_thread; i++){
                id_ep_send_buf_[i].reserve(n_surface_for_comm_ * 2 / n_thread);
            }
            for(S32 i=0; i<n_thread; i++){
                epj_for_force_[i].reserve(n_surface_for_comm_ * 2 / n_thread);
                spj_for_force_[i].reserve(1);
            }
            n_ep_sp_send_ = n_ep_sp_recv_ = NULL;
        }
        if(typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER)){
            for(S32 ith=0; ith<n_thread; ith++){
                epj_neighbor_[ith].reserve(10);
        }
        }
        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(1)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
        tp_loc_.reserve( epi_org_.capacity() );
        tp_glb_.reserve( epj_org_.capacity() );
        tp_buf_.reserve( epj_org_.capacity() );
        tc_loc_.reserve( tp_loc_.capacity() / n_leaf_limit_ * N_CHILDREN );
        tc_glb_.reserve( tp_glb_.capacity() / n_leaf_limit_ * N_CHILDREN );
        ipg_.reserve( std::min(epi_org_.capacity()/n_group_limit_*4, epi_org_.capacity()) );

        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(2)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
        for(S32 i=0; i<n_thread; i++) id_proc_send_[i] = new S32[n_proc];
        //epj_send_.reserve(n_surface_for_comm_);
        //epj_recv_.reserve(n_surface_for_comm_);
        epj_send_.reserve(100);
        epj_recv_.reserve(100);

        n_ep_send_ = new S32[n_proc];
        n_ep_recv_ = new S32[n_proc];
        n_ep_send_disp_ = new S32[n_proc+1];
        n_ep_recv_disp_ = new S32[n_proc+1];

        force_org_.reserve(epi_org_.capacity());
        force_sorted_.reserve(epi_org_.capacity());

        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(3)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
        // new variables for commnuication of LET
        // for scatterEP

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        req_send_ = new MPI_Request[n_proc];
        req_recv_ = new MPI_Request[n_proc];
        status_   = new MPI_Status[n_proc];
#endif

        if( typeid(TSM) == typeid(SEARCH_MODE_SCATTER) 
            || typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY) 
            || typeid(TSM) == typeid(SEARCH_MODE_GATHER)){
            for(S32 i=0; i<n_thread; i++){
                shift_image_domain_[i].reserve(5*5*5);
            }
        }

        if( typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY)
            || typeid(TSM) == typeid(SEARCH_MODE_GATHER) ){
            n_epj_recv_1st_ = new S32[n_proc];
            n_epj_recv_disp_1st_ = new S32[n_proc+1];
            n_epj_recv_2nd_ = new S32[n_proc];
            n_epj_recv_disp_2nd_ = new S32[n_proc+1];
            id_proc_src_ = new S32[n_proc];
            id_proc_dest_ = new S32[n_proc];
            adr_tc_level_partition_recv_1st_ = new S32*[n_thread];
            for(S32 i=0; i<n_thread; i++){
                id_ptcl_send_[i].reserve(n_surface_for_comm_);
                shift_image_box_[i].reserve(5*5*5);
                ip_disp_[i].reserve(3*3*3);
                tp_scatter_[i].reserve(n_surface_for_comm_);
                tc_recv_1st_[i].reserve(n_surface_for_comm_*4);
                adr_tc_level_partition_recv_1st_[i] = new S32[TREE_LEVEL_LIMIT+2];
            }
        }

        if( typeid(TSM) == typeid(SEARCH_MODE_SCATTER) 
            || typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY) ){
            for(S32 i=0; i<n_thread; i++){
                ep_send_buf_for_scatter_[i].reserve(n_surface_for_comm_);
            }
        }

        // for symmetry
        if( typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY) ){
            // LET 2nd
            epj_recv_1st_buf_.reserve(n_surface_for_comm_);
            epj_recv_2nd_buf_.reserve(n_surface_for_comm_);
            for(S32 i=0; i<n_thread; i++){
                epj_send_buf_[i].reserve(n_surface_for_comm_);
                epj_recv_1st_sorted_[i].reserve(n_surface_for_comm_*4);
            }
        }

        if( typeid(TSM) == typeid(SEARCH_MODE_GATHER) ){
            epjr_sorted_.reserve(epj_org_.capacity());
            epjr_send_.reserve(n_surface_for_comm_);
            epjr_recv_.reserve(n_surface_for_comm_);
            epjr_recv_1st_buf_.reserve(n_surface_for_comm_);
            epjr_recv_2nd_buf_.reserve(n_surface_for_comm_);
            for(S32 i=0; i<n_thread; i++){
                epjr_send_buf_[i].reserve(n_surface_for_comm_);
                epjr_send_buf_for_scatter_[i].reserve(n_surface_for_comm_);
                epjr_recv_1st_sorted_[i].reserve(n_surface_for_comm_);
            }
        }
        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(3)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tpsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setParticleLocalTree(const Tpsys & psys,
                         const bool clear){
        const F64 time_offset = GetWtime();
        const S32 nloc = psys.getNumberOfParticleLocal();
        if(clear){ n_loc_tot_ = 0;}
        //        const S32 offset = 0;
        const S32 offset = n_loc_tot_;
        n_loc_tot_ += nloc;
        epi_org_.resizeNoInitialize(n_loc_tot_);
        epj_org_.resizeNoInitialize(n_loc_tot_);
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
        if(!reuse){
            tp_loc_.resizeNoInitialize(n_loc_tot_);
            tp_buf_.resizeNoInitialize(n_loc_tot_);
            MortonKey::initialize( length_ * 0.5, center_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                tp_loc_[i].setFromEP(epj_org_[i], i);
            }
#ifdef USE_STD_SORT
            std::sort(tp_loc_.getPointer(), tp_loc_.getPointer()+n_loc_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );

#else
            rs_.lsdSort(tp_loc_.getPointer(), tp_buf_.getPointer(), 0, n_loc_tot_-1);
#endif
        } // end of if() no reuse
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = tp_loc_[i].adr_ptcl_;
            epi_sorted_[i] = epi_org_[adr];
            epj_sorted_[i] = epj_org_[adr];
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tp_loc_.size()="<<tp_loc_.size()<<" tp_buf_.size()="<<tp_buf_.size()<<std::endl;
        std::cout<<"epi_sorted_.size()="<<epi_sorted_.size()<<" epj_sorted_.size()="<<epj_sorted_.size()<<std::endl;
#endif
        time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset;
        time_profile_.make_local_tree += GetWtime() - wtime_offset;
    }

    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellLocalTreeOnly(){
        const F64 time_offset = GetWtime();
        //std::cerr<<"start link"<<std::endl;
        LinkCell(tc_loc_,  adr_tc_level_partition_loc_, tp_loc_.getPointer(), 
                 lev_max_loc_, n_loc_tot_, n_leaf_limit_);
        //std::cerr<<"finish link"<<std::endl;
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

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTree(const DomainInfo & dinfo){

        if(typeid(TSM) == typeid(SEARCH_MODE_LONG)
           && dinfo.getBoundaryCondition() != BOUNDARY_CONDITION_OPEN){
            PARTICLE_SIMULATOR_PRINT_ERROR("The forces w/o cutoff can be evaluated only under the open boundary condition");
            Abort(-1);
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(tc_loc_[0].n_ptcl_ <= 0){
            std::cout<<"The number of particles of this process is 0."<<std::endl;
            std::cout<<"tc_loc_[0].n_ptcl_="<<tc_loc_[0].n_ptcl_<<std::endl;
        }
#endif
        exchangeLocalEssentialTreeImpl(typename TSM::search_type(), dinfo);
        time_profile_.exchange_LET_tot = time_profile_.make_LET_1st + time_profile_.exchange_LET_1st + time_profile_.make_LET_2nd + time_profile_.exchange_LET_2nd;
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_send_.size()="<<epj_send_.size()<<" spj_send_.size()="<<spj_send_.size()<<std::endl;
        std::cout<<"epj_recv_.size()="<<epj_recv_.size()<<" spj_recv_.size()="<<spj_recv_.size()<<std::endl;
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLong, const DomainInfo & dinfo){
        F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            id_ep_send_buf_[ith].reserve(1000);
            id_sp_send_buf_[ith].reserve(1000);
            id_ep_send_buf_[ith].resizeNoInitialize(0);
            id_sp_send_buf_[ith].resizeNoInitialize(0);
            F64 r_crit_sq = LARGE_FLOAT;
            if (theta_ > 0.0) r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc; ib++){
                n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                if(my_rank == ib) continue;
                const F64ort pos_target_domain = dinfo.getPosDomain(ib);
                const S32 n_ep_cum = id_ep_send_buf_[ith].size();
                const S32 n_sp_cum = id_sp_send_buf_[ith].size();
                if (theta_ > 0.0) {
                    if(!tc_loc_[0].isLeaf(n_leaf_limit_)){
                            SearchSendParticleLong<TreeCell<Tmomloc>, Tepj>
                                (tc_loc_,               adr_tc_tmp,    epj_sorted_,
                                 id_ep_send_buf_[ith],  id_sp_send_buf_[ith],
                                 pos_target_domain,  r_crit_sq,   n_leaf_limit_);
                    }
                    else{
                        F64vec pos_tmp = tc_loc_[0].mom_.getPos();
                        if( (pos_target_domain.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0) ){
                            const S32 n_tmp = tc_loc_[0].n_ptcl_;
                            S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                            for(S32 ip=0; ip<n_tmp; ip++){
                                id_ep_send_buf_[ith].push_back(adr_tmp++);
                            }
                        }
                        else{
                            id_sp_send_buf_[ith].push_back(0); // set root
                        }
                    }
                } else {
                   // theta_ = 0 case
                   const S32 n_tmp = tc_loc_[0].n_ptcl_;
                   S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                   for(S32 ip=0; ip<n_tmp; ip++){
                       id_ep_send_buf_[ith].push_back(adr_tmp++);
                   }
                }

                n_ep_send_[ib] = id_ep_send_buf_[ith].size() - n_ep_cum;
                n_sp_send_[ib] = id_sp_send_buf_[ith].size() - n_sp_cum;
                id_proc_send_[ith][n_proc_cum++] = ib;
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    const S32 id_ep = id_ep_send_buf_[ith][n_ep_cnt++];
                    epj_send_[adr_ep_tmp++] = epj_sorted_[id_ep];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    const S32 id_sp = id_sp_send_buf_[ith][n_sp_cnt++];
                    spj_send_[adr_sp_tmp++].copyFromMoment(tc_loc_[id_sp].mom_);
                }
            }
        } // omp parallel scope

        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        wtime_exlet_comm_ = GetWtime();

        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send_[2*i] = n_ep_send_[i];
            n_ep_sp_send_[2*i+1] = n_sp_send_[i];
        }
        wtime_exlet_a2a_ = GetWtime();
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_sp_send_, 2, n_ep_sp_recv_);
#else
        Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_); // TEST
#endif //FAST_ALL_TO_ALL_FOR_K


        wtime_exlet_a2a_ = GetWtime() - wtime_exlet_a2a_;

        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_sp_recv_[2*i];
            n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
        }
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );

        wtime_exlet_a2av_ = GetWtime();
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tepj, 2> comm_a2a_epj_2d;
        static CommForAllToAll<Tspj, 2> comm_a2a_spj_2d;
        comm_a2a_epj_2d.executeV(epj_send_, epj_recv_, n_ep_send_, n_ep_recv_);
        comm_a2a_spj_2d.executeV(spj_send_, spj_recv_, n_sp_send_, n_sp_recv_);
#else
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_);
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_,
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_);
#endif //FAST_ALL_TO_ALL_FOR_K
        wtime_exlet_a2av_ = GetWtime() - wtime_exlet_a2av_;
        wtime_exlet_comm_ = GetWtime() - wtime_exlet_comm_;
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongCutoff,
                                   const DomainInfo & dinfo){
        F64 wtime_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        F64 r_crit_sq = LARGE_FLOAT;
        if (theta_ > 0.0) r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
        static bool first = true;
        static ReallocatableArray<Tepj> * epj_send_buf;
        static ReallocatableArray<Tspj> * spj_send_buf;
        static ReallocatableArray<F64ort> outer_boundaries_of_trees;
        if(first){
            const S32 n_thread = Comm::getNumberOfThread();
            epj_send_buf = new ReallocatableArray<Tepj>[n_thread];
            spj_send_buf = new ReallocatableArray<Tspj>[n_thread];
            for(S32 i=0; i<n_thread; i++){
                epj_send_buf[i].reserve(1000);
                spj_send_buf[i].reserve(1000);
            }
            outer_boundaries_of_trees.reserve(n_proc);
            first = false;
        }
        // Compute outer_boundaries_of_trees
        F64ort outer_boundary_of_my_tree;
        outer_boundary_of_my_tree.init();
        for (S32 ip=0; ip < epj_sorted_.size(); ip++) {
            outer_boundary_of_my_tree.merge(epj_sorted_[ip].getPos(),
                                            epj_sorted_[ip].getRSearch());
        }
        Comm::allGather(&outer_boundary_of_my_tree,1,outer_boundaries_of_trees.getPointer());

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            epj_send_buf[ith].clearSize();
            spj_send_buf[ith].clearSize();
            S32 n_epj_cum = 0;
            S32 n_spj_cum = 0;
            const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
            const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
            ReallocatableArray<F64vec> shift_image_domain(5*5*5);
            bool periodic_axis[DIMENSION];
            dinfo.getPeriodicAxis(periodic_axis);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_send_[i] = n_sp_send_[i] = n_ep_recv_[i] = n_sp_recv_[i] = 0;
                shift_image_domain.clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_domain.push_back( F64vec(0.0) );
                }
                else{
                    CalcNumberAndShiftOfImageDomain
                        (shift_image_domain, dinfo.getPosRootDomain().getFullLength(),
                         outer_boundary_of_my_tree, outer_boundaries_of_trees[i], periodic_axis);
                }
                S32 n_image = shift_image_domain.size();
                for(S32 j = 0; j < n_image; j++) {
                    if(my_rank == i && j == 0) continue;
                    const F64ort pos_target_domain = outer_boundaries_of_trees[i].shift(shift_image_domain[j]);
                    if (theta_ > 0.0) {
                        const F64ort cell_box = pos_root_cell_;
                        if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
/*
                            SearchSendParticleLongCutoff<TreeCell<Tmomloc>, Tepj, Tspj>
                                (tc_loc_,       adr_tc_tmp,
                                 epj_sorted_,   epj_send_buf[ith],
                                 spj_send_buf[ith],   cell_box,
                                 pos_target_domain,   r_crit_sq,
                                 r_cut_sq,            n_leaf_limit_,
                                 - shift_image_domain[j]);
*/

                            SearchSendParticleLongCutoffWithRootCellCheck
                                <TreeCell<Tmomloc>, Tepj, Tspj>
                                (tc_loc_,       adr_tc_tmp,
                                 epj_sorted_,   epj_send_buf[ith],
                                 spj_send_buf[ith],   cell_box,
                                 pos_target_domain,   r_crit_sq,
                                 r_cut_sq,            n_leaf_limit_,
                                 pos_root_cell_,
                                 - shift_image_domain[j]);

                        }
                        else{
                            const F64 dis_sq_cut = pos_target_domain.getDistanceMinSQ(cell_box);
                            if (dis_sq_cut <= r_cut_sq) {
                                const F64vec pos_tmp = tc_loc_[0].mom_.getPos();
                                const F64 dis_sq_crit = pos_target_domain.getDistanceMinSQ(pos_tmp);
                                if (dis_sq_crit <= r_crit_sq*4.0) {
                                    const S32 n_tmp = tc_loc_[0].n_ptcl_;
                                    S32 adr_ptcl_tmp = tc_loc_[0].adr_ptcl_;
                                    for (S32 ip=0; ip<n_tmp; ip++) {
                                        const F64vec pos = epj_sorted_[adr_ptcl_tmp].getPos();
                                        const F64vec pos_new = pos - shift_image_domain[j];
                                        if( pos_target_domain.notContained( pos ) ) continue;  
                                        epj_send_buf[ith].push_back(epj_sorted_[adr_ptcl_tmp++]);
                                        epj_send_buf[ith].back().setPos(pos_new);
                                    }
                                }
                                else {
                                    const F64vec pos = tc_loc_[0].mom_.getPos();
                                    const F64vec pos_new = pos - shift_image_domain[j];
                                    if( pos_target_domain.notContained( pos ) ) continue;
                                    spj_send_buf[ith].increaseSize();
                                    spj_send_buf[ith].back().copyFromMoment(tc_loc_[0].mom_);
                                    spj_send_buf[ith].back().setPos(pos_new);
                                }
                            }
                        }
                    }
                    else {
                        // theta_ = 0 case
                        if( !tc_loc_[0].isLeaf(n_leaf_limit_) ) {
                            const F64ort cell_box = pos_root_cell_;
                            // [Note] cell_box MUST be pos_root_cell_ for the correct behavior
                            // of the function makeChildCellBox, which is used in the function
                            // below.
                            SearchSendParticleLongCutoffForZeroTheta
                                <TreeCell<Tmomloc>, Tepj>
                                (tc_loc_,       adr_tc_tmp,
                                 epj_sorted_,   epj_send_buf[ith],
                                 cell_box,      pos_target_domain,
                                 r_cut_sq,      n_leaf_limit_,
                                 pos_root_cell_,
                                 - shift_image_domain[j]);
                        } else {
                            const F64ort cell_box = outer_boundary_of_my_tree;
                            const F64 dis_sq_cut = pos_target_domain.getDistanceMinSQ(cell_box);
                            if (dis_sq_cut <= r_cut_sq){
                                const S32 n_tmp = tc_loc_[0].n_ptcl_;
                                S32 adr_ptcl_tmp = tc_loc_[0].adr_ptcl_;
                                for (S32 ip=0; ip < epj_sorted_.size(); ip++) {
                                    const F64vec pos = epj_sorted_[ip].getPos();
                                    const F64vec pos_new = pos - shift_image_domain[j];
                                    if ( pos_target_domain.notContained( pos ) ) continue;
                                    epj_send_buf[ith].push_back(epj_sorted_[ip]);
                                    epj_send_buf[ith].back().setPos(pos_new);
                                }
                            }
                        }
                    }
                }
                n_ep_send_[i] = epj_send_buf[ith].size() - n_epj_cum;
                n_sp_send_[i] = spj_send_buf[ith].size() - n_spj_cum;
                n_epj_cum = epj_send_buf[ith].size();
                n_spj_cum = spj_send_buf[ith].size();
                id_proc_send_[ith][n_proc_cum++] = i;
            } // end of for loop
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    epj_send_[adr_ep_tmp++] = epj_send_buf[ith][n_ep_cnt++];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    spj_send_[adr_sp_tmp++] = spj_send_buf[ith][n_sp_cnt++];
                }
            }
        } // omp parallel scope

        time_profile_.make_LET_1st += GetWtime() - wtime_offset;
        wtime_offset = GetWtime();

        Comm::allToAll(n_ep_send_, 1, n_ep_recv_); // TEST
        Comm::allToAll(n_sp_send_, 1, n_sp_recv_); // TEST
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_); // TEST
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_, 
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_); // TEST

        time_profile_.exchange_LET_1st += GetWtime() - wtime_offset;

        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
        //time_profile_.exchange_LET_tot += time_profile_.make_LET_1st + time_profile_.exchange_LET_1st;
    }


    ///////////////
    // FOR P^3T
    // original version
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongScatter, const DomainInfo & dinfo){
        F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            id_ep_send_buf_[ith].reserve(1000);
            id_sp_send_buf_[ith].reserve(1000);
            id_ep_send_buf_[ith].resizeNoInitialize(0);
            id_sp_send_buf_[ith].resizeNoInitialize(0);
//#ifdef DEBUG_1023
#ifdef PARTICLE_SIMULATOR_EXCHANGE_LET_ALL
            /// NEW MUST REMOVED
            const F64 r_crit_sq = 1e10;
#else //PARTICLE_SIMULATOR_EXCHANGE_LET_ALL
            F64 r_crit_sq = LARGE_FLOAT;
            if (theta_ > 0.0) r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
#endif //PARTICLE_SIMULATOR_EXCHANGE_LET_ALL
            const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc; ib++){
                n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                if(my_rank == ib) continue;
                const F64ort pos_target_domain = dinfo.getPosDomain(ib);
                const S32 n_ep_cum = id_ep_send_buf_[ith].size();
                const S32 n_sp_cum = id_sp_send_buf_[ith].size();
                if (theta_ > 0.0) {
                    if(!tc_loc_[0].isLeaf(n_leaf_limit_)){
                        SearchSendParticleLongScatter<TreeCell<Tmomloc>, Tepj>
                            (tc_loc_,               adr_tc_tmp,    epj_sorted_,
                             id_ep_send_buf_[ith],  id_sp_send_buf_[ith],
                             pos_target_domain,  r_crit_sq,   n_leaf_limit_);
                    }
                    else if(tc_loc_[0].n_ptcl_==0) continue; // whiout this, if n_loc=0, this node send empty sp!
                    else{
                        F64vec pos_tmp = tc_loc_[0].mom_.getPos();
                        if( (pos_target_domain.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0) ){
                            const S32 n_tmp = tc_loc_[0].n_ptcl_;
                            S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                            for(S32 ip=0; ip<n_tmp; ip++){
                                id_ep_send_buf_[ith].push_back(adr_tmp++);
                            }
                        }
                        else{
                            id_sp_send_buf_[ith].push_back(0); // set root
                        }
                    }
                } else {
                    // theta_ = 0 case
                    const S32 n_tmp = tc_loc_[0].n_ptcl_;
                    S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                    for(S32 ip=0; ip<n_tmp; ip++){
                        id_ep_send_buf_[ith].push_back(adr_tmp++);
                    }
                }
                n_ep_send_[ib] = id_ep_send_buf_[ith].size() - n_ep_cum;
                n_sp_send_[ib] = id_sp_send_buf_[ith].size() - n_sp_cum;
                id_proc_send_[ith][n_proc_cum++] = ib;
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    const S32 id_ep = id_ep_send_buf_[ith][n_ep_cnt++];
                    epj_send_[adr_ep_tmp++] = epj_sorted_[id_ep];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    const S32 id_sp = id_sp_send_buf_[ith][n_sp_cnt++];
                    spj_send_[adr_sp_tmp++].copyFromMoment(tc_loc_[id_sp].mom_);
                }
            }
        } // omp parallel scope

        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send_[2*i] = n_ep_send_[i];
            n_ep_sp_send_[2*i+1] = n_sp_send_[i];
        }

        F64 wtime_offset_tmp = GetWtime();
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_sp_send_, 2, n_ep_sp_recv_);
#else
        Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_); // TEST
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st__a2a_n += GetWtime() - wtime_offset_tmp;

        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_sp_recv_[2*i];
            n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
        }
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );


#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tepj, 2> comm_a2a_epj_2d;
        static CommForAllToAll<Tspj, 2> comm_a2a_spj_2d;

        wtime_offset_tmp = GetWtime();
        comm_a2a_epj_2d.executeV(epj_send_, epj_recv_, n_ep_send_, n_ep_recv_);
        time_profile_.exchange_LET_1st__a2a_ep += GetWtime() - wtime_offset_tmp;

        wtime_offset_tmp = GetWtime();
        comm_a2a_spj_2d.executeV(spj_send_, spj_recv_, n_sp_send_, n_sp_recv_);
        time_profile_.exchange_LET_1st__a2a_sp += GetWtime() - wtime_offset_tmp;
#else
        wtime_offset_tmp = GetWtime();
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_);
        time_profile_.exchange_LET_1st__a2a_ep += GetWtime() - wtime_offset_tmp;

        wtime_offset_tmp = GetWtime();
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_,
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_);
        time_profile_.exchange_LET_1st__a2a_sp += GetWtime() - wtime_offset_tmp;
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;


        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
        //time_profile_.exchange_LET_tot += time_profile_.make_LET_1st + time_profile_.exchange_LET_1st;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongSymmetry, const DomainInfo & dinfo){
        F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        static ReallocatableArray< std::pair<F64ort, F64ort> > tree_top_boundary_pos;
        tree_top_boundary_pos.resizeNoInitialize(n_proc);
        std::pair<F64ort, F64ort> boundary_pos;
        boundary_pos.first  = tc_loc_[0].mom_.getVertexIn();
        boundary_pos.second = tc_loc_[0].mom_.getVertexOut();
        Comm::allGather(&boundary_pos, 1, tree_top_boundary_pos.getPointer());
        tree_top_inner_pos_.resizeNoInitialize(n_proc);
        tree_top_outer_pos_.resizeNoInitialize(n_proc);
        for(S32 i=0; i<n_proc; i++){
            tree_top_inner_pos_[i] = tree_top_boundary_pos[i].first;
            tree_top_outer_pos_[i] = tree_top_boundary_pos[i].second;
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            id_ep_send_buf_[ith].resizeNoInitialize(0);
            id_sp_send_buf_[ith].resizeNoInitialize(0);
            F64 r_crit_sq = LARGE_FLOAT;
            if (theta_ > 0.0) r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc; ib++){
                n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                if(my_rank == ib) continue;
                //const F64ort pos_target_domain = dinfo.getPosDomain(ib);
                const F64ort pos_target_box_in = tree_top_inner_pos_[ib];
                const F64ort pos_target_box_out = tree_top_outer_pos_[ib];
                const S32 n_ep_cum = id_ep_send_buf_[ith].size();
                const S32 n_sp_cum = id_sp_send_buf_[ith].size();
                if (theta_ > 0.0) {
                    if(!tc_loc_[0].isLeaf(n_leaf_limit_)){
                        SearchSendParticleLongSymmetry<TreeCell<Tmomloc>, Tepj>
                            (tc_loc_,               adr_tc_tmp,    epj_sorted_,
                             id_ep_send_buf_[ith],  id_sp_send_buf_[ith],
                             pos_target_box_in,  pos_target_box_out,
                             r_crit_sq,   n_leaf_limit_);
                    }
                    else if(tc_loc_[0].n_ptcl_==0) continue; // whiout this, if n_loc=0, this node send empty sp!
                    else{
                        F64vec pos_tmp = tc_loc_[0].mom_.getPos();
                        //if( (pos_target_domain.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0) ){
                        if( (pos_target_box_in.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0)
                            || (pos_target_box_out.contained(pos_tmp)) ){
                            const S32 n_tmp = tc_loc_[0].n_ptcl_;
                            S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                            for(S32 ip=0; ip<n_tmp; ip++){
                                id_ep_send_buf_[ith].push_back(adr_tmp++);
                            }
                        }
                        else{
                            id_sp_send_buf_[ith].push_back(0); // set root
                        }
                    }
                } else {
                    // theta_ = 0.0 case
                    const S32 n_tmp = tc_loc_[0].n_ptcl_;
                    S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                    for(S32 ip=0; ip<n_tmp; ip++){
                        id_ep_send_buf_[ith].push_back(adr_tmp++);
                    }
                }
                n_ep_send_[ib] = id_ep_send_buf_[ith].size() - n_ep_cum;
                n_sp_send_[ib] = id_sp_send_buf_[ith].size() - n_sp_cum;
                id_proc_send_[ith][n_proc_cum++] = ib;
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    const S32 id_ep = id_ep_send_buf_[ith][n_ep_cnt++];
                    epj_send_[adr_ep_tmp++] = epj_sorted_[id_ep];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    const S32 id_sp = id_sp_send_buf_[ith][n_sp_cnt++];
#if 0
                    if(tc_loc_[id_sp].mom_.mass < 1e-10){
                        std::cerr<<"check aaa"<<std::endl;
                    }
#endif
                    spj_send_[adr_sp_tmp++].copyFromMoment(tc_loc_[id_sp].mom_);
                }
            }
        } // omp parallel scope

        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send_[2*i] = n_ep_send_[i];
            n_ep_sp_send_[2*i+1] = n_sp_send_[i];
        }

        F64 wtime_offset_tmp = GetWtime();
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_sp_send_, 2, n_ep_sp_recv_);
#else
        Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_); // TEST
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st__a2a_n += GetWtime() - wtime_offset_tmp;

        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_sp_recv_[2*i];
            n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
        }
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );


#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tepj, 2> comm_a2a_epj_2d;
        static CommForAllToAll<Tspj, 2> comm_a2a_spj_2d;

        wtime_offset_tmp = GetWtime();
        comm_a2a_epj_2d.executeV(epj_send_, epj_recv_, n_ep_send_, n_ep_recv_);
        time_profile_.exchange_LET_1st__a2a_ep += GetWtime() - wtime_offset_tmp;

        wtime_offset_tmp = GetWtime();
        comm_a2a_spj_2d.executeV(spj_send_, spj_recv_, n_sp_send_, n_sp_recv_);
        time_profile_.exchange_LET_1st__a2a_sp += GetWtime() - wtime_offset_tmp;
#else
        wtime_offset_tmp = GetWtime();
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_);
        time_profile_.exchange_LET_1st__a2a_ep += GetWtime() - wtime_offset_tmp;

        wtime_offset_tmp = GetWtime();
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_,
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_);
        time_profile_.exchange_LET_1st__a2a_sp += GetWtime() - wtime_offset_tmp;
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;


        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
        //time_profile_.exchange_LET_tot += time_profile_.make_LET_1st + time_profile_.exchange_LET_1st;
    }


    


#if 0
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongSymmetry, const DomainInfo & dinfo){
        F64 wtime_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        static ReallocatableArray< std::pair<F64ort, F64ort> > tree_top_boundary_pos;
        tree_top_boundary_pos.resizeNoInitialize(n_proc);
        std::pair<F64ort, F64ort> boundary_pos;
        boundary_pos.first  = tc_loc_[0].mom_.getVertexIn();
        boundary_pos.second = tc_loc_[0].mom_.getVertexOut();
        Comm::allGather(&boundary_pos, 1, tree_top_boundary_pos.getPointer());
        tree_top_inner_pos_.resizeNoInitialize(n_proc);
        tree_top_outer_pos_.resizeNoInitialize(n_proc);
        for(S32 i=0; i<n_proc; i++){
            tree_top_inner_pos_[i] = boundary_pos.first;
            tree_top_outer_pos_[i] = boundary_pos.second;
        }
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            id_ep_send_buf_[ith].resizeNoInitialize(0);
            id_sp_send_buf_[ith].resizeNoInitialize(0);
            ReallocatableArray<Tepj> ep_list_dummy;
            ReallocatableArray<Tspj> sp_list_dummy, sp_list_dummy2;
            S32 n_ep_cum = 0;
            S32 n_sp_cum = 0;
            const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
            //const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc; ib++){
                n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                if(my_rank == ib) continue;
                MakeListUsingTreeRecursiveTop<TSM, MAKE_LIST_MODE_LET, LIST_CONTENT_ID, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                    (tc_loc_, 0, tp_loc_, epj_sorted_, ep_list_dummy,  id_ep_send_buf_[ith], sp_list_dummy, sp_list_dummy2,
                     id_sp_send_buf_[ith], tree_top_inner_pos_[ib], tree_top_outer_pos_[ib], r_crit_sq, n_leaf_limit_);
                n_ep_send_[ib] = id_ep_send_buf_[ith].size() - n_ep_cum;
                n_sp_send_[ib] = id_sp_send_buf_[ith].size() - n_sp_cum;
                n_ep_cum = id_ep_send_buf_[ith].size();
                n_sp_cum = id_sp_send_buf_[ith].size();
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );                
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    const S32 id_ep = id_ep_send_buf_[ith][n_ep_cnt++];
                    epj_send_[adr_ep_tmp++] = epj_sorted_[id_ep];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    const S32 id_sp = id_sp_send_buf_[ith][n_sp_cnt++];
                    spj_send_[adr_sp_tmp++].copyFromMoment(tc_loc_[id_sp].mom_);
                }
            }
        } // omp parallel scope

        time_profile_.make_LET_1st += GetWtime() - wtime_offset;
        wtime_offset = GetWtime();

        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send_[2*i] = n_ep_send_[i];
            n_ep_sp_send_[2*i+1] = n_sp_send_[i];
        }

        F64 wtime_offset_tmp = GetWtime();
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_sp_send_, 2, n_ep_sp_recv_);
#else
        Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_); // TEST
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st__a2a_n += GetWtime() - wtime_offset_tmp;

        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_sp_recv_[2*i];
            n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
        }
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );


#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tepj, 2> comm_a2a_epj_2d;
        static CommForAllToAll<Tspj, 2> comm_a2a_spj_2d;

        wtime_offset_tmp = GetWtime();
        comm_a2a_epj_2d.executeV(epj_send_, epj_recv_, n_ep_send_, n_ep_recv_);
        time_profile_.exchange_LET_1st__a2a_ep += GetWtime() - wtime_offset_tmp;

        wtime_offset_tmp = GetWtime();
        comm_a2a_spj_2d.executeV(spj_send_, spj_recv_, n_sp_send_, n_sp_recv_);
        time_profile_.exchange_LET_1st__a2a_sp += GetWtime() - wtime_offset_tmp;
#else
        wtime_offset_tmp = GetWtime();
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_);
        time_profile_.exchange_LET_1st__a2a_ep += GetWtime() - wtime_offset_tmp;

        wtime_offset_tmp = GetWtime();
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_,
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_);
        time_profile_.exchange_LET_1st__a2a_sp += GetWtime() - wtime_offset_tmp;
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st += GetWtime() - wtime_offset;
        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
        time_profile_.exchange_LET_tot += GetWtime() - wtime_offset;
    }
#endif
    











    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortScatter, const DomainInfo & dinfo){
        scatterEP(n_ep_send_,  n_ep_send_disp_,
                  n_ep_recv_,  n_ep_recv_disp_, 
                  epj_send_,   epj_recv_,
                  ep_send_buf_for_scatter_,
                  epj_sorted_, dinfo);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_send_.size()="<<epj_send_.size()<<" spj_send_.size()="<<spj_send_.size()<<std::endl;
        std::cout<<"epj_recv_.size()="<<epj_recv_.size()<<" spj_recv_.size()="<<spj_recv_.size()<<std::endl;
#endif
        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortGather, const DomainInfo & dinfo){
        exchangeLocalEssentialTreeGatherImpl(typename HasRSearch<Tepj>::type(), dinfo);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeGatherImpl(TagRSearch, const DomainInfo & dinfo){
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"EPJ has RSearch"<<std::endl;
#endif
        exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry(), dinfo);
    }



////////////////
/// gather mode (EPI has no getRSearch())
#if 1
// ver 1.1
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeGatherImpl(TagNoRSearch, const DomainInfo & dinfo){
        const S32 n_proc = Comm::getNumberOfProc();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const S32 my_rank = Comm::getRank();
#endif
        S32 n_proc_src_1st = 0;
        S32 n_proc_dest_1st = 0;
        epjr_sorted_.resizeNoInitialize(epj_sorted_.size());
        for(S32 i=0; i<epj_sorted_.size(); i++){
            epjr_sorted_[i].copyFromEPJ( epj_sorted_[i] );
            epjr_sorted_[i].setRSearch( epi_sorted_[i].getRSearch() );
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        for(S32 i=0; i<5; i++){
            std::cout<<"epjr_sorted_[i].getPos()="<<epjr_sorted_[i].getPos()<<std::endl;
            std::cout<<"epjr_sorted_[i].getRSearch()="<<epjr_sorted_[i].getRSearch()<<std::endl;
        }
#endif
        ////////////
        // 1st STEP (send j particles)
        // GATHER
        // FAST VERSION
        scatterEP(n_ep_send_,  n_ep_send_disp_,
                  n_epj_recv_1st_,     n_epj_recv_disp_1st_,
                  epjr_send_,   epjr_recv_1st_buf_,
                  epjr_send_buf_for_scatter_,
                  epjr_sorted_, dinfo);
        assert(epjr_recv_1st_buf_.size() == n_epj_recv_disp_1st_[n_proc]);
        //TexLET0_ = GetWtime() - TexLET0_;
        F64 time_offset = GetWtime();

        n_let_ep_send_1st_ += (CountT)epjr_send_.size();
        n_let_ep_recv_1st_ += (CountT)epjr_recv_1st_buf_.size();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"STEP 1"<<std::endl;
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_recv_1st_buf_.size()="<<epj_recv_1st_buf_.size()<<" epj_send_size()="<<epj_send_.size()<<std::endl;
#endif
        //TexLET1_ = GetWtime();



        ////////////
        // 2nd step (send j particles)
        // GATHER
        // FAST VERSION

        n_proc_src_1st = n_proc_dest_1st = 0;
        for(S32 i=0; i<n_proc; i++){
            if(n_ep_send_[i] > 0){
                id_proc_dest_[n_proc_dest_1st++] = i;
            }
            if(n_epj_recv_1st_[i] > 0){
                id_proc_src_[n_proc_src_1st++] = i;
            }
        }
        for(S32 i=0; i<n_proc; i++) n_ep_send_[i] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_src_2nd = 0;
            bool pa[DIMENSION];
            dinfo.getPeriodicAxis(pa);
            epjr_send_buf_[ith].clearSize();
            S32 n_ep_send_cum_old = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 ib=0; ib<n_proc_src_1st; ib++){
                const S32 id_proc_tmp = id_proc_src_[ib];
                const S32 i_head = n_epj_recv_disp_1st_[id_proc_tmp];
                const S32 i_tail = n_epj_recv_disp_1st_[id_proc_tmp+1];
                const F64ort pos_root_domain = dinfo.getPosRootDomain();
                const F64vec size_root_domain = pos_root_domain.getFullLength();
                S32vec id_image_new;
                S32vec id_image_old = -9999;
                shift_image_box_[ith].clearSize();
                ip_disp_[ith].clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_box_[ith].pushBackNoCheck( F64vec(0.0) );
                    ip_disp_[ith].pushBackNoCheck(i_head);
                    ip_disp_[ith].pushBackNoCheck(i_tail);
                }
                else{
                    for(S32 ip=i_head; ip<i_tail; ip++){
                        const F64vec pos_target = epjr_recv_1st_buf_[ip].getPos();
                        id_image_new = CalcIDOfImageDomain(pos_root_domain, pos_target, pa);
                        if(id_image_old != id_image_new){
                            ip_disp_[ith].push_back(ip);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
                            shift_image_box_[ith].push_back(F64vec(id_image_new.x*size_root_domain.x, id_image_new.y*size_root_domain.y));
#else
                            shift_image_box_[ith].push_back(F64vec(id_image_new.x*size_root_domain.x, id_image_new.y*size_root_domain.y, id_image_new.z*size_root_domain.z));
#endif
                            id_image_old = id_image_new;
                        }
                    }
                    ip_disp_[ith].push_back(i_tail);
                }
                const S32 n_image = shift_image_box_[ith].size();
                //////////////////
                // version 1.1
                for(S32 ii=0; ii<n_image; ii++){
                    //std::cout<<"check 0"<<std::endl;
                    const F64vec shift = shift_image_box_[ith][ii];
                    //const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
                    const F64ort pos_domain = dinfo.getPosDomain(id_proc_tmp).shift(shift);
                    //std::cout<<"check 1"<<std::endl;
                    ///////////////
                    // MAKE TREE A

                    S32 ncnt = 0;
                    tp_scatter_[ith].resizeNoInitialize(ip_disp_[ith][ii+1] - ip_disp_[ith][ii]);
                    //std::cout<<"ip_disp[ith][ii+1]="<<ip_disp[ith][ii+1]<<std::endl;
                    //std::cout<<"ip_disp[ith][ii]="<<ip_disp[ith][ii]<<std::endl;
                    for(S32 ip=ip_disp_[ith][ii]; ip<ip_disp_[ith][ii+1]; ip++, ncnt++){
                        //tp_scatter[ith][ncnt].setFromEP(epj_recv_1st_buf[ip], ncnt);
                        tp_scatter_[ith][ncnt].setFromEP(epjr_recv_1st_buf_[ip], ip);
                    }
                    std::sort(tp_scatter_[ith].getPointer(), tp_scatter_[ith].getPointer()+ncnt, LessOPKEY());

                    epjr_recv_1st_sorted_[ith].resizeNoInitialize(ncnt);
                    for(S32 ip=0; ip<ncnt; ip++){
                        const S32 adr = tp_scatter_[ith][ip].adr_ptcl_;
                        epjr_recv_1st_sorted_[ith][ip] = epjr_recv_1st_buf_[adr];
                        tp_scatter_[ith][ip].adr_ptcl_ = ip;
                    }
                    // OK

                    const S32 n_leaf_limit_A = 1;
                    S32 lev_max_A = 0;
                    //std::cout<<"check A"<<std::endl;
                    LinkCellST(tc_recv_1st_[ith], adr_tc_level_partition_recv_1st_[ith], 
                               tp_scatter_[ith].getPointer(), lev_max_A, ncnt, n_leaf_limit_A);
                    // OK

                    //std::cout<<"check B"<<std::endl;
                    //exit(1);
                    // probably this function is wrong!

                    CalcMomentST(adr_tc_level_partition_recv_1st_[ith], tc_recv_1st_[ith].getPointer(), 
                                 epjr_recv_1st_sorted_[ith].getPointer(), lev_max_A, n_leaf_limit_A);
                    // OK

                    //std::cout<<"check C"<<std::endl;
                    //exit(1);
                    id_ptcl_send_[ith].clearSize();
                    // OUT
                    //std::cout<<"check D"<<std::endl;
                    //exit(1);

                    MakeLETListByDoubleWalk(tc_recv_1st_[ith].getPointer(), tc_loc_.getPointer(), 
                                            epjr_sorted_.getPointer(),    pos_domain, 
                                            n_leaf_limit_A, n_leaf_limit_,
                                            id_ptcl_send_[ith]);

                    const S32 n_ep_per_image = id_ptcl_send_[ith].size();
                    for(S32 jp=0; jp<n_ep_per_image; jp++){
                        const S32 id_j = id_ptcl_send_[ith][jp];
                        const F64vec pos_j = epjr_sorted_[id_j].getPos();
                        //if( pos_root_cell_.notOverlapped(pos_j-shift) ) continue; // added by M.I. 2016/03/12
                        if( pos_root_cell_.notContained(pos_j-shift) ) continue; // added by M.I. 2016/09/06
                        epjr_send_buf_[ith].push_back(epjr_sorted_[id_j]);
                        epjr_send_buf_[ith].back().setPos(pos_j-shift);
                    }
                } // end of for ii=0 to n_image
                n_ep_send_[id_proc_tmp] = epjr_send_buf_[ith].size() - n_ep_send_cum_old;
                n_ep_send_cum_old = epjr_send_buf_[ith].size();
                if(n_ep_send_[id_proc_tmp] > 0){
                    id_proc_send_[ith][n_proc_src_2nd++] = id_proc_tmp;
                }
            } // end of pragma omp for
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                //std::cout<<"CHECK SINGLE"<<std::endl;
                n_ep_send_disp_[0] = 0;
                for(S32 i=0; i<n_proc; i++){
                    n_ep_send_disp_[i+1] = n_ep_send_disp_[i] + n_ep_send_[i];
                }
                epjr_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            for(S32 ib=0; ib<n_proc_src_2nd; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    epjr_send_[adr_ep_tmp++] = epjr_send_buf_[ith][n_ep_cnt++];
                }
            }
        }// omp parallel scope
        //TexLET1_ = GetWtime() - TexLET1_;
        //wtime_walk_LET_2nd_ = TexLET1_;
        time_profile_.make_LET_2nd += GetWtime() - time_offset;
        time_offset = GetWtime();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"STEP 2"<<std::endl;
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epjr_send_.size()="<<epjr_send_.size()<<std::endl;
#endif


        //TexLET2_ = GetWtime();

        /////////////////////
        // 3rd STEP
        // exchange # of particles

        for(S32 i=0; i<n_proc; i++){
            // probably not needed
            n_epj_recv_2nd_[i] = 0;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        for(S32 i=0; i<n_proc_src_1st; i++){
            S32 id_proc = id_proc_src_[i];
            S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
            MPI_Isend(n_ep_send_+id_proc, 1, GetDataType<S32>(), id_proc, tag, MPI_COMM_WORLD, &req_send_[i]);
        }
        for(S32 i=0; i<n_proc_dest_1st; i++){
            S32 id_proc = id_proc_dest_[i];
            S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
            MPI_Irecv(n_epj_recv_2nd_+id_proc, 1, GetDataType<S32>(), id_proc, tag, MPI_COMM_WORLD, &req_recv_[i]);
        }
        MPI_Waitall(n_proc_src_1st,  req_send_, status_);
        MPI_Waitall(n_proc_dest_1st, req_recv_, status_);
#else
        n_epj_recv_2nd_[0] = n_ep_send_[0];
#endif
        n_epj_recv_disp_2nd_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_epj_recv_disp_2nd_[i+1] = n_epj_recv_disp_2nd_[i] + n_epj_recv_2nd_[i];
        }
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_epj_recv_1st_[i] + n_epj_recv_2nd_[i];
            n_ep_recv_disp_[i] = n_epj_recv_disp_1st_[i] + n_epj_recv_disp_2nd_[i];
        }
        n_ep_recv_disp_[n_proc] = n_epj_recv_disp_1st_[n_proc] + n_epj_recv_disp_2nd_[n_proc];

        //TexLET2_ = GetWtime() - TexLET2_;


#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"STEP 3"<<std::endl;
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_ep_recv_disp_[n_proc]="<<n_ep_recv_disp_[n_proc]<<std::endl;
#endif


        //TexLET3_ = GetWtime();


        /////////////////////
        // 4th STEP
        // exchange EPJ

        epjr_recv_2nd_buf_.resizeNoInitialize( n_epj_recv_disp_2nd_[n_proc] );
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        S32 n_cnt_send = 0;
        S32 n_cnt_recv = 0;
        for(S32 i=0; i<n_proc_src_1st; i++){
            const S32 id_proc = id_proc_src_[i];
            if(n_ep_send_[id_proc] > 0){
                S32 adr = n_ep_send_disp_[id_proc];
                S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
                MPI_Isend(epjr_send_.getPointer()+adr, n_ep_send_[id_proc],
                          GetDataType<EPJWithR>(), id_proc, tag, MPI_COMM_WORLD, &req_send_[n_cnt_send++]);
            }
        }
        for(S32 i=0; i<n_proc_dest_1st; i++){
            const S32 id_proc = id_proc_dest_[i];
            if(n_epj_recv_2nd_[id_proc] > 0){
                S32 adr = n_epj_recv_disp_2nd_[id_proc];
                S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
                MPI_Irecv(epjr_recv_2nd_buf_.getPointer()+adr, n_epj_recv_2nd_[id_proc],
                          GetDataType<EPJWithR>(), id_proc, tag, MPI_COMM_WORLD, &req_recv_[n_cnt_recv++]);
            }
        }
        MPI_Waitall(n_cnt_send, req_send_, status_);
        MPI_Waitall(n_cnt_recv, req_recv_, status_);
#else
        S32 adr_send = n_ep_send_disp_[0];
        S32 adr_recv = n_epj_recv_disp_2nd_[0];
        for(S32 i=0; i<n_ep_send_[0]; i++){
            epjr_recv_2nd_buf_[adr_recv++] = epjr_send_[adr_send++];
        }
#endif
        n_let_ep_send_2nd_ += epj_send_.size();
        n_let_ep_recv_2nd_ += epj_recv_2nd_buf_.size();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"STEP 4"<<std::endl;
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epjr_recv_2nd_buf_.size()="<<epjr_recv_2nd_buf_.size()<<std::endl;
#endif

        /////////////////////
        // 5th STEP
        // set epj_recv_
        epjr_recv_.resizeNoInitialize(n_epj_recv_disp_1st_[n_proc]+n_epj_recv_disp_2nd_[n_proc]);
        epj_recv_.resizeNoInitialize(n_epj_recv_disp_1st_[n_proc]+n_epj_recv_disp_2nd_[n_proc]);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_proc; i++){
            S32 n_cnt = n_ep_recv_disp_[i];
            const S32 j_head_1st = n_epj_recv_disp_1st_[i];
            const S32 j_tail_1st = n_epj_recv_disp_1st_[i+1];
            for(S32 j=j_head_1st; j<j_tail_1st; j++){
                epj_recv_[n_cnt] = epjr_recv_1st_buf_[j].getEPJ();
                epjr_recv_[n_cnt++] = epjr_recv_1st_buf_[j];
            }
            const S32 j_head_2nd = n_epj_recv_disp_2nd_[i];
            const S32 j_tail_2nd = n_epj_recv_disp_2nd_[i+1];
            for(S32 j=j_head_2nd; j<j_tail_2nd; j++){
                epj_recv_[n_cnt] = epjr_recv_2nd_buf_[j].getEPJ();
                epjr_recv_[n_cnt++] = epjr_recv_2nd_buf_[j];
            }
        }

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"STEP 5"<<std::endl;
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif

        //TexLET4_ = GetWtime() - TexLET4_;
        time_profile_.exchange_LET_2nd += GetWtime() - time_offset;
        //time_profile_.exchange_LET_tot += time_profile_.make_LET_2nd + time_profile_.exchange_LET_2nd;
    }

#else // gather search
// ver 1.0
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeGatherImpl(TagNoRSearch, const DomainInfo & dinfo){
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"EPJ dose not have RSearch"<<std::endl;
#endif
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<EPXROnly> ep_x_r_send;
        static ReallocatableArray<EPXROnly> ep_x_r_recv;
        static ReallocatableArray<Tepj> * epj_send_buf; // for 1st communication
        static ReallocatableArray<Tepj> epj_recv_buf; // for 2nd communication
        static S32 * n_ep_recv_1st;
        static S32 * n_ep_recv_disp_1st;
        static S32 * n_ep_recv_2nd;
        static S32 * n_ep_recv_disp_2nd;
        static S32 * id_proc_src;
        static S32 * id_proc_dest;
        static ReallocatableArray<S32> * id_ptcl_send;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const S32 my_rank = Comm::getRank();
        static MPI_Request * req_send;
        static MPI_Request * req_recv;
        static MPI_Status  * status;
#endif
        static bool first = true;
        if(first){
            const S32 n_thread = Comm::getNumberOfThread();
            //ep_x_r_send.reserve(1000);
            //ep_x_r_recv.reserve(1000);
            ep_x_r_send.reserve(n_surface_for_comm_);
            ep_x_r_recv.reserve(n_surface_for_comm_);
            n_ep_recv_1st = new S32[n_proc];
            n_ep_recv_disp_1st = new S32[n_proc+1];
            n_ep_recv_2nd = new S32[n_proc];
            n_ep_recv_disp_2nd = new S32[n_proc+1];
            id_proc_src = new S32[n_proc];
            id_proc_dest = new S32[n_proc];
            id_ptcl_send = new ReallocatableArray<S32>[n_thread];
            epj_send_buf = new ReallocatableArray<Tepj>[n_thread];
            for(S32 i=0; i<n_thread; i++){
                //id_ptcl_send[i].reserve(1000);
                //epj_send_buf[i].reserve(1000);
                id_ptcl_send[i].reserve(n_surface_for_comm_ * 2 / n_thread);
                epj_send_buf[i].reserve(n_surface_for_comm_ * 2 / n_thread);
            }
            //epj_recv_buf.reserve(1000);
            epj_recv_buf.reserve(n_surface_for_comm_);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            req_send = new MPI_Request[n_proc];
            req_recv = new MPI_Request[n_proc];
            status   = new MPI_Status[n_proc];
#endif
            first = false;
        }
        S32 n_proc_src_1st = 0;
        S32 n_proc_dest_1st = 0;
        ////////////
        // 1st STEP
        // GATHER
        // NORMAL VERSION

        //TexLET0_ = GetWtime();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        scatterEPForGather(n_ep_send_,        n_ep_send_disp_,
                           n_ep_recv_1st,     n_ep_recv_disp_1st,
                           ep_x_r_send,       ep_x_r_recv,
                           epi_sorted_, dinfo);
        assert(ep_x_r_recv.size() == n_ep_recv_disp_1st[n_proc]);

        //TexLET0_ = GetWtime() - TexLET0_;

        n_let_ep_send_1st_ += (CountT)ep_x_r_send.size();
        n_let_ep_recv_1st_ += (CountT)ep_x_r_recv.size();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ep_x_r_send.size()="<<ep_x_r_send.size()
                 <<" ep_x_r_recv.size()="<<ep_x_r_recv.size()<<std::endl;
#endif

        //TexLET1_ = GetWtime();

        ////////////
        // 2nd step (send j particles)
        // GATHER
        // NORMAL VERSION
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        n_proc_src_1st = n_proc_dest_1st = 0;
        for(S32 i=0; i<n_proc; i++){
            if(n_ep_send_[i] > 0){
                id_proc_dest[n_proc_dest_1st++] = i;
            }
            if(n_ep_recv_1st[i] > 0){
                id_proc_src[n_proc_src_1st++] = i;
            }
        }
        for(S32 i=0; i<n_proc; i++) n_ep_send_[i] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_src_2nd = 0;
            //ReallocatableArray<F64ort> pos_image_box(5*5*5);
            ReallocatableArray<F64vec> shift_image_box(5*5*5);
            ReallocatableArray<S32> ip_disp(3*3*3);
            bool pa[DIMENSION];
            dinfo.getPeriodicAxis(pa);
            epj_send_buf[ith].clearSize();
            S32 n_ep_send_cum_old = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc_src_1st; ib++){
                const S32 id_proc_tmp = id_proc_src[ib];
                const S32 i_head = n_ep_recv_disp_1st[id_proc_tmp];
                const S32 i_tail = n_ep_recv_disp_1st[id_proc_tmp+1];
                const F64ort pos_root_domain = dinfo.getPosRootDomain();
                const F64vec size_root_domain = pos_root_domain.getFullLength();
                S32vec id_image_new;
                S32vec id_image_old = -9999;
                //pos_image_box.clearSize();
                shift_image_box.clearSize();
                ip_disp.clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_box.push_back( F64vec(0.0) );
                    ip_disp.push_back(i_head);
                    ip_disp.push_back(i_tail);
                    //pos_image_box.push_back( F64ort(9999.9, -9999.9) );
                    /*
                    for(S32 ip=i_head; ip<i_tail; ip++){
                        const F64vec pos_target = epj_recv_1st_buf[ip].getPos();
                        const F64 len_target = epj_recv_1st_buf[ip].getRSearch();
                        pos_image_box.back().merge(pos_target, len_target);
                    }
                    */
                }
                else{
                    for(S32 ip=i_head; ip<i_tail; ip++){
                        const F64vec pos_target = ep_x_r_recv[ip].getPos();
                        id_image_new = CalcIDOfImageDomain(pos_root_domain, pos_target, pa);
                        if(id_image_old != id_image_new){
                            ip_disp.push_back(ip);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
                            shift_image_box.push_back(F64vec(id_image_new.x*size_root_domain.x, id_image_new.y*size_root_domain.y));
#else
                            shift_image_box.push_back(F64vec(id_image_new.x*size_root_domain.x, id_image_new.y*size_root_domain.y, id_image_new.z*size_root_domain.z));
#endif
                            //pos_image_box.push_back( F64ort(9999.9, -9999.9) );
                            id_image_old = id_image_new;
                        }
                        //const F64 len_target = epj_recv_1st_buf[ip].getRSearch();
                        //pos_image_box.back().merge(pos_target, len_target);
                    }
                    ip_disp.push_back(i_tail);
                }
                const S32 n_image = shift_image_box.size();
                for(S32 ii=0; ii<n_image; ii++){
                    const F64vec shift = shift_image_box[ii];
                    const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#if 1
#ifdef UNORDERED_SET 
                    std::unordered_set<S32> id;
                    typedef std::unordered_set<S32>::iterator SetItr;
#else
                    std::set<S32> id;
                    typedef std::set<S32>::iterator SetItr;
#endif
                    for(S32 ip=ip_disp[ii]; ip<ip_disp[ii+1]; ip++){
                        const F64vec pos_target = ep_x_r_recv[ip].getPos();
                        const F64 len_target = ep_x_r_recv[ip].getRSearch();
                        const F64ort particle_box(pos_target, len_target);
                        id_ptcl_send[ith].clearSize();
                        if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
                            MakeListUsingInnerBoundaryForGatherModeNormalMode
                                (tc_loc_.getPointer(),      adr_tc_tmp,
                                 epi_sorted_.getPointer(),  id_ptcl_send[ith],
                                 particle_box,              n_leaf_limit_);
                        }
                        else{
                            S32 n_ptcl_tmp = tc_loc_[0].n_ptcl_;
                            S32 adr_ptcl_tmp = tc_loc_[0].adr_ptcl_;
                            for(S32 ip=0; ip<n_ptcl_tmp; ip++, adr_ptcl_tmp++){
                                const F64vec pos_tmp = epi_sorted_[adr_ptcl_tmp].getPos();
                                //if(particle_box.overlapped(pos_tmp)){
                                if(particle_box.contained(pos_tmp)){
                                    id_ptcl_send[ith].push_back(adr_ptcl_tmp);
                                }
                            }
                        }
                        for(S32 i2=0; i2<id_ptcl_send[ith].size(); i2++) id.insert(id_ptcl_send[ith][i2]);
                    }
                    id_ptcl_send[ith].reserveAtLeast( id.size() ); // TODO: fix here! 
                    id_ptcl_send[ith].clearSize();
                    //for(std::set<S32>::iterator itr = id.begin(); itr != id.end(); ++itr){
                    //for(std::unordered_set<S32>::iterator itr = id.begin(); itr != id.end(); ++itr){
                    for(SetItr itr = id.begin(); itr != id.end(); ++itr){
                        id_ptcl_send[ith].pushBackNoCheck(*itr);
                    }
                    assert( id_ptcl_send[ith].size() == (S32)id.size());
#elif 0
                    // sort version
                    id_ptcl_send[ith].clearSize();
                    for(S32 ip=ip_disp[ii]; ip<ip_disp[ii+1]; ip++){
                        const F64vec pos_target = ep_x_r_recv[ip].getPos();
                        const F64 len_target = ep_x_r_recv[ip].getRSearch();
                        const F64ort particle_box(pos_target, len_target);

                        if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
                            MakeListUsingInnerBoundaryForGatherModeNormalMode
                                (tc_loc_.getPointer(),      adr_tc_tmp,
                                 epi_sorted_.getPointer(),  id_ptcl_send[ith],
                                 particle_box,              n_leaf_limit_);
                        }
                        else{
                            S32 n_ptcl_tmp = tc_loc_[0].n_ptcl_;
                            S32 adr_ptcl_tmp = tc_loc_[0].adr_ptcl_;
                            for(S32 ip=0; ip<n_ptcl_tmp; ip++, adr_ptcl_tmp++){
                                const F64vec pos_tmp = epi_sorted_[adr_ptcl_tmp].getPos();
                                //if(particle_box.overlapped(pos_tmp)){
                                if(particle_box.contained(pos_tmp)){
                                    id_ptcl_send[ith].push_back(adr_ptcl_tmp);
                                }
                            }
                        }
                    }
                    std::sort(id_ptcl_send[ith].getPointer(), id_ptcl_send[ith].getPointer(id_ptcl_send[ith].size()));
                    S32 * adr_end = std::unique(id_ptcl_send[ith].getPointer(), id_ptcl_send[ith].getPointer(id_ptcl_send[ith].size()));
                    id_ptcl_send[ith].resizeNoInitialize( adr_end - id_ptcl_send[ith].getPointer() );
#else
                    id_ptcl_send[ith].clearSize();
                    MakeListUsingInnerBoundaryForGatherModeNormalMode
                        (tc_loc_.getPointer(),     adr_tc_tmp,
                         epi_sorted_.getPointer(), id_ptcl_send[ith],
                         pos_image_box[ii],      n_leaf_limit_);
#endif
                    const S32 n_ep_per_image = id_ptcl_send[ith].size();
                    for(S32 jp=0; jp<n_ep_per_image; jp++){
                        const S32 id_j = id_ptcl_send[ith][jp];
                        const F64vec pos_j = epi_sorted_[id_j].getPos(); // NOTE: not have to shift, because epj_recv have shifted vale.
                        for(S32 ip=ip_disp[ii]; ip<ip_disp[ii+1]; ip++){
                            const F64vec pos_i = ep_x_r_recv[ip].getPos();
                            const F64 len_sq_i = ep_x_r_recv[ip].getRSearch() * ep_x_r_recv[ip].getRSearch();
                            if(pos_j.getDistanceSQ(pos_i) <= len_sq_i){
                                epj_send_buf[ith].push_back( epj_sorted_[id_j] );
                                epj_send_buf[ith].back().setPos(pos_j-shift);
                                break;
                            }
                        }
                    }
                }
                n_ep_send_[id_proc_tmp] = epj_send_buf[ith].size() - n_ep_send_cum_old;
                n_ep_send_cum_old = epj_send_buf[ith].size();
                if(n_ep_send_[id_proc_tmp] > 0){
                    id_proc_send_[ith][n_proc_src_2nd++] = id_proc_tmp;
                }
            }// end of pragma omp for
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = 0;
                for(S32 i=0; i<n_proc; i++){
                    n_ep_send_disp_[i+1] = n_ep_send_disp_[i] + n_ep_send_[i];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            for(S32 ib=0; ib<n_proc_src_2nd; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    epj_send_[adr_ep_tmp++] = epj_send_buf[ith][n_ep_cnt++];
                }
            }
        }// omp parallel scope

        //TexLET1_ = GetWtime() - TexLET1_;
        wtime_walk_LET_2nd_ = TexLET1_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_send_.size()="<<epj_send_.size()<<std::endl;
#endif

        //TexLET2_ = GetWtime();

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL 
        /////////////////////
        // 3rd STEP
        // exchange # of particles
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        for(S32 i=0; i<n_proc; i++){ n_ep_recv_2nd[i] = 0; }
        for(S32 i=0; i<n_proc_src_1st; i++){
            S32 id_proc = id_proc_src[i];
            S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
            MPI_Isend(n_ep_send_+id_proc, 1, GetDataType<S32>(), id_proc, tag, MPI_COMM_WORLD, &req_send[i]);
        }
        for(S32 i=0; i<n_proc_dest_1st; i++){
            S32 id_proc = id_proc_dest[i];
            S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
            MPI_Irecv( n_ep_recv_2nd+id_proc, 1, GetDataType<S32>(), id_proc, tag, MPI_COMM_WORLD, &req_recv[i]);
        }
        MPI_Waitall(n_proc_src_1st,  req_send, status);
        MPI_Waitall(n_proc_dest_1st, req_recv, status);
#else
        n_ep_recv_2nd[0] = n_ep_send_[0];
#endif
        n_ep_recv_disp_2nd[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_2nd[i+1] = n_ep_recv_disp_2nd[i] + n_ep_recv_2nd[i];
        }
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_recv_2nd[i];
            n_ep_recv_disp_[i] = n_ep_recv_disp_2nd[i];
        }
        n_ep_recv_disp_[n_proc] = n_ep_recv_disp_2nd[n_proc];

        //TexLET2_ = GetWtime() - TexLET2_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_ep_recv_disp_[n_proc]="<<n_ep_recv_disp_[n_proc]<<std::endl;
#endif

        //TexLET3_ = GetWtime();

        /////////////////////
        // 4th STEP
        // exchange EPJ
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        epj_recv_buf.resizeNoInitialize( n_ep_recv_disp_2nd[n_proc] );
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL 
        S32 n_cnt_send = 0;
        S32 n_cnt_recv = 0;
        for(S32 i=0; i<n_proc_src_1st; i++){
            const S32 id_proc = id_proc_src[i];
            if(n_ep_send_[id_proc] > 0){
                S32 adr = n_ep_send_disp_[id_proc];
                S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
                MPI_Isend(epj_send_.getPointer()+adr, n_ep_send_[id_proc],
                          GetDataType<Tepj>(), id_proc, tag, MPI_COMM_WORLD, &req_send[n_cnt_send++]);
            }
        }
        for(S32 i=0; i<n_proc_dest_1st; i++){
            const S32 id_proc = id_proc_dest[i];
            if(n_ep_recv_2nd[id_proc] > 0){
                S32 adr = n_ep_recv_disp_2nd[id_proc];
                S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
                MPI_Irecv(epj_recv_buf.getPointer()+adr, n_ep_recv_2nd[id_proc],
                          GetDataType<Tepj>(), id_proc, tag, MPI_COMM_WORLD, &req_recv[n_cnt_recv++]);
            }
        }
        MPI_Waitall(n_cnt_send, req_send, status);
        MPI_Waitall(n_cnt_recv, req_recv, status);
#else
        S32 adr_send = n_ep_send_disp_[0];
        S32 adr_recv = n_ep_recv_disp_2nd[0];
        for(S32 i=0; i<n_ep_send_[0]; i++){
            epj_recv_buf[adr_recv++] = epj_send_[adr_send++];
        }
#endif
        //TexLET3_ = GetWtime() - TexLET3_;

        n_let_ep_send_2nd_ += (CountT)epj_send_.size();
        n_let_ep_recv_2nd_ += (CountT)epj_recv_2nd_buf_.size();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif

        //TexLET4_ = GetWtime();

        /////////////////////
        // 5th STEP
        // set epj_recv_
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_2nd[n_proc] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_proc; i++){
            S32 n_cnt = n_ep_recv_disp_[i];
            const S32 j_head_2nd = n_ep_recv_disp_2nd[i];
            const S32 j_tail_2nd = n_ep_recv_disp_2nd[i+1];
            for(S32 j=j_head_2nd; j<j_tail_2nd; j++){
                epj_recv_[n_cnt++] = epj_recv_buf[j];
            }
        }
        //TexLET4_ = GetWtime() - TexLET4_;
    }
#endif // gather search

    //////////////////////////////
    //  SYMMETRY MODE (OR IF EPJ HAS RSEARCH)
#if 1
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry, const DomainInfo & dinfo){
        //std::cout<<"TagSearchShortSymmetry"<<std::endl;
        const S32 n_proc = Comm::getNumberOfProc();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const S32 my_rank = Comm::getRank();
#endif
        S32 n_proc_src_1st = 0;
        S32 n_proc_dest_1st = 0;
        //TexLET0_ = GetWtime();
        ////////////
        // 1st STEP (send j particles)
        // SYMMETRY 
        // FAST VERSION
        scatterEP(n_ep_send_,  n_ep_send_disp_, 
                  n_epj_recv_1st_,     n_epj_recv_disp_1st_,
                  epj_send_,   epj_recv_1st_buf_,
                  ep_send_buf_for_scatter_,
                  epj_sorted_, dinfo);
        F64 time_offset = GetWtime();
        assert(epj_recv_1st_buf_.size() == n_epj_recv_disp_1st_[n_proc]);

        //TexLET0_ = GetWtime() - TexLET0_;

        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_1st_buf_.size();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_recv_1st_buf_.size()="<<epj_recv_1st_buf_.size()<<" epj_send_size()="<<epj_send_.size()<<std::endl;
#endif

        //TexLET1_ = GetWtime();

        ////////////
        // 2nd step (send j particles)
        // SYMMETRY
        // FAST VERSION
        n_proc_src_1st = n_proc_dest_1st = 0;
        for(S32 i=0; i<n_proc; i++){
            if(n_ep_send_[i] > 0){
                id_proc_dest_[n_proc_dest_1st++] = i;
            }
            if(n_epj_recv_1st_[i] > 0){
                id_proc_src_[n_proc_src_1st++] = i;
            }
        }
        for(S32 i=0; i<n_proc; i++) n_ep_send_[i] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_src_2nd = 0;
            bool pa[DIMENSION];
            dinfo.getPeriodicAxis(pa);
            epj_send_buf_[ith].clearSize();
            S32 n_ep_send_cum_old = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 ib=0; ib<n_proc_src_1st; ib++){
                const S32 id_proc_tmp = id_proc_src_[ib];
                const S32 i_head = n_epj_recv_disp_1st_[id_proc_tmp];
                const S32 i_tail = n_epj_recv_disp_1st_[id_proc_tmp+1];
                const F64ort pos_root_domain = dinfo.getPosRootDomain();
                const F64vec size_root_domain = pos_root_domain.getFullLength();
                S32vec id_image_new;
                S32vec id_image_old = -9999;
                shift_image_box_[ith].clearSize();
                ip_disp_[ith].clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_box_[ith].pushBackNoCheck( F64vec(0.0) );
                    ip_disp_[ith].pushBackNoCheck(i_head);
                    ip_disp_[ith].pushBackNoCheck(i_tail);
                }
                else{
                    for(S32 ip=i_head; ip<i_tail; ip++){
                        const F64vec pos_target = epj_recv_1st_buf_[ip].getPos();
                        id_image_new = CalcIDOfImageDomain(pos_root_domain, pos_target, pa);
                        if(id_image_old != id_image_new){
                            ip_disp_[ith].push_back(ip);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
                            shift_image_box_[ith].push_back(F64vec(id_image_new.x*size_root_domain.x, id_image_new.y*size_root_domain.y));
#else
                            shift_image_box_[ith].push_back(F64vec(id_image_new.x*size_root_domain.x, id_image_new.y*size_root_domain.y, id_image_new.z*size_root_domain.z));
#endif
                            id_image_old = id_image_new;
                        }
                    }
                    ip_disp_[ith].push_back(i_tail);
                }
                const S32 n_image = shift_image_box_[ith].size();
                //////////////////
                // version 1.1
                for(S32 ii=0; ii<n_image; ii++){
                    //std::cout<<"check 0"<<std::endl;
                    const F64vec shift = shift_image_box_[ith][ii];
                    //const S32 adr_tc_tmp = tc_loc_[0].adr_tc_; // commented out 2016/08/31
                    const F64ort pos_domain = dinfo.getPosDomain(id_proc_tmp).shift(shift);
                    //std::cout<<"check 1"<<std::endl;
                    ///////////////
                    // MAKE TREE A

                    S32 ncnt = 0;
                    tp_scatter_[ith].resizeNoInitialize(ip_disp_[ith][ii+1] - ip_disp_[ith][ii]);
                    //std::cout<<"ip_disp[ith][ii+1]="<<ip_disp[ith][ii+1]<<std::endl;
                    //std::cout<<"ip_disp[ith][ii]="<<ip_disp[ith][ii]<<std::endl;
                    for(S32 ip=ip_disp_[ith][ii]; ip<ip_disp_[ith][ii+1]; ip++, ncnt++){
                        //tp_scatter[ith][ncnt].setFromEP(epj_recv_1st_buf[ip], ncnt);
                        tp_scatter_[ith][ncnt].setFromEP(epj_recv_1st_buf_[ip], ip);
                    }
                    std::sort(tp_scatter_[ith].getPointer(), tp_scatter_[ith].getPointer()+ncnt, LessOPKEY());

                    epj_recv_1st_sorted_[ith].resizeNoInitialize(ncnt);
                    for(S32 ip=0; ip<ncnt; ip++){
                        const S32 adr = tp_scatter_[ith][ip].adr_ptcl_;
                        epj_recv_1st_sorted_[ith][ip] = epj_recv_1st_buf_[adr];
                        tp_scatter_[ith][ip].adr_ptcl_ = ip;
                    }
                    // OK

                    const S32 n_leaf_limit_A = 1;
                    S32 lev_max_A = 0;
                    //std::cout<<"check A"<<std::endl;
                    LinkCellST(tc_recv_1st_[ith], adr_tc_level_partition_recv_1st_[ith], 
                               tp_scatter_[ith].getPointer(), lev_max_A, ncnt, n_leaf_limit_A);
                    // OK

                    //std::cout<<"check B"<<std::endl;
                    //exit(1);
                    // probably this function is wrong!
                    /*
                      std::cout<<"1:Comm::getRank()="<<Comm::getRank()<<std::endl;
                      std::cout<<"1:tc_recv_1st[ith].n_ptcl_="<<tc_recv_1st[ith][0].n_ptcl_<<std::endl;
                      std::cout<<"1:tc_recv_1st[ith].mom_.getVertexOut()="<<tc_recv_1st[ith][0].mom_.getVertexOut()<<std::endl;
                      std::cout<<"1:tc_recv_1st[ith].mom_.getVertexIn()="<<tc_recv_1st[ith][0].mom_.getVertexIn()<<std::endl;
                    */
                    CalcMomentST(adr_tc_level_partition_recv_1st_[ith], tc_recv_1st_[ith].getPointer(), 
                                 epj_recv_1st_sorted_[ith].getPointer(), lev_max_A, n_leaf_limit_A);
                    // OK
                    /*
                      std::cout<<"2:Comm::getRank()="<<Comm::getRank()<<std::endl;
                      std::cout<<"2:tc_recv_1st[ith].n_ptcl_="<<tc_recv_1st[ith][0].n_ptcl_<<std::endl;
                      std::cout<<"2:tc_recv_1st[ith].mom_.getVertexOut()="<<tc_recv_1st[ith][0].mom_.getVertexOut()<<std::endl;
                      std::cout<<"2:tc_recv_1st[ith].mom_.getVertexIn()="<<tc_recv_1st[ith][0].mom_.getVertexIn()<<std::endl;
                    */
                    //std::cout<<"check C"<<std::endl;
                    //exit(1);
                    id_ptcl_send_[ith].clearSize();
                    // OUT
                    //std::cout<<"check D"<<std::endl;
                    //exit(1);
                    
                    MakeLETListByDoubleWalk(tc_recv_1st_[ith].getPointer(), tc_loc_.getPointer(), 
                                            epj_sorted_.getPointer(),    pos_domain, 
                                            n_leaf_limit_A, n_leaf_limit_,
                                            id_ptcl_send_[ith]);

                    const S32 n_ep_per_image = id_ptcl_send_[ith].size();
                    for(S32 jp=0; jp<n_ep_per_image; jp++){
                        const S32 id_j = id_ptcl_send_[ith][jp];
                        const F64vec pos_j = epj_sorted_[id_j].getPos();
                        //if( pos_root_cell_.notOverlapped(pos_j-shift) ) continue; // added by M.I. 2016/03/12
                        if( pos_root_cell_.notContained(pos_j-shift) ) continue; // added by M.I. 2016/03/12
                        epj_send_buf_[ith].push_back(epj_sorted_[id_j]);
                        epj_send_buf_[ith].back().setPos(pos_j-shift);
                    }
                } // end of for ii=0 to n_image
                n_ep_send_[id_proc_tmp] = epj_send_buf_[ith].size() - n_ep_send_cum_old;
                n_ep_send_cum_old = epj_send_buf_[ith].size();
                if(n_ep_send_[id_proc_tmp] > 0){
                    id_proc_send_[ith][n_proc_src_2nd++] = id_proc_tmp;
                }
            } // end of pragma omp for
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                //std::cout<<"CHECK SINGLE"<<std::endl;
                n_ep_send_disp_[0] = 0;
                for(S32 i=0; i<n_proc; i++){
                    n_ep_send_disp_[i+1] = n_ep_send_disp_[i] + n_ep_send_[i];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            for(S32 ib=0; ib<n_proc_src_2nd; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    epj_send_[adr_ep_tmp++] = epj_send_buf_[ith][n_ep_cnt++];
                }
            }
        }// omp parallel scope
        //TexLET1_ = GetWtime() - TexLET1_;
        //wtime_walk_LET_2nd_ = TexLET1_;

        time_profile_.make_LET_2nd += GetWtime() - time_offset;
        time_offset = GetWtime();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_send_.size()="<<epj_send_.size()<<std::endl;
#endif

        //TexLET2_ = GetWtime();

        /////////////////////
        // 3rd STEP
        // exchange # of particles
        for(S32 i=0; i<n_proc; i++){
            // probably not needed
            n_epj_recv_2nd_[i] = 0;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        for(S32 i=0; i<n_proc_src_1st; i++){
            S32 id_proc = id_proc_src_[i];
            S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
            MPI_Isend(n_ep_send_+id_proc, 1, GetDataType<S32>(), id_proc, tag, MPI_COMM_WORLD, &req_send_[i]);
        }
        for(S32 i=0; i<n_proc_dest_1st; i++){
            S32 id_proc = id_proc_dest_[i];
            S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
            MPI_Irecv(n_epj_recv_2nd_+id_proc, 1, GetDataType<S32>(), id_proc, tag, MPI_COMM_WORLD, &req_recv_[i]);
        }
        MPI_Waitall(n_proc_src_1st,  req_send_, status_);
        MPI_Waitall(n_proc_dest_1st, req_recv_, status_);
#else
        n_epj_recv_2nd_[0] = n_ep_send_[0];
#endif
        n_epj_recv_disp_2nd_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_epj_recv_disp_2nd_[i+1] = n_epj_recv_disp_2nd_[i] + n_epj_recv_2nd_[i];
        }
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_epj_recv_1st_[i] + n_epj_recv_2nd_[i];
            n_ep_recv_disp_[i] = n_epj_recv_disp_1st_[i] + n_epj_recv_disp_2nd_[i];
        }
        n_ep_recv_disp_[n_proc] = n_epj_recv_disp_1st_[n_proc] + n_epj_recv_disp_2nd_[n_proc];

        //TexLET2_ = GetWtime() - TexLET2_;


#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_ep_recv_disp_[n_proc]="<<n_ep_recv_disp_[n_proc]<<std::endl;
#endif

        //TexLET3_ = GetWtime();

        /////////////////////
        // 4th STEP
        // exchange EPJ
        epj_recv_2nd_buf_.resizeNoInitialize( n_epj_recv_disp_2nd_[n_proc] );
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        S32 n_cnt_send = 0;
        S32 n_cnt_recv = 0;
        for(S32 i=0; i<n_proc_src_1st; i++){
            const S32 id_proc = id_proc_src_[i];
            if(n_ep_send_[id_proc] > 0){
                S32 adr = n_ep_send_disp_[id_proc];
                S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
                MPI_Isend(epj_send_.getPointer()+adr, n_ep_send_[id_proc],
                          GetDataType<Tepj>(), id_proc, tag, MPI_COMM_WORLD, &req_send_[n_cnt_send++]);
            }
        }
        for(S32 i=0; i<n_proc_dest_1st; i++){
            const S32 id_proc = id_proc_dest_[i];
            if(n_epj_recv_2nd_[id_proc] > 0){
                S32 adr = n_epj_recv_disp_2nd_[id_proc];
                S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
                MPI_Irecv(epj_recv_2nd_buf_.getPointer()+adr, n_epj_recv_2nd_[id_proc],
                          GetDataType<Tepj>(), id_proc, tag, MPI_COMM_WORLD, &req_recv_[n_cnt_recv++]);
            }
        }
        MPI_Waitall(n_cnt_send, req_send_, status_);
        MPI_Waitall(n_cnt_recv, req_recv_, status_);
#else
        S32 adr_send = n_ep_send_disp_[0];
        S32 adr_recv = n_epj_recv_disp_2nd_[0];
        for(S32 i=0; i<n_ep_send_[0]; i++){
            epj_recv_2nd_buf_[adr_recv++] = epj_send_[adr_send++];
        }
#endif

        //TexLET3_ = GetWtime() - TexLET3_;

        n_let_ep_send_2nd_ += (CountT)epj_send_.size();
        n_let_ep_recv_2nd_ += (CountT)epj_recv_2nd_buf_.size();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_recv_2nd_buf_.size()="<<epj_recv_2nd_buf_.size()<<std::endl;
#endif

        //TexLET4_ = GetWtime();

        /////////////////////
        // 5th STEP
        // set epj_recv_
        epj_recv_.resizeNoInitialize(n_epj_recv_disp_1st_[n_proc]+n_epj_recv_disp_2nd_[n_proc]);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_proc; i++){
            S32 n_cnt = n_ep_recv_disp_[i];
            const S32 j_head_1st = n_epj_recv_disp_1st_[i];
            const S32 j_tail_1st = n_epj_recv_disp_1st_[i+1];
            for(S32 j=j_head_1st; j<j_tail_1st; j++){
                epj_recv_[n_cnt++] = epj_recv_1st_buf_[j];
            }
            const S32 j_head_2nd = n_epj_recv_disp_2nd_[i];
            const S32 j_tail_2nd = n_epj_recv_disp_2nd_[i+1];
            for(S32 j=j_head_2nd; j<j_tail_2nd; j++){
                epj_recv_[n_cnt++] = epj_recv_2nd_buf_[j];
            }
        }

        //TexLET4_ = GetWtime() - TexLET4_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_epj_recv_disp_1st_[n_proc]="<<n_epj_recv_disp_1st_[n_proc]<<" n_epj_recv_disp_2nd_[n_proc]"<<std::endl;
        std::cout<<"epj_recv_.size()="<<epj_recv_.size()<<std::endl;
#endif

        time_profile_.exchange_LET_2nd += GetWtime() - time_offset;
        //time_profile_.exchange_LET_tot += time_profile_.make_LET_2nd + time_profile_.exchange_LET_2nd;
    }
#else
    // original verion
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry, const DomainInfo & dinfo){
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<Tepj> epj_recv_1st_buf;
        static ReallocatableArray<Tepj> epj_recv_2nd_buf;
        static ReallocatableArray<Tepj> * epj_send_buf; // for 1st communication
        static S32 * n_ep_recv_1st;
        static S32 * n_ep_recv_disp_1st;
        static S32 * n_ep_recv_2nd;
        static S32 * n_ep_recv_disp_2nd;
        static S32 * id_proc_src;
        static S32 * id_proc_dest;
        static ReallocatableArray<S32> * id_ptcl_send;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        static MPI_Request * req_send;
        static MPI_Request * req_recv;
        static MPI_Status  * status;;
        const S32 my_rank = Comm::getRank();
#endif
        static bool first = true;
        if(first){
            const S32 n_thread = Comm::getNumberOfThread();
            //epj_recv_1st_buf.reserve(1000);
            //epj_recv_2nd_buf.reserve(1000);
            epj_recv_1st_buf.reserve(n_surface_for_comm_);
            epj_recv_2nd_buf.reserve(n_surface_for_comm_);
            n_ep_recv_1st = new S32[n_proc];
            n_ep_recv_disp_1st = new S32[n_proc+1];
            n_ep_recv_2nd = new S32[n_proc];
            n_ep_recv_disp_2nd = new S32[n_proc+1];
            id_proc_src = new S32[n_proc];
            id_proc_dest = new S32[n_proc];
            epj_send_buf = new ReallocatableArray<Tepj>[n_thread];
            id_ptcl_send = new ReallocatableArray<S32>[n_thread];
            for(S32 i=0; i<n_thread; i++){
                //id_ptcl_send[i].reserve(1000);
                //epj_send_buf[i].reserve(1000);
                id_ptcl_send[i].reserve(n_surface_for_comm_);
                epj_send_buf[i].reserve(n_surface_for_comm_);
            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            req_send = new MPI_Request[n_proc];
            req_recv = new MPI_Request[n_proc];
            status   = new MPI_Status[n_proc];
#endif
            first = false;
        }
        S32 n_proc_src_1st = 0;
        S32 n_proc_dest_1st = 0;

        ////////////
        // 1st STEP (send j particles)
        // SYMMETRY 
        // FAST VERSION
        scatterEP(n_ep_send_,  n_ep_send_disp_, 
                  n_ep_recv_1st,     n_ep_recv_disp_1st,
                  epj_send_,   epj_recv_1st_buf,
                  ep_send_buf_for_scatter_,
                  epj_sorted_, dinfo);
        assert(epj_recv_1st_buf.size() == n_ep_recv_disp_1st[n_proc]);

        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_1st_buf.size();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_recv_1st_buf.size()="<<epj_recv_1st_buf.size()<<" epj_send_size()="<<epj_send_.size()<<std::endl;
#endif

        ////////////
        // 2nd step (send j particles)
        // SYMMETRY
        // FAST VERSION
        n_proc_src_1st = n_proc_dest_1st = 0;
        for(S32 i=0; i<n_proc; i++){
            if(n_ep_send_[i] > 0){
                id_proc_dest[n_proc_dest_1st++] = i;
            }
            if(n_ep_recv_1st[i] > 0){
                id_proc_src[n_proc_src_1st++] = i;
            }
        }
        for(S32 i=0; i<n_proc; i++) n_ep_send_[i] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_src_2nd = 0;
            //ReallocatableArray<F64ort> pos_image_box(5*5*5);
            ReallocatableArray<F64vec> shift_image_box(5*5*5);
            ReallocatableArray<S32> ip_disp(3*3*3);
            bool pa[DIMENSION];
            dinfo.getPeriodicAxis(pa);
            epj_send_buf[ith].clearSize();
            S32 n_ep_send_cum_old = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc_src_1st; ib++){
                const S32 id_proc_tmp = id_proc_src[ib];
                const S32 i_head = n_ep_recv_disp_1st[id_proc_tmp];
                const S32 i_tail = n_ep_recv_disp_1st[id_proc_tmp+1];
                const F64ort pos_root_domain = dinfo.getPosRootDomain();
                const F64vec size_root_domain = pos_root_domain.getFullLength();
                S32vec id_image_new;
                S32vec id_image_old = -9999;
                //pos_image_box.clearSize();
                shift_image_box.clearSize();
                ip_disp.clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_box.push_back( F64vec(0.0) );
                    //pos_image_box.push_back( F64ort(9999.9, -9999.9) );
                    ip_disp.push_back(i_head);
                    ip_disp.push_back(i_tail);
/*
                    for(S32 ip=i_head; ip<i_tail; ip++){
                        const F64vec pos_target = epj_recv_1st_buf[ip].getPos();
                        const F64 len_target = epj_recv_1st_buf[ip].getRSearch();
                        pos_image_box.back().merge(pos_target, len_target);
                    }
*/
                }
                else{
                    for(S32 ip=i_head; ip<i_tail; ip++){
                        const F64vec pos_target = epj_recv_1st_buf[ip].getPos();
                        id_image_new = CalcIDOfImageDomain(pos_root_domain, pos_target, pa);
                        if(id_image_old != id_image_new){
                            ip_disp.push_back(ip);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
                            shift_image_box.push_back(F64vec(id_image_new.x*size_root_domain.x, id_image_new.y*size_root_domain.y));
#else
                            shift_image_box.push_back(F64vec(id_image_new.x*size_root_domain.x, id_image_new.y*size_root_domain.y, id_image_new.z*size_root_domain.z));
#endif
                            //pos_image_box.push_back( F64ort(9999.9, -9999.9) );
                            id_image_old = id_image_new;
                        }
                        //const F64 len_target = epj_recv_1st_buf[ip].getRSearch();
                        //pos_image_box.back().merge(pos_target, len_target);
                    }
                    ip_disp.push_back(i_tail);
                }

                const S32 n_image = shift_image_box.size();
                for(S32 ii=0; ii<n_image; ii++){
                    const F64vec shift = shift_image_box[ii];
                    const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
                    const F64ort pos_domain = dinfo.getPosDomain(id_proc_tmp).shift(shift);
#if 1
                    // use std::set
#ifdef UNORDERED_SET 
                    std::unordered_set<S32> id;
                    typedef std::unordered_set<S32>::iterator SetItr;
#else
                    std::set<S32> id;
                    typedef std::set<S32>::iterator SetItr;
#endif
                    for(S32 ip=ip_disp[ii]; ip<ip_disp[ii+1]; ip++){
                        const F64vec pos_target = epj_recv_1st_buf[ip].getPos();
                        const F64 len_target = epj_recv_1st_buf[ip].getRSearch();
                        const F64ort particle_box(pos_target, len_target);
                        id_ptcl_send[ith].clearSize();

                        if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
                            MakeListUsingInnerBoundaryForSymmetryExclusive
                                (tc_loc_.getPointer(),     adr_tc_tmp,
                                 epj_sorted_.getPointer(), id_ptcl_send[ith],
                                 particle_box,                   pos_domain,
                                 n_leaf_limit_);
                        }
                        else{
                            const S32 n_tmp = tc_loc_[0].n_ptcl_;
                            //id_ptcl_send[ith].reserve( id_ptcl_send[ith].size() + n_tmp );
                            id_ptcl_send[ith].reserveEmptyAreaAtLeast( n_tmp );
                            S32 adr_ptcl_tmp = tc_loc_[0].adr_ptcl_;
                            for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                                // NOTE: need to be concistent with MakeListUsingOuterBoundary()
                                const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                                const F64 len_sq = epj_sorted_[adr_ptcl_tmp].getRSearch() * epj_sorted_[adr_ptcl_tmp].getRSearch();
                                const F64 dis_sq0 = particle_box.getDistanceMinSQ(pos_tmp);
                                const F64 dis_sq1 = pos_domain.getDistanceMinSQ(pos_tmp);
                                if(dis_sq0 <= len_sq && dis_sq1 > len_sq){
                                    id_ptcl_send[ith].pushBackNoCheck(adr_ptcl_tmp);
                                }
                            }
                        }
                        for(S32 i2=0; i2<id_ptcl_send[ith].size(); i2++){
                            id.insert(id_ptcl_send[ith][i2]);
                        }
                    }
                    id_ptcl_send[ith].reserveAtLeast( id.size() );
                    id_ptcl_send[ith].clearSize();
                    //for(std::set<S32>::iterator itr = id.begin(); itr != id.end(); ++itr){
                    //for(std::unordered_set<S32>::iterator itr = id.begin(); itr != id.end(); ++itr){
                    for(SetItr itr = id.begin(); itr != id.end(); ++itr){
                        id_ptcl_send[ith].pushBackNoCheck(*itr);
                    }
                    assert( id_ptcl_send[ith].size() == (S32)id.size());
#elif 0
                    // sort version
                    id_ptcl_send[ith].clearSize();
                    for(S32 ip=ip_disp[ii]; ip<ip_disp[ii+1]; ip++){
                        const F64vec pos_target = epj_recv_1st_buf[ip].getPos();
                        const F64 len_target = epj_recv_1st_buf[ip].getRSearch();
                        const F64ort particle_box(pos_target, len_target);

                        if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
                            MakeListUsingInnerBoundaryForSymmetryExclusive
                                (tc_loc_.getPointer(),     adr_tc_tmp,
                                 epj_sorted_.getPointer(), id_ptcl_send[ith],
                                 particle_box,                   pos_domain,
                                 n_leaf_limit_);
                        }
                        else{
                            const S32 n_tmp = tc_loc_[0].n_ptcl_;
                            //id_ptcl_send[ith].reserve( id_ptcl_send[ith].size() + n_tmp );
                            id_ptcl_send[ith].reserveEmptyAreaAtLeast( n_tmp );
                            S32 adr_ptcl_tmp = tc_loc_[0].adr_ptcl_;
                            for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                                // NOTE: need to be concistent with MakeListUsingOuterBoundary()
                                const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                                const F64 len_sq = epj_sorted_[adr_ptcl_tmp].getRSearch() * epj_sorted_[adr_ptcl_tmp].getRSearch();
                                const F64 dis_sq0 = particle_box.getDistanceMinSQ(pos_tmp);
                                const F64 dis_sq1 = pos_domain.getDistanceMinSQ(pos_tmp);
                                if(dis_sq0 <= len_sq && dis_sq1 > len_sq){
                                    id_ptcl_send[ith].pushBackNoCheck(adr_ptcl_tmp);
                                }
                            }
                        }
                    }
                    std::sort(id_ptcl_send[ith].getPointer(), id_ptcl_send[ith].getPointer(id_ptcl_send[ith].size()));
                    S32 * adr_end = std::unique(id_ptcl_send[ith].getPointer(), id_ptcl_send[ith].getPointer(id_ptcl_send[ith].size()));
                    id_ptcl_send[ith].resizeNoInitialize( adr_end - id_ptcl_send[ith].getPointer() );
#else
                    id_ptcl_send[ith].clearSize();
                    MakeListUsingInnerBoundaryForSymmetryExclusive
                        (tc_loc_.getPointer(),     adr_tc_tmp,
                         epj_sorted_.getPointer(), id_ptcl_send[ith],
                         pos_image_box[ii],      pos_domain,
                         n_leaf_limit_);
#endif
                    const S32 n_ep_per_image = id_ptcl_send[ith].size();
                    for(S32 jp=0; jp<n_ep_per_image; jp++){
                        const S32 id_j = id_ptcl_send[ith][jp];
                        const F64vec pos_j = epj_sorted_[id_j].getPos(); // NOTE: not have to shift, because epj_recv have shifted vale.
                        for(S32 ip=ip_disp[ii]; ip<ip_disp[ii+1]; ip++){
                            const F64vec pos_i = epj_recv_1st_buf[ip].getPos();
                            const F64 len_sq_i = epj_recv_1st_buf[ip].getRSearch() * epj_recv_1st_buf[ip].getRSearch();
                            if(pos_j.getDistanceSQ(pos_i) <= len_sq_i){
                                epj_send_buf[ith].push_back(epj_sorted_[id_j]);
                                epj_send_buf[ith].back().setPos(pos_j-shift);
                                break;
                            }
                        }
                    }
                }
                n_ep_send_[id_proc_tmp] = epj_send_buf[ith].size() - n_ep_send_cum_old;
                n_ep_send_cum_old = epj_send_buf[ith].size();
                if(n_ep_send_[id_proc_tmp] > 0){
                    id_proc_send_[ith][n_proc_src_2nd++] = id_proc_tmp;
                }
            } // end of pragma omp for
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = 0;
                for(S32 i=0; i<n_proc; i++){
                    n_ep_send_disp_[i+1] = n_ep_send_disp_[i] + n_ep_send_[i];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            for(S32 ib=0; ib<n_proc_src_2nd; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    epj_send_[adr_ep_tmp++] = epj_send_buf[ith][n_ep_cnt++];
                }
            }
        }// omp parallel scope

        //TexLET1_ = GetWtime() - TexLET1_;

        //wtime_walk_LET_2nd_ = TexLET1_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_send_.size()="<<epj_send_.size()<<std::endl;
#endif

        //TexLET2_ = GetWtime();

        /////////////////////
        // 3rd STEP
        // exchange # of particles
        for(S32 i=0; i<n_proc; i++){
            // probably not needed
            n_ep_recv_2nd[i] = 0;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        for(S32 i=0; i<n_proc_src_1st; i++){
            S32 id_proc = id_proc_src[i];
            S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
            MPI_Isend(n_ep_send_+id_proc, 1, GetDataType<S32>(), id_proc, tag,MPI_COMM_WORLD, &req_send[i]);
        }
        for(S32 i=0; i<n_proc_dest_1st; i++){
            S32 id_proc = id_proc_dest[i];
            S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
            MPI_Irecv(n_ep_recv_2nd+id_proc, 1, GetDataType<S32>(), id_proc, tag,MPI_COMM_WORLD, &req_recv[i]);
        }
        MPI_Waitall(n_proc_src_1st, req_send, status);
        MPI_Waitall(n_proc_dest_1st, req_recv, status);
#else
        n_ep_recv_2nd[0] = n_ep_send_[0];
#endif
        n_ep_recv_disp_2nd[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_2nd[i+1] = n_ep_recv_disp_2nd[i] + n_ep_recv_2nd[i];
        }
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_recv_1st[i] + n_ep_recv_2nd[i];
            n_ep_recv_disp_[i] = n_ep_recv_disp_1st[i] + n_ep_recv_disp_2nd[i];
        }
        n_ep_recv_disp_[n_proc] = n_ep_recv_disp_1st[n_proc] + n_ep_recv_disp_2nd[n_proc];

        //TexLET2_ = GetWtime() - TexLET2_;


#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_ep_recv_disp_[n_proc]="<<n_ep_recv_disp_[n_proc]<<std::endl;
#endif

        //TexLET3_ = GetWtime();

        /////////////////////
        // 4th STEP
        // exchange EPJ
        epj_recv_2nd_buf.resizeNoInitialize( n_ep_recv_disp_2nd[n_proc] );
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        S32 n_cnt_send = 0;
        S32 n_cnt_recv = 0;
        for(S32 i=0; i<n_proc_src_1st; i++){
            const S32 id_proc = id_proc_src[i];
            if(n_ep_send_[id_proc] > 0){
                S32 adr = n_ep_send_disp_[id_proc];
                S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
                MPI_Isend(epj_send_.getPointer()+adr, n_ep_send_[id_proc],
                          GetDataType<Tepj>(), id_proc, tag, MPI_COMM_WORLD, &req_send[n_cnt_send++]);
            }
        }
        for(S32 i=0; i<n_proc_dest_1st; i++){
            const S32 id_proc = id_proc_dest[i];
            if(n_ep_recv_2nd[id_proc] > 0){
                S32 adr = n_ep_recv_disp_2nd[id_proc];
                S32 tag = (my_rank < id_proc) ? my_rank : id_proc;
                MPI_Irecv(epj_recv_2nd_buf.getPointer()+adr, n_ep_recv_2nd[id_proc],
                          GetDataType<Tepj>(), id_proc, tag, MPI_COMM_WORLD, &req_recv[n_cnt_recv++]);
            }
        }
        MPI_Waitall(n_cnt_send, req_send, status);
        MPI_Waitall(n_cnt_recv, req_recv, status);
#else
        S32 adr_send = n_ep_send_disp_[0];
        S32 adr_recv = n_ep_recv_disp_2nd[0];
        for(S32 i=0; i<n_ep_send_[0]; i++){
            epj_recv_2nd_buf[adr_recv++] = epj_send_[adr_send++];
        }
#endif

        //TexLET3_ = GetWtime() - TexLET3_;

        n_let_ep_send_2nd_ += (CountT)epj_send_.size();
        n_let_ep_recv_2nd_ += (CountT)epj_recv_2nd_buf.size();


#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_recv_2nd_buf.size()="<<epj_recv_2nd_buf.size()<<std::endl;
#endif

        //TexLET4_ = GetWtime();

        /////////////////////
        // 5th STEP
        // set epj_recv_
        epj_recv_.resizeNoInitialize(n_ep_recv_disp_1st[n_proc]+n_ep_recv_disp_2nd[n_proc]);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_proc; i++){
            S32 n_cnt = n_ep_recv_disp_[i];
            const S32 j_head_1st = n_ep_recv_disp_1st[i];
            const S32 j_tail_1st = n_ep_recv_disp_1st[i+1];
            for(S32 j=j_head_1st; j<j_tail_1st; j++){
                epj_recv_[n_cnt++] = epj_recv_1st_buf[j];
            }
            const S32 j_head_2nd = n_ep_recv_disp_2nd[i];
            const S32 j_tail_2nd = n_ep_recv_disp_2nd[i+1];
            for(S32 j=j_head_2nd; j<j_tail_2nd; j++){
                epj_recv_[n_cnt++] = epj_recv_2nd_buf[j];
            }
        }
        //TexLET4_ = GetWtime() - TexLET4_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_ep_recv_disp_1st[n_proc]="<<n_ep_recv_disp_1st[n_proc]<<" n_ep_recv_disp_2nd[n_proc]"<<std::endl;
        std::cout<<"epj_recv_.size()="<<epj_recv_.size()<<std::endl;
#endif
    }

#endif // END OF EXCHANGE LET FOR SYMMETRY

    
    //////////////////////////////
    /// SET LET TO GLOBAL TREE ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTree(){
        F64 time_offset = GetWtime();
        setLocalEssentialTreeToGlobalTreeImpl(typename TSM::force_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_loc_tot_="<<n_loc_tot_<<" n_glb_tot_="<<n_glb_tot_<<std::endl;
        std::cout<<"epj_org_.size()="<<epj_org_.size()<<" spj_org_.size()="<<spj_org_.size()<<" tp_glb_.size()="<<tp_glb_.size()<<std::endl;
#endif
        time_profile_.set_particle_global_tree += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTreeImpl(TagForceShort){
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 offset = this->n_loc_tot_;
        const S32 n_ep_add = this->n_ep_recv_disp_[n_proc];
        this->n_glb_tot_ = this->n_loc_tot_ + n_ep_add;
        this->tp_glb_.resizeNoInitialize( this->n_glb_tot_ );
        this->epj_org_.resizeNoInitialize( this->n_glb_tot_ );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<offset; i++){
                this->tp_glb_[i] = this->tp_loc_[i]; // NOTE: need to keep tp_loc_[]?
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<n_ep_add; i++){
                this->epj_org_[offset+i] = this->epj_recv_[i];
                this->tp_glb_[offset+i].setFromEP(this->epj_recv_[i], offset+i);
            }
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTreeImpl(TagForceLong){
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 offset = this->n_loc_tot_;
        const S32 n_ep_add = this->n_ep_recv_disp_[n_proc];
        const S32 n_sp_add = this->n_sp_recv_disp_[n_proc];
        const S32 offset2 = this->n_loc_tot_ + n_ep_add;
        this->n_glb_tot_ = this->n_loc_tot_ + n_ep_add + n_sp_add;
        this->tp_glb_.resizeNoInitialize( this->n_glb_tot_ );
        this->epj_org_.resizeNoInitialize( offset2 );
        this->spj_org_.resizeNoInitialize( this->n_glb_tot_ );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<offset; i++){
                this->tp_glb_[i] = this->tp_loc_[i]; // NOTE: need to keep tp_loc_[]?
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<n_ep_add; i++){
                this->epj_org_[offset+i] = this->epj_recv_[i];
                this->tp_glb_[offset+i].setFromEP(this->epj_recv_[i], offset+i);
                //this->tp_glb_[offset+i].key_ = ClearMSB(this->tp_glb_[offset+i].key_);
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<n_sp_add; i++){
                this->spj_org_[offset2+i] = this->spj_recv_[i];
                this->tp_glb_[offset2+i].setFromSP(this->spj_recv_[i], offset2+i);
            }
        }

    }

    ///////////////////////////////
    /// morton sort global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortGlobalTreeOnly(const bool reuse){
        F64 time_offset = GetWtime();
        if(!map_id_to_epj_.empty()){
            map_id_to_epj_.clear();
        }
        assert(map_id_to_epj_.empty());

        tp_glb_.resizeNoInitialize(n_glb_tot_);
        epj_sorted_.resizeNoInitialize( n_glb_tot_ );
        if(!reuse){
            tp_buf_.resizeNoInitialize(n_glb_tot_);
#ifdef USE_STD_SORT
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );
#else
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf_.getPointer(), 0, n_glb_tot_-1);
#endif
        }
        if( typeid(TSM) == typeid(SEARCH_MODE_LONG)
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) 
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER) 
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF_SCATTER)
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_SYMMETRY) ){
            spj_sorted_.resizeNoInitialize( spj_org_.size() );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                if( GetMSB(adr) == 0){
                    epj_sorted_[i] = epj_org_[adr];
                }
                else{
                    spj_sorted_[i] = spj_org_[ClearMSB(adr)];
                }
            }
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                epj_sorted_[i] = epj_org_[adr];
            }
        }

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tp_glb_.size()="<<tp_glb_.size()<<" tp_buf_.size()="<<tp_buf_.size()<<std::endl;
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
        /*
        // remove 2016 1/5
        // new 2015 Aug 06
        for(S32 ith=0; ith<Comm::getNumberOfThread(); ith++){
            epj_neighbor_[ith].clearSize();
        }
        */
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

    /////////////////////////////
    // new for multiwalk (send index)
    /////////////////////////
    /// add moment as spj ///
    /*
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Ttreecell>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpImpl(TagForceLong,
                      const ReallocatableArray<Ttreecell> & _tc,
                      ReallocatableArray<Tspj> & _spj){
        S32 n_spj_prev = _spj_.size();
        _spj.resizeNoInitialize(n_spj_prev+_tc.size());
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif //PARTICLE_SIMULATOR_THREAD_PARALLEL
        for(S32 i=0; i<_tc.size(); i++){
            _spj[n_spj_prev+i].copyFromMoment(_tc[i].mom_);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Ttreecell>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpImpl(TagForceShort,
                      const ReallocatableArray<Ttreecell> & _tc,
                      ReallocatableArray<Tspj> & _spj){    
        // do nothing
    }
    */
    
    
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
            ipg_[i].vertex_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
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
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionList(const S32 adr_ipg, const bool clear){
        makeInteractionListImpl(typename TSM::search_type(), adr_ipg, clear);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLong, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
       if (clear){
           epj_for_force_[ith].clearSize();
           spj_for_force_[ith].clearSize();
       }
#if 0
        const S32 adr_root_cell = 0;
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        //const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 4.0;
        n_cell_open_[ith] += MakeInteractionListLongEPSPIteration
            (tc_glb_,                adr_root_cell,
             tp_glb_,                epj_sorted_, 
             epj_for_force_[ith],    spj_sorted_, 
             spj_for_force_[ith],    pos_target_box, 
             r_crit_sq,              n_leaf_limit_);
#else
        if (theta_ > 0.0) {
            const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
                MakeInteractionListLongEPSP
                    (tc_glb_, tc_glb_[0].adr_tc_, 
                     tp_glb_, epj_sorted_, 
                     epj_for_force_[ith],
                     spj_sorted_, spj_for_force_[ith],
                     pos_target_box, r_crit_sq, n_leaf_limit_);
            }
            else{
                const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
                if( pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq*4.0 ){
                    const S32 n_tmp = tc_glb_[0].n_ptcl_;
                    S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                    epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                    spj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                    for(S32 ip=0; ip<n_tmp; ip++){
                        if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                            epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                        }
                        else{
                            spj_for_force_[ith].pushBackNoCheck(spj_sorted_[adr_ptcl_tmp++]);
                        }
                    }
                }
            }
        } else {
            // theta_ = 0 case
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
#endif
    }



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLongCutoff, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
        const F64ort cell_box = pos_root_cell_;
        if (clear) {
            epj_for_force_[ith].clearSize();
            spj_for_force_[ith].clearSize();
        }
        if (theta_ > 0.0){
            const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
#if 0
            // 2017.10.31
            if(pos_target_box.getDistanceMinSQ(cell_box) <= r_cut_sq
               && tc_glb_[0].n_ptcl_ > 0){
                const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
                if( pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0){
                    if( tc_glb_[0].isLeaf(n_leaf_limit_) ){
                        const S32 n_tmp = tc_glb_[0].n_ptcl_;
                        S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                        epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                        spj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                        for(S32 ip=0; ip<n_tmp; ip++){
                            if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                                epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                            }
                            else{
                                spj_for_force_[ith].pushBackNoCheck(spj_sorted_[adr_ptcl_tmp++]);
                            }
                        }
                    }
                    else{
                        MakeInteractionListLongCutoffEPSP
                            (tc_glb_, tc_glb_[0].adr_tc_, tp_glb_, 
                             epj_sorted_, epj_for_force_[ith],
                             spj_sorted_, spj_for_force_[ith],
                             cell_box,
                             pos_target_box, r_crit_sq, r_cut_sq, n_leaf_limit_); 
                    }
                }
                else{
                    spj_for_force_[ith].increaseSize();
                    spj_for_force_[ith].back().copyFromMoment(tc_glb_[0].mom_);
                }
            }        
#else
            // original
            if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
                MakeInteractionListLongCutoffEPSP
                    (tc_glb_, tc_glb_[0].adr_tc_, tp_glb_, 
                     epj_sorted_, epj_for_force_[ith],
                     spj_sorted_, spj_for_force_[ith],
                     cell_box,
                     pos_target_box, r_crit_sq, r_cut_sq, n_leaf_limit_); 
            }
            else{
                //std::cerr<<"check b"<<std::endl;
                if(pos_target_box.getDistanceMinSQ(cell_box) <= r_cut_sq){
                    const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
                    if( pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0){
                        const S32 n_tmp = tc_glb_[0].n_ptcl_;
                        S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                        epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                        spj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                        for(S32 ip=0; ip<n_tmp; ip++){
                            if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                                epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                            }
                            else{
                                spj_for_force_[ith].pushBackNoCheck(spj_sorted_[adr_ptcl_tmp++]);
                            }
                        }
                    }
                    else{
                        spj_for_force_[ith].increaseSize();
                        spj_for_force_[ith].back().copyFromMoment(tc_glb_[0].mom_);
                    }
                }
            }
#endif
        } else {
            // theta_ = 0 case
            const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
            if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
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
    }


    // FOR P^3T
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLongScatter, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
        F64 r_crit_sq = LARGE_FLOAT;
        if(clear){
            epj_for_force_[ith].clearSize();
            spj_for_force_[ith].clearSize();
        }
        if (theta_ > 0.0) {
            r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
                MakeInteractionListLongScatterEPSP
                    (tc_glb_, tc_glb_[0].adr_tc_, 
                     tp_glb_, epj_sorted_, 
                     epj_for_force_[ith],
                     spj_sorted_, spj_for_force_[ith],
                     pos_target_box, r_crit_sq, n_leaf_limit_);
            }
            else{
                const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
                if( pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq*4.0 ){
                    const S32 n_tmp = tc_glb_[0].n_ptcl_;
                    S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                    epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                    spj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                    for(S32 ip=0; ip<n_tmp; ip++){
                        if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                            epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                        }
                        else{
                            spj_for_force_[ith].pushBackNoCheck(spj_sorted_[adr_ptcl_tmp++]);
                        }
                    }
                }
                else{
                    spj_for_force_[ith].increaseSize();
                    spj_for_force_[ith].back().copyFromMoment(tc_glb_[0].mom_);
                }
            }
        } else {
            // theta_ = 0.0 case
            const S32 n_tmp = tc_glb_[0].n_ptcl_;
            S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
            epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
            for (S32 ip=0; ip<n_tmp; ip++){
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
    makeInteractionListImpl(TagSearchLongSymmetry, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box_in = ipg_[adr_ipg].vertex_;
        const F64ort pos_target_box_out = ipg_[adr_ipg].vertex_out_;
        F64 r_crit_sq = LARGE_FLOAT;
        if(clear){
            epj_for_force_[ith].clearSize();
            spj_for_force_[ith].clearSize();
        }
        if (theta_ > 0.0) {
            r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
                MakeInteractionListLongSymmetryEPSP
                    (tc_glb_, tc_glb_[0].adr_tc_, 
                     tp_glb_, epj_sorted_, 
                     epj_for_force_[ith],
                     spj_sorted_, spj_for_force_[ith],
                     pos_target_box_in, pos_target_box_out, r_crit_sq, n_leaf_limit_);
            }
            else{
                const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
                if( pos_target_box_in.getDistanceMinSQ(pos_tmp) <= r_crit_sq*4.0 ){
                    const S32 n_tmp = tc_glb_[0].n_ptcl_;
                    S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                    epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                    spj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                    for(S32 ip=0; ip<n_tmp; ip++){
                        if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                            epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                        }
                        else{
                            spj_for_force_[ith].pushBackNoCheck(spj_sorted_[adr_ptcl_tmp++]);
                        }
                    }
                }
                else{
                    spj_for_force_[ith].increaseSize();
                    spj_for_force_[ith].back().copyFromMoment(tc_glb_[0].mom_);
                }
            }
        }
        else {
            // theta_ = 0 case
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

#if 0    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLongSymmetry, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        if(clear){
            epj_for_force_[ith].clearSize();
            spj_for_force_[ith].clearSize();
        }
        ReallocatableArray<Tspj> sp_list_dummy, sp_list_dummy2;
        MakeListUsingTreeRecursiveTop<TSM, MAKE_LIST_MODE_INTERACTION, LIST_CONTENT_ID, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
            (tc_glb_, 0, tp_glb_, epj_sorted_, epj_for_force_[ith],  id_ep_send_buf_[ith], spj_sorted_, spj_for_force_[ith],
             id_sp_send_buf_[ith], tree_top_inner_pos_[ib], tree_top_outer_pos_[ib], r_crit_sq, n_leaf_limit_);        
    }
#endif


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortScatter, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
        if (clear){
            epj_for_force_[ith].clearSize();
        }
#if 0
        const S32 adr_root_cell = 0;
        n_cell_open_[ith] += MakeListUsingOuterBoundaryIteration
            (tc_glb_.getPointer(),     adr_root_cell,
             epj_sorted_.getPointer(), epj_for_force_[ith], 
             pos_target_box,           n_leaf_limit_);
#else
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeListUsingOuterBoundary
                (tc_glb_.getPointer(),     tc_glb_[0].adr_tc_,
                 epj_sorted_.getPointer(), epj_for_force_[ith], 
                 pos_target_box,   n_leaf_limit_);
        }
        else{
            //if( pos_target_box.overlapped( tc_glb_[0].mom_.getVertexOut()) ){
            if( pos_target_box.contained( tc_glb_[0].mom_.getVertexOut()) ){
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                //epj_for_force_[ith].reserve( epj_for_force_[ith].size() + n_tmp );
                epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                    const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                    const F64 size_tmp = epj_sorted_[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq_tmp = pos_target_box.getDistanceMinSQ(pos_tmp);
                    if(dis_sq_tmp > size_tmp*size_tmp) continue;
                    //epj_for_force_[ith].increaseSize();
                    //epj_for_force_[ith].back() = epj_sorted_[adr_ptcl_tmp];
                    //epj_for_force_[ith].push_back(epj_sorted_[adr_ptcl_tmp]);
                    epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp]);
                    const F64vec pos_new = epj_for_force_[ith].back().getPos();
                    epj_for_force_[ith].back().setPos(pos_new);
                }
            }
        }
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortGather, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
       if(clear){
           epj_for_force_[ith].clearSize();
       }
#if 0
        const S32 adr_root_cell = 0;
        n_cell_open_[ith] += MakeListUsingInnerBoundaryIteration
            (tc_glb_.getPointer(),     adr_root_cell,
             epj_sorted_.getPointer(), epj_for_force_[ith],
             pos_target_box,           n_leaf_limit_);
#else
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeListUsingInnerBoundary
                (tc_glb_.getPointer(),     tc_glb_[0].adr_tc_,
                 epj_sorted_.getPointer(), epj_for_force_[ith],
                 pos_target_box,                 n_leaf_limit_);
        }
        else{
            if( pos_target_box.overlapped( tc_glb_[0].mom_.getVertexIn()) ){
                //if( pos_target_box.contained( tc_glb_[0].mom_.getVertexIn()) ){
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                //epj_for_force_[ith].reserve( epj_for_force_[ith].size() + n_tmp );
                epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                    const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                    if( pos_target_box.overlapped( pos_tmp) ){
                        //if( pos_target_box.contained( pos_tmp) ){
                        //epj_for_force_[ith].push_back(epj_sorted_[adr_ptcl_tmp]);
                        epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp]);
                        const F64vec pos_new = epj_for_force_[ith].back().getPos();
                        epj_for_force_[ith].back().setPos(pos_new);
                    }
                }
            }
        }
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortSymmetry, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box_out = (ipg_[adr_ipg]).vertex_;
        const F64ort pos_target_box_in = (ipg_[adr_ipg]).vertex_in;
        if(clear){
            epj_for_force_[ith].clearSize();
        }
#if 1
        const S32 adr_root_cell = 0;
        n_cell_open_[ith] += MakeListUsingOuterBoundaryAndInnerBoundaryIteration
            (tc_glb_.getPointer(),     adr_root_cell,
             epj_sorted_.getPointer(), epj_for_force_[ith],
             pos_target_box_out,       pos_target_box_in, 
             n_leaf_limit_);
#else
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeListUsingOuterBoundaryAndInnerBoundary
                (tc_glb_.getPointer(),     tc_glb_[0].adr_tc_,
                 epj_sorted_.getPointer(), epj_for_force_[ith],
                 pos_target_box_out,       pos_target_box_in,
                 n_leaf_limit_);
        }
        else{
            /*
            if( pos_target_box_out.overlapped(tc_glb_[0].mom_.getVertexIn()) 
                || pos_target_box_in.overlapped(tc_glb_[9].mom_.getVertexOut()) ){
            */
            if( pos_target_box_out.contained(tc_glb_[0].mom_.getVertexIn()) 
                || pos_target_box_in.contained(tc_glb_[9].mom_.getVertexOut()) ){
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                epj_for_force_[ith].reserveAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                    const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                    const F64 size_tmp = epj_sorted_[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq_tmp = pos_target_box_in.getDistanceMinSQ(pos_tmp);
                    //if( pos_target_box_out.notOverlapped(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                    if( pos_target_box_out.notContained(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                    epj_for_force_[ith].pushBackNoCheck( epj_sorted_[adr_ptcl_tmp] );
                    const F64vec pos_new = epj_for_force_[ith].back().getPos();
                    epj_for_force_[ith].back().setPos(pos_new);
                }
            }
        }
#endif
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

    /////////////////
    // CALC FORCE ///
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
    copyForceOriginalOrder(){
        force_org_.resizeNoInitialize(n_loc_tot_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = ClearMSB(tp_loc_[i].adr_ptcl_);
            force_org_[adr] = force_sorted_[i];
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForce(Tfunc_ep_ep pfunc_ep_ep,
              const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
        //PROFILE::Start(profile.calc_force);
        F64 offset_walk_tree,offset_dispatch;
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                offset_walk_tree = GetWtime();
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
        //PROFILE::Stop(profile.calc_force);
        copyForceOriginalOrder();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
        std::cout<<"force_sorted_.size()="<<force_sorted_.size()<<std::endl;
        std::cout<<"force_org_.size()="<<force_org_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceWalkOnly(Tfunc_ep_ep pfunc_ep_ep,
                      const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
        //PROFILE::Start(profile.calc_force);
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                makeInteractionList(i);
                ni_tmp += ipg_[i].n_ptcl_;
                nj_tmp += epj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * epj_for_force_[Comm::getThreadNum()].size();
                //calcForceOnly( pfunc_ep_ep, i, clear);
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
        //PROFILE::Stop(profile.calc_force);
        copyForceOriginalOrder();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
        std::cout<<"force_sorted_.size()="<<force_sorted_.size()<<std::endl;
        std::cout<<"force_org_.size()="<<force_org_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForce(Tfunc_ep_ep pfunc_ep_ep,
              Tfunc_ep_sp pfunc_ep_sp,
              const bool clear){
        const F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
        //PROFILE::Start(profile.calc_force);
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp, n_interaction_ep_sp_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                makeInteractionList(i);
#if 0
                if(Comm::getRank()==0){
                    std::cout<<"Multi:no, index:no i= "<<i
                             <<" vertex= "<<ipg_[i].vertex_
                             <<" n_epi= "<<ipg_[i].n_ptcl_
                             <<" n_epj= "<<epj_for_force_[Comm::getThreadNum()].size()
                             <<" n_spj= "<<spj_for_force_[Comm::getThreadNum()].size()
                             <<" tc_glb_[0].mom_.vertex_out_= "<<tc_glb_[0].mom_.vertex_out_
                             <<" tc_glb_[0].n_ptcl_= "<<tc_glb_[0].n_ptcl_
                             <<std::endl;
                }
#endif
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
        std::cout<<"force_org_.size()="<<force_org_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
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
        //std::cerr<<"epi_org_.size()= "<<epi_org_.size()<<" n_epj_tmp= "<<n_epj_tmp<<std::endl;
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
#if 0
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchLongScatter, const Tptcl & ptcl, S32 & nnp){
        const F64vec pos_target = ptcl.getPos();
        const S32 id_thread = Comm::getThreadNum();

        const S32 adr = 0;
        const S32 size_old = epj_neighbor_[id_thread].size();
/*
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),       adr, 
                                             epj_sorted_,   epj_neighbor_, n_leaf_limit_);
*/
/*
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),       
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
*/
        bool error = false;
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),       
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_,
                                             error);
        if (error) { nnp = -1; }
        else {
            nnp = epj_neighbor_[id_thread].size() - size_old;
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        const S32 head = epj_neighbor_[id_thread].size();
        //std::cerr<<"head="<<head<<std::endl;
        S32 nnp = 0;
        getNeighborListOneParticleImpl(typename TSM::search_type(), ptcl, nnp);
        epj = epj_neighbor_[id_thread].getPointer(head);
        if(nnp == -1){
            epj_neighbor_[id_thread].clearSize();
        }
        return nnp;
    }
#else
    
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
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        return getNeighborListOneParticleImpl(typename TSM::search_type(), ptcl, epj);
    }
    
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
#endif

    
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
    freeMem(){
        tp_buf_.freeMem();
        tp_loc_.freeMem();
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
        epj_recv_.freeMem();
        spj_send_.freeMem();
        spj_recv_.freeMem();
        force_sorted_.freeMem();
        force_org_.freeMem();
        epjr_sorted_.freeMem();
        epjr_send_.freeMem();
        epjr_recv_.freeMem();
        epjr_recv_1st_buf_.freeMem();
        epjr_recv_2nd_buf_.freeMem();
        epj_recv_1st_buf_.freeMem();
        epj_recv_2nd_buf_.freeMem();
        const S32 n_thread = Comm::getNumberOfThread();
        for(S32 i=0; i<n_thread; i++){
            id_ep_send_buf_[i].freeMem();
            id_sp_send_buf_[i].freeMem();
            epj_for_force_[i].freeMem();
            spj_for_force_[i].freeMem();
            ep_send_buf_for_scatter_[i].freeMem();
            shift_image_domain_[i].freeMem();
            epjr_send_buf_[i].freeMem();
            epjr_send_buf_for_scatter_[i].freeMem();
            epjr_recv_1st_sorted_[i].freeMem();
            epj_send_buf_[i].freeMem();
            id_ptcl_send_[i].freeMem();
            shift_image_box_[i].freeMem();
            ip_disp_[i].freeMem();
            tp_scatter_[i].freeMem();
            tc_recv_1st_[i].freeMem();
            epj_recv_1st_sorted_[i].freeMem();
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    reallocMem(){
        tp_buf_.reallocMem();
        tp_loc_.reallocMem();
        tp_glb_.reallocMem();
        tc_loc_.reallocMem();
        tc_glb_.reallocMem();
        epi_sorted_.reallocMem();
        epi_org_.reallocMem();
        epj_sorted_.reallocMem();
        epj_org_.reallocMem();
        spj_sorted_.reallocMem();
        spj_org_.reallocMem();
        ipg_.reallocMem();
        epj_send_.reallocMem();
        epj_recv_.reallocMem();
        spj_send_.reallocMem();
        spj_recv_.reallocMem();
        force_sorted_.reallocMem();
        force_org_.reallocMem();
        epjr_sorted_.reallocMem();
        epjr_send_.reallocMem();
        epjr_recv_.reallocMem();
        epjr_recv_1st_buf_.reallocMem();
        epjr_recv_2nd_buf_.reallocMem();
        epj_recv_1st_buf_.reallocMem();
        epj_recv_2nd_buf_.reallocMem();
        const S32 n_thread = Comm::getNumberOfThread();
        for(S32 i=0; i<n_thread; i++){
            id_ep_send_buf_[i].reallocMem();
            id_sp_send_buf_[i].reallocMem();
            epj_for_force_[i].reallocMem();
            spj_for_force_[i].reallocMem();
            ep_send_buf_for_scatter_[i].reallocMem();
            shift_image_domain_[i].reallocMem();
            epjr_send_buf_[i].reallocMem();
            epjr_send_buf_for_scatter_[i].reallocMem();
            epjr_recv_1st_sorted_[i].reallocMem();
            epj_send_buf_[i].reallocMem();
            id_ptcl_send_[i].reallocMem();
            shift_image_box_[i].reallocMem();
            ip_disp_[i].reallocMem();
            tp_scatter_[i].reallocMem();
            tc_recv_1st_[i].reallocMem();
            epj_recv_1st_sorted_[i].reallocMem();
        }
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    clearSizeOfArray(){
        tp_buf_.clearSize();
        tp_loc_.clearSize();
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
        epj_recv_.clearSize();
        spj_send_.clearSize();
        spj_recv_.clearSize();
        force_sorted_.clearSize();
        force_org_.clearSize();
        epjr_sorted_.clearSize();
        epjr_send_.clearSize();
        epjr_recv_.clearSize();
        epjr_recv_1st_buf_.clearSize();
        epjr_recv_2nd_buf_.clearSize();
        epj_recv_1st_buf_.clearSize();
        epj_recv_2nd_buf_.clearSize();
        const S32 n_thread = Comm::getNumberOfThread();
        for(S32 i=0; i<n_thread; i++){
            id_ep_send_buf_[i].clearSize();
            id_sp_send_buf_[i].clearSize();
            epj_for_force_[i].clearSize();
            spj_for_force_[i].clearSize();
            ep_send_buf_for_scatter_[i].clearSize();
            shift_image_domain_[i].clearSize();
            epjr_send_buf_[i].clearSize();
            epjr_send_buf_for_scatter_[i].clearSize();
            epjr_recv_1st_sorted_[i].clearSize();
            epj_send_buf_[i].clearSize();
            id_ptcl_send_[i].clearSize();
            shift_image_box_[i].clearSize();
            ip_disp_[i].clearSize();
            tp_scatter_[i].clearSize();
            tc_recv_1st_[i].clearSize();
            epj_recv_1st_sorted_[i].clearSize();
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
        exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo,flag_reuse);
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
            //comm_table_.clearSize();
            FindScatterParticle<TSM, TreeCell<Tmomloc>, TreeParticle,
                                Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_loc_,
                 epj_sorted_, 
                 comm_table_.n_ep_send_,   comm_table_.adr_ep_send_, 
                 dinfo,          n_leaf_limit_,
                 spj_sorted_loc_,
                 comm_table_.n_sp_send_,   comm_table_.adr_sp_send_, 
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.n_sp_per_image_,
                 r_crit_sq);
            ExchangeNumber(comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                           comm_table_.n_sp_send_, comm_table_.n_sp_recv_);
        }
        ExchangeLet<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_,
                                     comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                     comm_table_.adr_ep_send_, epj_recv_,
                                     spj_sorted_loc_, comm_table_.n_sp_send_,
                                     comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
                                     comm_table_.adr_sp_send_, spj_recv_,
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
                (tc_loc_, tp_loc_, epj_sorted_,
                 comm_table_.n_ep_send_,  comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_);
            ExchangeNumber(comm_table_.n_ep_send_, comm_table_.n_ep_recv_);
        }
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_, epj_recv_,
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
                (tc_loc_, tp_loc_, epj_sorted_,
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
            
#if 0
            std::vector< std::vector<S32> > adr_per_proc(n_proc, std::vector<S32>(0) );
            if(Comm::getRank()==0){
                for(S32 i=0; i<n_proc; i++){
                    /*
                    std::cerr<<"rank= "<<i
                             <<" n_ep_send_per_proc_1st[i]= "<<n_ep_send_per_proc_1st[i]
                             <<" n_ep_send_per_proc_2nd[i]= "<<n_ep_send_per_proc_2nd[i]
                             <<std::endl;
                    */
                    for(S32 j=n_disp_ep_send_per_proc_1st[i]; j<n_disp_ep_send_per_proc_1st[i+1]; j++){
                        //std::cerr<<"adr_ep_send_1st[j]= "<<adr_ep_send_1st[j]<<std::endl;
                        adr_per_proc[i].push_back(adr_ep_send_1st[j]);
                    }
                    for(S32 j=n_disp_ep_send_per_proc_2nd[i]; j<n_disp_ep_send_per_proc_2nd[i+1]; j++){
                        //std::cerr<<"adr_ep_send_2nd[j]= "<<adr_ep_send_2nd[j]<<std::endl;
                        adr_per_proc[i].push_back(adr_ep_send_2nd[j]);
                    }
                    if(n_ep_send_per_proc_1st[i] > 0 && n_ep_send_per_proc_2nd[i] > 0){
                        std::sort(adr_per_proc[i].begin(), adr_per_proc[i].end());
                        for(S32 j=1; j<adr_per_proc[i].size(); j++){
                            assert(adr_per_proc[i][j-1] != adr_per_proc[i][j]);
                        }
                        //for(auto x : adr_per_proc[i]) std::cerr<<"x= "<<x<<std::endl;
                    }
                }
            }
#endif
        } // end of reuse
        /*
        if(my_rank == 0){
            std::cerr<<"Before exchange LET"<<std::endl;
            std::cerr<<"comm_table_.n_ep_send_.size()= "<<comm_table_.n_ep_send_.size()<<std::endl;
            std::cerr<<"comm_table_.n_ep_recv_.size()= "<<comm_table_.n_ep_recv_.size()<<std::endl;
            std::cerr<<"comm_table_.n_image_per_proc_.size()= "<<comm_table_.n_image_per_proc_.size()
                     <<std::endl;
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"i= "<<i
                         <<" n_ep_send_[i]= "<<comm_table_.n_ep_send_[i]
                         <<" n_ep_recv_[i]= "<<comm_table_.n_ep_recv_[i]
                         <<" n_image_per_proc_[i]= "<<comm_table_.n_image_per_proc_[i]
                         <<std::endl;
            }
            std::cerr<<"comm_table_.n_ep_per_image_.size()= "<<comm_table_.n_ep_per_image_.size()
                     <<std::endl;
            std::cerr<<"comm_table_.shift_per_image_.size()= "<<comm_table_.shift_per_image_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.n_ep_per_image_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" n_ep_per_image_[i]= "<<comm_table_.n_ep_per_image_[i]
                         <<" shift_per_image_[i]= "<<comm_table_.shift_per_image_[i]
                         <<std::endl;
            }
            std::cerr<<"comm_table_.adr_ep_send_.size()= "<<comm_table_.adr_ep_send_.size()
                     <<std::endl;
        }
        */
        //Comm::barrier();
        //exit(1);
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_, epj_recv_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);
    }


    ////////////////
    // GATHER MODE
#if 1
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchShortGather,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        const S32 n_proc = Comm::getNumberOfProc();
        //const S32 my_rank = Comm::getRank();
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
            /*
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_loc_, epj_sorted_,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            */
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, EssentialParticleBase, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_loc_, epi_base_sorted,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            ExchangeNumber(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st);
            // exchange epi_base particles
            /*
            ExchangeParticle(epi_sorted_, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            */
            ExchangeParticle(epi_base_sorted, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            
            ////////////
            // 2nd STEP (find j particle)
            /*
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, Tepi>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epi_sorted_,
                 center_, length_);
            */
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
            
#if 0
            std::vector< std::vector<S32> > adr_per_proc(n_proc, std::vector<S32>(0) );
            if(Comm::getRank()==0){
                for(S32 i=0; i<n_proc; i++){
                    for(S32 j=n_disp_ep_send_per_proc_1st[i]; j<n_disp_ep_send_per_proc_1st[i+1]; j++){
                        //std::cerr<<"adr_ep_send_1st[j]= "<<adr_ep_send_1st[j]<<std::endl;
                        adr_per_proc[i].push_back(adr_ep_send_1st[j]);
                    }
                    for(S32 j=n_disp_ep_send_per_proc_2nd[i]; j<n_disp_ep_send_per_proc_2nd[i+1]; j++){
                        //std::cerr<<"adr_ep_send_2nd[j]= "<<adr_ep_send_2nd[j]<<std::endl;
                        adr_per_proc[i].push_back(adr_ep_send_2nd[j]);
                    }
                    if(n_ep_send_per_proc_1st[i] > 0 && n_ep_send_per_proc_2nd[i] > 0){
                        std::sort(adr_per_proc[i].begin(), adr_per_proc[i].end());
                        for(S32 j=1; j<adr_per_proc[i].size(); j++){
                            assert(adr_per_proc[i][j-1] != adr_per_proc[i][j]);
                        }
                        //for(auto x : adr_per_proc[i]) std::cerr<<"x= "<<x<<std::endl;
                    }
                }
            }
#endif
        } // end of reuse flag
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_, epj_recv_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);
    }
#else
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchShortGather,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
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
        static ReallocatableArray<Tepi> ep_recv_1st;
        if(!flag_reuse){
            ////////////
            // 1st STEP (send j particles)
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_loc_, epj_sorted_,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            ExchangeNumber(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st);
            /*
            ExchangeParticle(epj_sorted_, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            */
            // exchange i particles
            ExchangeParticle(epi_sorted_, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            
            ////////////
            // 2nd STEP (send j particles)
            /*
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, Tepj>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epj_sorted_,
                 center_, length_);
            */
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, Tepi>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epi_sorted_,
                 center_, length_);
            /*
            if(Comm::getRank()==1){
                std::cerr<<"n_ep_send_per_proc_2nd.size()= "<<n_ep_send_per_proc_2nd.size()<<std::endl;
                for(S32 i=0; i<n_proc; i++){
                    std::cerr<<"i= "<<i
                             <<" n_ep_send_per_proc_2nd[i]= "<<n_ep_send_per_proc_2nd[i]
                             <<std::endl;
                }
                S32 n_image_tmp = shift_per_image_2nd.size();
                std::cerr<<"n_image_tmp= "<<n_image_tmp<<std::endl;
                for(S32 i=0; i<n_image_tmp; i++){
                    std::cerr<<"i= "<<i
                             <<" shift_per_image_2nd[i]= "<<shift_per_image_2nd[i]
                             <<std::endl;
                }
            }
            Comm::barrier();
            exit(1);
            */

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
            
#if 0
            std::vector< std::vector<S32> > adr_per_proc(n_proc, std::vector<S32>(0) );
            if(Comm::getRank()==0){
                for(S32 i=0; i<n_proc; i++){
                    /*
                    std::cerr<<"rank= "<<i
                             <<" n_ep_send_per_proc_1st[i]= "<<n_ep_send_per_proc_1st[i]
                             <<" n_ep_send_per_proc_2nd[i]= "<<n_ep_send_per_proc_2nd[i]
                             <<std::endl;
                    */
                    for(S32 j=n_disp_ep_send_per_proc_1st[i]; j<n_disp_ep_send_per_proc_1st[i+1]; j++){
                        //std::cerr<<"adr_ep_send_1st[j]= "<<adr_ep_send_1st[j]<<std::endl;
                        adr_per_proc[i].push_back(adr_ep_send_1st[j]);
                    }
                    for(S32 j=n_disp_ep_send_per_proc_2nd[i]; j<n_disp_ep_send_per_proc_2nd[i+1]; j++){
                        //std::cerr<<"adr_ep_send_2nd[j]= "<<adr_ep_send_2nd[j]<<std::endl;
                        adr_per_proc[i].push_back(adr_ep_send_2nd[j]);
                    }
                    if(n_ep_send_per_proc_1st[i] > 0 && n_ep_send_per_proc_2nd[i] > 0){
                        std::sort(adr_per_proc[i].begin(), adr_per_proc[i].end());
                        for(S32 j=1; j<adr_per_proc[i].size(); j++){
                            assert(adr_per_proc[i][j-1] != adr_per_proc[i][j]);
                        }
                        //for(auto x : adr_per_proc[i]) std::cerr<<"x= "<<x<<std::endl;
                    }
                }
            }
#endif
        } // end of reuse flag
        /*
        if(my_rank == 0){
            std::cerr<<"Before exchange LET"<<std::endl;
            std::cerr<<"comm_table_.n_ep_send_.size()= "<<comm_table_.n_ep_send_.size()<<std::endl;
            std::cerr<<"comm_table_.n_ep_recv_.size()= "<<comm_table_.n_ep_recv_.size()<<std::endl;
            std::cerr<<"comm_table_.n_image_per_proc_.size()= "<<comm_table_.n_image_per_proc_.size()
                     <<std::endl;
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"i= "<<i
                         <<" n_ep_send_[i]= "<<comm_table_.n_ep_send_[i]
                         <<" n_ep_recv_[i]= "<<comm_table_.n_ep_recv_[i]
                         <<" n_image_per_proc_[i]= "<<comm_table_.n_image_per_proc_[i]
                         <<std::endl;
            }
            std::cerr<<"comm_table_.n_ep_per_image_.size()= "<<comm_table_.n_ep_per_image_.size()
                     <<std::endl;
            std::cerr<<"comm_table_.shift_per_image_.size()= "<<comm_table_.shift_per_image_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.n_ep_per_image_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" n_ep_per_image_[i]= "<<comm_table_.n_ep_per_image_[i]
                         <<" shift_per_image_[i]= "<<comm_table_.shift_per_image_[i]
                         <<std::endl;
            }
            std::cerr<<"comm_table_.adr_ep_send_.size()= "<<comm_table_.adr_ep_send_.size()
                     <<std::endl;
        }
        Comm::barrier();
        exit(1);
        */
        /*
        const S32 n_loc = epj_sorted_.size();
        static ReallocatableArray<EPJWithR> epjr_sorted;
        epjr_sorted.resizeNoInitialize(n_loc);
        for(S32 i=0; i<n_loc; i++){
            epjr_sorted[i].copyFromEPJ(epj_sorted_[i]);
            epjr_sorted[i].r_search = epi_sorted_[i].getRSearch();
        }
        ExchangeLet(epjr_sorted, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_, epjr_recv_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);
        //Comm::barrier();
        //exit(1);
        */
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_, epj_recv_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);
    }
#endif
    
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
        const S32 adr_tree_sp_first = 0;
        const F64 r_crit_sq = 999.9;
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
                //const F64ort pos_target_box = GetIpgBoxForInteractionList<TSM>(ipg_[i]);
                TargetBox<TSM> target_box;
                GetTargetBox<TSM>(ipg_[i], target_box);
                /*
                if(Comm::getRank()==0){
                    std::cout<<"i= "<<i
                             <<" vertex_out= "<<target_box.vertex_out_
                             <<" vertex_in= "<<target_box.vertex_in_
                             <<std::endl;
                }
                */
                MakeListUsingTreeRecursiveTop
                    <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                     WALK_MODE_NORMAL, TagChopLeafTrue>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     spj_dummy,   adr_spj,
                     //pos_target_box,
                     target_box,
                     r_crit_sq, n_leaf_limit_,
                     adr_tree_sp_first, F64vec(0.0));
                interaction_list_.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                n_ep_cum_prev = adr_epj_tmp[ith].size();
                /*
                if(Comm::getRank()==0){
                    for(S32 j=0; j<n_ep_cum_prev; j++){
                        std::cerr<<"j= "<<j
                                 <<" adr_epj_tmp[ith][j]= "<<adr_epj_tmp[ith][j]<<std::endl;
                    }
                }
                */
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
        const S32 adr_tree_sp_first = spj_org_.size();
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
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
                TargetBox<TSM> target_box;
                GetTargetBox<TSM>(ipg_[i], target_box);
                /*
                if(Comm::getRank()==0){
                    std::cout<<"i= "<<i
                             <<" vertex_out= "<<target_box.vertex_out_
                             <<" vertex_in= "<<target_box.vertex_in_
                             <<std::endl;
                }
                */
                /*
                if(Comm::getRank()==0){
                    std::cerr<<"i= "<<i
                             <<" tc_glb[0].mom_.getVertexOut()= "<<tc_glb_[0].mom_.getVertexOut()
                             <<" tc_glb[0].mom_.getVertexIn()= "<<tc_glb_[0].mom_.getVertexIn()
                             <<" target_box.vertex_out_= "<<target_box.vertex_out_
                             <<" target_box.vertex_in_= "<<target_box.vertex_in_
                             <<std::endl;
                }
                */
                MakeListUsingTreeRecursiveTop
                    <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                     WALK_MODE_NORMAL, TagChopLeafTrue>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     spj_sorted_, adr_spj_tmp[ith],
                     target_box,
                     r_crit_sq, n_leaf_limit_,
                     adr_tree_sp_first, F64vec(0.0));
                /*
                if(Comm::getRank()==0){
                    std::cout<<"i= "<<i
                             <<" n_epi= "<<ipg_[i].n_ptcl_
                             <<" n_epj= "<<adr_epj_tmp[ith].size()-n_ep_cum_prev
                             <<" n_spj= "<<adr_spj_tmp[ith].size()-n_sp_cum_prev
                             <<std::endl;
                }
                */
#if 0
                if(Comm::getRank()==0){
                    std::cerr<<"n_ep_cum_prev= "<<n_ep_cum_prev
                             <<" n_sp_cum_prev= "<<n_sp_cum_prev
                             <<std::endl;
                    F64 mass_tmp = 0.0;
                    for(S32 j=n_ep_cum_prev; j<adr_epj_tmp[ith].size(); j++){
                        /*
                        std::cerr<<"j= "<<j
                                 <<" adr_epj_tmp[ith][j]= "<<adr_epj_tmp[ith][j]
                                 <<std::endl;
                        */
                        mass_tmp += epj_sorted_[adr_epj_tmp[ith][j]].getCharge();
                    }
                    for(S32 j=n_sp_cum_prev; j<adr_spj_tmp[ith].size(); j++){
                        /*
                        std::cerr<<"j= "<<j
                                 <<" adr_spj_tmp[ith][j]= "<<adr_spj_tmp[ith][j]
                                 <<std::endl;
                        */
                        mass_tmp += spj_sorted_[adr_spj_tmp[ith][j]].getCharge();
                    }
                    std::cerr<<"mass_tmp= "<<mass_tmp<<std::endl;
                }
#endif
                
                interaction_list_.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                interaction_list_.n_sp_[i] = adr_spj_tmp[ith].size() - n_sp_cum_prev;
                n_ep_cum_prev = adr_epj_tmp[ith].size();
                n_sp_cum_prev = adr_spj_tmp[ith].size();
            }
            n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
            n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
            /*
            std::cerr<<"A) ith= "<<ith
                     <<" n_disp_epj_tmp[ith].size()= "<<n_disp_epj_tmp[ith].size()
                     <<" n_disp_spj_tmp[ith].size()= "<<n_disp_spj_tmp[ith].size()
                         <<std::endl;
            */
        } // end of OMP

        interaction_list_.n_disp_ep_[0] = 0;
        interaction_list_.n_disp_sp_[0] = 0;
        for(S32 i=0; i<n_ipg; i++){
            interaction_list_.n_disp_ep_[i+1] = interaction_list_.n_disp_ep_[i] + interaction_list_.n_ep_[i];
            interaction_list_.n_disp_sp_[i+1] = interaction_list_.n_disp_sp_[i] + interaction_list_.n_sp_[i];
        }
        interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[n_ipg] );
        interaction_list_.adr_sp_.resizeNoInitialize( interaction_list_.n_disp_sp_[n_ipg] );
        /*
        if(Comm::getRank()==0){
            std::cerr<<"interaction_list_.adr_ep_.size()= "<<interaction_list_.adr_ep_.size()
                     <<" interaction_list_.adr_sp_.size()= "<<interaction_list_.adr_sp_.size()
                     <<std::endl;
        }
        */
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_thread; i++){
                std::cerr<<"i= "<<i
                         <<" n_disp_epj_tmp[i].size()= "<<n_disp_epj_tmp[i].size()
                         <<" n_disp_spj_tmp[i].size()= "<<n_disp_spj_tmp[i].size()
                         <<std::endl;
                for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                    const S32 adr_ipg = adr_ipg_tmp[i][j];
                    S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                    const S32 k_ep_h = n_disp_epj_tmp[i][j];
                    const S32 k_ep_e = n_disp_epj_tmp[i][j+1];
                    S32 adr_sp = interaction_list_.n_disp_sp_[adr_ipg];
                    const S32 k_sp_h = n_disp_spj_tmp[i][j];
                    const S32 k_sp_e = n_disp_spj_tmp[i][j+1];
                    std::cerr<<"k_ep_h= "<<k_ep_h<<" k_ep_e= "<<k_ep_e
                             <<" k_sp_h= "<<k_sp_h<<" k_sp_e= "<<k_sp_e
                             <<std::endl;
                }
            }
        }
        */
        //Comm::barrier();
        //exit(1);        
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
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<ipg_.size(); i++){
                for(S32 j=interaction_list_.n_disp_sp_[i]; j<interaction_list_.n_disp_sp_[i+1]; j++){
                    std::cerr<<"adr_sp_[j]= "<<interaction_list_.adr_sp_[j]<<std::endl;
                }
            }
        }
        Comm::barrier();
        exit(1);
        */
    }
    

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                    const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
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
                /*
                if(Comm::getRank() == 0){
                    std::cerr<<"adr_epi_head= "<<adr_epi_head
                             <<" n_epj= "<<n_epj
                             <<" force_sorted_[adr_epi_head].n_ngb= "<<force_sorted_[adr_epi_head].n_ngb<<std::endl;
                    for(S32 j=adr_epj_head; j<adr_epj_end; j++){
                        const S32 adr_epj = interaction_list_.adr_ep_[j];
                        std::cerr<<"j= "<<j
                                 <<" epj_sorted_[adr_epj].pos= "<<epj_sorted_[adr_epj].pos<<std::endl;
                    }
                }
                */
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
        force_org_.resizeNoInitialize(n_loc_tot_);
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
                //F64 mass_tmp = 0.0;
                for(S32 j=adr_epj_head; j<adr_epj_end; j++, n_ep_cnt++){
                    const S32 adr_epj = interaction_list_.adr_ep_[j];
                    epj_for_force_[ith][n_ep_cnt] = epj_sorted_[adr_epj];
                    //mass_tmp += epj_sorted_[adr_epj].mass;
                }
                pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            epj_for_force_[ith].getPointer(),   n_epj,
                            force_sorted_.getPointer(adr_epi_head));
                S32 n_sp_cnt = 0;
                for(S32 j=adr_spj_head; j<adr_spj_end; j++, n_sp_cnt++){
                    const S32 adr_spj = interaction_list_.adr_sp_[j];
                    spj_for_force_[ith][n_sp_cnt] = spj_sorted_[adr_spj];
                    //mass_tmp += spj_sorted_[adr_spj].mass;
                }
                //std::cerr<<"mass_tmp= "<<mass_tmp<<std::endl;
                pfunc_ep_sp(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            spj_for_force_[ith].getPointer(),   n_spj,
                            force_sorted_.getPointer(adr_epi_head));
                /*
                if(Comm::getRank() == 0){
                    std::cerr<<"adr_epi_head= "<<adr_epi_head
                             <<" n_epj= "<<n_epj
                             <<" n_spj= "<<n_spj
                             <<std::endl;
                    for(S32 j=adr_epj_head; j<adr_epj_end; j++){
                        const S32 adr_epj = interaction_list_.adr_ep_[j];
                        std::cerr<<"j= "<<j
                                 <<" epj_sorted_[adr_epj].pos= "<<epj_sorted_[adr_epj].pos<<std::endl;
                    }
                    for(S32 j=adr_spj_head; j<adr_spj_end; j++){
                        const S32 adr_spj = interaction_list_.adr_sp_[j];
                        std::cerr<<"j= "<<j
                                 <<" spj_sorted_[adr_spj].pos= "<<spj_sorted_[adr_spj].pos<<std::endl;
                    }
                }
                */
            }
        }
        n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
        copyForceOriginalOrder();
        time_profile_.calc_force += GetWtime() - time_offset;
    }




    

}

#include"tree_for_force_impl_force.hpp"
