/////////////////////////////////////////////////////////
// implementaion of check API methods of TreeForForce ///

namespace ParticleSimulator{
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMortonSortLocalTreeOnly(std::ostream & fout){
        CheckMortonSort(n_loc_tot_, tp_loc_.getPointer(), epj_sorted_.getPointer(), fout);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMortonSortGlobalTreeOnly(std::ostream & fout){
        checkMortonSortGlobalTreeOnlyImpl(typename TSM::force_type(), fout);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMortonSortGlobalTreeOnlyImpl(TagForceLong, std::ostream & fout){
        CheckMortonSort(n_glb_tot_, tp_glb_.getPointer(),
                        epj_sorted_.getPointer(), spj_sorted_.getPointer(), fout);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMortonSortGlobalTreeOnlyImpl(TagForceShort, std::ostream & fout){
        CheckMortonSort(n_glb_tot_, tp_glb_.getPointer(), epj_sorted_.getPointer(), fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeLocalTree(const F64 tolerance, std::ostream & fout){
        S32 err = 0;
        F64vec center = center_;
        tc_loc_.getPointer()->checkTree(epj_sorted_.getPointer(), tc_loc_.getPointer(),
                                        center, length_*0.5, n_leaf_limit_,
                                        tolerance, err, fout);
        if(!err) fout<<"makeLocalTree test: PASS"<<std::endl;
        else  fout<<"makeLocalTree test: FAIL err="<<err<<std::endl;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeGlobalTree(const F64 tolerance, std::ostream & fout){
        S32 err = 0;
        F64vec center = center_;
        checkMakeGlobalTreeImpl(typename TSM::force_type(),
                                err, center, tolerance, fout);
        if(!err) fout<<"makeGlobalTree test: PASS"<<std::endl;
        else  fout<<"makeGlobalTree test: FAIL err="<<err<<std::endl;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeGlobalTreeImpl(TagForceShort,
                            S32 & err,
                            const F64vec & center,
                            const F64 tolerance,
                            std::ostream & fout){
        fout<<"CheckMakeGlobalTreeDummy (SCATTER, GATHER, SYMMETRY)"<<std::endl;
        tc_glb_[0].checkTree(epj_sorted_.getPointer(),
                             tc_glb_.getPointer(), 
                             center, length_*0.5,
                             n_leaf_limit_, tolerance, err, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeGlobalTreeImpl(TagForceLong,
                            S32 & err,
                            const F64vec & center,
                            const F64 tolerance,
                            std::ostream & fout){
        std::cout<<"CheckMakeGlobalTreeDummy (LONG, LONGCUTOFF)"<<std::endl;
        tc_glb_.getPointer()->checkTreeLongGlobalTree
            (epj_sorted_.getPointer(), spj_sorted_.getPointer(), 
             tp_glb_.getPointer(),     tc_glb_.getPointer(), 
             center,  length_*0.5, n_leaf_limit_, tolerance, err, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkExchangeLocalEssentialTree(const DomainInfo & dinfo, 
                                    const F64 tolerance,
                                    std::ostream & fout){
        fout<<"checkExchangeLocalEssentialTree..."<<std::endl;
        checkExchangeLocalEssentialTreeImpl(typename TSM::force_type(), dinfo, tolerance, fout);
    }
    

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkExchangeLocalEssentialTreeImpl(TagForceLong,
                                        const DomainInfo & dinfo,
                                        const F64 tolerance,
                                        std::ostream & fout){
        checkExchangeLocalEssentialTreeForLongImpl(typename TSM::search_type(), dinfo, tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkExchangeLocalEssentialTreeForLongImpl(TagSearchLong,
                                               const DomainInfo & dinfo,
                                               const F64 tolerance, 
                                               std::ostream & fout){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	    
        std::cout<<"CheckExchangeLocalEssentialTreeDummy long"<<std::endl;
        F64 mass_cm_loc_direct = 0.0;
        F64vec pos_cm_loc_direct = 0.0;
        for(S32 i=0; i<this->n_loc_tot_; i++){
            mass_cm_loc_direct += this->epj_sorted_[i].getCharge();
            pos_cm_loc_direct += this->epj_sorted_[i].getCharge() * this->epj_sorted_[i].getPos();
        }
        F64 mass_cm_glb_direct = Comm::getSum(mass_cm_loc_direct);
        F64vec pos_cm_glb_direct = Comm::getSum(pos_cm_loc_direct);
        pos_cm_glb_direct /= mass_cm_glb_direct;
        F64 mass_cm_glb_tree = this->tc_loc_[0].mom_.mass;
        F64vec pos_cm_glb_tree = this->tc_loc_[0].mom_.pos * this->tc_loc_[0].mom_.mass;
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 n_ep = this->n_ep_recv_disp_[n_proc];
        const S32 n_sp = this->n_sp_recv_disp_[n_proc];
        fout<<"this->n_ep_recv_disp_[n_proc]="<<this->n_ep_recv_disp_[n_proc]<<std::endl;
        fout<<"this->n_sp_recv_disp_[n_proc]="<<this->n_sp_recv_disp_[n_proc]<<std::endl;
        for(S32 i=0; i<n_ep; i++){
            mass_cm_glb_tree += this->epj_recv_[i].getCharge();
            pos_cm_glb_tree += this->epj_recv_[i].getCharge() * this->epj_recv_[i].getPos();
        }
        for(S32 i=0; i<n_sp; i++){
            mass_cm_glb_tree += this->spj_recv_[i].getCharge();
            pos_cm_glb_tree += this->spj_recv_[i].getCharge() * this->spj_recv_[i].getPos();
        }
        pos_cm_glb_tree /= mass_cm_glb_tree;
        F64 dm = std::abs(mass_cm_glb_direct - mass_cm_glb_tree);
        F64vec dx = pos_cm_glb_direct - pos_cm_glb_tree;
        F64 dx_max = dx.applyEach( Abs<F64>() ).getMax();
        if( dm >= tolerance || dx_max >= tolerance){
            fout<<"CheckExchangeLocalEssentialTree (long): FAIL"<<std::endl;
        }
        else{
            fout<<"CheckExchangeLocalEssentialTree (long): PASS"<<std::endl;
        }
        fout<<"mass_cm_glb_direct="<<mass_cm_glb_direct<<" mass_cm_glb_tree="<<mass_cm_glb_tree<<std::endl;
        fout<<"pos_cm_glb_direct="<<pos_cm_glb_direct<<" pos_cm_glb_tree="<<pos_cm_glb_tree<<std::endl;
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL	    
    }

    // compare neighbouring mass(include own mass).
    // this function has NOT been checked in the case of PERIODIC BOUNDARY CONDITION
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkExchangeLocalEssentialTreeForLongImpl(TagSearchLongCutoff,
                                                        const DomainInfo & dinfo,
                                                        const F64 tolerance, 
                                                        std::ostream & fout){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL 	    
        fout<<"checkExchangeLocalEssentialTreeForLongImpl (for long w/ cutoff) is under construction"<<std::endl;
        fout<<"check if the total mass of the neighbour particle by using the tree is SIMILAR (NOT have to be the same) to that evaluated directly"<<std::endl;
        const S32 n_proc = Comm::getNumberOfProc();
        Tepj * ep_tmp;
        S32 * n_ptcl;
        S32 * n_ptcl_disp;
        const S32 n_loc = this->n_loc_tot_;
        AllGatherParticle(ep_tmp, n_ptcl, n_ptcl_disp, this->epj_org_.getPointer(), n_loc);
        // for PERIODIC
        ReallocatableArray<F64vec> shift_image_domain;
        shift_image_domain.resizeNoInitialize(0);
        if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
            shift_image_domain.push_back( F64vec(0.0) );
        }
        else{
            bool pa[DIMENSION];
            dinfo.getPeriodicAxis(pa);
            shift_image_domain.reserve(5*5*5);
            F64ort outer_boundary_of_my_tree = this->tc_loc_[0].mom_.vertex_out_;
            CalcNumberAndShiftOfImageDomain
                (shift_image_domain, dinfo.getPosRootDomain().getFullLength(),
                 outer_boundary_of_my_tree, dinfo.getPosRootDomain(),
                 pa);
        }
        const S32 n_image = shift_image_domain.size();
        F64 * mass_tree = new F64[n_loc];
        F64 r_search_sq = epj_org_[0].getRSearch() * epj_org_[0].getRSearch();
        for(S32 i=0; i<n_loc; i++){
            mass_tree[i] = 0.0;
            const S32 n_ep_recv = n_ep_recv_disp_[n_proc];
            for(S32 j=0; j<n_ep_recv; j++){
                F64vec dr = epj_org_[i].getPos() - epj_recv_[j].getPos();
                if( dr*dr < r_search_sq){
                    mass_tree[i] += epj_recv_[j].getCharge();
                }
            }
            for(S32 j=0; j<n_loc; j++){
                F64vec dr = epj_org_[i].getPos() - epj_org_[j].getPos();
                if( dr*dr < r_search_sq){
                    mass_tree[i] += epj_org_[j].getCharge();
                }
            }
            S32 n_sp = n_sp_recv_disp_[n_proc];
            for(S32 j=0; j<n_sp; j++){
                F64vec dr = epj_org_[i].getPos() - spj_recv_[j].getPos();
                if( dr*dr < r_search_sq){
                    mass_tree[i] += spj_recv_[j].getCharge();
                }
            }
        }
        F64 * mass_direct = new F64[n_loc];
        for(S32 i=0; i<n_loc; i++){
            mass_direct[i] = 0.0;
            const S32 n_ep = n_ptcl_disp[n_proc];
            for(S32 j=0; j<n_ep; j++){
                for(S32 ii=0; ii<n_image; ii++){
                    F64vec dr = epj_org_[i].getPos() - ep_tmp[j].getPos() - shift_image_domain[ii];
                    if( dr*dr < r_search_sq){
                        mass_direct[i] += ep_tmp[j].getCharge();
                    }
                }
            }
        }
        bool err = false;
        S32 n_err = 0;
        for(S32 i=0; i<n_loc; i++){
            if( fabs(mass_tree[i]-mass_direct[i]) > tolerance){
                err = true;
                fout<<"mass_tree[i-1]="<<mass_tree[i-1]<<" mass_direct[i-1]="<<mass_direct[i-1]<<std::endl;
                fout<<"mass_tree[i]="<<mass_tree[i]<<" mass_direct[i]="<<mass_direct[i]<<std::endl;
                fout<<"mass_tree[i+1]="<<mass_tree[i+1]<<" mass_direct[i+1]="<<mass_direct[i+1]<<std::endl;
                fout<<"fabs(mass_tree[i]-mass_direct[i])="<<fabs(mass_tree[i]-mass_direct[i])<<std::endl;
                fout<<std::endl;
                n_err++;
            }
        }
        if(err == false) fout<<"CheckExchangeLocalEssentialTreeDummy(LONGCUTOFF): PASS"<<std::endl;
        else {
            fout<<"CheckExchangeLocalEssentialTreeDummy(LONGCUTOFF): NOT PASS. BUT It dose not mean error."<<std::endl;
            fout<<" If n_err/n_loc is small enough. Exchange is acceseptable."<<std::endl;
        }
        fout<<"tolerance= "<<tolerance<<std::endl;
        fout<<"this->n_ep_recv_disp_[n_proc]= "<<this->n_ep_recv_disp_[n_proc]<<std::endl;
        fout<<"n_err="<<n_err<<" n_err/n_loc="<<(double)(n_err)/n_loc<<std::endl;
        delete [] ep_tmp;
        delete [] n_ptcl;
        delete [] n_ptcl_disp;
        delete [] mass_tree;
        delete [] mass_direct;
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL 	    
    }

    
    /////////////////
    /// FOR SHORT ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkExchangeLocalEssentialTreeImpl(TagForceShort,
                                        const DomainInfo & dinfo,
                                        const F64 tolerance,
                                        std::ostream & fout){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL 
        const S32 n_proc = Comm::getNumberOfProc();
        Tepj * ep_tmp;
        S32 * n_ptcl;
        S32 * n_ptcl_disp;
        const S32 n_loc = this->n_loc_tot_;
        AllGatherParticle(ep_tmp, n_ptcl, n_ptcl_disp, this->epj_org_.getPointer(), n_loc);
        S32 * n_recv_direct = new S32[n_proc];
        ReallocatableArray<F64vec> pos_direct;
        bool pa[DIMENSION];
        ReallocatableArray<F64vec> shift_image_domain;
        dinfo.getPeriodicAxis(pa);
        F64ort * outer_boundary_of_tree = new F64ort[n_proc];
        F64ort outer_boundary_of_my_tree = this->tc_loc_[0].mom_.vertex_out_;
        Comm::allGather(&outer_boundary_of_my_tree, 1, outer_boundary_of_tree);
        S32 n_cum = 0;
        bool err = false;
        for(S32 ib=0; ib<n_proc; ib++){
            S32 n_recv_per_proc = 0;
            n_recv_direct[ib] = 0;
            shift_image_domain.clearSize();
            if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                shift_image_domain.push_back( F64vec(0.0) );
            }
            else{
                shift_image_domain.reserve(5*5*5);
                CalcNumberAndShiftOfImageDomain
                    (shift_image_domain, dinfo.getPosRootDomain().getFullLength(),
                     outer_boundary_of_my_tree, outer_boundary_of_tree[ib], pa);
            }
            const S32 n_image_per_proc = shift_image_domain.size();
            //////////
            // this scope (for loop) is the only differenc from other check functions
            checkExchangeLocalEssentialTreeForShortImpl
                (typename TSM::search_type(), dinfo, ep_tmp, n_ptcl_disp[ib], 
                 n_ptcl_disp[ib+1],  ib, n_image_per_proc, shift_image_domain, 
                 n_recv_per_proc, pos_direct);
            // this scope (for loop) is the only differenc from other check functions
            //////////
	    
            bool err_tmp = false;
            n_recv_direct[ib] = n_recv_per_proc;
            if(this->n_ep_recv_[ib] > n_recv_direct[ib]){
                fout<<"It dose not mean error"<<std::endl;
                fout<<"Comm::getRank()= "<<Comm::getRank()<<" ib= "<<ib<<std::endl;
                fout<<"this->n_ep_recv_[ib]= "<<this->n_ep_recv_[ib]<<" n_recv_direct[ib]="<<n_recv_direct[ib]<<std::endl;
            }
            else if(this->n_ep_recv_[ib] < n_recv_direct[ib]){
                fout<<"error 1: this->n_ep_recv_[ib] is too small"<<std::endl;
                fout<<"Comm::getRank()= "<<Comm::getRank()<<" ib= "<<ib<<std::endl;
                fout<<"this->n_ep_recv_[ib]= "<<this->n_ep_recv_[ib]<<" n_recv_direct[ib]="<<n_recv_direct[ib]<<std::endl;
                err_tmp = true;
                err = true;
            }

            std::vector<F64vec> pos_tree_per_proc;
            pos_tree_per_proc.clear();
            pos_tree_per_proc.reserve(10000);
            for(S32 ip=this->n_ep_recv_disp_[ib]; ip<this->n_ep_recv_disp_[ib+1]; ip++){
                pos_tree_per_proc.push_back(this->epj_recv_[ip].getPos());
            }
            std::sort(pos_tree_per_proc.begin(), pos_tree_per_proc.end(), LessOPForVecX());
            std::vector<F64vec> pos_direct_per_proc;
            pos_direct_per_proc.clear();
            pos_direct_per_proc.reserve(1000);
            for(S32 ip=0; ip<n_recv_direct[ib]; ip++){
                pos_direct_per_proc.push_back(pos_direct[ip+n_cum]);
            }
            std::sort(pos_direct_per_proc.begin(), pos_direct_per_proc.end(), LessOPForVecX());
            const S32 n_tmp = std::min(pos_direct_per_proc.size(), pos_tree_per_proc.size());
            if(err_tmp == false){
                // check if epj by direct is a subset of epj by tree
                S32 jp1_loc = 0;
                bool found = false;
                for(S32 jp0=0; jp0<pos_direct_per_proc.size(); jp0++){
                    F64vec ref = pos_direct_per_proc[jp0];
                    for(S32 jp1=jp1_loc; jp1<pos_tree_per_proc.size(); jp1++){
                        if(pos_direct_per_proc[jp0] == pos_tree_per_proc[jp1]){
                            jp1_loc = jp1 + 1;
                            found = true;
                            break;
                        }
                    }
                    if(found){ found = false; }
                    else{ 
                        fout<<"error 2: ep by tree search is not a super set of ep by direct search"<<std::endl;
                        fout<<"Comm::getRank()= "<<Comm::getRank()<<" ib= "<<ib<<std::endl;
                        fout<<"jp1_loc="<<jp1_loc<<std::endl;
                        fout<<"pos_direct_per_proc[jp0]= "<<pos_direct_per_proc[jp0]<<" pos_tree_per_proc[jp1_loc]="<<pos_tree_per_proc[jp1_loc]<<std::endl;
                        err = true;    
                        break;
                    }
                }

            }
            n_cum += n_recv_direct[ib];
        }
        if(err == false) fout<<"CheckExchangeLocalEssentialTreeDummy(SHORT): PASS"<<std::endl;
        else fout<<"CheckExchangeLocalEssentialTreeDummy(SHORT): FAIL"<<std::endl;
        fout<<"this->n_ep_recv_disp_[n_proc]= "<<this->n_ep_recv_disp_[n_proc]<<std::endl;
        delete [] ep_tmp;
        delete [] n_ptcl;
        delete [] n_ptcl_disp;
        delete [] n_recv_direct;
        delete [] outer_boundary_of_tree;
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkExchangeLocalEssentialTreeForShortImpl
    (TagSearchShortScatter,
     const DomainInfo & dinfo,
     const Tep2 ep_tmp[],
     const S32 jp_head,
     const S32 jp_tail,
     const S32 rank_target,
     const S32 n_image_per_proc,
     const ReallocatableArray<F64vec> & shift_image_domain,
     S32 & n_recv_per_proc,
     ReallocatableArray<F64vec> & pos_direct){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL 
        const S32 my_rank = Comm::getRank();
        //for(S32 jp=n_ptcl_disp[ib]; jp<n_ptcl_disp[ib+1]; jp++){
        for(S32 jp=jp_head; jp<jp_tail; jp++){
            const F64 r_search_sq = ep_tmp[jp].getRSearch() * ep_tmp[jp].getRSearch();
            for(S32 k=0; k<n_image_per_proc; k++){
                if(rank_target == my_rank && k == 0) continue;
                F64 r_sq = dinfo.getPosDomain(my_rank).getDistanceMinSQ(ep_tmp[jp].getPos()+shift_image_domain[k]);
                if(r_sq < r_search_sq){
                    pos_direct.push_back(ep_tmp[jp].getPos()+shift_image_domain[k]);
                    n_recv_per_proc++;
                }
            }
        }
#endif	    
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkExchangeLocalEssentialTreeForShortImpl
    (TagSearchShortGather,
     const DomainInfo & dinfo,
     const Tep2 ep_tmp[],
     const S32 jp_head,
     const S32 jp_tail,
     const S32 rank_target,
     const S32 n_image_per_proc,
     const ReallocatableArray<F64vec> & shift_image_domain,
     S32 & n_recv_per_proc,
     ReallocatableArray<F64vec> & pos_direct){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const S32 my_rank = Comm::getRank();
        //for(S32 jp=n_ptcl_disp[ib]; jp<n_ptcl_disp[ib+1]; jp++){
        for(S32 jp=jp_head; jp<jp_tail; jp++){
            const F64vec pos_j = ep_tmp[jp].getPos();
            for(S32 ii=0; ii<n_image_per_proc; ii++){
                if(rank_target == my_rank && ii == 0) continue;
                const F64vec shift = shift_image_domain[ii];
                for(S32 ip=0; ip<this->n_loc_tot_; ip++){
                    const F64vec pos_i = this->epi_org_[ip].getPos();
                    const F64 len_sq_i = this->epi_org_[ip].getRSearch() * this->epi_org_[ip].getRSearch();
                    const F64 dis_sq_pp = pos_i.getDistanceSQ(pos_j+shift);
                    if(dis_sq_pp <= len_sq_i){
                        pos_direct.push_back(pos_j+shift);
                        n_recv_per_proc++;
                        break;
                    }
                }
            }
        }
#endif	    
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkExchangeLocalEssentialTreeForShortImpl
    (TagSearchShortSymmetry,
     const DomainInfo & dinfo,
     const Tep2 ep_tmp[],
     const S32 jp_head,
     const S32 jp_tail,
     const S32 rank_target,
     const S32 n_image_per_proc,
     const ReallocatableArray<F64vec> & shift_image_domain,
     S32 & n_recv_per_proc,
     ReallocatableArray<F64vec> & pos_direct){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const S32 my_rank = Comm::getRank();
        for(S32 jp=jp_head; jp<jp_tail; jp++){
            const F64vec pos_j = ep_tmp[jp].getPos();
            const F64 len_sq_j = ep_tmp[jp].getRSearch() * ep_tmp[jp].getRSearch();
            for(S32 ii=0; ii<n_image_per_proc; ii++){
                if(rank_target == my_rank && ii == 0) continue;
                const F64vec shift = shift_image_domain[ii];
                const F64 dis_sq = dinfo.getPosDomain(my_rank).getDistanceMinSQ(pos_j+shift);
                for(S32 ip=0; ip<this->n_loc_tot_; ip++){
                    const F64vec pos_i = this->epi_org_[ip].getPos();
                    const F64 len_sq_i = this->epi_org_[ip].getRSearch() * this->epi_org_[ip].getRSearch();
                    const F64 dis_sq_pp = pos_i.getDistanceSQ(pos_j+shift);
                    if(dis_sq <= len_sq_j || dis_sq_pp <= len_sq_i){
                        pos_direct.push_back(pos_j+shift);
                        n_recv_per_proc++;
                        break;
                    }
                }
            }
        }
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL	    
    }


    /////////////////
    // CALC MOMENT //
    // LOCAL 
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentLocalTree(const F64 tolerance, std::ostream & fout){
	checkCalcMomentLocalTreeImpl(typename TSM::search_type(), tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentLocalTreeImpl(TagSearchLong, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentLongLocalTree(epj_sorted_.getPointer(), n_loc_tot_,
				     tc_loc_.getPointer(), tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentLocalTreeImpl(TagSearchLongCutoff, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentOutOnly(epj_sorted_.getPointer(), n_loc_tot_,
			       tc_loc_.getPointer(), tolerance, fout);	
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentLocalTreeImpl(TagSearchShortScatter, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentShortInAndOut(epj_sorted_.getPointer(), n_loc_tot_,
				     tc_loc_.getPointer(), tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentLocalTreeImpl(TagSearchShortGather, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentShortInAndOut(epi_sorted_.getPointer(), n_loc_tot_,
				     tc_loc_.getPointer(), tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentLocalTreeImpl(TagSearchShortSymmetry, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentShortInAndOut(epi_sorted_.getPointer(), n_loc_tot_,
				     tc_loc_.getPointer(), tolerance, fout);
	CheckCalcMomentShortInAndOut(epj_sorted_.getPointer(), n_loc_tot_,
				     tc_loc_.getPointer(), tolerance, fout);
    }

    /////////////////
    // CALC MOMENT //
    // GLOBAL
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentGlobalTree(const F64 tolerance, std::ostream & fout){
	checkCalcMomentGlobalTreeImpl(typename TSM::search_type(), tolerance, fout);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentGlobalTreeImpl(TagSearchLong, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentLongGlobalTree(epj_sorted_.getPointer(), spj_sorted_.getPointer(), 
				      n_glb_tot_, tc_glb_.getPointer(),
				      tp_glb_.getPointer(), tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentGlobalTreeImpl(TagSearchLongCutoff, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentLongGlobalTree(epj_sorted_.getPointer(), spj_sorted_.getPointer(), 
				      n_glb_tot_, tc_glb_.getPointer(),
				      tp_glb_.getPointer(), tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentGlobalTreeImpl(TagSearchShortScatter, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentShortInAndOut(epj_sorted_.getPointer(), n_glb_tot_, 
				     tc_glb_.getPointer(),     tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentGlobalTreeImpl(TagSearchShortGather, const F64 tolerance, std::ostream & fout){
	CheckCalcMomentShortInOnly(epj_sorted_.getPointer(), n_glb_tot_, 
				   tc_glb_.getPointer(),     tolerance, fout);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkCalcMomentGlobalTreeImpl(TagSearchShortSymmetry, const F64 tolerance,
				  std::ostream & fout){
	CheckCalcMomentShortInAndOut(epj_sorted_.getPointer(), n_glb_tot_,
				     tc_glb_.getPointer(),     tolerance, fout);
    }

    //////////////////
    // MAKE IPGROUP //
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeInteractionListImpl(TagSearchLong,
                                 const DomainInfo & dinfo,
                                 const S32 adr_ipg,
                                 const S32 ith,
                                 const F64 tolerance,
                                 std::ostream & fout){
        F64 mass_cm = 0.0;
        F64vec pos_cm = 0.0;
        const S32 n_ep = epj_for_force_[ith].size();
        const S32 n_sp = spj_for_force_[ith].size();
        for(S32 i=0; i<n_ep; i++){
            mass_cm += epj_for_force_[ith][i].getCharge();
            pos_cm += epj_for_force_[ith][i].getCharge() * epj_for_force_[ith][i].getPos();
        }
        for(S32 i=0; i<n_sp; i++){
            mass_cm += spj_for_force_[ith][i].getCharge();
            pos_cm += spj_for_force_[ith][i].getCharge() * spj_for_force_[ith][i].getPos();
        }
        pos_cm /= mass_cm;
        F64 dm = std::abs(tc_glb_[0].mom_.mass - mass_cm);
        F64vec dx = tc_glb_[0].mom_.pos - pos_cm;
        F64 dx_max = dx.applyEach( Abs<F64>() ).getMax();
        if( dm >= tolerance || dx_max >= tolerance){
            fout<<"CheckMakeInteractionList: FAIL"<<std::endl;
            fout<<"mass_cm="<<mass_cm<<", tc_glb_[0].mom_.mass="<<tc_glb_[0].mom_.mass<<std::endl;
            fout<<"pos_cm="<<pos_cm<<", tc_glb_[0].mom_.pos="<<tc_glb_[0].mom_.pos<<std::endl;
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeInteractionListImpl(TagSearchLongCutoff,
                                 const DomainInfo & dinfo,
                                 const S32 adr_ipg,
                                 const S32 ith,
                                 const F64 tolerance,
                                 std::ostream & fout){
        fout<<"CheckMakeInteractionList(SEARCH_MODE_LONG_CUTOFF)"<<std::endl;
        fout<<"Under construction"<<std::endl;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeInteractionListImpl(TagSearchShortScatter,
                                 const DomainInfo & dinfo,
                                 const S32 adr_ipg,
                                 const S32 ith,
                                 const F64 tolerance,
                                 std::ostream & fout){
        //for SCATTER
        bool err = false;
	fout<<"CheckMakeInteractionListDummy dafult (SCATTER)"<<std::endl;
	const S32 adr_ptcl = ipg_[adr_ipg].adr_ptcl_;
	const S32 n_i = ipg_[adr_ipg].n_ptcl_;
	const S32 n_ep = epj_for_force_[ith].size();
	S32 * n_ngb = new S32[n_i];
	ReallocatableArray<F64vec> pos_tree;
	pos_tree.clearSize();
	ReallocatableArray<F64vec> pos_direct;
	pos_direct.clearSize();
	for(S32 i=0; i<n_i; i++){
	    n_ngb[i] = 0;
	    for(S32 j=0; j<n_ep; j++){
		const F64vec dr = epi_sorted_[i+adr_ptcl].getPos() - epj_for_force_[ith][j].getPos();
		const F64 r_crit_sq = epj_for_force_[ith][j].getRSearch() * epj_for_force_[ith][j].getRSearch();
		if( dr*dr <= r_crit_sq) n_ngb[i]++;
	    }
	}
	Tepj * epj_tmp;
	S32 n_tot_tmp;
	bool pa[DIMENSION];
	dinfo.getPeriodicAxis(pa);
	AllGatherParticle(epj_tmp, n_tot_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa);
	S32 * n_ngb_tmp = new S32[n_i];
	for(S32 i=0; i<n_i; i++){
	    n_ngb_tmp[i] = 0;
	    for(S32 j=0; j<n_tot_tmp; j++){
		const F64vec dr = epi_sorted_[i+adr_ptcl].getPos() - epj_tmp[j].getPos();
		const F64 r_crit_sq = epj_tmp[j].getRSearch() * epj_tmp[j].getRSearch();
		if( dr*dr <= r_crit_sq) n_ngb_tmp[i]++;
	    }
	    if(n_ngb_tmp[i] != n_ngb[i]){
		fout<<"CheckMakeInteractionList(SCATTER): FAIL"<<std::endl;
		fout<<"n_ngb[i]="<<n_ngb[i]<<" n_ngb_direct[i]="<<n_ngb_tmp[i]<<std::endl;
		err = true;
		pos_tree.clearSize();
		for(S32 j=0; j<epj_for_force_[ith].size(); j++){
		    const F64vec dr = epi_sorted_[i+adr_ptcl].getPos() - epj_for_force_[ith][j].getPos();
		    const F64 r_crit_sq = epj_for_force_[ith][j].getRSearch() * epj_for_force_[ith][j].getRSearch();
		    if( dr*dr <= r_crit_sq){
			pos_tree.push_back(epj_for_force_[ith][j].getPos());
		    }
		}
		pos_direct.clearSize();
		for(S32 j=0; j<n_tot_tmp; j++){
		    const F64vec dr = epi_sorted_[i+adr_ptcl].getPos() - epj_tmp[j].getPos();
		    const F64 r_crit_sq = epj_tmp[j].getRSearch() * epj_tmp[j].getRSearch();
		    if( dr*dr <= r_crit_sq){
			pos_direct.push_back(epj_tmp[j].getPos());
		    }
		}
		std::sort(pos_tree.getPointer(),   pos_tree.getPointer( pos_tree.size()+1 ), LessOPForVecX());
		std::sort(pos_direct.getPointer(), pos_direct.getPointer(  pos_direct.size()+1 ), LessOPForVecX());
		for(S32 ii=0; ii<n_ngb_tmp[ii]; ii++){
		    if(pos_tree[ii] != pos_direct[ii] || pos_tree.size() <= 0 || pos_direct.size() <= 0){
			fout<<"ii="<<ii<<std::endl;
			fout<<"pos_tree[ii-1]="<<pos_tree[ii-1]<<" pos_direct[ii-1]="<<pos_direct[ii-1]<<std::endl;
			fout<<"pos_tree[ii]="<<pos_tree[ii]<<" pos_direct[ii]="<<pos_direct[ii]<<std::endl;
			fout<<"pos_tree[ii+1]="<<pos_tree[ii+1]<<" pos_direct[ii+1]="<<pos_direct[ii+1]<<std::endl;
		    }
		}
	    }
	}
	if(!err){
	    fout<<"CheckMakeInteractionList(SCATTER): PASS"<<std::endl;
	}
	delete [] epj_tmp;
	delete [] n_ngb_tmp;
	delete [] n_ngb;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeInteractionListImpl(TagSearchShortGather,
				 const DomainInfo & dinfo,
				 const S32 adr_ipg,
				 const S32 ith,
				 const F64 tolerance,
				 std::ostream & fout){
	std::cout<<"under construction"<<std::endl;
	//for GATHER
	bool err = false;
	//fout<<"CheckMakeInteractionListDummy (GATHER)"<<std::endl;
	const S32 adr_ptcl = ipg_[adr_ipg].adr_ptcl_;
	const S32 n_i = ipg_[adr_ipg].n_ptcl_;
	const S32 n_ep = epj_for_force_[ith].size();
	S32 * n_ngb = new S32[n_i];
	for(S32 i=0; i<n_i; i++){
	    n_ngb[i] = 0;
	    const F64 r_crit_sq = epi_sorted_[i+adr_ptcl].getRSearch() * epi_sorted_[i+adr_ptcl].getRSearch();
	    for(S32 j=0; j<n_ep; j++){
		const F64vec dr = epi_sorted_[i+adr_ptcl].getPos() - epj_for_force_[ith][j].getPos();
		if( dr*dr <= r_crit_sq) n_ngb[i]++;
	    }
	}
	Tepj * epj_tmp;
	S32 n_tot_tmp;
	bool pa[DIMENSION];
	dinfo.getPeriodicAxis(pa);
	AllGatherParticle(epj_tmp, n_tot_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa);
	S32 * n_ngb_tmp = new S32[n_i];
	for(S32 i=0; i<n_i; i++){
	    n_ngb_tmp[i] = 0;
	    const F64 r_crit_sq = epi_sorted_[i+adr_ptcl].getRSearch() * epi_sorted_[i+adr_ptcl].getRSearch();
	    for(S32 j=0; j<n_tot_tmp; j++){
		const F64vec dr = epi_sorted_[i+adr_ptcl].getPos() - epj_tmp[j].getPos();
		if( dr*dr <= r_crit_sq) n_ngb_tmp[i]++;
	    }
	    if(n_ngb_tmp[i] != n_ngb[i]){
		fout<<"CheckMakeInteractionList(GATHER): FAIL"<<std::endl;
		fout<<"n_ngb[i]="<<n_ngb[i]<<" n_ngb_direct[i]="<<n_ngb_tmp[i]<<std::endl;
		err = true;
	    }
	}
	if(!err){
	    fout<<"CheckMakeInteractionList(GATHER): PASS"<<std::endl;
	}
	delete [] epj_tmp;
	delete [] n_ngb_tmp;
	delete [] n_ngb;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeInteractionListImpl(TagSearchShortSymmetry,
                                 const DomainInfo & dinfo,
                                 const S32 adr_ipg,
                                 const S32 ith,
                                 const F64 tolerance,
                                 std::ostream & fout){
        std::cout<<"under construction"<<std::endl;
        //for SYMMETRY
        bool err = false;
        const S32 adr_ptcl = ipg_[adr_ipg].adr_ptcl_;
        const S32 n_i = ipg_[adr_ipg].n_ptcl_;
        const S32 n_ep = epj_for_force_[ith].size();
        S32 * n_ngb = new S32[n_i];
        for(S32 i=0; i<n_i; i++){
            n_ngb[i] = 0;
            const F64 r_crit_sq = epi_sorted_[i+adr_ptcl].getRSearch() * epi_sorted_[i+adr_ptcl].getRSearch();
            for(S32 j=0; j<n_ep; j++){
                const F64vec dr = epi_sorted_[i+adr_ptcl].getPos() - epj_for_force_[ith][j].getPos();
                if( dr*dr <= r_crit_sq) n_ngb[i]++;
            }
        }
        Tepj * epj_tmp;
        S32 n_tot_tmp;
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        AllGatherParticle(epj_tmp, n_tot_tmp, epj_org_.getPointer(),
                          n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa);
        S32 * n_ngb_tmp = new S32[n_i];
        for(S32 i=0; i<n_i; i++){
            n_ngb_tmp[i] = 0;
            const F64 r_crit_sq = epi_sorted_[i+adr_ptcl].getRSearch() * epi_sorted_[i+adr_ptcl].getRSearch();
            for(S32 j=0; j<n_tot_tmp; j++){
                const F64vec dr = epi_sorted_[i+adr_ptcl].getPos() - epj_tmp[j].getPos();
                if( dr*dr <= r_crit_sq) n_ngb_tmp[i]++;
            }
            if(n_ngb_tmp[i] != n_ngb[i]){
                fout<<"CheckMakeInteractionList(SYMMETRY): FAIL"<<std::endl;
                fout<<"n_ngb[i]="<<n_ngb[i]<<" n_ngb_direct[i]="<<n_ngb_tmp[i]<<std::endl;
                err = true;
            }
        }
        if(!err){
            fout<<"CheckMakeInteractionList(SYMMETRY): PASS"<<std::endl;
        }
        delete [] epj_tmp;
        delete [] n_ngb_tmp;
        delete [] n_ngb;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkMakeIPGroup(const F64 tolerance,  std::ostream & fout){
        fout<<"checkMakeIPGroup"<<std::endl;
        bool err_overflow = false;
        bool err_wrongadr = false;
        bool err_outofbox = false;
        S32 n_cnt = 0;
        S32 n_ipg = ipg_.size();
        for(S32 i=0; i<n_ipg; i++){
            if(ipg_[i].n_ptcl_ > n_group_limit_){
                ipg_[i].dump();
                err_overflow = true;
            }
            if( n_cnt != ipg_[i].adr_ptcl_){
                err_wrongadr = true;
            }
            for(S32 j=0; j<ipg_[i].n_ptcl_; j++){
                if( ipg_[i].vertex_.getDistanceMinSQ( epi_sorted_[(ipg_[i].adr_ptcl_)+j].getPos() )  > tolerance * tolerance){
                    err_outofbox = true;
                }
            }
            n_cnt += ipg_[i].n_ptcl_;
        }
        if(err_wrongadr || err_overflow || err_outofbox){
            fout<<"checkMakeIPGroup: FAIL err_overflow="<<err_overflow
                <<"err_wrongadr="<<err_wrongadr
                <<"err_outofbox="<<err_outofbox<<std::endl;
        }
        else{
            fout<<"checkMakeIPGroup: PASS"<<std::endl;
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_compare>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    checkForce(Tfunc_ep_ep pfunc_ep_ep,
               Tfunc_compare func_compare,
               const DomainInfo & dinfo,
               std::ostream & fout){
        Tforce * force_direct = new Tforce[n_loc_tot_];
        for(S32 i=0; i<n_loc_tot_; i++) force_direct[i].clear();
        calcForceDirect(pfunc_ep_ep, force_direct, dinfo, true);
        func_compare(force_org_.getPointer(), force_direct, n_loc_tot_, fout );
        delete [] force_direct;
    }
} // end of namespace
