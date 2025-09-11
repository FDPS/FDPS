#pragma once

namespace ParticleSimulator {
// from bottom to top
template <class Ttc, class Tepj, typename Tmorton_key>
inline void CalcMomentLongLocalTree(S32 adr_tc_level_partition[], Ttc tc[], Tepj epj[], const S32 lev_max, const S32 n_leaf_limit,
                                    const S32 lev_leaf_limit, const Tmorton_key& morton_key) {
    for (S32 i = lev_max; i >= 0; --i) {
        const S32 head = adr_tc_level_partition[i];
        const S32 next = adr_tc_level_partition[i + 1];
        const F64 len = morton_key.getCorrespondingFullLength(i);
        PS_OMP_PARALLEL_FOR
        for (S32 j = head; j < next; j++) {
            Ttc* tc_tmp = tc + j;
            const int n_tmp = tc_tmp->n_ptcl_;
            tc_tmp->clearMoment();
            tc_tmp->geo_.setSize(len);
            if (n_tmp == 0)
                continue;
            else if (tc_tmp->isLeafPMM(n_leaf_limit, lev_leaf_limit)) {
                const S32 adr = tc_tmp->adr_ptcl_;
                tc_tmp->geo_.addFirstEpj(epj[adr]);
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                    tc_tmp->geo_.accumulateAtLeaf(epj[k]);
                }
                tc_tmp->mom_.set();
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
                }
            } else {
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ != 0) {
                        tc_tmp->geo_.copyBox(tc_tmp_tmp->geo_);
                        break;
                    }
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate(tc_tmp_tmp->mom_);
                    tc_tmp->geo_.accumulate(tc_tmp_tmp->geo_);
                }
                tc_tmp->mom_.set();
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate2(tc_tmp_tmp->mom_);
                }
            }
        }
    }
    // tc[0].dump();
}

// for short search
template <class Ttc, class Tepj>
inline void CalcMoment(S32 adr_tc_level_partition[], Ttc tc[], Tepj epj[], const S32 lev_max, const S32 n_leaf_limit) {
    for (S32 i = lev_max; i >= 0; --i) {
        const S32 head = adr_tc_level_partition[i];
        const S32 next = adr_tc_level_partition[i + 1];
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for (S32 j = head; j < next; j++) {
            Ttc* tc_tmp = tc + j;
            const int n_tmp = tc_tmp->n_ptcl_;
            // tc_tmp->mom_.init();
            tc_tmp->clearMoment();
            if (n_tmp == 0)
                continue;
            else if (tc_tmp->isLeaf(n_leaf_limit)) {
                const S32 adr = tc_tmp->adr_ptcl_;
                tc_tmp->geo_.addFirstEpj(epj[adr]);
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    // tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                    tc_tmp->geo_.accumulateAtLeaf(epj[k]);
                }
                // tc_tmp->mom_.set();
                // for(S32 k=adr; k<adr+n_tmp; k++){
                //     tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
                // }
            } else {
                /*
                for(S32 k=0; k<N_CHILDREN; k++){
                    Ttc * tc_tmp_tmp = tc + ((tc_tmp->adr_tc_)+k);
                    if(tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate( tc_tmp_tmp->mom_ );
                }
                tc_tmp->mom_.set();
                for(S32 k=0; k<N_CHILDREN; k++){
                    Ttc * tc_tmp_tmp = tc + ((tc_tmp->adr_tc_)+k);
                    if(tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate2( tc_tmp_tmp->mom_ );
                }
                */
                S32 k2 = 0;
                for (k2 = 0; k2 < N_CHILDREN; k2++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k2);
                    if (tc_tmp_tmp->n_ptcl_ != 0) {
                        tc_tmp->geo_ = tc_tmp_tmp->geo_;
                        break;
                    }
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ != 0) {
                        tc_tmp->geo_.copyBox(tc_tmp_tmp->geo_);
                        break;
                    }
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->geo_.accumulate(tc_tmp_tmp->geo_);
                }
            }
        }
    }
}

template <class Ttc, class Tepj>
inline void CalcMomentST(S32 adr_tc_level_partition[], Ttc tc[], Tepj epj[], const S32 lev_max, const S32 n_leaf_limit) {
    for (S32 i = lev_max; i >= 0; --i) {
        const S32 head = adr_tc_level_partition[i];
        const S32 next = adr_tc_level_partition[i + 1];
        for (S32 j = head; j < next; j++) {
            Ttc* tc_tmp = tc + j;
            const int n_tmp = tc_tmp->n_ptcl_;
            // tc_tmp->mom_.init();
            tc_tmp->clearMoment();
            if (n_tmp == 0)
                continue;
            else if (tc_tmp->isLeaf(n_leaf_limit)) {
                const S32 adr = tc_tmp->adr_ptcl_;
                tc_tmp->geo_.addFirstEpj(epj[adr]);
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                }
                tc_tmp->mom_.set();
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
                }
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    tc_tmp->geo_.accumulateAtLeaf(epj[k]);
                }
            } else {
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate(tc_tmp_tmp->mom_);
                }
                tc_tmp->mom_.set();
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate2(tc_tmp_tmp->mom_);
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ != 0) {
                        tc_tmp->geo_.copyBox(tc_tmp_tmp->geo_);
                        break;
                    }
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->geo_.accumulate(tc_tmp_tmp->geo_);
                }
            }
        }
    }
}

#if 1
// PMM
template <class Ttc, class Tepj, class Tspj, typename Tmorton_key>
inline void CalcMomentLongGlobalTree(S32 adr_tc_level_partition[], Ttc tc[], TreeParticle tp[], Tepj epj[], Tspj spj[], const S32 lev_max,
                                     const S32 n_leaf_limit,
                                     const S32 lev_leaf_limit,  // added one arg
                                     const Tmorton_key& morton_key) {
    for (S32 i = lev_max; i >= 0; --i) {
        const S32 head = adr_tc_level_partition[i];
        const S32 next = adr_tc_level_partition[i + 1];
        const F64 len = morton_key.getCorrespondingFullLength(i);
        /*
        std::cout<<"i= "<<i<<" lev_max= "<<lev_max<<std::endl;
        std::cout<<"head= "<<head
                 <<" next= "<<next
                 <<" len= "<<len
                 <<std::endl;
        */
        PS_OMP_PARALLEL_FOR
        for (S32 j = head; j < next; j++) {
            Ttc* tc_tmp = tc + j;
            const int n_tmp = tc_tmp->n_ptcl_;
            // tc_tmp->mom_.init();
            tc_tmp->clearMoment();
            tc_tmp->geo_.setSize(len);
            // std::cerr<<"j= "<<j<<std::endl;
            /*
            if(tc_tmp->adr_ptcl_ > 1000000){
                tc_tmp->dump(std::cerr);
                assert(0);
            }
            */
            // tc_tmp->dump(std::cerr);
            if (n_tmp == 0) {
                // std::cerr<<"check A"<<std::endl;
                continue;
            } else if (tc_tmp->isLeafPMM(n_leaf_limit, lev_leaf_limit)) {  // added one arg
                // std::cerr<<"check B"<<std::endl;
                const S32 adr = tc_tmp->adr_ptcl_;
                if (GetMSB(tp[adr].adr_ptcl_) == 0) {
                    const U32 adr_ep = tp[adr].adr_ptcl_;
                    tc_tmp->geo_.addFirstEpj(epj[adr_ep]);
                    // std::cout<<"adr_ep= "<<adr_ep
                    //          <<" epj[adr_ep].getPos()= "
                    //          <<epj[adr_ep].getPos()
                    //          <<std::endl;
                } else {
                    const U32 adr_sp = ClearMSB(tp[adr].adr_ptcl_);
                    tc_tmp->geo_.addFirstSpj(spj[adr_sp]);
                    // std::cout<<"adr_sp= "<<adr_sp
                    //          <<" spj[adr_sp].getPos()= "
                    //          <<spj[adr_sp].getPos()
                    //          <<std::endl;
                }
                // if(Comm::getRank()==0){
                //     std::cout<<"A) j= "<<j<<std::endl;
                //     tc_tmp->geo_.dump(std::cout);
                // }
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    if (GetMSB(tp[k].adr_ptcl_) == 0) {
                        const U32 adr_ep = tp[k].adr_ptcl_;
                        tc_tmp->mom_.accumulateAtLeaf(epj[adr_ep]);
                        tc_tmp->geo_.accumulateAtLeaf(epj[adr_ep]);
                    } else {
                        const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                        tc_tmp->mom_.accumulate(spj[adr_sp].convertToMoment());
                    }
                }
                // if(Comm::getRank()==0){
                //     std::cout<<"B) j= "<<j<<std::endl;
                //     tc_tmp->geo_.dump(std::cout);
                // }
                tc_tmp->mom_.set();
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    if (GetMSB(tp[k].adr_ptcl_) == 0) {
                        const U32 adr_ep = tp[k].adr_ptcl_;
                        tc_tmp->mom_.accumulateAtLeaf2(epj[adr_ep]);
                    } else {
                        const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                        tc_tmp->mom_.accumulate2(spj[adr_sp].convertToMoment());
                    }
                }
                // std::cerr<<"check B2"<<std::endl;
            } else {
                // std::cerr<<"check C"<<std::endl;
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ != 0) {
                        tc_tmp->geo_.copyBox(tc_tmp_tmp->geo_);
                        break;
                    }
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate(tc_tmp_tmp->mom_);
                    tc_tmp->geo_.accumulate(tc_tmp_tmp->geo_);
                }
                tc_tmp->mom_.set();
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate2(tc_tmp_tmp->mom_);
                }
                // std::cerr<<"check C2"<<std::endl;
            }
        }
        // std::cerr<<"check D"<<std::endl;
    }
}
#else
// original (not siginificant different)
template <class Ttc, class Tepj, class Tspj, typename Tmorton_key>
inline void CalcMomentLongGlobalTree(S32 adr_tc_level_partition[], Ttc tc[], TreeParticle tp[], Tepj epj[], Tspj spj[], const S32 lev_max,
                                     const S32 n_leaf_limit, const Tmorton_key& morton_key) {
    for (S32 i = lev_max; i >= 0; --i) {
        const S32 head = adr_tc_level_partition[i];
        const S32 next = adr_tc_level_partition[i + 1];
        const F64 len = morton_key.getCorrespondingFullLength(i);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for (S32 j = head; j < next; j++) {
            Ttc* tc_tmp = tc + j;
            const int n_tmp = tc_tmp->n_ptcl_;
            // tc_tmp->mom_.init();
            tc_tmp->clearMoment();
            tc_tmp->geo_.setSize(len);
            if (n_tmp == 0)
                continue;
            else if (tc_tmp->isLeaf(n_leaf_limit)) {
                const S32 adr = tc_tmp->adr_ptcl_;
                if (GetMSB(tp[adr].adr_ptcl_) == 0) {
                    const U32 adr_ep = tp[adr].adr_ptcl_;
                    tc_tmp->geo_.addFirstEpj(epj[adr_ep]);
                    // std::cout<<"adr_ep= "<<adr_ep
                    //          <<" epj[adr_ep].getPos()= "
                    //          <<epj[adr_ep].getPos()
                    //          <<std::endl;
                } else {
                    const U32 adr_sp = ClearMSB(tp[adr].adr_ptcl_);
                    tc_tmp->geo_.addFirstSpj(spj[adr_sp]);
                    // std::cout<<"adr_sp= "<<adr_sp
                    //          <<" spj[adr_sp].getPos()= "
                    //          <<spj[adr_sp].getPos()
                    //          <<std::endl;
                }
                // if(Comm::getRank()==0){
                //     std::cout<<"A) j= "<<j<<std::endl;
                //     tc_tmp->geo_.dump(std::cout);
                // }
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    if (GetMSB(tp[k].adr_ptcl_) == 0) {
                        const U32 adr_ep = tp[k].adr_ptcl_;
                        tc_tmp->mom_.accumulateAtLeaf(epj[adr_ep]);
                        tc_tmp->geo_.accumulateAtLeaf(epj[adr_ep]);
                    } else {
                        const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                        tc_tmp->mom_.accumulate(spj[adr_sp].convertToMoment());
                    }
                }
                // if(Comm::getRank()==0){
                //     std::cout<<"B) j= "<<j<<std::endl;
                //     tc_tmp->geo_.dump(std::cout);
                // }
                tc_tmp->mom_.set();
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    if (GetMSB(tp[k].adr_ptcl_) == 0) {
                        const U32 adr_ep = tp[k].adr_ptcl_;
                        tc_tmp->mom_.accumulateAtLeaf2(epj[adr_ep]);
                    } else {
                        const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                        tc_tmp->mom_.accumulate2(spj[adr_sp].convertToMoment());
                    }
                }
            } else {
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ != 0) {
                        tc_tmp->geo_.copyBox(tc_tmp_tmp->geo_);
                        break;
                    }
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate(tc_tmp_tmp->mom_);
                    tc_tmp->geo_.accumulate(tc_tmp_tmp->geo_);
                }
                tc_tmp->mom_.set();
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate2(tc_tmp_tmp->mom_);
                }
            }
        }
    }
}
#endif

#if 0
    // for PMMM
    template<class Ttc, class Tepj, class Tspj, typename Tmorton_key>
    inline void CalcMomentLongGlobalTreeP2P(S32 adr_tc_level_partition[],
                                            Ttc tc[],
                                            TreeParticle tp[],
                                            Tepj epj[],
                                            Tspj spj[],
                                            const S32 lev_max,
                                            const S32 n_leaf_limit,
                                            const CommTable & comm_table,
                                            const F64ort & pos_root_cell,
                                            const U32 adr_org_from_adr_sorted[],
                                            const Tmorton_key & morton_key){
        assert(0);
    }
#elif (1)
// new version befor PMMM
template <class Ttc, class Tepj, class Tspj, typename Tmorton_key>
inline void CalcMomentLongGlobalTreeP2P(S32 adr_tc_level_partition[], Ttc tc[], TreeParticle tp[], Tepj epj[], Tspj spj[], const S32 lev_max,
                                        const S32 n_leaf_limit, const CommTable& comm_table, const F64ort& pos_root_cell,
                                        const U32 adr_org_from_adr_sorted[], const Tmorton_key& morton_key) {
    for (S32 i = lev_max; i >= 0; --i) {
        const S32 head = adr_tc_level_partition[i];
        const S32 next = adr_tc_level_partition[i + 1];
        const F64 len = morton_key.getCorrespondingFullLength(i);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for (S32 j = head; j < next; j++) {
            Ttc* tc_tmp = tc + j;
            const int n_tmp = tc_tmp->n_ptcl_;
            tc_tmp->clearMoment();
            tc_tmp->geo_.setSize(len);
            if (n_tmp == 0)
                continue;
            else if (tc_tmp->isLeaf(n_leaf_limit)) {
                const S32 adr = tc_tmp->adr_ptcl_;
                if (GetMSB(tp[adr].adr_ptcl_) == 0) {
                    const U32 adr_ep = tp[adr].adr_ptcl_;
                    tc_tmp->geo_.addFirstEpj(epj[adr_ep]);
                } else {
                    const U32 adr_sp = ClearMSB(tp[adr].adr_ptcl_);
                    tc_tmp->geo_.addFirstSpj(spj[adr_sp]);
                }
                F64ort cell = morton_key.getPosTreeCell(i, tp[adr].key_);
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    if (GetMSB(tp[k].adr_ptcl_) == 0) {
                        const U32 adr_ep = tp[k].adr_ptcl_;
                        tc_tmp->mom_.accumulateAtLeaf(epj[adr_ep]);
                        tc_tmp->geo_.accumulateAtLeaf(epj[adr_ep]);
                    } else {
                        const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                        tc_tmp->mom_.accumulate(spj[adr_sp].convertToMoment());
                        const S32 adr_sp_org = ClearMSB(adr_org_from_adr_sorted[k]);
                        const S32 adr_sp_tmp = adr_sp_org - comm_table.n_sp_recv_tot_;
                        if (adr_sp_tmp >= 0) {
                            // assert(comm_table.pos_domain_allgather_.capacity() >= adr_sp_tmp);
                            const F64ort pos_domain = comm_table.pos_domain_allgather_[adr_sp_tmp];
                            const F64vec center = cell.getCenter();
                            const F64vec dxh = pos_domain.high_ - center;
                            const F64vec dxl = center - pos_domain.low_;
                            F64 half_len_cell = len * 0.5;
                            const F64 half_len = std::max(std::max(dxh.getMax(), dxl.getMax()), half_len_cell);
                            tc_tmp->geo_.setSize(std::max(tc_tmp->geo_.getSize(), half_len * 2.0));
                        }
                    }
                }
                tc_tmp->mom_.set();
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    if (GetMSB(tp[k].adr_ptcl_) == 0) {
                        const U32 adr_ep = tp[k].adr_ptcl_;
                        tc_tmp->mom_.accumulateAtLeaf2(epj[adr_ep]);
                    } else {
                        const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                        tc_tmp->mom_.accumulate2(spj[adr_sp].convertToMoment());
                    }
                }
            } else {
                const F64 child_full_len = len * 0.5;
                F64 dx_max = 0.0;
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    const Ttc& tc_tmp_tmp = tc[((tc_tmp->adr_tc_) + k)];
                    if (tc_tmp_tmp.n_ptcl_ != 0) {
                        tc_tmp->geo_.copyBox(tc_tmp_tmp.geo_);
                        break;
                    }
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    const Ttc& tc_tmp_tmp = tc[((tc_tmp->adr_tc_) + k)];
                    if (tc_tmp_tmp.n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate(tc_tmp_tmp.mom_);
                    tc_tmp->geo_.accumulate(tc_tmp_tmp.geo_);
                }
                tc_tmp->mom_.set();
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    const Ttc& tc_tmp_tmp = tc[((tc_tmp->adr_tc_) + k)];
                    dx_max = std::max(dx_max, tc_tmp_tmp.geo_.getSize() - child_full_len);
                    if (tc_tmp_tmp.n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate2(tc_tmp_tmp.mom_);
                }
                tc_tmp->geo_.setSize(len + dx_max);
            }
        }
    }
}
#else
// original
// #define MK_METHOD
template <class Ttc, class Tepj, class Tspj, typename Tmorton_key>
inline void CalcMomentLongGlobalTreeP2P(S32 adr_tc_level_partition[], Ttc tc[], TreeParticle tp[], Tepj epj[], Tspj spj[], const S32 lev_max,
                                        const S32 n_leaf_limit, const CommTable& comm_table, const F64ort& pos_root_cell,
                                        const U32 adr_org_from_adr_sorted[], const Tmorton_key& morton_key) {
    const S32 n_tc = adr_tc_level_partition[lev_max + 1];
    ReallocatableArray<F64ort> box(n_tc, n_tc, MemoryAllocMode::Pool);
    for (S32 i = lev_max; i >= 0; --i) {
        const S32 head = adr_tc_level_partition[i];
        const S32 next = adr_tc_level_partition[i + 1];
        const F64 len = morton_key.getCorrespondingFullLength(i);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for (S32 j = head; j < next; j++) {
            Ttc* tc_tmp = tc + j;
            F64ort& box_tmp = box[j];
            box_tmp.init();
            const int n_tmp = tc_tmp->n_ptcl_;
            // tc_tmp->mom_.init();
            tc_tmp->clearMoment();
            tc_tmp->geo_.setSize(len);
            if (n_tmp == 0)
                continue;
            else if (tc_tmp->isLeaf(n_leaf_limit)) {
                const S32 adr = tc_tmp->adr_ptcl_;
                if (GetMSB(tp[adr].adr_ptcl_) == 0) {
                    const U32 adr_ep = tp[adr].adr_ptcl_;
                    tc_tmp->geo_.addFirstEpj(epj[adr_ep]);
                } else {
                    const U32 adr_sp = ClearMSB(tp[adr].adr_ptcl_);
                    tc_tmp->geo_.addFirstSpj(spj[adr_sp]);
                }

                F64ort cell = morton_key.getPosTreeCell(i, tp[adr].key_);
                box_tmp = morton_key.getPosTreeCell(i, tp[adr].key_);
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    if (GetMSB(tp[k].adr_ptcl_) == 0) {
                        const U32 adr_ep = tp[k].adr_ptcl_;
                        tc_tmp->mom_.accumulateAtLeaf(epj[adr_ep]);
                        tc_tmp->geo_.accumulateAtLeaf(epj[adr_ep]);
                    } else {
                        const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                        tc_tmp->mom_.accumulate(spj[adr_sp].convertToMoment());
                        const S32 adr_sp_org = ClearMSB(adr_org_from_adr_sorted[k]);
                        const S32 adr_sp_tmp = adr_sp_org - comm_table.n_sp_recv_tot_;
                        if (adr_sp_tmp >= 0) {
#ifdef MK_METHOD
                            const F64ort pos_domain = comm_table.pos_domain_allgather_[adr_sp_tmp];
                            const F64vec center = cell.getCenter();
                            const F64vec dxh = pos_domain.high_ - center;
                            const F64vec dxl = center - pos_domain.low_;
                            F64 half_len_cell = len * 0.5;
                            const F64 half_len = std::max(std::max(dxh.getMax(), dxl.getMax()), half_len_cell);
                            tc_tmp->geo_.setSize(std::max(tc_tmp->geo_.getSize(), half_len * 2.0));
#elif 0
                            const F64vec center = cell.getCenter();
                            const F64ort pos_domain = comm_table.pos_domain_allgather_[adr_sp_tmp];
                            const F64vec dxh(fabs(box.high_.x - center.x), fabs(box.high_.y - center.y), fabs(box.high_.z - center.z));
                            const F64vec dxl(fabs(box.low_.x - center.x), fabs(box.low_.y - center.y), fabs(box.low_.z - center.z));
                            F64 half_len_cell = cell_tmp.getHalfLength()[0];
                            const F64 half_len = std::max(dxh.getMax(), dxl.getMax());
                            while (half_len_cell < half_len) {
                                half_len_cell *= 2.0;
                            }
                            const F64ort pos_domain(F64vec(center.x - half_len_cell, center.y - half_len_cell, center.z - half_len_cell),
                                                    F64vec(center.x + half_len_cell, center.y + half_len_cell, center.z + half_len_cell));
                            box_tmp.merge(pos_domain);
                            assert(pos_domain.contained(spj[adr_sp].getPos()));
#else
                            const F64ort pos_domain = comm_table.pos_domain_allgather_[adr_sp_tmp];
                            box_tmp.merge(pos_domain);
                            assert(pos_domain.contained(spj[adr_sp].getPos()));
#endif
                        }
                    }
                }
#ifndef MK_METHOD
                tc_tmp->geo_.setSize(box_tmp.getFullLength().getMax());
#endif
                tc_tmp->mom_.set();
                for (S32 k = adr; k < adr + n_tmp; k++) {
                    if (GetMSB(tp[k].adr_ptcl_) == 0) {
                        const U32 adr_ep = tp[k].adr_ptcl_;
                        tc_tmp->mom_.accumulateAtLeaf2(epj[adr_ep]);
                    } else {
                        const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                        tc_tmp->mom_.accumulate2(spj[adr_sp].convertToMoment());
                    }
                }
            } else {
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    F64ort& box_tmp_tmp = box[(tc_tmp->adr_tc_) + k];
                    box_tmp.merge(box_tmp_tmp);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate(tc_tmp_tmp->mom_);
                    // tc_tmp->geo_.accumulate( tc_tmp_tmp->geo_ );
                }
                tc_tmp->mom_.set();
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->mom_.accumulate2(tc_tmp_tmp->mom_);
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ != 0) {
                        tc_tmp->geo_.copyBox(tc_tmp_tmp->geo_);
                        break;
                    }
                }
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc* tc_tmp_tmp = tc + ((tc_tmp->adr_tc_) + k);
                    if (tc_tmp_tmp->n_ptcl_ == 0) continue;
                    tc_tmp->geo_.accumulate(tc_tmp_tmp->geo_);
                }
#ifdef MK_METHOD
                const F64 child_full_len = len * 0.5;
                F64 dx_max = 0.0;
                for (S32 k = 0; k < N_CHILDREN; k++) {
                    Ttc& tc_tmp_tmp = tc[(tc_tmp->adr_tc_) + k];
                    dx_max = std::max(dx_max, tc_tmp_tmp.geo_.getSize() - child_full_len);
                }
                tc_tmp->geo_.setSize(len + dx_max);
#else
                tc_tmp->geo_.setSize(box_tmp.getFullLength().getMax());
#endif
            }
        }
    }
}
#endif

//////////////////////////
//// CHECK CALC MOMENT ///
template <class Tep, class Ttc>
void CheckCalcMomentLongLocalTree(Tep ep[], const S32 n_ptcl, Ttc tc[], const F64 tolerance = 1e-5, std::ostream& fout = std::cout) {
    F64 mass_cm = 0.0;
    F64vec pos_cm = 0.0;
    for (S32 i = 0; i < n_ptcl; i++) {
        mass_cm += ep[i].getCharge();
        pos_cm += ep[i].getCharge() * ep[i].getPos();
    }
    pos_cm /= mass_cm;
    F64 dm = std::abs(tc[0].mom_.mass - mass_cm);
    F64vec dx = tc[0].mom_.pos - pos_cm;
    F64 dx_max = dx.applyEach(Abs<F64>()).getMax();
    if (dm >= tolerance || dx_max >= tolerance) {
        fout << "CalcMomentLong test: FAIL" << std::endl;
    } else {
        fout << "CalcMomentLong test: PASS" << std::endl;
    }
    if (Comm::getRank() == 0) {
        fout << "mass_cm=" << mass_cm << ", tc[0].mom_.mass=" << tc[0].mom_.mass << std::endl;
        fout << "pos_cm=" << pos_cm << ", tc[0].mom_.pos=" << tc[0].mom_.pos << std::endl;
    }
}

template <class Tep, class Tsp, class Ttc>
void CheckCalcMomentLongGlobalTree(const Tep ep[], const Tsp sp[], const S32 n_ptcl, const Ttc tc[], const TreeParticle tp[],
                                   const F64 tolerance = 1e-5, std::ostream& fout = std::cout) {
    F64 mass_cm_glb = 0.0;
    F64vec pos_cm_glb = 0.0;
    for (S32 i = 0; i < n_ptcl; i++) {
        if (GetMSB(tp[i].adr_ptcl_) == 0) {
            mass_cm_glb += ep[i].getCharge();
            pos_cm_glb += ep[i].getCharge() * ep[i].getPos();
        } else {
            mass_cm_glb += sp[i].getCharge();
            pos_cm_glb += sp[i].getCharge() * sp[i].getPos();
        }
    }
    pos_cm_glb /= mass_cm_glb;
    F64 dm = std::abs(tc[0].mom_.mass - mass_cm_glb);
    F64vec dx = tc[0].mom_.pos - pos_cm_glb;
    F64 dx_max = dx.applyEach(Abs<F64>()).getMax();
    if (dm >= tolerance || dx_max >= tolerance) {
        fout << "CheckCalcMomentLongGlobalTree: FAIL" << std::endl;
    } else {
        fout << "CheckCalcMomentLongGlobalTree: PASS" << std::endl;
    }
    if (Comm::getRank() == 0) {
        fout << "mass_cm_glb=" << mass_cm_glb << ", tc[0].mom_.mass=" << tc[0].mom_.mass << std::endl;
        fout << "pos_cm_glb=" << pos_cm_glb << ", tc[0].mom_.pos=" << tc[0].mom_.pos << std::endl;
    }
}

template <class Tep, class Ttc>
void CheckCalcMomentShortInAndOut(Tep ep[], const S32 n_ptcl, Ttc tc[], const F64 tolerance = 1e-5, std::ostream& fout = std::cout) {
    F64ort vertex_out_tmp;
    F64ort vertex_in_tmp;
    vertex_out_tmp.init();
    vertex_in_tmp.init();
    for (S32 i = 0; i < n_ptcl; i++) {
        vertex_out_tmp.merge(ep[i].getPos(), ep[i].getRSearch());
        vertex_in_tmp.merge(ep[i].getPos());
    }
    /*
    F64 d_out_high = (tc[0].mom_.vertex_out_.high_ - vertex_out_tmp.high_).applyEach( Abs<F64>() ).getMax();
    F64 d_out_low = (tc[0].mom_.vertex_out_.low_ - vertex_out_tmp.low_).applyEach( Abs<F64>() ).getMax();
    F64 d_in_high = (tc[0].mom_.vertex_in_.high_ - vertex_in_tmp.high_).applyEach( Abs<F64>() ).getMax();
    F64 d_in_low = (tc[0].mom_.vertex_in_.low_ - vertex_in_tmp.low_).applyEach( Abs<F64>() ).getMax();
    if( d_out_high >= tolerance || d_out_low >= tolerance || d_in_high >= tolerance || d_in_low >= tolerance){
        fout<<"CalcMomentShort test: FAIL"<<std::endl;
        fout<<"vertex_out_tmp.high_="<<vertex_out_tmp.high_<<" tc[0].mom_.vertex_out_.high_="<<tc[0].mom_.vertex_out_.high_<<std::endl;
        fout<<"vertex_out_tmp.low_="<<vertex_out_tmp.low_<<" tc[0].mom_.vertex_out_.low_="<<tc[0].mom_.vertex_out_.low_<<std::endl;
        fout<<"vertex_in_tmp.high_="<<vertex_in_tmp.high_<<" tc[0].mom_.vertex_in_.high_="<<tc[0].mom_.vertex_in_.high_<<std::endl;
        fout<<"vertex_in_tmp.low_="<<vertex_in_tmp.low_<<" tc[0].mom_.vertex_in_.low_="<<tc[0].mom_.vertex_in_.low_<<std::endl;
    }
    */
    F64 d_out_high = (tc[0].geo_.vertex_out_.high_ - vertex_out_tmp.high_).applyEach(Abs<F64>()).getMax();
    F64 d_out_low = (tc[0].geo_.vertex_out_.low_ - vertex_out_tmp.low_).applyEach(Abs<F64>()).getMax();
    F64 d_in_high = (tc[0].geo_.vertex_in_.high_ - vertex_in_tmp.high_).applyEach(Abs<F64>()).getMax();
    F64 d_in_low = (tc[0].geo_.vertex_in_.low_ - vertex_in_tmp.low_).applyEach(Abs<F64>()).getMax();
    if (d_out_high >= tolerance || d_out_low >= tolerance || d_in_high >= tolerance || d_in_low >= tolerance) {
        fout << "CalcMomentShort test: FAIL" << std::endl;
        fout << "vertex_out_tmp.high_=" << vertex_out_tmp.high_ << " tc[0].geo_.vertex_out_.high_=" << tc[0].geo_.vertex_out_.high_ << std::endl;
        fout << "vertex_out_tmp.low_=" << vertex_out_tmp.low_ << " tc[0].geo_.vertex_out_.low_=" << tc[0].geo_.vertex_out_.low_ << std::endl;
        fout << "vertex_in_tmp.high_=" << vertex_in_tmp.high_ << " tc[0].geo_.vertex_in_.high_=" << tc[0].geo_.vertex_in_.high_ << std::endl;
        fout << "vertex_in_tmp.low_=" << vertex_in_tmp.low_ << " tc[0].geo_.vertex_in_.low_=" << tc[0].geo_.vertex_in_.low_ << std::endl;
    } else {
        fout << "CalcMomentShort test: PASS" << std::endl;
    }
}

template <class Tep, class Ttc>
void CheckCalcMomentShortInOnly(Tep ep[], const S32 n_ptcl, Ttc tc[], const F64 tolerance = 1e-5, std::ostream& fout = std::cout) {
    F64ort vertex_in_tmp;
    vertex_in_tmp.init();
    for (S32 i = 0; i < n_ptcl; i++) {
        vertex_in_tmp.merge(ep[i].getPos());
    }
    // F64 d_in_high = (tc[0].mom_.vertex_in_.high_ - vertex_in_tmp.high_).applyEach( Abs<F64>() ).getMax();
    // F64 d_in_low = (tc[0].mom_.vertex_in_.low_ - vertex_in_tmp.low_).applyEach( Abs<F64>() ).getMax();
    F64 d_in_high = (tc[0].geo_.vertex_in_.high_ - vertex_in_tmp.high_).applyEach(Abs<F64>()).getMax();
    F64 d_in_low = (tc[0].geo_.vertex_in_.low_ - vertex_in_tmp.low_).applyEach(Abs<F64>()).getMax();
    if (d_in_high >= tolerance || d_in_low >= tolerance) {
        fout << "CalcMomentShort test: FAIL" << std::endl;
        // fout<<"vertex_in_tmp.high_="<<vertex_in_tmp.high_<<" tc[0].mom_.vertex_in_.high_="<<tc[0].mom_.vertex_in_.high_<<std::endl;
        // fout<<"vertex_in_tmp.low_="<<vertex_in_tmp.low_<<" tc[0].mom_.vertex_in_.low_="<<tc[0].mom_.vertex_in_.low_<<std::endl;
        fout << "vertex_in_tmp.high_=" << vertex_in_tmp.high_ << " tc[0].geo_.vertex_in_.high_=" << tc[0].geo_.vertex_in_.high_ << std::endl;
        fout << "vertex_in_tmp.low_=" << vertex_in_tmp.low_ << " tc[0].geo_.vertex_in_.low_=" << tc[0].geo_.vertex_in_.low_ << std::endl;
    } else {
        fout << "CalcMomentShort test: PASS" << std::endl;
    }
}

template <class Tep, class Ttc>
void CheckCalcMomentOutOnly(Tep ep[], const S32 n_ptcl, Ttc tc[], const F64 tolerance = 1e-5, std::ostream& fout = std::cout) {
    F64ort vertex_out_tmp;
    vertex_out_tmp.init();
    for (S32 i = 0; i < n_ptcl; i++) {
        vertex_out_tmp.merge(ep[i].getPos(), ep[i].getRSearch());
    }
    // F64 d_out_high = (tc[0].mom_.vertex_out_.high_ - vertex_out_tmp.high_).applyEach( Abs<F64>() ).getMax();
    // F64 d_out_low = (tc[0].mom_.vertex_out_.low_ - vertex_out_tmp.low_).applyEach( Abs<F64>() ).getMax();
    F64 d_out_high = (tc[0].geo_.vertex_out_.high_ - vertex_out_tmp.high_).applyEach(Abs<F64>()).getMax();
    F64 d_out_low = (tc[0].geo_.vertex_out_.low_ - vertex_out_tmp.low_).applyEach(Abs<F64>()).getMax();
    if (d_out_high >= tolerance || d_out_low >= tolerance) {
        fout << "CalcMomentShort test: FAIL" << std::endl;
        // fout<<"vertex_out_tmp.high_="<<vertex_out_tmp.high_<<" tc[0].mom_.vertex_out_.high_="<<tc[0].mom_.vertex_out_.high_<<std::endl;
        // fout<<"vertex_out_tmp.low_="<<vertex_out_tmp.low_<<" tc[0].mom_.vertex_out_.low_="<<tc[0].mom_.vertex_out_.low_<<std::endl;
        fout << "vertex_out_tmp.high_=" << vertex_out_tmp.high_ << " tc[0].geo_.vertex_out_.high_=" << tc[0].geo_.vertex_out_.high_ << std::endl;
        fout << "vertex_out_tmp.low_=" << vertex_out_tmp.low_ << " tc[0].geo_.vertex_out_.low_=" << tc[0].geo_.vertex_out_.low_ << std::endl;
    } else {
        fout << "CalcMomentOut test: PASS" << std::endl;
    }
}
}  // namespace ParticleSimulator
