#define ORIGINAL_SCATTER_MODE

#ifdef __HPC_ACE__
#include<emmintrin.h>
#endif

namespace ParticleSimulator{
    inline void CalcNumberAndShiftOfImageDomain
    (ReallocatableArray<F64vec> & shift_image_domain,
     const F64vec & size_root_domain,
     const F64ort & pos_my_domain,
     const F64ort & pos_target_domain,
     const bool periodic_axis[]){
        shift_image_domain.clearSize();
        if(pos_my_domain.overlapped(pos_target_domain)){
            shift_image_domain.push_back(F64vec(0.0)); // NOTE: sign is plus
        }
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        for(S32 lev=1;; lev++){
            S32 n_image_per_level = 0;
            S32 lev_x = 0;
            S32 lev_y = 0;
            if(periodic_axis[0] == true) lev_x = lev;
            if(periodic_axis[1] == true) lev_y = lev;
            for(S32 ix=-lev_x; ix<=lev_x; ix++){
                for(S32 iy=-lev_y; iy<=lev_y; iy++){
                    if( (std::abs(ix) !=lev_x || periodic_axis[0] == false) && (std::abs(iy) != lev_y  || periodic_axis[1] == false)) continue;
                    const F64vec shift_tmp(ix*size_root_domain.x, iy*size_root_domain.y);
                    const F64ort pos_image_tmp = pos_target_domain.shift(shift_tmp);
                    if(pos_my_domain.overlapped(pos_image_tmp)){
                        shift_image_domain.push_back(shift_tmp); // NOTE: sign is plus
                        n_image_per_level++;
                    }
                }
            }
            if(n_image_per_level == 0) break;
        }
#else
        for(S32 lev=1;; lev++){
            S32 n_image_per_level = 0;
            S32 lev_x = 0;
            S32 lev_y = 0;
            S32 lev_z = 0;
            if(periodic_axis[0] == true) lev_x = lev;
            if(periodic_axis[1] == true) lev_y = lev;
            if(periodic_axis[2] == true) lev_z = lev;
            for(S32 ix=-lev_x; ix<=lev_x; ix++){
                for(S32 iy=-lev_y; iy<=lev_y; iy++){
                    for(S32 iz=-lev_z; iz<=lev_z; iz++){
                        if( (std::abs(ix) !=lev_x || periodic_axis[0] == false) && (std::abs(iy) != lev_y  || periodic_axis[1] == false) && (std::abs(iz) != lev_z || periodic_axis[2] == false) ) continue;
                        const F64vec shift_tmp(ix*size_root_domain.x, iy*size_root_domain.y, iz*size_root_domain.z);
                        const F64ort pos_image_tmp = pos_target_domain.shift(shift_tmp);
                        if(pos_my_domain.overlapped(pos_image_tmp)){
                            shift_image_domain.push_back(shift_tmp); // NOTE: sign is plus
                            n_image_per_level++;
                        }
                    }
                }
            }
            if(n_image_per_level == 0) break;
        }
#endif
    }

    
    inline S32vec CalcIDOfImageDomain(const F64ort & pos_root_domain,
                                      const F64vec & pos_target,
                                      const bool periodic_axis[]){
        if( pos_root_domain.overlapped(pos_target) ){
            return S32vec(0.0);
        }
        const F64vec size_root_domain = pos_root_domain.getFullLength();
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        for(S32 lev=1;; lev++){
            S32 lev_x = 0;
            S32 lev_y = 0;
            if(periodic_axis[0] == true) lev_x = lev;
            if(periodic_axis[1] == true) lev_y = lev;
            for(S32 ix=-lev_x; ix<=lev_x; ix++){
                for(S32 iy=-lev_y; iy<=lev_y; iy++){
                    if( (std::abs(ix) !=lev_x || periodic_axis[0] == false) 
                        && (std::abs(iy) != lev_y  || periodic_axis[1] == false) ) continue;
                    const F64vec shift_tmp(ix*size_root_domain.x, iy*size_root_domain.y);
                    if( pos_root_domain.overlapped(pos_target + shift_tmp) ){
                        return S32vec(-ix, -iy);
                    }
                }
            }
        }
#else
        for(S32 lev=1;; lev++){
            S32 lev_x = 0;
            S32 lev_y = 0;
            S32 lev_z = 0;
            if(periodic_axis[0] == true) lev_x = lev;
            if(periodic_axis[1] == true) lev_y = lev;
            if(periodic_axis[2] == true) lev_z = lev;
            for(S32 ix=-lev_x; ix<=lev_x; ix++){
                for(S32 iy=-lev_y; iy<=lev_y; iy++){
                    for(S32 iz=-lev_z; iz<=lev_z; iz++){
                        if( (std::abs(ix) !=lev_x || periodic_axis[0] == false) 
                            && (std::abs(iy) != lev_y  || periodic_axis[1] == false) 
                            && (std::abs(iz) != lev_z || periodic_axis[2] == false) ) continue;
                        const F64vec shift_tmp(ix*size_root_domain.x, iy*size_root_domain.y, iz*size_root_domain.z);
                        if( pos_root_domain.overlapped(pos_target + shift_tmp) ){
                            return S32vec(-ix, -iy, -iz);
                        }
                    }
                }
            }
        }
#endif
    }

    template<class Tptcl>
    inline void AllGatherParticle(Tptcl *& ptcl,
                                  S32 *&n_ptcl,
                                  S32 *&n_ptcl_disp,
                                  const Tptcl ptcl_org[], 
                                  const S32 n_loc){
        n_ptcl = new S32[Comm::getNumberOfProc()];
        n_ptcl_disp = new S32[Comm::getNumberOfProc()+1];
        Comm::allGather(&n_loc, 1, n_ptcl); // TEST
        n_ptcl_disp[0] = 0;
        for(S32 i=0; i<Comm::getNumberOfProc(); i++){
            n_ptcl_disp[i+1] = n_ptcl_disp[i] + n_ptcl[i];
        }
        const S32 n_tot = n_ptcl_disp[Comm::getNumberOfProc()];
        ptcl = new Tptcl[ n_tot ];
        Comm::allGatherV(ptcl_org, n_loc, ptcl, n_ptcl, n_ptcl_disp); // TEST
    }

    template<class Tptcl>
    inline void AllGatherParticle(Tptcl *& ptcl, 
                                  S32 & n_tot, 
                                  const Tptcl ptcl_org[], 
                                  const S32 n_loc){
        //S32 * n_tmp = new S32[Comm::getNumberOfProc()];
        //S32 * n_disp_tmp = new S32[Comm::getNumberOfProc()+1];
        S32 * n_tmp;
        S32 * n_disp_tmp;
        AllGatherParticle(ptcl, n_tmp, n_disp_tmp, ptcl_org, n_loc);
        n_tot = n_disp_tmp[Comm::getNumberOfProc()];
        delete [] n_tmp;
        delete [] n_disp_tmp;
    }


    template<class Tptcl>
    inline void AllGatherParticle(Tptcl *& ptcl,
                                  S32 *&n_ptcl,
                                  S32 *&n_ptcl_disp,
                                  const Tptcl ptcl_org[],
                                  const S32 n_loc,
                                  const F64vec & domain_size,
                                  const F64ort & boundary,
                                  const bool periodic_axis[]){
        const S32 n_proc = Comm::getNumberOfProc();
        //std::cerr<<"n_proc="<<n_proc<<std::endl;
        n_ptcl = new S32[n_proc];
        n_ptcl_disp = new S32[n_proc+1];
        //std::cerr<<"n_loc="<<n_loc<<std::endl;
        Comm::allGather(&n_loc, 1, n_ptcl); // TEST
        //std::cerr<<"n_ptcl[0]="<<n_ptcl[0]<<std::endl;
        n_ptcl_disp[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ptcl_disp[i+1] = n_ptcl_disp[i] + n_ptcl[i];
        }
        const S32 n_tot = n_ptcl_disp[n_proc];
        //std::cerr<<"n_tot="<<n_tot<<std::endl;
        Tptcl * ptcl_tmp = new Tptcl[ n_tot ];
        //std::cerr<<"n_tot="<<n_tot<<std::endl;
        Comm::allGatherV(ptcl_org, n_loc, ptcl_tmp, n_ptcl, n_ptcl_disp); // TEST
        ReallocatableArray<F64vec> shift;
        CalcNumberAndShiftOfImageDomain(shift, domain_size, 
                                        boundary, boundary, periodic_axis);
        const S32 n_image = shift.size();
        //std::cerr<<"n_image="<<n_image<<std::endl;
        ptcl = new Tptcl[ n_tot*n_image ];
        //std::cerr<<"n_tot*n_image="<<n_tot*n_image<<std::endl;
        for(S32 i=0; i<n_tot; i++){
            for(S32 ii=0; ii<n_image; ii++){
                ptcl[i*n_image+ii] = ptcl_tmp[i];
                const F64vec pos_new = ptcl_tmp[i].getPos() + shift[ii];
                ptcl[i*n_image+ii].setPos(pos_new);
            }
        }
        for(S32 i=0; i<n_proc; i++){
            n_ptcl[i] *= n_image;
        }
        n_ptcl_disp[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ptcl_disp[i+1] = n_ptcl_disp[i] + n_ptcl[i];
        }

        delete [] ptcl_tmp;
    }

    template<class Tptcl>
    inline void AllGatherParticle(Tptcl *& ptcl, 
                                  S32 & n_tot, 
                                  const Tptcl ptcl_org[], 
                                  const S32 n_loc,
                                  const F64vec & domain_size,
                                  const F64ort & boundary,
                                  const bool periodic_axis[]){
        S32 * n_tmp;
        S32 * n_disp_tmp;
        AllGatherParticle(ptcl, n_tmp, n_disp_tmp, ptcl_org, n_loc, domain_size, boundary, periodic_axis);
        n_tot = n_disp_tmp[Comm::getNumberOfProc()];
        delete [] n_tmp;
        delete [] n_disp_tmp;
    }


#if 0
    // point to the first address
    // liner search
    template<class Ttp>
    inline S32 GetPartitionID(const Ttp tp[],  
                              const S32 left,
                              const S32 right, 
                              const S32 ref,
                              const S32 lev){
        for(S32 i=left; i<right+1; i++){
            //S32 tmp = MortonKey<DIMENSION>::getCellID(lev, tp[i].getKey());
	    S32 tmp = MortonKey::getCellID(lev, tp[i].getKey());
            if(ref == tmp) return i;
        }
        return -1;
    }
#else
    // binary search
    template<class Ttp>
    inline S32 GetPartitionID(const Ttp tp[],
                              const S32 left,
                              const S32 right,
                              const S32 ref,
                              const S32 lev){
        //if( MortonKey<DIMENSION>::getCellID(lev, tp[right].getKey()) < ref) return -1;
        //if( MortonKey<DIMENSION>::getCellID(lev, tp[left].getKey()) > ref) return -1;
	if( MortonKey::getCellID(lev, tp[right].getKey()) < ref) return -1;
        if( MortonKey::getCellID(lev, tp[left].getKey()) > ref) return -1;
        S32 l = left;
        S32 r = right;
        while(l < r){
            S32 cen = (l+r) / 2;
            //S32 val_cen = MortonKey<DIMENSION>::getCellID(lev, tp[cen].getKey());
	    S32 val_cen = MortonKey::getCellID(lev, tp[cen].getKey());
            if(ref > val_cen){
                l = cen+1;
            }
            else{
                r = cen;
            }
            //if( MortonKey<DIMENSION>::getCellID(lev, tp[l].getKey()) == ref ){ return l;
	    if( MortonKey::getCellID(lev, tp[l].getKey()) == ref ) return l;
	}
	return -1;
    }
#endif
  
    template<class Ttc>
    inline void LinkCell(ReallocatableArray<Ttc> & tc_array,
                         S32 adr_tc_level_partition[],
                         const TreeParticle tp[],
                         S32 & lev_max,
                         const S32 n_tot,
                         const S32 n_leaf_limit){
        tc_array.resizeNoInitialize(N_CHILDREN*2);
        tc_array[0].n_ptcl_ = n_tot;
        tc_array[0].adr_tc_ = N_CHILDREN;
        tc_array[0].adr_ptcl_ = 0;
        adr_tc_level_partition[0] = 0;
        adr_tc_level_partition[1] = N_CHILDREN;
        lev_max = 0;
        tc_array[0].level_ = lev_max;
        for(S32 k=1; k<N_CHILDREN; k++){
            tc_array[k].level_ = 0;
            tc_array[k].n_ptcl_ = 0;
            tc_array[k].adr_tc_ = SetMSB( (U32)(0) );
            tc_array[k].adr_ptcl_ = SetMSB( (U32)(0) );
        }
        if(tc_array[0].n_ptcl_ <= n_leaf_limit)return;
        S32 id_cell_left = 0;
        S32 id_cell_right = N_CHILDREN-1;
        while(1){
            S32 n_cell_new = 0;
            // assign particles to child cells and count # of particles in child cells
            // but loop over parent cells because they have indexes of particles
#pragma omp parallel for reduction(+: n_cell_new)
            for(S32 i=id_cell_left; i<id_cell_right+1; i++){
                const S32 n_ptcl_tmp = tc_array[i].n_ptcl_;
                if(n_ptcl_tmp <= n_leaf_limit) continue;
                const S32 adr_ptcl_tmp = tc_array[i].adr_ptcl_;
                S32 adr[N_CHILDREN];
                S32 n_cnt = 0;
                for(S32 j=0; j<N_CHILDREN; j++){
                    const S32 adr_tc_tmp = tc_array[i].adr_tc_ + j;
                    tc_array[adr_tc_tmp].adr_ptcl_ = GetPartitionID(tp, adr_ptcl_tmp, adr_ptcl_tmp+n_ptcl_tmp-1, j, lev_max+1);
                    tc_array[adr_tc_tmp].level_ = lev_max+1;
                    if(GetMSB(tc_array[adr_tc_tmp].adr_ptcl_)){
// GetPartitionID returnes -1, if there is no appropriate tree cell.
                        tc_array[adr_tc_tmp].n_ptcl_ = 0;
                        //tc[adr_tc_tmp].adr_tc_ = -1;
                        tc_array[adr_tc_tmp].adr_tc_ = SetMSB( (U32)(0) );
                    }
                    else{
                        adr[n_cnt] = adr_tc_tmp;
                        n_cnt++;
                    }
                }
                S32 n_ptcl_cum = 0;
                for(S32 j=0; j<n_cnt-1; j++){
                    tc_array[adr[j]].n_ptcl_ = tc_array[adr[j+1]].adr_ptcl_ - tc_array[adr[j]].adr_ptcl_;
                    n_ptcl_cum += tc_array[adr[j]].n_ptcl_;
                }
                tc_array[adr[n_cnt-1]].n_ptcl_ = n_ptcl_tmp - n_ptcl_cum;
                n_cell_new += N_CHILDREN;
            } // end of omp parallel scope
            if( n_cell_new == 0) break;
            // go deeper
            id_cell_left = id_cell_right + 1;
            id_cell_right += n_cell_new;
            lev_max++;
            adr_tc_level_partition[lev_max+1] = id_cell_right + 1;
            if(lev_max == TREE_LEVEL_LIMIT) break;
#if 0
            // probably need prefix sum
            // under construction
            //#pragma omp parallel for
            for(S32 i=id_cell_left; i<id_cell_right+1; i++){
                const S32 id_tmp = i;
            }
#else
            S32 offset = id_cell_right+1;
            for(S32 i=id_cell_left; i<id_cell_right+1; i++){
                const S32 id_tmp = i;
                if(!tc_array[id_tmp].isLeaf(n_leaf_limit)){
                    tc_array[id_tmp].adr_tc_ = offset;
                    offset += N_CHILDREN;
                }
            }
#endif
            tc_array.resizeNoInitialize(offset);
        }
    }



    template<class Ttc>
    inline void LinkCellST(ReallocatableArray<Ttc> & tc_array,
                           S32 adr_tc_level_partition[],
                           const TreeParticle tp[],
                           S32 & lev_max,
                           const S32 n_tot,
                           const S32 n_leaf_limit){
        tc_array.resizeNoInitialize(N_CHILDREN*2);
        tc_array[0].n_ptcl_ = n_tot;
        tc_array[0].adr_tc_ = N_CHILDREN;
        tc_array[0].adr_ptcl_ = 0;
        adr_tc_level_partition[0] = 0;
        adr_tc_level_partition[1] = N_CHILDREN;
        lev_max = 0;
        tc_array[0].level_ = lev_max;
        for(S32 k=1; k<N_CHILDREN; k++){
            tc_array[k].level_ = 0;
            tc_array[k].n_ptcl_ = 0;
            tc_array[k].adr_tc_ = SetMSB( (U32)(0) );
            tc_array[k].adr_ptcl_ = SetMSB( (U32)(0) );
        }
        if(tc_array[0].n_ptcl_ <= n_leaf_limit)return;
        S32 id_cell_left = 0;
        S32 id_cell_right = N_CHILDREN-1;
        while(1){
            S32 n_cell_new = 0;
            // assign particles to child cells and count # of particles in child cells
            // but loop over parent cells because they have indexes of particles
            for(S32 i=id_cell_left; i<id_cell_right+1; i++){
                const S32 n_ptcl_tmp = tc_array[i].n_ptcl_;
                if(n_ptcl_tmp <= n_leaf_limit) continue;
                const S32 adr_ptcl_tmp = tc_array[i].adr_ptcl_;
                S32 adr[N_CHILDREN];
                S32 n_cnt = 0;
                for(S32 j=0; j<N_CHILDREN; j++){
                    const S32 adr_tc_tmp = tc_array[i].adr_tc_ + j;
                    tc_array[adr_tc_tmp].adr_ptcl_ = GetPartitionID(tp, adr_ptcl_tmp, adr_ptcl_tmp+n_ptcl_tmp-1, j, lev_max+1);
                    tc_array[adr_tc_tmp].level_ = lev_max+1;
                    if(GetMSB(tc_array[adr_tc_tmp].adr_ptcl_)){
                        tc_array[adr_tc_tmp].n_ptcl_ = 0;
                        //tc[adr_tc_tmp].adr_tc_ = -1;
                        tc_array[adr_tc_tmp].adr_tc_ = SetMSB( (U32)(0) );
                    }
                    else{
                        adr[n_cnt] = adr_tc_tmp;
                        n_cnt++;
                    }
                }
                S32 n_ptcl_cum = 0;
                for(S32 j=0; j<n_cnt-1; j++){
                    tc_array[adr[j]].n_ptcl_ = tc_array[adr[j+1]].adr_ptcl_ - tc_array[adr[j]].adr_ptcl_;
                    n_ptcl_cum += tc_array[adr[j]].n_ptcl_;
                }
                tc_array[adr[n_cnt-1]].n_ptcl_ = n_ptcl_tmp - n_ptcl_cum;
                n_cell_new += N_CHILDREN;
            } // end of omp parallel scope
            if( n_cell_new == 0) break;
            // go deeper
            id_cell_left = id_cell_right + 1;
            id_cell_right += n_cell_new;
            lev_max++;
            adr_tc_level_partition[lev_max+1] = id_cell_right + 1;
            if(lev_max == TREE_LEVEL_LIMIT) break;
            S32 offset = id_cell_right+1;
            for(S32 i=id_cell_left; i<id_cell_right+1; i++){
                const S32 id_tmp = i;
                if(!tc_array[id_tmp].isLeaf(n_leaf_limit)){
                    tc_array[id_tmp].adr_tc_ = offset;
                    offset += N_CHILDREN;
                }
            }
            tc_array.resizeNoInitialize(offset);
        }
    }


    /////////////////////
    //// CALC MOMENT ////
    // from bottom to top
    template<class Ttc, class Tepj>
    inline void CalcMoment(S32 adr_tc_level_partition[], 
                           Ttc tc[],
                           Tepj epj[],
                           const S32 lev_max,
                           const S32 n_leaf_limit){
        for(S32 i=lev_max; i>=0; --i){
            const S32 head = adr_tc_level_partition[i];
            const S32 next = adr_tc_level_partition[i+1];
#pragma omp parallel for
            for(S32 j=head; j<next; j++){
                Ttc * tc_tmp = tc + j;
                const int n_tmp = tc_tmp->n_ptcl_;
                tc_tmp->mom_.init();
                if( n_tmp == 0) continue;
                else if(tc_tmp->isLeaf(n_leaf_limit)){
                    const S32 adr = tc_tmp->adr_ptcl_;
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                    }
                    tc_tmp->mom_.set();
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
                    }
                }
                else{
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
                }
            }
        }
    }

    template<class Ttc, class Tepj>
    inline void CalcMomentST(S32 adr_tc_level_partition[], 
			     Ttc tc[],
			     Tepj epj[],
			     const S32 lev_max,
			     const S32 n_leaf_limit){
        for(S32 i=lev_max; i>=0; --i){
            const S32 head = adr_tc_level_partition[i];
            const S32 next = adr_tc_level_partition[i+1];
            for(S32 j=head; j<next; j++){
                Ttc * tc_tmp = tc + j;
                const int n_tmp = tc_tmp->n_ptcl_;
                tc_tmp->mom_.init();
                if( n_tmp == 0) continue;
                else if(tc_tmp->isLeaf(n_leaf_limit)){
                    const S32 adr = tc_tmp->adr_ptcl_;
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                    }
                    tc_tmp->mom_.set();
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
                    }
                }
                else{
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
                }
            }
        }
    }

    template<class Ttc, class Tepj, class Tspj>
    inline void CalcMomentLongGlobalTree(S32 adr_tc_level_partition[], 
                                         Ttc tc[],
                                         TreeParticle tp[],
                                         Tepj epj[],
                                         Tspj spj[],
                                         const S32 lev_max,
                                         const S32 n_leaf_limit){
        for(S32 i=lev_max; i>=0; --i){
            const S32 head = adr_tc_level_partition[i];
            const S32 next = adr_tc_level_partition[i+1];
#pragma omp parallel for
            for(S32 j=head; j<next; j++){
                Ttc * tc_tmp = tc + j;
                const int n_tmp = tc_tmp->n_ptcl_;
                tc_tmp->mom_.init();
                if( n_tmp == 0) continue;
                else if(tc_tmp->isLeaf(n_leaf_limit)){
                    const S32 adr = tc_tmp->adr_ptcl_;
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        if( GetMSB(tp[k].adr_ptcl_) == 0 ){
                            tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                        }
                        else{
                            tc_tmp->mom_.accumulate(spj[k].convertToMoment());
                        }
                    }
                    tc_tmp->mom_.set();
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        if( GetMSB(tp[k].adr_ptcl_) == 0 ){
                            tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
                        }
                        else{
                            tc_tmp->mom_.accumulate2(spj[k].convertToMoment());
                        }
                    }
                }
                else{
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
                }
            }
        }
    }

/*
// the same as CalcMomentLongGlobalTree
    template<class Ttc, class Tepj, class Tspj>
    inline void CalcMomentLongCutoffGlobalTree(S32 adr_tc_level_partition[], 
                                               Ttc tc[],
                                               TreeParticle tp[],
                                               Tepj epj[],
                                               Tspj spj[],
                                               const S32 lev_max,
                                               const S32 n_leaf_limit){
        for(S32 i=lev_max; i>=0; --i){
            const S32 head = adr_tc_level_partition[i];
            const S32 next = adr_tc_level_partition[i+1];
#pragma omp parallel for
            for(S32 j=head; j<next; j++){
                Ttc * tc_tmp = tc + j;
                const int n_tmp = tc_tmp->n_ptcl_;
                tc_tmp->mom_.init();
                if( n_tmp == 0) continue;
                else if(tc_tmp->isLeaf(n_leaf_limit)){
                    const S32 adr = tc_tmp->adr_ptcl_;
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        if( GetMSB(tp[k].adr_ptcl_) == 0 ){
                            tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                        }
                        else{
                            tc_tmp->mom_.accumulate(spj[k].convertToMoment());
                        }
                    }
                    tc_tmp->mom_.set();
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        if( GetMSB(tp[k].adr_ptcl_) == 0 ){
                            tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
                        }
                        else{
                            tc_tmp->mom_.accumulate2(spj[k].convertToMoment());
                        }
                    }
                }
                else{
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
                }
            }
        }
    }
*/
    
    template<class Tipg, class Ttc, class Tepi>
    inline void MakeIPGroupLong(ReallocatableArray<Tipg> & ipg_first,
                                const ReallocatableArray<Ttc> & tc_first,
                                const ReallocatableArray<Tepi> & epi_first,
                                const S32 adr_tc,
                                const S32 n_grp_limit){
        const Ttc * tc_tmp = tc_first.getPointer(adr_tc);
        const S32 n_tmp = tc_tmp->n_ptcl_;
        if(n_tmp == 0) return;
        else if( tc_tmp->isLeaf(n_grp_limit) ){
            ipg_first.increaseSize();
            ipg_first.back().copyFromTC(*tc_tmp);
            ipg_first.back().vertex_ = GetMinBoxSingleThread(epi_first.data()+(tc_tmp->adr_ptcl_), n_tmp);
            return;
        }
        else{
            S32 adr_tc_tmp = tc_tmp->adr_tc_;
            for(S32 i=0; i<N_CHILDREN; i++){
                MakeIPGroupLong<Tipg, Ttc, Tepi>
                    (ipg_first, tc_first, epi_first, adr_tc_tmp+i, n_grp_limit);
            }
        }
    }

  // NOTE: This FUnction must be used after calc moment of global tree
    template<class Tipg, class Ttcloc, class Ttcglb, class Tepi>
    inline void MakeIPGroupUseGLBTreeLong(ReallocatableArray<Tipg> & ipg_first,
					  const ReallocatableArray<Ttcloc> & tc_loc_first,
					  const ReallocatableArray<Ttcglb> & tc_glb_first,
					  const ReallocatableArray<Tepi> & epi_first,
					  const S32 adr_tc_loc,
					  const S32 adr_tc_glb,
					  const S32 n_grp_limit,
					  const S32 n_leaf_limit){
        const Ttcloc * tc_loc_tmp = tc_loc_first.getPointer(adr_tc_loc);
        const Ttcglb * tc_glb_tmp = tc_glb_first.getPointer(adr_tc_glb);
        const S32 n_loc_tmp = tc_loc_tmp->n_ptcl_;
        const S32 n_glb_tmp = tc_glb_tmp->n_ptcl_;
	/*
	if(n_glb_tmp < n_loc_tmp){
	  std::cout<<"n_loc_tmp="<<n_loc_tmp<<std::endl;
	  std::cout<<"tc_loc_tmp.level_="<<tc_loc_tmp->level_<<std::endl;
	  std::cout<<"n_glb_tmp="<<n_glb_tmp<<std::endl;
	  std::cout<<"tc_glb_tmp.level_="<<tc_glb_tmp->level_<<std::endl;
	}
	assert(n_glb_tmp >= n_loc_tmp);
	*/
        if(n_loc_tmp == 0 || n_glb_tmp == 0) return;
        else if( tc_glb_tmp->isLeaf(n_grp_limit) || tc_loc_tmp->isLeaf(n_leaf_limit)){
            ipg_first.increaseSize();
            ipg_first.back().copyFromTC(*tc_loc_tmp);
	    const S32 adr_tmp = tc_loc_tmp->adr_ptcl_;
	    ipg_first.back().vertex_ = GetMinBoxSingleThread(epi_first.getPointer(adr_tmp), n_loc_tmp);
            return;
        }
        else{
            S32 adr_tc_loc_tmp = tc_loc_tmp->adr_tc_;
            S32 adr_tc_glb_tmp = tc_glb_tmp->adr_tc_;
            for(S32 i=0; i<N_CHILDREN; i++){
                MakeIPGroupUseGLBTreeLong
		  (ipg_first, tc_loc_first, tc_glb_first, epi_first, 
		   adr_tc_loc_tmp+i, adr_tc_glb_tmp+i, n_grp_limit, n_leaf_limit);
            }
        }
    }

    template<class Tipg, class Ttc, class Tepi>
    inline void MakeIPGroupShort(ReallocatableArray<Tipg> & ipg_first,
                                 const ReallocatableArray<Ttc> & tc_first,
                                 const ReallocatableArray<Tepi> & epi_first,
                                 const S32 adr_tc,
                                 const S32 n_grp_limit){
        const Ttc * tc_tmp = tc_first.getPointer() + adr_tc;
        const S32 n_tmp = tc_tmp->n_ptcl_;
        if(n_tmp == 0) return;
        else if( tc_tmp->isLeaf(n_grp_limit) ){
            ipg_first.increaseSize();
            ipg_first.back().copyFromTC(*tc_tmp);
            return;
        }
        else{
            S32 adr_tc_tmp = tc_tmp->adr_tc_;
            for(S32 i=0; i<N_CHILDREN; i++){
                MakeIPGroupShort<Tipg, Ttc, Tepi>
                    (ipg_first, tc_first, epi_first, adr_tc_tmp+i, n_grp_limit);
            }
        }
    }


    template<class Tep>
    void CheckMortonSort(const S32 n,
                         TreeParticle tp[],
                         Tep ep[],
                         std::ostream & fout = std::cout){
        S32 err_comp = 0;
        S32 err_order = 0;
        U64 * key = new U64[n];
        bool * flag = new bool[n];
        for(S32 i=0; i<n; i++){
            flag[i] = false;
	    key[i] = MortonKey::getKey(ep[i].getPos());
        }
        for(S32 i=0; i<n; i++) flag[tp[i].adr_ptcl_] = true;
        for(S32 i=0; i<n; i++){
            if(flag[i] == false){
                err_comp++;
            }
        }
        for(S32 i=1; i<n; i++){
            if( key[i] < key[i-1] || tp[i].getKey() < tp[i-1].getKey() ){
                err_order++;
            }
        }
        delete [] flag;
        delete [] key;
        if( err_comp || err_order){
            fout<<"checkMortonSort: FAIL err_comp="<<err_comp<<" err_order="<<err_order<<std::endl;
        }
        else{
            fout<<"checkMortonSort: PASS"<<std::endl;
        }
    }

    template<class Tep, class Tsp>
    void CheckMortonSort(const S32 n,
                         TreeParticle tp[],
                         Tep ep[],
                         Tsp sp[],
                         std::ostream & fout = std::cout){
        S32 err_comp = 0;
        S32 err_order = 0;
        F64 mass_cm_tmp = 0.0;
        U64 * key = new U64[n];
        bool * flag = new bool[n];
        for(S32 i=0; i<n; i++){
            const U32 adr = tp[i].adr_ptcl_;
            if( GetMSB(adr) ){
                key[i] = MortonKey::getKey(sp[i].getPos());
                mass_cm_tmp += sp[i].getCharge();
            }
            else{
                key[i] = MortonKey::getKey(ep[i].getPos());
                mass_cm_tmp += ep[i].getCharge();
            }
            flag[i] = false;
        }
        for(S32 i=0; i<n; i++){
            flag[ClearMSB(tp[i].adr_ptcl_)] = true;
        }
        for(S32 i=0; i<n; i++){
            if(flag[i] == false){
                err_comp++;
            }
        }
        for(S32 i=1; i<n; i++){
            if( key[i] < key[i-1] || tp[i].getKey() < tp[i-1].getKey() ){
                err_order++;
            }
        }
        delete [] flag;
        delete [] key;
        if( err_comp || err_order){
            fout<<"checkMortonSort: FAIL err_comp="<<err_comp<<" err_order="<<err_order<<std::endl;
        }
        else{
            fout<<"checkMortonSort: PASS"<<std::endl;
        }
        fout<<"mass_cm_tmp="<<mass_cm_tmp<<std::endl;
    }



    //////////////////////////
    //// CHECK CALC MOMENT ///
    template<class Tep, class Ttc>
    void CheckCalcMomentLongLocalTree(Tep ep[],
                                      const S32 n_ptcl,
                                      Ttc tc[],
                                      const F64 tolerance = 1e-5,
                                      std::ostream & fout = std::cout){
        F64 mass_cm = 0.0;
        F64vec pos_cm = 0.0;
        for(S32 i=0; i<n_ptcl; i++){
            mass_cm += ep[i].getCharge();
            pos_cm += ep[i].getCharge() * ep[i].getPos();
        }
        pos_cm /= mass_cm;
        F64 dm = std::abs(tc[0].mom_.mass - mass_cm);
        F64vec dx = tc[0].mom_.pos - pos_cm;
        F64 dx_max = dx.applyEach( Abs<F64>() ).getMax();
        if( dm >= tolerance || dx_max >= tolerance){
            fout<<"CalcMomentLong test: FAIL"<<std::endl;
        }
        else{
            fout<<"CalcMomentLong test: PASS"<<std::endl;
        }
        if(Comm::getRank() == 0){
            fout<<"mass_cm="<<mass_cm<<", tc[0].mom_.mass="<<tc[0].mom_.mass<<std::endl;
            fout<<"pos_cm="<<pos_cm<<", tc[0].mom_.pos="<<tc[0].mom_.pos<<std::endl;
        }
    }

    template<class Tep, class Tsp, class Ttc>
    void CheckCalcMomentLongGlobalTree(const Tep ep[],
                                       const Tsp sp[],
                                       const S32 n_ptcl,
                                       const Ttc tc[],
                                       const TreeParticle tp[],
                                       const F64 tolerance = 1e-5,
                                       std::ostream & fout = std::cout){
        F64 mass_cm_glb = 0.0;
        F64vec pos_cm_glb = 0.0;
        for(S32 i=0; i<n_ptcl; i++){
            if( GetMSB(tp[i].adr_ptcl_) == 0){
                mass_cm_glb += ep[i].getCharge();
                pos_cm_glb += ep[i].getCharge() * ep[i].getPos();
            }
            else{
                mass_cm_glb += sp[i].getCharge();
                pos_cm_glb += sp[i].getCharge() * sp[i].getPos();
            }
        }
        pos_cm_glb /= mass_cm_glb;
        F64 dm = std::abs(tc[0].mom_.mass - mass_cm_glb);
        F64vec dx = tc[0].mom_.pos - pos_cm_glb;
        F64 dx_max = dx.applyEach( Abs<F64>() ).getMax();
        if( dm >= tolerance || dx_max >= tolerance){
            fout<<"CheckCalcMomentLongGlobalTree: FAIL"<<std::endl;
        }
        else{
            fout<<"CheckCalcMomentLongGlobalTree: PASS"<<std::endl;
        }
        if(Comm::getRank() == 0){
            fout<<"mass_cm_glb="<<mass_cm_glb<<", tc[0].mom_.mass="<<tc[0].mom_.mass<<std::endl;
            fout<<"pos_cm_glb="<<pos_cm_glb<<", tc[0].mom_.pos="<<tc[0].mom_.pos<<std::endl;
        }
    }

    template<class Tep, class Ttc>
    void CheckCalcMomentShortInAndOut(Tep ep[],
                                      const S32 n_ptcl,
                                      Ttc tc[],
                                      const F64 tolerance = 1e-5,
                                      std::ostream & fout = std::cout){

        F64ort vertex_out_tmp;
        F64ort vertex_in_tmp;
        vertex_out_tmp.init();
        vertex_in_tmp.init();
        for(S32 i=0; i<n_ptcl; i++){
            vertex_out_tmp.merge(ep[i].getPos(), ep[i].getRSearch());
            vertex_in_tmp.merge(ep[i].getPos());
        }
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
        else{
            fout<<"CalcMomentShort test: PASS"<<std::endl;
        }
    }

    template<class Tep, class Ttc>
    void CheckCalcMomentShortInOnly(Tep ep[], 
                                    const S32 n_ptcl, 
                                    Ttc tc[], 
                                    const F64 tolerance = 1e-5, 
                                    std::ostream & fout = std::cout){
        F64ort vertex_in_tmp;
        vertex_in_tmp.init();
        for(S32 i=0; i<n_ptcl; i++){
            vertex_in_tmp.merge(ep[i].getPos());
        }
        F64 d_in_high = (tc[0].mom_.vertex_in_.high_ - vertex_in_tmp.high_).applyEach( Abs<F64>() ).getMax();
        F64 d_in_low = (tc[0].mom_.vertex_in_.low_ - vertex_in_tmp.low_).applyEach( Abs<F64>() ).getMax();
        if( d_in_high >= tolerance || d_in_low >= tolerance){
            fout<<"CalcMomentShort test: FAIL"<<std::endl;
            fout<<"vertex_in_tmp.high_="<<vertex_in_tmp.high_<<" tc[0].mom_.vertex_in_.high_="<<tc[0].mom_.vertex_in_.high_<<std::endl;
            fout<<"vertex_in_tmp.low_="<<vertex_in_tmp.low_<<" tc[0].mom_.vertex_in_.low_="<<tc[0].mom_.vertex_in_.low_<<std::endl;
        }
        else{
            fout<<"CalcMomentShort test: PASS"<<std::endl;
        }
    }

    template<class Tep, class Ttc>
    void CheckCalcMomentOutOnly(Tep ep[],
				const S32 n_ptcl,
				Ttc tc[],
				const F64 tolerance = 1e-5,
				std::ostream & fout = std::cout){
        F64ort vertex_out_tmp;
        vertex_out_tmp.init();
        for(S32 i=0; i<n_ptcl; i++){
            vertex_out_tmp.merge(ep[i].getPos(), ep[i].getRSearch());
        }
        F64 d_out_high = (tc[0].mom_.vertex_out_.high_ - vertex_out_tmp.high_).applyEach( Abs<F64>() ).getMax();
        F64 d_out_low = (tc[0].mom_.vertex_out_.low_ - vertex_out_tmp.low_).applyEach( Abs<F64>() ).getMax();
        if( d_out_high >= tolerance || d_out_low >= tolerance ){
            fout<<"CalcMomentShort test: FAIL"<<std::endl;
            fout<<"vertex_out_tmp.high_="<<vertex_out_tmp.high_<<" tc[0].mom_.vertex_out_.high_="<<tc[0].mom_.vertex_out_.high_<<std::endl;
            fout<<"vertex_out_tmp.low_="<<vertex_out_tmp.low_<<" tc[0].mom_.vertex_out_.low_="<<tc[0].mom_.vertex_out_.low_<<std::endl;
        }
        else{
            fout<<"CalcMomentOut test: PASS"<<std::endl;
        }
    }

    template<class Ttc, class Tep>
    inline void SearchSendParticleLong(const ReallocatableArray<Ttc> & tc_first,
                                       const S32 adr_tc,
                                       const ReallocatableArray<Tep> & ep_first,
                                       ReallocatableArray<S32> & id_ep_send,
                                       ReallocatableArray<S32> & id_sp_send,
                                       const F64ort & pos_target_box, // position of domain
                                       const F64 r_crit_sq,
                                       const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            //Ttc * tc_child = tc_first + adr_tc_child;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    SearchSendParticleLong<Ttc, Tep>
                        (tc_first,   tc_first[adr_tc_child].adr_tc_, ep_first,
                         id_ep_send, id_sp_send, pos_target_box,
                         r_crit_sq*0.25, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    //id_ep_send.reserve( id_ep_send.size() + n_child );
                    id_ep_send.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        id_ep_send.pushBackNoCheck(adr_ptcl_tmp++);
                    }
                }
            }
            else{
                id_sp_send.push_back(adr_tc_child);
            }
        }
    }

#if 1
    // original version
    // FOR P^3T
    template<class Ttc, class Tep>
    inline void SearchSendParticleLongScatter(const ReallocatableArray<Ttc> & tc_first,
                                              const S32 adr_tc,
                                              const ReallocatableArray<Tep> & ep_first,
                                              ReallocatableArray<S32> & id_ep_send,
                                              ReallocatableArray<S32> & id_sp_send,
                                              const F64ort & pos_target_box, // position of domain
                                              const F64 r_crit_sq,
                                              const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) 
                             || (pos_target_box.overlapped( tc_first[adr_tc+i].mom_.getVertexOut())) ) 
                           << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    SearchSendParticleLongScatter<Ttc, Tep>
                        (tc_first,   tc_first[adr_tc_child].adr_tc_, ep_first,
                         id_ep_send, id_sp_send, pos_target_box,
                         r_crit_sq*0.25, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    id_ep_send.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        id_ep_send.pushBackNoCheck(adr_ptcl_tmp++);
                    }
                }
            }
            else{
                id_sp_send.push_back(adr_tc_child);
            }
        }
    }
#else
    // cut-off verseion
    template<class Ttc, class Tep>
    inline void SearchSendParticleLongScatter(const ReallocatableArray<Ttc> & tc_first,
                                              const S32 adr_tc,
                                              const ReallocatableArray<Tep> & ep_first,
                                              ReallocatableArray<S32> & id_ep_send,
                                              ReallocatableArray<S32> & id_sp_send,
                                              const F64ort & pos_target_box, // position of domain
                                              const F64 r_crit_sq,
                                              const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) 
                             || (pos_target_box.overlapped( tc_first[adr_tc+i].mom_.getVertexOut())) ) 
                           << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    SearchSendParticleLongScatter<Ttc, Tep>
                        (tc_first,   tc_first[adr_tc_child].adr_tc_, ep_first,
                         id_ep_send, id_sp_send, pos_target_box,
                         r_crit_sq*0.25, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    id_ep_send.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        id_ep_send.pushBackNoCheck(adr_ptcl_tmp++);
                    }
                }
            }
            else{
                id_sp_send.push_back(adr_tc_child);
            }
        }
    }
#endif

    inline F64ort makeChildCellBox(const S32 i,
                                   const F64ort & cell_box) {
        const F64 half_length = 0.5 * (cell_box.high_.x - cell_box.low_.x);
        F64vec child_cell_high, child_cell_low;
        for(S32 k = 0; k < DIMENSION; k++) {
            child_cell_high[k] = cell_box.high_[k]
                - (((i >> (DIMENSION - 1 - k)) & 0x1) ^ 0x1) * half_length;
            child_cell_low[k]  = cell_box.low_[k]
                + ((i >> (DIMENSION - 1 - k)) & 0x1) * half_length;
        }
        //const F64ort child_cell_box(child_cell_high, child_cell_low);
        const F64ort child_cell_box(child_cell_low, child_cell_high);
        return child_cell_box;
    }

    template<class Ttc, class Tep, class Tsp>
    inline void SearchSendParticleLongCutoff
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<Tep> & ep_send,
     ReallocatableArray<Tsp> & sp_send,
     const F64ort & cell_box,
     const F64ort & pos_target_domain,
     const F64 r_crit_sq,
     const F64 r_cut_sq,
     const S32 n_leaf_limit, 
     const F64vec & shift = 0.){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( (pos_target_domain.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i); // using opening criterion
#if 0
            // vertex_out is not used. 
            const F64ort child_cell_box = makeChildCellBox(i, cell_box);
            open_bits |= ( (pos_target_domain.getDistanceMinSQ(child_cell_box) <= r_cut_sq) << (i + N_CHILDREN) ); // using cutoff criterion
#else
            // cell_box is not neede
            open_bits |= (pos_target_domain.overlapped( tc_first[adr_tc+i].mom_.getVertexOut() ) << (i + N_CHILDREN) );
#endif
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if( ( (open_bits >> (i+N_CHILDREN)) & 0x1 ) ^ 0x1 ) continue;
            const S32 adr_tc_child = adr_tc + i;
            Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    const F64ort child_cell_box = makeChildCellBox(i, cell_box);
                    SearchSendParticleLongCutoff<Ttc, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_,
                         ep_first, ep_send, sp_send, child_cell_box,
                         pos_target_domain, r_crit_sq*0.25, r_cut_sq, n_leaf_limit,
                         shift);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    //ep_send.reserve( ep_send.size()+n_child );
                    ep_send.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        ep_send.pushBackNoCheck(ep_first[adr_ptcl_tmp++]);
                        const F64vec pos_new = ep_send.back().getPos() + shift;
                        ep_send.back().setPos(pos_new);
                    }
                }
            }
            else{
                //sp_send.resizeNoInitialize(sp_send.size()+1);
                sp_send.increaseSize();
                sp_send.back().copyFromMoment(tc_child->mom_);
                const F64vec pos_new = sp_send.back().getPos() + shift;
                sp_send.back().setPos(pos_new);
            }
        }
    }
    
    /////////////////////////////
    /// MAKE INTERACTION LIST ///
#if 1
    template<class Ttc, class Ttp, class Tep, class Tsp>
    void MakeInteractionListLongEPSP(const ReallocatableArray<Ttc> & tc_first,
                                     const S32 adr_tc,
                                     const ReallocatableArray<Ttp> & tp_first,
                                     const ReallocatableArray<Tep> & ep_first,
                                     ReallocatableArray<Tep> & ep_list,
                                     const ReallocatableArray<Tsp> & sp_first,
                                     ReallocatableArray<Tsp> & sp_list,
                                     const F64ort & pos_target_box,
                                     const F64 r_crit_sq,
                                     const S32 n_leaf_limit){
        U32 open_bits = 0;
//#ifdef __HPC_ACE__
#ifdef FAST_WALK_K
        static const _fjsp_v2r8 ZERO = _fjsp_setzero_v2r8();
        const F64 x_h = pos_target_box.high_.x;
        const F64 x_l = pos_target_box.low_.x;
        const F64 y_h = pos_target_box.high_.y;
        const F64 y_l = pos_target_box.low_.y;
        const F64 z_h = pos_target_box.high_.z;
        const F64 z_l = pos_target_box.low_.z;
        const F64 cen_x = (x_h + x_l) * 0.5;
        const F64 cen_y = (y_h + y_l) * 0.5;
        const F64 cen_z = (z_h + z_l) * 0.5;
        const F64 len_x = (x_h - x_l) * 0.5;
        const F64 len_y = (y_h - y_l) * 0.5;
        const F64 len_z = (z_h - z_l) * 0.5;
        const _fjsp_v2r8 ix = _fjsp_set_v2r8(cen_x, cen_x);
        const _fjsp_v2r8 iy = _fjsp_set_v2r8(cen_y, cen_y);
        const _fjsp_v2r8 iz = _fjsp_set_v2r8(cen_z, cen_z);
        const _fjsp_v2r8 lx = _fjsp_set_v2r8(len_x, len_x);
        const _fjsp_v2r8 ly = _fjsp_set_v2r8(len_y, len_y);
        const _fjsp_v2r8 lz = _fjsp_set_v2r8(len_z, len_z);
        F64 dr_sq[N_CHILDREN];
        for(S32 i=0; i<N_CHILDREN/2; i++){
            const S32 i0 = i*2;
            const S32 i1 = i*2+1;
            const F64vec & mom0 = tc_first[adr_tc+i0].mom_.getPos();
            const F64  mom0_x = mom0.x;
            const F64  mom0_y = mom0.y;
            const F64  mom0_z = mom0.z;
            const F64vec & mom1 = tc_first[adr_tc+i1].mom_.getPos();
            const F64  mom1_x = mom1.x;
            const F64  mom1_y = mom1.y;
            const F64  mom1_z = mom1.z;
            const _fjsp_v2r8 jx = _fjsp_set_v2r8(mom0_x, mom1_x);
            _fjsp_v2r8 dx = _fjsp_abs_v2r8( _fjsp_sub_v2r8(ix, jx) );
            dx = _fjsp_sub_v2r8( dx, lx );
            dx = _fjsp_and_v2r8( dx, _fjsp_cmplt_v2r8(ZERO, dx) );

            const _fjsp_v2r8 jy = _fjsp_set_v2r8(mom0_y, mom1_y);
            _fjsp_v2r8 dy = _fjsp_abs_v2r8( _fjsp_sub_v2r8(iy, jy) );
            dy = _fjsp_sub_v2r8( dy, ly );
            dy = _fjsp_and_v2r8( dy, _fjsp_cmplt_v2r8(ZERO, dy) );

            const _fjsp_v2r8 jz = _fjsp_set_v2r8(mom0_z, mom1_z);
            _fjsp_v2r8 dz = _fjsp_abs_v2r8( _fjsp_sub_v2r8(iz, jz) );
            dz = _fjsp_sub_v2r8( dz, lz );
            dz = _fjsp_and_v2r8( dz, _fjsp_cmplt_v2r8(ZERO, dz) );

            _fjsp_v2r8 dr_sq_tmp = _fjsp_madd_v2r8(dz, dz,
                                               _fjsp_madd_v2r8(dy, dy,
                                                               _fjsp_mul_v2r8(dx, dx)));
            _fjsp_storeh_v2r8(dr_sq+i0, dr_sq_tmp);
            _fjsp_storel_v2r8(dr_sq+i1, dr_sq_tmp);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (( dr_sq[i] <= r_crit_sq) << i);
        }
#else //FAST_WALK_K
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i);
        }
#endif //FAST_WALK_K
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeInteractionListLongEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         ep_list,   sp_first, sp_list, 
                         pos_target_box, r_crit_sq*0.25, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    sp_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_) == 0){
                            ep_list.pushBackNoCheck(ep_first[adr_ptcl_tmp++]);
                        }
                        else{
                            sp_list.pushBackNoCheck(sp_first[adr_ptcl_tmp++]);
                        }
                    }
                }
            }
            else{
                sp_list.increaseSize();
                sp_list.back().copyFromMoment(tc_child->mom_);
            }
        }
    }
#else
    template<class Ttc, class Ttp, class Tep, class Tsp>
    void MakeInteractionListLongEPSP(const ReallocatableArray<Ttc> & tc_first,
                                     const S32 adr_tc,
                                     const ReallocatableArray<Ttp> & tp_first,
                                     const ReallocatableArray<Tep> & ep_first,
                                     ReallocatableArray<Tep> & ep_list,
                                     const ReallocatableArray<Tsp> & sp_first,
                                     ReallocatableArray<Tsp> & sp_list,
                                     const F64ort & pos_target_box,
                                     const F64 r_crit_sq,
                                     const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            //const Ttc * tc_child = tc_first.data() + adr_tc_child;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeInteractionListLongEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         ep_list,   sp_first, sp_list, 
                         pos_target_box, r_crit_sq*0.25, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    //ep_list.reserve( ep_list.size()+n_child );
                    //sp_list.reserve( sp_list.size()+n_child );
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    sp_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_) == 0){
                            //ep_list[n_ep] = ep_first[adr_ptcl_tmp++];
                            ep_list.pushBackNoCheck(ep_first[adr_ptcl_tmp++]);
                        }
                        else{
                            sp_list.pushBackNoCheck(sp_first[adr_ptcl_tmp++]);
                        }
                    }
                }
            }
            else{
                sp_list.increaseSize();
                sp_list.back().copyFromMoment(tc_child->mom_);
            }
        }
    }
#endif
    template<class Ttc, class Ttp, class Tep, class Tsp>
    void MakeInteractionListLongScatterEPSP(const ReallocatableArray<Ttc> & tc_first,
                                            const S32 adr_tc,
                                            const ReallocatableArray<Ttp> & tp_first,
                                            const ReallocatableArray<Tep> & ep_first,
                                            ReallocatableArray<Tep> & ep_list,
                                            const ReallocatableArray<Tsp> & sp_first,
                                            ReallocatableArray<Tsp> & sp_list,
                                            const F64ort & pos_target_box,
                                            const F64 r_crit_sq,
                                            const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            const F64ort box_tmp = tc_first[adr_tc+i].mom_.getVertexOut();
            //open_bits |= ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i);
            open_bits |= ( 
                ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) 
                  || ( pos_target_box.overlapped(box_tmp) ) )
                << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeInteractionListLongScatterEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         ep_list,   sp_first, sp_list, 
                         pos_target_box, r_crit_sq*0.25, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    sp_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_) == 0){
                            ep_list.pushBackNoCheck(ep_first[adr_ptcl_tmp++]);
                        }
                        else{
                            sp_list.pushBackNoCheck(sp_first[adr_ptcl_tmp++]);
                        }
                    }
                }
            }
            else{
                sp_list.increaseSize();
                sp_list.back().copyFromMoment(tc_child->mom_);
            }
        }
    }


    template<class Ttc, class Ttp, class Tep, class Tsp>
    void MakeInteractionListLongCutoffEPSP(const ReallocatableArray<Ttc> & tc_first,
					   const S32 adr_tc,
					   const ReallocatableArray<Ttp> & tp_first,
					   const ReallocatableArray<Tep> & ep_first,
					   ReallocatableArray<Tep> & ep_list,
					   const ReallocatableArray<Tsp> & sp_first,
					   ReallocatableArray<Tsp> & sp_list,
					   const F64ort & cell_box,
					   const F64ort & pos_target_box,
					   const F64 r_crit_sq,
					   const F64 r_cut_sq,
					   const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i);
#if 1
	    // vertex_out is not used
            const F64ort child_cell_box = makeChildCellBox(i, cell_box);
            open_bits |= ( (pos_target_box.getDistanceMinSQ(child_cell_box) <= r_cut_sq) << (i + N_CHILDREN));
#else
	    // cell_box is not used. maybe faster than above one.
	    // In this case, SPJ must have r_cutoff
	    // becuase vertex_out of global tree depends on r_cutoff of SPJ.
	    // It is not easy to tell r_cutoff to SPJ defined by users.
	    open_bits |= (pos_target_box.overlapped( tc_first[adr_tc+i].mom_.getVertexOut() ) << (i + N_CHILDREN) );
#endif
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if( ( (open_bits >> (i+N_CHILDREN)) & 0x1 ) ^ 0x1 ) continue;
            const S32 adr_tc_child = adr_tc + i;
            Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    const F64ort child_cell_box = makeChildCellBox(i, cell_box);
                    MakeInteractionListLongCutoffEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first,
                         ep_first, ep_list, sp_first, sp_list, child_cell_box,
                         pos_target_box, r_crit_sq*0.25, r_cut_sq, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    //ep_list.reserve(ep_list.size()+n_child);
                    //sp_list.reserve(sp_list.size()+n_child);
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    sp_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_) == 0){
                            ep_list.pushBackNoCheck(ep_first[adr_ptcl_tmp++]);
                        }
                        else{
                            sp_list.pushBackNoCheck(sp_first[adr_ptcl_tmp++]);
                        }
                    }
                }
            }
            else{
                //sp_list.resizeNoInitialize( sp_list.size()+1 );
                sp_list.increaseSize();
                sp_list.back().copyFromMoment(tc_child->mom_);
            }
        }
    }


    template<class Ttc, class Tep2, class Tep3>
    inline void MakeListUsingOuterBoundary(const Ttc * tc_first,
                                           const S32 adr_tc,
                                           const Tep2 * ep_first,
                                           ReallocatableArray<Tep3> & ep_list,
                                           const F64ort & pos_target_box, // position of domain
                                           const S32 n_leaf_limit,
                                           const F64vec & shift = F64vec(0.0) ){
        U32 open_bits = 0;
        //std::cout<<"CHECK 0 pos_target_box="<<pos_target_box<<std::endl;
        //std::cout<<"CHECK 0 adr_tc="<<adr_tc<<std::endl;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (pos_target_box.overlapped( tc_first[adr_tc+i].mom_.getVertexOut() ) << i);
            //std::cout<<"CHECK 0 tc_first[adr_tc+i].mom_.getVertexOut()="<<tc_first[adr_tc+i].mom_.getVertexOut()<<std::endl;
        }

        //std::cout<<"CHECK 0 open_bits="<<open_bits<<std::endl;
        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                //std::cout<<"CHECK A i="<<i<<std::endl;
                const S32 adr_tc_child = adr_tc + i;
                const Ttc * tc_child = tc_first + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    //std::cout<<"CHECK B n_child="<<n_child<<std::endl;
                    MakeListUsingOuterBoundary<Ttc, Tep2, Tep3>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, ep_list,
                         pos_target_box, n_leaf_limit, shift);
                }
                else{
                    //std::cout<<"CHECK C tc_child->adr_ptcl_="<<tc_child->adr_ptcl_<<std::endl;
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    //ep_list.reserve( ep_list.size()+n_child );
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
#ifdef ORIGINAL_SCATTER_MODE
                        const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                        const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq_tmp = pos_target_box.getDistanceMinSQ(pos_tmp);
                        //std::cout<<"CHECK D pos_tmp="<<pos_tmp<<std::endl;
                        //std::cout<<"CHECK D size_tmp="<<size_tmp<<std::endl;
                        //std::cout<<"CHECK D pos_target_box="<<pos_target_box<<std::endl;
                        //std::cout<<"CHECK D dis_sq_tmp="<<dis_sq_tmp<<std::endl;
                        if(dis_sq_tmp > size_tmp*size_tmp) continue;
#endif
                        //std::cout<<"CHECK E"<<std::endl;
                        ep_list.increaseSize();
                        ep_list.back() = ep_first[adr_ptcl_tmp];
                        const F64vec pos_new = ep_list.back().getPos() + shift;
                        ep_list.back().setPos(pos_new);
                    }
                }
            }
        }
    }


/*
    template<class Ttc, class Tep2, class Tep3, class Tep4>
    inline void MakeListUsingOuterBoundaryForGather(const Ttc * tc_first,
                                                    const S32 adr_tc,
                                                    const Tep2 * ep_first,
                                                    ReallocatableArray<Tep3> & ep_list,
                                                    ReallocatableArray<EPXROnly> & ep_x_r_list,
                                                    const F64ort & pos_target_box, // position of domain
                                                    const S32 n_leaf_limit,
                                                    const Tep4 * epi_first,
                                                    const F64vec & shift = F64vec(0.0) ){
        U32 open_bits = 0;
        //std::cout<<"CHECK 0 pos_target_box="<<pos_target_box<<std::endl;
        //std::cout<<"CHECK 0 adr_tc="<<adr_tc<<std::endl;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (pos_target_box.overlapped( tc_first[adr_tc+i].mom_.getVertexOut() ) << i);
            //std::cout<<"CHECK 0 tc_first[adr_tc+i].mom_.getVertexOut()="<<tc_first[adr_tc+i].mom_.getVertexOut()<<std::endl;
        }

        //std::cout<<"CHECK 0 open_bits="<<open_bits<<std::endl;
        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                //std::cout<<"CHECK A i="<<i<<std::endl;
                const S32 adr_tc_child = adr_tc + i;
                const Ttc * tc_child = tc_first + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    //std::cout<<"CHECK B n_child="<<n_child<<std::endl;
                    MakeListUsingOuterBoundaryForGather<Ttc, Tep2, Tep3>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, ep_list, ep_x_r_list,
                         pos_target_box, n_leaf_limit, epi_first, shift);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    ep_x_r_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                        const F64vec pos_tmp = epi_first[adr_ptcl_tmp].getPos();
                        const F64 size_tmp = epi_first[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq_tmp = pos_target_box.getDistanceMinSQ(pos_tmp);
                        //std::cout<<"CHECK D pos_tmp="<<pos_tmp<<std::endl;
                        //std::cout<<"CHECK D size_tmp="<<size_tmp<<std::endl;
                        //std::cout<<"CHECK D pos_target_box="<<pos_target_box<<std::endl;
                        //std::cout<<"CHECK D dis_sq_tmp="<<dis_sq_tmp<<std::endl;
                        if(dis_sq_tmp > size_tmp*size_tmp) continue;
                        //std::cout<<"CHECK E"<<std::endl;
                        ep_list.increaseSize();
                        ep_list.back() = ep_first[adr_ptcl_tmp];
                        const F64vec pos_new = ep_list.back().getPos() + shift;
                        ep_list.back().setPos(pos_new);
                        ep_x_r_list.increaseSize();
                        ep_x_r_list.back().copyFromEP(epi_first[adr_ptcl_tmp]);
                        ep_x_r_list.back().setPos(pos_new);
                    }
                }
            }
        }
    }
*/

    // NEW functions for 
    template<class Ttc, class Tep2, class Tep3>
    inline void MakeListUsingOuterBoundaryAndInnerBoundary(const Ttc * tc_first,
                                                           const S32 adr_tc,
                                                           const Tep2 * ep_first,
                                                           ReallocatableArray<Tep3> & ep_list,
                                                           const F64ort & pos_target_box_out, // position of domain
                                                           const F64ort & pos_target_box_in, // position of domain
                                                           const S32 n_leaf_limit,
                                                           const F64vec & shift = F64vec(0.0) ){
        U32 open_bits = 0;
        asm("# MakeListUsingOuterBoundaryAndInnerBoundary");
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= ( (pos_target_box_out.overlapped( tc_first[adr_tc+i].mom_.getVertexIn() ) 
                            || pos_target_box_in.overlapped( tc_first[adr_tc+i].mom_.getVertexOut() ) )
                           << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                const S32 adr_tc_child = adr_tc + i;
                const Ttc * tc_child = tc_first + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeListUsingOuterBoundaryAndInnerBoundary<Ttc, Tep2, Tep3>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, ep_list,
                         pos_target_box_out, pos_target_box_in, n_leaf_limit, shift);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    //ep_list.reserve( ep_list.size()+n_child );
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                        const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                        const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq_tmp = pos_target_box_in.getDistanceMinSQ(pos_tmp);
                        if( pos_target_box_out.notOverlapped(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                        ep_list.increaseSize();
                        ep_list.back() = ep_first[adr_ptcl_tmp];
                        const F64vec pos_new = ep_list.back().getPos() + shift;
                        ep_list.back().setPos(pos_new);
                    }
                }
            }
        }
    }

    template<class Ttc, class Tep2, class Tep3>
    inline void MakeListUsingInnerBoundary(const Ttc * tc_first,
                                           const S32 adr_tc,
                                           const Tep2 * ep_first,
                                           ReallocatableArray<Tep3> & ep_list,
                                           const F64ort & pos_target_box, // position of domain
                                           const S32 n_leaf_limit,
                                           const F64vec & shift = F64vec(0.0) ){
        U32 open_bits = 0;
        asm("# MakeListUsingInnerBoundary");
#if 0
//const F64ort vertex0 = 
        __builtin_v2r8 top_x, top_y, top_z, bot_x, bot_y, bot_z;
        __builtin_fj_storel_v2r8(&pos_target_box.low_.x, bot_x);
        __builtin_fj_storeh_v2r8(&pos_target_box.low_.x, bot_x);
        __builtin_fj_storel_v2r8(&pos_target_box.low_.y, bot_y);
        __builtin_fj_storeh_v2r8(&pos_target_box.low_.y, bot_y);
        __builtin_fj_storel_v2r8(&pos_target_box.low_.z, bot_z);
        __builtin_fj_storeh_v2r8(&pos_target_box.low_.z, bot_z);
        __builtin_fj_storel_v2r8(&pos_target_box.high_.x, top_x);
        __builtin_fj_storeh_v2r8(&pos_target_box.high_.x, top_x);
        __builtin_fj_storel_v2r8(&pos_target_box.high_.y, top_y);
        __builtin_fj_storeh_v2r8(&pos_target_box.high_.y, top_y);
        __builtin_fj_storel_v2r8(&pos_target_box.high_.z, top_z);
        __builtin_fj_storeh_v2r8(&pos_target_box.high_.z, top_z);
        __builtin_v2r8 top_x0, top_y0, top_z0, bot_x0, bot_y0, bot_z0;
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+0].mom_.getVertexIn().low_.x, bot_x0);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+0].mom_.getVertexIn().low_.y, bot_y0);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+0].mom_.getVertexIn().low_.z, bot_z0);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+0].mom_.getVertexIn().high_.x, top_x0);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+0].mom_.getVertexIn().high_.y, top_y0);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+0].mom_.getVertexIn().high_.z, top_z0);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+1].mom_.getVertexIn().low_.x, bot_x1);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+1].mom_.getVertexIn().low_.y, bot_y1);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+1].mom_.getVertexIn().low_.z, bot_z1);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+1].mom_.getVertexIn().high_.x, top_x1);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+1].mom_.getVertexIn().high_.y, top_y1);
        __builtin_fj_storel_v2r8(&tc_first[adr_tc+1].mom_.getVertexIn().high_.z, top_z1);

#else
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (pos_target_box.overlapped( tc_first[adr_tc+i].mom_.getVertexIn() ) << i);
        }
#endif
        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                const S32 adr_tc_child = adr_tc + i;
                const Ttc * tc_child = tc_first + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeListUsingInnerBoundary<Ttc, Tep2, Tep3>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, ep_list,
                         pos_target_box, n_leaf_limit, shift);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    //ep_list.reserve( ep_list.size()+n_child );
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        ep_list.pushBackNoCheck(ep_first[adr_ptcl_tmp]);
                        const F64vec pos_new = ep_first[adr_ptcl_tmp++].getPos() + shift; // for periodic mode
                        ep_list.back().setPos(pos_new);
                    }
                }
            }
         }
    }

/*
    template<class Ttc, class Tep>
    inline void MakeListUsingInnerBoundaryForGatherModeNormalMode
    (const Ttc * tc_first,
     const S32 adr_tc,
     const Tep * ep_first,
     ReallocatableArray<S32> & id_send,
     const F64ort & pos_box_target,
     const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (tc_first[adr_tc+i].mom_.getVertexIn().overlapped(pos_box_target) << i);
        }
         for(S32 i=0; i<N_CHILDREN; i++){
             if( (open_bits>>i) & 0x1){
                 const S32 adr_tc_child = adr_tc + i;
                 const Ttc * tc_child = tc_first + adr_tc_child;
                 const S32 n_child = tc_child->n_ptcl_;
                 if(n_child == 0) continue;
                 else{
                     if( !(tc_child->isLeaf(n_leaf_limit)) ){
                         MakeListUsingInnerBoundaryForGatherModeNormalMode<Ttc, Tep>
                             (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, id_send,
                              pos_box_target, n_leaf_limit);
                     }
                     else{
                         //id_send.reserve( id_send.size()+n_child );
                         id_send.reserveEmptyAreaAtLeast( n_child );
                         S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                         for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                             const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                             if(pos_box_target.overlapped(pos_tmp)){
                                 id_send.pushBackNoCheck(adr_ptcl_tmp);
                             }
                         }
                     }
                 }
             }
         }
    }
*/
/*
    template<class Ttc, class Tep>
    inline void MakeListUsingInnerBoundaryForSymmetryExclusive
    (const Ttc * tc_first,
     const S32 adr_tc,
     const Tep * ep_first,
     ReallocatableArray<S32> & id_send,
     const F64ort & pos_box_target,
     const F64ort & pos_domain,
     const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (tc_first[adr_tc+i].mom_.getVertexIn().overlapped(pos_box_target) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                const S32 adr_tc_child = adr_tc + i;
                const Ttc * tc_child = tc_first + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeListUsingInnerBoundaryForSymmetryExclusive<Ttc, Tep>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, id_send, 
                         pos_box_target, pos_domain, n_leaf_limit);
                }
                else{
                    //id_send.reserve( id_send.size()+n_child );
                    id_send.reserveEmptyAreaAtLeast( n_child );
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                        // NOTE: need to be concistent with MakeListUsingOuterBoundary()
                        const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                        
                        const F64 len_sq = ep_first[adr_ptcl_tmp].getRSearch() * ep_first[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq0 = pos_box_target.getDistanceMinSQ(pos_tmp);
                        const F64 dis_sq1 = pos_domain.getDistanceMinSQ(pos_tmp);
                        if(dis_sq0 <= len_sq && dis_sq1 > len_sq){ id_send.pushBackNoCheck(adr_ptcl_tmp); }
                    }
                }
            }
        }
    }
*/
/*
    template<class Ttc, class Tep>
    inline void MakeListUsingInnerBoundaryForSymmetryExclusiveSetFlag
    (const Ttc * tc_first,
     const S32 adr_tc,
     const Tep * ep_first,
     ReallocatableArray<bool> & flag_send,
     const F64ort & pos_box_target,
     const F64ort & pos_domain,
     const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (tc_first[adr_tc+i].mom_.getVertexIn().overlapped(pos_box_target) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                const S32 adr_tc_child = adr_tc + i;
                const Ttc * tc_child = tc_first + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeListUsingInnerBoundaryForSymmetryExclusiveSetFlag<Ttc, Tep>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, flag_send, 
                         pos_box_target, pos_domain, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                        // NOTE: need to be concistent with MakeListUsingOuterBoundary()
                        const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                        const F64 len_sq = ep_first[adr_ptcl_tmp].getRSearch() * ep_first[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq0 = pos_box_target.getDistanceMinSQ(pos_tmp);
                        const F64 dis_sq1 = pos_domain.getDistanceMinSQ(pos_tmp);
                        if(dis_sq0 <= len_sq && dis_sq1 > len_sq){
                            flag_send[adr_ptcl_tmp] = true;
                        }
                    }
                }
            }
        }
    }
*/
    template<class Ttc>
    inline void IsOverlappedUsingOuterBoundaryOfTreeImpl(const Ttc * tc_first,
							 const S32 adr_tc,
							 const F64ort & pos_box,
							 const S32 n_leaf_limit,
							 U32 & is_overlapped){

        if(is_overlapped) return;
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (tc_first[adr_tc+i].mom_.getVertexOut().overlapped(pos_box) << i);
        }

        for(S32 i=0; i<N_CHILDREN; i++){
            if(is_overlapped) return;
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first + adr_tc_child;
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            if( (open_bits>>i) & 0x1){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    IsOverlappedUsingOuterBoundaryOfTreeImpl(tc_first, tc_first[adr_tc_child].adr_tc_,
                                                             pos_box, n_leaf_limit, is_overlapped);
                }
                else{
                    // leaf and overlapped
                    is_overlapped = 1;
                    return;
                }
            }
        }
    }

    template<class Ttc>
    inline U32 IsOverlappedUsingOuterBoundaryOfTree(const Ttc * tc_first,
						    const F64ort & pos_box,
						    const S32 n_leaf_limit){
        U32 is_overlapped = 0;
        S32 adr_tc = N_CHILDREN;
        if( !tc_first->isLeaf(n_leaf_limit) ){
            IsOverlappedUsingOuterBoundaryOfTreeImpl(tc_first, adr_tc, pos_box, n_leaf_limit, is_overlapped);
        }
        else{
            is_overlapped = tc_first->mom_.getVertexOut().overlapped(pos_box);
        }
        return is_overlapped;
    }

    template<class Ttca, class Ttcb, class Tepb>
    inline void MakeLETListByDoubleWalkImpl(const Ttca * tc_first_A,
					    const Ttcb * tc_first_B,
					    const S32 adr_tc_B,
					    const Tepb * ep_first_B,
					    const F64ort & pos_domain,
					    const S32 n_leaf_limit_A,
					    const S32 n_leaf_limit_B,
					    ReallocatableArray<S32> & id_ptcl_send){
        U32 open_bits = 0;
        const U32 ONE = 1;
        //asm("# MakeLETListByDoubleWalkImpl");
        for(S32 i=0; i<N_CHILDREN; i++){
            if( tc_first_B[adr_tc_B+i].mom_.getVertexOut().overlapped(pos_domain) 
                || IsOverlappedUsingOuterBoundaryOfTree(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexIn(), n_leaf_limit_A)){
                open_bits |= ONE << i;
            }
/*
            open_bits |= (tc_first_B[adr_tc_B+i].mom_.getVertexOut().overlapped(pos_domain) << i);
            if( (~(open_bits>>i)) & 0x1 ){
                open_bits |= IsOverlappedUsingOuterBoundaryOfTree(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexIn(), n_leaf_limit_A) << i;
                //open_bits |= IsOverlappedUsingOuterBoundaryOfTree(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexOut(), n_leaf_limit_A) << i;
                //open_bits |= tc_first_A->mom_.getVertexOut().overlapped(tc_first_B[adr_tc_B+i].mom_.getVertexIn()) << i; // out
                //open_bits |= tc_first_A->mom_.getVertexOut().overlapped(tc_first_B[adr_tc_B+i].mom_.getVertexOut()) << i; // OK
                //open_bits |= tc_first_B[adr_tc_B+i].mom_.getVertexIn().overlapped(tc_first_A->mom_.getVertexOut()) << i; // out
            }
*/
        }


        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                const S32 adr_tc_child = adr_tc_B + i;
                const Ttcb * tc_child = tc_first_B + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                else if( !(tc_child->isLeaf(n_leaf_limit_B)) ){
                    MakeLETListByDoubleWalkImpl(tc_first_A, tc_first_B, tc_first_B[adr_tc_child].adr_tc_,
                                                ep_first_B, pos_domain, n_leaf_limit_A, n_leaf_limit_B,
                                                id_ptcl_send);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
#ifdef ORIGINAL_SCATTER_MODE
                        // NOTE: need to be concistent with MakeListUsingOuterBoundary()
                        const F64vec pos_tmp = ep_first_B[adr_ptcl_tmp].getPos();
                        const F64 len_sq = ep_first_B[adr_ptcl_tmp].getRSearch() * ep_first_B[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq = pos_domain.getDistanceMinSQ(pos_tmp);
                        if(dis_sq <= len_sq){
                            //std::cerr<<"original version"<<std::endl;
                            continue;
                        }
#else
                        if( pos_domain.overlapped(tc_first_B[adr_tc_child].mom_.getVertexOut()) ){
                            //std::cerr<<"new version"<<std::endl;
                            continue;
                        }
#endif
                        //if( pos_domain.overlapped(tc_first_B[adr_tc_child].mom_.getVertexOut()) ) continue;
                        id_ptcl_send.push_back(adr_ptcl_tmp);
                        //std::cerr<<"set ptcl"<<std::endl;
                    }
                }
            }
        }
    }

    template<class Ttca, class Ttcb, class Tepb>
    inline void MakeLETListByDoubleWalk(const Ttca * tc_first_A,
					const Ttcb * tc_first_B,
					const Tepb * ep_first_B,
					const F64ort & pos_domain,
					const S32 n_leaf_limit_A,
					const S32 n_leaf_limit_B,
					ReallocatableArray<S32> & id_ptcl_send){
        const S32 adr_tc_B = N_CHILDREN;
        if( !tc_first_B->isLeaf(n_leaf_limit_B) ){
            MakeLETListByDoubleWalkImpl(tc_first_A, tc_first_B, adr_tc_B,
                                        ep_first_B, pos_domain, n_leaf_limit_A, n_leaf_limit_B,
                                        id_ptcl_send);
        }
        else{
            if( tc_first_B->mom_.getVertexOut().overlapped(pos_domain) || IsOverlappedUsingOuterBoundaryOfTree(tc_first_A, tc_first_B->mom_.getVertexIn(), n_leaf_limit_A)){
                S32 adr_ptcl_tmp = tc_first_B->adr_ptcl_;
                const S32 n_child = tc_first_B->n_ptcl_;
                for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                    // NOTE: need to be concistent with MakeListUsingOuterBoundary()
#ifdef ORIGINAL_SCATTER_MODE
                    const F64vec pos_tmp = ep_first_B[adr_ptcl_tmp].getPos();
                    const F64 len_sq = ep_first_B[adr_ptcl_tmp].getRSearch() * ep_first_B[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq = pos_domain.getDistanceMinSQ(pos_tmp);
                    if(dis_sq <= len_sq) continue;
#else
                    if( tc_first_B->mom_.getVertexOut().overlapped(pos_domain) ) return;
#endif
                    id_ptcl_send.push_back(adr_ptcl_tmp);
                }
            }
        }
    }

#if 1
// for neighbour search
    template<class Tvec, class Ttc, class Ttp, class Tep>
    inline void SearchNeighborListOneParticleScatter(const Tvec & pos_target,
                                                     const Ttc * tc_first,
                                                     const Ttp * tp_first,
                                                     const S32 adr_tc,
                                                     const ReallocatableArray<Tep> & ep_first,
                                                     ReallocatableArray<Tep> & ep_list,
                                                     const S32 n_leaf_limit,
                                                     bool & error){
	//static S32 err = 0;
        if( tc_first[adr_tc].isLeaf(n_leaf_limit) ){
            const S32 n_tmp = tc_first[adr_tc].n_ptcl_;
            S32 adr_ptcl_tmp = tc_first[adr_tc].adr_ptcl_;
            //ep_list.reserveEmptyAreaAtLeast(n_tmp);
            for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                const F64vec dr = ep_first[adr_ptcl_tmp].getPos() - pos_target;
                const F64 r_search = ep_first[adr_ptcl_tmp].getRSearch();
                //if( dr*dr < r_search*r_search){
                if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_)==0 && dr*dr < r_search*r_search){
                    //ep_list.pushBackNoCheck( ep_first[adr_ptcl_tmp++] );
		    Tep * adr_old = ep_list.getPointer();
                    ep_list.push_back( ep_first[adr_ptcl_tmp] );
		    Tep * adr_new = ep_list.getPointer();
		    if(adr_old != adr_new){
			error = true;
			//std::cerr<<"ep_list.size()="<<ep_list.size()<<std::endl;
		    }
		    /*
		    if(adr_old != adr_new){
			std::cerr<<"adr_new="<<adr_new<<std::endl;
			std::cerr<<"adr_old="<<adr_old<<std::endl;
			err++;
			if(err >= 2){
			    std::cerr<<"adr_new->id="<<adr_new->id<<std::endl;
			    std::cerr<<"adr_old->id="<<adr_old->id<<std::endl;
			}
		    }
		    */


		    //assert(err < 2);
                }
            }
        }
        else{
            const S32 adr_tc_child = tc_first[adr_tc].adr_tc_;
            U32 open_bits = 0;
            for(S32 i=0; i<N_CHILDREN; i++){
                const F64ort vertex_cell = tc_first[adr_tc_child+i].mom_.getVertexOut();
                open_bits |= vertex_cell.overlapped(pos_target) << i;
            }
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_child = tc_first[adr_tc_child+i].n_ptcl_;
                if( n_child == 0 ) continue;
                else if( (open_bits >> i) & 0x1 ){
                    SearchNeighborListOneParticleScatter
                        (pos_target, tc_first, tp_first,    adr_tc_child+i,
                         ep_first,   ep_list,  n_leaf_limit, error);
                }
            }
        }
    }
#endif

#if 0
// for neighbour search
    template<class Tort, class Ttc, class Tep>
    inline void SearchNeighborListOneIPGroupScatter(const Tort & pos_box,
                                                    const Ttc * tc_first,
                                                    const S32 adr_tc,
                                                    const ReallocatableArray<Tep> & ep_first,
                                                    ReallocatableArray<Tep> & ep_list,
                                                    const S32 n_leaf_limit){
        if( tc_first[adr_tc].isLeaf(n_leaf_limit) ){
            const S32 n_tmp = tc_first[adr_tc].n_ptcl_;
            S32 adr_ptcl_tmp = tc_first[adr_tc].adr_ptcl_;
            for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                ep_list.push_back( ep_first[adr_ptcl_tmp] );
            }
        }
        else{
            const S32 adr_tc_child = tc_first[adr_tc].adr_tc_;
            U32 open_bits = 0;
            for(S32 i=0; i<N_CHILDREN; i++){
                const F64ort vertex_cell = tc_first[adr_tc_child+i].mom_.getVertexOut();
                open_bits |= vertex_cell.overlapped(pos_box) << i;
            }
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_child = tc_first[adr_tc_child+i].n_ptcl_;
                if( n_child == 0 ) continue;
                else if( (open_bits >> i) & 0x1 ){
                    SearchNeighborListOneIPGroupScatter
                        (pos_box, tc_first, adr_tc_child+i,
                         ep_first,   ep_list,  n_leaf_limit);
                }
            }
        }
    }
#endif
}
