#define ORIGINAL_SCATTER_MODE

#ifdef __HPC_ACE__
#include<emmintrin.h>
#endif

#ifdef FAST_WALK_AVX2
#include<immintrin.h>

/*
extern __inline __m256 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm256_undefined_ps (void)
{
  __m256 __Y = __Y;
  return __Y;
}
*/

static inline void transpose_4ymm(__m256 &v0, __m256 &v1, __m256 &v2, __m256 &v3){
    // transpose from
    // [ |w4|z4|y4|x4||w0|z0|y0|x0|, |w5|z5|y5|x5||w1|z1|y1|x1|, |w6|z6|y6|x6||w2|z2|y2|x2|, |w7|z7|y7|x7||w3|z3|y3|x3| ]
    // to
    // [ |x7|x6|x5|x4||x3|x2|x1|x0|, |y7|y6|y5|y4||y3|y2|y1|y0|, |z7|z6|z5|z4||z3|z2|z1|z0|, |w7|w6|w5|w4||w3|w2|w1|w0| ] 

    const __m256 y2y0x2x0 = _mm256_unpacklo_ps(v0, v2);
    const __m256 w2w0z2z0 = _mm256_unpackhi_ps(v0, v2);
    const __m256 y3y1x3x1 = _mm256_unpacklo_ps(v1, v3);
    const __m256 w3w1z3z1 = _mm256_unpackhi_ps(v1, v3);
    const __m256 xxxx = _mm256_unpacklo_ps(y2y0x2x0, y3y1x3x1);
    const __m256 yyyy = _mm256_unpackhi_ps(y2y0x2x0, y3y1x3x1);
    const __m256 zzzz = _mm256_unpacklo_ps(w2w0z2z0, w3w1z3z1);
    const __m256 wwww = _mm256_unpackhi_ps(w2w0z2z0, w3w1z3z1);

    //const v8sf y2y0x2x0 = __builtin_ia32_unpcklps256(v0, v2);
    //const v8sf w2w0z2z0 = __builtin_ia32_unpckhps256(v0, v2);
    //const v8sf y3y1x3x1 = __builtin_ia32_unpcklps256(v1, v3);
    //const v8sf w3w1z3z1 = __builtin_ia32_unpckhps256(v1, v3);
    //const v8sf xxxx = __builtin_ia32_unpcklps256(y2y0x2x0, y3y1x3x1);
    //const v8sf yyyy = __builtin_ia32_unpckhps256(y2y0x2x0, y3y1x3x1);
    //const v8sf zzzz = __builtin_ia32_unpcklps256(w2w0z2z0, w3w1z3z1);
    //const v8sf wwww = __builtin_ia32_unpckhps256(w2w0z2z0, w3w1z3z1);

    v0 = xxxx;
    v1 = yyyy;
    v2 = zzzz;
    v3 = wwww;
}

static inline __m256 _mm256_set_m128_forgcc(__m128 hi, __m128 lo){
    __m256 ret = {};
    //__m256 ret = _mm256_undefined_ps();
    //__m256 ret = ret;
    ret = _mm256_insertf128_ps(ret, lo, 0);
    ret = _mm256_insertf128_ps(ret, hi, 1);
    return ret;
}
#endif

#include"stack.hpp"

namespace ParticleSimulator{
    inline void CalcNumberAndShiftOfImageDomain
    (ReallocatableArray<F64vec> & shift_image_domain,
     const F64vec & size_root_domain,
     const F64ort & pos_my_domain,
     const F64ort & pos_target_domain,
     const bool periodic_axis[],
     const bool clear = true){
        if(clear) shift_image_domain.clearSize();
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        if( periodic_axis[0] == false && periodic_axis[1] == false){
            shift_image_domain.push_back(F64vec(0.0)); // NOTE: sign is plus
            return;
        }
#else
        if( periodic_axis[0] == false
            && periodic_axis[1] == false
            && periodic_axis[2] == false ){
            shift_image_domain.push_back(F64vec(0.0)); // NOTE: sign is plus
            return;
        }
#endif
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
                    if( (std::abs(ix) !=lev_x || periodic_axis[0] == false)
                        && (std::abs(iy) != lev_y  || periodic_axis[1] == false)) continue;
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
                        if( (std::abs(ix) !=lev_x || periodic_axis[0] == false)
                            && (std::abs(iy) != lev_y  || periodic_axis[1] == false)
                            && (std::abs(iz) != lev_z || periodic_axis[2] == false) ) continue;
                        const F64vec shift_tmp(ix*size_root_domain.x,
                                               iy*size_root_domain.y,
                                               iz*size_root_domain.z);
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
        if( pos_root_domain.contained(pos_target) ){
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
                    if( pos_root_domain.contained(pos_target + shift_tmp) ){
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
                        if( pos_root_domain.contained(pos_target + shift_tmp) ){
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
                                  const F64ort & pos_root_domain,
                                  const F64ort & pos_root_tree,
                                  const bool periodic_axis[]){
        const S32 n_proc = Comm::getNumberOfProc();
        n_ptcl = new S32[n_proc];
        n_ptcl_disp = new S32[n_proc+1];
        Comm::allGather(&n_loc, 1, n_ptcl); // TEST
        n_ptcl_disp[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ptcl_disp[i+1] = n_ptcl_disp[i] + n_ptcl[i];
        }
        const S32 n_tot = n_ptcl_disp[n_proc];
        Tptcl * ptcl_tmp = new Tptcl[ n_tot ];
        Comm::allGatherV(ptcl_org, n_loc, ptcl_tmp, n_ptcl, n_ptcl_disp); // TEST
        ReallocatableArray<F64vec> shift;
        F64vec domain_size = pos_root_domain.getFullLength();
        CalcNumberAndShiftOfImageDomain(shift, domain_size, 
                                        pos_root_tree, pos_root_domain, periodic_axis);
        const S32 n_image = shift.size();
        /*
        std::cerr<<"n_image= "<<n_image
                 <<" domain_size= "<<domain_size
                 <<" pos_root_tree= "<<pos_root_tree
                 <<" pos_root_domain= "<<pos_root_domain
                 <<" periodic_axis= "<<periodic_axis[0]<<" "<<periodic_axis[1]<<" "<<periodic_axis[2]
                 <<std::endl;
        */
        ptcl = new Tptcl[ n_tot*n_image ];
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
                                  const F64ort & pos_root_domain,
                                  const F64ort & pos_root_cell,
                                  const bool periodic_axis[]){
        S32 * n_tmp;
        S32 * n_disp_tmp;
        AllGatherParticle(ptcl, n_tmp, n_disp_tmp, ptcl_org, n_loc, pos_root_domain, pos_root_cell, periodic_axis);
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
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+: n_cell_new)
#endif
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
                //if (n_cnt == 0) {
                //    std::cout << "n_cnt = " << n_cnt << " "
                //              << "n_ptcl_tmp = " << n_ptcl_tmp << " "
                //              << "adr_ptcl_tmp = " << adr_ptcl_tmp << " "
                //              << "thread id = " << omp_get_thread_num() << " "
                //              << "rank = " << Comm::getRank() << " "
                //              << "lev_max = " << lev_max << " "
                //              << std::endl;
                //    //for (S32 j=0; j<N_CHILDREN; j++) {
                //    //    const S32 adr_tc_tmp = tc_array[i].adr_tc_ + j;
                //    //    std::cout << "j = " << j << " "
                //    //              << "tc_array[].adr_ptcl_ = " << tc_array[adr_tc_tmp].adr_ptcl_ << " "
                //    //              << "tc_array[].level_ = " << tc_array[adr_tc_tmp].level_ << " "
                //    //              << std::endl;
                //    //}
                //    for (S32 ip=adr_ptcl_tmp; ip < adr_ptcl_tmp+n_ptcl_tmp; ip++) 
                //        std::cout << "ip = " << ip << " "
                //                  << "tp[ip].key_ = " << tp[ip].key_ << " "
                //                  << "(" << GetBinString(tp[ip].key_) << ") "
                //                  << "val[0] = " << MortonKey::getCellID(0,tp[ip].key_) << " "
                //                  << "val[1] = " << MortonKey::getCellID(1,tp[ip].key_) << " "
                //                  << "val[2] = " << MortonKey::getCellID(2,tp[ip].key_) << " "
                //                  << "val[3] = " << MortonKey::getCellID(3,tp[ip].key_) << " "
                //                  << "val    = " << MortonKey::getCellID(lev_max+1,tp[ip].key_) << " "
                //                  << "tp[ip].adr_ptcl_ = " << tp[ip].adr_ptcl_ << " "
                //                  << std::endl;
                //     struct timespec ts;
                //     ts.tv_sec  = 0;
                //     ts.tv_nsec = 10000;
                //     clock_nanosleep(CLOCK_MONOTONIC, 0, &ts, NULL);
                //}
                for(S32 j=0; j<n_cnt-1; j++){
                    tc_array[adr[j]].n_ptcl_ = tc_array[adr[j+1]].adr_ptcl_ - tc_array[adr[j]].adr_ptcl_;
                    n_ptcl_cum += tc_array[adr[j]].n_ptcl_;
                }
                assert(n_cnt > 0);
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
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
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
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 j=head; j<next; j++){
                Ttc * tc_tmp = tc + j;
                const int n_tmp = tc_tmp->n_ptcl_;
                tc_tmp->mom_.init();
                if( n_tmp == 0) continue;
                else if(tc_tmp->isLeaf(n_leaf_limit)){
                    const S32 adr = tc_tmp->adr_ptcl_;
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        if( GetMSB(tp[k].adr_ptcl_) == 0 ){
                            const U32 adr_ep = tp[k].adr_ptcl_;
                            tc_tmp->mom_.accumulateAtLeaf(epj[adr_ep]);
                        }
                        else{
                            const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                            tc_tmp->mom_.accumulate(spj[adr_sp].convertToMoment());
                        }
                    }
                    tc_tmp->mom_.set();
                    for(S32 k=adr; k<adr+n_tmp; k++){
                        if( GetMSB(tp[k].adr_ptcl_) == 0 ){
                            const U32 adr_ep = tp[k].adr_ptcl_;
                            tc_tmp->mom_.accumulateAtLeaf2(epj[adr_ep]);
                        }
                        else{
                            const U32 adr_sp = ClearMSB(tp[k].adr_ptcl_);
                            tc_tmp->mom_.accumulate2(spj[adr_sp].convertToMoment());
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
            //ipg_first.back().vertex_ = GetMinBoxSingleThread(epi_first.data()+(tc_tmp->adr_ptcl_), n_tmp);
            return;
        }
        else{
            S32 adr_tc_tmp = tc_tmp->adr_tc_;
            for(S32 i=0; i<N_CHILDREN; i++){
                MakeIPGroupLong(ipg_first, tc_first, epi_first, adr_tc_tmp+i, n_grp_limit);                
                //MakeIPGroupLong<Tipg, Ttc, Tepi, Ti0, Ti1, Ti2>
                //    (ipg_first, tc_first, epi_first, adr_tc_tmp+i, n_grp_limit);
            }
        }
    }


    // NOTE: This FUnction mimic sirial code 
    template<class Tipg, class Ttcloc, class Ttcglb, class Tepj>
    inline void MakeIPGroupLongGLBTreeCellAsIPGBox(ReallocatableArray<Tipg> & ipg_first,
                                                   const ReallocatableArray<Ttcloc> & tc_loc_first,
                                                   const ReallocatableArray<Ttcglb> & tc_glb_first,
                                                   const ReallocatableArray<Tepj> & epj_first,
                                                   const S32 adr_tc_loc,
                                                   const S32 adr_tc_glb,
                                                   const S32 n_grp_limit,
                                                   const S32 n_leaf_limit){
        const Ttcloc * tc_loc_tmp = tc_loc_first.getPointer(adr_tc_loc);
        const Ttcglb * tc_glb_tmp = tc_glb_first.getPointer(adr_tc_glb);
        const S32 n_loc_tmp = tc_loc_tmp->n_ptcl_;
        const S32 n_glb_tmp = tc_glb_tmp->n_ptcl_;
        if(n_loc_tmp == 0 || n_glb_tmp == 0) return;
        else if( tc_glb_tmp->isLeaf(n_grp_limit) || tc_loc_tmp->isLeaf(n_leaf_limit)){
            ipg_first.increaseSize();
            ipg_first.back().copyFromTC(*tc_loc_tmp);
            const S32 adr_tmp = tc_glb_tmp->adr_ptcl_;
            ipg_first.back().vertex_ = GetMinBoxSingleThread(epj_first.getPointer(adr_tmp), n_glb_tmp);
            return;
        }
        else{
            S32 adr_tc_loc_tmp = tc_loc_tmp->adr_tc_;
            S32 adr_tc_glb_tmp = tc_glb_tmp->adr_tc_;
            for(S32 i=0; i<N_CHILDREN; i++){
                MakeIPGroupLongGLBTreeCellAsIPGBox(ipg_first,
                                                   tc_loc_first, tc_glb_first,
                                                   epj_first, 
                                                   adr_tc_loc_tmp+i, adr_tc_glb_tmp+i,
                                                   n_grp_limit, n_leaf_limit);
            }
        }
    }
    //#endif

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
        if(n_loc_tmp == 0 || n_glb_tmp == 0) return;
        else if( tc_glb_tmp->isLeaf(n_grp_limit) || tc_loc_tmp->isLeaf(n_leaf_limit)){
            ipg_first.increaseSize();
            ipg_first.back().copyFromTC(*tc_loc_tmp);
            //const S32 adr_tmp = tc_loc_tmp->adr_ptcl_;
            //ipg_first.back().vertex_ = GetMinBoxSingleThread(epi_first.getPointer(adr_tmp), n_loc_tmp);
            return;
        }
        else{
            S32 adr_tc_loc_tmp = tc_loc_tmp->adr_tc_;
            S32 adr_tc_glb_tmp = tc_glb_tmp->adr_tc_;
            for(S32 i=0; i<N_CHILDREN; i++){
                MakeIPGroupUseGLBTreeLong(ipg_first,
                                          tc_loc_first, tc_glb_first,
                                          epi_first, 
                                          adr_tc_loc_tmp+i, adr_tc_glb_tmp+i,
                                          n_grp_limit, n_leaf_limit);
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
                    //#ifdef DEBUG_1028
#ifdef PARTICLE_SIMULATOR_EXCHANGE_LET_ALL
                    SearchSendParticleLong<Ttc, Tep>
                        (tc_first,   tc_first[adr_tc_child].adr_tc_, ep_first,
                         id_ep_send, id_sp_send, pos_target_box,
                         r_crit_sq, n_leaf_limit);
#else
                    SearchSendParticleLong<Ttc, Tep>
                        (tc_first,   tc_first[adr_tc_child].adr_tc_, ep_first,
                         id_ep_send, id_sp_send, pos_target_box,
                         r_crit_sq*0.25, n_leaf_limit);
#endif
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
                             || (pos_target_box.contained( tc_first[adr_tc+i].mom_.getVertexOut())) ) 
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

    template<class Ttc, class Tep>
    inline void SearchSendParticleLongSymmetry(const ReallocatableArray<Ttc> & tc_first,
                                               const S32 adr_tc,
                                               const ReallocatableArray<Tep> & ep_first,
                                               ReallocatableArray<S32> & id_ep_send,
                                               ReallocatableArray<S32> & id_sp_send,
                                               const F64ort & pos_target_box_in,
                                               const F64ort & pos_target_box_out,
                                               const F64 r_crit_sq,
                                               const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            /*
            open_bits |= ( ( (pos_target_box_in.getDistanceMinSQ(pos_tmp) <= r_crit_sq) 
                             || (pos_target_box_in.contained( tc_first[adr_tc+i].mom_.getVertexOut())) 
                             || (pos_target_box_out.contained( tc_first[adr_tc+i].mom_.getVertexIn())) )
                           << i);
            */
            open_bits |= ( ( (pos_target_box_in.getDistanceMinSQ(pos_tmp) <= r_crit_sq) 
                             || (pos_target_box_in.overlapped( tc_first[adr_tc+i].mom_.getVertexOut())) 
                             || (pos_target_box_out.overlapped( tc_first[adr_tc+i].mom_.getVertexIn())) )
                           << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    SearchSendParticleLongSymmetry<Ttc, Tep>
                        (tc_first,   tc_first[adr_tc_child].adr_tc_, ep_first,
                         id_ep_send, id_sp_send, pos_target_box_in, pos_target_box_out,
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
                #if 0
                id_sp_send.push_back(adr_tc_child);
                if(tc_child->mom_.mass < 1e-10){
                    std::cerr<<"check bbb"<<std::endl;
                }
                #else
                // original
                id_sp_send.push_back(adr_tc_child);
                #endif
            }
        }
    }


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
            open_bits |= (pos_target_domain.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) << (i + N_CHILDREN) );
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
                sp_send.increaseSize();
                sp_send.back().copyFromMoment(tc_child->mom_);
                const F64vec pos_new = sp_send.back().getPos() + shift;
                sp_send.back().setPos(pos_new);
            }
        }
    }



    template<class Ttc, class Tep, class Tsp>
    inline void SearchSendParticleLongCutoffWithRootCellCheck
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
     const F64ort pos_root_cell,
     const F64vec & shift = 0.){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( (pos_target_domain.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i); // using opening criterion
#if 1
            // vertex_out is not used. 
            const F64ort child_cell_box = makeChildCellBox(i, cell_box);
            open_bits |= ( (pos_target_domain.getDistanceMinSQ(child_cell_box) <= r_cut_sq) << (i + N_CHILDREN) ); // using cutoff criterion
#else
            // cell_box is not neede
            open_bits |= (pos_target_domain.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) << (i + N_CHILDREN) );
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
                    SearchSendParticleLongCutoffWithRootCellCheck<Ttc, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_,
                         ep_first, ep_send, sp_send, child_cell_box,
                         pos_target_domain, r_crit_sq*0.25, r_cut_sq, n_leaf_limit,
                         pos_root_cell, shift);
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
                sp_send.increaseSize();
                sp_send.back().copyFromMoment(tc_child->mom_);
                const F64vec pos_new = sp_send.back().getPos() + shift;
                sp_send.back().setPos(pos_new);
            }
        }
    } 
   
    template<class Ttc, class Tep>
    inline void SearchSendParticleLongCutoffForZeroTheta
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<Tep> & ep_send,
     const F64ort & cell_box,
     const F64ort & pos_target_domain,
     const F64 r_cut_sq,
     const S32 n_leaf_limit,
     const F64ort pos_root_cell,
     const F64vec & shift = 0.){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64ort child_cell_box = makeChildCellBox(i, cell_box);
            open_bits |= ( (pos_target_domain.getDistanceMinSQ(child_cell_box) <= 1.01*r_cut_sq) << i );
            // open_bits is calculated using the cutoff radius only.
            // The factor 1.01 is introduce for the safety.
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if (n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    const F64ort child_cell_box = makeChildCellBox(i, cell_box);
                    SearchSendParticleLongCutoffForZeroTheta<Ttc, Tep>
                        (tc_first, tc_first[adr_tc_child].adr_tc_,
                         ep_first, ep_send, child_cell_box,
                         pos_target_domain, r_cut_sq, n_leaf_limit,
                         pos_root_cell, shift);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    ep_send.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        const S32 adr = adr_ptcl_tmp + ip;
                        const F64vec pos = ep_first[adr].getPos();
                        if (pos_target_domain.notContained(pos)) continue;
                        ep_send.pushBackNoCheck(ep_first[adr]);
                        const F64vec pos_new = ep_send.back().getPos() + shift;
                        ep_send.back().setPos(pos_new);
                    }
                }
            }
        }
    }    
    
    /////////////////////////////
    /// MAKE INTERACTION LIST ///
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
#elif FAST_WALK_AVX2
        __asm__("#tree walk calc distance");
        // single precision
        const __m256 r_crit_sq_v = _mm256_set1_ps((float)r_crit_sq);
        const __m256d half = _mm256_set1_pd(0.5);
        const __m256 mask_abs = (__m256)_mm256_set1_epi32(0x7fffffff);
        const __m256 zero = _mm256_set1_ps(0.0f);

        const __m256d high = _mm256_loadu_pd(&pos_target_box.high_.x);
        const __m256d low  = _mm256_loadu_pd(&pos_target_box.low_.x);
        const __m256d cen_d = _mm256_mul_pd(_mm256_add_pd(high, low), half);
        const __m256d len_d = _mm256_mul_pd(_mm256_sub_pd(high, low), half);

        const __m128 cen = _mm256_cvtpd_ps(cen_d);
        const __m128 len = _mm256_cvtpd_ps(len_d);

        const __m256 cen2 = _mm256_set_m128_forgcc(cen, cen);
        const __m256 len2 = _mm256_set_m128_forgcc(len, len);

        const __m256 cx = _mm256_shuffle_ps(cen2, cen2, 0x00);
        const __m256 cy = _mm256_shuffle_ps(cen2, cen2, 0x55);
        const __m256 cz = _mm256_shuffle_ps(cen2, cen2, 0xaa);

        const __m256 lx = _mm256_shuffle_ps(len2, len2, 0x00);
        const __m256 ly = _mm256_shuffle_ps(len2, len2, 0x55);
        const __m256 lz = _mm256_shuffle_ps(len2, len2, 0xaa);

        const F32vec & mom0 = tc_first[adr_tc+0].mom_.getPos();
        const F32vec & mom1 = tc_first[adr_tc+1].mom_.getPos();
        const F32vec & mom2 = tc_first[adr_tc+2].mom_.getPos();
        const F32vec & mom3 = tc_first[adr_tc+3].mom_.getPos();
        const F32vec & mom4 = tc_first[adr_tc+4].mom_.getPos();
        const F32vec & mom5 = tc_first[adr_tc+5].mom_.getPos();
        const F32vec & mom6 = tc_first[adr_tc+6].mom_.getPos();
        const F32vec & mom7 = tc_first[adr_tc+7].mom_.getPos();

        __m256 jx = _mm256_set_m128_forgcc(_mm_loadu_ps(&mom4.x), 
                                           _mm_loadu_ps(&mom0.x));
        __m256 jy = _mm256_set_m128_forgcc(_mm_loadu_ps(&mom5.x), 
                                           _mm_loadu_ps(&mom1.x));
        __m256 jz = _mm256_set_m128_forgcc(_mm_loadu_ps(&mom6.x), 
                                           _mm_loadu_ps(&mom2.x));
        __m256 jw = _mm256_set_m128_forgcc(_mm_loadu_ps(&mom7.x), 
                                           _mm_loadu_ps(&mom3.x));

        transpose_4ymm(jx, jy, jz, jw);

        __m256 dx = _mm256_and_ps(mask_abs, (cx - jx));
        dx = dx - lx;
        dx = _mm256_max_ps(dx, zero);
        //__m256 dy = _mm256_andnot_ps(minus_zero, (cy - jy));
        __m256 dy = _mm256_and_ps(mask_abs, (cy - jy));
        dy = dy - ly;
        dy = _mm256_max_ps(dy, zero);
        //__m256 dz = _mm256_andnot_ps(minus_zero, (cz - jz));
        __m256 dz = _mm256_and_ps(mask_abs, (cz - jz));
        dz = dz - lz;
        dz = _mm256_max_ps(dz, zero);
        __m256 rsq = _mm256_fmadd_ps(dz, dz, _mm256_fmadd_ps(dy, dy, (dx*dx)));

        open_bits = _mm256_movemask_ps(_mm256_cmp_ps(rsq, r_crit_sq_v, 0x02));
#else //FAST_WALK_K
        /*
        const F64 x_h = pos_target_box.high_.x;
        const F64 x_l = pos_target_box.low_.x;
        const F64 y_h = pos_target_box.high_.y;
        const F64 y_l = pos_target_box.low_.y;
        const F64 z_h = pos_target_box.high_.z;
        const F64 z_l = pos_target_box.low_.z;
        const F64 ix = (x_h + x_l) * 0.5;
        const F64 iy = (y_h + y_l) * 0.5;
        const F64 iz = (z_h + z_l) * 0.5;
        const F64 lx = (x_h - x_l) * 0.5;
        const F64 ly = (y_h - y_l) * 0.5;
        const F64 lz = (z_h - z_l) * 0.5;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            F64 dx = fabs(ix - pos_tmp.x);
            dx = dx - lx;
            dx = (dx > 0.0) ? dx : 0.0;
            F64 dy = fabs(iy - pos_tmp.y);
            dy = dy - ly;
            dy = (dy > 0.0) ? dy : 0.0;
            F64 dz = fabs(iz - pos_tmp.z);
            dz = dz - lz;
            dz = (dz > 0.0) ? dz : 0.0;
            F64 rsq = dx*dx + dy*dy + dz*dz;
            open_bits |= (( rsq <= r_crit_sq) << i);
        }
        */
        // original hozon
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i);
        }
#endif //FAST_WALK_K
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            __asm__("#tree walk  before branch");
            if(n_child == 0){
                __asm__("#tree walk  n_child == 0");
                continue;
            }
            else if( (open_bits>>i) & 0x1 ){
                __asm__("#tree walk  cell open");
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    __asm__("#tree walk cell is not leaf");
#ifdef PARTICLE_SIMULATOR_INTERACTION_LIST_ALL
                    MakeInteractionListLongEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         ep_list,   sp_first, sp_list, 
                         pos_target_box, r_crit_sq, n_leaf_limit);
#else
                    MakeInteractionListLongEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         ep_list,   sp_first, sp_list, 
                         pos_target_box, r_crit_sq*0.25, n_leaf_limit);
#endif
                }
                else{
                    __asm__("#tree walk cell is leaf");
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    sp_list.reserveEmptyAreaAtLeast( n_child );
#if 0
                    S32 iloc = ep_list.size();
                    // debug
                    for(S32 ip=0; ip<n_child; ip++){
                        ep_list[iloc].mass = ep_first[adr_ptcl_tmp].mass;
                        ep_list[iloc++].pos = ep_first[adr_ptcl_tmp++].pos;
                    }
                    ep_list.resizeNoInitialize(iloc);
#else
                    for(S32 ip=0; ip<n_child; ip++){
                        if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_) == 0){
                            Tep ep_tmp = ep_first[adr_ptcl_tmp++];
                            ep_list.pushBackNoCheck(ep_tmp);
                        }
                        else{
                            Tsp sp_tmp = sp_first[adr_ptcl_tmp++];
                            sp_list.pushBackNoCheck(sp_tmp);
                        }
                    }
                    /*
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                        U32 msb = GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_);
                        Tep ep_tmp = ep_first[adr_ptcl_tmp];
                        Tsp sp_tmp = sp_first[adr_ptcl_tmp];
                        ep_list.pushBackNoCheck(ep_tmp);
                        ep_list.resizeNoInitialize(ep_list.size()-msb);
                        sp_list.pushBackNoCheck(sp_tmp);
                        sp_list.resizeNoInitialize(sp_list.size()-(msb^0x1));
                    }
                    */
#endif
                }
            }
            else{
                __asm__("#tree walk substitute SP");
#if 0
                // debug
                ep_list.increaseSize();
                ep_list.back().mass = tc_child->mom_.mass;
                ep_list.back().pos = tc_child->mom_.pos;
#else
                sp_list.increaseSize();
                sp_list.back().copyFromMoment(tc_child->mom_);
#endif
            }
        }
    }

    /////////////////////////////
    // new for multiwalk (send index)
#if 0
    //underconstruction
    template<class Ttc, class Ttp, class Tep, class Tsp>
    void MakeInteractionListLongEPSPIndex(const ReallocatableArray<Ttc> & tc_first,
                                          const S32 adr_tc,
                                          const ReallocatableArray<Ttp> & tp_first,
                                          const ReallocatableArray<Tep> & ep_first,
                                          ReallocatableArray<S32> & id_ep_list,
                                          const ReallocatableArray<Tsp> & sp_first,
                                          ReallocatableArray<S32> & id_sp_list,
                                          const F64ort & pos_target_box,
                                          const F64 r_crit_sq,
                                          const S32 n_leaf_limit){
        if(tc_first[adr_tc].n_ptcl_ == 0) continue;


        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            open_bits |= ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            __asm__("#tree walk  before branch");
            if(n_child == 0){
                __asm__("#tree walk  n_child == 0");
                continue;
            }
            else if( (open_bits>>i) & 0x1 ){
                __asm__("#tree walk  cell open");
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    __asm__("#tree walk cell is not leaf");
                    MakeInteractionListLongEPSPIndex<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         id_ep_list,   sp_first, id_sp_list, 
                         pos_target_box, r_crit_sq*0.25, n_leaf_limit);
                }
                else{
                    __asm__("#tree walk cell is leaf");
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    id_ep_list.reserveEmptyAreaAtLeast( n_child );
                    id_sp_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++){
                        if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_) == 0){
                            id_ep_list.pushBackNoCheck(adr_ptcl_tmp++);
                        }
                        else{
                            id_sp_list.pushBackNoCheck(adr_ptcl_tmp++);
                        }
                    }
                }
            }
            else{
                __asm__("#tree walk substitute SP");
                S32 offset = ep_first.size();
                id_sp_list.push_back(adr_tc_child+offset);
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
                  || ( pos_target_box.contained(box_tmp) ) )
                << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    //#ifdef DEBUG_1023_2
#ifdef PARTICLE_SIMULATOR_INTERACTION_LIST_ALL
                    MakeInteractionListLongScatterEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         ep_list,   sp_first, sp_list, 
                         pos_target_box, r_crit_sq, n_leaf_limit);
#else
                    MakeInteractionListLongScatterEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         ep_list,   sp_first, sp_list, 
                         pos_target_box, r_crit_sq*0.25, n_leaf_limit);
#endif
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
    void MakeInteractionListLongSymmetryEPSP(const ReallocatableArray<Ttc> & tc_first,
                                             const S32 adr_tc,
                                             const ReallocatableArray<Ttp> & tp_first,
                                             const ReallocatableArray<Tep> & ep_first,
                                             ReallocatableArray<Tep> & ep_list,
                                             const ReallocatableArray<Tsp> & sp_first,
                                             ReallocatableArray<Tsp> & sp_list,
                                             const F64ort & pos_target_box_in,
                                             const F64ort & pos_target_box_out,
                                             const F64 r_crit_sq,
                                             const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos_tmp = tc_first[adr_tc+i].mom_.getPos();
            const F64ort box_tmp_out = tc_first[adr_tc+i].mom_.getVertexOut();
            const F64ort box_tmp_in = tc_first[adr_tc+i].mom_.getVertexIn();
            open_bits |= ( 
                ( (pos_target_box_in.getDistanceMinSQ(pos_tmp) <= r_crit_sq) 
                  || ( pos_target_box_out.contained(box_tmp_in))
                  || ( pos_target_box_in.contained(box_tmp_out)) )
                << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeInteractionListLongSymmetryEPSP<Ttc, Ttp, Tep, Tsp>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first, ep_first, 
                         ep_list,   sp_first, sp_list, 
                         pos_target_box_in, pos_target_box_out, r_crit_sq*0.25, n_leaf_limit);
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
                            //if(sp_first[adr_ptcl_tmp].mass == 0.0) { std::cerr<<"check a"<<std::endl; } // another node sent sp with mass of 0
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
#ifdef DEBUG_2017_10_30
            open_bits |= (pos_target_box.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) << (i + N_CHILDREN) );
#else
    #if 0
            // vertex_out is not used
            const F64ort child_cell_box = makeChildCellBox(i, cell_box);
            open_bits |= ( (pos_target_box.getDistanceMinSQ(child_cell_box) <= r_cut_sq) << (i + N_CHILDREN));
    #else
            // cell_box is not used. maybe faster than above one.
            // In this case, SPJ must have r_cutoff
            // becuase vertex_out of global tree depends on r_cutoff of SPJ.
            // It is not easy to tell r_cutoff to SPJ defined by users.
            open_bits |= (pos_target_box.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) << (i + N_CHILDREN) );
    #endif
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
                sp_list.increaseSize();
                sp_list.back().copyFromMoment(tc_child->mom_);
            }
        }
    }

    template<class Ttc, class Ttp, class Tep>
    void MakeInteractionListLongCutoffEPForZeroTheta(const ReallocatableArray<Ttc> & tc_first,
                                                     const S32 adr_tc,
                                                     const ReallocatableArray<Ttp> & tp_first,
                                                     const ReallocatableArray<Tep> & ep_first,
                                                     ReallocatableArray<Tep> & ep_list,
                                                     const F64ort & cell_box,
                                                     const F64ort & pos_target_box,
                                                     const F64 r_cut_sq,
                                                     const S32 n_leaf_limit){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
#if 0
            const F64ort child_cell_box = makeChildCellBox(i, cell_box);
            open_bits |= ( (pos_target_box.getDistanceMinSQ(child_cell_box) <= r_cut_sq) << i);
#else
            open_bits |= (pos_target_box.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) << i);
#endif
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            const S32 adr_tc_child = adr_tc + i;
            Ttc * tc_child = tc_first.getPointer(adr_tc_child);
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            else if( (open_bits>>i) & 0x1 ){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    const F64ort child_cell_box = makeChildCellBox(i, cell_box);
                    MakeInteractionListLongCutoffEPForZeroTheta<Ttc, Ttp, Tep>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, tp_first,
                         ep_first, ep_list, child_cell_box,
                         pos_target_box, r_cut_sq, n_leaf_limit);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    for (S32 ip=0; ip<n_child; ip++){
                        if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_) == 0){
                            ep_list.pushBackNoCheck(ep_first[adr_ptcl_tmp++]);
                        }
                    }
                }
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
        for(S32 i=0; i<N_CHILDREN; i++){
            //open_bits |= (pos_target_box.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) << i);
            open_bits |= (pos_target_box.overlapped( tc_first[adr_tc+i].mom_.getVertexOut() ) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                const S32 adr_tc_child = adr_tc + i;
                const Ttc * tc_child = tc_first + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeListUsingOuterBoundary<Ttc, Tep2, Tep3>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, ep_list,
                         pos_target_box, n_leaf_limit, shift);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
#ifdef ORIGINAL_SCATTER_MODE
                        const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                        const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq_tmp = pos_target_box.getDistanceMinSQ(pos_tmp);
                        if(dis_sq_tmp > size_tmp*size_tmp) continue;
#endif
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
    inline void MakeListUsingOuterBoundaryWithRootCellCheck
    (const Ttc * tc_first,
     const S32 adr_tc,
     const Tep2 * ep_first,
     ReallocatableArray<Tep3> & ep_list,
     const F64ort & pos_target_box, // position of domain
     const S32 n_leaf_limit,
     const F64ort pos_root_cell,
     const F64vec & shift = F64vec(0.0) ){
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (pos_target_box.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if( (open_bits>>i) & 0x1){
                const S32 adr_tc_child = adr_tc + i;
                const Ttc * tc_child = tc_first + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    MakeListUsingOuterBoundaryWithRootCellCheck<Ttc, Tep2, Tep3>
                        (tc_first, tc_first[adr_tc_child].adr_tc_, ep_first, ep_list,
                         pos_target_box, n_leaf_limit, pos_root_cell, shift);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                        const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
#ifdef ORIGINAL_SCATTER_MODE
                        const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq_tmp = pos_target_box.getDistanceMinSQ(pos_tmp);
                        if(dis_sq_tmp > size_tmp*size_tmp) continue;
#endif
                        if( pos_root_cell.notContained(pos_tmp+shift) ) continue;  // add M.I. 2016/03/12
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
    inline void MakeListUsingOuterBoundaryAndInnerBoundary(const Ttc * tc_first,
                                                           const S32 adr_tc,
                                                           const Tep2 * ep_first,
                                                           ReallocatableArray<Tep3> & ep_list,
                                                           const F64ort & pos_target_box_out, // position of domain
                                                           const F64ort & pos_target_box_in, // position of domain
                                                           const S32 n_leaf_limit,
                                                           const F64vec & shift = F64vec(0.0) ){
        //asm("# MakeListUsingOuterBoundaryAndInnerBoundary");
#ifdef FAST_WALK_K
        static const _fjsp_v2r8 ONE = _fjsp_set_v2r8(1.0, 1.0);
        F64 open_mask[N_CHILDREN];
        const _fjsp_v2r8 xi_h_in = _fjsp_set_v2r8(pos_target_box_in.high_.x, pos_target_box_in.high_.x);
        const _fjsp_v2r8 xi_l_in = _fjsp_set_v2r8(pos_target_box_in.low_.x, pos_target_box_in.low_.x);
        const _fjsp_v2r8 yi_h_in = _fjsp_set_v2r8(pos_target_box_in.high_.y, pos_target_box_in.high_.y);
        const _fjsp_v2r8 yi_l_in = _fjsp_set_v2r8(pos_target_box_in.low_.y, pos_target_box_in.low_.y);
        const _fjsp_v2r8 zi_h_in = _fjsp_set_v2r8(pos_target_box_in.high_.z, pos_target_box_in.high_.z);
        const _fjsp_v2r8 zi_l_in = _fjsp_set_v2r8(pos_target_box_in.low_.z, pos_target_box_in.low_.z);

        const _fjsp_v2r8 xi_h_out = _fjsp_set_v2r8(pos_target_box_out.high_.x, pos_target_box_out.high_.x);
        const _fjsp_v2r8 xi_l_out = _fjsp_set_v2r8(pos_target_box_out.low_.x, pos_target_box_out.low_.x);
        const _fjsp_v2r8 yi_h_out = _fjsp_set_v2r8(pos_target_box_out.high_.y, pos_target_box_out.high_.y);
        const _fjsp_v2r8 yi_l_out = _fjsp_set_v2r8(pos_target_box_out.low_.y, pos_target_box_out.low_.y);
        const _fjsp_v2r8 zi_h_out = _fjsp_set_v2r8(pos_target_box_out.high_.z, pos_target_box_out.high_.z);
        const _fjsp_v2r8 zi_l_out = _fjsp_set_v2r8(pos_target_box_out.low_.z, pos_target_box_out.low_.z);
        for(S32 i=0; i<N_CHILDREN/4; i++){
            const S32 i0 = i*4;
            const S32 i1 = i*4+1;
            const S32 i2 = i*4+2;
            const S32 i3 = i*4+3;
            const F64ort mom0_in = tc_first[adr_tc+i0].mom_.getVertexIn();
            const F64ort mom0_out = tc_first[adr_tc+i0].mom_.getVertexOut();
            const F64ort mom1_in = tc_first[adr_tc+i1].mom_.getVertexIn();
            const F64ort mom1_out = tc_first[adr_tc+i1].mom_.getVertexOut();
            const F64ort mom2_in = tc_first[adr_tc+i2].mom_.getVertexIn();
            const F64ort mom2_out = tc_first[adr_tc+i2].mom_.getVertexOut();
            const F64ort mom3_in = tc_first[adr_tc+i3].mom_.getVertexIn();
            const F64ort mom3_out = tc_first[adr_tc+i3].mom_.getVertexOut();

            const _fjsp_v2r8 xj_h_in_a = _fjsp_set_v2r8(mom0_in.high_.x,   mom1_in.high_.x);
            const _fjsp_v2r8 xj_l_in_a = _fjsp_set_v2r8(mom0_in.low_.x,    mom1_in.low_.x);
            const _fjsp_v2r8 xj_h_in_b = _fjsp_set_v2r8(mom2_in.high_.x,   mom3_in.high_.x);
            const _fjsp_v2r8 xj_l_in_b = _fjsp_set_v2r8(mom2_in.low_.x,    mom3_in.low_.x);
            const _fjsp_v2r8 xj_h_out_a = _fjsp_set_v2r8(mom0_out.high_.x, mom1_out.high_.x);
            const _fjsp_v2r8 xj_l_out_a = _fjsp_set_v2r8(mom0_out.low_.x,  mom1_out.low_.x);
            const _fjsp_v2r8 xj_h_out_b = _fjsp_set_v2r8(mom2_out.high_.x, mom3_out.high_.x);
            const _fjsp_v2r8 xj_l_out_b = _fjsp_set_v2r8(mom2_out.low_.x,  mom3_out.low_.x);



            const _fjsp_v2r8 yj_h_in_a = _fjsp_set_v2r8(mom0_in.high_.y,   mom1_in.high_.y);
            const _fjsp_v2r8 yj_l_in_a = _fjsp_set_v2r8(mom0_in.low_.y,    mom1_in.low_.y);
            const _fjsp_v2r8 yj_h_in_b = _fjsp_set_v2r8(mom2_in.high_.y,   mom3_in.high_.y);
            const _fjsp_v2r8 yj_l_in_b = _fjsp_set_v2r8(mom2_in.low_.y,    mom3_in.low_.y);
            const _fjsp_v2r8 yj_h_out_a = _fjsp_set_v2r8(mom0_out.high_.y, mom1_out.high_.y);
            const _fjsp_v2r8 yj_l_out_a = _fjsp_set_v2r8(mom0_out.low_.y,  mom1_out.low_.y);
            const _fjsp_v2r8 yj_h_out_b = _fjsp_set_v2r8(mom2_out.high_.y, mom3_out.high_.y);
            const _fjsp_v2r8 yj_l_out_b = _fjsp_set_v2r8(mom2_out.low_.y,  mom3_out.low_.y);

            const _fjsp_v2r8 zj_h_in_a = _fjsp_set_v2r8(mom0_in.high_.z,   mom1_in.high_.z);
            const _fjsp_v2r8 zj_l_in_a = _fjsp_set_v2r8(mom0_in.low_.z,    mom1_in.low_.z);
            const _fjsp_v2r8 zj_h_in_b = _fjsp_set_v2r8(mom2_in.high_.z,   mom3_in.high_.z);
            const _fjsp_v2r8 zj_l_in_b = _fjsp_set_v2r8(mom2_in.low_.z,    mom3_in.low_.z);
            const _fjsp_v2r8 zj_h_out_a = _fjsp_set_v2r8(mom0_out.high_.z, mom1_out.high_.z);
            const _fjsp_v2r8 zj_l_out_a = _fjsp_set_v2r8(mom0_out.low_.z,  mom1_out.low_.z);
            const _fjsp_v2r8 zj_h_out_b = _fjsp_set_v2r8(mom2_out.high_.z, mom3_out.high_.z);
            const _fjsp_v2r8 zj_l_out_b = _fjsp_set_v2r8(mom2_out.low_.z,  mom3_out.low_.z);


            const _fjsp_v2r8 mask_x_io_ji_a = _fjsp_and_v2r8( _fjsp_cmple_v2r8(xi_l_out, xj_h_in_a), _fjsp_cmple_v2r8(xj_l_in_a, xi_h_out) );  // i:out, j:in
            const _fjsp_v2r8 mask_x_io_ji_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(xi_l_out, xj_h_in_b), _fjsp_cmple_v2r8(xj_l_in_b, xi_h_out) ); 
            const _fjsp_v2r8 mask_x_ii_jo_a = _fjsp_and_v2r8( _fjsp_cmple_v2r8(xi_l_in, xj_h_out_a), _fjsp_cmple_v2r8(xj_l_out_a, xi_h_in) );  // i:in, j:out
            const _fjsp_v2r8 mask_x_ii_jo_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(xi_l_in, xj_h_out_b), _fjsp_cmple_v2r8(xj_l_out_b, xi_h_in) ); 

            const _fjsp_v2r8 mask_y_io_ji_a = _fjsp_and_v2r8( _fjsp_cmple_v2r8(yi_l_out, yj_h_in_a), _fjsp_cmple_v2r8(yj_l_in_a, yi_h_out) );  // i:out, j:in
            const _fjsp_v2r8 mask_y_io_ji_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(yi_l_out, yj_h_in_b), _fjsp_cmple_v2r8(yj_l_in_b, yi_h_out) ); 
            const _fjsp_v2r8 mask_y_ii_jo_a = _fjsp_and_v2r8( _fjsp_cmple_v2r8(yi_l_in, yj_h_out_a), _fjsp_cmple_v2r8(yj_l_out_a, yi_h_in) );  // i:in, j:out
            const _fjsp_v2r8 mask_y_ii_jo_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(yi_l_in, yj_h_out_b), _fjsp_cmple_v2r8(yj_l_out_b, yi_h_in) ); 

            const _fjsp_v2r8 mask_z_io_ji_a = _fjsp_and_v2r8( _fjsp_cmple_v2r8(zi_l_out, zj_h_in_a), _fjsp_cmple_v2r8(zj_l_in_a, zi_h_out) );  // i:out, j:in
            const _fjsp_v2r8 mask_z_io_ji_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(zi_l_out, zj_h_in_b), _fjsp_cmple_v2r8(zj_l_in_b, zi_h_out) ); 
            const _fjsp_v2r8 mask_z_ii_jo_a = _fjsp_and_v2r8( _fjsp_cmple_v2r8(zi_l_in, zj_h_out_a), _fjsp_cmple_v2r8(zj_l_out_a, zi_h_in) );  // i:in, j:out
            const _fjsp_v2r8 mask_z_ii_jo_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(zi_l_in, zj_h_out_b), _fjsp_cmple_v2r8(zj_l_out_b, zi_h_in) ); 

            const _fjsp_v2r8 mask_tmp_a = 
            _fjsp_and_v2r8( 
                _fjsp_or_v2r8( 
                    _fjsp_and_v2r8(mask_x_io_ji_a, _fjsp_and_v2r8(mask_y_io_ji_a, mask_z_io_ji_a)),
                    _fjsp_and_v2r8(mask_x_ii_jo_a, _fjsp_and_v2r8(mask_y_ii_jo_a, mask_z_ii_jo_a))
                    )
                , ONE);

            const _fjsp_v2r8 mask_tmp_b = 
            _fjsp_and_v2r8( 
                _fjsp_or_v2r8( 
                    _fjsp_and_v2r8(mask_x_io_ji_b, _fjsp_and_v2r8(mask_y_io_ji_b, mask_z_io_ji_b)),
                    _fjsp_and_v2r8(mask_x_ii_jo_b, _fjsp_and_v2r8(mask_y_ii_jo_b, mask_z_ii_jo_b))
                    )
                , ONE);
            _fjsp_storeh_v2r8(open_mask+i0, mask_tmp_a);
            _fjsp_storel_v2r8(open_mask+i1, mask_tmp_a);
            _fjsp_storeh_v2r8(open_mask+i2, mask_tmp_b);
            _fjsp_storel_v2r8(open_mask+i3, mask_tmp_b);
        }
#else //FAST_WALK_K
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= ( (pos_target_box_out.contained( tc_first[adr_tc+i].mom_.getVertexIn() ) 
                            || pos_target_box_in.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) )
                           << i);
        }
#endif //FAST_WALK_K
        Ttc tc_child[N_CHILDREN];
        for(S32 i=0; i<N_CHILDREN; i++){
            tc_child[i] = tc_first[adr_tc+i];
        }
        for(S32 i=0; i<N_CHILDREN; i++){
#ifdef FAST_WALK_K
            if( open_mask[i] == 1.0){
#else
            if( (open_bits>>i) & 0x1){
#endif
                const S32 n_child = tc_child[i].n_ptcl_;
                if(n_child == 0) continue;
                if( !(tc_child[i].isLeaf(n_leaf_limit)) ){
                    MakeListUsingOuterBoundaryAndInnerBoundary<Ttc, Tep2, Tep3>
                        (tc_first, tc_child[i].adr_tc_, ep_first, ep_list,
                         pos_target_box_out, pos_target_box_in, n_leaf_limit, shift);
                }
                else{
                    S32 adr_ptcl_tmp = tc_child[i].adr_ptcl_;
                    ep_list.reserveEmptyAreaAtLeast( n_child );
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
#if 1
                        ep_list.increaseSize();
                        ep_list.back() = ep_first[adr_ptcl_tmp];
                        const F64vec pos_new = ep_list.back().getPos() + shift;
                        ep_list.back().setPos(pos_new);
#else
                        const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                        const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq_tmp = pos_target_box_in.getDistanceMinSQ(pos_tmp);
                        if( pos_target_box_out.notContained(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                        ep_list.increaseSize();
                        ep_list.back() = ep_first[adr_ptcl_tmp];
                        const F64vec pos_new = ep_list.back().getPos() + shift;
                        ep_list.back().setPos(pos_new);
#endif
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

        //asm("# MakeListUsingInnerBoundary");
#ifdef FAST_WALK_K
        static const _fjsp_v2r8 ONE = _fjsp_set_v2r8(1.0, 1.0);
        F64 open_mask[N_CHILDREN];
        const _fjsp_v2r8 xi_h = _fjsp_set_v2r8(pos_target_box.high_.x, pos_target_box.high_.x);
        const _fjsp_v2r8 xi_l = _fjsp_set_v2r8(pos_target_box.low_.x, pos_target_box.low_.x);
        const _fjsp_v2r8 yi_h = _fjsp_set_v2r8(pos_target_box.high_.y, pos_target_box.high_.y);
        const _fjsp_v2r8 yi_l = _fjsp_set_v2r8(pos_target_box.low_.y, pos_target_box.low_.y);
        const _fjsp_v2r8 zi_h = _fjsp_set_v2r8(pos_target_box.high_.z, pos_target_box.high_.z);
        const _fjsp_v2r8 zi_l = _fjsp_set_v2r8(pos_target_box.low_.z, pos_target_box.low_.z);
        for(S32 i=0; i<N_CHILDREN/4; i++){
            const S32 i0 = i*4;
            const F64ort mom0 = tc_first[adr_tc+i0].mom_.getVertexIn();
            const S32 i1 = i*4+1;
            const F64ort mom1 = tc_first[adr_tc+i1].mom_.getVertexIn();
            const S32 i2 = i*4+2;
            const F64ort mom2 = tc_first[adr_tc+i2].mom_.getVertexIn();
            const S32 i3 = i*4+3;
            const F64ort mom3 = tc_first[adr_tc+i3].mom_.getVertexIn();

            const _fjsp_v2r8 xj_h = _fjsp_set_v2r8(mom0.high_.x, mom1.high_.x);
            const _fjsp_v2r8 xj_l = _fjsp_set_v2r8(mom0.low_.x, mom1.low_.x);
            const _fjsp_v2r8 xj_h_b = _fjsp_set_v2r8(mom2.high_.x, mom3.high_.x);
            const _fjsp_v2r8 xj_l_b = _fjsp_set_v2r8(mom2.low_.x, mom3.low_.x);
            const _fjsp_v2r8 mask_x = _fjsp_and_v2r8( _fjsp_cmple_v2r8(xi_l, xj_h), _fjsp_cmple_v2r8(xj_l, xi_h) ); 
            const _fjsp_v2r8 mask_x_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(xi_l, xj_h_b), _fjsp_cmple_v2r8(xj_l_b, xi_h) ); 

            const _fjsp_v2r8 yj_h = _fjsp_set_v2r8(mom0.high_.y, mom1.high_.y);
            const _fjsp_v2r8 yj_l = _fjsp_set_v2r8(mom0.low_.y, mom1.low_.y);
            const _fjsp_v2r8 yj_h_b = _fjsp_set_v2r8(mom2.high_.y, mom3.high_.y);
            const _fjsp_v2r8 yj_l_b = _fjsp_set_v2r8(mom2.low_.y, mom3.low_.y);
            const _fjsp_v2r8 mask_y = _fjsp_and_v2r8( _fjsp_cmple_v2r8(yi_l, yj_h), _fjsp_cmple_v2r8(yj_l, yi_h) );
            const _fjsp_v2r8 mask_y_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(yi_l, yj_h_b), _fjsp_cmple_v2r8(yj_l_b, yi_h) );

            const _fjsp_v2r8 zj_h = _fjsp_set_v2r8(mom0.high_.z, mom1.high_.z);
            const _fjsp_v2r8 zj_l = _fjsp_set_v2r8(mom0.low_.z, mom1.low_.z);
            const _fjsp_v2r8 zj_h_b = _fjsp_set_v2r8(mom2.high_.z, mom3.high_.z);
            const _fjsp_v2r8 zj_l_b = _fjsp_set_v2r8(mom2.low_.z, mom3.low_.z);
            const _fjsp_v2r8 mask_z = _fjsp_and_v2r8( _fjsp_cmple_v2r8(zi_l, zj_h), _fjsp_cmple_v2r8(zj_l, zi_h) );
            const _fjsp_v2r8 mask_z_b = _fjsp_and_v2r8( _fjsp_cmple_v2r8(zi_l, zj_h_b), _fjsp_cmple_v2r8(zj_l_b, zi_h) );

            // 1;open(go deeper), 0:not open
            const _fjsp_v2r8 mask_tmp_a = _fjsp_and_v2r8( _fjsp_and_v2r8(mask_x, _fjsp_and_v2r8(mask_y, mask_z)), ONE);
            const _fjsp_v2r8 mask_tmp_b = _fjsp_and_v2r8( _fjsp_and_v2r8(mask_x_b, _fjsp_and_v2r8(mask_y_b, mask_z_b)), ONE);
            _fjsp_storeh_v2r8(open_mask+i0, mask_tmp_a);
            _fjsp_storel_v2r8(open_mask+i1, mask_tmp_a);
            _fjsp_storeh_v2r8(open_mask+i2, mask_tmp_b);
            _fjsp_storel_v2r8(open_mask+i3, mask_tmp_b);

        }

#else //FAST_WALK_K
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (pos_target_box.contained( tc_first[adr_tc+i].mom_.getVertexIn() ) << i);
        }
#endif //FAST_WALK_K
        for(S32 i=0; i<N_CHILDREN; i++){
#ifdef FAST_WALK_K
            if( open_mask[i] == 1.0){
#else
            if( (open_bits>>i) & 0x1){
#endif
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

    template<class Ttc>
    inline void IsOverlappedUsingOuterBoundaryOfTreeImpl(const Ttc * tc_first,
                                                         const S32 adr_tc,
                                                         const F64ort & pos_box,
                                                         const S32 n_leaf_limit,
                                                         U32 & is_overlapped){

        if(is_overlapped) return;
        U32 open_bits = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bits |= (tc_first[adr_tc+i].mom_.getVertexOut().contained(pos_box) << i);
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
            is_overlapped = tc_first->mom_.getVertexOut().contained(pos_box);
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
            if( tc_first_B[adr_tc_B+i].mom_.getVertexOut().contained(pos_domain) 
                || IsOverlappedUsingOuterBoundaryOfTree(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexIn(), n_leaf_limit_A)){
                open_bits |= ONE << i;
            }
/*
            open_bits |= (tc_first_B[adr_tc_B+i].mom_.getVertexOut().contained(pos_domain) << i);
            if( (~(open_bits>>i)) & 0x1 ){
                open_bits |= IsOverlappedUsingOuterBoundaryOfTree(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexIn(), n_leaf_limit_A) << i;
                //open_bits |= IsOverlappedUsingOuterBoundaryOfTree(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexOut(), n_leaf_limit_A) << i;
                //open_bits |= tc_first_A->mom_.getVertexOut().contained(tc_first_B[adr_tc_B+i].mom_.getVertexIn()) << i; // out
                //open_bits |= tc_first_A->mom_.getVertexOut().contained(tc_first_B[adr_tc_B+i].mom_.getVertexOut()) << i; // OK
                //open_bits |= tc_first_B[adr_tc_B+i].mom_.getVertexIn().contained(tc_first_A->mom_.getVertexOut()) << i; // out
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
                        // NOTE: need to be consistent with MakeListUsingOuterBoundary()
                        const F64vec pos_tmp = ep_first_B[adr_ptcl_tmp].getPos();
                        const F64 len_sq = ep_first_B[adr_ptcl_tmp].getRSearch() * ep_first_B[adr_ptcl_tmp].getRSearch();
                        const F64 dis_sq = pos_domain.getDistanceMinSQ(pos_tmp);
                        if(dis_sq <= len_sq){
                            //std::cerr<<"original version"<<std::endl;
                            continue;
                        }
#else
                        if( pos_domain.contained(tc_first_B[adr_tc_child].mom_.getVertexOut()) ){
                            //std::cerr<<"new version"<<std::endl;
                            continue;
                        }
#endif
                        //if( pos_domain.contained(tc_first_B[adr_tc_child].mom_.getVertexOut()) ) continue;
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
            if( tc_first_B->mom_.getVertexOut().contained(pos_domain) || IsOverlappedUsingOuterBoundaryOfTree(tc_first_A, tc_first_B->mom_.getVertexIn(), n_leaf_limit_A)){
                S32 adr_ptcl_tmp = tc_first_B->adr_ptcl_;
                const S32 n_child = tc_first_B->n_ptcl_;
                for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                    // NOTE: need to be consistent with MakeListUsingOuterBoundary()
#ifdef ORIGINAL_SCATTER_MODE
                    const F64vec pos_tmp = ep_first_B[adr_ptcl_tmp].getPos();
                    const F64 len_sq = ep_first_B[adr_ptcl_tmp].getRSearch() * ep_first_B[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq = pos_domain.getDistanceMinSQ(pos_tmp);
                    if(dis_sq <= len_sq) continue;
#else
                    if( tc_first_B->mom_.getVertexOut().contained(pos_domain) ) return;
#endif
                    id_ptcl_send.push_back(adr_ptcl_tmp);
                }
            }
        }
    }

#if 0
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
                open_bits |= vertex_cell.contained(pos_target) << i;
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
#else
    // add 2016 2/5
    // for neighbour search
    template<class Tvec, class Ttc, class Ttp, class Tep>
    inline void SearchNeighborListOneParticleScatter(const Tvec & pos_target,
                                                     const Ttc * tc_first,
                                                     const Ttp * tp_first,
                                                     const S32 adr_tc,
                                                     const ReallocatableArray<Tep> & ep_first,
                                                     ReallocatableArray<Tep> & ep_list,
                                                     const S32 n_leaf_limit){
        if( tc_first[adr_tc].isLeaf(n_leaf_limit) ){
            const S32 n_tmp = tc_first[adr_tc].n_ptcl_;
            S32 adr_ptcl_tmp = tc_first[adr_tc].adr_ptcl_;
            for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                const F64vec dr = ep_first[adr_ptcl_tmp].getPos() - pos_target;
                const F64 r_search = ep_first[adr_ptcl_tmp].getRSearch();
#if 0 // debug
                if(Comm::getRank() == 0 && pos_target.x > 0.931 && pos_target.x < 0.932){
                    std::cerr<<"pos_target="<<pos_target<<std::endl;
                    std::cerr<<"dr="<<dr<<" dr*dr="<<dr*dr<<std::endl;
                    std::cerr<<"r_search="<<r_search<<" r_search*r_search="<<r_search*r_search<<std::endl;
                    std::cerr<<"GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_)="<<GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_)<<std::endl;
                    std::cerr<<"dr*dr <= r_search*r_search="<< ((dr*dr) <= (r_search*r_search)) <<std::endl;
                }
                if( dr*dr <= r_search*r_search){
#else // original
                if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_)==0 && dr*dr <= r_search*r_search){
#endif
#if 0 // debug
                    if(Comm::getRank() == 0 && pos_target.x > 0.931 && pos_target.x < 0.932){
                        std::cerr<<"regist ep_list"<<std::endl;
                        std::cerr<<"0: ep_list.size()="<<ep_list.size()<<std::endl;
                    }
#endif
                    ep_list.push_back( ep_first[adr_ptcl_tmp] );
#if 0 // debug
                    if(Comm::getRank() == 0 && pos_target.x > 0.931 && pos_target.x < 0.932){
                        std::cerr<<"1:ep_list.size()="<<ep_list.size()<<std::endl;
                    }
#endif

                }
            }
        }
        else{
            const S32 adr_tc_child = tc_first[adr_tc].adr_tc_;
            U32 open_bits = 0;
            for(S32 i=0; i<N_CHILDREN; i++){
                const F64ort vertex_cell = tc_first[adr_tc_child+i].mom_.getVertexOut();
                //open_bits |= vertex_cell.contained(pos_target) << i;
                open_bits |= vertex_cell.overlapped(pos_target) << i;
            }
#if 0 // debug
            /*
            if(Comm::getRank() == 0 && pos_target.x > 0.931 && pos_target.x < 0.932){
                std::cerr<<"open_bits="<<open_bits<<std::endl;
            }
            */
            open_bits = 255;
#endif
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_child = tc_first[adr_tc_child+i].n_ptcl_;
                if( n_child == 0 ) continue;
                else if( (open_bits >> i) & 0x1 ){
                    SearchNeighborListOneParticleScatter
                        (pos_target, tc_first, tp_first,    adr_tc_child+i,
                         ep_first,   ep_list,  n_leaf_limit);
                }
            }
        }
        /*
        if(Comm::getRank() == 0 && pos_target.x > 0.6701 && pos_target.x < 0.6702){
            std::cerr<<"&ep_list="<<&ep_list<<std::endl;
            std::cerr<<"2:ep_list.size()="<<ep_list.size()<<std::endl;
        }
        */
    }
#endif

    template<class Tvec, class Ttc, class Ttp, class Tep>
    inline void SearchNeighborListOneParticleGather(const Tvec & pos_target,
                                                    const F64  & r_search_sq,
                                                    const Ttc * tc_first,
                                                    const Ttp * tp_first,
                                                    const S32 adr_tc,
                                                    const ReallocatableArray<Tep> & ep_first,
                                                    ReallocatableArray<Tep> & ep_list,
                                                    const S32 n_leaf_limit){
        if( tc_first[adr_tc].isLeaf(n_leaf_limit) ){
            const S32 n_tmp = tc_first[adr_tc].n_ptcl_;
            S32 adr_ptcl_tmp = tc_first[adr_tc].adr_ptcl_;
            for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                const F64vec dr = ep_first[adr_ptcl_tmp].getPos() - pos_target;
                if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_)==0 && dr*dr <= r_search_sq){
                    ep_list.push_back( ep_first[adr_ptcl_tmp] );
                }
            }
        }
        else{
            const S32 adr_tc_child = tc_first[adr_tc].adr_tc_;
            U32 open_bits = 0;
            for(S32 i=0; i<N_CHILDREN; i++){
                const F64ort vertex_cell = tc_first[adr_tc_child+i].mom_.getVertexIn();
                open_bits |= (vertex_cell.getDistanceMinSQ(pos_target) <= r_search_sq)<< i;
            }
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_child = tc_first[adr_tc_child+i].n_ptcl_;
                if( n_child == 0 ) continue;
                else if( (open_bits >> i) & 0x1 ){
                    SearchNeighborListOneParticleGather
                        (pos_target, r_search_sq, tc_first, tp_first,    adr_tc_child+i,
                         ep_first,   ep_list,     n_leaf_limit);
                }
            }
        }
    }


    template<class Tvec, class Ttc, class Ttp, class Tep>
    inline void SearchNeighborListOneParticleSymmetry(const Tvec & pos_target,
                                                      const F64  & r_search_i_sq,
                                                      const Ttc * tc_first,
                                                      const Ttp * tp_first,
                                                      const S32 adr_tc,
                                                      const ReallocatableArray<Tep> & ep_first,
                                                      ReallocatableArray<Tep> & ep_list,
                                                      const S32 n_leaf_limit){
        if( tc_first[adr_tc].isLeaf(n_leaf_limit) ){
            const S32 n_tmp = tc_first[adr_tc].n_ptcl_;
            S32 adr_ptcl_tmp = tc_first[adr_tc].adr_ptcl_;
            for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
            const F64 r_search_j_sq = ep_first[adr_ptcl_tmp].getRSearch()*ep_first[adr_ptcl_tmp].getRSearch();
            const F64 r_search_sq =  (r_search_i_sq > r_search_j_sq) ?  r_search_i_sq : r_search_j_sq;
                const F64vec dr = ep_first[adr_ptcl_tmp].getPos() - pos_target;
                if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_)==0 && dr*dr <= r_search_sq){
                    ep_list.push_back( ep_first[adr_ptcl_tmp] );
                }
            }
        }
        else{
            const S32 adr_tc_child = tc_first[adr_tc].adr_tc_;
            U32 open_bits = 0;
            for(S32 i=0; i<N_CHILDREN; i++){
                const F64ort vertex_cell_in  = tc_first[adr_tc_child+i].mom_.getVertexIn();
                const F64ort vertex_cell_out = tc_first[adr_tc_child+i].mom_.getVertexOut();
                open_bits |= ( (vertex_cell_in.getDistanceMinSQ(pos_target) <= r_search_i_sq) || vertex_cell_out.overlapped(pos_target))<< i;
            }
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_child = tc_first[adr_tc_child+i].n_ptcl_;
                if( n_child == 0 ) continue;
                else if( (open_bits >> i) & 0x1 ){
                    SearchNeighborListOneParticleSymmetry
                        (pos_target, r_search_i_sq, tc_first, tp_first,    adr_tc_child+i,
                         ep_first,   ep_list,     n_leaf_limit);
                }
            }
        }
    }    

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
                open_bits |= vertex_cell.contained(pos_box) << i;
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

#if 0
    template<class Ttc, class Tep2, class Tep3>
    inline CountT MakeListUsingOuterBoundaryAndInnerBoundaryIteration
        (const Ttc * tc_first,
         const S32 adr_tc,
         const Tep2 * ep_first,
         ReallocatableArray<Tep3> & ep_list,
         const F64ort & pos_target_box_out, // position of domain
         const F64ort & pos_target_box_in, // position of domain
         const S32 n_leaf_limit,
         const F64vec & shift = F64vec(0.0) ){
        CountT n_cell_open = 0;

        if( tc_first[adr_tc].n_ptcl_ == 0) return n_cell_open;
        else if( tc_first[adr_tc].isLeaf(n_leaf_limit) ){
            S32 adr_ptcl_tmp = tc_first[adr_tc].adr_ptcl_;
            const S32 n_child = tc_first[adr_tc].n_ptcl_;
            for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
#if 1
                ep_list.increaseSize();
                ep_list.back() = ep_first[adr_ptcl_tmp];
                const F64vec pos_new = ep_list.back().getPos() + shift;
                ep_list.back().setPos(pos_new);
#else
                const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                const F64 dis_sq_tmp = pos_target_box_in.getDistanceMinSQ(pos_tmp);
                if( pos_target_box_out.notContained(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                ep_list.increaseSize();
                ep_list.back() = ep_first[adr_ptcl_tmp];
                const F64vec pos_new = ep_list.back().getPos() + shift;
                ep_list.back().setPos(pos_new);
#endif
            }
            return n_cell_open;
        }

        //static __thread Stack<S32, 1000> adr_tc_stack;
        Stack<S32> adr_tc_stack;
        adr_tc_stack.init();
        adr_tc_stack.push(tc_first[adr_tc].adr_tc_);  // push first child
        while( !adr_tc_stack.empty() ){
            n_cell_open++;
            S32 adr_tc_curr = adr_tc_stack.pop();
            U32 open_bits = 0;
            for(S32 i=0; i<N_CHILDREN; i++){
                open_bits |= ( (pos_target_box_out.contained( tc_first[adr_tc_curr+i].mom_.getVertexIn() ) 
                                || pos_target_box_in.contained( tc_first[adr_tc_curr+i].mom_.getVertexOut() ) )
                               << i);
            }
            for(S32 i=0; i<N_CHILDREN; i++){
                //const Ttc * tc_curr = tc_first.getPointer(adr_tc_curr+i);
                const Ttc * tc_curr = tc_first+(adr_tc_curr+i);
                const S32 n_child = tc_curr->n_ptcl_;
                if( (open_bits>>i) & 0x1){
                    if(n_child == 0) continue;
                    else if( !(tc_curr->isLeaf(n_leaf_limit)) ){
                        adr_tc_stack.push( tc_curr->adr_tc_ );
                    }
                    else{
                        S32 adr_ptcl_tmp = tc_curr->adr_ptcl_;
                        ep_list.reserveEmptyAreaAtLeast( n_child );
                        for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
#if 1
                            ep_list.increaseSize();
                            ep_list.back() = ep_first[adr_ptcl_tmp];
                            const F64vec pos_new = ep_list.back().getPos() + shift;
                            ep_list.back().setPos(pos_new);
#else
                            const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                            const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                            const F64 dis_sq_tmp = pos_target_box_in.getDistanceMinSQ(pos_tmp);
                            if( pos_target_box_out.notContained(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                            ep_list.increaseSize();
                            ep_list.back() = ep_first[adr_ptcl_tmp];
                            const F64vec pos_new = ep_list.back().getPos() + shift;
                            ep_list.back().setPos(pos_new);
#endif
                        }
                    }
                }
            }
        }
        return n_cell_open;
    }

#else
    template<class Ttc, class Tep2, class Tep3>
    inline CountT MakeListUsingOuterBoundaryAndInnerBoundaryIteration
        (const Ttc * tc_first,
         const S32 adr_tc,
         const Tep2 * ep_first,
         ReallocatableArray<Tep3> & ep_list,
         const F64ort & pos_target_box_out, // position of domain
         const F64ort & pos_target_box_in, // position of domain
         const S32 n_leaf_limit,
         const F64vec & shift = F64vec(0.0) ){
        CountT n_cell_open = 0;
        //static __thread Stack<S32, 1000> adr_tc_stack;
        Stack<S32> adr_tc_stack;
        adr_tc_stack.init();
        adr_tc_stack.push(adr_tc);  // push root cell
        while( !adr_tc_stack.empty() ){
            n_cell_open++;
            S32 adr_tc_curr = adr_tc_stack.pop();
            if( tc_first[adr_tc_curr].isLeaf(n_leaf_limit) ){
                // leaf
                const S32 n_child = tc_first[adr_tc_curr].n_ptcl_;
                S32 adr_ptcl_tmp = tc_first[adr_tc_curr].adr_ptcl_;
                for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
#if 1
                    ep_list.increaseSize();
                    ep_list.back() = ep_first[adr_ptcl_tmp];
                    const F64vec pos_new = ep_list.back().getPos() + shift;
                    ep_list.back().setPos(pos_new);
#else
                    const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                    const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq_tmp = pos_target_box_in.getDistanceMinSQ(pos_tmp);
                    if( pos_target_box_out.notContained(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                    ep_list.increaseSize();
                    ep_list.back() = ep_first[adr_ptcl_tmp];
                    const F64vec pos_new = ep_list.back().getPos() + shift;
                    ep_list.back().setPos(pos_new);
#endif
                }
            }
            else{
                // not leaf
                S32 adr_tc_child = tc_first[adr_tc_curr].adr_tc_;
                for(S32 i=0; i<N_CHILDREN; i++, adr_tc_child++){
                    if( pos_target_box_out.contained(tc_first[adr_tc_child].mom_.getVertexIn()) 
                        || pos_target_box_in.contained(tc_first[adr_tc_child].mom_.getVertexOut()) ){
                        adr_tc_stack.push(adr_tc_child);
                    }
                }
            }
        }
        return n_cell_open;
    }
#endif

    template<class Ttc, class Tep2, class Tep3>
    inline CountT MakeListUsingInnerBoundaryIteration
        (const Ttc * tc_first,
         const S32 adr_tc,
         const Tep2 * ep_first,
         ReallocatableArray<Tep3> & ep_list,
         const F64ort & pos_target_box, // position of domain
         const S32 n_leaf_limit,
         const F64vec & shift = F64vec(0.0) ){
        //asm("# MakeListUsingInnerBoundaryIteration");
        CountT n_cell_open = 0;
        //static __thread Stack<S32, 1000> adr_tc_stack;
        Stack<S32> adr_tc_stack;
        adr_tc_stack.init();
        adr_tc_stack.push(adr_tc);  // push root cell
        while( !adr_tc_stack.empty() ){
            n_cell_open++;
            S32 adr_tc_curr = adr_tc_stack.pop();
            if( tc_first[adr_tc_curr].isLeaf(n_leaf_limit) ){
                const S32 n_child = tc_first[adr_tc_curr].n_ptcl_;
                ep_list.reserveEmptyAreaAtLeast( n_child );
                U32 adr_ptcl_tmp = tc_first[adr_tc_curr].adr_ptcl_;
                for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                    ep_list.pushBackNoCheck(ep_first[adr_ptcl_tmp]);
                    const F64vec pos_new = ep_first[adr_ptcl_tmp].getPos() + shift; // for periodic mode
                    ep_list.back().setPos(pos_new);
                }
            }
            else{
                S32 adr_tc_child = tc_first[adr_tc_curr].adr_tc_;
                for(S32 i=0; i<N_CHILDREN; i++, adr_tc_child++){
                    if( pos_target_box.contained( tc_first[adr_tc_child].mom_.getVertexIn()) ){
                        adr_tc_stack.push(adr_tc_child);
                    }
                }
            }
        }
        return n_cell_open;
    }


    template<class Ttc, class Tep2, class Tep3>
    inline void MakeListUsingOuterBoundaryIteration
        (const Ttc * tc_first,
         const S32 adr_tc,
         const Tep2 * ep_first,
         ReallocatableArray<Tep3> & ep_list,
         const F64ort & pos_target_box, // position of domain
         const S32 n_leaf_limit,
         const F64vec & shift = F64vec(0.0) ){
        CountT n_cell_open = 0;
        //static __thread Stack<S32, 1000> adr_tc_stack;
        Stack<S32> adr_tc_stack;
        adr_tc_stack.init();
        adr_tc_stack.push(adr_tc);  // push root cell
        while( !adr_tc_stack.empty() ){
            n_cell_open++;
            S32 adr_tc_curr = adr_tc_stack.pop();
            if( tc_first[adr_tc_curr].isLeaf(n_leaf_limit) ){
                const S32 n_child = tc_first[adr_tc_curr].n_ptcl_;
                ep_list.reserveEmptyAreaAtLeast( n_child );
                U32 adr_ptcl_tmp = tc_first[adr_tc_curr].adr_ptcl_;
                for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
#ifdef ORIGINAL_SCATTER_MODE
                    const F64vec pos_tmp = ep_first[adr_ptcl_tmp].getPos();
                    const F64 size_tmp = ep_first[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq_tmp = pos_target_box.getDistanceMinSQ(pos_tmp);
                    if(dis_sq_tmp > size_tmp*size_tmp) continue;
#endif
                    ep_list.increaseSize();
                    ep_list.back() = ep_first[adr_ptcl_tmp];
                    const F64vec pos_new = ep_list.back().getPos() + shift;
                    ep_list.back().setPos(pos_new);
                }
            }
            else{
                S32 adr_tc_child = tc_first[adr_tc_curr].adr_tc_;
                for(S32 i=0; i<N_CHILDREN; i++, adr_tc_child++){
                    if( pos_target_box.contained( tc_first[adr_tc_child].mom_.getVertexOut()) ){
                        adr_tc_stack.push(adr_tc_child);
                    }
                }
            }
        }
    }

#if 1
    template<class Ttc, class Ttp, class Tep, class Tsp>
    inline CountT MakeInteractionListLongEPSPIteration
        (const ReallocatableArray<Ttc> & tc_first,
         const S32 adr_tc,
         const ReallocatableArray<Ttp> & tp_first,
         const ReallocatableArray<Tep> & ep_first,
         ReallocatableArray<Tep> & ep_list,
         const ReallocatableArray<Tsp> & sp_first,
         ReallocatableArray<Tsp> & sp_list,
         const F64ort & pos_target_box,
         const F64 r_crit_sq,
         const S32 n_leaf_limit){
        CountT n_cell_open = 0;
        if( tc_first[adr_tc].n_ptcl_ == 0) return n_cell_open;
        else if( tc_first[adr_tc].isLeaf(n_leaf_limit) ){
            S32 adr_ptcl_tmp = tc_first[adr_tc].adr_ptcl_;
            const S32 n_child = tc_first[adr_tc].n_ptcl_;
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
            return n_cell_open;
        }
        //static __thread Stack<S32, 1000> adr_tc_stack;
        Stack<S32> adr_tc_stack;
        //static __thread Stack<F64, 1000> r_crit_sq_stack;
        Stack<F64> r_crit_sq_stack;
        adr_tc_stack.init();
        r_crit_sq_stack.init();
        adr_tc_stack.push(tc_first[adr_tc].adr_tc_);  // push first child
        r_crit_sq_stack.push(r_crit_sq*0.25);  // push firt child
        while( !adr_tc_stack.empty() ){
            // root cell has children
            n_cell_open++;
            U32 open_bits = 0;
            const S32 adr_tc_curr = adr_tc_stack.pop();
            const F64 r_crit_sq_curr = r_crit_sq_stack.pop();
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
                const F64vec & mom0 = tc_first[adr_tc_curr+i0].mom_.getPos();
                const F64  mom0_x = mom0.x;
                const F64  mom0_y = mom0.y;
                const F64  mom0_z = mom0.z;
                const F64vec & mom1 = tc_first[adr_tc_curr+i1].mom_.getPos();
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
                open_bits |= (( dr_sq[i] <= r_crit_sq_curr) << i);
            }
#else //FAST_WALK_K
            for(S32 i=0; i<N_CHILDREN; i++){
                const F64vec pos_tmp = tc_first[adr_tc_curr+i].mom_.getPos();
                open_bits |= ( (pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq_curr) << i);
            }
#endif
            for(S32 i=0; i<N_CHILDREN; i++){
                const Ttc * tc_curr = tc_first.getPointer(adr_tc_curr+i);
                const S32 n_child = tc_curr->n_ptcl_;
                if(n_child == 0) continue;
                else if( (open_bits>>i) & 0x1 ){
                    if( !(tc_curr->isLeaf(n_leaf_limit)) ){
                        adr_tc_stack.push( tc_curr->adr_tc_ );
                        r_crit_sq_stack.push(r_crit_sq_curr*0.25);
                    }
                    else{
                        S32 adr_ptcl_tmp = tc_curr->adr_ptcl_;
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
                    sp_list.back().copyFromMoment(tc_curr->mom_);
                }
            }
        }
        return n_cell_open;
    }
#else
    template<class Ttc, class Ttp, class Tep, class Tsp>
    inline CountT MakeInteractionListLongEPSPIteration
        (const ReallocatableArray<Ttc> & tc_first,
         const S32 adr_tc,
         const ReallocatableArray<Ttp> & tp_first,
         const ReallocatableArray<Tep> & ep_first,
         ReallocatableArray<Tep> & ep_list,
         const ReallocatableArray<Tsp> & sp_first,
         ReallocatableArray<Tsp> & sp_list,
         const F64ort & pos_target_box,
         const F64 r_crit_sq,
         const S32 n_leaf_limit){

        CountT n_cell_open = 0;
        //static __thread Stack<S32, 1000> adr_tc_stack;
        Stack<S32> adr_tc_stack;
        //static __thread Stack<F64, 1000> r_crit_sq_stack;
        Stack<F64> r_crit_sq_stack;
        adr_tc_stack.init();
        r_crit_sq_stack.init();
        adr_tc_stack.push(adr_tc);  // push root cell
        r_crit_sq_stack.push(r_crit_sq);  // push root cell
        if( tc_first[adr_tc].n_ptcl_ == 0) return n_cell_open;
        while( !adr_tc_stack.empty() ){
            n_cell_open++;
            S32 adr_tc_curr = adr_tc_stack.pop();
            F64 r_crit_sq_curr = r_crit_sq_stack.pop();
            if( tc_first[adr_tc_curr].isLeaf(n_leaf_limit) ){
                const S32 n_child = tc_first[adr_tc_curr].n_ptcl_;
                S32 adr_ptcl_tmp = tc_first[adr_tc_curr].adr_ptcl_;
                ep_list.reserveEmptyAreaAtLeast( n_child );
                sp_list.reserveEmptyAreaAtLeast( n_child );
                for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                    if( GetMSB(tp_first[adr_ptcl_tmp].adr_ptcl_) == 0){
                        ep_list.pushBackNoCheck(ep_first[adr_ptcl_tmp]);
                    }
                    else{
                        sp_list.pushBackNoCheck(sp_first[adr_ptcl_tmp]);
                    }
                }
            }
            else{
                S32 adr_tc_child = tc_first[adr_tc_curr].adr_tc_;
                const F64 r_crit_sq_child = r_crit_sq_curr * 0.25;
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
                    const F64vec & mom0 = tc_first[adr_tc_child+i0].mom_.getPos();
                    const F64  mom0_x = mom0.x;
                    const F64  mom0_y = mom0.y;
                    const F64  mom0_z = mom0.z;
                    const F64vec & mom1 = tc_first[adr_tc_child+i1].mom_.getPos();
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
                for(S32 i=0; i<N_CHILDREN; i++, adr_tc_child++){
                    if(tc_first[adr_tc_child].n_ptcl_ == 0) continue;
                    else if( dr_sq[i] <= r_crit_sq_child ){
                        adr_tc_stack.push(adr_tc_child);
                        r_crit_sq_stack.push(r_crit_sq_child);
                    }
                    else{
                        sp_list.increaseSize();
                        sp_list.back().copyFromMoment(tc_first[adr_tc_child].mom_);
                    }
                }
#else //FAST_WALK_K
                for(S32 i=0; i<N_CHILDREN; i++, adr_tc_child++){
                    if(tc_first[adr_tc_child].n_ptcl_ == 0) continue;
                    const F64vec pos_tmp = tc_first[adr_tc_child].mom_.getPos();
                    if( pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq_child ){
                        adr_tc_stack.push(adr_tc_child);
                        r_crit_sq_stack.push(r_crit_sq_child);
                    }
                    else{
                        sp_list.increaseSize();
                        sp_list.back().copyFromMoment(tc_first[adr_tc_child].mom_);
                    }
                }
#endif //FAST_WALK_K
            }
        }
        return n_cell_open;
    }
#endif
}

#include"tree_for_force_utils2.hpp"
    
