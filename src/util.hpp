
namespace  ParticleSimulator{

    ////////////////////
    ///   calc distance
    //////////////////
    // between points
    PS_INLINE F64 GetDistanceMin1DImpl(const F64 a,
				       const F64 b){
        return std::fabs(a-b);
    }

    PS_INLINE F64 GetDistanceMin1DPeriImpl(const F64 a,
					   const F64 b,
					   const F64 len){
        const auto h0 = std::fabs(a-b);
	const auto h1 = len - h0;
	return std::min(h0, h1);
    }    

    template<int bit>  
    PS_INLINE F64 GetDistanceMin1D(const F64 a,
				   const F64 b,
				   const F64 len){
        Abort();
	return -1;
    }
    template<>  
    PS_INLINE F64 GetDistanceMin1D<0>(const F64 a,
				      const F64 b,
				      const F64 len){
        return GetDistanceMin1DImpl(a, b);
    }

    template<>  
    PS_INLINE F64 GetDistanceMin1D<1>(const F64 a,
				      const F64 b,
				      const F64 len){
      return GetDistanceMin1DPeriImpl(a, b, len);
    }
    
    //////////////////
    // between lines
    PS_INLINE F64 GetDistanceMin1DImpl(const F64 al,
				       const F64 ah,
				       const F64 bl,
				       const F64 bh){
        const auto dx0 = al - bh;
        const auto dx1 = bl - ah;
        const auto dx = (dx0*dx1 >= 0.0) ? 0.0 : ((dx0>0.0) ? dx0 : dx1);
        return dx;
    }

    PS_INLINE F64 GetDistanceMin1DPeriImpl(const F64 al,
					   const F64 ah,
					   const F64 bl,
					   const F64 bh,
					   const F64 len){
        const auto dx0 = al - bh;
        const auto dx1 = bl - ah;
        const auto dis0 = std::max(dx0, dx1);
	const auto dis1 = std::max( (len + std::min(dx0, dx1)), 0.0 );
	const auto dis  = (dx0*dx1 >= 0.0) ? 0.0 : std::min(dis0, dis1);
	return dis;
      /*
        const auto dx0 = GetDistanceMin1DImpl(al+len, ah+len, bl, bh);
        const auto dx1 = GetDistanceMin1DImpl(al-len, ah-len, bl, bh);
        const auto dx2 = GetDistanceMin1DImpl(al, ah, bl, bh);
        return std::min(std::min(dx0, dx1), dx2);
      */
    }

    template<int bit>
    PS_INLINE F64 GetDistanceMin1D(const F64 al,
				   const F64 ah,
				   const F64 bl,
				   const F64 bh,
				   const F64 len){
        Abort();
	return -1;
    }
    template<>
    PS_INLINE F64 GetDistanceMin1D<0>(const F64 al,
				      const F64 ah,
				      const F64 bl,
				      const F64 bh,
				      const F64 len){
      return GetDistanceMin1DImpl(al, ah, bl, bh);
    }
    template<>
    PS_INLINE F64 GetDistanceMin1D<1>(const F64 al,
				      const F64 ah,
				      const F64 bl,
				      const F64 bh,
				      const F64 len){
      return GetDistanceMin1DPeriImpl(al, ah, bl, bh, len);
    }

    //////////////////
    // between point and line
    PS_INLINE F64 GetDistanceMin1DImpl(const F64 al,
				       const F64 ah,
				       const F64 b){
        const auto dx0 = b - ah;
        const auto dx1 = al - b;
        const auto dx = (dx0*dx1 >= 0.0) ? 0.0 : ((dx0>0.0) ? dx0 : dx1);
        return dx;
    }

    PS_INLINE F64 GetDistanceMin1DPeriImpl(const F64 al,
					   const F64 ah,
					   const F64 b,
					   const F64 len){
        const auto dx0 = b - ah;
        const auto dx1 = al - b;
	const auto dis0 = std::max(dx0, dx1);
	//const auto dis1 = len + std::min(dx0, dx1);
	const auto dis1 = std::max((len + std::min(dx0, dx1)), 0.0);
	const auto dis  = (dx0*dx1 >= 0.0) ? 0.0 : std::min(dis0, dis1);
	return dis;
        /*
	const auto h0 = GetDistanceMin1DPeriImpl(al, b, len);
	const auto h1 = GetDistanceMin1DPeriImpl(ah, b, len);
	return std::min(h0, h1);
        */
    }

    template<int bit>
    PS_INLINE F64 GetDistanceMin1D(const F64 al,
				   const F64 ah,
				   const F64 b,
				   const F64 len){
        Abort();
	return -1;
    }
    template<>
    PS_INLINE F64 GetDistanceMin1D<0>(const F64 al,
				      const F64 ah,
				      const F64 b,
				      const F64 len){
        return GetDistanceMin1DImpl(al, ah, b);
    }
    template<>
    PS_INLINE F64 GetDistanceMin1D<1>(const F64 al,
				      const F64 ah,
				      const F64 b,
				      const F64 len){
        return GetDistanceMin1DPeriImpl(al, ah, b, len);
    }
  
    //////////////////
    // between F64vec
    PS_INLINE F64 GetDistanceMinSq(const F64vec & pos0,
				   const F64vec & pos1){
      const auto dpos = pos0 - pos1;
      return dpos*dpos;
    }

    template<int bits>
    PS_INLINE F64 GetDistanceMinSq(const F64vec & pos0,
				   const F64vec & pos1,
				   const F64vec & len_peri){
        const auto dx = GetDistanceMin1D<((bits>>0)&0x1)>(pos0.x, pos1.x, len_peri.x);
	const auto dy = GetDistanceMin1D<((bits>>1)&0x1)>(pos0.y, pos1.y, len_peri.y);
        auto dis = dx*dx + dy*dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	const auto dz = GetDistanceMin1D<((bits>>2)&0x1)>(pos0.z, pos1.z, len_peri.z);
        dis += dz*dz;
#endif
        return dis;
    }
  
    //////////////////
    // between F64ort
    PS_INLINE F64 GetDistanceMinSq(const F64ort & pos0,
				   const F64ort & pos1){
        const auto dx = GetDistanceMin1DImpl(pos0.low_.x, pos0.high_.x, pos1.low_.x, pos1.high_.x);
        const auto dy = GetDistanceMin1DImpl(pos0.low_.y, pos0.high_.y, pos1.low_.y, pos1.high_.y);
        auto dis = dx*dx + dy*dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        const auto dz = GetDistanceMin1DImpl(pos0.low_.z, pos0.high_.z, pos1.low_.z, pos1.high_.z);
        dis += dz*dz;
#endif
        return dis;
    }

    template<int bits=0>
    PS_INLINE F64 GetDistanceMinSq(const F64ort & pos0,
				   const F64ort & pos1,
				   const F64vec & len_peri){
        const auto dx = GetDistanceMin1D<((bits>>0)&0x1)>(pos0.low_.x, pos0.high_.x, pos1.low_.x, pos1.high_.x, len_peri.x);
        const auto dy = GetDistanceMin1D<((bits>>1)&0x1)>(pos0.low_.y, pos0.high_.y, pos1.low_.y, pos1.high_.y, len_peri.y);
        auto dis = dx*dx + dy*dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        const auto dz = GetDistanceMin1D<((bits>>2)&0x1)>(pos0.low_.z, pos0.high_.z, pos1.low_.z, pos1.high_.z, len_peri.z);
        dis += dz*dz;
#endif
        return dis;
    }

    //////////////////
    // between F64ort and F64vec
    PS_INLINE F64 GetDistanceMinSq(const F64ort & pos0,
				   const F64vec & pos1){
        const auto dx = GetDistanceMin1DImpl(pos0.low_.x, pos0.high_.x, pos1.x);
        const auto dy = GetDistanceMin1DImpl(pos0.low_.y, pos0.high_.y, pos1.y);
        auto dis = dx*dx + dy*dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        const F64 dz = GetDistanceMin1DImpl(pos0.low_.z, pos0.high_.z, pos1.z);
        dis += dz*dz;
#endif
        return dis;
    }

    PS_INLINE F64 GetDistanceMinSq(const F64vec & pos0,
				   const F64ort & pos1){
        return GetDistanceMinSq(pos1, pos0);
    }

    template<int bits>
    PS_INLINE F64 GetDistanceMinSq(const F64ort & pos0,
				   const F64vec & pos1,
				   const F64vec & len_peri){
        const auto dx = GetDistanceMin1D<((bits>>0)&0x1)>(pos0.low_.x, pos0.high_.x, pos1.x, len_peri.x);
        const auto dy = GetDistanceMin1D<((bits>>1)&0x1)>(pos0.low_.y, pos0.high_.y, pos1.y, len_peri.y);
        auto dis = dx*dx + dy*dy;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        const F64 dz = GetDistanceMin1D<((bits>>2)&0x1)>(pos0.low_.z, pos0.high_.z, pos1.z, len_peri.z);
        dis += dz*dz;
#endif
        return dis;
    }

    template<int bits>
    PS_INLINE F64 GetDistanceMinSq(const F64vec & pos0,
				   const F64ort & pos1,
				   const F64vec & len_peri){
        return GetDistanceMinSq<bits>(pos1, pos0, len_peri);
    }


    //////////////////
    // over lapped functions
#if __cplusplus >= PS_CPLUSPLUS_17
    template<typename ...>
    constexpr bool FALSE_V = false;
    template<typename T0, typename T1>
    PS_INLINE bool IsOverlapped(const T0 & pos0,
                                const T1 & pos1){
        if constexpr( (std::is_same_v<T0, F64vec> && std::is_same_v<T1, F64vec>)
                      || (std::is_same_v<T0, F64ort> && std::is_same_v<T1, F64ort>)
                      || (std::is_same_v<T0, F64vec> && std::is_same_v<T1, F64ort>)
                      || (std::is_same_v<T0, F64ort> && std::is_same_v<T1, F64vec>)){
                //asm("#IF TRUE");
            return GetDistanceMinSq(pos0, pos1) <= 0.0;
        }
        else{
            //asm("#IF FALSE");
            static_assert(FALSE_V<T0, T1>);
        }
    }
    template<int bits, typename T0, typename T1>
    PS_INLINE bool IsOverlapped(const T0 & pos0,
                                const T1 & pos1,
                                const F64vec & len_peri){
        if constexpr( (std::is_same_v<T0, F64vec> && std::is_same_v<T1, F64vec>)
                      || (std::is_same_v<T0, F64ort> && std::is_same_v<T1, F64ort>)
                      || (std::is_same_v<T0, F64vec> && std::is_same_v<T1, F64ort>)
                      || (std::is_same_v<T0, F64ort> && std::is_same_v<T1, F64vec>)){
                return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
        }
        else{
            static_assert(FALSE_V<T0, T1>);
        }
    }
#else
    PS_INLINE bool IsOverlapped(const F64ort & pos0,
                                const F64ort & pos1){
        return GetDistanceMinSq(pos0, pos1) <= 0.0;
    }
    PS_INLINE bool IsOverlapped(const F64vec & pos0,
                                const F64vec & pos1){
        return GetDistanceMinSq(pos0, pos1) <= 0.0;
    }
    PS_INLINE bool IsOverlapped(const F64ort & pos0,
                                const F64vec & pos1){
        return GetDistanceMinSq(pos0, pos1) <= 0.0;
    }
    PS_INLINE bool IsOverlapped(const F64vec & pos0,
                                const F64ort & pos1){
        return GetDistanceMinSq(pos0, pos1) <= 0.0;
    }

    template<int bits>
    PS_INLINE bool IsOverlapped(const F64ort & pos0,
                                const F64ort & pos1,
                                const F64vec & len_peri){
        return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
    }
    template<int bits>
    PS_INLINE bool IsOverlapped(const F64vec & pos0,
                                const F64vec & pos1,
                                const F64vec & len_peri){
        return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
    }
    template<int bits>
    PS_INLINE bool IsOverlapped(const F64ort & pos0,
                                const F64vec & pos1,
                                const F64vec & len_peri){
        return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
    }
    template<int bits>
    PS_INLINE bool IsOverlapped(const F64vec & pos0,
                                const F64ort & pos1,
                                const F64vec & len_peri){
        return GetDistanceMinSq<bits>(pos0, pos1, len_peri) <= 0.0;
    }
#endif
    
    struct CopyPos{
        template<typename Tsrc, typename Tdst> 
        void operator()(Tsrc & src, Tdst & dst){
            dst = src.getPos();
        }
    };

    
    struct Cmpvec{
	F64 F64vec::*loc;
	Cmpvec(F64 F64vec::*_loc) : loc(_loc) {}
	bool operator()(const F64vec &lhs, const F64vec &rhs){
	    return (lhs.*loc < rhs.*loc);
	}
    };
#ifdef ENABLE_AVX2
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
    
    static inline __m256 GetDistanceMinSQV8R4(const __m256 & ix,
					      const __m256 & iy,
					      const __m256 & iz,
					      const __m256 & jx,
					      const __m256 & jy,
					      const __m256 & jz,
					      const __m256 & lx,
					      const __m256 & ly,
					      const __m256 & lz){
        static const __m256 mask_abs = (__m256)_mm256_set1_epi32(0x7fffffff);
        static const __m256 zero = _mm256_set1_ps(0.0f);
	
        __m256 dx = _mm256_and_ps(mask_abs, (ix - jx));
        dx = dx - lx;
        dx = _mm256_max_ps(dx, zero);
        __m256 dy = _mm256_and_ps(mask_abs, (iy - jy));
        dy = dy - ly;
        dy = _mm256_max_ps(dy, zero);
        __m256 dz = _mm256_and_ps(mask_abs, (iz - jz));
        dz = dz - lz;
        dz = _mm256_max_ps(dz, zero);
	return _mm256_fmadd_ps(dz, dz, _mm256_fmadd_ps(dy, dy, (dx*dx)));
        //__m256 rsq = _mm256_fmadd_ps(dz, dz, _mm256_fmadd_ps(dy, dy, (dx*dx)));	
    }

    static inline __m256 GetDistanceMinSQV8R4(const F32ort & box, const F32vec pos[8]){
        const __m256 half = _mm256_set1_ps(0.5);
        const __m256 high = _mm256_loadu_ps(&box.high_.x);
        const __m256 low  = _mm256_loadu_ps(&box.low_.x);
        const __m256 cen = _mm256_mul_ps(_mm256_add_ps(high, low), half);
        const __m256 len = _mm256_mul_ps(_mm256_sub_ps(high, low), half);
        const __m256 ix = _mm256_shuffle_ps(cen, cen, 0x00);
        const __m256 iy = _mm256_shuffle_ps(cen, cen, 0x55);
        const __m256 iz = _mm256_shuffle_ps(cen, cen, 0xaa);
        const __m256 lx = _mm256_shuffle_ps(len, len, 0x00);
        const __m256 ly = _mm256_shuffle_ps(len, len, 0x55);
        const __m256 lz = _mm256_shuffle_ps(len, len, 0xaa);
        __m256 jx = _mm256_set_m128_forgcc(_mm_loadu_ps(&(pos[4].x)),
                                           _mm_loadu_ps(&(pos[0].x)));
        __m256 jy = _mm256_set_m128_forgcc(_mm_loadu_ps(&(pos[5].x)), 
                                           _mm_loadu_ps(&(pos[1].x)));
        __m256 jz = _mm256_set_m128_forgcc(_mm_loadu_ps(&(pos[6].x)), 
                                           _mm_loadu_ps(&(pos[2].x)));
        __m256 jw = _mm256_set_m128_forgcc(_mm_loadu_ps(&(pos[7].x)), 
                                           _mm_loadu_ps(&(pos[3].x)));
        transpose_4ymm(jx, jy, jz, jw);
	return GetDistanceMinSQV8VF(ix, iy, iz, jx, jy, jz, lx, ly, lz);
    }
#endif

#ifdef ENABLE_HPC_ACE
    static inline __m256 GetDistanceMinSQV2R8(const _fjsp_v2r8 & ix,
					      const _fjsp_v2r8 & iy,
					      const _fjsp_v2r8 & iz,
					      const _fjsp_v2r8 & jx,
					      const _fjsp_v2r8 & jy,
					      const _fjsp_v2r8 & jz,
					      const _fjsp_v2r8 & lx,
					      const _fjsp_v2r8 & ly,
					      const _fjsp_v2r8 & lz){
        static const _fjsp_v2r8 zero = _fjsp_setzero_v2r8();
	
        _fjsp_v2r8 dx = _mm256_and_ps(mask_abs, (ix - jx));
        dx = dx - lx;
        dx = _mm256_max_ps(dx, zero);
        _fjsp_v2r8 dy = _mm256_and_ps(mask_abs, (iy - jy));
        dy = dy - ly;
        dy = _mm256_max_ps(dy, zero);
        _fjsp_v2r8 dz = _mm256_and_ps(mask_abs, (iz - jz));
        dz = dz - lz;
        dz = _mm256_max_ps(dz, zero);
	return _mm256_fmadd_ps(dz, dz, _mm256_fmadd_ps(dy, dy, (dx*dx)));
        //__m256 rsq = _mm256_fmadd_ps(dz, dz, _mm256_fmadd_ps(dy, dy, (dx*dx)));	
    }
#endif

    static inline void CalcNumberAndShiftOfImageDomain
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

    static inline S32vec CalcIDOfImageDomain(const F64ort & pos_root_domain,
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

    template<typename T>
    static inline void PrefixSum(const ReallocatableArray<T> & n, ReallocatableArray<T> & n_disp){
        n_disp.resizeNoInitialize(n.size()+1);
        n_disp[0] = 0;
        for(int i=0; i<n.size(); i++){
            n_disp[i+1] = n_disp[i] + n[i];
        }
    }
    
    template<typename T>
    static inline void DeleteArray(T * array){
        if(array != nullptr) delete [] array;
        array = nullptr;
    }

    static PS_INLINE S32vec GetRankVec(const S32 n_domain[DIMENSION_LIMIT],
				const S32 rank_glb,
				const S32 dim){
	S32vec rank_1d;
	S32 rank_tmp = rank_glb;
	S32 div = 1;
	for(auto i=dim-1; i>=0; i--){
	    rank_1d[i] = rank_tmp % n_domain[i];
	    div *= n_domain[i];
	    rank_tmp = rank_glb / div;
	}
	return rank_1d;
    }

    static PS_INLINE S32 GetRankGlb(const S32vec rank_vec,
			     const S32 n_domain[DIMENSION_LIMIT]){
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
	return rank_vec[0]*n_domain[1] + rank_vec[1];
#else
	return rank_vec[0]*n_domain[1]*n_domain[2] + rank_vec[1]*n_domain[2] + rank_vec[2];
#endif
    }

    static PS_INLINE S32 GetRank1d(const S32 rank_glb,
				   const S32 n_domain[DIMENSION_LIMIT],
				   const S32 cid){
	S32 ret = -1;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
	if(cid==0)
	    ret = rank_glb / n_domain[1];
	else if(cid==1)
	    ret = rank_glb % n_domain[1];
#else
	if(cid==0)
	    ret = (rank_glb / n_domain[2]) / n_domain[1];
	else if(cid==1)
	    ret = (rank_glb / n_domain[2]) % n_domain[1];
	else if(cid==2)
	    ret = rank_glb % n_domain[2];
#endif
	return ret;
    }

    static inline void GetRankVec(S32 * rank_vec,
			   const S32 * n_domain,
			   const S32 rank_glb,
			   const S32 dim){
	S32 rank_tmp = rank_glb;
	S32 div = 1;
	for(auto i=dim-1; i>=0; i--){
	    rank_vec[i] = rank_tmp % n_domain[i];
	    div *= n_domain[i];
	    rank_tmp = rank_glb / div;
	}
    }
    
    static inline S32 GetRankGlb(const S32 * rank_vec,
			  const S32 * n_domain,
			  const S32 dim){
	S32 rank_new = 0;
	S32 radix = 1;
	for(auto i=dim-1; i>=0; i--){
	    rank_new += radix*rank_vec[i];
	    radix *= n_domain[i];
	}
	return rank_new;
    }
    
    static inline S32 GetRank1d(const S32 rank_glb,
			 const S32 * n_domain,
			 const S32 dim,
			 const S32 cid){
	S32 rank_new = rank_glb;
	for(auto i=dim-1; i>cid; i--){
	    rank_new /= n_domain[i];
	}
	return rank_new % n_domain[cid];
    }    

    static PS_INLINE void CalcRank1D(S32 & rank_1d, S32 & dnp_1d, F64 & dis_domain_1d, const F64vec & pos, const S32 rank_glb_org, const S32 n_domain[DIMENSION_LIMIT], const F64ort domain[], const F64 len_peri, const bool pa, const int cid){
	// rank_glb_org don't have to be my_rank. It is just a start rank to search new rank.
	const auto pos_1d = pos[cid];
	const auto rank_vec = GetRankVec(n_domain, rank_glb_org, DIMENSION);
	const auto rank_1d_org = GetRank1d(rank_glb_org, n_domain, cid);
	auto rank_glb_tmp = rank_glb_org;
	rank_1d           = rank_1d_org;
	S32 shift = 1;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
	if(cid == 0) shift = n_domain[1];
#else
	if(cid == 1) shift = n_domain[2];
	else if(cid == 0) shift = n_domain[1] * n_domain[2];
#endif
	if(pa == true){
	    // periodic
	    auto rank_vec_tmp = rank_vec;
	    //const auto rank_glb_max = GetRankGlb(rank_vec, n_domain);
	    if( domain[rank_glb_org].low_[cid] > pos_1d && (domain[rank_glb_org].low_[cid] - pos_1d) > 0.5*len_peri){
		rank_vec_tmp[cid] = 0;
	    }
	    else if( domain[rank_glb_org].high_[cid] <= pos_1d && pos_1d-domain[rank_glb_org].high_[cid] > 0.5*len_peri){
		rank_vec_tmp[cid] = n_domain[cid]-1;
	    }
	    rank_glb_tmp = GetRankGlb(rank_vec_tmp, n_domain);
	}

	while( pos_1d < domain[rank_glb_tmp].low_[cid] ){
	    rank_glb_tmp -= shift;
	}
	while( pos_1d >= domain[rank_glb_tmp].high_[cid] ){
	    rank_glb_tmp += shift;
	}

	if(pa == true)
	    dis_domain_1d = GetDistanceMin1DPeriImpl(domain[rank_glb_tmp].low_[cid], domain[rank_glb_tmp].high_[cid], pos_1d, len_peri);
	else
	    dis_domain_1d = GetDistanceMin1DImpl(domain[rank_glb_tmp].low_[cid], domain[rank_glb_tmp].high_[cid], pos_1d);

	rank_1d = GetRank1d(rank_glb_tmp, n_domain, cid);
	dnp_1d = std::abs(rank_1d-rank_1d_org);
	auto dnp_1d2  = n_domain[cid]-dnp_1d;
	if (dnp_1d > dnp_1d2) dnp_1d = dnp_1d2;
    }

    static inline int GetColorForCommSplit(const int my_rank_glb, const int * n_domain, const int dim, const int cid){
	const auto my_rank_1d = GetRank1d(my_rank_glb, n_domain, cid);
	int tmp = 1;
	for(auto i=dim-1; i>cid; i--){
	    tmp *= n_domain[i];
	}
	return (my_rank_glb - (my_rank_1d * tmp));
    }



#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
    template<typename T>
    static void MyAlltoall(T * send_buf, int cnt, T * recv_buf, CommInfo & comm_info, S32 n_domain[], S32 dim){
	const auto n_proc_glb = comm_info.getNumberOfProc();
	const auto my_rank_glb = comm_info.getRank();
	const auto n_recv_tot = cnt * n_proc_glb;
	ReallocatableArray<T> send_buf_tmp(n_recv_tot, n_recv_tot, MemoryAllocMode::Stack);
	for(auto i=0; i<n_recv_tot; i++){
	    recv_buf[i] = send_buf[i];
	}
	for(auto d=dim-1; d>=0; d--){
	    const auto radix = n_proc_glb / n_domain[d];
	    const auto color = GetColorForCommSplit(my_rank_glb, n_domain, dim, d);
	    const auto comm_1d  = comm_info.split(color, my_rank_glb);	    
	    for(auto i=0; i<n_proc_glb; i++){
		const int id_send = cnt * ( (i % radix) * n_domain[d] + i / radix );
		const int offset = i * cnt;
		for(auto j=0; j<cnt; j++){
		    send_buf_tmp[offset + j] = recv_buf[id_send + j];
		}
	    }
	    MPI_Alltoall(send_buf_tmp.getPointer(), cnt*radix, GetDataType<T>(),
			 recv_buf, cnt*radix, GetDataType<T>(), comm_1d.getCommunicator());
	}
    }

    template<typename T>
    static void MyAlltoallv(const ReallocatableArray<T> & send_buf, const int * n_send,
			    ReallocatableArray<T> & recv_buf, int * n_recv, 
			    CommInfo & comm_info, S32 n_domain[], S32 dim){
	const auto n_proc_glb  = comm_info.getNumberOfProc();
	const auto my_rank_glb = comm_info.getRank();
	ReallocatableArray<int> n_send_tmp(n_proc_glb, n_proc_glb, MemoryAllocMode::Default);
	ReallocatableArray<int> n_disp_recv_tmp(n_proc_glb+1, n_proc_glb+1, MemoryAllocMode::Default);
	n_disp_recv_tmp[0] = 0;
	for(auto i=0; i<n_proc_glb; i++){
	    n_recv[i] = n_send[i];
	    n_disp_recv_tmp[i+1] = n_disp_recv_tmp[i] + n_recv[i];
	}
	recv_buf.resizeNoInitialize(n_disp_recv_tmp[n_proc_glb]);
	for(auto i=0; i<n_disp_recv_tmp[n_proc_glb]; i++){
	    recv_buf[i] = send_buf[i];
	}
	// main loop
	for(auto d=dim-1; d>=0; d--){
	    const auto radix = n_proc_glb / n_domain[d];
	    const auto color = GetColorForCommSplit(my_rank_glb, n_domain, dim, d);
	    const auto comm_1d  = comm_info.split(color, my_rank_glb);
	    
	    n_disp_recv_tmp[0] = 0;
	    for(auto i=0; i<n_proc_glb; i++){
		n_disp_recv_tmp[i+1] = n_disp_recv_tmp[i] + n_recv[i];
	    }
	    recv_buf.resizeNoInitialize(n_disp_recv_tmp[n_proc_glb]);
	    
	    ReallocatableArray<T> send_buf_tmp(n_disp_recv_tmp[n_proc_glb], 0, MemoryAllocMode::Stack);
	    for(auto i=0; i<n_proc_glb; i++){
		const int id_send = ( (i % radix) * n_domain[d] + i / radix );
		//if(my_rank_glb==0){
		//    std::cerr<<"i= "<<i<<" id_send= "<<id_send<<" n_disp_recv_tmp[id_send]= "<<n_disp_recv_tmp[id_send]<<std::endl;
		//}
		n_send_tmp[i] = n_recv[id_send];
		for(auto j=0; j<n_recv[id_send]; j++){
		    send_buf_tmp.pushBackNoCheck(recv_buf[n_disp_recv_tmp[id_send]+j]);
		}
	    }
	    if(my_rank_glb==0){
		//std::cerr<<"send_buf_tmp.size()= "<<send_buf_tmp.size()<<std::endl;
		//for(auto i=0; i<send_buf_tmp.size(); i++){
		//    std::cerr<<"i= "<<i<<" send_buf_tmp[i].id= "<<send_buf_tmp[i].id<<std::endl;
		//}
	    }
	    
	    MPI_Alltoall(n_send_tmp.getPointer(), radix, GetDataType<int>(),
			 n_recv, radix, GetDataType<int>(), comm_1d.getCommunicator());

	    ReallocatableArray<int> n_recv_1d(n_domain[d], n_domain[d], MemoryAllocMode::Stack);
	    ReallocatableArray<int> n_send_1d(n_domain[d], n_domain[d], MemoryAllocMode::Stack);
	    ReallocatableArray<int> n_disp_recv_1d(n_domain[d]+1, n_domain[d]+1, MemoryAllocMode::Stack);
	    ReallocatableArray<int> n_disp_send_1d(n_domain[d]+1, n_domain[d]+1, MemoryAllocMode::Stack);
	    n_disp_recv_1d[0] = n_disp_send_1d[0] = 0;
	    for(auto i=0; i<n_domain[d]; i++){
		n_send_1d[i] = n_recv_1d[i] = 0;
		for(auto j=0; j<radix; j++){
		    n_send_1d[i] += n_send_tmp[i*radix+j];
		    n_recv_1d[i] += n_recv[i*radix+j];
		}
		n_disp_send_1d[i+1] = n_disp_send_1d[i] + n_send_1d[i];
		n_disp_recv_1d[i+1] = n_disp_recv_1d[i] + n_recv_1d[i];
	    }
	    recv_buf.resizeNoInitialize( n_disp_recv_1d[n_domain[d]] );
	    //if(my_rank_glb==0){
	    //	std::cerr<<"n_domain[d]= "<<n_domain[d]<<std::endl;
	    //	for(auto i=0; i<n_domain[d]; i++){
	    //	    std::cerr<<"i= "<<i<<" n_send_1d[i]= "<<n_send_1d[i]<<" n_recv_1d[i]= "<<n_recv_1d[i]<<std::endl;
	    //	}
	    //}
	    MPI_Alltoallv(send_buf_tmp.getPointer(), n_send_1d.getPointer(), n_disp_send_1d.getPointer(), GetDataType<T>(),
			  recv_buf.getPointer(), n_recv_1d.getPointer(), n_disp_recv_1d.getPointer(), GetDataType<T>(),
			  comm_1d.getCommunicator());
	    //if(my_rank_glb==0){
	    //	std::cerr<<"recv_buf.size()= "<<recv_buf.size()<<std::endl;
	    //	for(auto i=0; i<recv_buf.size(); i++){
	    //	    std::cerr<<"i= "<<i<<" recv_buf[i].id= "<<recv_buf[i].id<<std::endl;
	    //	}
	    //}
	}
    }

    static inline void Factorize(const S32 num, const S32 dim, S32 fact[]){
	S32 num_tmp = num;
	std::vector<S32> num_vec(dim);
        for(S32 d=dim; d > 0; d--){
	    S32 tmp = (S32)pow((F64)num_tmp+0.000001, (1.0/d)*1.000001 );
	    while(num_tmp%tmp!=0){
		tmp--;
	    }
	    num_vec[d-1] = tmp;
            num_tmp /= num_vec[d-1];
	}
	std::sort(num_vec.begin(), num_vec.end(), std::greater<S32>());
	for(auto k=0; k<dim; k++){
	    fact[k] = num_vec[k];
	}
	S32 num_check = fact[0];
	for(S32 d=1; d < dim; d++){
	    num_check *= fact[d];
	    assert(fact[d-1] >= fact[d]);
	}
	assert(num_check == num);
    }
    
    template<typename T>
    static inline void MyAlltoallReuse(const ReallocatableArray<T> & send_buf, int cnt, ReallocatableArray<T> * send_buf_multi_dim, ReallocatableArray<T> * recv_buf_multi_dim, const S32 n_domain[], const S32 dim, const CommInfo & comm_info){
	const auto n_proc_glb = comm_info.getNumberOfProc();
	const auto my_rank_glb = comm_info.getRank();
	//const auto n_recv_tot = cnt * n_proc_glb;
	for(auto d=dim-1; d>=0; d--){
	    const auto radix = n_proc_glb / n_domain[d];
	    const auto color = GetColorForCommSplit(my_rank_glb, n_domain, dim, d);
	    const auto comm_1d  = comm_info.split(color, my_rank_glb);
	    for(auto i=0; i<n_proc_glb; i++){
		const auto id_send = cnt * ( (i % radix) * n_domain[d] + i / radix );
		const auto offset = i * cnt;
		if(d==dim-1){
		    for(auto j=0; j<cnt; j++){
			send_buf_multi_dim[d][offset + j] = send_buf[id_send + j];
		    }
		} else {
		    for(auto j=0; j<cnt; j++){
			send_buf_multi_dim[d][offset + j] = recv_buf_multi_dim[d+1][id_send + j];
		    }
		}
	    }
	    MPI_Alltoall(send_buf_multi_dim[d].getPointer(), cnt*radix, GetDataType<T>(),
	    recv_buf_multi_dim[d].getPointer(), cnt*radix, GetDataType<T>(), comm_1d.getCommunicator());
	}
    }
    template<typename T>
    static inline void MyAlltoallReuse(const ReallocatableArray<T> & send_buf, int cnt, ReallocatableArray<T> * send_buf_multi_dim, ReallocatableArray<T> * recv_buf_multi_dim, const S32 n_domain[], const S32 dim){
	auto ci = Comm::getCommInfo();
	MyAlltoallReuse(send_buf, cnt, send_buf_multi_dim, recv_buf_multi_dim, n_domain, dim, ci);
    }

    template<typename T>
    static inline void MyAlltoallVReuse(const ReallocatableArray<T> & send_buf, ReallocatableArray<T> & recv_buf,
					const ReallocatableArray<int> * n_send_multi_dim, const ReallocatableArray<int> * n_recv_multi_dim, const S32 n_domain[], const S32 dim, const CommInfo & comm_info){
	const auto n_proc_glb  = comm_info.getNumberOfProc();
	const auto my_rank_glb = comm_info.getRank();
	ReallocatableArray<int> n_recv_tmp(n_proc_glb, n_proc_glb, MemoryAllocMode::Stack);
	ReallocatableArray<int> n_disp_recv_tmp(n_proc_glb+1, n_proc_glb+1, MemoryAllocMode::Stack);
	for(auto i=0; i<n_proc_glb; i++){
	    const auto radix = n_proc_glb / n_domain[dim-1];
	    const auto id_send = ( (i % radix) * n_domain[dim-1] + i / radix );
	    n_recv_tmp[id_send] = n_send_multi_dim[dim-1][i];
	}
	recv_buf.resizeNoInitialize(send_buf.size());
	for(auto i=0; i<send_buf.size(); i++){
	    recv_buf[i] = send_buf[i];
	}
	// main loop
	for(auto d=dim-1; d>=0; d--){
	    const auto radix = n_proc_glb / n_domain[d];
	    const auto color = GetColorForCommSplit(my_rank_glb, n_domain, dim, d);
	    const auto comm_1d  = comm_info.split(color, my_rank_glb);
	    n_disp_recv_tmp[0] = 0;
	    for(auto i=0; i<n_proc_glb; i++){
		n_disp_recv_tmp[i+1] = n_disp_recv_tmp[i] + n_recv_tmp[i];
	    }	    
	    recv_buf.resizeNoInitialize(n_disp_recv_tmp[n_proc_glb]);
	    ReallocatableArray<T> send_buf_tmp(n_disp_recv_tmp[n_proc_glb], 0, MemoryAllocMode::Stack);
	    for(auto i=0; i<n_proc_glb; i++){
		const auto id_send = ( (i % radix) * n_domain[d] + i / radix );
		for(auto j=0; j<n_recv_tmp[id_send]; j++){
		    send_buf_tmp.pushBackNoCheck(recv_buf[n_disp_recv_tmp[id_send]+j]);
		}
	    }
	    for(auto i=0; i<n_proc_glb; i++){
		n_recv_tmp[i] = n_recv_multi_dim[d][i];
	    }
	    ReallocatableArray<int> n_recv_1d(n_domain[d], n_domain[d], MemoryAllocMode::Stack);
	    ReallocatableArray<int> n_send_1d(n_domain[d], n_domain[d], MemoryAllocMode::Stack);
	    ReallocatableArray<int> n_disp_recv_1d(n_domain[d]+1, n_domain[d]+1, MemoryAllocMode::Stack);
	    ReallocatableArray<int> n_disp_send_1d(n_domain[d]+1, n_domain[d]+1, MemoryAllocMode::Stack);
	    n_disp_recv_1d[0] = n_disp_send_1d[0] = 0;
	    for(auto i=0; i<n_domain[d]; i++){
		n_send_1d[i] = n_recv_1d[i] = 0;
		for(auto j=0; j<radix; j++){
		    n_send_1d[i] += n_send_multi_dim[d][i*radix+j];
		    n_recv_1d[i] += n_recv_multi_dim[d][i*radix+j];
		}
		n_disp_send_1d[i+1] = n_disp_send_1d[i] + n_send_1d[i];
		n_disp_recv_1d[i+1] = n_disp_recv_1d[i] + n_recv_1d[i];
	    }
	    recv_buf.resizeNoInitialize( n_disp_recv_1d[n_domain[d]] );
	    MPI_Alltoallv(send_buf_tmp.getPointer(), n_send_1d.getPointer(), n_disp_send_1d.getPointer(), GetDataType<T>(),
			  recv_buf.getPointer(), n_recv_1d.getPointer(), n_disp_recv_1d.getPointer(), GetDataType<T>(),
			  comm_1d.getCommunicator());
	}
    }

    
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
  
}
