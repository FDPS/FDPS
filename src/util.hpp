
namespace  ParticleSimulator{
    struct Cmpvec{
	F64 F64vec::*loc;
	Cmpvec(F64 F64vec::*_loc) : loc(_loc) {}
	bool operator()(const F64vec &lhs, const F64vec &rhs){
	    return (lhs.*loc < rhs.*loc);
	}
    };
    /*
    inline F64 GetDistanceMinSQF64 (const F64vec & pos0, const F64vec & pos1){
	const F64vec r01 = pos0 - pos1;
	return r01 * r01;
    }
    
    inline F64 GetDistanceMinSQF64 (const F64vec & pos, const F64ort & box){
	const F64vec b_high = this->high_;
	const F64vec b_low = this->low_;
	const F64 dx = (vec.x > b_high.x) ? (vec.x - b_high.x) : ( (vec.x < b_low.x) ? (b_low.x - vec.x) : T(0) );
	const F64 dy = (vec.y > b_high.y) ? (vec.y - b_high.y) : ( (vec.y < b_low.y) ? (b_low.y - vec.y) : T(0) );
	const F64 dz = (vec.z > b_high.z) ? (vec.z - b_high.z) : ( (vec.z < b_low.z) ? (b_low.z - vec.z) : T(0) );
	return dx*dx + dy*dy + dz*dz;	
    }
    
    inline F64 GetDistanceMinSQF64 (const F64ort & box, const F64vec & pos){
	return GetDistanceMinSQF64 (pos, box);
    }
    
    inline F64 GetDistanceMinSQF64 (const F64ort & box0, const F64ort & box1){
            const F64vec a_high = ort.high_;
            const F64vec a_low  = ort.low_;
            const F64vec b_high = this->high_;
            const F64vec b_low  = this->low_;
            const F64 x0 = b_low.x - a_high.x;
            const F64 x1 = a_low.x - b_high.x;
            const F64 y0 = b_low.y - a_high.y;
            const F64 y1 = a_low.y - b_high.y;
            const F64 z0 = b_low.z - a_high.z;
            const F64 z1 = a_low.z - b_high.z;
            F64 dx = (x0 > x1) ? x0 : x1;
            F64 dy = (y0 > y1) ? y0 : y1;
            F64 dz = (z0 > z1) ? z0 : z1;
            dx = (x0*x1 < 0) ? dx : T(0);
            dy = (y0*y1 < 0) ? dy : T(0);
            dz = (z0*z1 < 0) ? dz : T(0);
            return dx*dx + dy*dy + dz*dz;	
    }
    */
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
    
}
