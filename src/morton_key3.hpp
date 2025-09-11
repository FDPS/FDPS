#pragma once

#include <iostream>
#include <iomanip>
#include "vector3.hpp"
#include "ps_defs.hpp"

namespace ParticleSimulator {

/*
  MSB is always 0.
  next 3bits represent octant index of 8 cells with level 1
 */
template <typename T>
class MortonKey {};

template <>
class MortonKey<TagMortonKeyNormal> {
   private:
    enum {
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        kLevMax = 21,
        kLevMaxHi = 21,
        kLevMaxLo = 0,
#elif defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        kLevMaxHi = 21,
        kLevMaxLo = 10,
        kLevMax = 31,
#else
        kLevMaxHi = 21,
        kLevMaxLo = 21,
        kLevMax = 42,
#endif
    };
    F64 half_len_;
    F64vec center_;
    F64 normalized_factor_;
    U64 separateBit(const U64 _s_in) const {
        U64 _s = _s_in;
        _s = (_s | _s << 32) & 0xffff00000000ffff;   // 11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
        _s = (_s | _s << 16) & 0x00ff0000ff0000ff;   // 00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
        _s = (_s | _s << 8) & 0xf00f00f00f00f00f;    // 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
        _s = (_s | _s << 4) & 0x30c30c30c30c30c3;    // 11 00 00 11 00 11 00 00 11
        return (_s | _s << 2) & 0x9249249249249249;  // 1 0 0 1 0 0 1 0 0 1 0 0 1
    }
    U32 separateBit32(const U32 _s_in) const {
        U32 _s = _s_in;
        _s = (_s | _s << 16) & 0xff0000ff;   // 11111111 00000000 00000000 11111111
        _s = (_s | _s << 8) & 0x0f00f00f;    // 1111 0000 0000 1111 0000 0000 1111
        _s = (_s | _s << 4) & 0xc30c30c3;    // 11 00 00 11 00 11 00 00 11
        return (_s | _s << 2) & 0x49249249;  // 1 0 0 1 0 0 1 0 0 1 0 0 1
    }

   public:
    // MortonKey(){};
    //~MortonKey(){};
    // MortonKey(const MortonKey &);
    // MortonKey & operator = (const MortonKey &);

    static inline S32 keylevelmaxhigh() { return S32(kLevMaxHi); }
    static inline S32 keylevelmaxlow() { return S32(kLevMaxLo); }
    void initialize(const F64 half_len, const F64vec& center = 0.0) {
        PS_COMPILE_TIME_MESSAGE("TagMortonKeyNormal");
        half_len_ = half_len;
        center_ = center;
        normalized_factor_ = (1.0 / (half_len * 2.0)) * (1UL << kLevMax);
    }
    KeyT getKey(F64vec pos) const {
        const F64vec cen = center_;
        const F64 hlen = half_len_;
        const F64 nfactor = normalized_factor_;
        U64 nx = (U64)((pos.x - cen.x + hlen) * nfactor);
        U64 ny = (U64)((pos.y - cen.y + hlen) * nfactor);
        U64 nz = (U64)((pos.z - cen.z + hlen) * nfactor);

        // avoid to overflow
        nx = (((nx >> 63) & 0x1) == 1) ? 0 : nx;
        ny = (((ny >> 63) & 0x1) == 1) ? 0 : ny;
        nz = (((nz >> 63) & 0x1) == 1) ? 0 : nz;
        U64 nmax = (1UL << kLevMax) - 1;
        nx = (nx > nmax) ? nmax : nx;
        ny = (ny > nmax) ? nmax : ny;
        nz = (nz > nmax) ? nmax : nz;
        // std::cout<<"A) nx= "<<ny<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
        KeyT ret;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        ret.hi_ = separateBit(nx) << 2 | separateBit(ny) << 1 | separateBit(nz);
#else
        // 96 or 128 bits
        U64 nx_hi = nx >> kLevMaxLo;
        U64 ny_hi = ny >> kLevMaxLo;
        U64 nz_hi = nz >> kLevMaxLo;
        ret.hi_ = (separateBit(nx_hi) << 2 | separateBit(ny_hi) << 1 | separateBit(nz_hi));
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        U32 nx_lo = (U32)(nx & 0x3ff);  // mask 10 bits
        U32 ny_lo = (U32)(ny & 0x3ff);
        U32 nz_lo = (U32)(nz & 0x3ff);
        ret.lo_ = (separateBit32(nx_lo) << 2 | separateBit32(ny_lo) << 1 | separateBit32(nz_lo));
#else
        U64 nx_lo = (U64)(nx & 0x1fffff);  // mask 21bits
        U64 ny_lo = (U64)(ny & 0x1fffff);
        U64 nz_lo = (U64)(nz & 0x1fffff);
        ret.lo_ = (separateBit(nx_lo) << 2 | separateBit(ny_lo) << 1 | separateBit(nz_lo));
#endif
#endif
        return ret;
    }

    template <typename Tkey>
    S32 getCellID(const S32 lev, const Tkey& mkey) const {
        U64 s;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        s = mkey.hi_ >> ((kLevMax - lev) * 3);
#else
        if (lev <= kLevMaxHi) {
            s = mkey.hi_ >> ((kLevMaxHi - lev) * 3);
        } else {
            s = mkey.lo_ >> ((kLevMaxLo - (lev - kLevMaxHi)) * 3);
        }
#endif
        return (S32)(s & 0x7);
    }

    template <typename Tkey>
    S32 getCellIDlow(const S32 lev, const Tkey& mkey) const {
        U64 s;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        s = mkey.hi_ >> ((kLevMax - lev) * 3);
#else
        //	    s = mkey.lo_ >> ( (kLevMaxLo - (lev-kLevMaxHi)) * 3 );
        s = mkey.lo_ >> lev;
#endif
        return (S32)(s & 0x7);
    }

    template <typename Tkey>
    S32 getCellIDhigh(const S32 lev, const Tkey& mkey) const {
        U64 s;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        s = 0;
#else
        //	    s = mkey.hi_ >> ( (kLevMaxHi - lev) * 3 );
        s = mkey.hi_ >> lev;
#endif
        return (S32)(s & 0x7);
    }

    F64ort getCorrespondingTreeCell(const F64ort& box) const {
        const auto low = getKey(box.low_);
        const auto high = getKey(box.high_);
        const auto tmp = high ^ low;
        S32 lev = 0;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        for (S32 i = TREE_LEVEL_LIMIT - 1; i >= 0; i--) {
            if ((((tmp.hi_ >> i * 3)) & 0x7) != 0) {
                lev = TREE_LEVEL_LIMIT - i - 1;
                break;
            }
        }
#else
        if (tmp.hi_ != 0x0) {
            for (S32 i = kLevMaxHi - 1; i >= 0; i--) {
                if ((((tmp.hi_ >> i * 3)) & 0x7) != 0) {
                    lev = kLevMaxHi - i - 1;
                    break;
                }
            }
        } else {
            for (S32 i = kLevMaxLo - 1; i >= 0; i--) {
                if ((((tmp.lo_ >> i * 3)) & 0x7) != 0) {
                    lev = kLevMaxHi + kLevMaxLo - i - 1;
                    break;
                }
            }
        }
#endif
        F64ort ret;
        const F64vec cen = center_;
        const F64 hlen = half_len_;
        if (lev == 0) {
            ret.low_.x = cen.x - hlen;
            ret.low_.y = cen.y - hlen;
            ret.low_.z = cen.z - hlen;
            ret.high_.x = cen.x + hlen;
            ret.high_.y = cen.y + hlen;
            ret.high_.z = cen.z + hlen;
        } else {
            F64vec cen_new = cen;
            F64 hlen_new = hlen;
            for (S32 i = 0; i < lev; i++) {
                S32 id = getCellID(i + 1, low);
                cen_new += hlen_new * SHIFT_CENTER[id];
                hlen_new *= 0.5;
            }
            ret.low_.x = cen_new.x - hlen_new;
            ret.low_.y = cen_new.y - hlen_new;
            ret.low_.z = cen_new.z - hlen_new;
            ret.high_.x = cen_new.x + hlen_new;
            ret.high_.y = cen_new.y + hlen_new;
            ret.high_.z = cen_new.z + hlen_new;
        }
        return ret;
    }
    F64 getCorrespondingFullLength(const S32 lev) const {
        F64 len = half_len_ * 2.0;
        for (S32 i = 0; i < lev; i++) {
            len *= 0.5;
        }
        return len;
    }
    template <typename Tkey>
    F64ort getPosTreeCell(const S32 lev, const Tkey& key) const {
        F64vec cen = center_;
        F64 hlen = half_len_;
        for (S32 i = 0; i < lev; i++) {
            S32 id = getCellID(i + 1, key);
            cen += hlen * SHIFT_CENTER[id];
            hlen *= 0.5;
        }
        F64ort ret;
        ret.low_.x = cen.x - hlen;
        ret.low_.y = cen.y - hlen;
        ret.low_.z = cen.z - hlen;
        ret.high_.x = cen.x + hlen;
        ret.high_.y = cen.y + hlen;
        ret.high_.z = cen.z + hlen;
        return ret;
    }
};

#if 0
    template <>
    class MortonKey<TagMortonKeyMeshBased>{
    private:

        enum{
#if defined(PARTICLE_SIMULATOR_USE_128BIT_KEY)
            kLevMaxHi = 21,
            kLevMaxLo = 21,
            kLevMax = 42,
#elif defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
            kLevMaxHi = 21,
            kLevMaxLo = 10,
            kLevMax = 31,
#else
            kLevMax = 21,
#endif
        };
        F64vec half_len_;
        F64vec center_;
        F64vec pos_ref_; // used only in getKey
        F64vec wid_ref_; // used only in getKey (PM cell width)
        F64vec wid_btm_; // used only in getKey (minimum tree cell width)
        S32vec idx_ref_;
        S32    lev_ref_;
        U64 separateBit(const U64 _s_in) const {
            U64 _s = _s_in;
            _s = (_s | _s<<32) & 0xffff00000000ffff; //11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
            _s = (_s | _s<<16) & 0x00ff0000ff0000ff; //00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
            _s = (_s | _s<<8) & 0xf00f00f00f00f00f;  //1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
            _s = (_s | _s<<4) & 0x30c30c30c30c30c3;  //11 00 00 11 00 11 00 00 11
            return (_s | _s<<2) & 0x9249249249249249;  //1 0 0 1 0 0 1 0 0 1 0 0 1
        }
        U32 separateBit32(const U32 _s_in) const {
            U32 _s = _s_in;
            _s = (_s | _s<<16)  & 0xff0000ff;  //11111111 00000000 00000000 11111111
            _s = (_s | _s<<8)   & 0x0f00f00f;  //1111 0000 0000 1111 0000 0000 1111
            _s = (_s | _s<<4)   & 0xc30c30c3;  //11 00 00 11 00 11 00 00 11
            return (_s | _s<<2) & 0x49249249;  //1 0 0 1 0 0 1 0 0 1 0 0 1
        }
        U64 shrinkBit(const U64 _s_in) const {
            U64 _s = _s_in;
            _s = (_s | _s>>2)  & 0x30c30c30c30c30c3ULL; // Repeat of pattern 00 00 11
            _s = (_s | _s>>4)  & 0xf00f00f00f00f00fULL; // Repeat of pattern 0000 0000 1111
            _s = (_s | _s>>8)  & 0x00ff0000ff0000ffULL; // Repeat of pattern 00000000 00000000 11111111
            _s = (_s | _s>>16) & 0xffff00000000ffffULL; // Repeat of pattern 0^{16} 0^{16} 1^{16},
                                                        // where x^{n} reprensets that x continues n times.
            _s = (_s | _s>>32) & 0x00000000ffffffffULL; // Repaat of pattern 0^{32} 0^{32} 1^{32}
            return _s;
        }
        U32 shrinkBit32(const U32 _s_in) const {
            U32 _s = _s_in;
            _s = (_s | _s>>2)  & 0xc30c30c3; // Repeat of pattern 00 00 11
            _s = (_s | _s>>4)  & 0x0f00f00f; // Repeat of pattern 0000 0000 1111
            _s = (_s | _s>>8)  & 0xff0000ff; // Repeat of pattern 00000000 00000000 11111111
            _s = (_s | _s>>16) & 0x0000ffff; // Repeat of pattern 0^{16} 0^{16} 1^{16},
                                             // where x^{n} reprensets that x continues n times.
            return _s;
        }
    public:
       
        void initialize(const F64 half_len, const F64vec & center) {
	    PS_COMPILE_TIME_MESSAGE("TagMortonKeyMeshBased");
            /* Do nothing */
        }
        void initialize(const F64vec half_len,
                        const F64vec & center,
                        const F64vec & pos_ref, // pos_unit_cell.low_
                        const S64vec & idx_ref, // idx_unit_cell_
                        const S32 lev_ref, // lev_pm_cell_
                        const F64vec & wid_ref){ // width_pm_cell_
	    PS_COMPILE_TIME_MESSAGE("TagMortonKeyMeshBased");
            half_len_ = half_len;
            center_ = center;
            pos_ref_ = pos_ref;
            idx_ref_ = idx_ref;
            wid_ref_ = wid_ref;
            lev_ref_ = lev_ref;
            wid_btm_ = wid_ref / (double)(1UL<<(kLevMax - lev_ref)); // minimum tree cell width
        }
        KeyT getKey(F64vec pos) const {
            // Compute cell index
            // idx_pm is an index of a PM cell, whose width is wid_ref_, to witch pos belongs. This index of (0, 0, 0) corresponds to the PM cell at the lower limit of the unit cell (computational box).
	    //std::cout<<"A) pos= "<<pos<<std::endl;
            S64vec idx_pm = GetCellIDMeasuredInUnitCell(pos_ref_, wid_ref_, pos);
            //std::cout<<"idx_pm= "<<idx_pm<<std::endl;
            //std::cout<<"pos_ref_= "<<pos_ref_<<std::endl;
            //std::cout<<"wid_ref_= "<<wid_ref_<<std::endl;
            // idx_btm is an index of a minimum tree cell, whose width is wid_btm_, to witch pos belongs. This index of (0, 0, 0) corresponds to the tree cell at the lower limit of the cell whose PM cell index is idx_pm
            S64vec idx_btm;
            if (kLevMax > lev_ref_) {
                // Compute sub-cell index in cell idx_pm.
                F64vec pos_cell;
                pos_cell.x = pos_ref_.x + wid_ref_.x * idx_pm.x;
                pos_cell.y = pos_ref_.y + wid_ref_.y * idx_pm.y;
                pos_cell.z = pos_ref_.z + wid_ref_.z * idx_pm.z;
                // [Notes]
                //    In principal, we can choose pos_cell and wid_btm
                //    different from those we are currently adopting
                //    because it is guaranteed that the particle of interest
                //    belongs to the cell to which it should belong.
                //    The only condition imposed on sub-cell indices is that
                //    the order of the particles belonging to that cell is
                //    preserved.
		//std::cout<<"B) pos= "<<pos<<std::endl;
                idx_btm = GetCellIDMeasuredInUnitCell(pos_cell, wid_btm_, pos);
                //std::cout<<"idx_btm= "<<idx_btm<<std::endl;
                //std::cout<<"pos_cell= "<<pos_cell<<std::endl;
                //std::cout<<"wid_btm_= "<<wid_btm_<<std::endl;
                const S64 maxval = (1L<<(kLevMax - lev_ref_)) - 1;
                if (idx_btm.x < 0) idx_btm.x = 0;
                if (idx_btm.x > maxval) idx_btm.x = maxval;
                if (idx_btm.y < 0) idx_btm.y = 0;
                if (idx_btm.y > maxval) idx_btm.y = maxval;
                if (idx_btm.z < 0) idx_btm.z = 0;
                if (idx_btm.z > maxval) idx_btm.z = maxval;
            } else {
                idx_btm.x = idx_btm.y = idx_btm.z = 0;
            }
            U64 nx = (U64)((idx_pm.x + idx_ref_.x) * (1L<<(kLevMax - lev_ref_)) + idx_btm.x);
            U64 ny = (U64)((idx_pm.y + idx_ref_.y) * (1L<<(kLevMax - lev_ref_)) + idx_btm.y);
            U64 nz = (U64)((idx_pm.z + idx_ref_.z) * (1L<<(kLevMax - lev_ref_)) + idx_btm.z);
	    /*
            std::cout<<"idx_ref_= "<<idx_ref_<<std::endl;
            std::cout<<"kLevMax= "<<kLevMax<<" lev_ref_= "<<lev_ref_<<std::endl;
            std::cout<<"(idx_pm.x + idx_ref_.x) * (1L<<(kLevMax - lev_ref_))= "<<(idx_pm.x + idx_ref_.x) * (1L<<(kLevMax - lev_ref_))<<std::endl;
            std::cout<<"(idx_pm.y + idx_ref_.y) * (1L<<(kLevMax - lev_ref_))= "<<(idx_pm.y + idx_ref_.y) * (1L<<(kLevMax - lev_ref_))<<std::endl;
            std::cout<<"(idx_pm.x + idx_ref_.x) * (1L<<(kLevMax - lev_ref_))+idx_btm.x= "<<(idx_pm.x + idx_ref_.x) * (1L<<(kLevMax - lev_ref_))+idx_btm.x<<std::endl;
            std::cout<<"(idx_pm.y + idx_ref_.y) * (1L<<(kLevMax - lev_ref_))+idx_btm.y= "<<(idx_pm.y + idx_ref_.y) * (1L<<(kLevMax - lev_ref_))+idx_btm.y<<std::endl;
            std::cout<<"idx_btm.x= "<<idx_btm.x<<std::endl;
            std::cout<<"idx_btm.y= "<<idx_btm.y<<std::endl;
            std::cout<<"nx= "<<ny<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
	    */
            // avoid to overflow
            nx = (((nx>>63)&0x1)==1) ? 0 : nx;
            ny = (((ny>>63)&0x1)==1) ? 0 : ny;
            nz = (((nz>>63)&0x1)==1) ? 0 : nz;
            U64 nmax = (1UL<<kLevMax)-1;
            nx = (nx > nmax) ? nmax : nx;
            ny = (ny > nmax) ? nmax : ny;
            nz = (nz > nmax) ? nmax : nz;
            //std::cout<<"nx= "<<ny<<" ny= "<<ny<<" nz= "<<nz<<std::endl;            
            KeyT ret;
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY) || defined(PARTICLE_SIMULATOR_USE_128BIT_KEY)
            // for 128bit and 96bit
            U64 nx_hi = nx>>kLevMaxLo;
            U64 ny_hi = ny>>kLevMaxLo;
            U64 nz_hi = nz>>kLevMaxLo;
            ret.hi_ = (separateBit(nx_hi)<<2
                       | separateBit(ny_hi)<<1
                       | separateBit(nz_hi));
#if defined(PARTICLE_SIMULATOR_USE_128BIT_KEY)
            // for 128bit
            U64 nx_lo = (U64)(nx & 0x1fffff); // mask 21bits
            U64 ny_lo = (U64)(ny & 0x1fffff);
            U64 nz_lo = (U64)(nz & 0x1fffff);
            ret.lo_ = ( separateBit(nx_lo)<<2
                        | separateBit(ny_lo)<<1
                        | separateBit(nz_lo) );
#else
            // for 96bit
            U32 nx_lo = (U32)(nx & 0x3ff); // mask 10 bits
            U32 ny_lo = (U32)(ny & 0x3ff);
            U32 nz_lo = (U32)(nz & 0x3ff);
            ret.lo_ = ( separateBit32(nx_lo)<<2
                        | separateBit32(ny_lo)<<1
                        | separateBit32(nz_lo) );
#endif  // end of (USE_96BIT_KEY) || defined (USE_128BIT_KEY)
#else
            // for 64bit
            ret.hi_ = separateBit(nx)<<2 | separateBit(ny)<<1 | separateBit(nz);
#endif            
            /*
#if defined (USE_96BIT_KEY) || defined (USE_128BIT_KEY)
            // for 128bit and 96bit
            U64 nx_hi = nx>>kLevMaxLo;
            U64 ny_hi = ny>>kLevMaxLo;
            U64 nz_hi = nz>>kLevMaxLo;
            ret.hi_ = (separateBit(nx_hi)<<2
                       | separateBit(ny_hi)<<1
                       | separateBit(nz_hi));
#if defined (USE_128BIT_KEY)
            // for 128bit
            U64 nx_lo = (U64)(nx & 0x1fffff); // mask 21bits
            U64 ny_lo = (U64)(ny & 0x1fffff);
            U64 nz_lo = (U64)(nz & 0x1fffff);
            ret.lo_ = ( separateBit(nx_lo)<<2
                        | separateBit(ny_lo)<<1
                        | separateBit(nz_lo) );
#else
            // for 96bit
            U32 nx_lo = (U32)(nx & 0x3ff); // mask 10 bits
            U32 ny_lo = (U32)(ny & 0x3ff);
            U32 nz_lo = (U32)(nz & 0x3ff);
            ret.lo_ = ( separateBit32(nx_lo)<<2
                        | separateBit32(ny_lo)<<1
                        | separateBit32(nz_lo) );
#endif // end of (USE_96BIT_KEY) || defined (USE_128BIT_KEY)
#else
            // for 64bit
            ret.hi_ = separateBit(nx)<<2 | separateBit(ny)<<1 | separateBit(nz);
#endif
            */
            return ret;
        }

        template<typename Tkey>
        S32 getCellID(const S32 lev, const Tkey & mkey) const {
            U64 s;
#if defined(PARTICLE_SIMULATOR_USE_128BIT_KEY) || defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
            if(lev <= kLevMaxHi){
                s = mkey.hi_ >> ( (kLevMaxHi - lev) * 3 );
            }
            else{
                s = mkey.lo_ >> ( (kLevMaxLo - (lev-kLevMaxHi)) * 3 );
            }
#else
            s = mkey.hi_ >> ( (kLevMax - lev) * 3 );
#endif            
            /*
#if defined (USE_128BIT_KEY) || defined (USE_96BIT_KEY)
            if(lev <= kLevMaxHi){
                s = mkey.hi_ >> ( (kLevMaxHi - lev) * 3 );
            }
            else{
                s = mkey.lo_ >> ( (kLevMaxLo - (lev-kLevMaxHi)) * 3 );
            }
#else
            s = mkey.hi_ >> ( (kLevMax - lev) * 3 );
#endif
            */
            return (S32)(s & 0x7);
        }
        
        F64ort getCorrespondingTreeCell(const F64ort & box) const {
            const auto low  = getKey(box.low_);
            const auto high = getKey(box.high_);
            const auto tmp = high^low;
            S32 lev = 0;
#if defined(PARTICLE_SIMULATOR_USE_128BIT_KEY) || defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
            if(tmp.hi_ != 0x0){
                for(S32 i=kLevMaxHi-1; i>=0; i--){
                    if( (((tmp.hi_ >> i*3)) & 0x7) != 0){
                        lev = kLevMaxHi-i-1;
                        break;
                    }
                }
            }
            else{
                for(S32 i=kLevMaxLo-1; i>=0; i--){
                    if( (((tmp.lo_ >> i*3)) & 0x7) != 0){
                        lev = kLevMaxHi+kLevMaxLo-i-1;
                        break;
                    }
                }
            }
#else
            for(S32 i=TREE_LEVEL_LIMIT-1; i>=0; i--){
                if( (((tmp.hi_ >> i*3)) & 0x7) != 0){
                    lev = TREE_LEVEL_LIMIT-i-1;
                    break;
                }
            }
#endif
            F64ort ret;
            const F64vec cen = center_;
            const F64vec hlen = half_len_;
            if(lev==0){
                ret.low_.x  = cen.x - hlen.x;
                ret.low_.y  = cen.y - hlen.y;
                ret.low_.z  = cen.z - hlen.z;
                ret.high_.x = cen.x + hlen.x;
                ret.high_.y = cen.y + hlen.y;
                ret.high_.z = cen.z + hlen.z;
            }
            else{
                F64vec cen_new = cen;
                F64vec hlen_new = hlen;
                for(S32 i=0; i<lev; i++){
                    S32 id = getCellID(i+1, low);
                    cen_new +=  hlen_new*SHIFT_CENTER[id];
                    hlen_new *= 0.5;
                }
                ret.low_.x  = cen_new.x - hlen_new.x;
                ret.low_.y  = cen_new.y - hlen_new.y;
                ret.low_.z  = cen_new.z - hlen_new.z;
                ret.high_.x = cen_new.x + hlen_new.x;
                ret.high_.y = cen_new.y + hlen_new.y;
                ret.high_.z = cen_new.z + hlen_new.z;
            }
            return ret;
        }
        F64 getCorrespondingFullLength(const S32 lev) const {
            F64vec len = half_len_*2.0;
            for(S32 i=0; i<lev; i++){
                len *= 0.5;
            }
            return std::max(std::max(len.x, len.y), len.z);
        }
        template<typename Tkey>
        F64ort getPosTreeCell(const S32 lev, const Tkey & key) const {
            F64vec cen = center_;
            F64vec hlen = half_len_;
            for(S32 i=0; i<lev; i++){
                S32 id = getCellID(i+1, key);
                cen +=  hlen*SHIFT_CENTER[id];
                hlen *= 0.5;
            }
            F64ort ret;
            ret.low_.x  = cen.x - hlen.x;
            ret.low_.y  = cen.y - hlen.y;
            ret.low_.z  = cen.z - hlen.z;
            ret.high_.x = cen.x + hlen.x;
            ret.high_.y = cen.y + hlen.y;
            ret.high_.z = cen.z + hlen.z;
            return ret;
        }
        // necessary functions for the PMMM feature
        S32vec getPMCellID(const F64vec pos) const {
            const KeyT mkey = getKey(pos);
            return getPMCellID(mkey);
        }
        template <typename TKey>
        S32vec getPMCellID(const TKey mkey) const {
	    PS_COMPILE_TIME_MESSAGE("getPMCellID(Tkey)@TagMortonKeyMeshBased");
            S32 ix,iy,iz;
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY) || defined(PARTICLE_SIMULATOR_USE_128BIT_KEY)
            if (lev_ref_ <= kLevMaxHi) {
                const U64 x = (mkey.hi_ & 0x4924924924924924ULL)>>2; // 0 100^{21}
                const U64 y = (mkey.hi_ & 0x2492492492492492ULL)>>1; // 0 010^{21}
                const U64 z =  mkey.hi_ & 0x1249249249249249ULL;     // 0 001^{21}
                ix = (S32)(shrinkBit(x)>>(kLevMaxHi - lev_ref_)) - idx_ref_.x;
                iy = (S32)(shrinkBit(y)>>(kLevMaxHi - lev_ref_)) - idx_ref_.y;
                iz = (S32)(shrinkBit(z)>>(kLevMaxHi - lev_ref_)) - idx_ref_.z;
            } else {
                const U64 x_hi = (mkey.hi_ & 0x4924924924924924ULL)>>2;
                const U64 y_hi = (mkey.hi_ & 0x2492492492492492ULL)>>1;
                const U64 z_hi =  mkey.hi_ & 0x1249249249249249ULL;
                U64 ix_hi = shrinkBit(x_hi)<<kLevMaxLo; // operation << is to make room
                U64 iy_hi = shrinkBit(y_hi)<<kLevMaxLo;
                U64 iz_hi = shrinkBit(z_hi)<<kLevMaxLo;
#if defined(PARTICLE_SIMULATOR_USE_128BIT_KEY)
                const U64 x_lo = (mkey.lo_ & 0x4924924924924924ULL)>>2;
                const U64 y_lo = (mkey.lo_ & 0x2492492492492492ULL)>>1;
                const U64 z_lo =  mkey.lo_ & 0x1249249249249249ULL;
                U64 ix_lo = (shrinkBit(x_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
                U64 iy_lo = (shrinkBit(y_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
                U64 iz_lo = (shrinkBit(z_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
#else   // defined(USE_128BIT_KEY)
                const U32 x_lo = (mkey.lo_ & 0x24924924)>>2; // 00 100^{10}
                const U32 y_lo = (mkey.lo_ & 0x12492492)>>1; // 00 010^{10}
                const U32 z_lo =  mkey.lo_ & 0x09249249;     // 00 001^{10}
                U64 ix_lo = (shrinkBit32(x_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
                U64 iy_lo = (shrinkBit32(y_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
                U64 iz_lo = (shrinkBit32(z_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
#endif  // defined(USE_128BIT_KEY)
                ix = (S32)(ix_hi | ix_lo) - idx_ref_.x;
                iy = (S32)(iy_hi | iy_lo) - idx_ref_.y;
                iz = (S32)(iz_hi | iz_lo) - idx_ref_.z;
            }
#else   // defined(USE_96BIT_KEY) || defined(USE_128BIT_KEY)
            const U64 x = (mkey.hi_ & 0x4924924924924924ULL)>>2;
            const U64 y = (mkey.hi_ & 0x2492492492492492ULL)>>1;
            const U64 z =  mkey.hi_ & 0x1249249249249249ULL;
            ix = (S32)(shrinkBit(x)>>(kLevMax - lev_ref_)) - idx_ref_.x;
            iy = (S32)(shrinkBit(y)>>(kLevMax - lev_ref_)) - idx_ref_.y;
            iz = (S32)(shrinkBit(z)>>(kLevMax - lev_ref_)) - idx_ref_.z;
#endif  // defined(USE_96BIT_KEY) || defined(USE_128BIT_KEY)  
            return S32vec(ix, iy, iz);
        }
    };
#endif

template <>
class MortonKey<TagMortonKeyMeshBased2> {
   private:
    enum {
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        kLevMax = 21,
        kLevMaxHi = 21,
        kLevMaxLo = 0,
#elif defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        kLevMaxHi = 21,
        kLevMaxLo = 10,
        kLevMax = 31,
#else
        kLevMaxHi = 21,
        kLevMaxLo = 21,
        kLevMax = 42,
#endif
    };
    F64vec half_len_;
    F64vec center_;
    F64vec pos_ref_;  // used only in getKey
    F64vec wid_ref_;  // used only in getKey (PM cell width)
    F64vec wid_btm_;  // used only in getKey (minimum tree cell width)
    S32vec idx_ref_;
    S32 lev_ref_;
    F64vec normalized_factor_;
    U64 separateBit(const U64 _s_in) const {
        U64 _s = _s_in;
        _s = (_s | _s << 32) & 0xffff00000000ffff;   // 11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
        _s = (_s | _s << 16) & 0x00ff0000ff0000ff;   // 00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
        _s = (_s | _s << 8) & 0xf00f00f00f00f00f;    // 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
        _s = (_s | _s << 4) & 0x30c30c30c30c30c3;    // 11 00 00 11 00 11 00 00 11
        return (_s | _s << 2) & 0x9249249249249249;  // 1 0 0 1 0 0 1 0 0 1 0 0 1
    }
    U32 separateBit32(const U32 _s_in) const {
        U32 _s = _s_in;
        _s = (_s | _s << 16) & 0xff0000ff;   // 11111111 00000000 00000000 11111111
        _s = (_s | _s << 8) & 0x0f00f00f;    // 1111 0000 0000 1111 0000 0000 1111
        _s = (_s | _s << 4) & 0xc30c30c3;    // 11 00 00 11 00 11 00 00 11
        return (_s | _s << 2) & 0x49249249;  // 1 0 0 1 0 0 1 0 0 1 0 0 1
    }
    U64 shrinkBit(const U64 _s_in) const {
        U64 _s = _s_in;
        _s = (_s | _s >> 2) & 0x30c30c30c30c30c3ULL;   // Repeat of pattern 00 00 11
        _s = (_s | _s >> 4) & 0xf00f00f00f00f00fULL;   // Repeat of pattern 0000 0000 1111
        _s = (_s | _s >> 8) & 0x00ff0000ff0000ffULL;   // Repeat of pattern 00000000 00000000 11111111
        _s = (_s | _s >> 16) & 0xffff00000000ffffULL;  // Repeat of pattern 0^{16} 0^{16} 1^{16},
                                                       // where x^{n} reprensets that x continues n times.
        _s = (_s | _s >> 32) & 0x00000000ffffffffULL;  // Repaat of pattern 0^{32} 0^{32} 1^{32}
        return _s;
    }
    U32 shrinkBit32(const U32 _s_in) const {
        U32 _s = _s_in;
        _s = (_s | _s >> 2) & 0xc30c30c3;   // Repeat of pattern 00 00 11
        _s = (_s | _s >> 4) & 0x0f00f00f;   // Repeat of pattern 0000 0000 1111
        _s = (_s | _s >> 8) & 0xff0000ff;   // Repeat of pattern 00000000 00000000 11111111
        _s = (_s | _s >> 16) & 0x0000ffff;  // Repeat of pattern 0^{16} 0^{16} 1^{16},
                                            // where x^{n} reprensets that x continues n times.
        return _s;
    }

   public:
    void initialize(const F64vec half_len, const F64vec& center) {
        PS_COMPILE_TIME_MESSAGE("TagMortonKeyMeshBased2");
        half_len_ = half_len;
        center_ = center;
        normalized_factor_.x = (1.0 / (half_len.x * 2.0)) * (1UL << kLevMax);
        normalized_factor_.y = (1.0 / (half_len.y * 2.0)) * (1UL << kLevMax);
        normalized_factor_.z = (1.0 / (half_len.z * 2.0)) * (1UL << kLevMax);
    }

    void initialize(const F64vec half_len, const F64vec& center,
                    const F64vec& pos_ref,    // pos_unit_cell.low_
                    const S32vec& idx_ref,    // idx_unit_cell_
                    const S32 lev_ref,        // lev_pm_cell_
                    const F64vec& wid_ref) {  // width_pm_cell_
        PS_COMPILE_TIME_MESSAGE("TagMortonKeyMeshBased2");
        half_len_ = half_len;
        center_ = center;
        pos_ref_ = pos_ref;
        idx_ref_ = idx_ref;
        wid_ref_ = wid_ref;
        lev_ref_ = lev_ref;
        wid_btm_ = wid_ref / (double)(1UL << (kLevMax - lev_ref));  // minimum tree cell width
        normalized_factor_.x = (1.0 / (half_len.x * 2.0)) * (1UL << kLevMax);
        normalized_factor_.y = (1.0 / (half_len.y * 2.0)) * (1UL << kLevMax);
        normalized_factor_.z = (1.0 / (half_len.z * 2.0)) * (1UL << kLevMax);
    }

    KeyT getKey(F64vec pos) const {
#if 0
// original version
        const F64vec cen = center_;
        const F64vec hlen = half_len_;
        const F64vec nfactor = normalized_factor_;
        // std::cout<<"pos= "<<pos<<std::endl;
        // std::cout<<"cen= "<<cen<<std::endl;
        // std::cout<<"hlen= "<<hlen<<std::endl;
        // std::cout<<"cen-hlen= "<<cen-hlen<<std::endl;
        // std::cout<<"cen+hlen= "<<cen+hlen<<std::endl;
        // std::cout<<"nfactor= "<<nfactor<<std::endl;
        U64 nx = (U64)((pos.x - cen.x + hlen.x) * nfactor.x);
        U64 ny = (U64)((pos.y - cen.y + hlen.y) * nfactor.y);
        U64 nz = (U64)((pos.z - cen.z + hlen.z) * nfactor.z);

        std::cout<<"pos= "<<pos<<" nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
        std::cout<<"cen= "<<cen<<" hlen= "<<hlen<<" nfactor= "<<nfactor<<std::endl;
        std::cout<<std::setprecision(15)<<"(pos.x - cen.x + hlen.x)= "<<(pos.x - cen.x + hlen.x)<<std::endl;

        // avoid to overflow
        nx = (((nx >> 63) & 0x1) == 1) ? 0 : nx;
        ny = (((ny >> 63) & 0x1) == 1) ? 0 : ny;
        nz = (((nz >> 63) & 0x1) == 1) ? 0 : nz;
        U64 nmax = (1UL << kLevMax) - 1;
        nx = (nx > nmax) ? nmax : nx;
        ny = (ny > nmax) ? nmax : ny;
        nz = (nz > nmax) ? nmax : nz;

        // std::cout<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<" nmax= "<<nmax<<std::endl;

        KeyT ret;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        ret.hi_ = separateBit(nx) << 2 | separateBit(ny) << 1 | separateBit(nz);
#else
        U64 nx_hi = nx >> kLevMaxLo;
        U64 ny_hi = ny >> kLevMaxLo;
        U64 nz_hi = nz >> kLevMaxLo;
        ret.hi_ = (separateBit(nx_hi) << 2 | separateBit(ny_hi) << 1 | separateBit(nz_hi));
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        U32 nx_lo = (U32)(nx & 0x3ff);  // mask 10 bits
        U32 ny_lo = (U32)(ny & 0x3ff);
        U32 nz_lo = (U32)(nz & 0x3ff);
        ret.lo_ = (separateBit32(nx_lo) << 2 | separateBit32(ny_lo) << 1 | separateBit32(nz_lo));
#else   // PARTICLE_SIMULATOR_USE_96BIT_KEY
        U64 nx_lo = (U64)(nx & 0x1fffff);  // mask 21bits
        U64 ny_lo = (U64)(ny & 0x1fffff);
        U64 nz_lo = (U64)(nz & 0x1fffff);
        ret.lo_ = (separateBit(nx_lo) << 2 | separateBit(ny_lo) << 1 | separateBit(nz_lo));
#endif  // PARTICLE_SIMULATOR_USE_96BIT_KEY
#endif  // PARTICLE_SIMULATOR_USE_64BIT_KEY
        return ret;
#else
        // original version
        // Compute cell index
        // idx_pm is an index of a PM cell, whose width is wid_ref_, to witch pos belongs. 
        // This index of (0, 0, 0) corresponds to the PM cell at the
        // lower limit of the unit cell (computational box).
        S32vec idx_pm = GetCellIDMeasuredInUnitCell(pos_ref_, wid_ref_, pos);
        // idx_btm is an index of a minimum tree cell, whose width is wid_btm_, to witch pos belongs. 
        // This index of (0, 0, 0) corresponds to the tree
        // cell at the lower limit of the cell whose PM cell index is idx_pm
        //S32vec idx_btm;
        S64vec idx_btm;

        //std::cout<<"*******"<<std::endl;
        //std::cout<<"pos= "<<pos<<std::endl;

        if (kLevMax > lev_ref_) {
            // Compute sub-cell index in cell idx_pm.
            F64vec pos_cell;
            pos_cell.x = pos_ref_.x + wid_ref_.x * idx_pm.x;
            pos_cell.y = pos_ref_.y + wid_ref_.y * idx_pm.y;
            pos_cell.z = pos_ref_.z + wid_ref_.z * idx_pm.z;
            // [Notes]
            //    In principal, we can choose pos_cell and wid_btm
            //    different from those we are currently adopting
            //    because it is guaranteed that the particle of interest
            //    belongs to the cell to which it should belong.
            //    The only condition imposed on sub-cell indices is that
            //    the order of the particles belonging to that cell is
            //    preserved.
            idx_btm = GetCellIDMeasuredInUnitCell(pos_cell, wid_btm_, pos);
            //const S32 maxval = (1 << (kLevMax - lev_ref_)) - 1;
            const S64 maxval = (1UL << (kLevMax - lev_ref_)) - 1;
            if (idx_btm.x < 0) idx_btm.x = 0;
            if (idx_btm.x > maxval) idx_btm.x = maxval;
            if (idx_btm.y < 0) idx_btm.y = 0;
            if (idx_btm.y > maxval) idx_btm.y = maxval;
            if (idx_btm.z < 0) idx_btm.z = 0;
            if (idx_btm.z > maxval) idx_btm.z = maxval;
            //std::cout<<"pos_ref_= "<<pos_ref_<<std::endl;
            //std::cout<<"pos_cell= "<<pos_cell<<std::endl;
            //std::cout<<"wid_ref_= "<<wid_ref_<<std::endl;
            //std::cout<<"wid_btm_= "<<wid_btm_<<std::endl;
        } else {
            idx_btm.x = idx_btm.y = idx_btm.z = 0;
        }
        //U64 nx = (U64)((idx_pm.x + idx_ref_.x) * (1 << (kLevMax - lev_ref_)) + idx_btm.x);
        U64 ONE = 1;
        U64 nx = (U64)((idx_pm.x + idx_ref_.x) * (ONE << (kLevMax - lev_ref_)) + idx_btm.x);
        U64 ny = (U64)((idx_pm.y + idx_ref_.y) * (ONE << (kLevMax - lev_ref_)) + idx_btm.y);
        U64 nz = (U64)((idx_pm.z + idx_ref_.z) * (ONE << (kLevMax - lev_ref_)) + idx_btm.z);


        //std::cout << "idx_pm= " << idx_pm << " idx_btm= " << idx_btm << " idx_ref_= " << idx_ref_ << std::endl;
        //std::cout<<"kLevMax= "<<kLevMax<<" lev_ref_= "<<lev_ref_<<" (ONE << (kLevMax - lev_ref_))= "<<(ONE << (kLevMax - lev_ref_))<<std::endl;
        //std::cout<<"((idx_pm.x + idx_ref_.x) * (1 << (kLevMax - lev_ref_)) + idx_btm.x)= "<<((idx_pm.x + idx_ref_.x) * (1 << (kLevMax - lev_ref_)) + idx_btm.x)<<std::endl;
        //std::cout << std::oct;
        //std::cout << "nx= " << nx << " ny= " << ny << " nz= " << nz << std::endl;
        //std::cout << std::dec;

        // avoid to overflow
        nx = (((nx >> 63) & 0x1) == 1) ? 0 : nx;
        ny = (((ny >> 63) & 0x1) == 1) ? 0 : ny;
        nz = (((nz >> 63) & 0x1) == 1) ? 0 : nz;

        //std::cout << "B) nx= " << nx << " ny= " << ny << " nz= " << nz << std::endl;

        U64 nmax = (1UL << kLevMax) - 1;
        nx = (nx > nmax) ? nmax : nx;
        ny = (ny > nmax) ? nmax : ny;
        nz = (nz > nmax) ? nmax : nz;

        //std::cout << "C) nx= " << nx << " ny= " << ny << " nz= " << nz << std::endl;

        KeyT ret;
        /*
#if defined (PARTICLE_SIMULATOR_USE_96BIT_KEY) || defined (PARTICLE_SIMULATOR_USE_128BIT_KEY)
        // for 128bit and 96bit
        U64 nx_hi = nx>>kLevMaxLo;
        U64 ny_hi = ny>>kLevMaxLo;
        U64 nz_hi = nz>>kLevMaxLo;
        ret.hi_ = (separateBit(nx_hi)<<2
                   | separateBit(ny_hi)<<1
                   | separateBit(nz_hi));
#if defined (PARTICLE_SIMULATOR_USE_128BIT_KEY)
        // for 128bit
        U64 nx_lo = (U64)(nx & 0x1fffff); // mask 21bits
        U64 ny_lo = (U64)(ny & 0x1fffff);
        U64 nz_lo = (U64)(nz & 0x1fffff);
        ret.lo_ = ( separateBit(nx_lo)<<2
                    | separateBit(ny_lo)<<1
                    | separateBit(nz_lo) );
#else
        // for 96bit
        U32 nx_lo = (U32)(nx & 0x3ff); // mask 10 bits
        U32 ny_lo = (U32)(ny & 0x3ff);
        U32 nz_lo = (U32)(nz & 0x3ff);
        ret.lo_ = ( separateBit32(nx_lo)<<2
                    | separateBit32(ny_lo)<<1
                    | separateBit32(nz_lo) );
#endif // end of (USE_96BIT_KEY) || defined (USE_128BIT_KEY)
#else
        // for 64bit
        ret.hi_ = separateBit(nx)<<2 | separateBit(ny)<<1 | separateBit(nz);
#endif
        */
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        ret.hi_ = separateBit(nx) << 2 | separateBit(ny) << 1 | separateBit(nz);
        std::cout<<std::oct;
        std::cout<<"nx= "<<nx<<std::endl;
        std::cout<<"ret.hi_= "<<ret.hi_<<std::endl;
        std::cout<<std::dec;                
#else
        U64 nx_hi = nx >> kLevMaxLo;
        U64 ny_hi = ny >> kLevMaxLo;
        U64 nz_hi = nz >> kLevMaxLo;
        ret.hi_ = (separateBit(nx_hi) << 2 | separateBit(ny_hi) << 1 | separateBit(nz_hi));
        //std::cout<<std::oct;
        //std::cout<<"nx= "<<nx<<" nx_hi= "<<nx_hi<<std::endl;
        //std::cout<<"ret.hi_= "<<ret.hi_<<std::endl;
        //std::cout<<std::dec;
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        U32 nx_lo = (U32)(nx & 0x3ff);  // mask 10 bits
        U32 ny_lo = (U32)(ny & 0x3ff);
        U32 nz_lo = (U32)(nz & 0x3ff);
        ret.lo_ = (separateBit32(nx_lo) << 2 | separateBit32(ny_lo) << 1 | separateBit32(nz_lo));
#else   // PARTICLE_SIMULATOR_USE_96BIT_KEY
        // 128bit
        U64 nx_lo = (U64)(nx & 0x1fffff);  // mask 21bits
        U64 ny_lo = (U64)(ny & 0x1fffff);
        U64 nz_lo = (U64)(nz & 0x1fffff);
        ret.lo_ = (separateBit(nx_lo) << 2 | separateBit(ny_lo) << 1 | separateBit(nz_lo));
#endif  // PARTICLE_SIMULATOR_USE_96BIT_KEY
#endif  // PARTICLE_SIMULATOR_USE_64BIT_KEY
        return ret;
#endif  // #if 1
    }

    template <typename Tkey>
    S32 getCellID(const S32 lev, const Tkey& mkey) const {
        U64 s;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        s = mkey.hi_ >> ((kLevMax - lev) * 3);
#else
        if (lev <= kLevMaxHi) {
            s = mkey.hi_ >> ((kLevMaxHi - lev) * 3);
        } else {
            s = mkey.lo_ >> ((kLevMaxLo - (lev - kLevMaxHi)) * 3);
        }
#endif
        return (S32)(s & 0x7);
    }
    template <typename Tkey>
    S32 getCellIDlow(const S32 lev, const Tkey& mkey) const {
        U64 s;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        s = mkey.hi_ >> ((kLevMax - lev) * 3);
#else
        //	    s = mkey.lo_ >> ( (kLevMaxLo - (lev-kLevMaxHi)) * 3 );
        s = mkey.lo_ >> lev;
#endif
        return (S32)(s & 0x7);
    }
    F64ort getCorrespondingTreeCell(const F64ort& box) const {
        const auto low = getKey(box.low_);
        const auto high = getKey(box.high_);
        const auto tmp = high ^ low;
        S32 lev = 0;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        for (S32 i = TREE_LEVEL_LIMIT - 1; i >= 0; i--) {
            if ((((tmp.hi_ >> i * 3)) & 0x7) != 0) {
                lev = TREE_LEVEL_LIMIT - i - 1;
                break;
            }
        }
#else
        if (tmp.hi_ != 0x0) {
            for (S32 i = kLevMaxHi - 1; i >= 0; i--) {
                if ((((tmp.hi_ >> i * 3)) & 0x7) != 0) {
                    lev = kLevMaxHi - i - 1;
                    break;
                }
            }
        } else {
            for (S32 i = kLevMaxLo - 1; i >= 0; i--) {
                if ((((tmp.lo_ >> i * 3)) & 0x7) != 0) {
                    lev = kLevMaxHi + kLevMaxLo - i - 1;
                    break;
                }
            }
        }
#endif
        F64ort ret;
        const F64vec cen = center_;
        const F64vec hlen = half_len_;
        if (lev == 0) {
            ret.low_.x = cen.x - hlen.x;
            ret.low_.y = cen.y - hlen.y;
            ret.low_.z = cen.z - hlen.z;
            ret.high_.x = cen.x + hlen.x;
            ret.high_.y = cen.y + hlen.y;
            ret.high_.z = cen.z + hlen.z;
        } else {
            F64vec cen_new = cen;
            F64vec hlen_new = hlen;
            for (S32 i = 0; i < lev; i++) {
                S32 id = getCellID(i + 1, low);
                cen_new.x += hlen_new.x * SHIFT_CENTER[id].x;
                cen_new.y += hlen_new.y * SHIFT_CENTER[id].y;
                cen_new.z += hlen_new.z * SHIFT_CENTER[id].z;
                hlen_new *= 0.5;
            }
            ret.low_.x = cen_new.x - hlen_new.x;
            ret.low_.y = cen_new.y - hlen_new.y;
            ret.low_.z = cen_new.z - hlen_new.z;
            ret.high_.x = cen_new.x + hlen_new.x;
            ret.high_.y = cen_new.y + hlen_new.y;
            ret.high_.z = cen_new.z + hlen_new.z;
        }
        return ret;
    }
    F64 getCorrespondingFullLength(const S32 lev) const {
        F64vec len = half_len_ * 2.0;
        for (S32 i = 0; i < lev; i++) {
            len *= 0.5;
        }
        return std::max(std::max(len.x, len.y), len.z);
    }
    template <typename Tkey>
    F64ort getPosTreeCell(const S32 lev, const Tkey& key) const {
        F64vec cen = center_;
        F64vec hlen = half_len_;
        for (S32 i = 0; i < lev; i++) {
            S32 id = getCellID(i + 1, key);
            cen += hlen * SHIFT_CENTER[id];
            hlen *= 0.5;
        }
        F64ort ret;
        ret.low_.x = cen.x - hlen.x;
        ret.low_.y = cen.y - hlen.y;
        ret.low_.z = cen.z - hlen.z;
        ret.high_.x = cen.x + hlen.x;
        ret.high_.y = cen.y + hlen.y;
        ret.high_.z = cen.z + hlen.z;
        return ret;
    }
    // necessary functions for the PMMM feature
    S32vec getPMCellID(const F64vec pos) const {
        const KeyT mkey = getKey(pos);
        //std::cout << "pos= " << pos << std::endl;
        //std::cout << "mkey= " << mkey << std::endl;
        return getPMCellID(mkey);
    }
    template <typename TKey>
    S32vec getPMCellID(const TKey mkey) const {
        PS_COMPILE_TIME_MESSAGE("getPMCellID(Tkey)@TagMortonKeyMeshBased2");
        S32 ix, iy, iz;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        PS_COMPILE_TIME_MESSAGE("key is 64 bit getPMCellID(Tkey)@TagMortonKeyMeshBased2");
        // 64 bit
        const U64 x = (mkey.hi_ & 0x4924924924924924ULL) >> 2;
        const U64 y = (mkey.hi_ & 0x2492492492492492ULL) >> 1;
        const U64 z = mkey.hi_ & 0x1249249249249249ULL;
        ix = (S32)(shrinkBit(x) >> (kLevMax - lev_ref_)) - idx_ref_.x;
        iy = (S32)(shrinkBit(y) >> (kLevMax - lev_ref_)) - idx_ref_.y;
        iz = (S32)(shrinkBit(z) >> (kLevMax - lev_ref_)) - idx_ref_.z;
#else
        // 128 or 96 bit
        // std::cout<<"lev_ref_= "<<lev_ref_<<std::endl;
        // std::cout<<"idx_ref_= "<<idx_ref_<<std::endl;
        if (lev_ref_ <= kLevMaxHi) {
            const U64 x = (mkey.hi_ & 0x4924924924924924ULL) >> 2;  // 0 100^{21}
            const U64 y = (mkey.hi_ & 0x2492492492492492ULL) >> 1;  // 0 010^{21}
            const U64 z = mkey.hi_ & 0x1249249249249249ULL;         // 0 001^{21}
            ix = (S32)(shrinkBit(x) >> (kLevMaxHi - lev_ref_)) - idx_ref_.x;
            iy = (S32)(shrinkBit(y) >> (kLevMaxHi - lev_ref_)) - idx_ref_.y;
            iz = (S32)(shrinkBit(z) >> (kLevMaxHi - lev_ref_)) - idx_ref_.z;
        } else {
            const U64 x_hi = (mkey.hi_ & 0x4924924924924924ULL) >> 2;
            const U64 y_hi = (mkey.hi_ & 0x2492492492492492ULL) >> 1;
            const U64 z_hi = mkey.hi_ & 0x1249249249249249ULL;
            U64 ix_hi = shrinkBit(x_hi) << kLevMaxLo;  // operation << is to make room
            U64 iy_hi = shrinkBit(y_hi) << kLevMaxLo;
            U64 iz_hi = shrinkBit(z_hi) << kLevMaxLo;
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
            PS_COMPILE_TIME_MESSAGE("key is 96 bit getPMCellID(Tkey)@TagMortonKeyMeshBased2");
            // 96bit
            const U32 x_lo = (mkey.lo_ & 0x24924924) >> 2;  // 00 100^{10}
            const U32 y_lo = (mkey.lo_ & 0x12492492) >> 1;  // 00 010^{10}
            const U32 z_lo = mkey.lo_ & 0x09249249;         // 00 001^{10}
            U64 ix_lo = (shrinkBit32(x_lo) >> (kLevMaxLo - (lev_ref_ - kLevMaxHi)));
            U64 iy_lo = (shrinkBit32(y_lo) >> (kLevMaxLo - (lev_ref_ - kLevMaxHi)));
            U64 iz_lo = (shrinkBit32(z_lo) >> (kLevMaxLo - (lev_ref_ - kLevMaxHi)));
#else
            PS_COMPILE_TIME_MESSAGE("key is 128 bit getPMCellID(Tkey)@TagMortonKeyMeshBased2");
            // 128bit
            const U64 x_lo = (mkey.lo_ & 0x4924924924924924ULL) >> 2;
            const U64 y_lo = (mkey.lo_ & 0x2492492492492492ULL) >> 1;
            const U64 z_lo = mkey.lo_ & 0x1249249249249249ULL;
            U64 ix_lo = (shrinkBit(x_lo) >> (kLevMaxLo - (lev_ref_ - kLevMaxHi)));
            U64 iy_lo = (shrinkBit(y_lo) >> (kLevMaxLo - (lev_ref_ - kLevMaxHi)));
            U64 iz_lo = (shrinkBit(z_lo) >> (kLevMaxLo - (lev_ref_ - kLevMaxHi)));
#endif
            // std::cout<<"ix_hi= "<<ix_hi<<" ix_lo= "<<ix_lo<<" idx_ref_.x= "<<idx_ref_.x<<std::endl;
            ix = (S32)(ix_hi | ix_lo) - idx_ref_.x;
            iy = (S32)(iy_hi | iy_lo) - idx_ref_.y;
            iz = (S32)(iz_hi | iz_lo) - idx_ref_.z;
        }
#endif
        return S32vec(ix, iy, iz);
        /*
        S32 ix,iy,iz;
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY) || defined(PARTICLE_SIMULATOR_USE_128BIT_KEY)
        if (lev_ref_ <= kLevMaxHi) {
            const U64 x = (mkey.hi_ & 0x4924924924924924ULL)>>2; // 0 100^{21}
            const U64 y = (mkey.hi_ & 0x2492492492492492ULL)>>1; // 0 010^{21}
            const U64 z =  mkey.hi_ & 0x1249249249249249ULL;     // 0 001^{21}
            ix = (S32)(shrinkBit(x)>>(kLevMaxHi - lev_ref_)) - idx_ref_.x;
            iy = (S32)(shrinkBit(y)>>(kLevMaxHi - lev_ref_)) - idx_ref_.y;
            iz = (S32)(shrinkBit(z)>>(kLevMaxHi - lev_ref_)) - idx_ref_.z;
        } else {
            const U64 x_hi = (mkey.hi_ & 0x4924924924924924ULL)>>2;
            const U64 y_hi = (mkey.hi_ & 0x2492492492492492ULL)>>1;
            const U64 z_hi =  mkey.hi_ & 0x1249249249249249ULL;
            U64 ix_hi = shrinkBit(x_hi)<<kLevMaxLo; // operation << is to make room
            U64 iy_hi = shrinkBit(y_hi)<<kLevMaxLo;
            U64 iz_hi = shrinkBit(z_hi)<<kLevMaxLo;
#if defined(PARTICLE_SIMULATOR_USE_128BIT_KEY)
            const U64 x_lo = (mkey.lo_ & 0x4924924924924924ULL)>>2;
            const U64 y_lo = (mkey.lo_ & 0x2492492492492492ULL)>>1;
            const U64 z_lo =  mkey.lo_ & 0x1249249249249249ULL;
            U64 ix_lo = (shrinkBit(x_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
            U64 iy_lo = (shrinkBit(y_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
            U64 iz_lo = (shrinkBit(z_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
#else // defined(USE_128BIT_KEY)
            const U32 x_lo = (mkey.lo_ & 0x24924924)>>2; // 00 100^{10}
            const U32 y_lo = (mkey.lo_ & 0x12492492)>>1; // 00 010^{10}
            const U32 z_lo =  mkey.lo_ & 0x09249249;     // 00 001^{10}
            U64 ix_lo = (shrinkBit32(x_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
            U64 iy_lo = (shrinkBit32(y_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
            U64 iz_lo = (shrinkBit32(z_lo)>>(kLevMaxLo - (lev_ref_ - kLevMaxHi)));
#endif // defined(USE_128BIT_KEY)
            ix = (S32)(ix_hi | ix_lo) - idx_ref_.x;
            iy = (S32)(iy_hi | iy_lo) - idx_ref_.y;
            iz = (S32)(iz_hi | iz_lo) - idx_ref_.z;
        }
#else // defined(USE_96BIT_KEY) || defined(USE_128BIT_KEY)
        const U64 x = (mkey.hi_ & 0x4924924924924924ULL)>>2;
        const U64 y = (mkey.hi_ & 0x2492492492492492ULL)>>1;
        const U64 z =  mkey.hi_ & 0x1249249249249249ULL;
        ix = (S32)(shrinkBit(x)>>(kLevMax - lev_ref_)) - idx_ref_.x;
        iy = (S32)(shrinkBit(y)>>(kLevMax - lev_ref_)) - idx_ref_.y;
        iz = (S32)(shrinkBit(z)>>(kLevMax - lev_ref_)) - idx_ref_.z;
#endif //defined(USE_96BIT_KEY) || defined(USE_128BIT_KEY)
        return S32vec(ix, iy, iz);
        */
    }
};
}  // namespace ParticleSimulator
