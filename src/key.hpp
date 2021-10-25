#pragma once

#include<iostream>
#include<iomanip>
#include"vector3.hpp"
#include"ps_defs.hpp"

namespace ParticleSimulator{

    /*
      MSB is always 0.
      next 3bits represent octant index of 8 cells with level 1
     */
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    class MortonKey{
    private:
        enum{
            kLevMax = 30,
        };
        F64 half_len_;
        F64vec center_;
        F64 normalized_factor_;
        U64 separateBit(const U64 _s_in) const {
            U64 _s = _s_in;
            _s = (_s | _s<<32) &  0x00000000ffffffff;  //0000 0000 0000 0000 0000 0000 0000 0000 1111 1111 1111 1111 1111 1111 1111 1111 
            _s = (_s | _s<<16) &  0x0000ffff0000ffff;  //0000 0000 0000 0000 1111 1111 1111 1111 0000 0000 0000 0000 1111 1111 1111 1111 
            _s = (_s | _s<<8) &  0x00ff00ff00ff00ff;  //0000 0000 1111 1111 0000 0000 1111 1111
            _s = (_s | _s<<4) &   0x0f0f0f0f0f0f0f0f;  //0000 1111 0000 1111 0000 1111
            _s = (_s | _s<<2) &   0x3333333333333333;  //00 11 00 11 00 11
            return (_s | _s<<1) & 0x5555555555555555;  //0101 0101
        }
    public:
        //MortonKey(){};
        //~MortonKey(){};
        //MortonKey(const MortonKey &){};
        //MortonKey & operator = (const MortonKey &);
        void initialize(const F64 half_len,
                        const F64vec & center=0.0){
            half_len_ = half_len;
            center_ = center;
            normalized_factor_ = (1.0/(half_len*2.0))*(1<<kLevMax);
        }

        KeyT getKey(const F64vec & pos) const {
            const F64vec cen = center_;
            const F64 hlen = half_len_;
            const F64 nfactor = normalized_factor_;
            const U64 nx = (U64)( (pos.x - cen.x + hlen) * nfactor);
            const U64 ny = (U64)( (pos.y - cen.y + hlen) * nfactor);
            return ( separateBit(nx)<<1 | separateBit(ny) );
        }
        template<typename Tkey>
        S32 getCellID(const S32 lev, const Tkey mkey) const {
            const auto s = mkey >> ( (kLevMax - lev) * 2 );
            return (S32)((s.hi_ & 0x3));
        }
        F64ort getCorrespondingTreeCell(const F64ort & box) const {
            const auto low  = getKey(box.low_);
            const auto high = getKey(box.high_);
            const auto tmp = high^low;
            S32 lev = 0;
            for(S32 i=TREE_LEVEL_LIMIT-1; i>=0; i--){
                if( (((tmp >> i*2).hi_) & 0x3) != 0){
                    lev = TREE_LEVEL_LIMIT-i-1;
                    break;
                }
            }
            F64ort ret;
            const F64vec cen = center_;
            const F64 hlen = half_len_;
            if(lev==0){
                ret.low_.x  = cen.x - hlen;
                ret.low_.y  = cen.y - hlen;
                ret.high_.x = cen.x + hlen;
                ret.high_.y = cen.y + hlen;
            }
            else{
                F64vec cen_new = cen;
                F64 hlen_new = hlen;
                for(S32 i=0; i<lev; i++){
                    S32 id = getCellID(i+1, low);
                    cen_new +=  hlen_new*SHIFT_CENTER[id];
                    hlen_new *= 0.5;
                }
                ret.low_.x  = cen_new.x - hlen_new;
                ret.low_.y  = cen_new.y - hlen_new;
                ret.high_.x = cen_new.x + hlen_new;
                ret.high_.y = cen_new.y + hlen_new;
            }
            return ret;
        }
        F64 getCorrespondingFullLength(const S32 lev) const {
             F64 len = half_len_*2.0;
             for(S32 i=0; i<lev; i++){
                 len *= 0.5;
             }
             return len;
        }
        template<typename Tkey>
        F64ort getPosTreeCell(const S32 lev, const Tkey & key) const {
            F64vec cen = center_;
            F64 hlen = half_len_;
            for(S32 i=0; i<lev; i++){
                S32 id = getCellID(i+1, key);
                cen +=  hlen*SHIFT_CENTER[id];
                hlen *= 0.5;
            }
            F64ort ret;
            ret.low_.x  = cen.x - hlen;
            ret.low_.y  = cen.y - hlen;
            ret.high_.x = cen.x + hlen;
            ret.high_.y = cen.y + hlen;
            return ret;
        }
    };
#else
    // 3D
    class MortonKey{
    private:
        enum{
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
            kLevMax = 21,
            kLevMaxHi = 21,
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
    public:
        //MortonKey(){};
        //~MortonKey(){};
        //MortonKey(const MortonKey &);
        //MortonKey & operator = (const MortonKey &);

	static inline S32 keylevelmaxhigh(){return S32 (kLevMaxHi);}
	static inline S32 keylevelmaxlow(){return S32 (kLevMaxLo);}
        void initialize(const F64 half_len,
                        const F64vec & center=0.0){
            half_len_ = half_len;
            center_ = center;
            normalized_factor_ = (1.0/(half_len*2.0))*(1UL<<kLevMax);
        }
        KeyT getKey(F64vec pos) const {
            const F64vec cen = center_;
            const F64 hlen = half_len_;
            const F64 nfactor = normalized_factor_;
            U64 nx = (U64)( (pos.x - cen.x + hlen) * nfactor);
            U64 ny = (U64)( (pos.y - cen.y + hlen) * nfactor);
            U64 nz = (U64)( (pos.z - cen.z + hlen) * nfactor);
            
            // avoid to overflow
            nx = (((nx>>63)&0x1)==1) ? 0 : nx;
            ny = (((ny>>63)&0x1)==1) ? 0 : ny;
            nz = (((nz>>63)&0x1)==1) ? 0 : nz;
            U64 nmax = (1UL<<kLevMax)-1;
            nx = (nx > nmax) ? nmax : nx;
            ny = (ny > nmax) ? nmax : ny;
            nz = (nz > nmax) ? nmax : nz;
            
            KeyT ret;
#if defined (PARTICLE_SIMULATOR_USE_64BIT_KEY)
            ret.hi_ = separateBit(nx)<<2 | separateBit(ny)<<1 | separateBit(nz);
#else
            U64 nx_hi = nx>>kLevMaxLo;
            U64 ny_hi = ny>>kLevMaxLo;
            U64 nz_hi = nz>>kLevMaxLo;
            ret.hi_ = (separateBit(nx_hi)<<2
                       | separateBit(ny_hi)<<1
                       | separateBit(nz_hi));
#if defined (PARTICLE_SIMULATOR_USE_96BIT_KEY)
            U32 nx_lo = (U32)(nx & 0x3ff); // mask 10 bits
            U32 ny_lo = (U32)(ny & 0x3ff);
            U32 nz_lo = (U32)(nz & 0x3ff);
            ret.lo_ = ( separateBit32(nx_lo)<<2
                        | separateBit32(ny_lo)<<1
                        | separateBit32(nz_lo) );
#else
            U64 nx_lo = (U64)(nx & 0x1fffff); // mask 21bits
            U64 ny_lo = (U64)(ny & 0x1fffff);
            U64 nz_lo = (U64)(nz & 0x1fffff);
            ret.lo_ = ( separateBit(nx_lo)<<2
                        | separateBit(ny_lo)<<1
                        | separateBit(nz_lo) );
#endif
#endif
            return ret;
        }

        template<typename Tkey>
        S32 getCellID(const S32 lev, const Tkey & mkey) const {
            U64 s;
#if defined (PARTICLE_SIMULATOR_USE_64BIT_KEY)
            s = mkey.hi_ >> ( (kLevMax - lev) * 3 );
#else
            if(lev <= kLevMaxHi){
                s = mkey.hi_ >> ( (kLevMaxHi - lev) * 3 );
            }
            else{
                s = mkey.lo_ >> ( (kLevMaxLo - (lev-kLevMaxHi)) * 3 );
            }
#endif
            return (S32)(s & 0x7);
        }
        
        template<typename Tkey>
        S32 getCellIDlow(const S32 lev, const Tkey & mkey) const {
            U64 s;
#if defined (PARTICLE_SIMULATOR_USE_64BIT_KEY)
            s = mkey.hi_ >> ( (kLevMax - lev) * 3 );
#else
	    //	    s = mkey.lo_ >> ( (kLevMaxLo - (lev-kLevMaxHi)) * 3 );
	    s = mkey.lo_ >> lev;
#endif
            return (S32)(s & 0x7);
        }
        
        template<typename Tkey>
        S32 getCellIDhigh(const S32 lev, const Tkey & mkey) const {
            U64 s;
#if defined (PARTICLE_SIMULATOR_USE_64BIT_KEY)
            s = 0;
#else
	    //	    s = mkey.hi_ >> ( (kLevMaxHi - lev) * 3 );
	    	    s = mkey.hi_ >> lev;
#endif
            return (S32)(s & 0x7);
        }
        
        F64ort getCorrespondingTreeCell(const F64ort & box) const {
            const auto low  = getKey(box.low_);
            const auto high = getKey(box.high_);
            const auto tmp = high^low;
            S32 lev = 0;	
#if defined (PARTICLE_SIMULATOR_USE_64BIT_KEY)
            for(S32 i=TREE_LEVEL_LIMIT-1; i>=0; i--){
                if( (((tmp.hi_ >> i*3)) & 0x7) != 0){
                    lev = TREE_LEVEL_LIMIT-i-1;
                    break;
                }
            }
#else
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
#endif
            F64ort ret;
            const F64vec cen = center_;
            const F64 hlen = half_len_;
            if(lev==0){
                ret.low_.x  = cen.x - hlen;
                ret.low_.y  = cen.y - hlen;
                ret.low_.z  = cen.z - hlen;
                ret.high_.x = cen.x + hlen;
                ret.high_.y = cen.y + hlen;
                ret.high_.z = cen.z + hlen;
            }
            else{
                F64vec cen_new = cen;
                F64 hlen_new = hlen;
                for(S32 i=0; i<lev; i++){
                    S32 id = getCellID(i+1, low);
                    cen_new +=  hlen_new*SHIFT_CENTER[id];
                    hlen_new *= 0.5;
                }
                ret.low_.x  = cen_new.x - hlen_new;
                ret.low_.y  = cen_new.y - hlen_new;
                ret.low_.z  = cen_new.z - hlen_new;
                ret.high_.x = cen_new.x + hlen_new;
                ret.high_.y = cen_new.y + hlen_new;
                ret.high_.z = cen_new.z + hlen_new;
            }
            return ret;
        }
        F64 getCorrespondingFullLength(const S32 lev) const {
            F64 len = half_len_*2.0;
            for(S32 i=0; i<lev; i++){
                len *= 0.5;
            }
            return len;
        }
        template<typename Tkey>
        F64ort getPosTreeCell(const S32 lev, const Tkey & key) const {
            F64vec cen = center_;
            F64 hlen = half_len_;
            for(S32 i=0; i<lev; i++){
                S32 id = getCellID(i+1, key);
                cen +=  hlen*SHIFT_CENTER[id];
                hlen *= 0.5;
            }
            F64ort ret;
            ret.low_.x  = cen.x - hlen;
            ret.low_.y  = cen.y - hlen;
            ret.low_.z  = cen.z - hlen;
            ret.high_.x = cen.x + hlen;
            ret.high_.y = cen.y + hlen;
            ret.high_.z = cen.z + hlen;
            return ret;
        }
    };
    
#endif
    
}

