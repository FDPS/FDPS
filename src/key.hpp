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
        MortonKey(){};
        ~MortonKey(){};
        MortonKey(const MortonKey &){};
        MortonKey & operator = (const MortonKey &);
        static MortonKey & getInstance(){
            static MortonKey inst;
            return inst;
        }
        static U64 separateBit(const U64 _s_in){
            U64 _s = _s_in;
            _s = (_s | _s<<32) &  0x00000000ffffffff;  //0000 0000 0000 0000 0000 0000 0000 0000 1111 1111 1111 1111 1111 1111 1111 1111 
            _s = (_s | _s<<16) &  0x0000ffff0000ffff;  //0000 0000 0000 0000 1111 1111 1111 1111 0000 0000 0000 0000 1111 1111 1111 1111 
            _s = (_s | _s<<8) &  0x00ff00ff00ff00ff;  //0000 0000 1111 1111 0000 0000 1111 1111
            _s = (_s | _s<<4) &   0x0f0f0f0f0f0f0f0f;  //0000 1111 0000 1111 0000 1111
            _s = (_s | _s<<2) &   0x3333333333333333;  //00 11 00 11 00 11
            return (_s | _s<<1) & 0x5555555555555555;  //0101 0101
        }
    public:
        static void initialize(const F64 half_len,
                               const F64vec & center=0.0){
            getInstance().half_len_ = half_len;
            getInstance().center_ = center;
            getInstance().normalized_factor_ = (1.0/(half_len*2.0))*(1<<kLevMax);
        }
        static U64 getKey(const F64vec & pos){
            const F64vec cen = getInstance().center_;
            const F64 hlen = getInstance().half_len_;
            const F64 nfactor = getInstance().normalized_factor_;
            const U64 nx = (U64)( (pos.x - cen.x + hlen) * nfactor);
            const U64 ny = (U64)( (pos.y - cen.y + hlen) * nfactor);
            //std::cerr<<"cen="<<cen<<" hlen="<<hlen<<std::endl;
            //std::cerr<<"nx="<<nx<<" ny="<<ny<<std::endl;
            return ( getInstance().separateBit(nx)<<1 | getInstance().separateBit(ny) );
        }
        static S32 getCellID(const S32 lev, const U64 mkey){
            U64 s = mkey >> ( (kLevMax - lev) * 2 );
            return (s & 0x3);
        }
    };


#else
    class MortonKey{
    private:
        enum{
            kLevMax = 21,
        };
        F64 half_len_;
        F64vec center_;
        F64 normalized_factor_;
        MortonKey(){};
        ~MortonKey(){};
        MortonKey(const MortonKey &);
        MortonKey & operator = (const MortonKey &);
        static MortonKey & getInstance(){
            static MortonKey inst;
            return inst;
        }
        static U64 separateBit(const U64 _s_in){
            U64 _s = _s_in;
            _s = (_s | _s<<32) & 0xffff00000000ffff; //11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
            _s = (_s | _s<<16) & 0x00ff0000ff0000ff; //00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
            _s = (_s | _s<<8) & 0xf00f00f00f00f00f;  //1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
            _s = (_s | _s<<4) & 0x30c30c30c30c30c3;  //11 00 00 11 00 11 00 00 11
            return (_s | _s<<2) & 0x9249249249249249;  //1 0 0 1 0 0 1 0 0 1 0 0 1
        }
    public:
        static void initialize(const F64 half_len,
                               const F64vec & center=0.0){
            getInstance().half_len_ = half_len;
            getInstance().center_ = center;
            getInstance().normalized_factor_ = (1.0/(half_len*2.0))*(1<<kLevMax);
        }
        static U64 getKey(F64vec pos){
            const F64vec cen = getInstance().center_;
            const F64 hlen = getInstance().half_len_;
            const F64 nfactor = getInstance().normalized_factor_;
            U64 nx = (U64)( (pos.x - cen.x + hlen) * nfactor);
            U64 ny = (U64)( (pos.y - cen.y + hlen) * nfactor);
            U64 nz = (U64)( (pos.z - cen.z + hlen) * nfactor);
            return ( getInstance().separateBit(nx)<<2 | getInstance().separateBit(ny)<<1 | getInstance().separateBit(nz) );
        }
        static S32 getCellID(const S32 lev, const U64 mkey){
            U64 s = mkey >> ( (kLevMax - lev) * 3 );
            return (s & 0x7);
        }
    };
    
#endif
    
}

