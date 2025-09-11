#include<cstdint>
#include"vector2.hpp"
#include"vector3.hpp"
#include"orthotope2.hpp"
#include"orthotope3.hpp"
#include"orthotope2i.hpp"
#include"orthotope3i.hpp"
#include"matrix_sym2.hpp"
#include"matrix_sym3.hpp"
#include"matrix2.hpp"

namespace ParticleSimulator{
    //typedef int S32;
    //typedef unsigned int U32;
    using S32 = std::int32_t;
    using U32 = std::uint32_t;
#if defined(PARTICLE_SIMULATOR_ALL_64BIT_PRECISION)
    typedef double F32;
#else
    typedef float F32;
#endif
    //typedef long long int S64;
    using S64 = std::int64_t;
    //typedef unsigned long long int U64;
    using U64 = std::uint64_t;
    typedef double F64;
    typedef Vector2<S32> S32vec2;
    typedef Vector3<S32> S32vec3;
    typedef Vector2<S64> S64vec2;
    typedef Vector3<S64> S64vec3;
    typedef Vector2<U64> U64vec2;
    typedef Vector3<U64> U64vec3;
    typedef Vector2<F32> F32vec2;
    typedef Vector3<F32> F32vec3;
    typedef Vector2<F64> F64vec2;
    typedef Vector3<F64> F64vec3;
    typedef MatrixSym2<F32> F32mat2;
    typedef MatrixSym3<F32> F32mat3;
    typedef MatrixSym2<F64> F64mat2;
    typedef MatrixSym3<F64> F64mat3;
    typedef Orthotope2i<S32> S32ort2;
    typedef Orthotope3i<S32> S32ort3;
    typedef Orthotope2<F32> F32ort2;
    typedef Orthotope3<F32> F32ort3;
    typedef Orthotope2<F64> F64ort2;
    typedef Orthotope3<F64> F64ort3;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    typedef S32vec2 S32vec;
    typedef S64vec2 S64vec;
    typedef F32vec2 F32vec;
    typedef F64vec2 F64vec;
    typedef F32mat2 F32mat;
    typedef F64mat2 F64mat;
    typedef S32ort2 S32ort;
    typedef F32ort2 F32ort;
    typedef F64ort2 F64ort;
#else
    typedef S32vec3 S32vec;
    typedef S64vec3 S64vec;
    typedef F32vec3 F32vec;
    typedef F64vec3 F64vec;
    typedef F32mat3 F32mat;
    typedef F64mat3 F64mat;
    typedef S32ort3 S32ort;
    typedef F32ort3 F32ort;
    typedef F64ort3 F64ort;
#endif

#if defined(PARTICLE_SIMULATOR_SPMOM_F32)
    typedef S32    SSP;
    typedef F32    FSP;
    typedef F32vec FSPvec;
    typedef F32mat FSPmat;
#else
    typedef S64    SSP;
    typedef F64    FSP;
    typedef F64vec FSPvec;
    typedef F64mat FSPmat;
#endif
    typedef U64 CountT;
    typedef CountT Count_t;
}
namespace PS = ParticleSimulator;
