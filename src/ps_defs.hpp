#pragma once

#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <cassert>
#include <typeinfo>
#include <cstdio>
#include <cstring>
#include <map>
#include <random>

#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
#include "mpi.h"
#endif

#define PS_CPLUSPLUS_20 202002L
#define PS_CPLUSPLUS_17 201703L
#define PS_CPLUSPLUS_14 201402L
#define PS_CPLUSPLUS_11 201103L

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include <omp.h>
#define PS_OMP(...) _Pragma(#__VA_ARGS__)
#define PS_OMP_PARALLEL_FOR _Pragma("omp parallel for")
#define PS_OMP_PARALLEL_FOR_G _Pragma("omp parallel for schedule(guided)")
#define PS_OMP_PARALLEL_FOR_D _Pragma("omp parallel for schedule(dynamic,1)")
#define PS_OMP_PARALLEL_FOR_D4 _Pragma("omp parallel for schedule(dynamic,4)")
#define PS_OMP_PARALLEL_FOR_S _Pragma("omp parallel for schedule(static)")

#define PS_OMP_PARALLEL _Pragma("omp parallel")
#define PS_OMP_FOR _Pragma("omp for")
#define PS_OMP_FOR_D _Pragma("omp for schedule(dynamic, 1)")
#define PS_OMP_FOR_D4 _Pragma("omp for schedule(dynamic, 4)")
#define PS_OMP_FOR_NOWAIT _Pragma("omp for nowait")
#define PS_OMP_CRITICAL _Pragma("omp critical")
#define PS_OMP_BARRIER _Pragma("omp barrier")
#define PS_OMP_SINGLE _Pragma("omp single")
#else
#define PS_OMP(...)
#define PS_OMP_PARALLEL_FOR
#define PS_OMP_PARALLEL_FOR_D4
#define PS_OMP_PARALLEL_FOR_D
#define PS_OMP_PARALLEL_FOR_G
#define PS_OMP_PARALLEL_FOR_S
#define PS_OMP_PARALLEL
#define PS_OMP_FOR
#define PS_OMP_FOR_D
#define PS_OMP_FOR_D4
#define PS_OMP_FOR_NOWAIT
#define PS_OMP_CRITICAL
#define PS_OMP_BARRIER
#define PS_OMP_SINGLE
#endif

#ifdef __GNUC__
#define PS_INLINE inline __attribute__((always_inline))
#else
#define PS_INLINE inline
#endif

#define SET_VAR_NAME(var)  \
    do {                   \
        var.setName(#var); \
    } while (0);

#define PS_PRINT_VARIABLE(var)                           \
    do {                                                 \
        std::cout << #var << "= " << (var) << std::endl; \
    } while (0);

#define PS_PRINT_VARIABLE2(var, var2)                                                      \
    do {                                                                                   \
        std::cout << #var << "= " << (var) << " " << #var2 << "= " << (var2) << std::endl; \
    } while (0);

#define PS_PRINT_VARIABLE3(var, var2, var3)                                                                                  \
    do {                                                                                                                     \
        std::cout << #var << "= " << (var) << " " << #var2 << "= " << (var2) << " " << #var3 << "= " << (var3) << std::endl; \
    } while (0);

#define PS_PRINT_VARIABLE4(var, var2, var3, var4)                                                                                                \
    do {                                                                                                                                         \
        std::cout << #var << "= " << (var) << " " << #var2 << "= " << (var2) << " " << #var3 << "= " << (var3) << " " << #var4 << "= " << (var4) \
                  << std::endl;                                                                                                                  \
    } while (0);

#define PS_PRINT_VARIABLE5(var, var2, var3, var4, var5)                                                                                          \
    do {                                                                                                                                         \
        std::cout << #var << "= " << (var) << " " << #var2 << "= " << (var2) << " " << #var3 << "= " << (var3) << " " << #var4 << "= " << (var4) \
                  << " " << #var5 << "= " << (var5) << std::endl;                                                                                \
    } while (0);

#define PS_PRINT_VARIABLE6(var, var2, var3, var4, var5, var6)                                                                                    \
    do {                                                                                                                                         \
        std::cout << #var << "= " << (var) << " " << #var2 << "= " << (var2) << " " << #var3 << "= " << (var3) << " " << #var4 << "= " << (var4) \
                  << " " << #var5 << "= " << (var5) << " " << #var6 << "= " << (var6) << std::endl;                                              \
    } while (0);

#if defined(PARTICLE_SIMULATOR_DEBUG_MODE)
#define PS_COMPILE_TIME_MESSAGE(msg) \
    do {                             \
        asm(".print " #msg);         \
    } while (0);

#define PS_SET_CALLER_INFO(obj)                              \
    do {                                                     \
        obj.setCallerInfo(__LINE__, __FUNCTION__, __FILE__); \
    } while (0);

#else

#define PS_COMPILE_TIME_MESSAGE(msg) \
    do {                             \
    } while (0);

#define PS_SET_CALLER_INFO(obj)                              \
    do {                                                     \
    } while (0);
    
#endif


#include "memory_pool.hpp"
#include "ps_types.hpp"

#define PS_DEBUG_CALL(func)                                                                                                                    \
    do {                                                                                                                                       \
        try {                                                                                                                                  \
            func;                                                                                                                              \
            std::cout << "[FDPS msg] " << #func << ": " << getMemSizeUsed() << ", " << Comm::getRank() << ", " << typeid(TSM).name() << ". "   \
                      << std::endl;                                                                                                            \
        } catch (std::bad_alloc & e) {                                                                                                         \
            std::cout << "[FDPS error] " << #func << ": " << getMemSizeUsed() << ", " << Comm::getRank() << ", " << typeid(TSM).name() << ". " \
                      << std::endl;                                                                                                            \
            MPI_Abort(MPI_COMM_WORLD, 9);                                                                                                      \
            std::exit(1);                                                                                                                      \
        } catch (...) {                                                                                                                        \
            std::cout << "[FDPS unknown error] " << #func << ": " << getMemSizeUsed() << ", " << Comm::getRank() << ", " << typeid(TSM).name() \
                      << ". " << std::endl;                                                                                                    \
            MPI_Abort(MPI_COMM_WORLD, 9);                                                                                                      \
            std::exit(1);                                                                                                                      \
        }                                                                                                                                      \
        MPI_Barrier(MPI_COMM_WORLD);                                                                                                           \
        if (Comm::getRank() == 0) std::cout << #func << " passed." << std::endl;                                                               \
    } while (0);

#define PARTICLE_SIMULATOR_PRINT_ERROR(msg) \
    { std::cout << "PS_ERROR: " << msg << " \n" << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl; }

#define PARTICLE_SIMULATOR_PRINT_LINE_INFO() \
    { std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl; }

#define PARTICLE_SIMULATOR_PRINT_MSG(msg, rank)                                                                                                   \
    do {                                                                                                                                          \
        if (Comm::getRank() == (rank)) {                                                                                                          \
            std::cout << "# " << msg << " function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << " ***" << std::endl; \
        }                                                                                                                                         \
    } while (0);

// This function is the same as PARTICLE_SIMULATOR_PRINT_MSG
#define PS_PRINT_MSG(msg, rank)                                                                                                                   \
    do {                                                                                                                                          \
        if (PS::Comm::getRank() == (rank)) {                                                                                                      \
            std::cout << "# " << msg << " function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << " ***" << std::endl; \
        }                                                                                                                                         \
    } while (0);

#define PARTICLE_SIMULATOR_PRINT_MSG_B(msg, rank)                                                                                                 \
    do {                                                                                                                                          \
        Comm::barrier();                                                                                                                          \
        if (Comm::getRank() == (rank)) {                                                                                                          \
            std::cout << "# " << msg << " function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << " ***" << std::endl; \
        }                                                                                                                                         \
        Comm::barrier();                                                                                                                          \
    } while (0);

template <typename Tcomm>
void DEBUG_PRINT_MAKING_TREE(Tcomm &comm_info) {
#if defined(DEBUG_PRINT_MAKING_TREE)
    comm_info.barrier();
    if (comm_info.getRank() == 0) PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
}

namespace ParticleSimulator {

static constexpr long long int LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1ll << 30;
static inline void Abort(const int err = -1) {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
    MPI_Abort(MPI_COMM_WORLD, err);
#else
    exit(err);
#endif
}
static constexpr S32 DIMENSION_LIMIT = 3;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
static constexpr S32 DIMENSION = 2;
static constexpr S32 N_CHILDREN = 4;
static constexpr S32 TREE_LEVEL_LIMIT = 30;
static constexpr F64vec SHIFT_CENTER[N_CHILDREN] = {F64vec(-0.5, -0.5), F64vec(-0.5, 0.5), F64vec(0.5, -0.5), F64vec(0.5, 0.5)};
#else
static constexpr S32 DIMENSION = 3;
static constexpr S32 N_CHILDREN = 8;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
static constexpr S32 TREE_LEVEL_LIMIT = 21;  // 21*3=63bit
#elif defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
static constexpr S32 TREE_LEVEL_LIMIT = 31;
static constexpr S32 KEY_LEVEL_MAX_HI = 21;  // 21*3=63bit
static constexpr S32 KEY_LEVEL_MAX_LO = 10;  // 10*3=30bit
#else
static constexpr S32 TREE_LEVEL_LIMIT = 42;
static constexpr S32 KEY_LEVEL_MAX_HI = 21;  // 21*3=63bit
static constexpr S32 KEY_LEVEL_MAX_LO = 21;  // 21*3=63bit
#endif
static constexpr F64vec SHIFT_CENTER[N_CHILDREN] = {F64vec(-0.5, -0.5, -0.5), F64vec(-0.5, -0.5, 0.5), F64vec(-0.5, 0.5, -0.5),
                                                    F64vec(-0.5, 0.5, 0.5),   F64vec(0.5, -0.5, -0.5), F64vec(0.5, -0.5, 0.5),
                                                    F64vec(0.5, 0.5, -0.5),   F64vec(0.5, 0.5, 0.5)};
#endif
}  // namespace ParticleSimulator

#include "reallocatable_array.hpp"
#include "comm_info.hpp"

namespace ParticleSimulator {

typedef U64 CountT;
typedef CountT Count_t;
static const F64 LARGE_DOUBLE = std::numeric_limits<F64>::max() * 0.0625;
static const F64 LARGE_FLOAT = std::numeric_limits<F32>::max() * 0.0625;
static const S64 LARGE_INT = std::numeric_limits<S32>::max() * 0.0625;

///// A.Tanikawa modified from
// In the upper line, the right-hand side is interpreted to be a 32bit-integer.
// static const S64 LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1<<31;
// static const S64 LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1ll<<30;
///// A.Tanikawa modified to

//////////////////
/// enum
enum EXCHANGE_LET_MODE {
    EXCHANGE_LET_A2A,
    EXCHANGE_LET_P2P_EXACT,
    EXCHANGE_LET_P2P_FAST,
};

enum INTERACTION_LIST_MODE {
    MAKE_LIST,
    MAKE_LIST_FOR_REUSE,
    REUSE_LIST,
};

inline void Check_INTERACTION_LIST_MODE(const INTERACTION_LIST_MODE mode) {
    if (mode != MAKE_LIST && mode != MAKE_LIST_FOR_REUSE && mode != REUSE_LIST) {
        PARTICLE_SIMULATOR_PRINT_ERROR("INTERACTION_LIST_MODE is invalid");
        Abort(-1);
    }
}

enum SET_ROOT_CELL_MODE {
    SET_ROOT_CELL_AUTO,
    SET_ROOT_CELL_CENTER_LENGHT,
};

enum TREE_TYPE {
    LOCAL_TREE = 0,
    LET_TREE,
    GLOBAL_TREE,
    N_TREE_TYPE,
};

enum SEARCH_MODE {
    LONG_NO_CUTOFF,
    LONG_CUTOFF,
    LONG_SCATTER,         // new for P^3T
    LONG_CUTOFF_SCATTER,  // new for P^3T + PM
    LONG_SYMMETRY,
    SHORT_GATHER,
    SHORT_SCATTER,
    SHORT_SYMMETRY,
};

enum FORCE_TYPE {
    FORCE_TYPE_LONG,
    FORCE_TYPE_SHORT,
};

enum CALC_DISTANCE_TYPE {
    CALC_DISTANCE_TYPE_NORMAL = 0,
    CALC_DISTANCE_TYPE_NEAREST_X = 1,
    CALC_DISTANCE_TYPE_NEAREST_Y = 2,
    CALC_DISTANCE_TYPE_NEAREST_XY = 3,
    CALC_DISTANCE_TYPE_NEAREST_Z = 4,
    CALC_DISTANCE_TYPE_NEAREST_XZ = 5,
    CALC_DISTANCE_TYPE_NEAREST_YZ = 6,
    CALC_DISTANCE_TYPE_NEAREST_XYZ = 7,
};

struct TagClassTypeParticleSystem {};
struct TagClassTypeDomainInfo {};
struct TagClassTypeTreeForForce {};
struct TagClassTypeParticleMesh {};
struct TagClassTypeParticleMeshMultipole {};

struct TagForceLong {
    enum {
        force_type = FORCE_TYPE_LONG,
    };
};
struct TagForceShort {
    enum {
        force_type = FORCE_TYPE_SHORT,
    };
};

// GEOMETRY OF TREE CELL
struct TagGeometrySize {};
struct TagGeometrySizeIn {};  // for PMMM
struct TagGeometrySizeOut {};
struct TagGeometrySizeInOut {};
struct TagGeometryIn {};
struct TagGeometryOut {};
struct TagGeometryInAndOut {};

// LOCAL TREE MOMENT TYPE
struct TagCalcMomLongEpjLt {};
struct TagCalcMomShortEpjLt {};
struct TagCalcMomShortEpiLt {};

// MORTON KEY GENERATION TYPE
struct TagMortonKeyNormal {};
struct TagMortonKeyMeshBased {};
struct TagMortonKeyMeshBased2 {};  // new

// IS PERIODIC BOUNDARY CONDITION AVAILABLE
struct TagSearchBoundaryConditionOpenOnly {};
struct TagSearchBoundaryConditionOpenPeriodic {};

// SEARCH TYPE
struct TagSearchLong {};
struct TagSearchLongCutoff {};
struct TagSearchLongScatter {};
struct TagSearchLongSymmetry {};
struct TagSearchLongCutoffScatter {};
struct TagSearchShortGather {};
struct TagSearchShortScatter {};
struct TagSearchShortSymmetry {};
struct TagSearchLongParticleMeshMultipole {};  // for PMMM

// IP-GROUP TYPE
struct TagIpgLongNormal {};
struct TagIpgIn {};
struct TagIpgInAndOut {};
struct TagIpgOut {};

struct TagWithoutCutoff {};
struct TagWithCutoff {};

// NEIGHBOR SEARCH TYPE
struct TagNeighborSearchSymmetry {};
struct TagNeighborSearchGather {};
struct TagNeighborSearchScatter {};
struct TagNeighborSearchNo {};

// TO calculate center and length of root cell
struct TagCalcCenterAndLengthUseEPI {};
struct TagCalcCenterAndLengthUseEPJ {};

// reposition_pos_root_cell_type
enum class RepositionPosRootCellType { No, PMMM };

enum class LET_USAGE_MODE { CONSTRUCT_GLOBAL_TREE, CONSTRUCT_LET_TREE };

struct SEARCH_MODE_LONG {
    typedef TagForceLong force_type;
    using mkey_type = TagMortonKeyNormal;
    typedef TagSearchLong search_type;
    typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
    typedef TagIpgLongNormal ipg_type;
    typedef TagNeighborSearchNo neighbor_search_type;
    using tree_cell_loc_geometry_type = TagGeometrySize;
    using tree_cell_glb_geometry_type = TagGeometrySize;
    using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
    using calc_center_and_length_type = TagCalcCenterAndLengthUseEPJ;
    static const inline auto reposition_pos_root_cell_type = RepositionPosRootCellType::No;
    // using reposition_pos_root_cell_type = TagRepositionPosRootCellNo;
};

struct SEARCH_MODE_LONG_CUTOFF {
    typedef TagForceLong force_type;
    using mkey_type = TagMortonKeyNormal;
    typedef TagSearchLongCutoff search_type;
    typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
    typedef TagIpgLongNormal ipg_type;
    typedef TagNeighborSearchNo neighbor_search_type;
    using tree_cell_loc_geometry_type = TagGeometrySizeOut;
    using tree_cell_glb_geometry_type = TagGeometrySizeOut;
    using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
    using calc_center_and_length_type = TagCalcCenterAndLengthUseEPJ;
    static const inline auto reposition_pos_root_cell_type = RepositionPosRootCellType::No;
    // inline RepositionPosRootCellType reposition_pos_root_cell_type = No;
    // using reposition_pos_root_cell_type = TagRepositionPosRootCellNo;
};
struct SEARCH_MODE_LONG_SCATTER {
    typedef TagForceLong force_type;
    using mkey_type = TagMortonKeyNormal;
    typedef TagSearchLongScatter search_type;
    typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
    typedef TagIpgIn ipg_type;
    typedef TagNeighborSearchScatter neighbor_search_type;
    using tree_cell_loc_geometry_type = TagGeometrySizeInOut;
    using tree_cell_glb_geometry_type = TagGeometrySizeInOut;
    using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
    using calc_center_and_length_type = TagCalcCenterAndLengthUseEPJ;
    static const inline auto reposition_pos_root_cell_type = RepositionPosRootCellType::No;
    // inline RepositionPosRootCellType reposition_pos_root_cell_type = No;
    // using reposition_pos_root_cell_type = TagRepositionPosRootCellNo;
};
struct SEARCH_MODE_LONG_SYMMETRY {
    typedef TagForceLong force_type;
    using mkey_type = TagMortonKeyNormal;
    typedef TagSearchLongSymmetry search_type;
    typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
    typedef TagIpgInAndOut ipg_type;
    typedef TagNeighborSearchSymmetry neighbor_search_type;
    using tree_cell_loc_geometry_type = TagGeometrySizeInOut;
    using tree_cell_glb_geometry_type = TagGeometrySizeInOut;
    using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
    using calc_center_and_length_type = TagCalcCenterAndLengthUseEPJ;
    static const inline auto reposition_pos_root_cell_type = RepositionPosRootCellType::No;
    // inline RepositionPosRootCellType reposition_pos_root_cell_type = No;
    // using reposition_pos_root_cell_type = TagRepositionPosRootCellNo;
};

struct SEARCH_MODE_GATHER {
    typedef TagForceShort force_type;
    using mkey_type = TagMortonKeyNormal;
    typedef TagSearchShortGather search_type;
    typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
    typedef TagIpgOut ipg_type;
    // typedef TagNeighborSearchGather neighbor_search_type;
    using neighbor_search_type = TagNeighborSearchGather;
    using tree_cell_loc_geometry_type = TagGeometryInAndOut;
    using tree_cell_glb_geometry_type = TagGeometryIn;
    using calc_moment_local_tree_type = TagCalcMomShortEpiLt;
    using calc_center_and_length_type = TagCalcCenterAndLengthUseEPI;
    static const inline auto reposition_pos_root_cell_type = RepositionPosRootCellType::No;
    // inline RepositionPosRootCellType reposition_pos_root_cell_type = No;
    // using reposition_pos_root_cell_type = TagRepositionPosRootCellNo;
};
struct SEARCH_MODE_SCATTER {
    typedef TagForceShort force_type;
    using mkey_type = TagMortonKeyNormal;
    typedef TagSearchShortScatter search_type;
    typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
    typedef TagIpgIn ipg_type;
    // typedef TagNeighborSearchScatter neighbor_search_type;
    using neighbor_search_type = TagNeighborSearchScatter;
    using tree_cell_loc_geometry_type = TagGeometryInAndOut;
    using tree_cell_glb_geometry_type = TagGeometryInAndOut;
    using calc_moment_local_tree_type = TagCalcMomShortEpjLt;
    using calc_center_and_length_type = TagCalcCenterAndLengthUseEPJ;
    static const inline auto reposition_pos_root_cell_type = RepositionPosRootCellType::No;
    // inline RepositionPosRootCellType reposition_pos_root_cell_type = No;
    // using reposition_pos_root_cell_type = TagRepositionPosRootCellNo;
};
struct SEARCH_MODE_SYMMETRY {
    typedef TagForceShort force_type;
    using mkey_type = TagMortonKeyNormal;
    typedef TagSearchShortSymmetry search_type;
    typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
    typedef TagIpgInAndOut ipg_type;
    typedef TagNeighborSearchSymmetry neighbor_search_type;
    using tree_cell_loc_geometry_type = TagGeometryInAndOut;
    using tree_cell_glb_geometry_type = TagGeometryInAndOut;
    using calc_moment_local_tree_type = TagCalcMomShortEpjLt;
    using calc_center_and_length_type = TagCalcCenterAndLengthUseEPJ;
    static const inline auto reposition_pos_root_cell_type = RepositionPosRootCellType::No;
    // inline RepositionPosRootCellType reposition_pos_root_cell_type = No;
    // using reposition_pos_root_cell_type = TagRepositionPosRootCellNo;
};
struct SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE {
    using force_type = TagForceLong;
    // using mkey_type  = TagMortonKeyMeshBased;
    using mkey_type = TagMortonKeyMeshBased2;
    using search_type = TagSearchLongParticleMeshMultipole;
    using search_boundary_type = TagSearchBoundaryConditionOpenPeriodic;
    using ipg_type = TagIpgLongNormal;
    using neighbor_search_type = TagNeighborSearchNo;
    using tree_cell_loc_geometry_type = TagGeometrySizeIn;
    using tree_cell_glb_geometry_type = TagGeometrySizeIn;
    using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
    using calc_center_and_length_type = TagCalcCenterAndLengthUseEPJ;
    static const inline auto reposition_pos_root_cell_type = RepositionPosRootCellType::PMMM;
    // inline RepositionPosRootCellType reposition_pos_root_cell_type = PMMM;
    // using reposition_pos_root_cell_type = TagRepositionPosRootCellPMM;
};

template <class T>
class ValueTypeReduction;
template <>
class ValueTypeReduction<float> {};
template <>
class ValueTypeReduction<double> {};
template <>
class ValueTypeReduction<int> {};
template <>
class ValueTypeReduction<long> {};

enum BOUNDARY_CONDITION {
    BOUNDARY_CONDITION_OPEN,
    BOUNDARY_CONDITION_PERIODIC_X,
    BOUNDARY_CONDITION_PERIODIC_Y,
    BOUNDARY_CONDITION_PERIODIC_Z,
    BOUNDARY_CONDITION_PERIODIC_XY,
    BOUNDARY_CONDITION_PERIODIC_XZ,
    BOUNDARY_CONDITION_PERIODIC_YZ,
    BOUNDARY_CONDITION_PERIODIC_XYZ,
    BOUNDARY_CONDITION_SHEARING_BOX,
    BOUNDARY_CONDITION_USER_DEFINED,
};

static inline void Initialize(int &argc, char **&argv, const S64 mpool_size = 1000000000) {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
    MPI_Init(&argc, &argv);
#endif
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
    omp_init_lock(&PS_OMP_LOCK);
#endif
    MemoryPool::initialize(mpool_size);
#if defined(MONAR)
    bool flag_monar = false;
    bool flag_MONAR = false;
    for (S32 i = 0; i < argc; i++) {
        if (strcmp(argv[i], "monar") == 0) flag_monar = true;
        if (strcmp(argv[i], "MONAR") == 0) flag_MONAR = true;
    }
#endif

    if (Comm::getRank() == 0) {
        std::cerr << "     //==================================\\\\" << std::endl;
        std::cerr << "     ||                                  ||" << std::endl;
        std::cerr << "     || ::::::: ::::::. ::::::. .::::::. ||" << std::endl;
        std::cerr << "     || ::      ::    : ::    : ::       ||" << std::endl;
        std::cerr << "     || ::::::  ::    : ::::::'  `:::::. ||" << std::endl;
        std::cerr << "     || ::      ::::::' ::      `......' ||" << std::endl;
        std::cerr << "     ||     Framework for Developing     ||" << std::endl;
        std::cerr << "     ||        Particle Simulator        ||" << std::endl;
        std::cerr << "     ||     Version 8.0 (2025/09)       ||" << std::endl;
        std::cerr << "     \\\\==================================//" << std::endl;
        std::cerr << "" << std::endl;
        std::cerr << "       Home   : https://github.com/fdps/fdps " << std::endl;
        std::cerr << "       E-mail : fdps-support@mail.jmlab.jp" << std::endl;
        std::cerr << "       Licence: MIT (see, https://github.com/FDPS/FDPS/blob/master/LICENSE)" << std::endl;
        std::cerr << "       Note   : Please cite the following papers." << std::endl;
        std::cerr << "                - Iwasawa et al. (2016, Publications of the Astronomical Society of Japan, 68, 54)" << std::endl;
        std::cerr << "                - Namekata et al. (2018, Publications of the Astronomical Society of Japan, 70, 70)" << std::endl;
        std::cerr << "                - Iwasawa et al. (2020, Publications of the Astronomical Society of Japan, 72, 13)" << std::endl;
        std::cerr << "" << std::endl;
        std::cerr << "       Copyright (C) 2015 " << std::endl;
        std::cerr << "         Masaki Iwasawa, Ataru Tanikawa, Natsuki Hosono," << std::endl;
        std::cerr << "         Keigo Nitadori, Takayuki Muranushi, Daisuke Namekata," << std::endl;
        std::cerr << "         Kentaro Nomura, Junichiro Makino and many others" << std::endl;
#ifdef MONAR
        if (flag_monar) {
            std::cerr << "　　 ^__^　 ／‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾" << std::endl;
            std::cerr << "　　( ´∀｀)＜******** FDPS has successfully begun. ********" << std::endl;
            std::cerr << "　　(     ) ＼＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿" << std::endl;
            std::cerr << "　　|  | |" << std::endl;
            std::cerr << "　　(__)_)" << std::endl;
        } else if (flag_MONAR) {
            std::cerr << "        ∧_∧   ／‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾" << std::endl;
            std::cerr << "       (´Д`) <  ******** FDPS has successfully begun. ********" << std::endl;
            std::cerr << "       ／_ /  ＼" << std::endl;
            std::cerr << "      (ぃ９｜  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾" << std::endl;
            std::cerr << "      /　　/、" << std::endl;
            std::cerr << "     /　　∧_二つ" << std::endl;
            std::cerr << "     ｜　　＼ " << std::endl;
            std::cerr << "     /　/~＼ ＼" << std::endl;
            std::cerr << "    /　/　　>　)" << std::endl;
            std::cerr << "   ／ ノ　　/ ／" << std::endl;
            std::cerr << "  / ／　　 / /" << std::endl;
            std::cerr << "`/ /　　　( ヽ" << std::endl;
            std::cerr << "(＿)　　　 ＼_)" << std::endl;
        } else {
            fprintf(stderr, "******** FDPS has successfully begun. ********\n");
        }
#else   // MONAR
        fprintf(stderr, "******** FDPS has successfully begun. ********\n");
#endif  // MONAR
    }
}

static inline void Finalize() {
    MemoryPool::finalize();
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
    MPI_Finalize();
#endif
    bool flag_monar = false;
    if (Comm::getRank() == 0) {
        if (flag_monar) {
        } else {
            fprintf(stderr, "******** FDPS has successfully finished. ********\n");
        }
    }
}

static const S32 N_SMP_PTCL_TOT_PER_PSYS_DEFAULT = 1000000;

////////
// util
inline F64 GetWtime() {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
#ifdef PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
    Comm::barrier();
#endif  // PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
    return MPI_Wtime();
#elif defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
    // PARTICLE_SIMULATOR_THREAD_PARALLEL
    return omp_get_wtime();
#else
    return (F64)clock() / CLOCKS_PER_SEC;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
}

inline F64 GetWtimeNoBarrier() {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
    return MPI_Wtime();
#elif defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
    return omp_get_wtime();
#else
    return clock() / CLOCKS_PER_SEC;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
}

struct TagRSearch {};
struct TagNoRSearch {};

template <bool T>
struct HasRSearchInner {
    typedef TagRSearch type;
};
template <>
struct HasRSearchInner<false> {
    typedef TagNoRSearch type;
};

template <class T>
class HasRSearch {
   private:
    typedef char One[1];
    typedef char Two[2];

    template <class T2, T2>
    class Check {};

    template <typename T3>
    static One &Func(Check<double (T3::*)(void), &T3::getRSearch> *);
    template <typename T3>
    static One &Func(Check<float (T3::*)(void), &T3::getRSearch> *);
    template <typename T3>
    static One &Func(Check<double (T3::*)(void) const, &T3::getRSearch> *);
    template <typename T3>
    static One &Func(Check<float (T3::*)(void) const, &T3::getRSearch> *);

    template <typename T3>
    static Two &Func(...);

   public:
    typedef typename HasRSearchInner<sizeof(Func<T>(NULL)) == 1>::type type;
};

template <typename Tptcl>
struct HasgetRSearchMethod {
#ifndef PARTICLE_SIMULATOR_ALL_64BIT_PRECISION
    template <typename U, F32 (U::*)()>
    struct SFINAE0 {};
    template <typename U, F32 (U::*)() const>
    struct SFINAE1 {};
#endif
    template <typename U, F64 (U::*)()>
    struct SFINAE2 {};
    template <typename U, F64 (U::*)() const>
    struct SFINAE3 {};
#ifndef PARTICLE_SIMULATOR_ALL_64BIT_PRECISION
    template <typename U>
    static char Test(SFINAE0<U, &U::getRSearch> *);
    template <typename U>
    static char Test(SFINAE1<U, &U::getRSearch> *);
#endif
    template <typename U>
    static char Test(SFINAE2<U, &U::getRSearch> *);
    template <typename U>
    static char Test(SFINAE3<U, &U::getRSearch> *);
    template <typename U>
    static int Test(...);
    static const bool value = sizeof(Test<Tptcl>(0)) == sizeof(char);
};
template <class Tptcl>
F64 GetMyRSearch(Tptcl ptcl, std::true_type) {
    return ptcl.getRSearch();
}
template <class Tptcl>
F64 GetMyRSearch(Tptcl ptcl, std::false_type) {
    return 0.0;
}
template <class Tptcl>
F64 GetMyRSearch(Tptcl ptcl) {
    return GetMyRSearch(ptcl, std::integral_constant<bool, HasgetRSearchMethod<Tptcl>::value>());
}

// (4) getDipole
template <typename Tptcl>
struct HasgetDipoleMethod {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    template <typename U, Vector2<float> (U::*)()>
    struct SFINAE0 {};
    template <typename U, Vector2<float> (U::*)() const>
    struct SFINAE1 {};
    template <typename U, Vector2<double> (U::*)()>
    struct SFINAE2 {};
    template <typename U, Vector2<double> (U::*)() const>
    struct SFINAE3 {};
#else
    template <typename U, Vector3<float> (U::*)()>
    struct SFINAE0 {};
    template <typename U, Vector3<float> (U::*)() const>
    struct SFINAE1 {};
    template <typename U, Vector3<double> (U::*)()>
    struct SFINAE2 {};
    template <typename U, Vector3<double> (U::*)() const>
    struct SFINAE3 {};
#endif
    template <typename U>
    static char Test(SFINAE0<U, &U::getDipole> *);
    template <typename U>
    static char Test(SFINAE1<U, &U::getDipole> *);
    template <typename U>
    static char Test(SFINAE2<U, &U::getDipole> *);
    template <typename U>
    static char Test(SFINAE3<U, &U::getDipole> *);
    template <typename U>
    static int Test(...);
    static const bool value = sizeof(Test<Tptcl>(0)) == sizeof(char);
};
template <class Tptcl>
F64vec GetMyDipole(Tptcl ptcl, std::true_type) {
    return ptcl.getDipole();
}
template <class Tptcl>
F64vec GetMyDipole(Tptcl ptcl, std::false_type) {
    return F64vec(0.0);
}
template <class Tptcl>
F64vec GetMyDipole(Tptcl ptcl) {
    return GetMyDipole(ptcl, std::integral_constant<bool, HasgetDipoleMethod<Tptcl>::value>());
}
// (5) getQuadrupole
template <typename Tptcl>
struct HasgetQuadrupoleMethod {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    template <typename U, MatrixSym2<float> (U::*)()>
    struct SFINAE0 {};
    template <typename U, MatrixSym2<float> (U::*)() const>
    struct SFINAE1 {};
    template <typename U, MatrixSym2<double> (U::*)()>
    struct SFINAE2 {};
    template <typename U, MatrixSym2<double> (U::*)() const>
    struct SFINAE3 {};
#else
    template <typename U, MatrixSym3<float> (U::*)()>
    struct SFINAE0 {};
    template <typename U, MatrixSym3<float> (U::*)() const>
    struct SFINAE1 {};
    template <typename U, MatrixSym3<double> (U::*)()>
    struct SFINAE2 {};
    template <typename U, MatrixSym3<double> (U::*)() const>
    struct SFINAE3 {};
#endif
    template <typename U>
    static char Test(SFINAE0<U, &U::getQuadrupole> *);
    template <typename U>
    static char Test(SFINAE1<U, &U::getQuadrupole> *);
    template <typename U>
    static char Test(SFINAE2<U, &U::getQuadrupole> *);
    template <typename U>
    static char Test(SFINAE3<U, &U::getQuadrupole> *);
    template <typename U>
    static int Test(...);
    static const bool value = sizeof(Test<Tptcl>(0)) == sizeof(char);
};
template <class Tptcl>
F64mat GetMyQuadrupole(Tptcl ptcl, std::true_type) {
    return ptcl.getQuadrupole();
}
template <class Tptcl>
F64mat GetMyQuadrupole(Tptcl ptcl, std::false_type) {
    return F64mat(0.0);
}
template <class Tptcl>
F64mat GetMyQuadrupole(Tptcl ptcl) {
    return GetMyQuadrupole(ptcl, std::integral_constant<bool, HasgetQuadrupoleMethod<Tptcl>::value>());
}

inline std::string GetBinString(const U32 val) {
    if (!val) return std::string("00000000000000000000000000000000");
    U32 tmp = val;
    std::string str;
    for (S32 i = 0; i < 32; i++) {
        if ((tmp & 1) == 0)
            str.insert(str.begin(), '0');
        else
            str.insert(str.begin(), '1');
        tmp >>= 1;
    }
    return str;
}
inline std::string GetBinString(const U32 *val) {
    if (!*val) return std::string("00000000000000000000000000000000");
    U32 tmp = *val;
    std::string str;
    for (S32 i = 0; i < 32; i++) {
        if ((tmp & 1) == 0)
            str.insert(str.begin(), '0');
        else
            str.insert(str.begin(), '1');
        tmp >>= 1;
    }
    return str;
}
inline std::string GetBinString(const U64 val) {
    if (!val) return std::string("0000000000000000000000000000000000000000000000000000000000000000");
    U64 tmp = val;
    std::string str;
    for (S32 i = 0; i < 64; i++) {
        if ((tmp & 1) == 0)
            str.insert(str.begin(), '0');
        else
            str.insert(str.begin(), '1');
        tmp >>= 1;
    }
    return str;
}
inline std::string GetBinString(const U64 *val) {
    if (!*val) return std::string("0000000000000000000000000000000000000000000000000000000000000000");
    U64 tmp = *val;
    std::string str;
    for (S32 i = 0; i < 64; i++) {
        if ((tmp & 1) == 0)
            str.insert(str.begin(), '0');
        else
            str.insert(str.begin(), '1');
        tmp >>= 1;
    }
    return str;
}

#ifdef TEST_VARIADIC_TEMPLATE
template <class T>
class IsParticleSystem {
    template <class T2>
    static auto check(T2 x) -> decltype(x.ParticleSystemDummyFunc(), std::true_type());
    static auto check(...) -> decltype(std::false_type());

   public:
    typedef decltype(check(std::declval<T>())) value;
};
template <class T>
class IsDomainInfo {
    template <class T2>
    static auto check(T2 x) -> decltype(x.DomainInfoDummyFunc(), std::true_type());
    static auto check(...) -> decltype(std::false_type());

   public:
    typedef decltype(check(std::declval<T>())) value;
};
#endif

#ifdef LOOP_TREE
static const S32 ADR_TREE_CELL_NULL = 1;
static const S32 SIMD_VEC_LEN = 8;
static const S32 N_THREAD_LIMIT = 64;
static const S32 BUF_SIZE_INTERACTION_LIST = 100000;
#endif

template <typename T>
struct my_is_char : std::false_type {};
template <>
struct my_is_char<char *> : std::true_type {};
template <>
struct my_is_char<const char *> : std::true_type {};
template <>
struct my_is_char<const char *const> : std::true_type {};

template <typename T>
struct my_is_string : std::false_type {};
template <>
struct my_is_string<std::string> : std::true_type {};
template <>
struct my_is_string<const std::string> : std::true_type {};

}  // namespace ParticleSimulator
// #include"comm_info.hpp"
#include "key_type.hpp"
#include "time_profile.hpp"
#include "utils.hpp"
#include "timer.hpp"
#include "utils_pmmm.hpp"
//#include "particle_mesh_multipole/fmm.hpp"
//#include "particle_mesh_multipole/cutoff.hpp"
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
#include "morton_key2.hpp"
#else
#include "morton_key3.hpp"
#endif
