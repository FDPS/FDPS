#pragma once

#include<cmath>
#include<vector>
#include<functional>
#include<algorithm>
#include<exception>
#include<stdexcept>
#include<cassert>
#include<typeinfo>
#include<cstdio>
#include<cstring>
#include<map>
#include<random>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include"mpi.h"
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include<omp.h>
#define PS_OMP(...) _Pragma(#__VA_ARGS__)
#define PS_OMP_PARALLEL_FOR _Pragma("omp parallel for")
#define PS_OMP_PARALLEL _Pragma("omp parallel")
#define PS_OMP_FOR _Pragma("omp for")
#define PS_OMP_CRITICAL _Pragma("omp critical")
#define PS_OMP_BARRIER _Pragma("omp barrier")
#define PS_OMP_SINGLE _Pragma("omp single")
#else
#define PS_OMP(...)
#define PS_OMP_PARALLEL_FOR
#define PS_OMP_PARALLEL
#define PS_OMP_FOR
#define PS_OMP_CRITICAL
#define PS_OMP_BARRIER
#define PS_OMP_SINGLE
#endif

#ifdef __GNUC__
#define PS_INLINE inline __attribute__((always_inline))
#else
#define PS_INLINE inline
#endif

#define PS_SET_CALLER_INFO(obj)					\
    do {							\
	obj.setCallerInfo(__LINE__, __FUNCTION__, __FILE__);	\
    } while(0);							\

#include"memory_pool.hpp"
#include"vector2.hpp"
#include"vector3.hpp"
#include"orthotope2.hpp"
#include"orthotope3.hpp"
#include"orthotope2i.hpp"
#include"orthotope3i.hpp"
#include"matrix_sym2.hpp"
#include"matrix_sym3.hpp"
#include"matrix2.hpp"

#define PS_DEBUG_CALL(func) \
    do {                                              \
        try{                                          \
           func;                                      \
           std::cout << "[FDPS msg] "                 \
                     << #func << ": "                 \
                     << getMemSizeUsed() << ", "      \
                     << Comm::getRank()  << ", "      \
                     << typeid(TSM).name() << ". "    \
                     << std::endl;                    \
        } catch (std::bad_alloc& e) {                 \
           std::cout << "[FDPS error] "               \
                     << #func << ": "                 \
                     << getMemSizeUsed() << ", "      \
                     << Comm::getRank()  << ", "      \
                     << typeid(TSM).name() << ". "    \
                     << std::endl;                    \
           MPI_Abort(MPI_COMM_WORLD,9);		      \
           std::exit(1);                              \
        } catch (...) {                               \
           std::cout << "[FDPS unknown error] "       \
                     << #func << ": "                 \
                     << getMemSizeUsed() << ", "      \
                     << Comm::getRank()  << ", "      \
                     << typeid(TSM).name() << ". "    \
                     << std::endl;                    \
           MPI_Abort(MPI_COMM_WORLD,9);		      \
           std::exit(1);                              \
        }                                             \
        MPI_Barrier(MPI_COMM_WORLD);		      \
        if (Comm::getRank() == 0)                     \
           std::cout << #func                         \
                     << " passed." << std::endl;      \
    } while(0);


#define PARTICLE_SIMULATOR_PRINT_ERROR(msg) \
    { std::cout<<"PS_ERROR: "<<msg<<" \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl; }

#define PARTICLE_SIMULATOR_PRINT_LINE_INFO() \
    { std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl; }

template<typename Tcomm>
void DEBUG_PRINT_MAKING_TREE(Tcomm & comm_info){
#if defined(DEBUG_PRINT_MAKING_TREE)
    comm_info.barrier();
    if(comm_info.getRank()==0)
	PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
}

namespace ParticleSimulator{
    static const long long int LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1ll<<30;
    static inline void Abort(const int err = -1){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
      MPI_Abort(MPI_COMM_WORLD, err);
#else
        exit(err);
#endif
    }
}

#include"reallocatable_array.hpp"

namespace ParticleSimulator{
    typedef int S32;
    typedef unsigned int U32;
#ifdef PARTICLE_SIMULATOR_ALL_64BIT_PRECISION
    typedef double F32;
#else
    typedef float F32;
#endif
    typedef long long int S64;
    typedef unsigned long long int U64;
    typedef double F64;
    typedef Vector2<S32> S32vec2;
    typedef Vector3<S32> S32vec3;
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

    static const S32 DIMENSION_LIMIT = 3;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    typedef S32vec2 S32vec;
    typedef F32vec2 F32vec;
    typedef F64vec2 F64vec;
    typedef F32mat2 F32mat;
    typedef F64mat2 F64mat;
    typedef S32ort2 S32ort;
    typedef F32ort2 F32ort;
    typedef F64ort2 F64ort;
    static const S32 DIMENSION = 2;
    static const S32 N_CHILDREN = 4;
    static const S32 TREE_LEVEL_LIMIT = 30;
    static const F64vec SHIFT_CENTER[N_CHILDREN] = 
        { F64vec(-0.5, -0.5), F64vec(-0.5, 0.5),
          F64vec( 0.5, -0.5), F64vec( 0.5, 0.5) };
#else
    typedef S32vec3 S32vec;
    typedef F32vec3 F32vec;
    typedef F64vec3 F64vec;
    typedef F32mat3 F32mat;
    typedef F64mat3 F64mat;
    typedef S32ort3 S32ort;
    typedef F32ort3 F32ort;
    typedef F64ort3 F64ort;
    static const S32 DIMENSION = 3;
    static const S32 N_CHILDREN = 8;
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
    static const S32 TREE_LEVEL_LIMIT = 21;
#elif defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
    static const S32 TREE_LEVEL_LIMIT = 31;
    static const S32 KEY_LEVEL_MAX_HI = 21;
    static const S32 KEY_LEVEL_MAX_LO = 10;
#else
    static const S32 TREE_LEVEL_LIMIT = 42;
    static const S32 KEY_LEVEL_MAX_HI = 21;
    static const S32 KEY_LEVEL_MAX_LO = 21;
#endif
    static const F64vec SHIFT_CENTER[N_CHILDREN] = 
        { F64vec(-0.5, -0.5, -0.5), F64vec(-0.5, -0.5, 0.5),
          F64vec(-0.5,  0.5, -0.5), F64vec(-0.5,  0.5,  0.5),
          F64vec( 0.5, -0.5, -0.5), F64vec( 0.5, -0.5,  0.5),
          F64vec( 0.5,  0.5, -0.5), F64vec( 0.5,  0.5,  0.5) };
#endif
    
#ifdef PARTICLE_SIMULATOR_SPMOM_F32
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
    static const F64 LARGE_DOUBLE = std::numeric_limits<F64>::max()*0.0625;
    static const F64 LARGE_FLOAT  = std::numeric_limits<F32>::max()*0.0625;
    static const S64 LARGE_INT    = std::numeric_limits<S32>::max()*0.0625;

    ///// A.Tanikawa modified from
    // In the upper line, the right-hand side is interpreted to be a 32bit-integer.
    // static const S64 LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1<<31;
    // static const S64 LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1ll<<30;
    ///// A.Tanikawa modified to

    //////////////////
    /// enum
    enum EXCHANGE_LET_MODE{
        EXCHANGE_LET_A2A,
        EXCHANGE_LET_P2P_EXACT,
        EXCHANGE_LET_P2P_FAST,
    };
    
    enum INTERACTION_LIST_MODE{
        MAKE_LIST,
        MAKE_LIST_FOR_REUSE,
        REUSE_LIST,
    };
    
    enum SEARCH_MODE{
        LONG_NO_CUTOFF,
        LONG_CUTOFF,
        LONG_SCATTER, // new for P^3T
        LONG_CUTOFF_SCATTER, // new for P^3T + PM
        LONG_SYMMETRY,
        SHORT_GATHER,
        SHORT_SCATTER,
        SHORT_SYMMETRY,
    };
    
    enum FORCE_TYPE{
        FORCE_TYPE_LONG,
        FORCE_TYPE_SHORT,
    };

    enum CALC_DISTANCE_TYPE{
        CALC_DISTANCE_TYPE_NORMAL = 0,
        CALC_DISTANCE_TYPE_NEAREST_X = 1,
	CALC_DISTANCE_TYPE_NEAREST_Y = 2,
        CALC_DISTANCE_TYPE_NEAREST_XY = 3,
	CALC_DISTANCE_TYPE_NEAREST_Z = 4,
	CALC_DISTANCE_TYPE_NEAREST_XZ = 5,
	CALC_DISTANCE_TYPE_NEAREST_YZ = 6,
	CALC_DISTANCE_TYPE_NEAREST_XYZ = 7,
    };

    struct TagForceLong{
        enum{
            force_type = FORCE_TYPE_LONG,
        };
    };
    struct TagForceShort{
        enum{
            force_type = FORCE_TYPE_SHORT,
        };
    };

    // GEOMETRY OF TREE CELL
    struct TagGeometrySize{};
    struct TagGeometrySizeOut{};
    struct TagGeometrySizeInOut{};
    struct TagGeometryIn{};
    struct TagGeometryOut{};
    struct TagGeometryInAndOut{};

    // LOCAL TREE MOMENT TYPE
    struct TagCalcMomLongEpjLt{};
    struct TagCalcMomShortEpjLt{};
    struct TagCalcMomShortEpiLt{};
    
    // IS PERIODIC BOUNDARY CONDITION AVAILABLE
    struct TagSearchBoundaryConditionOpenOnly{};
    struct TagSearchBoundaryConditionOpenPeriodic{};

    // SEARCH TYPE
    struct TagSearchLong{};
    struct TagSearchLongCutoff{};
    struct TagSearchLongScatter{};
    struct TagSearchLongSymmetry{};
    struct TagSearchLongCutoffScatter{};
    struct TagSearchShortGather{};
    struct TagSearchShortScatter{};
    struct TagSearchShortSymmetry{};

    // IP-GROUP TYPE
    struct TagIpgLongNormal{};
    struct TagIpgIn{};
    struct TagIpgInAndOut{};
    struct TagIpgOut{};

    struct TagWithoutCutoff{};
    struct TagWithCutoff{};

    // NEIGHBOR SEARCH TYPE
    struct TagNeighborSearchSymmetry{};
    struct TagNeighborSearchGather{};
    struct TagNeighborSearchScatter{};
    struct TagNeighborSearchNo{};

    struct SEARCH_MODE_LONG{
        typedef TagForceLong force_type;
        typedef TagSearchLong search_type;
        typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
        typedef TagIpgLongNormal ipg_type;
        typedef TagNeighborSearchNo neighbor_search_type;
        using tree_cell_loc_geometry_type = TagGeometrySize;
        using tree_cell_glb_geometry_type = TagGeometrySize;
        using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
      //enum{
      //    search_type_id = LONG_NO_CUTOFF,
      //};
    };
    
    struct SEARCH_MODE_LONG_CUTOFF{
        typedef TagForceLong force_type;
        typedef TagSearchLongCutoff search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgLongNormal ipg_type;
        typedef TagNeighborSearchNo neighbor_search_type;
        using tree_cell_loc_geometry_type = TagGeometrySizeOut;
        using tree_cell_glb_geometry_type = TagGeometrySizeOut;
        using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
      //enum{
      //    search_type_id = LONG_CUTOFF,
      //};
    };
    struct SEARCH_MODE_LONG_SCATTER{
        typedef TagForceLong force_type;
        typedef TagSearchLongScatter search_type;
        typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
        typedef TagIpgIn ipg_type;
        typedef TagNeighborSearchScatter neighbor_search_type;
        using tree_cell_loc_geometry_type = TagGeometrySizeInOut;
        using tree_cell_glb_geometry_type = TagGeometrySizeInOut;
        using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
      //enum{
      //    search_type_id = LONG_SCATTER,
      //};
    };
    struct SEARCH_MODE_LONG_SYMMETRY{
        typedef TagForceLong force_type;
        typedef TagSearchLongSymmetry search_type;
        typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
        typedef TagIpgInAndOut ipg_type;
        typedef TagNeighborSearchSymmetry neighbor_search_type;
        using tree_cell_loc_geometry_type = TagGeometrySizeInOut;
        using tree_cell_glb_geometry_type = TagGeometrySizeInOut;
        using calc_moment_local_tree_type = TagCalcMomLongEpjLt;
      //enum{
      //    search_type_id = LONG_SYMMETRY,
      //};
    };
    
    struct SEARCH_MODE_GATHER{
        typedef TagForceShort force_type;
        typedef TagSearchShortGather search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgOut ipg_type;
        typedef TagNeighborSearchGather neighbor_search_type;
        using tree_cell_loc_geometry_type = TagGeometryInAndOut;
        using tree_cell_glb_geometry_type = TagGeometryIn;
        using calc_moment_local_tree_type = TagCalcMomShortEpiLt;
      //enum{
      //    search_type_id = SHORT_GATHER,
      //};
    };
    struct SEARCH_MODE_SCATTER{
        typedef TagForceShort force_type;
        typedef TagSearchShortScatter search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgIn ipg_type;
        typedef TagNeighborSearchScatter neighbor_search_type;
        using tree_cell_loc_geometry_type = TagGeometryInAndOut;
        using tree_cell_glb_geometry_type = TagGeometryInAndOut;
        using calc_moment_local_tree_type = TagCalcMomShortEpjLt;
      //enum{
      //    search_type_id = SHORT_SCATTER,
      //};
    };
    struct SEARCH_MODE_SYMMETRY{
        typedef TagForceShort force_type;
        typedef TagSearchShortSymmetry search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgInAndOut ipg_type;
        typedef TagNeighborSearchSymmetry neighbor_search_type;
        using tree_cell_loc_geometry_type = TagGeometryInAndOut;
        using tree_cell_glb_geometry_type = TagGeometryInAndOut;
        using calc_moment_local_tree_type = TagCalcMomShortEpjLt;
      //enum{
      //    search_type_id = SHORT_SYMMETRY,
      //};
    };

    template<class T> class ValueTypeReduction;
    template<> class ValueTypeReduction<float>{};
    template<> class ValueTypeReduction<double>{};
    template<> class ValueTypeReduction<int>{};
    template<> class ValueTypeReduction<long>{};

    enum BOUNDARY_CONDITION{
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

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    template<class T> 
    inline MPI_Datatype GetDataType(){
      static MPI_Datatype type = MPI_DATATYPE_NULL;
      if( type == MPI_DATATYPE_NULL ){
          MPI_Type_contiguous(sizeof(T),MPI_BYTE,&type);
          MPI_Type_commit(&type);
      }
      return type;
    }
    template<class T> 
        inline MPI_Datatype GetDataType(const T &){
        return GetDataType<T>();
    }
    template<> inline MPI_Datatype GetDataType<int>(){return MPI_INT;}
    template<> inline MPI_Datatype GetDataType<long>(){return MPI_LONG;}
    template<> inline MPI_Datatype GetDataType<long long int>(){return MPI_LONG_LONG_INT;}
    template<> inline MPI_Datatype GetDataType<unsigned int>(){return MPI_UNSIGNED;}
    template<> inline MPI_Datatype GetDataType<unsigned long>(){return MPI_UNSIGNED_LONG;}
    template<> inline MPI_Datatype GetDataType<unsigned long long int>(){return MPI_UNSIGNED_LONG_LONG;}
    template<> inline MPI_Datatype GetDataType<float>(){return MPI_FLOAT;}
    template<> inline MPI_Datatype GetDataType<double>(){return MPI_DOUBLE;}

    template<class Tfloat, class Tint> 
    inline MPI_Datatype GetDataType();
    template<> inline MPI_Datatype GetDataType<int, int>(){return MPI_2INT;}
    template<> inline MPI_Datatype GetDataType<long, int>(){return MPI_LONG_INT;}
    template<> inline MPI_Datatype GetDataType<float, int>(){return MPI_FLOAT_INT;}
    template<> inline MPI_Datatype GetDataType<double, int>(){return MPI_DOUBLE_INT;}
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    template<class T, int DIM_COMM=2>
    class CommForAllToAll{
    private:
        int rank_glb_;
        int n_proc_glb_;
        MPI_Comm comm_glb_;
        int rank_1d_[DIM_COMM];
        int n_proc_1d_[DIM_COMM];
        MPI_Comm comm_1d_[DIM_COMM];
        ReallocatableArray<T> val_send_glb_;
        int * n_recv_disp_glb_;
        int * n_send_disp_glb_;
        int * n_send_glb_;
        int * n_recv_1d_;
        int * n_send_1d_;
        int * n_recv_disp_1d_;
        int * n_send_disp_1d_;

        template<int DIM>
        void divideProc(int np[],
                        int rank[],
                        const int nproc,
                        const int myrank){
            std::vector<int> npv;
            npv.resize(DIM);
            int np_tmp = nproc;
            for(int d=DIM, cid=0; cid<DIM-1; d--, cid++){
                int tmp = (int)pow(np_tmp+0.000001, (1.0/d)*1.000001 );
                while(np_tmp%tmp){
                    tmp--;
                }
                npv[cid] = tmp;
                np_tmp /= npv[cid];
            }
            npv[DIM-1] = np_tmp;
            int rank_tmp = myrank;
            std::sort(npv.begin(), npv.end(), std::greater<int>());
            for(int i=DIM-1; i>=0; i--){
                np[i] = npv[i];
                rank[i] = rank_tmp % np[i];
                rank_tmp /= np[i];
            }
        }

    public:
        void dumpRank(){
            for(int i=0; i<DIM_COMM; i++){
                int n_tmp = 0;
                MPI_Comm_size(comm_1d_[i], &n_tmp);
                std::cout<<"n_proc_1d_["<<i<<"]="<<n_proc_1d_[i]<<" n_tmp="<<n_tmp<<std::endl;
            }
            for(int i=0; i<DIM_COMM; i++){
                int r_tmp = 0;
                MPI_Comm_rank(comm_1d_[i], &r_tmp);
                std::cout<<"rank_1d_["<<i<<"]="<<rank_1d_[i]<<" r_tmp="<<r_tmp<<std::endl;
            }
        }

        CommForAllToAll(MPI_Comm comm = MPI_COMM_WORLD){
            comm_glb_ = comm;
            MPI_Comm_rank(comm_glb_, &rank_glb_);
            MPI_Comm_size(comm_glb_, &n_proc_glb_);
            n_recv_disp_glb_ = new int[n_proc_glb_ + 1];
            n_send_disp_glb_ = new int[n_proc_glb_ + 1];
            n_send_glb_ = new int[n_proc_glb_];
            n_recv_1d_ = new int[n_proc_glb_];
            n_send_1d_ = new int[n_proc_glb_];
            n_recv_disp_1d_ = new int[n_proc_glb_ + 1];
            n_send_disp_1d_ = new int[n_proc_glb_ + 1];
            divideProc<DIM_COMM>(n_proc_1d_, rank_1d_, n_proc_glb_, rank_glb_);
            int dim_max = -1;
            for(int i=0; i<DIM_COMM; i++){
                if(dim_max < n_proc_1d_[i]){
                    dim_max = n_proc_1d_[i];
                }
            }
            int split_color = 0;
            int factor = 1;
            for(int d=0; d<DIM_COMM; d++){
                split_color += rank_1d_[d] * factor;
                factor *= dim_max;
            }
            for(int d=DIM_COMM-1; d>=0; d--){
                factor = rank_1d_[d];
                for(int d0=0; d0<d; d0++){
                    factor *= dim_max;
                }
                MPI_Comm_split(comm_glb_, split_color-factor, rank_glb_, comm_1d_+d);
            }
        } // Constructor
	
        // alltoall
        void execute(const ReallocatableArray<T> & val_send,
                     const int cnt,
                     ReallocatableArray<T> & val_recv){
            const int n_recv_tot = cnt * n_proc_glb_;
            val_recv.resizeNoInitialize( n_recv_tot );
            val_send_glb_.resizeNoInitialize( n_recv_tot );
            for(int i=0; i<n_recv_tot; i++){
                val_recv[i] = val_send[i];
            }
            for(int d=DIM_COMM-1; d>=0; d--){
                const int radix = n_proc_glb_ / n_proc_1d_[d];
                for(int ib=0; ib<n_proc_glb_; ib++){
                    const int id_send = cnt * ( (ib % radix) * n_proc_1d_[d] + ib / radix );
                    const int offset = ib * cnt;
                    for(int i=0; i<cnt; i++){
                        val_send_glb_[offset + i] = val_recv[id_send + i];
                    }
                }
                MPI_Alltoall(val_send_glb_.getPointer(), cnt*radix, GetDataType<T>(),
                             val_recv.getPointer(), cnt*radix, GetDataType<T>(), comm_1d_[d]);
            }
        }

        void execute(const T val_send[],
                     const int cnt,
                     T val_recv[]){
            const int n_recv_tot = cnt * n_proc_glb_;
            val_send_glb_.resizeNoInitialize( n_recv_tot );
            for(int i=0; i<n_recv_tot; i++){
                val_recv[i] = val_send[i];
            }
            for(int d=DIM_COMM-1; d>=0; d--){
                const int radix = n_proc_glb_ / n_proc_1d_[d];
                for(int ib=0; ib<n_proc_glb_; ib++){
                    const int id_send = cnt * ( (ib % radix) * n_proc_1d_[d] + ib / radix );
                    const int offset = ib * cnt;
                    for(int i=0; i<cnt; i++){
                        val_send_glb_[offset + i] = val_recv[id_send + i];
                    }
                }
                MPI_Alltoall(val_send_glb_.getPointer(), cnt*radix, GetDataType<T>(),
                             val_recv, cnt*radix, GetDataType<T>(), comm_1d_[d]);
            }
        }

        void executeV(const ReallocatableArray<T> & val_send,
                      ReallocatableArray<T> & val_recv,
                      const int n_send[],
                      int n_recv[],
                      const int n_send_offset=0,
                      const int n_recv_offset=0){
            int cnt = 0;
            val_recv.reserveAtLeast(n_recv_offset+val_send.size()-n_send_offset);
            val_recv.resizeNoInitialize(n_recv_offset);
            for(int ib=0; ib<n_proc_glb_; ib++){
                n_recv[ib] = n_send[ib];
                for(int ip=0; ip<n_recv[ib]; ip++, cnt++){
                    val_recv.pushBackNoCheck(val_send[cnt+n_send_offset]);
                }
            }
            for(int d=DIM_COMM-1; d>=0; d--){
                int radix = n_proc_glb_ / n_proc_1d_[d];
                n_recv_disp_glb_[0] = 0;
                for(int i=0; i<n_proc_glb_; i++){
                    n_recv_disp_glb_[i+1] = n_recv_disp_glb_[i] + n_recv[i];
                }
                val_send_glb_.reserveAtLeast( val_recv.size()-n_recv_offset );
                val_send_glb_.clearSize();
                for(int ib0=0; ib0<n_proc_glb_; ib0++){
                    int id_send = (ib0 % radix) * n_proc_1d_[d] + ib0 / radix;
                    n_send_glb_[ib0] = n_recv[id_send];
                    int offset = n_recv_disp_glb_[id_send];
                    for(int ib1=0; ib1<n_send_glb_[ib0]; ib1++){
                        val_send_glb_.pushBackNoCheck(val_recv[ib1 + offset + n_recv_offset]);
                    }
                }
                MPI_Alltoall(n_send_glb_, radix, MPI_INT,
                             n_recv, radix, MPI_INT, comm_1d_[d]);
                n_send_disp_1d_[0] = n_recv_disp_1d_[0] = 0;
                for(int ib0=0; ib0<n_proc_1d_[d]; ib0++){
                    n_send_1d_[ib0] = n_recv_1d_[ib0] = 0;
                    int offset = ib0 * radix;
                    for(int ib1=0; ib1<radix; ib1++){
                        n_send_1d_[ib0] += n_send_glb_[offset + ib1];
                        n_recv_1d_[ib0] += n_recv[offset + ib1];
                    }
                    n_send_disp_1d_[ib0+1] = n_send_disp_1d_[ib0] + n_send_1d_[ib0];
                    n_recv_disp_1d_[ib0+1] = n_recv_disp_1d_[ib0] + n_recv_1d_[ib0];
                }
                val_recv.resizeNoInitialize(n_recv_disp_1d_[n_proc_1d_[d]] + n_recv_offset);
                MPI_Alltoallv(val_send_glb_.getPointer(), n_send_1d_, n_send_disp_1d_, GetDataType<T>(),
                              val_recv.getPointer(n_recv_offset), n_recv_1d_, n_recv_disp_1d_, GetDataType<T>(), comm_1d_[d]);
            }
        }
    }; //CommForAllToAll
#endif


    template <int DIM>
    inline void SetNumberOfDomainMultiDimension(const int n_proc, const int my_rank,
                                                int np[], int rank[]){
        for (int i=0; i<DIMENSION_LIMIT; i++){
            np[i] = 1;
            rank[i] = 1;
        }
        std::vector<S32> npv;
        npv.resize(DIM);
        auto np_tmp = n_proc;
        for(int d=DIM, cid=0; cid<DIM-1; d--, cid++){
            int tmp = (S32)pow(np_tmp+0.000001, (1.0/d)*1.000001 );
            while(np_tmp%tmp){
                tmp--;
            }
            npv[cid] = tmp;
            np_tmp /= npv[cid];
        }
        npv[DIM-1] = np_tmp;
        int rank_tmp = my_rank;
        std::sort(npv.begin(), npv.end(), std::greater<S32>());
        for(int i=DIM-1; i>=0; i--){
            np[i] = npv[i];
            rank[i] = rank_tmp % np[i];
            rank_tmp /= np[i];
        }
    }

    class CommInfo{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Comm comm_;
#endif
        int rank_;
        int n_proc_;
        S32 rank_multi_dim_[DIMENSION];
        S32 n_proc_multi_dim_[DIMENSION];
        template<class T>
        T allreduceMin(const T & val) const {
            T ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            T val_tmp = val;
            MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_MIN, comm_);
#endif
            return ret;
        }
        template<class T>
        Vector2<T> allreduceMin(const Vector2<T> & val) const {
            Vector2<T> ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            Vector2<T> val_tmp = val;
            MPI_Allreduce((T*)&val_tmp.x, (T*)&ret.x, 2, GetDataType<T>(), MPI_MIN, comm_);
#endif
            return ret;
        }
        template<class T>
        Vector3<T> allreduceMin(const Vector3<T> & val) const {
            Vector3<T> ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            Vector3<T> val_tmp = val;
            MPI_Allreduce((T*)&val_tmp.x, (T*)&ret.x, 3, GetDataType<T>(), MPI_MIN, comm_);
#endif
            return ret;
        }
        
        template<class T>
        T allreduceMax(const T & val) const {
            T ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            T val_tmp = val;
            MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_MAX, comm_);
#endif
            return ret;
        }
        template<class T>
        Vector2<T> allreduceMax(const Vector2<T> & val) const {
            Vector2<T> ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            Vector2<T> val_tmp = val;
            MPI_Allreduce((T*)&val_tmp.x, (T*)&ret.x, 2, GetDataType<T>(), MPI_MAX, comm_);
#endif
            return ret;
        }
        template<class T>
        Vector3<T> allreduceMax(const Vector3<T> & val) const {
            Vector3<T> ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            Vector3<T> val_tmp = val;
            MPI_Allreduce((T*)&val_tmp.x, (T*)&ret.x, 3, GetDataType<T>(), MPI_MAX, comm_);
#endif
            return ret;
        }
        template<class T>
        T allreduceSum(const T & val) const {
            T ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            T val_tmp = val;
            MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_SUM, comm_);
#endif
            return ret;
        }

        template<class T>
        Vector2<T> allreduceSum(const Vector2<T> & val) const {
            Vector2<T> ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            Vector2<T> val_tmp = val;
            MPI_Allreduce((T*)&val_tmp.x, (T*)&ret.x, 2, GetDataType<T>(), MPI_SUM, comm_);
#endif
            return ret;
        }
        template<class T>
        Vector3<T> allreduceSum(const Vector3<T> & val) const {
            Vector3<T> ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            Vector3<T> val_tmp = val;
            MPI_Allreduce((T*)&val_tmp.x, (T*)&ret.x, 3, GetDataType<T>(), MPI_SUM, comm_);
#endif
            return ret;
        }
        
        template<class Tfloat>
        void allreduceMin(const Tfloat & f_in,
                          const int & i_in,
                          Tfloat & f_out,
                          int & i_out){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            struct{
                Tfloat x;
                int y;
            } loc, glb;
            loc.x = f_in;
            loc.y = i_in;
            MPI_Allreduce(&loc, &glb, 1, GetDataType<Tfloat, int>(), MPI_MINLOC, comm_);
            f_out = glb.x;
            i_out = glb.y;
#else
            f_out = f_in;
            i_out = i_in;
#endif
        }

        template<class Tfloat>
        void allreduceMax(const Tfloat & f_in,
                          const int & i_in,
                          Tfloat & f_out,
                          int & i_out){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            struct{
                Tfloat x;
                int y;
            } loc, glb;
            loc.x = f_in;
            loc.y = i_in;
            MPI_Allreduce(&loc, &glb, 1, GetDataType<Tfloat, int>(), MPI_MAXLOC, comm_);
            f_out = glb.x;
            i_out = glb.y;
#else
            f_out = f_in;
            i_out = i_in;
#endif
        }
    public:
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        CommInfo(const MPI_Comm & c = MPI_COMM_NULL){
            setCommunicator(c);
        }
#else
        CommInfo(){
            setCommunicator();
        }
#endif
        void free(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            if(comm_ != MPI_COMM_NULL){
                MPI_Comm_free(&comm_);
            }
#endif
        }
        
        CommInfo create(const int n, const int rank[]) const {
            CommInfo comm_info_new = *this;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Group grp_old;
            MPI_Group grp_new;
            MPI_Comm comm_new;
            MPI_Comm_group(comm_, &grp_old);
            MPI_Group_incl(grp_old, n, rank, &grp_new);
            MPI_Comm_create(comm_, grp_new, &comm_new);
            comm_info_new.setCommunicator(comm_new);
#endif
            return comm_info_new;
        }

        ~CommInfo(){
            // DON'T call free in destructor, because all free function in the same communicater must be called at the same time.
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Comm getCommunicator() const {
            return comm_;
        }
        void setCommunicator(const MPI_Comm & c = MPI_COMM_WORLD){
            //std::cerr<<"A) setCommunicator"<<std::endl;
            //std::cerr<<"c= "<<c
            //         <<" MPI_COMM_WORLD= "<<MPI_COMM_WORLD
            //         <<" comm_= "<<comm_
            //         <<std::endl;
            if(c != MPI_COMM_NULL){
                //std::cerr<<"b1) setCommunicator"<<std::endl;
                MPI_Comm_dup(c, &comm_);
                MPI_Comm_rank(comm_, &rank_);
                MPI_Comm_size(comm_, &n_proc_);
            }
            else{
                //std::cerr<<"b2) setCommunicator"<<std::endl;
                comm_ = c;
                rank_ = -1;
                n_proc_ = -1;
            }
            //std::cerr<<"C) setCommunicator"<<std::endl;
        }
#else
        void setCommunicator(){
            rank_ = 0;
            n_proc_ = 1;
        }
#endif
        int getRank() const {
            //std::cerr<<"CommInfo::getRank()"<<std::endl;
            return rank_;
        }
        void setNumberOfProcMultiDim(const S32 id, const S32 n) {
            n_proc_multi_dim_[id] = n;
        }
        void setRankMultiDim(const S32 id, const S32 r) {
            rank_multi_dim_[id] = r;
        }
        S32 getRankMultiDim(const S32 id) {
            return rank_multi_dim_[id];
        }
        S32 getNumberOfProcMultiDim(const S32 id) {
            return n_proc_multi_dim_[id];
        }


        int getNumberOfProc() const {
            return n_proc_;
        }

        bool isCommNull() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            return this->comm_ == MPI_COMM_NULL;
#else
            return true;
#endif
        }
        bool isNotCommNull() const {
            return !isCommNull();
        }

        CommInfo split(int color, int key) const {
            CommInfo comm_info_new = *this;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Comm comm_new;
            MPI_Comm_split(comm_, color, key, &comm_new);
            comm_info_new.setCommunicator(comm_new);
#endif
            return comm_info_new;
        }

        const CommInfo & operator = (const CommInfo & rhs){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Comm c = rhs.getCommunicator();
            if(c != MPI_COMM_NULL){
                MPI_Comm_dup(c, &comm_);
            }
            else{
                comm_ = MPI_COMM_NULL;
            }
#endif
            rank_ = rhs.rank_;
            n_proc_ = rhs.n_proc_;
            return (*this);
        }
        
        void barrier(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Barrier(comm_);
#endif
        }
        bool synchronizeConditionalBranchAND(const bool & local){
            bool global = local;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            bool local_tmp = local;
            MPI_Allreduce(&local_tmp, &global, 1, MPI_C_BOOL, MPI_LAND, comm_);
#endif
            return global;
        }
        bool synchronizeConditionalBranchOR(const bool & local){
            bool global = local;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            bool local_tmp = local;
            MPI_Allreduce(&local_tmp, &global, 1, MPI_C_BOOL, MPI_LOR, comm_);
#endif
            return global;
        }
        
        ///////////////////////
        // MPI BCAST WRAPPER //
        template<class T>
        inline void broadcast(T * val, const int n, const int src=0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
            MPI_Bcast(val, n, GetDataType<T>(), src, comm_);
#endif
        }

        ////////////////////////
        // MPI GATHER WRAPPER //
        template<class T>
        inline void gather(T * val_send, // in
                           int n, // in
                           T * val_recv, // out
                           int dst=0 // in
                           ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Gather(val_send, n, GetDataType<T>(),
                       val_recv, n, GetDataType<T>(), dst, comm_);
#else
            for(int i=0; i<n; i++) val_recv[i] = val_send[i];
#endif
        }

        template<class T> 
        inline void gatherV(T * val_send, // in
                            int n_send,   // in
                            T * val_recv, // out
                            int * n_recv, // in
                            int * n_recv_disp, // in
                            int dst=0
                            ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Gatherv(val_send, n_send, GetDataType<T>(),
                        val_recv, n_recv, n_recv_disp, GetDataType<T>(), dst, comm_);
#else
            for(int i=0; i<n_send; i++) val_recv[i] = val_send[i];
#endif
        }

        template<class T>
        inline void gatherVAll(ReallocatableArray<T> & val_send,
                               int n_send,
                               ReallocatableArray<T> & val_recv,
                               int dst=0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            //const auto n_proc = this->getNumberOfProc();
            ReallocatableArray<int> n_recv(n_proc_, n_proc_, MemoryAllocMode::Pool);
            ReallocatableArray<int> n_disp_recv(n_proc_, n_proc_, MemoryAllocMode::Pool);
            this->gather(&n_send, 1, &n_recv[0], dst);
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc_; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            int n_recv_tot = 0;
            if(rank_ == dst){
                //if(this->getRank() == dst){
                n_recv_tot = n_disp_recv[n_proc_];
            }
            val_recv.resizeNoInitialize(n_recv_tot);
            this->gatherV(val_send.getPointer(), n_send, val_recv.getPointer(), n_recv.getPointer(), n_disp_recv.getPointer(), dst);
#else
            const int n = n_send;
            for(int i=0; i<n; i++)val_recv[i] = val_send[i];
#endif
        }


        /////////////////////////
        // MPI SCATTER WRAPPER //
        template<class T>
        inline void scatter(T * val_send,
                            int n,
                            T * val_recv,
                            int src = 0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Scatter(val_send, n, GetDataType<T>(),
                        val_recv, n, GetDataType<T>(),
                        src, comm_);
#else
            for(int i=0; i<n; i++)
                val_recv[i] = val_send[i];
#endif
        }
        template<class T>
        inline void scatterV(T * val_send,
                             int * n_send,
                             int * n_send_disp,
                             T * val_recv, // output
                             int n_recv,
                             int src = 0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Scatterv(val_send, n_send, n_send_disp, GetDataType<T>(),
                         val_recv, n_recv, GetDataType<T>(), 
                         src, comm_);
#else
            //const auto n_proc = this->getNumberOfProc();
            for(int i=0; i<n_send_disp[n_proc_]; i++) val_recv[i] = val_send[i];
#endif
        }
        template<class T>
        inline void scatterVAll(ReallocatableArray<T> & val_send,
                                int * n_send,
                                ReallocatableArray<T> & val_recv, // output
                                int src = 0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            int n_recv = 0;
            this->scatter(n_send, 1, &n_recv, src);
            //const auto n_proc = this->getNumberOfProc();
            ReallocatableArray<int> n_disp_send(n_proc_+1, n_proc_+1, MemoryAllocMode::Pool);
            n_disp_send[0] = 0;
            if(rank_ == src){
                //if(this->getRank() == src){
                for(int i=0; i<n_proc_; i++){
                    n_disp_send[i+1] = n_disp_send[i] + n_send[i];
                }
            }
            else{
                for(int i=0; i<n_proc_; i++){
                    n_send[i] = 0;
                    n_disp_send[i+1] = n_disp_send[i] + n_send[i];
                }
            }
            val_recv.resizeNoInitialize(n_recv);
            this->scatterV(val_send.getPointer(), n_send, n_disp_send.getPointer(),
                           val_recv.getPointer(), n_recv, src);
#else
            for(int i=0; i<n_send[0]; i++) val_recv[i] = val_send[i];
#endif
        }

        ///////////////////////////
        // MPI ALLGATHER WRAPPER //
        template<class T>
        inline void allGather(const T * val_send, // in
                              const int n,  // in
                              T * val_recv  // out
                              ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Allgather(val_send, n, GetDataType<T>(),
                          val_recv, n, GetDataType<T>(), 
                          comm_);
#else
            for(int i=0; i<n; i++)val_recv[i] = val_send[i];
#endif
        }



        template<class T>
        inline void allGatherV(T * val_send, // in
                               int n_send,   // in
                               T * val_recv, // out
                               int * n_recv, // in
                               int * n_recv_disp //in
                               ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Allgatherv(val_send, n_send, GetDataType<T>(),
                           val_recv, n_recv, n_recv_disp, GetDataType<T>(), 
                           comm_);
#else
            for(int i=0; i<n_send; i++) val_recv[i] = val_send[i];
            n_recv[0] = n_recv_disp[1] = n_send;
            n_recv_disp[0] = 0;
#endif
        }

        template<class T>
        inline void allGatherVAll(ReallocatableArray<T> & val_send, // in
                                  int n_send,   // in
                                  ReallocatableArray<T> & val_recv // out
                                  ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            //const auto n_proc = this->getNumberOfProc();
            ReallocatableArray<int> n_recv(n_proc_, n_proc_, MemoryAllocMode::Pool);
            this->allGather(&n_send, 1, n_recv.getPointer());
            ReallocatableArray<int> n_disp_recv(n_proc_+1, n_proc_+1, MemoryAllocMode::Pool);
            n_disp_recv[0] = 0;
            for(int i=0; i<n_proc_; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            const auto n_recv_tot = n_disp_recv[n_proc_];
            val_recv.resizeNoInitialize(n_recv_tot);
            this->allGatherV(val_send.getPointer(), n_send,
                             val_recv.getPointer(), n_recv.getPointer(),
                             n_disp_recv.getPointer());
#else
            for(int i=0; i<n_send; i++){ val_recv[i] = val_send[i]; }
#endif
        }
        
        ///////////////////////////
        // MPI ALLTOALL WRAPPER //
        template<class T>
        inline void allToAll(T * val_send, // in 
                             const int n, // in
                             T * val_recv // out
                             ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Alltoall(val_send, n, GetDataType<T>(), val_recv, n, GetDataType<T>(), comm_);
#else
            for(int i=0; i<n; i++) val_recv[i] = val_send[i];
#endif
        }


        template<class T>
        inline void allToAllV(T * val_send, // in
                              int * n_send, // in
                              int * n_send_disp, // in
                              T * val_recv, // out
                              int * n_recv, // in
                              int * n_recv_disp // in
                              ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Alltoallv(val_send, n_send, n_send_disp, GetDataType<T>(),
                          val_recv, n_recv, n_recv_disp, GetDataType<T>(),
                          comm_);
#else
            for(int i=0; i<n_send[0]; i++) val_recv[i] = val_send[i];
            n_recv[0] = n_send[0];
            n_recv_disp[0] = n_send_disp[0];
            n_recv_disp[1] = n_send_disp[1];
#endif
        }

        template<class T>
        inline void allToAllVAll(ReallocatableArray<T> & val_send, // in
                                 int * n_send, // in
                                 ReallocatableArray<T> & val_recv // out
                                 ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            //const auto n_proc = this->getNumberOfProc();
            ReallocatableArray<int> n_recv(n_proc_, n_proc_, MemoryAllocMode::Pool);
            this->allToAll(n_send, 1, n_recv.getPointer());
            ReallocatableArray<int> n_disp_recv(n_proc_+1, n_proc_+1, MemoryAllocMode::Pool);
            ReallocatableArray<int> n_disp_send(n_proc_+1, n_proc_+1, MemoryAllocMode::Pool);
            n_disp_recv[0] = n_disp_send[0] = 0;
            for(int i=0; i<n_proc_; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
                n_disp_send[i+1] = n_disp_send[i] + n_send[i];
            }
            const auto n_recv_tot = n_disp_recv[n_proc_];
            val_recv.resizeNoInitialize(n_recv_tot);
            this->allToAllV(val_send.getPointer(), n_send, n_disp_send.getPointer(), val_recv.getPointer(), n_recv.getPointer(), n_disp_recv.getPointer());
#else
            for(int i=0; i<n_send[0]; i++) val_recv[i] = val_send[i];
#endif
        }
        
        ///////////////////////////////////////
        // MPI ALLTOALL WRAPPER (send+irecv) //
        // under construction

        // val_send[n_rank_send*n_send]
        // val_recv[n_rank_recv*n_recv]
        template<class T>
        inline void sendIrecv(T * val_send, const int * rank_send, const int n_send, const int n_rank_send,
                              T * val_recv, const int * rank_recv, const int n_recv, const int n_rank_recv) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const int tag = 0;
            ReallocatableArray<MPI_Request> req_recv(n_rank_recv, n_rank_recv, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Status> status_recv(n_rank_recv, n_rank_recv, MemoryAllocMode::Pool);
            for(int i=0; i<n_rank_recv; i++){
                S32 rank = rank_recv[i];
                MPI_Irecv(&val_recv[n_recv*i], n_recv, GetDataType<T>(), rank, tag, comm_, req_recv.getPointer(i));
            }
            for(int i=0; i<n_rank_send; i++){
                S32 rank = rank_send[i];
                MPI_Send(&val_send[n_send*i], n_send, GetDataType<T>(), rank, tag, comm_);
            }
            MPI_Waitall(n_rank_recv, req_recv.getPointer(), status_recv.getPointer());
#else
            for(int i=0; i<n_send; i++) val_recv[i] = val_send[i];
#endif
        }

        template<class T>
        inline void sendIrecvV(T * val_send, const int * rank_send, const int * n_send, const int * n_disp_send, const int n_rank_send,
                               T * val_recv, const int * rank_recv, const int * n_recv, const int * n_disp_recv, const int n_rank_recv) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const int tag = 0;
            ReallocatableArray<MPI_Request> req_recv(n_rank_recv, n_rank_recv, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Status> status_recv(n_rank_recv, n_rank_recv, MemoryAllocMode::Pool);
            for(int i=0; i<n_rank_recv; i++){
                S32 rank = rank_recv[i];
                MPI_Irecv(&val_recv[n_disp_recv[i]], n_recv[i], GetDataType<T>(), rank, tag, comm_, &req_recv[i]);
            }
            for(int i=0; i<n_rank_send; i++){
                S32 rank = rank_send[i];
                MPI_Send(&val_send[n_disp_send[i]], n_send[i], GetDataType<T>(), rank, tag, comm_);
            }
            MPI_Waitall(n_rank_recv, req_recv.getPointer(), status_recv.getPointer());
#else
            for(int i=0; i<n_send[0]; i++) val_recv[i] = val_send[i];
#endif
        }
        template<class T0, class T1>
        inline void sendIrecvV(T0 * val_send_0, const int * rank_send_0, const int * n_send_0, const int * n_disp_send_0, const int n_rank_send_0,
                               T1 * val_send_1, const int * rank_send_1, const int * n_send_1, const int * n_disp_send_1, const int n_rank_send_1,
                               T0 * val_recv_0, const int * rank_recv_0, const int * n_recv_0, const int * n_disp_recv_0, const int n_rank_recv_0,
                               T1 * val_recv_1, const int * rank_recv_1, const int * n_recv_1, const int * n_disp_recv_1, const int n_rank_recv_1) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const int tag_0 = 0;
            const int tag_1 = 0;
            ReallocatableArray<MPI_Request> req_recv(n_rank_recv_0+n_rank_recv_1, n_rank_recv_0+n_rank_recv_1, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Status> status_recv(n_rank_recv_0+n_rank_recv_1, n_rank_recv_0+n_rank_recv_1, MemoryAllocMode::Pool);
            for(int i=0; i<n_rank_recv_0; i++){
                S32 rank = rank_recv_0[i];
                MPI_Irecv(&val_recv_0[n_disp_recv_0[i]], n_recv_0[i], GetDataType<T0>(), rank, tag_0, comm_, &req_recv[i]);
            }
            for(int i=0; i<n_rank_recv_1; i++){
                S32 rank = rank_recv_1[i];
                MPI_Irecv(&val_recv_1[n_disp_recv_1[i]], n_recv_1[i], GetDataType<T1>(), rank, tag_1, comm_, &req_recv[i+n_rank_recv_0]);
            }
            for(int i=0; i<n_rank_send_0; i++){
                S32 rank = rank_send_0[i];
                MPI_Send(&val_send_0[n_disp_send_0[i]], n_send_0[i], GetDataType<T0>(), rank, tag_0, comm_);
            }
            for(int i=0; i<n_rank_send_1; i++){
                S32 rank = rank_send_1[i];
                MPI_Send(&val_send_1[n_disp_send_1[i]], n_send_1[i], GetDataType<T1>(), rank, tag_1, comm_);
            }
            MPI_Waitall(n_rank_recv_0+n_rank_recv_1, req_recv.getPointer(), status_recv.getPointer());
#else
            for(int i=0; i<n_send_0[0]; i++) val_recv_0[i] = val_send_0[i];
            for(int i=0; i<n_send_1[0]; i++) val_recv_1[i] = val_send_1[i];
#endif
        }

        ///////////////////////////
        // MPI SEND/RECV WRAPPER //
        template<typename T>
        inline int send(const T * val, const int n, const int dst=0, const int tag=0) const {
            int ret = 0;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            ret = MPI_Send(val, n, GetDataType<T>(), dst, tag, comm_);
#endif
            return ret;
        }
        template<typename T>
        inline int recv(T * val, const int n, const int src=0, const int tag=0) const {
            int ret = 0;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Status * stat;
            ret = MPI_Recv(val, n, GetDataType<T>(), src, tag, comm_, stat);
#endif
            return ret;
        }
        
        ///////////////////////////
        // MPI ALLREDUCE WRAPPER //
        template<typename T>
        inline T getMinValue(const T & val) const {
            return allreduceMin(val);
        }
        template<typename T>
        inline T getMaxValue(const T & val) const {
            return allreduceMax(val);
        }
        template<class Tfloat>
        inline void getMinValue(const Tfloat & f_in, const int & i_in, Tfloat & f_out, int & i_out) {
            static_assert(std::is_floating_point<Tfloat>::value == true,
                          "Tfloat must be floating point type.");
            allreduceMin(f_in, i_in, f_out, i_out);
        }
        template<class Tfloat>
        inline void getMaxValue(const Tfloat & f_in, const int & i_in, Tfloat & f_out, int & i_out) {
            static_assert(std::is_floating_point<Tfloat>::value == true,
                          "Tfloat must be floating point type.");
            allreduceMax(f_in, i_in, f_out, i_out);
        }
        template<typename T>
        inline T getSum(const T & val) const {
            return allreduceSum(val);
        }
    };

    
    class Comm{
    private:
        Comm(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            comm_info_.setCommunicator(MPI_COMM_WORLD);
#else
            comm_info_.setCommunicator();
#endif
            rank_ = comm_info_.getRank();
            n_proc_ = comm_info_.getNumberOfProc();
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
            n_thread_ = omp_get_max_threads();
#else
            n_thread_ = 1;
#endif

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            constexpr int DIM = 2;
#else
            constexpr int DIM = 3;
#endif
            SetNumberOfDomainMultiDimension<DIM>(n_proc_, rank_,
                                                 n_proc_multi_dim_,
                                                 rank_multi_dim_);
        }
        ~Comm(){}
        Comm(const Comm & c){}
        Comm & operator=(const Comm& c);
        CommInfo comm_info_;
        S32 rank_;
        S32 n_proc_;
        S32 n_thread_;
        S32 rank_multi_dim_[DIMENSION];
        S32 n_proc_multi_dim_[DIMENSION];
        static Comm & getInstance(){
            static Comm inst;
            return inst;
        }
    public:
        static CommInfo create(const int n, const int rank[]) {
            return getInstance().comm_info_.create(n, rank);
        }
        static CommInfo getCommInfo() {
            CommInfo ret = getInstance().comm_info_;
            return ret;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        static MPI_Comm getCommunicator(){
            return getInstance().comm_info_.getCommunicator();
        }
#endif
        static void setNumberOfProcMultiDim(const S32 id, const S32 n) {
            getInstance().n_proc_multi_dim_[id] = n;
        }
        static void setRankMultiDim(const S32 id, const S32 r) {
            getInstance().rank_multi_dim_[id] = r;
        }
        static S32 getRank() {
            return getInstance().comm_info_.getRank();
        }
        static S32 getNumberOfProc() {
            return getInstance().comm_info_.getNumberOfProc();
        }
        static S32 getRankMultiDim(const S32 id) {
            return getInstance().rank_multi_dim_[id];
        }
        static S32 getNumberOfProcMultiDim(const S32 id) {
            return getInstance().n_proc_multi_dim_[id];
        }
        static CommInfo split(int color, int key){
            return getInstance().comm_info_.split(color, key);
        }

        static void barrier(){
            getInstance().comm_info_.barrier();
        }
        static bool synchronizeConditionalBranchAND(const bool & local){
            return getInstance().comm_info_.synchronizeConditionalBranchAND(local);
        }
        static bool synchronizeConditionalBranchOR(const bool & local){
            return getInstance().comm_info_.synchronizeConditionalBranchOR(local);
        }

        ///////////////////////
        // MPI BCAST WRAPPER //
        template<class T>
        static inline void broadcast(T * val,
                                     const int n,
                                     const int src=0){
            getInstance().comm_info_.broadcast(val, n, src);
        }
        
        ////////////////////////
        // MPI GATHER WRAPPER //
        template<class T>
        static inline void gather(T * val_send, // in
                           int n, // in
                           T * val_recv, // out
                           int dst=0 // in
                           ){
            getInstance().comm_info_.gather(val_send, n, val_recv, dst);
        }

        template<class T> 
        static inline void gatherV(T * val_send, // in
                            int n_send,   // in
                            T * val_recv, // out
                            int * n_recv, // in
                            int * n_disp_recv, // in
                            int dst=0
                            ){
            getInstance().comm_info_.gatherV(val_send, n_send, val_recv, n_recv, n_disp_recv, dst);
        }

        template<class T>
        static inline void gatherVAll(ReallocatableArray<T> & val_send,
                               int n_send,
                               ReallocatableArray<T> & val_recv,
                               int dst=0){
            getInstance().comm_info_.gatherVAll(val_send, n_send, val_recv, dst);
        }
        
        /////////////////////////
        // MPI SCATTER WRAPPER //
        template<class T>
        static inline void scatter(T * val_send,
                                   int n,
                                   T * val_recv,
                                   int src = 0){
            getInstance().comm_info_.scatter(val_send, n, val_recv, src);
        }
        template<class T>
        static inline void scatterV(T * val_send,
                                    int * n_send,
                                    int * n_disp_send,
                                    T * val_recv, // output
                                    int n_recv,
                                    int src = 0){
            getInstance().comm_info_.scatterV(val_send, n_send, n_disp_send, val_recv, n_recv, src);
        }
        template<class T>
        static inline void scatterVAll(ReallocatableArray<T> & val_send,
                                       int * n_send,
                                       ReallocatableArray<T> & val_recv, // output
                                       int src = 0){
            getInstance().comm_info_.scatterVAll(val_send, n_send, val_recv, src);
        }

        ///////////////////////////
        // MPI ALLGATHER WRAPPER //
        template<class T>
        static inline void allGather(const T * val_send, // in
                                     const int n,  // in
                                     T * val_recv  // out
                                     ){
            getInstance().comm_info_.allGather(val_send, n, val_recv);
        }
        template<class T>
        static inline void allGatherV(T * val_send, // in
                                      int n_send,   // in
                                      T * val_recv, // out
                                      int * n_recv, // in
                                      int * n_disp_recv //in
                                      ){
            getInstance().comm_info_.allGatherV(val_send, n_send, val_recv, n_recv, n_disp_recv);
        }
        template<class T>
        static inline void allGatherVAll(ReallocatableArray<T> & val_send, // in
                                         int n_send,   // in
                                         ReallocatableArray<T> & val_recv // out
                                         ){
            getInstance().comm_info_.allGatherVAll(val_send, n_send, val_recv);
        }
        
        ///////////////////////////
        // MPI ALLTOALL WRAPPER //
        template<class T>
        static inline void allToAll(T * val_send, // in 
                                    const int n, // in
                                    T * val_recv // out
                                    ){
            getInstance().comm_info_.allToAll(val_send, n, val_recv);
        }
        template<class T>
        static inline void allToAllV(T * val_send, // in
                                     int * n_send, // in
                                     int * n_disp_send, // in
                                     T * val_recv, // out
                                     int * n_recv, // in
                                     int * n_disp_recv // in
                                     ){
            getInstance().comm_info_.allToAllV(val_send, n_send, n_disp_send, val_recv, n_recv, n_disp_recv);
        }
        template<class T>
        static inline void allToAllVAll(ReallocatableArray<T> & val_send, // in
                                        int * n_send, // in
                                        ReallocatableArray<T> & val_recv // out
                                        ){
            getInstance().comm_info_.allToAllVAll(val_send, n_send, val_recv);
        }
        
        ///////////////////////////////////////
        // MPI ALLTOALL WRAPPER (send+irecv) //
        // under construction
        template<class T>
        static inline void sendIrecv(T * val_send, const int * rank_send, const int n_send, const int n_rank_send,
                                     T * val_recv, const int * rank_recv, const int n_recv, const int n_rank_recv){
            getInstance().comm_info_.sendIrecv(val_send, rank_send, n_send, n_rank_send,
                                               val_recv, rank_recv, n_recv, n_rank_recv);
        }

        template<class T>
        static inline void sendIrecvV(T * val_send, const int * rank_send, const int * n_send, const int * n_disp_send, const int n_rank_send,
                                      T * val_recv, const int * rank_recv, const int * n_recv, const int * n_disp_recv, const int n_rank_recv){
            getInstance().comm_info_.sendIrecvV(val_send, rank_send, n_send, n_disp_send, n_rank_send,
                                                val_recv, rank_recv, n_recv, n_disp_recv, n_rank_recv);
        }

        template<class T0, class T1>
        static inline void sendIrecvV(T0 * val_send_0, const int * rank_send_0, const int * n_send_0, const int * n_disp_send_0, const int n_rank_send_0,
                                      T1 * val_send_1, const int * rank_send_1, const int * n_send_1, const int * n_disp_send_1, const int n_rank_send_1,
                                      T0 * val_recv_0, const int * rank_recv_0, const int * n_recv_0, const int * n_disp_recv_0, const int n_rank_recv_0,
                                      T1 * val_recv_1, const int * rank_recv_1, const int * n_recv_1, const int * n_disp_recv_1, const int n_rank_recv_1){
            getInstance().comm_info_.sendIrecvV(val_send_0, rank_send_0, n_send_0, n_disp_send_0, n_rank_send_0,
                                                val_send_1, rank_send_1, n_send_1, n_disp_send_1, n_rank_send_1,
                                                val_recv_0, rank_recv_0, n_recv_0, n_disp_recv_0, n_rank_recv_0,
                                                val_recv_1, rank_recv_1, n_recv_1, n_disp_recv_1, n_rank_recv_1);
        }

        ///////////////////////////
        // MPI SEND/RECV WRAPPER //
        template<typename T>
        static inline int send(const T * val, const int n, const int dst=0, const int tag=0) {
            return getInstance().comm_info_.send(val, n, dst, tag);
        }
        template<typename T>
        static inline int recv(T * val, const int n, const int src=0, const int tag=0) {
            return getInstance().comm_info_.recv(val, n, src, tag);
        }
        
        
        
        ///////////////////////////
        // MPI ALLREDUCE WRAPPER //
        template<typename T>
        static inline T getMinValue(const T & val){
            return getInstance().comm_info_.getMinValue(val);
        }
        template<typename T>
        static inline T getMaxValue(const T & val){
            return getInstance().comm_info_.getMaxValue(val);
        }
        template<class Tfloat>
        static inline void getMinValue(const Tfloat & f_in, const int & i_in, Tfloat & f_out, int & i_out) {
            static_assert(std::is_floating_point<Tfloat>::value == true,
                          "Tfloat must be floating point type.");
            return getInstance().comm_info_.getMinValue(f_in, i_in, f_out, i_out);
        }
        template<class Tfloat>
        static inline void getMaxValue(const Tfloat & f_in, const int & i_in, Tfloat & f_out, int & i_out) {
            static_assert(std::is_floating_point<Tfloat>::value == true,
                          "Tfloat must be floating point type.");
            return getInstance().comm_info_.getMaxValue(f_in, i_in, f_out, i_out);
        }
        template<typename T>
        static inline T getSum(const T & val){
            return getInstance().comm_info_.getSum(val);
        }
        
        static S32 getNumberOfThread() { return getInstance().n_thread_; }
        static S32 getNumThreads() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
            return omp_get_num_threads();
#else
            return 1;
#endif
        }
        static S32 getThreadNum(){
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
            return omp_get_thread_num();
#else
            return 0;
#endif
        }
    }; // END OF Comm

    static inline void Initialize(int &argc,
                                  char **&argv,
                                  const S64 mpool_size=100000000){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Init(&argc, &argv);
#endif
        MemoryPool::initialize(mpool_size);
#ifdef MONAR
        bool flag_monar = false;
        bool flag_MONAR = false;
        for(S32 i=0; i<argc; i++){
            if(strcmp(argv[i], "monar") == 0) flag_monar = true;
            if(strcmp(argv[i], "MONAR") == 0) flag_MONAR = true;
        }
#endif

        if (Comm::getRank() == 0) {
            std::cerr << "     //==================================\\\\" << std::endl;
            std::cerr << "     ||                                  ||"   << std::endl;
            std::cerr << "     || ::::::: ::::::. ::::::. .::::::. ||"   << std::endl;
            std::cerr << "     || ::      ::    : ::    : ::       ||"   << std::endl;
            std::cerr << "     || ::::::  ::    : ::::::'  `:::::. ||"   << std::endl;
            std::cerr << "     || ::      ::::::' ::      `......' ||"   << std::endl;
            std::cerr << "     ||     Framework for Developing     ||"   << std::endl;
            std::cerr << "     ||        Particle Simulator        ||"   << std::endl;
            std::cerr << "     ||     Version 7.0 (2021/08)        ||"   << std::endl;
            std::cerr << "     \\\\==================================//" << std::endl;
            std::cerr << "" << std::endl;
            std::cerr << "       Home   : https://github.com/fdps/fdps " << std::endl;
            std::cerr << "       E-mail : fdps-support@mail.jmlab.jp" << std::endl;
            std::cerr << "       Licence: MIT (see, https://github.com/FDPS/FDPS/blob/master/LICENSE)" << std::endl;
            std::cerr << "       Note   : Please cite the following papers." << std::endl;
            std::cerr << "                - Iwasawa et al. (2016, Publications of the Astronomical Society of Japan, 68, 54)" << std::endl;
            std::cerr << "                - Namekata et al. (2018, Publications of the Astronomical Society of Japan, 70, 70)" << std::endl;
            std::cerr << "" << std::endl;
            std::cerr << "       Copyright (C) 2015 " << std::endl;
            std::cerr << "         Masaki Iwasawa, Ataru Tanikawa, Natsuki Hosono," << std::endl;
            std::cerr << "         Keigo Nitadori, Takayuki Muranushi, Daisuke Namekata," << std::endl;
            std::cerr << "         Kentaro Nomura, Junichiro Makino and many others" << std::endl;
#ifdef MONAR
            if(flag_monar){
                std::cerr<<"　　 ^__^　 ／‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"<<std::endl;
                std::cerr<<"　　( ´∀｀)＜******** FDPS has successfully begun. ********"<<std::endl;
                std::cerr<<"　　(     ) ＼＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿"<<std::endl;
                std::cerr<<"　　|  | |"<<std::endl;
                std::cerr<<"　　(__)_)"<<std::endl;
            }
            else if(flag_MONAR){
                std::cerr<<"        ∧_∧   ／‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"<<std::endl;
                std::cerr<<"       (´Д`) <  ******** FDPS has successfully begun. ********"<<std::endl;
                std::cerr<<"       ／_ /  ＼"<<std::endl;
                std::cerr<<"      (ぃ９｜  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"<<std::endl;
                std::cerr<<"      /　　/、"<<std::endl;
                std::cerr<<"     /　　∧_二つ"<<std::endl;
                std::cerr<<"     ｜　　＼ "<<std::endl;
                std::cerr<<"     /　/~＼ ＼"<<std::endl;
                std::cerr<<"    /　/　　>　)"<<std::endl;
                std::cerr<<"   ／ ノ　　/ ／"<<std::endl;
                std::cerr<<"  / ／　　 / /"<<std::endl;
                std::cerr<<"`/ /　　　( ヽ"<<std::endl;
                std::cerr<<"(＿)　　　 ＼_)"<<std::endl;
            }
            else{
                fprintf(stderr, "******** FDPS has successfully begun. ********\n");
            }
#else //MONAR
        fprintf(stderr, "******** FDPS has successfully begun. ********\n");
#endif //MONAR
        }
    }

    static inline void Finalize(){
        MemoryPool::finalize();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Finalize();
#endif  
        bool flag_monar = false;
        if(Comm::getRank() == 0) {
            if(flag_monar){
            }
            else{
                fprintf(stderr, "******** FDPS has successfully finished. ********\n");
            }
        }
    }
    
    static const S32 N_SMP_PTCL_TOT_PER_PSYS_DEFAULT = 1000000;

    ////////
    // util 
    inline F64 CalcSeparationSQPointToBox(const F64vec & point, 
                                          const F64vec & center, 
                                          const F64vec & half_length ){
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        auto dx = fabs(point.x - center.x);
        auto dy = fabs(point.y - center.y);
        dx = ( half_length.x < dx ) ? (dx - half_length.x) : 0.0;
        dy = ( half_length.y < dy ) ? (dy - half_length.y) : 0.0;
        return dx*dx + dy*dy;
#else
        auto dx = fabs(point.x - center.x);
        auto dy = fabs(point.y - center.y);
        auto dz = fabs(point.z - center.z);
        dx = ( half_length.x < dx ) ? (dx - half_length.x) : 0.0;
        dy = ( half_length.y < dy ) ? (dy - half_length.y) : 0.0;
        dz = ( half_length.z < dz ) ? (dz - half_length.z) : 0.0;
        return dx*dx + dy*dy + dz*dz;
#endif
    }

    // for check function
    inline bool IsInBox(const F64vec & pos,
                        const F64vec & center,
                        const F64 half_length,
                        const F64 tolerance=1e-6){
        const F64 tol = -std::fabs(tolerance);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        return ((pos.x - (center.x-half_length)) >= tol) && (((center.x+half_length) - pos.x) >= tol) &&
            ((pos.y - (center.y-half_length)) >= tol) && (((center.y+half_length) - pos.y) >= tol);
#else
        return ((pos.x - (center.x - half_length)) >= tol) && (((center.x+half_length) - pos.x) >= tol) &&
            ((pos.y - (center.y - half_length)) >= tol) && (((center.y+half_length) - pos.y) >= tol) &&
            ((pos.z - (center.z - half_length)) >= tol) &&  (((center.z+half_length) - pos.z) >= tol);
#endif
    }

    template<class T>
    class Abs{
    public:
        T operator () (const T val){
            return std::abs(val);
        }
    };
    template<class T>
    inline S32 Unique(T val[], const S32 n){
        S32 ret = 0;
        if( n > 0){
            ret = 1;
            T ref = val[0];
            for(S32 i=1; i<n; i++){
                if(val[i] > ref){
                    ref = val[i];
                    val[ret] = ref;
                    ret++;
                }
            }
        }
        return ret;
    }

    template<class T>
    inline T GetMSB(const T val);
    template<>
    inline U64 GetMSB(const U64 val){
        return (val>>63) & 0x1;
    }
    template<>
    inline U32 GetMSB(const U32 val){
        return (val>>31) & 0x1;
    }

    template<class T>
    inline T ClearMSB(const T val);
    template<>
    inline U64 ClearMSB(const U64 val){
        return val & 0x7fffffffffffffff;
    }
    template<>
    inline U32 ClearMSB(const U32 val){
        return val & 0x7fffffff;
    }

    template<class T>
    inline T SetMSB(const T val);
    template<>
    inline U64 SetMSB(const U64 val){
        return val | 0x8000000000000000;
    }
    template<>
    inline U32 SetMSB(const U32 val){
        return val | 0x80000000;
    }



    template<class Tp>
    inline F64ort GetMinBoxSingleThread(const Tp ptcl[], const S32 n){
        F64ort pos_box(ptcl[0].getPos(), ptcl[0].getPos());
        for(S32 i=1; i<n; i++){
            pos_box.merge(ptcl[i].getPos());
        }
        return pos_box;
    }
    
    template<class Tp>
    inline F64ort GetMinBox(const Tp ptcl[], const S32 n){
        F64ort box_loc;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif //PARTICLE_SIMULATOR_THREAD_PARALLEL
        {
            F64ort box_loc_tmp;
            box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
            for(S32 ip=0; ip<n; ip++){
                box_loc_tmp.merge(ptcl[ip].getPos());
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
            {
                box_loc.merge(box_loc_tmp);
            }
        }
        const F64vec xlow_loc = box_loc.low_;
        const F64vec xhigh_loc = box_loc.high_;
        const F64vec xlow_glb = Comm::getMinValue(xlow_loc);
        const F64vec xhigh_glb = Comm::getMaxValue(xhigh_loc);
        return F64ort(xlow_glb, xhigh_glb);
    }
    
    
    template<class Tp>
    inline F64ort GetMinBoxWithMargen(const Tp ptcl[], const S32 n){
        F64ort box_loc;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            F64ort box_loc_tmp;
            box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
            for(S32 ip=0; ip<n; ip++){
                box_loc_tmp.merge(ptcl[ip].getPos(), ptcl[ip].getRSearch());
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
            {
                box_loc.merge(box_loc_tmp);
            }
        }
        const F64vec xlow_loc = box_loc.low_;
        const F64vec xhigh_loc = box_loc.high_;
        const F64vec xlow_glb = Comm::getMinValue(xlow_loc);
        const F64vec xhigh_glb = Comm::getMaxValue(xhigh_loc);
        return F64ort(xlow_glb, xhigh_glb);
    }

    inline F64 GetWtime(){
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
#ifdef PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
        Comm::barrier();
#endif //PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
        return MPI_Wtime();
#elif defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
        //PARTICLE_SIMULATOR_THREAD_PARALLEL
        return omp_get_wtime();
#else
        return (F64)clock() / CLOCKS_PER_SEC;
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
    }


    inline F64 GetWtimeNoBarrier(){
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        return MPI_Wtime();
#elif defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
        return omp_get_wtime();
#else
        return clock() / CLOCKS_PER_SEC;
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
    }


    struct LessOPForVecX{
        bool operator() (const F64vec & left, const F64vec & right) const {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            return (left.x == right.x) ? ((left.y == right.y) ? true : right.y < right.y ) : left.x < right.x;
#else
            return (left.x == right.x) ? ((left.y == right.y) ? ((left.z == right.z) ? true : left.z < right.z) : left.y < right.y ) : left.x < right.x;
#endif
        }
    };

    struct LessOPX{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.x < right.x;
        }
    };
    struct LessOPY{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.y < right.y;
        }
    };
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    struct LessOPZ{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.z < right.z;
        }
    };
#endif
    struct LessOPKEY{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.key_ < right.key_;
        }
    };


    struct TagRSearch{};
    struct TagNoRSearch{};

    template<bool T>
    struct HasRSearchInner{
        typedef TagRSearch type;
    };
    template<>
    struct HasRSearchInner<false>{
        typedef TagNoRSearch type;
    };

    template <class T>
    class HasRSearch{
    private:
        typedef char One[1];
        typedef char Two[2];
    
        template <class T2, T2>
        class Check{};

        template <typename T3>
        static One & Func( Check<double (T3::*)(void), &T3::getRSearch>* );
        template <typename T3>
        static One & Func( Check<float (T3::*)(void), &T3::getRSearch>* );
        template <typename T3>
        static One & Func( Check<double (T3::*)(void) const, &T3::getRSearch>* );
        template <typename T3>
        static One & Func( Check<float (T3::*)(void) const, &T3::getRSearch>* );

        template <typename T3>
        static Two &  Func(...);
    
    public:
        typedef typename HasRSearchInner< sizeof(Func<T>(NULL)) == 1 >::type type;
    };

    template<typename Tptcl>
    struct HasgetRSearchMethod
    {
#ifndef PARTICLE_SIMULATOR_ALL_64BIT_PRECISION
       template<typename U, F32(U::*)() > struct SFINAE0 {};
       template<typename U, F32(U::*)() const > struct SFINAE1 {};
#endif
       template<typename U, F64(U::*)() > struct SFINAE2 {};
       template<typename U, F64(U::*)() const > struct SFINAE3 {};
#ifndef PARTICLE_SIMULATOR_ALL_64BIT_PRECISION
       template<typename U> static char Test(SFINAE0<U, &U::getRSearch> *);
       template<typename U> static char Test(SFINAE1<U, &U::getRSearch> *);
#endif
       template<typename U> static char Test(SFINAE2<U, &U::getRSearch> *);
       template<typename U> static char Test(SFINAE3<U, &U::getRSearch> *);
       template<typename U> static int Test(...);
       static const bool value = sizeof(Test<Tptcl>(0)) == sizeof(char);
    };
    template<class Tptcl>
    F64 GetMyRSearch(Tptcl ptcl, std::true_type)
    {
       return ptcl.getRSearch();
    }
    template<class Tptcl>
    F64 GetMyRSearch(Tptcl ptcl, std::false_type)
    {
       return 0.0;
    }
    template <class Tptcl>
    F64 GetMyRSearch(Tptcl ptcl) {
       return GetMyRSearch(ptcl, std::integral_constant<bool, HasgetRSearchMethod<Tptcl>::value>());
    }
    
    class TimeProfile{
    public:
        F64 collect_sample_particle;
        F64 decompose_domain;
        F64 exchange_particle;
        F64 set_particle_local_tree;
        F64 set_particle_global_tree;
        F64 make_local_tree;
        F64 make_global_tree;
        F64 set_root_cell;
        F64 calc_force;
        F64 calc_moment_local_tree;
        F64 calc_moment_global_tree;
        F64 make_LET_1st;
        F64 make_LET_2nd;
        F64 exchange_LET_1st;
        F64 exchange_LET_2nd;
        F64 write_back; // new but default

        // not public
        F64 morton_sort_local_tree;
        F64 link_cell_local_tree;
        F64 morton_sort_global_tree;
        F64 link_cell_global_tree;

        F64 make_local_tree_tot; // make_local_tree + calc_moment_local_tree
        F64 make_global_tree_tot;
        F64 exchange_LET_tot; // make_LET_1st + make_LET_2nd + exchange_LET_1st + exchange_LET_2nd

        F64 calc_force__core__walk_tree;
        F64 calc_force__core__keep_list;
        F64 calc_force__core__copy_ep;
        F64 calc_force__core__dispatch;
        F64 calc_force__core__retrieve;

        F64 calc_force__make_ipgroup;
        F64 calc_force__core;
        F64 calc_force__copy_original_order;

        F64 exchange_particle__find_particle;
        F64 exchange_particle__find_particle_2;
        F64 exchange_particle__find_particle_3;
        F64 exchange_particle__exchange_particle;

        F64 decompose_domain__sort_particle_1st;
        F64 decompose_domain__sort_particle_2nd;
        F64 decompose_domain__sort_particle_3rd;
        F64 decompose_domain__gather_particle;

        F64 decompose_domain__setup;
        F64 decompose_domain__determine_coord_1st;
        F64 decompose_domain__migrae_particle_1st;
        F64 decompose_domain__determine_coord_2nd;
        F64 decompose_domain__determine_coord_3rd;
        F64 decompose_domain__exchange_pos_domain;

        F64 exchange_LET_1st__a2a_n;
        F64 exchange_LET_1st__icomm_ptcl;
        F64 exchange_LET_1st__a2a_sp;
        F64 exchange_LET_1st__a2a_ep;

        F64 add_moment_as_sp_local;
        F64 add_moment_as_sp_global;

        F64 make_LET_1st__find_particle;
        F64 make_LET_1st__exchange_n;
        F64 make_LET_1st__exchange_top_moment;

        F64 getTotalTime() const {
            return collect_sample_particle + decompose_domain + exchange_particle
                + set_particle_local_tree + set_particle_global_tree
                + set_root_cell
                + calc_force + calc_moment_local_tree + calc_moment_global_tree + make_LET_1st + make_LET_2nd + exchange_LET_1st + exchange_LET_2nd
                + morton_sort_local_tree + link_cell_local_tree 
                + morton_sort_global_tree + link_cell_global_tree
                + add_moment_as_sp_local + add_moment_as_sp_global
                + write_back;
        }
        
        void dump(std::ostream & fout=std::cout,
                  const S32 level=0) const {
            fout<<"total_time= "<<getTotalTime()<<std::endl;
            fout<<"  collect_sample_particle= "<<collect_sample_particle<<std::endl;
            fout<<"  decompose_domain= "<<decompose_domain<<std::endl;
            fout<<"  exchange_particle= "<<exchange_particle<<std::endl;
            fout<<"  set_particle_local_tree= "<<set_particle_local_tree<<std::endl;
            fout<<"  set_particle_global_tree= "<<set_particle_global_tree<<std::endl;
            fout<<"  set_root_cell= "<<set_root_cell<<std::endl;
            fout<<"  morton_sort_local_tree= "<<morton_sort_local_tree<<std::endl;
            fout<<"  link_cell_local_tree= "<<link_cell_local_tree<<std::endl;
            fout<<"  morton_sort_global_tree= "<<morton_sort_global_tree<<std::endl;
            fout<<"  link_cell_global_tree= "<<link_cell_global_tree<<std::endl;
            fout<<"  calc_force= "<<calc_force<<std::endl;
            fout<<"  calc_moment_local_tree= "<<calc_moment_local_tree<<std::endl;
            fout<<"  calc_moment_global_tree= "<<calc_moment_global_tree<<std::endl;
            fout<<"  add_moment_as_sp_global= "<<add_moment_as_sp_global<<std::endl;
            fout<<"  make_LET_1st= "<<make_LET_1st<<std::endl;
            if(level > 0){
                fout<<"    make_LET_1st__find_particle= "<<make_LET_1st__find_particle<<std::endl;
                fout<<"    make_LET_1st__exchange_n= "<<make_LET_1st__exchange_n<<std::endl;
                fout<<"    make_LET_1st__exchange_top_moment= "<<make_LET_1st__exchange_top_moment<<std::endl;
            }
            fout<<"  make_LET_2nd= "<<make_LET_2nd<<std::endl;
            fout<<"  exchange_LET_1st= "<<exchange_LET_1st<<std::endl;
            if(level > 0){
                fout<<"    exchange_LET_1st__icomm_ptcl= "<<exchange_LET_1st__icomm_ptcl<<std::endl;
            }
            fout<<"  exchange_LET_2nd= "<<exchange_LET_2nd<<std::endl;
            fout<<"  write_back= "<<write_back<<std::endl;
        }

        TimeProfile () {
            collect_sample_particle = decompose_domain = exchange_particle = set_particle_local_tree = set_particle_global_tree = make_local_tree = make_global_tree = set_root_cell
                = calc_force = calc_moment_local_tree = calc_moment_global_tree = make_LET_1st = make_LET_2nd 
                = exchange_LET_1st = exchange_LET_2nd = 0.0;
            morton_sort_local_tree = link_cell_local_tree
                = morton_sort_global_tree = link_cell_global_tree = 0.0;
            make_local_tree_tot = make_global_tree_tot = exchange_LET_tot = 0.0;
            calc_force__make_ipgroup = calc_force__core = calc_force__copy_original_order = 0.0;

            exchange_particle__find_particle = exchange_particle__exchange_particle = 0.0;

	    exchange_particle__find_particle_2 = exchange_particle__find_particle_3 = 0.0;

            decompose_domain__sort_particle_1st = decompose_domain__sort_particle_2nd = decompose_domain__sort_particle_3rd = decompose_domain__gather_particle = 0.0;
            decompose_domain__setup = decompose_domain__determine_coord_1st = decompose_domain__migrae_particle_1st = decompose_domain__determine_coord_2nd 
                = decompose_domain__determine_coord_3rd = decompose_domain__exchange_pos_domain = 0.0;
            exchange_LET_1st__a2a_n = exchange_LET_1st__a2a_sp = exchange_LET_1st__icomm_ptcl = exchange_LET_1st__a2a_ep = 0.0;

            add_moment_as_sp_local = add_moment_as_sp_global = 0.0;
            write_back = 0.0;
            make_LET_1st__find_particle = make_LET_1st__exchange_n = make_LET_1st__exchange_top_moment = 0.0;
        }

        TimeProfile operator + (const TimeProfile & rhs) const{
            TimeProfile ret;
            ret.collect_sample_particle = this->collect_sample_particle + rhs.collect_sample_particle;
            ret.decompose_domain = this->decompose_domain + rhs.decompose_domain;
            ret.exchange_particle = this->exchange_particle + rhs.exchange_particle;
            ret.set_particle_local_tree = this->set_particle_local_tree + rhs.set_particle_local_tree;
            ret.set_particle_global_tree = this->set_particle_global_tree + rhs.set_particle_global_tree;
            ret.set_root_cell = this->set_root_cell + rhs.set_root_cell;
            ret.make_local_tree = this->make_local_tree + rhs.make_local_tree;
            ret.make_global_tree = this->make_global_tree + rhs.make_global_tree;
            ret.calc_force = this->calc_force + rhs.calc_force;
            ret.calc_moment_local_tree = this->calc_moment_local_tree + rhs.calc_moment_local_tree;
            ret.calc_moment_global_tree = this->calc_moment_global_tree + rhs.calc_moment_global_tree;
            ret.make_LET_1st = this->make_LET_1st + rhs.make_LET_1st;
            ret.make_LET_2nd = this->make_LET_2nd + rhs.make_LET_2nd;
            ret.exchange_LET_1st = this->exchange_LET_1st + rhs.exchange_LET_1st;
            ret.exchange_LET_2nd = this->exchange_LET_2nd + rhs.exchange_LET_2nd;

            ret.morton_sort_local_tree = this->morton_sort_local_tree + rhs.morton_sort_local_tree;
            ret.link_cell_local_tree = this->link_cell_local_tree + rhs.link_cell_local_tree;
            ret.morton_sort_global_tree = this->morton_sort_global_tree + rhs.morton_sort_global_tree;
            ret.link_cell_global_tree = this->link_cell_global_tree + rhs.link_cell_global_tree;

            ret.make_local_tree_tot = this->make_local_tree_tot + rhs.make_local_tree_tot;
            ret.make_global_tree_tot = this->make_global_tree_tot + rhs.make_global_tree_tot;
            ret.exchange_LET_tot = this->exchange_LET_tot + rhs.exchange_LET_tot;

            ret.calc_force__core__walk_tree = this->calc_force__core__walk_tree + rhs.calc_force__core__walk_tree;
            ret.calc_force__core__keep_list = this->calc_force__core__keep_list + rhs.calc_force__core__keep_list;
            ret.calc_force__core__dispatch = this->calc_force__core__dispatch + rhs.calc_force__core__dispatch;
            ret.calc_force__core__copy_ep = this->calc_force__core__copy_ep + rhs.calc_force__core__copy_ep;
            ret.calc_force__core__retrieve = this->calc_force__core__retrieve + rhs.calc_force__core__retrieve;

            ret.calc_force__make_ipgroup = this->calc_force__make_ipgroup + rhs.calc_force__make_ipgroup;
            ret.calc_force__core = this->calc_force__core + rhs.calc_force__core;
            ret.calc_force__copy_original_order = this->calc_force__copy_original_order + rhs.calc_force__copy_original_order;

            ret.exchange_particle__find_particle     = this->exchange_particle__find_particle     + rhs.exchange_particle__find_particle;  
            ret.exchange_particle__exchange_particle = this->exchange_particle__exchange_particle + rhs.exchange_particle__exchange_particle;


            ret.decompose_domain__sort_particle_1st = this->decompose_domain__sort_particle_1st + rhs.decompose_domain__sort_particle_1st;
            ret.decompose_domain__sort_particle_2nd = this->decompose_domain__sort_particle_2nd + rhs.decompose_domain__sort_particle_2nd;
            ret.decompose_domain__sort_particle_3rd = this->decompose_domain__sort_particle_3rd + rhs.decompose_domain__sort_particle_3rd;
            ret.decompose_domain__gather_particle = this->decompose_domain__gather_particle + rhs.decompose_domain__gather_particle;

            ret.decompose_domain__setup = this->decompose_domain__setup + rhs.decompose_domain__setup;
            ret.decompose_domain__determine_coord_1st = this->decompose_domain__determine_coord_1st + rhs.decompose_domain__determine_coord_1st;
            ret.decompose_domain__migrae_particle_1st = this->decompose_domain__migrae_particle_1st + rhs.decompose_domain__migrae_particle_1st;
            ret.decompose_domain__determine_coord_2nd = this->decompose_domain__determine_coord_2nd + rhs.decompose_domain__determine_coord_2nd;
            ret.decompose_domain__determine_coord_3rd = this->decompose_domain__determine_coord_3rd + rhs.decompose_domain__determine_coord_3rd;
            ret.decompose_domain__exchange_pos_domain = this->decompose_domain__exchange_pos_domain + rhs.decompose_domain__exchange_pos_domain;

            ret.exchange_LET_1st__a2a_n    = this->exchange_LET_1st__a2a_n    + rhs.exchange_LET_1st__a2a_n;
            ret.exchange_LET_1st__icomm_ptcl   = this->exchange_LET_1st__icomm_ptcl   + rhs.exchange_LET_1st__icomm_ptcl;
            ret.exchange_LET_1st__a2a_sp   = this->exchange_LET_1st__a2a_sp   + rhs.exchange_LET_1st__a2a_sp;
            ret.exchange_LET_1st__a2a_ep   = this->exchange_LET_1st__a2a_ep   + rhs.exchange_LET_1st__a2a_ep;

            ret.add_moment_as_sp_local = this->add_moment_as_sp_local + rhs.add_moment_as_sp_local;
            ret.add_moment_as_sp_global = this->add_moment_as_sp_global + rhs.add_moment_as_sp_global;

            ret.write_back = this->write_back + rhs.write_back;

            ret.make_LET_1st__find_particle = this->make_LET_1st__find_particle + rhs.make_LET_1st__find_particle;
            ret.make_LET_1st__exchange_n = this->make_LET_1st__exchange_n + rhs.make_LET_1st__exchange_n;
            ret.make_LET_1st__exchange_top_moment = this->make_LET_1st__exchange_top_moment + rhs.make_LET_1st__exchange_top_moment;
            
            return ret;
        }

        const TimeProfile & operator += (const TimeProfile & rhs) {
            (*this) = (*this) + rhs;
            return (*this);
        }

        void clear(){
            collect_sample_particle = decompose_domain = exchange_particle = make_local_tree = make_global_tree = set_particle_local_tree = set_particle_global_tree = set_root_cell
                = calc_force = calc_moment_local_tree = calc_moment_global_tree = make_LET_1st = make_LET_2nd = exchange_LET_1st = exchange_LET_2nd = 0.0;
            morton_sort_local_tree = link_cell_local_tree 
                = morton_sort_global_tree = link_cell_global_tree = 0.0;
            make_local_tree_tot = make_global_tree_tot = exchange_LET_tot = 0.0;
            calc_force__core__walk_tree = 0.0;
            calc_force__core__keep_list = 0.0;
            calc_force__core__copy_ep   = 0.0;
            calc_force__core__dispatch  = 0.0;
            calc_force__core__retrieve  = 0.0;
            calc_force__make_ipgroup = calc_force__core = calc_force__copy_original_order = 0.0;
            exchange_particle__find_particle = exchange_particle__exchange_particle = 0.0;

	    exchange_particle__find_particle_2 = exchange_particle__find_particle_3 = 0.0;

            decompose_domain__sort_particle_1st = decompose_domain__sort_particle_2nd = decompose_domain__sort_particle_3rd = decompose_domain__gather_particle = 0.0;
            decompose_domain__setup = decompose_domain__determine_coord_1st = decompose_domain__migrae_particle_1st = decompose_domain__determine_coord_2nd 
                = decompose_domain__determine_coord_3rd = decompose_domain__exchange_pos_domain = 0.0;
            exchange_LET_1st__a2a_n = exchange_LET_1st__a2a_sp = exchange_LET_1st__icomm_ptcl = exchange_LET_1st__a2a_ep = 0.0;

            add_moment_as_sp_local = add_moment_as_sp_global = 0.0;
            write_back = 0.0;
            make_LET_1st__find_particle = make_LET_1st__exchange_n = make_LET_1st__exchange_top_moment = 0.0;
        }
    };

    inline std::string GetBinString(const U32 val)
    {
        if( !val ) return std::string("00000000000000000000000000000000");
        U32 tmp = val;
        std::string str;
        for (S32 i=0; i<32; i++) {
            if ( (tmp & 1) == 0 ) str.insert(str.begin(), '0');
            else str.insert(str.begin(), '1');
            tmp >>= 1;
        }
        return str;
    }
    inline std::string GetBinString(const U32 *val)
    {
        if( !*val ) return std::string("00000000000000000000000000000000");
        U32 tmp = *val;
        std::string str;
        for (S32 i=0; i<32; i++) {
            if ( (tmp & 1) == 0 ) str.insert(str.begin(), '0');
            else str.insert(str.begin(), '1');
            tmp >>= 1;
        }
        return str;
    }
    inline std::string GetBinString(const U64 val)
    {
        if( !val ) return std::string("0000000000000000000000000000000000000000000000000000000000000000");
        U64 tmp = val;
        std::string str;
        for (S32 i=0; i<64; i++) {
            if ( (tmp & 1) == 0 ) str.insert(str.begin(), '0');
            else str.insert(str.begin(), '1');
            tmp >>= 1;
        }
        return str;
    }
    inline std::string GetBinString(const U64 *val)
    {
        if( !*val ) return std::string("0000000000000000000000000000000000000000000000000000000000000000");
        U64 tmp = *val;
        std::string str;
        for (S32 i=0; i<64; i++) {
            if ( (tmp & 1) == 0 ) str.insert(str.begin(), '0');
            else str.insert(str.begin(), '1');
            tmp >>= 1;
        }
        return str;
    }

#ifdef TEST_VARIADIC_TEMPLATE
    template<class T>
    class IsParticleSystem{
        template<class T2>
        static auto check(T2 x) -> decltype(x.ParticleSystemDummyFunc(), std::true_type());
        static auto check(...)  -> decltype(std::false_type());
    public:
        typedef decltype(check(std::declval<T>())) value;
    };
    template<class T>
    class IsDomainInfo{
        template<class T2>
        static auto check(T2 x) -> decltype(x.DomainInfoDummyFunc(), std::true_type());
        static auto check(...)  -> decltype(std::false_type());
    public:
        typedef decltype(check(std::declval<T>())) value;
    };
#endif

    template<typename T>
    inline void PackData(const ReallocatableArray<T> * data_in,
			 const S32 n_buf,
			 ReallocatableArray<T> & data_out,
			 const S32 offset = 0){
        S32 size_data = offset;
PS_OMP(omp parallel for reduction(+:size_data))
        for(S32 i=0; i<n_buf; i++){
            size_data += data_in[i].size();
        }
        data_out.resizeNoInitialize(size_data);
PS_OMP(omp parallel for)
        for(S32 i=0; i<n_buf; i++){
            S32 head_adr = offset;
            for(S32 j=0; j<i; j++){
                head_adr += data_in[j].size();
            }
            for(S32 j=0; j<data_in[i].size(); j++){
                data_out[j+head_adr] = data_in[i][j];
            }
        }        
    }

    template<typename T>
    inline void PackDataInOmp(const ReallocatableArray<T> * data_in,
			      ReallocatableArray<T> & data_out,
			      const S32 nth,
			      const S32 ith,
			      const S32 offset = 0){
	if(ith == 0){
	    S32 size_data = offset;
	    for(S32 i=0; i<nth; i++){
		size_data += data_in[i].size();
	    }
	    data_out.resizeNoInitialize(size_data);
	}
PS_OMP_BARRIER
	S32 head_adr = offset;
	for(S32 j=0; j<ith; j++){
	    head_adr += data_in[j].size();
	}
	for(S32 j=0; j<data_in[ith].size(); j++){
	    data_out[j+head_adr] = data_in[ith][j];
	}
    }



    inline void CalcAdrToSplitData(S32 & head,
                                   S32 & end,
                                   const S32 i_th,
                                   const S32 n_div,
                                   const S32 n_tot){
        head = (n_tot/n_div)*i_th + std::min(n_tot%n_div, i_th);
        end  = (n_tot/n_div)*(i_th+1) + std::min(n_tot%n_div, (i_th+1));
    }

    template<typename Tsrc, typename Tdst, typename Tcopy>
    inline void RandomSampling(Tsrc * val_src,
                               Tdst * val_dst,
                               const S32 n_src,
                               const S32 n_dst,
                               Tcopy copy){
        thread_local std::mt19937 mt(Comm::getRank()*Comm::getNumberOfThread()+Comm::getThreadNum());
        if(n_src == 0) return;
	std::uniform_int_distribution<S32> dist(0, n_src-1);
        ReallocatableArray<S32> record(n_dst, n_dst, MemoryAllocMode::Pool);
        for(S32 i = 0; i < n_dst; i++) {
            S32 j = dist(mt);
            Tsrc hold = val_src[j];
            val_src[j]   = val_src[i];
            val_src[i]   = hold;
            record[i]  = j;
        }
        for(S32 i=0; i<n_dst; i++) {
            copy(val_src[i], val_dst[i]);
        }
        for(S32 i = n_dst - 1; i >= 0; i--) {
            S32 j = record[i];
            Tsrc hold = val_src[j];
            val_src[j]   = val_src[i];
            val_src[i]   = hold;
        }
    }

    template<typename Tsrc, typename Tdst, typename Tcopy>
    inline void RandomSampling(Tsrc * val_src,
                               ReallocatableArray<Tdst> & val_dst,
                               const S32 n_src,
                               Tcopy copy){
        const auto n_dst = val_dst.size();
        thread_local std::mt19937 mt(Comm::getRank()*Comm::getNumberOfThread()+Comm::getThreadNum());
        if(n_src == 0) return;
	std::uniform_int_distribution<S32> dist(0, n_src-1);
        ReallocatableArray<S32> record(n_dst, n_dst, MemoryAllocMode::Pool);
        for(S32 i = 0; i < n_dst; i++) {
            S32 j = dist(mt);
            Tsrc hold = val_src[j];
            val_src[j]   = val_src[i];
            val_src[i]   = hold;
            record[i]  = j;
        }
        for(S32 i=0; i<n_dst; i++) {
            copy(val_src[i], val_dst[i]);
        }
        for(S32 i = n_dst - 1; i >= 0; i--) {
            S32 j = record[i];
            Tsrc hold = val_src[j];
            val_src[j]   = val_src[i];
            val_src[i]   = hold;
        }
    }

    
    template<typename Tsrc, typename Tdst, typename Tcopy>
    inline void RandomSampling(ReallocatableArray<Tsrc> & val_src_ar,
                               ReallocatableArray<Tdst> & val_dst_ar,
                               Tcopy copy,
                               const S32 offset_src=0,
                               const S32 offset_dst=0){
        const auto n_src = val_src_ar.size() - offset_src;
        const auto n_dst = val_dst_ar.size() - offset_dst;
        thread_local std::mt19937 mt(Comm::getRank()*Comm::getNumberOfThread()+Comm::getThreadNum());
        //thread_local std::mt19937 mt(0);
        if(n_src == 0) return;
	std::uniform_int_distribution<S32> dist(0, n_src-1);
        ReallocatableArray<S32> record(n_dst, n_dst, MemoryAllocMode::Pool);
        Tsrc * val_src = val_src_ar.getPointer(offset_src);
        Tdst * val_dst = val_dst_ar.getPointer(offset_dst);
        for(S32 i = 0; i < n_dst; i++) {
            S32 j = dist(mt);
            Tsrc hold = val_src[j];
            val_src[j]   = val_src[i];
            val_src[i]   = hold;
            record[i]  = j;
        }
        for(S32 i=0; i<n_dst; i++) {
            copy(val_src[i], val_dst[i]);
        }
        for(S32 i = n_dst - 1; i >= 0; i--) {
            S32 j = record[i];
            Tsrc hold = val_src[j];
            val_src[j]   = val_src[i];
            val_src[i]   = hold;
        }
    }
    
    inline void CountNObjInBucket(const S32 n_bucket,
                                  S32 * n_obj_in_bucket,
                                  const S32 n_obj,
                                  const S32 * id_bucket){
        for(S32 i=0; i<n_bucket; i++){
            n_obj_in_bucket[i] = 0;
        }
        for(S32 i=0; i<n_obj; i++){
            n_obj_in_bucket[id_bucket[id_bucket[i]]]++;
        }
    }
    
    inline void CountNObjInBucketOMP(const S32 n_bucket,
                                     S32 * n_obj_in_bucket,
                                     const S32 n_obj,
                                     const S32 * id_bucket){
        const auto n_thread = Comm::getNumberOfThread();
        std::vector< ReallocatableArray<S32> > n_obj_in_bucket_tmp(n_thread);
        for(S32 i=0; i<n_thread; i++){
            n_obj_in_bucket_tmp[i].setAllocMode(MemoryAllocMode::Pool);
            n_obj_in_bucket_tmp[i].resizeNoInitialize(n_bucket);
        }
PS_OMP_PARALLEL
        {
            const auto ith = Comm::getThreadNum();
            S32 head, end;
            CalcAdrToSplitData(head, end, ith, n_thread, n_obj);
            CountNObjInBucket(n_bucket, &n_obj_in_bucket_tmp[ith][0],
                              end-head, &id_bucket[head]);
        }
        for(S32 i=0; i<n_bucket; i++){
            n_obj_in_bucket[i] = 0;
        }
        for(S32 i=0; i<n_bucket; i++){
            for(S32 j=0; j<n_thread; j++){
                n_obj_in_bucket[i] += n_obj_in_bucket_tmp[j][i];
            }
        }
    }
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
    class KeyT{
    public:
        U64 hi_;
        KeyT() : hi_(0){}
        KeyT(const U64 hi) : hi_(hi){}
        void init(){
            hi_ = 0;
        }
        bool operator == (const KeyT & rhs) const {
            return (hi_ == rhs.hi_);
        }
        bool operator != (const KeyT & rhs) const {
            return !(*this == rhs);
        }
        bool operator < (const KeyT & rhs) const {
            return (hi_ < rhs.hi_);
        }
        bool operator <= (const KeyT & rhs) const {
            return (hi_ <= rhs.hi_);
        }
        bool operator > (const KeyT & rhs) const {
            return !(*this <= rhs);
        }
        bool operator >= (const KeyT & rhs) const {
            return !(*this < rhs);
        }
        
        KeyT operator & (const KeyT & rhs) const {
            return KeyT(hi_ & rhs.hi_);
        }
        KeyT operator ^ (const KeyT & rhs) const {
            return KeyT(hi_ ^ rhs.hi_);
        }
        KeyT operator | (const KeyT & rhs) const {
            return KeyT(hi_ | rhs.hi_);
        }
        const KeyT & operator |= (const KeyT & rhs){
            (this->hi_) |= rhs.hi_;
            return (*this);
        }
        // cast
        operator int() const {return (int)hi_;}
        operator unsigned int() const {return (unsigned int)hi_;}
        operator long() const {return (long)hi_;}
        operator unsigned long() const {return (unsigned long)hi_;}

        template<typename T>
        KeyT operator << (const T & s) const {
            return KeyT(hi_ << s);
        }
        template<typename T>
        KeyT operator >> (const T & s) const {
            return KeyT(hi_ >> s);
        }

        void dump(std::ostream & fout=std::cerr) const {
            fout<<std::setw(TREE_LEVEL_LIMIT)<<std::setfill('0')<<std::oct<<hi_<<std::endl;
            fout<<std::dec;
        }
        
        friend std::ostream & operator << (std::ostream & c, const KeyT & k){
            c<<k.hi_;
            return c;
        }
    };
#else // for 96bit or 128bit
    class KeyT{
    public:
        U64 hi_;
#ifdef PARTICLE_SIMULATOR_USE_96BIT_KEY
        U32 lo_;
        KeyT(const U64 hi, const U32 lo) : hi_(hi), lo_(lo){}
#else
        U64 lo_;
        KeyT(const U64 hi, const U64 lo) : hi_(hi), lo_(lo){}
#endif
        KeyT() : hi_(0), lo_(0){}
        template<typename Tint>
        KeyT(const Tint lo) : hi_(0), lo_(lo){}
        void init(){
            hi_ = 0;
            lo_ = 0;
        }
        bool operator == (const KeyT & rhs) const {
            return (hi_ == rhs.hi_) && (lo_ == rhs.lo_);
        }
        bool operator != (const KeyT & rhs) const {
            return !(*this == rhs);
        }
        bool operator < (const KeyT & rhs) const {
            return (hi_ < rhs.hi_) || ((hi_ == rhs.hi_) && (lo_ < rhs.lo_));
        }
        bool operator <= (const KeyT & rhs) const {
            return (hi_ < rhs.hi_) || ((hi_ == rhs.hi_) && (lo_ <= rhs.lo_));
        }
        bool operator > (const KeyT & rhs) const {
            return !(*this <= rhs);
        }
        bool operator >= (const KeyT & rhs) const {
            return !(*this < rhs);
        }
        KeyT operator & (const KeyT & rhs) const {
            return KeyT(hi_&rhs.hi_, lo_&rhs.lo_);
        }
        KeyT operator ^ (const KeyT & rhs) const {
            return KeyT(hi_^rhs.hi_, lo_^rhs.lo_);
        }
        KeyT operator | (const KeyT & rhs) const {
            return KeyT(hi_|rhs.hi_, lo_|rhs.lo_);
        }
        const KeyT & operator |= (const KeyT & rhs){
            (this->hi_) |= rhs.hi_;
            (this->lo_) |= rhs.lo_;
            return (*this);
        }
        // cast
        operator int() const {return (int)lo_;}
        operator unsigned int() const {return (unsigned int)lo_;}
        operator long() const {return (long)lo_;}
        operator unsigned long() const {return (unsigned long)lo_;}

        // shift
        template<typename T>
        KeyT operator << (const T & s) const {
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
            U32 lo = ((U64)s < sizeof(lo_)*8) ? lo_<<s : 0;
            U64 lo64 = (U64)lo_;
            const U64 rem = ((U64)s < sizeof(hi_)*8) ? hi_<<s : 0;
            const U64 inc = ((U64)s < sizeof(lo_)*8) ? (lo64>>(sizeof(lo_)*8-s)) : lo64<<(s-sizeof(lo_)*8);
            U64 hi = rem | inc;
#else
            U64 lo = ((U64)s < sizeof(lo_)*8) ? lo_<<s : 0;
            const U64 rem = (s==0) ? 0 : (lo_>>(sizeof(hi_)*8-s));
            U64 hi = ((U64)s < sizeof(lo_)*8) ? ((hi_<<s) | rem) : lo_ << (s-sizeof(lo_)*8);
#endif
            return KeyT(hi, lo);
        }
        template<typename T>
        KeyT operator >> (const T & s) const {
#if defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
            U64 hi = ((U64)s < sizeof(hi_)*8) ? hi_>>s : 0;
            const U32 rem = ((U64)s < sizeof(lo_)*8) ? lo_ >> s : 0;
            const U32 inc = (s==0) ? 0 : ((U64)s < sizeof(hi_)*8) ? (U32)((hi_<<(sizeof(hi_)*8-s))>>(sizeof(lo_)*8)) : (U32)(hi_>>(s-sizeof(lo_)*8));
            U32 lo = rem | inc;
#else
            //assert(s < 128);
            U64 hi = ((U64)s < sizeof(hi_)*8) ? hi_>>s : 0;
            const U64 rem = (s==0) ? 0 : (hi_<<(sizeof(hi_)*8-s));
            U64 lo = ((U64)s < sizeof(hi_)*8) ? (lo_>>s | rem) : hi_>>(s-sizeof(hi_)*8);	    
#endif
            return KeyT(hi, lo);
        }

        void dump(std::ostream & fout=std::cerr) const {
            fout<<std::setw(KEY_LEVEL_MAX_HI)<<std::setfill('0')<<std::oct<<hi_
                <<std::setw(KEY_LEVEL_MAX_LO)<<std::setfill('0')<<lo_<<std::endl;
            fout<<std::dec;
        }
        
        friend std::ostream & operator <<(std::ostream & c, const KeyT & k){
            c<<k.hi_<<"   "<<k.lo_;
            return c;
        }
    };
    
#endif
    
#ifdef LOOP_TREE
    static const S32 ADR_TREE_CELL_NULL = 1;
    static const S32 SIMD_VEC_LEN = 8;
    static const S32 N_THREAD_LIMIT = 12;
    static const S32 BUF_SIZE_INTERACTION_LIST = 100000;
#endif


    template<typename T>
    struct my_is_char: std::false_type{};
    template<>  struct my_is_char<char*>: std::true_type{};
    template<>  struct my_is_char<const char*>: std::true_type{};
    template<>  struct my_is_char<const char* const>: std::true_type{};

    template<typename T>
    struct my_is_string: std::false_type{};
    template<>  struct my_is_string<std::string>: std::true_type{};
    template<>  struct my_is_string<const std::string>: std::true_type{};


	
    
}

#include"util.hpp"
#include"timer.hpp"
