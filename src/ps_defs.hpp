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
#include "time.h"

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include"mpi.h"
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include<omp.h>
#endif

#include"memory_pool.hpp"
#include"vector2.hpp"
#include"vector3.hpp"
#include"orthotope2.hpp"
#include"orthotope3.hpp"
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

namespace ParticleSimulator{
    static const long long int LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1ll<<30;
    static inline void Abort(const int err = -1){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
      MPI_Abort(MPI_COMM_WORLD,err);
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
    //typedef long long int S32;
    //typedef unsigned long long int U32;
#else
    typedef float F32;
    //typedef int S32;
    //typedef unsigned int U32;
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
    typedef F32ort2 F32ort;
    typedef F64ort2 F64ort;
    static const S32 DIMENSION = 2;
    static const S32 N_CHILDREN = 4;
    static const S32 TREE_LEVEL_LIMIT = 30;
    static const F32vec SHIFT_CENTER[N_CHILDREN] = 
        { F32vec(-0.5, -0.5), F32vec(-0.5, 0.5),
          F32vec( 0.5, -0.5), F32vec( 0.5, 0.5) };
#else
    typedef S32vec3 S32vec;
    typedef F32vec3 F32vec;
    typedef F64vec3 F64vec;
    typedef F32mat3 F32mat;
    typedef F64mat3 F64mat;
    typedef F32ort3 F32ort;
    typedef F64ort3 F64ort;
    static const S32 DIMENSION = 3;
    static const S32 N_CHILDREN = 8;
    static const S32 TREE_LEVEL_LIMIT = 21;
    static const F32vec SHIFT_CENTER[N_CHILDREN] = 
        { F32vec(-0.5, -0.5, -0.5), F32vec(-0.5, -0.5, 0.5),
          F32vec(-0.5,  0.5, -0.5), F32vec(-0.5,  0.5,  0.5),
          F32vec( 0.5, -0.5, -0.5), F32vec( 0.5, -0.5,  0.5),
          F32vec( 0.5,  0.5, -0.5), F32vec( 0.5,  0.5,  0.5) };
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

    struct TagSearchBoundaryConditionOpenOnly{};
    struct TagSearchBoundaryConditionOpenPeriodic{};
    
    struct TagSearchLong{};
    struct TagSearchLongCutoff{};
    struct TagSearchLongScatter{};
    struct TagSearchLongSymmetry{};
    struct TagSearchLongCutoffScatter{};
    struct TagSearchShortGather{};
    struct TagSearchShortScatter{};
    struct TagSearchShortSymmetry{};

    struct TagIpgLongNormal{};
    struct TagIpgIn{};
    struct TagIpgInAndOut{};
    struct TagIpgOut{};

    struct TagWithoutCutoff{};
    struct TagWithCutoff{};
    
    struct SEARCH_MODE_LONG{
        typedef TagForceLong force_type;
        typedef TagSearchLong search_type;
        typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
        typedef TagIpgLongNormal ipg_type;
        enum{
            search_type_id = LONG_NO_CUTOFF,
        };
    };
    struct SEARCH_MODE_LONG_CUTOFF{
        typedef TagForceLong force_type;
        typedef TagSearchLongCutoff search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgLongNormal ipg_type;
        enum{
            search_type_id = LONG_CUTOFF,
        };
    };
    struct SEARCH_MODE_GATHER{
        typedef TagForceShort force_type;
        typedef TagSearchShortGather search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgOut ipg_type;
        enum{
            search_type_id = SHORT_GATHER,
        };
    };
    struct SEARCH_MODE_SCATTER{
        typedef TagForceShort force_type;
        typedef TagSearchShortScatter search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgIn ipg_type;
        enum{
            search_type_id = SHORT_SCATTER,
        };
    };
    struct SEARCH_MODE_SYMMETRY{
        typedef TagForceShort force_type;
        typedef TagSearchShortSymmetry search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgInAndOut ipg_type;
        enum{
            search_type_id = SHORT_SYMMETRY,
        };
    };

    // new TAG for P^3T
    struct SEARCH_MODE_LONG_SCATTER{
        typedef TagForceLong force_type;
        typedef TagSearchLongScatter search_type;
        typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
        typedef TagIpgIn ipg_type;
        enum{
            search_type_id = LONG_SCATTER,
        };
    };
    struct SEARCH_MODE_LONG_SYMMETRY{
        typedef TagForceLong force_type;
        typedef TagSearchLongSymmetry search_type;
        typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
        typedef TagIpgInAndOut ipg_type;
        enum{
            search_type_id = LONG_SYMMETRY,
        };
    };
    /*
    struct SEARCH_MODE_LONG_CUTOFF_SCATTER{
        typedef TagForceLong force_type;
        typedef TagSearchLongCutoffScatter search_type;
        typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
        typedef TagIpgInAndOut ipg_type;
        enum{
            search_type_id = LONG_CUTOFF_SCATTER,
        };
    };
    */

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
    typedef MPI_Request MpiRequest;
#else
    typedef int MpiRequest;
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    template<class T> 
    inline MPI_Datatype GetDataType(){
      static MPI_Datatype type = MPI_DATATYPE_NULL;
      if( type == MPI_DATATYPE_NULL ){
          MPI_Type_contiguous(sizeof(T),MPI_BYTE,&type);
          MPI_Type_commit(&type);
      }
      return type;
    };
    template<class T> 
        inline MPI_Datatype GetDataType(const T &){
        return GetDataType<T>();
    }
    template<> inline MPI_Datatype GetDataType<int>(){return MPI_INT;}
    template<> inline MPI_Datatype GetDataType<long>(){return MPI_LONG;}
    //template<> inline MPI_Datatype GetDataType<long long int>(){return MPI_LONG_LONG_INT;}
    template<> inline MPI_Datatype GetDataType<long long int>(){return MPI_LONG_LONG_INT;}
    template<> inline MPI_Datatype GetDataType<unsigned int>(){return MPI_UNSIGNED;}
    template<> inline MPI_Datatype GetDataType<unsigned long>(){return MPI_UNSIGNED_LONG;}
    //template<> inline MPI_Datatype GetDataType<unsigned long long int>(){return MPI_UNSIGNED_LONG;}
    template<> inline MPI_Datatype GetDataType<unsigned long long int>(){return MPI_UNSIGNED_LONG_LONG;}
    template<> inline MPI_Datatype GetDataType<float>(){return MPI_FLOAT;}
    template<> inline MPI_Datatype GetDataType<double>(){return MPI_DOUBLE;}

    template<class Tfloat, class Tint> 
    inline MPI_Datatype GetDataType();
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

        // alltoallv
        void executeV(const ReallocatableArray<T> & val_send,
                      ReallocatableArray<T> & val_recv,
                      const int n_send[],
                      int n_recv[]){
            int cnt = 0;
            //val_recv.reserveAtLeast( val_send.capacity() );
            val_recv.reserveAtLeast( val_send.size() );
            val_recv.clearSize();
            for(int ib=0; ib<n_proc_glb_; ib++){
                n_recv[ib] = n_send[ib];
                for(int ip=0; ip<n_recv[ib]; ip++, cnt++){
                    val_recv.pushBackNoCheck(val_send[cnt]);
                }
            }
            for(int d=DIM_COMM-1; d>=0; d--){
                int radix = n_proc_glb_ / n_proc_1d_[d];
                n_recv_disp_glb_[0] = 0;
                for(int i=0; i<n_proc_glb_; i++){
                    n_recv_disp_glb_[i+1] = n_recv_disp_glb_[i] + n_recv[i];
                }
                val_send_glb_.reserveAtLeast( val_recv.size() );
                val_send_glb_.clearSize();
                for(int ib0=0; ib0<n_proc_glb_; ib0++){
                    int id_send = (ib0 % radix) * n_proc_1d_[d] + ib0 / radix;
                    n_send_glb_[ib0] = n_recv[id_send];
                    int offset = n_recv_disp_glb_[id_send];
                    //val_send_glb_.reserveEmptyAreaAtLeast(n_send_glb_[ib0]);
                    for(int ib1=0; ib1<n_send_glb_[ib0]; ib1++){
                        val_send_glb_.pushBackNoCheck(val_recv[ib1 + offset]);
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
                val_recv.resizeNoInitialize( n_recv_disp_1d_[n_proc_1d_[d]] );
                MPI_Alltoallv(val_send_glb_.getPointer(), n_send_1d_, n_send_disp_1d_, GetDataType<T>(),
                              val_recv.getPointer(), n_recv_1d_, n_recv_disp_1d_, GetDataType<T>(), comm_1d_[d]);
            }
        }
    }; //CommForAllToAll
#endif
    
    class Comm{
    private:
        Comm(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
          MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
          MPI_Comm_size(MPI_COMM_WORLD,&n_proc_);
#else
          rank_ = 0;
          n_proc_ = 1;
#endif
          //#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
          n_thread_ = omp_get_max_threads();
#else
          n_thread_ = 1;
#endif
        }
        ~Comm(){}
        Comm(const Comm & c);
        Comm & operator=(const Comm& c);
        S32 rank_;
        S32 n_proc_;
        S32 n_thread_;
        S32 rank_multi_dim_[DIMENSION];
        S32 n_proc_multi_dim_[DIMENSION];
        //S32 boundary_condition_;
        static Comm & getInstance(){
            static Comm inst;
            return inst;
        }

    template<class T>
    static T allreduceMin(const T & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        T ret;
        T val_tmp = val;
        MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_MIN, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
}

    template<class T>
    static T allreduceMax(const T & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        T ret;
        T val_tmp = val;
        MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_MAX, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<class T>
    static T allreduceSum(const T & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        T ret;
        T val_tmp = val;
        MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_SUM, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }

/*
    template<class Tfloat, class Tint>
    static void allreduceMin(const Tfloat & f_in, const Tint & i_in, Tfloat & f_out, Tint & i_out){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        struct{
            Tfloat x;
            Tint y;
        } loc, glb;
        loc.x = f_in;
        loc.y = i_in;
        MPI_Allreduce(&loc, &glb, 1, GetDataType<Tfloat, Tint>(), MPI_MINLOC, MPI_COMM_WORLD);
        f_out = glb.x;
        i_out = glb.y;
#else
        f_out = f_in;
        i_out = i_in;
#endif
    }
*/
    template<class Tfloat>
    static void allreduceMin(const Tfloat & f_in, const int & i_in, Tfloat & f_out, int & i_out){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        struct{
            Tfloat x;
            int y;
        } loc, glb;
        loc.x = f_in;
        loc.y = i_in;
        MPI_Allreduce(&loc, &glb, 1, GetDataType<Tfloat, int>(), MPI_MINLOC, MPI_COMM_WORLD);
        f_out = glb.x;
        i_out = glb.y;
#else
        f_out = f_in;
        i_out = i_in;
#endif
    }

/*
    template<class Tfloat, class Tint>
    static void allreduceMax(const Tfloat & f_in, const Tint & i_in, Tfloat & f_out, Tint & i_out){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        struct{
            Tfloat x;
            Tint y;
        } loc, glb;
        loc.x = f_in;
        loc.y = i_in;
        MPI_Allreduce(&loc, &glb, 1, GetDataType<Tfloat, Tint>(), MPI_MAXLOC, MPI_COMM_WORLD);
        f_out = glb.x;
        i_out = glb.y;
        //if(Comm::getRank() == 0){std::cout<<"glb.x="<<glb.x<<" glb.y="<<glb.y<<std::endl;}
#else
        f_out = f_in;
        i_out = i_in;
#endif
    }
*/
    template<class Tfloat>
    static void allreduceMax(const Tfloat & f_in, const int & i_in, Tfloat & f_out, int & i_out){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        struct{
            Tfloat x;
            int y;
        } loc, glb;
        loc.x = f_in;
        loc.y = i_in;
        MPI_Allreduce(&loc, &glb, 1, GetDataType<Tfloat, int>(), MPI_MAXLOC, MPI_COMM_WORLD);
        f_out = glb.x;
        i_out = glb.y;
        //if(Comm::getRank() == 0){std::cout<<"glb.x="<<glb.x<<" glb.y="<<glb.y<<std::endl;}
#else
        f_out = f_in;
        i_out = i_in;
#endif
    }
    public:
    static void setNumberOfProcMultiDim(const S32 id, const S32 n) {
        getInstance().n_proc_multi_dim_[id] = n;
    }
    static void setRankMultiDim(const S32 id, const S32 r) {
        getInstance().rank_multi_dim_[id] = r;
    }
        static S32 getRank() { return getInstance().rank_; }
        static S32 getNumberOfProc() { return getInstance().n_proc_; }
        static S32 getRankMultiDim(const S32 id) { return getInstance().rank_multi_dim_[id]; }
        static S32 getNumberOfProcMultiDim(const S32 id) { return getInstance().n_proc_multi_dim_[id]; }
        static S32 getNumberOfThread() { return getInstance().n_thread_; }
        static S32 getThreadNum(){
        //#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
            return omp_get_thread_num();
#else
            return 0;
#endif
        }

        static void barrier(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        }

        static bool synchronizeConditionalBranchAND(const bool & local){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            bool global;
            bool local_tmp = local;
            MPI_Allreduce(&local_tmp, &global, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
            return global;
#else
            return local;
#endif
        }
        static bool synchronizeConditionalBranchOR(const bool & local){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            bool global;
            bool local_tmp = local;
            //std::cerr<<"rank="<<Comm::getRank()<<" local="<<local<<std::endl;
            MPI_Allreduce(&local_tmp, &global, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
            //std::cerr<<"rank="<<Comm::getRank()<<" global="<<global<<std::endl;
            return global;
#else
            return local;
#endif
        }

        ///////////////////////////
        // MPI ALLREDUCE WRAPPER //
        template<class T> static inline T getMinValue(const T & val);
        template<class T> static inline T getMaxValue(const T & val);
        //template<class Tfloat, class Tint> static inline void getMinValue(const Tfloat & f_in, const Tint & i_in, Tfloat & f_out, Tint & i_out);
        //template<class Tfloat, class Tint> static inline void getMaxValue(const Tfloat & f_in, const Tint & i_in, Tfloat & f_out, Tint & i_out);
        template<class Tfloat> static inline void getMinValue(const Tfloat & f_in, const int & i_in, Tfloat & f_out, int & i_out);
        template<class Tfloat> static inline void getMaxValue(const Tfloat & f_in, const int & i_in, Tfloat & f_out, int & i_out);
        template<class T> static inline T getSum(const T & val);

        ///////////////////////
        // MPI BCAST WRAPPER //
        // new functions 10 Feb 2015
        template<class T>
        static inline void broadcast(T * val, const int n, const int src=0){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
            MPI_Bcast(val, n, GetDataType<T>(), src, MPI_COMM_WORLD);
#else
            // NOP
#endif
        }

        ///////////////////////////
        // MPI GATHER WRAPPER //
        template<class T> static inline void gather(T * val_send,
                                                    int n,
                                                    T * val_recv){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Gather(val_send, n, GetDataType<T>(),
                       val_recv, n, GetDataType<T>(), 0, MPI_COMM_WORLD);
#else
            for(int i=0; i<n; i++)val_recv[i] = val_send[i];
#endif
        }

        template<class T> static inline void gatherV(T * val_send,
                                                     int n_send,
                                                     T * val_recv){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
            const int n_proc = Comm::getNumberOfProc();
            int * n_recv      = new int[n_proc];
            int * n_disp_recv = new int[n_proc+1];
            gather(&n_send, 1, n_recv);
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            MPI_Gatherv(val_send, n_send, GetDataType<T>(),
                        val_recv, n_recv, n_disp_recv, GetDataType<T>(), 0, MPI_COMM_WORLD);
            delete[] n_recv;
            delete[] n_disp_recv;
#else
            int n = n_send;
            for(int i=0; i<n; i++)val_recv[i] = val_send[i];
#endif
        }

        template<class T> static inline void gatherV(T * val_send, // in
                                                     int n_send,   // in
                                                     T * val_recv, // out
                                                     int * n_recv, // in
                                                     int * n_recv_disp // in
                                                     ){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Gatherv(val_send, n_send, GetDataType<T>(),
                        val_recv, n_recv, n_recv_disp, GetDataType<T>(), 0, MPI_COMM_WORLD);
#else
            for(int i=0; i<n_send; i++) val_recv[i] = val_send[i];
            //n_recv[0] = n_recv_disp[1] = n_send; //not needed ?
            //n_recv_disp[0] = 0; //not needed ?
#endif
        }

        ///////////////////////////
        // MPI SCATTER WRAPPER //
        template<class T> static inline void scatterV(T * val_send,
                                                      int * n_send,
                                                      int * n_send_disp,
                                                      T * val_recv, // output
                                                      int n_recv,
                                                      int root=0){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Scatterv(val_send, n_send, n_send_disp, GetDataType<T>(),
                         val_recv, n_recv, GetDataType<T>(), 
                         root, MPI_COMM_WORLD);
#else
            S32 n_proc = Comm::getNumberOfProc();
            for(int i=0; i<n_send_disp[n_proc]; i++) val_recv[i] = val_send[i];
#endif
        }



        ///////////////////////////
        // MPI ALLGATHER WRAPPER //
        template<class T> static inline void allGather(const T * val_send, // in
                                                       const int n,        // in
                                                       T * val_recv  // out
                                                       ){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            T * val_send_tmp = const_cast<T*>(val_send);
            MPI_Allgather(val_send_tmp, n, GetDataType<T>(),
                          val_recv, n, GetDataType<T>(), 
                          MPI_COMM_WORLD);
#else
            for(int i=0; i<n; i++)val_recv[i] = val_send[i];
#endif
        }

        template<class T> static inline void allGatherV(const T * val_send, // in
                                                        const int n_send,   // in
                                                        T * val_recv, // out
                                                        int * n_recv, // in
                                                        int * n_recv_disp //in
                                                        ){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            T * val_send_tmp = const_cast<T*>(val_send);
            int * n_recv_tmp = const_cast<int*>(n_recv);
            int * n_recv_disp_tmp = const_cast<int*>(n_recv_disp);
            MPI_Allgatherv(val_send_tmp, n_send, GetDataType<T>(),
                           val_recv, n_recv_tmp, n_recv_disp_tmp, GetDataType<T>(), 
                           MPI_COMM_WORLD);
#else
            for(int i=0; i<n_send; i++) val_recv[i] = val_send[i];
            n_recv[0] = n_recv_disp[1] = n_send;
            n_recv_disp[0] = 0;
#endif
        }



        ///////////////////////////
        // MPI ALLTOALL WRAPPER //
        template<class T>
        static inline void allToAll(const T * val_send, const int n, T * val_recv){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            T * val_send_tmp = const_cast<T*>(val_send);
            MPI_Alltoall(val_send_tmp, n, GetDataType<T>(),
                         val_recv, n, GetDataType<T>(), MPI_COMM_WORLD);
#else
            for(int i=0; i<n; i++) val_recv[i] = val_send[i];
#endif
        }
        template<class T>
        static inline void allToAllV(const T * val_send,
                                     const int * n_send,
                                     const int * n_send_disp,
                                     T * val_recv,
                                     int * n_recv,
                                     int * n_recv_disp){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            T * val_send_tmp = const_cast<T*>(val_send);
            int * n_send_tmp = const_cast<int*>(n_send);
            int * n_send_disp_tmp = const_cast<int*>(n_send_disp);
            int * n_recv_tmp = const_cast<int*>(n_recv);
            int * n_recv_disp_tmp = const_cast<int*>(n_recv_disp);
            MPI_Alltoallv(val_send_tmp, n_send_tmp, n_send_disp_tmp, GetDataType<T>(),
                          val_recv,     n_recv_tmp, n_recv_disp_tmp, GetDataType<T>(),
                          MPI_COMM_WORLD);
#else
            for(int i=0; i<n_send[0]; i++) val_recv[i] = val_send[i];
            n_recv[0] = n_send[0];
            n_recv_disp[0] = n_send_disp[0];
            n_recv_disp[1] = n_send_disp[1];
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
            std::cerr << "     ||     Version 5.0 (2018/11)        ||"   << std::endl;
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




    template<> inline float Comm::getMinValue<float>(const float & val){
        return allreduceMin(val);
    }
    template<> inline double Comm::getMinValue<double>(const double & val){
        return allreduceMin(val);
    }
    template<> inline int Comm::getMinValue<int>(const int & val){
        return allreduceMin(val);
    }
    template<> inline long Comm::getMinValue<long>(const long & val){
        return allreduceMin(val);
    }
    template<> inline long long Comm::getMinValue<long long>(const long long & val){
        return allreduceMin(val);
    }
    template<> inline unsigned int Comm::getMinValue<unsigned int>(const unsigned int & val){
        return allreduceMin(val);
    }    
    template<> inline unsigned long Comm::getMinValue<unsigned long>(const unsigned long & val){
        return allreduceMin(val);
    }    
    template<> inline unsigned long long Comm::getMinValue<unsigned long long>(const unsigned long long & val){
        return allreduceMin(val);
    }    
    template<> inline
    Vector2<float> Comm::getMinValue< Vector2<float> >(const Vector2<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector2<float> ret;
        MPI_Allreduce((float*)&val.x, (float*)&ret.x, 2, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline Vector2<double> Comm::getMinValue< Vector2<double> >(const Vector2<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector2<double> ret;
        MPI_Allreduce((double*)&val.x, (double*)&ret.x, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector3<float> Comm::getMinValue< Vector3<float> >(const Vector3<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector3<float> ret;
        MPI_Allreduce((float*)&val.x, (float*)&ret.x, 3, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector3<double> Comm::getMinValue< Vector3<double> >(const Vector3<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector3<double> ret;
        MPI_Allreduce((double*)&val.x, (double*)&ret.x, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }

    template<> inline float Comm::getMaxValue<float>(const float & val){
        return allreduceMax(val);
    }
    template<> inline double Comm::getMaxValue<double>(const double & val){
        return allreduceMax(val);
    }
    template<> inline int Comm::getMaxValue<int>(const int & val){
        return allreduceMax(val);
    }
    template<> inline long Comm::getMaxValue<long>(const long & val){
        return allreduceMax(val);
    }
    template<> inline long long Comm::getMaxValue<long long>(const long long & val){
        return allreduceMax(val);
    }
    template<> inline unsigned int Comm::getMaxValue<unsigned int>(const unsigned int & val){
        return allreduceMax(val);
    }    
    template<> inline unsigned long Comm::getMaxValue<unsigned long>(const unsigned long & val){
        return allreduceMax(val);
    }    
    template<> inline unsigned long long Comm::getMaxValue<unsigned long long>(const unsigned long long & val){
        return allreduceMax(val);
    }    
    template<> inline
    Vector2<float> Comm::getMaxValue< Vector2<float> >(const Vector2<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector2<float> ret;
        MPI_Allreduce((float*)&val.x, (float*)&ret.x, 2, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector2<double> Comm::getMaxValue< Vector2<double> >(const Vector2<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector2<double> ret;
        MPI_Allreduce((double*)&val.x, (double*)&ret.x, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector3<float> Comm::getMaxValue< Vector3<float> >(const Vector3<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector3<float> ret;
        MPI_Allreduce((float*)&val.x, (float*)&ret.x, 3, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector3<double> Comm::getMaxValue< Vector3<double> >(const Vector3<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector3<double> ret;
        MPI_Allreduce((double*)&val.x, (double*)&ret.x, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    
    template<> inline float Comm::getSum(const float & val){
        return allreduceSum(val);
    }
    template<> inline double Comm::getSum(const double & val){
        return allreduceSum(val);
    }
    template<> inline int Comm::getSum(const int & val){
        return allreduceSum(val);
    }
    template<> inline unsigned int Comm::getSum(const unsigned int & val){
        return allreduceSum(val);
    }
    template<> inline long Comm::getSum(const long & val){
        return allreduceSum(val);
    }
    template<> inline unsigned long Comm::getSum(const unsigned long & val){
        return allreduceSum(val);
    }
    template<> inline long long int Comm::getSum(const long long int & val){
        return allreduceSum(val);
    }
    template<> inline unsigned long long int Comm::getSum(const unsigned long long int & val){
        return allreduceSum(val);
    }
    template<> inline Vector3<float> Comm::getSum(const Vector3<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        Vector3<float> ret;
        MPI_Allreduce((float*)&val, (float*)&ret, 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline Vector3<double> Comm::getSum(const Vector3<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        Vector3<double> ret;
        MPI_Allreduce((double*)&val, (double*)&ret, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline Vector2<float> Comm::getSum(const Vector2<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        Vector2<float> ret;
        MPI_Allreduce((float*)&val, (float*)&ret, 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }
    template<> inline Vector2<double> Comm::getSum(const Vector2<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        Vector2<double> ret;
        MPI_Allreduce((double*)&val, (double*)&ret, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return ret;
#else
        return val;
#endif
    }

/*
    template<> inline void Comm::getMinValue<float, int>(const float & f_in, const int & i_in, 
                                                         float & f_out, int & i_out){
        allreduceMin(f_in, i_in, f_out, i_out);
    }
    template<> inline void Comm::getMinValue<double, int>(const double & f_in, const int & i_in, 
                                                          double & f_out, int & i_out){
        allreduceMin(f_in, i_in, f_out, i_out);
    }
*/
    template<> inline void Comm::getMinValue<float>(const float & f_in, const int & i_in, 
                                                    float & f_out, int & i_out){
        allreduceMin(f_in, i_in, f_out, i_out);
    }
    template<> inline void Comm::getMinValue<double>(const double & f_in, const int & i_in, 
                                                     double & f_out, int & i_out){
        allreduceMin(f_in, i_in, f_out, i_out);
    }
/*
    template<> inline
    void Comm::getMaxValue<float, int>(const float & f_in, const int & i_in, 
                                       float & f_out, int & i_out){
        allreduceMax(f_in, i_in, f_out, i_out);
    }

    template<> inline
    void Comm::getMaxValue<double, int>(const double & f_in, const int & i_in, 
                                        double & f_out, int & i_out){
        allreduceMax(f_in, i_in, f_out, i_out);
    }
*/
    template<> inline
    void Comm::getMaxValue<float>(const float & f_in, const int & i_in, 
                                  float & f_out, int & i_out){
        allreduceMax(f_in, i_in, f_out, i_out);
    }

    template<> inline
    void Comm::getMaxValue<double>(const double & f_in, const int & i_in, 
                                   double & f_out, int & i_out){
        allreduceMax(f_in, i_in, f_out, i_out);
    }

    static const S32 N_SMP_PTCL_TOT_PER_PSYS_DEFAULT = 1000000;

    ////////
    // util 
    inline F32 CalcSeparationSQPointToBox(const F32vec & point, 
                                          const F32vec & center, 
                                          const F32vec & half_length ){
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        F32 dx = fabs(point.x - center.x);
        F32 dy = fabs(point.y - center.y);
        dx = ( half_length.x < dx ) ? (dx - half_length.x) : 0.0;
        dy = ( half_length.y < dy ) ? (dy - half_length.y) : 0.0;
        return dx*dx + dy*dy;
#else
        F32 dx = fabs(point.x - center.x);
        F32 dy = fabs(point.y - center.y);
        F32 dz = fabs(point.z - center.z);
        dx = ( half_length.x < dx ) ? (dx - half_length.x) : 0.0;
        dy = ( half_length.y < dy ) ? (dy - half_length.y) : 0.0;
        dz = ( half_length.z < dz ) ? (dz - half_length.z) : 0.0;
        return dx*dx + dy*dy + dz*dz;
#endif
    }

    // for check function
    inline bool IsInBox(const F32vec & pos,
                        const F32vec & center,
                        const F32 half_length,
                        const F32 tolerance=1e-6){
        const F32 tol = -std::abs(tolerance);
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

/*
    template<class T>
    inline T GetMSB(const T val);
    template<>
    inline unsigned long long int GetMSB(const unsigned long long int val){
        return (val>>63) & 0x1;
    }
    template<>
    inline unsigned int GetMSB(const unsigned int val){
        return (val>>31) & 0x1;
    }
*/

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
/*
    template<class T>
    inline T ClearMSB(const T val);
    template<>
    inline unsigned long long int ClearMSB(const unsigned long long int val){
        return val & 0x7fffffffffffffff;
    }
    template<>
    inline unsigned int ClearMSB(const unsigned int val){
        return val & 0x7fffffff;
    }
*/

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

/*
    template<class T>
    inline T SetMSB(const T val);
    template<>
    inline unsigned long long int SetMSB(const unsigned long long int val){
        return val | 0x8000000000000000;
    }
    template<>
    inline unsigned int SetMSB(const unsigned int val){
        return val | 0x80000000;
    }
*/

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
        //std::cerr<<"box_loc="<<box_loc<<std::endl;
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
        F64 exchange_LET_1st__icomm_sp;
        F64 exchange_LET_1st__a2a_sp;
        F64 exchange_LET_1st__icomm_ep;
        F64 exchange_LET_1st__a2a_ep;

        F64 add_moment_as_sp_local;
        F64 add_moment_as_sp_global;

        void dump(std::ostream & fout=std::cout) const {
            fout<<"collect_sample_particle= "<<collect_sample_particle<<std::endl;
            fout<<"decompose_domain= "<<decompose_domain<<std::endl;
            fout<<"exchange_particle= "<<exchange_particle<<std::endl;
            fout<<"set_particle_local_tree= "<<set_particle_local_tree<<std::endl;
            fout<<"set_particle_global_tree= "<<set_particle_global_tree<<std::endl;
            fout<<"make_local_tree= "<<make_local_tree<<std::endl;
            fout<<"make_global_tree= "<<make_global_tree<<std::endl;
            fout<<"set_root_cell= "<<set_root_cell<<std::endl;
            fout<<"calc_force= "<<calc_force<<std::endl;
            fout<<"calc_moment_local_tree= "<<calc_moment_local_tree<<std::endl;
            fout<<"calc_moment_global_tree= "<<calc_moment_global_tree<<std::endl;
            fout<<"make_LET_1st= "<<make_LET_1st<<std::endl;
            fout<<"make_LET_2nd= "<<make_LET_2nd<<std::endl;
            fout<<"exchange_LET_1st= "<<exchange_LET_1st<<std::endl;
            fout<<"exchange_LET_2nd= "<<exchange_LET_2nd<<std::endl;
            fout<<"write_back= "<<write_back<<std::endl;
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

            decompose_domain__sort_particle_1st = decompose_domain__sort_particle_2nd = decompose_domain__sort_particle_3rd = decompose_domain__gather_particle = 0.0;
            decompose_domain__setup = decompose_domain__determine_coord_1st = decompose_domain__migrae_particle_1st = decompose_domain__determine_coord_2nd 
                = decompose_domain__determine_coord_3rd = decompose_domain__exchange_pos_domain = 0.0;
            exchange_LET_1st__a2a_n = exchange_LET_1st__a2a_sp = exchange_LET_1st__icomm_ep = exchange_LET_1st__icomm_sp = exchange_LET_1st__a2a_ep = 0.0;

            add_moment_as_sp_local = add_moment_as_sp_global = 0.0;
            write_back = 0.0;
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
            ret.exchange_LET_1st__icomm_ep   = this->exchange_LET_1st__icomm_ep   + rhs.exchange_LET_1st__icomm_ep;
            ret.exchange_LET_1st__a2a_sp   = this->exchange_LET_1st__a2a_sp   + rhs.exchange_LET_1st__a2a_sp;
            ret.exchange_LET_1st__icomm_sp = this->exchange_LET_1st__icomm_sp + rhs.exchange_LET_1st__icomm_sp;
            ret.exchange_LET_1st__a2a_ep   = this->exchange_LET_1st__a2a_ep   + rhs.exchange_LET_1st__a2a_ep;

            ret.add_moment_as_sp_local = this->add_moment_as_sp_local + rhs.add_moment_as_sp_local;
            ret.add_moment_as_sp_global = this->add_moment_as_sp_global + rhs.add_moment_as_sp_global;

            ret.write_back = this->write_back + rhs.write_back;
            return ret;
        }
        F64 getTotalTime() const {
            /*
            return collect_sample_particle + decompose_domain + exchange_particle + set_particle_local_tree + set_particle_global_tree + make_local_tree + make_global_tree + set_root_cell
                + calc_force + calc_moment_local_tree + calc_moment_global_tree + make_LET_1st + make_LET_2nd + exchange_LET_1st + exchange_LET_2nd;
            */
            return collect_sample_particle + decompose_domain + exchange_particle + set_particle_local_tree + set_particle_global_tree + make_local_tree + make_global_tree + set_root_cell
                + calc_force + calc_moment_local_tree + calc_moment_global_tree + make_LET_1st + make_LET_2nd + exchange_LET_1st + exchange_LET_2nd
                + morton_sort_local_tree + link_cell_local_tree 
                + morton_sort_global_tree + link_cell_global_tree
                + add_moment_as_sp_local + add_moment_as_sp_global
                + write_back;
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


            decompose_domain__sort_particle_1st = decompose_domain__sort_particle_2nd = decompose_domain__sort_particle_3rd = decompose_domain__gather_particle = 0.0;
            decompose_domain__setup = decompose_domain__determine_coord_1st = decompose_domain__migrae_particle_1st = decompose_domain__determine_coord_2nd 
                = decompose_domain__determine_coord_3rd = decompose_domain__exchange_pos_domain = 0.0;
            exchange_LET_1st__a2a_n = exchange_LET_1st__a2a_sp = exchange_LET_1st__icomm_ep = exchange_LET_1st__icomm_sp = exchange_LET_1st__a2a_ep = 0.0;

            add_moment_as_sp_local = add_moment_as_sp_global = 0.0;
            write_back = 0.0;
        }
    };

    inline F64 GetDistanceMinSq(const F64ort & pos0,
                                const F64ort & pos1,
                                const F64vec & len_peri){
        const F64vec cen0 = pos0.getCenter();
        const F64vec cen1 = pos1.getCenter();
        const F64vec len = pos0.getHalfLength() + pos1.getHalfLength();
        F64 dis = 0.0;
        
        F64 dx = fabs(cen0.x-cen1.x);
        dx = (dx < (len_peri.x-dx)) ? dx : (len_peri.x-dx);
        dx = (len.x < dx) ? dx-len.x : 0.0;

        F64 dy = fabs(cen0.y-cen1.y);
        dy = (dy < (len_peri.y-dy)) ? dy : (len_peri.y-dy);
        dy = (len.y < dy) ? dy-len.y : 0.0;
        dis = dx*dx + dy*dy;

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        F64 dz = fabs(cen0.z-cen1.z);
        dz = (dz < (len_peri.z-dz)) ? dz : (len_peri.z-dz);
        dz = (len.z < dz) ? dz-len.z : 0.0;
        dis += dz*dz;
#endif
        
        return dis;
    }
    
    inline F64 GetDistanceMinSq(const F64ort & pos0,
                                const F64vec & pos1,
                                const F64vec & len_peri){
        const F64vec cen0 = pos0.getCenter();
        const F64vec len  = pos0.getHalfLength();
        F64 dis = 0.0;
        
        F64 dx = fabs(cen0.x-pos1.x);
        dx = (dx < (len_peri.x-dx)) ? dx : (len_peri.x-dx);
        dx = (len.x < dx) ? dx-len.x : 0.0;
        
        F64 dy = fabs(cen0.y-pos1.y);
        dy = (dy < (len_peri.y-dy)) ? dy : (len_peri.y-dy);
        dy = (len.y < dy) ? dy-len.y : 0.0;
        
        dis = dx*dx + dy*dy;
        
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        F64 dz = fabs(cen0.z-pos1.z);
        dz = (dz < (len_peri.z-dz)) ? dz : (len_peri.z-dz);
        dz = (len.z < dz) ? dz-len.z : 0.0;
        dis += dz*dz;
#endif
        return dis;
    }


    inline F64 GetDistanceMinSq(const F64ort & pos0,
                                const F64ort & pos1){
        const F64vec cen0 = pos0.getCenter();
        const F64vec cen1 = pos1.getCenter();
        const F64vec len = pos0.getHalfLength() + pos1.getHalfLength();
        F64 dis = 0.0;
        
        F64 dx = fabs(cen0.x-cen1.x);
        dx = (len.x < dx) ? dx-len.x : 0.0;

        F64 dy = fabs(cen0.y-cen1.y);
        dy = (len.y < dy) ? dy-len.y : 0.0;
        dis = dx*dx + dy*dy;

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        F64 dz = fabs(cen0.z-cen1.z);
        dz = (len.z < dz) ? dz-len.z : 0.0;
        dis += dz*dz;
#endif
        
        return dis;
    }
    
    inline F64 GetDistanceMinSq(const F64ort & pos0,
                                const F64vec & pos1){
        const F64vec cen0 = pos0.getCenter();
        const F64vec len  = pos0.getHalfLength();
        F64 dis = 0.0;
        
        F64 dx = fabs(cen0.x-pos1.x);
        dx = (len.x < dx) ? dx-len.x : 0.0;
        
        F64 dy = fabs(cen0.y-pos1.y);
        dy = (len.y < dy) ? dy-len.y : 0.0;
        
        dis = dx*dx + dy*dy;
        
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        F64 dz = fabs(cen0.z-pos1.z);
        dz = (len.z < dz) ? dz-len.z : 0.0;
        dis += dz*dz;
#endif
        return dis;
    }

    template<class T>
    inline std::string GetBinString(const T val);
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
}

#include"util.hpp"

#include"timer.hpp"
