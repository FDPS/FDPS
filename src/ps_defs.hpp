#pragma once

#include<cmath>
#include<vector>
#include<algorithm>
#include<exception>
#include<stdexcept>
#include<cassert>
#include<typeinfo>
#include<cstdio>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include"mpi.h"
#endif
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include<omp.h>
#endif

namespace ParticleSimulator{
  //static const S64 LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1<<31;
}

//#include"reallocatable_array.hpp"
#include"vector2.hpp"
#include"vector3.hpp"
#include"orthotope2.hpp"
#include"orthotope3.hpp"
#include"matrix_sym2.hpp"
#include"matrix_sym3.hpp"


template<bool> struct CompileTimeError;
template<> struct CompileTimeError<true> {};
#define STATIC_ASSERT(expr, msg) \
{ CompileTimeError<expr> PS_ERROR_##msg; (void)ERROR_##msg; }

#define PARTICLE_SIMULATOR_STATIC_ASSERT(expr) \
{ CompileTimeError<expr> ERROR; }

#define PARTICLE_SIMULATOR_PRINT_ERROR(msg) \
    { std::cout<<"PS_ERROR: "<<msg<<" \n"<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl; }

#define PARTICLE_SIMULATOR_PRINT_LINE_INFO() \
    { std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl; }

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

    static const S64 LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1<<31;

    enum SEARCH_MODE{
        LONG_NO_CUTOFF,
        LONG_CUTOFF,
        SHORT_GATHER,
        SHORT_SCATTER,
        SHORT_SYMMETRY,
    };
    struct TagForceLong{};
    struct TagForceShort{};
    struct TagSearchLong{};
    struct TagSearchLongCutoff{};
    struct TagSearchShortGather{};
    struct TagSearchShortScatter{};
    struct TagSearchShortSymmetry{};
    struct SEARCH_MODE_LONG{
        typedef TagForceLong force_type;
        typedef TagSearchLong search_type;
        enum{
            search_type_id = LONG_NO_CUTOFF,
        };
    };
    struct SEARCH_MODE_LONG_CUTOFF{
        typedef TagForceLong force_type;
        typedef TagSearchLongCutoff search_type;
        enum{
            search_type_id = LONG_CUTOFF,
        };
    };
    struct SEARCH_MODE_GATHER{
        typedef TagForceShort force_type;
        typedef TagSearchShortGather search_type;
        enum{
            search_type_id = SHORT_GATHER,
        };
    };
    struct SEARCH_MODE_SCATTER{
        typedef TagForceShort force_type;
        typedef TagSearchShortScatter search_type;
        enum{
            search_type_id = SHORT_SCATTER,
        };
    };
    struct SEARCH_MODE_SYMMETRY{
        typedef TagForceShort force_type;
        typedef TagSearchShortSymmetry search_type;
        enum{
            search_type_id = SHORT_SYMMETRY,
        };
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
    typedef MPI::Request MpiRequest;
#else
    typedef int MpiRequest;
#endif
    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    template<class T> 
    inline MPI::Datatype GetDataType(){
        static MPI::Datatype type;
        if( type == MPI::Datatype() ){
            type = MPI::BYTE.Create_contiguous(sizeof(T));
            type.Commit();
        }
        return type;
    };
    template<class T> 
    inline MPI::Datatype GetDataType(const T &){
        return GetDataType<T>();
    }
    template<> inline MPI::Datatype GetDataType<int>(){return MPI::INT;}
    template<> inline MPI::Datatype GetDataType<long>(){return MPI::LONG;}
    template<> inline MPI::Datatype GetDataType<long long int>(){return MPI::LONG_LONG_INT;}
    template<> inline MPI::Datatype GetDataType<unsigned int>(){return MPI::UNSIGNED;}
    template<> inline MPI::Datatype GetDataType<unsigned long>(){return MPI::UNSIGNED_LONG;}
    template<> inline MPI::Datatype GetDataType<float>(){return MPI::FLOAT;}
    template<> inline MPI::Datatype GetDataType<double>(){return MPI::DOUBLE;}

    template<class Tfloat, class Tint> 
    inline MPI::Datatype GetDataType();
    template<> inline MPI::Datatype GetDataType<float, int>(){return MPI::FLOAT_INT;}
    template<> inline MPI::Datatype GetDataType<double, int>(){return MPI::DOUBLE_INT;}
#endif


    class Comm{
    private:
        Comm(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            rank_ = MPI::COMM_WORLD.Get_rank();
            n_proc_ = MPI::COMM_WORLD.Get_size();
#else
            rank_ = 0;
            n_proc_ = 1;
#endif
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
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
	    MPI::COMM_WORLD.Allreduce(&val, &ret, 1, GetDataType<T>(), MPI::MIN);
	    return ret;
#else
	    return val;
#endif
	}

	template<class T>
	static T allreduceMax(const T & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    T ret;
	    MPI::COMM_WORLD.Allreduce(&val, &ret, 1, GetDataType<T>(), MPI::MAX);
	    return ret;
#else
	    return val;
#endif
	}
	template<class T>
	static T allreduceSum(const T & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    T ret;
	    MPI::COMM_WORLD.Allreduce(&val, &ret, 1, GetDataType<T>(), MPI::SUM);
	    return ret;
#else
	    return val;
#endif
	}

	template<class Tfloat, class Tint>
	static void allreduceMin(const Tfloat & f_in, const Tint & i_in, Tfloat & f_out, Tint & i_out){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    struct{
            Tfloat x;
            Tint y;
	    } loc, glb;
	    loc.x = f_in;
	    loc.y = i_in;
	    MPI::COMM_WORLD.Allreduce(&loc, &glb, 1, GetDataType<Tfloat, Tint>(), MPI::MINLOC);
	    f_out = glb.x;
	    i_out = glb.y;
#else
	    f_out = f_in;
	    i_out = i_in;
#endif
	}
	template<class Tfloat, class Tint>
	static void allreduceMax(const Tfloat & f_in, const Tint & i_in, Tfloat & f_out, Tint & i_out){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    struct{
            Tfloat x;
            Tint y;
	    } loc, glb;
	    loc.x = f_in;
	    loc.y = i_in;
	    MPI::COMM_WORLD.Allreduce(&loc, &glb, 1, GetDataType<Tfloat, Tint>(), MPI::MAXLOC);
	    f_out = glb.x;
	    i_out = glb.y;
#else
	    f_out = f_in;
	    i_out = i_in;
#endif
	}

    public:
        static S32 getRank() { return getInstance().rank_; }
        static S32 getNumberOfProc() { return getInstance().n_proc_; }
        static S32 getRankMultiDim(const S32 id) { return getInstance().rank_multi_dim_[id]; }
        static S32 getNumberOfProcMultiDim(const S32 id) { return getInstance().n_proc_multi_dim_[id]; }
        static S32 getNumberOfThread() { return getInstance().n_thread_; }
        static S32 getThreadNum(){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            return omp_get_thread_num();
#else
            return 0;
#endif
        }

        static void barrier(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI::COMM_WORLD.Barrier();
#endif
        }
	
        static bool synchronizeConditionalBranchAND(const bool & local){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            bool global;
            MPI::COMM_WORLD.Allreduce(&local, &global, 1, MPI::BOOL, MPI::LAND);
            return global;
#else
            return local;
#endif
        }
        static bool synchronizeConditionalBranchOR(const bool & local){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            bool global;
            //std::cerr<<"rank="<<Comm::getRank()<<" local="<<local<<std::endl;
            MPI::COMM_WORLD.Allreduce(&local, &global, 1, MPI::BOOL, MPI::LOR);
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
        template<class Tfloat, class Tint> static inline void getMinValue(const Tfloat & f_in, const Tint & i_in, Tfloat & f_out, Tint & i_out);
        template<class Tfloat, class Tint> static inline void getMaxValue(const Tfloat & f_in, const Tint & i_in, Tfloat & f_out, Tint & i_out);
        template<class T> static inline T getSum(const T & val);

        ///////////////////////
        // MPI BCAST WRAPPER //
        // new functions 10 Feb 2015
        template<class T>
        static inline void broadcast(T * val, const int n, const int src=0){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
            MPI::COMM_WORLD.Bcast(val, n, GetDataType<T>(), src);
#else
            // NOP
#endif
        }

        ///////////////////////////
        // MPI ALLGATHER WRAPPER //
        template<class T> static inline void allGather(const T * val_send,
                                                       const int n,
                                                       T * val_recv){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
            MPI::COMM_WORLD.Allgather(val_send, n, GetDataType<T>(),
                                      val_recv, n, GetDataType<T>());
#else
            for(int i=0; i<n; i++)val_recv[i] = val_send[i];
#endif	    
        }

        template<class T> static inline void allGatherV(const T * val_send,
                                                        const int n_send,
                                                        T * val_recv,
                                                        int * n_recv,
                                                        int * n_recv_disp){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
            MPI::COMM_WORLD.Allgatherv(val_send, n_send, GetDataType<T>(),
                                       val_recv, n_recv, n_recv_disp, GetDataType<T>());
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
            MPI::COMM_WORLD.Alltoall(val_send, n, GetDataType<T>(),
                                     val_recv, n, GetDataType<T>());
#else
            for(int i=0; i<n; i++) val_recv[i] = val_send[i];
#endif
        }
        template<class T>
        static inline void allToAllV(const T * val_send, const int * n_send, const int * n_send_disp,
                                     T * val_recv,       int * n_recv,       int * n_recv_disp){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI::COMM_WORLD.Alltoallv(val_send, n_send, n_send_disp, GetDataType<T>(),
                                      val_recv, n_recv, n_recv_disp, GetDataType<T>());
#else
            for(int i=0; i<n_send[0]; i++) val_recv[i] = val_send[i];
            n_recv[0] = n_send[0];
            n_recv_disp[0] = n_send_disp[0];
            n_recv_disp[1] = n_send_disp[1];
#endif
        }
	
    };

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
    template<> inline
    Vector2<float> Comm::getMinValue< Vector2<float> >(const Vector2<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector2<float> ret;
        MPI::COMM_WORLD.Allreduce((float*)&val[0], (float*)&ret[0], 2, MPI::FLOAT, MPI::MIN);
        return ret;
#else
        return val;
#endif
    }
    template<> inline Vector2<double> Comm::getMinValue< Vector2<double> >(const Vector2<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector2<double> ret;
        MPI::COMM_WORLD.Allreduce((double*)&val[0], (double*)&ret[0], 2, MPI::DOUBLE, MPI::MIN);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector3<float> Comm::getMinValue< Vector3<float> >(const Vector3<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector3<float> ret;
        MPI::COMM_WORLD.Allreduce((float*)&val[0], (float*)&ret[0], 3, MPI::FLOAT, MPI::MIN);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector3<double> Comm::getMinValue< Vector3<double> >(const Vector3<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector3<double> ret;
        MPI::COMM_WORLD.Allreduce((double*)&val[0], (double*)&ret[0], 3, MPI::DOUBLE, MPI::MIN);
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

    template<> inline
    Vector2<float> Comm::getMaxValue< Vector2<float> >(const Vector2<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector2<float> ret;
        MPI::COMM_WORLD.Allreduce((float*)&val[0], (float*)&ret[0], 2, MPI::FLOAT, MPI::MAX);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector2<double> Comm::getMaxValue< Vector2<double> >(const Vector2<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector2<double> ret;
        MPI::COMM_WORLD.Allreduce((double*)&val[0], (double*)&ret[0], 2, MPI::DOUBLE, MPI::MAX);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector3<float> Comm::getMaxValue< Vector3<float> >(const Vector3<float> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector3<float> ret;
        MPI::COMM_WORLD.Allreduce((float*)&val[0], (float*)&ret[0], 3, MPI::FLOAT, MPI::MAX);
        return ret;
#else
        return val;
#endif
    }
    template<> inline
    Vector3<double> Comm::getMaxValue< Vector3<double> >(const Vector3<double> & val){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL	
        Vector3<double> ret;
        MPI::COMM_WORLD.Allreduce((double*)&val[0], (double*)&ret[0], 3, MPI::DOUBLE, MPI::MAX);
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
    template<> inline long Comm::getSum(const long & val){
        return allreduceSum(val);
    }
    template<> inline long long int Comm::getSum(const long long int & val){
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


    template<> inline void Comm::getMinValue<float, int>(const float & f_in, const int & i_in, 
                                                         float & f_out, int & i_out){
        allreduceMin(f_in, i_in, f_out, i_out);
    }
    template<> inline void Comm::getMinValue<double, int>(const double & f_in, const int & i_in, 
                                                          double & f_out, int & i_out){
        allreduceMin(f_in, i_in, f_out, i_out);
    }

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

    static inline void Initialize(int &argc, char **&argv){

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#ifndef MPI_VERSION
      CompileTimeError<false> PLEASE_COMPILE_USING_MPI;
#endif //MPI_VERSION
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#ifndef _OPENMP
        CompileTimeError<false> OPENMP_FLAG_IMCOMPATIBLE;
#endif //_OPENMP
#else  //PARTICLE_SIMULATOR_THREAD_PARALLEL
#ifdef _OPENMP
        CompileTimeError<false> OPENMP_FLAG_IMCOMPATIBLE;
#endif //_OPENMP
#endif //PARTICLE_SIMULATOR_THREAD_PARALLEL

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI::Init(argc, argv);
#endif
    }

    static inline void Finalize(){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if(Comm::getRank() == 0) {
            fprintf(stderr, "******** FDPS has successfully finished. ********\n");
        }
        MPI::Finalize();
#else
        fprintf(stderr, "******** FDPS has successfully finished. ********\n");
#endif
    }

    static inline void Abort(const S32 err = -1){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI::COMM_WORLD.Abort(err);
#else
        exit(err);
#endif
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
#pragma omp parallel
	{
	    F64ort box_loc_tmp;
	    box_loc_tmp.init();
#pragma omp for nowait
	    for(S32 ip=0; ip<n; ip++){
		box_loc_tmp.merge(ptcl[ip].getPos());
	    }
#pragma omp critical
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
#pragma omp parallel
        {
            F64ort box_loc_tmp;
            box_loc_tmp.init();
#pragma omp for nowait
            for(S32 ip=0; ip<n; ip++){
                box_loc_tmp.merge(ptcl[ip].getPos(), ptcl[ip].getRSearch());
            }
#pragma omp critical
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
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        return MPI::Wtime();
#else
	return clock();
#endif
    }
    struct LessOPForVecX{
        bool operator() (const F64vec & left, const F64vec & right) const {
            //return (left.x == right.x) ? left.y < right.y : left.x < right.x;
	    return (left.x == right.x) ? ((left.y == right.y) ? ((left.z == right.z) ? left.z : left.z < right.z) : left.y < right.y ) : left.x < right.x;
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
        static One & Func( Check<F64 (T3::*)(void), &T3::getRSearch>* );
        template <typename T3>
        static One & Func( Check<F32 (T3::*)(void), &T3::getRSearch>* );
        template <typename T3>
        static One & Func( Check<F64 (T3::*)(void) const, &T3::getRSearch>* );
        template <typename T3>
        static One & Func( Check<F32 (T3::*)(void) const, &T3::getRSearch>* );
    
        template <typename T3>
        static Two &  Func(...);
    
    public:
        typedef typename HasRSearchInner< sizeof(Func<T>(NULL)) == 1 >::type type;
    };
}

#include"reallocatable_array.hpp"

#include"timer.hpp"
