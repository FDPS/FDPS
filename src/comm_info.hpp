namespace ParticleSimulator {

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
template <class T>
inline MPI_Datatype GetDataType() {
    static MPI_Datatype type = MPI_DATATYPE_NULL;
    if (type == MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(sizeof(T), MPI_BYTE, &type);
        MPI_Type_commit(&type);
    }
    return type;
}
template <class T>
inline MPI_Datatype GetDataType(const T &) {
    return GetDataType<T>();
}
template <>
inline MPI_Datatype GetDataType<int>() {
    return MPI_INT;
}
template <>
inline MPI_Datatype GetDataType<long>() {
    return MPI_LONG;
}
template <>
inline MPI_Datatype GetDataType<long long int>() {
    return MPI_LONG_LONG_INT;
}
template <>
inline MPI_Datatype GetDataType<unsigned int>() {
    return MPI_UNSIGNED;
}
template <>
inline MPI_Datatype GetDataType<unsigned long>() {
    return MPI_UNSIGNED_LONG;
}
template <>
inline MPI_Datatype GetDataType<unsigned long long int>() {
    return MPI_UNSIGNED_LONG_LONG;
}
template <>
inline MPI_Datatype GetDataType<float>() {
    return MPI_FLOAT;
}
template <>
inline MPI_Datatype GetDataType<double>() {
    return MPI_DOUBLE;
}

template <class Tfloat, class Tint>
inline MPI_Datatype GetDataType();
template <>
inline MPI_Datatype GetDataType<int, int>() {
    return MPI_2INT;
}
template <>
inline MPI_Datatype GetDataType<long, int>() {
    return MPI_LONG_INT;
}
template <>
inline MPI_Datatype GetDataType<float, int>() {
    return MPI_FLOAT_INT;
}
template <>
inline MPI_Datatype GetDataType<double, int>() {
    return MPI_DOUBLE_INT;
}
#endif

template <typename Tcomm>
void BarrierForMeasure(const Tcomm &comm) {
#if defined(BARRIER_FOR_MEASURE)
    comm.barrier();
#endif
}

template <int DIM>
inline void SetNumberOfDomainMultiDimension(const int n_proc, const int my_rank, int np[], int rank[]) {
    for (int i = 0; i < DIMENSION_LIMIT; i++) {
        np[i] = 1;
        rank[i] = 1;
    }
    std::vector<S32> npv;
    npv.resize(DIM);
    auto np_tmp = n_proc;
    for (int d = DIM, cid = 0; cid < DIM - 1; d--, cid++) {
        int tmp = (S32)pow(np_tmp + 0.000001, (1.0 / d) * 1.000001);
        while (np_tmp % tmp) {
            tmp--;
        }
        npv[cid] = tmp;
        np_tmp /= npv[cid];
    }
    npv[DIM - 1] = np_tmp;
    int rank_tmp = my_rank;
    std::sort(npv.begin(), npv.end(), std::greater<S32>());
    for (int i = DIM - 1; i >= 0; i--) {
        np[i] = npv[i];
        rank[i] = rank_tmp % np[i];
        rank_tmp /= np[i];
    }
}

#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
template <class T, int DIM_COMM = 2>
class CommForAllToAll {
   private:
    int rank_glb_;
    int n_proc_glb_;
    MPI_Comm comm_glb_;
    int rank_1d_[DIM_COMM];
    int n_proc_1d_[DIM_COMM];
    MPI_Comm comm_1d_[DIM_COMM];
    ReallocatableArray<T> val_send_glb_;
    int *n_recv_disp_glb_;
    int *n_send_disp_glb_;
    int *n_send_glb_;
    int *n_recv_1d_;
    int *n_send_1d_;
    int *n_recv_disp_1d_;
    int *n_send_disp_1d_;

    template <int DIM>
    void divideProc(int np[], int rank[], const int nproc, const int myrank) {
        std::vector<int> npv;
        npv.resize(DIM);
        int np_tmp = nproc;
        for (int d = DIM, cid = 0; cid < DIM - 1; d--, cid++) {
            int tmp = (int)pow(np_tmp + 0.000001, (1.0 / d) * 1.000001);
            while (np_tmp % tmp) {
                tmp--;
            }
            npv[cid] = tmp;
            np_tmp /= npv[cid];
        }
        npv[DIM - 1] = np_tmp;
        int rank_tmp = myrank;
        std::sort(npv.begin(), npv.end(), std::greater<int>());
        for (int i = DIM - 1; i >= 0; i--) {
            np[i] = npv[i];
            rank[i] = rank_tmp % np[i];
            rank_tmp /= np[i];
        }
    }

   public:
    void dumpRank() {
        for (int i = 0; i < DIM_COMM; i++) {
            int n_tmp = 0;
            MPI_Comm_size(comm_1d_[i], &n_tmp);
            std::cout << "n_proc_1d_[" << i << "]=" << n_proc_1d_[i] << " n_tmp=" << n_tmp << std::endl;
        }
        for (int i = 0; i < DIM_COMM; i++) {
            int r_tmp = 0;
            MPI_Comm_rank(comm_1d_[i], &r_tmp);
            std::cout << "rank_1d_[" << i << "]=" << rank_1d_[i] << " r_tmp=" << r_tmp << std::endl;
        }
    }

    CommForAllToAll(MPI_Comm comm = MPI_COMM_WORLD) {
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
        for (int i = 0; i < DIM_COMM; i++) {
            if (dim_max < n_proc_1d_[i]) {
                dim_max = n_proc_1d_[i];
            }
        }
        int split_color = 0;
        int factor = 1;
        for (int d = 0; d < DIM_COMM; d++) {
            split_color += rank_1d_[d] * factor;
            factor *= dim_max;
        }
        for (int d = DIM_COMM - 1; d >= 0; d--) {
            factor = rank_1d_[d];
            for (int d0 = 0; d0 < d; d0++) {
                factor *= dim_max;
            }
            MPI_Comm_split(comm_glb_, split_color - factor, rank_glb_, comm_1d_ + d);
        }
    }  // Constructor

    // alltoall
    void execute(const ReallocatableArray<T> &val_send, const int cnt, ReallocatableArray<T> &val_recv) {
        const int n_recv_tot = cnt * n_proc_glb_;
        val_recv.resizeNoInitialize(n_recv_tot);
        val_send_glb_.resizeNoInitialize(n_recv_tot);
        for (int i = 0; i < n_recv_tot; i++) {
            val_recv[i] = val_send[i];
        }
        for (int d = DIM_COMM - 1; d >= 0; d--) {
            const int radix = n_proc_glb_ / n_proc_1d_[d];
            for (int ib = 0; ib < n_proc_glb_; ib++) {
                const int id_send = cnt * ((ib % radix) * n_proc_1d_[d] + ib / radix);
                const int offset = ib * cnt;
                for (int i = 0; i < cnt; i++) {
                    val_send_glb_[offset + i] = val_recv[id_send + i];
                }
            }
            MPI_Alltoall(val_send_glb_.getPointer(), cnt * radix, GetDataType<T>(), val_recv.getPointer(), cnt * radix, GetDataType<T>(),
                         comm_1d_[d]);
        }
    }

    void execute(const T val_send[], const int cnt, T val_recv[]) {
        const int n_recv_tot = cnt * n_proc_glb_;
        val_send_glb_.resizeNoInitialize(n_recv_tot);
        for (int i = 0; i < n_recv_tot; i++) {
            val_recv[i] = val_send[i];
        }
        for (int d = DIM_COMM - 1; d >= 0; d--) {
            const int radix = n_proc_glb_ / n_proc_1d_[d];
            for (int ib = 0; ib < n_proc_glb_; ib++) {
                const int id_send = cnt * ((ib % radix) * n_proc_1d_[d] + ib / radix);
                const int offset = ib * cnt;
                for (int i = 0; i < cnt; i++) {
                    val_send_glb_[offset + i] = val_recv[id_send + i];
                }
            }
            MPI_Alltoall(val_send_glb_.getPointer(), cnt * radix, GetDataType<T>(), val_recv, cnt * radix, GetDataType<T>(), comm_1d_[d]);
        }
    }

    void executeV(const ReallocatableArray<T> &val_send, ReallocatableArray<T> &val_recv, const int n_send[], int n_recv[],
                  const int n_send_offset = 0, const int n_recv_offset = 0) {
        int cnt = 0;
        val_recv.reserveAtLeast(n_recv_offset + val_send.size() - n_send_offset);
        val_recv.resizeNoInitialize(n_recv_offset);
        for (int ib = 0; ib < n_proc_glb_; ib++) {
            n_recv[ib] = n_send[ib];
            for (int ip = 0; ip < n_recv[ib]; ip++, cnt++) {
                val_recv.pushBackNoCheck(val_send[cnt + n_send_offset]);
            }
        }
        for (int d = DIM_COMM - 1; d >= 0; d--) {
            int radix = n_proc_glb_ / n_proc_1d_[d];
            n_recv_disp_glb_[0] = 0;
            for (int i = 0; i < n_proc_glb_; i++) {
                n_recv_disp_glb_[i + 1] = n_recv_disp_glb_[i] + n_recv[i];
            }
            val_send_glb_.reserveAtLeast(val_recv.size() - n_recv_offset);
            val_send_glb_.clearSize();
            for (int ib0 = 0; ib0 < n_proc_glb_; ib0++) {
                int id_send = (ib0 % radix) * n_proc_1d_[d] + ib0 / radix;
                n_send_glb_[ib0] = n_recv[id_send];
                int offset = n_recv_disp_glb_[id_send];
                for (int ib1 = 0; ib1 < n_send_glb_[ib0]; ib1++) {
                    val_send_glb_.pushBackNoCheck(val_recv[ib1 + offset + n_recv_offset]);
                }
            }
            MPI_Alltoall(n_send_glb_, radix, MPI_INT, n_recv, radix, MPI_INT, comm_1d_[d]);
            n_send_disp_1d_[0] = n_recv_disp_1d_[0] = 0;
            for (int ib0 = 0; ib0 < n_proc_1d_[d]; ib0++) {
                n_send_1d_[ib0] = n_recv_1d_[ib0] = 0;
                int offset = ib0 * radix;
                for (int ib1 = 0; ib1 < radix; ib1++) {
                    n_send_1d_[ib0] += n_send_glb_[offset + ib1];
                    n_recv_1d_[ib0] += n_recv[offset + ib1];
                }
                n_send_disp_1d_[ib0 + 1] = n_send_disp_1d_[ib0] + n_send_1d_[ib0];
                n_recv_disp_1d_[ib0 + 1] = n_recv_disp_1d_[ib0] + n_recv_1d_[ib0];
            }
            val_recv.resizeNoInitialize(n_recv_disp_1d_[n_proc_1d_[d]] + n_recv_offset);
            MPI_Alltoallv(val_send_glb_.getPointer(), n_send_1d_, n_send_disp_1d_, GetDataType<T>(), val_recv.getPointer(n_recv_offset), n_recv_1d_,
                          n_recv_disp_1d_, GetDataType<T>(), comm_1d_[d]);
        }
    }
};  // CommForAllToAll
#endif

class CommInfo {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
    MPI_Comm comm_;
#endif
    int rank_;
    int n_proc_;
    int n_thread_;
    int rank_multi_dim_[DIMENSION];
    int n_proc_multi_dim_[DIMENSION];
    template <class T>
    T allreduceMin(const T &val) const {
        T ret = val;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        T val_tmp = val;
        MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_MIN, comm_);
#endif
        return ret;
    }
    template <class T>
    Vector2<T> allreduceMin(const Vector2<T> &val) const {
        Vector2<T> ret = val;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        Vector2<T> val_tmp = val;
        MPI_Allreduce((T *)&val_tmp.x, (T *)&ret.x, 2, GetDataType<T>(), MPI_MIN, comm_);
#endif
        return ret;
    }
    template <class T>
    Vector3<T> allreduceMin(const Vector3<T> &val) const {
        Vector3<T> ret = val;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        Vector3<T> val_tmp = val;
        MPI_Allreduce((T *)&val_tmp.x, (T *)&ret.x, 3, GetDataType<T>(), MPI_MIN, comm_);
#endif
        return ret;
    }

    template <class T>
    T allreduceMax(const T &val) const {
        T ret = val;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        T val_tmp = val;
        MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_MAX, comm_);
#endif
        return ret;
    }
    template <class T>
    Vector2<T> allreduceMax(const Vector2<T> &val) const {
        Vector2<T> ret = val;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        Vector2<T> val_tmp = val;
        MPI_Allreduce((T *)&val_tmp.x, (T *)&ret.x, 2, GetDataType<T>(), MPI_MAX, comm_);
#endif
        return ret;
    }
    template <class T>
    Vector3<T> allreduceMax(const Vector3<T> &val) const {
        Vector3<T> ret = val;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        Vector3<T> val_tmp = val;
        MPI_Allreduce((T *)&val_tmp.x, (T *)&ret.x, 3, GetDataType<T>(), MPI_MAX, comm_);
#endif
        return ret;
    }

    template <class T>
    void allreduceSum(T dst[], const T src[], const S32 n) const {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        MPI_Allreduce(src, dst, n, GetDataType<T>(), MPI_SUM, comm_);
#else
        for (S32 i = 0; i < n; i++) dst[i] = src[i];
#endif
    }

    template <class T>
    T allreduceSum(const T &val) const {
        T ret = val;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        T val_tmp = val;
        MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_SUM, comm_);
#endif
        return ret;
    }

    template <class T>
    Vector2<T> allreduceSum(const Vector2<T> &val) const {
        Vector2<T> ret = val;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        Vector2<T> val_tmp = val;
        MPI_Allreduce((T *)&val_tmp.x, (T *)&ret.x, 2, GetDataType<T>(), MPI_SUM, comm_);
#endif
        return ret;
    }
    template <class T>
    Vector3<T> allreduceSum(const Vector3<T> &val) const {
        Vector3<T> ret = val;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        Vector3<T> val_tmp = val;
        MPI_Allreduce((T *)&val_tmp.x, (T *)&ret.x, 3, GetDataType<T>(), MPI_SUM, comm_);
#endif
        return ret;
    }

    template <class Tfloat>
    void allreduceMin(const Tfloat &f_in, const int &i_in, Tfloat &f_out, int &i_out) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        struct {
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

    template <class Tfloat>
    void allreduceMax(const Tfloat &f_in, const int &i_in, Tfloat &f_out, int &i_out) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        struct {
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
    CommInfo(const MPI_Comm &c = MPI_COMM_NULL) {
        setCommunicator(c);
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        n_thread_ = omp_get_max_threads();
#else
        n_thread_ = 1;
#endif
    }
#else
    CommInfo() {
        setCommunicator();
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        n_thread_ = omp_get_max_threads();
#else
        n_thread_ = 1;
#endif
    }
#endif
    void free() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if (comm_ != MPI_COMM_NULL) {
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

    ~CommInfo() {
        // DON'T call free in destructor, because all free function in the same communicater must be called at the same time.
    }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Comm getCommunicator() const { return comm_; }
    void setCommunicator(const MPI_Comm &c = MPI_COMM_WORLD) {
        // std::cerr<<"A) setCommunicator"<<std::endl;
        // std::cerr<<"c= "<<c
        //          <<" MPI_COMM_WORLD= "<<MPI_COMM_WORLD
        //          <<" comm_= "<<comm_
        //          <<std::endl;
        if (c != MPI_COMM_NULL) {
            // std::cerr<<"b1) setCommunicator"<<std::endl;
            MPI_Comm_dup(c, &comm_);
            MPI_Comm_rank(comm_, &rank_);
            MPI_Comm_size(comm_, &n_proc_);
        } else {
            // std::cerr<<"b2) setCommunicator"<<std::endl;
            comm_ = c;
            rank_ = -1;
            n_proc_ = -1;
        }
        // std::cerr<<"C) setCommunicator"<<std::endl;
    }
#else
    void setCommunicator() {
        rank_ = 0;
        n_proc_ = 1;
    }
#endif
    int getRank() const {
        // std::cerr<<"CommInfo::getRank()"<<std::endl;
        return rank_;
    }
    void setNumberOfProcMultiDim(const int id, const int n) { n_proc_multi_dim_[id] = n; }
    void setRankMultiDim(const int id, const int r) { rank_multi_dim_[id] = r; }
    int getRankMultiDim(const int id) { return rank_multi_dim_[id]; }
    int getNumberOfProcMultiDim(const int id) { return n_proc_multi_dim_[id]; }

    int getNumberOfProc() const { return n_proc_; }

    bool isCommNull() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        return this->comm_ == MPI_COMM_NULL;
#else
        return true;
#endif
    }
    bool isNotCommNull() const { return !isCommNull(); }

    CommInfo split(int color, int key) const {
        CommInfo comm_info_new = *this;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Comm comm_new;
        MPI_Comm_split(comm_, color, key, &comm_new);
        comm_info_new.setCommunicator(comm_new);
#endif
        return comm_info_new;
    }

    const CommInfo &operator=(const CommInfo &rhs) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Comm c = rhs.getCommunicator();
        if (c != MPI_COMM_NULL) {
            MPI_Comm_dup(c, &comm_);
        } else {
            comm_ = MPI_COMM_NULL;
        }
#endif
        rank_ = rhs.rank_;
        n_proc_ = rhs.n_proc_;
        return (*this);
    }

    void barrier() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Barrier(comm_);
#endif
    }
    bool synchronizeConditionalBranchAND(const bool &local) {
        bool global = local;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        bool local_tmp = local;
        MPI_Allreduce(&local_tmp, &global, 1, MPI_C_BOOL, MPI_LAND, comm_);
#endif
        return global;
    }
    bool synchronizeConditionalBranchOR(const bool &local) {
        bool global = local;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        bool local_tmp = local;
        MPI_Allreduce(&local_tmp, &global, 1, MPI_C_BOOL, MPI_LOR, comm_);
#endif
        return global;
    }

    ///////////////////////
    // MPI BCAST WRAPPER //
    template <class T>
    inline void broadcast(T *val, const int n, const int src = 0) const {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        constexpr S64 dlimit = std::numeric_limits<int>::max() / 4;
        MPI_Datatype dtype = MPI_BYTE;
        int factor = 1;
        if (sizeof(T) % sizeof(long long int) == 0) {
            factor = sizeof(T) / sizeof(long long int);
            dtype = MPI_LONG_LONG_INT;
        } else if (sizeof(T) % sizeof(long) == 0) {
            factor = sizeof(T) / sizeof(long);
            dtype = MPI_LONG;
        } else if (sizeof(T) % sizeof(int) == 0) {
            factor = sizeof(T) / sizeof(int);
            dtype = MPI_INT;
        }
        S64 dsize = (S64)n * (S64)factor;
        if (dsize < dlimit) {
            MPI_Bcast(val, (int)dsize, dtype, src, comm_);
        } else {
            MPI_Bcast(val, n, GetDataType<T>(), src, comm_);
        }
#endif
    }

    ////////////////////////
    // MPI GATHER WRAPPER //
    template <class T>
    inline void gather(T *val_send,  // in
                       int n,        // in
                       T *val_recv,  // out
                       int dst = 0   // in
    ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Gather(val_send, n, GetDataType<T>(), val_recv, n, GetDataType<T>(), dst, comm_);
#else
        for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
    }

    template <class T>
    inline void gatherV(T *val_send,       // in
                        int n_send,        // in
                        T *val_recv,       // out
                        int *n_recv,       // in
                        int *n_recv_disp,  // in
                        int dst = 0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Gatherv(val_send, n_send, GetDataType<T>(), val_recv, n_recv, n_recv_disp, GetDataType<T>(), dst, comm_);
#else
        for (int i = 0; i < n_send; i++) val_recv[i] = val_send[i];
#endif
    }

    template <class T>
    inline void gatherVAll(ReallocatableArray<T> &val_send, int n_send, ReallocatableArray<T> &val_recv, int dst = 0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // const auto n_proc = this->getNumberOfProc();
        ReallocatableArray<int> n_recv(n_proc_, n_proc_, MemoryAllocMode::Pool);
        ReallocatableArray<int> n_disp_recv(n_proc_, n_proc_, MemoryAllocMode::Pool);
        this->gather(&n_send, 1, &n_recv[0], dst);
        n_disp_recv[0] = 0;
        for (int i = 0; i < n_proc_; i++) {
            n_disp_recv[i + 1] = n_disp_recv[i] + n_recv[i];
        }
        int n_recv_tot = 0;
        if (rank_ == dst) {
            // if(this->getRank() == dst){
            n_recv_tot = n_disp_recv[n_proc_];
        }
        val_recv.resizeNoInitialize(n_recv_tot);
        this->gatherV(val_send.getPointer(), n_send, val_recv.getPointer(), n_recv.getPointer(), n_disp_recv.getPointer(), dst);
#else
        const int n = n_send;
        for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
    }

    /////////////////////////
    // MPI SCATTER WRAPPER //
    template <class T>
    inline void scatter(T *val_send, int n, T *val_recv, int src = 0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Scatter(val_send, n, GetDataType<T>(), val_recv, n, GetDataType<T>(), src, comm_);
#else
        for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
    }
    template <class T>
    inline void scatterV(T *val_send, int *n_send, int *n_send_disp,
                         T *val_recv,  // output
                         int n_recv, int src = 0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Scatterv(val_send, n_send, n_send_disp, GetDataType<T>(), val_recv, n_recv, GetDataType<T>(), src, comm_);
#else
        // const auto n_proc = this->getNumberOfProc();
        for (int i = 0; i < n_send_disp[n_proc_]; i++) val_recv[i] = val_send[i];
#endif
    }
    template <class T>
    inline void scatterVAll(ReallocatableArray<T> &val_send, int *n_send,
                            ReallocatableArray<T> &val_recv,  // output
                            int src = 0) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        int n_recv = 0;
        this->scatter(n_send, 1, &n_recv, src);
        // const auto n_proc = this->getNumberOfProc();
        ReallocatableArray<int> n_disp_send(n_proc_ + 1, n_proc_ + 1, MemoryAllocMode::Pool);
        n_disp_send[0] = 0;
        if (rank_ == src) {
            // if(this->getRank() == src){
            for (int i = 0; i < n_proc_; i++) {
                n_disp_send[i + 1] = n_disp_send[i] + n_send[i];
            }
        } else {
            for (int i = 0; i < n_proc_; i++) {
                n_send[i] = 0;
                n_disp_send[i + 1] = n_disp_send[i] + n_send[i];
            }
        }
        val_recv.resizeNoInitialize(n_recv);
        this->scatterV(val_send.getPointer(), n_send, n_disp_send.getPointer(), val_recv.getPointer(), n_recv, src);
#else
        for (int i = 0; i < n_send[0]; i++) val_recv[i] = val_send[i];
#endif
    }

    ///////////////////////////
    // MPI ALLGATHER WRAPPER //
    template <class T>
    inline void allGather(const T *val_send,  // in
                          const int n,        // in
                          T *val_recv         // out
    ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Allgather(val_send, n, GetDataType<T>(), val_recv, n, GetDataType<T>(), comm_);
#else
        for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
    }

    template <class T>
    inline void allGatherV(const T *val_send,  // in
                           const int n_send,   // in
                           T *val_recv,        // out
                           int *n_recv,        // in
                           int *n_recv_disp    // in
    ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Allgatherv(val_send, n_send, GetDataType<T>(), val_recv, n_recv, n_recv_disp, GetDataType<T>(), comm_);
#else
        for (int i = 0; i < n_send; i++) val_recv[i] = val_send[i];
        n_recv[0] = n_recv_disp[1] = n_send;
        n_recv_disp[0] = 0;
#endif
    }

    template <class T>
    inline void allGatherVAll(ReallocatableArray<T> &val_send,  // in
                              int n_send,                       // in
                              ReallocatableArray<T> &val_recv   // out
    ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // const auto n_proc = this->getNumberOfProc();
        ReallocatableArray<int> n_recv(n_proc_, n_proc_, MemoryAllocMode::Pool);
        this->allGather(&n_send, 1, n_recv.getPointer());
        ReallocatableArray<int> n_disp_recv(n_proc_ + 1, n_proc_ + 1, MemoryAllocMode::Pool);
        n_disp_recv[0] = 0;
        for (int i = 0; i < n_proc_; i++) {
            n_disp_recv[i + 1] = n_disp_recv[i] + n_recv[i];
        }
        const auto n_recv_tot = n_disp_recv[n_proc_];
        val_recv.resizeNoInitialize(n_recv_tot);
        this->allGatherV(val_send.getPointer(), n_send, val_recv.getPointer(), n_recv.getPointer(), n_disp_recv.getPointer());
#else
        for (int i = 0; i < n_send; i++) {
            val_recv[i] = val_send[i];
        }
#endif
    }

    ///////////////////////////
    // MPI ALLTOALL WRAPPER //
    template <class T>
    inline void allToAll(T *val_send,  // in
                         const int n,  // in
                         T *val_recv   // out
    ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Alltoall(val_send, n, GetDataType<T>(), val_recv, n, GetDataType<T>(), comm_);
#else
        for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
    }

    template <class T>
    inline void allToAllV(T *val_send,       // in
                          int *n_send,       // in
                          int *n_send_disp,  // in
                          T *val_recv,       // out
                          int *n_recv,       // in
                          int *n_recv_disp   // in
    ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Alltoallv(val_send, n_send, n_send_disp, GetDataType<T>(), val_recv, n_recv, n_recv_disp, GetDataType<T>(), comm_);
#else
        for (int i = 0; i < n_send[0]; i++) val_recv[i] = val_send[i];
        n_recv[0] = n_send[0];
        n_recv_disp[0] = n_send_disp[0];
        n_recv_disp[1] = n_send_disp[1];
#endif
    }

    template <class T>
    inline void allToAllVAll(ReallocatableArray<T> &val_send,  // in
                             int *n_send,                      // in
                             ReallocatableArray<T> &val_recv   // out
    ) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // const auto n_proc = this->getNumberOfProc();
        ReallocatableArray<int> n_recv(n_proc_, n_proc_, MemoryAllocMode::Pool);
        this->allToAll(n_send, 1, n_recv.getPointer());
        ReallocatableArray<int> n_disp_recv(n_proc_ + 1, n_proc_ + 1, MemoryAllocMode::Pool);
        ReallocatableArray<int> n_disp_send(n_proc_ + 1, n_proc_ + 1, MemoryAllocMode::Pool);
        n_disp_recv[0] = n_disp_send[0] = 0;
        for (int i = 0; i < n_proc_; i++) {
            n_disp_recv[i + 1] = n_disp_recv[i] + n_recv[i];
            n_disp_send[i + 1] = n_disp_send[i] + n_send[i];
        }
        const auto n_recv_tot = n_disp_recv[n_proc_];
        val_recv.resizeNoInitialize(n_recv_tot);
        this->allToAllV(val_send.getPointer(), n_send, n_disp_send.getPointer(), val_recv.getPointer(), n_recv.getPointer(),
                        n_disp_recv.getPointer());
#else
        for (int i = 0; i < n_send[0]; i++) val_recv[i] = val_send[i];
#endif
    }

    ///////////////////////////////////////
    // MPI ALLTOALL WRAPPER (send+irecv) //
    // under construction

    // val_send[n_rank_send*n_send]
    // val_recv[n_rank_recv*n_recv]
    template <class T>
    inline void sendIrecv(T *val_send, const int *rank_send, const int n_send, const int n_rank_send, T *val_recv, const int *rank_recv,
                          const int n_recv, const int n_rank_recv) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const int tag = 0;
        ReallocatableArray<MPI_Request> req_recv(n_rank_recv, n_rank_recv, MemoryAllocMode::Pool);
        ReallocatableArray<MPI_Status> status_recv(n_rank_recv, n_rank_recv, MemoryAllocMode::Pool);
        for (int i = 0; i < n_rank_recv; i++) {
            int rank = rank_recv[i];
            MPI_Irecv(&val_recv[n_recv * i], n_recv, GetDataType<T>(), rank, tag, comm_, req_recv.getPointer(i));
        }
        for (int i = 0; i < n_rank_send; i++) {
            int rank = rank_send[i];
            MPI_Send(&val_send[n_send * i], n_send, GetDataType<T>(), rank, tag, comm_);
        }
        MPI_Waitall(n_rank_recv, req_recv.getPointer(), status_recv.getPointer());
#else
        for (int i = 0; i < n_send; i++) val_recv[i] = val_send[i];
#endif
    }

    template <class T>
    inline void sendIrecvV(T *val_send, const int *rank_send, const int *n_send, const int *n_disp_send, const int n_rank_send, T *val_recv,
                           const int *rank_recv, const int *n_recv, const int *n_disp_recv, const int n_rank_recv) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const int tag = 0;
        ReallocatableArray<MPI_Request> req_recv(n_rank_recv, n_rank_recv, MemoryAllocMode::Pool);
        ReallocatableArray<MPI_Status> status_recv(n_rank_recv, n_rank_recv, MemoryAllocMode::Pool);
        for (int i = 0; i < n_rank_recv; i++) {
            int rank = rank_recv[i];
            MPI_Irecv(&val_recv[n_disp_recv[i]], n_recv[i], GetDataType<T>(), rank, tag, comm_, &req_recv[i]);
        }
        for (int i = 0; i < n_rank_send; i++) {
            int rank = rank_send[i];
            MPI_Send(&val_send[n_disp_send[i]], n_send[i], GetDataType<T>(), rank, tag, comm_);
        }
        MPI_Waitall(n_rank_recv, req_recv.getPointer(), status_recv.getPointer());
#else
        for (int i = 0; i < n_send[0]; i++) val_recv[i] = val_send[i];
#endif
    }
    template <class T0, class T1>
    inline void sendIrecvV(T0 *val_send_0, const int *rank_send_0, const int *n_send_0, const int *n_disp_send_0, const int n_rank_send_0,
                           T1 *val_send_1, const int *rank_send_1, const int *n_send_1, const int *n_disp_send_1, const int n_rank_send_1,
                           T0 *val_recv_0, const int *rank_recv_0, const int *n_recv_0, const int *n_disp_recv_0, const int n_rank_recv_0,
                           T1 *val_recv_1, const int *rank_recv_1, const int *n_recv_1, const int *n_disp_recv_1, const int n_rank_recv_1) const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const int tag_0 = 0;
        const int tag_1 = 0;
        ReallocatableArray<MPI_Request> req_recv(n_rank_recv_0 + n_rank_recv_1, n_rank_recv_0 + n_rank_recv_1, MemoryAllocMode::Pool);
        ReallocatableArray<MPI_Status> status_recv(n_rank_recv_0 + n_rank_recv_1, n_rank_recv_0 + n_rank_recv_1, MemoryAllocMode::Pool);
        for (int i = 0; i < n_rank_recv_0; i++) {
            int rank = rank_recv_0[i];
            MPI_Irecv(&val_recv_0[n_disp_recv_0[i]], n_recv_0[i], GetDataType<T0>(), rank, tag_0, comm_, &req_recv[i]);
        }
        for (int i = 0; i < n_rank_recv_1; i++) {
            int rank = rank_recv_1[i];
            MPI_Irecv(&val_recv_1[n_disp_recv_1[i]], n_recv_1[i], GetDataType<T1>(), rank, tag_1, comm_, &req_recv[i + n_rank_recv_0]);
        }
        for (int i = 0; i < n_rank_send_0; i++) {
            int rank = rank_send_0[i];
            MPI_Send(&val_send_0[n_disp_send_0[i]], n_send_0[i], GetDataType<T0>(), rank, tag_0, comm_);
        }
        for (int i = 0; i < n_rank_send_1; i++) {
            int rank = rank_send_1[i];
            MPI_Send(&val_send_1[n_disp_send_1[i]], n_send_1[i], GetDataType<T1>(), rank, tag_1, comm_);
        }
        MPI_Waitall(n_rank_recv_0 + n_rank_recv_1, req_recv.getPointer(), status_recv.getPointer());
#else
        for (int i = 0; i < n_send_0[0]; i++) val_recv_0[i] = val_send_0[i];
        for (int i = 0; i < n_send_1[0]; i++) val_recv_1[i] = val_send_1[i];
#endif
    }

    ///////////////////////////
    // MPI SEND/RECV WRAPPER //
    template <typename T>
    inline int send(const T *val, const int n, const int dst = 0, const int tag = 0) const {
        int ret = 0;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        ret = MPI_Send(val, n, GetDataType<T>(), dst, tag, comm_);
#endif
        return ret;
    }
    template <typename T>
    inline int recv(T *val, const int n, const int src = 0, const int tag = 0) const {
        int ret = 0;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Status *stat;
        ret = MPI_Recv(val, n, GetDataType<T>(), src, tag, comm_, stat);
#endif
        return ret;
    }

    ///////////////////////////
    // MPI ALLREDUCE WRAPPER //
    template <typename T>
    inline T getMinValue(const T &val) const {
        return allreduceMin(val);
    }
    template <typename T>
    inline T getMaxValue(const T &val) const {
        return allreduceMax(val);
    }
    template <class Tfloat>
    inline void getMinValue(const Tfloat &f_in, const int &i_in, Tfloat &f_out, int &i_out) {
        static_assert(std::is_floating_point<Tfloat>::value == true, "Tfloat must be floating point type.");
        allreduceMin(f_in, i_in, f_out, i_out);
    }
    template <class Tfloat>
    inline void getMaxValue(const Tfloat &f_in, const int &i_in, Tfloat &f_out, int &i_out) {
        static_assert(std::is_floating_point<Tfloat>::value == true, "Tfloat must be floating point type.");
        allreduceMax(f_in, i_in, f_out, i_out);
    }
    template <typename T>
    inline void getSum(T dst[], const T src[], const S32 n) const {
        allreduceSum(dst, src, n);
    }    
    template <typename T>
    inline T getSum(const T &val) const {
        return allreduceSum(val);
    }

    int getNumberOfThread() const { return n_thread_; }
    int getNumThreads() const {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        return omp_get_num_threads();
#else
        return 1;
#endif
    }
    int getMaxThreads() const {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        return omp_get_max_threads();
#else
        return 1;
#endif
    }
    int getThreadNum() const {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        return omp_get_thread_num();
#else
        return 0;
#endif
    }
};

class Comm {
   private:
    Comm() {
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
        SetNumberOfDomainMultiDimension<DIM>(n_proc_, rank_, n_proc_multi_dim_, rank_multi_dim_);
    }
    ~Comm() {}
    Comm(const Comm &c) {}
    Comm &operator=(const Comm &c);
    CommInfo comm_info_;
    int rank_;
    int n_proc_;
    int n_thread_;
    int rank_multi_dim_[DIMENSION];
    int n_proc_multi_dim_[DIMENSION];
    static Comm &getInstance() {
        static Comm inst;
        return inst;
    }

   public:
    static CommInfo create(const int n, const int rank[]) { return getInstance().comm_info_.create(n, rank); }
    static CommInfo getCommInfo() {
        CommInfo ret = getInstance().comm_info_;
        return ret;
    }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    static MPI_Comm getCommunicator() { return getInstance().comm_info_.getCommunicator(); }
#endif
    static void setNumberOfProcMultiDim(const int id, const int n) { getInstance().n_proc_multi_dim_[id] = n; }
    static void setRankMultiDim(const int id, const int r) { getInstance().rank_multi_dim_[id] = r; }
    static int getRank() { return getInstance().comm_info_.getRank(); }
    static int getNumberOfProc() { return getInstance().comm_info_.getNumberOfProc(); }
    static int getRankMultiDim(const int id) { return getInstance().rank_multi_dim_[id]; }
    static int getNumberOfProcMultiDim(const int id) { return getInstance().n_proc_multi_dim_[id]; }
    static CommInfo split(int color, int key) { return getInstance().comm_info_.split(color, key); }

    static void barrier() { getInstance().comm_info_.barrier(); }
    static bool synchronizeConditionalBranchAND(const bool &local) { return getInstance().comm_info_.synchronizeConditionalBranchAND(local); }
    static bool synchronizeConditionalBranchOR(const bool &local) { return getInstance().comm_info_.synchronizeConditionalBranchOR(local); }

    ///////////////////////
    // MPI BCAST WRAPPER //
    template <class T>
    static inline void broadcast(T *val, const int n, const int src = 0) {
        getInstance().comm_info_.broadcast(val, n, src);
    }

    ////////////////////////
    // MPI GATHER WRAPPER //
    template <class T>
    static inline void gather(T *val_send,  // in
                              int n,        // in
                              T *val_recv,  // out
                              int dst = 0   // in
    ) {
        getInstance().comm_info_.gather(val_send, n, val_recv, dst);
    }

    template <class T>
    static inline void gatherV(T *val_send,       // in
                               int n_send,        // in
                               T *val_recv,       // out
                               int *n_recv,       // in
                               int *n_disp_recv,  // in
                               int dst = 0) {
        getInstance().comm_info_.gatherV(val_send, n_send, val_recv, n_recv, n_disp_recv, dst);
    }

    template <class T>
    static inline void gatherVAll(ReallocatableArray<T> &val_send, int n_send, ReallocatableArray<T> &val_recv, int dst = 0) {
        getInstance().comm_info_.gatherVAll(val_send, n_send, val_recv, dst);
    }

    /////////////////////////
    // MPI SCATTER WRAPPER //
    template <class T>
    static inline void scatter(T *val_send, int n, T *val_recv, int src = 0) {
        getInstance().comm_info_.scatter(val_send, n, val_recv, src);
    }
    template <class T>
    static inline void scatterV(T *val_send, int *n_send, int *n_disp_send,
                                T *val_recv,  // output
                                int n_recv, int src = 0) {
        getInstance().comm_info_.scatterV(val_send, n_send, n_disp_send, val_recv, n_recv, src);
    }
    template <class T>
    static inline void scatterVAll(ReallocatableArray<T> &val_send, int *n_send,
                                   ReallocatableArray<T> &val_recv,  // output
                                   int src = 0) {
        getInstance().comm_info_.scatterVAll(val_send, n_send, val_recv, src);
    }

    ///////////////////////////
    // MPI ALLGATHER WRAPPER //
    template <class T>
    static inline void allGather(const T *val_send,  // in
                                 const int n,        // in
                                 T *val_recv         // out
    ) {
        getInstance().comm_info_.allGather(val_send, n, val_recv);
    }
    template <class T>
    static inline void allGatherV(T *val_send,      // in
                                  int n_send,       // in
                                  T *val_recv,      // out
                                  int *n_recv,      // in
                                  int *n_disp_recv  // in
    ) {
        getInstance().comm_info_.allGatherV(val_send, n_send, val_recv, n_recv, n_disp_recv);
    }
    template <class T>
    static inline void allGatherVAll(ReallocatableArray<T> &val_send,  // in
                                     int n_send,                       // in
                                     ReallocatableArray<T> &val_recv   // out
    ) {
        getInstance().comm_info_.allGatherVAll(val_send, n_send, val_recv);
    }

    ///////////////////////////
    // MPI ALLTOALL WRAPPER //
    template <class T>
    static inline void allToAll(T *val_send,  // in
                                const int n,  // in
                                T *val_recv   // out
    ) {
        getInstance().comm_info_.allToAll(val_send, n, val_recv);
    }
    template <class T>
    static inline void allToAllV(T *val_send,       // in
                                 int *n_send,       // in
                                 int *n_disp_send,  // in
                                 T *val_recv,       // out
                                 int *n_recv,       // in
                                 int *n_disp_recv   // in
    ) {
        getInstance().comm_info_.allToAllV(val_send, n_send, n_disp_send, val_recv, n_recv, n_disp_recv);
    }
    template <class T>
    static inline void allToAllVAll(ReallocatableArray<T> &val_send,  // in
                                    int *n_send,                      // in
                                    ReallocatableArray<T> &val_recv   // out
    ) {
        getInstance().comm_info_.allToAllVAll(val_send, n_send, val_recv);
    }

    ///////////////////////////////////////
    // MPI ALLTOALL WRAPPER (send+irecv) //
    // under construction
    template <class T>
    static inline void sendIrecv(T *val_send, const int *rank_send, const int n_send, const int n_rank_send, T *val_recv, const int *rank_recv,
                                 const int n_recv, const int n_rank_recv) {
        getInstance().comm_info_.sendIrecv(val_send, rank_send, n_send, n_rank_send, val_recv, rank_recv, n_recv, n_rank_recv);
    }

    template <class T>
    static inline void sendIrecvV(T *val_send, const int *rank_send, const int *n_send, const int *n_disp_send, const int n_rank_send, T *val_recv,
                                  const int *rank_recv, const int *n_recv, const int *n_disp_recv, const int n_rank_recv) {
        getInstance().comm_info_.sendIrecvV(val_send, rank_send, n_send, n_disp_send, n_rank_send, val_recv, rank_recv, n_recv, n_disp_recv,
                                            n_rank_recv);
    }

    template <class T0, class T1>
    static inline void sendIrecvV(T0 *val_send_0, const int *rank_send_0, const int *n_send_0, const int *n_disp_send_0, const int n_rank_send_0,
                                  T1 *val_send_1, const int *rank_send_1, const int *n_send_1, const int *n_disp_send_1, const int n_rank_send_1,
                                  T0 *val_recv_0, const int *rank_recv_0, const int *n_recv_0, const int *n_disp_recv_0, const int n_rank_recv_0,
                                  T1 *val_recv_1, const int *rank_recv_1, const int *n_recv_1, const int *n_disp_recv_1, const int n_rank_recv_1) {
        getInstance().comm_info_.sendIrecvV(val_send_0, rank_send_0, n_send_0, n_disp_send_0, n_rank_send_0, val_send_1, rank_send_1, n_send_1,
                                            n_disp_send_1, n_rank_send_1, val_recv_0, rank_recv_0, n_recv_0, n_disp_recv_0, n_rank_recv_0, val_recv_1,
                                            rank_recv_1, n_recv_1, n_disp_recv_1, n_rank_recv_1);
    }

    ///////////////////////////
    // MPI SEND/RECV WRAPPER //
    template <typename T>
    static inline int send(const T *val, const int n, const int dst = 0, const int tag = 0) {
        return getInstance().comm_info_.send(val, n, dst, tag);
    }
    template <typename T>
    static inline int recv(T *val, const int n, const int src = 0, const int tag = 0) {
        return getInstance().comm_info_.recv(val, n, src, tag);
    }

    ///////////////////////////
    // MPI ALLREDUCE WRAPPER //
    template <typename T>
    static inline T getMinValue(const T &val) {
        return getInstance().comm_info_.getMinValue(val);
    }
    template <typename T>
    static inline T getMaxValue(const T &val) {
        return getInstance().comm_info_.getMaxValue(val);
    }
    template <class Tfloat>
    static inline void getMinValue(const Tfloat &f_in, const int &i_in, Tfloat &f_out, int &i_out) {
        static_assert(std::is_floating_point<Tfloat>::value == true, "Tfloat must be floating point type.");
        return getInstance().comm_info_.getMinValue(f_in, i_in, f_out, i_out);
    }
    template <class Tfloat>
    static inline void getMaxValue(const Tfloat &f_in, const int &i_in, Tfloat &f_out, int &i_out) {
        static_assert(std::is_floating_point<Tfloat>::value == true, "Tfloat must be floating point type.");
        return getInstance().comm_info_.getMaxValue(f_in, i_in, f_out, i_out);
    }

    template <class T>
    static inline void getSum(T dst[], const T src[], const S32 n) {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        getInstance().comm_info_.getSum(dst, src, n);
#else
        for (S32 i = 0; i < n; i++) dst[i] = src[i];
#endif
    }


    template <typename T>
    static inline T getSum(const T &val) {
        return getInstance().comm_info_.getSum(val);
    }

    static int getNumberOfThread() { return getInstance().n_thread_; }

    static int getNumThreads() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        return omp_get_num_threads();
#else
        return 1;
#endif
    }

    static int getMaxThreads() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        return omp_get_max_threads();
#else
        return 1;
#endif
    }

    static int getThreadNum() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        return omp_get_thread_num();
#else
        return 0;
#endif
    }
};  // END OF Comm

}  // namespace ParticleSimulator