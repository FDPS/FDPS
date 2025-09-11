namespace ParticleSimulator {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)

static inline int GetColorForCommSplit(const int my_rank_glb, const int *n_domain, const int dim, const int cid) {
    const auto my_rank_1d = GetRank1d(my_rank_glb, n_domain, cid);
    int tmp = 1;
    for (auto i = dim - 1; i > cid; i--) {
        tmp *= n_domain[i];
    }
    return (my_rank_glb - (my_rank_1d * tmp));
}

template <typename T>
static void MyAlltoall(T *send_buf, int cnt, T *recv_buf, CommInfo &comm_info, S32 n_domain[], S32 dim) {
    const auto n_proc_glb = comm_info.getNumberOfProc();
    const auto my_rank_glb = comm_info.getRank();
    const auto n_recv_tot = cnt * n_proc_glb;
    ReallocatableArray<T> send_buf_tmp(n_recv_tot, n_recv_tot, MemoryAllocMode::Stack);
    for (auto i = 0; i < n_recv_tot; i++) {
        recv_buf[i] = send_buf[i];
    }
    for (auto d = dim - 1; d >= 0; d--) {
        const auto radix = n_proc_glb / n_domain[d];
        const auto color = GetColorForCommSplit(my_rank_glb, n_domain, dim, d);
        const auto comm_1d = comm_info.split(color, my_rank_glb);
        for (auto i = 0; i < n_proc_glb; i++) {
            const int id_send = cnt * ((i % radix) * n_domain[d] + i / radix);
            const int offset = i * cnt;
            for (auto j = 0; j < cnt; j++) {
                send_buf_tmp[offset + j] = recv_buf[id_send + j];
            }
        }
        MPI_Alltoall(send_buf_tmp.getPointer(), cnt * radix, GetDataType<T>(), recv_buf, cnt * radix, GetDataType<T>(), comm_1d.getCommunicator());
    }
}

template <typename T>
static void MyAlltoallv(const ReallocatableArray<T> &send_buf, const int *n_send, ReallocatableArray<T> &recv_buf, int *n_recv, CommInfo &comm_info,
                        S32 n_domain[], S32 dim) {
    const auto n_proc_glb = comm_info.getNumberOfProc();
    const auto my_rank_glb = comm_info.getRank();
    ReallocatableArray<int> n_send_tmp(n_proc_glb, n_proc_glb, MemoryAllocMode::Default);
    ReallocatableArray<int> n_disp_recv_tmp(n_proc_glb + 1, n_proc_glb + 1, MemoryAllocMode::Default);
    n_disp_recv_tmp[0] = 0;
    for (auto i = 0; i < n_proc_glb; i++) {
        n_recv[i] = n_send[i];
        n_disp_recv_tmp[i + 1] = n_disp_recv_tmp[i] + n_recv[i];
    }
    recv_buf.resizeNoInitialize(n_disp_recv_tmp[n_proc_glb]);
    for (auto i = 0; i < n_disp_recv_tmp[n_proc_glb]; i++) {
        recv_buf[i] = send_buf[i];
    }
    // main loop
    for (auto d = dim - 1; d >= 0; d--) {
        const auto radix = n_proc_glb / n_domain[d];
        const auto color = GetColorForCommSplit(my_rank_glb, n_domain, dim, d);
        const auto comm_1d = comm_info.split(color, my_rank_glb);

        n_disp_recv_tmp[0] = 0;
        for (auto i = 0; i < n_proc_glb; i++) {
            n_disp_recv_tmp[i + 1] = n_disp_recv_tmp[i] + n_recv[i];
        }
        recv_buf.resizeNoInitialize(n_disp_recv_tmp[n_proc_glb]);

        ReallocatableArray<T> send_buf_tmp(n_disp_recv_tmp[n_proc_glb], 0, MemoryAllocMode::Stack);
        for (auto i = 0; i < n_proc_glb; i++) {
            const int id_send = ((i % radix) * n_domain[d] + i / radix);
            // if(my_rank_glb==0){
            //     std::cerr<<"i= "<<i<<" id_send= "<<id_send<<" n_disp_recv_tmp[id_send]= "<<n_disp_recv_tmp[id_send]<<std::endl;
            // }
            n_send_tmp[i] = n_recv[id_send];
            for (auto j = 0; j < n_recv[id_send]; j++) {
                send_buf_tmp.pushBackNoCheck(recv_buf[n_disp_recv_tmp[id_send] + j]);
            }
        }
        if (my_rank_glb == 0) {
            // std::cerr<<"send_buf_tmp.size()= "<<send_buf_tmp.size()<<std::endl;
            // for(auto i=0; i<send_buf_tmp.size(); i++){
            //     std::cerr<<"i= "<<i<<" send_buf_tmp[i].id= "<<send_buf_tmp[i].id<<std::endl;
            // }
        }

        MPI_Alltoall(n_send_tmp.getPointer(), radix, GetDataType<int>(), n_recv, radix, GetDataType<int>(), comm_1d.getCommunicator());

        ReallocatableArray<int> n_recv_1d(n_domain[d], n_domain[d], MemoryAllocMode::Stack);
        ReallocatableArray<int> n_send_1d(n_domain[d], n_domain[d], MemoryAllocMode::Stack);
        ReallocatableArray<int> n_disp_recv_1d(n_domain[d] + 1, n_domain[d] + 1, MemoryAllocMode::Stack);
        ReallocatableArray<int> n_disp_send_1d(n_domain[d] + 1, n_domain[d] + 1, MemoryAllocMode::Stack);
        n_disp_recv_1d[0] = n_disp_send_1d[0] = 0;
        for (auto i = 0; i < n_domain[d]; i++) {
            n_send_1d[i] = n_recv_1d[i] = 0;
            for (auto j = 0; j < radix; j++) {
                n_send_1d[i] += n_send_tmp[i * radix + j];
                n_recv_1d[i] += n_recv[i * radix + j];
            }
            n_disp_send_1d[i + 1] = n_disp_send_1d[i] + n_send_1d[i];
            n_disp_recv_1d[i + 1] = n_disp_recv_1d[i] + n_recv_1d[i];
        }
        recv_buf.resizeNoInitialize(n_disp_recv_1d[n_domain[d]]);
        // if(my_rank_glb==0){
        //	std::cerr<<"n_domain[d]= "<<n_domain[d]<<std::endl;
        //	for(auto i=0; i<n_domain[d]; i++){
        //	    std::cerr<<"i= "<<i<<" n_send_1d[i]= "<<n_send_1d[i]<<" n_recv_1d[i]= "<<n_recv_1d[i]<<std::endl;
        //	}
        // }
        MPI_Alltoallv(send_buf_tmp.getPointer(), n_send_1d.getPointer(), n_disp_send_1d.getPointer(), GetDataType<T>(), recv_buf.getPointer(),
                      n_recv_1d.getPointer(), n_disp_recv_1d.getPointer(), GetDataType<T>(), comm_1d.getCommunicator());
        // if(my_rank_glb==0){
        //	std::cerr<<"recv_buf.size()= "<<recv_buf.size()<<std::endl;
        //	for(auto i=0; i<recv_buf.size(); i++){
        //	    std::cerr<<"i= "<<i<<" recv_buf[i].id= "<<recv_buf[i].id<<std::endl;
        //	}
        // }
    }
}

template <typename T>
static inline void MyAlltoallReuse(const ReallocatableArray<T> &send_buf, int cnt, ReallocatableArray<T> *send_buf_multi_dim,
                                   ReallocatableArray<T> *recv_buf_multi_dim, const S32 n_domain[], const S32 dim, const CommInfo &comm_info) {
    const auto n_proc_glb = comm_info.getNumberOfProc();
    const auto my_rank_glb = comm_info.getRank();
    // const auto n_recv_tot = cnt * n_proc_glb;
    for (auto d = dim - 1; d >= 0; d--) {
        const auto radix = n_proc_glb / n_domain[d];
        const auto color = GetColorForCommSplit(my_rank_glb, n_domain, dim, d);
        const auto comm_1d = comm_info.split(color, my_rank_glb);
        for (auto i = 0; i < n_proc_glb; i++) {
            const auto id_send = cnt * ((i % radix) * n_domain[d] + i / radix);
            const auto offset = i * cnt;
            if (d == dim - 1) {
                for (auto j = 0; j < cnt; j++) {
                    send_buf_multi_dim[d][offset + j] = send_buf[id_send + j];
                }
            } else {
                for (auto j = 0; j < cnt; j++) {
                    send_buf_multi_dim[d][offset + j] = recv_buf_multi_dim[d + 1][id_send + j];
                }
            }
        }
        MPI_Alltoall(send_buf_multi_dim[d].getPointer(), cnt * radix, GetDataType<T>(), recv_buf_multi_dim[d].getPointer(), cnt * radix,
                     GetDataType<T>(), comm_1d.getCommunicator());
    }
}
template <typename T>
static inline void MyAlltoallReuse(const ReallocatableArray<T> &send_buf, int cnt, ReallocatableArray<T> *send_buf_multi_dim,
                                   ReallocatableArray<T> *recv_buf_multi_dim, const S32 n_domain[], const S32 dim) {
    auto ci = Comm::getCommInfo();
    MyAlltoallReuse(send_buf, cnt, send_buf_multi_dim, recv_buf_multi_dim, n_domain, dim, ci);
}

template <typename T>
static inline void MyAlltoallVReuse(const ReallocatableArray<T> &send_buf, ReallocatableArray<T> &recv_buf,
                                    const ReallocatableArray<int> *n_send_multi_dim, const ReallocatableArray<int> *n_recv_multi_dim,
                                    const S32 n_domain[], const S32 dim, const CommInfo &comm_info) {
    const auto n_proc_glb = comm_info.getNumberOfProc();
    const auto my_rank_glb = comm_info.getRank();
    ReallocatableArray<int> n_recv_tmp(n_proc_glb, n_proc_glb, MemoryAllocMode::Stack);
    ReallocatableArray<int> n_disp_recv_tmp(n_proc_glb + 1, n_proc_glb + 1, MemoryAllocMode::Stack);
    for (auto i = 0; i < n_proc_glb; i++) {
        const auto radix = n_proc_glb / n_domain[dim - 1];
        const auto id_send = ((i % radix) * n_domain[dim - 1] + i / radix);
        n_recv_tmp[id_send] = n_send_multi_dim[dim - 1][i];
    }
    recv_buf.resizeNoInitialize(send_buf.size());
    for (auto i = 0; i < send_buf.size(); i++) {
        recv_buf[i] = send_buf[i];
    }
    // main loop
    for (auto d = dim - 1; d >= 0; d--) {
        const auto radix = n_proc_glb / n_domain[d];
        const auto color = GetColorForCommSplit(my_rank_glb, n_domain, dim, d);
        const auto comm_1d = comm_info.split(color, my_rank_glb);
        n_disp_recv_tmp[0] = 0;
        for (auto i = 0; i < n_proc_glb; i++) {
            n_disp_recv_tmp[i + 1] = n_disp_recv_tmp[i] + n_recv_tmp[i];
        }
        recv_buf.resizeNoInitialize(n_disp_recv_tmp[n_proc_glb]);
        ReallocatableArray<T> send_buf_tmp(n_disp_recv_tmp[n_proc_glb], 0, MemoryAllocMode::Stack);
        for (auto i = 0; i < n_proc_glb; i++) {
            const auto id_send = ((i % radix) * n_domain[d] + i / radix);
            for (auto j = 0; j < n_recv_tmp[id_send]; j++) {
                send_buf_tmp.pushBackNoCheck(recv_buf[n_disp_recv_tmp[id_send] + j]);
            }
        }
        for (auto i = 0; i < n_proc_glb; i++) {
            n_recv_tmp[i] = n_recv_multi_dim[d][i];
        }
        ReallocatableArray<int> n_recv_1d(n_domain[d], n_domain[d], MemoryAllocMode::Stack);
        ReallocatableArray<int> n_send_1d(n_domain[d], n_domain[d], MemoryAllocMode::Stack);
        ReallocatableArray<int> n_disp_recv_1d(n_domain[d] + 1, n_domain[d] + 1, MemoryAllocMode::Stack);
        ReallocatableArray<int> n_disp_send_1d(n_domain[d] + 1, n_domain[d] + 1, MemoryAllocMode::Stack);
        n_disp_recv_1d[0] = n_disp_send_1d[0] = 0;
        for (auto i = 0; i < n_domain[d]; i++) {
            n_send_1d[i] = n_recv_1d[i] = 0;
            for (auto j = 0; j < radix; j++) {
                n_send_1d[i] += n_send_multi_dim[d][i * radix + j];
                n_recv_1d[i] += n_recv_multi_dim[d][i * radix + j];
            }
            n_disp_send_1d[i + 1] = n_disp_send_1d[i] + n_send_1d[i];
            n_disp_recv_1d[i + 1] = n_disp_recv_1d[i] + n_recv_1d[i];
        }
        recv_buf.resizeNoInitialize(n_disp_recv_1d[n_domain[d]]);
        MPI_Alltoallv(send_buf_tmp.getPointer(), n_send_1d.getPointer(), n_disp_send_1d.getPointer(), GetDataType<T>(), recv_buf.getPointer(),
                      n_recv_1d.getPointer(), n_disp_recv_1d.getPointer(), GetDataType<T>(), comm_1d.getCommunicator());
    }
}


#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL

}  // namespace ParticleSimulator