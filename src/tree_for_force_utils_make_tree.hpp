#pragma once
namespace ParticleSimulator {
inline bool IsLeafCell(const S32 n_cur, const S32 n_leaf_limit, const S32 lev_cur, const S32 lev_leaf_limit) {
    return (n_cur == 0) || (n_cur <= n_leaf_limit && lev_cur >= lev_leaf_limit);
}
template <class Ttc, typename Tmorton_key>
inline void LinkCell(ReallocatableArray<Ttc> &tc_array, S32 adr_tc_level_partition[], const TreeParticle tp[], S32 &lev_max, const S32 n_tot,
                     const S32 n_leaf_limit, const S32 lev_leaf_limit, const Tmorton_key &morton_key) {
                        //std::cout<<"n_leaf_limit= "<<n_leaf_limit<<" lev_leaf_limit= "<<lev_leaf_limit<<std::endl;
#if defined(PS_DEBUG_PRINT_LinkCell)
    PARTICLE_SIMULATOR_PRINT_MSG("--START LinkCell --", 0);
#endif
    tc_array.resizeNoInitialize(N_CHILDREN * 2);
    tc_array[0].n_ptcl_ = n_tot;
    tc_array[0].adr_tc_ = N_CHILDREN;
    tc_array[0].adr_ptcl_ = 0;
    adr_tc_level_partition[0] = 0;
    adr_tc_level_partition[1] = N_CHILDREN;
    lev_max = 0;
    tc_array[0].level_ = lev_max;
    for (S32 k = 1; k < N_CHILDREN; k++) {
        tc_array[k].level_ = 0;
        tc_array[k].n_ptcl_ = 0;
        tc_array[k].adr_tc_ = SetMSB((U32)(0));
        tc_array[k].adr_ptcl_ = SetMSB((U32)(0));
    }
    // if(tc_array[0].n_ptcl_ <= n_leaf_limit)return;
    if (IsLeafCell(tc_array[0].n_ptcl_, n_leaf_limit, lev_max, lev_leaf_limit)) return;
    S32 id_cell_left = 0;
    S32 id_cell_right = N_CHILDREN - 1;
    while (1) {
        // std::cerr<<"id_cell_left= "<<id_cell_left<<std::endl;
        S32 n_cell_new = 0;
        // assign particles to child cells and count # of particles in child cells
        // but loop over parent cells because they have indexes of particles
        // std::cerr<<"before loop"<<std::endl;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
#pragma omp parallel for reduction(+ : n_cell_new) schedule(static)
#endif
        for (S32 i = id_cell_left; i < id_cell_right + 1; i++) {
            const S32 n_ptcl_tmp = tc_array[i].n_ptcl_;
            /*
            if(i==id_cell_left){
                std::cout<<"n_ptcl_tmp= "<<n_ptcl_tmp<<std::endl;
                std::cout<<"n_leaf_limit= "<<n_leaf_limit<<std::endl;
                std::cout<<"lev_max= "<<lev_max<<std::endl;
                std::cout<<"lev_leaf_limit= "<<lev_leaf_limit<<std::endl;
                std::cout<<"IsLeafCell(n_ptcl_tmp, n_leaf_limit, lev_max, lev_leaf_limit)= "
                         <<IsLeafCell(n_ptcl_tmp, n_leaf_limit, lev_max, lev_leaf_limit)<<std::endl;
                         std::cout<<"tc_array[i].adr_ptcl_= "<<tc_array[i].adr_ptcl_<<std::endl;
            }
                         */
            // if(n_ptcl_tmp <= n_leaf_limit) continue;
            /*
            std::cerr<<"i= "<<i<<std::endl;
            std::cerr<<"n_ptcl_tmp= "<<n_ptcl_tmp<<std::endl;
            std::cerr<<"n_leaf_limit= "<<n_leaf_limit<<std::endl;
            std::cerr<<"lev_max= "<<lev_max<<std::endl;
            std::cerr<<"lev_leaf_limit= "<<lev_leaf_limit<<std::endl;
            std::cerr<<"IsLeafCell(n_ptcl_tmp, n_leaf_limit, lev_max, lev_leaf_limit)= "
                     <<IsLeafCell(n_ptcl_tmp, n_leaf_limit, lev_max, lev_leaf_limit)<<std::endl;
            */
            if (IsLeafCell(n_ptcl_tmp, n_leaf_limit, lev_max, lev_leaf_limit)) continue;
            const S32 adr_ptcl_tmp = tc_array[i].adr_ptcl_;
            S32 adr[N_CHILDREN];
            S32 n_cnt = 0;
            for (S32 j = 0; j < N_CHILDREN; j++) {
                const S32 adr_tc_tmp = tc_array[i].adr_tc_ + j;
                tc_array[adr_tc_tmp].adr_ptcl_ = GetPartitionID(tp, adr_ptcl_tmp, adr_ptcl_tmp + n_ptcl_tmp - 1, j, lev_max + 1, morton_key);
                tc_array[adr_tc_tmp].level_ = lev_max + 1;
                if (GetMSB(tc_array[adr_tc_tmp].adr_ptcl_)) {
                    // GetPartitionID returnes -1, if there is no appropriate tree cell.
                    tc_array[adr_tc_tmp].n_ptcl_ = 0;
                    // tc[adr_tc_tmp].adr_tc_ = -1;
                    tc_array[adr_tc_tmp].adr_tc_ = SetMSB((U32)(0));
                } else {
                    adr[n_cnt] = adr_tc_tmp;
                    n_cnt++;
                }
            }
            S32 n_ptcl_cum = 0;
            for (S32 j = 0; j < n_cnt - 1; j++) {
                tc_array[adr[j]].n_ptcl_ = tc_array[adr[j + 1]].adr_ptcl_ - tc_array[adr[j]].adr_ptcl_;
                n_ptcl_cum += tc_array[adr[j]].n_ptcl_;
            }
            assert(n_cnt > 0);
            tc_array[adr[n_cnt - 1]].n_ptcl_ = n_ptcl_tmp - n_ptcl_cum;
            n_cell_new += N_CHILDREN;
        }  // end of omp parallel scope
        //std::cout<<"n_cell_new= "<<n_cell_new<<std::endl;
        if (n_cell_new == 0) break;
        // go deeper
        id_cell_left = id_cell_right + 1;
        id_cell_right += n_cell_new;
        lev_max++;
        adr_tc_level_partition[lev_max + 1] = id_cell_right + 1;
        if (lev_max == TREE_LEVEL_LIMIT) break;
        S32 offset = id_cell_right + 1;
        auto n_thread = Comm::getNumberOfThread();
        // std::cerr<<"n_thread= "<<n_thread<<std::endl;
        if (n_thread == 1 || (id_cell_right + 1 - id_cell_left) < 10000) {
        //if (1) {
            // std::cerr<<"check a"<<std::endl;
            for (S32 i = id_cell_left; i < id_cell_right + 1; i++) {
                const S32 id_tmp = i;
                // if(!tc_array[id_tmp].isLeaf(n_leaf_limit)){
                if (!tc_array[id_tmp].isLeafPMM(n_leaf_limit, lev_leaf_limit)) {
                    tc_array[id_tmp].adr_tc_ = offset;
                    offset += N_CHILDREN;
                }
            }
        } else {
#if 1
            //S32 offset_tmp[id_cell_right - id_cell_left + 1];
            std::vector<S32> offset_tmp(id_cell_right - id_cell_left + 1);
            S32 th_offset[n_thread + 1];
            PS_OMP_PARALLEL
            {
                auto id_thread = Comm::getThreadNum();
                auto [mystart, myend] = CalcAdrToSplitData(id_thread, n_thread, id_cell_right - id_cell_left + 1);
                mystart += id_cell_left;
                myend += id_cell_left;
                auto localoffset = 0;
                /*
                for(int i = 0; i < n_thread + 1; i++) {
                    if(i == id_thread && Comm::getRank() == 0) {
                        std::cout<<"id_thread= "<<id_thread<<" mystart= "<<mystart<<" myend= "<<myend<<std::endl;
                        std::cout<<lev_max<<" lev_max"<<std::endl;
                    }
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                }
            */
                for (auto i = mystart; i < myend; i++) {
                    // if(!tc_array[i].isLeaf(n_leaf_limit)){
                    if (!tc_array[i].isLeafPMM(n_leaf_limit, lev_leaf_limit)) {
                        offset_tmp[i - id_cell_left] = localoffset;
                        localoffset += N_CHILDREN;
                    } else {
                        offset_tmp[i - id_cell_left] = -1;
                    }
                }
                th_offset[id_thread + 1] = localoffset;
PS_OMP_BARRIER
                if (id_thread == 0) {
                    th_offset[0] = offset;
                    for (auto i = 0; i < n_thread; i++) {
                        th_offset[i + 1] += th_offset[i];
                    }
                    offset = th_offset[n_thread];
                    // for(auto i = 0;i<n_thread; i++){
                    //     printf("thoffset[%d]=%d\n", i,th_offset[i]);
                    // }
                    // printf("offset=%d\n", offset);
                }
PS_OMP_BARRIER                
                for (auto i = mystart; i < myend; i++) {
                    if (offset_tmp[i - id_cell_left] >= 0) {
                        tc_array[i].adr_tc_ = th_offset[id_thread] + offset_tmp[i - id_cell_left];
                    }
                }
            }
#else
            // std::cerr<<"check b"<<std::endl;
            // std::cerr<<"id_cell_left= "<<id_cell_left<<" id_cell_right= "<<id_cell_right<<std::endl;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
            // OpenMP test
            // std::cerr<<"check b1"<<std::endl;
            // S32 offset_tmp[(id_cell_right + 128)];
            //std::vector<S32> offset_tmp(id_cell_right + 128);
            //S32 offset_tmp[id_cell_right - id_cell_left + 128];
            std::vector<S32> offset_tmp(id_cell_right - id_cell_left + 128);
            // std::cerr<<"check b2"<<std::endl;
            S32 th_offset[n_thread + 1];
            // std::cerr<<"check b3"<<std::endl;
#pragma omp parallel
            {
                // std::cerr<<"check c"<<std::endl;
                auto id_thread = Comm::getThreadNum();
                auto mywork = (id_cell_right + 1 + n_thread - id_cell_left) / n_thread;
                auto mystart = id_cell_left + id_thread * mywork;
                auto myend = mystart + mywork;
                if (myend > id_cell_right + 1) myend = id_cell_right + 1;
                // printf("idlr= %d %d, start, end=%d %d\n",
                // 	   id_cell_left, id_cell_right, mystart, myend);
                auto localoffset = 0;
                // std::cerr<<"mystart= "<<mystart<<" myend= "<<myend<<std::endl;
                for (auto i = mystart; i < myend; i++) {
                    // if(!tc_array[i].isLeaf(n_leaf_limit)){
                    if (!tc_array[i].isLeafPMM(n_leaf_limit, lev_leaf_limit)) {
                        offset_tmp[i - id_cell_left] = localoffset;
                        localoffset += N_CHILDREN;
                    } else {
                        offset_tmp[i - id_cell_left] = -1;
                    }
                }
                th_offset[id_thread + 1] = localoffset;
#pragma omp barrier
                if (id_thread == 0) {
                    th_offset[0] = offset;
                    for (auto i = 0; i < n_thread; i++) {
                        th_offset[i + 1] += th_offset[i];
                    }
                    offset = th_offset[n_thread];
                    // for(auto i = 0;i<n_thread; i++){
                    //     printf("thoffset[%d]=%d\n", i,th_offset[i]);
                    // }
                    // printf("offset=%d\n", offset);
                }
#pragma omp barrier
                for (auto i = mystart; i < myend; i++) {
                    if (offset_tmp[i - id_cell_left] >= 0) {
                        tc_array[i].adr_tc_ = th_offset[id_thread] + offset_tmp[i - id_cell_left];
                    }
                }
            } // end of omp parallel
#endif  // #if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
#endif  // #if 1
        }
        // std::cerr<<"offset= "<<offset<<std::endl;
        tc_array.resizeNoInitialize(offset);
    }

#ifdef LOOP_TREE
    auto LinkNext = [&](const U32 parent) {
        U32 parent_next = tc_array[parent].adr_tc_next_;
        U32 first_child = tc_array[parent].adr_tc_;
        tc_array[first_child + N_CHILDREN - 1].adr_tc_next_ = parent_next;  // set most left cell. point to root
        for (U32 k = first_child; k < first_child + N_CHILDREN - 1; k++) {
            tc_array[k].adr_tc_next_ = k + 1;
        }
    };

    const U32 adr_root = 0;
    const U32 adr_null = ADR_TREE_CELL_NULL;
    // set next pointer of root cell
    tc_array[adr_root].adr_tc_next_ = adr_null;

    // set next and more pointers of null cell
    tc_array[adr_null].adr_tc_next_ = adr_null;
    tc_array[adr_null].adr_tc_ = adr_null;

    LinkNext(adr_root);
    /*
    U32 parent = adr_root;
    U32 parent_next = tc_array[parent].adr_tc_next_;
    U32 first_child = tc_array[parent].adr_tc_;
    tc_array[first_child+N_CHILDREN-1].adr_tc_next_ = parent_next; // set most left cell. point to root
    for(U32 k=first_child; k<first_child+N_CHILDREN-1; k++){
        tc_array[k].adr_tc_next_ = k+1;
    }
    */
    for (U32 i = 1; i < lev_max + 1; i++) {
#pragma omp parallel for schedule(dynamic, 4)
        for (U32 j = adr_tc_level_partition[i]; j < adr_tc_level_partition[i + 1]; j++) {
            if (tc_array[j].isLeaf(n_leaf_limit)) continue;
            LinkNext(j);
            /*
            parent = j;
            parent_next = tc_array[parent].adr_tc_next_;
            first_child = tc_array[parent].adr_tc_;
            tc_array[first_child+N_CHILDREN-1].adr_tc_next_ = parent_next; // set most left cell. point to root
            for(U32 k=first_child; k<first_child+N_CHILDREN-1; k++){
                tc_array[k].adr_tc_next_ = k+1;
            }
            */
        }
    }

    // std::cerr<<"tc_array.size()= "<<tc_array.size()<<std::endl;
    // for(S32 i=0; i<tc_array.size(); i++){
    //     std::cerr<<"i= "<<i<<" more_= "<<tc_array[i].adr_tc_<<" next= "<<tc_array[i].adr_tc_next_<<std::endl;
    // }
    // Comm::barrier();
    // Abort();
#endif
#if defined(PS_DEBUG_PRINT_LinkCell)
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE2(tc_array.size(), lev_max);
#if PS_DEBUG_PRINT_LinkCell >= 2
        for (int i = 0; i < lev_max + 1; i++) {
            PS_PRINT_VARIABLE2(i, adr_tc_level_partition[i]);
        }
#endif
    }
    PARTICLE_SIMULATOR_PRINT_MSG("--END LinkCell --", 0);
#endif
}

template <class Ttc, typename Tmorton_key>
inline void LinkCellST(ReallocatableArray<Ttc> &tc_array, S32 adr_tc_level_partition[], const TreeParticle tp[], S32 &lev_max, const S32 n_tot,
                       const S32 n_leaf_limit, const Tmorton_key &morton_key) {
    tc_array.resizeNoInitialize(N_CHILDREN * 2);
    tc_array[0].n_ptcl_ = n_tot;
    tc_array[0].adr_tc_ = N_CHILDREN;
    tc_array[0].adr_ptcl_ = 0;
    adr_tc_level_partition[0] = 0;
    adr_tc_level_partition[1] = N_CHILDREN;
    lev_max = 0;
    tc_array[0].level_ = lev_max;
    for (S32 k = 1; k < N_CHILDREN; k++) {
        tc_array[k].level_ = 0;
        tc_array[k].n_ptcl_ = 0;
        tc_array[k].adr_tc_ = SetMSB((U32)(0));
        tc_array[k].adr_ptcl_ = SetMSB((U32)(0));
    }
    if (tc_array[0].n_ptcl_ <= n_leaf_limit) return;
    S32 id_cell_left = 0;
    S32 id_cell_right = N_CHILDREN - 1;
    while (1) {
        S32 n_cell_new = 0;
        // assign particles to child cells and count # of particles in child cells
        // but loop over parent cells because they have indexes of particles
        for (S32 i = id_cell_left; i < id_cell_right + 1; i++) {
            const S32 n_ptcl_tmp = tc_array[i].n_ptcl_;
            if (n_ptcl_tmp <= n_leaf_limit) continue;
            const S32 adr_ptcl_tmp = tc_array[i].adr_ptcl_;
            S32 adr[N_CHILDREN];
            S32 n_cnt = 0;
            for (S32 j = 0; j < N_CHILDREN; j++) {
                const S32 adr_tc_tmp = tc_array[i].adr_tc_ + j;
                tc_array[adr_tc_tmp].adr_ptcl_ = GetPartitionID(tp, adr_ptcl_tmp, adr_ptcl_tmp + n_ptcl_tmp - 1, j, lev_max + 1, morton_key);
                tc_array[adr_tc_tmp].level_ = lev_max + 1;
                if (GetMSB(tc_array[adr_tc_tmp].adr_ptcl_)) {
                    tc_array[adr_tc_tmp].n_ptcl_ = 0;
                    // tc[adr_tc_tmp].adr_tc_ = -1;
                    tc_array[adr_tc_tmp].adr_tc_ = SetMSB((U32)(0));
                } else {
                    adr[n_cnt] = adr_tc_tmp;
                    n_cnt++;
                }
            }
            S32 n_ptcl_cum = 0;
            for (S32 j = 0; j < n_cnt - 1; j++) {
                tc_array[adr[j]].n_ptcl_ = tc_array[adr[j + 1]].adr_ptcl_ - tc_array[adr[j]].adr_ptcl_;
                n_ptcl_cum += tc_array[adr[j]].n_ptcl_;
            }
            tc_array[adr[n_cnt - 1]].n_ptcl_ = n_ptcl_tmp - n_ptcl_cum;
            n_cell_new += N_CHILDREN;
        }  // end of omp parallel scope
        if (n_cell_new == 0) break;
        // go deeper
        id_cell_left = id_cell_right + 1;
        id_cell_right += n_cell_new;
        lev_max++;
        adr_tc_level_partition[lev_max + 1] = id_cell_right + 1;
        if (lev_max == TREE_LEVEL_LIMIT) break;
        S32 offset = id_cell_right + 1;
        for (S32 i = id_cell_left; i < id_cell_right + 1; i++) {
            const S32 id_tmp = i;
            if (!tc_array[id_tmp].isLeaf(n_leaf_limit)) {
                tc_array[id_tmp].adr_tc_ = offset;
                offset += N_CHILDREN;
            }
        }
        tc_array.resizeNoInitialize(offset);
    }
}
}  // namespace ParticleSimulator