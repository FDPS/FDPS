namespace ParticleSimulator {
template <class Ttc, class Ttp, class Tepj, typename Tmorton_key>
inline void FindExchangeParticleDoubleWalk(const ReallocatableArray<Tepj>& epj_A,  // received particles
                                           const ReallocatableArray<Ttc>& tc_first_B, const ReallocatableArray<S32>& n_epj_src_per_proc,
                                           const ReallocatableArray<S32>& n_image_send_per_proc_irnai,  // not needed
                                           const DomainInfo& dinfo, const S32 n_leaf_limit_B, ReallocatableArray<S32>& n_epj_send_per_proc,
                                           ReallocatableArray<S32>& n_epj_send_per_image, ReallocatableArray<S32>& n_image_send_per_proc,
                                           ReallocatableArray<S32>& adr_ep_send, ReallocatableArray<F64vec>& shift_per_image,
                                           const ReallocatableArray<Tepj>& epj_B,  // assigned
                                           const F64vec& center_tree, const F64& full_len_tree, const Tmorton_key& morton_key) {
    const auto n_proc = dinfo.getCommInfo().getNumberOfProc();
    const S32 n_thread = Comm::getNumberOfThread();
    const F64ort pos_root_domain = dinfo.getPosRootDomain();
    const F64vec len_root_domain = pos_root_domain.getFullLength();
    n_image_send_per_proc.resizeNoInitialize(n_proc);
    n_epj_send_per_proc.resizeNoInitialize(n_proc);
    S32 n_disp_epj_per_proc[n_proc + 1];
    // ReallocatableArray<S32> n_disp_epj_per_proc(n_proc+1, n_proc+1, MemoryAllocMode::Pool);
    // n_disp_epj_per_proc.resizeNoInitialize(n_proc+1);
    n_disp_epj_per_proc[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_epj_per_proc[i + 1] = n_disp_epj_per_proc[i] + n_epj_src_per_proc[i];
    }
    for (S32 i = 0; i < n_proc; i++) {
        n_epj_send_per_proc[i] = 0;
    }

    // std::vector<ReallocatableArray<Tepj>> epj_sorted_tmp(n_thread, ReallocatableArray<Tepj>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<F64vec>> shift_per_image_tmp(n_thread, ReallocatableArray<F64vec>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<Ttp>> tp_tmp(n_thread, ReallocatableArray<Ttp>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<Ttc>> tc_tmp(n_thread, ReallocatableArray<Ttc>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<S32>> adr_tc_level_partition_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<S32>> adr_ptcl_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<S32>> rank_dst_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<S32>> adr_epj_src_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<S32>> n_epj_src_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));
    // std::vector<ReallocatableArray<S32>> n_epj_send_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Pool));

    // std::vector<ReallocatableArray<Tepj>> epj_sorted_tmp(n_thread, ReallocatableArray<Tepj>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<F64vec>> shift_per_image_tmp(n_thread, ReallocatableArray<F64vec>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<Ttp>> tp_tmp(n_thread, ReallocatableArray<Ttp>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<Ttc>> tc_tmp(n_thread, ReallocatableArray<Ttc>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<S32>> adr_tc_level_partition_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<S32>> adr_ptcl_send_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<S32>> rank_dst_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<S32>> adr_epj_src_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<S32>> n_epj_src_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Default));
    // std::vector<ReallocatableArray<S32>> n_epj_send_per_image_tmp(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Default));
    ReallocatableArray<Tepj> epj_sorted_tmp[n_thread];
    ReallocatableArray<F64vec> shift_per_image_tmp[n_thread];
    // ReallocatableArray<Ttp> tp_tmp[n_thread];
    ReallocatableArray<Ttc> tc_tmp[n_thread];
    ReallocatableArray<S32> adr_tc_level_partition_tmp[n_thread];
    ReallocatableArray<S32> adr_ptcl_send_tmp[n_thread];
    ReallocatableArray<S32> rank_dst_tmp[n_thread];
    ReallocatableArray<S32> adr_epj_src_per_image_tmp[n_thread];
    ReallocatableArray<S32> n_epj_src_per_image_tmp[n_thread];
    ReallocatableArray<S32> n_epj_send_per_image_tmp[n_thread];

    for (int i = 0; i < n_thread; i++) {
        epj_sorted_tmp[i].setAllocMode(MemoryAllocMode::Pool);

        // shift_per_image_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        // shift_per_image_tmp[i].setAllocMode(MemoryAllocMode::Default);

        // tp_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        // tp_tmp[i].setAllocMode(MemoryAllocMode::Default);

        // tc_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        tc_tmp[i].setAllocMode(MemoryAllocMode::Default);

        // adr_tc_level_partition_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        adr_tc_level_partition_tmp[i].setAllocMode(MemoryAllocMode::Default);

        // adr_ptcl_send_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        adr_ptcl_send_tmp[i].setAllocMode(MemoryAllocMode::Default);

        rank_dst_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        rank_dst_tmp[i].reserve(n_proc);

        // adr_epj_src_per_image_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        adr_epj_src_per_image_tmp[i].setAllocMode(MemoryAllocMode::Default);

        // n_epj_src_per_image_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        n_epj_src_per_image_tmp[i].setAllocMode(MemoryAllocMode::Default);

        // n_epj_send_per_image_tmp[i].setAllocMode(MemoryAllocMode::Pool);
        n_epj_send_per_image_tmp[i].setAllocMode(MemoryAllocMode::Default);
    }

    for (S32 i = 0; i < n_thread; i++) {
        adr_tc_level_partition_tmp[i].resizeNoInitialize(TREE_LEVEL_LIMIT + 2);
    }
    ReallocatableArray<S32> rank_src(n_proc, 0, MemoryAllocMode::Pool);
    for (S32 i = 0; i < n_proc; i++) {
        n_image_send_per_proc[i] = 0;
        if (n_epj_src_per_proc[i] > 0) rank_src.push_back(i);
    }
    // ReallocatableArray<F64> len_peri(DIMENSION, DIMENSION, MemoryAllocMode::Pool);
    bool pa[DIMENSION];
    dinfo.getPeriodicAxis(pa);
    // for(S32 i=0; i<DIMENSION; i++){
    //     if(pa[i]==false) len_peri[i] = 0.0;
    // }
    PS_OMP_PARALLEL {
        const S32 ith = Comm::getThreadNum();
        epj_sorted_tmp[ith].clearSize();
        shift_per_image_tmp[ith].clearSize();
        adr_ptcl_send_tmp[ith].clearSize();
        rank_dst_tmp[ith].clearSize();
        tc_tmp[ith].clearSize();
        // tp_tmp[ith].clearSize();
        adr_epj_src_per_image_tmp[ith].clearSize();
        n_epj_src_per_image_tmp[ith].clearSize();
        n_epj_send_per_image_tmp[ith].clearSize();
        // S32 n_ep_send_cum_old = 0;
        // bool first_loop = true;
        //std::cerr<<"rank_src.size()= "<<rank_src.size()<<std::endl;
        PS_OMP_FOR_D4
        for (S32 ib = 0; ib < rank_src.size(); ib++) {
            const S32 rank_tmp = rank_src[ib];
            rank_dst_tmp[ith].push_back(rank_tmp);
            if (n_epj_src_per_proc[rank_tmp] <= 0) continue;
            const S32 adr_ptcl_head = n_disp_epj_per_proc[rank_tmp];
            const S32 adr_ptcl_end = n_disp_epj_per_proc[rank_tmp + 1];
            S32vec id_image_old = -9999;
            S32 n_image = 0;
            for (S32 ip = adr_ptcl_head; ip < adr_ptcl_end; ip++) {
                const F64vec pos_target = epj_A[ip].getPos();
                const S32vec id_image_new = CalcIDOfImageDomain(pos_root_domain, pos_target, pa);
                if (id_image_old != id_image_new) {
                    adr_epj_src_per_image_tmp[ith].push_back(ip);
                    n_epj_src_per_image_tmp[ith].push_back(0);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
                    shift_per_image_tmp[ith].push_back(F64vec(id_image_new.x * len_root_domain.x, id_image_new.y * len_root_domain.y));
#else
                    shift_per_image_tmp[ith].push_back(
                        F64vec(id_image_new.x * len_root_domain.x, id_image_new.y * len_root_domain.y, id_image_new.z * len_root_domain.z));
#endif
                    id_image_old = id_image_new;
                    n_image++;
                }
                n_epj_src_per_image_tmp[ith].back()++;
            }
            /*
              if(Comm::getRank()==0 && rank_tmp==1){
              std::cerr<<"adr_ptcl_head= "<<adr_ptcl_head
              <<" adr_ptcl_end= "<<adr_ptcl_end
              <<std::endl;
              }
            */
            /*
              if(Comm::getRank()==0){
              std::cerr<<"rank_tmp= "<<rank_tmp
              <<" n_image= "<<n_image
              <<" n_disp_epj_per_image_tmp[ith].size()= "<<n_disp_epj_per_image_tmp[ith].size()
              <<std::endl;
              for(S32 ip=adr_ptcl_head; ip<adr_ptcl_end; ip++){
              std::cerr<<"ip, pos= "<<ip<<", "<<epj_A[ip].getPos()<<std::endl;
              }
              }
            */
            n_image_send_per_proc[rank_tmp] = n_image;
            const S32 adr_image_end = n_epj_src_per_image_tmp[ith].size();
            const S32 adr_image_head = adr_image_end - n_image;
            const S32 n_epj_send_cum_prev = adr_ptcl_send_tmp[ith].size();
            for (S32 ii = adr_image_head; ii < adr_image_end; ii++) {
                const F64vec shift = shift_per_image_tmp[ith][ii];
                /*
                  if(Comm::getRank()==0){
                  std::cerr<<"rank_tmp= "<<rank_tmp
                  <<" ii= "<<ii
                  <<" shift= "<<shift
                  <<" n_epj_src_per_image_tmp[ith][ii]= "<<n_epj_src_per_image_tmp[ith][ii]
                  <<std::endl;
                  }
                */
                const F64ort pos_domain = dinfo.getPosDomain(rank_tmp).shift(shift);
                // const F64ort pos_domain = dinfo.getPosDomain(rank_tmp).shift(-shift);
                /*
                  if(Comm::getRank() == 0 && rank_tmp==1){
                  std::cerr<<"rank_tmp= "<<rank_tmp
                  <<" shift= "<<shift
                  <<" pos_domain= "<<pos_domain
                  <<" n_disp_epj_per_image_tmp[ith][ii]= "
                  <<n_disp_epj_per_image_tmp[ith][ii]
                  <<" n_disp_epj_per_image_tmp[ith][ii+1]= "
                  <<n_disp_epj_per_image_tmp[ith][ii+1]
                  <<std::endl;
                  }
                */
                ///////////////
                // MAKE TREE A
                S32 n_cnt = 0;
                // tp_tmp[ith].resizeNoInitialize(n_epj_src_per_image_tmp[ith][ii]);
                Ttp tp_tmp[n_epj_src_per_image_tmp[ith][ii]];
                for (S32 ip = adr_epj_src_per_image_tmp[ith][ii]; ip < adr_epj_src_per_image_tmp[ith][ii] + n_epj_src_per_image_tmp[ith][ii];
                     ip++, n_cnt++) {
                    // tp_tmp[ith][n_cnt].setFromEP(epj_A[ip], ip, morton_key);
                    tp_tmp[n_cnt].setFromEP(epj_A[ip], ip, morton_key);
                    // if(Comm::getRank()==0) std::cerr<<"epj_A[ip].pos= "<<epj_A[ip].pos<<std::endl;
                }
                // std::sort(tp_tmp[ith].getPointer(), tp_tmp[ith].getPointer()+n_cnt, LessOPKEY());
                //std::sort(tp_tmp, tp_tmp+n_cnt, LessOPKEY());
                std::sort(tp_tmp, tp_tmp + n_cnt, [](const auto& lhs, const auto& rhs) { return lhs.getKey() < rhs.getKey(); });
                /*
                bool flag = false;
                PS_OMP_CRITICAL
                    {
                        if(n_cnt > epj_sorted_tmp[ith].capacity_){
                            std::cout<<"A"<<std::endl;
                            PS_PRINT_VARIABLE3(ith, n_cnt, epj_sorted_tmp[ith].capacity_);
                            flag = true;
                        }
                    }
                */

                epj_sorted_tmp[ith].resizeNoInitialize(n_cnt);

                /*
                PS_OMP_CRITICAL
                    {
                        if(flag){
                            std::cout<<"B"<<std::endl;
                            PS_PRINT_VARIABLE3(ith, n_cnt, epj_sorted_tmp[ith].capacity_);
                        }
                    }
                */

                for (S32 ip = 0; ip < n_cnt; ip++) {
                    // const S32 adr = tp_tmp[ith][ip].adr_ptcl_;
                    const S32 adr = tp_tmp[ip].adr_ptcl_;
                    epj_sorted_tmp[ith][ip] = epj_A[adr];
                    // tp_tmp[ith][ip].adr_ptcl_ = ip;
                    tp_tmp[ip].adr_ptcl_ = ip;
                    /*
                      if(Comm::getRank()==0 && rank_tmp == 1){
                      std::cout<<"ip= "<<ip<<" epj_A[adr].pos= "<<epj_A[adr].pos<<std::endl;
                      }
                    */
                }
                /*
                  if(Comm::getRank() == 0){
                  for(S32 ip=0; ip<n_cnt; ip++){
                  std::cerr<<"epj_sorted_tmp[ith][ip].pos= "<<epj_sorted_tmp[ith][ip].pos<<std::endl;
                  }
                  }
                */
                const S32 n_leaf_limit_A = 1;
                S32 lev_max_A = 0;
                // LinkCellST(tc_tmp[ith], adr_tc_level_partition_tmp[ith].getPointer(), tp_tmp[ith].getPointer(), lev_max_A, n_cnt, n_leaf_limit_A,
                // morton_key);
                LinkCellST(tc_tmp[ith], adr_tc_level_partition_tmp[ith].getPointer(), tp_tmp, lev_max_A, n_cnt, n_leaf_limit_A, morton_key);
                CalcMomentST(adr_tc_level_partition_tmp[ith].getPointer(), tc_tmp[ith].getPointer(), epj_sorted_tmp[ith].getPointer(), lev_max_A,
                             n_leaf_limit_A);
                /*
                  if(Comm::getRank() == 0 && rank_tmp == 1){
                  std::cerr<<"tc_tmp[ith][0].mom_.getVertexOut()= "<<tc_tmp[ith][0].mom_.getVertexOut()
                  <<"tc_tmp[ith][0].mom_.getVertexIn()= "<<tc_tmp[ith][0].mom_.getVertexIn()
                  <<" epj_sorted_tmp[ith][0].pos= "<<epj_sorted_tmp[ith][0].pos
                  <<std::endl;
                  }
                */
                /*
                  if(Comm::getRank() == 0){
                  S32 err = 0;
                  tc_tmp[ith].getPointer()->checkTree(epj_sorted_tmp[ith].getPointer(),
                  tc_tmp[ith].getPointer(),
                  center_tree, full_len_tree*0.5,
                  n_leaf_limit_A, 1e-4,
                  err);
                  tc_tmp[ith].getPointer()->dumpTree(epj_sorted_tmp[ith].getPointer(),
                  tc_tmp[ith].getPointer(),
                  center_tree, full_len_tree*0.5,
                  n_leaf_limit_A);
                  }
                */
                const S32 n_epj_send_per_image_prev = adr_ptcl_send_tmp[ith].size();
                MakeListDoubleWalkTop(tc_tmp[ith], tc_first_B, epj_A, pos_domain, n_leaf_limit_A, n_leaf_limit_B, adr_ptcl_send_tmp[ith]);
                n_epj_send_per_image_tmp[ith].push_back(adr_ptcl_send_tmp[ith].size() - n_epj_send_per_image_prev);

                if (Comm::getRank() == 0 && rank_tmp == 1) {
                    /*
                      std::cout<<"rank_tmp= "<<rank_tmp
                      <<" shift= "<<shift
                      <<" pos_domain= "<<pos_domain
                      <<" tc_tmp[ith].size()= "<<tc_tmp[ith].size()
                      <<" tc_first_B.size()= "<<tc_first_B.size()
                      <<" n_epj_send_per_image_tmp[ith].back()= "<<n_epj_send_per_image_tmp[ith].back()
                      <<" adr_ptcl_send_tmp[ith].size()= "<<adr_ptcl_send_tmp[ith].size()
                      <<std::endl;
                    */
                    /*
                      for(S32 iii=n_epj_send_per_image_prev; iii<adr_ptcl_send_tmp[ith].size(); iii++){
                      std::cerr<<"adr_ptcl_send_tmp[ith][iii]= "<<adr_ptcl_send_tmp[ith][iii]
                      <<" epj_B[adr_ptcl_send_tmp[ith][iii]].id= "<<epj_B[adr_ptcl_send_tmp[ith][iii]].id
                      <<" pos= "<<epj_B[adr_ptcl_send_tmp[ith][iii]].pos
                      <<std::endl;
                      }
                    */
                }

                /*
                  const S32 n_ep_per_image = adr_ptcl_send_tmp[ith].size();
                  n_ep_per_image_tmp[ith].push_back(n_ep_per_image);
                  for(S32 jp=0; jp<n_ep_per_image; jp++){
                  const S32 adr = adr_ptcl_send_tmp[ith][jp];
                  //const F64vec pos_j = epj[adr].getPos();
                  //if( pos_root_cell_.notOverlapped(pos_j-shift) ) continue;
                  }
                */
            }  // end of for ii=0 to n_image
            n_epj_send_per_proc[rank_tmp] = adr_ptcl_send_tmp[ith].size() - n_epj_send_cum_prev;
        }  // end of OMP for
    }  // end of OMP scope
    /*
    if(Comm::getRank()==0){
        for(S32 i=0; i<n_proc; i++){
            std::cerr<<"n_image_send_per_proc[i]= "<<n_image_send_per_proc[i]
                     <<" n_epj_send_per_proc[i]= "<<n_epj_send_per_proc[i]
                     <<std::endl;
        }
    }
    */
    /*
    if(Comm::getRank()==0){
        for(S32 i=0; i<n_thread; i++){
            std::cerr<<"rank_send_tmp[i].size()= "<<rank_send_tmp[i].size()<<std::endl;
            for(S32 j=0; j<rank_send_tmp[i].size(); j++){
                std::cerr<<"rank_send_tmp[i][j]= "<<rank_send_tmp[i][j]<<std::endl;
            }
        }
    }
    */

    ReallocatableArray<S32> n_disp_image_send_per_proc(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
    n_disp_image_send_per_proc.resizeNoInitialize(n_proc + 1);
    n_disp_image_send_per_proc[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
        n_disp_image_send_per_proc[i + 1] = n_disp_image_send_per_proc[i] + n_image_send_per_proc[i];
    }
    /*
    if(Comm::getRank()==0){
        for(S32 i=0; i<n_proc; i++){
            std::cerr<<"n_image_send_per_proc[i]= "<<n_image_send_per_proc[i]
                     <<" n_disp_image_send_per_proc[i]= "<<n_disp_image_send_per_proc[i]
                     <<std::endl;
        }
    }
    */

    /*
    S32 n_image_send_tot = 0;
    for(S32 i=0; i<n_thread; i++){
        n_image_send_tot += shift_per_image_tmp[i].size();
    }
    if(Comm::getRank()==0){
        std::cerr<<"n_image_send_tot= "<<n_image_send_tot<<std::endl;
    }
    */

    const S32 n_image_send_tot = n_disp_image_send_per_proc[n_proc];
    n_epj_send_per_image.resizeNoInitialize(n_image_send_tot);
    shift_per_image.resizeNoInitialize(n_image_send_tot);
    for (S32 i = 0; i < n_thread; i++) {
        S32 n_cnt_image = 0;
        // S32 n_cnt_ep = 0;
        for (S32 j = 0; j < rank_dst_tmp[i].size(); j++) {
            const S32 rank = rank_dst_tmp[i][j];
            const S32 adr_image_head = n_disp_image_send_per_proc[rank];
            const S32 adr_image_end = n_disp_image_send_per_proc[rank + 1];
            /*
            if(Comm::getRank()==0){
                std::cerr<<"rank= "<<rank
                         <<" adr_image_head= "<<adr_image_head
                         <<" adr_image_end= "<<adr_image_end
                         <<std::endl;
            }
            */
            for (S32 k = adr_image_head; k < adr_image_end; k++, n_cnt_image++) {
                n_epj_send_per_image[k] = n_epj_send_per_image_tmp[i][n_cnt_image];
                shift_per_image[k] = shift_per_image_tmp[i][n_cnt_image];
                /*
                if(Comm::getRank()==0){
                    std::cerr<<"k= "<<k
                             <<" n_epj_send_per_image[k]= "<<n_epj_send_per_image[k]
                             <<std::endl;
                }
                */
            }
        }
    }
    ReallocatableArray<S32> n_disp_epj_send_per_image(n_image_send_tot + 1, n_image_send_tot + 1, MemoryAllocMode::Pool);
    // n_disp_epj_send_per_image.resizeNoInitialize(n_image_send_tot+1);
    n_disp_epj_send_per_image[0] = 0;
    for (S32 i = 0; i < n_image_send_tot; i++) {
        n_disp_epj_send_per_image[i + 1] = n_disp_epj_send_per_image[i] + n_epj_send_per_image[i];
    }
    S32 n_epj_send_tot = 0;
    for (S32 i = 0; i < n_proc; i++) n_epj_send_tot += n_epj_send_per_proc[i];
    adr_ep_send.resizeNoInitialize(n_epj_send_tot);

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
    for (S32 i = 0; i < n_thread; i++) {
        S32 n_cnt_image = 0;
        S32 n_cnt_ep = 0;
        for (S32 j = 0; j < rank_dst_tmp[i].size(); j++) {
            const S32 rank = rank_dst_tmp[i][j];
            const S32 adr_image_head = n_disp_image_send_per_proc[rank];
            const S32 adr_image_end = n_disp_image_send_per_proc[rank + 1];
            for (S32 k = adr_image_head; k < adr_image_end; k++, n_cnt_image++) {
                n_epj_send_per_image[k] = n_epj_send_per_image_tmp[i][n_cnt_image];
                shift_per_image[k] = shift_per_image_tmp[i][n_cnt_image];
                const S32 adr_epj_head = n_disp_epj_send_per_image[k];
                const S32 adr_epj_end = n_disp_epj_send_per_image[k + 1];
                for (S32 l = adr_epj_head; l < adr_epj_end; l++, n_cnt_ep++) {
                    adr_ep_send[l] = adr_ptcl_send_tmp[i][n_cnt_ep];
                    /*
                    if(Comm::getRank()==0){
                        std::cerr<<"l= "<<l
                                 <<" adr_ep_send[l]= "<<adr_ep_send[l]
                                 <<std::endl;
                    }
                    */
                }
            }
        }
    }
    /*
    if(Comm::getRank()==0){
        for(S32 i=0; i<n_proc; i++){
            std::cout<<"rank= "<<i<<std::endl;
            const S32 adr_image_head = n_disp_image_send_per_proc[i];
            const S32 adr_image_end  = n_disp_image_send_per_proc[i+1];
            for(S32 j=adr_image_head; j<adr_image_end; j++){
                std::cout<<"image= "<<j
                         <<" shift_per_image[j]= "<<shift_per_image[j]
                         <<std::endl;
                const S32 adr_epj_head = n_disp_epj_send_per_image[j];
                const S32 adr_epj_end  = n_disp_epj_send_per_image[j+1];
                for(S32 k=adr_epj_head; k<adr_epj_end; k++){
                    std::cout<<"k= "<<k
                             <<" adr_ep_send[k]= "<<adr_ep_send[k]
                             <<std::endl;
                }
            }
        }
    }
    */
    // Comm::barrier();
    // exit(1);
}
}  // namespace ParticleSimulator
