#pragma once

#include<reallocatable_array.hpp>

namespace ParticleSimulator{
    template<class T>
    inline size_t MyGetClassSize(){
        return sizeof(T);
    }
    
    template<>
    inline size_t MyGetClassSize<KeyT>(){
#if defined(PARTICLE_SIMULATOR_USE_64BIT_KEY)
        return (size_t)8;
#elif defined(PARTICLE_SIMULATOR_USE_96BIT_KEY)
        return (size_t)12;
#else
        return (size_t)16;
#endif
    }
    
    // T is unsigned intger only
    template<class T, int NBIT=8>
    class RadixSort{
    private:
        int n_thread_;
        int n_bucket_;
        T mask_;
        int ** bucket_size_;
        int ** prefix_sum_;
    public:
        RadixSort(){
            //static const unsigned long int one = 1;
            const T one(1);
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
            n_thread_ = omp_get_max_threads();
#else
            n_thread_ = 1;
#endif
            n_bucket_ = 1<<NBIT;
            mask_ = 0;
            for(int i=0; i<NBIT; i++) mask_ |= one<<i;
            //std::cerr<<"mask= ";
            //mask_.dump(std::cerr);
            bucket_size_ = new int * [n_thread_];
            prefix_sum_ = new int * [n_thread_];
            for(int ith=0; ith<n_thread_; ith++){
                bucket_size_[ith] = new int [n_bucket_];
                prefix_sum_[ith] = new int [n_bucket_];
            }
        }
        ~RadixSort(){
            for(int ith=0; ith<n_thread_; ith++){
                delete [] bucket_size_[ith];
                delete [] prefix_sum_[ith];
            }
            delete [] prefix_sum_;
            delete [] bucket_size_;
        }

        template<class Tobj>
        void lsdSort(Tobj * val,
                     Tobj * val_buf,
                     const int first,
                     const int last){
            const int n_tot = last - first + 1;
            
            //int n_loop = (sizeof(T) * 8) / NBIT;
            //if((sizeof(T) * 8) % NBIT != 0) n_loop++;

            int n_loop = (MyGetClassSize<T>() * 8) / NBIT;
            if((MyGetClassSize<T>() * 8) % NBIT != 0) n_loop++;

            //std::cerr<<"sizeof(T)= "<<sizeof(T)
            //         <<" MyGetClassSize<T>()= "<<MyGetClassSize<T>()
            //         <<" NBIT= "<<NBIT
            //         <<" n_loop= "<<n_loop<<std::endl;
            
            for(int loop=0, shift=0; ; loop++, shift += NBIT){
                //std::cerr<<std::dec<<"loop= "<<loop<<" n_loop= "<<n_loop<<std::endl;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL		
#pragma omp parallel
#endif		
                {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)		    
                    int id_th = omp_get_thread_num();
#else
                    int id_th = 0;
#endif
                    for(int ibkt=0; ibkt<n_bucket_; ibkt++){
                        bucket_size_[id_th][ibkt] = 0;
                    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
                    for(int i=0; i<n_tot; i++){
                        S32 radix = (val[i].getKey() >> shift) & mask_;
                        //std::cerr<<std::hex<<"radix= "<<radix<<std::endl;
                        //std::cerr<<"shift= "<<shift
                        //         <<" val[i].getKey()= "
                        //         <<std::hex<<val[i].getKey()
                        //         <<" (val[i].getKey() >> shift)= "
                        //         <<(val[i].getKey() >> shift)
                        //         <<std::endl;
                        bucket_size_[id_th][radix]++;
                    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
                    {
                        prefix_sum_[0][0] = 0;
                        for(int j=1; j<n_thread_; j++){
                            prefix_sum_[j][0] = prefix_sum_[j-1][0] + bucket_size_[j-1][0];
                        }
                        for(int i=1; i<n_bucket_; i++){
                            prefix_sum_[0][i] = prefix_sum_[n_thread_-1][i-1] + bucket_size_[n_thread_-1][i-1];
                            for(int j=1; j<n_thread_; j++){
                                prefix_sum_[j][i] = prefix_sum_[j-1][i] + bucket_size_[j-1][i];
                            }
                        }
                    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL		    
#pragma omp for
#endif
                    for(int i=0; i<n_tot; i++){
                        S32 radix = (val[i].getKey() >> shift) & mask_;
                        const int offset = prefix_sum_[id_th][radix];
                        val_buf[offset] = val[i];
                        prefix_sum_[id_th][radix]++;
                    }
                } // OMP parallel scope
		if( loop+1 == n_loop) break;
                Tobj * Ptmp = val;
                val = val_buf;
                val_buf = Ptmp;
            }
            if(n_loop % 2 == 1){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL		
#pragma omp parallel for
#endif
                for(int i=0; i<n_tot; i++){
                    val[i] = val_buf[i];
                }
            }
        }
    };

    /*
    template<typename T, typename Tpred>
    void SampleSortOmp(T * first,
                       T * last,
                       const int n,
                       Tpred pred){
    */
    template<typename T, typename Tpred>
    inline void SampleSortOmp(ReallocatableArray<T> & val,
			      const int n0,
			      const int n1,
			      Tpred pred){
        assert(n1>=n0);
#ifndef PARTICLE_SIMULATOR_THREAD_PARALLEL
        std::sort(val.getPointer(n0), val.getPointer(n1), pred);
#else
        constexpr int N_PARALLEL_LIMIT = 500;
        const int n = n1-n0;
        auto round_up_div = [](const int a, const int b){ return (a+(b-1)) / b;};
        const int n_thread = omp_get_max_threads();
        if(n_thread==1 || n < N_PARALLEL_LIMIT){
            std::sort(val.getPointer(n0), val.getPointer(n1), pred);
            return;
        }

        //const int cash_size_limit_in_byte = 1u << 21; // 4MiB per thread
        const int cash_size_limit_in_byte = 1u << 18; // 512 kiB per thread
        //const int cash_size_limit_in_byte = 1u << 16; // 128 kiB per thread
        //const int cash_size_limit_in_byte = 1u << 10; // 1 kiB per thread
        const int n_bucket_per_thread = round_up_div(round_up_div(n*sizeof(T), cash_size_limit_in_byte), n_thread);

        //const int n_bucket_per_thread = 1;
        const int n_bucket = n_bucket_per_thread * n_thread;
        const int n_smp = std::max(((n > n_bucket*2) ? n/(n_bucket*2) : n), n_bucket);

        // set partition
        ReallocatableArray<T> smp(n_smp, n_smp, MemoryAllocMode::Pool);
        ReallocatableArray<T> partition(n_bucket-1, n_bucket-1, MemoryAllocMode::Pool);
        ReallocatableArray<int> id_bucket(n, n, MemoryAllocMode::Pool);
        std::vector<ReallocatableArray<int>> bucket_size_per_thread(n_thread, ReallocatableArray<int>(MemoryAllocMode::Pool));
        for(int i=0; i<n_thread; i++)
            bucket_size_per_thread[i].resizeNoInitialize(n_bucket);        
        ReallocatableArray<int> bucket_size(n_bucket, n_bucket, MemoryAllocMode::Pool);
        ReallocatableArray<int> n_disp_val_in_bucket(n_bucket+1, n_bucket+1, MemoryAllocMode::Pool);
        std::vector<ReallocatableArray<int>> n_disp_val_per_thread(n_thread, ReallocatableArray<int>(MemoryAllocMode::Pool));
        for(int i=0; i<n_thread; i++)
            n_disp_val_per_thread[i].resizeNoInitialize(n_bucket+1);
        ReallocatableArray<T> val_tmp(n, n, MemoryAllocMode::Pool);
        std::vector<ReallocatableArray<int>> n_cnt(n_thread, ReallocatableArray<int>(MemoryAllocMode::Pool));
        for(int i=0; i<n_thread; i++)
            n_cnt[i].resizeNoInitialize(n_bucket);

        //const double t0 = GetWtime();
        
        RandomSampling(val, smp, [](const T & src, T & dst){dst = src;});
        std::sort(smp.getPointer(0), smp.getPointer(n_smp), pred);
        for(int i=0; i<n_bucket-1; i++){
            const int id = (n_smp/n_bucket)*(i+1) + std::min(n_smp%n_bucket, i);
            partition[i] = smp[id];
        }

        //const double t1 = GetWtime() - t0;
        
        // calc bucked id and count bucket size
#pragma omp parallel
        {

            const auto ith = omp_get_thread_num();
            for(int i=0; i<n_bucket; i++){
                bucket_size_per_thread[ith][i] = 0;
            }
            int head, end;
            CalcAdrToSplitData(head, end, ith, n_thread, n);

            for(int i=head; i<end; i++){
                int id = -1;
                //auto * itr = std::lower_bound(partition.getPointer(0), partition.getPointer(n_bucket-1), val[i], pred);
                auto * itr = std::lower_bound(partition.getPointer(0), partition.getPointer(n_bucket-1), val[i+n0], pred);
                if(itr != partition.getPointer(n_bucket-1)){
                    id = itr - partition.getPointer(0);
                }
                else{
                    id = n_bucket-1;
                }
                id_bucket[i] = id;
                bucket_size_per_thread[ith][id]++;
            }
#pragma omp barrier
#pragma omp for
            for(int ib=0; ib<n_bucket; ib++){
                bucket_size[ib] = 0;
                for(int ith=0; ith<n_thread; ith++){
                    bucket_size[ib] += bucket_size_per_thread[ith][ib];
                }
            }
        }

        //const double t2 = GetWtime() - t0;
        
        n_disp_val_in_bucket[0] = 0;
        for(int i=0; i<n_bucket; i++){
            n_disp_val_in_bucket[i+1] = n_disp_val_in_bucket[i] + bucket_size[i];
        }
        
        for(int ith=0; ith<n_thread; ith++)
            n_disp_val_per_thread[ith][0] = 0;
        
        for(int ib=0; ib<n_bucket; ib++){
            n_disp_val_per_thread[0][ib] = n_disp_val_in_bucket[ib];
            for(int ith=0; ith<n_thread-1; ith++){
                n_disp_val_per_thread[ith+1][ib] = bucket_size_per_thread[ith][ib] + n_disp_val_per_thread[ith][ib];
            }
        }
        
        //const double t3 = GetWtime() - t0;
        
#pragma omp parallel
        {
            const auto ith = omp_get_thread_num();
            for(int i=0; i<n_bucket; i++)
                n_cnt[ith][i] = 0;
            int head, end;
            CalcAdrToSplitData(head, end, ith, n_thread, n);
            for(int i=head; i<end; i++){
                auto id = id_bucket[i];
                //val_tmp[n_disp_val_per_thread[ith][id]+n_cnt[ith][id]] = val[i];
                val_tmp[n_disp_val_per_thread[ith][id]+n_cnt[ith][id]] = val[i+n0];
                n_cnt[ith][id]++;
            }
        }

        //const double t4 = GetWtime() - t0;
        
#pragma omp parallel
        {
            const auto ith = omp_get_thread_num();
#pragma omp barrier
            for(int ib=n_bucket_per_thread*ith; ib<n_bucket_per_thread*(ith+1); ib++){
                const auto ip_head = n_disp_val_in_bucket[ib];
                const auto ip_end  = n_disp_val_in_bucket[ib+1];
                std::sort(val_tmp.getPointer(ip_head), val_tmp.getPointer(ip_end), pred);
            }
#pragma omp barrier
#pragma omp for
            for(int i=0; i<n; i++){
                //val[i] = val_tmp[i];
                val[i+n0] = val_tmp[i];
            }
        }
        //const double t5 = GetWtime() - t0;
        //std::cerr<<"t1= "<<t1<<" t2= "<<t2<<" t3= "<<t3<<" t4= "<<t4<<" t5= "<<t5<<std::endl;
#endif        
    }

    

    template<typename T, typename Tpred>
    inline void FindSeparator(int & id_big, int & id_sml, const T * big, const T * sml, const int n_big, const int n_sml, const int id_mrg, Tpred pred){
        const int id_start_big = std::min(id_mrg, n_big);
        const int id_start_sml = std::max(0, id_mrg-n_big);
        const int id_end_big   = std::max(0, id_mrg-n_sml);
        if(id_start_big==0 && id_start_sml==0){
            id_big = id_start_big;
            id_sml = id_start_sml;
            return;
        }
#if 0
        //const int id_end_sml   = std::min(id_mrg, n_sml);
        // linear search
        id_big = id_start_big;
        id_sml = id_start_sml;
        for(int i=id_start_big; i>id_end_big; i--){
            if(big[id_big-1] <= sml[id_sml]){
                break;
            }
            id_big--;
            id_sml++;
        }
#else
        // binary search
        id_big = id_start_big;
        id_sml = id_start_sml;        
        int top = -1; // top
        int bot = id_start_big-id_end_big+1; // bottom
        int loop = 0;
        int shift = 0;
        while( (bot-top) > 1){
            shift = (top+bot)/2;
            if( (id_sml+shift >= n_sml) || (id_big-shift-1<0) || pred(big[id_big-shift-1], sml[id_sml+shift]) ){
                bot = shift;
            }
            else{
                top = shift;
            }
            assert(loop<50);
            loop++;
        }
        id_big -= bot;
        id_sml += bot;
#endif
    }


    /*
      template<typename T, typename Tpred>
    void MergeSortOmpImpl(ReallocatableArray<T> & val0, const int nb0, const int ne0,
                          ReallocatableArray<T> & val1, const int nb1, const int be1,
                          ReallocatableArray<T> & res, Tpred pred){    
    */
    template<typename T, typename Tpred>
    inline void MergeSortOmpImpl(T * first0, T * last0, T * first1, T * last1, T * res, Tpred pred){
        //const auto n0 = ne0 - nb0;
        //const auto n1 = ne1 - nb1;

#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
        const auto n0 = last0 - first0;
        const auto n1 = last1 - first1;
        const auto n  = n0+n1;
        //T * big = (n0 >= n1) ? val0.getPointer(nb0) : val1.getPointer(nb1);
        //T * sml = (n0 < n1)  ? val0.getPointer(nb0) : val1.getPointer(nb1);
        T * big = (n0 >= n1) ? first0 : first1;
        T * sml = (n0 < n1)  ? first0 : first1;
        const int n_big = std::max(n0, n1);
        const int n_sml = std::min(n0, n1);
        const int n_thread = omp_get_max_threads();
        /*
        // is it fast ?
        if((n_big+n_sml)/n_thread <= 10){
            std::merge(first0, last0, first1, last1, res, pred);
            return;
        }
        */
        ReallocatableArray<int> head_r(n_thread, n_thread, MemoryAllocMode::Pool);
        ReallocatableArray<int> end_r(n_thread, n_thread, MemoryAllocMode::Pool);
        ReallocatableArray<int> id_big(n_thread-1, n_thread-1, MemoryAllocMode::Pool);
        ReallocatableArray<int> id_sml(n_thread-1, n_thread-1, MemoryAllocMode::Pool);
        for(int i=0; i<n_thread; i++){
            CalcAdrToSplitData(head_r[i], end_r[i], i, n_thread, n);
        }
#pragma omp parallel
        {
            const int ith = omp_get_thread_num();
            // FIND SEPARATOR
            if(ith < n_thread-1){
                const auto id_mrg = end_r[ith];
                FindSeparator(id_big[ith], id_sml[ith], big, sml, n_big, n_sml, id_mrg, pred);
            }
#pragma omp barrier
            int head_b, head_s, end_b, end_s;
            if(ith==0){
                head_b = head_s = 0;
            }
            else{
                head_b = id_big[ith-1];
                head_s = id_sml[ith-1];
            }
            if(ith==n_thread-1){
                end_b = n_big;
                end_s = n_sml;
            }
            else{
                end_b = id_big[ith];
                end_s = id_sml[ith];
            }
            std::merge(big+head_b, big+end_b, sml+head_s, sml+end_s, res+head_r[ith], pred);
        }
#else
	std::merge(first0, last0, first1, last1, res, pred);
#endif
    }

    template<typename T, typename Tpred>
    /*
    void MergeSortOmp(T * first,
                      T * last,
                      const int n,
                      Tpred pred){
    */
    inline void MergeSortOmp(ReallocatableArray<T> & val,
			     const int n0,
			     const int n1,
			     Tpred pred){
#ifndef PARTICLE_SIMULATOR_THREAD_PARALLEL
        std::sort(val.getPointer(n0), val.getPointer(n1), pred);
#else
        const int n_thread = omp_get_max_threads();
        const int n = n1-n0;
        //const int n_block = n_thread;
        const int n_block_per_thread = 8;
        const int n_block = n_thread*n_block_per_thread;
        ReallocatableArray<T>   tmp(n, n, MemoryAllocMode::Pool);
        ReallocatableArray<int> head(n_block, n_block, MemoryAllocMode::Pool);
        ReallocatableArray<int> end(n_block, n_block, MemoryAllocMode::Pool);
        for(int i=0; i<n_block; i++){
            CalcAdrToSplitData(head[i], end[i], i, n_block, n);
        }
        //const double t0 = GetWtime();
#pragma omp parallel
        {
            const int ith = omp_get_thread_num();
            for(int i=0; i<n_block_per_thread; i++){
                const int id_block = ith*n_block_per_thread + i;
                std::sort(val.getPointer(head[id_block]+n0), val.getPointer(end[id_block]+n0), pred);
            }
            //std::sort(val.getPointer(head[ith]+n0), val.getPointer(end[ith]+n0), pred);
        }
        //const double t1 = GetWtime() - t0;
        
        //int n_sort = n_thread;
        int n_sort = n_block;
        int loop = 0;
        int amari = 0;
        while(1){
            amari = n_sort % 2;
            n_sort /= 2;
            //T * src = val.getPointer(n0);
            //T * dst = tmp.getPointer();
            auto * src = &val;
            auto * dst = &tmp;
            int offset_src = n0;
            int offset_dst = 0;
            if(loop % 2 == 1){
                //src = tmp.getPointer();
                //dst = val.getPointer(n0);
                src = &tmp;
                dst = &val;
                offset_src = 0;
                offset_dst = n0;
            }
            for(int i=0; i<n_sort; i++){
                MergeSortOmpImpl(src->getPointer(head[i*2+0]+offset_src), src->getPointer(end[i*2+0]+offset_src),
                                 src->getPointer(head[i*2+1]+offset_src), src->getPointer(end[i*2+1]+offset_src),
                                 dst->getPointer(head[i*2+0]+offset_dst), pred);
                head[i] = head[i*2+0];
                end[i]  = end[i*2+1];
            }
            if(amari){
                head[n_sort] = head[n_sort*2];
                end[n_sort]  = end[n_sort*2];
#pragma omp parallel for
                for(int i=head[n_sort]; i<end[n_sort]; i++){
                    (*dst)[i+offset_dst] = (*src)[i+offset_src];
                }
            }
            n_sort += amari;
            if(n_sort == 1) break;
            loop++;
        }
        if(loop % 2 == 0){
#pragma omp parallel for
            for(int i=0; i<n; i++){
                val[i+n0] = tmp[i];
            }
        }
        //const double t2 = GetWtime() - t0;
        //std::cerr<<"t1= "<<t1<<" t2= "<<t2<<std::endl;
#endif
    }
}



