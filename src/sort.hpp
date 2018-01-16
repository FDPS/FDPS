#pragma once

namespace ParticleSimulator{

    /*
    template<class T>
    class GetKeyForSort{
    public:
        U64 operator () (const T & val) const {
            return val.getKey();
        }
    };
    */

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
            static const unsigned long int one = 1;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
            n_thread_ = omp_get_max_threads();
#else
            n_thread_ = 1;
#endif
            //std::cerr<<"n_thread_="<<n_thread_<<std::endl;
            n_bucket_ = 1<<NBIT;
            //std::cerr<<"n_bucket_="<<n_bucket_<<std::endl;
            mask_ = 0;
            for(int i=0; i<NBIT; i++) mask_ |= one<<i;
            //std::cerr<<std::hex<<"mask_="<<mask_<<std::endl;
            //std::cerr<<std::dec;
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
            //GetKeyForSort<Tobj> gk;
            //Tgk gk;
            const int n_tot = last - first + 1;
            int n_loop = (sizeof(T) * 8) / NBIT;
            if((sizeof(T) * 8) % NBIT != 0) n_loop++;

            for(int loop=0, shift=0; ; loop++, shift += NBIT){
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
                        T radix = (val[i].getKey() >> shift) & mask_;
                        //T radix = ( gk(val[i]) >> shift) & mask_;
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
                        T radix = (val[i].getKey() >> shift) & mask_;
                        //T radix = ( gk(val[i]) >> shift) & mask_;
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
}
