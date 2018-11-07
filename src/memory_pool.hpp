#include<iostream>
#include<cstdlib>
#include<cassert>

namespace ParticleSimulator{
    class MemoryPool{
    private:
        enum{
            ALLIGE_SIZE       = 8,
            N_SEGMENT_LIMIT   = 10000,
        };
        MemoryPool(){}
        ~MemoryPool(){}
        MemoryPool(const MemoryPool & mem);
        MemoryPool & operator = (const MemoryPool & mem);
        void * bottom_;
        void * top_;
        size_t cap_;
        size_t size_;
        size_t n_segment_;
        size_t cap_per_seg_[N_SEGMENT_LIMIT];
        bool   used_per_seg_[N_SEGMENT_LIMIT];
        void * ptr_data_per_seg_[N_SEGMENT_LIMIT];
        static MemoryPool & getInstance(){
            static MemoryPool inst;
            return inst;
        }
        static size_t getAllignSize(const size_t _size){
            return (((_size-1)/ALLIGE_SIZE)+1)*ALLIGE_SIZE;
        }
    public:

        static size_t getNSegment(){
            return getInstance().n_segment_;
        }

        static size_t getSize(){
            return getInstance().size_;
        }
        
        static void initialize(const size_t _cap){
            getInstance().cap_ = _cap;
            getInstance().size_ = 0;
            getInstance().bottom_ = malloc(getInstance().cap_);
            getInstance().top_ = getInstance().bottom_;
            getInstance().n_segment_ = 0;
            for(size_t i=0; i<N_SEGMENT_LIMIT; i++){
                getInstance().cap_per_seg_[i] = 0;
                getInstance().used_per_seg_[i] = false;
            }
        }
        static void reInitialize(const size_t size_add){
//#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
//#pragma omp critical
//#endif
            {
                const size_t cap_old = getInstance().cap_;
                const size_t cap_new = getInstance().getAllignSize( cap_old+size_add );
                void * bottom_old = getInstance().bottom_;
                void * bottom_new = malloc(cap_new);
                getInstance().bottom_ = bottom_new;
                const size_t diff = (size_t)getInstance().top_ - (size_t)bottom_old;
                getInstance().top_  = (void*)((char*)bottom_new + diff);
                getInstance().cap_  = cap_new;
                memcpy(bottom_new, bottom_old, cap_old);
                size_t cap_cum = 0;
                for(size_t i=0; i<getInstance().n_segment_; i++){
                    void * p_new = (void*)((char*)bottom_new + cap_cum);
                    memcpy(getInstance().ptr_data_per_seg_[i], &p_new, sizeof(void*));
                    cap_cum += getInstance().cap_per_seg_[i];
                }
                free(bottom_old);
            }
        }

        static void alloc(const size_t _size, int & _id_mpool, void * ptr_data, void *& ret){
//#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
//#pragma omp critical
//#endif
            {
                const size_t size_allign = getAllignSize(_size);
                size_t cap_cum = 0;
                bool flag_break = false;
                //std::cerr<<"A) _id_mpool= "<<_id_mpool
                //         <<" getInstance().n_segment_= "<<getInstance().n_segment_
                //         <<" getInstance().size_= "<<getInstance().size_
                //         <<" getInstance().top_= "<<getInstance().top_
                //         <<std::endl;
                //for(size_t i=0; i<getInstance().n_segment_; i++){
                //    std::cerr<<"i= "<<i<<" getInstance().cap_per_seg_[i]= "<<getInstance().cap_per_seg_[i]<<std::endl;
                //}
                for(size_t i=0; i<getInstance().n_segment_; i++){
                    if( !getInstance().used_per_seg_[i] && getInstance().cap_per_seg_[i] >= size_allign){
                        // insert to middle 
                        getInstance().used_per_seg_[i] = true;
                        getInstance().ptr_data_per_seg_[i] = ptr_data;
                        _id_mpool = i;
                        ret = (void*)((char*)getInstance().bottom_ + cap_cum);
                        flag_break = true;
                        break;
                    }
                    cap_cum += getInstance().cap_per_seg_[i];
                }
                if(!flag_break){
                    // add to tail
                    assert(N_SEGMENT_LIMIT > getInstance().n_segment_+1);

                    if((size_t)_id_mpool == getInstance().n_segment_-1 && getInstance().n_segment_ > 0){
                        //std::cerr<<"B) _id_mpool= "<<_id_mpool
                        //         <<" getInstance().n_segment_= "<<getInstance().n_segment_
                        //         <<" getInstance().size_= "<<getInstance().size_
                        //         <<" getInstance().top_= "<<getInstance().top_
                        //         <<" getInstance().cap_per_seg_[getInstance().n_segment_-1]= "<<getInstance().cap_per_seg_[getInstance().n_segment_-1]
                        //         <<std::endl;
                        getInstance().n_segment_--;
                        getInstance().top_  = ((char*)getInstance().top_) - getInstance().cap_per_seg_[getInstance().n_segment_];
                        getInstance().size_ = getInstance().size_ - getInstance().cap_per_seg_[getInstance().n_segment_];
                        //std::cerr<<"C) _id_mpool= "<<_id_mpool
                        //         <<" getInstance().n_segment_= "<<getInstance().n_segment_
                        //         <<" getInstance().size_= "<<getInstance().size_
                        //         <<" getInstance().top_= "<<getInstance().top_
                        //         <<" getInstance().cap_per_seg_[getInstance().n_segment_-1]= "<<getInstance().cap_per_seg_[getInstance().n_segment_-1]
                        //         <<std::endl;
                    }

                    if(getInstance().cap_ < getInstance().size_ + size_allign){
                        reInitialize(size_allign);
                    }
                    void * top_prev = getInstance().top_;
                    getInstance().top_   = ((char*)getInstance().top_) + size_allign;
                    getInstance().size_ += size_allign;
                    getInstance().cap_per_seg_[getInstance().n_segment_] = size_allign;
                    getInstance().used_per_seg_[getInstance().n_segment_] = true;
                    getInstance().ptr_data_per_seg_[getInstance().n_segment_] = ptr_data;
                    _id_mpool = getInstance().n_segment_;
                    getInstance().n_segment_++;
                    ret = top_prev;
                    //std::cerr<<"D) _id_mpool= "<<_id_mpool
                    //         <<" getInstance().n_segment_= "<<getInstance().n_segment_
                    //         <<" getInstance().size_= "<<getInstance().size_
                    //         <<" getInstance().top_= "<<getInstance().top_
                    //         <<" getInstance().cap_per_seg_[getInstance().n_segment_-1]= "<<getInstance().cap_per_seg_[getInstance().n_segment_-1]
                    //         <<std::endl;
                }
                //std::cerr<<"E) _id_mpool= "<<_id_mpool
                //         <<" getInstance().n_segment_= "<<getInstance().n_segment_
                //         <<" getInstance().size_= "<<getInstance().size_
                //         <<" getInstance().top_= "<<getInstance().top_
                //         <<" getInstance().cap_per_seg_[getInstance().n_segment_-1]= "<<getInstance().cap_per_seg_[getInstance().n_segment_-1]
                //         <<std::endl;
            }
        }

        static void freeMem(const int id_seg){
//#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
//#pragma omp critical
//#endif
            {
                getInstance().used_per_seg_[id_seg] = false;
                if((size_t)id_seg == getInstance().n_segment_-1){
                    for(int i=id_seg; i>=0; i--){
                        if(getInstance().used_per_seg_[i] == true) break;
                        getInstance().size_ -= getInstance().cap_per_seg_[i];
                        getInstance().cap_per_seg_[i] = 0;
                        getInstance().ptr_data_per_seg_[getInstance().n_segment_] = NULL;
                        getInstance().n_segment_--;
                    }
                }
                getInstance().top_ = ((char*)getInstance().bottom_) + getInstance().size_;
            }
        }


        static void dump(){
            std::cerr<<"bottom_= "<<getInstance().bottom_<<std::endl;
            std::cerr<<"top_= "<<getInstance().top_<<std::endl;
            std::cerr<<"cap_= "<<getInstance().cap_<<std::endl;
            std::cerr<<"size_= "<<getInstance().size_<<std::endl;
            std::cerr<<"n_segment= "<<getInstance().n_segment_<<std::endl;
            for(size_t i=0; i<getInstance().n_segment_; i++){
                std::cerr<<"i= "<<i
                         <<" cap= "<<getInstance().cap_per_seg_[i]
                         <<" used= "<<getInstance().used_per_seg_[i]
                         <<std::endl;
            }
        }
    };
}
