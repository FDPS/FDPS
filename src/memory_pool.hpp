#pragma once
#include<iostream>
#include<cstdlib>
#include<cassert>
#include<vector>
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
#include<omp.h>
#endif

extern int N_STEP;

namespace ParticleSimulator{

    enum class MemoryAllocMode : int {Default=0, Pool=1, Stack=2};
    
    class MemoryPool{
    private:
        enum{
            ALIGN_SIZE        = 8,
            //N_SEGMENT_LIMIT   = 3000000,
            N_SEGMENT_LIMIT   = 1000,
            NAME_SIZE = 32, // This value must be same as that of corresponding ReallcatableArray variable.
        };
        typedef struct {
            void * data;
            size_t cap;
            bool used;
            void * ptr_data;
            void * ptr_id_mpool;
        } EmergencyBuffer;

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
        //char   name_[NAME_SIZE][N_SEGMENT_LIMIT];
        char   name_[N_SEGMENT_LIMIT][NAME_SIZE];
        std::vector<EmergencyBuffer> emerg_bufs_;

        static MemoryPool & getInstance(){
            static MemoryPool inst;
            return inst;
        }
        static size_t getAlignSize(const size_t _size){
            return (((_size-1)/ALIGN_SIZE)+1)*ALIGN_SIZE;
        }
        static bool isLastSegment(const int id_seg) {
            if (id_seg < N_SEGMENT_LIMIT &&
                getInstance().n_segment_ > 0 && 
                (size_t)id_seg == getInstance().n_segment_-1) {
                return true;
            } else {
                return false;
            }
        }
        static bool inParallelRegion() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
            const int n_thread = omp_get_num_threads();
#else
            const int n_thread = 1;
#endif
            if (n_thread > 1) return true;
            else return false;
        }
    public:
        static void checkEmpty(){
            if(getInstance().top_ != getInstance().bottom_
               || getInstance().n_segment_ != 0
               || getInstance().size_ != 0){
                MemoryPool::dump();
            }
            assert(getInstance().top_ == getInstance().bottom_);
            assert(getInstance().n_segment_==0);
            assert(getInstance().size_==0);
        }
        static size_t getNSegment(){
            return getInstance().n_segment_;
        }

        static size_t getSize(){
            return getInstance().size_;
        }

        static size_t getCapacity(){
            return getInstance().cap_;
        }
        
        static void initialize(const size_t _cap){
            getInstance().cap_  = _cap;
            getInstance().size_ = 0;
            getInstance().bottom_ = malloc(getInstance().cap_);
            if (getInstance().bottom_ != NULL) {
                getInstance().top_ = getInstance().bottom_;
                getInstance().n_segment_ = 0;
                for(size_t i=0; i<N_SEGMENT_LIMIT; i++){
                    getInstance().cap_per_seg_[i] = 0;
                    getInstance().used_per_seg_[i] = false;
                }
            } else {
                std::cerr << "PS_ERROR: malloc failed. (function: " << __func__
                          << ", line: " << __LINE__ << ", file: " << __FILE__ 
                          << ")" <<std::endl;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                MPI_Abort(MPI_COMM_WORLD,-1);
#endif
                std::exit(-1);
            }
            //std::cerr<<"getInstance().cap_= "<<getInstance().cap_<<std::endl;
        }
        static void reInitialize(const size_t size_add){
            const size_t cap_old = getInstance().cap_;
            const size_t cap_new = getInstance().getAlignSize( cap_old+size_add );
            void * bottom_old = getInstance().bottom_;
            void * bottom_new = malloc(cap_new);
            if (bottom_new != NULL) {
                getInstance().bottom_ = bottom_new;
                const size_t diff = (size_t)getInstance().top_ - (size_t)bottom_old;
                getInstance().top_  = (void*)((char*)bottom_new + diff);
                getInstance().cap_  = cap_new;
                memcpy(bottom_new, bottom_old, cap_old);
                size_t cap_cum = 0;
                for(size_t i=0; i<getInstance().n_segment_; i++){
                    void * p_new = (void*)((char*)bottom_new + cap_cum);
                    if (getInstance().used_per_seg_[i]) {
                        memcpy(getInstance().ptr_data_per_seg_[i], &p_new, sizeof(void*));
                    }
                    cap_cum += getInstance().cap_per_seg_[i];
                }
                if (bottom_old != NULL) free(bottom_old);
            } else {
                std::cerr << "PS_ERROR: malloc failed. (function: " << __func__
                          << ", line: " << __LINE__ << ", file: " << __FILE__ 
                          << ")" <<std::endl;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                MPI_Abort(MPI_COMM_WORLD,-1);
#endif
                std::exit(-1);
            }
        }

        static void unifyMem() {
            const int n_emerg_bufs = getInstance().emerg_bufs_.size();
            if (!inParallelRegion() && n_emerg_bufs > 0) {
               
                size_t size_add = 0;
                for (int i=0; i<n_emerg_bufs; i++) {
                    if (getInstance().emerg_bufs_[i].used)
                        size_add += getInstance().emerg_bufs_[i].cap;
                }
                if (getInstance().cap_ < getInstance().size_ + size_add) {
                    reInitialize(size_add);
                }
                for (int i=0; i<n_emerg_bufs; i++) {
                    if (!getInstance().emerg_bufs_[i].used) continue;
                    void * data = getInstance().emerg_bufs_[i].data;
                    const size_t cap = getInstance().emerg_bufs_[i].cap;
                    void * ptr_data = getInstance().emerg_bufs_[i].ptr_data;
                    void * ptr_id_mpool = getInstance().emerg_bufs_[i].ptr_id_mpool;
                    void * top_prev = getInstance().top_;
                    // Copy a single emergency buffer to the segment at the tail
                    getInstance().top_ = (void *)((char*)getInstance().top_ + cap);
                    getInstance().size_ += cap;
                    getInstance().cap_per_seg_[getInstance().n_segment_] = cap;
                    getInstance().used_per_seg_[getInstance().n_segment_] = true;
                    getInstance().ptr_data_per_seg_[getInstance().n_segment_] = ptr_data;
                    const int id_mpool = getInstance().n_segment_;
                    getInstance().n_segment_++;
                    memcpy(top_prev, data, cap);
                    memcpy(ptr_data, &top_prev, sizeof(void *));
                    memcpy(ptr_id_mpool, &id_mpool, sizeof(int));
                }
                // Release emergency buffers
                for (int i=0; i<n_emerg_bufs; i++) {
                    if (getInstance().emerg_bufs_[i].data != NULL)
                        free(getInstance().emerg_bufs_[i].data);
                }
                getInstance().emerg_bufs_.clear();
            } 
        }

        static void alloc(const size_t _size, int & _id_mpool, void * ptr_data, void *& ret, MemoryAllocMode mode=MemoryAllocMode::Pool,
                          const char nm[] = "nil"){

#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
            int rank = -1;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            //if(N_STEP == 1431 && rank == 14){
            //    std::cout<<"alloc 0"<<std::endl;
            //}
            
            //std::cerr<<"A) getInstance().cap_= "<<getInstance().cap_<<std::endl;
            if(_size <= 0) return;

            //if(N_STEP == 1431 && rank == 14){
            //    std::cout<<"alloc 1"<<std::endl;
            //}
            
            unifyMem();

            //if(N_STEP == 1431 && rank == 14){
            //    std::cout<<"alloc 2"<<std::endl;
            //}
            
            const size_t size_align = getAlignSize(_size);

            //if(N_STEP == 1431 && rank == 14){
            //    std::cout<<"alloc 3 size_align= "<<size_align<<std::endl;
            //}
            
            size_t cap_cum = 0;
            bool flag_break = false;
            //std::cerr<<"B) getInstance().cap_= "<<getInstance().cap_<<std::endl;
            if(mode == MemoryAllocMode::Pool){
                for(size_t i=0; i<getInstance().n_segment_; i++){
                    if( !getInstance().used_per_seg_[i] && getInstance().cap_per_seg_[i] >= size_align){
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
            }
            //std::cerr<<"C) getInstance().cap_= "<<getInstance().cap_<<std::endl;

            //if(N_STEP == 1431 && rank == 14){
            //    std::cout<<"alloc 4"<<std::endl;
            //}
            
            if(!flag_break){
                // In this case, we add a new segment to the tail of the memory pool
                assert(N_SEGMENT_LIMIT > getInstance().n_segment_+1);
                // Delete the last segment first if _id_mpool points to the last segment.
                if (isLastSegment(_id_mpool)) {

                    //if(N_STEP == 1431 && rank == 14){
                    //    std::cout<<"alloc 5"<<std::endl;
                    //}
                    
                    getInstance().n_segment_--;
                    getInstance().top_  = ((char*)getInstance().top_) - getInstance().cap_per_seg_[getInstance().n_segment_];
                    getInstance().size_ = getInstance().size_ - getInstance().cap_per_seg_[getInstance().n_segment_];
                    getInstance().cap_per_seg_[getInstance().n_segment_] = 0;
                    getInstance().used_per_seg_[getInstance().n_segment_] = false;
                    getInstance().ptr_data_per_seg_[getInstance().n_segment_] = NULL;
                }
                // Choose an operation mode
                bool flag_realloc = false;
                //std::cerr<<"getInstance().cap_= "<<getInstance().cap_
                //	 <<" getInstance().size_= "<<getInstance().size_
                //	 <<" size_align= "<<size_align
                //	 <<std::endl;
                if (getInstance().cap_ < getInstance().size_ + size_align) flag_realloc = true;
                bool flag_use_emerg_bufs = false;
                if (flag_realloc && inParallelRegion()) flag_use_emerg_bufs = true;
                // Add a new segment to the tail of the memory pool.
                if (!flag_use_emerg_bufs) {


                    //if(N_STEP == 1431 && rank == 14){
                    //    std::cout<<"alloc 6"<<std::endl;
                    //}
                    
                    if(flag_realloc) reInitialize(size_align);
                    void * top_prev = getInstance().top_;


                    //if(N_STEP == 1431 && rank == 14){
                    //    std::cout<<"alloc 7"<<std::endl;
                    //}
                    
                    getInstance().top_   = ((char*)getInstance().top_) + size_align;
                    getInstance().size_ += size_align;
                    getInstance().cap_per_seg_[getInstance().n_segment_] = size_align;
                    getInstance().used_per_seg_[getInstance().n_segment_] = true;
                    getInstance().ptr_data_per_seg_[getInstance().n_segment_] = ptr_data;


                    //if(N_STEP == 1431 && rank == 14){
                    //    std::cout<<"alloc 8"<<std::endl;
                    //    std::cout<<"nm= "<<nm<<std::endl;
                    //    std::cout<<"getInstance().n_segment_= "<<getInstance().n_segment_<<std::endl;
                    //    std::cout<<"getInstance().name_[getInstance().n_segment_]= "<<getInstance().name_[getInstance().n_segment_]<<std::endl;
                    //}
                    
                    strcpy(getInstance().name_[getInstance().n_segment_], nm);


                    //if(N_STEP == 1431 && rank == 14){
                    //    std::cout<<"alloc 9"<<std::endl;
                    //}
                    
                    _id_mpool = getInstance().n_segment_;

                    //if(N_STEP == 1431 && rank == 14){
                    //    std::cout<<"alloc 10"<<std::endl;
                    //}
                    
                    getInstance().n_segment_++;
                    ret = top_prev;

                    //if(N_STEP == 1431 && rank == 14){
                    //    std::cout<<"alloc 11"<<std::endl;
                    //}
                    
                } else {
                    int n_emerg_bufs = getInstance().emerg_bufs_.size();
                    // Newly create an emergency buffer
                    const int idx = n_emerg_bufs;
                    n_emerg_bufs++;
                    getInstance().emerg_bufs_.reserve(n_emerg_bufs);
                    getInstance().emerg_bufs_.resize(n_emerg_bufs);
                    void * data = malloc(size_align);
                    if (data != NULL) {
                        getInstance().emerg_bufs_[idx].data = data;
                        getInstance().emerg_bufs_[idx].cap = size_align;
                        getInstance().emerg_bufs_[idx].used = true;
                        getInstance().emerg_bufs_[idx].ptr_data = ptr_data;
                        getInstance().emerg_bufs_[idx].ptr_id_mpool = &_id_mpool;
                        _id_mpool = idx + N_SEGMENT_LIMIT;
                        ret = getInstance().emerg_bufs_[idx].data;
                    } else {
                        std::cerr << "PS_ERROR: malloc failed. (function: " << __func__
                                  << ", line: " << __LINE__ << ", file: " << __FILE__ 
                                  << ")" <<std::endl;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                        MPI_Abort(MPI_COMM_WORLD,-1);
#endif
                        std::exit(-1);
                    }
                }
            }
        }

        static void freeMem(const int id_seg){
            if (id_seg >= 0) {
                // This condition is needed to prevent out-of-range accesses                        
                // to the arrays cap_per_seg_[], used_per_seg_[], etc.,                             
                // which occurs when id_seg = -1, which in turn occurs when
                // the member function freeMem of class ReallocatableArray
                // is executed before the destructor of the class is executed.                      
                if (id_seg < N_SEGMENT_LIMIT) {
                    if (getInstance().cap_per_seg_[id_seg] <= 0) return;
                    getInstance().used_per_seg_[id_seg] = false;
                    if((size_t)id_seg == getInstance().n_segment_-1){                               
                        for(int i=id_seg; i>=0; i--){
                            if(getInstance().used_per_seg_[i] == true) break;                       
                            getInstance().size_ -= getInstance().cap_per_seg_[i];                   
                            getInstance().cap_per_seg_[i] = 0;
                            getInstance().ptr_data_per_seg_[i] = NULL;                              
                            getInstance().n_segment_--;                                             
                        }                                                                           
                    }
                    getInstance().top_ = ((char*)getInstance().bottom_) + getInstance().size_;      
                } else {
                    const int idx = id_seg - N_SEGMENT_LIMIT;                                       
                    getInstance().emerg_bufs_[idx].used = false;                                    
                }
                unifyMem(); 
            }
        }

        static void finalize() {
            void * bottom_old = getInstance().bottom_;
            if (bottom_old != NULL) {
                const size_t n_segment = getInstance().n_segment_;                                  
                for (size_t i=0; i<n_segment; i++) {
                    //if (getInstance().used_per_seg_[i]) {
                    if (getInstance().cap_per_seg_[i] > 0 || getInstance().used_per_seg_[i]) {
                        getInstance().used_per_seg_[i] = false;
                        getInstance().cap_per_seg_[i] = 0;
                        void * p_new = NULL;
                        memcpy(getInstance().ptr_data_per_seg_[i], &p_new, sizeof(void*));
                        getInstance().ptr_data_per_seg_[i] = NULL;
                        getInstance().n_segment_--;
                    }
                }
                assert(getInstance().n_segment_ == (size_t)0);
                getInstance().cap_ = 0;
                getInstance().size_ = 0;
                free(bottom_old);
                getInstance().bottom_ = NULL;
                getInstance().top_ = NULL;
                getInstance().n_segment_ = 0;
                int my_rank;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
#else
                my_rank = 0;
#endif
                if (my_rank == 0) {
                    std::cout << "MemoryPool::finalize() is completed!" << std::endl;
                }
            }
        }

        static void clearAll(){
            void * bottom_old = getInstance().bottom_;
            if (bottom_old != NULL) {
                const size_t n_segment = getInstance().n_segment_;
                for (size_t i=0; i<n_segment; i++) {
                    if (getInstance().cap_per_seg_[i] > 0 || getInstance().used_per_seg_[i]) {
                        getInstance().used_per_seg_[i] = false;
                        getInstance().cap_per_seg_[i] = 0;
                        void * p_new = NULL;
                        memcpy(getInstance().ptr_data_per_seg_[i], &p_new, sizeof(void*));
                        getInstance().ptr_data_per_seg_[i] = NULL;
                        getInstance().n_segment_--;
                    }
                }
                assert(getInstance().n_segment_ == (size_t)0);
                //getInstance().cap_ = 0;
                getInstance().size_ = 0;
                //getInstance().bottom_ = NULL;
                getInstance().top_ = NULL;
                getInstance().n_segment_ = 0;
            }
        }

        static void dump(std::ostream & fout=std::cerr){
            fout<<"bottom_= "<<getInstance().bottom_<<std::endl;
            fout<<"top_= "<<getInstance().top_<<std::endl;
            fout<<"cap_= "<<getInstance().cap_<<std::endl;
            fout<<"size_= "<<getInstance().size_<<std::endl;
            fout<<"n_segment= "<<getInstance().n_segment_<<std::endl;
            for(size_t i=0; i<getInstance().n_segment_; i++){
                fout<<"i= "<<i
                    <<" cap= "<<getInstance().cap_per_seg_[i]
                    <<" used= "<<getInstance().used_per_seg_[i]
                    //<<" name= "<<(ReallocatableArray<char>*)(getInstance().ptr_data_per_seg_[i]).name_
                    <<" name= "<<getInstance().name_[i]
                    <<std::endl;
            }
        }

        static size_t getRemainSize(){
            return getInstance().cap_ - getInstance().size_;
        }
        
        static void setName(const int id, const char * str){
            strcpy(getInstance().name_[id], str);
        }
        
    };
}
