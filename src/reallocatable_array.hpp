#include<typeinfo>
#include<cassert>
#include<iostream>

namespace  ParticleSimulator{
    template<class T>
        class ReallocatableArray{
    private:
        ReallocatableArray(const ReallocatableArray &);
        ReallocatableArray & operator = (const ReallocatableArray &);
        T * data_;
        int size_;
        int capacity_;
        int capacity_org_;
        int id_mpool_;
        int alloc_mode_;

        void checkCapacityOverInConstructor(){
            if(capacity_ >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
		int rank_tmp;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank_tmp);
                std::cerr<<"rank="<<rank_tmp<<std::endl;
#else
                std::cerr<<"rank=0"<<std::endl;
#endif
                std::cerr<<"typeid(T).name():"<<typeid(T).name()<<std::endl;
                std::cerr<<"capacity_="<<capacity_<<std::endl;
                std::cerr<<"LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE="<<LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE<<std::endl;
                Abort(-1);
            }            
        }
        void constructorFunction(){
            checkCapacityOverInConstructor();
            if(alloc_mode_ == 0){
                data_ = new T[capacity_];
            }
            else if(alloc_mode_ == 1) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                {
                    void * ret = NULL;
                    MemoryPool::alloc(sizeof(T)*capacity_, id_mpool_, &data_, ret);
                    data_ = (T*)ret;
                }
            }
        }
        
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        void dumpImpl() const {
            std::cerr<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            std::cerr<<"data_= "<<data_<<std::endl;
            std::cerr<<"size_="<<size_<<std::endl;
            std::cerr<<"capacity_="<<capacity_<<std::endl;
            std::cerr<<"capacity_org_="<<capacity_org_<<std::endl;
            std::cerr<<"id_mpool_="<<id_mpool_<<std::endl;
            std::cerr<<"n_expand_="<<n_expand_<<std::endl;
        }
        int n_expand_;
        void increaseNExpand(const int n_input){
            std::cerr<<"expand capacity"<<std::endl;
            std::cerr<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            std::cerr<<"n_input="<<n_input<<std::endl;
            std::cerr<<"size_="<<size_<<std::endl;
            std::cerr<<"capacity_="<<capacity_<<std::endl;
            std::cerr<<"n_expand_="<<n_expand_<<std::endl;
            n_expand_++;
        }
#else
        void dumpImpl() const {
            std::cerr<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            std::cerr<<"data_= "<<data_<<std::endl;
            std::cerr<<"size_="<<size_<<std::endl;
            std::cerr<<"capacity_="<<capacity_<<std::endl;
            std::cerr<<"capacity_org_="<<capacity_org_<<std::endl;
            std::cerr<<"id_mpool_="<<id_mpool_<<std::endl;
            std::cerr<<"data_= "<<data_<<std::endl;
        }
#endif

#ifdef __GNUC__
        void __attribute__((always_inline)) ReallocInner(const int new_cap){
#else
        void ReallocInner(const int new_cap){
#endif
            capacity_ = new_cap;
            T * data_old = data_;
            if(alloc_mode_ == 0){
                T * data_new = new T[capacity_];
                for( int i=0; i<size_; i++){
                    data_new[i] = data_old[i];
                }
                data_ = data_new;
                delete[] data_old;
            }
            else if(alloc_mode_ == 1){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                {
                    capacity_ = new_cap;
                    int id_old = id_mpool_;
                    void * ret = NULL;
                    MemoryPool::alloc(sizeof(T)*capacity_, id_mpool_, &data_, ret);
                    data_ = (T*)ret;
                    for( int i=0; i<size_; i++){
                        data_[i] = data_old[i];
                    }
                    if(id_old != id_mpool_){
                        MemoryPool::freeMem(id_old);
                    }
                }
            }
        }
	
    public:
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        ReallocatableArray() : data_(NULL), size_(0), capacity_(0), capacity_org_(0), id_mpool_(-1), alloc_mode_(0), n_expand_(0) {}
        ReallocatableArray(int cap) : data_(NULL), size_(0), capacity_(cap), capacity_org_(0), id_mpool_(-1), alloc_mode_(0), n_expand_(0) {
            constructorFunction();
        }
        ReallocatableArray(int cap, int size) : data_(NULL), size_(size), capacity_(cap), capacity_org_(0), id_mpool_(-1), alloc_mode_(0), n_expand_(0) {
            constructorFunction();
        }
        ReallocatableArray(int cap, int size, int alloc_mode) : data_(NULL), size_(size), capacity_(cap), capacity_org_(0), id_mpool_(-1), alloc_mode_(alloc_mode), n_expand_(0) {
            constructorFunction();
        }
#else
        ReallocatableArray() : data_(NULL), size_(0), capacity_(0), capacity_org_(0), id_mpool_(-1), alloc_mode_(0) {}
        ReallocatableArray(int cap) : data_(NULL), size_(0), capacity_(cap), capacity_org_(0), id_mpool_(-1), alloc_mode_(0) {
            constructorFunction();
        }
        ReallocatableArray(int cap, int size) : data_(NULL), size_(size), capacity_(cap), capacity_org_(0), id_mpool_(-1), alloc_mode_(0) {
            constructorFunction();
        }
        ReallocatableArray(int cap, int size, int alloc_mode) : data_(NULL), size_(size), capacity_(cap), capacity_org_(0), id_mpool_(-1), alloc_mode_(alloc_mode) {
            constructorFunction();
        }
#endif
        ~ReallocatableArray(){
            if(alloc_mode_ == 0){
                if(capacity_ > 0) delete [] data_;
            }
            else if(alloc_mode_ == 1) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                {
                    if(capacity_ > 0) MemoryPool::freeMem(id_mpool_);
                }
            }
            data_ = NULL;
        }

        void setAllocMode(const int alloc_mode){
            alloc_mode_ = alloc_mode;
        }

        void initialize(int cap, int size, int alloc_mode){
            capacity_ = cap;
            size_ = size;
            alloc_mode_ = alloc_mode;
            constructorFunction();
        }
        
        void reserve(const int n){
            if( n > capacity_){
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 1
                increaseNExpand(n);
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
                capacity_ = n;
                if(capacity_ >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
		    int rank_tmp;
		    MPI_Comm_rank(MPI_COMM_WORLD,&rank_tmp);
                    std::cerr<<"rank="<<rank_tmp<<std::endl;
#else
                    std::cerr<<"rank=0"<<std::endl;
#endif
                    std::cerr<<"typeid(T).name():"<<typeid(T).name()<<std::endl;
                    std::cerr<<"capacity_="<<capacity_<<std::endl;
                    std::cerr<<"LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE="<<LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE<<std::endl;
                    ParticleSimulator::Abort(-1);
                }
                ReallocInner(n);
            }
            
        }

        int size() const {return size_; }

        int capacity() const {return capacity_; }

        const T & operator [] (const int i) const {
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
            if(i > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE || i < 0 || capacity_ <= i){
                dumpImpl();
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(i <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
            assert(i >= 0);
            assert(capacity_ > i);
#endif
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 1
            if(size_<= i ){
                dumpImpl();
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(size_ > i);
#endif
            return data_[i]; 
        }

        T & operator [] (const int i){
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
            if(i > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE || i < 0 || capacity_ <= i){
                dumpImpl();
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(i <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
            assert(i >= 0);
            assert(capacity_ > i);
#endif
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 1
            if(size_ <= i){
                dumpImpl();
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(size_ > i);
#endif
            return data_[i]; 
        }

        T & front(){ return data_[0]; }

        T & back(){ return data_[size_-1]; }

        T * data(){ return data_; }

        const T * data() const { return data_; }

        void push_back(const T & val){
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
            if(size_+1 > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                dumpImpl();
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(size_+1 <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
#endif
            const int size_new = size_ + 1;
            resizeNoInitialize(size_new);
            data_[size_-1] = val;
        }
#ifdef __GNUC__
        void __attribute__((always_inline)) resizeNoInitialize (const int n){
#else
        void resizeNoInitialize (const int n){
#endif
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
            if(n > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                dumpImpl();
                std::cout<<"n="<<n<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(n <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
#endif
            const int new_cap = std::max((n+(n+3)/3)+100, capacity_org_);
            if(data_ == NULL){
                if(alloc_mode_ == 0){data_ = new T[new_cap];}
                else if(alloc_mode_ == 1){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                    {
                        void * ret = NULL;
                        MemoryPool::alloc(sizeof(T)*new_cap, id_mpool_, &data_, ret);
                        data_ = (T*)ret;
                    }
                }
                capacity_ = new_cap;
                size_ = n;
            }
            if(n > capacity_){
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 1
                increaseNExpand(n);
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
                if(new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
		    int rank_tmp;
		    MPI_Comm_rank(MPI_COMM_WORLD,&rank_tmp);
                    std::cerr<<"rank="<<rank_tmp<<std::endl;
#else
                    std::cerr<<"rank=0"<<std::endl;
#endif
                    std::cerr<<"typeid(T).name():"<<typeid(T).name()<<std::endl;
                    std::cerr<<"capacity_="<<capacity_<<std::endl;
                    std::cerr<<"LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE="<<LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE<<std::endl;
                    Abort(-1);
                }
                ReallocInner(new_cap);
            }
            size_ = n;
        }

        void dump(const std::string str=""){
            std::cerr<<str<<std::endl;
            dumpImpl();
        }

        size_t getMemSize() const { return capacity_ * sizeof(T); }

        T * getPointer(const int i=0) const { return data_+i; }
	
        void pushBackNoCheck(const T & val){
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
            assert(size_ <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
            if(size_ >= capacity_){
                dumpImpl();
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(size_ < capacity_);
#endif
            data_[size_++] = val;
        }

        void clearSize() {size_ = 0;}

        void increaseSize(const int n=1){ resizeNoInitialize(size_+n); }
        void decreaseSize(const int n=1){ resizeNoInitialize(size_-n); } // no check
        void reserveAtLeast(const int n){
            if( n >= capacity_){
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 1
                increaseNExpand(n);
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
                const int new_cap = (n+(n+3)/3) + 100;
                if(new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
		    int rank_tmp;
		    MPI_Comm_rank(MPI_COMM_WORLD,&rank_tmp);
                    std::cerr<<"rank="<<rank_tmp<<std::endl;
#else
                    std::cerr<<"rank=0"<<std::endl;
#endif
                    std::cerr<<"typeid(T).name():"<<typeid(T).name()<<std::endl;
                    std::cerr<<"capacity_="<<capacity_<<std::endl;
                    std::cerr<<"LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE="<<LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE<<std::endl;
                    Abort(-1);
                }
                ReallocInner(new_cap);
            }
        }
#ifdef __GNUC__
        void __attribute__((always_inline)) reserveEmptyAreaAtLeast (const int n_add){
#else
        void reserveEmptyAreaAtLeast (const int n_add){
#endif
            const int n = n_add + size_;
            if( n >= capacity_){
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 1
                increaseNExpand(n);
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
                const int new_cap = (n+(n+3)/3) + 100;
                if(new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
		    int rank_tmp;
		    MPI_Comm_rank(MPI_COMM_WORLD,&rank_tmp);
                    std::cerr<<"rank="<<rank_tmp<<std::endl;
#else
                    std::cerr<<"rank=0"<<std::endl;
#endif
                    std::cerr<<"typeid(T).name():"<<typeid(T).name()<<std::endl;
                    std::cerr<<"capacity_="<<capacity_<<std::endl;
                    std::cerr<<"LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE="<<LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE<<std::endl;
                    Abort(-1);
                }
                ReallocInner(new_cap);
            }
        }

        void freeMem(const int free_alloc_mode = -1){
            if( alloc_mode_ == 0 && (free_alloc_mode == -1  || free_alloc_mode == 0) ){
                if(capacity_ > 0){
                    capacity_org_ = capacity_;
                    size_ = capacity_ = 0;
                    delete [] data_;
                    capacity_ = 10;
                    data_ = new T[capacity_];
                }
                else{
                    capacity_org_ = 0;
                }
            }
            else if( alloc_mode_ == 1 && (free_alloc_mode == -1  || free_alloc_mode == 1) ){
                if(data_ == NULL) return;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                {
                    MemoryPool::freeMem(id_mpool_);
                    capacity_org_ = capacity_;
                    capacity_ = 0;
                    size_     = 0;
                    data_     = NULL;
                    id_mpool_ = -1;
                }
            }
        }

        
    };
}

