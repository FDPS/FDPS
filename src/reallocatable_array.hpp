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
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        void dumpImpl() const {
            std::cout<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            std::cout<<"size_="<<size_<<std::endl;
            std::cout<<"capacity_="<<capacity_<<std::endl;
            std::cout<<"n_expand_="<<n_expand_<<std::endl;
        }
        int n_expand_;
        void increaseNExpand(const int n_input){
            std::cout<<"expand capacity"<<std::endl;
            std::cout<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            std::cout<<"n_input="<<n_input<<std::endl;
            std::cout<<"size_="<<size_<<std::endl;
            std::cout<<"capacity_="<<capacity_<<std::endl;
            std::cout<<"n_expand_="<<n_expand_<<std::endl;
            n_expand_++;
        }
#endif

#ifdef __GNUC__
        void __attribute__((always_inline)) ReallocInner(const int new_cap){
#else
        void ReallocInner(const int new_cap){
#endif
            capacity_ = new_cap;
            T * new_data = new T[capacity_];
            T * old_data = data_;
            for( int i=0; i<size_; i++){
                new_data[i] = old_data[i];
            }
            data_ = new_data;
            delete[] old_data;
        }
	
    public:
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        ReallocatableArray() : data_(NULL), size_(0), capacity_(0), capacity_org_(0), n_expand_(0) {}
        ReallocatableArray(int cap) : size_(0), capacity_(cap), capacity_org_(0), n_expand_(0) {
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
            data_ = new T[capacity_];
        }
#else
        ReallocatableArray() : data_(NULL), size_(0), capacity_(0), capacity_org_(0) {}
        ReallocatableArray(int cap) : size_(0), capacity_(cap), capacity_org_(0)  {
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
            data_ = new T[capacity_];
        }
#endif
        ~ReallocatableArray(){
            if(capacity_ > 0) delete [] data_;
            data_ = NULL;
        }
        void reserve(const int n){
            if( n > capacity_){
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
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
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
            if(n > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                dumpImpl();
                std::cout<<"n="<<n<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(n <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
#endif
            if(n > capacity_){
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
            size_ = n;
        }

        void dump(const std::string str=""){
//#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
            std::cout<<str<<std::endl;
            dumpImpl();
#endif
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

        void freeMem(){
            if(capacity_ > 0){
                capacity_org_ = capacity_;
                size_ = capacity_ = 0;
                delete [] data_;
                capacity_ = 100;
                data_ = new T[capacity_];
            }
            else{
                capacity_org_ = 0;
            }
        }

        void reallocMem(){
            if(capacity_org_ > 0){
                capacity_ = capacity_org_;
                data_ = new T[capacity_];
                capacity_org_ = 0;
            }
        }

        void setDataPointer(const void * _data){
            if(capacity_ > 0) delete [] data_;
            data_ = (T*)_data;
        }
        
    };
}

