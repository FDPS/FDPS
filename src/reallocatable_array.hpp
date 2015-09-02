#include<typeinfo>
#include<cassert>
#include<iostream>

namespace  ParticleSimulator{
    template<class T>
        class ReallocatableArray{
    private:
        T * data_;
        int size_;
        int capacity_;
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
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

        void dumpImpl() const {
            std::cout<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            std::cout<<"size_="<<size_<<std::endl;
            std::cout<<"capacity_="<<capacity_<<std::endl;
            std::cout<<"n_expand_="<<n_expand_<<std::endl;
        }
#endif

        void ReallocInner(const int new_cap){
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
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
        ReallocatableArray() : data_(NULL), size_(0), capacity_(0), n_expand_(0) {}
        ReallocatableArray(int cap) : size_(0), capacity_(cap), n_expand_(0) {
            if(capacity_ >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
                //std::cerr<<"rank="<<Comm::getRank()<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                std::cerr<<"rank="<<MPI::COMM_WORLD.Get_rank()<<std::endl;
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
        ReallocatableArray() : data_(NULL), size_(0), capacity_(0) {}
        ReallocatableArray(int cap) : size_(0), capacity_(cap) {
            if(capacity_ >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
                //std::cerr<<"rank="<<Comm::getRank()<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                std::cerr<<"rank="<<MPI::COMM_WORLD.Get_rank()<<std::endl;
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
            delete [] data_;
            data_ = NULL;
        }
        void reserve(const int n){
            if( n > capacity_){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
                increaseNExpand(n);
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
                capacity_ = n;
                if(capacity_ >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
                    //std::cerr<<"rank="<<Comm::getRank()<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    std::cerr<<"rank="<<MPI::COMM_WORLD.Get_rank()<<std::endl;
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
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
            if(i > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE || i < 0 || capacity_ <= i || size_<= i ){
                dumpImpl();
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(i <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
            assert(i >= 0);
            assert(capacity_ > i);
            assert(size_ > i);
#endif
            return data_[i]; 
        }

        T & operator [] (const int i){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
            if(i > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE || i < 0 || capacity_ <= i || size_<= i ){
                dumpImpl();
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(i <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
            assert(i >= 0);
            assert(capacity_ > i);
            if(size_ <= i){
                std::cerr<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
                std::cerr<<"i="<<i<<std::endl;
                std::cerr<<"size_="<<size_<<std::endl;
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
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
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

        void resizeNoInitialize(const int n){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
            if(n > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                dumpImpl();
                std::cout<<"n="<<n<<std::endl;
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
            }
            assert(n <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
#endif
            if(n > capacity_){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
                increaseNExpand(n);
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
                const int new_cap = n * 1.3 + 100;
                if(new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
                    //std::cerr<<"rank="<<Comm::getRank()<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    std::cerr<<"rank="<<MPI::COMM_WORLD.Get_rank()<<std::endl;
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
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
            std::cout<<str<<std::endl;
            dumpImpl();
#endif
        }

        size_t getMemSize() const { return capacity_ * sizeof(T); }

        T * getPointer(const int i=0) const { return data_+i; }

        void pushBackNoCheck(const T & val){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
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

        void increaseSize(){ resizeNoInitialize(size_+1); }

        void reserveAtLeast(const int n){
            if( n >= capacity_){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
                increaseNExpand(n);
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
                const int new_cap = n * 1.3 + 100;
                if(new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
                    //std::cerr<<"rank="<<Comm::getRank()<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    std::cerr<<"rank="<<MPI::COMM_WORLD.Get_rank()<<std::endl;
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

        void reserveEmptyAreaAtLeast(const int n_add){
            const int n = n_add + size_;
            if( n >= capacity_){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
                increaseNExpand(n);
                std::cout<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
#endif
                const int new_cap = n * 1.3 + 100;
                if(new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
                    //std::cerr<<"rank="<<Comm::getRank()<<std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    std::cerr<<"rank="<<MPI::COMM_WORLD.Get_rank()<<std::endl;
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

    };
}
