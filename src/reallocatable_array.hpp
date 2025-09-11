#pragma once

#include <typeinfo>
#include <cassert>
#include <iostream>
#include <fstream>

#define PARTICLE_SIMULATOR_ROUND_UP_CAPACITY

extern int N_STEP;

constexpr int NAME_SIZE = 32;

namespace ParticleSimulator {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
static omp_lock_t PS_OMP_LOCK;
#endif

template <class T>
class ReallocatableArray {
    //    private:
   public:
    // ReallocatableArray(const ReallocatableArray &){};
    ReallocatableArray &operator=(const ReallocatableArray &) {};
    T *data_;
    int size_;
    int capacity_;
    int capacity_org_;
    int id_mpool_;
    char name_[NAME_SIZE];
    MemoryAllocMode alloc_mode_;

    int my_omp_get_num_threads() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        const int n_thread = omp_get_num_threads();
#else
        const int n_thread = 1;
#endif
        // std::cout<<"n_thread= "<<n_thread<<std::endl;
        return n_thread;
    }

    int my_omp_get_thread_num() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        const int id_thread = omp_get_thread_num();
#else
        const int id_thread = 0;
#endif
        return id_thread;
    }

    class CallerInfo {
#if defined(SET_CALLER_INFO_REALLOCATABLE_ARRAY)
        int line_num_;
        char func_name_[256];
        char file_name_[256];
#endif
       public:
        void dump(std::ostream &fout = std::cerr) const {
#if defined(SET_CALLER_INFO_REALLOCATABLE_ARRAY)
            fout << "line_num_= " << line_num_ << " func_name= " << func_name_ << " file_name= " << file_name_ << std::endl;
#endif
        }
        void set(const int l = __LINE__, const char *func = __FUNCTION__, const char *file = __FILE__) {
#if defined(SET_CALLER_INFO_REALLOCATABLE_ARRAY)
            line_num_ = l;
            strcpy(func_name_, func);
            strcpy(file_name_, file);
#endif
        }
    };
    CallerInfo caller_info_;

    bool inParallelRegion() {
        /*
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
        const int n_thread = omp_get_num_threads();
#else
        const int n_thread = 1;
#endif
        */
        const int n_thread = my_omp_get_num_threads();
        if (n_thread > 1)
            return true;
        else
            return false;
    }

    void calcAdrToSplitData(int &head, int &end, const int i_th, const int n_div, const int n_tot) {
        head = (n_tot / n_div) * i_th + std::min(n_tot % n_div, i_th);
        end = (n_tot / n_div) * (i_th + 1) + std::min(n_tot % n_div, (i_th + 1));
    }

    int roundUp(const int cap_prev, const int factor = 8) { return ((cap_prev + factor - 1) / factor) * factor; }

    int getNewCapacity(const int cap_prev) {
#if defined(PARTICLE_SIMULATOR_ROUND_UP_CAPACITY)
        return roundUp(cap_prev);
#else
        return cap_prev;
#endif
    }
    void checkCapacityOverInConstructor() {
        if (capacity_ >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE) {
            PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
            int rank_tmp;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank_tmp);
            std::cerr << "rank=" << rank_tmp << std::endl;
#else
            std::cerr << "rank=0" << std::endl;
#endif
            std::cerr << "typeid(T).name():" << typeid(T).name() << std::endl;
            std::cerr << "capacity_=" << capacity_ << std::endl;
            std::cerr << "LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE=" << LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE << std::endl;
            Abort(-1);
        }
    }
    void constructorFunction(const int cap) {
        capacity_ = getNewCapacity(cap);
        checkCapacityOverInConstructor();
        strcpy(name_, "no_name");
        if (capacity_ == 0) {
            data_ = nullptr;
        } else if (alloc_mode_ == MemoryAllocMode::Default) {
            data_ = new T[capacity_];
        } else if (alloc_mode_ == MemoryAllocMode::Pool || alloc_mode_ == MemoryAllocMode::Stack) {
            // omp_set_lock(&PS_OMP_LOCK);
            // #pragma omp flush
            PS_OMP_CRITICAL {
                void *ret = nullptr;
                MemoryPool::alloc(sizeof(T) * capacity_, id_mpool_, &data_, ret, alloc_mode_);
                MemoryPool::setName(id_mpool_, name_);
                data_ = (T *)ret;
            }
            // #pragma omp flush
            // omp_unset_lock(&PS_OMP_LOCK);
        }
    }

#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
    void dumpImpl(std::ostream &fout = std::cerr) const {
        fout << "typeid(T).name()" << typeid(T).name() << std::endl;
        fout << "data_= " << data_ << std::endl;
        fout << "size_=" << size_ << std::endl;
        fout << "capacity_=" << capacity_ << std::endl;
        fout << "capacity_org_=" << capacity_org_ << std::endl;
        fout << "id_mpool_=" << id_mpool_ << std::endl;
        fout << "n_expand_=" << n_expand_ << std::endl;
        fout << "alloc_mode_= " << static_cast<int>(alloc_mode_) << std::endl;
        fout << "name_= " << name_ << std::endl;
        caller_info_.dump(fout);
    }
    int n_expand_;
    void increaseNExpand(const int n_input, std::ostream &fout = std::cerr) {
        fout << "expand capacity" << std::endl;
        fout << "typeid(T).name()" << typeid(T).name() << std::endl;
        fout << "n_input=" << n_input << std::endl;
        fout << "size_=" << size_ << std::endl;
        fout << "capacity_=" << capacity_ << std::endl;
        fout << "n_expand_=" << n_expand_ << std::endl;
        fout << "alloc_mode_= " << static_cast<int>(alloc_mode_) << std::endl;
        fout << "name_= " << name_ << std::endl;
        n_expand_++;
    }
#else
    void dumpImpl(std::ostream &fout = std::cerr) const {
        fout << "typeid(T).name()" << typeid(T).name() << std::endl;
        fout << "data_= " << data_ << std::endl;
        fout << "size_=" << size_ << std::endl;
        fout << "capacity_=" << capacity_ << std::endl;
        fout << "capacity_org_=" << capacity_org_ << std::endl;
        fout << "id_mpool_=" << id_mpool_ << std::endl;
        fout << "data_= " << data_ << std::endl;
        fout << "alloc_mode_= " << static_cast<int>(alloc_mode_) << std::endl;
        fout << "name_= " << name_ << std::endl;
        caller_info_.dump(fout);
    }
#endif

#if defined(__GNUC__)
    void __attribute__((always_inline)) ReallocInner(const int new_cap) {
#else
    void ReallocInner(const int new_cap) {
#endif

        // int rank = -1;
        // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (alloc_mode_ == MemoryAllocMode::Default) {
            capacity_ = getNewCapacity(new_cap);
            T *data_old = data_;
            T *data_new = new T[capacity_];
            if (inParallelRegion()) {
                for (int i = 0; i < size_; i++) {
                    data_new[i] = data_old[i];
                }
            } else {
                // PS_OMP_PARALLEL_FOR
                for (int i = 0; i < size_; i++) {
                    data_new[i] = data_old[i];
                }
            }
            data_ = data_new;
            // DELETE_ARRAY(data_old);
            if (data_old != nullptr) delete[] data_old;
        } else if (alloc_mode_ == MemoryAllocMode::Pool || alloc_mode_ == MemoryAllocMode::Stack) {
            if (inParallelRegion()) {
                // in parallel scope
                // omp_set_lock(&PS_OMP_LOCK);
                // #pragma omp flush
                PS_OMP_CRITICAL {
                    if (data_ != nullptr) {
                        // Save old data
                        T *data_old = new T[capacity_];
                        memcpy((void *)data_old, (void *)data_, (size_t)sizeof(T) * capacity_);
                        // Set new capacity
                        capacity_ = getNewCapacity(new_cap);
                        int id_old = id_mpool_;
                        void *ret = nullptr;
                        MemoryPool::alloc(sizeof(T) * capacity_, id_mpool_, &data_, ret, alloc_mode_);
                        MemoryPool::setName(id_mpool_, name_);
                        data_ = (T *)ret;
                        if (id_old != id_mpool_) {
                            MemoryPool::freeMem(id_old);
                        }
                        // Copy old data
                        memcpy((void *)data_, (void *)data_old, (size_t)sizeof(T) * size_);
                        // Release memory of temporal buffer
                        delete[] data_old;
                    } else {
                        // Set new capacity
                        capacity_ = getNewCapacity(new_cap);
                        void *ret = nullptr;
                        MemoryPool::alloc(sizeof(T) * capacity_, id_mpool_, &data_, ret, alloc_mode_);
                        MemoryPool::setName(id_mpool_, name_);
                        data_ = (T *)ret;
                    }
                }  // end of OMP_CRITICAL
                // #pragma omp flush
                // omp_unset_lock(&PS_OMP_LOCK);
            } else {
                // not in parallel scope
                if (data_ != nullptr) {
                    // Save old data
                    T *data_old = new T[capacity_];
                    PS_OMP_PARALLEL {
                        int head, end;
                        // int i_th = omp_get_thread_num();
                        // int n_th = omp_get_num_threads();
                        int i_th = my_omp_get_thread_num();
                        int n_th = my_omp_get_num_threads();
                        calcAdrToSplitData(head, end, i_th, n_th, capacity_);

                        // if(N_STEP == 1431 && rank == 14){
                        //     std::cout<<"ReallocInner 1"<<std::endl;
                        // }

                        memcpy((void *)(data_old + head), (void *)(data_ + head), (size_t)(sizeof(T) * (end - head)));

                        // if(N_STEP == 1431 && rank == 14){
                        //     std::cout<<"ReallocInner 2"<<std::endl;
                        // }
                    }
                    // Set new capacity
                    capacity_ = getNewCapacity(new_cap);
                    int id_old = id_mpool_;
                    void *ret = nullptr;

                    // if(N_STEP == 1431 && rank == 14){
                    //     std::cout<<"ReallocInner 3"<<std::endl;
                    //     MemoryPool::dump(std::cout);
                    // }

                    MemoryPool::alloc(sizeof(T) * capacity_, id_mpool_, &data_, ret, alloc_mode_);

                    // if(N_STEP == 1431 && rank == 14){
                    //     std::cout<<"ReallocInner 4"<<std::endl;
                    //     MemoryPool::dump(std::cout);
                    // }

                    MemoryPool::setName(id_mpool_, name_);

                    // if(N_STEP == 1431 && rank == 14){
                    //     std::cout<<"ReallocInner 5"<<std::endl;
                    //     MemoryPool::dump(std::cout);
                    // }

                    data_ = (T *)ret;

                    // if(N_STEP == 1431 && rank == 14){
                    //     std::cout<<"ReallocInner 6"<<std::endl;
                    //     MemoryPool::dump(std::cout);
                    // }

                    if (id_old != id_mpool_) {
                        // if(N_STEP == 1431 && rank == 14){
                        //     std::cout<<"ReallocInner 7"<<std::endl;
                        //     MemoryPool::dump(std::cout);
                        // }

                        MemoryPool::freeMem(id_old);

                        // if(N_STEP == 1431 && rank == 14){
                        //     std::cout<<"ReallocInner 8"<<std::endl;
                        //     MemoryPool::dump(std::cout);
                        // }
                    }
                    // Copy old data
                    PS_OMP_PARALLEL {
                        int head, end;
                        // int i_th = omp_get_thread_num();
                        // int n_th = omp_get_num_threads();
                        int i_th = my_omp_get_thread_num();
                        int n_th = my_omp_get_num_threads();
                        calcAdrToSplitData(head, end, i_th, n_th, size_);
                        memcpy((void *)(data_ + head), (void *)(data_old + head), (size_t)(sizeof(T) * (end - head)));
                    }
                    // Release memory of temporal buffer
                    delete[] data_old;
                } else {
                    // Set new capacity
                    capacity_ = getNewCapacity(new_cap);
                    // int rank = -1;
                    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                    // if(N_STEP == 1431 && rank == 14){
                    //     std::cout<<"memory pool A"<<std::endl;
                    //     MemoryPool::dump(std::cout);
                    // }
                    void *ret = nullptr;
                    MemoryPool::alloc(sizeof(T) * capacity_, id_mpool_, &data_, ret, alloc_mode_);

                    // if(N_STEP == 1431 && rank == 14){
                    //     std::cout<<"memory pool B"<<std::endl;
                    //     MemoryPool::dump(std::cout);
                    // }

                    MemoryPool::setName(id_mpool_, name_);
                    data_ = (T *)ret;
                }
            }
        }
    }

   public:
    ///////////////
    // constructor
    ReallocatableArray(int cap, int size, MemoryAllocMode alloc_mode)
        : data_(nullptr), size_(size), capacity_org_(0), id_mpool_(-1), alloc_mode_(alloc_mode) {
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        n_expand_ = 0;
#endif
#if defined(USE_ALLOC_MODE_NEW)
        alloc_mode_ = MemoryAllocMode::Default;
#endif
        constructorFunction(cap);
    }
    ReallocatableArray(MemoryAllocMode alloc_mode) : ReallocatableArray(0, 0, alloc_mode) {}
    ReallocatableArray() : ReallocatableArray(0, 0, MemoryAllocMode::Default) {}
    ~ReallocatableArray() {
        // std::cerr<<"****** destroctor ******"<<std::endl;
        // dumpImpl();
        if (alloc_mode_ == MemoryAllocMode::Default) {
            if (capacity_ > 0) delete[] data_;
        } else if (alloc_mode_ == MemoryAllocMode::Pool || alloc_mode_ == MemoryAllocMode::Stack) {
            // omp_set_lock(&PS_OMP_LOCK);
            // #pragma omp flush
            PS_OMP_CRITICAL { MemoryPool::freeMem(id_mpool_); }
            // #pragma omp flush
            // omp_unset_lock(&PS_OMP_LOCK);
        }
        data_ = nullptr;
    }

    void setName(const char *str) {
        if (strlen(str) < NAME_SIZE) {
            strcpy(name_, str);
            /*
            if(alloc_mode_ == MemoryAllocMode::Pool){
                MemoryPool::setName(id_mpool_, str);
            }
            */
        }
    }

    void setAllocMode(const MemoryAllocMode alloc_mode) {
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        assert(data_ == nullptr);
#endif
        alloc_mode_ = alloc_mode;
#if defined(USE_ALLOC_MODE_NEW)
        alloc_mode_ = MemoryAllocMode::Default;
#endif
    }
    void setCallerInfo(const int l = __LINE__, const char *func = __FUNCTION__, const char *file = __FILE__) { caller_info_.set(l, func, file); }

    void initialize(int cap, int size, MemoryAllocMode alloc_mode) {
        // std::cerr<<"MemoryPool::getCapacity()= "<<MemoryPool::getCapacity()<<std::endl;
        capacity_ = getNewCapacity(cap);
        size_ = size;
        alloc_mode_ = alloc_mode;
        constructorFunction(capacity_);
    }

    void reserve(const int n) {
        if (n > capacity_) {
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 2
            increaseNExpand(n);
            std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#endif
            capacity_ = getNewCapacity(n);
            if (capacity_ >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE) {
                PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                int rank_tmp;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank_tmp);
                std::cerr << "rank= " << rank_tmp << std::endl;
#else
                std::cerr << "rank= 0" << std::endl;
#endif
                std::cerr << "typeid(T).name():" << typeid(T).name() << std::endl;
                std::cerr << "capacity_=" << capacity_ << std::endl;
                std::cerr << "LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE=" << LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE << std::endl;
                ParticleSimulator::Abort(-1);
            }
            ReallocInner(n);
        }
    }

    int size() const { return size_; }

    int capacity() const { return capacity_; }

    const T &operator[](const int i) const {
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        if (i > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE || i < 0 || capacity_ <= i) {
            dumpImpl();
            std::cerr << "i=" << i << std::endl;
            std::cerr << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#if defined(THROW_EXCEPTION_REALLOCATABLE_ARRAY)
            std::stringstream msg;
            msg << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
            throw msg.str().c_str();
#endif
        }
        assert(i <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
        assert(i >= 0);
        assert(capacity_ > i);
#endif
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 1
        if (size_ <= i) {
            dumpImpl();
            std::cout << "i=" << i << std::endl;
            std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#if defined(THROW_EXCEPTION_REALLOCATABLE_ARRAY)
            std::stringstream msg;
            msg << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
            throw msg.str().c_str();
#endif
        }
        assert(size_ > i);
#endif
        return data_[i];
    }

    T &operator[](const int i) {
// #ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        if (i > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE || i < 0 || capacity_ <= i) {
            dumpImpl();
            std::cerr << "i=" << i << std::endl;
            std::cerr << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#if defined(THROW_EXCEPTION_REALLOCATABLE_ARRAY)
            std::stringstream msg;
            msg << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
            throw msg.str().c_str();
#endif
        }
        assert(i <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
        assert(i >= 0);
        assert(capacity_ > i);
#endif
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 1
        if (size_ <= i) {
            dumpImpl();
            std::cerr << "i=" << i << std::endl;
            std::cerr << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#if defined(THROW_EXCEPTION_REALLOCATABLE_ARRAY)
            std::stringstream msg;
            msg << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
            throw msg.str().c_str();
#endif
        }
        assert(size_ > i);
#endif
        return data_[i];
    }

    T &front() { return data_[0]; }

    T &back() { return data_[size_ - 1]; }

    T *frontP() { return data_; }

    T *backP() { return data_ + (size_ - 1); }

    T *endP() { return data_ + size_; }

    T *data() { return data_; }

    const T *data() const { return data_; }

    void push_back(const T &val) {
// #ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        if (size_ + 1 > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE) {
            dumpImpl();
            std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#if defined(THROW_EXCEPTION_REALLOCATABLE_ARRAY)
            std::stringstream msg;
            msg << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
            throw msg.str().c_str();
#endif
        }
        assert(size_ + 1 <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
#endif
        const int size_new = size_ + 1;
        resizeNoInitialize(size_new);
        data_[size_ - 1] = val;
    }

#if defined(__GNUC__)
    void __attribute__((always_inline)) resizeNoInitialize(const int n) {
#else
    void resizeNoInitialize(const int n) {
#endif
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        if (n > LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE) {
            dumpImpl();
            std::cout << "n=" << n << std::endl;
            std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#if defined(THROW_EXCEPTION_REALLOCATABLE_ARRAY)
            std::stringstream msg;
            msg << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
            throw msg.str().c_str();
#endif
        }
        assert(n <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
#endif
        const int new_cap = std::max(getNewCapacity((n + (n + 10) / 10) + 100), capacity_org_);
        // std::cout<<std::dec<<"new_cap= "<<new_cap<<" data_= "<<data_<<std::endl;
        if (data_ == nullptr) {
            if (alloc_mode_ == MemoryAllocMode::Default) {
                data_ = new T[new_cap];
            } else if (alloc_mode_ == MemoryAllocMode::Pool || alloc_mode_ == MemoryAllocMode::Stack) {
                // omp_set_lock(&PS_OMP_LOCK);
                // #pragma omp flush
                PS_OMP_CRITICAL {
                    void *ret = nullptr;
                    // std::cout<<"check 1"<<std::endl;
                    MemoryPool::alloc(sizeof(T) * new_cap, id_mpool_, &data_, ret, alloc_mode_);
                    // std::cout<<"check 2"<<std::endl;
                    MemoryPool::setName(id_mpool_, name_);
                    // std::cout<<"check 3"<<std::endl;
                    data_ = (T *)ret;
                }
                // #pragma omp flush
                // omp_unset_lock(&PS_OMP_LOCK);
            }
            capacity_ = new_cap;
            size_ = n;
        }
        if (n > capacity_) {
            // std::cerr<<"n= "<<n<<" capacity_= "<<capacity_<<std::endl;
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 2
            increaseNExpand(n);
            std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#endif
            if (new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE) {
                PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                int rank_tmp;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank_tmp);
                std::cerr << "rank=" << rank_tmp << std::endl;
#else
                std::cerr << "rank=0" << std::endl;
#endif
                std::cerr << "typeid(T).name():" << typeid(T).name() << std::endl;
                std::cerr << "capacity_=" << capacity_ << std::endl;
                std::cerr << "LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE=" << LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE << std::endl;
                Abort(-1);
            }
            ReallocInner(new_cap);
        }
        size_ = n;
    }

    void dump(const std::string str = "", std::ostream &fout = std::cerr) const {
        fout << str << std::endl;
        dumpImpl(fout);
    }

    size_t getMemSize() const { return capacity_ * sizeof(T); }

    T *getPointer(const int i = 0) const { return data_ + i; }

    void pushBackNoCheck(const T &val) {
// #ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 0
        assert(size_ <= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE);
        if (size_ >= capacity_) {
            dumpImpl();
            std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#if defined(THROW_EXCEPTION_REALLOCATABLE_ARRAY)
            std::stringstream msg;
            msg << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
            throw msg.str().c_str();
#endif
        }
        assert(size_ < capacity_);
#endif
        data_[size_++] = val;
    }
    void clearSize() { size_ = 0; }
    void increaseSize(const int n = 1) { resizeNoInitialize(size_ + n); }
    void decreaseSize(const int n = 1) { resizeNoInitialize(size_ - n); }  // no check
    void reserveAtLeast(const int n) {
        if (n >= capacity_) {
// #ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 2
            increaseNExpand(n);
            std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#endif
            const int new_cap = (n + (n + 10) / 10) + 100;
            if (new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE) {
                PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                int rank_tmp;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank_tmp);
                std::cerr << "rank=" << rank_tmp << std::endl;
#else
                std::cerr << "rank=0" << std::endl;
#endif
                std::cerr << "typeid(T).name():" << typeid(T).name() << std::endl;
                std::cerr << "capacity_=" << capacity_ << std::endl;
                std::cerr << "LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE=" << LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE << std::endl;
                Abort(-1);
            }
            ReallocInner(new_cap);
        }
    }
#if defined(__GNUC__)
    void __attribute__((always_inline)) reserveEmptyAreaAtLeast(const int n_add){
#else
    void reserveEmptyAreaAtLeast(const int n_add) {
#endif
        const int n = n_add + size_;
    if (n >= capacity_) {
#if SANITY_CHECK_REALLOCATABLE_ARRAY > 2
        increaseNExpand(n);
        std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ << ", file: " << __FILE__ << std::endl;
#endif
        const int new_cap = (n + (n + 10) / 10) + 100;
        if (new_cap >= LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE) {
            PARTICLE_SIMULATOR_PRINT_ERROR("The number of particles of this process is beyound the FDPS limit number");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            int rank_tmp;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank_tmp);
            std::cerr << "rank=" << rank_tmp << std::endl;
#else
                std::cerr << "rank=0" << std::endl;
#endif
            std::cerr << "typeid(T).name():" << typeid(T).name() << std::endl;
            std::cerr << "capacity_=" << capacity_ << std::endl;
            std::cerr << "LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE=" << LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE << std::endl;
            Abort(-1);
        }
        ReallocInner(new_cap);
    }
}
        
        void freeMem(const int free_alloc_mode = -1){
    if (alloc_mode_ == MemoryAllocMode::Default) {
        if (capacity_ > 0) {
            capacity_org_ = capacity_;
            size_ = capacity_ = 0;
            delete[] data_;
            data_ = nullptr;
        } else {
            capacity_org_ = 0;
        }
    } else if (alloc_mode_ == MemoryAllocMode::Pool || alloc_mode_ == MemoryAllocMode::Stack) {
        if (data_ == nullptr) return;
        // omp_set_lock(&PS_OMP_LOCK);
        // #pragma omp flush
        PS_OMP_CRITICAL {
            MemoryPool::freeMem(id_mpool_);
            capacity_org_ = capacity_;
            capacity_ = 0;
            size_ = 0;
            data_ = nullptr;
            id_mpool_ = -1;
        }
        // #pragma omp flush
        // omp_unset_lock(&PS_OMP_LOCK);
    }
}
};  // namespace ParticleSimulator

// 複数のオブジェクトのfreemem関数を呼ぶ．
template <class T, class... Args>
void FreeMemVA(T & obj, Args&&... args){
    obj.freeMem();
    FreeMemVA(args...);
}
template <class T>
void FreeMemVA(T & obj){
    obj.freeMem();
}

// 複数のオブジェクトのresizeNoInitialize関数を呼び，引数nを渡す．
template <class T, class... Args>
void ResizeNoInitializeVA(T &obj, int n, Args &&...args) {
    obj.resizeNoInitialize(n);
    ResizeNoInitializeVA(args...);
}
template <class T>
void ResizeNoInitializeVA(T &obj, int n) {
    obj.resizeNoInitialize(n);
}

// 複数のオブジェクトのclearSize関数を呼ぶ．
template <class T, class... Args>
void ClearSizeVA(T &obj, Args &&...args) {
    obj.clearSize();
    ClearSizeVA(args...);
}
template <class T>
void ClearSizeVA(T &obj) {
    obj.clearSize();
}

}  // namespace ParticleSimulator
