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
    public:
        ReallocatableArray() : data_(NULL), size_(0), capacity_(0) {}
        ReallocatableArray(int cap) : size_(0), capacity_(cap) {
            data_ = new T[capacity_];
        }
        ReallocatableArray(int size, const T & val) : size_(size), capacity_(size+1) {
            data_ = new T[capacity_];
            for(int i=0; i<capacity_; i++){
                data_[i] = val;
            }
        }
        ReallocatableArray(const ReallocatableArray & from){
            const int size_tmp = from.size_;
            for(int i=0; i<size_tmp; i++){
                data_[i] = from.data_[i];
            }
        }
        ~ReallocatableArray(){
            delete [] data_;
            data_ = NULL;
        }
        void reserve(const int n){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
	    assert(n <= 1000000);
#endif
            if( n > capacity_){
                capacity_ = n;
                T * new_data = new T[capacity_];
                T * old_data = data_;
                for( int i=0; i<size_; i++){
                    new_data[i] = old_data[i];
                }
                data_ = new_data;
                delete[] old_data;
            }
        }
        int size() const {return size_; }
        int capacity() const {return capacity_; }
        const T & operator [] (const int i) const {
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
	    /*
	    if(i > 1000000 || i < 0 || capacity_ <= i){
		std::cout<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"size_="<<size_<<std::endl;
                std::cout<<"capacity_="<<capacity_<<std::endl;
	    }
	    assert(i <= 1000000);
	    assert(i >= 0);
	    assert(capacity_ > i);
	    */

            if(i > 1000000 || i < 0 || capacity_ <= i || size_<= i ){
                std::cout<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"size_="<<size_<<std::endl;
                std::cout<<"capacity_="<<capacity_<<std::endl;
            }
	    assert(i <= 1000000);
	    assert(i >= 0);
	    assert(capacity_ > i);
	    assert(size_ > i);

#endif
            return data_[i]; 
        }
        T & operator [] (const int i){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
	    /*
	    if(i > 1000000 || i < 0 || capacity_ <= i){
		std::cout<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"size_="<<size_<<std::endl;
                std::cout<<"capacity_="<<capacity_<<std::endl;
	    }
	    assert(i <= 1000000);
	    assert(i >= 0);
	    assert(capacity_ > i);
	    */
            if(i > 1000000 || i < 0 || capacity_ <= i || size_<= i ){
                std::cout<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"size_="<<size_<<std::endl;
                std::cout<<"capacity_="<<capacity_<<std::endl;
            }
	    assert(i <= 1000000);
	    assert(i >= 0);
	    assert(capacity_ > i);
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
	    assert(size_+1 <= 1000000);
#endif
            const int size_new = size_ + 1;
            resizeNoInitialize(size_new);
            data_[size_-1] = val;
        }
	void resizeNoInitialize(const int n){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
	    assert(n <= 1000000);
#endif
            if(n > capacity_){
                capacity_ = n * 2;
                T * new_data = new T[capacity_];
                T * old_data = data_;
                for( int i=0; i<size_; i++){
                    new_data[i] = old_data[i];
                }
                data_ = new_data;
                delete[] old_data;
            }
            size_ = n;
        }
        void dump(const std::string str=""){
            std::cout<<str<<std::endl;
            std::cout<<"size_="<<size_<<std::endl;
            std::cout<<"capacity_="<<capacity_<<std::endl;
        }
        size_t getMemSize() const { return capacity_ * sizeof(T); }
        T * getPointer(const int i=0) const { return data_+i; }
        void pushBackNoCheck(const T & val){
#ifdef SANITY_CHECK_REALLOCATABLE_ARRAY
	    assert(size_ <= 1000000);
            if(size_ >= capacity_){
                std::cout<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
                std::cout<<"size_="<<size_<<std::endl;
                std::cout<<"capacity_="<<capacity_<<std::endl;
            }
            assert(size_ < capacity_);
#endif
	    data_[size_++] = val;
        }
	//void setSize(const int n){ size_ = n; }
	void clearSize() {size_ = 0;}
	void increaseSize(){
	    resizeNoInitialize(size_+1);
	}
    };
}
