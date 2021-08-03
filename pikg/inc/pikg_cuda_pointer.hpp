#ifndef PIKG_CUDA_POINTER
#define PIKG_CUDA_POINTER

#include <cassert>
#include <cstdlib>
#include <pikg_vector.hpp>

namespace PIKG{
  template <typename T>
  class CUDAPointer{
  public:
    T* hst = nullptr;
    T* dev = nullptr;
    U64 _size = 0;
    U64 _limit = 0;

    ~CUDAPointer(){
      free();
    }
    void free_host(){
      if(hst != nullptr){
	delete[] hst;
	hst = nullptr;
      }
    }
    void free_dev(){
      if(dev != nullptr){
	cudaFree(dev);
	dev = nullptr;
      }
    }
    void free(){
      free_host();
      free_dev();
    }

    void allocate(const U64 n){
      //free();
      hst = new T[n];
      cudaMalloc((void**)&dev,n*sizeof(T));
      _size = n;
      _limit = n;

      //assert(hst != nullptr);
      assert(dev != nullptr);
    }
    void resize(const U64 n){
      reserve(n);
      _size = n;
    }

    void reserve(const U64 n){
      if(_size==0 || _limit < n) allocate(n);
    }

    void h2d(const U64 n){
      cudaMemcpy(dev,hst,n*sizeof(T),cudaMemcpyHostToDevice);
    }
    void d2h(const U64 n) const {
      cudaMemcpy(hst,dev,n*sizeof(T),cudaMemcpyDeviceToHost);
    }

    const T& operator[](const int index) const {
      return hst[index];
    }
    T& operator[](const int index){
      return hst[index];
    }

    operator T*(){
      return dev;
    }
  };
};

#endif // PIKG_CUDA_POINTER
