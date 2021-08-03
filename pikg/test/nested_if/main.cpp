#include <pikg_vector.hpp>
struct Ep{
  PIKG::F32vec pos;
  PIKG::S32    id;
};
struct Force{
  PIKG::F32    r2min;
  PIKG::S32    count;
};

#include "kernel.hpp"

#include <random>
#include <limits>
#include <cassert>

int main(){
  constexpr PIKG::S32 N = 1024;
  Ep* ep = new Ep[N];

  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::uniform_real_distribution<> dist(-1.0, 1.0);
  for(int i=0;i<N;i++){
    ep[i].id = i;
    ep[i].pos.x = dist(engine);
    ep[i].pos.y = dist(engine);
    ep[i].pos.z = dist(engine);
  }

  Force* f_pikg = new Force[N];
  Force* f_orig = new Force[N];
  for(int i=0;i<N;i++){
    f_orig[i].r2min = std::numeric_limits<float>::max();
    f_orig[i].count = 0;
  }
  std::cout << "-- original kernel --" << std::endl;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      const PIKG::F32vec dr = ep[i].pos - ep[j].pos;
#if 0
      const PIKG::F32 r2 = dr*dr;
#else
      PIKG::F32 r2 = dr.y*dr.y;
      r2 = std::fma(dr.x,dr.x,r2);
      r2 = std::fma(dr.z,dr.z,r2);
#endif
      if(r2 < 1.0f){
	if(ep[i].id != ep[j].id){
	  f_orig[i].r2min = std::min(r2,f_orig[i].r2min);
	  f_orig[i].count += 1;
	}	
      }
    }
  }
  std::cout << "-- pikg kernel --" << std::endl;
  Kernel kernel;

  bool isOK = true;
#ifdef SIMD_KERNEL
  for(int k=1;k<=2;k++){
    std::cout << "k=" << k << std::endl;
#else
    int k = 1;
#endif
    for(int i=0;i<N;i++){
      f_pikg[i].r2min = std::numeric_limits<float>::max();
      f_pikg[i].count = 0;
    }
    kernel(ep,N,ep,N,f_pikg,k);
    for(int i=0;i<N;i++){
      if(f_pikg[i].r2min != f_orig[i].r2min || f_pikg[i].count != f_orig[i].count){
	isOK = false;
	std::cout << i << ": " << f_pikg[i].r2min << " != " << f_orig[i].r2min << std::endl;
	std::cout << i << ": " << f_pikg[i].count << " != " << f_orig[i].count << std::endl;
      }
    }
#ifdef SIMD_KERNEL
  }
#endif
  std::cout << "-- pikg kernel end --" << std::endl;
  if(isOK){
    std::cout << "Test passed" << std::endl;
    return 0;
  }else{
    std::cout << "Test failed" << std::endl;
    return -1;
  }
}
