#pragma once
#include <unistd.h>
//
// a library for parallel sample sort
//
// usage: #include "PATH/OF/THIS/FILE/samplesortlib.hpp"
//
//   Provide:
//   template<class T>
//   void samplesort_bodies(T * b, int n);
//   and
//   template<class T, class GetKey>
//   void samplesort_bodies(T * b,
//		       int n,
//                     GetKey getkey);
//   in namespace SampleSortLib
//
//  b is the pointer to array (could be vector?) of class T
//  n is the number of elements
//  Class T should provide a member function getsortkey() 
//  using which this function sort the array. (1st form)
//  Class GetKey is class for a function which returns the sort key value.
//  taking class T as its argument (2nd form)
//
//    Copyright   2021 Jun Makino
//    Licencse: MIT
//    https://github.com/jmakino/sortlib/blob/main/LICENSE
//

#ifdef _OPENMP
#define SAMPLESORTLIB_OMP_GET_MAX_THREADS omp_get_max_threads()
#define SAMPLESORTLIB_OMP_GET_THREAD_NUM omp_get_thread_num()
#else
#define SAMPLESORTLIB_OMP_GET_MAX_THREADS 1
#define SAMPLESORTLIB_OMP_GET_THREAD_NUM 0
#endif
namespace SampleSortLib{

    
    
    inline double GetWtime()
    {
	struct timespec ts;
	if (clock_gettime(CLOCK_MONOTONIC,&ts )){
	    printf("GetWtime Failed\n");
	}
	return ((double) ts.tv_sec)+ ts.tv_nsec*1e-9;
    }

    
    template<class T>
    void sort_bodies(T * b,
		     int n)
    {
	
	std::sort(b, b+n, 
		  [](T & l, T & r )
		  ->bool{return l.getsortkey() < r.getsortkey();} );
    }
    
    template<class T, class Compare>
    void sort_bodies(T * b,
		     int n,
		     Compare c)
    {
	
	std::sort(b, b+n, c);
    }
    
    //#define SORTLIB_MEASURE_TIME
    void showdt(const char *s,
		bool initialize=false)
    {
#ifdef SORTLIB_MEASURE_TIME
	static double t0;
    
	if (!initialize){
	    printf("%s dt=%g\n",s, GetWtime()-t0);
	}
	t0=GetWtime();
#endif    
    }


    template<class T, class GetKey>
    inline bool compare(T &a, T& b, GetKey getkey)
    {
	return getkey(a)<getkey(b);
    }
    
    
    const int randomize_offset = 48271;     // This is a prime number

    template<class T, class GetKey>
    void samplesort(T * b,
		    int n,
		    GetKey getkey)
    {
	typedef decltype(getkey(*b)) KeyType;
	class KeyValuePair{
	public:
	    int64_t key;
	    KeyType value;
	    inline auto getsortkey() {return value;}    
	};
    
	showdt("", true);
	auto nt = SAMPLESORTLIB_OMP_GET_MAX_THREADS;
	//    auto nt = 4;
	//	if (n < nt*5 || n < 2*nt*nt|| nt==1){
	if (n < nt*5 || n<1000|| nt==1){
	    //	printf("single thread sort called\n");
	    // n too small for parallization. Call single-thread sort
	    std::sort(b,b+n,
		      [&](T & l, T & r )
		      ->bool{return compare(l,r, getkey);} );
	    
	    return;
	}
	auto nwork0= (n+nt-1)/nt;
	auto nworkm= n/nt;
	auto nworkp= nworkm+1;
	auto ntplus = n%nt;
	auto nsample= nwork0/nt/10;
	if (nsample < 2) nsample = 2;
	auto nsampletotal = nsample * nt;
	KeyType samplearray[nsampletotal];
	KeyValuePair key[nt][nwork0];
	//    printf("nt=%d  nwork=%d nsample=%d \n", nt, nwork0, nsample);
	int isrcstart[nt][nt+1];
	int destsize[nt];
	int deststart[nt];
	auto maxdestsize=0;
#pragma omp parallel
	{
	    auto it=SAMPLESORTLIB_OMP_GET_THREAD_NUM;
	    if (it==0) showdt("Enter parallel region");
#if 0
	    auto mystartdest = it*nsample;
	    auto mystartlocal = it*nwork0;
	    auto myendlocal =mystartlocal+nwork0;
	    if (myendlocal > n) myendlocal=n;
	    auto myrange = myendlocal - mystartlocal;
#else
	    auto mystartdest = it*nsample;
	    int mystartlocal;
	    int myrange;
	    if (it < ntplus) {
		mystartlocal = it*nworkp;
		myrange=nworkp;
	    }else{
		mystartlocal = ntplus*nworkp + (it-ntplus)*nworkm;
		myrange=nworkm;
	    }
	    
	    auto myendlocal =mystartlocal + myrange;
#endif	    
	    //	    printf("it=%d  mystart=%d, myend=%d, myrange=%d\n",
	    //	   it, mystartlocal, myendlocal,myrange );
	    auto myoffset = randomize_offset ;
	    if (myoffset % myrange == 0 ||myrange % myoffset == 0){
		myoffset ++;
	    }
	    for(auto ilocal=0; ilocal<nsample; ilocal++){
		auto myi = ((ilocal + nt)*myoffset)%myrange;
		//	    printf("it= %d ilocal=%d myi=%d dest=%d source[%d]=%g\n",
		//		   it, ilocal, myi, mystartdest+ilocal,
		//		   mystartlocal+myi,b[mystartlocal+myi].getsortkey());
	    
		samplearray[mystartdest+ilocal]= getkey(b[mystartlocal+myi]);
	    }
#pragma omp barrier	
	    if (it==0) showdt("Sample made");
	    //	printf("sort size=%d %d\n", myrange,it);
	    for(auto i=0;i<myrange; i++){
		key[it][i].key =i+mystartlocal;
		key[it][i].value =getkey(b[i+mystartlocal]);
	    }
	    //	sort_bodies(b+mystartlocal, myrange);
	    if (it==0) showdt("Key made");
	    sort_bodies(key[it], myrange);
	    if (it==0) showdt("Initial sort end");
	    //	printf("sort end %d\n", it);
	    if (it==0){
		showdt("Intial sampling and internal sort");
		//		printf("sample sort size=%d\n", nsampletotal);
		std::sort(samplearray, samplearray+nsampletotal, 
			  [](KeyType  l, KeyType  r )   ->bool{return l < r;} );
		// for(auto i=0;i<nsampletotal; i+=nsample){
		// 	printf("sample[%d]=%g\n", i, samplearray[i]);
		// }
		showdt(" sample sort");
	    }
#pragma omp barrier	
#if 0	
	}
#pragma omp parallel for schedule(static)
	for(auto it=0; it<nt; it++){
	    auto mystartdest = it*nsample;
	    auto mystartlocal = it*nwork0;
	    auto myendlocal =mystartlocal+nwork0;
	    if (myendlocal > n) myendlocal=n;
	    auto myrange = myendlocal - mystartlocal;
#endif	
	    //	auto bstart =  b+mystartlocal;
	    auto bstart =  key[it];
	    if (myendlocal > n) myendlocal=n;
	    isrcstart[it][0]=0;
	    isrcstart[it][nt]=myrange;
	    // for(auto i=0; i<myrange; i++){
	    //     printf(" %g", bstart[i].getsortkey());
	    // }
	    // printf("\n");
	    for (auto itdest=1; itdest<nt; itdest++){
		
		// find the next start location using bisection
		auto nextval = samplearray[itdest*nsample];
		int start = isrcstart[it][itdest-1];
		//	    if (bstart[start].getsortkey() < nextval) start ++;
		int end = myrange-1;
		// printf("%d %d %d %d val=%g start,end = %g %g\n",
		// 	   it, itdest, start, end,
		// 	   nextval,
		// 	   bstart[start].getsortkey(),
		// 	   bstart[end].getsortkey()  );
		if (isrcstart[it][itdest-1]==myrange){
		    //		printf("end already reached\n");
		    isrcstart[it][itdest]=myrange;
		
		}else if (bstart[start].getsortkey() >= nextval){
		    //		printf("no element\n");
		    isrcstart[it][itdest]=isrcstart[it][itdest-1];
		}else if (bstart[end].getsortkey() <nextval){
		    //		printf("end reached\n");
		    isrcstart[it][itdest]=myrange;
		}else{
		    while (start < end-1){
			auto k= (start + end)/2;
			//		    printf("%d %d %d %d %d\n", it, itdest, start, end, k );
			if (bstart[k].getsortkey() < nextval){
			    start=k;
			}else{
			    end=k;
			}
		    }
		    //		printf("final: %d %d %d %d\n", it, itdest, start, end);
		    isrcstart[it][itdest]=end;
		}
	    }
	    if (nt == 0){
		showdt(" range determination");
	    
		// for(auto it=0; it<nt; it++){
		// 	auto mystartlocal = it*nwork0;
		// 	auto  bstart =  key[it];
		// 	//	printf("isrcstart[%d]:", it);
		// 	for (auto itdest=0; itdest<nt; itdest++){
		// 	       	    printf(" %d", isrcstart[it][itdest]);
		// 	    if (itdest!=0){
		// 	    	printf(" %g %g", bstart[isrcstart[it][itdest]-1].getsortkey(),
		// 	    	       bstart[isrcstart[it][itdest]].getsortkey()    );
		// 	    }
	    
		// 	}
		// 	   	printf("\n");
		// }

	    }
#pragma omp barrier	
	    auto mydestsize =0;
	    for(auto itsrc=0; itsrc < nt; itsrc++){
		mydestsize += isrcstart[itsrc][it+1]-isrcstart[itsrc][it];
	    }
	    destsize[it]=mydestsize;
#pragma omp barrier	
	    if (it==0){
		deststart[0]=0;
		for(auto i=1;i<nt;i++) deststart[i] = deststart[i-1]+destsize[i-1];
		// for (auto i=0;i<nt;i++){
		// 	printf("dest start, size[%d]= %d, %d\n", i, deststart[i], destsize[i]);
		// }
		for(auto i=0;i<nt; i++) {
		    if (maxdestsize < destsize[i]) maxdestsize = destsize[i];
		}
	    }
	}

	KeyValuePair localcopy[nt][maxdestsize];
	T bodylocalcopy[nt][maxdestsize];
	//    showdt(" range determination2");
#pragma omp parallel
	{
	    auto it=SAMPLESORTLIB_OMP_GET_THREAD_NUM;
	    auto mydestsize = destsize[it];
	    int idest=0;
	    for(auto itsrc=0; itsrc < nt; itsrc++){
		auto bstart = key[itsrc];
		for (auto isrc=isrcstart[itsrc][it];
		     isrc < isrcstart[itsrc][it+1]; isrc++){
		    //		printf("to local: %d src= %d dest=%d\n",isrc, idest, isrc);
		    localcopy[it][idest] = bstart[isrc];
		    idest++;
		}
	    }
	    //	printf("sort size=%d %d\n", mydestsize,it);
	    sort_bodies(localcopy[it], mydestsize);
	    //	printf("sort end %d\n", it);

	    if (it==0)showdt(" 2nd sort");
	    mydestsize = destsize[it];
	    auto offset = deststart[it];
	    for(auto idest=0; idest < mydestsize; idest++){
		bodylocalcopy[it][idest]=b[localcopy[it][idest].key];
	    }
	    if (it==0)    showdt(" Copyforward");
#pragma omp barrier	
	    //	auto mydestsize = destsize[it];
	    offset = deststart[it];
	    for(auto idest=0; idest < mydestsize; idest++){
		//	    printf("it, idest, offst=%d %d %d\n", it, idest, offset);
		b[offset+idest] = bodylocalcopy[it][idest];
	    }
	}
	showdt(" Copyback");
    }
    template<class T>
    void samplesort(T * b,
		    int n)
    {
	samplesort( b, n,
		    [&](T & l) ->auto{return l.getsortkey();});
    }

    // old samplesort_bodies APIs for backward compatibility
    template<class T, class GetKey>
    void samplesort_bodies(T * b,
			   int n,
			   GetKey getkey)
    {
	samplesort(b, n, getkey);
    }
    template<class T>
    void samplesort_bodies(T * b,
		    int n)
    {
	samplesort(b, n);
    }
}

