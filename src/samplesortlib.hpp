#pragma once
#include <unistd.h>
#include <cstring>
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

#include "simdsort.hpp"

#define USELOCALSAMPLESORT 1
#ifdef _OPENMP
#define SAMPLESORTLIB_OMP_GET_MAX_THREADS omp_get_max_threads()
#define SAMPLESORTLIB_OMP_GET_THREAD_NUM omp_get_thread_num()
#else
#define SAMPLESORTLIB_OMP_GET_MAX_THREADS 1
#define SAMPLESORTLIB_OMP_GET_THREAD_NUM 0
#endif
namespace SampleSortLib{

    
    constexpr int simdbufsize=32;    
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

    template <typename T>
    class Has_get_hi_key
    {
	typedef char one;
	typedef long two;
	template <typename C> static one test( decltype(&C::get_hi_key) ) ;
	template <typename C> static two test(...);
    public:
	enum { value = sizeof(test<T>(0)) == sizeof(char) };
    };
    
    template<class T, class GetKey>
    void simd_sort_if_possible(T* b, int n, GetKey getkey,
			       bool indexonly=false)
    {
	constexpr bool  simd_possible = Has_get_hi_key<decltype(getkey(*b))>::value;
	if constexpr (simd_possible)
		     //if constexpr (0)
	    {		
		//fprintf(stderr, "SIMD sort called\n"); 
		uint64_t hi[n+simdbufsize];
		uint64_t lo[n+simdbufsize];
		uint64_t index[n+simdbufsize];
		for (auto i=0;i<n;i++){
		    hi[i]=b[i].value.get_hi_key();
		    lo[i]=b[i].value.get_lo_key();
		    index[i]=b[i].key;
		    //  fprintf(stderr, "%d %lu %lu %lu\n",
		    //	    i, hi[i], lo[i], index[i]);
		    
		}
		//fprintf(stderr, "before sort\n");
		SIMDSortLib::simd_sort(hi, lo, index, n);
		//fprintf(stderr, "after sort\n");
		if (!indexonly){
		    for (auto i=0;i<n;i++){
			b[i].value.set_hi_key(hi[i]);
			b[i].value.set_lo_key(lo[i]);
			b[i].key=index[i];
			//  fprintf(stderr, "%d %lu %lu %lu\n",
			//	    i, hi[i], lo[i], index[i]);
		    }
		}else{
		    for (auto i=0;i<n;i++){
			b[i].key=index[i];
		    }
	        }
	    }else{
	    sort_bodies(b, n);
	}		
    }

#ifdef SVE    

    template<class T>
    inline void copy(T* dest, T* src)
    {
	constexpr int nsize = sizeof(T);
	auto s = (int8_t*) src;
	auto d = (int8_t*) dest;
	int nsve = svcntb();
	for(int i=0; i<nsize; i+= nsve){
	    svbool_t pmask= svwhilelt_b8(i,nsize);
	    svst1(pmask, d+i, svld1(pmask, s+i));
	}
	    
	    
    }
#else
    template<class T>
    inline void copy(T* dest, T* src)
    {
	*dest = *src;
    }
#endif    

    template<class T>
    void settree(T** partition,
		 T* treepartition,
		 int n, 
		 int nlevel,
		 int ilevel,
		 int locincurrentlevel)
    {
	int offset = (1<<ilevel)-1;
	//	fprintf(stderr, "settree called with %d %d %d %d offset=%d\n",
	//	n, nlevel, ilevel, locincurrentlevel, offset);
	int myindex = offset + locincurrentlevel;
	if (ilevel > nlevel || myindex>=n ) {
	    //	    fprintf(stderr, "settree do nothing\n");
	    return;
	}
	settree(partition, treepartition,  n, nlevel,
		ilevel+1,locincurrentlevel*2);
	treepartition[myindex]= **partition;
    //	fprintf(stderr, "settree for index %d = (%d)\n", myindex,
    //		**partition);
	(*partition) ++;
	settree(partition, treepartition,  n, nlevel,
		ilevel+1,locincurrentlevel*2+1);
    }

    template<class T>
    bool is_sorted(T* data,
		   int n)
    {
	for(auto i=1;i<n; i++){
	    if (data[i-1].value > data[i].value){
		return false;
	    }
	}
	return true;
    }
    
    
    template<class T>
    void dump_data(T* data,
		   int n,
		   std::string message)
    {
	fprintf(stderr, "%s:\n", &(message[0]));
	for(auto i=0;i<n; i++){
	    fprintf(stderr, " %d", data[i]);
	}
	fprintf(stderr, "\n");
    }

    template<class Value>
    inline  int find_partition(Value data,
				Value* tree,
				int n,
				int nlevel)
    {
	int ipart = 0;
	int ilevel;
	bool unfilled = false;
	//	fprintf(stderr, "find_partition called with %d %d %d\n", data, n, nlevel);
	for(ilevel=0;ilevel <nlevel; ilevel++){
	    //	    fprintf(stderr, "level=%d part=%d\n", ilevel, ipart);
	    if (ipart > n-1) {
		//		fprintf(stderr, "break at %d %d\n",ilevel,  ipart);
		unfilled = true;
	    }else{
		int inc = data <= tree[ipart]? 0:1;
		ipart = ipart*2+1+inc;
	    }
	}
	int offset = (1<<nlevel)-1;
	if (unfilled){
	    offset = offset -n -1;
	}
	//	fprintf(stderr, "loop out; %d %d offset=%d retval=%d  \n",
	//		ilevel,  ipart, offset,ipart - offset );
	return ipart - offset;
	
    }


    int blocksize = 16000;
    constexpr int nsampleperblock = 100;
    constexpr int minblocks = 4;
    const int randomize_offset = 48271;     // This is a prime number

    void setblocksize(int n){blocksize=n;}

    template<class T>
    void adjust_partition(int n,
			  int nprocessed,
			  int ioverflown,
			  int nblocks,
			  int * size,
			  T** bdest,
			  int * ndestlimit)
    {
	int bufftotalsize=0;
	for(auto i=0;i<nblocks; i++)bufftotalsize += ndestlimit[i];
	auto nremain = bufftotalsize - nprocessed;
	auto nothers = nremain/2/(nblocks - 1);
	auto noverflown = nremain - nothers*(nblocks-1);
	T copybuff[nprocessed];
	auto icopy=0;

	// fprintf(stderr, "adjust partiton nremain=%d nothers=%d nover=%d\n",
	// 	nremain, nothers, noverflown);
	// for(auto i=0;i<nblocks; i++){
	//     fprintf(stderr, "i=%d  %d %d\n", i, size[i], ndestlimit[i]);
	// }
	
	for(auto i=0;i<nblocks; i++){
	    for(auto j=0;j<size[i];j++){
		copybuff[icopy+j] = bdest[i][j];
	    }
	    icopy += size[i];
	}
	int newbdestindex = 0;
	for(auto i=0;i<nblocks;i++){
	    auto ninc = nothers;
	    if (i==ioverflown) ninc = noverflown;
	    ndestlimit[i]=size[i]+ninc;
	    newbdestindex += ndestlimit[i];
	    if (i< nblocks - 1){
		bdest[i+1] = bdest[0] + newbdestindex;
	    }
	}
	icopy=0;
	for(auto i=0;i<nblocks; i++){
	    for(auto j=0;j<size[i];j++){
		bdest[i][j] = copybuff[icopy+j];
	    }
	    icopy += size[i];
	}
	
	
    }	
	
    
    template<class T, class ValueType>
    bool multipartition(T * b,
			 int n,
			 int nblocks,
			ValueType * partition,
			 int * size,
			 T** bdist,
			 int ndestmax)
    {
	auto nparts = nblocks-1;
	if (nblocks < 2){
	    fprintf(stderr, "something broken: n, nb= %d %d\n", n, nblocks);
	}
	int ndestlimit[nblocks];
	for(auto i=0;i<nblocks; i++) ndestlimit[i] = ndestmax;
	// assumption: initially, bdest[i] point to
	//  bdest[0] + i*ndestmax 
	ValueType * pp = partition;
	ValueType treepartition[nparts];
	int nlevel = 0;
	for(int i=nblocks-1;i>0; i>>=1 )nlevel++;
	settree(&pp, treepartition, nparts, nlevel,0,0);
	//dump_data(treepartition, nparts, "tree partitions:");
	for(int i=0;i<nblocks; i++)size[i]=0;
	int lastoffset = nblocks - 1 - (1<<(nlevel-1));
	for(int i=0;i<n;i++){
	    auto  loc = find_partition(b[i].value, treepartition, nparts, nlevel);
	    bdist[loc][size[loc]]=b[i];
	    size[loc]++;
	    if (size[loc]>= ndestlimit[loc]) {
		adjust_partition(n, i+1, loc,  nblocks, size, bdist, ndestlimit);
	    }
	}
	for(int i=0;i<nblocks; i++){
	    if (size[i]> ndestlimit[i]) return false;
	}
	return true;
    }
    
    template<class T>
    bool samplepartition(T * b,
			 int n,
			 int nblocks,
			 int * size,
			 T** bdist,
			 int ndestmax)
    {
	typedef decltype(b->value) ValueType;
	auto nsample = nblocks* nsampleperblock;
	auto nsperblock = nsampleperblock;
	if(nsample > n/2) {
	    nsample = n/2;
	    nsperblock = (nsample + nblocks - 1)/nblocks;
	}
	auto myoffset = randomize_offset ;
	auto nparts = nblocks-1;
	if (myoffset % n == 0 ||n% myoffset == 0){
	    myoffset ++;
	}
	ValueType work[nsample];
	for(auto i=0; i<nsample; i++){
	    auto myi = (i*myoffset)%n;
	    work[i] = b[myi].value;
	}
	// fprintf(stderr, "samples (%d):\n", nsample);
	// for(auto i=0;i<nsample; i++){
	//     fprintf(stderr, " %d", work[i]);
	// }
	// fprintf(stderr, "\n");
	std::sort(work, work+nsample, 
		  [](ValueType l, ValueType r )   ->bool{return l < r;} );
	//dump_data(work, nsample, "samples after sort:");
	ValueType partition[nparts];
	ValueType * pp = partition;
	ValueType treepartition[nparts];
	for (auto i=0; i<nparts; i++){
	    partition[i] = work[(i+1)*nsperblock];
	}
	return  multipartition(b, n, nblocks, partition, size, bdist, ndestmax);
    }
    
    template<class T>
    bool samplepartition_old(T * b,
			 int n,
			 int nblocks,
			 int * size,
			 T** bdist,
			 int ndestmax)
    {
	typedef decltype(b->value) ValueType;
	auto nsample = nblocks* nsampleperblock;
	auto nsperblock = nsampleperblock;
	if(nsample > n/2) {
	    nsample = n/2;
	    nsperblock = (nsample + nblocks - 1)/nblocks;
	}
	auto myoffset = randomize_offset ;
	auto nparts = nblocks-1;
	if (myoffset % n == 0 ||n% myoffset == 0){
	    myoffset ++;
	}
	ValueType work[nsample];
	for(auto i=0; i<nsample; i++){
	    auto myi = (i*myoffset)%n;
	    work[i] = b[myi].value;
	}
	// fprintf(stderr, "samples (%d):\n", nsample);
	// for(auto i=0;i<nsample; i++){
	//     fprintf(stderr, " %d", work[i]);
	// }
	// fprintf(stderr, "\n");
	std::sort(work, work+nsample, 
		  [](ValueType l, ValueType r )   ->bool{return l < r;} );
	//dump_data(work, nsample, "samples after sort:");
	ValueType partition[nparts];
	ValueType * pp = partition;
	ValueType treepartition[nparts];
	for (auto i=0; i<nparts; i++){
	    partition[i] = work[(i+1)*nsperblock];
	}
	//dump_data(partition, nparts, "partitions:");
	int nlevel = 0;
	for(int i=nblocks-1;i>0; i>>=1 )nlevel++;
	settree(&pp, treepartition, nparts, nlevel,0,0);
	//dump_data(treepartition, nparts, "tree partitions:");
	for(int i=0;i<nblocks; i++)size[i]=0;
	int lastoffset = nblocks - 1 - (1<<(nlevel-1));
	for(int i=0;i<n;i++){
	    auto  loc = find_partition(b[i].value, treepartition, nparts, nlevel);
	    bdist[loc][size[loc]]=b[i];
	    size[loc]++;
	    if (size[loc]>= ndestmax) return false;
	}
	for(int i=0;i<nblocks; i++){
	    if (size[i]> ndestmax) return false;

	}
	return true;
    }
    
    
    template<class T>
    bool localsamplesort_try(T * b,
			     int n,
			     int nblocks,
			     int itry)
    {
	if (nblocks < 3){
	    simd_sort_if_possible(b, n, [](T& l)->auto{return l.value;});
	    return true;
	}

	int ndest = n*1.2/nblocks + sqrt(n)*4;
	if (itry) ndest *= (1<<itry);
	if (ndest > n) ndest=n;
	//	fprintf(stderr, "localsamplesort called %d %d  %d %d\n", n,
	//		nblocks, ndest,itry);
	int sizes[nblocks];
	T bdest[nblocks][ndest];
	T * pbdest[nblocks];
	for(auto i=0; i<nblocks; i++) pbdest[i]= bdest[i];
	if (samplepartition(b, n, nblocks, sizes,pbdest, ndest) == false){
	    return false;
	}
	// for(auto i=0;i<nblocks; i++){
	//     printf("block %d:", i);
	//     for(auto j=0;j<sizes[i];j++){
	// 	printf(" (%d %d)", pbdest[i][j].key, pbdest[i][j].value);
	//     }
	//     printf("\n\n");
	// }
	int start = 0;
	for(auto i=0; i<nblocks; i++){
	    //	    std::sort(pbdest[i], pbdest[i]+sizes[i], 
	    //		  [](T l, T r )   ->bool{return l.value < r.value;} );
	    simd_sort_if_possible(pbdest[i], sizes[i], 
				  [](T& l)->auto{return l.value;});

	    // printf("block %d after sort:", i);
	    // for(auto j=0;j<sizes[i];j++){
	    // 	printf(" (%d %d)", pbdest[i][j].key, pbdest[i][j].value);
	    // }
	    // printf("\n\n");
	    for(auto ii=0;ii<sizes[i]; ii++){
		b[start+ii]=pbdest[i][ii];
	    }
	    start += sizes[i];
	}
	return true;
    }

	    

    template<class T>
    bool localsamplesort(T * b,
			 int n,
			 int nblocks)
    {
	for(auto i=0;i<10; i++){
	    if (localsamplesort_try( b, n,nblocks,i)){
		// if (!is_sorted(b, n)){
		//     fprintf(stderr, "ERROR in localsamplesort %d %d\n",
		// 	    n, nblocks);
		// }
		return true;
	    }
	}
	return false;
    }
	

    template<class T, class GetKey>
    bool samplesort_try(T * b,
			int n,
			GetKey getkey,
			int itry)
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

	    return true;
	}
	auto nwork0= (n+nt-1)/nt;
	auto nworkm= n/nt;
	auto nworkp= nworkm+1;
	auto ntplus = n%nt;
	auto nsample= nwork0/nt/10;
	if (nsample < 2) nsample = 2;
	auto nsampletotal = nsample * nt;
	KeyType samplearray[nsampletotal];
	//    printf("nt=%d  nwork=%d nsample=%d \n", nt, nwork0, nsample);
	int isrcstart[nt][nt+1];
	int destsize[nt];
	int deststart[nt];
	auto maxdestsize=0;
	int  destarraysize = ((n*2.5)/nt+sqrt(n)*20+200)*(1<<itry);
	if (destarraysize > n) destarraysize = n;
	KeyValuePair localcopy[nt][destarraysize];
	T bodylocalcopy[nt][destarraysize];

	int ndest = (nwork0/nt*2 + sqrt(nwork0)*4+100)*(1<<itry);
	if (ndest > nwork0) ndest = nwork0;
	//int ndest = nwork0;
	
	KeyValuePair keyparted[nt][nt][ndest];
	int sizeparted[nt][nt];
	auto nparts = nt-1;
	KeyType partition[nparts];
	KeyType * pp = partition;
	KeyType treepartition[nparts];
	bool partition_succeeded[nt];
	bool partition_succeeded_global=true;
	KeyValuePair * bdest[nt][nt];
#pragma omp parallel
	{
	    auto it=SAMPLESORTLIB_OMP_GET_THREAD_NUM;
	    if (it==0) showdt("Enter parallel region");
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
	    if (it==0){
		showdt("Intial sampling");
		//		printf("sample sort size=%d\n", nsampletotal);
		std::sort(samplearray, samplearray+nsampletotal, 
			  [](KeyType  l, KeyType  r )   ->bool{return l < r;} );
		// for(auto i=0;i<nsampletotal; i+=nsample){
		// 	printf("sample[%d]=%g\n", i, samplearray[i]);
		// }
		for (auto itdest=1; itdest<nt; itdest++){
		    partition[itdest-1] = samplearray[itdest*nsample];
		}
		

		showdt(" sample sort");
	    }
	    //	printf("sort size=%d %d\n", myrange,it);
	    for(auto i=0;i<myrange; i++){
		localcopy[it][i].key =i+mystartlocal;
		localcopy[it][i].value =getkey(b[i+mystartlocal]);
	    }
#pragma omp barrier
	    for(auto itdest=0; itdest < nt; itdest++){
		bdest[it][itdest]=keyparted[it][itdest];
	    }
	    partition_succeeded[it]=  multipartition(localcopy[it],
						      myrange, nt, partition,
						      sizeparted[it], bdest[it],
						      ndest);
	    
	    if (it==0)showdt(" partition");
	    destsize[it]=0;
	    for(int itsrc=0;itsrc<nt; itsrc++){
		destsize[it]+= sizeparted[itsrc][it];
	    }
#pragma omp barrier
	    if (it == 0){
		for(int itsrc=0;itsrc<nt; itsrc++){
		    if (!partition_succeeded[itsrc]
			|| destsize[it] > destarraysize){
			partition_succeeded_global =false;
		    }
		}
	    }
	    
	    int idest=0;
	    for(auto itsrc=0; itsrc < nt; itsrc++){
		for (auto isrc=0; isrc < sizeparted[itsrc][it]; isrc++){
		    localcopy[it][idest] = bdest[itsrc][it][isrc];
		    idest++;
		}
	    }
	    auto mydestsize = idest;
	    destsize[it]=mydestsize;
#pragma omp barrier
	    if (it == 0){
		deststart[0]=0;
		for(auto i=1;i<nt;i++) deststart[i] = deststart[i-1]+destsize[i-1];
	    }
#pragma omp barrier
	    
	    //	printf("sort size=%d %d\n", mydestsize,it);
	    //sort_bodies(localcopy[it], mydestsize);
#ifndef USELOCALSAMPLESORT	    
	    simd_sort_if_possible(localcopy[it], mydestsize,
				  [](KeyValuePair& l)
				  ->auto{return l.value;}, true);
#else
	    if (!localsamplesort(localcopy[it], mydestsize,
				 mydestsize/blocksize)){
		fprintf(stderr, "local sample sort failed, exiting...\n");
		exit(-1);
	    }
#endif	    
	    //	printf("sort end %d\n", it);

	    if (it==0)showdt(" 2nd sort");
	    mydestsize = destsize[it];
	    auto offset = deststart[it];
#ifdef SVE
	    svbool_t ptrue =svptrue_b64();
#endif
	    for(auto idest=0; idest < mydestsize; idest++){
		bodylocalcopy[it][idest]=b[localcopy[it][idest].key];
		//copy(&(bodylocalcopy[it][idest]),
		//   b+localcopy[it][idest].key);

#ifdef SVE // improves performance of this part by 60% or more
		double * p = (double*) (b+localcopy[it][idest+16].key);
		svprfd(ptrue, p, SV_PLDL1STRM);
		if constexpr (sizeof(T) > 64){
			svprfd(ptrue, p+8, SV_PLDL1STRM);
		    }

		if constexpr (sizeof(T) > 128){
			svprfd(ptrue, p+16, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 192){
			svprfd(ptrue, p+24, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 256){
			svprfd(ptrue, p+32, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 320){
			svprfd(ptrue, p+40, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 384){
			svprfd(ptrue, p+48, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 448){
			svprfd(ptrue, p+56, SV_PLDL1STRM);
		    }
#endif		
#ifdef AVX512
		double * p = (double*) (b+localcopy[it][idest+16].key);
		_mm_prefetch((char*)p, _MM_HINT_T0);
		if constexpr (sizeof(T) > 64){
			_mm_prefetch((char*)(p+8), _MM_HINT_T0);
		    }

		if constexpr (sizeof(T) > 128){
			_mm_prefetch((char*)(p+16), _MM_HINT_T0);
		    }
		if constexpr (sizeof(T) > 192){
			_mm_prefetch((char*)(p+24), _MM_HINT_T0);
		    }
		if constexpr (sizeof(T) > 256){
			_mm_prefetch((char*)(p+32), _MM_HINT_T0);
		    }
		if constexpr (sizeof(T) > 320){
			_mm_prefetch((char*)(p+40), _MM_HINT_T0);
		    }
		if constexpr (sizeof(T) > 384){
			_mm_prefetch((char*)(p+48), _MM_HINT_T0);
		    }
		if constexpr (sizeof(T) > 448){
			_mm_prefetch((char*)(p+56), _MM_HINT_T0);
		    }
#endif		

	    }
	    if (it==0)    showdt(" Copyforward");
#pragma omp barrier	
	    //	auto mydestsize = destsize[it];
	    offset = deststart[it];
#ifndef SVE
#ifndef AVX512
	    for(auto idest=0; idest < mydestsize; idest++){
		//	    printf("it, idest, offst=%d %d %d\n", it, idest, offset);
		b[offset+idest] = bodylocalcopy[it][idest];
	    }
#else //AVX512
	    int nsize = sizeof(T) * mydestsize;
	    auto d = (int8_t*) (b+offset);
	    auto s = (int8_t*) bodylocalcopy[it];
	    int64_t doffset  = ((int64_t) d) & 63;
	    if (doffset){
		doffset = 64 - doffset;
		memcpy(d, s, doffset);
		d+= doffset;
		s+= doffset;
		nsize-=doffset;
	    }
		
	    
	    int nsve = 64;
	    auto nlimit = (nsize/nsve)*nsve;
	    //fprintf(stderr, "offset, remain = %ld %d\n", doffset, nsize-nlimit);
	    for(int i=0; i<nlimit; i+= nsve){
		_mm512_stream_pd((double*)(d+i), _mm512_loadu_pd((double*)(s+i)));
	    }
	    memcpy(d+nlimit, s+nlimit, nsize-nlimit);

	    
#endif	    
#else //SVE
	    int nsize = sizeof(T) * mydestsize;
	    auto d = (int8_t*) (b+offset);
	    auto s = (int8_t*) bodylocalcopy[it];
	    int nsve = svcntb();
#pragma fj zfill
#pragma fj unroll(8)
#pragma fj loop prefetch_sequential	    
	    for(int i=0; i<nsize; i+= nsve){
		svbool_t pmask= svwhilelt_b8(i,nsize);
		svst1(pmask, d+i, svld1(pmask, s+i));
	    }

	    //	    for(auto idest=0; idest < mydestsize; idest++){
	    //	    printf("it, idest, offst=%d %d %d\n", it, idest, offset);
	    //	b[offset+idest] = bodylocalcopy[it][idest];
	    //	    }
#endif	    
	}
	showdt(" Copyback");
	return true;
    }

    template<class T, class GetKey>
    bool fdps_samplesort_try(T * b,
			int n,
			GetKey getkey,
			int itry)
    {
	typedef decltype(getkey(*b)) KeyType;
	class KeyValuePair{
	public:
	    uint64_t key;
	    KeyType value;
	    inline auto getsortkey() {return value;}    
	};
	//        fprintf(stderr, "fdps_samplesort_try called with %d\n", n);
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

	    return true;
	}
	auto nwork0= (n+nt-1)/nt;
	auto nworkm= n/nt;
	auto nworkp= nworkm+1;
	auto ntplus = n%nt;
	auto nsample= nwork0/nt/10;
	if (nsample < 2) nsample = 2;
	auto nsampletotal = nsample * nt;
	KeyType samplearray[nsampletotal];
	//    printf("nt=%d  nwork=%d nsample=%d \n", nt, nwork0, nsample);
	int isrcstart[nt][nt+1];
	int destsize[nt];
	int deststart[nt];
	auto maxdestsize=0;
	int  destarraysize = ((n*2.5)/nt+sqrt(n)*20+200)*(1<<itry);
	if (destarraysize > n) destarraysize = n;
	KeyValuePair localcopy[nt][destarraysize];
	//	T bodylocalcopy[nt][destarraysize];

	int ndest = (nwork0/nt*2 + sqrt(nwork0)*4+100)*(1<<itry);
	if (ndest > nwork0) ndest = nwork0;
	//int ndest = nwork0;
	
	KeyValuePair keyparted[nt][nt][ndest];
	int sizeparted[nt][nt];
	auto nparts = nt-1;
	KeyType partition[nparts];
	KeyType * pp = partition;
	KeyType treepartition[nparts];
	bool partition_succeeded[nt];
	bool partition_succeeded_global=true;
	KeyValuePair * bdest[nt][nt];
#pragma omp parallel
	{
	    auto it=SAMPLESORTLIB_OMP_GET_THREAD_NUM;
	    if (it==0) showdt("Enter parallel region");
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
	    if (it==0){
		showdt("Intial sampling");
		//		printf("sample sort size=%d\n", nsampletotal);
		std::sort(samplearray, samplearray+nsampletotal, 
			  [](KeyType  l, KeyType  r )   ->bool{return l < r;} );
		// for(auto i=0;i<nsampletotal; i+=nsample){
		// 	printf("sample[%d]=%g\n", i, samplearray[i]);
		// }
		for (auto itdest=1; itdest<nt; itdest++){
		    partition[itdest-1] = samplearray[itdest*nsample];
		}
		

		showdt(" sample sort");
	    }
	    //	printf("sort size=%d %d\n", myrange,it);
	    for(auto i=0;i<myrange; i++){
		localcopy[it][i].key =b[i+mystartlocal].getData();
		localcopy[it][i].value =getkey(b[i+mystartlocal]);
	    }
#pragma omp barrier
	    for(auto itdest=0; itdest < nt; itdest++){
		bdest[it][itdest]=keyparted[it][itdest];
	    }
	    partition_succeeded[it]=  multipartition(localcopy[it],
						      myrange, nt, partition,
						      sizeparted[it], bdest[it],
						      ndest);
	    
	    if (it==0)showdt(" partition");
	    destsize[it]=0;
	    for(int itsrc=0;itsrc<nt; itsrc++){
		destsize[it]+= sizeparted[itsrc][it];
	    }
#pragma omp barrier
	    if (it == 0){
		for(int itsrc=0;itsrc<nt; itsrc++){
		    if (!partition_succeeded[itsrc]
			|| destsize[it] > destarraysize){
			partition_succeeded_global =false;
		    }
		}
	    }
	    
	    int idest=0;
	    for(auto itsrc=0; itsrc < nt; itsrc++){
		for (auto isrc=0; isrc < sizeparted[itsrc][it]; isrc++){
		    localcopy[it][idest] = bdest[itsrc][it][isrc];
		    idest++;
		}
	    }
	    auto mydestsize = idest;
	    destsize[it]=mydestsize;
#pragma omp barrier
	    if (it == 0){
		deststart[0]=0;
		for(auto i=1;i<nt;i++) deststart[i] = deststart[i-1]+destsize[i-1];
	    }
#pragma omp barrier
	    
	    //	printf("sort size=%d %d\n", mydestsize,it);
	    //sort_bodies(localcopy[it], mydestsize);
#ifndef USELOCALSAMPLESORT	    
	    simd_sort_if_possible(localcopy[it], mydestsize,
				  [](KeyValuePair& l)
				  ->auto{return l.value;}, true);
#else
	    if (!localsamplesort(localcopy[it], mydestsize,
				 mydestsize/blocksize)){
		fprintf(stderr, "local sample sort failed, exiting...\n");
		exit(-1);
	    }
#endif	    
	    //	printf("sort end %d\n", it);

	    if (it==0)showdt(" 2nd sort");
	    mydestsize = destsize[it];
	    auto offset = deststart[it];
	    for(auto idest=0; idest < mydestsize; idest++){
		//	    printf("it, idest, offst=%d %d %d\n", it, idest, offset);
		b[offset+idest].setKey(localcopy[it][idest].value);
		b[offset+idest].setData(localcopy[it][idest].key);
	    }
	}
	showdt(" Copyback");
	return true;
    }

    template<class T, class GetKey>
    void samplesort(T * b,
	       int n,
	       GetKey getkey)
    {
	for(auto itry=0;itry<10; itry++){
	    if (samplesort_try(b, n, getkey, itry)){
		return;
	    }
	}
	fprintf(stderr, "sample sort failed.\n");
	exit(-1);
    }
    
    template<class T, class GetKey>
    void fdps_samplesort(T * b,
	       int n,
	       GetKey getkey)
    {
	for(auto itry=0;itry<10; itry++){
	    if (fdps_samplesort_try(b, n, getkey, itry)){
		return;
	    }
	}
	fprintf(stderr, "sample sort failed.\n");
	exit(-1);
    }
    
    template<class T, class GetKey>
    void samplesort_old(T * b,
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
	int  destarraysize = (n*2.5)/nt+sqrt(n)*20+200;
	KeyValuePair localcopy[nt][destarraysize];
	T bodylocalcopy[nt][destarraysize];
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
	    //	printf("sort size=%d %d\n", myrange,it);
	    for(auto i=0;i<myrange; i++){
		key[it][i].key =i+mystartlocal;
		key[it][i].value =getkey(b[i+mystartlocal]);
	    }
	    //	sort_bodies(b+mystartlocal, myrange);
	    if (it==0) showdt("Key made");
#ifndef USELOCALSAMPLESORT	    
	    simd_sort_if_possible(key[it], myrange,
	    			  [](KeyValuePair& l)
				  ->auto{return l.value;});
#else
	    if(! localsamplesort(key[it], myrange, myrange/blocksize)){
		fprintf(stderr, "local sample sort failed, exiting...\n");
		exit(-1);
	    }
#endif	    
	    //    sort_bodies(key[it], myrange);
	    if (it==0) showdt("Initial sort end");
	    //	printf("sort end %d\n", it);
#pragma omp barrier	
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
		if (destarraysize < maxdestsize){
		    fprintf(stderr, "n, max, size = %d %d %d\n",
			    n, maxdestsize, destarraysize);
		    exit(-1);
		}
		
	    }
#pragma omp barrier	

	    mydestsize = destsize[it];
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
	    //sort_bodies(localcopy[it], mydestsize);
#ifndef USELOCALSAMPLESORT	    
	    simd_sort_if_possible(localcopy[it], mydestsize,
				   [](KeyValuePair& l)
				  ->auto{return l.value;}, true);
#else
	    if (!localsamplesort(localcopy[it], mydestsize,
				 mydestsize/blocksize)){
		fprintf(stderr, "local sample sort failed, exiting...\n");
		exit(-1);
	    }
#endif	    
	    //	printf("sort end %d\n", it);

	    if (it==0)showdt(" 2nd sort");
	    mydestsize = destsize[it];
	    auto offset = deststart[it];
#ifdef SVE
		svbool_t ptrue =svptrue_b64();
#endif
	    for(auto idest=0; idest < mydestsize; idest++){
		bodylocalcopy[it][idest]=b[localcopy[it][idest].key];
		//copy(&(bodylocalcopy[it][idest]),
		//   b+localcopy[it][idest].key);

#ifdef SVE // improves performance of this part by 60% or more
		double * p = (double*) (b+localcopy[it][idest+16].key);
		svprfd(ptrue, p, SV_PLDL1STRM);
		if constexpr (sizeof(T) > 64){
			svprfd(ptrue, p+8, SV_PLDL1STRM);
		    }

		if constexpr (sizeof(T) > 128){
			svprfd(ptrue, p+16, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 192){
			svprfd(ptrue, p+24, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 256){
			svprfd(ptrue, p+32, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 320){
			svprfd(ptrue, p+40, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 384){
			svprfd(ptrue, p+48, SV_PLDL1STRM);
		    }
		if constexpr (sizeof(T) > 448){
			svprfd(ptrue, p+56, SV_PLDL1STRM);
		    }
#endif		

	    }
	    if (it==0)    showdt(" Copyforward");
#pragma omp barrier	
	    //	auto mydestsize = destsize[it];
	    offset = deststart[it];
#ifndef SVE	    
	    for(auto idest=0; idest < mydestsize; idest++){
		//	    printf("it, idest, offst=%d %d %d\n", it, idest, offset);
		b[offset+idest] = bodylocalcopy[it][idest];
	    }
#else
	    int nsize = sizeof(T) * mydestsize;
	    auto d = (int8_t*) (b+offset);
	    auto s = (int8_t*) bodylocalcopy[it];
	    int nsve = svcntb();
#pragma fj zfill
#pragma fj unroll(8)
#pragma fj loop prefetch_sequential	    
	    for(int i=0; i<nsize; i+= nsve){
		svbool_t pmask= svwhilelt_b8(i,nsize);
		svst1(pmask, d+i, svld1(pmask, s+i));
	    }

	    //	    for(auto idest=0; idest < mydestsize; idest++){
		//	    printf("it, idest, offst=%d %d %d\n", it, idest, offset);
	    //	b[offset+idest] = bodylocalcopy[it][idest];
	    //	    }
#endif	    
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

