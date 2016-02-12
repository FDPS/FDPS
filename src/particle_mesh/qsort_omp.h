#include "treepm_header.h"



namespace ParticleSimulator{
    namespace ParticleMesh{

inline void particleSortUsingKey( long long int *key, 
				  pParticle particle, 
				  const int n){

  int i = 0;
  int j = n-1;
  long long int  pivot = key[n/2];

  while(1){
    while( key[i] < pivot){
      i++;
    }
    while( key[j] > pivot){
      j--;
    }
    if( i >= j){
       break;
    }

    long long int  tmp_key = key[j];
    key[j] = key[i];
    key[i] = tmp_key;
    Particle tmp_p = particle[j];
    particle[j] = particle[i];
    particle[i] = tmp_p;
    i++;
    j--;

  }

  if( 0 < i-1){
    particleSortUsingKey( key, particle, i);
  }
  if( j+1 < n-1){
    particleSortUsingKey( key+j+1, particle+j+1, n-j-1);
  }

}



inline void qsortKey( int *index, 
		      long long int *key, 
		      const int n){

  int i = 0;
  int j = n-1;
  long long int  pivot = key[n/2];

  while(1){
    while( key[i] < pivot){
      i++;
    }
    while( key[j] > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    long long int  tmp_key = key[j];
    key[j] = key[i];
    key[i] = tmp_key;
    int tmp_index = index[j];
    index[j] = index[i];
    index[i] = tmp_index;
    i++;
    j--;

  }

  if( 0 < i-1){
    qsortKey( index, key, i);
  }
  if( j+1 < n-1){
    qsortKey( index+j+1, key+j+1, n-j-1);
  }

}



inline void qsortKeyNonParallel( int *index, long long int *key, const int n, 
				 int level, const int max_level, 
				 int **start_index, long long int **start_key,
				 int *n_block, int *nsort){

  if( level == max_level){
    start_index[*nsort] = index;
    start_key[*nsort] = key;
    n_block[*nsort] = n;
    (*nsort) ++;
    return;
  }

  int i = 0;
  int j = n-1;
  long long int pivot = key[n/2];

  while(1){
    while( key[i] < pivot){
      i++;
    }
    while( key[j] > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    long long int  tmp_key = key[j];
    key[j] = key[i];
    key[i] = tmp_key;
    int tmp_index = index[j];
    index[j] = index[i];
    index[i] = tmp_index;
    i++;
    j--;

  }

  if( 0 < i-1){
    level ++;
    qsortKeyNonParallel( index, key, i, level, max_level,
			 start_index, start_key, n_block, nsort);
    level --;
  }
  if( j+1 < n-1){
    level ++;
    qsortKeyNonParallel( index+j+1, key+j+1, n-j-1, level, max_level,
			 start_index, start_key, n_block, nsort);
    level --;
  }

}



inline void qsortKey2( int *index, long long int *key, const int n){

  static int *start_index[NSORTMAX];
  static long long int *start_key[NSORTMAX];
  static int n_block[NSORTMAX];

  int nsort = 0;
  double nowtime = 0.0;
  getTime(&nowtime);
  qsortKeyNonParallel( index, key, n, 0, NUMBER_OF_LEVEL_QSORT_SERIAL,
		       start_index, start_key, n_block, &nsort);
  fprintf( stderr, "%lf\n", getTime(&nowtime));

#pragma omp parallel for CHUNK_QSORT
  for( int i=0; i<nsort; i++){
    qsortKey( start_index[i], start_key[i], n_block[i]);
  }

  fprintf( stderr, "%lf\n", getTime(&nowtime));

}



inline void qsortKeyWrapper( int *index, long long int *key, const int n){

#ifndef _OPENMP
  qsortKey( index, key, n);
#else
  qsortKey2( index, key, n);
#endif

}
///////////////////////////////////////////////////////////////////////////////


template<typename type_obj1, typename type_obj2, typename type_pivot>
inline void qsortPivotNonParallel( type_obj1 *obj1, 
				   type_obj2 *obj2, 
				   type_pivot *obj_pivot, 
				   const int n, 
				   int level, const int max_level, 
				   type_obj1 **start_obj1, 
				   type_obj2 **start_obj2, 
				   type_pivot **start_obj_pivot,
				   int *n_block, int *nsort){

  if( level == max_level){
    start_obj1[*nsort] = obj1;
    start_obj2[*nsort] = obj2;
    start_obj_pivot[*nsort] = obj_pivot;
    n_block[*nsort] = n;
    (*nsort) ++;
    return;
  }

  int i = 0;
  int j = n-1;
  type_pivot pivot = obj_pivot[n/2];

  while(1){
    while( obj_pivot[i] < pivot){
      i++;
    }
    while( obj_pivot[j] > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    type_pivot tmp_obj_pivot = obj_pivot[j];
    obj_pivot[j] = obj_pivot[i];
    obj_pivot[i] = tmp_obj_pivot;
    type_obj1 tmp_obj1 = obj1[j];
    obj1[j] = obj1[i];
    obj1[i] = tmp_obj1;
    type_obj2 tmp_obj2 = obj2[j];
    obj2[j] = obj2[i];
    obj2[i] = tmp_obj2;
    i++;
    j--;

  }

  if( 0 < i-1){
    level ++;
    qsortPivotNonParallel( obj1, obj2, obj_pivot, i, level, max_level,
			   start_obj1, start_obj2, start_obj_pivot, 
			   n_block, nsort);
    level --;
  }
  if( j+1 < n-1){
    level ++;
    qsortPivotNonParallel( obj1+j+1, obj2+j+1, obj_pivot+j+1, n-j-1, level, max_level,
			   start_obj1, start_obj2, start_obj_pivot, n_block, nsort);
    level --;
  }

}




template< typename type_obj1, typename type_obj2, typename type_pivot>
inline void qsortPivot( type_obj1 *obj1, 
			type_obj2 *obj2, 
			type_pivot *obj_pivot, 
			const int n){

  int i = 0;
  int j = n-1;
  type_pivot pivot = obj_pivot[n/2];

  while(1){
    while( obj_pivot[i] < pivot){
      i++;
    }
    while( obj_pivot[j] > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    type_pivot tmp_obj_pivot = obj_pivot[j];
    obj_pivot[j] = obj_pivot[i];
    obj_pivot[i] = tmp_obj_pivot;
    type_obj1 tmp_obj1 = obj1[j];
    obj1[j] = obj1[i];
    obj1[i] = tmp_obj1;
    type_obj2 tmp_obj2 = obj2[j];
    obj2[j] = obj2[i];
    obj2[i] = tmp_obj2;
    i++;
    j--;

  }

  if( 0 < i-1){
    qsortPivot( obj1, obj2, obj_pivot, i);
  }
  if( j+1 < n-1){
    qsortPivot( obj1+j+1, obj2+j+1, obj_pivot+j+1, n-j-1);
  }

}



template< typename type_obj1, typename type_obj2, typename type_pivot>
inline void qsortPivotOmp( type_obj1 *obj1, 
			   type_obj2 *obj2, 
			   type_pivot *obj_pivot, 
			   const int n,
			   const int level_serial, 
			   const int chunk_qsort){

  static type_obj1 *start_obj1[NSORTMAX];
  static type_obj2 *start_obj2[NSORTMAX];
  static type_pivot *start_obj_pivot[NSORTMAX];
  static int n_block[NSORTMAX];

  int nsort = 0;
  qsortPivotNonParallel( obj1, obj2, obj_pivot, n, 0, level_serial,
			 start_obj1, start_obj2, start_obj_pivot, n_block, &nsort);

#ifdef OMP_SCHDULE_DISABLE
#pragma omp parallel for
#else
#pragma omp parallel for schedule( dynamic, chunk_qsort)
#endif
  for( int i=0; i<nsort; i++){
    qsortPivot( start_obj1[i], start_obj2[i], start_obj_pivot[i], n_block[i]);
  }

}



template< typename type_obj1, typename type_obj2, typename type_pivot>
inline void qsortPivotWrapper( type_obj1 *obj1, 
			       type_obj2 *obj2, 
			       type_pivot *obj_pivot, 
			       const int n,
			       const int level_serial, 
			       const int chunk_qsort){

  if( n==0)  return;

#ifndef _OPENMP
  qsortPivot( obj1, obj2, obj_pivot, n);
#else
  qsortPivotOmp( obj1, obj2, obj_pivot, n, level_serial, chunk_qsort);
#endif

}




///////////////////////////////////////////////////////////////////////////////
template<typename type_obj, typename type_pivot>
inline void qsortPivotNonParallel( type_obj *obj, 
				   type_pivot *obj_pivot, 
				   const int n, 
				   int level, 
				   const int max_level, 
				   type_obj **start_obj, 
				   type_pivot **start_obj_pivot,
				   int *n_block, 
				   int *nsort){

#if 0
  if( n < 35000){
#else
  if( level == max_level){
#endif
    start_obj[*nsort] = obj;
    start_obj_pivot[*nsort] = obj_pivot;
    n_block[*nsort] = n;
    (*nsort) ++;
    return;
  }

  int i = 0;
  int j = n-1;
  type_pivot pivot = obj_pivot[n/2];

  while(1){
    while( obj_pivot[i] < pivot){
      i++;
    }
    while( obj_pivot[j] > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    type_pivot tmp_obj_pivot = obj_pivot[j];
    obj_pivot[j] = obj_pivot[i];
    obj_pivot[i] = tmp_obj_pivot;
    type_obj tmp_obj = obj[j];
    obj[j] = obj[i];
    obj[i] = tmp_obj;
    i++;
    j--;

  }

  if( 0 < i-1){
    level ++;
    qsortPivotNonParallel( obj, obj_pivot, i, level, max_level,
			   start_obj, start_obj_pivot, n_block, nsort);
    level --;
  }
  if( j+1 < n-1){
    level ++;
    qsortPivotNonParallel( obj+j+1, obj_pivot+j+1, n-j-1, level, max_level,
			   start_obj, start_obj_pivot, n_block, nsort);
    level --;
  }

}




template<typename type_obj, typename type_pivot>
inline void qsortPivot( type_obj *obj, 
			type_pivot *obj_pivot, 
			const int n){

  int i = 0;
  int j = n-1;
  type_pivot pivot = obj_pivot[n/2];

  while(1){
    while( obj_pivot[i] < pivot){
      i++;
    }
    while( obj_pivot[j] > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    type_pivot tmp_obj_pivot = obj_pivot[j];
    obj_pivot[j] = obj_pivot[i];
    obj_pivot[i] = tmp_obj_pivot;
    type_obj tmp_obj = obj[j];
    obj[j] = obj[i];
    obj[i] = tmp_obj;
    i++;
    j--;

  }

  if( 0 < i-1){
    qsortPivot( obj, obj_pivot, i);
  }
  if( j+1 < n-1){
    qsortPivot( obj+j+1, obj_pivot+j+1, n-j-1);
  }

}



template<typename type_obj, typename type_pivot>
inline void qsortPivotOmp( type_obj *obj, 
			   type_pivot *obj_pivot, 
			   const int n,
			   const int level_serial, 
			   const int chunk_qsort){

  static type_obj *start_obj[NSORTMAX];
  static type_pivot *start_obj_pivot[NSORTMAX];
  static int n_block[NSORTMAX];

  int nsort = 0;
  //  double nowtime = 0.0;
  //  getTime(&nowtime);
  qsortPivotNonParallel( obj, obj_pivot, n, 0, level_serial,
			 start_obj, start_obj_pivot, n_block, &nsort);
  //  fprintf( stderr, "%lf\n", getTime(&nowtime));

#ifdef OMP_SCHDULE_DISABLE
#pragma omp parallel for
#else
#pragma omp parallel for schedule( dynamic, chunk_qsort)
#endif
  for( int i=0; i<nsort; i++){
    qsortPivot( start_obj[i], start_obj_pivot[i], n_block[i]);
  }

  //  fprintf( stderr, "%lf\n", getTime(&nowtime));

}



template<typename type_obj, typename type_pivot>
inline void qsortPivotWrapper( type_obj *obj, 
			       type_pivot *obj_pivot, 
			       const int n,
			       const int level_serial, 
			       const int chunk_qsort){

  if( n==0)  return;

#ifndef _OPENMP
  qsortPivot( obj, obj_pivot, n);
#else
  qsortPivotOmp( obj, obj_pivot, n, level_serial, chunk_qsort);
#endif

}





///////////////////////////////////////////////////////////////////////////////

template<typename type_obj>
inline void qsortBasicNonParallel( type_obj *obj, 
				   const int n, 
				   int level, 
				   const int max_level, 
				   type_obj **start_obj, 
				   int *n_block, 
				   int *nsort){

  if( level == max_level){
    start_obj[*nsort] = obj;
    n_block[*nsort] = n;
    (*nsort) ++;
    return;
  }

  int i = 0;
  int j = n-1;
  type_obj pivot = obj[n/2];

  while(1){
    while( obj[i] < pivot){
      i++;
    }
    while( obj[j] > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    type_obj tmp_obj = obj[j];
    obj[j] = obj[i];
    obj[i] = tmp_obj;
    i++;
    j--;

  }

  if( 0 < i-1){
    level ++;
    qsortBasicNonParallel( obj, i, level, max_level,
			   start_obj, n_block, nsort);
    level --;
  }
  if( j+1 < n-1){
    level ++;
    qsortBasicNonParallel( obj+j+1, n-j-1, level, max_level,
			   start_obj, n_block, nsort);
    level --;
  }

}



template<typename type_obj>
inline void qsortBasic( type_obj *obj, 
			const int n){

  int i = 0;
  int j = n-1;
  type_obj pivot = obj[n/2];

  while(1){
    while( obj[i] < pivot){
      i++;
    }
    while( obj[j] > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    type_obj tmp_obj = obj[j];
    obj[j] = obj[i];
    obj[i] = tmp_obj;
    i++;
    j--;

  }

  if( 0 < i-1){
    qsortBasic( obj, i);
  }
  if( j+1 < n-1){
    qsortBasic( obj+j+1, n-j-1);
  }

}



template<typename type_obj>
inline void qsortBasicOmp( type_obj *obj, 
			   const int n,
			   const int level_serial, 
			   const int chunk_qsort){

  static type_obj *start_obj[NSORTMAX];
  static int n_block[NSORTMAX];

  int nsort = 0;
  qsortBasicNonParallel( obj, n, 0, level_serial,
			 start_obj, n_block, &nsort);

#ifdef OMP_SCHDULE_DISABLE
#pragma omp parallel for
#else
#pragma omp parallel for schedule( dynamic, chunk_qsort)
#endif
  for( int i=0; i<nsort; i++){
    qsortBasic( start_obj[i], n_block[i]);
  }

}



template<typename type_obj>
inline void qsortBasicWrapper( type_obj *obj, 
			       const int n,
			       const int level_serial, 
			       const int chunk_qsort){

  if( n==0)  return;

#ifndef _OPENMP
  qsortBasic( obj, n);
#else
  qsortBasicOmp( obj, n, level_serial, chunk_qsort);
#endif

}

///////////////////////////////////////////////////////////////////////////////
inline void particleSortUsingXposNonParallel( pParticle p, 
					      const int n, 
					      int level, 
					      const int max_level, 
					      Particle **start_p, 
					      int *n_block, 
					      int *nsort){

  if( level == max_level){
    start_p[*nsort] = p;
    n_block[*nsort] = n;
    (*nsort) ++;
    return;
  }

  int i = 0;
  int j = n-1;
  double pivot = p[n/2].xpos;

  while(1){
    while( p[i].xpos < pivot){
      i++;
    }
    while( p[j].xpos > pivot){
      j--;
    }
    if( i >= j){
      break;
    }

    Particle tmp_p = p[j];
    p[j] = p[i];
    p[i] = tmp_p;
    i++;
    j--;

  }

  if( 0 < i-1){
    level ++;
    particleSortUsingXposNonParallel( p, i, level, max_level,
				      start_p, n_block, nsort);
    level --;
  }
  if( j+1 < n-1){
    level ++;
    particleSortUsingXposNonParallel( p+j+1, n-j-1, level, max_level,
				      start_p, n_block, nsort);
    level --;
  }

}


inline void particleSortUsingXposOmp( pParticle p,
				      const int n,
				      const int level_serial, 
				      const int chunk_qsort){

  static pParticle start_p[NSORTMAX];
  static int n_block[NSORTMAX];

  int nsort = 0;
  particleSortUsingXposNonParallel( p, n, 0, level_serial,
				    start_p, n_block, &nsort);

#ifdef OMP_SCHDULE_DISABLE
#pragma omp parallel for
#else
#pragma omp parallel for schedule( dynamic, chunk_qsort)
#endif
  for( int i=0; i<nsort; i++){
    particleSortUsingXpos( start_p[i], n_block[i]);
  }

}


inline void qsortParticleUsingXpos( pParticle p,
				    const int n,
				    const int level_serial, 
				    const int chunk_qsort){

  if( n==0)  return;

#ifndef _OPENMP
  particleSortUsingXpos( p, n);
#else
  particleSortUsingXposOmp( p, n, level_serial, chunk_qsort);
#endif

}



    } // namespace ParticleMesh
}     // namespace ParticleSimulator
