#include "pp.h"

namespace ParticleSimulator{
    namespace ParticleMesh{


static char cbuf[CHARMAX];

double (*p_cache)[4];


void dumpAllParticle( const pParticle particle, const int n,
                      const int inode, const char *filename){

  sprintf( cbuf, "%s-%d", filename, inode);
  FILE *fout = fopen( cbuf, "w");
  dumpAllParticle( particle, n, fout);
  fclose( fout);

}



void dumpAllParticle( const pParticle particle, const int n, FILE *fout){

  for( int i=0; i<n; i++){
    fprintf( fout, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%lld\n",
             i, 
             particle[i].xpos, particle[i].ypos, particle[i].zpos,
             particle[i].xvel, particle[i].yvel, particle[i].zvel,
             particle[i].id);
  }

}



void dumpParticle( const Particle *particle, FILE *outstream){

  fprintf( outstream, "%e\t%e\t%e\t%e\t%e\t%e\t%lld\n",
	   particle->xpos, particle->ypos, particle->zpos,
	   particle->xvel, particle->yvel, particle->zvel,
	   particle->id);

}



void dumpTree( const TreeTmp *treetmp, FILE *outstream){

  for( int i=0; i<treetmp->clistmax; i++){
    //    if( treetmp->tree_cell[i].n == 0)  continue;
    fprintf( outstream, "%lld\t%d\t%d\t%d\t%e\t%e\t%e\t%e\n",
#ifdef TREE2
	     (long long int)i,
#else
	     treetmp->tree_cell[i].key,
#endif
	     treetmp->tree_cell[i].first_index,
	     treetmp->tree_cell[i].n,
	     treetmp->tree_cell[i].next,
	     treetmp->tree_cell[i].cm[0],
	     treetmp->tree_cell[i].cm[1],
	     treetmp->tree_cell[i].cm[2],
	     treetmp->tree_cell[i].m);
  }

}



/* constructer : struct TreeCell
 * initialize member variable
 */
void TreeCellTreeCell(pTreeCell tree_cell, const int n){

  int i;
  for( i=0; i<n; i++){
#ifndef TREE2
    tree_cell[i].key = 0;
    tree_cell[i].zeroflag = 0;
#endif
    tree_cell[i].first_index = 0;
    tree_cell[i].n = 0;
    tree_cell[i].next = -1;
    tree_cell[i].m = 0.0;
    tree_cell[i].cm[0] = 0.0;
    tree_cell[i].cm[1] = 0.0;
    tree_cell[i].cm[2] = 0.0;
  }

}



/* constructer : struct TreeTmp
   initialize member variable
   allocate tree cell
*/

void TreeTmpTreeTmp( pTreeTmp tree_tmp, const pTreeParam treeparam,
		     const int n){

#ifdef TREECONSTRUCTION_PARALLEL
  int clistmax = n;
#else
  int clistmax = (int)(250000.0 * (30.0/treeparam->tree_nleaf) * (n/1.0e6));
#endif
  tree_tmp->clistmax = clistmax;
  //fprintf( stderr, "#clistmax:%d\t", clistmax);

  int tree_clistmask = (int)(log((double)clistmax)/log(2.0));
  //fprintf( stderr, "bit-mask-h:%d\t", tree_clistmask);
  tree_clistmask = (1<<(tree_clistmask-1))-1;
  tree_tmp->tree_clistmask = tree_clistmask;
  //fprintf( stderr, "bit-mask:%d\n", tree_clistmask);

  tree_tmp->tree_cell = (pTreeCell) my_malloc( sizeof(TreeCell) * clistmax);
  TreeCellTreeCell(tree_tmp->tree_cell, clistmax);

  int walklist_size = (int)((double)n/(double)treeparam->tree_ncrit*10.0);
  tree_tmp->walklist_size = walklist_size;

#ifdef TREE2
  tree_tmp->walklist = (pICell) my_malloc( sizeof(ICell) * tree_tmp->walklist_size);
#else
  tree_tmp->walklist = (int *) my_malloc( sizeof(int) * tree_tmp->walklist_size);
#endif

  tree_tmp->key = (long long int *) my_malloc( sizeof(long long int) * n);
  tree_tmp->index = (int *) my_malloc( sizeof(int) * n);

  for( int i=0; i<n; i++){
    tree_tmp->index[i] = i;
  }

}


void TreeTmpTreeTmp( pTreeTmp tree_tmp, const int n, const int offset){

#ifndef TREE2
  pTreeCell tree_cell = tree_tmp->tree_cell;
  for( int i=0; i<tree_tmp->clistmax; i++){
    tree_cell[i].next = -1;
  }
#endif

  int *index = tree_tmp->index;;
  for( int i=offset; i<(n+offset); i++){
    index[i] = i;
  }
  
}


/* destructer : struct TreeTmp */
void TreeTmpDTreeTmp( pTreeTmp tree_tmp){

  free( tree_tmp->tree_cell);
  free( tree_tmp->walklist);
  free( tree_tmp->key);
  free( tree_tmp->index);

}



void TreeTmpPrintMemsize( pTreeTmp tree_tmp, const int n){

  int clistmax = tree_tmp->clistmax;
  int walklist_size = tree_tmp->walklist_size;

  double memsize = sizeof(TreeCell) * clistmax;
#ifdef TREE2
  memsize += sizeof(ICell) * walklist_size;
#else
  memsize += sizeof(int) * walklist_size;
#endif
  memsize += sizeof(long long int) * n;
  memsize += sizeof(int) * n;

  fprintf(stderr, "#MEMSIZE Tree Structure : %4.0lfMbyte %ld\n", 
	  memsize/1.0e6, sizeof(TreeCell));

}



void initTreeParam( pTreeParam treeparam, const double theta,
		   const int ncrit, const int nleaf, const double uniform_mass){

  treeparam->tree_theta = theta;
  treeparam->tree_theta2 = theta*theta;
  treeparam->tree_theta2_quarter = theta*theta*0.25;
  treeparam->tree_ncrit = ncrit;
  treeparam->tree_nleaf = nleaf;
  treeparam->uniform_mass = uniform_mass;

}



#ifndef TREE2
int keyToAdr( const long long int key, const pTreeCell tree_cell, 
	      const int tree_clistmask){

  int adr = tree_clistmask & key;

  if( tree_cell[adr].key != key){
    do{
      adr = tree_cell[adr].next;
    }while(tree_cell[adr].key != key);
  }

  return adr;

}
#endif


#ifndef TREE2
void getCellRange2( const pTreeCell tree_cell, double *ccell, double *length){

  int key_level = ((0x0000ff00) & tree_cell->zeroflag) >> 8;
  //int key_level = tree_cell->zeroflag >> 8;
  int level = (63 - key_level) / 3;
  double fac = pow(2.0, (double)level);
  length[0] = 0.5 / fac;
  length[1] = 0.5 / fac;
  length[2] = 0.5 / fac;
  double tmpscale = 1.0 / pow(2.0, 21.0);
  long long int tmpkey[3] = {0,0,0};
  long long int key2 = tree_cell->key << key_level;

  int ii;
  for( ii=0; ii<21; ii++){
    tmpkey[2] |= ((key2>>(3*ii))&1)<<ii;
    tmpkey[1] |= ((key2>>(3*ii+1))&1)<<ii;
    tmpkey[0] |= ((key2>>(3*ii+2))&1)<<ii;
  }

  ccell[0] = tmpkey[0]*tmpscale + length[0];
  ccell[1] = tmpkey[1]*tmpscale + length[1];
  ccell[2] = tmpkey[2]*tmpscale + length[2];

}
#endif



void getCellRange2( const pTreeCell tree_cell, double *ccell, double *length,
		    const int key_level, const long long int key, const double maxx){

  int level = (63 - key_level) / 3;
  double fac = pow(2.0, (double)level);
  length[0] = 0.5 / fac;
  length[1] = 0.5 / fac;
  length[2] = 0.5 / fac;
  double tmpscale = 1.0 / pow(2.0, 21.0);
  long long int tmpkey[3] = {0,0,0};
  long long int key2 = key << key_level;

  int ii;
  for( ii=0; ii<21; ii++){
    tmpkey[2] |= ((key2>>(3*ii))&1)<<ii;
    tmpkey[1] |= ((key2>>(3*ii+1))&1)<<ii;
    tmpkey[0] |= ((key2>>(3*ii+2))&1)<<ii;
  }
  ccell[0] = tmpkey[0]*tmpscale + length[0];
  ccell[1] = tmpkey[1]*tmpscale + length[1];
  ccell[2] = tmpkey[2]*tmpscale + length[2];

}



#ifndef TREE2
void getCellRange3( const pTreeCell tree_cell, double *ccell, double *length,
		    const pTreeTmp treetmp){

  int key_level = ((0x0000ff00) & tree_cell->zeroflag) >> 8;
  int level = (63 - key_level) / 3;
  double fac = pow(2.0, (double)level);

  length[0] = 0.5 * treetmp->maxx / fac;
  length[1] = length[0];
  length[2] = length[0];
  double tmpscale = treetmp->maxx  / (double)(1<<21);

  long long int tmpkey[3] = {0,0,0};
  long long int key2 = tree_cell->key << key_level;

  int ii;
  for( ii=0; ii<21; ii++){
    tmpkey[2] |= ((key2>>(3*ii))&1)<<ii;
    tmpkey[1] |= ((key2>>(3*ii+1))&1)<<ii;
    tmpkey[0] |= ((key2>>(3*ii+2))&1)<<ii;
  }

  ccell[0] = tmpkey[0]*tmpscale + length[0] + treetmp->bmin[0];
  ccell[1] = tmpkey[1]*tmpscale + length[1] + treetmp->bmin[1];
  ccell[2] = tmpkey[2]*tmpscale + length[2] + treetmp->bmin[2];

}
#endif



#ifndef TREE2
void makeInteractionList(const int ikey, const double cell_length2,
			 const double *cicell, const double *i_half_length,
			 const pTreeCell tree_cell,
			 pJList jlist, int *njlist,
			 const long long int current_key, const int tree_clistmask, 
			 const pParticle particle, const int *index, 
			 const pTreeParam treeparam){


  const int adr = keyToAdr(current_key, tree_cell, tree_clistmask);

  // calc theta
  double dr = 0.0;
  double peri[3] = {0.0, 0.0, 0.0};
  for( int j=0; j<3; j++){
    double drj = cicell[j] - tree_cell[adr].cm[j];
    if( drj > 0.5)   peri[j] =  1.0;
    if( drj < -0.5)  peri[j] = -1.0;
    drj = fabs(drj);
    if( drj > 0.5)  drj = 1.0 - drj;
    drj -= i_half_length[j];
    if( drj > 0)  dr += drj * drj;
  }
  const double theta2_dr2 = treeparam->tree_theta2 * dr;

  if(theta2_dr2 > cell_length2){   // theta > d/l
    if( dr < RADIUS_FOR_PP2){
#ifdef KCOMPUTER
      jlist->x[*njlist][3] = tree_cell[adr].m;
#else
      jlist->m[*njlist] = tree_cell[adr].m;
#endif
      jlist->x[*njlist][0] = tree_cell[adr].cm[0];
      jlist->x[*njlist][1] = tree_cell[adr].cm[1];
      jlist->x[*njlist][2] = tree_cell[adr].cm[2];
      jlist->x[*njlist][0] += peri[0];
      jlist->x[*njlist][1] += peri[1];
      jlist->x[*njlist][2] += peri[2];
      (*njlist) ++;
    }
  }
  else{
    if( tree_cell[adr].n > treeparam->tree_nleaf){  //have child
      for( int i=0; i<8; i++){
	const long long int tmpkey = ((current_key<<3) | i);
	if( tmpkey != tree_cell[ikey].key){  //not own cell
	  const int zeroflag = ( tree_cell[adr].zeroflag>>i) & 1;
	  if( zeroflag == 1)  continue;
	  makeInteractionList(ikey, cell_length2*0.25, cicell, i_half_length,
			      tree_cell, jlist, njlist,
			      tmpkey, tree_clistmask, particle, index, treeparam);
	}
	else{ // own cell -> i-particle copy to j-particle
	  const int c = keyToAdr(tmpkey, tree_cell, tree_clistmask);
	  const int start = tree_cell[c].first_index;
	  const int n = tree_cell[c].n;
	  const int nj = *njlist;
	  for( int j=0; j<n; j++){
	    int ji = nj + j;
	    int pi = start + j;
#ifdef TREE_PARTICLE_CACHE
	    jlist->x[ji][0] = p_cache[pi][0];
	    jlist->x[ji][1] = p_cache[pi][1];
	    jlist->x[ji][2] = p_cache[pi][2];
	    double mass = p_cache[pi][3];
#else
	    int ii = index[pi];
	    getPos( &particle[ii], jlist->x[ji]);
	    double mass = getMass( &particle[ii], treeparam->uniform_mass);
#endif
#ifdef KCOMPUTER
	    jlist->x[ji][3] = mass;
#else
	    jlist->m[ji] = mass;
#endif
	  }
	  (*njlist) += n;
	}
      }
    }
    else{  // lowest cell
      const int start = tree_cell[adr].first_index;
      const int n = tree_cell[adr].n;
      const int nj = *njlist;
      for( int j=0; j<n; j++){
	int ji = nj + j;
	int pi = start + j;
#ifdef TREE_PARTICLE_CACHE
	jlist->x[ji][0] = p_cache[pi][0];
	jlist->x[ji][1] = p_cache[pi][1];
	jlist->x[ji][2] = p_cache[pi][2];
	double mass = p_cache[pi][3];
#else
	int ii = index[pi];
	getPos( &particle[ii], jlist->x[ji]);
	double mass = getMass( &particle[ii], treeparam->uniform_mass);
#endif
#ifdef KCOMPUTER
	jlist->x[ji][3] = mass;
#else
	jlist->m[ji] = mass;
#endif
	jlist->x[ji][0] += peri[0];
	jlist->x[ji][1] += peri[1];
	jlist->x[ji][2] += peri[2];
      }
      (*njlist) += n;
    }
  }


}
#endif



#define _out_
void scan_keys( const long long int key[],
		const int    len,
		const int    rshift,
		_out_ int    off[],
		_out_ int    num[]){
  off[0] = 0;
  for(int i=0, ic=0; ic<8; i++){
    int c = (key[i] >> rshift) & 7;
    if(i == len) c = 8;
    // printf("%lx %d\n", key[i], c);
    while(ic != c){
      off[++ic] = i;
    }
  }
  assert(off[8] == len);
#if 0
  puts("check");
  for(int ic=0; ic<8; ic++){
    for(int i=off[ic]; i<off[ic+1]; i++){
      int c = (key[i] >> rshift) & 7;
      assert(c == ic);
    }
  }
#endif
  for(int ic=0; ic<8; ic++){
    num[ic] = off[ic+1] - off[ic];
  }
}



void scan_keys_brute( const long long int key[],
		      const int    len,
		      const int    rshift,
		      _out_ int    off[],
		      _out_ int    num[]){

  register int n0=0, n1=0, n2=0, n3=0, n4=0, n5=0, n6=0, n7=0;
#pragma omp parallel for reduction(+ : n0,n1,n2,n3,n4,n5,n6,n7)
  for(int i=0; i<len; i++){
    int c = (key[i] >> rshift) & 7;
    if(c&4){
      if(c&2){
	if(c&1) n7++;
	else    n6++;
      }else{
	if(c&1) n5++;
	else    n4++;
      }
    }else{
      if(c&2){
	if(c&1) n3++;
	else    n2++;
      }else{
	if(c&1) n1++;
	else    n0++;
      }
    }
  } // end for(i)
  num[0] = n0; num[1] = n1; num[2] = n2; num[3] = n3;
  num[4] = n4; num[5] = n5; num[6] = n6; num[7] = n7;

  off[0] = 0;
  for(int ic=0; ic<8; ic++){
    off[ic+1] = off[ic] + num[ic];
  }
}



void treeConstruction( pTreeTmp tree_tmp,
		       const int first_index, const int n,
		       long long int current_key, int key_level,
		       pICell walklist, int *nwalk, int ipflag,
		       const int parent_cell, const pTreeParam treeparam,
		       int *ncell_all, double (*p_cache)[4]){

  /* new version */

  int ncell[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int incell_first_index[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  long long int *key = tree_tmp->key;
  TreeCell *tree_cell = tree_tmp->tree_cell;
  key_level -= 3;

  // count n in child cell
#if 1
  scan_keys( key+first_index, n, key_level, 
	     incell_first_index, ncell);
#else
#if 1
  for( int i=first_index; i<(first_index+n); i++){
    int c = (key[i] >> key_level) & 0x7;
    ncell[c] ++;
  }
#else
  int ic = 0; 
  for(int i=first_index; i<(first_index+n); i++){
    int c = (key[i] >> key_level) & 0x7;
    while(ic != c){
      incell_first_index[++ic] = i;
    }
  }
  while(ic != 8){
    incell_first_index[ic++] = first_index + n;
  }
  for(int ic=0; ic<8; ic++){
    ncell[ic] = incell_first_index[ic+1] - incell_first_index[ic];
  }
#endif
#endif

  int first_cell = *ncell_all;
  tree_cell[parent_cell].next = first_cell;
  *ncell_all += 8;

  int first_index0 = first_index;
  for( int i=0; i<8; i++){
    int c = first_cell + i;
    tree_cell[c].n = ncell[i];
    tree_cell[c].first_index = first_index0;
    first_index0 += ncell[i];
  }

  double cmtmp[3] = {0.0, 0.0, 0.0};
  double mtmp = 0.0;

  for( int i=0; i<8; i++){
    if( ncell[i] == 0)  continue;

    int c = first_cell + i;
    int ipflag_this = 0;
    if(ipflag == 0){
      if( ncell[i] <= treeparam->tree_ncrit){
	int nw = *nwalk;
	walklist[nw].id = c;
	walklist[nw].key = ((current_key<<3) | i);
	walklist[nw].key_level = key_level;
	(*nwalk) ++;
	ipflag_this = 1;
      }
    }
    else{
      ipflag_this = 2;
    }

    if( ncell[i] > treeparam->tree_nleaf){  // construct tree recursively
      long long int tmp_key =  ((current_key<<3) | i);
      treeConstruction( tree_tmp, tree_cell[c].first_index, ncell[i],
		        tmp_key, key_level, 
			walklist, nwalk, ipflag_this, 
			c, treeparam, ncell_all, p_cache);
      double m = (double)tree_cell[c].m;
      cmtmp[0] += tree_cell[c].cm[0] * m;
      cmtmp[1] += tree_cell[c].cm[1] * m;
      cmtmp[2] += tree_cell[c].cm[2] * m;
      mtmp += m;
    }
    else{ // calculate gravity center
      tree_cell[c].next = -1;
      double cpos[3] = {0.0, 0.0, 0.0};
      double m = 0.0;
      int start = tree_cell[c].first_index;
      int end = start + ncell[i];
      for( int j=start; j<end; j++){
	double mk = p_cache[j][3];
	cpos[0] += p_cache[j][0] * mk;
	cpos[1] += p_cache[j][1] * mk;
	cpos[2] += p_cache[j][2] * mk;
	m += mk;
      }
      double minv = 1.0 / m;
      tree_cell[c].m = m;
      tree_cell[c].cm[0] = cpos[0] * minv;
      tree_cell[c].cm[1] = cpos[1] * minv;
      tree_cell[c].cm[2] = cpos[2] * minv;
      mtmp += m;
      cmtmp[0] += cpos[0];
      cmtmp[1] += cpos[1];
      cmtmp[2] += cpos[2];
    }
  }

  // child gravity center -> parent node
  double minv = 1.0 / (double)mtmp;
  tree_cell[parent_cell].m = mtmp;
  tree_cell[parent_cell].cm[0] = cmtmp[0] * minv;
  tree_cell[parent_cell].cm[1] = cmtmp[1] * minv;
  tree_cell[parent_cell].cm[2] = cmtmp[2] * minv;

}



#ifdef TREE2
static void treeConstructionWithoutGC( pTreeTmp tree_tmp,
				       const int first_index, const int n,
				       long long int current_key, int key_level,
				       const int nleaf, pICell leaf_list, int *ncell_leaf,
				       const int parent_cell,
				       int *ncell_all){

  /* new version */

  int ncell[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int incell_first_index[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  long long int *key = tree_tmp->key;
  TreeCell *tree_cell = tree_tmp->tree_cell;
  key_level -= 3;

  // count n in child cell
  scan_keys_brute( key+first_index, n, key_level, incell_first_index, ncell);

  int first_cell = *ncell_all;
  tree_cell[parent_cell].next = first_cell;
  *ncell_all += 8;

  int first_index0 = first_index;
  for( int i=0; i<8; i++){
    int c = first_cell + i;
    tree_cell[c].first_index = first_index0;
    tree_cell[c].n = ncell[i];
    first_index0 += ncell[i];
  }

  for( int i=0; i<8; i++){
    if( ncell[i] == 0)  continue;
    int c = first_cell + i;
    long long int tmp_key =  ((current_key<<3) | i);
    if( ncell[i] > nleaf){  // construct tree recursively
      treeConstructionWithoutGC( tree_tmp, tree_cell[c].first_index, ncell[i],
				 tmp_key, key_level, 
				 nleaf, leaf_list, ncell_leaf,
				 c, ncell_all);
    }
    else{
      tree_cell[c].next = -1;
      leaf_list[*ncell_leaf].id = c;
      leaf_list[*ncell_leaf].key = tmp_key;
      leaf_list[*ncell_leaf].key_level = key_level;
      (*ncell_leaf) ++;
    }
  }


}



static int getWalkList( const pTreeTmp tree_tmp,
			const int ncrit,
			int c,
			long long int current_key, 
			int key_level,
			pICell walklist){

  pTreeCell tree_cell = tree_tmp->tree_cell;

  int nwalk = 0;
  if( tree_cell[c].n <= ncrit){
    walklist[nwalk].key = current_key;
    walklist[nwalk].key_level = key_level;
    walklist[nwalk].id = c;
    nwalk ++;
  }
  else{
    key_level -= 3;
    for( int i=0; i<8; i++){
      int ic = tree_cell[c].next + i;
      if( tree_cell[ic].n == 0)  continue;
      long long int tmp_key =  ((current_key<<3) | i);
      nwalk += getWalkList( tree_tmp, ncrit, ic, tmp_key, key_level, walklist+nwalk);
    }
  }
  
  return nwalk;

}



static void calcGC( pTreeTmp tree_tmp,
		    const int ncell_all,
		    const int c,
		    const pTreeParam treeparam){

  TreeCell *tree_cell = tree_tmp->tree_cell;

#pragma omp parallel for
  for( int i=0; i<8; i++){
    int ic = tree_cell[c].next + i;
    if( tree_cell[ic].n <= treeparam->tree_nleaf)  continue;
    if( tree_cell[ic].next >= ncell_all){
      double cm[3] = {0.0, 0.0, 0.0};
      double m = 0.0;
      for( int i=0; i<8; i++){
	int ic2 = tree_cell[ic].next + i;
	double im = tree_cell[ic2].m;
	cm[0] += tree_cell[ic2].cm[0] * im;
	cm[1] += tree_cell[ic2].cm[1] * im;
	cm[2] += tree_cell[ic2].cm[2] * im;
	m += im;
      }
      double minv = 1.0 / m;
      tree_cell[ic].cm[0] = cm[0] * minv;
      tree_cell[ic].cm[1] = cm[1] * minv;
      tree_cell[ic].cm[2] = cm[2] * minv;
      tree_cell[ic].m = m;
    }
    else{
      calcGC( tree_tmp, ncell_all, ic, treeparam);
    }
  }

  double cm[3] = {0.0, 0.0, 0.0};
  double m = 0.0;
  for( int i=0; i<8; i++){
    int ic = tree_cell[c].next + i;
    double im = tree_cell[ic].m;
    cm[0] += tree_cell[ic].cm[0] * im;
    cm[1] += tree_cell[ic].cm[1] * im;
    cm[2] += tree_cell[ic].cm[2] * im;
    m += im;
  }
  double minv = 1.0 / m;
  tree_cell[c].cm[0] = cm[0] * minv;;
  tree_cell[c].cm[1] = cm[1] * minv;
  tree_cell[c].cm[2] = cm[2] * minv;
  tree_cell[c].m = m;

}



int treeConstructionParallel( pTreeTmp tree_tmp,
			      const int n,
			      pICell walklist,
			      const pTreeParam treeparam,
			      double (*p_cache)[4]){

  //  double nowtime = 0.0;
  //  getTime( &nowtime);

  pTreeCell tree_cell = tree_tmp->tree_cell;
  int key_level = 63;

  const int nleaf0 = NLEAF_PARALLEL_TREE_CONSTRUCTION;
#ifdef STATIC_ARRAY
  //  static int nleafmax = NUMBER_OF_PART * NLEAF / NLEAF_PARALLEL_TREE_CONSTRUCTION;
  static ICell leaf_list[nleafmax];
#else
  int nleafmax = tree_tmp->clistmax * treeparam->tree_nleaf / nleaf0;
  pICell leaf_list = new ICell[nleafmax];
#endif
  int ncell_leaf = 0;
  int ncell_all = 1;

  treeConstructionWithoutGC( tree_tmp, 0, n, 1, key_level, 
			     nleaf0, leaf_list, &ncell_leaf,
			     0, &ncell_all);

  int ncell_residual = tree_tmp->clistmax - ncell_all;
  int ntotal = 0;
  for( int i=0; i<ncell_leaf; i++){
    int c = leaf_list[i].id;
    if( n <= treeparam->tree_nleaf)  continue;
    ntotal += tree_cell[c].n;
  }

  static int start_leaf[NUMBER_OF_CHUNK_TREE_PARALLEL+1];
  static int ncell_leaf_threads[NUMBER_OF_CHUNK_TREE_PARALLEL+1];
  start_leaf[0] = 0;
  int nsum_thershold = ntotal / NUMBER_OF_CHUNK_TREE_PARALLEL;
  int nsum = 0;
  int it = 1;
  for( int i=0; i<ncell_leaf; i++){
    int c = leaf_list[i].id;
    int n = tree_cell[c].n;
    if( n <= treeparam->tree_nleaf)  continue;
    if( nsum > nsum_thershold){
      start_leaf[it] = i;
      it ++;
      //      fprintf( stderr, "%d\t%d\t%d\n", nsum, nsum_thershold, ntotal);
      nsum = 0;
    }
    nsum += n;
  }

  while( it <= NUMBER_OF_CHUNK_TREE_PARALLEL){
    start_leaf[it] = ncell_leaf;
    it ++;
  }

  for( int i=0; i<NUMBER_OF_CHUNK_TREE_PARALLEL; i++){
    ncell_leaf_threads[i] = start_leaf[i+1] - start_leaf[i];
    //fprintf( stderr, "%d\t%d\t%d\n", ncell_leaf_threads[i], start_leaf[i], ncell_leaf);
  }

  //  fprintf( stderr, "%e\n", getTime( &nowtime));

  int ncell_max = (int)( (double)ncell_residual * (double)nsum_thershold / (double)ntotal);
  //  cerr << ncell_max << endl;
#pragma omp parallel for
  for( int k=0; k<NUMBER_OF_CHUNK_TREE_PARALLEL; k++){
    if( ncell_leaf_threads[k] == 0)  continue;
    int ki = start_leaf[k];
    int cstart = leaf_list[ki].id;
    int k_ncell_all0 = ncell_all + ncell_max*k;
    int k_ncell_all = k_ncell_all0;
    tree_cell[cstart].next = k_ncell_all;
    for( int i=0; i<ncell_leaf_threads[k]; i++){
      int il = ki + i;
      int c = leaf_list[il].id;
      int n = tree_cell[c].n;
      if( n <= treeparam->tree_nleaf){
	if( n == 0)  continue;
	double cpos[3] = {0.0, 0.0, 0.0};
	double m = 0.0;
	int start = tree_cell[c].first_index;
	int end = start + n;
	for( int j=start; j<end; j++){
	  double mk = p_cache[j][3];
	  cpos[0] += p_cache[j][0] * mk;
	  cpos[1] += p_cache[j][1] * mk;
	  cpos[2] += p_cache[j][2] * mk;
	  m += mk;
	}
	double minv = 1.0 / m;
	tree_cell[c].m = m;
	tree_cell[c].cm[0] = cpos[0] * minv;
	tree_cell[c].cm[1] = cpos[1] * minv;
	tree_cell[c].cm[2] = cpos[2] * minv;
	continue;
      }
      int i_nwalk = 0;
      treeConstruction( tree_tmp, tree_cell[c].first_index, tree_cell[c].n,
			leaf_list[il].key, leaf_list[il].key_level,
			tree_tmp->walklist, &i_nwalk, 0,
			c, treeparam, &k_ncell_all, p_cache);
    }
    int ncell_use =  k_ncell_all - k_ncell_all0;
    if( ncell_use > ncell_max){
      //      fprintf( stderr, "%d\t%d\t%d\t%d\t%d\n", k_ncell_all0, k_ncell_all, 
      //	       ncell_use, n, ncell_max);
      cerr << ncell_use << "\t" << n << "\t" << ncell_max << endl;
    }
    assert( ncell_use < ncell_max);
  }

  //  fprintf( stderr, "%e\n", getTime( &nowtime));

  calcGC( tree_tmp, ncell_all, 0, treeparam);
  //  fprintf( stderr, "%e\n", getTime( &nowtime));

  int nwalk = getWalkList( tree_tmp, treeparam->tree_ncrit, 0, 
			   1, key_level, walklist);
  //  fprintf( stderr, "%e\n", getTime( &nowtime));

#ifndef STATIC_ARRAY
  delete [] leaf_list;
#endif

  return nwalk;


}
#endif


#ifndef TREE2
void treeConstruction( pTreeCell tree_cell, const long long int *key, 
		       const int *index, const pParticle particle, 
		       const int first_index, const int n,
		       const int tree_clistmask, int *col_cell_index,
		       long long int current_key, int key_level,
		       int *walklist, int *nwalk, int ipflag,
		       const int parent_cell, const pTreeParam treeparam){

  int i, j;
  int ncell[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int incell_first_index[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int cell_ipflag[8];
  int adrlist[8];

  key_level -= 3;
  int key_mask = 0x0000ff00;

  // count n in child cell
  for( i=first_index; i<(first_index+n); i++){
    int c = (key[i] >> key_level) & 0x7;
    if( ncell[c] == 0)  incell_first_index[c] = i;
    ncell[c] ++;
  }

  for( i=0; i<8; i++){
    if( ncell[i] == 0){
      tree_cell[parent_cell].zeroflag |= (1<<i);
    }
  }

  for( i=0; i<8; i++){
    if( ncell[i] == 0) continue;
    long long int tmp_key =  ((current_key<<3) | i);
    int c = tree_clistmask & tmp_key;
    int ncell_i = ncell[i];
    if(tree_cell[c].key != 0){ // collisions in the hash table
      int tmpc = c;
      c = *col_cell_index;
      int tmpnext = tree_cell[tmpc].next;
      tree_cell[tmpc].next = c;
      tree_cell[c].next = tmpnext;
      (*col_cell_index) ++;
    }

    tree_cell[c].key = tmp_key;
    tree_cell[c].first_index = incell_first_index[i];
    tree_cell[c].n = ncell_i;
    tree_cell[c].zeroflag |= (key_mask & ( key_level << 8));
    //tree_cell[c].zeroflag = key_level << 8;
    //tree_cell[c].next = -1;

    // if( ncell[i] == 0) continue;
    // determine i-particle
    if(ipflag == 0){
      if(ncell[i] <= treeparam->tree_ncrit){
	walklist[*nwalk] = c;
	(*nwalk) ++;
	cell_ipflag[i] = 1;
      }
      else{
	cell_ipflag[i] = 0;
      }
    }
    else{
      cell_ipflag[i] = 2;
    }
    adrlist[i] = c;
  }

  double cmtmp[3] = {0.0, 0.0, 0.0};
  double mtmp = 0.0;

  for( i=0; i<8; i++){
    if( ncell[i] == 0)  continue;

    long long int tmp_key =  ((current_key<<3) | i);
    int c = adrlist[i];

    if( ncell[i] > treeparam->tree_nleaf){  // construct tree recursively
      treeConstruction( tree_cell, key, index, particle, incell_first_index[i], ncell[i],
			tree_clistmask, col_cell_index, tmp_key, key_level, 
			walklist, nwalk, cell_ipflag[i], c, treeparam);
      double m = (double)tree_cell[c].m;
      cmtmp[0] += tree_cell[c].cm[0] * m;
      cmtmp[1] += tree_cell[c].cm[1] * m;
      cmtmp[2] += tree_cell[c].cm[2] * m;
      mtmp += m;
    }
    else{ // calculate gravity center
      double cpos[3] = {0.0, 0.0, 0.0};
      double m = 0.0;
      int start = incell_first_index[i];
      int end = start + ncell[i];
      for( j=start; j<end; j++){
	double posk[3];
#ifdef TREE_PARTICLE_CACHE
	posk[0] = p_cache[j][0];
	posk[1] = p_cache[j][1];
	posk[2] = p_cache[j][2];
	double mk = p_cache[j][3];
#else
	int k = index[j];
	getPos( &particle[k], posk);
	double mk = getMass( &particle[k], treeparam->uniform_mass);
#endif
	cpos[0] += posk[0] * mk;
	cpos[1] += posk[1] * mk;
	cpos[2] += posk[2] * mk;
	m += mk;
      }
      double minv = 1.0 / m;
      tree_cell[c].m = m;
      tree_cell[c].cm[0] = cpos[0] * minv;
      tree_cell[c].cm[1] = cpos[1] * minv;
      tree_cell[c].cm[2] = cpos[2] * minv;
      mtmp += m;
      cmtmp[0] += cpos[0];
      cmtmp[1] += cpos[1];
      cmtmp[2] += cpos[2];
    }
  }

  // child gravity center -> parent node
  double minv = 1.0 / (double)mtmp;
  tree_cell[parent_cell].m = mtmp;
  tree_cell[parent_cell].cm[0] = cmtmp[0] * minv;
  tree_cell[parent_cell].cm[1] = cmtmp[1] * minv;
  tree_cell[parent_cell].cm[2] = cmtmp[2] * minv;

}
#endif



#ifdef UNSTABLE
template <typename VEC, typename REAL>
struct morton_key{
	typedef unsigned long long key_t;
	key_t val;

	morton_key() : val(0) {}
	morton_key(VEC &vec, const REAL size=16.0){
		static key_t table[128] = {
			#include "key_table"
		};
#if 1
		int xi = int(vec[0]);
		int yi = int(vec[1]);
		int zi = int(vec[2]);
#else
		const REAL scale = (1<<20) / size;
		int xi = int((vec[0] + size) * scale);
		int yi = int((vec[1] + size) * scale);
		int zi = int((vec[2] + size) * scale);
#endif
		assert((xi >> 21) == 0);
		assert((yi >> 21) == 0);
		assert((zi >> 21) == 0);
		key_t xkey = (table[xi&127]) | (table[(xi>>7)&127] << 21) | (table[(xi>>14)&127] << 42);
		key_t ykey = (table[yi&127]) | (table[(yi>>7)&127] << 21) | (table[(yi>>14)&127] << 42);
		key_t zkey = (table[zi&127]) | (table[(zi>>7)&127] << 21) | (table[(zi>>14)&127] << 42);
		val = (xkey<<2) | (ykey<<1) | zkey;
	}
};
#endif



/* make morton key */
double makeKey( pParticle particle, long long int *key, const int n){

  double maxx = 0.0;
  double xscale = 0;
  maxx = 1.0;
  xscale = (double)(1<<21);
  long long int tmpone = ((long long int)1) << 63;

#pragma omp parallel for
  for( int i=0;i<n;i++){
    double x[3];
    getPos( &particle[i], x);
#ifdef UNSTABLE
    int tmpi0[3];
    tmpi0[0] = (int)(x[0] * xscale);
    tmpi0[1] = (int)(x[1] * xscale);
    tmpi0[2] = (int)(x[2] * xscale);
    morton_key<int [3], double> mk( tmpi0, 1.0);
    key[i] = tmpone | mk.val;
#else //UNSTABLE
    long long int tmpi[3];
    tmpi[0] = (int)(x[0] * xscale);
    tmpi[1] = (int)(x[1] * xscale);
    tmpi[2] = (int)(x[2] * xscale);
    key[i] = tmpone;
    key[i] |= ((long long int)(((tmpi[0]&0x1)<<2)|((tmpi[1]&0x1)<<1)|((tmpi[2]&0x1))));
    key[i] |= ((long long int)(((tmpi[0]&0x2)<<2)|((tmpi[1]&0x2)<<1)|((tmpi[2]&0x2)))<<2);
    key[i] |= ((long long int)(((tmpi[0]&0x4)<<2)|((tmpi[1]&0x4)<<1)|((tmpi[2]&0x4)))<<4);
    key[i] |= ((long long int)(((tmpi[0]&0x8)<<2)|((tmpi[1]&0x8)<<1)|((tmpi[2]&0x8)))<<6);
    key[i] |= ((long long int)(((tmpi[0]&0x10)<<2)|((tmpi[1]&0x10)<<1)|((tmpi[2]&0x10)))<<8);
    key[i] |= ((long long int)(((tmpi[0]&0x20)<<2)|((tmpi[1]&0x20)<<1)|((tmpi[2]&0x20)))<<10);
    key[i] |= ((long long int)(((tmpi[0]&0x40)<<2)|((tmpi[1]&0x40)<<1)|((tmpi[2]&0x40)))<<12);
    key[i] |= ((long long int)(((tmpi[0]&0x80)<<2)|((tmpi[1]&0x80)<<1)|((tmpi[2]&0x80)))<<14);
    key[i] |= ((long long int)(((tmpi[0]&0x100)<<2)|((tmpi[1]&0x100)<<1)|((tmpi[2]&0x100)))<<16);
    key[i] |= ((long long int)(((tmpi[0]&0x200)<<2)|((tmpi[1]&0x200)<<1)|((tmpi[2]&0x200)))<<18);
    key[i] |= ((long long int)(((tmpi[0]&0x400)<<2)|((tmpi[1]&0x400)<<1)|((tmpi[2]&0x400)))<<20);
    key[i] |= ((long long int)(((tmpi[0]&0x800)<<2)|((tmpi[1]&0x800)<<1)|((tmpi[2]&0x800)))<<22);
    key[i] |= ((long long int)(((tmpi[0]&0x1000)<<2)|((tmpi[1]&0x1000)<<1)|((tmpi[2]&0x1000)))<<24);
    key[i] |= ((long long int)(((tmpi[0]&0x2000)<<2)|((tmpi[1]&0x2000)<<1)|((tmpi[2]&0x2000)))<<26);
    key[i] |= ((long long int)(((tmpi[0]&0x4000)<<2)|((tmpi[1]&0x4000)<<1)|((tmpi[2]&0x4000)))<<28);
    key[i] |= ((long long int)(((tmpi[0]&0x8000)<<2)|((tmpi[1]&0x8000)<<1)|((tmpi[2]&0x8000)))<<30);
    key[i] |= ((long long int)(((tmpi[0]&0x10000)<<2)|((tmpi[1]&0x10000)<<1)|((tmpi[2]&0x10000)))<<32);
    key[i] |= ((long long int)(((tmpi[0]&0x20000)<<2)|((tmpi[1]&0x20000)<<1)|((tmpi[2]&0x20000)))<<34);
    key[i] |= ((long long int)(((tmpi[0]&0x40000)<<2)|((tmpi[1]&0x40000)<<1)|((tmpi[2]&0x40000)))<<36);
    key[i] |= ((long long int)(((tmpi[0]&0x80000)<<2)|((tmpi[1]&0x80000)<<1)|((tmpi[2]&0x80000)))<<38);
    key[i] |= ((long long int)(((tmpi[0]&0x100000)<<2)|((tmpi[1]&0x100000)<<1)|((tmpi[2]&0x100000)))<<40);
#endif // UNSTABLE
    //if( i==0)  cerr << key[i] << "\t" << key0 << endl;
  }

  return maxx;

}



#ifndef GRAPE_OFF
void calculateForce2( pParticle particle, const int *ilist, const int nilist,
		      const pJList jlist, const int njlist, float (*a)[3],
		      double *cicell){

  static double xi[NCRIT][3];
  static double atmp[NCRIT][3];
  static double ptmp[NCRIT];

  int npipes = g5_get_number_of_pipelines();
  int jmemsize = g5_get_jmemsize();

  /* offset */
  for( int i=0; i<nilist; i++){
    int iii = ilist[i];
#ifdef TREE_PARTICLE_CACHE
    xi[i][0] = p_cache[iii][0];
    xi[i][1] = p_cache[iii][1];
    xi[i][2] = p_cache[iii][2];
#else
    getPos( &particle[iii], xi[i]);
#endif
    xi[i][0] -= cicell[0];
    xi[i][1] -= cicell[1];
    xi[i][2] -= cicell[2];
  }
  for( int i=0; i<njlist; i++){
    jlist->x[i][0] -= cicell[0];
    jlist->x[i][1] -= cicell[1];
    jlist->x[i][2] -= cicell[2];
  }

  for( int j=0; j<njlist; j+=jmemsize){

    int njlist1 = jmemsize;
    if( (jmemsize+j)>njlist)  njlist1 = njlist - j;

    g5_set_xmj( 0, njlist1, &jlist->x[j], &jlist->m[j]);
    g5_set_n(njlist1);

    for( int i=0; i<nilist; i+=npipes){
      int nn = npipes;
      if((i+npipes)>nilist) nn = nilist-i;
      g5_set_xi(nn,xi+i);
      g5_run();
      g5_get_force(nn,atmp+i,ptmp+i);
    }
    for( int i=0; i<nilist; i++){
      a[i][0] += atmp[i][0];
      a[i][1] += atmp[i][1];
      a[i][2] += atmp[i][2];
#ifdef CALCPOT
      particle[i].pot -= ptmp[i];
#endif
    } 
  }
}
#endif



static inline double plummerCutoffFunc( double r){

  r /= (SFT_FOR_PM*0.50);

  double fc;

  if (r<1.0){
    fc = r*r*r*(224.0+r*r*(-224.0+r*(70.0+r*(48.0-r*21.0))))/140.0;
  }
  else if (r<2.0) {
    fc = (12.0+r*r*(-224.0+r*(896.0+r*(-840.0+r*(224.0+r*(70.0+r*(-48.0+r*7.0)))))))/140.0;
  }
  else{
    fc = 1.0;
  }

  fc = 1.0-fc;

  return fc;

}



static inline double plummerCutoffFunc(double r, double r2){

  r  *= 1.0/(SFT_FOR_PM*0.50);
  if( r>2){
    return 0;
  }
  r2 *= 1.0/(SFT_FOR_PM*SFT_FOR_PM*0.25);

  const double r3 = r * r2;
  const double s = (r-1.0) + fabs(r-1.0);
  const double s2 = s*s;
  const double s6 = s2*s2*s2;
  const double poly1 = 1. + r3*(-8./5. + r2*(8./5. + r*(-1./2. +
                                                        r*(-12./35. + r*(3./20.)))));
  const double poly2 = s6 * (3./35./64. + r*(18./35./64. + r*(1./5./64.)));
  return poly1 - poly2;

}



void calculateForceHost( pParticle particle, const int *ilist, const int nilist,
			 const pJList jlist, const int njlist, float (*a)[3], 
			 const double eps){

  double eps2 = eps * eps;

  int i,j;
  for( i=0; i<nilist; i++){
    int ii = ilist[i];
    double pos[3];
#ifdef TREE_PARTICLE_CACHE
    pos[0] = p_cache[ii][0];
    pos[1] = p_cache[ii][1];
    pos[2] = p_cache[ii][2];
#else
    getPos( &particle[ii], pos);
#endif
    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;
    for( j=0; j<njlist; j++){
      double dx = jlist->x[j][0] - pos[0];
      double dy = jlist->x[j][1] - pos[1];
      double dz = jlist->x[j][2] - pos[2];
      double r2 = dx*dx + dy*dy + dz*dz + eps2;
      double rinv = 1.0/sqrt(r2);
      double r2inv = rinv*rinv;
      double g = plummerCutoffFunc(r2*rinv, r2);
      double r3inv = r2inv * rinv;
#ifndef KCOMPUTER
      ax += g * jlist->m[j] * r3inv * dx;
      ay += g * jlist->m[j] * r3inv * dy;
      az += g * jlist->m[j] * r3inv * dz;
#else
      ax += g * jlist->x[j][3] * r3inv * dx;
      ay += g * jlist->x[j][3] * r3inv * dy;
      az += g * jlist->x[j][3] * r3inv * dz;
#endif
    }
    a[i][0] = ax;
    a[i][1] = ay;
    a[i][2] = az;
  } 


}



#define PRINT_TANIKAWA
#ifdef PRINT_TANIKAWA
static bool first_tanik = true;
static char file_pos[64], file_app[64], file_apm[64];
static FILE *fppos, *fpapp, *fpapm;
#endif

/* make interaction list 
 * and calculate PP(using Tree) and PM(using mesh potential) force
 * and velocity kick
 */
#ifndef _OPENMP
void calculateForceUsingTreeAndPMForce2(const pTreeCell tree_cell, const int ncell,
					const pParticle particle, const int *index,
					const pTreeTmp tree_tmp, const int nwalk,
					const double maxx, const int tree_clistmask,
					const double eps, const pTreeParam treeparam,
					PMForce &pm,
					const int ncalc, const int npart,
					const double vfac, 
					const double afac,
					pRunParam this_run){

#ifdef MAKE_LIST_PROF
  t_calc_theta = 0.0;
  t_calc_adr = 0.0;
  t_push_intr = 0.0;
#endif

  double nowtime = 0.0;
  static double time1 = 0.0;
  static double time2 = 0.0;
  static double time3 = 0.0;
  static double time4 = 0.0;

  int i, j;
  static int ilist[NCRIT];
  static JList jlist;
  static float a_ilist[NCRIT][3];

  int nisum = 0;
  int njsum = 0;
  long long int ninter = 0;
  int nimin = 1 << 30;
  int nimax = 0;
  int njmin = 1 << 30;
  int njmax = 0;

  double maxa2 = 0.0;
  double maxv2 = 0.0;
  double mina2 = 1.0e30;
  double minv2 = 1.0e30;
  double suma2 = 0.0;
  double sumv2 = 0.0;

#ifdef TREE2
  pICell walklist = tree_tmp->walklist;
#else
  int *walklist = tree_tmp->walklist;
#endif

#ifdef PRINT_TANIKAWA
    if(first_tanik) {
        int theta_int = (int)(treeparam->tree_theta * 10.);
        sprintf(file_pos, "greem_n32.theta%02d.pos.%04d", theta_int, this_run->inode);
        fppos = fopen(file_pos, "w");
        sprintf(file_app, "greem_n32.theta%02d.app.%04d", theta_int, this_run->inode);
        fpapp = fopen(file_app, "w");
        sprintf(file_apm, "greem_n32.theta%02d.apm.%04d", theta_int, this_run->inode);
        fpapm = fopen(file_apm, "w");
        if(this_run->inode == 0) {
            printf("print_tanikawa %02d\n", theta_int);
        }
    }
#endif

  time1 = time2 = time3 = time4 = 0.0;
  for( i=0; i<nwalk; i++){

#ifdef TREE2
    int c = walklist[i].id;
    long long int key = walklist[i].key;
    int key_level = walklist[i].key_level;
#else
    int c = walklist[i];
#endif

    // assign i-particle
    int istart = tree_cell[c].first_index;
    int nilist = tree_cell[c].n;
    int nilist2 = 0;
    for( j=0; j<nilist; j++){
      if( index[istart+j] >= ncalc)  continue;
#ifdef TREE_PARTICLE_CACHE
      ilist[nilist2] = istart + j;
#else
      ilist[nilist2] = index[istart+j];
#endif
      a_ilist[nilist2][0] = a_ilist[nilist2][1] = a_ilist[nilist2][2] = 0.0;
      nilist2 ++;
    }
    if( nilist2 == 0)  continue;

    // assign j-particle
    int njlist = 0;
    double cicell[3], ilength[3];
    double cell_length = 1.0;
    double cell_length2 = cell_length * cell_length;
    getTime(&nowtime);

#ifdef TREE2
    getCellRange2( &tree_cell[c], cicell, ilength,
		   key_level, key, tree_tmp->maxx);
#else
    getCellRange2( &tree_cell[c], cicell, ilength);
#endif
    /*
    fprintf( stdout, "%e\t%e\t%e\t%e\t%e\t%e\n",
	     cicell[0], cicell[1], cicell[2],
	     ilength[0], ilength[1], ilength[2]);
    */
#ifdef TREE2
    float peri[4] = {0.0, 0.0, 0.0, 0.0};
    makeInteractionList( tree_cell, c, cicell, ilength,
			 0, peri, cell_length2, &jlist, &njlist,
			 treeparam, p_cache);
    /*
    MPI_Barrier(MPI_COMM_WORLD);
    cerr << njlist << endl;
    exit(1);
    */
#else
    makeInteractionList(c, cell_length2, cicell, ilength, 
			tree_cell, &jlist, &njlist, 
			1, tree_clistmask, particle, index, treeparam);
#endif
    time1 += getTime(&nowtime);
#ifdef CALCPOT
    for( j=0; j<nilist2; j++){
      int ii = ilist[j];
      particle[ii].pot = 0.0;
    }
#else
#ifdef GRAPE_OFF
    calculateForceHost( particle, ilist, nilist2, &jlist, njlist, a_ilist, eps);
#else
    //    cerr << nilist2 << "\t" << njlist << endl;
    calculateForce2( particle, ilist, nilist2, &jlist, njlist, a_ilist, cicell);
#endif // GRAPE_OFF
#endif // CALCPOT
#ifdef CALCPOT
    calculateForcePotSSE( particle, ilist, nilist2, &jlist, njlist, eps);
#endif // CALCPOT
    time2 += getTime(&nowtime);

    for( j=0; j<nilist2; j++){
      int ii = ilist[j];
      float ia[3] = { 0.0, 0.0, 0.0};
#ifndef MULTI_TIMESTEP
      Particle p;
#ifdef TREE_PARTICLE_CACHE
      p.xpos = p_cache[ii][0];
      p.ypos = p_cache[ii][1];
      p.zpos = p_cache[ii][2];
      ii = index[ii];
      pm.forceInterpolation( &p, ia);
#else // TREE_PARTICLE_CACHE
      pm.forceInterpolation( &particle[ii], ia);
#endif // TREE_PARTICLE_CACHE
#endif //MULTI_TIMESTEP

#ifdef PRINT_TANIKAWA
      if(first_tanik){
          fprintf(fppos, "%8lld", particle[ii].id);
          fprintf(fppos, " %+.10e %+.10e %+.10e", particle[ii].xpos, particle[ii].ypos, particle[ii].zpos);
          fprintf(fppos, "\n");
          fprintf(fpapp, "%8lld", particle[ii].id);
          fprintf(fpapp, " %+.10e %+.10e %+.10e", a_ilist[j][0], a_ilist[j][1], a_ilist[j][2]);
          fprintf(fpapp, "\n");
          fprintf(fpapm, "%8lld", particle[ii].id);
          fprintf(fpapm, " %+.10e %+.10e %+.10e", ia[0], ia[1], ia[2]);
          fprintf(fpapm, "\n");
      }
#endif

#if 1
      ia[0] += a_ilist[j][0];
      ia[1] += a_ilist[j][1];
      ia[2] += a_ilist[j][2];
#endif

#ifdef NOACC
#ifdef MULTI_TIMESTEP
      kick( &particle[ii], ia, afac);
#else  
      kick( &particle[ii], ia, vfac, afac);
#endif
#else  // NOACC
#endif

      double asize = ia[0]*ia[0] + ia[1]*ia[1] + ia[2]*ia[2];
      double vsize = particle[ii].xvel*particle[ii].xvel +
                     particle[ii].yvel*particle[ii].yvel +
	             particle[ii].zvel*particle[ii].zvel;

#ifndef NOACC
      particle[ii].xacc = ia[0];
      particle[ii].yacc = ia[1];
      particle[ii].zacc = ia[2];
#endif

      if( maxa2 < asize)  maxa2 = asize;
      if( maxv2 < vsize)  maxv2 = vsize;
      if( mina2 > asize)  mina2 = asize;
      if( minv2 > vsize)  minv2 = vsize;
      suma2 += asize;
      sumv2 += vsize;
      //time3 += getTime(&nowtime);
    }

    nisum += nilist2;
    njsum += njlist;
    ninter += nilist2*njlist;
    if( nimin > nilist2)  nimin = nilist2;
    if( nimax < nilist2)  nimax = nilist2;
    if( njmin > njlist)  njmin = njlist;
    if( njmax < njlist)  njmax = njlist;
    time4 += getTime(&nowtime);

  }

#ifdef PRINT_TANIKAWA
  if(first_tanik) {
      fclose(fppos);
      fclose(fpapp);
      fclose(fpapm);
      first_tanik = false;
  }
#endif

  this_run->t.pp_make_intr_list = time1;
  this_run->t.pp_calc_force = time2;
#ifdef MULTI_TIMESTEP
  this_run->t.pp_vel_kick = time4;
#else
  this_run->t.pp_vel_kick = time3;
  this_run->t.pm_interpolation = time4;
#endif

  this_run->nwalk = nwalk;
  this_run->nisum = nisum;
  this_run->nimin = nimin;
  this_run->nimax = nimax;
  this_run->njsum = njsum;
  this_run->njmin = njmin;
  this_run->njmax = njmax;
  this_run->ninteraction = ninter;

  this_run->a2max = maxa2;
  this_run->v2max = maxv2;
  this_run->a2min = mina2;
  this_run->v2min = minv2;
  this_run->a2ave = suma2 / this_run->npart;
  this_run->v2ave = sumv2 / this_run->npart;
  this_run->ninteraction = ninter;

#ifdef MAKE_LIST_PROF
  fprintf( stderr, "makeInteractionList: %e\n", time1);
  fprintf( stderr, "makeInteractionList_push: %e\n", t_push_intr);
  fprintf( stderr, "makeInteractionList_adr: %e\n", t_calc_adr);
  fprintf( stderr, "makeInteractionList_theta: %e\n", t_calc_theta);
#endif


}
#endif



void copyPCache( const pParticle particle, 
		 const int *index, const int n, 
		 const double uniform_mass){

#ifdef TREE_PARTICLE_CACHE
#pragma omp parallel for
  for( int i=0; i<n; i++){
    int ii = index[i];
    p_cache[i][0] = particle[ii].xpos;
    p_cache[i][1] = particle[ii].ypos;
    p_cache[i][2] = particle[ii].zpos;
    p_cache[i][3] = getMass( &particle[ii], uniform_mass);
  }
#endif

}


/* construct tree and calculate gravity and velocity kick*/
void calculateForceUsingMakingTreeAndPMForce2( pTreeTmp tree_tmp, 
					       const pParticle particle, 
					       pRunParam this_run,
					       const double eps,  
					       const pTreeParam treeparam,
					       PMForce &pm,
					       const int ncalc,
					       const int n,
					       const double vfac,
					       const double afac){


  double nowtime;

  TreeCellTreeCell(tree_tmp->tree_cell, tree_tmp->clistmax);
  getTime( &nowtime);
  fprintf_verbose( stderr, "#PP makeKey start\n");
#ifdef BUFFER_FOR_TREE
  double maxx = makeKey( particle+ncalc, &(tree_tmp->key[ncalc]), n-ncalc); 
#else
  double maxx = makeKey( particle, tree_tmp->key, n); 
#endif
  tree_tmp->maxx = maxx;
  this_run->t.pp_make_key = getTime(&nowtime);

  fprintf_verbose( stderr, "#PP qsort start\n");
  qsortPivotWrapper(tree_tmp->index, tree_tmp->key, n, NLEVEL_KEY_SORT, NCHUNK_QSORT);
  this_run->t.pp_key_sort = getTime(&nowtime);

#ifdef TREE_PARTICLE_CACHE
  copyPCache( particle, tree_tmp->index, n, treeparam->uniform_mass);
#endif

  fprintf_verbose( stderr, "#PP treeConstruction start\n");
  int col_cell_index = tree_tmp->tree_clistmask + 1;  //for hash collision
  int nwalk = 0;         // for i-particle
  int key_level = 63;
  int ipflag = 0;   //1:i-particle 0:upper 2:lower
#ifdef TREE2
  tree_tmp->tree_cell[0].first_index = 0;
  tree_tmp->tree_cell[0].n = n;
  tree_tmp->tree_cell[0].next = 1;
#ifdef TREECONSTRUCTION_PARALLEL
  nwalk = treeConstructionParallel( tree_tmp,
				    n,
				    tree_tmp->walklist,
				    treeparam,
				    p_cache);
#else
  int ncell_all = 1;
  treeConstruction( tree_tmp, 0, n, 1, key_level, 
		    tree_tmp->walklist, &nwalk, ipflag,
		    0, treeparam, &ncell_all, p_cache);
#endif


#else
  int key_mask = 0x0000ff00;
  long long int current_key = 1;
  tree_tmp->tree_cell[current_key].n = n;
  tree_tmp->tree_cell[current_key].key = current_key;
  tree_tmp->tree_cell[current_key].zeroflag |= (key_mask & ( key_level << 8));
  treeConstruction( tree_tmp->tree_cell, tree_tmp->key, tree_tmp->index, particle, 0, n, 
		    tree_tmp->tree_clistmask, &col_cell_index, current_key, key_level,
		    tree_tmp->walklist, &nwalk, ipflag, 1, treeparam);
#endif

  this_run->t.pp_construct_tree = getTime(&nowtime);

  assert( col_cell_index < tree_tmp->clistmax);
  assert( nwalk < tree_tmp->walklist_size);

  fprintf_verbose( stderr, "#PP force loop start\n");
  calculateForceUsingTreeAndPMForce2( tree_tmp->tree_cell, col_cell_index, 
				      particle, tree_tmp->index, 
				      tree_tmp, nwalk, maxx, 
				      tree_tmp->tree_clistmask, 
				      eps, treeparam, pm,
				      ncalc, n,
				      vfac, afac, this_run);



}



void getNextDt( pRunParam this_run, const double dtmid, const float eps){

#ifdef CONSTANT_TIMESTEP
  this_run->tnow += dtmid;
  this_run->dtime = constant_timestep;
  update_now(this_run);
#else
  int step = this_run->nstep;
  double a2max = this_run->a2max;
  double v2max = this_run->v2max;
  double a2min = this_run->a2min;
  double v2min = this_run->v2min;
  double a2ave = this_run->a2ave * this_run->npart;
  double v2ave = this_run->v2ave * this_run->npart;
  double a2max2, v2max2, a2min2, v2min2;
  double a2ave2 = 0.0;
  double v2ave2 = 0.0;
  MPI_Allreduce( &a2max, &a2max2, 1, MPI_DOUBLE, MPI_MAX, this_run->MPI_COMM_INTERNAL);
  MPI_Allreduce( &v2max, &v2max2, 1, MPI_DOUBLE, MPI_MAX, this_run->MPI_COMM_INTERNAL);
  MPI_Allreduce( &a2min, &a2min2, 1, MPI_DOUBLE, MPI_MIN, this_run->MPI_COMM_INTERNAL);
  MPI_Allreduce( &v2min, &v2min2, 1, MPI_DOUBLE, MPI_MIN, this_run->MPI_COMM_INTERNAL);
  a2max = a2max2;
  v2max = v2max2;
  a2min = a2min2;
  v2min = v2min2;
  a2ave = a2ave2 / this_run->npart_total;
  v2ave = v2ave2 / this_run->npart_total;

  this_run->tnow += dtmid;
  update_now(this_run);
  this_run->dtprev = this_run->dtime;
  double dta = sqrt(eps/sqrt(a2max));
  double dtv = eps/sqrt(v2max);

  double znow = this_run->znow;
  double eta = ETA_TIMESTEP;
  if( znow > Z_SWITCH_ETA){
    eta = ETA_HIGHZ;
  }

  double anow = this_run->anow;
  double anow2 = anow * anow;
  dta *= anow*sqrt(anow);
  dtv *= anow2;
  double dttmp = dta;
  this_run->dtime = dttmp * eta;

  this_run->dta = dta*eta;
  this_run->dtv = dtv*eta;

  if(step == 0) this_run->dtime_ini = this_run->dtime;


#endif

  if( LOADBALANCE_METHOD == 0){
    long long int ninteraction = 0;
    MPI_Allreduce( &this_run->ninteraction, &ninteraction, 1,
		   MPI_LONG_LONG_INT, MPI_SUM, this_run->MPI_COMM_INTERNAL);
    this_run->nrate = (double)this_run->ninteraction / (double)ninteraction;
  }
  else if( LOADBALANCE_METHOD==1  ||  LOADBALANCE_METHOD==3  || LOADBALANCE_METHOD==4){
    double t_calc_sum = 0.0;
    double calc_time = this_run->t.pp + this_run->t.pm;
    if( LOADBALANCE_METHOD == 4){
      calc_time = this_run->t.pp + this_run->t.drift;
    }
    MPI_Allreduce( &calc_time, &t_calc_sum, 1,
		   MPI_DOUBLE, MPI_SUM, this_run->MPI_COMM_INTERNAL);
    this_run->nrate = calc_time / t_calc_sum;
    if( LOADBALANCE_METHOD == 3){
      //      fprintf( stderr, "#L3\n");
      double nrate2 = (double)this_run->npart / (double)this_run->npart_total;
      nrate2 /= SAMPLING_LOWER_LIMIT_FACTOR;
      fprintf( stderr, "%d %d %e %e\n", 
	       this_run->inode, this_run->npart, this_run->nrate, nrate2);
      if(this_run->nrate < nrate2){
        this_run->nrate = nrate2;
      }
      fprintf( stderr, "%d %d %e %e\n", 
	       this_run->inode, this_run->npart, this_run->nrate, nrate2);
      double nrate2_sum = 0.0;
      MPI_Allreduce( &this_run->nrate, &nrate2_sum, 1,
                     MPI_DOUBLE, MPI_SUM, this_run->MPI_COMM_INTERNAL);
      this_run->nrate /= nrate2_sum;
    }
  }
  else if( LOADBALANCE_METHOD == 2){
    this_run->nrate = (double)this_run->npart / (double)this_run->npart_total;
  }

}



void calcPMKick( pParticle particle, pRunParam this_run,
		 PMForce &pm, 
		 const double vfac,
		 const double afac){

  double nowtime;
  getTime(&nowtime);

#pragma omp parallel for
  for( int i=0; i<this_run->npart; i++){
    float a[3] = { 0.0, 0.0, 0.0};
    pm.forceInterpolation( particle+i, a);
    kick( particle+i, a, vfac, afac);
  }
  this_run->t.pm_interpolation = getTime(&nowtime);

}



void calcPPKick( pParticle particle, pRunParam this_run,
		 TreeParam treeparam, PMForce &pm, 
		 const double eps,
		 const double vfac,
		 const double afac){

  double nowtime, nowtime2;
  int n = this_run->npart;
  int ncalc = n;

#ifdef TREE_PARTICLE_CACHE
  static int first_call = 0;
  if( first_call == 0){
    p_cache = new double[NUMBER_OF_PART][4];
    first_call = 1;
  }
#endif

#ifdef BUFFER_FOR_TREE
  static TreeTmp tree_tmp;
  static int first_call2 = 0;
  if( first_call2 == 0){
    TreeTmpTreeTmp( &tree_tmp, &treeparam, NUMBER_OF_PART);
    first_call2 = 1;
  }
#endif

  getTime( &nowtime);
  getTime( &nowtime2);
  fprintf_verbose( stderr, "#PP getBoundaryParticle start\n");
  int bn = getBoundaryParticle( particle, 
				this_run, 
				&treeparam
#ifdef BUFFER_FOR_TREE
				,tree_tmp
#endif
				);
  this_run->t.pp_get_boundary = getTime(&nowtime2);
  this_run->bn = bn;

  if( (n+bn) > NUMBER_OF_PART){
    char *ctmp = (char *)"error.log";
    FILE *fout = fopen( ctmp, "aw");
    fprintf_verbose( stderr, "#n+bn > NUMBER_OF_PART\n");
    fprintf( fout, "n+bn > NUMBER_OF_PART\n");
    fclose(fout);
    exit(1);
  }

  correctBoundaryCondition( particle+n, bn);

  fprintf_verbose( stderr, "#PP TreeTmpTreeTmp start\n");
#ifdef BUFFER_FOR_TREE
  TreeTmpTreeTmp( &tree_tmp, bn, n);
  //TreeTmpTreeTmp( &tree_tmp, n+bn, 0);
#else
  TreeTmp tree_tmp;
  TreeTmpTreeTmp( &tree_tmp, &treeparam, n+bn);
#endif
  if( this_run->nstep == 0){
    TreeTmpPrintMemsize( &tree_tmp, n+bn);
  }
  fprintf_verbose( stderr, "#PP main start\n");
  calculateForceUsingMakingTreeAndPMForce2( &tree_tmp, particle, this_run, 
					   eps, &treeparam, 
					   pm, ncalc, ncalc+bn,
					   vfac, afac);
  this_run->t.pp = getTime(&nowtime);

#ifndef BUFFER_FOR_TREE
  TreeTmpDTreeTmp( &tree_tmp);  
#endif


}



void calcTreePMForce2( pParticle particle, pRunParam this_run,
		       TreeParam treeparam, PMForce &pm,
		       const float eps){


  int step = this_run->nstep;
  
  if(step==0) this_run->dtprev = this_run->dtime;
  double dtmid = (this_run->dtprev + this_run->dtime) * 0.50;
  double vfac, afac;
  getIntegralFactor( this_run, dtmid, &vfac, &afac);
  //  if( this_run->inode == 0)  cerr << "vfac,afac" << vfac << "\t" << afac << endl;

  int n = this_run->npart;
  int ncalc = n;

#ifdef TREE_PARTICLE_CACHE
  static int first_call = 0;
  if( first_call == 0){
    p_cache = new double[NUMBER_OF_PART][4];
    first_call = 1;
  }
#endif

#ifdef BUFFER_FOR_TREE
  static TreeTmp tree_tmp;
  static int first_call2 = 0;
  if( first_call2 == 0){
    TreeTmpTreeTmp( &tree_tmp, &treeparam, NUMBER_OF_PART);
    first_call2 = 1;
  }
#endif

  int bn = 0;
  double nowtime;
  getTime( &nowtime);
  if( this_run->nnode != 1){
    double nowtime2;
    getTime( &nowtime2);
    fprintf_verbose( stderr, "#PP getBoundaryParticle start\n");
    bn = getBoundaryParticle( particle, 
			      this_run, 
			      &treeparam
#ifdef BUFFER_FOR_TREE
			      ,tree_tmp
#endif
			      );
    this_run->t.pp_get_boundary = getTime(&nowtime2);
  }
  else{
    bn = 0;
  }
  this_run->bn = bn;

  if( (n+bn) > NUMBER_OF_PART){
    char *ctmp = (char *)"error.log";
    FILE *fout = fopen( ctmp, "aw");
    fprintf_verbose( stderr, "#n+bn > NUMBER_OF_PART\n");
    fprintf( fout, "n+bn > NUMBER_OF_PART\n");
    fclose(fout);
    exit(1);
  }

  correctBoundaryCondition( particle+n, bn);

  /* initialize tree */
  fprintf_verbose( stderr, "#PP TreeTmpTreeTmp start\n");

#ifdef BUFFER_FOR_TREE
  TreeTmpTreeTmp( &tree_tmp, bn, n);
  //TreeTmpTreeTmp( &tree_tmp, n+bn, 0);
#else
  TreeTmp tree_tmp;
  TreeTmpTreeTmp( &tree_tmp, &treeparam, n+bn);
#endif
  if( step == 0){
    TreeTmpPrintMemsize( &tree_tmp, n+bn);
  }
  fprintf_verbose( stderr, "#PP main start\n");
  calculateForceUsingMakingTreeAndPMForce2( &tree_tmp, particle, this_run, 
					   eps, &treeparam, 
					   pm, ncalc, ncalc+bn,
					   vfac, afac);
  this_run->t.pp = getTime(&nowtime);

  /* drift */
  getTime( &nowtime);
  drift( particle, n, this_run->dtime);
  this_run->t.drift = getTime(&nowtime);

  getNextDt( this_run, dtmid, eps);
  this_run->t.next_dt = getTime(&nowtime);

  //  fprintf( stderr, "%.15lf\n", this_run->dtime);

#ifndef BUFFER_FOR_TREE
  TreeTmpDTreeTmp( &tree_tmp);  
#endif

}


    } // namespace ParticleMesh
}     // namespace ParticleSimulator
