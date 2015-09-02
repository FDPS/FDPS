#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "avx.h"
#include "avx2.h"
#include "avx_type.h"
#include "gravity.h"

#define IPARA 2
#define JPARA 4
#define SEC 2.0
#define THD 3.0

#define JMEMSIZE 262144
#define ALIGN32  __attribute__ ((aligned(32)))
#define ALIGN128 __attribute__ ((aligned(128)))
#define ALIGN256 __attribute__ ((aligned(256)))

#define predict(dt, x, v, aby2, jby6) \
  ((x) + (v) * (dt) + (aby2) * (dt) * (dt) + (jby6) * (dt) * (dt) * (dt))

static float ALIGN32 three[NVECS] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
static float ALIGN32 threefourth[NVECS] = {0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75};
static float ALIGN32 flag[NVECS] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

static double time;
static int nblen[NPIPES];
static int nbl[NPIPES][MAXLEN];
static int nblerror;

static struct Ptcl_Mem{
  double pos[3];
  double vel[3];
  double acc[3];
  double jrk[3];
  double mss;
  double tim;
  int    idx;
  int    pad[3];
} ptcl_mem[JMEMSIZE] ALIGN128;

typedef struct Pred_Mem * pPred_Mem;
static struct Pred_Mem{
  double xpos[NVECD], ypos[NVECD], zpos[NVECD];
  float  indx[NVECS], mass[NVECS];
  float  xvel[NVECS], yvel[NVECS], zvel[NVECS];
} pred_mem[JMEMSIZE] ALIGN256;

typedef struct NeighbourList * pNeighbourList;
static struct NeighbourList{
  float flag[NVECS];
} (*neighbour)[JMEMSIZE];

typedef struct Iparticle * pIparticle;
struct Iparticle{
  double xpos0[NVECD], xpos1[NVECD];
  double ypos0[NVECD], ypos1[NVECD];
  double zpos0[NVECD], zpos1[NVECD];
  float  xvel01[NVECS], yvel01[NVECS], zvel01[NVECS];
  float  id01[NVECS], veps2[NVECS];
  double xacc[NVECD], yacc[NVECD], zacc[NVECD], pot[NVECD];
  float  xjrk[NVECS], yjrk[NVECS], zjrk[NVECS];
  float  rmin2[NVECS], in[NVECS];
  float  hinv[NVECS];
};
#define NVAR_IP 21

void avx_debugfunc(void)
{
  int j;

  for(j = 0; j < 1024; j++){
    printf("%4d %+.13E %+.13E\n", j, ptcl_mem[j].acc[0], ptcl_mem[j].jrk[0]);
  }

  return;
}

void avx_open(int nthread)
{
  int ret;
  
  ret = posix_memalign((void **)&neighbour, 32, sizeof(struct NeighbourList) * JMEMSIZE * nthread);
  assert(ret == 0);
  
  return;
}

void avx_close(void)
{
  free(neighbour);
  return;
}

void avx_set_j_particle(int padr, int pidx, double tim, double mss,
			double *pos, double *vel, double *acc, double *jrk)
{
  ptcl_mem[padr].pos[0] = pos[0];
  ptcl_mem[padr].pos[1] = pos[1];
  ptcl_mem[padr].pos[2] = pos[2];
  ptcl_mem[padr].vel[0] = vel[0];
  ptcl_mem[padr].vel[1] = vel[1];
  ptcl_mem[padr].vel[2] = vel[2];
  ptcl_mem[padr].acc[0] = acc[0];
  ptcl_mem[padr].acc[1] = acc[1];
  ptcl_mem[padr].acc[2] = acc[2];
  ptcl_mem[padr].jrk[0] = jrk[0];
  ptcl_mem[padr].jrk[1] = jrk[1];
  ptcl_mem[padr].jrk[2] = jrk[2];
  ptcl_mem[padr].mss    = mss;
  ptcl_mem[padr].tim    = tim;
  ptcl_mem[padr].idx    = pidx;

  return;
}

void avx_set_ti(double tim)
{
  time = tim;
  return;
}

void avx_initialize_neighbourlist(void)
{
  nblerror = 0;

  return;
}

int avx_get_neighbourlist_error(void)
{
  return nblerror;
}

int avx_get_neighbourlist(int ipipe, int maxlen, int *nblenfunc, int *nblfunc)
{
  int j;

  if(nblen[ipipe] > maxlen){
    return 1;
  }else{
    *nblenfunc = nblen[ipipe];
    for(j = 0; j < nblen[ipipe]; j++)
      nblfunc[j] = nbl[ipipe][j];
    return 0;
  }
}

void avx_predict_j_particle(int nj)
{
  int j, jmod;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(j = 0; j < nj; j += JPARA){
    int jmem = j  / JPARA;
    int jj;
    for(jj = 0; jj < JPARA; jj++){
      int jadr  = j  + jj;
      int j2    = jj + JPARA;
      double dt = time - ptcl_mem[jadr].tim;

      pred_mem[jmem].xpos[jj] = predict(dt, ptcl_mem[jadr].pos[0], ptcl_mem[jadr].vel[0], ptcl_mem[jadr].acc[0], ptcl_mem[jadr].jrk[0]);
      pred_mem[jmem].ypos[jj] = predict(dt, ptcl_mem[jadr].pos[1], ptcl_mem[jadr].vel[1], ptcl_mem[jadr].acc[1], ptcl_mem[jadr].jrk[1]);
      pred_mem[jmem].zpos[jj] = predict(dt, ptcl_mem[jadr].pos[2], ptcl_mem[jadr].vel[2], ptcl_mem[jadr].acc[2], ptcl_mem[jadr].jrk[2]);

      pred_mem[jmem].indx[jj] = pred_mem[jmem].indx[j2] = (float)ptcl_mem[jadr].idx;
      pred_mem[jmem].mass[jj] = pred_mem[jmem].mass[j2] = (float)ptcl_mem[jadr].mss;

      pred_mem[jmem].xvel[jj] = predict(dt, ptcl_mem[jadr].vel[0], SEC * ptcl_mem[jadr].acc[0], THD * ptcl_mem[jadr].jrk[0], 0.0);
      pred_mem[jmem].xvel[j2] = pred_mem[jmem].xvel[jj];
      pred_mem[jmem].yvel[jj] = predict(dt, ptcl_mem[jadr].vel[1], SEC * ptcl_mem[jadr].acc[1], THD * ptcl_mem[jadr].jrk[1], 0.0);
      pred_mem[jmem].yvel[j2] = pred_mem[jmem].yvel[jj];
      pred_mem[jmem].zvel[jj] = predict(dt, ptcl_mem[jadr].vel[2], SEC * ptcl_mem[jadr].acc[2], THD * ptcl_mem[jadr].jrk[2], 0.0);
      pred_mem[jmem].zvel[j2] = pred_mem[jmem].zvel[jj];
    }   
  }

  if((jmod = nj % JPARA) != 0){
    int jj;
    int jmem = nj / JPARA;
    for(jj = JPARA - 1; jj >= jmod; jj--){
      pred_mem[jmem].xpos[jj]   = 1e3;
      pred_mem[jmem].ypos[jj]   = 1e3;
      pred_mem[jmem].zpos[jj]   = 1e3;
      pred_mem[jmem].mass[jj]   = 0.0;
      pred_mem[jmem].mass[jj+4] = 0.0;
    }
  }

  return;
}

void gravity_kernels(int nj, double eps2, pPrdPosVel posvel, pNewAccJrk accjerk)
{
  int i, j, jj;
  float id;
  double  pxp, pyp, pzp;
  double  vxp, vyp, vzp;
  double  ax, ay, az;
  double  jx, jy, jz;
  double  pot;
  double  r2, rinv, rinv2, rinv3, rv;
  double  dpx, dpy, dpz;
  double  dvx, dvy, dvz;
  pPred_Mem jptr;

  for(i = 0; i < IPARA; i++){
    id  = posvel[i].id;
    pxp = posvel[i].xpos;
    pyp = posvel[i].ypos;
    pzp = posvel[i].zpos;
    vxp = posvel[i].xvel;
    vyp = posvel[i].yvel;
    vzp = posvel[i].zvel;
    ax = ay = az = jx = jy = jz = pot = 0.0;
    for(j = 0, jptr = pred_mem; j < nj; j += JPARA, jptr++){      
      for(jj = 0; jj < JPARA; jj++){
	if(jptr->indx[jj] == id)
	  continue;
	dpx = jptr->xpos[jj] - pxp;
	dpy = jptr->ypos[jj] - pyp;
	dpz = jptr->zpos[jj] - pzp;
	dvx = jptr->xvel[jj] - vxp;
	dvy = jptr->yvel[jj] - vyp;
	dvz = jptr->zvel[jj] - vzp;
	r2  = dpx * dpx + dpy * dpy + dpz * dpz + eps2;
	rv  = dpx * dvx + dpy * dvy + dpz * dvz;

	rinv2 = 1.0 / r2;
	rinv  = sqrt(rinv2);
	rv   *= rinv2 * 3.0;
	rinv *= jptr->mass[jj];
	rinv3 = rinv * rinv2;

	pot -= rinv;
	dpx *= rinv3; ax += dpx;
	dpy *= rinv3; ay += dpy;
	dpz *= rinv3; az += dpz;
	dvx *= rinv3; jx += dvx;
	dvy *= rinv3; jy += dvy;
	dvz *= rinv3; jz += dvz;
	dpx *= rv;    jx -= dpx;
	dpy *= rv;    jy -= dpy;
	dpz *= rv;    jz -= dpz;
      }
    }
    accjerk[i].xacc = ax;
    accjerk[i].yacc = ay;
    accjerk[i].zacc = az;
    accjerk[i].xjrk = jx;
    accjerk[i].yjrk = jy;
    accjerk[i].zjrk = jz;
    accjerk[i].pot  = pot;
  }
  return;
}

void gravity_kernel(int nj, pPrdPosVel posvel, pNewAccJrk accjerk)
{
  int ret;
  int j;
  pPred_Mem jptr = pred_mem;
  pIparticle iptr;

  ret = posix_memalign((void **)&iptr, 32, NVAR_IP * 32);
  assert(ret == 0);

  VBROADCASTSD(posvel[0].xpos, YMM00);
  VBROADCASTSD(posvel[0].ypos, YMM01);
  VBROADCASTSD(posvel[0].zpos, YMM02);

  VBROADCASTSD(posvel[1].xpos, YMM03);
  VBROADCASTSD(posvel[1].ypos, YMM04);
  VBROADCASTSD(posvel[1].zpos, YMM05);

  VBROADCASTSS(posvel[0].xvel, XMM06);
  VBROADCASTSS(posvel[1].xvel, XMM07);
  VMERGE(YMM06, YMM07, YMM06);

  VBROADCASTSS(posvel[0].yvel, XMM08);
  VBROADCASTSS(posvel[1].yvel, XMM09);
  VMERGE(YMM08, YMM09, YMM07);

  VBROADCASTSS(posvel[0].zvel, XMM10);
  VBROADCASTSS(posvel[1].zvel, XMM11);
  VMERGE(YMM10, YMM11, YMM08);

  VBROADCASTSS(posvel[0].id, XMM12);
  VBROADCASTSS(posvel[1].id, XMM13);
  VMERGE(YMM12, YMM13, YMM09);

  VBROADCASTSS(posvel[0].eps2, XMM14);
  VBROADCASTSS(posvel[1].eps2, XMM15);
  VMERGE(YMM14, YMM15, YMM10);

  VSTORPD(YMM00, iptr->xpos0[0]);
  VSTORPD(YMM01, iptr->ypos0[0]);
  VSTORPD(YMM02, iptr->zpos0[0]);
  VSTORPD(YMM03, iptr->xpos1[0]);
  VSTORPD(YMM04, iptr->ypos1[0]);
  VSTORPD(YMM05, iptr->zpos1[0]);
  VSTORPS(YMM06, iptr->xvel01[0]);
  VSTORPS(YMM07, iptr->yvel01[0]);
  VSTORPS(YMM08, iptr->zvel01[0]);
  VSTORPS(YMM09, iptr->id01[0]);
  VSTORPS(YMM10, iptr->veps2[0]);

  VZEROALL;
  for(j = 0; j < nj; j += JPARA, jptr++){ // if nj % 2 != 0 ATARU
    // dx -> YMM03
    VLOADPD(jptr->xpos[0], YMM00);
    VSUBPD_M(iptr->xpos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->xpos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM03);
    // dy -> YMM04
    VLOADPD(jptr->ypos[0], YMM00);
    VSUBPD_M(iptr->ypos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->ypos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM04);
    // dz -> YMM05
    VLOADPD(jptr->zpos[0], YMM00);
    VSUBPD_M(iptr->zpos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->zpos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM05);
    // dr^2
    VLOADPS(iptr->veps2[0], YMM01);
    VFMADDPS(YMM01, YMM03, YMM03);
    VFMADDPS(YMM01, YMM04, YMM04);
    VFMADDPS(YMM01, YMM05, YMM05);
    // - 2 / r -> YMM01
    VRSQRTPS(YMM01, YMM02);
    VMULPS(YMM02, YMM01, YMM01);
    VFMSUB213PS_M(three[0], YMM02, YMM01);
    VMULPS(YMM02, YMM01, YMM01);
    // mask
    VLOADPS(jptr->indx[0], YMM02);
    VLOADPS(iptr->id01[0], YMM00);
    VCMPNEQPS(YMM00, YMM02, YMM02);
    VANDPS(YMM02, YMM01, YMM01);    
    // potential
    VMULPS_M(jptr->mass[0], YMM01, YMM02);
    VCVTPS2PD(XMM02, YMM00);
    VUP2LOW(YMM02, XMM06);
    VCVTPS2PD(XMM06, YMM06);
    VHADDPD(YMM06, YMM00, YMM07);
    VADDPD(YMM07, YMM09, YMM09);
    // dvx, dvy, dvz (vj - vi)
    VLOADPS(jptr->xvel[0], YMM06);
    VSUBPS_M(iptr->xvel01[0], YMM06, YMM06);
    VLOADPS(jptr->yvel[0], YMM07);
    VSUBPS_M(iptr->yvel01[0], YMM07, YMM07);
    VLOADPS(jptr->zvel[0], YMM08);
    VSUBPS_M(iptr->zvel01[0], YMM08, YMM08);
    // xv -> YMM00
    VMULPS(YMM03, YMM06, YMM00);
    VFMADDPS(YMM00, YMM04, YMM07);
    VFMADDPS(YMM00, YMM05, YMM08);
    // YMM00: 3.0 * xv / r^2, YMM02: - m / r^3
    VMULPS_M(jptr->mass[0], YMM01, YMM02);
    VMULPS(YMM01, YMM01, YMM01);
    VMULPS(YMM01, YMM00, YMM00);
    VMULPS(YMM01, YMM02, YMM02);
    VMULPS_M(threefourth[0], YMM00, YMM00);
    // prefetch
    PREFETCH((jptr+1)->xpos[0]);
    PREFETCH((jptr+1)->zpos[0]);
    PREFETCH((jptr+1)->mass[0]);
    PREFETCH((jptr+1)->yvel[0]);
    // jx1, jy1, jz1
    VFMADDPS(YMM13, YMM02, YMM06);
    VFMADDPS(YMM14, YMM02, YMM07);
    VFMADDPS(YMM15, YMM02, YMM08);
    // ax
    VMULPS(YMM02, YMM03, YMM03);
    VCVTPS2PD(XMM03, YMM06);
    VUP2LOW(YMM03, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM10, YMM10);
    // ay
    VMULPS(YMM02, YMM04, YMM04);
    VCVTPS2PD(XMM04, YMM06);
    VUP2LOW(YMM04, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM11, YMM11);
    // az
    VMULPS(YMM02, YMM05, YMM05);
    VCVTPS2PD(XMM05, YMM06);
    VUP2LOW(YMM05, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM12, YMM12);
    // jx2, jy2, jz2
    VFNMADDPS(YMM13, YMM00, YMM03);
    VFNMADDPS(YMM14, YMM00, YMM04);
    VFNMADDPS(YMM15, YMM00, YMM05);
  }

  VSTORPD(YMM09, iptr->pot[0]);
  VSTORPD(YMM10, iptr->xacc[0]);
  VSTORPD(YMM11, iptr->yacc[0]);
  VSTORPD(YMM12, iptr->zacc[0]);
  VSTORPS(YMM13, iptr->xjrk[0]);
  VSTORPS(YMM14, iptr->yjrk[0]);
  VSTORPS(YMM15, iptr->zjrk[0]);

  VZEROUPPER;

  accjerk[0].xacc = iptr->xacc[0] + iptr->xacc[2];
  accjerk[0].yacc = iptr->yacc[0] + iptr->yacc[2];
  accjerk[0].zacc = iptr->zacc[0] + iptr->zacc[2];
  accjerk[0].pot  = iptr->pot[0]  + iptr->pot[2];
  accjerk[0].xjrk = iptr->xjrk[0] + iptr->xjrk[1] + iptr->xjrk[2] + iptr->xjrk[3];
  accjerk[0].yjrk = iptr->yjrk[0] + iptr->yjrk[1] + iptr->yjrk[2] + iptr->yjrk[3];
  accjerk[0].zjrk = iptr->zjrk[0] + iptr->zjrk[1] + iptr->zjrk[2] + iptr->zjrk[3];

  accjerk[1].xacc = iptr->xacc[1] + iptr->xacc[3];
  accjerk[1].yacc = iptr->yacc[1] + iptr->yacc[3];
  accjerk[1].zacc = iptr->zacc[1] + iptr->zacc[3];
  accjerk[1].pot  = iptr->pot[1]  + iptr->pot[3];
  accjerk[1].xjrk = iptr->xjrk[4] + iptr->xjrk[5] + iptr->xjrk[6] + iptr->xjrk[7];
  accjerk[1].yjrk = iptr->yjrk[4] + iptr->yjrk[5] + iptr->yjrk[6] + iptr->yjrk[7];
  accjerk[1].zjrk = iptr->zjrk[4] + iptr->zjrk[5] + iptr->zjrk[6] + iptr->zjrk[7];

  free(iptr);

  return;
}

void gravity_kernel2(int nj, pPrdPosVel posvel, pNewAccJrk accjerk)
{
  int ret;
  int j;
  double true_rmin2;
  pPred_Mem jptr = pred_mem;
  pIparticle iptr;
  float ten = 10.0, minusone = -1.0;

  ret = posix_memalign((void **)&iptr, 32, NVAR_IP * 32);
  assert(ret == 0);

  VBROADCASTSD(posvel[0].xpos, YMM00);
  VBROADCASTSD(posvel[0].ypos, YMM01);
  VBROADCASTSD(posvel[0].zpos, YMM02);

  VBROADCASTSD(posvel[1].xpos, YMM03);
  VBROADCASTSD(posvel[1].ypos, YMM04);
  VBROADCASTSD(posvel[1].zpos, YMM05);

  VBROADCASTSS(posvel[0].xvel, XMM06);
  VBROADCASTSS(posvel[1].xvel, XMM07);
  VMERGE(YMM06, YMM07, YMM06);

  VBROADCASTSS(posvel[0].yvel, XMM08);
  VBROADCASTSS(posvel[1].yvel, XMM09);
  VMERGE(YMM08, YMM09, YMM07);

  VBROADCASTSS(posvel[0].zvel, XMM10);
  VBROADCASTSS(posvel[1].zvel, XMM11);
  VMERGE(YMM10, YMM11, YMM08);

  VBROADCASTSS(posvel[0].id, XMM12);
  VBROADCASTSS(posvel[1].id, XMM13);
  VMERGE(YMM12, YMM13, YMM09);

  VBROADCASTSS(posvel[0].eps2, XMM14);
  VBROADCASTSS(posvel[1].eps2, XMM15);
  VMERGE(YMM14, YMM15, YMM10);

  VBROADCASTSS(ten, YMM11);
  VBROADCASTSS(minusone, YMM12);

  VSTORPD(YMM00, iptr->xpos0[0]);
  VSTORPD(YMM01, iptr->ypos0[0]);
  VSTORPD(YMM02, iptr->zpos0[0]);
  VSTORPD(YMM03, iptr->xpos1[0]);
  VSTORPD(YMM04, iptr->ypos1[0]);
  VSTORPD(YMM05, iptr->zpos1[0]);
  VSTORPS(YMM06, iptr->xvel01[0]);
  VSTORPS(YMM07, iptr->yvel01[0]);
  VSTORPS(YMM08, iptr->zvel01[0]);
  VSTORPS(YMM09, iptr->id01[0]);
  VSTORPS(YMM10, iptr->veps2[0]);
  VSTORPS(YMM11, iptr->rmin2[0]);
  VSTORPS(YMM12, iptr->in[0]);

  VZEROALL;
  for(j = 0; j < nj; j += JPARA, jptr++){ // if nj % 2 != 0 ATARU
    // dx -> YMM03
    VLOADPD(jptr->xpos[0], YMM00);
    VSUBPD_M(iptr->xpos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->xpos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM03);
    // dy -> YMM04
    VLOADPD(jptr->ypos[0], YMM00);
    VSUBPD_M(iptr->ypos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->ypos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM04);
    // dz -> YMM05
    VLOADPD(jptr->zpos[0], YMM00);
    VSUBPD_M(iptr->zpos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->zpos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM05);
    // dr^2
    VLOADPS(iptr->veps2[0], YMM01);
    VFMADDPS(YMM01, YMM03, YMM03);
    VFMADDPS(YMM01, YMM04, YMM04);
    VFMADDPS(YMM01, YMM05, YMM05);
    // - 2 / r -> YMM01
    VRSQRTPS(YMM01, YMM02);
    VMULPS(YMM02, YMM01, YMM01);
    VFMSUB213PS_M(three[0], YMM02, YMM01);
    VMULPS(YMM02, YMM01, YMM01);
    // mask
    VLOADPS(jptr->indx[0], YMM02);
    VLOADPS(iptr->id01[0], YMM00);
    VCMPNEQPS(YMM00, YMM02, YMM02);
    VANDPS(YMM02, YMM01, YMM01);    
    // nearest neighbour (free: YMM00, YMM02, YMM06, YMM07, YMM08)
    VLOADPS(iptr->rmin2[0], YMM00);
    VMINPS(YMM01, YMM00, YMM02);
    VSTORPS(YMM02, iptr->rmin2[0]);
    VCMPPS(YMM01, YMM00, YMM02, GT);
    VLOADPS(jptr->indx[0], YMM06);
    VANDPS(YMM02, YMM06, YMM07);
    VCMPPS(YMM01, YMM00, YMM08, LE);
    VANDPS_M(iptr->in[0], YMM08, YMM08);
    VADDPS(YMM08, YMM07, YMM07);
    VSTORPS(YMM07, iptr->in[0]);
    // potential
    VMULPS_M(jptr->mass[0], YMM01, YMM02);
    VCVTPS2PD(XMM02, YMM00);
    VUP2LOW(YMM02, XMM06);
    VCVTPS2PD(XMM06, YMM06);
    VHADDPD(YMM06, YMM00, YMM07);
    VADDPD(YMM07, YMM09, YMM09);
    // dvx, dvy, dvz (vj - vi)
    VLOADPS(jptr->xvel[0], YMM06);
    VSUBPS_M(iptr->xvel01[0], YMM06, YMM06);
    VLOADPS(jptr->yvel[0], YMM07);
    VSUBPS_M(iptr->yvel01[0], YMM07, YMM07);
    VLOADPS(jptr->zvel[0], YMM08);
    VSUBPS_M(iptr->zvel01[0], YMM08, YMM08);
    // xv -> YMM00
    VMULPS(YMM03, YMM06, YMM00);
    VFMADDPS(YMM00, YMM04, YMM07);
    VFMADDPS(YMM00, YMM05, YMM08);
    // YMM00: 3.0 * xv / r^2, YMM02: - m / r^3
    VMULPS_M(jptr->mass[0], YMM01, YMM02);
    VMULPS(YMM01, YMM01, YMM01);
    VMULPS(YMM01, YMM00, YMM00);
    VMULPS(YMM01, YMM02, YMM02);
    VMULPS_M(threefourth[0], YMM00, YMM00);
    // prefetch
    PREFETCH((jptr+1)->xpos[0]);
    PREFETCH((jptr+1)->zpos[0]);
    PREFETCH((jptr+1)->mass[0]);
    PREFETCH((jptr+1)->yvel[0]);
    // jx1, jy1, jz1
    VFMADDPS(YMM13, YMM02, YMM06);
    VFMADDPS(YMM14, YMM02, YMM07);
    VFMADDPS(YMM15, YMM02, YMM08);
    // ax
    VMULPS(YMM02, YMM03, YMM03);
    VCVTPS2PD(XMM03, YMM06);
    VUP2LOW(YMM03, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM10, YMM10);
    // ay
    VMULPS(YMM02, YMM04, YMM04);
    VCVTPS2PD(XMM04, YMM06);
    VUP2LOW(YMM04, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM11, YMM11);
    // az
    VMULPS(YMM02, YMM05, YMM05);
    VCVTPS2PD(XMM05, YMM06);
    VUP2LOW(YMM05, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM12, YMM12);
    // jx2, jy2, jz2
    VFNMADDPS(YMM13, YMM00, YMM03);
    VFNMADDPS(YMM14, YMM00, YMM04);
    VFNMADDPS(YMM15, YMM00, YMM05);
  }

  VSTORPD(YMM09, iptr->pot[0]);
  VSTORPD(YMM10, iptr->xacc[0]);
  VSTORPD(YMM11, iptr->yacc[0]);
  VSTORPD(YMM12, iptr->zacc[0]);
  VSTORPS(YMM13, iptr->xjrk[0]);
  VSTORPS(YMM14, iptr->yjrk[0]);
  VSTORPS(YMM15, iptr->zjrk[0]);

  accjerk[0].xacc = iptr->xacc[0] + iptr->xacc[2];
  accjerk[0].yacc = iptr->yacc[0] + iptr->yacc[2];
  accjerk[0].zacc = iptr->zacc[0] + iptr->zacc[2];
  accjerk[0].pot  = iptr->pot[0]  + iptr->pot[2];
  accjerk[0].xjrk = iptr->xjrk[0] + iptr->xjrk[1] + iptr->xjrk[2] + iptr->xjrk[3];
  accjerk[0].yjrk = iptr->yjrk[0] + iptr->yjrk[1] + iptr->yjrk[2] + iptr->yjrk[3];
  accjerk[0].zjrk = iptr->zjrk[0] + iptr->zjrk[1] + iptr->zjrk[2] + iptr->zjrk[3];
  for(true_rmin2 = 1e30, j = 0; j < JPARA; j++){
    if(iptr->rmin2[j] < true_rmin2){
      true_rmin2    = iptr->rmin2[j];
      accjerk[0].rnnb = - 2.0 / true_rmin2;
      accjerk[0].nnb  = (int)iptr->in[j];
    }
  }

  accjerk[1].xacc = iptr->xacc[1] + iptr->xacc[3];
  accjerk[1].yacc = iptr->yacc[1] + iptr->yacc[3];
  accjerk[1].zacc = iptr->zacc[1] + iptr->zacc[3];
  accjerk[1].pot  = iptr->pot[1]  + iptr->pot[3];
  accjerk[1].xjrk = iptr->xjrk[4] + iptr->xjrk[5] + iptr->xjrk[6] + iptr->xjrk[7];
  accjerk[1].yjrk = iptr->yjrk[4] + iptr->yjrk[5] + iptr->yjrk[6] + iptr->yjrk[7];
  accjerk[1].zjrk = iptr->zjrk[4] + iptr->zjrk[5] + iptr->zjrk[6] + iptr->zjrk[7];
  for(true_rmin2 = 1e30, j = 4; j < 4 + JPARA; j++){
    if(iptr->rmin2[j] < true_rmin2){
      true_rmin2    = iptr->rmin2[j];
      accjerk[1].rnnb = - 2.0 / true_rmin2;
      accjerk[1].nnb = (int)iptr->in[j];
    }
  }

  free(iptr);

  return;
}

void gravity_kerneln(int nj, pPrdPosVel posvel, pNewAccJrk accjerk, int i, int ithread)
{
  int ret;
  int j;
  float hinv0, hinv1;
  pPred_Mem jptr = pred_mem;
  pIparticle iptr;
  pNeighbourList nbptr, nbptr0 = neighbour[ithread];

  if(posvel[0].h2 == 0.0)
    hinv0 = - 1e10;
  else
    hinv0 = - 2.0 / sqrt(posvel[0].h2);
  if(posvel[1].h2 == 0.0)
    hinv1 = - 1e10;
  else
    hinv1 = - 2.0 / sqrt(posvel[1].h2);

  ret = posix_memalign((void **)&iptr, 32, NVAR_IP * 32);
  assert(ret == 0);

  VBROADCASTSD(posvel[0].xpos, YMM00);
  VBROADCASTSD(posvel[0].ypos, YMM01);
  VBROADCASTSD(posvel[0].zpos, YMM02);

  VBROADCASTSD(posvel[1].xpos, YMM03);
  VBROADCASTSD(posvel[1].ypos, YMM04);
  VBROADCASTSD(posvel[1].zpos, YMM05);

  VBROADCASTSS(posvel[0].xvel, XMM06);
  VBROADCASTSS(posvel[1].xvel, XMM07);
  VMERGE(YMM06, YMM07, YMM06);

  VBROADCASTSS(posvel[0].yvel, XMM08);
  VBROADCASTSS(posvel[1].yvel, XMM09);
  VMERGE(YMM08, YMM09, YMM07);

  VBROADCASTSS(posvel[0].zvel, XMM10);
  VBROADCASTSS(posvel[1].zvel, XMM11);
  VMERGE(YMM10, YMM11, YMM08);

  VBROADCASTSS(posvel[0].id, XMM12);
  VBROADCASTSS(posvel[1].id, XMM13);
  VMERGE(YMM12, YMM13, YMM09);

  VBROADCASTSS(posvel[0].eps2, XMM14);
  VBROADCASTSS(posvel[1].eps2, XMM15);
  VMERGE(YMM14, YMM15, YMM10);

  VBROADCASTSS(hinv0, XMM11);
  VBROADCASTSS(hinv1, XMM12);
  VMERGE(YMM11, YMM12, YMM11);

  VSTORPD(YMM00, iptr->xpos0[0]);
  VSTORPD(YMM01, iptr->ypos0[0]);
  VSTORPD(YMM02, iptr->zpos0[0]);
  VSTORPD(YMM03, iptr->xpos1[0]);
  VSTORPD(YMM04, iptr->ypos1[0]);
  VSTORPD(YMM05, iptr->zpos1[0]);
  VSTORPS(YMM06, iptr->xvel01[0]);
  VSTORPS(YMM07, iptr->yvel01[0]);
  VSTORPS(YMM08, iptr->zvel01[0]);
  VSTORPS(YMM09, iptr->id01[0]);
  VSTORPS(YMM10, iptr->veps2[0]);
  VSTORPS(YMM11, iptr->hinv[0]);

  VZEROALL;
  for(j = 0, nbptr = nbptr0; j < nj; j += JPARA, jptr++, nbptr++){ // if nj % 2 != 0 ATARU
    // dx -> YMM03
    VLOADPD(jptr->xpos[0], YMM00);
    VSUBPD_M(iptr->xpos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->xpos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM03);
    // dy -> YMM04
    VLOADPD(jptr->ypos[0], YMM00);
    VSUBPD_M(iptr->ypos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->ypos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM04);
    // dz -> YMM05
    VLOADPD(jptr->zpos[0], YMM00);
    VSUBPD_M(iptr->zpos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->zpos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM05);
    // dr^2
    VLOADPS(iptr->veps2[0], YMM01);
    VFMADDPS(YMM01, YMM03, YMM03);
    VFMADDPS(YMM01, YMM04, YMM04);
    VFMADDPS(YMM01, YMM05, YMM05);
    // - 2 / r -> YMM01
    VRSQRTPS(YMM01, YMM02);
    VMULPS(YMM02, YMM01, YMM01);
    VFMSUB213PS_M(three[0], YMM02, YMM01);
    VMULPS(YMM02, YMM01, YMM01);
    // mask
    VLOADPS(jptr->indx[0], YMM02);
    VLOADPS(iptr->id01[0], YMM00);
    VCMPNEQPS(YMM00, YMM02, YMM02);
    VANDPS(YMM02, YMM01, YMM01);    
    // neighbour list
    VLOADPS(iptr->hinv[0], YMM00);
    VCMPPS(YMM00, YMM01, YMM00, LE);
    VLOADPS(flag[0], YMM02);
    VANDPS(YMM02, YMM00, YMM00);
    VSTORPS(YMM00, nbptr->flag[0]);
    // potential
    VMULPS_M(jptr->mass[0], YMM01, YMM02);
    VCVTPS2PD(XMM02, YMM00);
    VUP2LOW(YMM02, XMM06);
    VCVTPS2PD(XMM06, YMM06);
    VHADDPD(YMM06, YMM00, YMM07);
    VADDPD(YMM07, YMM09, YMM09);
    // dvx, dvy, dvz (vj - vi)
    VLOADPS(jptr->xvel[0], YMM06);
    VSUBPS_M(iptr->xvel01[0], YMM06, YMM06);
    VLOADPS(jptr->yvel[0], YMM07);
    VSUBPS_M(iptr->yvel01[0], YMM07, YMM07);
    VLOADPS(jptr->zvel[0], YMM08);
    VSUBPS_M(iptr->zvel01[0], YMM08, YMM08);
    // xv -> YMM00
    VMULPS(YMM03, YMM06, YMM00);
    VFMADDPS(YMM00, YMM04, YMM07);
    VFMADDPS(YMM00, YMM05, YMM08);
    // YMM00: 3.0 * xv / r^2, YMM02: - m / r^3
    VMULPS_M(jptr->mass[0], YMM01, YMM02);
    VMULPS(YMM01, YMM01, YMM01);
    VMULPS(YMM01, YMM00, YMM00);
    VMULPS(YMM01, YMM02, YMM02);
    VMULPS_M(threefourth[0], YMM00, YMM00);
    // prefetch
    PREFETCH((jptr+1)->xpos[0]);
    PREFETCH((jptr+1)->zpos[0]);
    PREFETCH((jptr+1)->mass[0]);
    PREFETCH((jptr+1)->yvel[0]);
    // jx1, jy1, jz1
    VFMADDPS(YMM13, YMM02, YMM06);
    VFMADDPS(YMM14, YMM02, YMM07);
    VFMADDPS(YMM15, YMM02, YMM08);
    // ax
    VMULPS(YMM02, YMM03, YMM03);
    VCVTPS2PD(XMM03, YMM06);
    VUP2LOW(YMM03, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM10, YMM10);
    // ay
    VMULPS(YMM02, YMM04, YMM04);
    VCVTPS2PD(XMM04, YMM06);
    VUP2LOW(YMM04, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM11, YMM11);
    // az
    VMULPS(YMM02, YMM05, YMM05);
    VCVTPS2PD(XMM05, YMM06);
    VUP2LOW(YMM05, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM12, YMM12);
    // jx2, jy2, jz2
    VFNMADDPS(YMM13, YMM00, YMM03);
    VFNMADDPS(YMM14, YMM00, YMM04);
    VFNMADDPS(YMM15, YMM00, YMM05);
  }

  VSTORPD(YMM09, iptr->pot[0]);
  VSTORPD(YMM10, iptr->xacc[0]);
  VSTORPD(YMM11, iptr->yacc[0]);
  VSTORPD(YMM12, iptr->zacc[0]);
  VSTORPS(YMM13, iptr->xjrk[0]);
  VSTORPS(YMM14, iptr->yjrk[0]);
  VSTORPS(YMM15, iptr->zjrk[0]);

  accjerk[0].xacc = iptr->xacc[0] + iptr->xacc[2];
  accjerk[0].yacc = iptr->yacc[0] + iptr->yacc[2];
  accjerk[0].zacc = iptr->zacc[0] + iptr->zacc[2];
  accjerk[0].pot  = iptr->pot[0]  + iptr->pot[2];
  accjerk[0].xjrk = iptr->xjrk[0] + iptr->xjrk[1] + iptr->xjrk[2] + iptr->xjrk[3];
  accjerk[0].yjrk = iptr->yjrk[0] + iptr->yjrk[1] + iptr->yjrk[2] + iptr->yjrk[3];
  accjerk[0].zjrk = iptr->zjrk[0] + iptr->zjrk[1] + iptr->zjrk[2] + iptr->zjrk[3];

  accjerk[1].xacc = iptr->xacc[1] + iptr->xacc[3];
  accjerk[1].yacc = iptr->yacc[1] + iptr->yacc[3];
  accjerk[1].zacc = iptr->zacc[1] + iptr->zacc[3];
  accjerk[1].pot  = iptr->pot[1]  + iptr->pot[3];
  accjerk[1].xjrk = iptr->xjrk[4] + iptr->xjrk[5] + iptr->xjrk[6] + iptr->xjrk[7];
  accjerk[1].yjrk = iptr->yjrk[4] + iptr->yjrk[5] + iptr->yjrk[6] + iptr->yjrk[7];
  accjerk[1].zjrk = iptr->zjrk[4] + iptr->zjrk[5] + iptr->zjrk[6] + iptr->zjrk[7];

  int jj;
  int nn0, nn1;
  for(nn0 = nn1 = 0, j = 0, jptr = pred_mem, nbptr = nbptr0; j < nj; j += JPARA, jptr++, nbptr++){
    for(jj = 0; jj < JPARA; jj++)
      if(nbptr->flag[jj] == 1.0){
	nbl[i][nn0] = (int)jptr->indx[jj];
	++nn0;
      }
    for(jj = 4; jj < JPARA + 4; jj++)
      if(nbptr->flag[jj] == 1.0){
	nbl[i+1][nn1] = (int)jptr->indx[jj];
	++nn1;
      }
  }

  if(nn0 > MAXLEN || nn1 > MAXLEN)
    nblerror = 1;

  nblen[i]   = nn0;
  nblen[i+1] = nn1;

  free(iptr);

  return;
}

void gravity_kernel2n(int nj, pPrdPosVel posvel, pNewAccJrk accjerk, int i, int ithread)
{
  int ret;
  int j;
  double true_rmin2;
  float hinv0, hinv1;
  pPred_Mem jptr = pred_mem;
  pIparticle iptr;
  pNeighbourList nbptr, nbptr0 = neighbour[ithread];
  float ten = 10.0, minusone = -1.0;

  if(posvel[0].h2 == 0.0)
    hinv0 = - 1e10;
  else
    hinv0 = - 2.0 / sqrt(posvel[0].h2);
  if(posvel[1].h2 == 0.0)
    hinv1 = - 1e10;
  else
    hinv1 = - 2.0 / sqrt(posvel[1].h2);

  ret = posix_memalign((void **)&iptr, 32, NVAR_IP * 32);
  assert(ret == 0);

  VBROADCASTSD(posvel[0].xpos, YMM00);
  VBROADCASTSD(posvel[0].ypos, YMM01);
  VBROADCASTSD(posvel[0].zpos, YMM02);

  VBROADCASTSD(posvel[1].xpos, YMM03);
  VBROADCASTSD(posvel[1].ypos, YMM04);
  VBROADCASTSD(posvel[1].zpos, YMM05);

  VBROADCASTSS(posvel[0].xvel, XMM06);
  VBROADCASTSS(posvel[1].xvel, XMM07);
  VMERGE(YMM06, YMM07, YMM06);

  VBROADCASTSS(posvel[0].yvel, XMM08);
  VBROADCASTSS(posvel[1].yvel, XMM09);
  VMERGE(YMM08, YMM09, YMM07);

  VBROADCASTSS(posvel[0].zvel, XMM10);
  VBROADCASTSS(posvel[1].zvel, XMM11);
  VMERGE(YMM10, YMM11, YMM08);

  VBROADCASTSS(posvel[0].id, XMM12);
  VBROADCASTSS(posvel[1].id, XMM13);
  VMERGE(YMM12, YMM13, YMM09);

  VBROADCASTSS(posvel[0].eps2, XMM14);
  VBROADCASTSS(posvel[1].eps2, XMM15);
  VMERGE(YMM14, YMM15, YMM10);

  VBROADCASTSS(hinv0, XMM11);
  VBROADCASTSS(hinv1, XMM12);
  VMERGE(YMM11, YMM12, YMM11);

  VBROADCASTSS(ten, YMM12);
  VBROADCASTSS(minusone, YMM13);

  VSTORPD(YMM00, iptr->xpos0[0]);
  VSTORPD(YMM01, iptr->ypos0[0]);
  VSTORPD(YMM02, iptr->zpos0[0]);
  VSTORPD(YMM03, iptr->xpos1[0]);
  VSTORPD(YMM04, iptr->ypos1[0]);
  VSTORPD(YMM05, iptr->zpos1[0]);
  VSTORPS(YMM06, iptr->xvel01[0]);
  VSTORPS(YMM07, iptr->yvel01[0]);
  VSTORPS(YMM08, iptr->zvel01[0]);
  VSTORPS(YMM09, iptr->id01[0]);
  VSTORPS(YMM10, iptr->veps2[0]);
  VSTORPS(YMM11, iptr->hinv[0]);
  VSTORPS(YMM12, iptr->rmin2[0]);
  VSTORPS(YMM13, iptr->in[0]);

  VZEROALL;
  for(j = 0, nbptr = nbptr0; j < nj; j += JPARA, jptr++, nbptr++){ // if nj % 2 != 0 ATARU
    // dx -> YMM03
    VLOADPD(jptr->xpos[0], YMM00);
    VSUBPD_M(iptr->xpos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->xpos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM03);
    // dy -> YMM04
    VLOADPD(jptr->ypos[0], YMM00);
    VSUBPD_M(iptr->ypos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->ypos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM04);
    // dz -> YMM05
    VLOADPD(jptr->zpos[0], YMM00);
    VSUBPD_M(iptr->zpos0[0], YMM00, YMM01);
    VCVTPD2PS(YMM01, XMM01);
    VSUBPD_M(iptr->zpos1[0], YMM00, YMM02);
    VCVTPD2PS(YMM02, XMM02);
    VMERGE(YMM01, YMM02, YMM05);
    // dr^2
    VLOADPS(iptr->veps2[0], YMM01);
    VFMADDPS(YMM01, YMM03, YMM03);
    VFMADDPS(YMM01, YMM04, YMM04);
    VFMADDPS(YMM01, YMM05, YMM05);
    // - 2 / r -> YMM01
    VRSQRTPS(YMM01, YMM02);
    VMULPS(YMM02, YMM01, YMM01);
    VFMSUB213PS_M(three[0], YMM02, YMM01);
    VMULPS(YMM02, YMM01, YMM01);
    // mask
    VLOADPS(jptr->indx[0], YMM02);
    VLOADPS(iptr->id01[0], YMM00);
    VCMPNEQPS(YMM00, YMM02, YMM02);
    VANDPS(YMM02, YMM01, YMM01);    
    // nearest neighbour (free: YMM00, YMM02, YMM06, YMM07, YMM08)
    VLOADPS(iptr->rmin2[0], YMM00);
    VMINPS(YMM01, YMM00, YMM02);
    VSTORPS(YMM02, iptr->rmin2[0]);
    VCMPPS(YMM01, YMM00, YMM02, GT);
    VLOADPS(jptr->indx[0], YMM06);
    VANDPS(YMM02, YMM06, YMM07);
    VCMPPS(YMM01, YMM00, YMM08, LE);
    VANDPS_M(iptr->in[0], YMM08, YMM08);
    VADDPS(YMM08, YMM07, YMM07);
    VSTORPS(YMM07, iptr->in[0]);
    // neighbour list
    VLOADPS(iptr->hinv[0], YMM00);
    VCMPPS(YMM00, YMM01, YMM00, LE);
    VLOADPS(flag[0], YMM02);
    VANDPS(YMM02, YMM00, YMM00);
    VSTORPS(YMM00, nbptr->flag[0]);
    // potential
    VMULPS_M(jptr->mass[0], YMM01, YMM02);
    VCVTPS2PD(XMM02, YMM00);
    VUP2LOW(YMM02, XMM06);
    VCVTPS2PD(XMM06, YMM06);
    VHADDPD(YMM06, YMM00, YMM07);
    VADDPD(YMM07, YMM09, YMM09);
    // dvx, dvy, dvz (vj - vi)
    VLOADPS(jptr->xvel[0], YMM06);
    VSUBPS_M(iptr->xvel01[0], YMM06, YMM06);
    VLOADPS(jptr->yvel[0], YMM07);
    VSUBPS_M(iptr->yvel01[0], YMM07, YMM07);
    VLOADPS(jptr->zvel[0], YMM08);
    VSUBPS_M(iptr->zvel01[0], YMM08, YMM08);
    // xv -> YMM00
    VMULPS(YMM03, YMM06, YMM00);
    VFMADDPS(YMM00, YMM04, YMM07);
    VFMADDPS(YMM00, YMM05, YMM08);
    // YMM00: 3.0 * xv / r^2, YMM02: - m / r^3
    VMULPS_M(jptr->mass[0], YMM01, YMM02);
    VMULPS(YMM01, YMM01, YMM01);
    VMULPS(YMM01, YMM00, YMM00);
    VMULPS(YMM01, YMM02, YMM02);
    VMULPS_M(threefourth[0], YMM00, YMM00);
    // prefetch
    PREFETCH((jptr+1)->xpos[0]);
    PREFETCH((jptr+1)->zpos[0]);
    PREFETCH((jptr+1)->mass[0]);
    PREFETCH((jptr+1)->yvel[0]);
    // jx1, jy1, jz1
    VFMADDPS(YMM13, YMM02, YMM06);
    VFMADDPS(YMM14, YMM02, YMM07);
    VFMADDPS(YMM15, YMM02, YMM08);
    // ax
    VMULPS(YMM02, YMM03, YMM03);
    VCVTPS2PD(XMM03, YMM06);
    VUP2LOW(YMM03, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM10, YMM10);
    // ay
    VMULPS(YMM02, YMM04, YMM04);
    VCVTPS2PD(XMM04, YMM06);
    VUP2LOW(YMM04, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM11, YMM11);
    // az
    VMULPS(YMM02, YMM05, YMM05);
    VCVTPS2PD(XMM05, YMM06);
    VUP2LOW(YMM05, XMM07);
    VCVTPS2PD(XMM07, YMM07);
    VHADDPD(YMM07, YMM06, YMM06);
    VADDPD(YMM06, YMM12, YMM12);
    // jx2, jy2, jz2
    VFNMADDPS(YMM13, YMM00, YMM03);
    VFNMADDPS(YMM14, YMM00, YMM04);
    VFNMADDPS(YMM15, YMM00, YMM05);
  }

  VSTORPD(YMM09, iptr->pot[0]);
  VSTORPD(YMM10, iptr->xacc[0]);
  VSTORPD(YMM11, iptr->yacc[0]);
  VSTORPD(YMM12, iptr->zacc[0]);
  VSTORPS(YMM13, iptr->xjrk[0]);
  VSTORPS(YMM14, iptr->yjrk[0]);
  VSTORPS(YMM15, iptr->zjrk[0]);

  accjerk[0].xacc = iptr->xacc[0] + iptr->xacc[2];
  accjerk[0].yacc = iptr->yacc[0] + iptr->yacc[2];
  accjerk[0].zacc = iptr->zacc[0] + iptr->zacc[2];
  accjerk[0].pot  = iptr->pot[0]  + iptr->pot[2];
  accjerk[0].xjrk = iptr->xjrk[0] + iptr->xjrk[1] + iptr->xjrk[2] + iptr->xjrk[3];
  accjerk[0].yjrk = iptr->yjrk[0] + iptr->yjrk[1] + iptr->yjrk[2] + iptr->yjrk[3];
  accjerk[0].zjrk = iptr->zjrk[0] + iptr->zjrk[1] + iptr->zjrk[2] + iptr->zjrk[3];
  for(true_rmin2 = 1e30, j = 0; j < JPARA; j++){
    if(iptr->rmin2[j] < true_rmin2){
      true_rmin2    = iptr->rmin2[j];
      accjerk[0].rnnb = - 2.0 / true_rmin2;
      accjerk[0].nnb = (int)iptr->in[j];
    }
  }

  accjerk[1].xacc = iptr->xacc[1] + iptr->xacc[3];
  accjerk[1].yacc = iptr->yacc[1] + iptr->yacc[3];
  accjerk[1].zacc = iptr->zacc[1] + iptr->zacc[3];
  accjerk[1].pot  = iptr->pot[1]  + iptr->pot[3];
  accjerk[1].xjrk = iptr->xjrk[4] + iptr->xjrk[5] + iptr->xjrk[6] + iptr->xjrk[7];
  accjerk[1].yjrk = iptr->yjrk[4] + iptr->yjrk[5] + iptr->yjrk[6] + iptr->yjrk[7];
  accjerk[1].zjrk = iptr->zjrk[4] + iptr->zjrk[5] + iptr->zjrk[6] + iptr->zjrk[7];
  for(true_rmin2 = 1e30, j = 4; j < 4 + JPARA; j++){
    if(iptr->rmin2[j] < true_rmin2){
      true_rmin2    = iptr->rmin2[j];
      accjerk[1].rnnb = - 2.0 / true_rmin2;
      accjerk[1].nnb = (int)iptr->in[j];
    }
  }

  int jj;
  int nn0, nn1;
  for(nn0 = nn1 = 0, j = 0, jptr = pred_mem, nbptr = nbptr0; j < nj; j += JPARA, jptr++, nbptr++){
    for(jj = 0; jj < JPARA; jj++)
      if(nbptr->flag[jj] == 1.0){
	nbl[i][nn0] = (int)jptr->indx[jj];
	++nn0;
      }
    for(jj = 4; jj < JPARA + 4; jj++)
      if(nbptr->flag[jj] == 1.0){
	nbl[i+1][nn1] = (int)jptr->indx[jj];
	++nn1;
      }
  }

  if(nn0 > MAXLEN || nn1 > MAXLEN)
    nblerror = 1;

  nblen[i]   = nn0;
  nblen[i+1] = nn1;

  free(iptr);

  return;
}
