#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "sse.h"
#include "avx.h"
#include "sse_type.h"
#include "pg5_table.h"

#define NUNROLL 2

#define ALIGN32 __attribute__ ((aligned(32)))

#define ZI   BUF3
#define ZI_X BUF3_X

typedef struct ipdata_reg{
  float x[8];
  float y[8];
} Ipdata_reg, *pIpdata_reg;

void gravity_kernel(pIpdata ipdata, pJpdata jp, pFodata fodata, 
		int nj, float fcut[][2], v4sf r2cut, v4sf accscale)
{
  int j;
  int ALIGN32 idx[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  Ipdata_reg ALIGN32 ipdata_reg;
  static v4sf two = {2.0f, 2.0f, 2.0f, 2.0f};

  fcut -= (1<<(30-(23-FRC_BIT)));

  VZEROALL;
  VLOADPS(ipdata->x[0], X2_X);
  VLOADPS(ipdata->y[0], Y2_X);
  VLOADPS(ipdata->z[0], Z2_X);
  VLOADPS(r2cut, R2CUT_X);
  VLOADPS(two, TWO_X);
  VPERM2F128(X2, X2, X2, 0x00);
  VSTORPS(X2, ipdata_reg.x[0]);
  VPERM2F128(Y2, Y2, Y2, 0x00);
  VSTORPS(Y2, ipdata_reg.y[0]);
  VPERM2F128(Z2, Z2, ZI, 0x00);
  VPERM2F128(R2CUT, R2CUT, R2CUT, 0x00);
  VPERM2F128(TWO, TWO, TWO, 0x00);

  VLOADPS(*jp, MJ);
  jp += 2;
  VSHUFPS(MJ, MJ, X2, 0x00);
  VSHUFPS(MJ, MJ, Y2, 0x55);
  VSHUFPS(MJ, MJ, Z2, 0xaa);

  VSUBPS_M(*ipdata_reg.x, X2, DX);
  VMULPS(DX, DX, X2);
  VADDPS(TWO, X2, X2);

  VSUBPS_M(*ipdata_reg.y, Y2, DY);
  VMULPS(DY, DY, Y2);
  VADDPS(X2, Y2, Y2);

  VSUBPS(ZI, Z2, DZ);
  VMULPS(DZ, DZ, Z2);
  VADDPS(Y2, Z2, Y2);

  VSHUFPS(MJ, MJ, MJ, 0xff);
  VMULPS(MJ, DX, DX);
  VMULPS(MJ, DY, DY);
  VMULPS(MJ, DZ, DZ);
  
  VMINPS(R2CUT, Y2, Z2);
#if NUNROLL == 2
  for(j = 0; j < nj; j += 2){
    VLOADPS(*jp, MJ);
    jp += 2;

    VEXTRACTF128(Z2, Y2_X, 0x01);

    VPSRLD(23-FRC_BIT, Y2_X, Y2_X);
    VPSRLD(23-FRC_BIT, Z2_X, X2_X);

    MEMORYFLUSH;
    VSTORPS(X2_X, idx[0]);
    VSTORPS(Y2_X, idx[4]);

    VPSLLD(23-FRC_BIT, Y2_X, Y2_X);
    VPSLLD(23-FRC_BIT, X2_X, X2_X);

    VPERM2F128(Y2, X2, Y2, 0x02);
    VSUBPS(Y2, Z2, Z2);

    VSHUFPS(MJ, MJ, X2, 0x00);
    VSHUFPS(MJ, MJ, Y2, 0x55);
    VSUBPS_M(*ipdata_reg.x, X2, X2);
    VSUBPS_M(*ipdata_reg.y, Y2, Y2);

    VLOADLPS(*fcut[idx[4]], BUF0_X, BUF0_X);
    VLOADHPS(*fcut[idx[5]], BUF0_X, BUF0_X);

    VLOADLPS(*fcut[idx[0]], BUF1_X, BUF1_X);
    VLOADHPS(*fcut[idx[1]], BUF1_X, BUF1_X);
    VPERM2F128(BUF0, BUF1, BUF1, 0x02); // 55441100
    VLOADLPS(*fcut[idx[6]], BUF2_X, BUF2_X);
    VLOADHPS(*fcut[idx[7]], BUF2_X, BUF2_X);
    VLOADLPS(*fcut[idx[2]], BUF0_X, BUF0_X);
    VLOADHPS(*fcut[idx[3]], BUF0_X, BUF0_X);
    VPERM2F128(BUF2, BUF0, BUF2, 0x02); // 77663322
    VSHUFPS(BUF1, BUF2, BUF0, 0xdd); // derivative
    VSHUFPS(BUF1, BUF2, BUF2, 0x88); // main

    VMULPS(Z2, BUF0, BUF0);
    VSHUFPS(MJ, MJ, Z2, 0xaa);
    VSHUFPS(MJ, MJ, MJ, 0xff);
    VSUBPS(ZI, Z2, Z2);
    VADDPS(BUF0, BUF2, BUF2);

    VMULPS(BUF2, DX, DX);
    VMULPS(BUF2, DY, DY);
    VMULPS(BUF2, DZ, DZ);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);

    VMOVAPS(X2, DX);
    VMOVAPS(Y2, DY);
    VMOVAPS(Z2, DZ);

    VMULPS(X2, X2, X2);
    VMULPS(Y2, Y2, Y2);
    VMULPS(Z2, Z2, Z2);

    VADDPS(TWO, X2, X2);
    VADDPS(Z2, Y2, Y2);
    VADDPS(X2, Y2, Y2);

    VMULPS(MJ, DX, DX);
    VMULPS(MJ, DY, DY);
    VMULPS(MJ, DZ, DZ);
    VMINPS(R2CUT, Y2, Z2);

  }
#elif NUNROLL == 4
  for(j = 0; j < nj; j += 4){
    //    PREFETCH(*(jp+4));

    VLOADPS(*jp, MJ);
    jp += 2;

    VEXTRACTF128(Z2, Y2_X, 0x01);

    VPSRLD(23-FRC_BIT, Y2_X, Y2_X);
    VPSRLD(23-FRC_BIT, Z2_X, X2_X);

    VEXTINT0(Y2_X, idx[0]);
    VEXTINT1(Y2_X, idx[1]);
    VEXTINT0(X2_X, idx[4]);
    VEXTINT1(X2_X, idx[5]);
    VEXTINT2(Y2_X, idx[2]);
    VEXTINT3(Y2_X, idx[3]);
    VEXTINT2(X2_X, idx[6]);
    VEXTINT3(X2_X, idx[7]);

    VPSLLD(23-FRC_BIT, Y2_X, Y2_X);
    VPSLLD(23-FRC_BIT, X2_X, X2_X);

    VPERM2F128(X2, Y2, Y2, 0x20);
    VSUBPS(Y2, Z2, Z2);

    VSHUFPS(MJ, MJ, X2, 0x00);
    VSHUFPS(MJ, MJ, Y2, 0x55);
    VSUBPS_M(*ipdata_reg.x, X2, X2);
    VSUBPS_M(*ipdata_reg.y, Y2, Y2);

    VLOADLPS(*fcut[idx[0]], BUF0_X, BUF0_X);
    VLOADHPS(*fcut[idx[1]], BUF0_X, BUF0_X);
    VLOADLPS(*fcut[idx[4]], BUF1_X, BUF1_X);
    VLOADHPS(*fcut[idx[5]], BUF1_X, BUF1_X);
    VPERM2F128(BUF0, BUF1, BUF1, 0x02); // 55441100
    VLOADLPS(*fcut[idx[2]], BUF2_X, BUF2_X);
    VLOADHPS(*fcut[idx[3]], BUF2_X, BUF2_X);
    VLOADLPS(*fcut[idx[6]], BUF0_X, BUF0_X);
    VLOADHPS(*fcut[idx[7]], BUF0_X, BUF0_X);
    VPERM2F128(BUF2, BUF0, BUF2, 0x02); // 77663322
    VSHUFPS(BUF1, BUF2, BUF0, 0xdd); // derivative
    VSHUFPS(BUF1, BUF2, BUF2, 0x88); // main

    VMULPS(Z2, BUF0, BUF0);
    VSHUFPS(MJ, MJ, Z2, 0xaa);
    VSHUFPS(MJ, MJ, MJ, 0xff);
    VSUBPS(ZI, Z2, Z2);
    VADDPS(BUF0, BUF2, BUF2);
    
    VMULPS(BUF2, DX, DX);
    VMULPS(BUF2, DY, DY);
    VMULPS(BUF2, DZ, DZ);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);

    VPERM2F128(X2, X2, DX, 0x10);
    VPERM2F128(Y2, Y2, DY, 0x10);
    VPERM2F128(Z2, Z2, DZ, 0x10);

    VMULPS(X2, X2, X2);
    VMULPS(Y2, Y2, Y2);
    VMULPS(Z2, Z2, Z2);

    VADDPS(TWO, X2, X2);
    VADDPS(Z2, Y2, Y2);
    VADDPS(X2, Y2, Y2);

    VMULPS(MJ, DX, DX);
    VMULPS(MJ, DY, DY);
    VMULPS(MJ, DZ, DZ);
    VMINPS(R2CUT, Y2, Z2);

    // ------------------------

    VLOADPS(*jp, MJ);
    jp += 2;

    VEXTRACTF128(Z2, Y2_X, 0x01);

    VPSRLD(23-FRC_BIT, Y2_X, Y2_X);
    VPSRLD(23-FRC_BIT, Z2_X, X2_X);

    VEXTINT0(Y2_X, idx[0]);
    VEXTINT1(Y2_X, idx[1]);
    VEXTINT0(X2_X, idx[4]);
    VEXTINT1(X2_X, idx[5]);
    VEXTINT2(Y2_X, idx[2]);
    VEXTINT3(Y2_X, idx[3]);
    VEXTINT2(X2_X, idx[6]);
    VEXTINT3(X2_X, idx[7]);

    VPSLLD(23-FRC_BIT, Y2_X, Y2_X);
    VPSLLD(23-FRC_BIT, X2_X, X2_X);

    VPERM2F128(X2, Y2, Y2, 0x20);
    VSUBPS(Y2, Z2, Z2);

    VSHUFPS(MJ, MJ, X2, 0x00);
    VSHUFPS(MJ, MJ, Y2, 0x55);
    VSUBPS_M(*ipdata_reg.x, X2, X2);
    VSUBPS_M(*ipdata_reg.y, Y2, Y2);

    VLOADLPS(*fcut[idx[0]], BUF0_X, BUF0_X);
    VLOADHPS(*fcut[idx[1]], BUF0_X, BUF0_X);
    VLOADLPS(*fcut[idx[4]], BUF1_X, BUF1_X);
    VLOADHPS(*fcut[idx[5]], BUF1_X, BUF1_X);
    VPERM2F128(BUF0, BUF1, BUF1, 0x02); // 55441100
    VLOADLPS(*fcut[idx[2]], BUF2_X, BUF2_X);
    VLOADHPS(*fcut[idx[3]], BUF2_X, BUF2_X);
    VLOADLPS(*fcut[idx[6]], BUF0_X, BUF0_X);
    VLOADHPS(*fcut[idx[7]], BUF0_X, BUF0_X);
    VPERM2F128(BUF2, BUF0, BUF2, 0x02); // 77663322
    VSHUFPS(BUF1, BUF2, BUF0, 0xdd); // derivative
    VSHUFPS(BUF1, BUF2, BUF2, 0x88); // main

    VMULPS(Z2, BUF0, BUF0);
    VSHUFPS(MJ, MJ, Z2, 0xaa);
    VSHUFPS(MJ, MJ, MJ, 0xff);
    VSUBPS(ZI, Z2, Z2);
    VADDPS(BUF0, BUF2, BUF2);

    VMULPS(BUF2, DX, DX);
    VMULPS(BUF2, DY, DY);
    VMULPS(BUF2, DZ, DZ);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);

    VPERM2F128(X2, X2, DX, 0x10);
    VPERM2F128(Y2, Y2, DY, 0x10);
    VPERM2F128(Z2, Z2, DZ, 0x10);

    VMULPS(X2, X2, X2);
    VMULPS(Y2, Y2, Y2);
    VMULPS(Z2, Z2, Z2);

    VADDPS(TWO, X2, X2);
    VADDPS(Z2, Y2, Y2);
    VADDPS(X2, Y2, Y2);

    VMULPS(MJ, DX, DX);
    VMULPS(MJ, DY, DY);
    VMULPS(MJ, DZ, DZ);
    VMINPS(R2CUT, Y2, Z2);
  }
#else
#error
#endif

  VPERM2F128(AX, AX, X2, 0x01);
  VADDPS(AX, X2, AX);
  VPERM2F128(AY, AY, Y2, 0x01);
  VADDPS(AY, Y2, AY);
  VPERM2F128(AZ, AZ, Z2, 0x01);
  VADDPS(AZ, Z2, AZ);

  VMULPS_M(accscale, AX_X, AX_X);
  VMULPS_M(accscale, AY_X, AY_X);
  VMULPS_M(accscale, AZ_X, AZ_X);
  
  VSTORPS(AX_X, *fodata->ax);
  VSTORPS(AY_X, *fodata->ay);
  VSTORPS(AZ_X, *fodata->az);

}
