#include <assert.h>
#include "avx_type.h"

#include "sse.h"
#include "avx.h"

#define AX    YMM08
#define AY    YMM09
#define AZ    YMM10
#define PHI   YMM11

#define DX    YMM12
#define DY    YMM13
#define DZ    YMM14
#define MJ    YMM07

#define J1    YMM00
#define J2    YMM01
#define X2    YMM02
#define Y2    YMM03

#define XI    YMM04
#define YI    YMM05
#define ZI    YMM06
#define EPS2  YMM15

void GravityKernel(pIpdata ipdata, pFodata fodata, pJpdata jpdata, int nj)
{
  int j;

  PREFETCH(jpdata[0]);

  VZEROALL;
  VLOADPS(*ipdata->x, XMM04);
  VLOADPS(*ipdata->y, XMM05);
  VLOADPS(*ipdata->z, XMM06);
  VLOADPS(*ipdata->eps2, XMM15);
  VPERM2F128(XI, XI, XI, 0x00);
  VPERM2F128(YI, YI, YI, 0x00);
  VPERM2F128(ZI, ZI, ZI, 0x00);
  VPERM2F128(EPS2, EPS2, EPS2, 0x00);

#if (2 == NUNROLL)
  VLOADPS(*(jpdata), J1);
  jpdata += 2;

  VSHUFPS(J1, J1, X2, 0x00);
  VSHUFPS(J1, J1, J2, 0xaa);
  VSHUFPS(J1, J1, MJ, 0xff);
  VSHUFPS(J1, J1, Y2, 0x55);

  for(j = 0; j < nj; j += 2){

    VSUBPS(XI, X2, DX);
    VSUBPS(ZI, J2, DZ);
    VSUBPS(YI, Y2, DY);

    VLOADPS(*(jpdata), J1);
    jpdata += 2;
    
    VMULPS(DX, DX, X2);
    VMULPS(DZ, DZ, J2);
    VMULPS(DY, DY, Y2);

    VADDPS(X2, J2, J2);
    VADDPS(EPS2, Y2, Y2);
    VADDPS(J2, Y2, Y2);

    VRSQRTPS(Y2, X2);

    VMULPS(X2, MJ, MJ);
    VMULPS(X2, X2, Y2);

    VMULPS(MJ, Y2, Y2);
    VSUBPS(MJ, PHI, PHI);

    VMULPS(Y2, DX, DX);
    VMULPS(Y2, DY, DY);
    VMULPS(Y2, DZ, DZ);

    VSHUFPS(J1, J1, X2, 0x00);
    VSHUFPS(J1, J1, J2, 0xaa);
    VSHUFPS(J1, J1, MJ, 0xff);
    VSHUFPS(J1, J1, Y2, 0x55);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);
  }
#elif (4 == NUNROLL)
#if 1
  VLOADPS(*(jpdata), J1);
  VLOADPS(*(jpdata+2), J2);

  jpdata += 4;

  VSHUFPS(J1, J1, X2, 0x00);
  VSHUFPS(J1, J1, Y2, 0x55);
  VSHUFPS(J1, J1, MJ, 0xff);
  VSHUFPS(J1, J1, J1, 0xaa);

  for(j = 0 ; j < nj; j += 4) {

    VSUBPS(XI, X2, DX);
    VSUBPS(YI, Y2, DY);
    VSUBPS(ZI, J1, DZ);

    VMULPS(DX, DX, X2);
    VMULPS(DZ, DZ, J1);
    VMULPS(DY, DY, Y2);

    VADDPS(J1, X2, X2);
    VADDPS(EPS2, Y2, Y2);
    VADDPS(Y2, X2, Y2);

    VLOADPS(*(jpdata), J1);

    VRSQRTPS(Y2, X2);

    VMULPS(X2, MJ, MJ);
    VMULPS(X2, X2, Y2);

    VMULPS(MJ, Y2, Y2);
    VSUBPS(MJ, PHI, PHI);

    VMULPS(Y2, DX, DX);
    VMULPS(Y2, DY, DY);
    VMULPS(Y2, DZ, DZ);

    VSHUFPS(J2, J2, X2, 0x00);
    VSHUFPS(J2, J2, MJ, 0xff);
    VSHUFPS(J2, J2, Y2, 0x55);
    VSHUFPS(J2, J2, J2, 0xaa);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);

    VSUBPS(XI, X2, DX);
    VSUBPS(YI, Y2, DY);
    VSUBPS(ZI, J2, DZ);

    VMULPS(DX, DX, X2);
    VMULPS(DZ, DZ, J2);
    VMULPS(DY, DY, Y2);

    VADDPS(J2, X2, X2);
    VADDPS(EPS2, Y2, Y2);
    VADDPS(Y2, X2, Y2);

    VLOADPS(*(jpdata+2), J2);    

    VRSQRTPS(Y2, X2);

    VMULPS(X2, MJ, MJ);
    VMULPS(X2, X2, Y2);

    jpdata += 4;
    PREFETCH(*(jpdata));

    VMULPS(MJ, Y2, Y2);
    VSUBPS(MJ, PHI, PHI);

    VMULPS(Y2, DX, DX);
    VMULPS(Y2, DY, DY);
    VMULPS(Y2, DZ, DZ);

    VSHUFPS(J1, J1, X2, 0x00);
    VSHUFPS(J1, J1, MJ, 0xff);
    VSHUFPS(J1, J1, Y2, 0x55);
    VSHUFPS(J1, J1, J1, 0xaa);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);
  }
#else
  VLOADPS(*(jpdata), J1);
  VLOADPS(*(jpdata+2), J2);

  jpdata += 4;

  VSHUFPS(J1, J1, X2, 0x00);
  VSHUFPS(J1, J1, Y2, 0x55);
  VSHUFPS(J1, J1, MJ, 0xaa);
  VSHUFPS(J1, J1, J1, 0xff);

  for(j = 0 ; j < nj; j += 4) {

    VSUBPS(XI, X2, DX);
    VSUBPS(YI, Y2, DY);
    VSUBPS(ZI, MJ, DZ);

    VMULPS(DX, DX, X2);
    VMULPS(DY, DY, Y2);
    VMULPS(DZ, DZ, MJ);

    VADDPS(X2, Y2, Y2);
    VADDPS(EPS2, MJ, MJ);
    VADDPS(Y2, MJ, Y2);

    VRSQRTPS(Y2, X2);

    VMULPS(X2, J1, Y2);
    VMULPS(X2, X2, X2);

    VLOADPS(*(jpdata), J1);

    VSUBPS(Y2, PHI, PHI);
    VMULPS(X2, Y2, Y2);

    VMULPS(Y2, DX, DX);
    VMULPS(Y2, DY, DY);
    VMULPS(Y2, DZ, DZ);

    VSHUFPS(J2, J2, X2, 0x00);
    VSHUFPS(J2, J2, Y2, 0x55);
    VSHUFPS(J2, J2, MJ, 0xaa);
    VSHUFPS(J2, J2, J2, 0xff);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);

    VSUBPS(XI, X2, DX);
    VSUBPS(YI, Y2, DY);
    VSUBPS(ZI, MJ, DZ);

    VMULPS(DX, DX, X2);
    VMULPS(DY, DY, Y2);
    VMULPS(DZ, DZ, MJ);

    VADDPS(X2, Y2, Y2);
    VADDPS(EPS2, MJ, MJ);
    VADDPS(Y2, MJ, Y2);

    VRSQRTPS(Y2, X2);

    VMULPS(X2, J2, Y2);
    VMULPS(X2, X2, X2);

    VLOADPS(*(jpdata+2), J2);    

    jpdata += 4;
    PREFETCH(*(jpdata));

    VSUBPS(Y2, PHI, PHI);
    VMULPS(X2, Y2, Y2);

    VMULPS(Y2, DX, DX);
    VMULPS(Y2, DY, DY);
    VMULPS(Y2, DZ, DZ);

    VSHUFPS(J1, J1, X2, 0x00);
    VSHUFPS(J1, J1, Y2, 0x55);
    VSHUFPS(J1, J1, MJ, 0xaa);
    VSHUFPS(J1, J1, J1, 0xff);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);
  }
#endif
#else
#error
#endif

  VEXTRACTF128(AX, XMM00, 0x01);
  VADDPS(AX, YMM00, AX);
  VEXTRACTF128(AY, XMM01, 0x01);
  VADDPS(AY, YMM01, AY);
  VEXTRACTF128(AZ, XMM02, 0x01);
  VADDPS(AZ, YMM02, AZ);
  VEXTRACTF128(PHI, XMM03, 0x01);
  VADDPS(PHI, YMM03, PHI);

  VSTORPS(XMM08,  *fodata->ax);
  VSTORPS(XMM09,  *fodata->ay);
  VSTORPS(XMM10,  *fodata->az);
  VSTORPS(XMM11, *fodata->phi);

}

#define DEBUG 1
#define Z2 J1
#define EPSI2 EPS2
#define EPSJ2 J2
#if DEBUG == 0
#else
void GravityKernel0(pIpdata ipdata, pFodata fodata, pJpdata0 jpdata, int nj)
{
  int j;

  PREFETCH(jpdata[0]);

  VZEROALL;
  VLOADPS(*ipdata->x, XMM04);
  VLOADPS(*ipdata->y, XMM05);
  VLOADPS(*ipdata->z, XMM06);
  VLOADPS(*ipdata->eps2, XMM15);
  VPERM2F128(XI, XI, XI, 0x00);
  VPERM2F128(YI, YI, YI, 0x00);
  VPERM2F128(ZI, ZI, ZI, 0x00);
  VPERM2F128(EPSI2, EPSI2, EPSI2, 0x00);

  VLOADPS(jpdata->xm[0][0], Z2);
  VADDPS_M(jpdata->ep[0][0], EPSI2, EPSJ2);
  jpdata++;

  VSHUFPS(Z2, Z2, X2, 0x00);
  VSHUFPS(Z2, Z2, MJ, 0xff);
  VSHUFPS(Z2, Z2, Y2, 0x55);
  VSHUFPS(Z2, Z2, Z2, 0xaa);

  for(j = 0; j < nj; j += 2){

    VSUBPS(XI, X2, DX);
    VSUBPS(ZI, Z2, DZ);
    VSUBPS(YI, Y2, DY);

    VMULPS(DX, DX, X2);
    VMULPS(DZ, DZ, Z2);
    VMULPS(DY, DY, Y2);

    VADDPS(X2, Z2, X2);
    VADDPS(EPSJ2, Y2, Y2);
    VADDPS(X2, Y2, Y2);

    VLOADPS(jpdata->xm[0][0], Z2);
    VADDPS_M(jpdata->ep[0][0], EPSI2, EPSJ2);
    jpdata++;

    VRSQRTPS(Y2, X2);

    VMULPS(X2, MJ, MJ);
    VMULPS(X2, X2, Y2);

    VMULPS(MJ, Y2, Y2);
    VSUBPS(MJ, PHI, PHI);

    VMULPS(Y2, DX, DX);
    VMULPS(Y2, DY, DY);
    VMULPS(Y2, DZ, DZ);

    VSHUFPS(Z2, Z2, X2, 0x00);
    VSHUFPS(Z2, Z2, MJ, 0xff);
    VSHUFPS(Z2, Z2, Y2, 0x55);
    VSHUFPS(Z2, Z2, Z2, 0xaa);

    VADDPS(DX, AX, AX);
    VADDPS(DY, AY, AY);
    VADDPS(DZ, AZ, AZ);

  }

  VEXTRACTF128(AX, XMM00, 0x01);
  VADDPS(AX, YMM00, AX);
  VEXTRACTF128(AY, XMM01, 0x01);
  VADDPS(AY, YMM01, AY);
  VEXTRACTF128(AZ, XMM02, 0x01);
  VADDPS(AZ, YMM02, AZ);
  VEXTRACTF128(PHI, XMM03, 0x01);
  VADDPS(PHI, YMM03, PHI);

  VSTORPS(XMM08,  *fodata->ax);
  VSTORPS(XMM09,  *fodata->ay);
  VSTORPS(XMM10,  *fodata->az);
  VSTORPS(XMM11, *fodata->phi);

}
#endif
