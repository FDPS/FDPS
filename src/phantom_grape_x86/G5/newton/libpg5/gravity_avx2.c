#include <assert.h>
#include "avx_type.h"

#include "sse.h"
#include "avx.h"
#include "avx2.h"

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

    VMOVAPS(EPS2, Y2);
    VFMADDPS(Y2, DX, DX);
    VFMADDPS(Y2, DY, DY);
    VFMADDPS(Y2, DZ, DZ);

    VRSQRTPS(Y2, X2);

    VFNMADDPS(PHI, X2, MJ);
    VMULPS(X2, MJ, MJ);
    VMULPS(X2, X2, Y2);
    VMULPS(MJ, Y2, Y2);

    VFMADDPS(AX, DX, Y2);
    VFMADDPS(AY, DY, Y2);
    VFMADDPS(AZ, DZ, Y2);

    VSHUFPS(J1, J1, X2, 0x00);
    VSHUFPS(J1, J1, J2, 0xaa);
    VSHUFPS(J1, J1, MJ, 0xff);
    VSHUFPS(J1, J1, Y2, 0x55);

  }
#elif (4 == NUNROLL)
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

    VMOVAPS(EPS2, Y2);
    VFMADDPS(Y2, DX, DX);
    VFMADDPS(Y2, DY, DY);
    VFMADDPS(Y2, DZ, DZ);    

    VRSQRTPS(Y2, X2);

    VFNMADDPS(PHI, X2, MJ);
    VMULPS(X2, MJ, MJ);
    VMULPS(X2, X2, Y2);
    VMULPS(MJ, Y2, Y2);

    VFMADDPS(AX, DX, Y2);
    VFMADDPS(AY, DY, Y2);
    VFMADDPS(AZ, DZ, Y2);

    VLOADPS(*(jpdata), J1);

    VSHUFPS(J2, J2, X2, 0x00);
    VSHUFPS(J2, J2, MJ, 0xff);
    VSHUFPS(J2, J2, Y2, 0x55);
    VSHUFPS(J2, J2, J2, 0xaa);

    VSUBPS(XI, X2, DX);
    VSUBPS(YI, Y2, DY);
    VSUBPS(ZI, J2, DZ);

    VMOVAPS(EPS2, Y2);
    VFMADDPS(Y2, DX, DX);
    VFMADDPS(Y2, DY, DY);
    VFMADDPS(Y2, DZ, DZ);    

    VLOADPS(*(jpdata+2), J2);    

    VRSQRTPS(Y2, X2);

    VFNMADDPS(PHI, X2, MJ);
    VMULPS(X2, MJ, MJ);
    VMULPS(X2, X2, Y2);
    VMULPS(MJ, Y2, Y2);

    VFMADDPS(AX, Y2, DX);
    VFMADDPS(AY, Y2, DY);
    VFMADDPS(AZ, Y2, DZ);

    VSHUFPS(J1, J1, X2, 0x00);
    VSHUFPS(J1, J1, MJ, 0xff);
    VSHUFPS(J1, J1, Y2, 0x55);
    VSHUFPS(J1, J1, J1, 0xaa);

    jpdata += 4;
    PREFETCH(*(jpdata));

  }
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

#define Z2 J1
#define EPSI2 EPS2
#define EPSJ2 J2
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

    VMOVAPS(EPSJ2, Y2);
    VFMADDPS(Y2, DX, DX);
    VFMADDPS(Y2, DY, DY);
    VFMADDPS(Y2, DZ, DZ);

    VLOADPS(jpdata->xm[0][0], Z2);
    VADDPS_M(jpdata->ep[0][0], EPSI2, EPSJ2);
    jpdata++;

    VRSQRTPS(Y2, X2);

    VFNMADDPS(PHI, MJ, X2);
    VMULPS(X2, MJ, MJ);
    VMULPS(X2, X2, Y2);
    VMULPS(MJ, Y2, Y2);

    VFMADDPS(AX, Y2, DX);
    VFMADDPS(AY, Y2, DY);
    VFMADDPS(AZ, Y2, DZ);

    VSHUFPS(Z2, Z2, X2, 0x00);
    VSHUFPS(Z2, Z2, MJ, 0xff);
    VSHUFPS(Z2, Z2, Y2, 0x55);
    VSHUFPS(Z2, Z2, Z2, 0xaa);

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
