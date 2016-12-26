#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cosmology.h"

#define FUNC(x,y) ((*func)(x,y))
COSM_REAL cosm_midpnt(COSM_REAL (*func)(COSM_REAL, struct cosmology), 
		  struct cosmology cosm,
		  COSM_REAL a, COSM_REAL b, int n)
{
  COSM_REAL x, tnm, sum, del, ddel;
  static COSM_REAL s;
  int it, j;
  
  if (n == 1) {
    return (s=(b-a)*FUNC(0.5*(a+b),cosm));
  } else{
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC(x,cosm);
      x += ddel;
      sum += FUNC(x,cosm);
      x += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}

COSM_REAL cosm_trapzd(COSM_REAL (*func)(COSM_REAL,struct cosmology), 
                   struct cosmology cosm, COSM_REAL a, COSM_REAL b, int n)
{
  COSM_REAL x,tnm,sum,del;
  static COSM_REAL s;
  int it,j;

  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC(a,cosm)+FUNC(b,cosm)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x,cosm);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}
#undef FUNC

void cosm_polint(COSM_REAL *xa,COSM_REAL *ya,int n,COSM_REAL x,COSM_REAL *y,COSM_REAL *dy)
{
  int i,m,ns=1;
  COSM_REAL den,dif,dift,ho,hp,w;
  COSM_REAL *c,*d;
   
  COSM_REAL *ctmp = (COSM_REAL *)malloc(n*sizeof(COSM_REAL));
  COSM_REAL *dtmp = (COSM_REAL *)malloc(n*sizeof(COSM_REAL));
  c = ctmp-1;
  d = dtmp-1;
   
  dif=fabs(x-xa[1]);

  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[(ns--)];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      den=ho-hp;
      if (den == 0.0){
        fprintf(stderr,"Error in routine cosm_polint\n");
        exit(1);
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
   
  free((void *)ctmp);
  free((void *)dtmp);
}



#define NITER (7)

COSM_REAL ddplus(COSM_REAL anow, struct cosmology cosm)
{
  if(anow == 0.0) {
    return 0.0e0;
  }else{
    COSM_REAL eta, omega_k;
    omega_k = 1.0-cosm.omega_m-cosm.omega_v;

    eta = sqrt(cosm.omega_m/anow + cosm.omega_v*anow*anow + omega_k);
    return 2.5/(eta*eta*eta);
  }
}

COSM_REAL dplus(COSM_REAL anow, struct cosmology cosm)
{
  int it;
  COSM_REAL result;
  COSM_REAL eta, omega_k;

  for(it=0;it<NITER;it++) {
    result = cosm_trapzd(ddplus, cosm, 0.0, anow, it+1);
  }
  
  omega_k = 1.0-cosm.omega_m-cosm.omega_v;
    
  eta = sqrt(cosm.omega_m/anow + cosm.omega_v*anow*anow + omega_k);

  return result*eta/anow;
  
}

COSM_REAL fomega(COSM_REAL anow, struct cosmology cosm)
/* dlog(D+)/dlog(a) */
{
  if(cosm.omega_m==1.0 && cosm.omega_v==0.0) {
    return 1.0;
  }else{
    COSM_REAL omega_k,eta;
    COSM_REAL fo;

    omega_k = 1.0-cosm.omega_m-cosm.omega_v;
    eta = sqrt(cosm.omega_m/anow+cosm.omega_v*anow*anow+omega_k);
    fo = (2.5/dplus(anow, cosm)-1.5*cosm.omega_m/anow-omega_k)/(eta*eta);

    return fo;
  }
}

COSM_REAL dladt(COSM_REAL anow, struct cosmology cosm)
/* dlog(a)/dtau where tau is conformal time */
{
  COSM_REAL eta, omega_k;
  
  omega_k = 1.0-cosm.omega_m-cosm.omega_v;

  eta = sqrt(cosm.omega_m/anow + cosm.omega_v*anow*anow + omega_k);

  return (anow*eta);
}

COSM_REAL dtda(COSM_REAL anow, struct cosmology cosm)
{
  COSM_REAL eta;
  COSM_REAL om,ov;

  om = cosm.omega_m;
  ov = cosm.omega_v;

  eta = sqrt(anow/(om+(1.0-om-ov)*anow+ov*anow*anow*anow));
  
  return  (eta);
}

COSM_REAL atotime(COSM_REAL anow, struct cosmology cosm)
{
  COSM_REAL ss, dss;
  static COSM_REAL h[21],s[21];
  int i,k=5;
  COSM_REAL eps=1.0e-7;
  
  h[1]=1.e0;
  for(i=1;i<=20;i++) {
    s[i] = cosm_trapzd(dtda,cosm,0.0,anow,i);
    if(i>=k){
      cosm_polint(&h[i-k],&s[i-k],k,0.e0,&ss,&dss);
      if(fabs(dss)<=eps*fabs(ss)) return ss;
    }
    s[i+1]=s[i];
    h[i+1]=0.25*h[i];
  }
  fprintf(stderr,"too many steps in atotime..\n");

  exit(EXIT_FAILURE);
}

COSM_REAL ztotime(COSM_REAL znow, struct cosmology cosm)
{
  COSM_REAL ss, dss, anow;
  static COSM_REAL h[21],s[21];
  int i,k=5;
  COSM_REAL eps=1.0e-7;
  
  anow = 1.0/(1.0+znow);

  h[1]=1.e0;
  for(i=1;i<=20;i++) {
    s[i] = cosm_trapzd(dtda,cosm,0.0,anow,i);
    if(i>=k){
      cosm_polint(&h[i-k],&s[i-k],k,0.e0,&ss,&dss);
      if(fabs(dss)<=eps*fabs(ss)) return ss;
    }
    s[i+1]=s[i];
    h[i+1]=0.25*h[i];
  }
  fprintf(stderr,"too many steps in ztotime..\n");
  fprintf(stderr,"omega_m = %14.6e / omega_b = %14.6e / omega_v = %14.6e\n",
	  cosm.omega_m, cosm.omega_b, cosm.omega_v);

  exit(EXIT_FAILURE);
}

#define MAXIT 100
COSM_REAL timetoa(COSM_REAL tnow, struct cosmology cosm)
{
  int j;
  COSM_REAL anow;
  COSM_REAL xacc;

  COSM_REAL df,dx,dxold,f,fh,fl;
  COSM_REAL temp,xh,xl,rts;
  COSM_REAL x1,x2;

  xacc = 1.0e-7;

  x1 = 1.e-4;
  x2 = 2.0;

  fl = atotime(x1, cosm)-tnow;
  fh = atotime(x2, cosm)-tnow;

  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)){
    fprintf(stderr,"Root must be bracketed in timetoa at tnow = %14.6e\n", tnow);
    fprintf(stderr,"fl at a = %14.6e :: %14.6e\n", x1, fl);
    fprintf(stderr,"fh at a = %14.6e :: %14.6e\n", x2, fh);
    exit(EXIT_FAILURE);
  }
  
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  f  = atotime(rts, cosm)-tnow;
  anow = 1.0/(1.0+rts);
  df = -dtda(anow,cosm)/(1.0+rts)/(1.0+rts);
  for (j=1;j<=MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
	|| (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    f  = atotime(rts, cosm)-tnow;
    anow = 1.0/(1.0+rts);
    df = -dtda(anow,cosm)/(1.0+rts)/(1.0+rts);
    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
  }
  fprintf(stderr, "Maximum number of iterations exceeded in timetoa\n");
  exit(EXIT_FAILURE);

}

COSM_REAL timetoz(COSM_REAL tnow, struct cosmology cosm)
{
  COSM_REAL anow;

  anow = timetoa(tnow, cosm);

  return (1.0/anow-1.0);

}

#undef MAXIT

COSM_REAL drift(COSM_REAL time, struct cosmology cosm)
{
  COSM_REAL anow;
  anow = timetoa(time, cosm);

  return (1.0/anow);
}

COSM_REAL kick(COSM_REAL time, struct cosmology cosm)
{
  COSM_REAL anow;
  anow = timetoa(time, cosm);

  return (1.0/anow);
}

COSM_REAL drift_factor(COSM_REAL time1, COSM_REAL time2, struct cosmology cosm)
{


  COSM_REAL drift_fac;
  int it;

  for(it=0;it<NITER;it++) {
    drift_fac = cosm_midpnt(drift, cosm, time1, time2, it+1);
  }

  return (drift_fac);

}

COSM_REAL kick_factor(COSM_REAL time1, COSM_REAL time2, struct cosmology cosm)
{

  COSM_REAL ss, dss;
  static COSM_REAL h[21],s[21];
  int i,k=5;
  COSM_REAL eps=1.0e-6;

  h[1]=1.e0;
  for(i=1;i<=20;i++) {
    s[i] = cosm_trapzd(kick,cosm,time1,time2,i);
    if(i>=k){
      cosm_polint(&h[i-k],&s[i-k],k,0.e0,&ss,&dss);
      if(fabs(dss)<=eps*fabs(ss)) return ss;
    }
    s[i+1]=s[i];
    h[i+1]=0.25*h[i];
  }
  fprintf(stderr,"too many steps in kick_factor..\n");
  exit(EXIT_FAILURE);
  
}

#if 0
int main(void)
{
  struct cosmology cosm;
  COSM_REAL anow;

  cosm.omega_m = 0.3;
  cosm.omega_v = 0.7;
  cosm.hubble = 0.7;

  //for(anow=1.e-3;anow<1.0;anow*=1.01){
    //    printf("%14.6e %14.6e\n", anow, timetoa(atotime(anow,cosm),cosm));
  //    printf("%14.6e %14.6e\n", anow, fomega(anow, cosm));
  //  }

  for(cosm.omega_m = 0.05; cosm.omega_m < 1.0; cosm.omega_m += 0.05) {
    cosm.omega_v = 1.0-cosm.omega_m;
    printf("%14.6e %14.6e\n", cosm.omega_m, fomega(1.0, cosm));
  }
    

  //  COSM_REAL time = 0.001;
  //  COSM_REAL dtime = 0.0001;
  //  printf("%14.6e\n",timetoa(time,cosm));
  //  printf("%14.6e %14.6e\n",dtime/timetoa(time+0.5*dtime,cosm),drift_factor(time, time+dtime, cosm));

  //  printf("%14.6e\n",ztotime(0.5062900,cosm));
  //  printf("%14.6e\n",timetoz(0.6340764,cosm));

}
#endif
