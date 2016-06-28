#ifndef __COSMOLOGY_HPP__
#define __COSMOLOGY_HPP__

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace COSM {
  typedef float REAL;
}

namespace COSM {
  class cosmology {
  public:
    REAL omegam, omegav, omegab, omeganu;
    REAL hubble;
    
    REAL dtda(REAL anow)
    {
      REAL eta;
      REAL om,ov;
      
      om = omegam;
      ov = omegav;
      
      eta = sqrt(anow/(om+(1.0-om-ov)*anow+ov*anow*anow*anow));
      
      return  (eta);
    }
    
#define FUNC(x) ((this->*func)(x))
    REAL trapzd(REAL (cosmology::*func)(REAL),
		REAL &a, REAL &b, int &n)
    {
      REAL x,tnm,sum,del;
      static REAL s;
      int it,j;
      
      if (n == 1) {
	return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
      } else {
	for (it=1,j=1;j<n-1;j++) it <<= 1;
	tnm=it;
	del=(b-a)/tnm;
	x=a+0.5*del;
	for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
	s=0.5*(s+(b-a)*sum/tnm);
	return s;
      }
    }
#undef FUNC
    
    void polint(REAL *xa,REAL *ya,int n,REAL x,REAL *y,REAL *dy)
    {
      int i,m,ns=1;
      REAL den,dif,dift,ho,hp,w;
      REAL *c,*d;
      
      REAL *ctmp = (REAL *)malloc(n*sizeof(REAL));
      REAL *dtmp = (REAL *)malloc(n*sizeof(REAL));
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
	    fprintf(stderr,"Error in routine polint\n");
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
    
    REAL atotime(REAL anow) {
      REAL ss, dss, zero;
      static REAL h[21],s[21];
      int i,k=5;
      REAL eps=1.0e-7;
      zero = 0.0;
      
      h[1]=1.e0;
      for(i=1;i<=20;i++) {
	s[i] = trapzd(&cosmology::dtda,zero,anow,i);
	if(i>=k){
	  polint(&h[i-k],&s[i-k],k,zero,&ss,&dss);
	  if(fabs(dss)<=eps*fabs(ss)) return ss;
	}
	s[i+1]=s[i];
	h[i+1]=0.25*h[i];
      }
      fprintf(stderr,"too many steps in atotime..\n");
      exit(EXIT_FAILURE);
      
    }
    
#define MAXIT 100
    REAL timetoa(REAL tnow)
    {
      int j;
      REAL anow;
      REAL xacc;
      
      REAL df,dx,dxold,f,fh,fl;
      REAL temp,xh,xl,rts;
      REAL x1,x2;
      
      xacc = 1.0e-7;
      
      x1 = 1.e-4;
      x2 = 2.0;
      
      fl = atotime(x1)-tnow;
      fh = atotime(x2)-tnow;
      
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
      f  = atotime(rts)-tnow;
      anow = 1.0/(1.0+rts);
      df = -dtda(anow)/(1.0+rts)/(1.0+rts);
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
	f  = atotime(rts)-tnow;
	anow = 1.0/(1.0+rts);
	df = -dtda(anow)/(1.0+rts)/(1.0+rts);
	if (f < 0.0)
	  xl=rts;
	else
	  xh=rts;
      }
      fprintf(stderr, "Maximum number of iterations exceeded in timetoa\n");
      exit(EXIT_FAILURE);
      
    }
    
    REAL ztotime(REAL zred)
    {
      REAL anow = 1.0/(1.0+zred);
      
      return atotime(anow);
    }
    
    REAL timetoz(REAL tnow)
    {
      REAL anow = timetoa(tnow);
      
      return (1.0/anow - 1.0);
    }
    
  };
}

#endif /* __COSMOLOGY_HPP__ */
