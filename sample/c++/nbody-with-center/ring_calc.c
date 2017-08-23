//ring_calc.c
// create division
//
//  Version 0.0 2017/4/11 J. Makino 
//
//usage:
//ring_calc
//enter:width nx ny
// Sampple:  ring_calc
// 0.01 10 10
/*
  High-level function:
  This function assumes the ring with is dt and ring central radius is 1
  x and y divisions are nx and ny.
    one can get box dimensions  of box ix iy (starting fron 0) by calling
  void box(int nx,
	 int ny,
	 int ix,
	 int iy,
	 double dr,
	 double *px0,
	 double *px1,
	 double *py0,
	 double *py1 )

This function returns the range of angle for which a circle of radius
r interesect with a box defined by x0, x1, y0, y1 (left and right x
coordinates and bottom and top y coordinates. 0<=x0<x1 and 0<=y0<y1)

void  determine_range_first_q(double x0,
			     double x1,
			     double y0,
			     double y1,
			     double r,
			     double * theta0,
			     double * theta1)
	 
Low-level functions	 
  ring has the inner radius r0 and outer radius r1
  we divide [0,1] to n "equal are", by calling 
    calc_x0(r0,r1,n,x1)
  n times. We first need to set x1=r1, where x1 is the
  right boundary of the rightmost box
  this function returns x0, left boundary of the box
  then you can get the next box by calling the same fucnction
  with x0 instead of x0.

  The function 
  calc_y0(x0,x1,r0,r1,m,y1);
  is very similar, except that it does the y-direction division
  and m is total process count (in 1st quadrant)

The function  
int count_particles(double r0,
		     double r1,
		     int nr,
		     int nt,
		     double x0,
		     double x1,
		     double y0,
		     double y1)

shows the example usage of determine_range_first_q. It actually
generates particle locations

*/

#include <stdio.h>
#include <math.h>

double cut_fan_area(double r,double x)
{
    if ( x >= r){
	return 0.0;
    }else{
	double theta=acos(x/r);
	return    (r*r*theta -r*x*sin(theta))/2;
    }
}

double band_area(double x0,double x1, double r0, double r1)
{
    return    cut_fan_area(r1,x0)-cut_fan_area(r1,x1) -
	cut_fan_area(r0,x0)+cut_fan_area(r0,x1) ;
}

double  patch_dx(double y,double x0,double x1, double r0, double r1)
{
    
    double yx0r1 = sqrt(r1*r1-x0*x0);
    double  yx0r0=0;
    double  yx1r1=sqrt(r1*r1-x1*x1);
    double  yx1r0=0;
    if (r0>x0){
	yx0r0 = sqrt(r0*r0-x0*x0);
    }
    if (r0>x1){
	yx1r0 = sqrt(r0*r0-x1*x1) ;
    }
    double  xmin=x0;
    double      xmax=x1;
    if (y < yx0r0){
	xmin = sqrt(r0*r0-y*y);
    }
    if (y> yx1r1){
	xmax = sqrt(r1*r1-y*y);
    }
    double dx = xmax-xmin;
    if (dx<0.0)  dx=0.0;
    return dx;
}

double  patch_area(double y0,
		   double y1,
		   double x0,
		   double x1,
		   double r0,
		   double r1)
{
    int  n=20;
    double  s= (patch_dx(y0,x0,x1,r0,r1)+patch_dx(y1,x0,x1,r0,r1))/2;
    int i;
    for(i=0;i<n-1;i++){
	s+=patch_dx(y0+(i+1)*(y1-y0)/n,x0,x1,r0,r1);
    }
    return   s/n*(y1-y0);
}

double  calc_x0(double r0, double r1, int n, double x1)
{
    double s= M_PI*(r1*r1-r0*r0)/n/4;
    double  a=0;
    double  b=x1;
    double  x=b;
    int i;
    for (i=0;i<50;i++){
	x= (a+b)/2;
	if (band_area(x,x1,r0,r1)-s > 0.0){
	    a=x;
	}else{
	    b=x;
	}
    }
    return x;
}

double  calc_y0(double x0, double x1, double r0, double r1, int n, double y1)
{
    //    fprintf(stderr,"calc y0 called with %e %e %e %e  %d %e\n",
    //	    x0, x1, r0, r1, n, y1);
    double s= M_PI*(r1*r1-r0*r0)/n/4;
    double  yx1r0=0;
    if (r0>x1){
	yx1r0 = sqrt(r0*r0-x1*x1) ;
    }
    double a=yx1r0;
    double  b=y1;
    double   y=b;
    int i;
    for (i=0;i<50;i++){
	y= (a+b)/2;
	if (patch_area(y,y1,x0,x1,r0,r1)-s > 0.0){
	    a=y;
	}else{
	    b=y;
	}
    }
    return y;
}

void box(int nx,
	 int ny,
	 int ix,
	 int iy,
	 double dr,
	 double *px0,
	 double *px1,
	 double *py0,
	 double *py1 )
{
    int i,j;
    double r0=1-dr/2;
    double r1=1+dr/2;
    double x1=r1;
    double x1old;
    fprintf(stderr, "test routine with ix, iy=%d %d\n", ix, iy);
    double x0;
    for (i=nx;i>ix;i--){
	 x0=calc_x0(r0,r1,nx,x1);
	 x1old=x1;
	 x1=x0;
    }
    
    double y1 = sqrt(r1*r1-x0*x0);
    double y1old;
    double y0;
    for(j=ny;j>iy;j--){
	y0=calc_y0(x0,x1old,r0,r1,nx*ny,y1);
	//	fprintf(stderr,"y0, y1= %e %e\n", y0,y1);
	y1old=y1;
	y1=y0;
    }
    *px0=x0;
    *px1=x1old;
    *py0=y0;
    *py1=y1old;
}

void determine_range_first_q(double x0,
			     double x1,
			     double y0,
			     double y1,
			     double r,
			     double * theta0,
			     double * theta1)
{
    double tmin0=0;
    double tmax0=M_PI/2;
    double tmin1=0;
    double tmax1=M_PI/2;
    if (x1 <= r) tmin0 = acos(x1/r);
    if (x0<= r)  tmax0 = acos(x0/r);
    if ((x1>r)&&(x0>r)){
	tmin0=tmax0=0;
    }
    if (y1 <= r) tmax1= asin(y1/r);
    if (y0 <= r) tmin1 = asin(y0/r);
    if ((y1>r)&&(y0>r)){
	tmin1=tmax1=0;
    }

    *theta0= tmin0;
    if (tmin1 > tmin0) *theta0=tmin1;
    *theta1= tmax0;
    if (tmax1 < tmax0) *theta1=tmax1;
    if (*theta1 <= *theta0)  *theta1 = *theta0;
    //    fprintf(stderr,"range for %e %e %e %e %e: %e %e\n",
    //	    x0, x1, y0, y1,r, *theta0, *theta1);
}

    
    
#ifdef TESTMAIN
int main()
{
    fprintf(stderr,"Enter dr, n, m:");
    double dr;
    int n, m;
    scanf("%lf%d%d", &dr, &n, &m);
    fprintf(stderr,"dr, n, m= %e, %d %d\n", dr, n, m);
    double r0=1-dr/2;
    double r1=1+dr/2;
    double x1=r1;
    printf( "rmax= %e\n",r1);
    int i;
    for (i=0;i<n;i++){
	double x0=calc_x0(r0,r1,n,x1);
	double y1 = sqrt(r1*r1-x0*x0);
	int j;
	for(j=0;j<m;j++){
	    double y0=calc_y0(x0,x1,r0,r1,n*m,y1);
	    //	    fprintf(stderr,"Y0, Y1= %e %e\n", y0,y1);
	    //    printf("box %e %e %e %e\n", x0, x1, y0, y1);
	    fprintf(stderr,"box %e %e %e %e\n", x0, x1, y0, y1);
	    double xx0,xx1, yy0, yy1;
	    box( n,m, n-1-i, m-1-j, dr,
		&xx0, &xx1, &yy0, &yy1);
	    fprintf(stderr,"box %e %e %e %e\n", xx0, xx1, yy0, yy1);
	    printf("box %e %e %e %e\n", xx0, xx1, yy0, yy1);

	    y1=y0;
	}
    x1=x0;
    }
}
#endif

int count_particles(double r0,
		     double r1,
		     int nr,
		     int nt,
		     double x0,
		     double x1,
		     double y0,
		     double y1)
{
    int i;
    double dtheta = 2*M_PI/nt;
    int sum=0;
    for (i=0;i<nr;i++){
	double r = r0+ (r1-r0)*(i+0.5)/nr;
	double  theta0;
	double  theta1;
	determine_range_first_q(x0, x1, y0, y1, r,&theta0,& theta1);
	int i0= theta0/dtheta;
	int i1= (theta1/dtheta);
	//	fprintf(stderr,"i0,i1 = %d %d\n", i0, i1);
	sum += i1-i0;
#ifdef PLOT_PARTICLES	
	int ii;
	for(ii=i0;ii<i1;ii++){
	    double theta=ii*dtheta;
	    printf("dot %.8f %.8f\n", r*cos(theta), r*sin(theta));
	}
#endif
    }
    return sum;
}


#ifdef TESTPARTICLEGENERATION
int main()
{
    fprintf(stderr,"Enter dr, nx, ny, nr, nt:");
    double dr;
    int nx, ny, nr, nt;
    scanf("%lf%d%d%d%d", &dr, &nx, &ny, &nr, &nt);
    fprintf(stderr,"dr, nx, ny, nr, nt= %e, %d %d %d %d\n", dr, nx, ny, nr, nt);
    double r0=1-dr/2;
    double r1=1+dr/2;
    double x1=r1;
    printf( "rmax= %e\n",r1);
    int i;
    for (i=0;i<nx;i++){
	int j;
	for(j=0;j<ny;j++){
	    double x0,x1, y0, y1;
	    box( nx,ny, i, j, dr,
		&x0, &x1, &y0, &y1);
	    int np=count_particles(r0,r1,nr,nt,x0, x1, y0, y1);
	    printf("box %e %e %e %e %d\n", x0, x1, y0, y1,np);
	}
    }
}
#endif

  
  
