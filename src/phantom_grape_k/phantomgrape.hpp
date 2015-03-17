#pragma once

#include "emmintrin.h"

template<class T>
class Vec4{
public:
    T x,y,z,w;
    Vec4() : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {
        assert(sizeof(Vec4) == sizeof(T)*4);
    }
    Vec4(const T &_z, const T &_y, const T &_x, const T &_w) : x(_x), y(_y), z(_z), w(_w), {
        assert(sizeof(Vec4) == sizeof(T)*4);
    }
    Vec4(const T &_x){
        assert(sizeof(Vec4) == sizeof(T)*4);
        z = _x;
        y = _x;
        x = _x;
        w = _x;
    }
};

typedef Vec4<float> v4r4;
typedef Vec4<double> v4r8;

typedef __builtin_v2r8 V2R8;

#define V2R8_set    __builtin_fj_set_v2r8
#define V2R8_storel __builtin_fj_storel_v2r8
#define V2R8_storeh __builtin_fj_storeh_v2r8
#define V2R8_mul    __builtin_fj_mul_v2r8
#define V2R8_add    __builtin_fj_add_v2r8
#define V2R8_sub    __builtin_fj_sub_v2r8
#define V2R8_max    __builtin_fj_max_v2r8
#define V2R8_min    __builtin_fj_min_v2r8
#define V2R8_rcp    __builtin_fj_rcpa_v2r8
#define V2R8_rsqrt  __builtin_fj_rsqrta_v2r8
#define V2R8_madd   __builtin_fj_madd_v2r8
#define V2R8_msub   __builtin_fj_msub_v2r8
#define V2R8_nmadd  __builtin_fj_nmadd_v2r8
#define V2R8_nmsub  __builtin_fj_nmsub_v2r8
#define V2R8_abs    __builtin_fj_abs_v2r8
#define V2R8_and    __builtin_fj_and_v2r8
#define V2R8_cmplt  __builtin_fj_cmplt_v2r8
#define V2I8_slli __builtin_fj_slli_v2i8
#define V2R8_setzero __builtin_fj_setzero_v2r8
#define immr8(x)    __builtin_fj_set_v2r8(x, x)


//#warning COMPILING WITHOUT GNU C COMPATIBLE MODE
#define PREFETCH(addr) (void)(addr)
#define INLINE
#define NOINLINE

#include <cassert>

class PhantomGRAPE{
public:
	enum{
        PG_NIMAX = 4096,
		PG_NJMAX = 131072,
	};
    
private:
	double xscale, ascale, eps2, pad[5];
	V2R8 xjbuf[PG_NJMAX/2][4];
	V2R8 xibuf[PG_NIMAX][3];
	V2R8 aibuf[PG_NIMAX][3];
	V2R8 pibuf[PG_NIMAX];

public:
	PhantomGRAPE() : xscale(1.0), ascale(1.0), eps2(1.e-4) {}
	void set_eps2(const double _eps2){
		this->eps2 = _eps2;
	}

	INLINE
	void set_xj(const int nj, const double xj[][4]){ // x, y, z, m
		assert(nj <= PG_NJMAX);
		const V2R8 ss = immr8(this->xscale);
		for(int jj=0; jj<nj/2; jj++){
			const int j0  = 2*jj + 0;
			const int j1  = 2*jj + 1;
			const V2R8 x = V2R8_mul(ss, V2R8_set(xj[j0][0], xj[j1][0]));
			const V2R8 y = V2R8_mul(ss, V2R8_set(xj[j0][1], xj[j1][1]));
			const V2R8 z = V2R8_mul(ss, V2R8_set(xj[j0][2], xj[j1][2]));
			const V2R8 m =         (    V2R8_set(xj[j0][3], xj[j1][3]));
			xjbuf[jj][0] = x;
			xjbuf[jj][1] = y;
			xjbuf[jj][2] = z;
			xjbuf[jj][3] = m;
		}
		if(nj%2 == 1){
			const int jj = nj/2;
			const int j0  = 2*jj + 0;
			const V2R8 x = V2R8_mul(ss, V2R8_set(xj[j0][0], 1.0));
			const V2R8 y = V2R8_mul(ss, V2R8_set(xj[j0][1], 1.0));
			const V2R8 z = V2R8_mul(ss, V2R8_set(xj[j0][2], 1.0));
			const V2R8 m =         (    V2R8_set(xj[j0][3], 0.0));
			xjbuf[jj][0] = x;
			xjbuf[jj][1] = y;
			xjbuf[jj][2] = z;
			xjbuf[jj][3] = m;
		}
	}

	void set_xj(const int nj, 
                const v4r4 xj[]){
		assert(nj <= PG_NJMAX);
		for(int jj=0; jj<nj/2; jj++){
			const int j0  = 2*jj + 0;
			const int j1  = 2*jj + 1;
			const V2R8 x = V2R8_set((double)xj[j0].x, (double)xj[j1].x);
			const V2R8 y = V2R8_set((double)xj[j0].y, (double)xj[j1].y);
			const V2R8 z = V2R8_set((double)xj[j0].z, (double)xj[j1].z);
			const V2R8 m = V2R8_set((double)xj[j0].w, (double)xj[j1].w);
			xjbuf[jj][0] = x;
			xjbuf[jj][1] = y;
			xjbuf[jj][2] = z;
			xjbuf[jj][3] = m;
		}
		if(nj%2 == 1){
			const int jj = nj/2;
			const int j0  = 2*jj + 0;
			const V2R8 x = V2R8_set((double)xj[j0].x, 1.0);
			const V2R8 y = V2R8_set((double)xj[j0].y, 1.0);
			const V2R8 z = V2R8_set((double)xj[j0].z, 1.0);
			const V2R8 m = V2R8_set((double)xj[j0].w, 0.0);
			xjbuf[jj][0] = x;
			xjbuf[jj][1] = y;
			xjbuf[jj][2] = z;
			xjbuf[jj][3] = m;
		}
	}

    INLINE
    void set_xj_one(const unsigned long addr,
                    const double x,
                    const double y,
                    const double z,
                    const double m){
        const unsigned long ah = addr/2;
        const unsigned long al = addr%2;
        const double s = this->xscale;
        double *p = (double *)(xjbuf[ah]);
        p[al+0] = x * s;
        p[al+2] = y * s;
        p[al+4] = z * s;
        p[al+6] = m;
	}
		
	INLINE
	void set_xi(const int ni, const double xi[][3]){
		assert(ni <= PG_NIMAX);
		const double s = this->xscale;
		for(int i=0; i<ni; i++){
			const V2R8 x = immr8(s * xi[i][0]);
			const V2R8 y = immr8(s * xi[i][1]);
			const V2R8 z = immr8(s * xi[i][2]);
			xibuf[i][0] = x;
			xibuf[i][1] = y;
			xibuf[i][2] = z;
		}
	}

	void set_xi(const int ni, const v4r4 xi[]){
		assert(ni <= PG_NIMAX);
		for(int i=0; i<ni; i++){
			const V2R8 x = immr8((double)xi[i].x);
			const V2R8 y = immr8((double)xi[i].y);
			const V2R8 z = immr8((double)xi[i].z);
			xibuf[i][0] = x;
			xibuf[i][1] = y;
			xibuf[i][2] = z;
		}
	}

	INLINE
    void get_ai(const int ni, double ai[][3]) const {
		assert(ni <= PG_NIMAX);
		const double s = this->ascale;
		for(int i=0; i<ni; i++){
			const double *ptr = (const double *)aibuf[i];
			const double ax = s * (ptr[0] + ptr[1]);
			const double ay = s * (ptr[2] + ptr[3]);
			const double az = s * (ptr[4] + ptr[5]);
			ai[i][0] = ax;
			ai[i][1] = ay;
			ai[i][2] = az;
		}
	}

	// added by T.I. 20121128
	void get_ai(const int ni, double ai[][3], double *pi) const {
		assert(ni <= PG_NIMAX);
		const double s = this->ascale;
		for(int i=0; i<ni; i++){
			const double *ptr = (const double *)aibuf[i];
			const double ax = s * (ptr[0] + ptr[1]);
			const double ay = s * (ptr[2] + ptr[3]);
			const double az = s * (ptr[4] + ptr[5]);
			const double *ptr2 = (const double *)(pibuf+i);
			pi[i] = ptr2[0] + ptr2[1];
			ai[i][0] = ax;
			ai[i][1] = ay;
			ai[i][2] = az;
		}
	}

	void get_ai(const int ni, v4r4 ap[]) const {
		assert(ni <= PG_NIMAX);
		for(int i=0; i<ni; i++){
			const double *ptr = (const double *)aibuf[i];
			const double ax = ptr[0] + ptr[1];
			const double ay = ptr[2] + ptr[3];
			const double az = ptr[4] + ptr[5];
			const double *ptr2 = (const double *)(pibuf+i);
			ap[i].w = ptr2[0] + ptr2[1];
			ap[i].z = az;
			ap[i].y = ay;
			ap[i].x = ax;
		}
	}

	INLINE
	void run(const int ni, const int nj){
		assert(ni <= PG_NIMAX);
		assert(nj <= PG_NJMAX);
		this->gravity(ni, nj);
	}


private:

// 9op (3madd/msub + 3)
	INLINE
	static V2R8 rsqrt_x3(const V2R8 a){
		const V2R8 c0  = V2R8_set(1.0, 1.0);
		const V2R8 c1  = V2R8_set(0.5, 0.5);
		const V2R8 c2  = V2R8_set(0.375, 0.375);
		const V2R8 x0  = V2R8_rsqrt(a);
		const V2R8 ax0 = V2R8_mul(a, x0);
		const V2R8 h0  = V2R8_nmsub(ax0, x0, c0);
		const V2R8 hpoly = V2R8_madd(V2R8_madd(c2, h0, c1), h0, c0);
		return V2R8_mul(x0, hpoly);
	}

	static V2R8 rsqrt_x3_2(const V2R8 a){
		const V2R8 c0  = V2R8_set(1.0, 1.0);
		const V2R8 c1  = V2R8_set(0.5, 0.5);
		const V2R8 c2  = V2R8_set(0.375, 0.375);
		const V2R8 x0  = V2R8_rsqrt(a);
		const V2R8 ax0 = V2R8_mul(a, x0);
		const V2R8 h0  = V2R8_nmsub(ax0, x0, c0);
		const V2R8 hpoly0 = V2R8_madd(V2R8_madd(c2, h0, c1), h0, c0);
        const V2R8 x1 = V2R8_mul(x0, hpoly0); 
		const V2R8 ax1 = V2R8_mul(a, x1);
		const V2R8 h1  = V2R8_nmsub(ax1, x1, c0);
		const V2R8 hpoly1 = V2R8_madd(V2R8_madd(c2, h1, c1), h1, c0);
		return V2R8_mul(x1, hpoly1); 
	}

	// added by T.I 20121128
// 19op (6madd/msub + 7) + 9op(3madd/msub + 3)
	static void gravity_kernel(const V2R8 xi, const V2R8 yi, const V2R8 zi, const V2R8 eps2,
                               const V2R8 xj, const V2R8 yj, const V2R8 zj, const V2R8 mj,
                               V2R8 &ax, V2R8 &ay, V2R8 &az, V2R8 &p){
        const V2R8 dx = V2R8_sub(xj, xi);
        const V2R8 dy = V2R8_sub(yj, yi);
        const V2R8 dz = V2R8_sub(zj, zi);
        const V2R8 r2 = V2R8_madd(dz, dz, V2R8_madd(dy, dy, V2R8_madd(dx, dx, eps2)));
        const V2R8 ri = rsqrt_x3(r2);
        //const V2R8 ri = rsqrt_x3_2(r2);
		//const V2R8 ri  = V2R8_rsqrt(r2);
        p   = V2R8_nmsub( mj, ri, p); // count as 1op
        const V2R8 mri3 = V2R8_mul( V2R8_mul(mj, ri), V2R8_mul(ri, ri));
        ax  = V2R8_madd(mri3, dx, ax);
        ay  = V2R8_madd(mri3, dy, ay);
        az  = V2R8_madd(mri3, dz, az);
	}

	// NOINLINE
	void gravity(const int ni, const int nj){
		const int ni4 = ni/4 + (ni%4 ? 1 : 0);
		const int nj2 = nj/2 + (nj%2 ? 1 : 0);
        const V2R8 veps2 = immr8(eps2);
		for(int i4=0; i4<ni4; i4++){
			const int i = i4*4;
			const V2R8 xi0 = xibuf[i+0][0];
			const V2R8 yi0 = xibuf[i+0][1];
			const V2R8 zi0 = xibuf[i+0][2];
			V2R8 ax0 = immr8(0.0);
			V2R8 ay0 = immr8(0.0);
			V2R8 az0 = immr8(0.0);
			V2R8  p0 = immr8(0.0);
			
			const V2R8 xi1 = xibuf[i+1][0];
			const V2R8 yi1 = xibuf[i+1][1];
			const V2R8 zi1 = xibuf[i+1][2];
			V2R8 ax1 = immr8(0.0);
			V2R8 ay1 = immr8(0.0);
			V2R8 az1 = immr8(0.0);
			V2R8  p1 = immr8(0.0);

			const V2R8 xi2 = xibuf[i+2][0];
			const V2R8 yi2 = xibuf[i+2][1];
			const V2R8 zi2 = xibuf[i+2][2];
			V2R8 ax2 = immr8(0.0);
			V2R8 ay2 = immr8(0.0);
			V2R8 az2 = immr8(0.0);
			V2R8  p2 = immr8(0.0);
			
			const V2R8 xi3 = xibuf[i+3][0];
			const V2R8 yi3 = xibuf[i+3][1];
			const V2R8 zi3 = xibuf[i+3][2];
			V2R8 ax3 = immr8(0.0);
			V2R8 ay3 = immr8(0.0);
			V2R8 az3 = immr8(0.0);
			V2R8  p3 = immr8(0.0);
			
			for(int j=0; j<nj2; j++){
				PREFETCH(xjbuf[j+1]);
				const V2R8 xj = xjbuf[j][0];
				const V2R8 yj = xjbuf[j][1];
				const V2R8 zj = xjbuf[j][2];
				const V2R8 mj = xjbuf[j][3];

				gravity_kernel(xi0, yi0, zi0, veps2, xj, yj, zj, mj, ax0, ay0, az0, p0);
				gravity_kernel(xi1, yi1, zi1, veps2, xj, yj, zj, mj, ax1, ay1, az1, p1);
				gravity_kernel(xi2, yi2, zi2, veps2, xj, yj, zj, mj, ax2, ay2, az2, p2);
				gravity_kernel(xi3, yi3, zi3, veps2, xj, yj, zj, mj, ax3, ay3, az3, p3);
			}

			aibuf[i+0][0] = ax0;
			aibuf[i+0][1] = ay0;
			aibuf[i+0][2] = az0;
			pibuf[i+0]    = p0;
			
			aibuf[i+1][0] = ax1;
			aibuf[i+1][1] = ay1;
			aibuf[i+1][2] = az1;
			pibuf[i+1]    = p1;

			aibuf[i+2][0] = ax2;
			aibuf[i+2][1] = ay2;
			aibuf[i+2][2] = az2;
			pibuf[i+2]    = p2;

			aibuf[i+3][0] = ax3;
			aibuf[i+3][1] = ay3;
			aibuf[i+3][2] = az3;
			pibuf[i+3]    = p3;
		}
	}
                    
};
