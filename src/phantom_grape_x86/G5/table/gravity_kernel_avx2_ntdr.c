enum{
	EXP_BIT  = 4,
	FRC_BIT  = 5,
};

typedef struct iptdata{
	float x[4];
	float y[4];
	float z[4];
	float eps2[4]; // not used in this implementation
} Ipdata, *pIpdata;

typedef struct fodata{
	float ax[4];
	float ay[4];
	float az[4];
	float phi[4];
} Fodata, *pFodata;

typedef struct jpdata{
	float x, y, z, m;
} Jpdata, *pJpdata;

typedef float     v4sf __attribute__ ((vector_size(16)));
typedef double    v2df __attribute__ ((vector_size(16)));
typedef int       v4si __attribute__ ((vector_size(16)));
typedef long long v2di __attribute__ ((vector_size(16)));
typedef float     v8sf __attribute__ ((vector_size(32)));
typedef double    v4df __attribute__ ((vector_size(32)));
typedef int       v8si __attribute__ ((vector_size(32)));
typedef long long v4di __attribute__ ((vector_size(32)));

#define REP4(x) {x, x, x, x}
#define REP8(x) {x, x, x, x, x, x, x, x}

void gravity_kernel(pIpdata ipdata, pJpdata jp, pFodata fodata, 
		int nj, float fcut[][2], v4sf r2cut, v4sf accscale);

static inline void vec_extract(v8si vec, int idx[]){
#if 1
# if 0
	v2di lo = (v2di)__builtin_ia32_vextractf128_ps256(vec, 0);
	v2di hi = (v2di)__builtin_ia32_vextractf128_ps256(vec, 1);
	unsigned long long ll01 = __builtin_ia32_vec_ext_v2di (lo, 0);
	unsigned long long ll23 = __builtin_ia32_vec_ext_v2di (lo, 1);
	unsigned long long ll45 = __builtin_ia32_vec_ext_v2di (hi, 0);
	unsigned long long ll67 = __builtin_ia32_vec_ext_v2di (hi, 1);
# else
	unsigned long long buf[4] __attribute__((aligned(32)));
	*(v8si *)buf = vec;
	unsigned long long ll01 = buf[0];
	unsigned long long ll23 = buf[1];
	unsigned long long ll45 = buf[2];
	unsigned long long ll67 = buf[3];
# endif
	idx[0] = ll01;
	idx[1] = ll01 >> 32;
	idx[2] = ll23;
	idx[3] = ll23 >> 32;
	idx[4] = ll45;
	idx[5] = ll45 >> 32;
	idx[6] = ll67;
	idx[7] = ll67 >> 32;
#else
	v4si lo = (v4si)__builtin_ia32_vextractf128_ps256(vec, 0);
	v4si hi = (v4si)__builtin_ia32_vextractf128_ps256(vec, 1);
	idx[0] = __builtin_ia32_vec_ext_v4si(lo, 0);
	idx[1] = __builtin_ia32_vec_ext_v4si(lo, 1);
	idx[2] = __builtin_ia32_vec_ext_v4si(lo, 2);
	idx[3] = __builtin_ia32_vec_ext_v4si(lo, 3);
	idx[4] = __builtin_ia32_vec_ext_v4si(hi, 0);
	idx[5] = __builtin_ia32_vec_ext_v4si(hi, 1);
	idx[6] = __builtin_ia32_vec_ext_v4si(hi, 2);
	idx[7] = __builtin_ia32_vec_ext_v4si(hi, 3);
#endif
}

void gravity_kernel(
		pIpdata ipdata, 
		pJpdata jpdata, 
		pFodata fodata, 
		int     nj, 
		float   fcut[][2], 
		v4sf    v4_r2cut, 
		v4sf    accscale)
{
	fcut -= (1<<(30-(23-FRC_BIT)));

	v8sf xi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(ipdata->x));
	v8sf yi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(ipdata->y));
	v8sf zi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(ipdata->z));
	v8sf r2cut = REP8(0.0f);
	r2cut = __builtin_ia32_vinsertf128_ps256(r2cut, v4_r2cut, 0);
	r2cut = __builtin_ia32_vinsertf128_ps256(r2cut, v4_r2cut, 1);
	v8sf two = REP8(2.0f);
	v8sf ax = REP8(0.0f);
	v8sf ay = REP8(0.0f);
	v8sf az = REP8(0.0f);

	v8sf cc = REP8(0.0f);

	int j;
	v8sf jp = *(v8sf *)(jpdata + 0);
	v8sf xj = __builtin_ia32_shufps256(jp, jp, 0x00);
	v8sf yj = __builtin_ia32_shufps256(jp, jp, 0x55);
	v8sf zj = __builtin_ia32_shufps256(jp, jp, 0xaa);
	v8sf mj = __builtin_ia32_shufps256(jp, jp, 0xff);
	/*
	 * ~22 cycle per loop on Core-i7 4770 with TB disabled
	 */
	for(j=0; j<nj; j+=2){
		jp = *(v8sf *)(jpdata+=2);

		v8sf dx = xj - xi;
		v8sf dy = yj - yi;
		v8sf dz = zj - zi;
#if 1
		v8sf r2 = ((two + dx*dx) + dy*dy) + dz*dz;
#else
		v8sf r2 = (two + dx*dx) + (dy*dy + dz*dz);
#endif
#if 0
		v8sf mask_v8 = r2 < r2cut;
		v8si r2_sr = __builtin_ia32_psrldi256((v8si)r2, 23-FRC_BIT);
		v8si r2_sl = __builtin_ia32_pslldi256(r2_sr,    23-FRC_BIT);
		v8si izero = REP8(0);
		v4di idx_0145 = (v4di)__builtin_ia32_punpckldq256(r2_sr, izero);
		v4di idx_2367 = (v4di)__builtin_ia32_punpckhdq256(r2_sr, izero);
		v4di msk_0145 = (v4di)__builtin_ia32_unpcklps256(mask_v8, mask_v8);
		v4di msk_2367 = (v4di)__builtin_ia32_unpckhps256(mask_v8, mask_v8);
		v4di lzero = REP4(0);
		v4di tbl_0145 = __builtin_ia32_gatherdiv4di(lzero, (long long *)fcut, idx_0145, msk_0145, 8);
		v4di tbl_2367 = __builtin_ia32_gatherdiv4di(lzero, (long long *)fcut, idx_2367, msk_2367, 8);
#else
		r2 = __builtin_ia32_minps256(r2, r2cut);
		v8si r2_sr = __builtin_ia32_psrldi256((v8si)r2, 23-FRC_BIT);
		v8si r2_sl = __builtin_ia32_pslldi256(r2_sr,    23-FRC_BIT);
		unsigned int idx[8] __attribute__((aligned(32)));
# if 1
		*(v8si *)idx = r2_sr;
# else
		vec_extract(r2_sr, idx);
# endif
# if 1
		const long long *ptr = (long long *)fcut;
		v4di tbl_0145 = {ptr[idx[0]], ptr[idx[1]], ptr[idx[4]], ptr[idx[5]]};
		v4di tbl_2367 = {ptr[idx[2]], ptr[idx[3]], ptr[idx[6]], ptr[idx[7]]};
# else
		const double *ptr = (double *)fcut;
		v4df tbl_0145 = {ptr[idx[0]], ptr[idx[1]], ptr[idx[4]], ptr[idx[5]]};
		v4df tbl_2367 = {ptr[idx[2]], ptr[idx[3]], ptr[idx[6]], ptr[idx[7]]};
# endif
#endif

		v8sf ff = __builtin_ia32_shufps256((v8sf)tbl_0145, (v8sf)tbl_2367, 0x88);
		v8sf df = __builtin_ia32_shufps256((v8sf)tbl_0145, (v8sf)tbl_2367, 0xdd);
		v8sf dr2 = r2 - (v8sf)r2_sl;
		ff += dr2 * df;

		v8sf mf = mj * ff;
		xj = __builtin_ia32_shufps256(jp, jp, 0x00);
		yj = __builtin_ia32_shufps256(jp, jp, 0x55);
		zj = __builtin_ia32_shufps256(jp, jp, 0xaa);
		mj = __builtin_ia32_shufps256(jp, jp, 0xff);
		ax += mf * dx;
		ay += mf * dy;
		az += mf * dz;

		// cc += __builtin_ia32_andps256(mask_v8, (v8sf)REP8(1.0f));
	}
	v4sf ax4 = __builtin_ia32_vextractf128_ps256(ax, 0)
	         + __builtin_ia32_vextractf128_ps256(ax, 1);
	v4sf ay4 = __builtin_ia32_vextractf128_ps256(ay, 0)
	         + __builtin_ia32_vextractf128_ps256(ay, 1);
	v4sf az4 = __builtin_ia32_vextractf128_ps256(az, 0)
	         + __builtin_ia32_vextractf128_ps256(az, 1);
	*(v4sf *)(fodata->ax) = accscale * ax4;
	*(v4sf *)(fodata->ay) = accscale * ay4;
	*(v4sf *)(fodata->az) = accscale * az4;

	v4sf cc4 = __builtin_ia32_vextractf128_ps256(cc, 0)
	         + __builtin_ia32_vextractf128_ps256(cc, 1);
	*(v4sf *)(fodata->phi) = cc4;
}

