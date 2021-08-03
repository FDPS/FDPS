#pragma once

#define SIMD_WIDTH 64 // 512bit

struct Epi0{
    PS::F32vec pos;
    PS::F32vec vel;
    PS::F32    mass;
    PS::F32    r_coll;
    PS::S32    id;
};
struct Epj0{
    PS::F32vec pos;
    PS::F32vec vel;
    PS::F32    mass;
    PS::F32    r_coll;
    PS::S32    id;
};
struct Force0{
    PS::F32vec acc;
    PS::F32vec acc_dash;
    PS::F32    pot;
};
struct Epi1{
    PS::F32vec pos;
};
struct Epj1{
    PS::F32vec pos;
    PS::F32    mass;
};
struct Force1{
    PS::F32vec acc;
    PS::F32    pot;
};

struct Epi2{
    PS::F32vec pos;
};
struct Epj2{
    PS::F32vec pos;
    PS::F32    mass;
    PS::F32    quad_xx;
    PS::F32    quad_yy;
    PS::F32    quad_zz;
    PS::F32    quad_xy;
    PS::F32    quad_xz;
    PS::F32    quad_yz;

};
struct Force2{
    PS::F32vec acc;
    PS::F32    pot;
};

#include"kernel_ep.hpp"

template<typename Tpi, typename Tpj, typename Tfi>
struct CalcForceEp{
    void operator () (const Tpi * ep_i,
		      const PS::S32 n_ip,
		      const Tpj * ep_j,
		      const PS::S32 n_jp,
		      Tfi * force) {
	const PS::F32 eps2  = (PS::F32)(FP_t::eps*FP_t::eps);
	const PS::F32 kappa = (PS::F32)FP_t::kappa;
	const PS::F32 eta   = (PS::F32)FP_t::eta;  
	alignas(SIMD_WIDTH) Epi0 epi[n_ip];
	alignas(SIMD_WIDTH) Force0 f[n_ip];
	for(int i=0;i<n_ip;i++){
	    epi[i].pos = ep_i[i].pos_car - ep_i[0].pos_car;
	    epi[i].mass   = ep_i[i].mass;
	    epi[i].r_coll = ep_i[i].r_coll;
	    epi[i].vel    = ep_i[i].vel_full;
	    epi[i].id     = ep_i[i].id;

	    f[i].acc      = 0.0;
	    f[i].acc_dash = 0.0;
	    f[i].pot      = 0.0;
	    
	    force[i].acc      = 0.0;
	    force[i].acc_dash = 0.0;
	    force[i].pot      = 0.0;
	}
	alignas(SIMD_WIDTH) Epj0 epj[n_jp];
	for(int i=0;i<n_jp;i++){
	    epj[i].pos = ep_j[i].pos_car - ep_i[0].pos_car;
	    epj[i].mass   = ep_j[i].mass;
	    epj[i].r_coll = ep_j[i].r_coll;
	    epj[i].vel    = ep_j[i].vel_full;
	    epj[i].id     = ep_j[i].id;
	}
	CalcForceEpEpImpl{eps2, kappa, eta}(epi,n_ip,epj,n_jp,f);
	for(int i=0;i<n_ip;i++){
	    force[i].acc      = f[i].acc;
	    force[i].acc_dash = f[i].acc_dash;
	    force[i].pot      = f[i].pot;
	}
    }
};


#include"kernel_sp.hpp"

template<typename Tpi, typename Tpj, typename Tfi>
struct CalcForceSpMono{
    void operator () (const Tpi * ep_i,
		      const PS::S32 n_ip,
		      const Tpj * ep_j,
		      const PS::S32 n_jp,
		      Tfi * force) {
	PS::F32 eps2  = FP_t::eps*FP_t::eps;
	alignas(SIMD_WIDTH) Epi1 epi[n_ip];
	alignas(SIMD_WIDTH) Force1 f[n_ip];
	for(int i=0;i<n_ip;i++){
	    epi[i].pos = ep_i[i].pos_car - ep_i[0].pos_car;
	    f[i].acc = 0.0;
	    f[i].pot = 0.0;
	}
	alignas(SIMD_WIDTH) Epj1 epj[n_jp];
	for(int i=0;i<n_jp;i++){
	    epj[i].pos = ep_j[i].pos_car - ep_i[0].pos_car;
	    epj[i].mass = ep_j[i].mass;
	}
	CalcForceEpSpImpl{eps2}(epi,n_ip,epj,n_jp,f);
	for(int i=0;i<n_ip;i++){
	    force[i].acc += f[i].acc;
	    force[i].pot += f[i].pot;
	}
    }
};


#include"kernel_sp_quad.hpp"

template<typename Tpi, typename Tpj, typename Tfi>
struct CalcForceSpQuad{
    void operator () (const Tpi * ep_i,
		      const PS::S32 n_ip,
		      const Tpj * ep_j,
		      const PS::S32 n_jp,
		      Tfi * force) {
        PS::F32 eps2  = FP_t::eps*FP_t::eps;
	alignas(SIMD_WIDTH) Epi2 epi[n_ip];
	alignas(SIMD_WIDTH) Force2 f[n_ip];
	for(int i=0;i<n_ip;i++){
	    epi[i].pos = ep_i[i].pos_car - ep_i[0].pos_car;
	    f[i].acc = 0.0;
	    f[i].pot = 0.0;
	}
	alignas(SIMD_WIDTH) Epj2 epj[n_jp];
	for(int i=0;i<n_jp;i++){
	    epj[i].pos     = ep_j[i].pos_car - ep_i[0].pos_car;
	    epj[i].mass    = ep_j[i].mass;
	    epj[i].quad_xx = ep_j[i].quad.xx;
	    epj[i].quad_yy = ep_j[i].quad.yy;
	    epj[i].quad_zz = ep_j[i].quad.zz;
	    epj[i].quad_xy = ep_j[i].quad.xy;
	    epj[i].quad_xz = ep_j[i].quad.xz;
	    epj[i].quad_yz = ep_j[i].quad.yz;
	}
	CalcForceEpSpQuadImpl{eps2}(epi,n_ip,epj,n_jp,f);
	for(int i=0;i<n_ip;i++){
	    force[i].acc += f[i].acc;
	    force[i].pot += f[i].pot;
	}
    }
};

