//#include "class.hpp"
//#include "force.hpp"
#include<particle_simulator.hpp>
#include "cuda_pointer.h"
#include "force_gpu_cuda.hpp"

enum{
	N_THREAD_GPU = 32,
	N_WALK_LIMIT = 1000,
	NI_LIMIT     = N_WALK_LIMIT*1000,
	NJ_LIMIT     = N_WALK_LIMIT*10000,
};

struct EpiGPU{
	float3 pos;
	int    id_walk;
};

struct EpjGPU{
	float4 posm;
};

struct ForceGPU{
	float4 accp;
};

inline __device__ float4 dev_gravity(
		float  eps2,
		float3 ipos,
		float4 jposm,
		float4 accp)
{
	float dx = jposm.x - ipos.x;
	float dy = jposm.y - ipos.y;
	float dz = jposm.z - ipos.z;

	float r2   = eps2 + dx*dx + dy*dy + dz*dz;
	float rinv = rsqrtf(r2);
	float pij  = jposm.w * rinv;
	float mri3 = rinv*rinv * pij;

	accp.x += mri3 * dx;
	accp.y += mri3 * dy;
	accp.z += mri3 * dz;
	accp.w -= pij;

	return accp;
}

#if 0
__global__ void ForceKernel(
		const int2   * ij_disp,
		const EpiGPU * epi,
		const EpjGPU * epj, 
		ForceGPU     * force,
		const float    eps2)
{
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
	const float3 ipos = epi[tid].pos;
    const int j_head = ij_disp[epi[tid].id_walk  ].y;
    const int j_tail = ij_disp[epi[tid].id_walk+1].y;

	float4 accp = make_float4(0.f, 0.f, 0.f, 0.f);
    for(int j=j_head; j<j_tail; j++){
		float4 jposm = epj[j].posm;
		accp = dev_gravity(eps2, ipos, jposm, accp);
	}

	force[tid].accp = accp;
}
#else
__device__ float4 ForceKernel_1walk(
		float4       *jpsh,
		const float3  ipos,
		const int     id_walk,
		const int2   *ij_disp,
		const EpjGPU *epj, 
		float4        accp,
		const float   eps2)
{
    const int tid = threadIdx.x;
    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;

	for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
		// __syncthreads();
		jpsh[tid] = ((float4 *)(epj + j)) [tid];
		// __syncthreads();

		if(j_tail-j < N_THREAD_GPU){
			for(int jj=0; jj<j_tail-j; jj++){
				accp = dev_gravity(eps2, ipos, jpsh[jj], accp);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				accp = dev_gravity(eps2, ipos, jpsh[jj], accp);
			}
		}
	}
	
	return accp;
}

__device__ float4 ForceKernel_2walk(
		float4        jpsh[2][N_THREAD_GPU],
		const float3  ipos,
		const int     id_walk,
		const int     iwalk0,
		const int     iwalk1,
		const int2   *ij_disp,
		const EpjGPU *epj, 
		float4        accp,
		const float   eps2)
{
	const int jbeg0 = ij_disp[iwalk0].y;
	const int jbeg1 = ij_disp[iwalk1].y;
	const int jend0 = ij_disp[iwalk0 + 1].y;
	const int jend1 = ij_disp[iwalk1 + 1].y;
	const int nj0   = jend0 - jbeg0;
	const int nj1   = jend1 - jbeg1;

	const int nj_longer  = nj0 > nj1 ? nj0 : nj1;
	const int nj_shorter = nj0 > nj1 ? nj1 : nj0;
	const int walk_longer= nj0 > nj1 ? 0 : 1;
	const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;

	const int mywalk = id_walk==iwalk0 ? 0 : 1;

    const int tid = threadIdx.x;
	for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
		jpsh[0][tid] = ((float4 *)(epj + jbeg0 + j)) [tid];
		jpsh[1][tid] = ((float4 *)(epj + jbeg1 + j)) [tid];
		if(nj_shorter-j < N_THREAD_GPU){
			for(int jj=0; jj<nj_shorter-j; jj++){
				accp = dev_gravity(eps2, ipos, jpsh[mywalk][jj], accp);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				accp = dev_gravity(eps2, ipos, jpsh[mywalk][jj], accp);
			}
		}
	}
	for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
		jpsh[0][tid] = ((float4 *)(epj + jbeg_longer +  j)) [tid];
		int jrem = nj_longer - j;
		if(jrem < N_THREAD_GPU){
			for(int jj=0; jj<jrem; jj++){
				if(mywalk == walk_longer)
				accp = dev_gravity(eps2, ipos, jpsh[0][jj], accp);
			}
		}else{
#pragma unroll
			for(int jj=0; jj<N_THREAD_GPU; jj++){
				if(mywalk == walk_longer)
				accp = dev_gravity(eps2, ipos, jpsh[0][jj], accp);
			}
		}
	}

	return accp;
}

__device__ float4 ForceKernel_multiwalk(
		const float3  ipos,
		const int     id_walk,
		const int2   *ij_disp,
		const EpjGPU *epj, 
		float4        accp,
		const float   eps2)
{
    const int j_head = ij_disp[id_walk  ].y;
    const int j_tail = ij_disp[id_walk+1].y;

#if 1
    for(int j=j_head; j<j_tail; j++){
		float4 jposm = epj[j].posm;
		accp = dev_gravity(eps2, ipos, jposm, accp);
	}
#else
	int njmin = j_tail - j_head;
	njmin = min(njmin, __shfl_xor(njmin, 1));
	njmin = min(njmin, __shfl_xor(njmin, 2));
	njmin = min(njmin, __shfl_xor(njmin, 4));
	njmin = min(njmin, __shfl_xor(njmin, 8));
	njmin = min(njmin, __shfl_xor(njmin, 16));
	
	njmin &= 3;;
	for(int j=0; j<njmin; j+=4){
#pragma unroll 4
		for(int jj=0; jj<4; jj++){
			float4 jposm = epj[j_head + j + jj].posm;
			float4 jposm = jpf[jj];
			accp = dev_gravity(eps2, ipos, jposm, accp);
		}
	}
    for(int j=j_head+njmin; j<j_tail; j++){
		float4 jposm = epj[j].posm;
		accp = dev_gravity(eps2, ipos, jposm, accp);
	}
#endif
	return accp;
}

__global__ void ForceKernel(
		const int2   * ij_disp,
		const EpiGPU * epi,
		const EpjGPU * epj, 
		ForceGPU     * force,
		const float    eps2)
{
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
	float3 ipos    = epi[tid].pos;
	int    id_walk = epi[tid].id_walk;
	float4 accp    = make_float4(0.f, 0.f, 0.f, 0.f);


	int t_head = blockDim.x * blockIdx.x;
	int t_tail = t_head + N_THREAD_GPU - 1;
	int nwalk_in_block = 1 + (epi[t_tail].id_walk - epi[t_head].id_walk);

	__shared__ float4 jpsh[2][N_THREAD_GPU];

	if(1 == nwalk_in_block){
		accp = ForceKernel_1walk(jpsh[0], ipos, id_walk, ij_disp, epj, accp, eps2);
	} else if(2 == nwalk_in_block){
		// accp = ForceKernel_multiwalk(ipos, id_walk, ij_disp, epj, accp, eps2);
		int iwalk0 = epi[t_head].id_walk;
		int iwalk1 = epi[t_tail].id_walk;
		accp = ForceKernel_2walk(jpsh, ipos, id_walk, iwalk0, iwalk1, ij_disp, epj, accp, eps2);
	} else{
		accp = ForceKernel_multiwalk(ipos, id_walk, ij_disp, epj, accp, eps2);
	}
	force[tid].accp = accp;
}
#endif

static cudaPointer<EpiGPU>   dev_epi;
static cudaPointer<EpjGPU>   dev_epj;
static cudaPointer<ForceGPU> dev_force;
static cudaPointer<int2>     ij_disp;
static bool init_call = true;

PS::S32 DispatchKernelWithSP(
                             const PS::S32          tag,
                             const PS::S32          n_walk,
                             const FPGrav          *epi[],
                             const PS::S32          n_epi[],
                             const FPGrav          *epj[],
                             const PS::S32          n_epj[],
                             const PS::SPJMonopole *spj[],
                             const PS::S32          n_spj[]){
    assert(n_walk <= N_WALK_LIMIT);
    if(init_call){
		dev_epi  .allocate(NI_LIMIT);
		dev_epj  .allocate(NJ_LIMIT);
		dev_force.allocate(NI_LIMIT);
		ij_disp  .allocate(N_WALK_LIMIT+2);
		init_call = false;
    }
    const float eps2 = FPGrav::eps * FPGrav::eps;
    ij_disp[0].x = 0;
    ij_disp[0].y = 0;
    for(int k=0; k<n_walk; k++){
        ij_disp[k+1].x = ij_disp[k].x + n_epi[k];
        ij_disp[k+1].y = ij_disp[k].y + (n_epj[k] + n_spj[k]);
    }
    ij_disp[n_walk+1] = ij_disp[n_walk];

    assert(ij_disp[n_walk].x < NI_LIMIT);
    assert(ij_disp[n_walk].y < NJ_LIMIT);
    ij_disp.htod(n_walk + 2);

    int ni_tot_reg = ij_disp[n_walk].x;
    if(ni_tot_reg % N_THREAD_GPU){
        ni_tot_reg /= N_THREAD_GPU;
        ni_tot_reg++;
        ni_tot_reg *= N_THREAD_GPU;
    }

    int ni_tot = 0;
    int nj_tot = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<n_epi[iw]; i++){
            dev_epi[ni_tot].pos.x = epi[iw][i].pos.x;
            dev_epi[ni_tot].pos.y = epi[iw][i].pos.y;
            dev_epi[ni_tot].pos.z = epi[iw][i].pos.z;
            dev_epi[ni_tot].id_walk = iw;
            ni_tot++;
        }
        for(int j=0; j<n_epj[iw]; j++){
            dev_epj[nj_tot].posm.x  = epj[iw][j].pos.x;
            dev_epj[nj_tot].posm.y  = epj[iw][j].pos.y;
            dev_epj[nj_tot].posm.z  = epj[iw][j].pos.z;
            dev_epj[nj_tot].posm.w  = epj[iw][j].mass;
            nj_tot++;
        }
        for(int j=0; j<n_spj[iw]; j++){
            dev_epj[nj_tot].posm.x  = spj[iw][j].pos.x;
            dev_epj[nj_tot].posm.y  = spj[iw][j].pos.y;
            dev_epj[nj_tot].posm.z  = spj[iw][j].pos.z;
            dev_epj[nj_tot].posm.w  = spj[iw][j].getCharge();
            nj_tot++;
        }
    }
    for(int i=ni_tot; i<ni_tot_reg; i++){
        dev_epi[i].id_walk = n_walk;
    }

    dev_epi.htod(ni_tot_reg);
    dev_epj.htod(nj_tot);

    int nblocks  = ni_tot_reg / N_THREAD_GPU;
    int nthreads = N_THREAD_GPU;
    ForceKernel <<<nblocks, nthreads>>> (ij_disp, dev_epi, dev_epj, dev_force, eps2);

    return 0;
}

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       FPGrav    *force[])
{
    int ni_tot = 0;
    for(int k=0; k<n_walk; k++){
        ni_tot += ni[k];
    }
    dev_force.dtoh(ni_tot);

    int n_cnt = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<ni[iw]; i++){
            force[iw][i].acc.x = dev_force[n_cnt].accp.x;
            force[iw][i].acc.y = dev_force[n_cnt].accp.y;
            force[iw][i].acc.z = dev_force[n_cnt].accp.z;
            force[iw][i].pot   = dev_force[n_cnt].accp.w;
            n_cnt++;
        }
    }
    return 0;
}
