#include<pikg_vector.hpp>
#include<cmath>
#include<limits>
#include<chrono>

#include <pikg_avx2.hpp>
struct CalcForceEpEpImpl{
PIKG::F32 eps2;
PIKG::F32 kappa;
PIKG::F32 eta;
CalcForceEpEpImpl(){}
CalcForceEpEpImpl(PIKG::F32 eps2,PIKG::F32 kappa,PIKG::F32 eta):eps2(eps2),kappa(kappa),eta(eta){}
void initialize(PIKG::F32 eps2_,PIKG::F32 kappa_,PIKG::F32 eta_){
eps2 = eps2_;
kappa = kappa_;
eta = eta_;
}
int kernel_id = 0;
void operator()(const Epi0* __restrict__ epi,const int ni,const Epj0* __restrict__ epj,const int nj,Force0* __restrict__ force,const int kernel_select = 1){
static_assert(sizeof(Epi0) == 36,"check consistency of EPI member variable definition between PIKG source and original source");
static_assert(sizeof(Epj0) == 36,"check consistency of EPJ member variable definition between PIKG source and original source");
static_assert(sizeof(Force0) == 28,"check consistency of FORCE member variable definition between PIKG source and original source");
if(kernel_select>=0) kernel_id = kernel_select;
if(kernel_id == 0){
std::cout << "ni: " << ni << " nj:" << nj << std::endl;
Force0* force_tmp = new Force0[ni];
std::chrono::system_clock::time_point  start, end;
double min_time = std::numeric_limits<double>::max();
{ // test Kernel_I8_J1
for(int i=0;i<ni;i++) force_tmp[i] = force[i];
start = std::chrono::system_clock::now();
Kernel_I8_J1(epi,ni,epj,nj,force_tmp);
end = std::chrono::system_clock::now();
double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
std::cerr << "kerel 1: " << elapsed << " ns" << std::endl;
if(min_time > elapsed){
min_time = elapsed;
kernel_id = 1;
}
}
{ // test Kernel_I1_J8
for(int i=0;i<ni;i++) force_tmp[i] = force[i];
start = std::chrono::system_clock::now();
Kernel_I1_J8(epi,ni,epj,nj,force_tmp);
end = std::chrono::system_clock::now();
double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
std::cerr << "kerel 2: " << elapsed << " ns" << std::endl;
if(min_time > elapsed){
min_time = elapsed;
kernel_id = 2;
}
}
delete[] force_tmp;
} // if(kernel_id == 0)
if(kernel_id == 1) Kernel_I8_J1(epi,ni,epj,nj,force);
if(kernel_id == 2) Kernel_I1_J8(epi,ni,epj,nj,force);
} // operator() definition 
void Kernel_I8_J1(const Epi0* __restrict__ epi,const PIKG::S32 ni,const Epj0* __restrict__ epj,const PIKG::S32 nj,Force0* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
for(i = 0;i < (ni/8)*8;i += 8){
__m256i EPI_id;

int index_gather_load0[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load0 = _mm256_load_si256((const __m256i*)index_gather_load0);
EPI_id = _mm256_i32gather_epi32(((int*)&epi[i+0].id),vindex_gather_load0,4);
__m256 EPI_mass;

int index_gather_load1[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load1 = _mm256_load_si256((const __m256i*)index_gather_load1);
EPI_mass = _mm256_i32gather_ps(((float*)&epi[i+0].mass),vindex_gather_load1,4);
__m256x3 EPI_pos;

int index_gather_load2[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load2 = _mm256_load_si256((const __m256i*)index_gather_load2);
EPI_pos.v0 = _mm256_i32gather_ps(((float*)&epi[i+0].pos.x),vindex_gather_load2,4);
int index_gather_load3[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load3 = _mm256_load_si256((const __m256i*)index_gather_load3);
EPI_pos.v1 = _mm256_i32gather_ps(((float*)&epi[i+0].pos.y),vindex_gather_load3,4);
int index_gather_load4[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load4 = _mm256_load_si256((const __m256i*)index_gather_load4);
EPI_pos.v2 = _mm256_i32gather_ps(((float*)&epi[i+0].pos.z),vindex_gather_load4,4);
__m256 EPI_r_coll;

int index_gather_load5[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load5 = _mm256_load_si256((const __m256i*)index_gather_load5);
EPI_r_coll = _mm256_i32gather_ps(((float*)&epi[i+0].r_coll),vindex_gather_load5,4);
__m256x3 EPI_vel;

int index_gather_load6[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load6 = _mm256_load_si256((const __m256i*)index_gather_load6);
EPI_vel.v0 = _mm256_i32gather_ps(((float*)&epi[i+0].vel.x),vindex_gather_load6,4);
int index_gather_load7[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load7 = _mm256_load_si256((const __m256i*)index_gather_load7);
EPI_vel.v1 = _mm256_i32gather_ps(((float*)&epi[i+0].vel.y),vindex_gather_load7,4);
int index_gather_load8[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load8 = _mm256_load_si256((const __m256i*)index_gather_load8);
EPI_vel.v2 = _mm256_i32gather_ps(((float*)&epi[i+0].vel.z),vindex_gather_load8,4);
__m256x3 FORCE_acc;

FORCE_acc.v0 = _mm256_set1_ps(0.0f);
FORCE_acc.v1 = _mm256_set1_ps(0.0f);
FORCE_acc.v2 = _mm256_set1_ps(0.0f);
__m256x3 FORCE_acc_dash;

FORCE_acc_dash.v0 = _mm256_set1_ps(0.0f);
FORCE_acc_dash.v1 = _mm256_set1_ps(0.0f);
FORCE_acc_dash.v2 = _mm256_set1_ps(0.0f);
__m256 FORCE_pot;

FORCE_pot = _mm256_set1_ps(0.0f);
for(j = 0;j < (nj/1)*1;++j){
__m256i EPJ_id;

EPJ_id = _mm256_set1_epi32(epj[j+0].id);

__m256 EPJ_mass;

EPJ_mass = _mm256_set1_ps(epj[j+0].mass);

__m256x3 EPJ_pos;

EPJ_pos.v0 = _mm256_set1_ps(epj[j+0].pos.x);

EPJ_pos.v1 = _mm256_set1_ps(epj[j+0].pos.y);

EPJ_pos.v2 = _mm256_set1_ps(epj[j+0].pos.z);

__m256 EPJ_r_coll;

EPJ_r_coll = _mm256_set1_ps(epj[j+0].r_coll);

__m256x3 EPJ_vel;

EPJ_vel.v0 = _mm256_set1_ps(epj[j+0].vel.x);

EPJ_vel.v1 = _mm256_set1_ps(epj[j+0].vel.y);

EPJ_vel.v2 = _mm256_set1_ps(epj[j+0].vel.z);

__m256x3 acc_sprg_tmp;

__m256x3 acc_dash_tmp;

__m256 pot_sprg_tmp;

__m256 pot_grav_tmp;

__m256x3 rij;

__m256 r_coll_tmp;

__m256 r_coll_sq;

__m256 r_real_sq;

__m256 over_r_real;

__m256 over_r_real_sq;

__m256 r_coll_cu;

__m256 over_r_coll_cu;

__m256 __fkg_tmp0;

__m256x3 __fkg_tmp4;

__m256 pot_offset;

__m256 __fkg_tmp5;

__m256 m_red;

__m256 r_real;

__m256 dr;

__m256 __fkg_tmp1;

__m256x3 __fkg_tmp6;

__m256 __fkg_tmp7;

__m256x3 vij;

__m256 rv;

__m256 __fkg_tmp2;

__m256x3 __fkg_tmp8;

__m256x3 acc_grav_tmp;

__m256 m_over_r_real;

__m256 __fkg_tmp3;

__m256 __fkg_tmp10;

__m256 __fkg_tmp9;

__m256 __fkg_tmp12;

__m256 __fkg_tmp11;

__m256 __fkg_tmp14;

__m256 __fkg_tmp13;

__m256 __fkg_tmp15;

acc_sprg_tmp.v0 = _mm256_set1_ps(0.0f);
acc_sprg_tmp.v1 = _mm256_set1_ps(0.0f);
acc_sprg_tmp.v2 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v0 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v1 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v2 = _mm256_set1_ps(0.0f);
pot_sprg_tmp = _mm256_set1_ps(0.0f);
pot_grav_tmp = _mm256_set1_ps(0.0f);
{
__m256 pg1;
__m256 pg0;
pg1 = _mm256_cmp_ps(_mm256_castsi256_ps(EPI_id),_mm256_castsi256_ps(EPJ_id),_CMP_NEQ_OQ);
pg0 = pg1;

rij.v0 = _mm256_sub_ps(EPI_pos.v0,EPJ_pos.v0);
rij.v1 = _mm256_sub_ps(EPI_pos.v1,EPJ_pos.v1);
rij.v2 = _mm256_sub_ps(EPI_pos.v2,EPJ_pos.v2);
r_coll_tmp = _mm256_add_ps(EPI_r_coll,EPJ_r_coll);
r_coll_sq = _mm256_mul_ps(r_coll_tmp,r_coll_tmp);
r_real_sq = _mm256_fmadd_ps(rij.v2,rij.v2,_mm256_fmadd_ps(rij.v1,rij.v1,_mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2))));
over_r_real = rsqrt(r_real_sq);
over_r_real_sq = _mm256_mul_ps(over_r_real,over_r_real);
{
__m256 pg3;
__m256 pg2;
pg3 = _mm256_cmp_ps(r_coll_sq,r_real_sq,_CMP_GT_OS);
pg2 = pg3;
pg3 = _mm256_and_ps(pg3,pg1);

r_coll_cu = _mm256_mul_ps(r_coll_sq,r_coll_tmp);
over_r_coll_cu = _mm256_div_ps(_mm256_set1_ps(1.0f),r_coll_cu);
__fkg_tmp0 = _mm256_mul_ps(EPJ_mass,over_r_coll_cu);
__fkg_tmp4.v0 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v0);
__fkg_tmp4.v1 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v1);
__fkg_tmp4.v2 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v2);
pot_offset = _mm256_div_ps(_mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),_mm256_set1_ps(1.5f)),r_coll_tmp);
__fkg_tmp5 = _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.25f),EPJ_mass),_mm256_fmadd_ps(r_real_sq,over_r_coll_cu,pot_offset));
m_red = _mm256_div_ps(EPJ_mass,_mm256_add_ps(EPI_mass,EPJ_mass));
r_real = _mm256_mul_ps(r_real_sq,over_r_real);
dr = _mm256_sub_ps(r_coll_tmp,r_real);
__fkg_tmp1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real);
__fkg_tmp6.v0 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v0);
__fkg_tmp6.v1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v1);
__fkg_tmp6.v2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v2);
__fkg_tmp7 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.25f),_mm256_set1_ps(kappa)),m_red),dr),dr);
vij.v0 = _mm256_sub_ps(EPI_vel.v0,EPJ_vel.v0);
vij.v1 = _mm256_sub_ps(EPI_vel.v1,EPJ_vel.v1);
vij.v2 = _mm256_sub_ps(EPI_vel.v2,EPJ_vel.v2);
rv = _mm256_fmadd_ps(rij.v2,vij.v2,_mm256_fmadd_ps(rij.v0,vij.v0,_mm256_mul_ps(rij.v1,vij.v1)));
__fkg_tmp2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq);
__fkg_tmp8.v0 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v0);
__fkg_tmp8.v1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v1);
__fkg_tmp8.v2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v2);
acc_grav_tmp.v0 = _mm256_blendv_ps(acc_grav_tmp.v0,__fkg_tmp4.v0,pg3);;
acc_grav_tmp.v1 = _mm256_blendv_ps(acc_grav_tmp.v1,__fkg_tmp4.v1,pg3);;
acc_grav_tmp.v2 = _mm256_blendv_ps(acc_grav_tmp.v2,__fkg_tmp4.v2,pg3);;
pot_grav_tmp = _mm256_blendv_ps(pot_grav_tmp,__fkg_tmp5,pg3);;
acc_sprg_tmp.v0 = _mm256_blendv_ps(acc_sprg_tmp.v0,__fkg_tmp6.v0,pg3);;
acc_sprg_tmp.v1 = _mm256_blendv_ps(acc_sprg_tmp.v1,__fkg_tmp6.v1,pg3);;
acc_sprg_tmp.v2 = _mm256_blendv_ps(acc_sprg_tmp.v2,__fkg_tmp6.v2,pg3);;
pot_sprg_tmp = _mm256_blendv_ps(pot_sprg_tmp,__fkg_tmp7,pg3);;
acc_dash_tmp.v0 = _mm256_blendv_ps(acc_dash_tmp.v0,__fkg_tmp8.v0,pg3);;
acc_dash_tmp.v1 = _mm256_blendv_ps(acc_dash_tmp.v1,__fkg_tmp8.v1,pg3);;
acc_dash_tmp.v2 = _mm256_blendv_ps(acc_dash_tmp.v2,__fkg_tmp8.v2,pg3);;
pg3 = _mm256_xor_ps(pg2,_mm256_castsi256_ps(_mm256_set1_epi32(-1)));
pg3 = _mm256_and_ps(pg3,pg1);

m_over_r_real = _mm256_mul_ps(EPJ_mass,over_r_real);
__fkg_tmp3 = _mm256_mul_ps(m_over_r_real,over_r_real_sq);
__fkg_tmp4.v0 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v0);
__fkg_tmp4.v1 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v1);
__fkg_tmp4.v2 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v2);
__fkg_tmp5 = _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),_mm256_set1_ps(0.5f)),m_over_r_real);
acc_grav_tmp.v0 = _mm256_blendv_ps(acc_grav_tmp.v0,__fkg_tmp4.v0,pg3);;
acc_grav_tmp.v1 = _mm256_blendv_ps(acc_grav_tmp.v1,__fkg_tmp4.v1,pg3);;
acc_grav_tmp.v2 = _mm256_blendv_ps(acc_grav_tmp.v2,__fkg_tmp4.v2,pg3);;
pot_grav_tmp = _mm256_blendv_ps(pot_grav_tmp,__fkg_tmp5,pg3);;
}

}

__fkg_tmp10 = _mm256_add_ps(acc_grav_tmp.v0,acc_sprg_tmp.v0);
__fkg_tmp9 = _mm256_add_ps(__fkg_tmp10,acc_dash_tmp.v0);
FORCE_acc.v0 = _mm256_add_ps(FORCE_acc.v0,__fkg_tmp9);
__fkg_tmp12 = _mm256_add_ps(acc_grav_tmp.v1,acc_sprg_tmp.v1);
__fkg_tmp11 = _mm256_add_ps(__fkg_tmp12,acc_dash_tmp.v1);
FORCE_acc.v1 = _mm256_add_ps(FORCE_acc.v1,__fkg_tmp11);
__fkg_tmp14 = _mm256_add_ps(acc_grav_tmp.v2,acc_sprg_tmp.v2);
__fkg_tmp13 = _mm256_add_ps(__fkg_tmp14,acc_dash_tmp.v2);
FORCE_acc.v2 = _mm256_add_ps(FORCE_acc.v2,__fkg_tmp13);
FORCE_acc_dash.v0 = _mm256_add_ps(FORCE_acc_dash.v0,acc_dash_tmp.v0);
FORCE_acc_dash.v1 = _mm256_add_ps(FORCE_acc_dash.v1,acc_dash_tmp.v1);
FORCE_acc_dash.v2 = _mm256_add_ps(FORCE_acc_dash.v2,acc_dash_tmp.v2);
__fkg_tmp15 = _mm256_add_ps(pot_grav_tmp,pot_sprg_tmp);
FORCE_pot = _mm256_add_ps(FORCE_pot,__fkg_tmp15);
} // loop of j

{
__m256 __fkg_tmp_accum;
int index_gather_load9[8] = {0,7,14,21,28,35,42,49};
__m256i vindex_gather_load9 = _mm256_load_si256((const __m256i*)index_gather_load9);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.x),vindex_gather_load9,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,FORCE_acc.v0);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc.x)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc.x)[7] = __fkg_store_tmp[1];
((float*)&force[i+0].acc.x)[14] = __fkg_store_tmp[2];
((float*)&force[i+0].acc.x)[21] = __fkg_store_tmp[3];
((float*)&force[i+0].acc.x)[28] = __fkg_store_tmp[4];
((float*)&force[i+0].acc.x)[35] = __fkg_store_tmp[5];
((float*)&force[i+0].acc.x)[42] = __fkg_store_tmp[6];
((float*)&force[i+0].acc.x)[49] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
int index_gather_load10[8] = {0,7,14,21,28,35,42,49};
__m256i vindex_gather_load10 = _mm256_load_si256((const __m256i*)index_gather_load10);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.y),vindex_gather_load10,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,FORCE_acc.v1);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc.y)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc.y)[7] = __fkg_store_tmp[1];
((float*)&force[i+0].acc.y)[14] = __fkg_store_tmp[2];
((float*)&force[i+0].acc.y)[21] = __fkg_store_tmp[3];
((float*)&force[i+0].acc.y)[28] = __fkg_store_tmp[4];
((float*)&force[i+0].acc.y)[35] = __fkg_store_tmp[5];
((float*)&force[i+0].acc.y)[42] = __fkg_store_tmp[6];
((float*)&force[i+0].acc.y)[49] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
int index_gather_load11[8] = {0,7,14,21,28,35,42,49};
__m256i vindex_gather_load11 = _mm256_load_si256((const __m256i*)index_gather_load11);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc.z),vindex_gather_load11,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,FORCE_acc.v2);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc.z)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc.z)[7] = __fkg_store_tmp[1];
((float*)&force[i+0].acc.z)[14] = __fkg_store_tmp[2];
((float*)&force[i+0].acc.z)[21] = __fkg_store_tmp[3];
((float*)&force[i+0].acc.z)[28] = __fkg_store_tmp[4];
((float*)&force[i+0].acc.z)[35] = __fkg_store_tmp[5];
((float*)&force[i+0].acc.z)[42] = __fkg_store_tmp[6];
((float*)&force[i+0].acc.z)[49] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
int index_gather_load12[8] = {0,7,14,21,28,35,42,49};
__m256i vindex_gather_load12 = _mm256_load_si256((const __m256i*)index_gather_load12);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc_dash.x),vindex_gather_load12,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,FORCE_acc_dash.v0);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc_dash.x)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc_dash.x)[7] = __fkg_store_tmp[1];
((float*)&force[i+0].acc_dash.x)[14] = __fkg_store_tmp[2];
((float*)&force[i+0].acc_dash.x)[21] = __fkg_store_tmp[3];
((float*)&force[i+0].acc_dash.x)[28] = __fkg_store_tmp[4];
((float*)&force[i+0].acc_dash.x)[35] = __fkg_store_tmp[5];
((float*)&force[i+0].acc_dash.x)[42] = __fkg_store_tmp[6];
((float*)&force[i+0].acc_dash.x)[49] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
int index_gather_load13[8] = {0,7,14,21,28,35,42,49};
__m256i vindex_gather_load13 = _mm256_load_si256((const __m256i*)index_gather_load13);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc_dash.y),vindex_gather_load13,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,FORCE_acc_dash.v1);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc_dash.y)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc_dash.y)[7] = __fkg_store_tmp[1];
((float*)&force[i+0].acc_dash.y)[14] = __fkg_store_tmp[2];
((float*)&force[i+0].acc_dash.y)[21] = __fkg_store_tmp[3];
((float*)&force[i+0].acc_dash.y)[28] = __fkg_store_tmp[4];
((float*)&force[i+0].acc_dash.y)[35] = __fkg_store_tmp[5];
((float*)&force[i+0].acc_dash.y)[42] = __fkg_store_tmp[6];
((float*)&force[i+0].acc_dash.y)[49] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
int index_gather_load14[8] = {0,7,14,21,28,35,42,49};
__m256i vindex_gather_load14 = _mm256_load_si256((const __m256i*)index_gather_load14);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].acc_dash.z),vindex_gather_load14,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,FORCE_acc_dash.v2);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].acc_dash.z)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].acc_dash.z)[7] = __fkg_store_tmp[1];
((float*)&force[i+0].acc_dash.z)[14] = __fkg_store_tmp[2];
((float*)&force[i+0].acc_dash.z)[21] = __fkg_store_tmp[3];
((float*)&force[i+0].acc_dash.z)[28] = __fkg_store_tmp[4];
((float*)&force[i+0].acc_dash.z)[35] = __fkg_store_tmp[5];
((float*)&force[i+0].acc_dash.z)[42] = __fkg_store_tmp[6];
((float*)&force[i+0].acc_dash.z)[49] = __fkg_store_tmp[7];
}
}

{
__m256 __fkg_tmp_accum;
int index_gather_load15[8] = {0,7,14,21,28,35,42,49};
__m256i vindex_gather_load15 = _mm256_load_si256((const __m256i*)index_gather_load15);
__fkg_tmp_accum = _mm256_i32gather_ps(((float*)&force[i+0].pot),vindex_gather_load15,4);
__fkg_tmp_accum = _mm256_add_ps(__fkg_tmp_accum,FORCE_pot);
{
PIKG::F32 __fkg_store_tmp[8];
_mm256_storeu_ps(__fkg_store_tmp,__fkg_tmp_accum);
((float*)&force[i+0].pot)[0] = __fkg_store_tmp[0];
((float*)&force[i+0].pot)[7] = __fkg_store_tmp[1];
((float*)&force[i+0].pot)[14] = __fkg_store_tmp[2];
((float*)&force[i+0].pot)[21] = __fkg_store_tmp[3];
((float*)&force[i+0].pot)[28] = __fkg_store_tmp[4];
((float*)&force[i+0].pot)[35] = __fkg_store_tmp[5];
((float*)&force[i+0].pot)[42] = __fkg_store_tmp[6];
((float*)&force[i+0].pot)[49] = __fkg_store_tmp[7];
}
}

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
PIKG::S32 EPI_id;

EPI_id = epi[i+0].id;
PIKG::F32 EPI_mass;

EPI_mass = epi[i+0].mass;
PIKG::F32vec EPI_pos;

EPI_pos.x = epi[i+0].pos.x;
EPI_pos.y = epi[i+0].pos.y;
EPI_pos.z = epi[i+0].pos.z;
PIKG::F32 EPI_r_coll;

EPI_r_coll = epi[i+0].r_coll;
PIKG::F32vec EPI_vel;

EPI_vel.x = epi[i+0].vel.x;
EPI_vel.y = epi[i+0].vel.y;
EPI_vel.z = epi[i+0].vel.z;
PIKG::F32vec FORCE_acc;

FORCE_acc.x = 0.0f;
FORCE_acc.y = 0.0f;
FORCE_acc.z = 0.0f;
PIKG::F32vec FORCE_acc_dash;

FORCE_acc_dash.x = 0.0f;
FORCE_acc_dash.y = 0.0f;
FORCE_acc_dash.z = 0.0f;
PIKG::F32 FORCE_pot;

FORCE_pot = 0.0f;
for(j = 0;j < nj;++j){
PIKG::S32 EPJ_id;

EPJ_id = epj[j+0].id;
PIKG::F32 EPJ_mass;

EPJ_mass = epj[j+0].mass;
PIKG::F32vec EPJ_pos;

EPJ_pos.x = epj[j+0].pos.x;
EPJ_pos.y = epj[j+0].pos.y;
EPJ_pos.z = epj[j+0].pos.z;
PIKG::F32 EPJ_r_coll;

EPJ_r_coll = epj[j+0].r_coll;
PIKG::F32vec EPJ_vel;

EPJ_vel.x = epj[j+0].vel.x;
EPJ_vel.y = epj[j+0].vel.y;
EPJ_vel.z = epj[j+0].vel.z;
PIKG::F32vec acc_sprg_tmp;

PIKG::F32vec acc_dash_tmp;

PIKG::F32 pot_sprg_tmp;

PIKG::F32 pot_grav_tmp;

PIKG::F32vec rij;

PIKG::F32 r_coll_tmp;

PIKG::F32 r_coll_sq;

PIKG::F32 r_real_sq;

PIKG::F32 over_r_real;

PIKG::F32 over_r_real_sq;

PIKG::F32 r_coll_cu;

PIKG::F32 over_r_coll_cu;

PIKG::F32 __fkg_tmp0;

PIKG::F32vec __fkg_tmp4;

PIKG::F32 pot_offset;

PIKG::F32 __fkg_tmp5;

PIKG::F32 m_red;

PIKG::F32 r_real;

PIKG::F32 dr;

PIKG::F32 __fkg_tmp1;

PIKG::F32vec __fkg_tmp6;

PIKG::F32 __fkg_tmp7;

PIKG::F32vec vij;

PIKG::F32 rv;

PIKG::F32 __fkg_tmp2;

PIKG::F32vec __fkg_tmp8;

PIKG::F32vec acc_grav_tmp;

PIKG::F32 m_over_r_real;

PIKG::F32 __fkg_tmp3;

PIKG::F32 __fkg_tmp10;

PIKG::F32 __fkg_tmp9;

PIKG::F32 __fkg_tmp12;

PIKG::F32 __fkg_tmp11;

PIKG::F32 __fkg_tmp14;

PIKG::F32 __fkg_tmp13;

PIKG::F32 __fkg_tmp15;

acc_sprg_tmp.x = 0.0f;
acc_sprg_tmp.y = 0.0f;
acc_sprg_tmp.z = 0.0f;
acc_dash_tmp.x = 0.0f;
acc_dash_tmp.y = 0.0f;
acc_dash_tmp.z = 0.0f;
pot_sprg_tmp = 0.0f;
pot_grav_tmp = 0.0f;
if((EPI_id!=EPJ_id)){
rij.x = (EPI_pos.x-EPJ_pos.x);
rij.y = (EPI_pos.y-EPJ_pos.y);
rij.z = (EPI_pos.z-EPJ_pos.z);
r_coll_tmp = (EPI_r_coll+EPJ_r_coll);
r_coll_sq = (r_coll_tmp*r_coll_tmp);
r_real_sq = (rij.z*rij.z+(rij.y*rij.y+(rij.x*rij.x+eps2)));
over_r_real = rsqrt(r_real_sq);
over_r_real_sq = (over_r_real*over_r_real);
if((r_coll_sq>r_real_sq)){
r_coll_cu = (r_coll_sq*r_coll_tmp);
over_r_coll_cu = (1.0f/r_coll_cu);
__fkg_tmp0 = (EPJ_mass*over_r_coll_cu);
__fkg_tmp4.x = ((EPJ_mass*over_r_coll_cu)*rij.x);
__fkg_tmp4.y = ((EPJ_mass*over_r_coll_cu)*rij.y);
__fkg_tmp4.z = ((EPJ_mass*over_r_coll_cu)*rij.z);
pot_offset = (-(1.5f)/r_coll_tmp);
__fkg_tmp5 = ((0.25f*EPJ_mass)*(r_real_sq*over_r_coll_cu+pot_offset));
m_red = (EPJ_mass/(EPI_mass+EPJ_mass));
r_real = (r_real_sq*over_r_real);
dr = (r_coll_tmp-r_real);
__fkg_tmp1 = (((kappa*m_red)*dr)*over_r_real);
__fkg_tmp6.x = ((((kappa*m_red)*dr)*over_r_real)*rij.x);
__fkg_tmp6.y = ((((kappa*m_red)*dr)*over_r_real)*rij.y);
__fkg_tmp6.z = ((((kappa*m_red)*dr)*over_r_real)*rij.z);
__fkg_tmp7 = ((((0.25f*kappa)*m_red)*dr)*dr);
vij.x = (EPI_vel.x-EPJ_vel.x);
vij.y = (EPI_vel.y-EPJ_vel.y);
vij.z = (EPI_vel.z-EPJ_vel.z);
rv = (rij.z*vij.z+(rij.x*vij.x+(rij.y*vij.y)));
__fkg_tmp2 = (((eta*m_red)*rv)*over_r_real_sq);
__fkg_tmp8.x = ((((eta*m_red)*rv)*over_r_real_sq)*rij.x);
__fkg_tmp8.y = ((((eta*m_red)*rv)*over_r_real_sq)*rij.y);
__fkg_tmp8.z = ((((eta*m_red)*rv)*over_r_real_sq)*rij.z);
acc_grav_tmp.x = __fkg_tmp4.x;
acc_grav_tmp.y = __fkg_tmp4.y;
acc_grav_tmp.z = __fkg_tmp4.z;
pot_grav_tmp = __fkg_tmp5;
acc_sprg_tmp.x = __fkg_tmp6.x;
acc_sprg_tmp.y = __fkg_tmp6.y;
acc_sprg_tmp.z = __fkg_tmp6.z;
pot_sprg_tmp = __fkg_tmp7;
acc_dash_tmp.x = __fkg_tmp8.x;
acc_dash_tmp.y = __fkg_tmp8.y;
acc_dash_tmp.z = __fkg_tmp8.z;
}else{
m_over_r_real = (EPJ_mass*over_r_real);
__fkg_tmp3 = (m_over_r_real*over_r_real_sq);
__fkg_tmp4.x = ((m_over_r_real*over_r_real_sq)*rij.x);
__fkg_tmp4.y = ((m_over_r_real*over_r_real_sq)*rij.y);
__fkg_tmp4.z = ((m_over_r_real*over_r_real_sq)*rij.z);
__fkg_tmp5 = (-(0.5f)*m_over_r_real);
acc_grav_tmp.x = __fkg_tmp4.x;
acc_grav_tmp.y = __fkg_tmp4.y;
acc_grav_tmp.z = __fkg_tmp4.z;
pot_grav_tmp = __fkg_tmp5;
}
}
__fkg_tmp10 = (acc_grav_tmp.x+acc_sprg_tmp.x);
__fkg_tmp9 = (__fkg_tmp10+acc_dash_tmp.x);
FORCE_acc.x = (FORCE_acc.x+__fkg_tmp9);
__fkg_tmp12 = (acc_grav_tmp.y+acc_sprg_tmp.y);
__fkg_tmp11 = (__fkg_tmp12+acc_dash_tmp.y);
FORCE_acc.y = (FORCE_acc.y+__fkg_tmp11);
__fkg_tmp14 = (acc_grav_tmp.z+acc_sprg_tmp.z);
__fkg_tmp13 = (__fkg_tmp14+acc_dash_tmp.z);
FORCE_acc.z = (FORCE_acc.z+__fkg_tmp13);
FORCE_acc_dash.x = (FORCE_acc_dash.x+acc_dash_tmp.x);
FORCE_acc_dash.y = (FORCE_acc_dash.y+acc_dash_tmp.y);
FORCE_acc_dash.z = (FORCE_acc_dash.z+acc_dash_tmp.z);
__fkg_tmp15 = (pot_grav_tmp+pot_sprg_tmp);
FORCE_pot = (FORCE_pot+__fkg_tmp15);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+FORCE_acc.x);
force[i+0].acc.y = (force[i+0].acc.y+FORCE_acc.y);
force[i+0].acc.z = (force[i+0].acc.z+FORCE_acc.z);
force[i+0].acc_dash.x = (force[i+0].acc_dash.x+FORCE_acc_dash.x);
force[i+0].acc_dash.y = (force[i+0].acc_dash.y+FORCE_acc_dash.y);
force[i+0].acc_dash.z = (force[i+0].acc_dash.z+FORCE_acc_dash.z);
force[i+0].pot = (force[i+0].pot+FORCE_pot);
} // loop of i
} // end loop of reference 
} // Kernel_I8_J1 definition 
void Kernel_I1_J8(const Epi0* __restrict__ epi,const PIKG::S32 ni,const Epj0* __restrict__ epj,const PIKG::S32 nj,Force0* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
for(i = 0;i < (ni/1)*1;++i){
__m256i EPI_id;

EPI_id = _mm256_set1_epi32(epi[i+0].id);

__m256 EPI_mass;

EPI_mass = _mm256_set1_ps(epi[i+0].mass);

__m256x3 EPI_pos;

EPI_pos.v0 = _mm256_set1_ps(epi[i+0].pos.x);

EPI_pos.v1 = _mm256_set1_ps(epi[i+0].pos.y);

EPI_pos.v2 = _mm256_set1_ps(epi[i+0].pos.z);

__m256 EPI_r_coll;

EPI_r_coll = _mm256_set1_ps(epi[i+0].r_coll);

__m256x3 EPI_vel;

EPI_vel.v0 = _mm256_set1_ps(epi[i+0].vel.x);

EPI_vel.v1 = _mm256_set1_ps(epi[i+0].vel.y);

EPI_vel.v2 = _mm256_set1_ps(epi[i+0].vel.z);

__m256x3 FORCE_acc;

FORCE_acc.v0 = _mm256_set1_ps(0.0f);
FORCE_acc.v1 = _mm256_set1_ps(0.0f);
FORCE_acc.v2 = _mm256_set1_ps(0.0f);
__m256x3 FORCE_acc_dash;

FORCE_acc_dash.v0 = _mm256_set1_ps(0.0f);
FORCE_acc_dash.v1 = _mm256_set1_ps(0.0f);
FORCE_acc_dash.v2 = _mm256_set1_ps(0.0f);
__m256 FORCE_pot;

FORCE_pot = _mm256_set1_ps(0.0f);
for(j = 0;j < (nj/8)*8;j += 8){
__m256i EPJ_id;

int index_gather_load16[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load16 = _mm256_load_si256((const __m256i*)index_gather_load16);
EPJ_id = _mm256_i32gather_epi32(((int*)&epj[j+0].id),vindex_gather_load16,4);
__m256 EPJ_mass;

int index_gather_load17[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load17 = _mm256_load_si256((const __m256i*)index_gather_load17);
EPJ_mass = _mm256_i32gather_ps(((float*)&epj[j+0].mass),vindex_gather_load17,4);
__m256x3 EPJ_pos;

int index_gather_load18[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load18 = _mm256_load_si256((const __m256i*)index_gather_load18);
EPJ_pos.v0 = _mm256_i32gather_ps(((float*)&epj[j+0].pos.x),vindex_gather_load18,4);
int index_gather_load19[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load19 = _mm256_load_si256((const __m256i*)index_gather_load19);
EPJ_pos.v1 = _mm256_i32gather_ps(((float*)&epj[j+0].pos.y),vindex_gather_load19,4);
int index_gather_load20[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load20 = _mm256_load_si256((const __m256i*)index_gather_load20);
EPJ_pos.v2 = _mm256_i32gather_ps(((float*)&epj[j+0].pos.z),vindex_gather_load20,4);
__m256 EPJ_r_coll;

int index_gather_load21[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load21 = _mm256_load_si256((const __m256i*)index_gather_load21);
EPJ_r_coll = _mm256_i32gather_ps(((float*)&epj[j+0].r_coll),vindex_gather_load21,4);
__m256x3 EPJ_vel;

int index_gather_load22[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load22 = _mm256_load_si256((const __m256i*)index_gather_load22);
EPJ_vel.v0 = _mm256_i32gather_ps(((float*)&epj[j+0].vel.x),vindex_gather_load22,4);
int index_gather_load23[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load23 = _mm256_load_si256((const __m256i*)index_gather_load23);
EPJ_vel.v1 = _mm256_i32gather_ps(((float*)&epj[j+0].vel.y),vindex_gather_load23,4);
int index_gather_load24[8] = {0,9,18,27,36,45,54,63};
__m256i vindex_gather_load24 = _mm256_load_si256((const __m256i*)index_gather_load24);
EPJ_vel.v2 = _mm256_i32gather_ps(((float*)&epj[j+0].vel.z),vindex_gather_load24,4);
__m256x3 acc_sprg_tmp;

__m256x3 acc_dash_tmp;

__m256 pot_sprg_tmp;

__m256 pot_grav_tmp;

__m256x3 rij;

__m256 r_coll_tmp;

__m256 r_coll_sq;

__m256 r_real_sq;

__m256 over_r_real;

__m256 over_r_real_sq;

__m256 r_coll_cu;

__m256 over_r_coll_cu;

__m256 __fkg_tmp0;

__m256x3 __fkg_tmp4;

__m256 pot_offset;

__m256 __fkg_tmp5;

__m256 m_red;

__m256 r_real;

__m256 dr;

__m256 __fkg_tmp1;

__m256x3 __fkg_tmp6;

__m256 __fkg_tmp7;

__m256x3 vij;

__m256 rv;

__m256 __fkg_tmp2;

__m256x3 __fkg_tmp8;

__m256x3 acc_grav_tmp;

__m256 m_over_r_real;

__m256 __fkg_tmp3;

__m256 __fkg_tmp10;

__m256 __fkg_tmp9;

__m256 __fkg_tmp12;

__m256 __fkg_tmp11;

__m256 __fkg_tmp14;

__m256 __fkg_tmp13;

__m256 __fkg_tmp15;

acc_sprg_tmp.v0 = _mm256_set1_ps(0.0f);
acc_sprg_tmp.v1 = _mm256_set1_ps(0.0f);
acc_sprg_tmp.v2 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v0 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v1 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v2 = _mm256_set1_ps(0.0f);
pot_sprg_tmp = _mm256_set1_ps(0.0f);
pot_grav_tmp = _mm256_set1_ps(0.0f);
{
__m256 pg1;
__m256 pg0;
pg1 = _mm256_cmp_ps(_mm256_castsi256_ps(EPI_id),_mm256_castsi256_ps(EPJ_id),_CMP_NEQ_OQ);
pg0 = pg1;

rij.v0 = _mm256_sub_ps(EPI_pos.v0,EPJ_pos.v0);
rij.v1 = _mm256_sub_ps(EPI_pos.v1,EPJ_pos.v1);
rij.v2 = _mm256_sub_ps(EPI_pos.v2,EPJ_pos.v2);
r_coll_tmp = _mm256_add_ps(EPI_r_coll,EPJ_r_coll);
r_coll_sq = _mm256_mul_ps(r_coll_tmp,r_coll_tmp);
r_real_sq = _mm256_fmadd_ps(rij.v2,rij.v2,_mm256_fmadd_ps(rij.v1,rij.v1,_mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2))));
over_r_real = rsqrt(r_real_sq);
over_r_real_sq = _mm256_mul_ps(over_r_real,over_r_real);
{
__m256 pg3;
__m256 pg2;
pg3 = _mm256_cmp_ps(r_coll_sq,r_real_sq,_CMP_GT_OS);
pg2 = pg3;
pg3 = _mm256_and_ps(pg3,pg1);

r_coll_cu = _mm256_mul_ps(r_coll_sq,r_coll_tmp);
over_r_coll_cu = _mm256_div_ps(_mm256_set1_ps(1.0f),r_coll_cu);
__fkg_tmp0 = _mm256_mul_ps(EPJ_mass,over_r_coll_cu);
__fkg_tmp4.v0 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v0);
__fkg_tmp4.v1 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v1);
__fkg_tmp4.v2 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v2);
pot_offset = _mm256_div_ps(_mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),_mm256_set1_ps(1.5f)),r_coll_tmp);
__fkg_tmp5 = _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.25f),EPJ_mass),_mm256_fmadd_ps(r_real_sq,over_r_coll_cu,pot_offset));
m_red = _mm256_div_ps(EPJ_mass,_mm256_add_ps(EPI_mass,EPJ_mass));
r_real = _mm256_mul_ps(r_real_sq,over_r_real);
dr = _mm256_sub_ps(r_coll_tmp,r_real);
__fkg_tmp1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real);
__fkg_tmp6.v0 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v0);
__fkg_tmp6.v1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v1);
__fkg_tmp6.v2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v2);
__fkg_tmp7 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.25f),_mm256_set1_ps(kappa)),m_red),dr),dr);
vij.v0 = _mm256_sub_ps(EPI_vel.v0,EPJ_vel.v0);
vij.v1 = _mm256_sub_ps(EPI_vel.v1,EPJ_vel.v1);
vij.v2 = _mm256_sub_ps(EPI_vel.v2,EPJ_vel.v2);
rv = _mm256_fmadd_ps(rij.v2,vij.v2,_mm256_fmadd_ps(rij.v0,vij.v0,_mm256_mul_ps(rij.v1,vij.v1)));
__fkg_tmp2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq);
__fkg_tmp8.v0 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v0);
__fkg_tmp8.v1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v1);
__fkg_tmp8.v2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v2);
acc_grav_tmp.v0 = _mm256_blendv_ps(acc_grav_tmp.v0,__fkg_tmp4.v0,pg3);;
acc_grav_tmp.v1 = _mm256_blendv_ps(acc_grav_tmp.v1,__fkg_tmp4.v1,pg3);;
acc_grav_tmp.v2 = _mm256_blendv_ps(acc_grav_tmp.v2,__fkg_tmp4.v2,pg3);;
pot_grav_tmp = _mm256_blendv_ps(pot_grav_tmp,__fkg_tmp5,pg3);;
acc_sprg_tmp.v0 = _mm256_blendv_ps(acc_sprg_tmp.v0,__fkg_tmp6.v0,pg3);;
acc_sprg_tmp.v1 = _mm256_blendv_ps(acc_sprg_tmp.v1,__fkg_tmp6.v1,pg3);;
acc_sprg_tmp.v2 = _mm256_blendv_ps(acc_sprg_tmp.v2,__fkg_tmp6.v2,pg3);;
pot_sprg_tmp = _mm256_blendv_ps(pot_sprg_tmp,__fkg_tmp7,pg3);;
acc_dash_tmp.v0 = _mm256_blendv_ps(acc_dash_tmp.v0,__fkg_tmp8.v0,pg3);;
acc_dash_tmp.v1 = _mm256_blendv_ps(acc_dash_tmp.v1,__fkg_tmp8.v1,pg3);;
acc_dash_tmp.v2 = _mm256_blendv_ps(acc_dash_tmp.v2,__fkg_tmp8.v2,pg3);;
pg3 = _mm256_xor_ps(pg2,_mm256_castsi256_ps(_mm256_set1_epi32(-1)));
pg3 = _mm256_and_ps(pg3,pg1);

m_over_r_real = _mm256_mul_ps(EPJ_mass,over_r_real);
__fkg_tmp3 = _mm256_mul_ps(m_over_r_real,over_r_real_sq);
__fkg_tmp4.v0 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v0);
__fkg_tmp4.v1 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v1);
__fkg_tmp4.v2 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v2);
__fkg_tmp5 = _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),_mm256_set1_ps(0.5f)),m_over_r_real);
acc_grav_tmp.v0 = _mm256_blendv_ps(acc_grav_tmp.v0,__fkg_tmp4.v0,pg3);;
acc_grav_tmp.v1 = _mm256_blendv_ps(acc_grav_tmp.v1,__fkg_tmp4.v1,pg3);;
acc_grav_tmp.v2 = _mm256_blendv_ps(acc_grav_tmp.v2,__fkg_tmp4.v2,pg3);;
pot_grav_tmp = _mm256_blendv_ps(pot_grav_tmp,__fkg_tmp5,pg3);;
}

}

__fkg_tmp10 = _mm256_add_ps(acc_grav_tmp.v0,acc_sprg_tmp.v0);
__fkg_tmp9 = _mm256_add_ps(__fkg_tmp10,acc_dash_tmp.v0);
FORCE_acc.v0 = _mm256_add_ps(FORCE_acc.v0,__fkg_tmp9);
__fkg_tmp12 = _mm256_add_ps(acc_grav_tmp.v1,acc_sprg_tmp.v1);
__fkg_tmp11 = _mm256_add_ps(__fkg_tmp12,acc_dash_tmp.v1);
FORCE_acc.v1 = _mm256_add_ps(FORCE_acc.v1,__fkg_tmp11);
__fkg_tmp14 = _mm256_add_ps(acc_grav_tmp.v2,acc_sprg_tmp.v2);
__fkg_tmp13 = _mm256_add_ps(__fkg_tmp14,acc_dash_tmp.v2);
FORCE_acc.v2 = _mm256_add_ps(FORCE_acc.v2,__fkg_tmp13);
FORCE_acc_dash.v0 = _mm256_add_ps(FORCE_acc_dash.v0,acc_dash_tmp.v0);
FORCE_acc_dash.v1 = _mm256_add_ps(FORCE_acc_dash.v1,acc_dash_tmp.v1);
FORCE_acc_dash.v2 = _mm256_add_ps(FORCE_acc_dash.v2,acc_dash_tmp.v2);
__fkg_tmp15 = _mm256_add_ps(pot_grav_tmp,pot_sprg_tmp);
FORCE_pot = _mm256_add_ps(FORCE_pot,__fkg_tmp15);
} // loop of j

if(j<nj){ // tail j loop
__m256x3 __fkg_tmp16;

__fkg_tmp16.v0 = FORCE_acc.v0;
__fkg_tmp16.v1 = FORCE_acc.v1;
__fkg_tmp16.v2 = FORCE_acc.v2;
__m256x3 __fkg_tmp17;

__fkg_tmp17.v0 = FORCE_acc_dash.v0;
__fkg_tmp17.v1 = FORCE_acc_dash.v1;
__fkg_tmp17.v2 = FORCE_acc_dash.v2;
__m256 __fkg_tmp18;

__fkg_tmp18 = FORCE_pot;
for(;j < nj;++j){
__m256i EPJ_id;

EPJ_id = _mm256_set1_epi32(epj[j+0].id);

__m256 EPJ_mass;

EPJ_mass = _mm256_set1_ps(epj[j+0].mass);

__m256x3 EPJ_pos;

EPJ_pos.v0 = _mm256_set1_ps(epj[j+0].pos.x);

EPJ_pos.v1 = _mm256_set1_ps(epj[j+0].pos.y);

EPJ_pos.v2 = _mm256_set1_ps(epj[j+0].pos.z);

__m256 EPJ_r_coll;

EPJ_r_coll = _mm256_set1_ps(epj[j+0].r_coll);

__m256x3 EPJ_vel;

EPJ_vel.v0 = _mm256_set1_ps(epj[j+0].vel.x);

EPJ_vel.v1 = _mm256_set1_ps(epj[j+0].vel.y);

EPJ_vel.v2 = _mm256_set1_ps(epj[j+0].vel.z);

__m256x3 acc_sprg_tmp;

__m256x3 acc_dash_tmp;

__m256 pot_sprg_tmp;

__m256 pot_grav_tmp;

__m256x3 rij;

__m256 r_coll_tmp;

__m256 r_coll_sq;

__m256 r_real_sq;

__m256 over_r_real;

__m256 over_r_real_sq;

__m256 r_coll_cu;

__m256 over_r_coll_cu;

__m256 __fkg_tmp0;

__m256x3 __fkg_tmp4;

__m256 pot_offset;

__m256 __fkg_tmp5;

__m256 m_red;

__m256 r_real;

__m256 dr;

__m256 __fkg_tmp1;

__m256x3 __fkg_tmp6;

__m256 __fkg_tmp7;

__m256x3 vij;

__m256 rv;

__m256 __fkg_tmp2;

__m256x3 __fkg_tmp8;

__m256x3 acc_grav_tmp;

__m256 m_over_r_real;

__m256 __fkg_tmp3;

__m256 __fkg_tmp10;

__m256 __fkg_tmp9;

__m256 __fkg_tmp12;

__m256 __fkg_tmp11;

__m256 __fkg_tmp14;

__m256 __fkg_tmp13;

__m256 __fkg_tmp15;

acc_sprg_tmp.v0 = _mm256_set1_ps(0.0f);
acc_sprg_tmp.v1 = _mm256_set1_ps(0.0f);
acc_sprg_tmp.v2 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v0 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v1 = _mm256_set1_ps(0.0f);
acc_dash_tmp.v2 = _mm256_set1_ps(0.0f);
pot_sprg_tmp = _mm256_set1_ps(0.0f);
pot_grav_tmp = _mm256_set1_ps(0.0f);
{
__m256 pg1;
__m256 pg0;
pg1 = _mm256_cmp_ps(_mm256_castsi256_ps(EPI_id),_mm256_castsi256_ps(EPJ_id),_CMP_NEQ_OQ);
pg0 = pg1;

rij.v0 = _mm256_sub_ps(EPI_pos.v0,EPJ_pos.v0);
rij.v1 = _mm256_sub_ps(EPI_pos.v1,EPJ_pos.v1);
rij.v2 = _mm256_sub_ps(EPI_pos.v2,EPJ_pos.v2);
r_coll_tmp = _mm256_add_ps(EPI_r_coll,EPJ_r_coll);
r_coll_sq = _mm256_mul_ps(r_coll_tmp,r_coll_tmp);
r_real_sq = _mm256_fmadd_ps(rij.v2,rij.v2,_mm256_fmadd_ps(rij.v1,rij.v1,_mm256_fmadd_ps(rij.v0,rij.v0,_mm256_set1_ps(eps2))));
over_r_real = rsqrt(r_real_sq);
over_r_real_sq = _mm256_mul_ps(over_r_real,over_r_real);
{
__m256 pg3;
__m256 pg2;
pg3 = _mm256_cmp_ps(r_coll_sq,r_real_sq,_CMP_GT_OS);
pg2 = pg3;
pg3 = _mm256_and_ps(pg3,pg1);

r_coll_cu = _mm256_mul_ps(r_coll_sq,r_coll_tmp);
over_r_coll_cu = _mm256_div_ps(_mm256_set1_ps(1.0f),r_coll_cu);
__fkg_tmp0 = _mm256_mul_ps(EPJ_mass,over_r_coll_cu);
__fkg_tmp4.v0 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v0);
__fkg_tmp4.v1 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v1);
__fkg_tmp4.v2 = _mm256_mul_ps(_mm256_mul_ps(EPJ_mass,over_r_coll_cu),rij.v2);
pot_offset = _mm256_div_ps(_mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),_mm256_set1_ps(1.5f)),r_coll_tmp);
__fkg_tmp5 = _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.25f),EPJ_mass),_mm256_fmadd_ps(r_real_sq,over_r_coll_cu,pot_offset));
m_red = _mm256_div_ps(EPJ_mass,_mm256_add_ps(EPI_mass,EPJ_mass));
r_real = _mm256_mul_ps(r_real_sq,over_r_real);
dr = _mm256_sub_ps(r_coll_tmp,r_real);
__fkg_tmp1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real);
__fkg_tmp6.v0 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v0);
__fkg_tmp6.v1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v1);
__fkg_tmp6.v2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(kappa),m_red),dr),over_r_real),rij.v2);
__fkg_tmp7 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.25f),_mm256_set1_ps(kappa)),m_red),dr),dr);
vij.v0 = _mm256_sub_ps(EPI_vel.v0,EPJ_vel.v0);
vij.v1 = _mm256_sub_ps(EPI_vel.v1,EPJ_vel.v1);
vij.v2 = _mm256_sub_ps(EPI_vel.v2,EPJ_vel.v2);
rv = _mm256_fmadd_ps(rij.v2,vij.v2,_mm256_fmadd_ps(rij.v0,vij.v0,_mm256_mul_ps(rij.v1,vij.v1)));
__fkg_tmp2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq);
__fkg_tmp8.v0 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v0);
__fkg_tmp8.v1 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v1);
__fkg_tmp8.v2 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(eta),m_red),rv),over_r_real_sq),rij.v2);
acc_grav_tmp.v0 = _mm256_blendv_ps(acc_grav_tmp.v0,__fkg_tmp4.v0,pg3);;
acc_grav_tmp.v1 = _mm256_blendv_ps(acc_grav_tmp.v1,__fkg_tmp4.v1,pg3);;
acc_grav_tmp.v2 = _mm256_blendv_ps(acc_grav_tmp.v2,__fkg_tmp4.v2,pg3);;
pot_grav_tmp = _mm256_blendv_ps(pot_grav_tmp,__fkg_tmp5,pg3);;
acc_sprg_tmp.v0 = _mm256_blendv_ps(acc_sprg_tmp.v0,__fkg_tmp6.v0,pg3);;
acc_sprg_tmp.v1 = _mm256_blendv_ps(acc_sprg_tmp.v1,__fkg_tmp6.v1,pg3);;
acc_sprg_tmp.v2 = _mm256_blendv_ps(acc_sprg_tmp.v2,__fkg_tmp6.v2,pg3);;
pot_sprg_tmp = _mm256_blendv_ps(pot_sprg_tmp,__fkg_tmp7,pg3);;
acc_dash_tmp.v0 = _mm256_blendv_ps(acc_dash_tmp.v0,__fkg_tmp8.v0,pg3);;
acc_dash_tmp.v1 = _mm256_blendv_ps(acc_dash_tmp.v1,__fkg_tmp8.v1,pg3);;
acc_dash_tmp.v2 = _mm256_blendv_ps(acc_dash_tmp.v2,__fkg_tmp8.v2,pg3);;
pg3 = _mm256_xor_ps(pg2,_mm256_castsi256_ps(_mm256_set1_epi32(-1)));
pg3 = _mm256_and_ps(pg3,pg1);

m_over_r_real = _mm256_mul_ps(EPJ_mass,over_r_real);
__fkg_tmp3 = _mm256_mul_ps(m_over_r_real,over_r_real_sq);
__fkg_tmp4.v0 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v0);
__fkg_tmp4.v1 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v1);
__fkg_tmp4.v2 = _mm256_mul_ps(_mm256_mul_ps(m_over_r_real,over_r_real_sq),rij.v2);
__fkg_tmp5 = _mm256_mul_ps(_mm256_sub_ps(_mm256_set1_ps((PIKG::F32)0.0),_mm256_set1_ps(0.5f)),m_over_r_real);
acc_grav_tmp.v0 = _mm256_blendv_ps(acc_grav_tmp.v0,__fkg_tmp4.v0,pg3);;
acc_grav_tmp.v1 = _mm256_blendv_ps(acc_grav_tmp.v1,__fkg_tmp4.v1,pg3);;
acc_grav_tmp.v2 = _mm256_blendv_ps(acc_grav_tmp.v2,__fkg_tmp4.v2,pg3);;
pot_grav_tmp = _mm256_blendv_ps(pot_grav_tmp,__fkg_tmp5,pg3);;
}

}

__fkg_tmp10 = _mm256_add_ps(acc_grav_tmp.v0,acc_sprg_tmp.v0);
__fkg_tmp9 = _mm256_add_ps(__fkg_tmp10,acc_dash_tmp.v0);
FORCE_acc.v0 = _mm256_add_ps(FORCE_acc.v0,__fkg_tmp9);
__fkg_tmp12 = _mm256_add_ps(acc_grav_tmp.v1,acc_sprg_tmp.v1);
__fkg_tmp11 = _mm256_add_ps(__fkg_tmp12,acc_dash_tmp.v1);
FORCE_acc.v1 = _mm256_add_ps(FORCE_acc.v1,__fkg_tmp11);
__fkg_tmp14 = _mm256_add_ps(acc_grav_tmp.v2,acc_sprg_tmp.v2);
__fkg_tmp13 = _mm256_add_ps(__fkg_tmp14,acc_dash_tmp.v2);
FORCE_acc.v2 = _mm256_add_ps(FORCE_acc.v2,__fkg_tmp13);
FORCE_acc_dash.v0 = _mm256_add_ps(FORCE_acc_dash.v0,acc_dash_tmp.v0);
FORCE_acc_dash.v1 = _mm256_add_ps(FORCE_acc_dash.v1,acc_dash_tmp.v1);
FORCE_acc_dash.v2 = _mm256_add_ps(FORCE_acc_dash.v2,acc_dash_tmp.v2);
__fkg_tmp15 = _mm256_add_ps(pot_grav_tmp,pot_sprg_tmp);
FORCE_pot = _mm256_add_ps(FORCE_pot,__fkg_tmp15);
} // loop of j
FORCE_acc.v0 = _mm256_blend_ps(__fkg_tmp16.v0,FORCE_acc.v0,0b00000001);
FORCE_acc.v1 = _mm256_blend_ps(__fkg_tmp16.v1,FORCE_acc.v1,0b00000001);
FORCE_acc.v2 = _mm256_blend_ps(__fkg_tmp16.v2,FORCE_acc.v2,0b00000001);
FORCE_acc_dash.v0 = _mm256_blend_ps(__fkg_tmp17.v0,FORCE_acc_dash.v0,0b00000001);
FORCE_acc_dash.v1 = _mm256_blend_ps(__fkg_tmp17.v1,FORCE_acc_dash.v1,0b00000001);
FORCE_acc_dash.v2 = _mm256_blend_ps(__fkg_tmp17.v2,FORCE_acc_dash.v2,0b00000001);
FORCE_pot = _mm256_blend_ps(__fkg_tmp18,FORCE_pot,0b00000001);
} // if of j tail loop

{
FORCE_acc.v0 = _mm256_add_ps(FORCE_acc.v0,_mm256_shuffle_ps(FORCE_acc.v0,FORCE_acc.v0,0xb1));
FORCE_acc.v0 = _mm256_add_ps(FORCE_acc.v0,_mm256_shuffle_ps(FORCE_acc.v0,FORCE_acc.v0,0xee));
FORCE_acc.v0 = _mm256_add_ps(FORCE_acc.v0,_mm256_castps128_ps256(_mm256_extractf128_ps(FORCE_acc.v0,1)));
((float*)&force[i+0].acc.x)[0] = (((float*)&force[i+0].acc.x)[0]+FORCE_acc.v0[0]);
}

{
FORCE_acc.v1 = _mm256_add_ps(FORCE_acc.v1,_mm256_shuffle_ps(FORCE_acc.v1,FORCE_acc.v1,0xb1));
FORCE_acc.v1 = _mm256_add_ps(FORCE_acc.v1,_mm256_shuffle_ps(FORCE_acc.v1,FORCE_acc.v1,0xee));
FORCE_acc.v1 = _mm256_add_ps(FORCE_acc.v1,_mm256_castps128_ps256(_mm256_extractf128_ps(FORCE_acc.v1,1)));
((float*)&force[i+0].acc.y)[0] = (((float*)&force[i+0].acc.y)[0]+FORCE_acc.v1[0]);
}

{
FORCE_acc.v2 = _mm256_add_ps(FORCE_acc.v2,_mm256_shuffle_ps(FORCE_acc.v2,FORCE_acc.v2,0xb1));
FORCE_acc.v2 = _mm256_add_ps(FORCE_acc.v2,_mm256_shuffle_ps(FORCE_acc.v2,FORCE_acc.v2,0xee));
FORCE_acc.v2 = _mm256_add_ps(FORCE_acc.v2,_mm256_castps128_ps256(_mm256_extractf128_ps(FORCE_acc.v2,1)));
((float*)&force[i+0].acc.z)[0] = (((float*)&force[i+0].acc.z)[0]+FORCE_acc.v2[0]);
}

{
FORCE_acc_dash.v0 = _mm256_add_ps(FORCE_acc_dash.v0,_mm256_shuffle_ps(FORCE_acc_dash.v0,FORCE_acc_dash.v0,0xb1));
FORCE_acc_dash.v0 = _mm256_add_ps(FORCE_acc_dash.v0,_mm256_shuffle_ps(FORCE_acc_dash.v0,FORCE_acc_dash.v0,0xee));
FORCE_acc_dash.v0 = _mm256_add_ps(FORCE_acc_dash.v0,_mm256_castps128_ps256(_mm256_extractf128_ps(FORCE_acc_dash.v0,1)));
((float*)&force[i+0].acc_dash.x)[0] = (((float*)&force[i+0].acc_dash.x)[0]+FORCE_acc_dash.v0[0]);
}

{
FORCE_acc_dash.v1 = _mm256_add_ps(FORCE_acc_dash.v1,_mm256_shuffle_ps(FORCE_acc_dash.v1,FORCE_acc_dash.v1,0xb1));
FORCE_acc_dash.v1 = _mm256_add_ps(FORCE_acc_dash.v1,_mm256_shuffle_ps(FORCE_acc_dash.v1,FORCE_acc_dash.v1,0xee));
FORCE_acc_dash.v1 = _mm256_add_ps(FORCE_acc_dash.v1,_mm256_castps128_ps256(_mm256_extractf128_ps(FORCE_acc_dash.v1,1)));
((float*)&force[i+0].acc_dash.y)[0] = (((float*)&force[i+0].acc_dash.y)[0]+FORCE_acc_dash.v1[0]);
}

{
FORCE_acc_dash.v2 = _mm256_add_ps(FORCE_acc_dash.v2,_mm256_shuffle_ps(FORCE_acc_dash.v2,FORCE_acc_dash.v2,0xb1));
FORCE_acc_dash.v2 = _mm256_add_ps(FORCE_acc_dash.v2,_mm256_shuffle_ps(FORCE_acc_dash.v2,FORCE_acc_dash.v2,0xee));
FORCE_acc_dash.v2 = _mm256_add_ps(FORCE_acc_dash.v2,_mm256_castps128_ps256(_mm256_extractf128_ps(FORCE_acc_dash.v2,1)));
((float*)&force[i+0].acc_dash.z)[0] = (((float*)&force[i+0].acc_dash.z)[0]+FORCE_acc_dash.v2[0]);
}

{
FORCE_pot = _mm256_add_ps(FORCE_pot,_mm256_shuffle_ps(FORCE_pot,FORCE_pot,0xb1));
FORCE_pot = _mm256_add_ps(FORCE_pot,_mm256_shuffle_ps(FORCE_pot,FORCE_pot,0xee));
FORCE_pot = _mm256_add_ps(FORCE_pot,_mm256_castps128_ps256(_mm256_extractf128_ps(FORCE_pot,1)));
((float*)&force[i+0].pot)[0] = (((float*)&force[i+0].pot)[0]+FORCE_pot[0]);
}

} // loop of i
{ // tail loop of reference 
for(;i < ni;++i){
PIKG::S32 EPI_id;

EPI_id = epi[i+0].id;
PIKG::F32 EPI_mass;

EPI_mass = epi[i+0].mass;
PIKG::F32vec EPI_pos;

EPI_pos.x = epi[i+0].pos.x;
EPI_pos.y = epi[i+0].pos.y;
EPI_pos.z = epi[i+0].pos.z;
PIKG::F32 EPI_r_coll;

EPI_r_coll = epi[i+0].r_coll;
PIKG::F32vec EPI_vel;

EPI_vel.x = epi[i+0].vel.x;
EPI_vel.y = epi[i+0].vel.y;
EPI_vel.z = epi[i+0].vel.z;
PIKG::F32vec FORCE_acc;

FORCE_acc.x = 0.0f;
FORCE_acc.y = 0.0f;
FORCE_acc.z = 0.0f;
PIKG::F32vec FORCE_acc_dash;

FORCE_acc_dash.x = 0.0f;
FORCE_acc_dash.y = 0.0f;
FORCE_acc_dash.z = 0.0f;
PIKG::F32 FORCE_pot;

FORCE_pot = 0.0f;
for(j = 0;j < nj;++j){
PIKG::S32 EPJ_id;

EPJ_id = epj[j+0].id;
PIKG::F32 EPJ_mass;

EPJ_mass = epj[j+0].mass;
PIKG::F32vec EPJ_pos;

EPJ_pos.x = epj[j+0].pos.x;
EPJ_pos.y = epj[j+0].pos.y;
EPJ_pos.z = epj[j+0].pos.z;
PIKG::F32 EPJ_r_coll;

EPJ_r_coll = epj[j+0].r_coll;
PIKG::F32vec EPJ_vel;

EPJ_vel.x = epj[j+0].vel.x;
EPJ_vel.y = epj[j+0].vel.y;
EPJ_vel.z = epj[j+0].vel.z;
PIKG::F32vec acc_sprg_tmp;

PIKG::F32vec acc_dash_tmp;

PIKG::F32 pot_sprg_tmp;

PIKG::F32 pot_grav_tmp;

PIKG::F32vec rij;

PIKG::F32 r_coll_tmp;

PIKG::F32 r_coll_sq;

PIKG::F32 r_real_sq;

PIKG::F32 over_r_real;

PIKG::F32 over_r_real_sq;

PIKG::F32 r_coll_cu;

PIKG::F32 over_r_coll_cu;

PIKG::F32 __fkg_tmp0;

PIKG::F32vec __fkg_tmp4;

PIKG::F32 pot_offset;

PIKG::F32 __fkg_tmp5;

PIKG::F32 m_red;

PIKG::F32 r_real;

PIKG::F32 dr;

PIKG::F32 __fkg_tmp1;

PIKG::F32vec __fkg_tmp6;

PIKG::F32 __fkg_tmp7;

PIKG::F32vec vij;

PIKG::F32 rv;

PIKG::F32 __fkg_tmp2;

PIKG::F32vec __fkg_tmp8;

PIKG::F32vec acc_grav_tmp;

PIKG::F32 m_over_r_real;

PIKG::F32 __fkg_tmp3;

PIKG::F32 __fkg_tmp10;

PIKG::F32 __fkg_tmp9;

PIKG::F32 __fkg_tmp12;

PIKG::F32 __fkg_tmp11;

PIKG::F32 __fkg_tmp14;

PIKG::F32 __fkg_tmp13;

PIKG::F32 __fkg_tmp15;

acc_sprg_tmp.x = 0.0f;
acc_sprg_tmp.y = 0.0f;
acc_sprg_tmp.z = 0.0f;
acc_dash_tmp.x = 0.0f;
acc_dash_tmp.y = 0.0f;
acc_dash_tmp.z = 0.0f;
pot_sprg_tmp = 0.0f;
pot_grav_tmp = 0.0f;
if((EPI_id!=EPJ_id)){
rij.x = (EPI_pos.x-EPJ_pos.x);
rij.y = (EPI_pos.y-EPJ_pos.y);
rij.z = (EPI_pos.z-EPJ_pos.z);
r_coll_tmp = (EPI_r_coll+EPJ_r_coll);
r_coll_sq = (r_coll_tmp*r_coll_tmp);
r_real_sq = (rij.z*rij.z+(rij.y*rij.y+(rij.x*rij.x+eps2)));
over_r_real = rsqrt(r_real_sq);
over_r_real_sq = (over_r_real*over_r_real);
if((r_coll_sq>r_real_sq)){
r_coll_cu = (r_coll_sq*r_coll_tmp);
over_r_coll_cu = (1.0f/r_coll_cu);
__fkg_tmp0 = (EPJ_mass*over_r_coll_cu);
__fkg_tmp4.x = ((EPJ_mass*over_r_coll_cu)*rij.x);
__fkg_tmp4.y = ((EPJ_mass*over_r_coll_cu)*rij.y);
__fkg_tmp4.z = ((EPJ_mass*over_r_coll_cu)*rij.z);
pot_offset = (-(1.5f)/r_coll_tmp);
__fkg_tmp5 = ((0.25f*EPJ_mass)*(r_real_sq*over_r_coll_cu+pot_offset));
m_red = (EPJ_mass/(EPI_mass+EPJ_mass));
r_real = (r_real_sq*over_r_real);
dr = (r_coll_tmp-r_real);
__fkg_tmp1 = (((kappa*m_red)*dr)*over_r_real);
__fkg_tmp6.x = ((((kappa*m_red)*dr)*over_r_real)*rij.x);
__fkg_tmp6.y = ((((kappa*m_red)*dr)*over_r_real)*rij.y);
__fkg_tmp6.z = ((((kappa*m_red)*dr)*over_r_real)*rij.z);
__fkg_tmp7 = ((((0.25f*kappa)*m_red)*dr)*dr);
vij.x = (EPI_vel.x-EPJ_vel.x);
vij.y = (EPI_vel.y-EPJ_vel.y);
vij.z = (EPI_vel.z-EPJ_vel.z);
rv = (rij.z*vij.z+(rij.x*vij.x+(rij.y*vij.y)));
__fkg_tmp2 = (((eta*m_red)*rv)*over_r_real_sq);
__fkg_tmp8.x = ((((eta*m_red)*rv)*over_r_real_sq)*rij.x);
__fkg_tmp8.y = ((((eta*m_red)*rv)*over_r_real_sq)*rij.y);
__fkg_tmp8.z = ((((eta*m_red)*rv)*over_r_real_sq)*rij.z);
acc_grav_tmp.x = __fkg_tmp4.x;
acc_grav_tmp.y = __fkg_tmp4.y;
acc_grav_tmp.z = __fkg_tmp4.z;
pot_grav_tmp = __fkg_tmp5;
acc_sprg_tmp.x = __fkg_tmp6.x;
acc_sprg_tmp.y = __fkg_tmp6.y;
acc_sprg_tmp.z = __fkg_tmp6.z;
pot_sprg_tmp = __fkg_tmp7;
acc_dash_tmp.x = __fkg_tmp8.x;
acc_dash_tmp.y = __fkg_tmp8.y;
acc_dash_tmp.z = __fkg_tmp8.z;
}else{
m_over_r_real = (EPJ_mass*over_r_real);
__fkg_tmp3 = (m_over_r_real*over_r_real_sq);
__fkg_tmp4.x = ((m_over_r_real*over_r_real_sq)*rij.x);
__fkg_tmp4.y = ((m_over_r_real*over_r_real_sq)*rij.y);
__fkg_tmp4.z = ((m_over_r_real*over_r_real_sq)*rij.z);
__fkg_tmp5 = (-(0.5f)*m_over_r_real);
acc_grav_tmp.x = __fkg_tmp4.x;
acc_grav_tmp.y = __fkg_tmp4.y;
acc_grav_tmp.z = __fkg_tmp4.z;
pot_grav_tmp = __fkg_tmp5;
}
}
__fkg_tmp10 = (acc_grav_tmp.x+acc_sprg_tmp.x);
__fkg_tmp9 = (__fkg_tmp10+acc_dash_tmp.x);
FORCE_acc.x = (FORCE_acc.x+__fkg_tmp9);
__fkg_tmp12 = (acc_grav_tmp.y+acc_sprg_tmp.y);
__fkg_tmp11 = (__fkg_tmp12+acc_dash_tmp.y);
FORCE_acc.y = (FORCE_acc.y+__fkg_tmp11);
__fkg_tmp14 = (acc_grav_tmp.z+acc_sprg_tmp.z);
__fkg_tmp13 = (__fkg_tmp14+acc_dash_tmp.z);
FORCE_acc.z = (FORCE_acc.z+__fkg_tmp13);
FORCE_acc_dash.x = (FORCE_acc_dash.x+acc_dash_tmp.x);
FORCE_acc_dash.y = (FORCE_acc_dash.y+acc_dash_tmp.y);
FORCE_acc_dash.z = (FORCE_acc_dash.z+acc_dash_tmp.z);
__fkg_tmp15 = (pot_grav_tmp+pot_sprg_tmp);
FORCE_pot = (FORCE_pot+__fkg_tmp15);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+FORCE_acc.x);
force[i+0].acc.y = (force[i+0].acc.y+FORCE_acc.y);
force[i+0].acc.z = (force[i+0].acc.z+FORCE_acc.z);
force[i+0].acc_dash.x = (force[i+0].acc_dash.x+FORCE_acc_dash.x);
force[i+0].acc_dash.y = (force[i+0].acc_dash.y+FORCE_acc_dash.y);
force[i+0].acc_dash.z = (force[i+0].acc_dash.z+FORCE_acc_dash.z);
force[i+0].pot = (force[i+0].pot+FORCE_pot);
} // loop of i
} // end loop of reference 
} // Kernel_I1_J8 definition 
PIKG::F64 rsqrt(PIKG::F64 op){ return 1.0/std::sqrt(op); }
PIKG::F64 sqrt(PIKG::F64 op){ return std::sqrt(op); }
PIKG::F64 inv(PIKG::F64 op){ return 1.0/op; }
PIKG::F64 max(PIKG::F64 a,PIKG::F64 b){ return std::max(a,b);}
PIKG::F64 min(PIKG::F64 a,PIKG::F64 b){ return std::min(a,b);}
PIKG::F32 rsqrt(PIKG::F32 op){ return 1.f/std::sqrt(op); }
PIKG::F32 sqrt(PIKG::F32 op){ return std::sqrt(op); }
PIKG::F32 inv(PIKG::F32 op){ return 1.f/op; }
PIKG::S64 max(PIKG::S64 a,PIKG::S64 b){ return std::max(a,b);}
PIKG::S64 min(PIKG::S64 a,PIKG::S64 b){ return std::min(a,b);}
PIKG::S32 max(PIKG::S32 a,PIKG::S32 b){ return std::max(a,b);}
PIKG::S32 min(PIKG::S32 a,PIKG::S32 b){ return std::min(a,b);}
PIKG::F64 table(PIKG::F64 tab[],PIKG::S64 i){ return tab[i]; }
PIKG::F32 table(PIKG::F32 tab[],PIKG::S32 i){ return tab[i]; }
PIKG::F64 to_float(PIKG::U64 op){return (PIKG::F64)op;}
PIKG::F32 to_float(PIKG::U32 op){return (PIKG::F32)op;}
PIKG::F64 to_float(PIKG::S64 op){return (PIKG::F64)op;}
PIKG::F32 to_float(PIKG::S32 op){return (PIKG::F32)op;}
PIKG::S64   to_int(PIKG::F64 op){return (PIKG::S64)op;}
PIKG::S32   to_int(PIKG::F32 op){return (PIKG::S32)op;}
PIKG::U64  to_uint(PIKG::F64 op){return (PIKG::U64)op;}
PIKG::U32  to_uint(PIKG::F32 op){return (PIKG::U32)op;}
template<typename T> PIKG::F64 to_f64(const T& op){return (PIKG::F64)op;}
template<typename T> PIKG::F32 to_f32(const T& op){return (PIKG::F32)op;}
template<typename T> PIKG::S64 to_s64(const T& op){return (PIKG::S64)op;}
template<typename T> PIKG::S32 to_s32(const T& op){return (PIKG::S32)op;}
template<typename T> PIKG::U64 to_u64(const T& op){return (PIKG::U64)op;}
template<typename T> PIKG::U32 to_u32(const T& op){return (PIKG::U32)op;}
__m256 rsqrt(__m256 op){
  return _mm256_rsqrt_ps(op);
}
__m256 sqrt(__m256 op){ return _mm256_sqrt_ps(op); }
__m256 inv(__m256 op){
return _mm256_rcp_ps(op);
}
__m256d rsqrt(__m256d op){
  __m256d y = _mm256_castsi256_pd(_mm256_sub_epi64(_mm256_set1_epi64x(0x5fe6eb50c7b537a9LL),_mm256_srlv_epi64(_mm256_castpd_si256(op),_mm256_set1_epi64x(1))));
  __m256d h = _mm256_mul_pd(op,y);
  h = _mm256_fnmadd_pd(h,y,_mm256_set1_pd(1.0));
  __m256d poly = _mm256_fmadd_pd(h,_mm256_set1_pd(0.375),_mm256_set1_pd(0.5));
  poly = _mm256_mul_pd(poly,h);
  y = _mm256_fmadd_pd(y,poly,y);
  return y;
}__m256d sqrt(__m256d op){
  return _mm256_sqrt_pd(op);
}__m256d inv(__m256d op){
  __m256 x = _mm256_castps128_ps256(_mm256_cvtpd_ps(op));
  x = inv(x);
  return _mm256_cvtps_pd(_mm256_castps256_ps128(x));
}__m256d max(__m256d a,__m256d b){ return _mm256_max_pd(a,b);}
__m256d min(__m256d a,__m256d b){ return _mm256_min_pd(a,b);}
__m256  max(__m256  a,__m256  b){ return _mm256_max_ps(a,b);}
__m256  min(__m256  a,__m256  b){ return _mm256_min_ps(a,b);}
__m256i max(__m256i a,__m256i b){ return _mm256_max_epi32(a,b);}
__m256i min(__m256i a,__m256i b){ return _mm256_min_epi32(a,b);}
__m256d table(__m256d tab,__m256i index){ return _mm256_permutexvar_pd(index,tab);}
__m256  table(__m256  tab,__m256i index){ return _mm256_permutexvar_ps(index,tab);}
__m256  to_float(__m256i op){ return _mm256_cvtepi32_ps(op);}
__m256i  to_int(__m256  op){ return _mm256_cvtps_epi32(op);}
void transpose2x2_pd(__m256d& a, __m256d& b){
  __m256d tmp = _mm256_unpacklo_pd(a,b);
  b = _mm256_unpackhi_pd(a,b);
  a = tmp;
}
void transpose4x4_pd(__m256d& a, __m256d& b,__m256d& c,__m256d& d){
  __m256d tmp0 = _mm256_unpacklo_pd(a, b);
  __m256d tmp1 = _mm256_unpackhi_pd(a, b);
  __m256d tmp2 = _mm256_unpacklo_pd(c, d);
  __m256d tmp3 = _mm256_unpackhi_pd(c, d);
  a = _mm256_permute2f128_pd(tmp0, tmp2, 0|(2<<4));
  b = _mm256_permute2f128_pd(tmp1, tmp3, 0|(2<<4));
  c = _mm256_permute2f128_pd(tmp0, tmp2, 1|(3<<4));
  d = _mm256_permute2f128_pd(tmp1, tmp3, 1|(3<<4));
}
void unpack2x2_pd(__m256d& a,__m256d& b){
  transpose2x2_pd(a,b);
}
void pack2x2_pd(__m256d& a,__m256d& b){
  transpose2x2_pd(a,b);
}
void unpack4x4_pd(__m256d& a,__m256d& b,__m256d& c,__m256d& d){
  transpose4x4_pd(a,b,c,d);
}
void pack4x4_pd(__m256d& a,__m256d& b,__m256d& c,__m256d& d){
  transpose4x4_pd(a,b,c,d);
}
void unpack2x2_ps(__m256& a, __m256& b){
  __m256 tmp = _mm256_shuffle_ps(a,b,0xdd);
  b = _mm256_shuffle_ps(a,b,0x88);
  a = tmp;
}
void pack2x2_ps(__m256& a, __m256& b){
  __m256 tmp = _mm256_unpackhi_ps(a,b);
  b = _mm256_unpacklo_ps(a,b);
  a = tmp;
}
void gather4_ps(__m256& a, __m256& b,__m256& c,__m256& d){
  __m256 src0 = _mm256_permute2f128_ps(a,c,0x20); // x0 y0 z0 w0 x4 y4 z4 w4
  __m256 src1 = _mm256_permute2f128_ps(a,c,0x31); // x1 y1 z1 w1 x5 y5 z5 w5
  __m256 src2 = _mm256_permute2f128_ps(b,d,0x20);
  __m256 src3 = _mm256_permute2f128_ps(b,d,0x31);
  __m256d tmp0 = _mm256_castps_pd(_mm256_unpacklo_ps(src0, src1)); // x0 x2 y0 y2 x1 y1 x3 y3
  __m256d tmp1 = _mm256_castps_pd(_mm256_unpackhi_ps(src0, src1)); // z0 z2 w0 w2 z1 w1 z3 w3
  __m256d tmp2 = _mm256_castps_pd(_mm256_unpacklo_ps(src2, src3)); // x4 x6 y4 y6 x5 x7 y5 y7
  __m256d tmp3 = _mm256_castps_pd(_mm256_unpackhi_ps(src2, src3)); // z4 z6 w4 w6 z5 z7 w5 w7
  a = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp0, tmp2)); // x0 x2 x4 x6 x1 x3 x5 x7
  b = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp0, tmp2)); // y0 y2 y4 y6 y1 y3 y5 y7
  c = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp1, tmp3)); // z0 z2 z4 z6 z1 z3 z5 z7
  d = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp1, tmp3)); // w0 w2 w4 w6 w1 w3 w5 w7
}
void scatter4_ps(__m256& a, __m256& b,__m256& c,__m256& d){
  __m256d tmp0 = _mm256_castps_pd(_mm256_unpacklo_ps(a, b)); // x0 x2 y0 y2 x1 y1 x3 y3
  __m256d tmp1 = _mm256_castps_pd(_mm256_unpackhi_ps(a, b)); // z0 z2 w0 w2 z1 w1 z3 w3
  __m256d tmp2 = _mm256_castps_pd(_mm256_unpacklo_ps(c, d)); // x4 x6 y4 y6 x5 x7 y5 y7
  __m256d tmp3 = _mm256_castps_pd(_mm256_unpackhi_ps(c, d)); // z4 z6 w4 w6 z5 z7 w5 w7
  __m256 dst0 = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp0, tmp2)); // x0 x2 x4 x6 x1 x3 x5 x7
  __m256 dst1 = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp0, tmp2)); // y0 y2 y4 y6 y1 y3 y5 y7
  __m256 dst2 = _mm256_castpd_ps(_mm256_unpacklo_pd(tmp1, tmp3)); // z0 z2 z4 z6 z1 z3 z5 z7
  __m256 dst3 = _mm256_castpd_ps(_mm256_unpackhi_pd(tmp1, tmp3)); // w0 w2 w4 w6 w1 w3 w5 w7
  a = _mm256_permute2f128_ps(dst0,dst1,0x20);
  b = _mm256_permute2f128_ps(dst2,dst3,0x20);
  c = _mm256_permute2f128_ps(dst0,dst1,0x31);
  d = _mm256_permute2f128_ps(dst2,dst3,0x31);
}
void unpack4x4_ps(__m256& a, __m256& b,__m256& c,__m256& d){
  gather4_ps(a,b,c,d);
}
void pack4x4_ps(__m256& a, __m256& b,__m256& c,__m256& d){
  scatter4_ps(a,b,c,d);
}
void unpack4x4_ps(__m256x4& v){
  unpack4x4_ps(v.v0,v.v1,v.v2,v.v3);
}
void pack4x4_ps(__m256x4& v){
  unpack4x4_ps(v.v0,v.v1,v.v2,v.v3);
}
void unpack4x4_pd(__m256dx4& v){
  unpack4x4_pd(v.v0,v.v1,v.v2,v.v3);
}
void pack4x4_pd(__m256dx4& v){
  unpack4x4_pd(v.v0,v.v1,v.v2,v.v3);
}
void print_ps(const __m256 v){ for(int i=0;i<8;i++) printf(" %f",v[i]); printf("\n");}
void print_ps(const __m256x4 v){
  print_ps(v.v0);
  print_ps(v.v1);
  print_ps(v.v2);
  print_ps(v.v3);
}
};// kernel functor definition 
