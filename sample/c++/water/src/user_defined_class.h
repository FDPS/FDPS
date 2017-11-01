#ifndef H_USER_DEFINED_CLASS
#define H_USER_DEFINED_CLASS

#include "matrix3.hpp"


// =================
// Force class
// =================
class Force{
public:
  PS::F64vec acc;
  PS::F64    pot;
  void clear(){
    acc = 0.0;
    pot = 0.0;
  }
};

// =========================================
// Full particle and essential particle class
// =========================================
class FP{
public:
  PS::S64 id;
  PS::S64 type;

  PS::F64vec3 gpos;
  PS::F64vec3 pos, vel, acc;

  PS::F64 mass,pot;
  PS::F64 sigma, epsilon, charge;

  PS::F64 search_radius;

  PS::F64vec prev_pos; // for RATTLE

  PS::F64 getRsearch() const {
    return search_radius;
  }
  PS::F64vec getPos() const {
    return gpos;
  }
  void copyFromForce(const Force& force){
    acc = force.acc;
    pot = force.pot;
  }

  // for constraint and integration
  void KeepCurrentPos(){
    prev_pos = pos;
  }
  void IntegratePos(const PS::F64 dt,const PS::F64vec cell_size){
    pos += vel*dt;
  }
  void IntegrateVel(const PS::F64 dth){
    vel += acc/mass*dth;
  }
};

class EP{
public:
  PS::S64 id;
  PS::F64vec pos;

  PS::F64 sigma;
  PS::F64 epsilon;
  PS::F64 charge;

  PS::F64 search_radius;
  PS::F64 getRSearch() const {
    return search_radius;
  }
  PS::F64vec getPos() const { return pos;}
  PS::F64    getCharge() const {return charge;}
  void copyFromFP(const FP & fp){
    pos     = fp.pos;
    id      = fp.id;
    sigma   = fp.sigma;
    epsilon = fp.epsilon;
    charge  = fp.charge;
    search_radius = fp.search_radius;
  }
  void setPos(const PS::F64vec3 _pos){
    pos = _pos;
  }
};

// =================================
// Moment and Superparticle class for IPS/Tree
//==================================

class MomentUnit{
 public:
  PS::F64vec pos;
  PS::F64    mass;
  PS::F64vec di;
  PS::F64mat quad;

  MomentUnit(){
    init();
  }
  MomentUnit(const MomentUnit& m){
    pos  = m.pos;
    mass = m.mass;
    di   = m.di;
    quad = m.quad;
  }

  void init(){
    pos  = 0.0;
    mass = 0.0;
    di   = 0.0;
    quad = 0.0;
  }
  void set(){
    if(mass == 0.0) pos = 0.0;
    else            pos = pos / mass; // gravity center as a center of multipole expansion
  }

  PS::F64vec getPos() const {
    return pos;
  }
  //PS::F64ort getVertexOut() const { return vertex_out_; }

  void accumulateAtLeaf(const EP& epj){
    mass += epj.charge;
    pos  += epj.charge*epj.pos;
  }
  void accumulateAtLeaf2(const EP& epj){
    const PS::F64vec r = epj.pos - pos;
    di += epj.charge * r;
    quad.xx += epj.charge * r.x * r.x;
    quad.yy += epj.charge * r.y * r.y;
    quad.zz += epj.charge * r.z * r.z;
    quad.xy += epj.charge * r.x * r.y;
    quad.xz += epj.charge * r.x * r.z;
    quad.yz += epj.charge * r.y * r.z;
  }
  void accumulate(const MomentUnit& m){
    pos   += m.mass * m.pos;
    mass  += m.mass;
  }
  void accumulate2(const MomentUnit& m){
    const PS::F64vec r = m.pos - pos;
    di += m.di - m.mass * r;
    quad.xx += m.quad.xx - 2.0 * m.di.x * r.x + m.mass * r.x * r.x;
    quad.yy += m.quad.yy - 2.0 * m.di.y * r.y + m.mass * r.y * r.y;
    quad.zz += m.quad.zz - 2.0 * m.di.z * r.z + m.mass * r.z * r.z;
    quad.xy += m.quad.xy - (m.di.x*r.y + m.di.y*r.x) + m.mass * r.x * r.y;
    quad.xz += m.quad.xz - (m.di.x*r.z + m.di.z*r.x) + m.mass * r.x * r.z;
    quad.yz += m.quad.yz - (m.di.y*r.z + m.di.z*r.y) + m.mass * r.y * r.z;
  }
};


const PS::F64mat3asym mul(const PS::F64mat3asym& lhs,const PS::F64mat3asym& rhs){
  return PS::F64mat3asym(PS::F64vec(lhs.xx,lhs.xy,lhs.xz) * PS::F64vec(rhs.xx,rhs.yx,rhs.zx),
			 PS::F64vec(lhs.xx,lhs.xy,lhs.xz) * PS::F64vec(rhs.xy,rhs.yy,rhs.zy),
			 PS::F64vec(lhs.xx,lhs.xy,lhs.xz) * PS::F64vec(rhs.xz,rhs.yz,rhs.zz),

			 PS::F64vec(lhs.yx,lhs.yy,lhs.yz) * PS::F64vec(rhs.xx,rhs.yx,rhs.zx),
			 PS::F64vec(lhs.yx,lhs.yy,lhs.yz) * PS::F64vec(rhs.xy,rhs.yy,rhs.zy),
			 PS::F64vec(lhs.yx,lhs.yy,lhs.yz) * PS::F64vec(rhs.xz,rhs.yz,rhs.zz),

			 PS::F64vec(lhs.zx,lhs.zy,lhs.zz) * PS::F64vec(rhs.xx,rhs.yx,rhs.zx),
			 PS::F64vec(lhs.zx,lhs.zy,lhs.zz) * PS::F64vec(rhs.xy,rhs.yy,rhs.zy),
			 PS::F64vec(lhs.zx,lhs.zy,lhs.zz) * PS::F64vec(rhs.xz,rhs.yz,rhs.zz));
}

const PS::F64mat3asym trans(const PS::F64mat3asym& m){
  return PS::F64mat3asym(m.xx,m.yx,m.zx,
			 m.xy,m.yy,m.zy,
			 m.xz,m.yz,m.zz);
}

void Jacobi3x3(const PS::F64mat &A, PS::F64 lambda[3], PS::F64 U[9]){
  PS::F64mat3asym A0(A.xx,A.xy,A.xz,
		     A.xy,A.yy,A.yz,
		     A.xz,A.yz,A.zz);

  PS::F64mat3asym U0(1.0,0.0,0.0,
			  0.0,1.0,0.0,
			  0.0,0.0,1.0);
  PS::F64mat3asym R;
  PS::F64 theta,cs,sn;
  const PS::F64 eps = 1e-6;
  PS::F64 max = 1.0;
  PS::S32 count = 0;
  while(max > eps && count < 100){
    if(fabs(A0.xy)>eps){
      if(fabs(A0.xx-A0.yy)<eps) theta = 0.25*M_PI;
      else                      theta = 0.5*atan(2.0*A0.xy / (A0.xx - A0.yy));
      cs = cos(theta);
      sn = sin(theta);
      R = PS::F64mat3asym( cs, -sn, 0.0,
			   sn,  cs, 0.0,
			  0.0, 0.0, 1.0);
      A0 = mul(mul(trans(R),A0),R);
      U0 = mul(U0,R);
    }
    if(fabs(A0.yz)>eps){
      if(fabs(A0.yy-A0.zz)<eps) theta = 0.25*M_PI;
      else                      theta = 0.5*atan(2.0*A0.yz / (A0.yy - A0.zz));
      cs = cos(theta);
      sn = sin(theta);
      R = PS::F64mat3asym(1.0, 0.0, 0.0,
			  0.0,  cs, -sn,
			  0.0,  sn,  cs);
      A0 = mul(mul(trans(R),A0),R);
      U0 = mul(U0,R);
    }
    if(fabs(A0.xz)>eps){
      if(fabs(A0.xx-A0.zz)<eps) theta = 0.25*M_PI;
      else                      theta = 0.5*atan(2.0*A0.xz / (A0.xx - A0.zz));
      cs = cos(theta);
      sn = sin(theta);
      R = PS::F64mat3asym( cs, 0.0, -sn,
			  0.0, 1.0, 0.0,
			   sn, 0.0,  cs);
      A0 = mul(mul(trans(R),A0),R);
      U0 = mul(U0,R);
    }
    max = std::max(fabs(A0.xy),std::max(fabs(A0.yz),fabs(A0.zx)));
    count++;
  }
  // normalize eigenvector
  PS::F64 evl[3] = {A0.xx,A0.yy,A0.zz};
  PS::F64vec evc[3];
  evc[0] = PS::F64vec(U0.xx,U0.yx,U0.zx);
  evc[1] = PS::F64vec(U0.xy,U0.yy,U0.zy);
  evc[2] = PS::F64vec(U0.xz,U0.yz,U0.zz);
  evc[0] = evc[0] / sqrt(evc[0]*evc[0]);
  evc[1] = evc[1] / sqrt(evc[1]*evc[1]);
  evc[2] = evc[2] / sqrt(evc[2]*evc[2]);
  // rearrange eigenvalues in ascending order
  if(evl[0]>evl[1]){
    std::swap(evl[0],evl[1]);
    std::swap(evc[0],evc[1]);
  }
  if(evl[0]>evl[2]){
    std::swap(evl[0],evl[2]);
    std::swap(evc[0],evc[2]);
  }
  if(evl[1]>evl[2]){
    std::swap(evl[1],evl[2]);
    std::swap(evc[1],evc[2]);
  }

  for(int i=0;i<3;i++){
    lambda[i] = evl[i];
    for(int j=0;j<3;j++)
      U[3*i+j] = evc[i][j];
  }
}

PS::F64mat CalcQuadrupoleMoment(const PS::F64mat &quad){
  PS::F64mat qmom;

  PS::F64 tmp = quad.xx + quad.yy + quad.zz;
  qmom.xx = 1.5*quad.xx - 0.5*tmp;
  qmom.yy = 1.5*quad.yy - 0.5*tmp;
  qmom.zz = 1.5*quad.zz - 0.5*tmp;
  qmom.xy = 1.5*quad.xy;
  qmom.xz = 1.5*quad.xz;
  qmom.yz = 1.5*quad.yz;

  return qmom;
}

PS::F64mat CalcQuadrupole(const EP pp[3]){
  PS::F64mat quad = 0.0;
  for(int i=0;i<3;i++){
    quad.xx += pp[i].charge * pp[i].pos.x * pp[i].pos.x;
    quad.yy += pp[i].charge * pp[i].pos.y * pp[i].pos.y;
    quad.zz += pp[i].charge * pp[i].pos.z * pp[i].pos.z;
    quad.xy += pp[i].charge * pp[i].pos.x * pp[i].pos.y;
    quad.xz += pp[i].charge * pp[i].pos.x * pp[i].pos.z;
    quad.yz += pp[i].charge * pp[i].pos.y * pp[i].pos.z;
  }
  return quad;
}

void GeneratePseudoParticles(EP pp[3],
			     const MomentUnit& mom){
  if(mom.mass == 0.0){
    pp[0].charge = pp[1].charge = pp[2].charge = 0.0;
    pp[0].pos = pp[1].pos = pp[2].pos = 0.0;
    return;
  }
  PS::F64mat quad = CalcQuadrupoleMoment(mom.quad);
  const PS::F64 Unorm = sqrt(quad.xy*quad.xy + quad.xz*quad.xz + quad.yz*quad.yz);
  const PS::F64 Dnorm = sqrt(quad.xx*quad.xx + quad.yy*quad.yy + quad.zz*quad.zz);
  if(Unorm < 1e-8 && Dnorm < 1e-8){
    pp[0].charge = mom.mass; pp[1].charge = pp[2].charge = 0.0;
    pp[0].pos = pp[1].pos = pp[2].pos = 0.0;
    return;
  }
  if(mom.mass < 0.0){
    quad = PS::F64mat(-quad.xx,-quad.yy,-quad.zz,-quad.xy,-quad.xz,-quad.yz);
  }
  double w[3],A[9];
  Jacobi3x3(quad,w,A);

  const PS::F64 mass = fabs(mom.mass);
  assert(w[2]>=w[1] && w[1]>=w[0]);
  assert(fabs(w[2]+w[1]+w[0]) < 1e-5);
  const PS::F64 alpha = sqrt((2.0*w[2] + w[1]) / mass);
  const PS::F64 beta  =
    (w[2]+2.0*w[1] <= 0.0) ? 0.0 : sqrt((w[2] + 2.0*w[1]) / (3.0*mass));

  // center of PPs must be (0,0,0)
  pp[0].pos = PS::F64vec(   0.0, 2.0*beta, 0.0);
  pp[1].pos = PS::F64vec( alpha,    -beta, 0.0);
  pp[2].pos = PS::F64vec(-alpha,    -beta, 0.0);
  pp[0].charge = pp[1].charge = pp[2].charge = mom.mass/3.0;
  // convert to original coordinate
  const PS::F64vec e2(A[0],A[1],A[2]);
  const PS::F64vec e1(A[3],A[4],A[5]);
  const PS::F64vec e0(A[6],A[7],A[8]);
  for(int k=0;k<3;k++){
    const PS::F64vec tmp = pp[k].pos;
    pp[k].pos = tmp.x*e0 + tmp.y*e1 + tmp.z*e2;
  }
}

class MomentQuadrupole{
public:
  PS::F64    mass;
  PS::F64vec pos;
  MomentUnit positive;
  MomentUnit negative;
  PS::F64ort vertex_out_;

  MomentQuadrupole(){
    mass = 0.0;
    pos  = 0.0;
    positive.init();
    negative.init();
    vertex_out_.init();
  }
  MomentQuadrupole(const PS::F64 & _m,
		   const PS::F64vec & _p,
		   const MomentUnit & _pos,
		   const MomentUnit & _neg,
		   const PS::F64ort & v_out){
    mass = _m;
    pos  = _p;
    positive = _pos;
    negative = _neg;
    vertex_out_ = v_out;
  }
  void init(){
    mass = 0.0;
    pos  = 0.0;
    positive.init();
    negative.init();
    vertex_out_.init();
  }

  PS::F64vec getPos() const {
    return pos;
  }
  PS::F64ort getVertexOut() const { return vertex_out_; }

  void accumulateAtLeaf(EP & epj){
    if(epj.charge > 0.0){
      mass += epj.charge;
      pos  += epj.charge*epj.pos;
      positive.accumulateAtLeaf(epj);
    }else if(epj.charge < 0.0){
      mass -= epj.charge;
      pos  -= epj.charge*epj.pos;
      negative.accumulateAtLeaf(epj);
    }else{
      assert(false && "error: no charge particle exist!");
    }
    vertex_out_.merge(epj.pos, epj.getRSearch());
  }
  void accumulateAtLeaf2(EP & epj){
    if(epj.getCharge() > 0.0){
      positive.accumulateAtLeaf2(epj);
    }else if(epj.getCharge() < 0.0){
      negative.accumulateAtLeaf2(epj);
    }
  }
  void set(){
    pos = pos / mass;
    positive.set();
    negative.set();
  }
  void accumulate(const MomentQuadrupole & mom){
    const PS::F64 m = fabs(mom.mass);
    mass += m;
    pos  += m * mom.pos;
    positive.accumulate(mom.positive);
    negative.accumulate(mom.negative);
    vertex_out_.merge(mom.vertex_out_);
  }
  void accumulate2(const MomentQuadrupole & mom){
    positive.accumulate2(mom.positive);
    negative.accumulate2(mom.negative);
  }
};


class SPJQuadrupole{
public:
  EP pp[6]; // pseudo particle of positive(3) and negative(3) charge
  void copyFromMoment(const MomentQuadrupole & mom){
    GeneratePseudoParticles(pp  ,mom.positive);
    GeneratePseudoParticles(pp+3,mom.negative);
#ifdef IPS_TREE_PSEUDOPARTICLE_QUADRUPOLE_MOMENT_CHECK
    // check generated pseudo particle reproduce physical quadrupole moment
    {
      const PS::F64mat qm_ppp = CalcQuadrupoleMoment(CalcQuadrupole(pp));
      const PS::F64mat qm_mom = CalcQuadrupoleMoment(mom.positive.quad);
#if 0
      printf("pseudo:   %e %e %e %e %e %e\n",
	     qm_ppp.xx, qm_ppp.yy, qm_ppp.zz,
	     qm_ppp.xy, qm_ppp.xz, qm_ppp.yz);
      printf("original: %e %e %e %e %e %e\n",
  	     qm_mom.xx,qm_mom.yy,qm_mom.zz,
	     qm_mom.xy,qm_mom.xz,qm_mom.yz);
#endif
#if 0
      assert(fabs(qm_mom.xx) < 1e-6 || fabs((qm_mom.xx-qm_ppp.xx)/qm_mom.xx) < 1e-6);
      assert(fabs(qm_mom.yy) < 1e-6 || fabs((qm_mom.yy-qm_ppp.yy)/qm_mom.yy) < 1e-6);
      assert(fabs(qm_mom.zz) < 1e-6 || fabs((qm_mom.zz-qm_ppp.zz)/qm_mom.zz) < 1e-6);
      assert(fabs(qm_mom.xy) < 1e-6 || fabs((qm_mom.xy-qm_ppp.xy)/qm_mom.xy) < 1e-6);
      assert(fabs(qm_mom.xz) < 1e-6 || fabs((qm_mom.xz-qm_ppp.xz)/qm_mom.xz) < 1e-6);
      assert(fabs(qm_mom.yz) < 1e-6 || fabs((qm_mom.yz-qm_ppp.yz)/qm_mom.yz) < 1e-6);
#else
      assert(fabs(qm_mom.xx) < 1e-6 || fabs(qm_mom.xx-qm_ppp.xx) < 1e-5);
      assert(fabs(qm_mom.yy) < 1e-6 || fabs(qm_mom.yy-qm_ppp.yy) < 1e-5);
      assert(fabs(qm_mom.zz) < 1e-6 || fabs(qm_mom.zz-qm_ppp.zz) < 1e-5);
      assert(fabs(qm_mom.xy) < 1e-6 || fabs(qm_mom.xy-qm_ppp.xy) < 1e-5);
      assert(fabs(qm_mom.xz) < 1e-6 || fabs(qm_mom.xz-qm_ppp.xz) < 1e-5);
      assert(fabs(qm_mom.yz) < 1e-6 || fabs(qm_mom.yz-qm_ppp.yz) < 1e-5);
#endif
    }
    {
      const PS::F64mat qm_ppn = CalcQuadrupoleMoment(CalcQuadrupole(pp+3));
      const PS::F64mat qm_mom = CalcQuadrupoleMoment(mom.negative.quad);
#if 0
      fprintf(stderr,
	      "pseudo:   %e %e %e %e %e %e\n",
	      qm_ppn.xx, qm_ppn.yy, qm_ppn.zz,
	      qm_ppn.xy, qm_ppn.xz, qm_ppn.yz);
      fprintf(stderr,
	      "original: %e %e %e %e %e %e\n",
	      qm_mom.xx,qm_mom.yy,qm_mom.zz,
	      qm_mom.xy,qm_mom.xz,qm_mom.yz);
#endif
#if 0
      assert(fabs(qm_mom.xx) < 1e-6 || fabs((qm_mom.xx-qm_ppn.xx)/qm_mom.xx) < 1e-6);
      assert(fabs(qm_mom.yy) < 1e-6 || fabs((qm_mom.yy-qm_ppn.yy)/qm_mom.yy) < 1e-6);
      assert(fabs(qm_mom.zz) < 1e-6 || fabs((qm_mom.zz-qm_ppn.zz)/qm_mom.zz) < 1e-6);
      assert(fabs(qm_mom.xy) < 1e-6 || fabs((qm_mom.xy-qm_ppn.xy)/qm_mom.xy) < 1e-6);
      assert(fabs(qm_mom.xz) < 1e-6 || fabs((qm_mom.xz-qm_ppn.xz)/qm_mom.xz) < 1e-6);
      assert(fabs(qm_mom.yz) < 1e-6 || fabs((qm_mom.yz-qm_ppn.yz)/qm_mom.yz) < 1e-6);
#else
      assert(fabs(qm_mom.xx) < 1e-6 || fabs(qm_mom.xx-qm_ppn.xx) < 1e-5);
      assert(fabs(qm_mom.yy) < 1e-6 || fabs(qm_mom.yy-qm_ppn.yy) < 1e-5);
      assert(fabs(qm_mom.zz) < 1e-6 || fabs(qm_mom.zz-qm_ppn.zz) < 1e-5);
      assert(fabs(qm_mom.xy) < 1e-6 || fabs(qm_mom.xy-qm_ppn.xy) < 1e-5);
      assert(fabs(qm_mom.xz) < 1e-6 || fabs(qm_mom.xz-qm_ppn.xz) < 1e-5);
      assert(fabs(qm_mom.yz) < 1e-6 || fabs(qm_mom.yz-qm_ppn.yz) < 1e-5);
#endif
    }
#endif
    for(int i=0;i<3;i++){
      pp[i  ].pos += mom.positive.pos;
      pp[i+3].pos += mom.negative.pos;
    }
  }
  void clear(){
    for(int i=0;i<6;i++){
      pp[i].charge = 0.0; pp[i].pos = 0.0;
    }
  }

  PS::F32vec getPos() const {
    PS::F64    mass = 0.0;
    PS::F64vec pos  = 0.0;
    for(int i=0;i<6;i++){
      mass += fabs(pp[i].charge);
      pos  += fabs(pp[i].charge)*pp[i].pos;
    }
    return pos/mass;
  }
  void setPos(const PS::F64vec & pos_new) {
    const PS::F64vec diff = pos_new - getPos();
    for(int i=0;i<6;i++){
      if(pp[i].charge != 0.0) pp[i].pos += diff;
    }
  }
  MomentQuadrupole convertToMoment() const {
    PS::F64    mass = 0.0;
    for(int i=0;i<6;i++) mass += fabs(pp[i].getCharge());
    PS::F64vec pos = getPos();
    MomentUnit positive,negative;
    for(int i=0;i<3;i++){
      positive.accumulateAtLeaf(pp[i]);
      negative.accumulateAtLeaf(pp[i+3]);
    }
    positive.set();
    negative.set();

    PS::F64ort vtx; // dummy vertex
    vtx.init();
    return MomentQuadrupole(mass,pos,positive,negative,vtx);
  }
};
#endif
