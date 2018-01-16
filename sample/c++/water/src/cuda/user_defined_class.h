#ifndef H_USER_DEFINED_CLASS
#define H_USER_DEFINED_CLASS

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

// ====================
// Full particle class
// ====================
class FP{
public:
  PS::S64 id;
  PS::S64 type;
  PS::S64 cid[3]; // constraint id of H1, O, H2

  PS::F64vec pos;
  PS::F64vec vel;
  PS::F64vec acc;

  PS::F64 mass;
  PS::F64 sigma;
  PS::F64 epsilon;
  PS::F64 charge;

  PS::F64 pot;

  PS::F64 search_radius;

  PS::F64 getRsearch() const {
    return search_radius;
  }
  PS::F64vec getPos() const {
    return pos;
  }
  void copyFromForce(const Force& force){
    acc = force.acc;
    pot = force.pot;
  }

  void IntegratePos(const PS::F64 dt,const PS::F64vec cell_size){
    pos += vel*dt;
  }
  void IntegrateVel(const PS::F64 dth){
    vel += acc/mass*dth;
  }
};

// ====================
// Essential particle class
// ====================
class EP{
public:
  PS::S64 id;
  PS::F64vec pos;

  PS::F64 sigma;
  PS::F64 epsilon;
  PS::F64 charge;

  PS::F64 search_radius;
  PS::F64 getRSearch() const {return search_radius;}
  PS::F64vec getPos()  const {return pos;}
  PS::F64 getCharge()  const {return charge;}
  PS::S64 getId()      const {return id;}
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

#endif
