struct CalcForceEpEpImpl{
  PS::F32 eps2;
  PS::F32 kappa;
  PS::F32 eta;
  CalcForceEpEpImpl(){}
  CalcForceEpEpImpl(PIKG::F32 eps2,PIKG::F32 kappa,PIKG::F32 eta):eps2(eps2),kappa(kappa),eta(eta){}
  void initialize(PIKG::F32 eps2_,PIKG::F32 kappa_,PIKG::F32 eta_){
    eps2 = eps2_;
    kappa = kappa_;
    eta = eta_;
  }
  template<typename Tepi, typename Tepj, typename Tforce>
  void operator()(const Epi0* __restrict__ epi,const int ni,const Epj0* __restrict__ epj,const int nj,Force0* __restrict__ force,const int kernel_select = 1){  
  }
  
};
