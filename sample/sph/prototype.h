#pragma once

void DisplayInfo(void);
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64*, boundary*);

void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system);

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system);
