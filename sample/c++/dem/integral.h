#pragma once

template <class ThisPtcl> double GetGlobalTimestep(const PS::ParticleSystem<ThisPtcl>& ptcl){
	PS::F64 dt = ptcl[0].dt;
	for(PS::S64 i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
		//timestep determined by contact
		dt = std::min(ptcl[i].dt, dt);
		//timestep determined by velocity
		//dt = std::min(1.0 / 50.0 * ptcl[i].rad / sqrt(ptcl[i].vel * ptcl[i].vel), dt);
		//timestep determined by acc
		//dt = std::min(1.0 / 50.0 * sqrt(ptcl[i].rad / sqrt(ptcl[i].acc * ptcl[i].acc)), dt);
	}
	return PS::Comm::getMinValue(dt);
}

