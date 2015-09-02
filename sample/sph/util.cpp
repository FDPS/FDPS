#include "header.h"

void DisplayInfo(void){
	if(PS::Comm::getRank() == 0){
		std::cout << "//==================================\\\\" << std::endl;
		std::cout << "||                                  ||" << std::endl;
		std::cout << "|| ::::::: ::::::. ::::::. .::::::. ||" << std::endl;
		std::cout << "|| ::      ::    : ::    : ::       ||" << std::endl;
		std::cout << "|| ::::::  ::    : ::::::'  `:::::. ||" << std::endl;
		std::cout << "|| ::      ::::::' ::      `......' ||" << std::endl;
		std::cout << "||     Framework for Developing     ||" << std::endl;
		std::cout << "||        Particle Simulator        ||" << std::endl;
		std::cout << "\\\\==================================//" << std::endl;
		std::cout << "//=====================================" << std::endl;
		std::cout << "This is a sample program of Smoothed Particle Hydrodynamics on FDPS!" << std::endl;
		std::cout << "# of proc is   " << PS::Comm::getNumberOfProc() << std::endl;
		#ifdef _OPENMP
//		std::cout << "# of thread is " << omp_get_num_threads() << std::endl;
		std::cout << "# of thread is " << PS::Comm::getNumberOfThread() << std::endl;
		#else 
		std::cout << "No OpenMP Mode... " << std::endl;
		#endif
		std::cout << "//=====================================" << std::endl;
	}
	return ;
}

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system){
	PS::F64vec Mom;//total momentum
	PS::F64    Eng;//total enegry
	Mom = 0;
	Eng = 0;
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		Mom += sph_system[i].vel * sph_system[i].mass;
		Eng += (sph_system[i].eng + 0.5 * sph_system[i].vel * sph_system[i].vel)  * sph_system[i].mass;
	}
	Eng = PS::Comm::getSum(Eng);
	Mom = PS::Comm::getSum(Mom);
    if(PS::Comm::getRank() == 0){
        printf("%.16e\n", Eng);
        printf("%.16e\n", Mom.x);
        printf("%.16e\n", Mom.y);
        printf("%.16e\n", Mom.z);
    }
}
