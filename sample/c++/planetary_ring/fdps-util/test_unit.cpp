#include<iostream>
#include"particle_simulator.hpp"
#include"fdps-util.hpp"
#include"kepler.hpp"
using namespace MY_LIB::LITERAL;
int main(){
    const auto time = MY_LIB::CalcTimeUnit(1.0_au2cm, 1.0_msun2g);
    std::cerr<<"time= "<<time<<" 2.0*3.1415*time= "<<2.0*3.1416*time<<std::endl;
    
}
