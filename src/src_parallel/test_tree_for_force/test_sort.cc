#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>

int main(int argc, char *argv[]){

    if(argc < 3){
        std::cerr<<"too few args.; arg1: problem size, arg2: repeat count"<<std::endl;
        return 1;
    }

    PS::RadixSort<PS::U64, 8> RS;

    PS::S32 n_size = std::atoi(argv[1]);
    PS::S32 n_repeat = std::atoi(argv[2]);
    PS::TreeParticle * data = new PS::TreeParticle[n_size];
    PS::TreeParticle * data_buf = new PS::TreeParticle[n_size];
    bool * flag = new bool[n_size];
    PS::S32 err_comp = 0;
    PS::S32 err_order = 0;

    for(PS::S32 i=0; i<n_repeat; i++){
        for(int j=0; j<n_size; j++){
            data[j].setKey(PS::U64(abs(rand()))<<32 | PS::U64(abs(rand())));
            data[j].adr_ptcl_ = j;
            flag[j] = false;
        }
        RS.lsdSort(data, data_buf, 0, n_size-1);
        for(PS::S32 j=0; j<n_size; j++) flag[data[j].adr_ptcl_] = true;
        for(PS::S32 j=0; j<n_size; j++){
            if(flag[j] == false){
                err_comp++;
            }
        }
        for(PS::S32 j=1; j<n_size; j++){
            if(data[j].getKey() < data[j-1].getKey() ){
                err_order++;
            }
        }
    }

    if( err_comp || err_order){
        std::cout<<"FAIL err_comp="<<err_comp<<" err_order="<<err_order<<std::endl;
    }
    else{
        std::cout<<"PASS"<<std::endl;
    }

    return 0;
}
