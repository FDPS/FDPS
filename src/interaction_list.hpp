#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{
    class InteractionList{
    public:
        ReallocatableArray<S32> n_ep_;
        ReallocatableArray<S32> n_sp_;
        ReallocatableArray<S32> n_disp_ep_;
        ReallocatableArray<S32> n_disp_sp_;
        ReallocatableArray<S32> adr_ep_;
        ReallocatableArray<S32> adr_sp_;
        void clearSize(){
            n_ep_.clearSize();
            n_sp_.clearSize();
            n_disp_ep_.clearSize();
            n_disp_sp_.clearSize();
            adr_ep_.clearSize();
            adr_sp_.clearSize();
        }
    };
}
