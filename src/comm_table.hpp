#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{

    class CommTable{
    public:
        //CommInfo comm_info_;
        S32 n_ep_send_tot_, n_ep_recv_tot_, n_sp_send_tot_, n_sp_recv_tot_;
        ReallocatableArray<S32> n_ep_send_;
        ReallocatableArray<S32> n_ep_recv_;
        ReallocatableArray<S32> adr_ep_send_;
        
        ReallocatableArray<S32> n_sp_send_;
        ReallocatableArray<S32> n_sp_recv_;
        ReallocatableArray<S32> adr_sp_send_;

        ReallocatableArray<F64vec> shift_per_image_;
        ReallocatableArray<S32> n_image_per_proc_;

        ReallocatableArray<S32> n_ep_per_image_;
        ReallocatableArray<S32> n_sp_per_image_;

        ReallocatableArray<S32> rank_send_;
        ReallocatableArray<S32> rank_recv_;

        ReallocatableArray<S32> rank_recv_allgather_;
        ReallocatableArray<F64ort> pos_domain_allgather_;
        
    public:
        /*
        CommTable(){
            comm_info_.setCommunicator();
        }
        void setCommInfo(const CommInfo & c){
            comm_info_ = c;
        }
        */
        void setPosDomainAllgather(const DomainInfo & dinfo,
                                   const F64ort & pos_root_cell){
            pos_domain_allgather_.resizeNoInitialize(rank_recv_allgather_.size());
            for(S32 i=0; i<pos_domain_allgather_.size(); i++){
                pos_domain_allgather_[i] = dinfo.getPosDomain(pos_root_cell, rank_recv_allgather_[i]);
            }
        }
        void clearSize(){
            adr_ep_send_.clearSize();
            adr_sp_send_.clearSize();
            shift_per_image_.clearSize();
            n_ep_per_image_.clearSize();
            n_sp_per_image_.clearSize();
        }
        void dump(std::ostream &fout=std::cerr){
            fout<<"n_ep_send_tot_="<<n_ep_send_tot_<<std::endl;
            fout<<"n_ep_recv_tot_="<<n_ep_recv_tot_<<std::endl;
            fout<<"n_sp_send_tot_="<<n_sp_send_tot_<<std::endl;
            fout<<"n_sp_recv_tot_="<<n_sp_recv_tot_<<std::endl;
            fout<<"rank_send_.size()="<<rank_send_.size()<<std::endl;
            fout<<"rank_recv_.size()="<<rank_recv_.size()<<std::endl;
            fout<<"rank_recv_allgather_.size()="<<rank_recv_allgather_.size()<<std::endl;
            fout<<"pos_domain_allgather_.size()="<<pos_domain_allgather_.size()<<std::endl;

        }
    };
}
