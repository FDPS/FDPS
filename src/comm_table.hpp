#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{
    template<class Tepj, class Tspj>
    class CommTable{
    public:
        S32 n_ep_send_tot_, n_ep_recv_tot_, n_sp_send_tot_, n_sp_recv_tot_;
        ReallocatableArray<S32> n_ep_send_;
        ReallocatableArray<S32> n_ep_recv_;
        ReallocatableArray<S32> adr_ep_send_;
        //ReallocatableArray<S32> adr_ep_recv_;
        
        ReallocatableArray<S32> n_sp_send_;
        ReallocatableArray<S32> n_sp_recv_;
        ReallocatableArray<S32> adr_sp_send_;
        //ReallocatableArray<S32> adr_sp_recv_;

        ReallocatableArray<F64vec> shift_per_image_;
        ReallocatableArray<S32> n_image_per_proc_;

        ReallocatableArray<S32> n_ep_per_image_;
        ReallocatableArray<S32> n_sp_per_image_;
        
    public:
        void initialize(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_ep_send_.resizeNoInitialize(n_proc);
            n_ep_recv_.resizeNoInitialize(n_proc);
            n_sp_send_.resizeNoInitialize(n_proc);
            n_sp_recv_.resizeNoInitialize(n_proc);
            n_ep_send_tot_ = n_ep_recv_tot_ = n_sp_send_tot_ = n_sp_recv_tot_ = 0;
        }

        void clearSize(){
            adr_ep_send_.clearSize();
            //adr_ep_recv_.clearSize();
            adr_sp_send_.clearSize();
            //adr_sp_recv_.clearSize();
            shift_per_image_.clearSize();
            n_ep_per_image_.clearSize();
            n_sp_per_image_.clearSize();
        }
#if 0
        void setDispSend(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_disp_ep_send_[0] = n_disp_sp_send_[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_[i+1] = n_disp_ep_send_[i] + n_ep_send_[i];
                n_disp_sp_send_[i+1] = n_disp_sp_send_[i] + n_sp_send_[i];
            }
        }

        void setDispRecv(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_disp_ep_recv_[0] = n_disp_sp_recv_[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_recv_[i+1] = n_disp_ep_recv_[i] + n_ep_recv_[i];
                n_disp_sp_recv_[i+1] = n_disp_sp_recv_[i] + n_sp_recv_[i];
            }
        }

        void setSizeSendBuf(){
            const S32 n_proc = Comm::getNumberOfProc();
            ep_send_.resizeNoInitialize( n_disp_ep_send_[n_proc] );
            sp_send_.resizeNoInitialize( n_disp_sp_send_[n_proc] );
        }

        void setSizeRecvBuf(){
            const S32 n_proc = Comm::getNumberOfProc();
            ep_recv_.resizeNoInitialize( n_disp_ep_recv_[n_proc] );
            sp_recv_.resizeNoInitialize( n_disp_sp_recv_[n_proc] );
        }

        void exchangeLet(const bool reuse = false){
            if(!reuse){
                const S32 n_proc = Comm::getNumberOfProc();
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i] = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    n_ep_recv_[i] = n_ep_sp_recv_[2*i];
                    n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
                }
                setDispRecv();
                setSizeRecvBuf();
            }
            Comm::allToAllV(ep_send_.getPointer(), n_ep_send_, n_disp_ep_send_,
                            ep_recv_.getPointer(), n_ep_recv_, n_disp_ep_recv_);
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
        }

        void setSPTop(const Tspj & _sp_top){
            sp_top_ = _sp_top;
        }

        void allgahterTreeOuterPos(const F64ort & _tree_outer_pos){
            Comm::allGather(&_tree_outer_pos, 1, tree_outer_pos_.getPointer());
        }

        void dumpMemSizeUsed(std::ostream & fout){
            fout<<"ep_send_.getMemSize()= "<<ep_send_.getMemSize()<<std::endl;
            fout<<"sp_send_.getMemSize()= "<<sp_send_.getMemSize()<<std::endl;
            fout<<"ep_recv_.getMemSize()= "<<ep_recv_.getMemSize()<<std::endl;
            fout<<"sp_recv_.getMemSize()= "<<sp_recv_.getMemSize()<<std::endl;
            fout<<"rank_sp_isend_.getMemSize()= "<<rank_sp_isend_.getMemSize()<<std::endl;
            fout<<"rank_sp_irecv_.getMemSize()= "<<rank_sp_irecv_.getMemSize()<<std::endl;
            fout<<"rank_sp_top_recv_.getMemSize()= "<<rank_sp_top_recv_.getMemSize()<<std::endl;
            fout<<"rank_ep_send_.getMemSize()= "<<rank_ep_send_.getMemSize()<<std::endl;
            fout<<"rank_ep_recv_.getMemSize()= "<<rank_ep_recv_.getMemSize()<<std::endl;
            fout<<"sum= "<<
                (double)(ep_send_.getMemSize()
                         +sp_send_.getMemSize()
                         +ep_recv_.getMemSize()
                         +sp_recv_.getMemSize()
                         +rank_sp_isend_.getMemSize()
                         +rank_sp_irecv_.getMemSize()
                         +rank_sp_top_recv_.getMemSize()
                         +rank_ep_send_.getMemSize()
                         +rank_ep_recv_.getMemSize()) / 1e9
                <<" [GB]"
                <<std::endl;
        }
#endif
    };
}
