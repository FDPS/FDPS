#pragma once

//#define UNORDERED_SET 

#if __cplusplus <= 199711L
#include<map>
#else
#include<unordered_map>
#endif

#include<sort.hpp>
#include<tree.hpp>
#include<comm_table.hpp>
#include<interaction_list.hpp>
#include<tree_walk.hpp>
#include<tree_for_force_utils.hpp>



namespace ParticleSimulator{

    ///////////////////////////////////////
    /// TREE FOR FORCE CLASS DEFINITION ///
    template<
        class TSM, // search mode
        class Tforce, // USER def
        class Tepi, // USER def
        class Tepj, // USER def
        class Tmomloc, // PS or USER def
        class Tmomglb, // PS or USER def
        class Tspj // PS or USER def
        >
    class TreeForForce{

    public:
        F64 length_; // length of a side of the root cell
        F64vec center_; // new member (not used)
        //F64ort pos_root_cell_;

        //private:
    public:

        TimeProfile time_profile_;
#if __cplusplus <= 199711L
        typedef std::map<S64, Tepj*> MyMap;
#else
        typedef std::unordered_map<S64, Tepj*> MyMap;
#endif

        MyMap map_id_to_epj_;
        
        CountT n_interaction_ep_ep_local_, n_interaction_ep_sp_local_, n_walk_local_;
        CountT n_let_ep_send_1st_, n_let_ep_recv_1st_, n_let_sp_send_1st_, n_let_sp_recv_1st_,
            n_let_ep_send_2nd_, n_let_ep_recv_2nd_;
        CountT * n_cell_open_;
        CountT n_proc_send_exchange_LET_1st__icomm_sp_, n_proc_recv_exchange_LET_1st__icomm_sp_, 
            n_proc_send_exchange_LET_1st__icomm_ep_, n_proc_recv_exchange_LET_1st__icomm_ep_; 

        //F64 Tcomm_tmp_;
        F64 wtime_exlet_comm_;
        F64 wtime_exlet_a2a_;
        F64 wtime_exlet_a2av_;
        F64 Tcomm_scatterEP_tmp_;
        //F64 TexLET0_, TexLET1_, TexLET2_, TexLET3_, TexLET4_;
        //S64 n_ep_send_1st_, n_ep_recv_1st_, n_ep_send_2nd_, n_ep_recv_2nd_;
        F64 wtime_walk_LET_1st_, wtime_walk_LET_2nd_;

        bool is_initialized_;

        //S64 n_interaction_;
        S64 n_interaction_ep_ep_;
        S64 n_interaction_ep_sp_;
        S32 ni_ave_;
        S32 nj_ave_;

        RadixSort<U64, 8> rs_;
        S32 n_loc_tot_; // # of all kinds of assigned particles in local proc
        S64 n_glb_tot_; // n_loc_tot_ + LETs
        S32 n_leaf_limit_;
        S32 n_group_limit_;

        //S32 adr_tc_level_partition_[TREE_LEVEL_LIMIT+2];
        //S32 lev_max_;
        F64 theta_;

        //F64 length_; // length of a side of the root cell
        //F64vec center_; // new member (not used)
        F64ort pos_root_cell_;

        ReallocatableArray< TreeParticle > tp_buf_, tp_loc_, tp_glb_;
        ReallocatableArray< TreeCell< Tmomloc > > tc_loc_;
        ReallocatableArray< TreeCell< Tmomglb > > tc_glb_;
        ReallocatableArray< Tepi > epi_sorted_, epi_org_;
        ReallocatableArray< Tepj > epj_sorted_, epj_org_;
        ReallocatableArray< Tspj > spj_sorted_, spj_org_;
        ReallocatableArray< Tspj > spj_sorted_loc_; // for reuse the list
        ReallocatableArray< IPGroup<TSM> > ipg_;

        ReallocatableArray<Tepj> epj_send_;
        ReallocatableArray<Tepj> epj_recv_;
        ReallocatableArray<Tspj> spj_send_;
        ReallocatableArray<Tspj> spj_recv_;

        S32 * n_ep_send_; // * n_proc
        S32 * n_sp_send_; // * n_proc
        S32 * n_ep_send_disp_; // * n_proc+1
        S32 * n_sp_send_disp_; // * n_proc+1
        S32 * n_ep_recv_; // * n_proc
        S32 * n_sp_recv_; // * n_proc
        S32 * n_ep_recv_disp_; // * n_proc+1
        S32 * n_sp_recv_disp_; // * n_proc+1

        S32 * n_ep_sp_send_; // 2 * n_proc: even id is # of EP and odd id is # of SP
        S32 * n_ep_sp_recv_; // 2 * n_proc: even id is # of EP and odd id is # of SP


        ReallocatableArray<S32> * id_ep_send_buf_;
        ReallocatableArray<S32> * id_sp_send_buf_;

        S32 ** id_proc_send_; // id_proc_send_[n_thread][n_proc]

        ReallocatableArray<Tforce> force_sorted_;
        ReallocatableArray<Tforce> force_org_;

        ReallocatableArray<Tepj> * epj_for_force_;
        ReallocatableArray<Tspj> * spj_for_force_;

        ReallocatableArray<S32> * id_epj_for_force_;
        ReallocatableArray<S32> * id_spj_for_force_;

        S32 n_surface_for_comm_;



        // new variables for commnuication of LET
        // for scatterEP
        ReallocatableArray<Tepj> * ep_send_buf_for_scatter_;
        ReallocatableArray<F64vec> * shift_image_domain_;


        //PROFILE::Profile profile;

        // for gather mode
        class EPJWithR{
            Tepj epj_;
            F64 r_search_;
        public:
            Tepj getEPJ() const { return epj_; }
            F64vec getPos() const { return epj_.getPos(); }
            F64 getCharge() const { return epj_.getCharge(); }
            F64 getRSearch() const { return r_search_; }
            void copyFromEPJ(const Tepj & epj){ epj_ = epj; }
            void setRSearch(const F64 r_search){ r_search_ = r_search; }
            void setPos(const F64vec & pos){ epj_.setPos(pos);}
        };
        ReallocatableArray<EPJWithR> epjr_sorted_; // cirectly copied from EPJ + RSearch
        ReallocatableArray<EPJWithR> epjr_send_;
        ReallocatableArray<EPJWithR> epjr_recv_;
        ReallocatableArray<EPJWithR> epjr_recv_1st_buf_;
        ReallocatableArray<EPJWithR> epjr_recv_2nd_buf_;
        ReallocatableArray<EPJWithR> * epjr_send_buf_; // for 1st communication
        ReallocatableArray<EPJWithR> * epjr_send_buf_for_scatter_;
        ReallocatableArray<EPJWithR> * epjr_recv_1st_sorted_;



        // for symmeteric mode
        ReallocatableArray<Tepj> epj_recv_1st_buf_;
        ReallocatableArray<Tepj> epj_recv_2nd_buf_;
        ReallocatableArray<Tepj> * epj_send_buf_; // for 1st communication
        S32 * n_epj_recv_1st_;
        S32 * n_epj_recv_disp_1st_;
        S32 * n_epj_recv_2nd_;
        S32 * n_epj_recv_disp_2nd_;
        S32 * id_proc_src_;
        S32 * id_proc_dest_;
        ReallocatableArray<S32> * id_ptcl_send_;
        ReallocatableArray<F64vec> * shift_image_box_;
        ReallocatableArray<S32> * ip_disp_;
        // new val
        ReallocatableArray< TreeParticle > * tp_scatter_;
        ReallocatableArray< TreeCell< Tmomloc > > * tc_recv_1st_;
        ReallocatableArray<Tepj> * epj_recv_1st_sorted_;
        S32 ** adr_tc_level_partition_recv_1st_;
        ReallocatableArray<F64ort> tree_top_outer_pos_;
        ReallocatableArray<F64ort> tree_top_inner_pos_;
        
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Request * req_send_;
        MPI_Request * req_recv_;
        MPI_Status  * status_;
#endif

        template<class Tep2, class Tep3>
        inline void scatterEP(S32 n_send[],
                              S32 n_send_disp[],
                              S32 n_recv[],
                              S32 n_recv_disp[],
                              ReallocatableArray<Tep2> & ep_send,  // send buffer
                              ReallocatableArray<Tep2> & ep_recv,  // recv buffer
                              ReallocatableArray<Tep2> * ep_end_buf,  // send buffer
                              const ReallocatableArray<Tep3> & ep_org, // original
                              const DomainInfo & dinfo);

        template<class Tep2, class Tep3>
        inline void scatterEPForGather(S32 n_send[],
                                       S32 n_send_disp[],
                                       S32 n_recv[],
                                       S32 n_recv_disp[],
                                       ReallocatableArray<Tep2> & ep_send,  // send buffer
                                       ReallocatableArray<Tep2> & ep_recv,  // recv buffer
                                       const ReallocatableArray<Tep3> & ep_org, // original
                                       const DomainInfo & dinfo);

        void calcMomentLocalTreeOnlyImpl(TagSearchLong);
        void calcMomentLocalTreeOnlyImpl(TagSearchLongCutoff);
        void calcMomentLocalTreeOnlyImpl(TagSearchLongScatter);
        void calcMomentLocalTreeOnlyImpl(TagSearchLongSymmetry);
        void calcMomentLocalTreeOnlyImpl(TagSearchLongCutoffScatter);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortScatter);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortGather);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortSymmetry);

        void exchangeLocalEssentialTreeImpl(TagSearchLong, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchLongCutoff, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchLongScatter, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchLongSymmetry, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchLongCutoffScatter, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortScatter, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortGather, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeGatherImpl(TagRSearch, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeGatherImpl(TagNoRSearch, const DomainInfo & dinfo);

        void setLocalEssentialTreeToGlobalTreeImpl(TagForceShort);
        void setLocalEssentialTreeToGlobalTreeImpl(TagForceLong);

// probably, calcMomentGlobalTreeOnlyImpl is classified depending on force_type.
        void calcMomentGlobalTreeOnlyImpl(TagSearchLong);
        void calcMomentGlobalTreeOnlyImpl(TagSearchLongCutoff);
        void calcMomentGlobalTreeOnlyImpl(TagSearchLongScatter);
        void calcMomentGlobalTreeOnlyImpl(TagSearchLongSymmetry);
        void calcMomentGlobalTreeOnlyImpl(TagSearchShortScatter);
        void calcMomentGlobalTreeOnlyImpl(TagSearchShortGather);
        void calcMomentGlobalTreeOnlyImpl(TagSearchShortSymmetry);

        void makeIPGroupImpl(TagForceLong);
        void makeIPGroupImpl(TagForceShort);

        void makeInteractionListImpl(TagSearchLong, const S32 adr_ipg, const bool clear);
        void makeInteractionListImpl(TagSearchLongCutoff, const S32 adr_ipg, const bool clear);
        void makeInteractionListImpl(TagSearchLongScatter, const S32 adr_ipg, const bool clear);
        void makeInteractionListImpl(TagSearchLongSymmetry, const S32 adr_ipg, const bool clear);
        void makeInteractionListImpl(TagSearchShortScatter, const S32 adr_ipg, const bool clear);
        void makeInteractionListImpl(TagSearchShortGather, const S32 adr_ipg, const bool clear);
        void makeInteractionListImpl(TagSearchShortSymmetry, const S32 adr_ipg, const bool clear);

        void makeInteractionListIndexImpl(TagSearchLong, const S32 adr_ipg, const bool clear);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkImpl(TagForceLong,
                                   Tfunc_dispatch pfunc_dispatch,
                                   Tfunc_retrieve pfunc_retrieve,
                                   const S32 tag_max,
                                   const S32 n_walk_limit,
                                   const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkImpl(TagForceShort,
                                   Tfunc_dispatch pfunc_dispatch,
                                   Tfunc_retrieve pfunc_retrieve,
                                   const S32 tag_max,
                                   const S32 n_walk_limit,
                                   const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndexImpl(TagForceLong,
                                        Tfunc_dispatch pfunc_dispatch,
                                        Tfunc_retrieve pfunc_retrieve,
                                        const S32 tag_max,
                                        const S32 n_walk_limit,
                                        const bool flag_keep_list,
                                        const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndexImpl(TagForceShort,
                                        Tfunc_dispatch pfunc_dispatch,
                                        Tfunc_retrieve pfunc_retrieve,
                                        const S32 tag_max,
                                        const S32 n_walk_limit,
                                        const bool flag_keep_list,
                                        const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndexNewImpl(TagForceLong,
                                           Tfunc_dispatch pfunc_dispatch,
                                           Tfunc_retrieve pfunc_retrieve,
                                           const S32 tag_max,
                                           const S32 n_walk_limit,
                                           const bool flag_keep_list,
                                           const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndexNewImpl(TagForceShort,
                                           Tfunc_dispatch pfunc_dispatch,
                                           Tfunc_retrieve pfunc_retrieve,
                                           const S32 tag_max,
                                           const S32 n_walk_limit,
                                           const bool flag_keep_list,
                                           const bool clear=true);


        void checkMakeGlobalTreeImpl(TagForceLong, S32 & err, const F64vec & center, const F64 tolerance, std::ostream & fout);
        void checkMakeGlobalTreeImpl(TagForceShort, S32 & err, const F64vec & center, const F64 tolerance, std::ostream & fout);

        void checkCalcMomentLocalTreeImpl(TagSearchLong, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentLocalTreeImpl(TagSearchLongCutoff, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentLocalTreeImpl(TagSearchShortScatter, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentLocalTreeImpl(TagSearchShortGather, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentLocalTreeImpl(TagSearchShortSymmetry, const F64 tolerance, std::ostream & fout);
        
        void checkCalcMomentGlobalTreeImpl(TagSearchLong, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentGlobalTreeImpl(TagSearchLongCutoff, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentGlobalTreeImpl(TagSearchShortScatter, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentGlobalTreeImpl(TagSearchShortGather, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentGlobalTreeImpl(TagSearchShortSymmetry, const F64 tolerance, std::ostream & fout);


        ///////////////////////
        ///// for open boundary
        // for P^3T
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLongScatter){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLongSymmetry){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        // for P^3T
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLongCutoffScatter){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchShortScatter){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchShortGather){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epi_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchShortSymmetry){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epi_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLongCutoff){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLong){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        template<class Tep2>
        void calcCenterAndLengthOfRootCellOpenNoMargenImpl(const Tep2 ep[]);
        template<class Tep2>
        void calcCenterAndLengthOfRootCellOpenWithMargenImpl(const Tep2 ep[]);

        /////////////
        //// PERIODIC
        // for P^3T
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLongCutoffScatter){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchShortScatter){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchShortGather){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epi_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchShortSymmetry){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epi_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLongCutoff){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epj_org_.getPointer());
        }
        ///////////////
        // for compile 
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLong){}
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLongScatter){}
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLongSymmetry){}

        template<class Tep2>
        void calcCenterAndLengthOfRootCellPeriodicImpl2(const Tep2 ep[]);

        void checkMortonSortGlobalTreeOnlyImpl(TagForceLong, std::ostream & fout);
        void checkMortonSortGlobalTreeOnlyImpl(TagForceShort, std::ostream & fout);

        void checkMakeInteractionListImpl(TagSearchLong,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkMakeInteractionListImpl(TagSearchLongCutoff,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkMakeInteractionListImpl(TagSearchShortScatter,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkMakeInteractionListImpl(TagSearchShortGather,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkMakeInteractionListImpl(TagSearchShortSymmetry,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkExchangeLocalEssentialTreeImpl(TagForceLong,
                                                 const DomainInfo & dinfo,
                                                 const F64 tolerance,
                                                 std::ostream & fout);
        void checkExchangeLocalEssentialTreeForLongImpl(TagSearchLong,
                                                        const DomainInfo & dinfo,
                                                        const F64 tolerance, 
                                                        std::ostream & fout);
        void checkExchangeLocalEssentialTreeForLongImpl(TagSearchLongCutoff,
                                                        const DomainInfo & dinfo,
                                                        const F64 tolerance, 
                                                        std::ostream & fout);
        void checkExchangeLocalEssentialTreeImpl(TagForceShort,
                                                 const DomainInfo & dinfo,
                                                 const F64 tolerance,
                                                 std::ostream & fout);
        template<class Tep2>
        void checkExchangeLocalEssentialTreeForShortImpl
        (TagSearchShortScatter,
         const DomainInfo & dinfo,
         const Tep2 ep_tmp[],
         const S32 jp_head,
         const S32 jp_tail,
         const S32 rank_target,
         const S32 n_image_per_proc,
         const ReallocatableArray<F64vec> & shift_image_domain,
         S32 & n_recv_per_proc,
         ReallocatableArray<F64vec> & pos_direct);
        
        template<class Tep2>
        void checkExchangeLocalEssentialTreeForShortImpl
        (TagSearchShortGather,
         const DomainInfo & dinfo,
         const Tep2 ep_tmp[],
         const S32 jp_head,
         const S32 jp_tail,
         const S32 rank_target,
         const S32 n_image_per_proc,
         const ReallocatableArray<F64vec> & shift_image_domain,
         S32 & n_recv_per_proc,
         ReallocatableArray<F64vec> & pos_direct);

        template<class Tep2>
        void checkExchangeLocalEssentialTreeForShortImpl
        (TagSearchShortSymmetry,
         const DomainInfo & dinfo,
         const Tep2 ep_tmp[],
         const S32 jp_head,
         const S32 jp_tail,
         const S32 rank_target,
         const S32 n_image_per_proc,
         const ReallocatableArray<F64vec> & shift_image_domain,
         S32 & n_recv_per_proc,
         ReallocatableArray<F64vec> & pos_direct);

    public:
// new 
        void setPrefixOfProfile(const char * str){
            //profile.setPrefix(str);
        }

        // for neighbour search
        ReallocatableArray<Tepj> * epj_neighbor_;
        TimeProfile getTimeProfile() const {return time_profile_;}
        void clearTimeProfile(){time_profile_.clear();}
        CountT getNumberOfWalkLocal() const { return n_walk_local_; }
        CountT getNumberOfInteractionEPEPLocal() const { return n_interaction_ep_ep_local_; }
        CountT getNumberOfInteractionEPSPLocal() const { return n_interaction_ep_sp_local_; }
        CountT getNumberOfWalkGlobal() const { return Comm::getSum(n_walk_local_); }
        CountT getNumberOfInteractionEPEPGlobal() const { return Comm::getSum(n_interaction_ep_ep_local_); }
        CountT getNumberOfInteractionEPSPGlobal() const { return Comm::getSum(n_interaction_ep_sp_local_); }

        CountT getNumberOfLETEPSend1stLocal() const {return n_let_ep_send_1st_;}
        CountT getNumberOfLETEPRecv1stLocal() const {return n_let_ep_recv_1st_;}
        CountT getNumberOfLETSPSend1stLocal() const {return n_let_sp_send_1st_;}
        CountT getNumberOfLETSPRecv1stLocal() const {return n_let_sp_recv_1st_;}
        CountT getNumberOfLETEPSend2ndLocal() const {return n_let_ep_send_2nd_;}
        CountT getNumberOfLETEPRecv2ndLocal() const {return n_let_ep_recv_2nd_;}

        CountT getNumberOfLETEPSend1stGlobal() const {return Comm::getSum(n_let_ep_send_1st_);}
        CountT getNumberOfLETEPRecv1stGlobal() const {return Comm::getSum(n_let_ep_recv_1st_);}
        CountT getNumberOfLETSPSend1stGlobal() const {return Comm::getSum(n_let_sp_send_1st_);}
        CountT getNumberOfLETSPRecv1stGlobal() const {return Comm::getSum(n_let_sp_recv_1st_);}
        CountT getNumberOfLETEPSend2ndGlobal() const {return Comm::getSum(n_let_ep_send_2nd_);}
        CountT getNumberOfLETEPRecv2ndGlobal() const {return Comm::getSum(n_let_ep_recv_2nd_);}

        CountT getNumberOfCellOpenLocal() const {return n_cell_open_[0]; }
        CountT getNumberOfCellOpenGlobal() const {return Comm::getSum(n_cell_open_[0]); }
        CountT getNumberOfCellGlobal() const {return tc_glb_.size(); }

        CountT getNumberOfProcSendLET1stICommSP() const {return n_proc_send_exchange_LET_1st__icomm_sp_;}
        CountT getNumberOfProcRecvLET1stICommSP() const {return n_proc_recv_exchange_LET_1st__icomm_sp_;}
        CountT getNumberOfProcSendLET1stICommEP() const {return n_proc_send_exchange_LET_1st__icomm_ep_;}
        CountT getNumberOfProcRecvLET1stICommEP() const {return n_proc_recv_exchange_LET_1st__icomm_ep_;}

        void clearNumberOfInteraction(){
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        }
        void clearCounterAll(){
            n_let_ep_send_1st_ = n_let_ep_recv_1st_ = n_let_sp_send_1st_ = n_let_sp_recv_1st_  = n_let_ep_send_2nd_ = n_let_ep_recv_2nd_ =  0;
            clearNumberOfInteraction();
            clearTimeProfile();
            //time_profile_.clear();
        }

        TreeForForce() : is_initialized_(false){}
        ~TreeForForce(){
            delete [] n_ep_send_;
            delete [] n_ep_recv_;
            delete [] n_cell_open_;
            delete [] n_ep_send_disp_;
            delete [] n_ep_recv_disp_;
            const S32 n_thread = Comm::getNumberOfThread();
            for(S32 i=0; i<n_thread; i++) delete [] id_proc_send_[i];
            delete [] id_proc_send_;
            delete [] id_sp_send_buf_;
            delete [] id_ptcl_send_;
            delete [] ip_disp_;
            delete [] tp_scatter_;
            delete [] tc_recv_1st_;
            delete [] shift_image_box_;
            delete [] epj_send_buf_;
            delete [] epj_recv_1st_sorted_;
            delete [] epjr_send_buf_;
            delete [] epjr_send_buf_for_scatter_;
            delete [] epjr_recv_1st_sorted_;
            delete [] epj_neighbor_;
            delete [] spj_for_force_;
            delete [] shift_image_domain_;
            delete [] id_ep_send_buf_;
            delete [] ep_send_buf_for_scatter_;
            delete [] epj_for_force_;
            delete [] id_epj_for_force_;
            delete [] id_spj_for_force_;
            if( typeid(TSM) == typeid(SEARCH_MODE_LONG) ||
                typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) ||
                typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER) ||
                typeid(TSM) == typeid(SEARCH_MODE_LONG_SYMMETRY)){
                delete [] n_sp_send_;
                delete [] n_sp_recv_;
                delete [] n_sp_send_disp_;
                delete [] n_sp_recv_disp_;
                delete [] n_ep_sp_send_;
                delete [] n_ep_sp_recv_;
            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            delete [] req_send_;
            delete [] req_recv_;
#endif
            if( typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY)
                || typeid(TSM) == typeid(SEARCH_MODE_GATHER) ){
                delete [] n_epj_recv_1st_;
                delete [] n_epj_recv_disp_1st_;
                delete [] n_epj_recv_2nd_;
                delete [] n_epj_recv_disp_2nd_;
                delete [] id_proc_src_;
                delete [] id_proc_dest_;
                for(S32 i=0; i<n_thread; i++){
                    delete [] adr_tc_level_partition_recv_1st_[i];
                }
                delete [] adr_tc_level_partition_recv_1st_;
            }
        }
        
        size_t getMemSizeUsed()const;

        void setNInteractionEPEP(const S64 n_ep_ep){
            n_interaction_ep_ep_ = n_ep_ep;
        }
        void setNInteractionEPSP(const S64 n_ep_sp){
            n_interaction_ep_sp_ = n_ep_sp;
        }

        S64 getNInteractionEPEP() const { return n_interaction_ep_ep_;}
        S64 getNInteractionEPSP() const { return n_interaction_ep_sp_;}


        void initialize(const U64 n_glb_tot,
                        const F64 theta=0.7,
                        const U32 n_leaf_limit=8,
                        const U32 n_group_limit=64);

        void reallocMem();
        void freeMem();
        void clearSizeOfArray();

        template<class Tpsys>
        void setParticleLocalTree(const Tpsys & psys, const bool clear=true);
        void setRootCell(const DomainInfo & dinfo);
        void setRootCell(const F64 l, const F64vec & c=F64vec(0.0));
        template<class Ttree>  void copyRootCell(const Ttree & tree);
        void mortonSortLocalTreeOnly(const bool reuse=false);
        void linkCellLocalTreeOnly();
        void linkCellGlobalTreeOnly();
        void calcMomentLocalTreeOnly();
        void calcMomentGlobalTreeOnly();
        //void addMomentAsSp(); // for multiwalk (send index)
        void makeIPGroup();
        void exchangeLocalEssentialTree(const DomainInfo & dinfo);
        void setLocalEssentialTreeToGlobalTree();
        void mortonSortGlobalTreeOnly(const bool reuse=false);
        S32 getNumberOfIPG() const { return ipg_.size();}
        void makeInteractionList(const S32 adr_ipg, const bool clear=true);
        template<class Tfunc_ep_ep>
        void calcForceOnly(Tfunc_ep_ep pfunc_ep_ep,
                           const S32 adr_ipg,
                           const bool clear = true);
        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForceOnly(Tfunc_ep_ep pfunc_ep_ep,
                           Tfunc_ep_sp pfunc_ep_sp,
                           const S32 adr_ipg,
                           const bool clear = true);
        void copyForceOriginalOrder();
        template<class Tfunc_ep_ep>
        void calcForce(Tfunc_ep_ep pfunc_ep_ep,
                       const bool clear=true);

        template<class Tfunc_ep_ep>
        void calcForceWalkOnly(Tfunc_ep_ep pfunc_ep_ep,
                               const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalk(Tfunc_dispatch pfunc_dispatch,
                               Tfunc_retrieve pfunc_retrieve,
                               const S32 tag_max,
                               const S32 n_walk_limit,
                               const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                    Tfunc_retrieve pfunc_retrieve,
                                    const S32 tag_max,
                                    const S32 n_walk_limit,
                                    const bool flag_keep_list,
                                    const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndexNew(Tfunc_dispatch pfunc_dispatch,
                                       Tfunc_retrieve pfunc_retrieve,
                                       const S32 tag_max,
                                       const S32 n_walk_limit,
                                       const bool flag_keep_list,
                                       const bool clear=true);

        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForce(Tfunc_ep_ep pfunc_ep_ep,
                       Tfunc_ep_sp pfunc_ep_sp,
                       const bool clear=true);

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                   Tpsys & psys,
                                   const bool clear=true){
            calcForce(pfunc_ep_ep, clear);
            const F64 time_offset = GetWtime();
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.calc_force += GetWtime() - time_offset;
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                   Tfunc_ep_sp pfunc_ep_sp,
                                   Tpsys & psys,
                                   const bool clear=true){
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear);
            const F64 time_offset = GetWtime();
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.calc_force += GetWtime() - time_offset;
        }



        Tforce getForce(const S32 i) const { return force_org_[i]; }

        ///////////////////////
        /// CHECK FUNCTIONS ///
        void checkMortonSortLocalTreeOnly(std::ostream & fout = std::cout);
        void checkMortonSortGlobalTreeOnly(std::ostream & fout = std::cout);
        void checkMakeLocalTree(const F64 tolerance = 1e-6, std::ostream & fout = std::cout);
        void checkMakeGlobalTree(const F64 tolerance = 1e-6, std::ostream & fout = std::cout);
        void checkCalcMomentLocalTree(const F64 tolerance = 1e-5, std::ostream & fout = std::cout);
        void checkCalcMomentGlobalTree(const F64 tolerance = 1e-5, std::ostream & fout = std::cout);
        void checkExchangeLocalEssentialTree(const DomainInfo & dinfo, 
                                             const F64 tolerance = 1e-5, 
                                             std::ostream & fout = std::cout);
        void checkMakeIPGroup(const F64 tolerance = 1e-5, std::ostream & fout = std::cout);
        void checkMakeInteractionList(const DomainInfo & dinfo,
                                      const S32 adr_ipg = 0, 
                                      const S32 ith = 0, 
                                      const F64 tolerance = 1e-5, 
                                      std::ostream & fout = std::cout){
            checkMakeInteractionListImpl(TSM::search_type(),  dinfo, adr_ipg, ith, tolerance, fout);
        }
        template<class Tfunc_ep_ep, class Tfunc_compare>
        void checkForce(Tfunc_ep_ep pfunc_ep_ep,
                        Tfunc_compare func_compare,
                        const DomainInfo & dinfo,
                        std::ostream & fout=std::cout);

        ////////////
        // for neighbour search APIs
        template<class Tptcl>
        S32 getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchShortScatter, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchShortGather, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchShortSymmetry, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchLongScatter, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchLongSymmetry, const Tptcl & ptcl, Tepj * & epj);


        F64ort getOuterBoundaryOfLocalTree(){
            return getOuterBoundaryOfLocalTreeImpl(typename TSM::search_type());
        }
        F64ort getInnerBoundaryOfLocalTree(){
            return getInnerBoundaryOfLocalTreeImpl(typename TSM::search_type());
        }
        F64ort getOuterBoundaryOfLocalTreeImpl(TagSearchLongSymmetry);
        F64ort getInnerBoundaryOfLocalTreeImpl(TagSearchLongSymmetry);
        F64ort getOuterBoundaryOfLocalTreeImpl(TagSearchShortSymmetry);
        F64ort getInnerBoundaryOfLocalTreeImpl(TagSearchShortSymmetry);        


        ///////////////
        // UTILS
        Tepj * getEpjFromId(const S64 id, const Tepj * epj_tmp=NULL);

        ///////////////
        // DEBUG
        S32 getNumberOfEpjSorted() const {return epj_sorted_.size();}

        /////////////////////////////////
        // FOR REUSING INTERACTION LIST
        template<class Ttreecell>
        void addMomentAsSpImpl(TagForceLong, ReallocatableArray<Ttreecell> & );
        template<class Ttreecell>
        void addMomentAsSpImpl(TagForceShort, ReallocatableArray<Ttreecell> & );

      void exchangeLocalEssentialTreeReuseList(const DomainInfo & dinfo,const bool flag_reuse = false);
        void exchangeLocalEssentialTreeReuseListImpl(TagSearchShortSymmetry,
                                                     const DomainInfo & dinfo,
                                                     const bool flag_reuse=false);
        void exchangeLocalEssentialTreeReuseListImpl(TagSearchShortScatter,
                                                     const DomainInfo & dinfo,
                                                     const bool flag_reuse=false);
        void exchangeLocalEssentialTreeReuseListImpl(TagSearchShortGather,
                                                     const DomainInfo & dinfo,
                                                     const bool flag_reuse=false);
        void exchangeLocalEssentialTreeReuseListImpl(TagSearchLong,
                                                     const DomainInfo & dinfo,
                                                     const bool flag_reuse=false);
        void exchangeLocalEssentialTreeReuseListImpl(TagSearchLongCutoff,
                                                     const DomainInfo & dinfo,
                                                     const bool flag_reuse=false);
        void exchangeLocalEssentialTreeReuseListImpl(TagSearchLongScatter,
                                                     const DomainInfo & dinfo,
                                                     const bool flag_reuse=false);
        void exchangeLocalEssentialTreeReuseListImpl(TagSearchLongSymmetry,
                                                     const DomainInfo & dinfo,
                                                     const bool flag_reuse=false);
        void exchangeLocalEssentialTreeReuseListLong(const DomainInfo & dinfo,
                                                     const bool flag_reuse=false);

        
        CommTable<Tepj, Tspj> comm_table_;
        InteractionList interaction_list_;

        void makeInteractionListIndexLong();
        void makeInteractionListIndexShort();
        void makeInteractionListIndex(TagForceLong){
            makeInteractionListIndexLong();
        }
        void makeInteractionListIndex(TagForceShort){
            makeInteractionListIndexShort();
        }
        template<class Tfunc_ep_ep>
        void calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                             const bool clear=true);
        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                             Tfunc_ep_sp pfunc_ep_sp,
                             const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkImpl(TagForceLong,
                                             Tfunc_dispatch pfunc_dispatch,
                                             Tfunc_retrieve pfunc_retrieve,
                                             const S32 n_walk_limit,
                                             const bool clear);
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkImpl(TagForceShort,
                                             Tfunc_dispatch pfunc_dispatch,
                                             Tfunc_retrieve pfunc_retrieve,
                                             const S32 n_walk_limit,
                                             const bool clear);        
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                         Tfunc_retrieve pfunc_retrieve,
                                         const S32 n_walk_limit,
                                         const bool clear=true){
            F64 wtime_offset = GetWtime();
            calcForceNoWalkForMultiWalkImpl(typename TSM::force_type(),
                                            pfunc_dispatch,
                                            pfunc_retrieve,
                                            n_walk_limit,
                                            clear);
            time_profile_.calc_force += GetWtime() - wtime_offset;
        }


        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkNewImpl(TagForceLong,
                                                Tfunc_dispatch pfunc_dispatch,
                                                Tfunc_retrieve pfunc_retrieve,
                                                const S32 n_walk_limit,
                                                const bool clear);
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkNewImpl(TagForceShort,
                                                Tfunc_dispatch pfunc_dispatch,
                                                Tfunc_retrieve pfunc_retrieve,
                                                const S32 n_walk_limit,
                                                const bool clear);        
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkNew(Tfunc_dispatch pfunc_dispatch,
                                            Tfunc_retrieve pfunc_retrieve,
                                            const S32 n_walk_limit,
                                            const bool clear=true){
            F64 wtime_offset = GetWtime();
            calcForceNoWalkForMultiWalkNewImpl(typename TSM::force_type(),
                                               pfunc_dispatch,
                                               pfunc_retrieve,
                                               n_walk_limit,
                                               clear);
            time_profile_.calc_force += GetWtime() - wtime_offset;
        }
        
        S32 adr_tc_level_partition_loc_[TREE_LEVEL_LIMIT+2];
        S32 lev_max_loc_;
        S32 adr_tc_level_partition_glb_[TREE_LEVEL_LIMIT+2];
        S32 lev_max_glb_;
        
        ////////////////////////////
        /// HIGH LEVEL FUNCTIONS ///

        template<class Tfunc_ep_ep>
        void calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep, 
                                 DomainInfo & dinfo,
                                 const bool clear_force=true,
                                 const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            if(list_mode == MAKE_LIST){
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                exchangeLocalEssentialTree(dinfo);
                setLocalEssentialTreeToGlobalTree();
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                makeIPGroup();
                calcForce(pfunc_ep_ep, clear_force);
            }
            else if(list_mode == MAKE_LIST_FOR_REUSE){
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                exchangeLocalEssentialTreeReuseList(dinfo);
                setLocalEssentialTreeToGlobalTree2();
                //this->n_glb_tot_ = tp_glb_.size();
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                makeIPGroup();
                makeInteractionListIndexShort();
                calcForceNoWalk(pfunc_ep_ep, clear_force);
            }
            else if(list_mode == REUSE_LIST){
                mortonSortLocalTreeOnly(true);
                exchangeLocalEssentialTreeReuseList(dinfo, true);
#if 1
                setLocalEssentialTreeToGlobalTree2(true);
                mortonSortGlobalTreeOnly(true);
#else
                // under construction
                CopyPtclLocPtclRecvToPtclSortedGlobalTree(typename TSM::force_type(),
                                                          epj_recv_, tp_loc_,
                                                          epj_org_, tp_glb_);
#endif
                calcForceNoWalk(pfunc_ep_ep, clear_force);
            }
            else{
                PARTICLE_SIMULATOR_PRINT_ERROR("INVALID INTERACTION_LIST_MODE.");
                std::cerr<<"INTERACTION_LIST_MODE: "<<list_mode<<std::endl;
                Abort(-1);
            }
        }
        
        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep,
                          Tpsys & psys,
                          DomainInfo & dinfo,
                          const bool clear_force=true,
                          const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            setParticleLocalTree(psys, true);
            calcForceMakingTree(pfunc_ep_ep, dinfo, clear_force, list_mode);
        }
        
        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllWalkOnly(Tfunc_ep_ep pfunc_ep_ep, 
                                  Tpsys & psys,
                                  DomainInfo & dinfo,
                                  const bool clear_force=true){
            setParticleLocalTree(psys);
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
            calcMomentGlobalTreeOnly();
            makeIPGroup();
            calcForceWalkOnly(pfunc_ep_ep, clear_force);
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
                                   Tpsys & psys,
                                   DomainInfo & dinfo,
                                   const bool clear_force=true){
            setParticleLocalTree(psys);
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            checkMortonSortLocalTreeOnly(); // check morton sort
            linkCellLocalTreeOnly();
            checkMakeLocalTree(); // check link cell
            calcMomentLocalTreeOnly();
            checkCalcMomentLocalTree(); // check calc moment
            exchangeLocalEssentialTree(dinfo);
            checkExchangeLocalEssentialTree(dinfo); // check ex let
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            checkMortonSortGlobalTreeOnly(); // check morton sort 
            linkCellGlobalTreeOnly();
            checkMakeGlobalTree(); // check link cell
            calcMomentGlobalTreeOnly();
            checkCalcMomentGlobalTree(); // check calc moment 
            makeIPGroup();
            checkMakeIPGroup(); // check  make ipg
            calcForce(pfunc_ep_ep, clear_force);
        }
        
        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                      Tpsys & psys,
                                      DomainInfo & dinfo,
                                      const bool clear_force = true,
                                      const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            calcForceAll(pfunc_ep_ep, psys, dinfo, clear_force, list_mode); 
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }


        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllWalkOnlyAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, 
                                              Tpsys & psys,
                                              DomainInfo & dinfo,
                                              const bool clear_force = true){
            calcForceAllWalkOnly(pfunc_ep_ep, psys, dinfo, clear_force); 
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }


        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllAndWriteBackWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               const bool clear_force = true){
            calcForceAllWithCheck(pfunc_ep_ep, psys, dinfo, clear_force);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }

        
        //////////////////
        // FOR LONG FORCE
        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep, 
                                 Tfunc_ep_sp pfunc_ep_sp,  
                                 DomainInfo & dinfo,
                                 const bool clear_force=true,
                                 const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            if(list_mode == MAKE_LIST){
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                exchangeLocalEssentialTree(dinfo);
                setLocalEssentialTreeToGlobalTree();
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                makeIPGroup();
                calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
            }
            else if(list_mode == MAKE_LIST_FOR_REUSE){
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                spj_sorted_loc_.clearSize();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_loc_, 0, spj_sorted_loc_);
                exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo);
                setLocalEssentialTreeToGlobalTree2(false);
                //setLocalEssentialTreeToGlobalTreeImpl2(typename TSM::force_type(), false);
                //this->n_glb_tot_ = tp_glb_.size();
                mortonSortGlobalTreeOnly();
                for(S32 i=0; i<tc_glb_.size(); i++){
                    tc_glb_[i].clear();
                }
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                SetOuterBoxGlobalTreeForLongCutoffTop(typename TSM::search_type(),
                                                      tc_glb_.getPointer(),
                                                      epj_org_.getPointer(), //to get rcut
                                                      n_leaf_limit_,
                                                      length_*0.5, center_);
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_,
                                  spj_sorted_.size(), spj_sorted_); // new
                makeIPGroup();
                makeInteractionListIndexLong();
                calcForceNoWalk(pfunc_ep_ep, pfunc_ep_sp, clear_force);
            }
            else if(list_mode == REUSE_LIST){
                mortonSortLocalTreeOnly(true);
                calcMomentLocalTreeOnly();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_loc_,
                                  0, spj_sorted_loc_);
                exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo, true);
                setLocalEssentialTreeToGlobalTreeImpl2(typename TSM::force_type(), true);
                mortonSortGlobalTreeOnly(true);
                calcMomentGlobalTreeOnly();
                S32 offset = epi_org_.size() + epj_recv_.size() + spj_recv_.size();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_,
                                  offset, spj_sorted_);
                calcForceNoWalk(pfunc_ep_ep, pfunc_ep_sp, clear_force);
            }
            else{
                PARTICLE_SIMULATOR_PRINT_ERROR("INVALID INTERACTION_LIST_MODE.");
                std::cerr<<"INTERACTION_LIST_MODE: "<<list_mode<<std::endl;
                Abort(-1);
            }
        }
        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep,
                          Tfunc_ep_sp pfunc_ep_sp,
                          Tpsys & psys,
                          DomainInfo & dinfo,
                          const bool clear_force=true,
                          const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            setParticleLocalTree(psys, true);
            calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, dinfo, clear_force, list_mode);
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
                                   Tfunc_ep_sp pfunc_ep_sp,  
                                   Tpsys & psys,
                                   DomainInfo & dinfo,
                                   const bool clear_force,
                                   std::ostream & fout){

            fout<<"setParticleLocalTree"<<std::endl;
            setParticleLocalTree(psys);
            fout<<"setRootCell"<<std::endl;
            setRootCell(dinfo);
            fout<<"mortonSortLocalTreeOnly"<<std::endl;
            mortonSortLocalTreeOnly();

            checkMortonSortLocalTreeOnly(fout); // check morton sort
            fout<<"linkCellLocalTreeOnly"<<std::endl;
            linkCellLocalTreeOnly();
            checkMakeLocalTree(1e-6, fout); // check link cell
            fout<<"calcMomentLocalTreeOnly"<<std::endl;
            calcMomentLocalTreeOnly();
            checkCalcMomentLocalTree(1e-6, fout); // check calc moment
            fout<<"exchangeLocalEssentialTree"<<std::endl;
            exchangeLocalEssentialTree(dinfo);
            checkExchangeLocalEssentialTree(dinfo); // check ex let
            fout<<"setLocalEssentialTreeToGlobalTree"<<std::endl;
            setLocalEssentialTreeToGlobalTree();
            fout<<"mortonSortGlobalTreeOnly"<<std::endl;
            mortonSortGlobalTreeOnly();
            checkMortonSortGlobalTreeOnly(fout); // check morton sort 
            fout<<"linkCellGlobalTreeOnly"<<std::endl;
            linkCellGlobalTreeOnly();
            checkMakeGlobalTree(1e-6, fout); // check link cell
            fout<<"calcMomentGlobalTreeOnly"<<std::endl;
            calcMomentGlobalTreeOnly();
            checkCalcMomentGlobalTree(1e-6, fout); // check calc moment 
            fout<<"makeIPGroup"<<std::endl;
            makeIPGroup();
            checkMakeIPGroup(1e-6, fout); // check  make ipg
            fout<<"calcForce"<<std::endl;
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                      Tfunc_ep_sp pfunc_ep_sp,
                                      Tpsys & psys,
                                      DomainInfo & dinfo,
                                      const bool clear_force=true,
                                      const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            if(list_mode != REUSE_LIST){ clearSizeOfArray(); }
            calcForceAll(pfunc_ep_ep, pfunc_ep_sp, psys, dinfo, clear_force, list_mode);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }





        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBackWithCheck(Tfunc_ep_ep pfunc_ep_ep,
                                               Tfunc_ep_sp pfunc_ep_sp,  
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               const bool clear_force,
                                               std::ostream & fout=std::cout){
            calcForceAllWithCheck(pfunc_ep_ep, pfunc_ep_sp, psys, dinfo, clear_force, fout);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }

        ////////////////
        // multiwalk
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMakingTreeMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                         Tfunc_retrieve pfunc_retrieve,
                                         const S32 tag_max,
                                         DomainInfo & dinfo,
                                         const S32 n_walk_limit,
                                         const bool clear=true,
                                         const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            S32 ret = 0;
            if(list_mode == MAKE_LIST){
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                exchangeLocalEssentialTree(dinfo);
                setLocalEssentialTreeToGlobalTree();
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                makeIPGroup();
                ret = calcForceMultiWalk(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, clear);
            }
            else if(list_mode == MAKE_LIST_FOR_REUSE){
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_loc_,
                                  0, spj_sorted_loc_);
                exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo);
                setLocalEssentialTreeToGlobalTree2(false);
                //setLocalEssentialTreeToGlobalTreeImpl2(typename TSM::force_type(), false);
                //this->n_glb_tot_ = tp_glb_.size();
                mortonSortGlobalTreeOnly();
                for(S32 i=0; i<tc_glb_.size(); i++){
                    tc_glb_[i].clear();
                }
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                SetOuterBoxGlobalTreeForLongCutoffTop(typename TSM::search_type(),
                                                      tc_glb_.getPointer(),
                                                      epj_org_.getPointer(), //to get rcut
                                                      n_leaf_limit_,
                                                      length_*0.5, center_);
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_,
                                  spj_sorted_.size(), spj_sorted_);
                makeIPGroup();
                ret = calcForceMultiWalkIndexNew(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, true, clear);
            }
            else if(list_mode == REUSE_LIST){
                mortonSortLocalTreeOnly(true);
                calcMomentLocalTreeOnly();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_loc_,
                                  0, spj_sorted_loc_);
                exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo, true);
                
                setLocalEssentialTreeToGlobalTree2(true);
                mortonSortGlobalTreeOnly(true);
                calcMomentGlobalTreeOnly();
                S32 offset = epi_org_.size() + epj_recv_.size() + spj_recv_.size();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_,
                                  offset, spj_sorted_);
                calcForceNoWalkForMultiWalkNew(pfunc_dispatch, pfunc_retrieve, n_walk_limit, clear);
            }
            return ret;
        }
        
        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                  Tfunc_retrieve pfunc_retrieve,
                                  const S32 tag_max,
                                  Tpsys & psys,
                                  DomainInfo & dinfo,
                                  const S32 n_walk_limit,
                                  const bool clear=true,
                                  const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            S32 ret = 0;
            setParticleLocalTree(psys, true);
            ret = calcForceMakingTreeMultiWalk(pfunc_dispatch,
                                         pfunc_retrieve,
                                         tag_max,
                                         dinfo,
                                         n_walk_limit,
                                         clear,
                                         list_mode);
            return ret;
        }

        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllAndWriteBackMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                              Tfunc_retrieve pfunc_retrieve,
                                              const S32 tag_max,
                                              Tpsys & psys,
                                              DomainInfo & dinfo,
                                              const S32 n_walk_limit,
                                              const bool clear=true,
                                              const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            S32 ret = 0;
            ret = calcForceAllMultiWalk(pfunc_dispatch, pfunc_retrieve,
                                        tag_max, psys, dinfo, n_walk_limit, clear, list_mode);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            return ret;
        }

        ///////
        // new
        void setLocalEssentialTreeToGlobalTree2(const bool flag_reuse = false){
            F64 time_offset = GetWtime();
            setLocalEssentialTreeToGlobalTreeImpl2(typename TSM::force_type(), flag_reuse);
            this->n_glb_tot_ = tp_glb_.size();
            time_profile_.set_particle_global_tree += GetWtime() - time_offset;
        } 
        void setLocalEssentialTreeToGlobalTreeImpl2(TagForceLong,
                                                    const bool flag_reuse = false){
            SetLocalEssentialTreeToGlobalTreeImpl(epj_recv_, spj_recv_, tp_loc_,
                                                  epj_org_,  spj_org_, tp_glb_,
                                                  flag_reuse);
        }
        void setLocalEssentialTreeToGlobalTreeImpl2(TagForceShort,
                                                    const bool flag_reuse = false){
            SetLocalEssentialTreeToGlobalTreeImpl(epj_recv_, tp_loc_,
                                                  epj_org_,  tp_glb_,
                                                  flag_reuse);
        }

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMakingTreeMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                              Tfunc_retrieve pfunc_retrieve,
                                              const S32 tag_max,
                                              DomainInfo & dinfo,
                                              const S32 n_walk_limit,
                                              const bool clear=true,
                                              const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            S32 ret = 0;
            if(list_mode == MAKE_LIST){
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                exchangeLocalEssentialTree(dinfo);
                setLocalEssentialTreeToGlobalTree();
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                SetOuterBoxGlobalTreeForLongCutoffTop(typename TSM::search_type(),
                                                      tc_glb_.getPointer(),
                                                      epj_org_.getPointer(), //to get rcut
                                                      n_leaf_limit_,
                                                      length_*0.5, center_);
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_,
                                  spj_sorted_.size(), spj_sorted_);
                makeIPGroup();
                ret = calcForceMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, 
                                              tag_max, n_walk_limit, false, clear);
            }
            else if(list_mode == MAKE_LIST_FOR_REUSE){
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_loc_,
                                  0, spj_sorted_loc_);
                exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo);
                setLocalEssentialTreeToGlobalTree2(false);
                //this->n_glb_tot_ = tp_glb_.size();
                mortonSortGlobalTreeOnly();
                for(S32 i=0; i<tc_glb_.size(); i++){
                    tc_glb_[i].clear();
                }
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                SetOuterBoxGlobalTreeForLongCutoffTop(typename TSM::search_type(),
                                                      tc_glb_.getPointer(),
                                                      epj_org_.getPointer(), //to get rcut
                                                      n_leaf_limit_,
                                                      length_*0.5, center_);
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_,
                                  spj_sorted_.size(), spj_sorted_);
                makeIPGroup();
#if 1
                ret = calcForceMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, 
                                              tag_max, n_walk_limit, true, clear);
#else
                makeInteractionListIndex(typename TSM::force_type());
                calcForceNoWalkForMultiWalk(pfunc_dispatch, pfunc_retrieve, n_walk_limit, clear);
#endif
            }
            else if(list_mode == REUSE_LIST){
                F64 wtime_0 = GetWtime();
                mortonSortLocalTreeOnly(true);
                F64 wtime_1 = GetWtime() - wtime_0;

                calcMomentLocalTreeOnly();
                F64 wtime_2 = GetWtime() - wtime_0;

                F64 wtime_tmp = GetWtime();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_loc_,
                                  0, spj_sorted_loc_);
                time_profile_.add_moment_as_sp_local += GetWtime() - wtime_tmp;
                F64 wtime_3 = GetWtime() - wtime_0;


                exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo, true);
                F64 wtime_4 = GetWtime() - wtime_0;

                setLocalEssentialTreeToGlobalTree2(true);
                F64 wtime_5 = GetWtime() - wtime_0;

                mortonSortGlobalTreeOnly(true);
                F64 wtime_6 = GetWtime() - wtime_0;

                calcMomentGlobalTreeOnly();
                F64 wtime_7 = GetWtime() - wtime_0;

                wtime_tmp = GetWtime();
                S32 offset = epi_org_.size() + epj_recv_.size() + spj_recv_.size();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_,
                                  offset, spj_sorted_);
                time_profile_.add_moment_as_sp_global += GetWtime() - wtime_tmp;
                F64 wtime_8 = GetWtime() - wtime_0;

                wtime_tmp = GetWtime();
                calcForceNoWalkForMultiWalk(pfunc_dispatch, pfunc_retrieve, n_walk_limit, clear);
                wtime_tmp = GetWtime() - wtime_tmp;
                F64 wtime_9 = GetWtime() - wtime_0;

                std::cerr<<"wtime_1= "<<wtime_1
                         <<" wtime_2= "<<wtime_2
                         <<" wtime_3= "<<wtime_3
                         <<" wtime_4= "<<wtime_4
                         <<" wtime_5= "<<wtime_5
                         <<" wtime_6= "<<wtime_6
                         <<" wtime_7= "<<wtime_7
                         <<" wtime_8= "<<wtime_8
                         <<" wtime_9= "<<wtime_9
                         <<std::endl;
            }
            else{
                PARTICLE_SIMULATOR_PRINT_ERROR("INVALID INTERACTION_LIST_MODE.");
                std::cerr<<"INTERACTION_LIST_MODE: "<<list_mode<<std::endl;
                Abort(-1);
            }
            return ret;
        }
        
        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                       Tfunc_retrieve pfunc_retrieve,
                                       const S32 tag_max,
                                       Tpsys & psys,
                                       DomainInfo & dinfo,
                                       const S32 n_walk_limit,
                                       const bool clear=true,
                                       const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            S32 ret = 0;
            setParticleLocalTree(psys, true);
            ret = calcForceMakingTreeMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, tag_max, dinfo, n_walk_limit, clear, list_mode);
            return ret;
        }

        
        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllAndWriteBackMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                                   Tfunc_retrieve pfunc_retrieve,
                                                   const S32 tag_max,
                                                   Tpsys & psys,
                                                   DomainInfo & dinfo,
                                                   const S32 n_walk_limit,
                                                   const bool clear=true,
                                                   const INTERACTION_LIST_MODE list_mode = MAKE_LIST){
            S32 ret = 0;
            ret = calcForceAllMultiWalkIndex(pfunc_dispatch, pfunc_retrieve,
                                             tag_max, psys, dinfo, n_walk_limit, clear, list_mode);
            F64 wtime_0 = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.write_back += GetWtime() - wtime_0;
            return ret;
        }


        template<class Tfunc_ep_ep>
        void calcForceDirect(Tfunc_ep_ep pfunc_ep_ep,
                             Tforce force[],
                             const DomainInfo & dinfo,
                             const bool clear=true);

        template<class Tfunc_ep_ep>
        void calcForceDirectAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                         const DomainInfo & dinfo,
                                         const bool clear=true);

        void dump(std::ostream & fout){
            fout<<"n_loc_tot_="<<n_loc_tot_<<std::endl;
            fout<<"n_glb_tot_="<<n_glb_tot_<<std::endl;
            fout<<"length_="<<length_<<std::endl;
            fout<<"center_="<<center_<<std::endl;
            fout<<"pos_root_cell_="<<pos_root_cell_<<std::endl;
        }

        void exchangeLocalEssentialTreeUsingCommTable(){}

        void dumpMemSizeUsed(std::ostream & fout){
            S32 n_thread = Comm::getNumberOfThread();
            if (Comm::getRank() == 0) {
                fout<<"tp_buf_.getMemSize()= "<<tp_buf_.getMemSize()<<std::endl;
                fout<<"tp_loc_.getMemSize()= "<<tp_loc_.getMemSize()<<std::endl;
                fout<<"tp_glb_.getMemSize()= "<<tp_glb_.getMemSize()<<std::endl;
                fout<<"tc_loc_.getMemSize()= "<<tc_loc_.getMemSize()<<std::endl;
                fout<<"tc_glb_.getMemSize()= "<<tc_glb_.getMemSize()<<std::endl;
                fout<<"epi_sorted_.getMemSize()= "<<epi_sorted_.getMemSize()<<std::endl;
                fout<<"epi_org_.getMemSize()= "<<epi_org_.getMemSize()<<std::endl;            
                fout<<"epj_sorted_.getMemSize()= "<<epj_sorted_.getMemSize()<<std::endl;
                fout<<"epj_org_.getMemSize()= "<<epj_org_.getMemSize()<<std::endl;
                fout<<"spj_sorted_.getMemSize()= "<<spj_sorted_.getMemSize()<<std::endl;
                fout<<"spj_org_.getMemSize()= "<<spj_org_.getMemSize()<<std::endl;
                //fout<<"epj_sorted_loc_.getMemSize()= "<<epj_sorted_loc_.getMemSize()<<std::endl;
                //fout<<"spj_sorted_loc_.getMemSize()= "<<spj_sorted_loc_.getMemSize()<<std::endl;
                fout<<"ipg_.getMemSize()= "<<ipg_.getMemSize()<<std::endl;
                fout<<"epj_send_.getMemSize()= "<<epj_send_.getMemSize()<<std::endl;
                fout<<"epj_recv_.getMemSize()= "<<epj_recv_.getMemSize()<<std::endl;
                fout<<"spj_send_.getMemSize()= "<<spj_send_.getMemSize()<<std::endl;
                fout<<"spj_recv_.getMemSize()= "<<spj_recv_.getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"id_ep_send_buf_["<<i<<"].getMemSize()= "<<id_ep_send_buf_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"id_sp_send_buf_["<<i<<"].getMemSize()= "<<id_sp_send_buf_[i].getMemSize()<<std::endl;
                //fout<<"id_proc_send_[0].getMemSize()= "<<id_proc_send_[0].getMemSize()<<std::endl;
                //fout<<"rank_proc_send_[0].getMemSize()= "<<rank_proc_send_[0].getMemSize()<<std::endl;
                fout<<"force_org_.getMemSize()= "<<force_org_.getMemSize()<<std::endl;
                fout<<"force_sorted_.getMemSize()= "<<force_sorted_.getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epj_for_force_["<<i<<"].getMemSize()= "<<epj_for_force_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"spj_for_force_["<<i<<"].getMemSize()= "<<spj_for_force_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"id_epj_for_force_["<<i<<"].getMemSize()= "<<id_epj_for_force_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"id_spj_for_force_["<<i<<"].getMemSize()= "<<id_spj_for_force_[i].getMemSize()<<std::endl;


                for(S32 i=0; i<n_thread; i++) fout<<"ep_send_buf_for_scatter_["<<i<<"].getMemSize()= "<<ep_send_buf_for_scatter_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"shift_image_domain_["<<i<<"].getMemSize()= "<<shift_image_domain_[i].getMemSize()<<std::endl;
                fout<<"epjr_sorted_.getMemSie()= "<<epjr_sorted_.getMemSize()<<std::endl;
                fout<<"epjr_send_.getMemSie()= "<<epjr_send_.getMemSize()<<std::endl;
                fout<<"epjr_recv_.getMemSie()= "<<epjr_recv_.getMemSize()<<std::endl;
                fout<<"epjr_recv_1st_buf_.getMemSie()= "<<epjr_recv_1st_buf_.getMemSize()<<std::endl;
                fout<<"epjr_recv_2nd_buf_.getMemSie()= "<<epjr_recv_2nd_buf_.getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epjr_send_buf_["<<i<<"].getMemSize()= "<<epjr_send_buf_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epjr_send_buf_for_scatter_["<<i<<"].getMemSize()= "<<epjr_send_buf_for_scatter_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epjr_recv_1st_sorted_["<<i<<"].getMemSize()= "<<epjr_recv_1st_sorted_[i].getMemSize()<<std::endl;
                fout<<"epj_recv_1st_buf_= "<<epj_recv_1st_buf_.getMemSize()<<std::endl;
                fout<<"epj_recv_2nd_buf_= "<<epj_recv_2nd_buf_.getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epj_send_buf_["<<i<<"].getMemSize()= "<<epj_send_buf_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"id_ptcl_send_["<<i<<"].getMemSize()= "<<id_ptcl_send_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"shift_image_box_["<<i<<"].getMemSize()= "<<shift_image_box_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"ip_disp_["<<i<<"].getMemSize()= "<<ip_disp_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"tp_scatter_["<<i<<"].getMemSize()= "<<tp_scatter_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"tc_recv_1st_["<<i<<"].getMemSize()= "<<tc_recv_1st_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epj_recv_1st_sorted_["<<i<<"].getMemSize()= "<<epj_recv_1st_sorted_[i].getMemSize()<<std::endl;
                fout<<"tree_top_outer_pos_.getMemSize()= "<<tree_top_outer_pos_.getMemSize()<<std::endl;
                fout<<"tree_top_inner_pos_.getMemSize()= "<<tree_top_inner_pos_.getMemSize()<<std::endl;
                
                //fout<<"id_epj_recorder_for_force_[0].getMemSize()= "<<id_epj_recorder_for_force_[0].getMemSize()<<std::endl;
                //fout<<"id_spj_recorder_for_force_[0].getMemSize()= "<<id_spj_recorder_for_force_[0].getMemSize()<<std::endl;
                //fout<<"n_epi_recorder_for_force_[0].getMemSize()= "<<n_epi_recorder_for_force_[0].getMemSize()<<std::endl;
                //fout<<"n_disp_epi_recorder_for_force_[0].getMemSize()= "<<n_disp_epi_recorder_for_force_[0].getMemSize()<<std::endl;
                //fout<<"n_epj_recorder_for_force_[0].getMemSize()= "<<n_epj_recorder_for_force_[0].getMemSize()<<std::endl;
                //fout<<"n_spj_recorder_for_force_[0].getMemSize()= "<<n_spj_recorder_for_force_[0].getMemSize()<<std::endl;
                //fout<<"n_disp_epj_recorder_for_force_[0].getMemSize()= "<<n_disp_epj_recorder_for_force_[0].getMemSize()<<std::endl;
                //fout<<"n_disp_spj_recorder_for_force_[0].getMemSize()= "<<n_disp_spj_recorder_for_force_[0].getMemSize()<<std::endl;
                //fout<<"adr_epj_loc2glb_.getMemSize()= "<<adr_epj_loc2glb_.getMemSize()<<std::endl;
                //fout<<"adr_epj_buf2glb_.getMemSize()= "<<adr_epj_buf2glb_.getMemSize()<<std::endl;
                //fout<<"adr_spj_buf2glb_.getMemSize()= "<<adr_spj_buf2glb_.getMemSize()<<std::endl;;
                //fout<<"adr_epj_org2glb_.getMemSize()= "<<adr_epj_org2glb_.getMemSize()<<std::endl;
                //fout<<"adr_ptcl_of_tp_loc_.getMemSize()= "<<adr_ptcl_of_tp_loc_.getMemSize()<<std::endl;
                //fout<<"id_recorder_for_interaction_.getMemSize()= "<<id_recorder_for_interaction_.getMemSize()<<std::endl;

                size_t size_id_ep_send_buf_tot = 0;
                size_t size_id_sp_send_buf_tot = 0;
                size_t size_epj_for_force_tot = 0;
                size_t size_spj_for_force_tot = 0;
                size_t size_id_epj_for_force_tot = 0;
                size_t size_id_spj_for_force_tot = 0;
                
                size_t size_ep_send_buf_for_scatter_tot = 0;
                size_t size_shift_image_domain_tot = 0;
                size_t size_epjr_send_buf_tot = 0;
                size_t size_epjr_send_buf_for_scatter_tot = 0;
                size_t size_epjr_recv_1st_sorted_tot = 0;
                size_t size_epj_send_buf_tot = 0;
                size_t size_id_ptcl_send_tot = 0;
                size_t size_shift_image_box_tot = 0;
                size_t size_ip_disp_tot = 0;
                size_t size_tp_scatter_tot = 0;
                size_t size_tc_recv_1st_tot = 0;
                size_t size_epj_recv_1st_sorted_tot = 0;
                for(S32 i=0; i<n_thread; i++){
                    size_id_ep_send_buf_tot += id_ep_send_buf_[i].getMemSize();
                    size_id_sp_send_buf_tot += id_sp_send_buf_[i].getMemSize();
                    size_epj_for_force_tot += epj_for_force_[i].getMemSize();
                    size_spj_for_force_tot += spj_for_force_[i].getMemSize();
                    size_id_epj_for_force_tot += id_epj_for_force_[i].getMemSize();
                    size_id_spj_for_force_tot += id_spj_for_force_[i].getMemSize();

                    size_ep_send_buf_for_scatter_tot += ep_send_buf_for_scatter_[i].getMemSize();
                    size_shift_image_domain_tot += shift_image_domain_[i].getMemSize();
                    size_epjr_send_buf_tot += epjr_send_buf_[i].getMemSize();
                    size_epjr_send_buf_for_scatter_tot += epjr_send_buf_for_scatter_[i].getMemSize();
                    size_epjr_recv_1st_sorted_tot += epjr_recv_1st_sorted_[i].getMemSize();
                    size_epj_send_buf_tot += epj_send_buf_[i].getMemSize();
                    size_id_ptcl_send_tot += id_ptcl_send_[i].getMemSize();
                    size_shift_image_box_tot += shift_image_box_[i].getMemSize();
                    size_ip_disp_tot += ip_disp_[i].getMemSize();
                    size_tp_scatter_tot += tp_scatter_[i].getMemSize();
                    size_tc_recv_1st_tot += tc_recv_1st_[i].getMemSize();
                    size_epj_recv_1st_sorted_tot += epj_recv_1st_sorted_[i].getMemSize();
                }
                
                fout<<"sum= "<<
                    (double)
                    (
                     tp_buf_.getMemSize()
                     +tp_loc_.getMemSize()
                     +tp_glb_.getMemSize()
                     +tc_loc_.getMemSize()
                     +tc_glb_.getMemSize()
                     +epi_sorted_.getMemSize()+epi_org_.getMemSize()+epj_sorted_.getMemSize()
                     +epj_org_.getMemSize()+spj_sorted_.getMemSize()+spj_org_.getMemSize()
                     //+epj_sorted_loc_.getMemSize()
                     //+spj_sorted_loc_.getMemSize()
                     +ipg_.getMemSize()
                     +epj_send_.getMemSize()
                     +epj_recv_.getMemSize()
                     +spj_send_.getMemSize()
                     +spj_recv_.getMemSize()
                     //+id_ep_send_buf_[0].getMemSize()
                     //+id_sp_send_buf_[0].getMemSize()
                     +size_id_ep_send_buf_tot
                     +size_id_sp_send_buf_tot
                     +force_org_.getMemSize()
                     +force_sorted_.getMemSize()
                     +size_epj_for_force_tot
                     +size_spj_for_force_tot
                     +size_id_epj_for_force_tot
                     +size_id_spj_for_force_tot

                     +size_ep_send_buf_for_scatter_tot
                     +size_shift_image_domain_tot
                     +epjr_sorted_.getMemSize()
                     +epjr_send_.getMemSize()
                     +epjr_recv_.getMemSize()
                     +epjr_recv_1st_buf_.getMemSize()
                     +epjr_recv_2nd_buf_.getMemSize()
                     +size_epjr_send_buf_tot
                     +size_epjr_send_buf_for_scatter_tot
                     +size_epjr_recv_1st_sorted_tot
                     +epj_recv_1st_buf_.getMemSize()
                     +epj_recv_2nd_buf_.getMemSize()
                     +size_epj_send_buf_tot
                     +size_id_ptcl_send_tot
                     +size_shift_image_box_tot
                     +size_ip_disp_tot
                     +size_tp_scatter_tot
                     +size_tc_recv_1st_tot
                     +size_epj_recv_1st_sorted_tot
                     +tree_top_outer_pos_.getMemSize()
                     +tree_top_inner_pos_.getMemSize()
                     
                     //+epj_for_force_[0].getMemSize()
                     //+spj_for_force_[0].getMemSize()
                     //+id_epj_for_force_[0].getMemSize()
                     //+id_spj_for_force_[0].getMemSize()
                     //+id_epj_recorder_for_force_[0].getMemSize()
                     //+id_spj_recorder_for_force_[0].getMemSize()
                     //+n_epi_recorder_for_force_[0].getMemSize()
                     //+n_disp_epi_recorder_for_force_[0].getMemSize()
                     //+n_epj_recorder_for_force_[0].getMemSize()
                     //+n_spj_recorder_for_force_[0].getMemSize()
                     //+n_disp_epj_recorder_for_force_[0].getMemSize()
                     //+n_disp_spj_recorder_for_force_[0].getMemSize()
                     //+adr_epj_loc2glb_.getMemSize()
                     //+adr_epj_buf2glb_.getMemSize()
                     //+adr_spj_buf2glb_.getMemSize()
                     //+adr_epj_org2glb_.getMemSize()
                     //+adr_ptcl_of_tp_loc_.getMemSize()
                     //+id_recorder_for_interaction_.getMemSize()
                     ) / 1e9
                    <<" [GB]"
                    <<std::endl;
                //comm_table_.dumpMemSizeUsed(fout);
            }
        }

        
    };

    template<class Tforce, class Tepi, class Tepj, class Tmom=void, class Tsp=void>
    class TreeForForceLong{
    public:
        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> Normal;

        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> WithCutoff;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SCATTER,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> WithScatterSearch; // for P^3T

        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF_SCATTER,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> WithCutoffScatterSearch; // for P^3T
    };

    template<class Tforce, class Tepi, class Tepj>
    class TreeForForceLong<Tforce, Tepi, Tepj, void, void>{
    public:
        // for P^3T
        typedef TreeForForce
        <SEARCH_MODE_LONG_SCATTER,
         Tforce, Tepi, Tepj,
         MomentMonopoleScatter,
         MomentMonopoleScatter,
         SPJMonopoleScatter> MonopoleWithScatterSearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SCATTER,
         Tforce, Tepi, Tepj,
         MomentQuadrupoleScatter,
         MomentQuadrupoleScatter,
         SPJQuadrupoleScatter> QuadrupoleWithScatterSearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SYMMETRY,
         Tforce, Tepi, Tepj,
         MomentMonopoleInAndOut,
         MomentMonopoleInAndOut,
         SPJMonopoleInAndOut> MonopoleWithSymmetrySearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SYMMETRY,
         Tforce, Tepi, Tepj,
         MomentQuadrupoleInAndOut,
         MomentQuadrupoleInAndOut,
         SPJQuadrupoleInAndOut> QuadrupoleWithSymmetrySearch;

        // for P^3T + PM
        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF_SCATTER,
         Tforce, Tepi, Tepj,
         MomentMonopoleCutoffScatter,
         MomentMonopoleCutoffScatter,
         SPJMonopoleCutoffScatter> MonopoleWithCutoffScatterSearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentMonopole,
         MomentMonopole,
         SPJMonopole> Monopole;

        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF,
         Tforce, Tepi, Tepj,
         MomentMonopoleCutoff,
         MomentMonopoleCutoff,
         SPJMonopoleCutoff> MonopoleWithCutoff;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentQuadrupole,
         MomentQuadrupole,
         SPJQuadrupole> Quadrupole;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentMonopoleGeometricCenter,
         MomentMonopoleGeometricCenter,
         SPJMonopoleGeometricCenter> MonopoleGeometricCenter;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentDipoleGeometricCenter,
         MomentDipoleGeometricCenter,
         SPJDipoleGeometricCenter> DipoleGeometricCenter;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentQuadrupoleGeometricCenter,
         MomentQuadrupoleGeometricCenter,
         SPJQuadrupoleGeometricCenter> QuadrupoleGeometricCenter;
    };

    template<class Tforce, class Tepi, class Tepj>
    class TreeForForceShort{
    public:

        typedef TreeForForce 
        <SEARCH_MODE_SYMMETRY,
         Tforce, Tepi, Tepj,
         MomentSearchInAndOut,
         MomentSearchInAndOut,
         SuperParticleBase> Symmetry;

        typedef TreeForForce 
        <SEARCH_MODE_GATHER,
         Tforce, Tepi, Tepj,
         MomentSearchInAndOut,
         MomentSearchInOnly,
         SuperParticleBase> Gather;

        // send_tree: out
        // recv_tree: in
        // loc_tree: in
        // glb_tree: out
        typedef TreeForForce 
        <SEARCH_MODE_SCATTER,
         Tforce, Tepi, Tepj,
         MomentSearchInAndOut,
         MomentSearchInAndOut,
         SuperParticleBase> Scatter;
    };


}
#include"tree_for_force_impl.hpp"
#include"tree_for_force_check_impl.hpp"





