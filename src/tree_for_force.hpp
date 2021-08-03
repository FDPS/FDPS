#pragma once

//#if __cplusplus <= 199711L
//#include<map>
//#else
//#include<unordered_map>
//#endif

#include<unordered_map>
#include<bitset>

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
        class Tspj, // PS or USER def
        enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE = CALC_DISTANCE_TYPE_NORMAL
      >
    class TreeForForce{

#if defined(DEV_CODE)
    public:
#else
    private:
#endif

	DomainInfo * p_domain_info_;
	
        CommInfo comm_info_;
        MortonKey morton_key_;
        F64ort inner_boundary_of_local_tree_;
        F64 length_; // length of a side of the root cell
        F64vec center_; // new member (not used)
        F64ort pos_root_cell_;
        TimeProfile time_profile_;
        using MyMap = std::unordered_map<S64, Tepj*>;
        
        MyMap map_id_to_epj_;

        enum EXCHANGE_LET_MODE exchange_let_mode_ = EXCHANGE_LET_A2A;
        
        CountT n_interaction_ep_ep_local_, n_interaction_ep_sp_local_, n_walk_local_;
        CountT * n_cell_open_;
        
        bool is_initialized_;

        S64 n_interaction_ep_ep_;
        S64 n_interaction_ep_sp_;
        S32 ni_ave_;
        S32 nj_ave_;

        RadixSort<KeyT, 8> rs_;
        S32 n_loc_tot_; // # of all kinds of assigned particles in local proc
        S64 n_let_sp_;
        S64 n_glb_tot_; // n_loc_tot_ + LETs
        S32 n_leaf_limit_;
        S32 n_group_limit_;
        F64 theta_;

        using TreeCellLoc = TreeCell< Tmomloc, Geometry<typename TSM::tree_cell_loc_geometry_type> >;
        using TreeCellGlb = TreeCell< Tmomglb, Geometry<typename TSM::tree_cell_glb_geometry_type> >;
        
        ReallocatableArray< TreeParticle> tp_glb_; // not removed for neighbour search
        ReallocatableArray< TreeCellLoc > tc_loc_; // not removed for reusing method (calc moment LT)
        ReallocatableArray< TreeCellGlb > tc_glb_; // not removed for neighbour search

        
        ReallocatableArray< Tepj > epj_sorted_; // not removed for neighbour search
        ReallocatableArray< Tepj > epj_send_; // not removed (sometimes needed (no API))
        ReallocatableArray< Tspj > spj_send_; // not removed (sometimes needed (no API))

        ReallocatableArray< IPGroup< typename TSM::ipg_type > > ipg_;  // not removed for reusing method
        ReallocatableArray< Tforce > force_org_;
        ReallocatableArray<S32> adr_org_from_adr_sorted_loc_;
        ReallocatableArray<U32> adr_org_from_adr_sorted_glb_;
        
        ReallocatableArray< Tepi > epi_sorted_; // msortLT --- final
        ReallocatableArray< Tspj > spj_sorted_; //  --- final
        ReallocatableArray< Tepi > epi_org_; // setPtclLT ---
        ReallocatableArray< Tspj > spj_org_; // insted of it, use spj_recv
        ReallocatableArray< Tepj > epj_org_;
        ReallocatableArray< Tforce > force_sorted_; // -- final

        ReallocatableArray<Tepj> * epj_for_force_;
        ReallocatableArray<Tspj> * spj_for_force_;
        ReallocatableArray<S32> * adr_epj_for_force_;
        ReallocatableArray<S32> * adr_spj_for_force_;
        ReallocatableArray<S32> * adr_ipg_for_force_;

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
        
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Request * req_send_;
        MPI_Request * req_recv_;
        MPI_Status  * status_;
#endif
        void calcMomentLocalTreeOnlyImpl(TagCalcMomLongEpjLt);
        void calcMomentLocalTreeOnlyImpl(TagCalcMomShortEpjLt);
        void calcMomentLocalTreeOnlyImpl(TagCalcMomShortEpiLt);
// probably, calcMomentGlobalTreeOnlyImpl is classified depending on force_type.
        void calcMomentGlobalTreeOnlyImpl(TagForceLong);
        void calcMomentGlobalTreeOnlyImpl(TagForceShort);

        void makeIPGroupImpl(TagForceLong);
        void makeIPGroupImpl(TagForceShort);

        void makeInteractionListLongForZeroTheta(TagWithoutCutoff, const S32 adr_ipg);
        void makeInteractionListLongForZeroTheta(TagWithCutoff, const S32 adr_ipg);
        void makeInteractionListImpl(TagForceLong, const S32 adr_ipg, const bool clear);
        void makeInteractionListImpl(TagForceShort, const S32 adr_ipg, const bool clear);
        void makeInteractionListIndexImpl(TagSearchLong, const S32 adr_ipg, const bool clear);

        template<typename Tpsys>
        void writeBack(Tpsys & psys){
            const F64 time_offset = GetWtime();
PS_OMP_PARALLEL_FOR
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.write_back += GetWtime() - time_offset;
        }
        
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
        S32 calcForceMultiWalkPtclImpl(TagForceLong,
                                            Tfunc_dispatch pfunc_dispatch,
                                            Tfunc_retrieve pfunc_retrieve,
                                            const S32 tag_max,
                                            const S32 n_walk_limit,
                                            const bool flag_keep_list,
                                            const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkPtclImpl(TagForceShort,
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

        // for short
        void calcCenterAndLengthOfRootCell(TagSearchShortGather,
                                           const DomainInfo & dinfo){
            calcCenterAndLengthOfRootCellImpl(epi_org_.getPointer(), dinfo);
        }
        void calcCenterAndLengthOfRootCell(TagSearchShortScatter,
                                           const DomainInfo & dinfo){
            calcCenterAndLengthOfRootCellImpl(epj_org_.getPointer(), dinfo);
        }
        void calcCenterAndLengthOfRootCell(TagSearchShortSymmetry,
                                           const DomainInfo & dinfo){
            calcCenterAndLengthOfRootCellImpl(epi_org_.getPointer(), dinfo);
        }
        // for long
        void calcCenterAndLengthOfRootCell(TagSearchLong,
                                           const DomainInfo & dinfo){
            calcCenterAndLengthOfRootCellImpl(epj_org_.getPointer(), dinfo);
        }
        void calcCenterAndLengthOfRootCell(TagSearchLongCutoff,
                                               const DomainInfo & dinfo){
            calcCenterAndLengthOfRootCellImpl(epj_org_.getPointer(), dinfo);
        }
        void calcCenterAndLengthOfRootCell(TagSearchLongScatter,
                                           const DomainInfo & dinfo){
            calcCenterAndLengthOfRootCellImpl(epj_org_.getPointer(), dinfo);
        }
        void calcCenterAndLengthOfRootCell(TagSearchLongSymmetry,
                                           const DomainInfo & dinfo){
            calcCenterAndLengthOfRootCellImpl(epj_org_.getPointer(), dinfo);
        }
        void calcCenterAndLengthOfRootCell(TagSearchLongCutoffScatter,
                                           const DomainInfo & dinfo){
            calcCenterAndLengthOfRootCellImpl(epj_org_.getPointer(), dinfo);
        }
        // for compile
        template<class Tep2>
        void calcCenterAndLengthOfRootCellImpl(const Tep2 ep[], const DomainInfo & dinfo);

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

      //void setDomain_info(const DomainInfo & dinfo){p_domain_info_ = &dinfo;}

	
        void freeObjectFromMemoryPool(){
#ifdef PARTICLE_SIMULATOR_USE_MEMORY_POOL
            epi_sorted_.freeMem(1);
            spj_sorted_.freeMem(1);
            force_sorted_.freeMem(1);
            for(int i=0; i<Comm::getNumberOfThread(); i++){
                epj_for_force_[i].freeMem(1);
                spj_for_force_[i].freeMem(1);
            }
            epi_org_.freeMem(1);
            spj_org_.freeMem(1);
            epj_org_.freeMem(1);
            epj_send_.freeMem(1);
            spj_send_.freeMem(1);
            assert(MemoryPool::getSize() == 0);
#else
            for(int i=0; i<Comm::getNumberOfThread(); i++){
                spj_for_force_[i].freeMem(1);
                epj_for_force_[i].freeMem(1);
            }
            spj_send_.freeMem(1);
            epj_send_.freeMem(1);
            epi_org_.freeMem(1);
#endif
        }
        
    public:
        using TypeEpi = Tepi;
        using TypeEpj = Tepj;

        void setCommInfo(const CommInfo & c);
        F64ort getPosRootCell() const {
            return pos_root_cell_;
        }
        
        void setPrefixOfProfile(const char * str){
            //profile.setPrefix(str);
        }

        void setExchangeLETMode(enum EXCHANGE_LET_MODE elm){
            exchange_let_mode_ = elm;
        }
        
        // for neighbour search
        ReallocatableArray<Tepj> * epj_neighbor_;
        
        TimeProfile getTimeProfile() const {return time_profile_;}
        void clearTimeProfile(){time_profile_.clear();}
        CountT getNumberOfWalkLocal() const { return n_walk_local_; }
        CountT getNumberOfInteractionEPEPLocal() const { return n_interaction_ep_ep_local_; }
        CountT getNumberOfInteractionEPSPLocal() const { return n_interaction_ep_sp_local_; }
        //CountT getNumberOfWalkGlobal() const { return Comm::getSum(n_walk_local_); }
        CountT getNumberOfWalkGlobal() const { return comm_info_.getSum(n_walk_local_); }
        //CountT getNumberOfInteractionEPEPGlobal() const { return Comm::getSum(n_interaction_ep_ep_local_); }
        CountT getNumberOfInteractionEPEPGlobal() const { return comm_info_.getSum(n_interaction_ep_ep_local_); }
        //CountT getNumberOfInteractionEPSPGlobal() const { return Comm::getSum(n_interaction_ep_sp_local_); }
        CountT getNumberOfInteractionEPSPGlobal() const { return comm_info_.getSum(n_interaction_ep_sp_local_); }
        
        CountT getNumberOfCellOpenLocal() const {return n_cell_open_[0]; }
        //CountT getNumberOfCellOpenGlobal() const {return Comm::getSum(n_cell_open_[0]); }
        CountT getNumberOfCellOpenGlobal() const {return comm_info_.getSum(n_cell_open_[0]); }
        CountT getNumberOfCellGlobal() const {return tc_glb_.size(); }
        
        CountT getNumberOfExLetAllgatherLocal() const {return comm_table_.rank_recv_allgather_.size();}
        //CountT getNumberOfExLetAllgatherGlobal() const {return Comm::getSum(comm_table_.rank_recv_allgather_.size());}
        CountT getNumberOfExLetAllgatherGlobal() const {return comm_info_.getSum(comm_table_.rank_recv_allgather_.size());}

        CountT getNumberOfExLetEpSendLocal() const {return comm_table_.n_ep_send_tot_;}
        CountT getNumberOfExLetSpSendLocal() const {return comm_table_.n_sp_send_tot_;}
        CountT getNumberOfExLetEpRecvLocal() const {return comm_table_.n_ep_recv_tot_;}
        CountT getNumberOfExLetSpRecvLocal() const {return comm_table_.n_sp_recv_tot_;}
        //CountT getNumberOfExLetEpGlobal() const {return Comm::getSum(comm_table_.n_ep_send_tot_);}
        CountT getNumberOfExLetEpGlobal() const {return comm_info_.getSum(comm_table_.n_ep_send_tot_);}
        //CountT getNumberOfExLetSpGlobal() const {return Comm::getSum(comm_table_.n_sp_send_tot_);}
        CountT getNumberOfExLetSpGlobal() const {return comm_info_.getSum(comm_table_.n_sp_send_tot_);}
        
        void clearNumberOfInteraction(){
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        }
        //void clearCommTable(){ comm_table_.clearCounter(); }
        void clearCounterAll(){
            //clearCommTable();
            clearNumberOfInteraction();
            clearTimeProfile();
        }

        TreeForForce() : is_initialized_(false){}
        ~TreeForForce(){
            delete [] n_cell_open_;
            delete [] epjr_send_buf_;
            delete [] epjr_send_buf_for_scatter_;
            delete [] epjr_recv_1st_sorted_;
            delete [] epj_neighbor_;
            delete [] spj_for_force_;
            delete [] epj_for_force_;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            delete [] req_send_;
            delete [] req_recv_;
#endif
        }
        
        size_t getMemSizeUsed()const;
        size_t getUsedMemorySize()const; //same function as getMemSizeUsed()

        void setNInteractionEPEP(const S64 n_ep_ep){
            n_interaction_ep_ep_ = n_ep_ep;
        }
        void setNInteractionEPSP(const S64 n_ep_sp){
            n_interaction_ep_sp_ = n_ep_sp;
        }

        S64 getNInteractionEPEP() const { return n_interaction_ep_ep_; }
        S64 getNInteractionEPSP() const { return n_interaction_ep_sp_; }

        void initialize(const U64 n_glb_tot,
                        const F64 theta=FDPS_DFLT_VAL_THETA,
                        const U32 n_leaf_limit=FDPS_DFLT_VAL_N_LEAF_LIMIT,
                        const U32 n_group_limit=FDPS_DFLT_VAL_N_GROUP_LIMIT);

        void reallocMem();
        void freeMem();
        void clearSizeOfArray();
        
        //#if __cplusplus <= 201112L
#ifdef TEST_VARIADIC_TEMPLATE
    #if 0
        template<class Tail>
        void setParticleLocalTreeRecursive(bool & first_set, Tail & tail){
        }
        template<class Head, class ... Tail>
        void setParticleLocalTreeRecursive(bool & first_set, Head & head, Tail & ... tail){
            bool clear = false;
            if(first_set){
                clear = true;
                first_set = false;
            }
            setParticleLocalTreeImpl(head, clear);
            setParticleLocalTreeRecursive(first_set, tail ...);
        }
        template<class... Args>
        void setParticleLocalTree(Args & ... args){
            std::cerr<<"check"<<std::endl;
            bool first_set = true;
            setParticleLocalTreeRecursive(first_set, args ...);
        }
    #endif
        // tuple version
        template<int N, class Ttpl>
        void setParticleLocalTreeTplRecursive(Ttpl & sys_tpl){}
        template<int N, class Ttpl, class Head, class ... Tail>
        void setParticleLocalTreeTplRecursive(Ttpl & sys_tpl){
            bool clear = false;
            if(N == 0) clear = true;
            setParticleLocalTreeImpl(std::get<N>(sys_tpl), clear);
            setParticleLocalTreeTplRecursive<N+1, Ttpl, Tail...>(sys_tpl);
        }
        template<class... Args>
        void setParticleLocalTree(std::tuple<Args...> & sys_tpl){
            setParticleLocalTreeTplRecursive<0, std::tuple<Args...>, Args...>(sys_tpl);
        }
        /*
        template<class Tpsys>
        void setParticleLocalTreeImpl(const Tpsys & psys, const bool clear);
        template<class Tpsys>
        void setParticleLocalTree(const Tpsys & psys, const bool clear=true){
            setParticleLocalTreeImpl(psys, clear);
        };
        */
#endif
        template<class Tpsys>
        void setParticleLocalTreeImpl(const Tpsys & psys, const bool clear);
        template<class Tpsys>
        void setParticleLocalTree(const Tpsys & psys, const bool clear=true){
            const int n_loc = psys.getNumberOfParticleLocal();
            epi_sorted_.resizeNoInitialize(n_loc);
            setParticleLocalTreeImpl(psys, clear);
        }

        void setDomainInfo(const DomainInfo & dinfo){p_domain_info_ = const_cast<DomainInfo*>(&dinfo);}
	void setDomainInfo(DomainInfo & dinfo){p_domain_info_ = &dinfo;}
        void setRootCell(const DomainInfo & dinfo);
        void setRootCell(const F64 l, const F64vec & c=F64vec(0.0));
        template<class Ttree>  void copyRootCell(const Ttree & tree);
        void mortonSortLocalTreeOnly(const bool reuse=false);
        void linkCellLocalTreeOnly();
        void linkCellGlobalTreeOnly();
        void calcMomentLocalTreeOnly();
        void calcMomentGlobalTreeOnly();
        void makeIPGroup();
        void mortonSortGlobalTreeOnly(const bool reuse=false);
        CountT getNumberOfIPG() const { return (CountT)ipg_.size();}
        CountT getNumberOfIPGLocal() const {
            return getNumberOfIPG();
        }
        CountT getNumberOfIPGGlobal() const {
            CountT n_loc = getNumberOfIPGLocal();
            //CountT n_glb = Comm::getSum(n_loc);
            CountT n_glb = comm_info_.getSum(n_loc);
            return n_glb;
        }
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

        template<class Tfunc_ep_ep>
        void calcForce(Tfunc_ep_ep pfunc_ep_ep,
                       const bool clear=true);

        template<class Tfunc_ep_ep>
        void calcForceWalkOnly(Tfunc_ep_ep pfunc_ep_ep,
                               const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                    Tfunc_retrieve pfunc_retrieve,
                                    const S32 tag_max,
                                    const S32 n_walk_limit,
                                    const bool flag_keep_list,
                                    const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkPtcl(Tfunc_dispatch pfunc_dispatch,
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
                                   const bool clear=true,
                                   const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            calcForce(pfunc_ep_ep, clear);
#if 1
            writeBack(psys);
#else
            const F64 time_offset = GetWtime();
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.calc_force += GetWtime() - time_offset;
#endif
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                   Tfunc_ep_sp pfunc_ep_sp,
                                   Tpsys & psys,
                                   const bool clear=true,
                                   const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear);
#if 1
            writeBack(psys);
#else
            const F64 time_offset = GetWtime();
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.calc_force += GetWtime() - time_offset;
#endif
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
        S32 getNeighborListOneParticleImpl(TagNeighborSearchScatter, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagNeighborSearchGather, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagNeighborSearchSymmetry, const Tptcl & ptcl, Tepj * & epj);
#ifdef PARTICLE_SIMULATOR_CHECK_SEARCH_MODE
        // for test_all_search_mode
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagNeighborSearchNo, const Tptcl & ptcl, Tepj * & epj);
#endif
        template<class Tptcl>
        S32 getNeighborListOneParticle(const Tptcl & ptcl, const S32 n_ngb, Tepj * & epj);
        

        F64ort getOuterBoundaryOfLocalTree(){
            return getOuterBoundaryOfLocalTreeImpl(typename TSM::search_type());
        }
        F64ort getInnerBoundaryOfLocalTree(){
            return inner_boundary_of_local_tree_;
        }
        F64ort getOuterBoundaryOfLocalTreeImpl(TagSearchLongSymmetry);
        F64ort getOuterBoundaryOfLocalTreeImpl(TagSearchShortSymmetry);


        ///////////////
        // UTILS
        Tepj * getEpjFromId(const S64 id, const Tepj * epj_tmp=NULL);

        ///////////////
        // DEBUG
        S32 getNumberOfEpiSorted() const {return epi_sorted_.size();}
        S32 getNumberOfEpjSorted() const {return epj_sorted_.size();}
        S32 getNumberOfSpjSorted() const {return spj_sorted_.size();}

        /////////////////////////////////
        // FOR REUSING INTERACTION LIST
        void addMomentAsSp(){
            F64 wtime_offset = GetWtime();
            AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_, n_let_sp_, spj_sorted_);
            time_profile_.add_moment_as_sp_global += GetWtime() - wtime_offset;
        }

        /////////////////////////////////
        // EXCHANGE LET
        void exchangeLocalEssentialTree
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        // LONG
        void exchangeLocalEssentialTreeImpl
        (TagSearchLong,
         const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeImpl
        (TagSearchLongCutoff,
         const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeImpl
        (TagSearchLongScatter,
         const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeImpl
        (TagSearchLongSymmetry,
         const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeImpl
        (TagSearchShortScatter,
         const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeImpl
        (TagSearchShortSymmetry,
         const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeImpl
        (TagSearchShortGather,
         const DomainInfo & dinfo,
         const bool flag_reuse=false);

        void exchangeLocalEssentialTreeLong
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeLongP2P
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeLongCutoffP2P
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeLongP2P2
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);        
        void exchangeLocalEssentialTreeShortScatter
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeShortScatterP2P
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeShortSymmetry
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeShortSymmetryP2P
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeShortGather
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        void exchangeLocalEssentialTreeShortGatherP2P
        (const DomainInfo & dinfo,
         const bool flag_reuse=false);
        
        
        CommTable comm_table_;
        InteractionList interaction_list_;

        void makeInteractionListIndexLong(DomainInfo & dinfo);
        void makeInteractionListIndexShort(DomainInfo & dinfo);
        void makeInteractionListIndex(TagForceLong, DomainInfo & dinfo){
            makeInteractionListIndexLong(dinfo);
        }
        void makeInteractionListIndex(TagForceShort, DomainInfo & dinfo){
            makeInteractionListIndexShort(dinfo);
        }
        template<class Tfunc_ep_ep>
        void calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                             const bool clear=true);
        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                             Tfunc_ep_sp pfunc_ep_sp,
                             const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkIndexImpl(TagForceLong,
                                                  Tfunc_dispatch pfunc_dispatch,
                                                  Tfunc_retrieve pfunc_retrieve,
                                                  const S32 n_walk_limit,
                                                  const bool clear);
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkIndexImpl(TagForceShort,
                                                  Tfunc_dispatch pfunc_dispatch,
                                                  Tfunc_retrieve pfunc_retrieve,
                                                  const S32 n_walk_limit,
                                                  const bool clear);
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                              Tfunc_retrieve pfunc_retrieve,
                                              const S32 n_walk_limit,
                                              const bool clear=true){
            F64 wtime_offset = GetWtime();
            calcForceNoWalkForMultiWalkIndexImpl(typename TSM::force_type(),
                                                 pfunc_dispatch,
                                                 pfunc_retrieve,
                                                 n_walk_limit,
                                                 clear);
            time_profile_.calc_force += GetWtime() - wtime_offset;
        }


        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkPtclImpl(TagForceLong,
                                                 Tfunc_dispatch pfunc_dispatch,
                                                 Tfunc_retrieve pfunc_retrieve,
                                                 const S32 n_walk_limit,
                                                 const bool clear);
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkPtclImpl(TagForceShort,
                                                 Tfunc_dispatch pfunc_dispatch,
                                                 Tfunc_retrieve pfunc_retrieve,
                                                 const S32 n_walk_limit,
                                                 const bool clear);
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        void calcForceNoWalkForMultiWalkPtcl(Tfunc_dispatch pfunc_dispatch,
                                             Tfunc_retrieve pfunc_retrieve,
                                             const S32 n_walk_limit,
                                             const bool clear=true){
            F64 wtime_offset = GetWtime();
            calcForceNoWalkForMultiWalkPtclImpl(typename TSM::force_type(),
                                                pfunc_dispatch,
                                                pfunc_retrieve,
                                                n_walk_limit,
                                                clear);
            time_profile_.calc_force += GetWtime() - wtime_offset;
        }

        void copyForceOriginalOrder();
        
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
                                 const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                                 const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
	    setDomainInfo(dinfo);
            if(list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE){
                setRootCell(dinfo);
                
                mortonSortLocalTreeOnly();

                epi_org_.freeMem(1);

                linkCellLocalTreeOnly();

                calcMomentLocalTreeOnly();

                exchangeLocalEssentialTree(dinfo, false);

                setLocalEssentialTreeToGlobalTree();


		
                epj_send_.freeMem(1);

                mortonSortGlobalTreeOnly();

                linkCellGlobalTreeOnly();

                spj_org_.freeMem(1);
                
                calcMomentGlobalTreeOnly();


		
                makeIPGroup();

                if(list_mode == MAKE_LIST){
                    calcForce(pfunc_ep_ep, clear_force);
                }
                else{
                    makeInteractionListIndexShort(dinfo);
                    calcForceNoWalk(pfunc_ep_ep, clear_force);
                }
                
            }
            else if(list_mode == REUSE_LIST){
                mortonSortLocalTreeOnly(true);

                epi_org_.freeMem(1);
                
                exchangeLocalEssentialTree(dinfo, true);
                setLocalEssentialTreeToGlobalTree(true);

                epj_send_.freeMem(1);
                
                mortonSortGlobalTreeOnly(true);

                spj_org_.freeMem(1);
                
                n_walk_local_ += ipg_.size();
                calcForceNoWalk(pfunc_ep_ep, clear_force);
            }
            else{
                PARTICLE_SIMULATOR_PRINT_ERROR("INVALID INTERACTION_LIST_MODE.");
                std::cerr<<"INTERACTION_LIST_MODE: "<<list_mode<<std::endl;
                Abort(-1);
            }
            freeObjectFromMemoryPool();
#ifdef PS_DEBUG_TREE_FORCE
            MemoryPool::checkEmpty();
#endif
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep,
                          Tpsys & psys,
                          DomainInfo & dinfo,
                          const bool clear_force=true,
                          const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                          const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            setParticleLocalTree(psys, true);
            calcForceMakingTree(pfunc_ep_ep, dinfo, clear_force, list_mode, flag_serialize);
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                      Tpsys & psys,
                                      DomainInfo & dinfo,
                                      const bool clear_force = true,
                                      const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                                      const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            calcForceAll(pfunc_ep_ep, psys, dinfo, clear_force, list_mode, flag_serialize);
            writeBack(psys);
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllWalkOnlyAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, 
                                              Tpsys & psys,
                                              DomainInfo & dinfo,
                                              const bool clear_force = true){
            calcForceAllWalkOnly(pfunc_ep_ep, psys, dinfo, clear_force);
            writeBack(psys);
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllAndWriteBackWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               const bool clear_force = true){
            calcForceAllWithCheck(pfunc_ep_ep, psys, dinfo, clear_force);
#if 1
            writeBack(psys);
#else
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
#endif
        }

        //////////////////
        // FOR LONG FORCE
        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep,
                                 Tfunc_ep_sp pfunc_ep_sp,  
                                 DomainInfo & dinfo,
                                 const bool clear_force=true,
                                 const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                                 const bool flag_serialize=false){
	    DEBUG_PRINT_MAKING_TREE(comm_info_);
            if (flag_serialize==true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
	    setDomainInfo(dinfo);
            if(list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE){
		DEBUG_PRINT_MAKING_TREE(comm_info_);
                setRootCell(dinfo);
		
		DEBUG_PRINT_MAKING_TREE(comm_info_);
                mortonSortLocalTreeOnly();

                epi_org_.freeMem();
		
		DEBUG_PRINT_MAKING_TREE(comm_info_);
                linkCellLocalTreeOnly();

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                calcMomentLocalTreeOnly();

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                exchangeLocalEssentialTree(dinfo, false);

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                setLocalEssentialTreeToGlobalTree(false);

                epj_send_.freeMem();
                spj_send_.freeMem();
		DEBUG_PRINT_MAKING_TREE(comm_info_);
                mortonSortGlobalTreeOnly();

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                linkCellGlobalTreeOnly();

                epj_org_.freeMem();
                spj_org_.freeMem();

		/*
                if(comm_info_.getRank()==0){
		  for(auto i=0; i<comm_info_.getNumberOfProc(); i++){
		    std::cerr<<"i= "<<i
			     <<" n_ep_send[i]= "<<comm_table_.n_ep_send_[i]
			     <<" n_sp_send[i]= "<<comm_table_.n_sp_send_[i]
			     <<" n_ep_recv[i]= "<<comm_table_.n_ep_recv_[i]
			     <<" n_sp_recv[i]= "<<comm_table_.n_sp_recv_[i]
			     <<std::endl;
		  }
		}
		*/
		/*
		std::cerr<<"Comm::getRank()= "<<Comm::getRank()
			 <<" comm_table_.rank_recv_allgather_.size()= "
			 <<comm_table_.rank_recv_allgather_.size()
			 <<" comm_table_.n_sp_recv_tot_= "<<comm_table_.n_sp_recv_tot_
			 <<" spj_sorted_.size()= "<<spj_sorted_.size()<<std::endl;
		*/
#ifdef DEBUG_EXLET_P2P
                //if(Comm::getRank()==0){
                if(comm_info_.getRank()==0){
                    std::cerr<<"comm_table_.rank_recv_allgather_.size()= "<<comm_table_.rank_recv_allgather_.size()<<std::endl;
                    std::cerr<<"comm_table_.n_sp_recv_tot_= "<<comm_table_.n_sp_recv_tot_<<std::endl;
                    std::cerr<<"spj_sorted_.size()= "<<spj_sorted_.size()<<std::endl;
                }
                for(S32 i=0; i<tp_glb_.size(); i++){
                    if(GetMSB(adr_org_from_adr_sorted_glb_[i])){
                        U32 adr_org = ClearMSB(adr_org_from_adr_sorted_glb_[i]);
                        S32 adr_rank = adr_org - (spj_sorted_.size()-comm_table_.rank_recv_allgather_.size());
                        U32 adr_sorted = ClearMSB(tp_glb_[i].adr_ptcl_);
                        if(adr_rank >= 0){
                            S32 rank_allgather = comm_table_.rank_recv_allgather_[adr_rank];
                            //if(Comm::getRank()==0){
                            if(comm_info_.getRank()==0){
                                std::cerr<<" rank_allgather= "
                                         <<rank_allgather
                                         <<" dinfo.getPosDomain(pos_root_cell_, rank_allgather)= "<<dinfo.getPosDomain(pos_root_cell_, rank_allgather)
                                         <<" spj_sorted_[adr_sorted].getPos()= "<<spj_sorted_[adr_sorted].getPos()
                                         <<std::endl;
                            }
                            //assert(dinfo.getPosDomain(pos_root_cell_, rank_allgather).contained(spj_sorted_[adr_sorted].getPos()));
                            assert(comm_table_.pos_domain_allgather_[adr_rank].contained(spj_sorted_[adr_sorted].getPos()));
                        }
                    }
                }
                //Comm::barrier();
                comm_info_.barrier();
                //exit(1);
#endif

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                calcMomentGlobalTreeOnly();

                if(list_mode == MAKE_LIST){

		    DEBUG_PRINT_MAKING_TREE(comm_info_);
                    makeIPGroup();

		    DEBUG_PRINT_MAKING_TREE(comm_info_);
                    calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);

                }
                else{

		    DEBUG_PRINT_MAKING_TREE(comm_info_);
                    addMomentAsSp();

		    DEBUG_PRINT_MAKING_TREE(comm_info_);
                    makeIPGroup();

		    DEBUG_PRINT_MAKING_TREE(comm_info_);
                    makeInteractionListIndexLong(dinfo);

		    DEBUG_PRINT_MAKING_TREE(comm_info_);
                    calcForceNoWalk(pfunc_ep_ep, pfunc_ep_sp, clear_force);
                }
            }
            else if(list_mode == REUSE_LIST){
		DEBUG_PRINT_MAKING_TREE(comm_info_);
                mortonSortLocalTreeOnly(true);

                epi_org_.freeMem();

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                calcMomentLocalTreeOnly();

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                exchangeLocalEssentialTree(dinfo, true);

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                setLocalEssentialTreeToGlobalTree(true);

                epj_send_.freeMem();
                spj_send_.freeMem();

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                mortonSortGlobalTreeOnly(true);

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                calcMomentGlobalTreeOnly();

                spj_org_.freeMem();

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_, spj_sorted_.size(), spj_sorted_);
                n_walk_local_ += ipg_.size();

		DEBUG_PRINT_MAKING_TREE(comm_info_);
                calcForceNoWalk(pfunc_ep_ep, pfunc_ep_sp, clear_force);
            }
            else{
                PARTICLE_SIMULATOR_PRINT_ERROR("INVALID INTERACTION_LIST_MODE.");
                std::cerr<<"INTERACTION_LIST_MODE: "<<list_mode<<std::endl;
                Abort(-1);
            }
            freeObjectFromMemoryPool();
	    DEBUG_PRINT_MAKING_TREE(comm_info_);
#ifdef PS_DEBUG_TREE_FORCE
            MemoryPool::checkEmpty();
#endif
        }

        
#ifdef TEST_VARIADIC_TEMPLATE
    #if 0 // variadic template version
        template<class Tfunc0, class Tfunc1, class Tdinfo>
        void calcForceAllImpl3l(std::true_type,
                                bool & flag_first_set,
                                Tfunc0 & func0,
                                Tfunc1 & func1,
                                Tdinfo & dinfo,
                                bool  clear_force=true,
                                INTERACTION_LIST_MODE list_mode=MAKE_LIST){
            calcForceMakingTree(func0, func1, dinfo,
                                clear_force, list_mode);
        }
        template<class Tfunc0, class Tfunc1, class Head0, class Head1,
                 class ... Tail>
        void calcForceAllImpl3l(std::false_type,
                                bool & flag_first_set,
                                Tfunc0 & func0,
                                Tfunc1 & func1,
                                Head0 & head0,
                                Head1 & head1,
                                Tail & ... tail){
            bool flag_clear = false;
            if(flag_first_set){
                flag_clear = true;
                flag_first_set = false;
            }
            setParticleLocalTree(head0, flag_clear);
            calcForceAllImpl3l(typename IsDomainInfo<Head1>::value(), flag_first_set, func0, func1, head1, tail ...);
        }
        template<class Head0, class Head1, class ... Tail>
        void calcForceAllImpl2(std::true_type,
                               Head0 & head0,
                               Head1 & head1,
                               Tail & ... tail){
            // head0 is a for function and head1 is a particle_system class
            std::cerr<<"head1 is a particle_system class"<<std::endl;
            //setParticleLocalTree(head1, true);
            //calcForceAllImpl3s(head0, tail ...);
        }
        template<class Head0, class Head1, class Head2, class ... Tail>
        void calcForceAllImpl2(std::false_type,
                               Head0 & head0,
                               Head1 & head1,
                               Head2 & head2,
                               Tail & ... tail){
            // head0 and head1 are force function
            std::cerr<<"head1 is a force function"<<std::endl;
            bool flag_first_set = true;
            calcForceAllImpl3l(typename IsDomainInfo<Head2>::value(), flag_first_set, head0, head1, head2, tail ...);
        }
        
        template<class Head0, class Head1, class ... Tail>
        void calcForceAllImpl(Head0 & head0,
                              Head1 & head1,
                              Tail & ... tail){
            calcForceAllImpl2(typename IsParticleSystem<Head1>::value(),
                              head0, head1, tail...);
        }
        template<class ... Args>
        void calcForceAll(Args && ... args){
            calcForceAllImpl(args ...);
        }
    #endif // variadic template version
        // tuple version
        template<class Tfunc_ep_ep, class ... Args >
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep,
                          std::tuple<Args ...> sys_tpl,
                          DomainInfo & dinfo,
                          const bool clear_force=true,
                          const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE){
            setParticleLocalTree(sys_tpl);
            calcForceMakingTree(pfunc_ep_ep, dinfo, clear_force, list_mode);
        }
        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class ... Args >
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep,
                          Tfunc_ep_sp pfunc_ep_sp,
                          std::tuple<Args ...> sys_tpl,
                          DomainInfo & dinfo,
                          const bool clear_force=true,
                          const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE){
            setParticleLocalTree(sys_tpl);
            calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, dinfo, clear_force, list_mode);
        }
        template<size_t N, class Ttpl>
        void setForceToParticleSystemImpl(S32 n_tot,
                                          Ttpl & sys_tpl){ }
        template<size_t N, class Ttpl, class Head, class ... Tail>
        void setForceToParticleSystemImpl(S32 n_tot,
                                          Ttpl & sys_tpl){
            S32 n_loc = (std::get<N>(sys_tpl)).getNumberOfParticleLocal();
            for(S32 i=0; i<n_loc; i++){
                const S32 i_dst = i + n_tot;
                (std::get<N>(sys_tpl))[i].copyFromForce(force_org_[i_dst]);
            }
            n_tot += n_loc;
            setForceToParticleSystemImpl<N+1, Ttpl, Tail ...>(n_tot, sys_tpl);
        }
        template<class ... Args>
        void setForceToParticleSystem(std::tuple<Args ... > sys_tpl){
            S32 n_tot = 0;
            setForceToParticleSystemImpl<0, std::tuple<Args ... >, Args ...>(n_tot, sys_tpl);
        }
        template<class Tfunc_ep_ep, class ... Args>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                      std::tuple<Args ...> sys_tpl,
                                      DomainInfo & dinfo,
                                      const bool clear_force=true,
                                      const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE){
            if(list_mode != REUSE_LIST){ clearSizeOfArray(); }
            calcForceAll(pfunc_ep_ep, sys_tpl, dinfo, clear_force, list_mode);
            setForceToParticleSystem(sys_tpl);
        }
        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class ... Args>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                      Tfunc_ep_sp pfunc_ep_sp,
                                      std::tuple<Args ...> sys_tpl,
                                      DomainInfo & dinfo,
                                      const bool clear_force=true,
                                      const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE){
            if(list_mode != REUSE_LIST){ clearSizeOfArray();}
            calcForceAll(pfunc_ep_ep, pfunc_ep_sp, sys_tpl, dinfo, clear_force, list_mode);
            setForceToParticleSystem(sys_tpl);
        }
#endif
        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep,
                          Tfunc_ep_sp pfunc_ep_sp,
                          Tpsys & psys,
                          DomainInfo & dinfo,
                          const bool clear_force=true,
                          const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                          const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            setParticleLocalTree(psys, true);
            calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, dinfo, clear_force, list_mode, flag_serialize);
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
				      Tfunc_ep_sp pfunc_ep_sp,
				      Tpsys & psys,
				      DomainInfo & dinfo,
				      const bool clear_force=true,
				      const INTERACTION_LIST_MODE list_mode   = FDPS_DFLT_VAL_LIST_MODE,
				      const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            if(list_mode != REUSE_LIST){ clearSizeOfArray(); }
            calcForceAll(pfunc_ep_ep, pfunc_ep_sp, psys, dinfo, clear_force, list_mode, flag_serialize);
            writeBack(psys);
        }





        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBackWithCheck(Tfunc_ep_ep pfunc_ep_ep,
                                               Tfunc_ep_sp pfunc_ep_sp,  
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               const bool clear_force,
                                               std::ostream & fout=std::cout){
            calcForceAllWithCheck(pfunc_ep_ep, pfunc_ep_sp, psys, dinfo, clear_force, fout);
#if 1
            writeBack(psys);
#else
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
#endif
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
                                         const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE){
            S32 ret = 0;
	    setDomainInfo(dinfo);
            if(list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE){
                bool flag_keep_list = false;
                if(list_mode == MAKE_LIST_FOR_REUSE){
                    flag_keep_list = true;
                }
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();

                epi_org_.freeMem(1);
                
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                exchangeLocalEssentialTree(dinfo, false);
                setLocalEssentialTreeToGlobalTree(false);

                epj_send_.freeMem(1);
                spj_send_.freeMem(1);
                
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();

                spj_org_.freeMem(1);
                
                calcMomentGlobalTreeOnly();

                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_, spj_sorted_.size(), spj_sorted_);
                makeIPGroup();
                ret = calcForceMultiWalkPtcl(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, flag_keep_list, clear);
            }
            else if(list_mode == REUSE_LIST){
                mortonSortLocalTreeOnly(true);

                epi_org_.freeMem(1);
                
                calcMomentLocalTreeOnly();
                exchangeLocalEssentialTree(dinfo, true);
                setLocalEssentialTreeToGlobalTree(true);

                epj_send_.freeMem(1);
                spj_send_.freeMem(1);
                
                mortonSortGlobalTreeOnly(true);

                spj_org_.freeMem(1);
                
                calcMomentGlobalTreeOnly();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_, spj_sorted_.size(), spj_sorted_);
                n_walk_local_ += ipg_.size();
                calcForceNoWalkForMultiWalkPtcl(pfunc_dispatch, pfunc_retrieve, n_walk_limit, clear);
            }
            freeObjectFromMemoryPool();
#ifdef PS_DEBUG_TREE_FORCE
            MemoryPool::checkEmpty();
#endif
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
                                  const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE){
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
                                              const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                                              const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            S32 ret = 0;
            ret = calcForceAllMultiWalk(pfunc_dispatch, pfunc_retrieve,
                                        tag_max, psys, dinfo, n_walk_limit, clear, list_mode);
#if 1
            writeBack(psys);
#else
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
#endif
            return ret;
        }

        ///////
        // SET LET
        void setLocalEssentialTreeToGlobalTree(const bool flag_reuse = false){
            F64 time_offset = GetWtime();
            setLocalEssentialTreeToGlobalTreeImpl(typename TSM::force_type(), flag_reuse);
            this->n_glb_tot_ = tp_glb_.size();
            time_profile_.set_particle_global_tree += GetWtime() - time_offset;
        }
        void setLocalEssentialTreeToGlobalTreeImpl(TagForceLong,
                                                   const bool flag_reuse = false){
            SetLocalEssentialTreeToGlobalTreeLong(epj_org_, spj_org_, n_loc_tot_, tp_glb_, morton_key_, flag_reuse);
            n_let_sp_ = spj_org_.size();
        }
        void setLocalEssentialTreeToGlobalTreeImpl(TagForceShort,
                                                    const bool flag_reuse = false){
            SetLocalEssentialTreeToGlobalTreeShort(epj_org_, n_loc_tot_, tp_glb_, morton_key_, flag_reuse);
        }

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMakingTreeMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                              Tfunc_retrieve pfunc_retrieve,
                                              const S32 tag_max,
                                              DomainInfo & dinfo,
                                              const S32 n_walk_limit,
                                              const bool clear=true,
                                              const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE){
            S32 ret = 0;
	    setDomainInfo(dinfo);
            if(list_mode == MAKE_LIST || list_mode == MAKE_LIST_FOR_REUSE){
                bool flag_keep_list = false;
                if(list_mode == MAKE_LIST_FOR_REUSE){
                    flag_keep_list = true;
                }
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();


                
                epi_org_.freeMem(1);
                
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();

                //std::cerr<<"A) epj_sorted_.size()= "<<epj_sorted_.size()
                //         <<" epj_org_.size()= "<<epj_org_.size()
                //         <<std::endl;
                
                exchangeLocalEssentialTree(dinfo, false);

                //std::cerr<<"B) epj_sorted_.size()= "<<epj_sorted_.size()
                //         <<" epj_org_.size()= "<<epj_org_.size()
                //         <<std::endl;
                
                setLocalEssentialTreeToGlobalTree(false);

                //std::cerr<<"C) epj_sorted_.size()= "<<epj_sorted_.size()
                //         <<" epj_org_.size()= "<<epj_org_.size()
                //         <<std::endl;
                
                epj_send_.freeMem(1);
                spj_send_.freeMem(1);
                
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();

                //epj_org_.freeMem(1);
                spj_org_.freeMem(1);
                
                calcMomentGlobalTreeOnly();
                /*
                SetOuterBoxGlobalTreeForLongCutoffTop(typename TSM::search_type(),
                                                      tc_glb_.getPointer(),
                                                      epj_sorted_.getPointer(), //to get rcut
                                                      n_leaf_limit_,
                                                      length_*0.5, center_);
                */
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_,
                                  spj_sorted_.size(), spj_sorted_);
                makeIPGroup();

                //tc_loc_.freeMem(1);

                //std::cerr<<"epj_sorted_.size()= "<<epj_sorted_.size()
                //         <<" spj_sorted_.size()= "<<spj_sorted_.size()
                //         <<std::endl;
                
                ret = calcForceMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, flag_keep_list, clear);
            }                
            else if(list_mode == REUSE_LIST){
                mortonSortLocalTreeOnly(true);
                
                epi_org_.freeMem(1);
                
                calcMomentLocalTreeOnly();
                exchangeLocalEssentialTree(dinfo, true);
                setLocalEssentialTreeToGlobalTree(true);

                epj_send_.freeMem(1);
                spj_send_.freeMem(1);
                
                mortonSortGlobalTreeOnly(true);

                //epj_org_.freeMem(1);
                spj_org_.freeMem(1);
                
                calcMomentGlobalTreeOnly();
                AddMomentAsSpImpl(typename TSM::force_type(), tc_glb_, spj_sorted_.size(), spj_sorted_);
                n_walk_local_ += ipg_.size();
                calcForceNoWalkForMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, n_walk_limit, clear);
            }
            else{
                PARTICLE_SIMULATOR_PRINT_ERROR("INVALID INTERACTION_LIST_MODE.");
                std::cerr<<"INTERACTION_LIST_MODE: "<<list_mode<<std::endl;
                Abort(-1);
            }
            freeObjectFromMemoryPool();
#ifdef PS_DEBUG_TREE_FORCE
            MemoryPool::checkEmpty();
#endif
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
                                       const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE){
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
                                                   const INTERACTION_LIST_MODE list_mode = FDPS_DFLT_VAL_LIST_MODE,
                                                   const bool flag_serialize=false){
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            S32 ret = 0;
            ret = calcForceAllMultiWalkIndex(pfunc_dispatch, pfunc_retrieve,
                                             tag_max, psys, dinfo, n_walk_limit, clear, list_mode);

#if 1
            writeBack(psys);
#else
            F64 wtime_0 = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.write_back += GetWtime() - wtime_0;
#endif
            return ret;
        }


        template<class Tfunc_ep_ep, class Tsys>
        void calcForceDirectParallel(Tfunc_ep_ep pfunc_ep_ep,
                                     Tforce force[],
                                     const Tsys & system,
                                     const DomainInfo & dinfo,
                                     const bool clear=true);
        
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
            //if (Comm::getRank() == 0) {
            if (comm_info_.getRank() == 0) {
                fout<<"tp_glb_.getMemSize()= "<<tp_glb_.getMemSize()<<std::endl;
                fout<<"tc_loc_.getMemSize()= "<<tc_loc_.getMemSize()<<std::endl;
                fout<<"tc_glb_.getMemSize()= "<<tc_glb_.getMemSize()<<std::endl;
                fout<<"epi_sorted_.getMemSize()= "<<epi_sorted_.getMemSize()<<std::endl;
                fout<<"epi_org_.getMemSize()= "<<epi_org_.getMemSize()<<std::endl;            
                fout<<"epj_sorted_.getMemSize()= "<<epj_sorted_.getMemSize()<<std::endl;
                fout<<"epj_org_.getMemSize()= "<<epj_org_.getMemSize()<<std::endl;
                fout<<"spj_sorted_.getMemSize()= "<<spj_sorted_.getMemSize()<<std::endl;
                fout<<"spj_org_.getMemSize()= "<<spj_org_.getMemSize()<<std::endl;
                fout<<"ipg_.getMemSize()= "<<ipg_.getMemSize()<<std::endl;
                fout<<"epj_send_.getMemSize()= "<<epj_send_.getMemSize()<<std::endl;
                fout<<"spj_send_.getMemSize()= "<<spj_send_.getMemSize()<<std::endl;
                fout<<"force_org_.getMemSize()= "<<force_org_.getMemSize()<<std::endl;
                fout<<"force_sorted_.getMemSize()= "<<force_sorted_.getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epj_for_force_["<<i<<"].getMemSize()= "<<epj_for_force_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"spj_for_force_["<<i<<"].getMemSize()= "<<spj_for_force_[i].getMemSize()<<std::endl;
                fout<<"epjr_sorted_.getMemSie()= "<<epjr_sorted_.getMemSize()<<std::endl;
                fout<<"epjr_send_.getMemSie()= "<<epjr_send_.getMemSize()<<std::endl;
                fout<<"epjr_recv_.getMemSie()= "<<epjr_recv_.getMemSize()<<std::endl;
                fout<<"epjr_recv_1st_buf_.getMemSie()= "<<epjr_recv_1st_buf_.getMemSize()<<std::endl;
                fout<<"epjr_recv_2nd_buf_.getMemSie()= "<<epjr_recv_2nd_buf_.getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epjr_send_buf_["<<i<<"].getMemSize()= "<<epjr_send_buf_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epjr_send_buf_for_scatter_["<<i<<"].getMemSize()= "<<epjr_send_buf_for_scatter_[i].getMemSize()<<std::endl;
                for(S32 i=0; i<n_thread; i++) fout<<"epjr_recv_1st_sorted_["<<i<<"].getMemSize()= "<<epjr_recv_1st_sorted_[i].getMemSize()<<std::endl;
                size_t size_epj_for_force_tot = 0;
                size_t size_spj_for_force_tot = 0;
                size_t size_epjr_send_buf_tot = 0;
                size_t size_epjr_send_buf_for_scatter_tot = 0;
                size_t size_epjr_recv_1st_sorted_tot = 0;
                for(S32 i=0; i<n_thread; i++){
                    size_epj_for_force_tot += epj_for_force_[i].getMemSize();
                    size_spj_for_force_tot += spj_for_force_[i].getMemSize();
                    size_epjr_send_buf_tot += epjr_send_buf_[i].getMemSize();
                    size_epjr_send_buf_for_scatter_tot += epjr_send_buf_for_scatter_[i].getMemSize();
                    size_epjr_recv_1st_sorted_tot += epjr_recv_1st_sorted_[i].getMemSize();
                }
                
                fout<<"sum= "<<
                    (double)
                    (
                     tp_glb_.getMemSize()
                     +tc_loc_.getMemSize()
                     +tc_glb_.getMemSize()
                     +epi_sorted_.getMemSize()+epi_org_.getMemSize()+epj_sorted_.getMemSize()
                     +epj_org_.getMemSize()+spj_sorted_.getMemSize()
                     +spj_org_.getMemSize()
                     +ipg_.getMemSize()
                     +epj_send_.getMemSize()
                     +spj_send_.getMemSize()
                     +force_org_.getMemSize()
                     +force_sorted_.getMemSize()
                     +size_epj_for_force_tot
                     +size_spj_for_force_tot
                     +epjr_sorted_.getMemSize()
                     +epjr_send_.getMemSize()
                     +epjr_recv_.getMemSize()
                     +epjr_recv_1st_buf_.getMemSize()
                     +epjr_recv_2nd_buf_.getMemSize()
                     +size_epjr_send_buf_tot
                     +size_epjr_send_buf_for_scatter_tot
                     +size_epjr_recv_1st_sorted_tot
                     ) / 1e9
                    <<" [GB]"
                    <<std::endl;
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
        <SEARCH_MODE_LONG_SYMMETRY,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> WithSymmetrySearch; // added by D.N. at 2018/11/28.
    };

    template<class Tforce, class Tepi, class Tepj>
    class TreeForForceLong<Tforce, Tepi, Tepj, void, void>{
    public:
        // for P^3T
        typedef TreeForForce
        <SEARCH_MODE_LONG_SCATTER,
         Tforce, Tepi, Tepj,
         MomentMonopole,
         MomentMonopole,
         SPJMonopole> MonopoleWithScatterSearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SCATTER,
         Tforce, Tepi, Tepj,
         MomentQuadrupole,
         MomentQuadrupole,
         SPJQuadrupole> QuadrupoleWithScatterSearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SYMMETRY,
         Tforce, Tepi, Tepj,
         MomentMonopole,
         MomentMonopole,
         SPJMonopole> MonopoleWithSymmetrySearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SYMMETRY,
         Tforce, Tepi, Tepj,
         MomentQuadrupole,
         MomentQuadrupole,
         SPJQuadrupole> QuadrupoleWithSymmetrySearch;

        // for P^3T + PM
        //typedef TreeForForce
        //<SEARCH_MODE_LONG_CUTOFF_SCATTER,
        // Tforce, Tepi, Tepj,
        // MomentMonopoleCutoffScatter,
        // MomentMonopoleCutoffScatter,
        // SPJMonopoleCutoffScatter> MonopoleWithCutoffScatterSearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentMonopole,
         MomentMonopole,
         SPJMonopole> Monopole;

        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF,
         Tforce, Tepi, Tepj,
         MomentMonopole,
         MomentMonopole,
         SPJMonopole> MonopoleWithCutoff;

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
         MomentShort,
         MomentShort,
         SuperParticleBase> Symmetry;


        typedef TreeForForce 
        <SEARCH_MODE_GATHER,
         Tforce, Tepi, Tepj,
         MomentShort,
         MomentShort,
         SuperParticleBase> Gather;

        typedef TreeForForce 
        <SEARCH_MODE_SCATTER,
         Tforce, Tepi, Tepj,
         MomentShort,
         MomentShort,
         SuperParticleBase> Scatter;
    };


}
#include"tree_for_force_impl.hpp"
//#include"tree_for_force_check_impl.hpp"





