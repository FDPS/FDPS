

#pragma once

//#define UNORDERED_SET

#ifdef UNORDERED_SET 
#include<unordered_set>
#else
#include<set>
#endif

#include<sort.hpp>
#include<tree.hpp>

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
        //F64 length_; // length of a side of the root cell
        //F64vec center_; // new member (not used)
        //F64ort pos_root_cell_;

    private:
        F64 Tcomm_tmp_;
        F64 Tcomm_scatterEP_tmp_;
        F64 TexLET0_, TexLET1_, TexLET2_, TexLET3_, TexLET4_;
        S64 n_ep_send_1st_, n_ep_recv_1st_, n_ep_send_2nd_, n_ep_recv_2nd_;
        F64 wtime_walk_LET_1st_, wtime_walk_LET_2nd_;

        bool is_initialized_;

        S64 n_interaction_;
        S32 ni_ave_;
        S32 nj_ave_;

        RadixSort<U64, 8> rs_;
        S32 n_loc_tot_; // # of all kinds of particles in local process
        S32 n_glb_tot_; // # of all kinds of particles in all processes
        S32 n_leaf_limit_;
        S32 n_group_limit_;

        S32 adr_tc_level_partition_[TREE_LEVEL_LIMIT+2];
        S32 lev_max_;
        F64 theta_;

        F64 length_; // length of a side of the root cell
        F64vec center_; // new member (not used)
        F64ort pos_root_cell_;

        ReallocatableArray< TreeParticle > tp_buf_, tp_loc_, tp_glb_;
        ReallocatableArray< TreeCell< Tmomloc > > tc_loc_;
        ReallocatableArray< TreeCell< Tmomglb > > tc_glb_;
        ReallocatableArray< Tepi > epi_sorted_, epi_org_;
        ReallocatableArray< Tepj > epj_sorted_, epj_org_;
        ReallocatableArray< Tspj > spj_sorted_, spj_org_;
        ReallocatableArray< IPGroup<TSM> > ipg_;
	
        S32 * n_ep_send_; // * n_proc
        S32 * n_sp_send_; // * n_proc
        S32 * n_ep_send_disp_; // * n_proc+1
        S32 * n_sp_send_disp_; // * n_proc+1
        S32 * n_ep_recv_; // * n_proc
        S32 * n_sp_recv_; // * n_proc
        S32 * n_ep_recv_disp_; // * n_proc+1
        S32 * n_sp_recv_disp_; // * n_proc+1

        ReallocatableArray<Tepj> epj_send_;
        ReallocatableArray<Tepj> epj_recv_;
        ReallocatableArray<Tspj> spj_send_;
        ReallocatableArray<Tspj> spj_recv_;

        ReallocatableArray<S32> * id_ep_send_buf_;
        ReallocatableArray<S32> * id_sp_send_buf_;

        S32 ** id_proc_send_; // id_proc_send_[n_thread][n_proc]

        ReallocatableArray<Tforce> force_sorted_;
        ReallocatableArray<Tforce> force_org_;

        ReallocatableArray<Tepj> * epj_for_force_;
        ReallocatableArray<Tspj> * spj_for_force_;

        S32 n_surface_for_comm_;

        template<class Tep2, class Tep3>
        inline void scatterEP(S32 n_send[],
                              S32 n_send_disp[],
                              S32 n_recv[],
                              S32 n_recv_disp[],
                              ReallocatableArray<Tep2> & ep_send,  // send buffer
                              ReallocatableArray<Tep2> & ep_recv,  // recv buffer
                              const ReallocatableArray<Tep3> & ep_org, // original
                              const DomainInfo & dinfo);
	
        void calcMomentLocalTreeOnlyImpl(TagSearchLong);
        void calcMomentLocalTreeOnlyImpl(TagSearchLongCutoff);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortScatter);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortGather);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortSymmetry);

        void exchangeLocalEssentialTreeImpl(TagSearchLong, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchLongCutoff, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortScatter, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortGather, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeGatherImpl(TagRSearch, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeGatherImpl(TagNoRSearch, const DomainInfo & dinfo);

        void setLocalEssentialTreeToGlobalTreeImpl(TagForceShort);
        void setLocalEssentialTreeToGlobalTreeImpl(TagForceLong);

        void calcMomentGlobalTreeOnlyImpl(TagSearchLong);
        void calcMomentGlobalTreeOnlyImpl(TagSearchLongCutoff);
        void calcMomentGlobalTreeOnlyImpl(TagSearchShortScatter);
        void calcMomentGlobalTreeOnlyImpl(TagSearchShortGather);
        void calcMomentGlobalTreeOnlyImpl(TagSearchShortSymmetry);

        void makeIPGroupImpl(TagForceLong);
        void makeIPGroupImpl(TagForceShort);

        void makeInteractionListImpl(TagSearchLong, const S32 adr_ipg);
        void makeInteractionListImpl(TagSearchLongCutoff, const S32 adr_ipg);
        void makeInteractionListImpl(TagSearchShortScatter, const S32 adr_ipg);
        void makeInteractionListImpl(TagSearchShortGather, const S32 adr_ipg);
        void makeInteractionListImpl(TagSearchShortSymmetry, const S32 adr_ipg); 

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
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLong){}

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

        TreeForForce() : is_initialized_(false){}

        size_t getMemSizeUsed()const;
	
        void initialize(const U64 n_glb_tot,
                        const F64 theta=0.7,
                        const U32 n_leaf_limit=8,
                        const U32 n_group_limit=64);

        template<class Tpsys>
        void setParticleLocalTree(const Tpsys & psys, const bool clear=true);
        void setRootCell(const DomainInfo & dinfo);
        void setRootCell(const F64 l, const F64vec & c=F64vec(0.0));
        template<class Ttree>  void copyRootCell(const Ttree & tree);
        void mortonSortLocalTreeOnly();
        void linkCellLocalTreeOnly();
        void linkCellGlobalTreeOnly();
        void calcMomentLocalTreeOnly();
        void calcMomentGlobalTreeOnly();
        void makeIPGroup();
        void exchangeLocalEssentialTree(const DomainInfo & dinfo);
        void setLocalEssentialTreeToGlobalTree();
        void mortonSortGlobalTreeOnly();
        S32 getNumberOfIPG() const { return ipg_.size();}
        void makeInteractionList(const S32 adr_ipg);
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
        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForce(Tfunc_ep_ep pfunc_ep_ep,
                       Tfunc_ep_sp pfunc_ep_sp,
                       const bool clear=true);
        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                   Tpsys & psys,
                                   const bool clear=true);
        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                   Tfunc_ep_sp pfunc_ep_sp,
                                   Tpsys & psys,
                                   const bool clear=true);
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


        //////////////////////////////
        /// MIDDLE LEVEL FUNCTIONS ///
        //////////////////////////////
        void makeLocalTree(DomainInfo & dinfo){
            setRootCell(dinfo);
	    mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
        }
        void makeLocalTree(const DomainInfo & dinfo){
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
        }
        void calcMomentGlobalTree(){
            calcMomentGlobalTreeOnly();
            makeIPGroup();
        }

        ////////////////////////////
        /// HIGH LEVEL FUNCTIONS ///
        //////////////////
        // FOR LONG FORCE
        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep, 
                          Tpsys & psys,
                          DomainInfo & dinfo,
                          const bool clear_force = true){
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
            calcForce(pfunc_ep_ep, clear_force);
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
                                   Tpsys & psys,
                                   DomainInfo & dinfo,
                                   const bool clear_force = true){

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
            //checkForce( pfunc_ep_ep, pfunc_compare_grav, dinfo); // check calc force
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, 
                                      Tpsys & psys,
                                      DomainInfo & dinfo,
                                      const bool clear_force = true){
            calcForceAll(pfunc_ep_ep, psys, dinfo, clear_force); 
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
        // FOR SHORT FORCE
        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep, 
                          Tfunc_ep_sp pfunc_ep_sp,  
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
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
        }

	template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
				   Tfunc_ep_sp pfunc_ep_sp,  
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
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
        }


        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllAndWriteBackWithTimer(Tfunc_ep_ep pfunc_ep_ep, 
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               Timer & tm,
                                               const bool clear_force=true){
            setParticleLocalTree(psys);
            tm.restart("setParticleLocalTree");
            setRootCell(dinfo);
            tm.restart("setRootCell");
            mortonSortLocalTreeOnly();
            tm.restart("mortonSortLocalTreeOnly");
            linkCellLocalTreeOnly();
            tm.restart("linkCellLocalTreeOnly");
            calcMomentLocalTreeOnly();
            tm.restart("calcMomentLocalTreeOnly");
            exchangeLocalEssentialTree(dinfo);
            tm.restart("exchangeLocalEssentialTree");
            setLocalEssentialTreeToGlobalTree();
            tm.restart("setLocalEssentialTreeToGlobalTree");
            mortonSortGlobalTreeOnly();
            tm.restart("mortonSortGlobalTreeOnly");
            linkCellGlobalTreeOnly();
            tm.restart("linkCellGlobalTreeOnly");
            calcMomentGlobalTreeOnly();
            tm.restart("calcMomentGlobalTreeOnly");
            makeIPGroup();
            tm.restart("makeIPGroup");
            calcForce(pfunc_ep_ep, clear_force);
            tm.restart("calcForce");
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            tm.restart("write back");
        }



// for debug
        template<class Tfunc_ep_ep, class Tpsys, class Ttree>
        void calcForceAllAndWriteBackWithTimer2(Tfunc_ep_ep pfunc_ep_ep, 
                                                Tpsys & psys,
                                                DomainInfo & dinfo,
                                                Timer & tm,
                                                const Ttree & tree,
                                                const bool clear_force=true){
            setParticleLocalTree(psys);
            tm.restart("setParticleLocalTree");
            //setRootCell(dinfo);
            copyRootCell(tree);
            tm.restart("setcopytree");
            mortonSortLocalTreeOnly();
            tm.restart("mortonSortLocalTreeOnly");
            linkCellLocalTreeOnly();
            tm.restart("linkCellLocalTreeOnly");
            calcMomentLocalTreeOnly();
            tm.restart("calcMomentLocalTreeOnly");
            exchangeLocalEssentialTree(dinfo);
            tm.restart("exchangeLocalEssentialTree");
            setLocalEssentialTreeToGlobalTree();
            tm.restart("setLocalEssentialTreeToGlobalTree");
            mortonSortGlobalTreeOnly();
            tm.restart("mortonSortGlobalTreeOnly");
            linkCellGlobalTreeOnly();
            tm.restart("linkCellGlobalTreeOnly");
            calcMomentGlobalTreeOnly();
            tm.restart("calcMomentGlobalTreeOnly");
            makeIPGroup();
            tm.restart("makeIPGroup");
            calcForce(pfunc_ep_ep, clear_force);
            tm.restart("calcForce");
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            tm.restart("write back");
        }



        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllWithTimer(Tfunc_ep_ep pfunc_ep_ep, 
                                   Tfunc_ep_sp pfunc_ep_sp,  
                                   Tpsys & psys,
                                   DomainInfo & dinfo,
                                   Timer & tm,
                                   const bool clear_force=true){
            setParticleLocalTree(psys);
	    tm.restart("setParticleLocalTree");
            setRootCell(dinfo);
	    tm.restart("setRootCell");
            mortonSortLocalTreeOnly();
	    tm.restart("mortonSortLocalTreeOnly");
            linkCellLocalTreeOnly();
	    tm.restart("linkCellLocalTreeOnly");
            calcMomentLocalTreeOnly();
	    tm.restart("calcMomentLocalTreeOnly");
            exchangeLocalEssentialTree(dinfo);
	    tm.restart("exchangeLocalEssentialTree");
            setLocalEssentialTreeToGlobalTree();
	    tm.restart("setLocalEssentialTreeToGlobalTree");
            mortonSortGlobalTreeOnly();
	    tm.restart("mortonSortGlobalTreeOnly");
            linkCellGlobalTreeOnly();
	    tm.restart("linkCellGlobalTreeOnly");
            calcMomentGlobalTreeOnly();
	    tm.restart("calcMomentGlobalTreeOnly");
            makeIPGroup();
	    tm.restart("makeIPGroup");
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
            tm.restart("calcForce");

        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, 
                                      Tfunc_ep_sp pfunc_ep_sp,  
                                      Tpsys & psys,
                                      DomainInfo & dinfo,
                                      const bool clear_force=true){
	    calcForceAll(pfunc_ep_ep, pfunc_ep_sp, psys, dinfo, clear_force);
	    for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBackWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
					       Tfunc_ep_sp pfunc_ep_sp,  
					       Tpsys & psys,
					       DomainInfo & dinfo,
					       const bool clear_force=true){
	    calcForceAll(pfunc_ep_ep, pfunc_ep_sp, psys, dinfo, clear_force);
	    for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }
	
        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBackWithTimer(Tfunc_ep_ep pfunc_ep_ep, 
                                               Tfunc_ep_sp pfunc_ep_sp,  
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               Timer & tm,
                                               const bool clear_force=true){
            setParticleLocalTree(psys);
	    tm.restart("setParticleLocalTree");
            setRootCell(dinfo);
	    tm.restart("setRootCell");
            mortonSortLocalTreeOnly();
	    tm.restart("mortonSortLocalTreeOnly");
            linkCellLocalTreeOnly();
	    tm.restart("linkCellLocalTreeOnly");
            calcMomentLocalTreeOnly();
	    tm.restart("calcMomentLocalTreeOnly");
            exchangeLocalEssentialTree(dinfo);
	    tm.restart("exchangeLocalEssentialTree");
            setLocalEssentialTreeToGlobalTree();
	    tm.restart("setLocalEssentialTreeToGlobalTree");
            mortonSortGlobalTreeOnly();
	    tm.restart("mortonSortGlobalTreeOnly");
            linkCellGlobalTreeOnly();
	    tm.restart("linkCellGlobalTreeOnly");
            calcMomentGlobalTreeOnly();
	    tm.restart("calcMomentGlobalTreeOnly");
            makeIPGroup();
            tm.restart("makeIPGroup");
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
            tm.restart("calcForce");
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);

            tm.restart("write back");
        }






        void dump_calc_cost(const double & tcal, std::ostream & fout){
            double speed_per_node = 128.0*1e9; // for K com
            double n_op_per_interaction = 28.0; // use the 3rd orde convergence with potential
            S64 n_interaction_tot = Comm::getSum(n_interaction_);
            fout<<"ni_ave_= "<<ni_ave_<<" nj_ave_= "<<nj_ave_<<" n_interaction_= "<<n_interaction_
                <<" speed: "<<(double)n_interaction_tot*n_op_per_interaction/tcal*1e-12<<" [Tflops] efficiency= "
                <<(double)n_interaction_tot*n_op_per_interaction/tcal/(speed_per_node*Comm::getNumberOfProc())<<" Tcomm_tmp_= "<<Tcomm_tmp_<<" Tcomm_scatterEP_tmp_= "<<Tcomm_scatterEP_tmp_<<std::endl;
            fout<<"TexLET0_= "<<TexLET0_<<"   TexLET1_= "<<TexLET1_<<"   TexLET2_= "<<TexLET2_<<"   TexLET3_= "<<TexLET3_<<"   TexLET4_= "<<TexLET4_<<std::endl;

            double max_wtime_walk_LET_1st;
            int rank_max_wtime_walk_LET_1st;
            Comm::getMaxValue(wtime_walk_LET_1st_, Comm::getRank(), max_wtime_walk_LET_1st, rank_max_wtime_walk_LET_1st);
            double max_wtime_walk_LET_2nd;
            int rank_max_wtime_walk_LET_2nd;
            Comm::getMaxValue(wtime_walk_LET_2nd_, Comm::getRank(), max_wtime_walk_LET_2nd, rank_max_wtime_walk_LET_2nd);
            fout<<"wtime_walk_LET_1st_= "<<wtime_walk_LET_1st_<<"   max_wtime_walk_LET_1st= "<<max_wtime_walk_LET_1st<<" (@"<<rank_max_wtime_walk_LET_1st<<")"
                <<"   wtime_walk_LET_2nd_= "<<wtime_walk_LET_2nd_<<"   max_wtime_walk_LET_2nd= "<<max_wtime_walk_LET_2nd<<" (@"<<rank_max_wtime_walk_LET_2nd<<")"<<std::endl;
            fout<<"n_ep_send_1st_= "<<n_ep_send_1st_<<"   n_ep_recv_1st_= "<<n_ep_recv_1st_
                <<"   n_ep_send_2nd_= "<<n_ep_send_2nd_<<"   n_ep_recv_2nd_= "<<n_ep_recv_2nd_<<std::endl;

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
    };

    template<class Tforce, class Tepi, class Tepj>
    class TreeForForceLong<Tforce, Tepi, Tepj, void, void>{
    public:
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





