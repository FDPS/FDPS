#pragma once

#include<particle_simulator.hpp>

extern int N_STEP;

namespace ParticleSimulator{


    ////////////
    // walk mode
    struct TagChopLeafTrue{};
    struct TagChopLeafFalse{};
    struct TagChopNonleafTrue{}; // PMMM
    struct TagChopNonleafFalse{}; // PMMM
    // walk mode
    ////////////

    ////////////
    // copy info close mode
    struct TagCopyInfoCloseNormal{};
    struct TagCopyInfoCloseNoSp{};
    struct TagCopyInfoCloseWithTpAdrptcl{};
    // copy info close mode
    ////////////
    
    ///////////
    // IS OPEN
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchLong,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){

        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen A ***************************");
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }

    // for consistency with PMM
    template<class TSM, class Ttc, class Tmorton_key> 
    inline bool IsOpen(TagSearchLong,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen B ***************************");
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }    
    
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongScatter,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen C ***************************");
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }

    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongSymmetry,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen D ***************************");
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }

    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongCutoff,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen E ***************************");
#if 1
        return ( IsOverlapped<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc].geo_.getVertexOut(), len_peri)
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#else
        return ( (target_box.vertex_.overlapped(tc_first[adr_tc].geo_.getVertexOut()))
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#endif
    }

    //////
    // PMMM
    template<class TSM, class Ttc, class Tmorton_key> 
    inline bool IsOpen(TagSearchLongParticleMeshMultipole,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen F ***************************");
        const S32 n_ptcl = tc_first[adr_tc].n_ptcl_;
        if (n_ptcl > 0) {
            const F64ort vertex_in = tc_first[adr_tc].geo_.getVertexIn();
            S32ort box;
            box.low_ = morton_key.getPMCellID(vertex_in.low_);
            box.high_ = morton_key.getPMCellID(vertex_in.high_);
            return (target_box.vertex_out_.overlapped(box));
        } else {
            return false;
        }
        // [Notes]
        // The above implementation is almost the same as the implementation
        // shown below, in which vertex_out_ is class F64ort. However,
        // the following implementation is not secure. 
        // return ( (target_box.vertex_out_.overlapped(tc_first[adr_tc].mom_.getVertexIn()))
        //          &&  (tc_first[adr_tc].n_ptcl_ > 0) );
    }

    // PMMM
    //////
    
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortScatter,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen G ***************************");
#if 1
        return ( IsOverlapped<CALC_DISTANCE_TYPE>(target_box.vertex_in_, tc_first[adr_tc].geo_.getVertexOut(), len_peri)
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#else
        return ( (target_box.vertex_in_.overlapped(tc_first[adr_tc].geo_.getVertexOut()))
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#endif
    }
    
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortSymmetry,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen H ***************************");
#if 1
        return ( (IsOverlapped<CALC_DISTANCE_TYPE>(target_box.vertex_in_, tc_first[adr_tc].geo_.getVertexOut(), len_peri)
                  || target_box.vertex_out_.overlapped(tc_first[adr_tc].geo_.getVertexIn()))
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#else
        return ( ( (target_box.vertex_in_.overlapped(tc_first[adr_tc].geo_.getVertexOut()))
                   || (target_box.vertex_out_.overlapped(tc_first[adr_tc].geo_.getVertexIn())) )
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#endif
    }
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortGather,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** IsOpen I ***************************");
#if 1
        return ( IsOverlapped<CALC_DISTANCE_TYPE>(target_box.vertex_out_, tc_first[adr_tc].geo_.getVertexIn(), len_peri)
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#else
        return ( target_box.vertex_out_.overlapped(tc_first[adr_tc].geo_.getVertexIn())
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#endif
    }    
    // IS OPEN
    ///////////
    
    ///////////
    // COPY INFO DISTANT

    template <class TSM, class Ttc, class Tmorton_key>
    inline void CopyInfoDistant2(const S32 adr_tc,
                                 const ReallocatableArray<Ttc> & tc_first,
                                 const S32 adr_sp,
                                 ReallocatableArray<S32> & adr_sp_list,
                                 const TargetBox<TSM> & target_box,
                                 const Tmorton_key & morton_key){

        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoDistant2 A ***************************");
        
        if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
            PS_COMPILE_TIME_MESSAGE("SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE@CopyInfoDistant2");
            if (target_box.contains(tc_first[adr_tc].mom_.getPos(), morton_key)) {
#ifdef NO_PUSH_BACK
                adr_sp_list.pushBackNoCheck(adr_sp);
#else
                adr_sp_list.push_back(adr_sp);
#endif
            }
        } else if constexpr ( std::is_same_v<typename TSM::force_type, TagForceLong>) {
	    PS_COMPILE_TIME_MESSAGE("TagForceLong@CopyInfoDistant2");
#ifdef NO_PUSH_BACK
            adr_sp_list.pushBackNoCheck(adr_sp);
#else
            adr_sp_list.push_back(adr_sp);
#endif            
        } else if constexpr ( std::is_same_v<typename TSM::force_type, TagForceShort>) {
  	    PS_COMPILE_TIME_MESSAGE("TagForceShort@CopyInfoDistant2");
            // do nothing
        } else {
            static_assert([]{return false;}());
        }
    }
    
    // PMMM
    template <class TSM, class Ttc, class Tmorton_key>
    inline void CopyInfoDistant(TagForceShort,
                                TagChopNonleafTrue,
                                const S32 adr_tc,
                                const ReallocatableArray<Ttc> & tc_first,
                                const S32 adr_sp,
                                ReallocatableArray<S32> & adr_sp_list,
                                const TargetBox<TSM> & target_box,
                                const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoDistant B ***************************");
        // do nothing
    }
    template <class TSM, class Ttc, class Tmorton_key>
    inline void CopyInfoDistant(TagForceShort,
                                TagChopNonleafFalse,
                                const S32 adr_tc,
                                const ReallocatableArray<Ttc> & tc_first,
                                const S32 adr_sp,
                                ReallocatableArray<S32> & adr_sp_list,
                                const TargetBox<TSM> & target_box,
                                const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoDistant C ***************************");
        // do nothing
    }
    template <class TSM, class Ttc, class Tmorton_key>
    inline void CopyInfoDistant(TagForceLong,
                                TagChopNonleafTrue,
                                const S32 adr_tc,
                                const ReallocatableArray<Ttc> & tc_first,
                                const S32 adr_sp,
                                ReallocatableArray<S32> & adr_sp_list,
                                const TargetBox<TSM> & target_box,
                                const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoDistant D ***************************");
        //DebugUtils::fouts(Comm::getThreadNum()) << "distant, begin" << std::endl;
        if (target_box.contains(tc_first[adr_tc].mom_.getPos(), morton_key)) {
#ifdef NO_PUSH_BACK
            adr_sp_list.pushBackNoCheck(adr_sp);
#else
            adr_sp_list.push_back(adr_sp);
#endif
        }
    }
    template <class TSM, class Ttc, class Tmorton_key>
    inline void CopyInfoDistant(TagForceLong,
                                TagChopNonleafFalse,
                                const S32 adr_tc,
                                const ReallocatableArray<Ttc> & tc_first,
                                const S32 adr_sp,
                                ReallocatableArray<S32> & adr_sp_list,
                                const TargetBox<TSM> & target_box,
                                const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoDistant E ***************************");
#ifdef NO_PUSH_BACK
      adr_sp_list.pushBackNoCheck(adr_sp);
#else
      adr_sp_list.push_back(adr_sp);
#endif
    }

    // original
    inline void CopyInfoDistant(TagForceShort,
                                const S32 adr_sp,
                                ReallocatableArray<S32> & adr_sp_list){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoDistant F ***************************");
        // do nothing
    }
    inline void CopyInfoDistant(TagForceLong,
                                const S32 adr_sp,
                                ReallocatableArray<S32> & adr_sp_list){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoDistant G ***************************");
        adr_sp_list.push_back(adr_sp);
    }

    // COPY INFO DISTANT
    ///////////

    ///////////
    // COPY INFO CLOSE
    template<class TSM, class Ttp, class Tep, class Tsp, class Tmorton_key, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClosePMM(TagChopLeafTrue,
                                 TagCopyInfoCloseNoSp,
                                 const ReallocatableArray<Ttp> & tp_first,
                                 const S32 adr_ptcl,
                                 const S32 n_ptcl,
                                 const ReallocatableArray<Tep> & ep_first,
                                 ReallocatableArray<S32> & adr_ep_list,
                                 const ReallocatableArray<Tsp> & sp_first,
                                 ReallocatableArray<S32> & adr_sp_list,
                                 const TargetBox<TSM> & target_box,
                                 const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClosePMM A ***************************");
        PS_COMPILE_TIME_MESSAGE("TagChopLeafTrue and TagCopyInfoCloseNoSp, @ CopyInfoClosePMM")
#if !defined(NO_PUSH_BACK)
            adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
#endif
        S32 cnt_adr_ptcl = adr_ptcl;        
        if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
                const KeyT key = tp_first[cnt_adr_ptcl].key_;
                if (target_box.contains(key, morton_key)) {
                    adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
                }
            }
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
                for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                    assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
                    adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
                }
            }
        }
    }

    
    template<class TSM, class Ttp, class Tep, class Tsp, class Tmorton_key, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClosePMM(TagChopLeafFalse,
                                 TagCopyInfoCloseNoSp,
                                 const ReallocatableArray<Ttp> & tp_first,
                                 const S32 adr_ptcl,
                                 const S32 n_ptcl,
                                 const ReallocatableArray<Tep> & ep_first,
                                 ReallocatableArray<S32> & adr_ep_list,
                                 const ReallocatableArray<Tsp> & sp_first,
                                 ReallocatableArray<S32> & adr_sp_list,
                                 const TargetBox<TSM> & target_box,
                                 const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClosePMM B ***************************");
        PS_COMPILE_TIME_MESSAGE("A) TagChopLeafFalse and TagCopyInfoCloseNoSp, @ CopyInfoClosePMM")
#if !defined(NO_PUSH_BACK)
        //if(N_STEP == 1431 && Comm::getRank() == 14){
        //    std::cout<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
        //    adr_ep_list.dump("check 1", std::cout);
        //}
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        //if(N_STEP == 1431 && Comm::getRank() == 14){
        //    std::cout<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
        //    adr_ep_list.dump("check 2", std::cout);
        //}
#endif
        S32 cnt_adr_ptcl = adr_ptcl;        
        if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                //if(N_STEP == 1431 && Comm::getRank() == 14){
                //    std::cout<<"ip= "<<ip<<" check 3"<<std::endl;
                //}
                assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
                const KeyT key = tp_first[cnt_adr_ptcl].key_;
                //adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);

                //if(N_STEP == 1431 && Comm::getRank() == 14){
                //    std::cout<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                //    adr_ep_list.dump("check 4", std::cout);
                //}
                
                adr_ep_list.push_back(cnt_adr_ptcl);

                //if(N_STEP == 1431 && Comm::getRank() == 14){
                //    std::cout<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                //    adr_ep_list.dump("check 5", std::cout);
                //}
                
            }
            
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
                adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            }
        }
    }


    template<class TSM, class Ttp, class Tep, class Tsp, class Tmorton_key, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClosePMM(TagChopLeafTrue,
                               TagCopyInfoCloseWithTpAdrptcl,
                               const ReallocatableArray<Ttp> & tp_first,
                               const S32 adr_ptcl,
                               const S32 n_ptcl,
                               const ReallocatableArray<Tep> & ep_first,
                               ReallocatableArray<S32> & adr_ep_list,
                               const ReallocatableArray<Tsp> & sp_first,
                               ReallocatableArray<S32> & adr_sp_list,
                               const TargetBox<TSM> & target_box,
                               const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClosePMM C ***************************");
        PS_COMPILE_TIME_MESSAGE("TagChopLeafTrue and TagCopyInfoCloseWithTpAdrptcl, @ CopyInfoClosePMM")
#if !defined(NO_PUSH_BACK)
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
#endif
        S32 cnt_adr_ptcl = adr_ptcl;        
        if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                    const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                    if (target_box.contains(ep_first[adr_ep].getPos(), morton_key)) {
                        adr_ep_list.pushBackNoCheck(adr_ep);
                    }
                }
                else{
                    const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                    if (target_box.contains(sp_first[adr_sp].getPos(), morton_key)) {
                        adr_sp_list.pushBackNoCheck(adr_sp);
                    }
                }
            }
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
	    // under construction
	    assert(0);
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
	    // under construction
	    assert(0);
            // to do need to chop particles 
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                    const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                    adr_ep_list.pushBackNoCheck(adr_ep);
                }
                else{
                    const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                    adr_sp_list.pushBackNoCheck(adr_sp);
                }
            }
        }
    }



    template<class TSM, class Ttp, class Tep, class Tsp, class Tmorton_key, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClosePMM(TagChopLeafFalse,
                               TagCopyInfoCloseWithTpAdrptcl,
                               const ReallocatableArray<Ttp> & tp_first,
                               const S32 adr_ptcl,
                               const S32 n_ptcl,
                               const ReallocatableArray<Tep> & ep_first,
                               ReallocatableArray<S32> & adr_ep_list,
                               const ReallocatableArray<Tsp> & sp_first,
                               ReallocatableArray<S32> & adr_sp_list,
                               const TargetBox<TSM> & target_box,
                               const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClosePMM D ***************************");
        PS_COMPILE_TIME_MESSAGE("B) TagChopLeafFalse and TagCopyInfoCloseWithTpAdrptcl, @ CopyInfoClosePMM")
#if !defined(NO_PUSH_BACK)
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
#endif
        S32 cnt_adr_ptcl = adr_ptcl;        
        if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                    const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                    //adr_ep_list.pushBackNoCheck(adr_ep);
                    adr_ep_list.push_back(adr_ep);
                }
                else{
                    const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                    adr_sp_list.pushBackNoCheck(adr_sp);
                }
            }
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceShort>) {
        } else if constexpr (std::is_same_v<typename TSM::force_type, TagForceLong>) {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                    const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                    adr_ep_list.pushBackNoCheck(adr_ep);
                }
                else{
                    const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                    adr_sp_list.pushBackNoCheck(adr_sp);
                }
            }
        }
    }

    
    
    // PMMM
    template<class TSM, class Ttp, class Tep, class Tsp, class Tmorton_key, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose E ***************************");
        PS_COMPILE_TIME_MESSAGE("TagForceLong, TagChopLeafTrue and TagCopyInfoCloseNoSp, @ CopyInfoClose")
#ifndef NO_PUSH_BACK
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
#endif
        S32 cnt_adr_ptcl = adr_ptcl;

        if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
                const KeyT key = tp_first[cnt_adr_ptcl].key_;
                if (target_box.contains(key, morton_key)) {
                    adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
                }
            }
        } else {
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
                if (target_box.contains(ep_first[cnt_adr_ptcl].getPos(), morton_key)) {
                    adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
                }
            }
        }
    }

    // original
    template<class TSM, class Ttp, class Tep, class Tsp, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose F ***************************");
        PS_COMPILE_TIME_MESSAGE("TagForceLong, TagChopLeafTrue and TagCopyInfoCloseNoSp, @ CopyInfoClose")
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
            adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
        }
    }





    // PMMM
    template<class TSM, class Ttp, class Tep, class Tsp, class Tmorton_key, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafFalse,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose G ***************************");
        PS_COMPILE_TIME_MESSAGE("TagForceLong, TagChopLeafFalse and TagCopyInfoCloseNoSp, @ CopyInfoClose")
#ifndef NO_PUSH_BACK
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
#endif
        S32 cnt_adr_ptcl = adr_ptcl;
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
            adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
        }
    }

    template<class TSM, class Ttp, class Tep, class Tsp, class Tmorton_key, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              TagCopyInfoCloseWithTpAdrptcl,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose H ***************************");
        PS_COMPILE_TIME_MESSAGE("TagForceLong, TagChopLeafTrue and TagCopyInfoCloseWithTpAdrptcl, @ CopyInfoClose")
#ifndef NO_PUSH_BACK
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
#endif
        S32 cnt_adr_ptcl = adr_ptcl;
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                if (target_box.contains(ep_first[adr_ep].getPos(), morton_key)) {
                    adr_ep_list.pushBackNoCheck(adr_ep);
                }
            }
            else{
                const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                if (target_box.contains(sp_first[adr_sp].getPos(), morton_key)) {
                    adr_sp_list.pushBackNoCheck(adr_sp);
                }
            }
        }
    }


    // old
    template<class TSM, class Ttp, class Tep, class Tsp, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafFalse,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose I ***************************");
        PS_COMPILE_TIME_MESSAGE("TagForceLong, TagChopLeafFalse and TagCopyInfoCloseNoSp, @ CopyInfoClose")
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
            adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
        }
    }

    // NORMAL LONG mode
    template<class TSM, class Ttp, class Tep, class Tsp, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              TagCopyInfoCloseWithTpAdrptcl,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose J ***************************");
        //std::cout<<"CopyInfoClose B"<<std::endl;
        //std::cout<<"n_ptcl= "<<n_ptcl<<std::endl;
        S32 cnt_adr_ptcl = adr_ptcl;
        //std::cout<<"cnt_adr_ptcl= "<<cnt_adr_ptcl<<std::endl;
        //adr_ep_list.dump("adr_ep_list", std::cout);
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        //adr_sp_list.dump("adr_sp_list", std::cout);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        //std::cout<<"befor copy"<<std::endl;
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                adr_ep_list.pushBackNoCheck(adr_ep);
            }
            else{
                const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                adr_sp_list.pushBackNoCheck(adr_sp);
            }
        }
    }
    template<class TSM, class Ttp, class Tep, class Tsp, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafFalse,
                              TagCopyInfoCloseWithTpAdrptcl,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose K ***************************");
        PS_COMPILE_TIME_MESSAGE("TagForceLong, TagChopLeafFalse and TagCopyInfoCloseWithTpAdrptcl, @ CopyInfoClose")
        //std::cout<<"CopyInfoClose C"<<std::endl;        
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                adr_ep_list.pushBackNoCheck(adr_ep);
            }
            else{
                const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                adr_sp_list.pushBackNoCheck(adr_sp);
            }
        }
    }

    template<class TSM, class Ttp, class Tep, class Tsp, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceShort,
                              TagChopLeafTrue,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose L ***************************");
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++){
            const S32 adr = adr_ptcl + ip;
            /*
#if 1
	    const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box, ep_first[adr].getPos(), len_peri);
            if( dis_sq <= 0.0 ){
#else
            if( target_box.isInEpList(ep_first, adr) ){
#endif
            */
            //if( target_box.isInEpList<Tep, CALC_DISTANCE_TYPE>(ep_first, adr, len_peri) ){
            if( target_box.template isInEpList<Tep, CALC_DISTANCE_TYPE>(ep_first, adr, len_peri) ){
                adr_ep_list.pushBackNoCheck(adr);
            }

        }
    }
    template<class TSM, class Ttp, class Tep, class Tsp, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void CopyInfoClose(TagForceShort,
                              TagChopLeafFalse,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box,
                              const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("*************************** CopyInfoClose M ***************************");
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++){
            const S32 adr = adr_ptcl + ip;
            adr_ep_list.pushBackNoCheck(adr);
        }
    }
    // COPY INFO CLOSE
    ///////////
    
    /////////////////
    // GET OPEN BITS

    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBit(TagSearchLong,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBit A ***************************");
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, pos, len_peri);
            const F64 size = tc_first[adr_tc+i].geo_.getSize();
            const F64 r_crit_sq = size*size*inv_theta_sq;
            open_bit |= (dis_sq <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
        return open_bit;
    }

    // for consistency with PMM
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc, class Tmorton_key>
    inline U32 GetOpenBit(TagSearchLong,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq,
                          const S32 lev_crit,
                          const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBit B ***************************");
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
	    const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, pos, len_peri);
            const F64 size = tc_first[adr_tc+i].geo_.getSize();
            const F64 r_crit_sq = size*size*inv_theta_sq;
            open_bit |= (dis_sq <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
        return open_bit;
    }
    
    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBit(TagSearchLongScatter,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBit C ***************************");
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            if(  tc_first[adr_tc+i].n_ptcl_ >0){
                const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
                const F64 size = tc_first[adr_tc+i].geo_.getSize();
                const F64 r_crit_sq = size*size*inv_theta_sq;
                const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, pos, len_peri);
                //if(target_box.vertex_.low_.x > 0.77815 && target_box.vertex_.low_.x < 0.77816){
                //if(target_box.vertex_.low_.x > 0.590872 && target_box.vertex_.low_.x < 0.590873){
                //if(Comm::getRank()==6){
                //if(Comm::getRank()==0){
                /*
                if(Comm::getRank()==7){
                    std::cout<<"??????????????????????????????"<<std::endl;
                    std::cout<<"sqrt(r_crit_sq)= "<<sqrt(r_crit_sq)<<std::endl;
                    std::cout<<"sqrt(dis_sq)= "<<sqrt(dis_sq)<<std::endl;
                    std::cout<<"pos= "<<pos<<std::endl;
                    std::cout<<"len_peri= "<<len_peri<<std::endl;
                    std::cout<<"target_box.vertex_= "<<target_box.vertex_<<std::endl;
                    std::cout<<"tc_first[adr_tc+i].dump()"<<std::endl;
                    tc_first[adr_tc+i].dump(std::cout);
                    std::cout<<"GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri)= "<<GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri)<<std::endl;
                    std::cout<<"??????????????????????????????"<<std::endl;
                }
                */
                
                open_bit |= ( GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri) <= 0.0 || dis_sq <= r_crit_sq) << i;
                //const F64 dis_sq = target_box.vertex_.getDistanceMinSq(pos);
                //open_bit |= ( target_box.vertex_.overlapped(tc_first[adr_tc+i].geo_.getVertexOut()) || dis_sq <= r_crit_sq) << i;
                //if(Comm::getRank()==0){
                //  std::cerr<<"target_box.vertex_= "<<target_box.vertex_<<" tc_first[adr_tc+i].geo_.getVertexOut()= "<<tc_first[adr_tc+i].geo_.getVertexOut()
                //	       <<" dis= "<<GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri)
                //	       <<std::endl;
                //}
                open_bit |= ( 1<< (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        return open_bit;
    }

    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBit(TagSearchLongSymmetry,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBit D ***************************");
        // long symmetry
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
#if 0
            const F64 dis_sq = target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].geo_.getVertexOut());
#else
            //const F64 dis_sq = target_box.vertex_in_.getDistanceMinSq(pos);
            const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_in_, pos, len_peri);
#endif
            const F64 size = tc_first[adr_tc+i].geo_.getSize();
            const F64 r_crit_sq = size*size*inv_theta_sq;
#if 0
            open_bit |= (dis_sq <= r_crit_sq) << i; // OK
            //open_bit |= (target_box.vertex_out_.overlapped(tc_first[adr_tc+i].geo_.getVertexIn()) || dis_sq <= r_crit_sq) << i;
#else
            //open_bit |= ( (target_box.vertex_in_.overlapped(tc_first[adr_tc+i].geo_.getVertexOut())
            //               || target_box.vertex_out_.overlapped(tc_first[adr_tc+i].geo_.getVertexIn()) )
            //              || dis_sq <= r_crit_sq) << i;
            /*
            open_bit |= ( GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_in_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri) <= 0.0
                          || GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_out_, tc_first[adr_tc+i].geo_.getVertexIn(), len_peri) <= 0.0
                          || dis_sq <= r_crit_sq) << i;
            */
            open_bit |= ( IsOverlapped<CALC_DISTANCE_TYPE>(target_box.vertex_in_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri)
                          || IsOverlapped<CALC_DISTANCE_TYPE>(target_box.vertex_out_, tc_first[adr_tc+i].geo_.getVertexIn(), len_peri)
                          || dis_sq <= r_crit_sq ) << i;
#endif
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBit(TagSearchLongCutoff,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBit E ***************************");
        U32 open_bit = 0;
        if (inv_theta_sq > 0.0) {
            // theta_ > 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
#if 1
                const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
                const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, pos, len_peri);
                const F64 size = tc_first[adr_tc+i].geo_.getSize();
                const F64 r_crit_sq = size*size*inv_theta_sq;
                open_bit |= (dis_sq <= r_crit_sq) << i;
                const F64 dis_sq_2 = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri);
                open_bit |= ((dis_sq_2 <= 0.0) && (tc_first[adr_tc+i].n_ptcl_ > 0)) << (i + N_CHILDREN); // if true, it should be checked
#else
                const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
                const F64 dis_sq = target_box.vertex_.getDistanceMinSq(pos);
                const F64 size = tc_first[adr_tc+i].geo_.getSize();
                const F64 r_crit_sq = size*size*inv_theta_sq;
                open_bit |= (dis_sq <= r_crit_sq) << i;
                //open_bit |= ( ( (target_box.vertex_.overlapped( tc_first[adr_tc+i].mom_.getVertexOut()))
                //                && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                //              << (i + N_CHILDREN) ); // if true, it should be checked
                open_bit |= ( ( (target_box.vertex_.overlapped( tc_first[adr_tc+i].geo_.getVertexOut()))
                                && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                              << (i + N_CHILDREN) ); // if true, it should be checked
#endif
            }
        }
        else {
            // theta_ = 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                open_bit |= 1 << i;
                //open_bit |= ( ( (target_box.vertex_.overlapped( tc_first[adr_tc+i].mom_.getVertexOut()))
                //                && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                //              << (i + N_CHILDREN) ); // if true, it should be checked
                open_bit |= ( ( (target_box.vertex_.overlapped( tc_first[adr_tc+i].geo_.getVertexOut()))
                                && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                              << (i + N_CHILDREN) ); // if true, it should be checked                
            }
        }
        return open_bit;
    }

    // PMMM
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc, class Tmorton_key>
    inline U32 GetOpenBitPMM(TagSearchLongParticleMeshMultipole,
			     const ReallocatableArray<Ttc> & tc_first,
			     const S32 adr_tc,
			     const TargetBox<TSM> & target_box,
			     const F64vec & len_peri,
			     const F64 inv_theta_sq,
			     const S32 lev_crit,
			     const Tmorton_key & morton_key){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBitPMM F ***************************");
        U32 open_bit = 0;
        if (inv_theta_sq >= 0.0) {
            // theta_ > 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_ptcl = tc_first[adr_tc+i].n_ptcl_;
                const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
                const F64 dis_sq = target_box.vertex_in_.getDistanceMinSq(pos); // type of vertex_in_ is F64ort
                const F64 size = tc_first[adr_tc+i].geo_.getSize();
                const F64 r_crit_sq = size*size*inv_theta_sq;
                const F64ort tc_vertex_in = tc_first[adr_tc+i].geo_.getVertexIn();
                if (n_ptcl > 0) {
                    S32ort box;
                    box.low_  = morton_key.getPMCellID(tc_vertex_in.low_);
                    box.high_ = morton_key.getPMCellID(tc_vertex_in.high_);
                    open_bit |= ( (dis_sq <= r_crit_sq) ||
                                  (target_box.vertex_in_.overlapped(tc_vertex_in)) ||
                                  (target_box.vertex_out_.doesNotContain(box)) ||  // type of vertex_out is S32ort. compare using integer
                                  (tc_first[adr_tc+i].level_ < lev_crit)) << i;
                    open_bit |= ( target_box.vertex_out_.overlapped( box )
                                  << (i + N_CHILDREN) ); // if true, it should be checked
		    // The above criterion is a little bit complex
		    // it means: the tree cell is not opend if target_box.vertex_out_ completely includes tc_vertex_in (i.e. target_box.vertex_out_.Contain(box) )
		    
#if 0
                    // for debug
                    if (f_debug_GetOpenBit && i == 0) {
                        std::cout << "target_box.vertex_in_.low_  = " << target_box.vertex_in_.low_ << std::endl;
                        std::cout << "target_box.vertex_in_.high_ = " << target_box.vertex_in_.high_ << std::endl;
                        std::cout << "tc_first[].vertex_in.low_   = " << tc_vertex_in.low_ << std::endl;
                        std::cout << "tc_first[].vertex_in.high_  = " << tc_vertex_in.high_ << std::endl;
                        std::cout << "tc_first[].level_ = " << tc_first[adr_tc+i].level_ << std::endl;
                        std::cout << "lev_crit          = " << lev_crit << std::endl;
                        std::cout << "tc_first[].pos = " << pos << std::endl; 
                        std::cout << "dis = " << dis << std::endl; 
                        std::cout << "r_crit_sq = " << r_crit_sq << std::endl;
                        std::cout << "box.low_  = " << box.low_ << std::endl; 
                        std::cout << "box.high_ = " << box.high_ << std::endl;
                        std::cout << "target_box.vertex_out_.low_  = " << target_box.vertex_out_.low_ << std::endl;
                        std::cout << "target_box.vertex_out_.high_ = " << target_box.vertex_out_.high_ << std::endl;
                        std::cout << "Cond. (1) = " << (dis <= r_crit_sq) << std::endl;
                        std::cout << "Cond. (2) = " << (target_box.vertex_in_.overlapped(tc_vertex_in)) << std::endl;
                        std::cout << "Cond. (3) = " << (target_box.vertex_out_.doesNotContain(box)) << std::endl;
                        std::cout << "Cond. (4) = " << (tc_first[adr_tc+i].level_ < lev_crit) << std::endl;
                    }
                    // for debug
#endif
                }
                // [Notes]
                // The above implementation is almost the same as the implementation
                // shown below, in which vertex_out_ is class F64ort. However,
                // the following implementation is not secure. 
                // open_bit |= ( (dis_sq <= r_crit_sq) ||
                //               (target_box.vertex_out_.doesNotContain(tc_vertex_in)) ||
                //               (tc_first[adr_tc+i].level_ < lev_crit)) << i;
                // open_bit |= ( ( (target_box.vertex_out_.overlapped( tc_vertex_in ))
                //                 && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                //               << (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        else {
            // theta_ = 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_ptcl = tc_first[adr_tc+i].n_ptcl_;
                if (n_ptcl > 0) {
                    const F64ort tc_vertex_in = tc_first[adr_tc+i].geo_.getVertexIn();
                    S32ort box;
                    box.low_ = morton_key.getPMCellID(tc_vertex_in.low_);
                    box.high_ = morton_key.getPMCellID(tc_vertex_in.high_);
                    // Note that `box` cannot be calculated correctly when n_ptcl_ = 0
                    // since tc_vertex_in is not defined for n_ptcl_ = 0.
                    open_bit |= 1 << i;
                    open_bit |= ( target_box.vertex_out_.overlapped( box ) 
                                  << (i + N_CHILDREN) ); // if true, it should be checked
                }
                // [Notes]
                // The above implementation is almost the same as the implementation
                // shown below, in which vertex_out_ is class F64ort. However,
                // the following implementation is not secure. 
                // open_bit |= 1 << i;
                // open_bit |= ( ( (target_box.vertex_out_.overlapped( tc_first[adr_tc+i].mom_.getVertexIn()))
                //                 && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                //               << (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        return open_bit;
    }
    // PMMM
    
    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBit(TagSearchShortScatter,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBit G ***************************");
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
#if 1
            const auto dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_in_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri);
            open_bit |= ( dis_sq <= 0.0 ) << i;
#else
            open_bit |= ( target_box.vertex_in_.overlapped(tc_first[adr_tc+i].geo_.getVertexOut()) ) << i;
#endif
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBit(TagSearchShortSymmetry,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBit H ***************************");
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
#if 1
            const auto dis_sq = std::min( GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_in_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri),
                                          GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_out_, tc_first[adr_tc+i].geo_.getVertexIn(), len_peri) );
            open_bit |= ( dis_sq <= 0.0 ) << i;
#else
            //open_bit |= ( target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut())
            //              || target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) ) << i;
            open_bit |= ( target_box.vertex_in_.overlapped(tc_first[adr_tc+i].geo_.getVertexOut())
                          || target_box.vertex_out_.overlapped(tc_first[adr_tc+i].geo_.getVertexIn()) ) << i;
#endif
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBit(TagSearchShortGather,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        PS_COMPILE_TIME_MESSAGE("*************************** GetOpenBit I ***************************");
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
#if 1
            const auto dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_out_, tc_first[adr_tc+i].geo_.getVertexIn(), len_peri);
            open_bit |= ( dis_sq <= 0.0 ) << i;
#else
            //open_bit |= ( target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) ) << i;
            open_bit |= ( target_box.vertex_out_.overlapped(tc_first[adr_tc+i].geo_.getVertexIn()) ) << i;
#endif
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }
    // GET OPEN BITS
    /////////////////

    // PMMM
    // for both long mode and short mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tmorton_key, class Tchopleaf, class Tchopnonleaf, class Tcopyinfoclose, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void MakeListUsingTreeRecursivePMM
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     const TargetBox<TSM> & target_box,
     const S32 n_leaf_limit,
     const S32 lev_leaf_limit,
     const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
     const F64vec & len_peri,
     const F64 theta,
     const S32 lev_crit,
     const Tmorton_key & morton_key){
        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl_child = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        const F64 inv_theta_sq = (theta > 0.0) ? 1.0 / (theta*theta) : -1.0;
        if( !(tc_cur->isLeafPMM(n_leaf_limit, lev_leaf_limit)) ){ // not leaf
            //U32 open_bit = GetOpenBit<CALC_DISTANCE_TYPE>(typename TSM::search_type(), tc_first, adr_tc_child, target_box, len_peri, inv_theta_sq, lev_crit, morton_key);
            U32 open_bit = GetOpenBitPMM<CALC_DISTANCE_TYPE>(typename TSM::search_type(), tc_first, adr_tc_child, target_box, len_peri, inv_theta_sq, lev_crit, morton_key);
            for(S32 i=0; i<N_CHILDREN; i++){
                if( !((open_bit>>(i+N_CHILDREN)) & 0x1) ) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    MakeListUsingTreeRecursivePMM
                        <TSM, Ttc, Ttp, Tep, Tsp, Tmorton_key, Tchopleaf, Tchopnonleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>
                        (tc_first, adr_tc_child+i, tp_first, ep_first, adr_ep_list, sp_first, adr_sp_list,
                         target_box, n_leaf_limit, lev_leaf_limit, adr_tree_sp_first, len_peri, theta, lev_crit,
                         morton_key);
                }
                else{ // far
                    CopyInfoDistant2(adr_tc_child+i, tc_first, adr_tc_child+adr_tree_sp_first+i, adr_sp_list, target_box, morton_key);
                }
            }
        }
        else{ //leaf
            CopyInfoClosePMM<TSM, Ttp, Tep, Tsp, Tmorton_key, CALC_DISTANCE_TYPE>(Tchopleaf(), Tcopyinfoclose(), tp_first, adr_ptcl_child, n_ptcl, ep_first, adr_ep_list, sp_first, adr_sp_list, target_box, morton_key);
        }
    }

    // original
    // new
    // for both long mode and short mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tchopleaf, class Tcopyinfoclose, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void MakeListUsingTreeRecursive
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     const TargetBox<TSM> & target_box,
     const S32 n_leaf_limit,
     const S32 adr_tree_sp_first, // address of the first sp from the (global) tree.
     const F64vec & len_peri,
     const F64 theta){
        //std::cout<<"MakeListUsingTreeRecursive"<<std::endl;
        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl_child = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        const F64 inv_theta_sq = (theta > 0.0) ? 1.0 / (theta*theta) : -1.0;
        if( !(tc_cur->isLeaf(n_leaf_limit)) ){ // not leaf
            //std::cout<<"not leaf"<<std::endl;
            U32 open_bit = GetOpenBit<TSM, Ttc, CALC_DISTANCE_TYPE>(typename TSM::search_type(), tc_first, adr_tc_child, target_box, len_peri, inv_theta_sq);
            //std::cout<<"open_bit= "<<std::hex<<open_bit<<std::endl;
            //std::cout<<std::dec;
            for(S32 i=0; i<N_CHILDREN; i++){
                if( !((open_bit>>(i+N_CHILDREN)) & 0x1) ) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp, Tchopleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>
                        (tc_first, adr_tc_child+i, tp_first, ep_first, adr_ep_list, sp_first, adr_sp_list,
                         target_box, n_leaf_limit, adr_tree_sp_first, len_peri, theta);
                } else { // far
                    CopyInfoDistant(typename TSM::force_type(), adr_tc_child+adr_tree_sp_first+i, adr_sp_list);
                }
            }
        }
        else{ //leaf
            //std::cout<<"leaf0 "<<std::endl;
            CopyInfoClose<TSM, Ttp, Tep, Tsp, CALC_DISTANCE_TYPE>
                (typename TSM::force_type(), Tchopleaf(), Tcopyinfoclose(), tp_first, adr_ptcl_child, n_ptcl,
                 ep_first, adr_ep_list, sp_first, adr_sp_list, target_box, len_peri);
            //std::cout<<"leaf1 "<<std::endl;
        }
    }


    // PMMM
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tmorton_key, class Tchopleaf, class Tchopnonleaf, class Tcopyinfoclose, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void MakeListUsingTreeRecursiveTopPMM
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     const TargetBox<TSM> & target_box,
     const S32 n_leaf_limit,
     const S32 lev_leaf_limit,
     const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
     const F64vec & len_peri,
     const F64 theta,
     const S32 lev_crit,
     const Tmorton_key & morton_key){
        if ((theta > 0.0) ||
            (typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF)) ||
            (typeid(TSM) == typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE))) {
            if( IsOpen(typename TSM::search_type(), tc_first, adr_tc, target_box, morton_key) ){
                MakeListUsingTreeRecursivePMM
                    <TSM, Ttc, Ttp, Tep, Tsp, Tmorton_key, Tchopleaf, Tchopnonleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>
                    (tc_first, adr_tc, tp_first, ep_first, adr_ep_list, sp_first,
                     adr_sp_list, target_box, n_leaf_limit, lev_leaf_limit,
                     adr_tree_sp_first, len_peri, theta, lev_crit, morton_key);
            }
        }
        else {
            const S32 n_ptcl = tc_first[0].n_ptcl_;
            S32 adr_ptcl = tc_first[0].adr_ptcl_;
            CopyInfoClosePMM<TSM, Ttp, Tep, Tsp, Tmorton_key, CALC_DISTANCE_TYPE>
                (Tchopleaf(), Tcopyinfoclose(),
                 tp_first, adr_ptcl, n_ptcl,
                 ep_first, adr_ep_list, sp_first, adr_sp_list,
                 target_box, morton_key);            
        }
    }    

    // for long mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tchopleaf, class Tcopyinfoclose, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void MakeListUsingTreeRecursiveTop
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     const TargetBox<TSM> & target_box,
     const S32 n_leaf_limit,
     const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
     const F64vec & len_peri,
     const F64 theta){
        PS_COMPILE_TIME_MESSAGE("MakeListUsingTreeRecursiveTop for Long mode and short mode only if exchangeLET");
        //std::cerr<<"CALC_DISTANCE_TYPE_NORMAL= "<<CALC_DISTANCE_TYPE_NORMAL<<std::endl;
        if ((theta > 0.0) || (typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF))) {
            // theta_ > 0 case or PS::SEARCH_MODE_LONG_CUTOFF
            if( IsOpen<CALC_DISTANCE_TYPE>(typename TSM::search_type(), tc_first, adr_tc, target_box, len_peri) ){
                MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp, Tchopleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>
                    (tc_first, adr_tc, tp_first, ep_first, adr_ep_list, sp_first,
                     adr_sp_list, target_box, n_leaf_limit,
                     adr_tree_sp_first, len_peri, theta);
            }
        }
        else {
            // theta_ = 0 case with the modes other than PS::SEARCH_MODE_LONG_CUTOFF
            const S32 n_ptcl = tc_first[0].n_ptcl_;
            S32 adr_ptcl = tc_first[0].adr_ptcl_;
            CopyInfoClose<TSM, Ttp, Tep, Tsp, CALC_DISTANCE_TYPE>
                (typename TSM::force_type(), Tchopleaf(), Tcopyinfoclose(), tp_first, adr_ptcl, n_ptcl,
                 ep_first, adr_ep_list, sp_first, adr_sp_list, target_box, len_peri);
        }
    }
    

    // for short mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Tchopleaf, class Tcopyinfoclose, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void MakeListUsingTreeShortRecursiveTop
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const TargetBox<TSM> & target_box,
     const S32 n_leaf_limit,
     const F64vec & len_peri){
        PS_COMPILE_TIME_MESSAGE("MakeListUsingTreeRecursiveTop for Short mode");
        const F64 theta_tmp = 1.0;
        S32 adr_tree_sp_first = 0;
        ReallocatableArray<SuperParticleBase> sp_first;
        ReallocatableArray<S32> adr_sp_list;

        MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, SuperParticleBase, Tchopleaf, TagCopyInfoCloseNoSp, CALC_DISTANCE_TYPE>
            (tc_first, adr_tc, tp_first, ep_first, adr_ep_list, sp_first,
             adr_sp_list, target_box, n_leaf_limit,
             adr_tree_sp_first, len_peri, theta_tmp);

    }
    
    ////////////////////
    // FOR DOUBLE WALK
    //////////////////////////
    // for short symmetry or gather mode
    template<class Ttc>
    inline void IsOverlapped(const Ttc * tc_first,
                             const S32 adr_tc,
                             const F64ort & pos_box,
                             const S32 n_leaf_limit,
                             U32 & is_overlapped){
        if(is_overlapped) return;
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            //open_bit |= (tc_first[adr_tc+i].mom_.getVertexOut().overlapped(pos_box) << i);
            open_bit |= (tc_first[adr_tc+i].geo_.getVertexOut().overlapped(pos_box) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if(is_overlapped) return;
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first + adr_tc_child;
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            if( (open_bit>>i) & 0x1){
                if( !(tc_child->isLeaf(n_leaf_limit)) ){
                    // not leaf
                    IsOverlapped(tc_first, tc_first[adr_tc_child].adr_tc_,
                                 pos_box, n_leaf_limit, is_overlapped);
                }
                else{
                    // leaf and overlapped
                    is_overlapped = 1;
                    return;
                }
            }
        }
    }

    template<class Ttc>
    inline U32 IsOverlappedTop(const ReallocatableArray<Ttc> & tc_first,
                               const F64ort & pos_box,
                               const S32 n_leaf_limit){
        U32 is_overlapped = 0;
        S32 adr_tc = N_CHILDREN;
        if( !tc_first[0].isLeaf(n_leaf_limit) ){
            IsOverlapped(tc_first.getPointer(), adr_tc, pos_box, n_leaf_limit, is_overlapped);
        }
        else{
            //is_overlapped = tc_first[0].mom_.getVertexOut().overlapped(pos_box);
            is_overlapped = tc_first[0].geo_.getVertexOut().overlapped(pos_box);
        }
        return is_overlapped;
    }
    
    template<class Ttca, class Ttcb>
    inline U32 GetOpenBitDoubleWalk(const ReallocatableArray<Ttca> & tc_first_A,
                                    const ReallocatableArray<Ttcb> & tc_first_B,
                                    const S32 adr_tc_B,
                                    const F64ort & pos_domain,
                                    const S32 n_leaf_limit_A){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= (tc_first_B[adr_tc_B+i].geo_.getVertexOut().overlapped(pos_domain)
                         || IsOverlappedTop(tc_first_A, tc_first_B[adr_tc_B+i].geo_.getVertexIn(), n_leaf_limit_A)) << i;
            open_bit |= (tc_first_B[adr_tc_B+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if false, just skip
        }
        return open_bit;
    }
    
    template<class Ttca, class Ttcb, class Tepb>
    inline void MakeListDoubleWalk(const ReallocatableArray<Ttca> & tc_first_A,
                                   const ReallocatableArray<Ttcb> & tc_first_B,
                                   const S32 adr_tc_B,
                                   const ReallocatableArray<Tepb> & ep_first_B,
                                   const F64ort & pos_domain,
                                   const S32 n_leaf_limit_A,
                                   const S32 n_leaf_limit_B,
                                   ReallocatableArray<S32> & adr_ptcl_send){
        U32 open_bit = GetOpenBitDoubleWalk(tc_first_A, tc_first_B, adr_tc_B, pos_domain, n_leaf_limit_A);
        for(S32 i=0; i<N_CHILDREN; i++){
            if( !((open_bit>>(i+N_CHILDREN)) & 0x1) ) continue;
            else if( (open_bit>>i) & 0x1){
                const S32 adr_tc_child = adr_tc_B + i;
                const Ttcb * tc_child = tc_first_B.getPointer() + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if( !(tc_child->isLeaf(n_leaf_limit_B)) ){
                    MakeListDoubleWalk(tc_first_A, tc_first_B,
                                       tc_first_B[adr_tc_child].adr_tc_,
                                       ep_first_B, pos_domain, n_leaf_limit_A, n_leaf_limit_B,
                                       adr_ptcl_send);
                }
                else{
                    //if( pos_domain.overlapped(tc_first_B[adr_tc_child].mom_.getVertexOut()) ){
                    if( pos_domain.overlapped(tc_first_B[adr_tc_child].geo_.getVertexOut()) ){
                        //if(Comm::getRank()==0) std::cerr<<"tc_first_B[adr_tc_child].mom_.getVertexOut()= "<<tc_first_B[adr_tc_child].mom_.getVertexOut()<<std::endl;
                        continue;
                    }
                    else{
                        S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                        for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                            adr_ptcl_send.push_back(adr_ptcl_tmp);
                        }
                    }
                }
            }
        }
    }

    
    template<class Ttca, class Ttcb, class Tepb>
    inline void MakeListDoubleWalkTop
    (const ReallocatableArray<Ttca> & tc_first_A, // tree of reciving ptcls
     const ReallocatableArray<Ttcb> & tc_first_B, // tree of assigned ptcls
     const ReallocatableArray<Tepb> & ep_first_B,
     const F64ort & pos_domain,
     const S32 n_leaf_limit_A,
     const S32 n_leaf_limit_B,
     ReallocatableArray<S32> & adr_ptcl_send){
        const S32 adr_tc_B = N_CHILDREN;
        if( !tc_first_B[0].isLeaf(n_leaf_limit_B) ){
            MakeListDoubleWalk(tc_first_A, tc_first_B, adr_tc_B,
                               ep_first_B, pos_domain, n_leaf_limit_A, n_leaf_limit_B,
                               adr_ptcl_send);
        }
        else{
            // original
            //if( tc_first_B[0].mom_.getVertexOut().overlapped(pos_domain) ){
            if( tc_first_B[0].geo_.getVertexOut().overlapped(pos_domain) ){
                // already sent
                return;
            }
            else{
                S32 adr_ptcl_tmp = tc_first_B[0].adr_ptcl_;
                const S32 n_child = tc_first_B[0].n_ptcl_;
                // not sent yet (might contain distant ptcls)
                // this is sufficient condition
                for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                    adr_ptcl_send.push_back(adr_ptcl_tmp);
                }
            }
        }
    }

#ifdef LOOP_TREE
    template<typename Ttc, typename Ttp>
    inline void CopyInfoClose(const Ttc & tc,
                              const ReallocatableArray<Ttp> & tp_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              ReallocatableArray<S32> & adr_sp_list){
        const S32 n_ptcl = tc.n_ptcl_;
        S32 adr_tp = tc.adr_ptcl_;
        for(S32 i=0; i<n_ptcl; i++, adr_tp++){
            U32 adr_ptcl = tp_first[adr_tp].adr_ptcl_;
            if( GetMSB(adr_ptcl) == 0 ){
                adr_ep_list.pushBackNoCheck(adr_ptcl);
            }
            else{
                adr_sp_list.pushBackNoCheck(ClearMSB(adr_ptcl));
            }
        }
    }
    
    template<class TSM, class Ttc, class Ttp>
    inline void MakeListUsingTreeLoop
    (const ReallocatableArray<Ttc> & tc_first,
     const ReallocatableArray<Ttp> & tp_first,
     const TargetBox<TSM> target_box[],
     ReallocatableArray<S32> adr_ep_list[],
     ReallocatableArray<S32> adr_sp_list[],
     const S32 n_ipg,
     const S32 n_leaf_limit,
     const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
     const F64 inv_theta_sq,
     const S32 buf_size_limit = 100000){
        S32 adr_tc[128];
        assert(n_ipg < 128);
        for(S32 i=0; i<n_ipg; i++){
            adr_tc[i] = 0;
        }
        S32 loop = 0;
        while(1){
            S32 n_ipg_rem = n_ipg;
            //std::cerr<<"loop= "<<loop<<std::endl;
            loop++;
            for(S32 i=0; i<n_ipg; i++){
                const U32 adr = adr_tc[i];
                //std::cerr<<"i= "<<i<<" loop= "<<loop<<" adr= "<<adr<<std::endl;
                const Ttc & tc_cur = tc_first[adr];
                const F64vec tc_pos = tc_cur.mom_.getPos();
                const F64 size = tc_cur.geo_.getSize();
                const F64 r_crit_sq = size*size*inv_theta_sq;
                const F64 dis_sq = target_box[i].vertex_.getDistanceMinSq(tc_pos);
                if(r_crit_sq < dis_sq){
                    //std::cerr<<"check A"<<std::endl;
                    // distant
                    adr_sp_list[i].pushBackNoCheck(adr+adr_tree_sp_first);
                    adr_tc[i] = tc_cur.adr_tc_next_;
                }
                else{
                    // close
                    if(tc_cur.isLeaf(n_leaf_limit)){
                        //std::cerr<<"check B: tc_cur.n_ptcl_= "<<tc_cur.n_ptcl_<<std::endl;
                        CopyInfoClose(tc_cur, tp_first, adr_ep_list[i], adr_sp_list[i]);
                        adr_tc[i] = tc_cur.adr_tc_next_;
                    }
                    else{
                        //std::cerr<<"check C"<<std::endl;
                        // not leaf. go deeper.
                        adr_tc[i] = tc_cur.adr_tc_;
                    }
                }
                if(&tc_cur-tc_first.getPointer() == 1){
                    n_ipg_rem--;
                }
            }
            if(n_ipg_rem == 0) break;
        }
    }


    
#endif
    
}
