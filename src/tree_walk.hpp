#pragma once

/*
リスト長がMAKE_LISTとMAKE_LIST_FOR_REUSEで異なる問題について。Leaf
Cellの粒子を全部送るかどうかが原因。TargetBoxにメンバ関数isInEpList(粒
子を追加するかを切り捨てるかを記述)を追加することで対処しようとしたが、
Scatterモードで粒子の切り捨てを行った場合、EXLETではGATHERでもSYMMETRY
でも第一段階通信でScatterモードでツリーウォークを行うので、第二段階通
信と
 */

#include<particle_simulator.hpp>

namespace ParticleSimulator{

    ////////////
    // walk mode
    struct TagWalkModeNormal{};
    struct TagWalkModeNearestImage{};
    struct WALK_MODE_NORMAL{
        typedef TagWalkModeNormal walk_type;
    };
    struct WALK_MODE_NEAREST_IMAGE{
        typedef TagWalkModeNearestImage walk_type;
    };
    // walk mode
    ////////////

    ////////////
    // walk mode
    struct TagChopLeafTrue{};
    struct TagChopLeafFalse{};
    // walk mode
    ////////////

    ////////////
    // targe box class
    template<class TSM>
    struct TargetBox{
        F64ort vertex_;
    };
    template<>
    struct TargetBox<SEARCH_MODE_LONG>{
        F64ort vertex_;
    };
    template<>
    struct TargetBox<SEARCH_MODE_LONG_SYMMETRY>{
        F64ort vertex_in_;
        F64ort vertex_out_;
    };
    template<>
    struct TargetBox<SEARCH_MODE_LONG_CUTOFF>{
        F64ort vertex_;
    };
    template<>
    struct TargetBox<SEARCH_MODE_SCATTER>{
        F64ort vertex_in_;
        template<class Tep>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr) const {
            const F64 dis_sq = vertex_in_.getDistanceMinSQ(ep_first[adr].getPos());
            const F64 r_crit_sq = ep_first[adr].getRSearch() * ep_first[adr].getRSearch();
            return dis_sq <= r_crit_sq;
        }
    };
    template<>
    struct TargetBox<SEARCH_MODE_SYMMETRY>{
        F64ort vertex_out_;
        F64ort vertex_in_;
        template<class Tep>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr) const {
            //std::cerr<<"symmetry"<<std::endl;
            const F64 dis_sq = vertex_in_.getDistanceMinSQ(ep_first[adr].getPos());
            const F64 r_crit_sq = ep_first[adr].getRSearch() * ep_first[adr].getRSearch();
            return (dis_sq <= r_crit_sq) || (vertex_out_.overlapped(ep_first[adr].getPos()));
        }
    };
    template<>
    struct TargetBox<SEARCH_MODE_GATHER>{
        F64ort vertex_out_;
        template<class Tep>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr) const {
            //std::cerr<<"gather"<<std::endl;
            //std::cout<<"ep_first[adr].mass= "<<ep_first[adr].mass<<std::endl;
            //std::cout<<"ep_first[adr].pos= "<<ep_first[adr].pos<<std::endl;
            return vertex_out_.overlapped(ep_first[adr].getPos());
            //return true;
        }
    };
    // targe box class
    ////////////


    /////////////
    // set target box for exching let
    template<class TSM>
    inline void SetTargetBoxExLet(TargetBox<TSM> & box,
                                  const F64ort & in,
                                  const F64ort & out){
        assert(0);
    }
    template<>
    inline void SetTargetBoxExLet<SEARCH_MODE_LONG_SYMMETRY>
    (TargetBox<SEARCH_MODE_LONG_SYMMETRY> & box,
     const F64ort & in,
     const F64ort & out){
        box.vertex_in_  = in;
        box.vertex_out_ = out;
    }
    template<>
    inline void SetTargetBoxExLet<SEARCH_MODE_LONG_SCATTER>
    (TargetBox<SEARCH_MODE_LONG_SCATTER> & box,
     const F64ort & in,
     const F64ort & out){
        box.vertex_ = in;
    }
    template<>
    inline void SetTargetBoxExLet<SEARCH_MODE_LONG>
    (TargetBox<SEARCH_MODE_LONG> & box,
     const F64ort & in,
     const F64ort & out){
        box.vertex_ = in;
    }
    template<>
    inline void SetTargetBoxExLet<SEARCH_MODE_LONG_CUTOFF>
    (TargetBox<SEARCH_MODE_LONG_CUTOFF> & box,
     const F64ort & in,
     const F64ort & out){
        box.vertex_ = in;
    }
    // set target box for exching let
    /////////////

    /////////////
    // set target box for interaction list
    template<class TSM>
    inline void GetTargetBox(const IPGroup<TSM> & ipg,
                      TargetBox<TSM> & target_box){
        assert(0);
    }
    template<>
    inline void GetTargetBox<SEARCH_MODE_SYMMETRY>(const IPGroup<SEARCH_MODE_SYMMETRY> & ipg,
                                                   TargetBox<SEARCH_MODE_SYMMETRY> & target_box){
        target_box.vertex_in_  = ipg.vertex_in;
        target_box.vertex_out_ = ipg.vertex_;
    }
    template<>
    inline void GetTargetBox<SEARCH_MODE_SCATTER>(const IPGroup<SEARCH_MODE_SCATTER> & ipg,
                                                  TargetBox<SEARCH_MODE_SCATTER> & target_box){
        target_box.vertex_in_ = ipg.vertex_;
    }
    template<>
    inline void GetTargetBox<SEARCH_MODE_GATHER>(const IPGroup<SEARCH_MODE_GATHER> & ipg,
                                                 TargetBox<SEARCH_MODE_GATHER> & target_box){
        target_box.vertex_out_ = ipg.vertex_;
    }
    template<>
    inline void GetTargetBox<SEARCH_MODE_LONG>(const IPGroup<SEARCH_MODE_LONG> & ipg,
                                               TargetBox<SEARCH_MODE_LONG> & target_box){
        target_box.vertex_ = ipg.vertex_;
    }
    template<>
    inline void GetTargetBox<SEARCH_MODE_LONG_CUTOFF>(const IPGroup<SEARCH_MODE_LONG_CUTOFF> & ipg,
                                                      TargetBox<SEARCH_MODE_LONG_CUTOFF> & target_box){
        target_box.vertex_ = ipg.vertex_;
    }
    template<>
    inline void GetTargetBox<SEARCH_MODE_LONG_SCATTER>(const IPGroup<SEARCH_MODE_LONG_SCATTER> & ipg,
                                                TargetBox<SEARCH_MODE_LONG_SCATTER> & target_box){
        target_box.vertex_ = ipg.vertex_;
    }
    template<>
    inline void GetTargetBox<SEARCH_MODE_LONG_SYMMETRY>(const IPGroup<SEARCH_MODE_LONG_SYMMETRY> & ipg,
                                                 TargetBox<SEARCH_MODE_LONG_SYMMETRY> & target_box){
        target_box.vertex_in_  = ipg.vertex_;
        target_box.vertex_out_ = ipg.vertex_out_;
    }
    // set target box for interaction list
    /////////////
    
    ///////////
    // IS OPEN
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchLong,
                const ReallocatableArray<Ttc> & tc_first,
                const S32 adr_tc,
                const TargetBox<TSM> & target_box,
                const F64 r_crit_sq,
                const F64vec & len_peri,
                TagWalkModeNormal){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }
    
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongScatter,
                const ReallocatableArray<Ttc> & tc_first,
                const S32 adr_tc,
                const TargetBox<TSM> & target_box,
                const F64 r_crit_sq,
                const F64vec & len_peri,
                TagWalkModeNormal){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }

    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongSymmetry,
                const ReallocatableArray<Ttc> & tc_first,
                const S32 adr_tc,
                const TargetBox<TSM> & target_box,
                const F64 r_crit_sq,
                const F64vec & len_peri,
                TagWalkModeNormal){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }
    
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongCutoff,
                const ReallocatableArray<Ttc> & tc_first,
                const S32 adr_tc,
                const TargetBox<TSM> & target_box,
                const F64 r_crit_sq,
                const F64vec & len_peri,
                TagWalkModeNormal){
        return ( (target_box.vertex_.overlapped(tc_first[adr_tc].mom_.getVertexOut()))
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
    }
    
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortScatter,
                const ReallocatableArray<Ttc> & tc_first,
                const S32 adr_tc,
                //const F64ort & pos_target_box,
                const TargetBox<TSM> & target_box,
                const F64 r_crit_sq,
                const F64vec & len_peri,
                TagWalkModeNormal){
        return ( (target_box.vertex_in_.overlapped(tc_first[adr_tc].mom_.getVertexOut()))
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
    }
    
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortSymmetry,
                const ReallocatableArray<Ttc> & tc_first,
                const S32 adr_tc,
                //const F64ort & pos_target_box,
                const TargetBox<TSM> & target_box,
                const F64 r_crit_sq,
                const F64vec & len_peri,
                TagWalkModeNormal){
        return ( ( (target_box.vertex_in_.overlapped(tc_first[adr_tc].mom_.getVertexOut()))
                 || (target_box.vertex_out_.overlapped(tc_first[adr_tc].mom_.getVertexIn())) )
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
        /*
        F64 dis = std::min(target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc].mom_.getVertexOut()),
                           target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc].mom_.getVertexIn())); 
        bool ret = ( dis <= 0.0 );
        return ret;
        */
    }
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortGather,
                const ReallocatableArray<Ttc> & tc_first,
                const S32 adr_tc,
                //const F64ort & pos_target_box,
                const TargetBox<TSM> & target_box,
                const F64 r_crit_sq,
                const F64vec & len_peri,
                TagWalkModeNormal){
        return ( target_box.vertex_out_.overlapped(tc_first[adr_tc].mom_.getVertexIn())
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
        /*
        F64 dis = target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc].mom_.getVertexIn());
        bool ret = ( dis <= 0.0 );
        return ret;
        */
    }    
    // IS OPEN
    ///////////
    
    ///////////
    // COPY INFO DISTANT
    inline void CopyInfoDistant(TagForceShort,
                                const S32 adr_sp,
                                ReallocatableArray<S32> & adr_sp_list){
        // do nothing
    }
    inline void CopyInfoDistant(TagForceLong,
                                const S32 adr_sp,
                                ReallocatableArray<S32> & adr_sp_list){
        adr_sp_list.push_back(adr_sp);
    }
    // COPY INFO DISTANT
    ///////////

    ///////////
    // COPY INFO CLOSE
    template<class TSM, class Ttp, class Tep, class Tsp>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            }
            else{
                adr_sp_list.pushBackNoCheck(cnt_adr_ptcl);
            }
        }
    }
    
    template<class TSM, class Ttp, class Tep, class Tsp>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafFalse,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            }
            else{
                adr_sp_list.pushBackNoCheck(cnt_adr_ptcl);
            }
        }
    }

    
    template<class TSM, class Ttp, class Tep, class Tsp>
    inline void CopyInfoClose(TagForceShort,
                              TagChopLeafTrue,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++){
            const S32 adr = adr_ptcl + ip;
            if( target_box.isInEpList(ep_first, adr) ){
                adr_ep_list.pushBackNoCheck(adr);
            }
        }
    }
    template<class TSM, class Ttp, class Tep, class Tsp>
    inline void CopyInfoClose(TagForceShort,
                              TagChopLeafFalse,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++){
            const S32 adr = adr_ptcl + ip;
            adr_ep_list.pushBackNoCheck(adr);
            //if( target_box.isInEpList(ep_first, adr) ) adr_ep_list.pushBackNoCheck(adr);
        }

    }
    // COPY INFO CLOSE
    ///////////
    
    /////////////////
    // GET OPEN BITS
    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLong,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis = target_box.vertex_.getDistanceMinSq(pos);
            open_bit |= (dis <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLongScatter,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            //const F64 dis0 = target_box.vertex_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut());
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis1 = target_box.vertex_.getDistanceMinSq(pos);
            open_bit |= ( target_box.vertex_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()) || dis1 <= r_crit_sq) << i;
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLongSymmetry,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          TagWalkModeNormal){
        U32 open_bit = 0;
#if 1
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis1 = target_box.vertex_in_.getDistanceMinSq(pos);
            open_bit |= ( (target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut())
                           || target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) )
                          || dis1 <= r_crit_sq) << i;
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
            /*
            const F64 dis0 = std::min( target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut()),
                                       target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexIn()) );
            if(Comm::getRank()==0){
                if((target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()))){
                    std::cerr<<"A) dis0= "<<target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut())<<std::endl;
                }
                if((target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()))){
                    std::cerr<<"B) dis1= "<<target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexIn())<<std::endl;
                }
            }
            */
#if 0
            if(Comm::getRank()==0){
                if(target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut()) <= 0.0){
                    
                    if(target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()) == 0){
                        std::cout<<"A) n_ptcl= "<<tc_first[adr_tc+i].n_ptcl_
                                 <<" adr_tc+i= "<<adr_tc+i
                                 <<" tc_first[adr_tc+i].mom_.getVertexOut()= "<<tc_first[adr_tc+i].mom_.getVertexOut()
                                 <<" tc_first[adr_tc+i].mom_.getVertexIn()= "<<tc_first[adr_tc+i].mom_.getVertexIn()
                                 <<std::endl;
                    }
                        /*
                    std::cerr<<"A) over= "<<target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut())
                        //<<" dis= "<<target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut())
                        //<<" target_box.vertex_in_= "<<target_box.vertex_in_
                        //<<" tc_first[adr_tc+i].mom_.getVertexOut()= "<<tc_first[adr_tc+i].mom_.getVertexOut()
                             <<" n_ptcl= "<<tc_first[adr_tc+i].n_ptcl_
                             <<std::endl;
                        */
                }
                /*
                if(target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexIn()) <= 0.0){
                    if(target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) == 0){
                        std::cerr<<"B) n_ptcl= "<<tc_first[adr_tc+i].n_ptcl_
                                 <<" tc_first[adr_tc+i].mom_.getVertexOut()= "<<tc_first[adr_tc+i].mom_.getVertexOut()
                                 <<std::endl;
                    }
                }
                */
            }
#endif
        }

        /*
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64 dis0 = std::min( target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut()),
                                       target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexIn()) );
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis1 = target_box.vertex_in_.getDistanceMinSq(pos);
            open_bit |= ( dis0 <= 0.0 || dis1 <= r_crit_sq) << i;
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
        }
        */
        /*
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis = target_box.vertex_in_.getDistanceMinSq(pos);
            open_bit |= (dis <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
        */
#else
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis = target_box.vertex_in_.getDistanceMinSq(pos);
            open_bit |= (dis <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
#endif
        return open_bit;
    }
    
    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLongCutoff,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        if (r_crit_sq >= 0.0) {
            // theta_ > 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
                const F64 dis = target_box.vertex_.getDistanceMinSq(pos);
                open_bit |= (dis <= r_crit_sq) << i;
                open_bit |= ( ( (target_box.vertex_.overlapped( tc_first[adr_tc+i].mom_.getVertexOut()))
                                && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                              << (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        else {
            // theta_ = 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                open_bit |= 1 << i;
                open_bit |= ( ( (target_box.vertex_.overlapped( tc_first[adr_tc+i].mom_.getVertexOut()))
                                && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                              << (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        return open_bit;
    }
    
    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchShortScatter,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          //const F64ort & pos_target_box,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
#if 0
            const F64 dis = target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut());
            open_bit |= ( dis <= 0.0 ) << i;
#else
            open_bit |= ( target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()) ) << i;
#endif
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchShortSymmetry,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          //const F64ort & pos_target_box,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= ( target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut())
                          || target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) ) << i;            
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchShortGather,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          //const F64ort & pos_target_box,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= ( target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) ) << i;
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }
    // GET OPEN BITS
    /////////////////

    // for both long mode and short mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Twalkmode, class Tchopleaf>
    inline void MakeListUsingTreeRecursive
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     //const F64ort & pos_target_box,
     const TargetBox<TSM> & target_box,
     const F64 r_crit_sq,
     const S32 n_leaf_limit,
     const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
     const F64vec & len_peri){
        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl_child = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        if( !(tc_cur->isLeaf(n_leaf_limit)) ){ // not leaf
            U32 open_bit = GetOpenBit(typename TSM::search_type(), tc_first, adr_tc_child,
                                      target_box, r_crit_sq*0.25, len_peri,
                                      typename Twalkmode::walk_type());
            for(S32 i=0; i<N_CHILDREN; i++){
                if( !((open_bit>>(i+N_CHILDREN)) & 0x1) ) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp, Twalkmode, Tchopleaf>
                        (tc_first, adr_tc_child+i, tp_first, ep_first, adr_ep_list, sp_first, adr_sp_list,
                         target_box, r_crit_sq*0.25, n_leaf_limit, adr_tree_sp_first, len_peri);
                }
                else{ // far
                    CopyInfoDistant(typename TSM::force_type(), adr_tc_child+adr_tree_sp_first+i, adr_sp_list);
                }
            }
        }
        else{ //leaf
#if 0
            CopyInfoClose(typename TSM::force_type(), tp_first, adr_ptcl_child, n_ptcl,
                          ep_first, adr_ep_list, sp_first, adr_sp_list);
#elif 0
            adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
            for(S32 ip=0; ip<n_ptcl; ip++){
                F64 dis_sq = target_box.vertex_in_.getDistanceMinSQ(ep_first[adr_ptcl_child+ip].pos);
                if(dis_sq > ep_first[adr_ptcl_child+ip].getRSearch()*ep_first[adr_ptcl_child+ip].getRSearch()) continue;
                adr_ep_list.pushBackNoCheck(adr_ptcl_child+ip);
            }
#else
            CopyInfoClose(typename TSM::force_type(), Tchopleaf(), tp_first, adr_ptcl_child, n_ptcl,
                          ep_first, adr_ep_list, sp_first, adr_sp_list, target_box);
#endif
        }
    }

    // for long mode 
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Twalkmode, class Tchopleaf>
    inline void MakeListUsingTreeRecursiveTop
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     const TargetBox<TSM> & target_box,
     const F64 r_crit_sq,
     const S32 n_leaf_limit,
     const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
     const F64vec & len_peri){
        if ((r_crit_sq >= 0.0) || (typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF))) {
            // theta_ > 0 case or PS::SEARCH_MODE_LONG_CUTOFF
            if( IsOpen(typename TSM::search_type(), tc_first, adr_tc, target_box,
                       r_crit_sq, len_peri, typename Twalkmode::walk_type()) ){
                MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp, Twalkmode, Tchopleaf>
                    (tc_first, adr_tc, tp_first, ep_first, adr_ep_list, sp_first,
                     adr_sp_list, target_box, r_crit_sq, n_leaf_limit,
                     adr_tree_sp_first, len_peri);
            } 
        }
        else {
            // theta_ = 0 case with the modes other than PS::SEARCH_MODE_LONG_CUTOFF
            const S32 n_ptcl = tc_first[0].n_ptcl_;
            S32 adr_ptcl = tc_first[0].adr_ptcl_;
            CopyInfoClose(typename TSM::force_type(), Tchopleaf(), tp_first, adr_ptcl, n_ptcl,
                          ep_first, adr_ep_list, sp_first, adr_sp_list, target_box);
        }
    }

    // for short mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Twalkmode, class Tchopleaf>
    inline void MakeListUsingTreeRecursiveTop
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const TargetBox<TSM> & target_box,
     const F64 r_crit_sq,
     const S32 n_leaf_limit,
     const F64vec & len_peri){
        S32 adr_tree_sp_first = 0;
        const ReallocatableArray<SuperParticleBase> sp_first;
        ReallocatableArray<S32> adr_sp_list;
        MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, SuperParticleBase, Twalkmode, Tchopleaf>
            (tc_first, adr_tc, tp_first, ep_first, adr_ep_list, sp_first,
             adr_sp_list, target_box, r_crit_sq, n_leaf_limit,
             adr_tree_sp_first, len_peri);
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
            open_bit |= (tc_first[adr_tc+i].mom_.getVertexOut().overlapped(pos_box) << i);
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
            is_overlapped = tc_first[0].mom_.getVertexOut().overlapped(pos_box);
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
            open_bit |= (tc_first_B[adr_tc_B+i].mom_.getVertexOut().overlapped(pos_domain)
                         || IsOverlappedTop(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexIn(), n_leaf_limit_A)) << i;
            open_bit |= (tc_first_B[adr_tc_B+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if false, just skip
        }
        /*
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= (tc_first_B[adr_tc_B+i].mom_.getVertexOut().overlapped(pos_domain)) << i;
            if( !(tc_first_B[adr_tc_B+i].mom_.getVertexOut().overlapped(pos_domain)) ){
                open_bit |= IsOverlappedTop(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexIn(), n_leaf_limit_A) << i;
            }
            open_bit |= (tc_first_B[adr_tc_B+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if false, just skip
        }
        */
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
                    if( pos_domain.overlapped(tc_first_B[adr_tc_child].mom_.getVertexOut()) ){
                        //if(Comm::getRank()==0) std::cerr<<"tc_first_B[adr_tc_child].mom_.getVertexOut()= "<<tc_first_B[adr_tc_child].mom_.getVertexOut()<<std::endl;
                        continue;
                    }
                    else{
                        /*
                        if(Comm::getRank()==0){
                            std::cerr<<"pos_domain= "<<pos_domain
                                     <<" tc_first_B[adr_tc_child].mom_.getVertexOut()= "
                                     <<tc_first_B[adr_tc_child].mom_.getVertexOut()
                                     <<" tc_first_B[adr_tc_child].mom_.getVertexIn()= "
                                     <<tc_first_B[adr_tc_child].mom_.getVertexIn()
                                     <<" overlapped= "
                                     <<IsOverlappedTop(tc_first_A, tc_first_B[adr_tc_child].mom_.getVertexIn(), n_leaf_limit_A)
                                     <<std::endl;
                        }
                        */                       
                        S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                        for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                            /*
                            if(Comm::getRank()==0){
                                std::cerr<<"adr_ptcl_tmp= "<<adr_ptcl_tmp<<std::endl;
                            }
                            */
                            adr_ptcl_send.push_back(adr_ptcl_tmp);
                        }
                    }
                    /*
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                        if( pos_domain.overlapped(tc_first_B[adr_tc_child].mom_.getVertexOut()) ) continue;
                        adr_ptcl_send.push_back(adr_ptcl_tmp);
                    }
                    */
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
#if 0
            S32 adr_ptcl_tmp = tc_first_B[0].adr_ptcl_;
            const S32 n_child = tc_first_B[0].n_ptcl_;
#else
            // original
            if( tc_first_B[0].mom_.getVertexOut().overlapped(pos_domain) ){
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
#endif
        }
    }
}
