#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{

    ////////////
    // walk mode
    struct TagChopLeafTrue{};
    struct TagChopLeafFalse{};
    // walk mode
    ////////////

    ////////////
    // copy info close mode
    struct TagCopyInfoCloseNormal{};
    struct TagCopyInfoCloseNoSp{};
    struct TagCopyInfoCloseWithTpAdrptcl{};
    // copy info close mode
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
        template<class Tipg>
        void set(const Tipg & ipg){
            vertex_ = ipg.vertex_in_;
        }
        template<typename T>
        void setForExLet(const T & top, const F64vec & s){
            vertex_ = top.vertex_in_.shift(s);
        }
        void dump(std::ostream & fout=std::cerr){
            fout<<"vertex_= "<<vertex_<<std::endl;
        }
    };
  
    template<>
    struct TargetBox<SEARCH_MODE_LONG_SYMMETRY>{
        F64ort vertex_in_;
        F64ort vertex_out_;
        template<class Tipg>
        void set(const Tipg & ipg){
            vertex_in_  = ipg.vertex_in_;
            vertex_out_ = ipg.vertex_out_;
        }
        template<typename T>
        void setForExLet(const T & top, const F64vec & s){
            vertex_in_  = top.vertex_in_.shift(s);
            vertex_out_ = top.vertex_out_.shift(s);
        }
    };

    template<>
    struct TargetBox<SEARCH_MODE_LONG_SCATTER>{
        F64ort vertex_;
        template<class Tipg>
        void set(const Tipg & ipg){
            vertex_ = ipg.vertex_in_;
        }
        template<typename T>
        void setForExLet(const T & top, const F64vec & s){
            vertex_ = top.vertex_in_.shift(s);
        }
        void dump(std::ostream & fout=std::cerr){
            fout<<"vertex_= "<<vertex_<<std::endl;
        }
    };
    
    template<>
    struct TargetBox<SEARCH_MODE_LONG_CUTOFF>{
        F64ort vertex_;
        template<class Tipg>
        void set(const Tipg & ipg){
            vertex_ = ipg.vertex_in_;
        }
        template<typename T>
        void setForExLet(const T & top, const F64vec & s){
            vertex_ = top.vertex_in_.shift(s);
        }
        void dump(std::ostream & fout=std::cerr){
            fout<<"vertex_= "<<vertex_<<std::endl;
        }
    };
    
    template<>
    struct TargetBox<SEARCH_MODE_SCATTER>{
        F64ort vertex_in_;
        template<class Tipg>
        void set(const Tipg & ipg){
            vertex_in_ = ipg.vertex_in_;
        }
        template<class Tep, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr,
                               const F64vec & len_peri) const {
            //const auto dis_sq = vertex_in_.getDistanceMinSQ(ep_first[adr].getPos());
            const F64vec ep_pos = ep_first[adr].getPos();
            //const auto dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_in_, ep_pos, len_peri);
            const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_in_, ep_pos, len_peri);
            const F64 r_crit_sq = ep_first[adr].getRSearch() * ep_first[adr].getRSearch();
            return dis_sq <= r_crit_sq;
        }
    };
    
    template<>
    struct TargetBox<SEARCH_MODE_SYMMETRY>{
        F64ort vertex_out_;
        F64ort vertex_in_;
        template<class Tipg>
        void set(const Tipg & ipg){
            vertex_out_ = ipg.vertex_out_;
            vertex_in_ = ipg.vertex_in_;
        }
        template<class Tep, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr,
                               const F64vec & len_peri) const {
            //const F64 dis_sq = vertex_in_.getDistanceMinSQ(ep_first[adr].getPos());
	    const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_in_, ep_first[adr].getPos(), len_peri);
            const F64 r_crit_sq = ep_first[adr].getRSearch() * ep_first[adr].getRSearch();
#if 1
            //const F64 dis_sq_2 = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_out_, ep_first[adr].getPos(), len_peri);
            return (dis_sq <= r_crit_sq) || IsOverlapped<CALC_DISTANCE_TYPE>(vertex_out_, ep_first[adr].getPos(), len_peri);
#else
            return (dis_sq <= r_crit_sq) || (vertex_out_.overlapped(ep_first[adr].getPos()));
#endif
        }
    };
  
    template<>
    struct TargetBox<SEARCH_MODE_GATHER>{
        F64ort vertex_out_;
        template<class Tipg>
        void set(const Tipg & ipg){
            vertex_out_ = ipg.vertex_out_;
        }
        template<class Tep, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr,
                               const F64vec & len_peri) const {
            return IsOverlapped<CALC_DISTANCE_TYPE>(vertex_out_, ep_first[adr].getPos(), len_peri);
	    //const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_out_, ep_first[adr].getPos(), len_peri);
            //return dis_sq;
            //return vertex_out_.overlapped(ep_first[adr].getPos());
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

    
    ///////////
    // IS OPEN
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchLong,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }
    
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongScatter,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }

    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongSymmetry,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }

    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongCutoff,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
#if 1
        return ( IsOverlapped<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc].geo_.getVertexOut(), len_peri)
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#else
        return ( (target_box.vertex_.overlapped(tc_first[adr_tc].geo_.getVertexOut()))
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
#endif
    }
    
    template<enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortScatter,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> target_box,
                       const F64vec & len_peri){
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
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
            adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
        }
    }
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
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
            adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
        }
    }

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
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            //const F64 dis_sq = target_box.vertex_.getDistanceMinSq(pos);
	    const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, pos, len_peri);
            const F64 size = tc_first[adr_tc+i].geo_.getSize();
            const F64 r_crit_sq = size*size*inv_theta_sq;
            open_bit |= (dis_sq <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBitOld(TagSearchLongScatter,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
	    const F64 size = tc_first[adr_tc+i].geo_.getSize();
            const F64 r_crit_sq = size*size*inv_theta_sq;

	    const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, pos, len_peri);
	    open_bit |= ( GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri) <= 0.0 || dis_sq <= r_crit_sq) << i;

	    
            //const F64 dis_sq = target_box.vertex_.getDistanceMinSq(pos);
            //open_bit |= ( target_box.vertex_.overlapped(tc_first[adr_tc+i].geo_.getVertexOut()) || dis_sq <= r_crit_sq) << i;
	    //if(Comm::getRank()==0){
	    //  std::cerr<<"target_box.vertex_= "<<target_box.vertex_<<" tc_first[adr_tc+i].geo_.getVertexOut()= "<<tc_first[adr_tc+i].geo_.getVertexOut()
	    //	       <<" dis= "<<GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri)
	    //	       <<std::endl;
	    //}
	    

            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
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
        U32 open_bit = 0;
	//	std::cerr << "GetOpenBit long scatter  called\n";
        for(S32 i=0; i<N_CHILDREN; i++){
	    if(  tc_first[adr_tc+i].n_ptcl_ >0){
		const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
		const F64 size = tc_first[adr_tc+i].geo_.getSize();
		const F64 r_crit_sq = size*size*inv_theta_sq;
		
		const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_, pos, len_peri);
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
    
    template<class TSM, class Ttc, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline U32 GetOpenBit(TagSearchShortScatter,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64vec & len_peri,
                          const F64 inv_theta_sq){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
#if 1
            const auto dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(target_box.vertex_in_, tc_first[adr_tc+i].geo_.getVertexOut(), len_peri);
            open_bit |= ( dis_sq <= 0.0 ) << i;
#else
            //open_bit |= ( target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()) ) << i;
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
      //if(Comm::getRank()==0){
      //    std::cerr<<"B) CALC_DISTANCE_TYPE= "<<CALC_DISTANCE_TYPE<<std::endl;
      //}
        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl_child = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        const F64 inv_theta_sq = (theta > 0.0) ? 1.0 / (theta*theta) : -1.0;
        if( !(tc_cur->isLeaf(n_leaf_limit)) ){ // not leaf
	  U32 open_bit = GetOpenBit<TSM, Ttc, CALC_DISTANCE_TYPE>(typename TSM::search_type(), tc_first, adr_tc_child,
						target_box, len_peri, inv_theta_sq);
            for(S32 i=0; i<N_CHILDREN; i++){
                if( !((open_bit>>(i+N_CHILDREN)) & 0x1) ) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
		    MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp, Tchopleaf, Tcopyinfoclose, CALC_DISTANCE_TYPE>
                        (tc_first, adr_tc_child+i, tp_first, ep_first, adr_ep_list, sp_first, adr_sp_list,
                         target_box, n_leaf_limit, adr_tree_sp_first, len_peri, theta);
                }
                else{ // far
                    CopyInfoDistant(typename TSM::force_type(), adr_tc_child+adr_tree_sp_first+i, adr_sp_list);
                }
            }
        }
        else{ //leaf
            CopyInfoClose<TSM, Ttp, Tep, Tsp, CALC_DISTANCE_TYPE>
                (typename TSM::force_type(), Tchopleaf(), Tcopyinfoclose(), tp_first, adr_ptcl_child, n_ptcl,
                 ep_first, adr_ep_list, sp_first, adr_sp_list, target_box, len_peri);
        }
    }

    
    // new
    // for long mode
    //template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tchopleaf, class Tcopyinfoclose, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE=CALC_DISTANCE_TYPE_NORMAL>
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
    template<class TSM, class Ttc, class Ttp, class Tep, class Tchopleaf, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline void MakeListUsingTreeRecursiveTop
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const TargetBox<TSM> & target_box,
     const S32 n_leaf_limit,
     const F64vec & len_peri){
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
