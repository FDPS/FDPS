#pragma once
#include<key.hpp>

namespace ParticleSimulator{

    class TreeParticle{
    public:
        U64 key_;
        U32 adr_ptcl_; // U32 because its MSB is used for distinguish if it is EP or SP
        TreeParticle() : key_(SetMSB( U64(0) )), adr_ptcl_(SetMSB( U32(0) )) {}

        template<class Tfp>
        void setFromFP(const Tfp & fp, const U32 adr){
            //key_ = MortonKey<DIMENSION>::getKey( fp.getPos() );
            key_ = MortonKey::getKey( fp.getPos() );
            adr_ptcl_ = adr;
        }
        template<class Tep>
        void setFromEP(const Tep & ep, const U32 adr){
            //key_ = MortonKey<DIMENSION>::getKey( ep.getPos() );
            key_ = MortonKey::getKey( ep.getPos() );
            adr_ptcl_ = adr;
        }
        template<class Tsp>
        void setFromSP(const Tsp & sp, const U32 adr){
            //key_ = MortonKey<DIMENSION>::getKey( sp.getPos() );
            key_ = MortonKey::getKey( sp.getPos() );
            adr_ptcl_ = SetMSB(adr);
        }

        U64 getKey() const {
            return key_;
        }
        // for DEBUG
        void setKey(const U64 key){
            key_ = key;
        }
        void dump(std::ostream & fout=std::cout) const {
            fout<<std::oct<<"key_="<<key_<<std::dec<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
        }
    };

    template<class Tmom>
    class TreeCell{
    public:
        S32 n_ptcl_;
        U32 adr_tc_;
        U32 adr_ptcl_;
        S32 level_;
        Tmom mom_;
        TreeCell(){
            n_ptcl_ = adr_tc_ = adr_ptcl_ = level_ = 0;
            mom_.init();
        }
        void clear(){
            n_ptcl_ = adr_tc_ = adr_ptcl_ = level_ = 0;
            mom_.init();
        }
        void clearMoment(){
            mom_.init();
        }
        bool isLeaf(const S32 n_leaf_limit) const {
            return ( n_ptcl_ <= n_leaf_limit || level_ == TREE_LEVEL_LIMIT);
        }
        
        // for DEBUG
        void dump(std::ostream & fout=std::cout) const {
            fout<<"n_ptcl_="<<n_ptcl_<<std::endl;
            fout<<"adr_tc_="<<adr_tc_<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
            mom_.dump(fout);
        }
        
        template<class Tep>
        void dumpTree(const Tep * const first_ep_ptr,
                      const TreeCell * const first_tc_ptr,
                      const F32vec & center,
                      const F32 half_length,
                      const S32 n_leaf_limit,
                      std::ostream & fout=std::cout) const {
            fout<<std::endl;
            fout<<"cell info"<<std::endl;
            dump(fout);
            if( !(this->isLeaf(n_leaf_limit)) ){
                fout<<"this cell is not a leaf"<<std::endl;
                fout<<"half_length="<<half_length<<std::endl;
                fout<<"center="<<center<<std::endl;
                fout<<"level="<<this->level_<<std::endl;
                const TreeCell * child = first_tc_ptr + adr_tc_;
                for(S32 ic=0; ic<N_CHILDREN; ic++){
                    if((child + ic)->n_ptcl_ <= 0) continue;
                    const Tep * ptcl = first_ep_ptr + adr_ptcl_;
                    for(S32 ip=0; ip<n_ptcl_; ip++, ptcl++){
                        if(!IsInBox(ptcl->getPos(), center, half_length)){
                            fout<<"out of box(Cell)"<<std::endl;
                            fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                            fout<<"center="<<center<<std::endl;
                            fout<<"half_length="<<half_length<<std::endl;
                        }
                        else{
                            fout<<"in box(Cell)"<<std::endl;
                            fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                            fout<<"center="<<center<<std::endl;
                            fout<<"half_length="<<half_length<<std::endl;
                        }
                    }
                    fout<<"octid="<<ic<<std::endl;
                    (child + ic)->dumpTree
                        (first_ep_ptr, first_tc_ptr, 
                         center+SHIFT_CENTER[ic]*half_length, half_length*0.5, 
                         n_leaf_limit, fout);
                }
            }
            else{
                fout<<"this cell is a leaf"<<std::endl;
                fout<<"half_length="<<half_length<<std::endl;
                fout<<"center="<<center<<std::endl;
                fout<<"level="<<this->level_<<std::endl;
                const Tep * ptcl = first_ep_ptr + adr_ptcl_;
                for(S32 ip=0; ip<n_ptcl_; ip++, ptcl++){
                    if(!IsInBox(ptcl->getPos(), center, half_length)){
                        fout<<"out of box(LeafCell)"<<std::endl;
                        fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                        fout<<"center="<<center<<std::endl;
                        fout<<"half_length="<<half_length<<std::endl;
                    }
                    else{
                        fout<<"in box(LeafCell)"<<std::endl;
                        fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                        fout<<"center="<<center<<std::endl;
                        fout<<"half_length="<<half_length<<std::endl;
                    }
                    //ptcl->dump(fout);
                    fout<<std::endl;
                }
                fout<<std::endl;
            }
        }
        
        template<class Tep>
        void checkTree(const Tep * const first_ep_ptr,
                       const TreeCell * const first_tc_ptr,
                       const F32vec & center,
                       const F32 half_length,
                       const S32 n_leaf_limit,
                       const F32 tolerance,
                       S32 & err,
                       std::ostream & fout=std::cout) const {
            if( !(this->isLeaf(n_leaf_limit)) ){
                //std::cerr<<"adr_tc_="<<adr_tc_<<std::endl;
                const TreeCell * child = first_tc_ptr + adr_tc_;
                for(S32 ic=0; ic<N_CHILDREN; ic++){
                    //std::cerr<<"(child + ic)->n_ptcl_="<<(child + ic)->n_ptcl_<<std::endl;
                    if((child + ic)->n_ptcl_ <= 0) continue;
                    const Tep * ptcl = first_ep_ptr + adr_ptcl_;
                    for(S32 ip=0; ip<n_ptcl_; ip++, ptcl++){
                        if(!IsInBox(ptcl->getPos(), center, half_length, tolerance)){
                            fout<<"out of box(Cell)"<<std::endl;
                            fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                            fout<<"center="<<center<<std::endl;
                            fout<<"half_length="<<half_length<<std::endl;
                            fout<<"(center+half_length)-ptcl->getPos()="<<(center+half_length)-ptcl->getPos()<<std::endl;
                            fout<<"(center-half_length)-ptcl->getPos()="<<(center-half_length)-ptcl->getPos()<<std::endl;
                            err++;
                        }
                    }
                    //std::cerr<<"ic="<<ic<<std::endl;
                    //std::cerr<<"SHIFT_CENTER[ic]="<<SHIFT_CENTER[ic]<<std::endl;
                    //std::cerr<<"center+SHIFT_CENTER[ic]*half_length="<<center+SHIFT_CENTER[ic]*half_length<<std::endl;
                    (child + ic)->checkTree
                        (first_ep_ptr, first_tc_ptr,
                         center+SHIFT_CENTER[ic]*half_length, half_length*0.5,
                         n_leaf_limit, tolerance, err, fout);
                }
            }
            else{
                const Tep * ptcl = first_ep_ptr + adr_ptcl_;
                for(S32 ip=0; ip<n_ptcl_; ip++, ptcl++){
                    if(!IsInBox(ptcl->getPos(), center, half_length, tolerance)){
                        fout<<"out of box(Leaf)"<<std::endl;
                        fout<<"center="<<center<<std::endl;
                        fout<<"half_length="<<half_length<<std::endl;
                        fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                        fout<<"(center+half_length)-ptcl->getPos()="<<(center+half_length)-ptcl->getPos()<<std::endl;
                        fout<<"(center-half_length)-ptcl->getPos()="<<(center-half_length)-ptcl->getPos()<<std::endl;
                        err++;
                    }
                }
            }
        }
        
        template<class Tep, class Tsp>
        void checkTreeLongGlobalTree(const Tep * const first_ep_ptr,
                                     const Tsp * const first_sp_ptr,
                                     const TreeParticle * const first_tp_ptr,
                                     const TreeCell * const first_tc_ptr,
                                     const F32vec & center,
                                     const F32 half_length,
                                     const S32 n_leaf_limit,
                                     const F32 tolerance,
                                     S32 & err,
                                     std::ostream & fout=std::cout) const {
            if( !(this->isLeaf(n_leaf_limit)) ){
                const TreeCell * child = first_tc_ptr + adr_tc_;
                for(S32 ic=0; ic<N_CHILDREN; ic++){
                    if((child + ic)->n_ptcl_ <= 0) continue;
                    const TreeParticle * tp_tmp = first_tp_ptr + adr_ptcl_;
                    for(S32 ip=0; ip<n_ptcl_; ip++, tp_tmp++){
                        F32vec pos_tmp;
                        /*
                          const U32 adr = tp_tmp->adr_ptcl_;
                          if( GetMSB(adr) ) pos_tmp = first_sp_ptr[ClearMSB(adr)].getPos();
                          else pos_tmp = first_ep_ptr[adr].getPos();
                        */
                        if(GetMSB(tp_tmp->adr_ptcl_)) pos_tmp = first_sp_ptr[adr_ptcl_+ip].getPos();
                        else pos_tmp = first_ep_ptr[adr_ptcl_+ip].getPos();
                        //if(Comm::getRank() == 0){
                        if(!IsInBox(pos_tmp, center, half_length, tolerance)){
                            fout<<"out of box(Cell)"<<std::endl;
                            fout<<"pos_tmp="<<pos_tmp<<std::endl;
                            fout<<"center="<<center<<std::endl;
                            fout<<"half_length="<<half_length<<std::endl;
                            //fout<<"adr="<<adr<<std::endl;
                            fout<<"adr_ptcl_+ip="<<adr_ptcl_+ip<<std::endl;
                            fout<<"(center+half_length)-pos_tmp="<<(center+half_length)-pos_tmp<<std::endl;
                            fout<<"(center-half_length)-pos_tmp="<<(center-half_length)-pos_tmp<<std::endl;
                            err++;
                            //}
                        }
                    }
                    /*
                    (child + ic)->checkTreeLongGlobalTree < Tep, Tsp >
                        (first_ep_ptr, first_sp_ptr, first_tp_ptr, first_tc_ptr,
                         center+SHIFT_CENTER[ic]*half_length, half_length*0.5,
                         n_leaf_limit, tolerance, err, fout);
                    */
                    (child + ic)->checkTreeLongGlobalTree
                        (first_ep_ptr, first_sp_ptr, first_tp_ptr, first_tc_ptr,
                         center+SHIFT_CENTER[ic]*half_length, half_length*0.5,
                         n_leaf_limit, tolerance, err, fout);
                }
            }
            else{
                const TreeParticle * tp_tmp = first_tp_ptr + adr_ptcl_;
                for(S32 ip=0; ip<n_ptcl_; ip++, tp_tmp++){
                    F32vec pos_tmp;
                    if(GetMSB(tp_tmp->adr_ptcl_)) pos_tmp = first_sp_ptr[adr_ptcl_+ip].getPos();
                    else pos_tmp = first_ep_ptr[adr_ptcl_+ip].getPos();
                    //if(Comm::getRank() == 0){
                    if(!IsInBox(pos_tmp, center, half_length, tolerance)){
                        fout<<"out of box(Leaf)"<<std::endl;
                        fout<<"center="<<center<<std::endl;
                        fout<<"half_length="<<half_length<<std::endl;
                        //fout<<"adr="<<adr<<std::endl;
                        fout<<"adr_ptcl_+ip="<<adr_ptcl_+ip<<std::endl;
                        fout<<"pos_tmp="<<pos_tmp<<std::endl;
                        fout<<"(center+half_length)-pos_tmp="<<(center+half_length)-pos_tmp<<std::endl;
                        fout<<"(center-half_length)-pos_tmp="<<(center-half_length)-pos_tmp<<std::endl;
                        err++;
                    }
                    // }
                }
            }
        }
    };

#if 1
    template<class Tipg>
    class IPGroup{};
    
    template<>
    class IPGroup<TagIpgLongNormal>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_in_;
        template<class Ttc>
        void copyFromTC(const Ttc & tc){
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
        }
    };

    template<>
    class IPGroup<TagIpgIn>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_in_;
        template<class Ttc>
        void copyFromTC(const Ttc & tc){
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_in_ = tc.mom_.vertex_in_;
        }
    };
    
    template<>
    class IPGroup<TagIpgInAndOut>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_in_;
        F64ort vertex_out_;
        template<class Ttc>
        void copyFromTC(const Ttc & tc){
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_in_ = tc.mom_.vertex_in_;
            vertex_out_ = tc.mom_.vertex_out_;
        }
    };
    
    template<>
    class IPGroup<TagIpgOut>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_out_;
        template<class Ttc>
        void copyFromTC(const Ttc & tc){
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_out_ = tc.mom_.vertex_out_;
        }
    };
#else
    template<class TSM>
    class IPGroup{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_;
        template<class Ttc> 
        void copyFromTC(const Ttc & tc){
            CopyFromTCDummy<TSM, Ttc>()(this, tc);
        }
        //
        template<class TSM2, class Ttc2, class Tdummy=void>
        struct CopyFromTCDummy{
            void operator () (IPGroup * ip, const Ttc2 & tc){
                // for SCATTER, FIXED
                ip->n_ptcl_ = tc.n_ptcl_;
                ip->adr_ptcl_ = tc.adr_ptcl_;
                ip->vertex_ = tc.mom_.vertex_in_;
            }
        };
        template<class Ttc2, class Tdummy>
        struct CopyFromTCDummy<SEARCH_MODE_GATHER, Ttc2, Tdummy>{
            void operator () (IPGroup * ip, const Ttc2 & tc){
                // for GAHTER
                ip->n_ptcl_ = tc.n_ptcl_;
                ip->adr_ptcl_ = tc.adr_ptcl_;
                ip->vertex_ = tc.mom_.vertex_out_;
            }
        };

        template<class Ttc2, class Tdummy>
        struct CopyFromTCDummy<SEARCH_MODE_LONG, Ttc2, Tdummy>{
            void operator () (IPGroup * ip, const Ttc2 & tc){
                ip->n_ptcl_ = tc.n_ptcl_;
                ip->adr_ptcl_ = tc.adr_ptcl_;
            }
        };
        template<class Ttc2, class Tdummy>
        struct CopyFromTCDummy<SEARCH_MODE_LONG_CUTOFF, Ttc2, Tdummy>{
            void operator () (IPGroup * ip, const Ttc2 & tc){
                ip->n_ptcl_ = tc.n_ptcl_;
                ip->adr_ptcl_ = tc.adr_ptcl_;
            }
        };
        // for DEBUG
        void dump(std::ostream & fout = std::cout){
            fout<<"n_ptcl_="<<n_ptcl_<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
            fout<<"vertex_="<<vertex_<<std::endl;
        }
    };

    template<>
    class IPGroup<SEARCH_MODE_SYMMETRY>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_;
        F64ort vertex_in_;
        template<class Ttc> 
        void copyFromTC(const Ttc & tc){
            // for SYMMETRY
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_ = tc.mom_.vertex_out_;
            vertex_in_ = tc.mom_.vertex_in_;
        }
        // for DEBUG
        void dump(std::ostream & fout = std::cout){
            fout<<"n_ptcl_="<<n_ptcl_<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
            fout<<"vertex_ing_="<<vertex_in_<<std::endl;
        }
    };

    template<>
    class IPGroup<SEARCH_MODE_LONG_SYMMETRY>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_;
        F64ort vertex_out_;
        template<class Ttc> 
        void copyFromTC(const Ttc & tc){
            // for SYMMETRY
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_ = tc.mom_.vertex_in_;
            vertex_out_ = tc.mom_.vertex_out_;
        }
        // for DEBUG
        void dump(std::ostream & fout = std::cout){
            fout<<"n_ptcl_="<<n_ptcl_<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
            fout<<"vertex_="<<vertex_<<std::endl;
            fout<<"vertex_out_="<<vertex_out_<<std::endl;
        }
    };    
#endif
    
    ////////////////////////
    // ESSENTIAL PARTICLE //
    class EPXROnly{
        F64 r_search;
        F64vec pos;
    public:
        F64 getRSearch() const { return r_search;}
        F64vec getPos() const { return pos;}
        void setPos(const F64vec & pos_in) { pos = pos_in;}
        template<class Tep>
        void copyFromEP(const Tep & ep){
            r_search = ep.getRSearch();
            pos = ep.getPos();
        }
        template<class Tptcl>
        const EPXROnly & operator = (const Tptcl & ptcl){
            r_search = ptcl.getRSearch();
            pos = ptcl.getPos();
            return (*this);
        }
    };


    //////////////
    /// Moment ///
    // new
    // these are the same as
    class MomentMonopoleInAndOut{
    public:
        FSP mass;
        FSPvec pos;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentMonopoleInAndOut(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentMonopoleInAndOut(const FSP m, const FSPvec & p){ 
            mass = m;
            pos = p;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopoleInAndOut & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentMonopoleInAndOut & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };
    typedef MomentMonopoleInAndOut MomentMonopoleSymmetry;


    class MomentQuadrupoleInAndOut{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentQuadrupoleInAndOut(){
            mass = 0.0;
            pos = 0.0;
	    quad = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentQuadrupoleInAndOut(const FSP m, const FSPvec & p, const FSPmat & q){ 
            mass = m;
            pos = p;
            quad = q;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quad.xx += cx * ptmp.x;
            this->quad.yy += cy * ptmp.y;
            this->quad.zz += cz * ptmp.z;
            this->quad.xy += cx * ptmp.y;
            this->quad.xz += cx * ptmp.z;
            this->quad.yz += cy * ptmp.z;
#else
	    // under construction
#endif
	}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentQuadrupoleInAndOut & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentQuadrupoleInAndOut & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 mtmp = mom.mass;
            F64vec ptmp = mom.pos - this->pos;
            F64 cx = mtmp * ptmp.x;
            F64 cy = mtmp * ptmp.y;
            F64 cz = mtmp * ptmp.z;
            this->quad.xx += cx * ptmp.x + mom.quad.xx;
            this->quad.yy += cy * ptmp.y + mom.quad.yy;
            this->quad.zz += cz * ptmp.z + mom.quad.zz;
            this->quad.xy += cx * ptmp.y + mom.quad.xy;
            this->quad.xz += cx * ptmp.z + mom.quad.xz;
            this->quad.yz += cy * ptmp.z + mom.quad.yz;
#else
	    // under construction
#endif
	}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };
    typedef MomentQuadrupoleInAndOut MomentQuadrupoleSymmetry;
    
    
    // for P^3T
    class MomentMonopoleScatter{
    public:
        FSP mass;
        FSPvec pos;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentMonopoleScatter(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentMonopoleScatter(const FSP m, const FSPvec & p){ 
            mass = m;
            pos = p;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopoleScatter & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentMonopoleScatter & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };


    class MomentQuadrupoleScatter{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentQuadrupoleScatter(){
            mass = 0.0;
            pos = 0.0;
	    quad = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentQuadrupoleScatter(const FSP m, const FSPvec & p, const FSPmat & q){ 
            mass = m;
            pos = p;
            quad = q;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quad.xx += cx * ptmp.x;
            this->quad.yy += cy * ptmp.y;
            this->quad.zz += cz * ptmp.z;
            this->quad.xy += cx * ptmp.y;
            this->quad.xz += cx * ptmp.z;
            this->quad.yz += cy * ptmp.z;
#else
	    // under construction
#endif
	}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentQuadrupoleScatter & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentQuadrupoleScatter & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 mtmp = mom.mass;
            F64vec ptmp = mom.pos - this->pos;
            F64 cx = mtmp * ptmp.x;
            F64 cy = mtmp * ptmp.y;
            F64 cz = mtmp * ptmp.z;
            this->quad.xx += cx * ptmp.x + mom.quad.xx;
            this->quad.yy += cy * ptmp.y + mom.quad.yy;
            this->quad.zz += cz * ptmp.z + mom.quad.zz;
            this->quad.xy += cx * ptmp.y + mom.quad.xy;
            this->quad.xz += cx * ptmp.z + mom.quad.xz;
            this->quad.yz += cy * ptmp.z + mom.quad.yz;
#else
	    // under construction
#endif
	}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };


    // for P^3T + PM
    class MomentMonopoleCutoffScatter{
    public:
        FSP mass;
        FSPvec pos;
        F64ort vertex_out_; // cutoff
        F64ort vertex_out2_; // search ep
        F64ort vertex_in_;
        MomentMonopoleCutoffScatter(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_out2_.init();
            vertex_in_.init();
        }
        MomentMonopoleCutoffScatter(const FSP m, const FSPvec & p){ 
            mass = m;
            pos = p;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexOut2() const { return vertex_out2_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_out2_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_out2_).merge(epj.getPos(), epj.getRSearch2());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopoleCutoffScatter & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_out2_).merge(mom.vertex_out2_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentMonopoleCutoffScatter & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentMonopole{
    public:
        FSP mass;
        FSPvec pos;
        MomentMonopole(){
            mass = 0.0;
            pos = 0.0;
        }
        MomentMonopole(const FSP m, const FSPvec & p){
            mass = m;
            pos = p;
        }
        void init(){
            mass = 0.0;
            pos = 0.0;
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            mass += epj.getCharge();
            pos += epj.getCharge() * epj.getPos();
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopole & mom){
            mass += mom.mass;
            pos += mom.mass * mom.pos;
        }
        void accumulate2(const MomentMonopole & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class MomentQuadrupole{
    public:
        FSPvec pos;
        FSP mass;
        FSPmat quad;
        void init(){
            pos = 0.0;
            mass = 0.0;
            quad = 0.0;
        }
        MomentQuadrupole(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
        }
        MomentQuadrupole(const FSP m, const FSPvec & p, const FSPmat & q){
            mass = m;
            pos = p;
            quad = q;
        }
        FSPvec getPos() const {
            return pos;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            mass += epj.getCharge();
            pos += epj.getCharge() * epj.getPos();
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quad.xx += cx * ptmp.x;
            this->quad.yy += cy * ptmp.y;
            this->quad.zz += cz * ptmp.z;
            this->quad.xy += cx * ptmp.y;
            this->quad.xz += cx * ptmp.z;
            this->quad.yz += cy * ptmp.z;
#else
	    // under construction
#endif
        }
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentQuadrupole & mom){
            mass += mom.mass;
            pos += mom.mass * mom.pos;
        }
        void accumulate2(const MomentQuadrupole & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 mtmp = mom.mass;
            F64vec ptmp = mom.pos - this->pos;
            F64 cx = mtmp * ptmp.x;
            F64 cy = mtmp * ptmp.y;
            F64 cz = mtmp * ptmp.z;
            this->quad.xx += cx * ptmp.x + mom.quad.xx;
            this->quad.yy += cy * ptmp.y + mom.quad.yy;
            this->quad.zz += cz * ptmp.z + mom.quad.zz;
            this->quad.xy += cx * ptmp.y + mom.quad.xy;
            this->quad.xz += cx * ptmp.z + mom.quad.xz;
            this->quad.yz += cy * ptmp.z + mom.quad.yz;
#else
	    // under construction
#endif
        }
    };

    class MomentMonopoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        MomentMonopoleGeometricCenter(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        MomentMonopoleGeometricCenter(const FSP c, const FSPvec & p, const SSP n){
            n_ptcl = n;
            charge = c;
            pos = p;
        }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            charge += epj.getCharge();
            pos += epj.getPos();
            n_ptcl++;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentMonopoleGeometricCenter & mom){
            charge += mom.charge;
            pos += mom.n_ptcl * mom.pos;
            n_ptcl += mom.n_ptcl;
        }
        void accumulate2(const MomentMonopoleGeometricCenter & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
        }
    };

    class MomentDipoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        MomentDipoleGeometricCenter(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
        }
        MomentDipoleGeometricCenter(const FSP c, const FSPvec & p, 
                                    const SSP n, const FSPvec & di){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
        }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            charge += epj.getCharge();
            pos += epj.getPos();
            n_ptcl++;
        }
        void accumulate(const MomentDipoleGeometricCenter & mom){
            charge += mom.charge;
            pos += mom.n_ptcl * mom.pos;
            n_ptcl += mom.n_ptcl;
        }
        void set(){
            pos = pos / n_ptcl;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
            this->dipole += epj.getCharge() * (epj.getPos() - this->pos);
        }
        void accumulate2(const MomentDipoleGeometricCenter & mom){
            //dipole += mom.charge * (mom.pos - this->pos);
            this->dipole += mom.charge * (mom.pos - this->pos) + mom.dipole;
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
        }
    };


    class MomentQuadrupoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        MomentQuadrupoleGeometricCenter(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
        }
        MomentQuadrupoleGeometricCenter(const FSP c, const FSPvec & p, 
                                        const SSP n, const FSPvec & di,
                                        const FSPmat & q){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
            quadrupole = q;
        }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
            quadrupole = 0.0;
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            charge += epj.getCharge();
            pos += epj.getPos();
            n_ptcl++;
        }
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentQuadrupoleGeometricCenter & mom){
            charge += mom.charge;
            pos += mom.n_ptcl * mom.pos;
            n_ptcl += mom.n_ptcl;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            this->dipole += ctmp * ptmp;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x;
            this->quadrupole.yy += cy * ptmp.y;
            this->quadrupole.zz += cz * ptmp.z;
            this->quadrupole.xy += cx * ptmp.y;
            this->quadrupole.xz += cx * ptmp.z;
            this->quadrupole.yz += cy * ptmp.z;
#else
	    // underconstruction
#endif
        }
        void accumulate2(const MomentQuadrupoleGeometricCenter & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = mom.charge;
            F64vec ptmp = mom.pos - this->pos;
            this->dipole += ctmp * ptmp + mom.dipole;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x + mom.quadrupole.xx;
            this->quadrupole.yy += cy * ptmp.y + mom.quadrupole.yy;
            this->quadrupole.zz += cz * ptmp.z + mom.quadrupole.zz;
            this->quadrupole.xy += cx * ptmp.y + mom.quadrupole.xy;
            this->quadrupole.xz += cx * ptmp.z + mom.quadrupole.xz;
            this->quadrupole.yz += cy * ptmp.z + mom.quadrupole.yz;
#else
	    // underconstruction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
            fout<<"quadrupole="<<quadrupole<<std::endl;
        }
    };
    
    class MomentMonopoleCutoff{
    public:
        FSP mass;
        FSPvec pos;
        F64ort vertex_out_;
        MomentMonopoleCutoff(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
        }
        MomentMonopoleCutoff(const FSP m,
			     const FSPvec & p,
			     const F64ort & v_out){
            mass = m;
            pos = p;
            vertex_out_ = v_out;
        }
        void init(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        F64ort getVertexOut() const { return vertex_out_; }

        template<class Tepj>
        void accumulateAtLeaf(Tepj & epj){
            mass += epj.getCharge();
            pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        }
        template<class Tepj>
        void accumulateAtLeaf2(Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopoleCutoff & mom){
            mass += mom.mass;
            pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
        }
        void accumulate2(const MomentMonopoleCutoff & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
	    fout<<"vertex_out_="<<vertex_out_<<std::endl;
        }
    };
    
    typedef MomentMonopoleCutoff MomentMonopolePeriodic;
    
    class MomentSearchInAndOut{
    public:
        F64ort vertex_out_;
        F64ort vertex_in_;
        void init(){
            vertex_out_.init();
            vertex_in_.init();
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        template<class Tep>
        void accumulateAtLeaf(const Tep & ep){
            (this->vertex_out_).merge(ep.getPos(), ep.getRSearch());
            (this->vertex_in_).merge(ep.getPos());
        }
        template<class Tep>
        void accumulateAtLeaf(Tep & ep){
            (this->vertex_out_).merge(ep.getPos(), ep.getRSearch());
            (this->vertex_in_).merge(ep.getPos());
        }
        void set(){}
        void accumulate(const MomentSearchInAndOut & _mom){
            (this->vertex_out_).merge(_mom.vertex_out_);
            (this->vertex_in_).merge(_mom.vertex_in_);
        }
        template<class Tep>
        void accumulateAtLeaf2(Tep & ep){}
        void accumulate2(const MomentSearchInAndOut & _mom){}
        void dump(std::ostream & fout=std::cout) const {
            fout<<"vertex_out_.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out_.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in_.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in_.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentSearchInOnly{
    public:
        //F32ort vertex_in_;
        F64ort vertex_in_;
        void init(){
            vertex_in_.init();
        }
        //F32ort getVertexIn() const { return vertex_in_; }
        F64ort getVertexIn() const { return vertex_in_; }
        template<class Tep>
        void accumulateAtLeaf(const Tep & ep){
            (this->vertex_in_).merge(ep.getPos());
        }
        void set(){}
        void accumulate(const MomentSearchInOnly & _mom){
            (this->vertex_in_).merge(_mom.vertex_in_);
        }
        template<class Tep>
        void accumulateAtLeaf2(const Tep & ep){}
        void accumulate2(const MomentSearchInOnly & _mom){}
        void dump(std::ostream & fout=std::cout) const {
            fout<<"vertex_in_.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in_.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    ///////////
    /// SPJ ///
    class SPJMonopoleInAndOut{
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            this->mass = mass;
            this->pos = pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleInAndOut convertToMoment() const {
            return MomentMonopoleInAndOut(mass, pos);
        }
    };
    typedef SPJMonopoleInAndOut SPJMonopoleSymmetry;

    class SPJQuadrupoleInAndOut{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            FSPmat quad = mom.quad;
            this->mass = mass;
            this->pos = pos;
            this->quad = quad;
        }
        MomentQuadrupoleInAndOut convertToMoment() const {
            return MomentQuadrupoleInAndOut(mass, pos, quad);
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
	    quad = 0.0;
        }
        void dump(std::ostream & fout=std::cout) const {
	    fout<<"mass="<<mass<<std::endl;
	    fout<<"pos="<<pos<<std::endl;
	}
    };
    typedef SPJQuadrupoleInAndOut SPJQuadrupoleSymmetry;
    
    
    class SPJMonopole{
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            this->mass = mass;
            this->pos = pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopole convertToMoment() const {
            return MomentMonopole(mass, pos);
        }
    };

    // the same as SPJMonopole
    // for P^3T
    class SPJMonopoleScatter{
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            this->mass = mass;
            this->pos = pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleScatter convertToMoment() const {
            return MomentMonopoleScatter(mass, pos);
        }
        void dump(std::ostream & fout=std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class SPJQuadrupoleScatter{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            FSPmat quad = mom.quad;
            this->mass = mass;
            this->pos = pos;
            this->quad = quad;
        }
        MomentQuadrupoleScatter convertToMoment() const {
            return MomentQuadrupoleScatter(mass, pos, quad);
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
	    quad = 0.0;
        }
        void dump(std::ostream & fout=std::cout) const {
	    fout<<"mass="<<mass<<std::endl;
	    fout<<"pos="<<pos<<std::endl;
	}
    };

    // for P^3T + PM
    class SPJMonopoleCutoffScatter{
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            mass = mom.mass;
            pos = mom.pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleCutoffScatter convertToMoment() const {
            return MomentMonopoleCutoffScatter(mass, pos);
        }
    };

    class SPJQuadrupole{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        void copyFromMoment(const MomentQuadrupole & mom){
            mass = mom.mass;
            pos = mom.pos;
            quad = mom.quad;
        }
        MomentQuadrupole convertToMoment() const {
            return MomentQuadrupole(mass, pos, quad);
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
        }
    };

    class SPJMonopoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        void copyFromMoment(const MomentMonopoleGeometricCenter & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleGeometricCenter convertToMoment() const {
            return MomentMonopoleGeometricCenter(charge, pos, n_ptcl);
        }
    };

    class SPJDipoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        void copyFromMoment(const MomentDipoleGeometricCenter & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentDipoleGeometricCenter convertToMoment() const {
            return MomentDipoleGeometricCenter(charge, pos, n_ptcl, dipole);
        }
    };

    class SPJQuadrupoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        void copyFromMoment(const MomentQuadrupoleGeometricCenter & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
            quadrupole = mom.quadrupole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentQuadrupoleGeometricCenter convertToMoment() const {
            return MomentQuadrupoleGeometricCenter(charge, pos, n_ptcl, dipole, quadrupole);
        }
    };

    class SPJMonopoleCutoff{
    public:
        FSP mass;
        FSPvec pos;
        //static FSP r_cutoff;
        void copyFromMoment(const MomentMonopoleCutoff & mom){
            mass = mom.mass;
            pos = mom.pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }        
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
	/*
        MomentMonopoleCutoff convertToMoment(const FSPvec & r_cutoff) const {
            //const F32ort v_out_dummy(0.0, 0.0);
	    // don't forget to expand box in calc moment global tree
	    const F32ort v_out_dummy(pos-r_cutoff, pos+r_cutoff); 
            return MomentMonopoleCutoff(mass, pos, v_out_dummy);
        }
	*/
        MomentMonopoleCutoff convertToMoment() const {
            //const F32ort v_out_dummy(0.0, 0.0);
            const F64ort v_out_dummy(0.0, 0.0);
            return MomentMonopoleCutoff(mass, pos, v_out_dummy);
        }
    };
    //F32 SPJMonopoleCutoff::r_cutoff;
    typedef SPJMonopoleCutoff SPJMonopolePeriodic;
    
    class SuperParticleBase{
    public:
        void clear(){}
    };

    class EssentialParticleBase{
    public:
        F64vec pos;
        F64    r_search;
        F64vec getPos() const {
            return pos;
        }
        F64 getRSearch() const {
            return r_search;
        }
        void setPos(const F64vec & _pos){
            pos = _pos;
        }
        void clear(){}
    };

}
