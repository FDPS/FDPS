#pragma once

namespace ParticleSimulator {

class TreeParticle {
   private:
   public:
    KeyT key_;
    U32 adr_ptcl_;  // U32 because its MSB is used for distinguish if it is EP or SP

    template <typename Tfp, typename Tmorton_key>
    void setFromFP(const Tfp& fp, const U32 adr, const Tmorton_key& morton_key) {
        key_ = morton_key.getKey(fp.getPos());
        adr_ptcl_ = adr;
    }
    template <class Tep, typename Tmorton_key>
    void setFromEP(const Tep& ep, const U32 adr, const Tmorton_key& morton_key) {
        key_ = morton_key.getKey(ep.getPos());
        adr_ptcl_ = adr;
    }
    template <class Tsp, typename Tmorton_key>
    void setFromSP(const Tsp& sp, const U32 adr, const Tmorton_key& morton_key) {
        key_ = morton_key.getKey(sp.getPos());
        adr_ptcl_ = SetMSB(adr);
    }

    KeyT getKey() const { return key_; }
    void dump(std::ostream& fout = std::cout) const {
        // fout<<std::oct<<"key_="<<key_<<std::dec<<std::endl;
        key_.dump(fout);
        fout << "adr_ptcl_=" << adr_ptcl_ << std::endl;
    }
    U32 getData() { return adr_ptcl_; }
    void setData(U32 val) { adr_ptcl_ = val; }
    void setKey(KeyT k) { key_ = k; }
};

template <typename Tmom, typename Tgeo>
class TreeCell {
   public:
    S32 n_ptcl_;
    U32 adr_tc_;
    U32 adr_ptcl_;
    S32 level_;
#ifdef LOOP_TREE
    U32 adr_tc_next_;
    bool loopFinish() const {
        bool ret = false;
        if (adr_tc_next_ == ADR_TREE_CELL_NULL) {
            ret = true;
        }
        return ret;
    }
#endif
    Tmom mom_;
    Tgeo geo_;
    TreeCell() {
        n_ptcl_ = adr_tc_ = adr_ptcl_ = level_ = 0;
        mom_.init();
    }
    void clear() {
        n_ptcl_ = adr_tc_ = adr_ptcl_ = level_ = 0;
        mom_.init();
        geo_.clear();
    }
    void clearMoment() {
        mom_.init();
        geo_.clear();
    }

    // PMMM
    bool isGroup(const S32 n_group_limit, const S32 lev_group_limit) const {
        return ((n_ptcl_ <= n_group_limit && level_ >= lev_group_limit) || (level_ == TREE_LEVEL_LIMIT));
    }
    bool isGroup(const S32 n_group_limit, const S32 lev_group_limit, const F64 size_limit) const {
        return ((n_ptcl_ <= n_group_limit && level_ >= lev_group_limit && geo_.getEffectiveSize() <= size_limit) || (level_ == TREE_LEVEL_LIMIT));
    }

    bool isLeafPMM(const S32 n_leaf_limit, const S32 lev_leaf_limit) const {
        return ((n_ptcl_ == 0) || (n_ptcl_ <= n_leaf_limit && level_ >= lev_leaf_limit) || (level_ == TREE_LEVEL_LIMIT));
    }

    bool isLeaf(const S32 n_leaf_limit) const { return (n_ptcl_ <= n_leaf_limit || level_ == TREE_LEVEL_LIMIT); }

    // for DEBUG
    void dump(std::ostream& fout = std::cout) const {
        fout << "n_ptcl_=" << n_ptcl_ << std::endl;
        fout << "adr_tc_=" << adr_tc_ << std::endl;
        fout << "adr_ptcl_=" << adr_ptcl_ << std::endl;
        mom_.dump(fout);
        geo_.dump(fout);
    }
    template <class Tep>
    void dumpTree(const Tep* const first_ep_ptr, const TreeCell* const first_tc_ptr, const F32vec& center, const F32 half_length,
                  const S32 n_leaf_limit, std::ostream& fout = std::cout) const {
        fout << std::endl;
        fout << "cell info" << std::endl;
        dump(fout);
        if (!(this->isLeaf(n_leaf_limit))) {
            fout << "this cell is not a leaf" << std::endl;
            fout << "half_length=" << half_length << std::endl;
            fout << "center=" << center << std::endl;
            fout << "level=" << this->level_ << std::endl;
            const TreeCell* child = first_tc_ptr + adr_tc_;
            for (S32 ic = 0; ic < N_CHILDREN; ic++) {
                if ((child + ic)->n_ptcl_ <= 0) continue;
                const Tep* ptcl = first_ep_ptr + adr_ptcl_;
                for (S32 ip = 0; ip < n_ptcl_; ip++, ptcl++) {
                    if (!IsInBox(ptcl->getPos(), center, half_length)) {
                        fout << "out of box(Cell)" << std::endl;
                        fout << "ptcl->getPos()=" << ptcl->getPos() << std::endl;
                        fout << "center=" << center << std::endl;
                        fout << "half_length=" << half_length << std::endl;
                    } else {
                        fout << "in box(Cell)" << std::endl;
                        fout << "ptcl->getPos()=" << ptcl->getPos() << std::endl;
                        fout << "center=" << center << std::endl;
                        fout << "half_length=" << half_length << std::endl;
                    }
                }
                fout << "octid=" << ic << std::endl;
                (child + ic)->dumpTree(first_ep_ptr, first_tc_ptr, center + SHIFT_CENTER[ic] * half_length, half_length * 0.5, n_leaf_limit, fout);
            }
        } else {
            fout << "this cell is a leaf" << std::endl;
            fout << "half_length=" << half_length << std::endl;
            fout << "center=" << center << std::endl;
            fout << "level=" << this->level_ << std::endl;
            const Tep* ptcl = first_ep_ptr + adr_ptcl_;
            for (S32 ip = 0; ip < n_ptcl_; ip++, ptcl++) {
                if (!IsInBox(ptcl->getPos(), center, half_length)) {
                    fout << "out of box(LeafCell)" << std::endl;
                    fout << "ptcl->getPos()=" << ptcl->getPos() << std::endl;
                    fout << "center=" << center << std::endl;
                    fout << "half_length=" << half_length << std::endl;
                } else {
                    fout << "in box(LeafCell)" << std::endl;
                    fout << "ptcl->getPos()=" << ptcl->getPos() << std::endl;
                    fout << "center=" << center << std::endl;
                    fout << "half_length=" << half_length << std::endl;
                }
                // ptcl->dump(fout);
                fout << std::endl;
            }
            fout << std::endl;
        }
    }

    template <class Tep>
    void checkTree(const Tep* const first_ep_ptr, const TreeCell* const first_tc_ptr, const F64vec& center, const F64 half_length,
                   const S64 n_leaf_limit, const S64 lev_leaf_limit, const F64 tolerance, S32& err, std::ostream& fout = std::cout) const {
        if (!(this->isLeaf(n_leaf_limit))) {
            const TreeCell* child = first_tc_ptr + adr_tc_;
            for (S32 ic = 0; ic < N_CHILDREN; ic++) {
                if ((child + ic)->n_ptcl_ <= 0) continue;
                const Tep* ptcl = first_ep_ptr + adr_ptcl_;
                for (S32 ip = 0; ip < n_ptcl_; ip++, ptcl++) {
                    if (!IsInBox(ptcl->getPos(), center, half_length, tolerance)) {
                        fout << "out of box(Cell)" << std::endl;
                        fout << "ptcl->getPos()=" << ptcl->getPos() << std::endl;
                        fout << "center=" << center << std::endl;
                        fout << "half_length=" << half_length << std::endl;
                        fout << "(center+half_length)-ptcl->getPos()=" << (center + half_length) - ptcl->getPos() << std::endl;
                        fout << "(center-half_length)-ptcl->getPos()=" << (center - half_length) - ptcl->getPos() << std::endl;
                        err++;
                    }
                }
                (child + ic)
                    ->checkTree(first_ep_ptr, first_tc_ptr, center + SHIFT_CENTER[ic] * half_length, half_length * 0.5, n_leaf_limit, tolerance, err,
                                fout);
            }
        } else {
            const Tep* ptcl = first_ep_ptr + adr_ptcl_;
            for (S32 ip = 0; ip < n_ptcl_; ip++, ptcl++) {
                if (!IsInBox(ptcl->getPos(), center, half_length, tolerance)) {
                    fout << "out of box(Leaf)" << std::endl;
                    fout << "center=" << center << std::endl;
                    fout << "half_length=" << half_length << std::endl;
                    fout << "ptcl->getPos()=" << ptcl->getPos() << std::endl;
                    fout << "(center+half_length)-ptcl->getPos()=" << (center + half_length) - ptcl->getPos() << std::endl;
                    fout << "(center-half_length)-ptcl->getPos()=" << (center - half_length) - ptcl->getPos() << std::endl;
                    err++;
                }
            }
        }
    }

    template <class Tep, class Tsp>
    void checkTreeLongGlobalTree(const Tep* const first_ep_ptr, const Tsp* const first_sp_ptr, const TreeParticle* const first_tp_ptr,
                                 const TreeCell* const first_tc_ptr, const F32vec& center, const F32 half_length, const S32 n_leaf_limit,
                                 const F32 tolerance, S32& err, std::ostream& fout = std::cout) const {
        if (!(this->isLeaf(n_leaf_limit))) {
            const TreeCell* child = first_tc_ptr + adr_tc_;
            for (S32 ic = 0; ic < N_CHILDREN; ic++) {
                if ((child + ic)->n_ptcl_ <= 0) continue;
                const TreeParticle* tp_tmp = first_tp_ptr + adr_ptcl_;
                for (S32 ip = 0; ip < n_ptcl_; ip++, tp_tmp++) {
                    F32vec pos_tmp;
                    /*
                      const U32 adr = tp_tmp->adr_ptcl_;
                      if( GetMSB(adr) ) pos_tmp = first_sp_ptr[ClearMSB(adr)].getPos();
                      else pos_tmp = first_ep_ptr[adr].getPos();
                    */
                    if (GetMSB(tp_tmp->adr_ptcl_))
                        pos_tmp = first_sp_ptr[adr_ptcl_ + ip].getPos();
                    else
                        pos_tmp = first_ep_ptr[adr_ptcl_ + ip].getPos();
                    // if(Comm::getRank() == 0){
                    if (!IsInBox(pos_tmp, center, half_length, tolerance)) {
                        fout << "out of box(Cell)" << std::endl;
                        fout << "pos_tmp=" << pos_tmp << std::endl;
                        fout << "center=" << center << std::endl;
                        fout << "half_length=" << half_length << std::endl;
                        // fout<<"adr="<<adr<<std::endl;
                        fout << "adr_ptcl_+ip=" << adr_ptcl_ + ip << std::endl;
                        fout << "(center+half_length)-pos_tmp=" << (center + half_length) - pos_tmp << std::endl;
                        fout << "(center-half_length)-pos_tmp=" << (center - half_length) - pos_tmp << std::endl;
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
                (child + ic)
                    ->checkTreeLongGlobalTree(first_ep_ptr, first_sp_ptr, first_tp_ptr, first_tc_ptr, center + SHIFT_CENTER[ic] * half_length,
                                              half_length * 0.5, n_leaf_limit, tolerance, err, fout);
            }
        } else {
            const TreeParticle* tp_tmp = first_tp_ptr + adr_ptcl_;
            for (S32 ip = 0; ip < n_ptcl_; ip++, tp_tmp++) {
                F32vec pos_tmp;
                if (GetMSB(tp_tmp->adr_ptcl_))
                    pos_tmp = first_sp_ptr[adr_ptcl_ + ip].getPos();
                else
                    pos_tmp = first_ep_ptr[adr_ptcl_ + ip].getPos();
                // if(Comm::getRank() == 0){
                if (!IsInBox(pos_tmp, center, half_length, tolerance)) {
                    fout << "out of box(Leaf)" << std::endl;
                    fout << "center=" << center << std::endl;
                    fout << "half_length=" << half_length << std::endl;
                    // fout<<"adr="<<adr<<std::endl;
                    fout << "adr_ptcl_+ip=" << adr_ptcl_ + ip << std::endl;
                    fout << "pos_tmp=" << pos_tmp << std::endl;
                    fout << "(center+half_length)-pos_tmp=" << (center + half_length) - pos_tmp << std::endl;
                    fout << "(center-half_length)-pos_tmp=" << (center - half_length) - pos_tmp << std::endl;
                    err++;
                }
                // }
            }
        }
    }
};

#if 1
// PMMM (just copy from original PMM version but added dump functions)
template <class Tipg>
class IPGroup {};

template <>
class IPGroup<TagIpgLongNormal> {
    // in the case of PMMM
   public:
    S32 n_ptcl_;
    S32 n_epi_;
    S32 adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
    S32 adr_tc_;
#endif
    F64ort vertex_in_;
    template <class Ttc>
    void copyFromTC(const Ttc& tc, const S32 adr_tc) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
        adr_tc_ = adr_tc;
#endif
    }
    void dump(std::ostream& fout = std::cout) const {
        fout << "n_ptcl_= " << n_ptcl_ << " n_epi_= " << n_epi_ << " adr_ptcl_= " << adr_ptcl_ << " vertex_in_= " << vertex_in_ << std::endl;
    }
};

template <>
class IPGroup<TagIpgIn> {
   public:
    S32 n_ptcl_;
    S32 n_epi_;
    S32 adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
    S32 adr_tc_;
#endif
    F64ort vertex_in_;
    template <class Ttc>
    void copyFromTC(const Ttc& tc, const S32 adr_tc) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
        adr_tc_ = adr_tc;
#endif
        // vertex_in_ = tc.mom_.vertex_in_;
        vertex_in_ = tc.geo_.vertex_in_;
    }
    template <class Ttc, class Tep>
    void copyFromTCAndRecalcVertex(const Ttc& tc, const S32 adr_tc, const ReallocatableArray<Tep>& ep_first, const U32 pop_id) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
        adr_tc_ = adr_tc;
#endif
        n_epi_ = 0;
        vertex_in_.init();
        for (S32 i = 0; i < n_ptcl_; i++) {
            Tep* ep = ep_first.getPointer(adr_ptcl_ + i);
            const S32 my_pop_id = GetMyPopId(*ep);
            if (my_pop_id == static_cast<S32>(pop_id)) {
                vertex_in_.merge(ep->getPos());
                n_epi_++;
            }
        }
    }
    void dump(std::ostream& fout = std::cout) const {
        fout << "n_ptcl_= " << n_ptcl_ << " n_epi_= " << n_epi_ << " adr_ptcl_= " << adr_ptcl_ << " vertex_in_= " << vertex_in_ << std::endl;
    }
};

template <>
class IPGroup<TagIpgInAndOut> {
   public:
    S32 n_ptcl_;
    S32 n_epi_;
    S32 adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
    S32 adr_tc_;
#endif
    F64ort vertex_in_;
    F64ort vertex_out_;
    template <class Ttc>
    void copyFromTC(const Ttc& tc, const S32 adr_tc) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
        adr_tc_ = adr_tc;
#endif
        // vertex_in_ = tc.mom_.vertex_in_;
        // vertex_out_ = tc.mom_.vertex_out_;
        vertex_in_ = tc.geo_.vertex_in_;
        vertex_out_ = tc.geo_.vertex_out_;
    }
    template <class Ttc, class Tep>
    void copyFromTCAndRecalcVertex(const Ttc& tc, const S32 adr_tc, const ReallocatableArray<Tep>& ep_first, const U32 pop_id) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
        adr_tc_ = adr_tc;
#endif
        n_epi_ = 0;
        vertex_in_.init();
        vertex_out_.init();
        for (S32 i = 0; i < n_ptcl_; i++) {
            Tep* ep = ep_first.getPointer(adr_ptcl_ + i);
            const S32 my_pop_id = GetMyPopId(*ep);
            if (my_pop_id == static_cast<S32>(pop_id)) {
                const F64vec pos = ep->getPos();
                const F64 r_srch = GetMyRSearch(*ep);
                vertex_in_.merge(pos);
                vertex_out_.merge(pos, r_srch);
                n_epi_++;
            }
        }
    }
    void dump(std::ostream& fout = std::cout) const {
        fout << "n_ptcl_= " << n_ptcl_ << " n_epi_= " << n_epi_ << " adr_ptcl_= " << adr_ptcl_ << " vertex_in_= " << vertex_in_
             << " vertex_out_= " << vertex_out_ << std::endl;
    }
};

template <>
class IPGroup<TagIpgOut> {
   public:
    S32 n_ptcl_;
    S32 n_epi_;
    S32 adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
    S32 adr_tc_;
#endif
    F64ort vertex_out_;
    template <class Ttc>
    void copyFromTC(const Ttc& tc, const S32 adr_tc) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
        adr_tc_ = adr_tc;
#endif
        // vertex_out_ = tc.mom_.vertex_out_;
        vertex_out_ = tc.geo_.vertex_out_;
    }
    template <class Ttc, class Tep>
    void copyFromTCAndRecalcVertex(const Ttc& tc, const S32 adr_tc, const ReallocatableArray<Tep>& ep_first, const U32 pop_id) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
#ifdef PARTICLE_SIMULATOR_MAKE_IPGROUP_HAVE_ADR_TC
        adr_tc_ = adr_tc;
#endif
        n_epi_ = 0;
        vertex_out_.init();
        for (S32 i = 0; i < n_ptcl_; i++) {
            Tep* ep = ep_first.getPointer(adr_ptcl_ + i);
            const S32 my_pop_id = GetMyPopId(*ep);
            if (my_pop_id == static_cast<S32>(pop_id)) {
                const F64vec pos = ep->getPos();
                const F64 r_srch = GetMyRSearch(*ep);
                vertex_out_.merge(pos, r_srch);
                n_epi_++;
            }
        }
    }
    void dump(std::ostream& fout = std::cout) const {
        fout << "n_ptcl_= " << n_ptcl_ << " n_epi_= " << n_epi_ << " adr_ptcl_= " << adr_ptcl_ << " vertex_out_= " << vertex_out_ << std::endl;
    }
};

#else
// befor PMMM
template <class Tipg>
class IPGroup {};

template <>
class IPGroup<TagIpgLongNormal> {
   public:
    S32 n_ptcl_;
    S32 adr_ptcl_;
    F64ort vertex_in_;
    template <class Ttc>
    void copyFromTC(const Ttc& tc) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
    }
};

template <>
class IPGroup<TagIpgIn> {
   public:
    S32 n_ptcl_;
    S32 adr_ptcl_;
    F64ort vertex_in_;
    template <class Ttc>
    void copyFromTC(const Ttc& tc) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
        // vertex_in_ = tc.mom_.vertex_in_;
        vertex_in_ = tc.geo_.vertex_in_;
    }
};

template <>
class IPGroup<TagIpgInAndOut> {
   public:
    S32 n_ptcl_;
    S32 adr_ptcl_;
    F64ort vertex_in_;
    F64ort vertex_out_;
    template <class Ttc>
    void copyFromTC(const Ttc& tc) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
        // vertex_in_ = tc.mom_.vertex_in_;
        // vertex_out_ = tc.mom_.vertex_out_;
        vertex_in_ = tc.geo_.vertex_in_;
        vertex_out_ = tc.geo_.vertex_out_;
    }
};

template <>
class IPGroup<TagIpgOut> {
   public:
    S32 n_ptcl_;
    S32 adr_ptcl_;
    F64ort vertex_out_;
    template <class Ttc>
    void copyFromTC(const Ttc& tc) {
        n_ptcl_ = tc.n_ptcl_;
        adr_ptcl_ = tc.adr_ptcl_;
        // vertex_out_ = tc.mom_.vertex_out_;
        vertex_out_ = tc.geo_.vertex_out_;
    }
};
#endif

////////////////////////
// ESSENTIAL PARTICLE //
class EPXROnly {
    F64 r_search;
    F64vec pos;

   public:
    F64 getRSearch() const { return r_search; }
    F64vec getPos() const { return pos; }
    void setPos(const F64vec& pos_in) { pos = pos_in; }
    template <class Tep>
    void copyFromEP(const Tep& ep) {
        r_search = ep.getRSearch();
        pos = ep.getPos();
    }
    template <class Tptcl>
    const EPXROnly& operator=(const Tptcl& ptcl) {
        r_search = ptcl.getRSearch();
        pos = ptcl.getPos();
        return (*this);
    }
};

//////////////
/// Moment ///
// new
// these are the same as
// for P^3T
// for P^3T + PM
class MomentMonopoleCutoffScatter {
   public:
    FSP mass;
    FSPvec pos;
    // F64ort vertex_out_; // cutoff
    // F64ort vertex_out2_; // search ep
    // F64ort vertex_in_;
    MomentMonopoleCutoffScatter() {
        mass = 0.0;
        pos = 0.0;
        // vertex_out_.init();
        // vertex_out2_.init();
        // vertex_in_.init();
    }
    MomentMonopoleCutoffScatter(const FSP m, const FSPvec& p) {
        mass = m;
        pos = p;
    }
    // F64ort getVertexOut() const { return vertex_out_; }
    // F64ort getVertexOut2() const { return vertex_out2_; }
    // F64ort getVertexIn() const { return vertex_in_; }
    void init() {
        mass = 0.0;
        pos = 0.0;
        // vertex_out_.init();
        // vertex_out2_.init();
        // vertex_in_.init();
    }
    FSPvec getPos() const { return pos; }
    FSP getCharge() const { return mass; }
    template <class Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        this->mass += epj.getCharge();
        this->pos += epj.getCharge() * epj.getPos();
        //(this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        //(this->vertex_out2_).merge(epj.getPos(), epj.getRSearch2());
        //(this->vertex_in_).merge(epj.getPos());
    }
    template <class Tepj>
    void accumulateAtLeaf2(const Tepj& epj) {}
    void set() { pos = pos / mass; }
    void accumulate(const MomentMonopoleCutoffScatter& mom) {
        this->mass += mom.mass;
        this->pos += mom.mass * mom.pos;
        //(this->vertex_out_).merge(mom.vertex_out_);
        //(this->vertex_out2_).merge(mom.vertex_out2_);
        //(this->vertex_in_).merge(mom.vertex_in_);
    }
    void accumulate2(const MomentMonopoleCutoffScatter& mom) {}
    FSP getMonopole() const { return mass; }
    // for DEBUG
    void dump(std::ostream& fout = std::cout) const {
        fout << "mass=" << mass << std::endl;
        fout << "pos=" << pos << std::endl;
        // fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
        // fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
        // fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
        // fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
    }
};

class MomentMonopole {
   public:
    FSP mass;
    FSPvec pos;
    MomentMonopole() {
        mass = 0.0;
        pos = 0.0;
    }
    MomentMonopole(const FSP m, const FSPvec& p) {
        mass = m;
        pos = p;
    }
    void init() {
        mass = 0.0;
        pos = 0.0;
    }
    FSPvec getPos() const { return pos; }
    FSP getCharge() const { return mass; }
    template <class Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        mass += epj.getCharge();
        pos += epj.getCharge() * epj.getPos();
    }
    template <class Tepj>
    void accumulateAtLeaf2(const Tepj& epj) {}
    void set() { pos = pos / mass; }
    void accumulate(const MomentMonopole& mom) {
        mass += mom.mass;
        pos += mom.mass * mom.pos;
    }
    void accumulate2(const MomentMonopole& mom) {}
    FSP getMonopole() const { return mass; }
    // for DEBUG
    void dump(std::ostream& fout = std::cout) const {
        fout << "mass=" << mass << std::endl;
        fout << "pos=" << pos << std::endl;
    }
};


class MomentMonopoleCutoff {
   public:
    FSP mass;
    FSPvec pos;
    F64ort vertex_out_;
    MomentMonopoleCutoff() {
        mass = 0.0;
        pos = 0.0;
        vertex_out_.init();
    }
    MomentMonopoleCutoff(const FSP m,
                         const FSPvec& p,
                         const F64ort& v_out) {
        mass = m;
        pos = p;
        vertex_out_ = v_out;
    }
    void init() {
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
    template <class Tepj>
    void accumulateAtLeaf(Tepj& epj) {
        mass += epj.getCharge();
        pos += epj.getCharge() * epj.getPos();
        (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
    }
    template <class Tepj>
    void accumulateAtLeaf2(Tepj& epj) {}
    void set() {
        pos = pos / mass;
    }
    void accumulate(const MomentMonopoleCutoff& mom) {
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        (this->vertex_out_).merge(mom.vertex_out_);
    }
    void accumulate2(const MomentMonopoleCutoff& mom) {}
    // for DEBUG
    void dump(std::ostream& fout = std::cout) const {
        fout << "mass=" << mass << std::endl;
        fout << "pos=" << pos << std::endl;
        fout << "vertex_out_=" << vertex_out_ << std::endl;
    }
};

class MomentQuadrupole {
   public:
    FSPvec pos;
    FSP mass;
    FSPmat quad;
    void init() {
        pos = 0.0;
        mass = 0.0;
        quad = 0.0;
    }
    MomentQuadrupole() {
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
    }
    MomentQuadrupole(const FSP m, const FSPvec& p, const FSPmat& q) {
        mass = m;
        pos = p;
        quad = q;
    }
    FSPvec getPos() const { return pos; }
    template <class Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        mass += epj.getCharge();
        pos += epj.getCharge() * epj.getPos();
    }
    template <class Tepj>
    void accumulateAtLeaf2(const Tepj& epj) {
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
    void set() { pos = pos / mass; }
    void accumulate(const MomentQuadrupole& mom) {
        mass += mom.mass;
        pos += mom.mass * mom.pos;
    }
    void accumulate2(const MomentQuadrupole& mom) {
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
    FSP getMonopole() const { return mass; }
    FSPmat getQuadrupole() const { return quad; }
    void dump(std::ostream& fout = std::cout) const {
        fout << "mass= " << mass << std::endl;
        fout << "pos= " << pos << std::endl;
        fout << "quad= " << quad << std::endl;
    }
};

class MomentMonopoleGeometricCenter {
   public:
    SSP n_ptcl;
    FSP charge;
    FSPvec pos;
    MomentMonopoleGeometricCenter() {
        n_ptcl = 0;
        charge = 0.0;
        pos = 0.0;
    }
    MomentMonopoleGeometricCenter(const FSP c, const FSPvec& p, const SSP n) {
        n_ptcl = n;
        charge = c;
        pos = p;
    }
    void init() {
        n_ptcl = 0;
        charge = 0.0;
        pos = 0.0;
    }
    FSPvec getPos() const { return pos; }
    FSP getCharge() const { return charge; }
    template <class Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        charge += epj.getCharge();
        pos += epj.getPos();
        n_ptcl++;
    }
    template <class Tepj>
    void accumulateAtLeaf2(const Tepj& epj) {}
    void set() { pos = pos / n_ptcl; }
    void accumulate(const MomentMonopoleGeometricCenter& mom) {
        charge += mom.charge;
        pos += mom.n_ptcl * mom.pos;
        n_ptcl += mom.n_ptcl;
    }
    void accumulate2(const MomentMonopoleGeometricCenter& mom) {}
    // for DEBUG
    void dump(std::ostream& fout = std::cout) const {
        fout << "charge=" << charge << std::endl;
        fout << "pos=" << pos << std::endl;
        fout << "n_ptcl=" << n_ptcl << std::endl;
    }
    FSP getMonopole() const { return charge; }
};

class MomentDipoleGeometricCenter {
   public:
    SSP n_ptcl;
    FSP charge;
    FSPvec pos;
    FSPvec dipole;
    MomentDipoleGeometricCenter() {
        n_ptcl = 0;
        charge = 0.0;
        pos = 0.0;
        dipole = 0.0;
    }
    MomentDipoleGeometricCenter(const FSP c, const FSPvec& p, const SSP n, const FSPvec& di) {
        charge = c;
        pos = p;
        n_ptcl = n;
        dipole = di;
    }
    void init() {
        n_ptcl = 0;
        charge = 0.0;
        pos = dipole = 0.0;
    }
    FSPvec getPos() const { return pos; }
    FSP getCharge() const { return charge; }
    template <class Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        charge += epj.getCharge();
        pos += epj.getPos();
        n_ptcl++;
    }
    void accumulate(const MomentDipoleGeometricCenter& mom) {
        charge += mom.charge;
        pos += mom.n_ptcl * mom.pos;
        n_ptcl += mom.n_ptcl;
    }
    void set() { pos = pos / n_ptcl; }
    template <class Tepj>
    void accumulateAtLeaf2(const Tepj& epj) {
        this->dipole += epj.getCharge() * (epj.getPos() - this->pos);
    }
    void accumulate2(const MomentDipoleGeometricCenter& mom) {
        // dipole += mom.charge * (mom.pos - this->pos);
        this->dipole += mom.charge * (mom.pos - this->pos) + mom.dipole;
    }
    FSP getMonopole() const { return charge; }
    FSPvec getDipole() const { return dipole; }
    // for DEBUG
    void dump(std::ostream& fout = std::cout) const {
        fout << "charge=" << charge << std::endl;
        fout << "pos=" << pos << std::endl;
        fout << "n_ptcl=" << n_ptcl << std::endl;
        fout << "dipole=" << dipole << std::endl;
    }
};

class MomentQuadrupoleGeometricCenter {
   public:
    SSP n_ptcl;
    FSP charge;
    FSPvec pos;
    FSPvec dipole;
    FSPmat quadrupole;
    MomentQuadrupoleGeometricCenter() {
        n_ptcl = 0;
        charge = 0.0;
        pos = 0.0;
        dipole = 0.0;
        quadrupole = 0.0;
    }
    MomentQuadrupoleGeometricCenter(const FSP c, const FSPvec& p, const SSP n, const FSPvec& di, const FSPmat& q) {
        charge = c;
        pos = p;
        n_ptcl = n;
        dipole = di;
        quadrupole = q;
    }
    void init() {
        n_ptcl = 0;
        charge = 0.0;
        pos = dipole = 0.0;
        quadrupole = 0.0;
    }
    FSPvec getPos() const { return pos; }
    FSP getCharge() const { return charge; }
    template <class Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        charge += epj.getCharge();
        pos += epj.getPos();
        n_ptcl++;
    }
    void set() { pos = pos / n_ptcl; }
    void accumulate(const MomentQuadrupoleGeometricCenter& mom) {
        charge += mom.charge;
        pos += mom.n_ptcl * mom.pos;
        n_ptcl += mom.n_ptcl;
    }
    template <class Tepj>
    void accumulateAtLeaf2(const Tepj& epj) {
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
    void accumulate2(const MomentQuadrupoleGeometricCenter& mom) {
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
#if 1
        auto ctmp = mom.charge;
        auto ptmp = mom.pos - this->pos;
        auto di_src = mom.dipole;
        this->dipole += di_src + ctmp * ptmp;
        F64 cx = ctmp * ptmp.x;
        F64 cy = ctmp * ptmp.y;
        F64 cz = ctmp * ptmp.z;
        this->quadrupole.xx += mom.quadrupole.xx + 2.0 * ptmp.x * di_src.x + cx * ptmp.x;
        this->quadrupole.yy += mom.quadrupole.yy + 2.0 * ptmp.y * di_src.y + cy * ptmp.y;
        this->quadrupole.zz += mom.quadrupole.zz + 2.0 * ptmp.z * di_src.z + cz * ptmp.z;
        this->quadrupole.xy += mom.quadrupole.xy + ptmp.x * di_src.y + ptmp.y * di_src.x + cx * ptmp.y;
        this->quadrupole.xz += mom.quadrupole.xz + ptmp.x * di_src.z + ptmp.z * di_src.x + cx * ptmp.z;
        this->quadrupole.yz += mom.quadrupole.yz + ptmp.y * di_src.z + ptmp.z * di_src.y + cy * ptmp.z;

        /*
        this->quadrupole.xx += cx * ptmp.x + mom.quadrupole.xx - (mom.dipole.x * ptmp.x + mom.dipole.x * ptmp.x);
        this->quadrupole.yy += cy * ptmp.y + mom.quadrupole.yy - (mom.dipole.y * ptmp.y + mom.dipole.y * ptmp.y);
        this->quadrupole.zz += cz * ptmp.z + mom.quadrupole.zz - (mom.dipole.z * ptmp.z + mom.dipole.z * ptmp.z);
        this->quadrupole.xy += cx * ptmp.y + mom.quadrupole.xy - (mom.dipole.x * ptmp.y + mom.dipole.y * ptmp.x);
        this->quadrupole.xz += cx * ptmp.z + mom.quadrupole.xz - (mom.dipole.x * ptmp.z + mom.dipole.z * ptmp.x);
        this->quadrupole.yz += cy * ptmp.z + mom.quadrupole.yz - (mom.dipole.y * ptmp.z + mom.dipole.z * ptmp.y);
        */
#else
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
#endif
#else
        // underconstruction
#endif
    }
    FSP getMonopole() const { return charge; }
    FSPvec getDipole() const { return dipole; }    
    FSPmat getQuadrupole() const { return quadrupole; }
    // for DEBUG
    void dump(std::ostream& fout = std::cout) const {
        fout << "charge=" << charge << std::endl;
        fout << "pos=" << pos << std::endl;
        fout << "n_ptcl=" << n_ptcl << std::endl;
        fout << "dipole=" << dipole << std::endl;
        fout << "quadrupole=" << quadrupole << std::endl;
    }
};

class MomentShort {
   public:
    void init() {}
    template <class Tepj>
    void accumulateAtLeaf(Tepj& epj) {}
    template <class Tepj>
    void accumulateAtLeaf2(Tepj& epj) {}
    void set() {}
    void accumulate(const MomentShort& mom) {}
    void accumulate2(const MomentShort& mom) {}
    // for DEBUG
    void dump(std::ostream& fout = std::cout) const {}
};

///////////
/// SPJ ///
class SPJMonopole {
   public:
    FSP mass;
    FSPvec pos;

    template <class Tmom>
    void copyFromMoment(const Tmom& mom) {
        FSP mass = mom.mass;
        FSPvec pos = mom.pos;
        this->mass = mass;
        this->pos = pos;
    }
    void clear() {
        mass = 0.0;
        pos = 0.0;
    }
    FSP getCharge() const { return mass; }
    FSPvec getPos() const { return pos; }
    void setPos(const FSPvec& pos_new) { pos = pos_new; }
    MomentMonopole convertToMoment() const { return MomentMonopole(mass, pos); }
    FSP getMonopole() const { return getCharge(); }
};

// for P^3T + PM
class SPJMonopoleCutoffScatter {
   public:
    FSP mass;
    FSPvec pos;
    template <class Tmom>
    void copyFromMoment(const Tmom& mom) {
        mass = mom.mass;
        pos = mom.pos;
    }
    void clear() {
        mass = 0.0;
        pos = 0.0;
    }
    FSP getCharge() const { return mass; }
    FSPvec getPos() const { return pos; }
    void setPos(const FSPvec& pos_new) { pos = pos_new; }
    MomentMonopoleCutoffScatter convertToMoment() const { return MomentMonopoleCutoffScatter(mass, pos); }
    FSP getMonopole() const { return getCharge(); }
};

/*
class SPJMonopoleCutoff {
   public:
    FSP mass;
    FSPvec pos;
    // static FSP r_cutoff;
    void copyFromMoment(const MomentMonopoleCutoff& mom) {
        mass = mom.mass;
        pos = mom.pos;
    }
    void clear() {
        mass = 0.0;
        pos = 0.0;
    }
    FSP getCharge() const {
        return mass;
    }
    FSPvec getPos() const {
        return pos;
    }
    void setPos(const FSPvec& pos_new) {
        pos = pos_new;
    }
    MomentMonopoleCutoff convertToMoment() const {
        // const F32ort v_out_dummy(0.0, 0.0);
        const F64ort v_out_dummy(0.0, 0.0);
        return MomentMonopoleCutoff(mass, pos, v_out_dummy);
    }
};
*/
  //using SPJMonopoleCutoff = SPJMonopoleCutoffScatter;
  using SPJMonopoleCutoff = SPJMonopole;

class SPJQuadrupole {
   public:
    FSP mass;
    FSPvec pos;
    FSPmat quad;
    FSP getCharge() const { return mass; }
    FSPvec getPos() const { return pos; }
    void setPos(const FSPvec& pos_new) { pos = pos_new; }
    void copyFromMoment(const MomentQuadrupole& mom) {
        mass = mom.mass;
        pos = mom.pos;
        quad = mom.quad;
    }
    MomentQuadrupole convertToMoment() const { return MomentQuadrupole(mass, pos, quad); }
    void clear() {
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
    }
    FSP getMonopole() const { return getCharge(); }
    FSPmat getQuadrupole() const { return quad; }
};

class SPJMonopoleGeometricCenter {
   public:
    SSP n_ptcl;
    FSP charge;
    FSPvec pos;
    void copyFromMoment(const MomentMonopoleGeometricCenter& mom) {
        n_ptcl = mom.n_ptcl;
        charge = mom.charge;
        pos = mom.pos;
    }
    void clear() {
        n_ptcl = 0;
        charge = 0.0;
        pos = 0.0;
    }
    FSP getCharge() const { return charge; }
    FSPvec getPos() const { return pos; }
    void setPos(const FSPvec& pos_new) { pos = pos_new; }
    MomentMonopoleGeometricCenter convertToMoment() const { return MomentMonopoleGeometricCenter(charge, pos, n_ptcl); }
    FSP getMonopole() const { return getCharge(); }
};

class SPJDipoleGeometricCenter {
   public:
    SSP n_ptcl;
    FSP charge;
    FSPvec pos;
    FSPvec dipole;
    void copyFromMoment(const MomentDipoleGeometricCenter& mom) {
        n_ptcl = mom.n_ptcl;
        charge = mom.charge;
        pos = mom.pos;
        dipole = mom.dipole;
    }
    void clear() {
        n_ptcl = 0;
        charge = 0.0;
        pos = 0.0;
        dipole = 0.0;
    }
    FSP getCharge() const { return charge; }
    FSPvec getPos() const { return pos; }
    void setPos(const FSPvec& pos_new) { pos = pos_new; }
    MomentDipoleGeometricCenter convertToMoment() const { return MomentDipoleGeometricCenter(charge, pos, n_ptcl, dipole); }
    FSP getMonopole() const { return getCharge(); }
    FSPvec getDipole() const { return dipole; }        
};

class SPJQuadrupoleGeometricCenter {
   public:
    SSP n_ptcl;
    FSP charge;
    FSPvec pos;
    FSPvec dipole;
    FSPmat quadrupole;
    void copyFromMoment(const MomentQuadrupoleGeometricCenter& mom) {
        n_ptcl = mom.n_ptcl;
        charge = mom.charge;
        pos = mom.pos;
        dipole = mom.dipole;
        quadrupole = mom.quadrupole;
    }
    void clear() {
        n_ptcl = 0;
        charge = 0.0;
        pos = 0.0;
        dipole = 0.0;
        quadrupole = 0.0;
    }
    FSP getCharge() const { return charge; }
    FSPvec getPos() const { return pos; }
    void setPos(const FSPvec& pos_new) { pos = pos_new; }
    MomentQuadrupoleGeometricCenter convertToMoment() const { return MomentQuadrupoleGeometricCenter(charge, pos, n_ptcl, dipole, quadrupole); }
    FSP getMonopole() const { return getCharge(); }
    FSPvec getDipole() const { return dipole; }    
    FSPmat getQuadrupole() const { return quadrupole; } 
    void dump(std::ostream& fout = std::cout){
        fout<<"pos= "<<pos<<std::endl;
        fout<<"charge= "<<charge<<std::endl;
        fout<<"dipole= "<<dipole<<std::endl;
        fout<<"quadrupole= "<<quadrupole<<std::endl;
    }
};

class SuperParticleBase {
   public:
    void clear() {}
};

class EssentialParticleBase {
   public:
    F64vec pos;
    F64 r_search;
    F64vec getPos() const { return pos; }
    F64 getRSearch() const { return r_search; }
    void setPos(const F64vec& _pos) { pos = _pos; }
    void clear() {}
};

/////////////////////////////////
// GEOMETRIC INFO OF TREE CELL //
template <typename T>
class Geometry {};

template <>
class Geometry<TagGeometrySize> {
   public:
    F64 size_;
    void setSize(const F64 s) { size_ = s; }
    F64 getSize() const { return size_; }
    void clear() { size_ = 0.0; }
    F64 getEffectiveSize() const { return size_; }
    template <typename Tepj>
    void addFirstEpj(const Tepj& epj) {}
    template <typename Tspj>
    void addFirstSpj(const Tspj& spj) {}
    void copyBox(const Geometry<TagGeometrySize>& box) {}
    template <typename Tepj>
    void accumulateAtLeaf(const Tepj& epj) {}
    void accumulate(const Geometry<TagGeometrySize>& g) { size_ = (size_ > g.size_) ? size_ : g.size_; }
    void dump(std::ostream& fout = std::cerr) const { fout << "size_= " << size_ << std::endl; }
};

//////
// for PMMM
template <>
class Geometry<TagGeometrySizeIn> {
   public:
    F64 size_;
    F64ort vertex_in_;
    void clear() {
        size_ = 0.0;
        vertex_in_.low_ = 1.0;
        vertex_in_.high_ = -1.0;
    }
    void setSize(const F64 s) { size_ = s; }
    F64 getSize() const { return size_; }
    F64ort getVertexIn() const { return vertex_in_; }
    F64 getEffectiveSize() const { return vertex_in_.getFullLength().getMax(); }
    template <typename Tepj>
    void addFirstEpj(const Tepj& epj) {
        vertex_in_.low_ = epj.getPos();
        vertex_in_.high_ = epj.getPos();
    }
    template <typename Tspj>
    void addFirstSpj(const Tspj& spj) {
        vertex_in_.low_ = vertex_in_.high_ = spj.getPos();
    }
    void copyBox(const Geometry<TagGeometrySizeIn>& box) { vertex_in_ = box.vertex_in_; }
    template <typename Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        vertex_in_.merge(epj.getPos());
    }
    void accumulate(const Geometry<TagGeometrySizeIn>& g) {
        vertex_in_.merge(g.vertex_in_);
        size_ = (size_ > g.size_) ? size_ : g.size_;
    }
    void dump(std::ostream& fout = std::cerr) const {
        fout << "size_= " << size_ << std::endl;
        fout << "vertex_in_= " << vertex_in_ << std::endl;
    }
};
// for PMMM
//////

template <>
class Geometry<TagGeometrySizeOut> {
   public:
    F64 size_;
    F64ort vertex_out_;
    void clear() {
        size_ = 0.0;
        vertex_out_.low_ = 1.0;
        vertex_out_.high_ = -1.0;
    }
    void setSize(const F64 s) { size_ = s; }
    F64 getSize() const { return size_; }
    F64 getEffectiveSize() const { return size_; }
    F64ort getVertexOut() const { return vertex_out_; }
    template <typename Tepj>
    void addFirstEpj(const Tepj& epj) {
        vertex_out_.low_ = epj.getPos() - epj.getRSearch();
        vertex_out_.high_ = epj.getPos() + epj.getRSearch();
    }
    template <typename Tspj>
    void addFirstSpj(const Tspj& spj) {
        vertex_out_.low_ = vertex_out_.high_ = spj.getPos();
    }
    void copyBox(const Geometry<TagGeometrySizeOut>& box) { vertex_out_ = box.vertex_out_; }
    template <typename Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        vertex_out_.merge(epj.getPos(), epj.getRSearch());
    }
    void accumulate(const Geometry<TagGeometrySizeOut>& g) {
        vertex_out_.merge(g.vertex_out_);
        size_ = (size_ > g.size_) ? size_ : g.size_;
    }
    void dump(std::ostream& fout = std::cerr) const {
        fout << "size_= " << size_ << std::endl;
        fout << "vertex_out_= " << vertex_out_ << std::endl;
    }
};

template <>
class Geometry<TagGeometrySizeInOut> {
   public:
    F64 size_;
    F64ort vertex_in_;
    F64ort vertex_out_;
    void clear() {
        size_ = 0.0;
        vertex_in_.low_ = vertex_out_.low_ = 1.0;
        vertex_in_.high_ = vertex_out_.high_ = -1.0;
    }
    void setSize(const F64 s) { size_ = s; }
    F64 getSize() const { return size_; }
    F64 getEffectiveSize() const {
        return vertex_in_.getFullLength().getMax();
        // return size_;
    }
    F64ort getVertexIn() const { return vertex_in_; }
    F64ort getVertexOut() const { return vertex_out_; }
    template <typename Tepj>
    void addFirstEpj(const Tepj& epj) {
        vertex_out_.low_ = epj.getPos() - epj.getRSearch();
        vertex_out_.high_ = epj.getPos() + epj.getRSearch();
        vertex_in_.low_ = epj.getPos();
        vertex_in_.high_ = epj.getPos();
    }
    template <typename Tspj>
    void addFirstSpj(const Tspj& spj) {
        vertex_out_.low_ = vertex_out_.high_ = spj.getPos();
        vertex_in_.low_ = vertex_in_.high_ = spj.getPos();
    }
    void copyBox(const Geometry<TagGeometrySizeInOut>& box) {
        vertex_out_ = box.vertex_out_;
        vertex_in_ = box.vertex_in_;
    }
    template <typename Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        vertex_in_.merge(epj.getPos());
        vertex_out_.merge(epj.getPos(), epj.getRSearch());
    }
    void accumulate(const Geometry<TagGeometrySizeInOut>& g) {
        vertex_out_.merge(g.vertex_out_);
        vertex_in_.merge(g.vertex_in_);
        size_ = (size_ > g.size_) ? size_ : g.size_;
    }
    void dump(std::ostream& fout = std::cerr) const {
        fout << "size_= " << size_ << std::endl;
        fout << "vertex_in_= " << vertex_in_ << std::endl;
        fout << "vertex_out_= " << vertex_out_ << std::endl;
    }
};

template <>
class Geometry<TagGeometryOut> {
   public:
    F64ort vertex_out_;
    void clear() {
        vertex_out_.low_ = 1.0;
        vertex_out_.high_ = -1.0;
    }
    F64ort getVertexOut() const { return vertex_out_; }
    /*
    F64 getEffectiveSize() const {
        return size_;
    }
    */
    template <typename Tepj>
    void addFirstEpj(const Tepj& epj) {
        vertex_out_.low_ = epj.getPos() - epj.getRSearch();
        vertex_out_.high_ = epj.getPos() + epj.getRSearch();
    }
    template <typename Tspj>
    void addFirstSpj(const Tspj& spj) {
        vertex_out_.low_ = vertex_out_.high_ = spj.getPos();
    }
    void copyBox(const Geometry<TagGeometryOut>& box) { vertex_out_ = box.vertex_out_; }
    template <typename Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        vertex_out_.merge(epj.getPos(), epj.getRSearch());
    }
    void accumulate(const Geometry<TagGeometryOut>& g) { vertex_out_.merge(g.vertex_out_); }
    void dump(std::ostream& fout = std::cerr) const { fout << "vertex_out_= " << vertex_out_ << std::endl; }
};

template <>
class Geometry<TagGeometryIn> {
   public:
    F64ort vertex_in_;
    void clear() {
        vertex_in_.low_ = 1.0;
        vertex_in_.high_ = -1.0;
    }
    F64ort getVertexIn() const { return vertex_in_; }
    F64 getEffectiveSize() const { return vertex_in_.getFullLength().getMax(); }
    template <typename Tepj>
    void addFirstEpj(const Tepj& epj) {
        vertex_in_.low_ = epj.getPos();
        vertex_in_.high_ = epj.getPos();
    }
    template <typename Tspj>
    void addFirstSpj(const Tspj& spj) {
        vertex_in_.low_ = vertex_in_.high_ = spj.getPos();
    }
    void copyBox(const Geometry<TagGeometryIn>& box) { vertex_in_ = box.vertex_in_; }
    template <typename Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        vertex_in_.merge(epj.getPos());
    }
    void accumulate(const Geometry<TagGeometryIn>& g) { vertex_in_.merge(g.vertex_in_); }
    void dump(std::ostream& fout = std::cerr) const { fout << "vertex_in_= " << vertex_in_ << std::endl; }
};

template <>
class Geometry<TagGeometryInAndOut> {
   public:
    F64ort vertex_in_;
    F64ort vertex_out_;
    void clear() {
        vertex_in_.low_ = vertex_out_.low_ = 1.0;
        vertex_in_.high_ = vertex_out_.high_ = -1.0;
    }
    F64 getEffectiveSize() const { return vertex_in_.getFullLength().getMax(); }
    F64ort getVertexIn() const { return vertex_in_; }
    F64ort getVertexOut() const { return vertex_out_; }
    template <typename Tepj>
    void addFirstEpj(const Tepj& epj) {
        vertex_out_.low_ = epj.getPos() - epj.getRSearch();
        vertex_out_.high_ = epj.getPos() + epj.getRSearch();
        vertex_in_.low_ = epj.getPos();
        vertex_in_.high_ = epj.getPos();
    }
    template <typename Tspj>
    void addFirstSpj(const Tspj& spj) {
        vertex_in_.low_ = vertex_in_.high_ = spj.getPos();
        vertex_out_.low_ = vertex_out_.high_ = spj.getPos();
    }
    void copyBox(const Geometry<TagGeometryInAndOut>& box) {
        vertex_in_ = box.vertex_in_;
        vertex_out_ = box.vertex_out_;
    }
    template <typename Tepj>
    void accumulateAtLeaf(const Tepj& epj) {
        vertex_in_.merge(epj.getPos());
        vertex_out_.merge(epj.getPos(), epj.getRSearch());
    }
    void accumulate(const Geometry<TagGeometryInAndOut>& g) {
        vertex_in_.merge(g.vertex_in_);
        vertex_out_.merge(g.vertex_out_);
    }
    void dump(std::ostream& fout = std::cerr) const {
        fout << "vertex_out_= " << vertex_out_ << std::endl;
        fout << "vertex_in_= " << vertex_in_ << std::endl;
    }
};

////////////
// targe box class
template <class TSM>
struct TargetBox {
    F64ort vertex_;
};

template <>
struct TargetBox<SEARCH_MODE_LONG> {
    F64ort vertex_;
    template <class Tipg>
    void set(const Tipg& ipg) {
        vertex_ = ipg.vertex_in_;
    }

    // for cinsistency of PMM mode
    template <class Tipg, class Tmorton_key>
    void set(const Tipg& ipg, const S32 icut,
             const Tmorton_key& morton_key) {  // used in MakeInteractionList*
        vertex_ = ipg.vertex_in_;
    }

    template <typename T>
    void setForExLet(const T& top, const F64vec& s) {
        vertex_ = top.vertex_in_.shift(s);
    }
    void dump(std::ostream& fout = std::cerr) { fout << "vertex_= " << vertex_ << std::endl; }
};

template <>
struct TargetBox<SEARCH_MODE_LONG_SYMMETRY> {
    F64ort vertex_in_;
    F64ort vertex_out_;
    template <class Tipg>
    void set(const Tipg& ipg) {
        vertex_in_ = ipg.vertex_in_;
        vertex_out_ = ipg.vertex_out_;
    }
    // for cinsistency of PMM mode
    template <class Tipg, class Tmorton_key>
    void set(const Tipg& ipg, const S32 icut,
             const Tmorton_key& morton_key) {  // used in MakeInteractionList*
        set(ipg);
        // vertex_in_ = ipg.vertex_in_;
    }
    template <typename T>
    void setForExLet(const T& top, const F64vec& s) {
        vertex_in_ = top.vertex_in_.shift(s);
        vertex_out_ = top.vertex_out_.shift(s);
    }
};

template <>
struct TargetBox<SEARCH_MODE_LONG_SCATTER> {
    F64ort vertex_;
    template <class Tipg>
    void set(const Tipg& ipg) {
        vertex_ = ipg.vertex_in_;
    }
    // for cinsistency of PMM mode
    template <class Tipg, class Tmorton_key>
    void set(const Tipg& ipg, const S32 icut,
             const Tmorton_key& morton_key) {  // used in MakeInteractionList*
        // vertex_ = ipg.vertex_in_;
        set(ipg);
    }
    template <typename T>
    void setForExLet(const T& top, const F64vec& s) {
        vertex_ = top.vertex_in_.shift(s);
    }
    void dump(std::ostream& fout = std::cerr) { fout << "vertex_= " << vertex_ << std::endl; }
};

template <>
struct TargetBox<SEARCH_MODE_LONG_CUTOFF> {
    F64ort vertex_;
    template <class Tipg>
    void set(const Tipg& ipg) {
        vertex_ = ipg.vertex_in_;
    }
    // for cinsistency of PMM mode
    template <class Tipg, class Tmorton_key>
    void set(const Tipg& ipg, const S32 icut,
             const Tmorton_key& morton_key) {  // used in MakeInteractionList*
        // vertex_ = ipg.vertex_in_;
        set(ipg);
    }
    template <typename T>
    void setForExLet(const T& top, const F64vec& s) {
        vertex_ = top.vertex_in_.shift(s);
    }
    void dump(std::ostream& fout = std::cerr) { fout << "vertex_= " << vertex_ << std::endl; }
};

// PMMM
template <>
struct TargetBox<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE> {
    F64ort vertex_in_;
    S32ort vertex_out_;
    template <class Tipg, class Tmorton_key>
    void set(const Tipg& ipg, const S32 icut,
             const Tmorton_key& morton_key) {  // used in MakeInteractionList*
        vertex_in_ = ipg.vertex_in_;
        vertex_out_.low_ = morton_key.getPMCellID(vertex_in_.low_) - S32vec(icut);
        vertex_out_.high_ = morton_key.getPMCellID(vertex_in_.high_) + S32vec(icut);
    }
    template <typename T, typename Tmorton_key>
    void setForExLet(const T& top, const F64vec& s, const S32 icut,
                     const Tmorton_key& morton_key) {  // used in FindScatterParticle*
        vertex_in_ = top.vertex_in_.shift(s);
        vertex_out_.low_ = morton_key.getPMCellID(vertex_in_.low_) - S32vec(icut);
        vertex_out_.high_ = morton_key.getPMCellID(vertex_in_.high_) + S32vec(icut);
    }
    template <typename Tmorton_key>
    unsigned int contains(const F64vec& pos, const Tmorton_key& morton_key) const {
        const S32vec pos_tmp = morton_key.getPMCellID(pos);
        return vertex_out_.contains(pos_tmp);
    }
    template <typename Tmorton_key>
    unsigned int contains(const KeyT mkey, const Tmorton_key& morton_key) const {
        // This is faster than the above `contains`. We can use this function
        // only for EPJ and SPJ having a Morton key.
        const S32vec pos_tmp = morton_key.getPMCellID(mkey);
        return vertex_out_.contains(pos_tmp);
    }
    void dump(std::ostream& fout = std::cerr) const {
        fout << "vertex_in_= " << vertex_in_ << std::endl;
        fout << "vertex_out_= " << vertex_out_ << std::endl;
    }
};

template <>
struct TargetBox<SEARCH_MODE_SCATTER> {
    F64ort vertex_in_;
    template <class Tipg>
    void set(const Tipg& ipg) {
        vertex_in_ = ipg.vertex_in_;
    }
    template <class Tipg, class Tmorton_key>
    void set(const Tipg& ipg, const S32 icut, const Tmorton_key& morton_key) {
        set(ipg);
    }
    template <class Tep, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline bool isInEpList(const ReallocatableArray<Tep>& ep_first, const S32 adr, const F64vec& len_peri) const {
        // const auto dis_sq = vertex_in_.getDistanceMinSQ(ep_first[adr].getPos());
        const F64vec ep_pos = ep_first[adr].getPos();
        // const auto dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_in_, ep_pos, len_peri);
        const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_in_, ep_pos, len_peri);
        const F64 r_crit_sq = ep_first[adr].getRSearch() * ep_first[adr].getRSearch();
        return dis_sq <= r_crit_sq;
    }
};

template <>
struct TargetBox<SEARCH_MODE_SYMMETRY> {
    F64ort vertex_out_;
    F64ort vertex_in_;
    template <class Tipg>
    void set(const Tipg& ipg) {
        vertex_out_ = ipg.vertex_out_;
        vertex_in_ = ipg.vertex_in_;
    }
    template <class Tipg, class Tmorton_key>
    void set(const Tipg& ipg, const S32 icut, const Tmorton_key& morton_key) {
        set(ipg);
    }
    template <class Tep, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline bool isInEpList(const ReallocatableArray<Tep>& ep_first, const S32 adr, const F64vec& len_peri) const {
        // const F64 dis_sq = vertex_in_.getDistanceMinSQ(ep_first[adr].getPos());
        const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_in_, ep_first[adr].getPos(), len_peri);
        const F64 r_crit_sq = ep_first[adr].getRSearch() * ep_first[adr].getRSearch();
        return (dis_sq <= r_crit_sq) || IsOverlapped<CALC_DISTANCE_TYPE>(vertex_out_, ep_first[adr].getPos(), len_peri);
    }
};

template <>
struct TargetBox<SEARCH_MODE_GATHER> {
    F64ort vertex_out_;
    template <class Tipg>
    void set(const Tipg& ipg) {
        vertex_out_ = ipg.vertex_out_;
    }
    template <class Tipg, class Tmorton_key>
    void set(const Tipg& ipg, const S32 icut, const Tmorton_key& morton_key) {
        set(ipg);
    }
    template <class Tep, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    inline bool isInEpList(const ReallocatableArray<Tep>& ep_first, const S32 adr, const F64vec& len_peri) const {
        return IsOverlapped<CALC_DISTANCE_TYPE>(vertex_out_, ep_first[adr].getPos(), len_peri);
        // const F64 dis_sq = GetDistanceMinSq<CALC_DISTANCE_TYPE>(vertex_out_, ep_first[adr].getPos(), len_peri);
        // return dis_sq;
        // return vertex_out_.overlapped(ep_first[adr].getPos());
    }
};
// targe box class
////////////

/////////////
// set target box for exching let
template <class TSM>
inline void SetTargetBoxExLet(TargetBox<TSM>& box, const F64ort& in, const F64ort& out) {
    assert(0);
}
template <>
inline void SetTargetBoxExLet<SEARCH_MODE_LONG_SYMMETRY>(TargetBox<SEARCH_MODE_LONG_SYMMETRY>& box, const F64ort& in, const F64ort& out) {
    box.vertex_in_ = in;
    box.vertex_out_ = out;
}
template <>
inline void SetTargetBoxExLet<SEARCH_MODE_LONG_SCATTER>(TargetBox<SEARCH_MODE_LONG_SCATTER>& box, const F64ort& in, const F64ort& out) {
    box.vertex_ = in;
}
template <>
inline void SetTargetBoxExLet<SEARCH_MODE_LONG>(TargetBox<SEARCH_MODE_LONG>& box, const F64ort& in, const F64ort& out) {
    box.vertex_ = in;
}
template <>
inline void SetTargetBoxExLet<SEARCH_MODE_LONG_CUTOFF>(TargetBox<SEARCH_MODE_LONG_CUTOFF>& box, const F64ort& in, const F64ort& out) {
    box.vertex_ = in;
}
// set target box for exching let
/////////////

}  // namespace ParticleSimulator

// #include<tree_unused.hpp>
