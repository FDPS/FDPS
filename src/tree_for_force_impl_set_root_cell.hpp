#pragma once
namespace ParticleSimulator {

/*
  setRootCell:  determine following parameters
  - pos_root_cell_ // position of the root cell
  - center_ // center of the root cell
  - lenght_ // (full) length of the root cell 
  - pos_unit_cell_ is a computational domain, which is the same as pos_root_cell_ only if OPEN
  BOUNDAYR, othewise (i.e. PERIODIC BOUNDARY) pos_root_domain.
  - idx_pm_cell_ : index of the unit_cell from pos_root_cell
  - lev_pm_cell_
  - width_pm_cell_
    root_cell includes unit_cell
 */

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setRootCell(const DomainInfo& dinfo) {
#if defined(PS_DEBUG_PRINT_setRootCell)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START setRootCell(const DomainInfo&)", 0);
#endif
    const F64 time_offset = GetWtime();
    // calcCenterAndLengthOfRootCell(typename TSM::search_type(), dinfo);
    calcCenterAndLengthOfRootCell(dinfo);
    if constexpr (std::is_same_v<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>) {
        pos_unit_cell_ = dinfo.getPosRootDomain();
    }

    /*
    if (dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN && dinfo.getFlagSetPosRootDomain() == false) {
        pos_unit_cell_ = pos_root_cell_;
    } else {
        pos_unit_cell_ = p_domain_info_->getPosRootDomain();
        // pos_unit_cell_ = root_domain;
    }
    */
#if defined(PS_DEBUG_PRINT_setRootCell)
    PARTICLE_SIMULATOR_PRINT_MSG("**** before repositionPosRootCellToEmbedParticleMesh", 0);
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE(pos_root_cell_);
        PS_PRINT_VARIABLE(pos_unit_cell_);
        PS_PRINT_VARIABLE2(center_, length_);
    }
#endif
    repositionPosRootCellToEmbedParticleMesh();
#if defined(PS_DEBUG_PRINT_setRootCell)
    PARTICLE_SIMULATOR_PRINT_MSG("**** after repositionPosRootCellToEmbedParticleMesh", 0);
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE(pos_root_cell_);
        PS_PRINT_VARIABLE(pos_unit_cell_);
        PS_PRINT_VARIABLE2(center_, length_);
    }
    // exit(1);
#endif
    time_profile_.set_root_cell += GetWtime() - time_offset;
#if defined(PS_DEBUG_PRINT_setRootCell)
    PARTICLE_SIMULATOR_PRINT_MSG("**** Root cell infomation", 0);
    if (comm_info_.getRank() == 0) {
        std::cout << "pos_unit_cell_ = " << pos_unit_cell_ << std::endl;
        std::cout << "idx_pm3_unit_cell_ = " << idx_pm3_unit_cell_ << std::endl;
        std::cout << "lev_pm_cell_     = " << lev_pm_cell_ << std::endl;
        std::cout << "width_pm_cell_   = " << width_pm_cell_ << std::endl;
        std::cout << "pos_root_cell_   = " << pos_root_cell_ << std::endl;
        std::cout << "pos_root_cell_.getFullLength() = " << pos_root_cell_.getFullLength() << std::endl;
        std::cout << "length_ " << length_ << std::endl;
        std::cout << "center_ " << center_ << std::endl;
    }
    PARTICLE_SIMULATOR_PRINT_MSG("**** END setRootCell(const DomainInfo&)", 0);
    // exit(1);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::setRootCell(const F64 l, const F64vec& c) {
#if defined(PS_DEBUG_PRINT_setRootCell)
    PARTICLE_SIMULATOR_PRINT_MSG("**** START setRootCell(const F64, const F64vec&)", 0);
#endif
    const F64 time_offset = GetWtime();
    center_ = c;
    length_ = l;
    pos_root_cell_.low_ = center_ - F64vec(length_ * 0.5);
    pos_root_cell_.high_ = center_ + F64vec(length_ * 0.5);
    if (p_domain_info_->getBoundaryCondition() == BOUNDARY_CONDITION_OPEN) {
        pos_unit_cell_ = pos_root_cell_;
    } else {
        pos_unit_cell_ = p_domain_info_->getPosRootDomain();
    }
    time_profile_.set_root_cell += GetWtime() - time_offset;
#if defined(PS_DEBUG_PRINT_setRootCell)
    if (comm_info_.getRank() == 0) {
        std::cout << "pos_unit_cell_ = " << pos_unit_cell_ << std::endl;
        std::cout << "idx_pm3_unit_cell_ = " << idx_pm3_unit_cell_ << std::endl;
        std::cout << "lev_pm_cell_     = " << lev_pm_cell_ << std::endl;
        std::cout << "width_pm_cell_   = " << width_pm_cell_ << std::endl;
        std::cout << "pos_root_cell_   = " << pos_root_cell_ << std::endl;
        std::cout << "pos_root_cell_.getFullLength() = " << pos_root_cell_.getFullLength() << std::endl;
        std::cout << "length_ " << length_ << std::endl;
        std::cout << "center_ " << center_ << std::endl;
    }
    PARTICLE_SIMULATOR_PRINT_MSG("**** END setRootCell(const F64, const F64vec&)", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcCenterAndLengthOfRootCell(const DomainInfo& dinfo) {
    if constexpr (std::is_same_v<typename TSM::calc_center_and_length_type(), TagCalcCenterAndLengthUseEPJ>) {
        calcCenterAndLengthOfRootCellImpl(epj_org_.getPointer(), dinfo);
    } else {
        calcCenterAndLengthOfRootCellImpl(epi_org_.getPointer(), dinfo);
    }
}
/*
template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcCenterAndLengthOfRootCell(const DomainInfo& dinfo,
                                                                                                                      const F64ort& outer_boundary) {
    if constexpr (std::is_same_v<TSM::calc_center_and_length_type(), TagCalcCenterAndLengthUseEPJ>) {
        calcCenterAndLengthOfRootCellImpl(epj_org_.getPointer(), dinfo, outer_boundary);
    } else {
        calcCenterAndLengthOfRootCellImpl(epi_org_.getPointer(), dinfo, outer_boundary);
    }
}
*/

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tep2>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcCenterAndLengthOfRootCellImpl(const Tep2 ep[],
                                                                                                                          const DomainInfo& dinfo) {
#if defined(PS_DEBUG_PRINT_calcCenterAndLengthOfRootCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START calcCenterAndLengthOfRootCellImpl", 0);
#endif
    // determine outer boundary
    F64ort box_loc;
    box_loc.initNegativeVolume();
    PS_OMP_PARALLEL {
        F64ort box_loc_tmp;
        box_loc_tmp.init();
        PS_OMP(omp for nowait)
        for (S32 ip = 0; ip < n_loc_tot_; ip++) {
            box_loc_tmp.merge(ep[ip].getPos(), GetMyRSearch(ep[ip]) * 1.000001);
        }
        PS_OMP(omp critical) { box_loc.merge(box_loc_tmp); }
    }
    F64vec xlow_loc = box_loc.low_;
    F64vec xhigh_loc = box_loc.high_;
    F64vec xlow_glb = comm_info_.getMinValue(xlow_loc);
    F64vec xhigh_glb = comm_info_.getMaxValue(xhigh_loc);
    const F64ort outer_boundary(xlow_glb, xhigh_glb);
#if defined(PS_DEBUG_PRINT_calcCenterAndLengthOfRootCellImpl)
    if (Comm::getRank() == 0) {
        PS_PRINT_VARIABLE(outer_boundary);
    }
    PARTICLE_SIMULATOR_PRINT_MSG("***** before calcCenterAndLengthOfRootCellImpl", 0);
#endif
    calcCenterAndLengthOfRootCellImpl(ep, dinfo, outer_boundary);
#if defined(PS_DEBUG_PRINT_calcCenterAndLengthOfRootCellImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** END calcCenterAndLengthOfRootCellImpl", 0);
#endif
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
template <typename Tep2>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::calcCenterAndLengthOfRootCellImpl(
    const Tep2 ep[], const DomainInfo& dinfo, const F64ort& outer_boundary) {
    const F64ort root_domain = dinfo.getPosRootDomain();
    const F64vec shift = root_domain.high_ - root_domain.low_;
    bool pa[DIMENSION_LIMIT];
    dinfo.getPeriodicAxis(pa);
    S32 num_shift_p[DIMENSION_LIMIT];
    S32 num_shift_n[DIMENSION_LIMIT];
    S32 num_shift_max[DIMENSION_LIMIT];
    for (S32 cid = 0; cid < DIMENSION; cid++) {
        if (pa[cid]) {
            num_shift_p[cid] = num_shift_n[cid] = 0;
            F64vec shift_tmp(0.0);
            shift_tmp[cid] = shift[cid];
            CalcNumShift(root_domain, outer_boundary,  shift_tmp, num_shift_p[cid]); // @util.hpp
            CalcNumShift(root_domain, outer_boundary, -shift_tmp, num_shift_n[cid]);
            num_shift_max[cid] = std::max(num_shift_p[cid], num_shift_n[cid]);
        } else {
            num_shift_max[cid] = 0;
        }
    }
    /*
    std::cout<<"root_domain= "<<root_domain<<std::endl;
    std::cout<<"outer_boundary= "<<outer_boundary<<std::endl;
    std::cout<<"shift= "<<shift<<std::endl;
    std::cout<<"num_shift_p="<<num_shift_p[0]<<" "<<num_shift_p[1]<<" "<<num_shift_p[2]<<std::endl;
    std::cout<<"num_shift_n="<<num_shift_n[0]<<" "<<num_shift_n[1]<<" "<<num_shift_n[2]<<std::endl;
    std::cout<<"num_shift_max="<<num_shift_max[0]<<" "<<num_shift_max[1]<<" "<<num_shift_max[2]<<std::endl;
    */
    length_ = 0.0;
    for (S32 cid = 0; cid < DIMENSION; cid++) {
        if (pa[cid]) {
            F64 length_tmp = (2 * num_shift_max[cid] + 1) * shift[cid];
            if (length_tmp > length_) length_ = length_tmp;
            center_[cid] = root_domain.getCenter()[cid];
        } else {
            F64 length_tmp = outer_boundary.getFullLength()[cid];
            if (length_tmp > length_) length_ = length_tmp;
            center_[cid] = outer_boundary.getCenter()[cid];
        }
    }
    length_ *= 1.000001;
    pos_root_cell_.low_ = center_ - F64vec(length_ * 0.5);
    pos_root_cell_.high_ = center_ + F64vec(length_ * 0.5);
    //Comm::barrier();
    //exit(1);
    //std::cout << "************* outer_boundary= " << outer_boundary << std::endl;
}





template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::repositionPosRootCellToEmbedParticleMesh() {
    if constexpr (TSM::reposition_pos_root_cell_type == RepositionPosRootCellType::PMMM) {
        PS_COMPILE_TIME_MESSAGE("RepositionPosRootCellType::PMMM@repositionPosRootCellToEmbedParticleMesh");
        repositionPosRootCellToEmbedParticleMeshImpl();
    }
}

template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE>::repositionPosRootCellToEmbedParticleMeshImpl() {
#if defined(PS_DEBUG_PRINT_repositionPosRootCellToEmbedParticleMeshImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START repositionPosRootCellToEmbedParticleMeshImpl", 0);
#endif
    // Compute cell widths
    width_pm_cell_ = GetWidthOfParticleMeshCell(pos_unit_cell_, n_cell_pm3_);
#if defined(PS_DEBUG_PRINT_repositionPosRootCellToEmbedParticleMeshImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** in repositionPosRootCellToEmbedParticleMeshImpl", 0);
    std::cout << "width_pm_cell_= " << width_pm_cell_ << std::endl;
    std::cout << "1) pos_root_cell_= " << pos_root_cell_ << std::endl;
#endif
    // Compute tree level
    S32vec n_cell_margin, n_cell_tot_reqmin, level;
    n_cell_margin.x = (i_cut_ / n_cell_pm3_.x + 1) * n_cell_pm3_.x;
    n_cell_margin.y = (i_cut_ / n_cell_pm3_.y + 1) * n_cell_pm3_.y;
    n_cell_margin.z = (i_cut_ / n_cell_pm3_.z + 1) * n_cell_pm3_.z;
    n_cell_tot_reqmin.x = n_cell_pm3_.x + 2 * n_cell_margin.x;
    n_cell_tot_reqmin.y = n_cell_pm3_.y + 2 * n_cell_margin.y;
    n_cell_tot_reqmin.z = n_cell_pm3_.z + 2 * n_cell_margin.z;

#if defined(PS_DEBUG_PRINT_repositionPosRootCellToEmbedParticleMeshImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** in repositionPosRootCellToEmbedParticleMeshImpl", 0);
    std::cout << "n_cell_margin= " << n_cell_margin << std::endl;
    std::cout << "n_cell_tot_reqmin= " << n_cell_tot_reqmin << std::endl;
#endif

    // [Notes]
    // (1) In order to get a correct PM cell ID in the exchange LET,
    //     we must require that a shifted domain is fully contained
    //     in pos_root_cell.
    // (2) The following is an old implementation:
    //     n_cell_tot_reqmin.x = n_cell_.x + 2*icut_;
    //     n_cell_tot_reqmin.y = n_cell_.y + 2*icut_;
    //     n_cell_tot_reqmin.z = n_cell_.z + 2*icut_;

    F64vec length;
    for (auto i = 0; i < DIMENSION; i++) {
        length[i] = width_pm_cell_[i] * n_cell_tot_reqmin[i];
        while (length[i] < pos_root_cell_.getFullLength()[i]) {
            n_cell_tot_reqmin[i]++;
            length[i] = width_pm_cell_[i] * n_cell_tot_reqmin[i];
        }
        //level[i] = getMinExpOfTwoGT(n_cell_tot_reqmin[i]);
        level[i] = GetMinExpOfTwoGT(n_cell_tot_reqmin[i]); // i.e. (log2(n_cell_tot_reqmin)) @util.hpp
    }

#if defined(PS_DEBUG_PRINT_repositionPosRootCellToEmbedParticleMeshImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** in repositionPosRootCellToEmbedParticleMeshImpl", 0);
    std::cout << "length= " << length << std::endl;
    std::cout << "level= " << level << std::endl;
#endif

    lev_pm_cell_ = MyMax(level.x, level.y, level.z);
    lev_leaf_limit_ = lev_group_limit_ = lev_pm_cell_;

#if defined(PS_DEBUG_PRINT_repositionPosRootCellToEmbedParticleMeshImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** in repositionPosRootCellToEmbedParticleMeshImpl", 0);
    std::cout << "lev_pm_cell= " << lev_pm_cell_ << std::endl;
    std::cout << "lev_leaf_limit= " << lev_leaf_limit_ << std::endl;
#endif

    S32vec n_cell_tot; // the QUANTIZED number of PM cells with margin
    n_cell_tot.x = 1 << lev_pm_cell_;  // 2^{lev_pm_cell_}
    n_cell_tot.y = 1 << lev_pm_cell_;
    n_cell_tot.z = 1 << lev_pm_cell_;
#if defined(PS_DEBUG_PRINT_repositionPosRootCellToEmbedParticleMeshImpl)
    PARTICLE_SIMULATOR_PRINT_MSG("***** in repositionPosRootCellToEmbedParticleMeshImpl", 0);
    std::cout << "n_cell_tot= " << n_cell_tot << std::endl;
#endif
    // Compute idx_unit_cell_, pos_root_cell_
    idx_pm3_unit_cell_.x = (n_cell_tot.x - n_cell_tot_reqmin.x) / 2 + n_cell_margin.x; // the index of the unit cell in the PM3 mesh
    idx_pm3_unit_cell_.y = (n_cell_tot.y - n_cell_tot_reqmin.y) / 2 + n_cell_margin.y;
    idx_pm3_unit_cell_.z = (n_cell_tot.z - n_cell_tot_reqmin.z) / 2 + n_cell_margin.z;
    assert(idx_pm3_unit_cell_.x > (i_cut_ - 1));  // this line assert there are at least i_cut_ cells to the left of the cell with idx_pm3_unit_cell_
    assert(idx_pm3_unit_cell_.y > (i_cut_ - 1));
    assert(idx_pm3_unit_cell_.z > (i_cut_ - 1));
    pos_root_cell_.low_.x = pos_unit_cell_.low_.x - width_pm_cell_.x * idx_pm3_unit_cell_.x; // the position of the root cell
    pos_root_cell_.low_.y = pos_unit_cell_.low_.y - width_pm_cell_.y * idx_pm3_unit_cell_.y;
    pos_root_cell_.low_.z = pos_unit_cell_.low_.z - width_pm_cell_.z * idx_pm3_unit_cell_.z;
    pos_root_cell_.high_.x = pos_root_cell_.low_.x + width_pm_cell_.x * n_cell_tot.x;
    pos_root_cell_.high_.y = pos_root_cell_.low_.y + width_pm_cell_.y * n_cell_tot.y;
    pos_root_cell_.high_.z = pos_root_cell_.low_.z + width_pm_cell_.z * n_cell_tot.z;

    // is it correct that non-cubic tree cell is enabled?
    // Recalculate length_ & center_
    const F64vec lenvec = pos_root_cell_.getFullLength();
    // length_ = std::max(std::max(lenvec.x, lenvec.y), lenvec.z);
    length_ = MyMax(lenvec.x, lenvec.y, lenvec.z);
    center_ = pos_root_cell_.getCenter();

#if defined(PS_DEBUG_PRINT_repositionPosRootCellToEmbedParticleMeshImpl)
    if (comm_info_.getRank() == 0) {
        std::cout << "pos_unit_cell_ = " << pos_unit_cell_ << std::endl;
        std::cout << "idx_pm3_unit_cell_ = " << idx_pm3_unit_cell_ << std::endl;
        std::cout << "(n_cell_tot      = " << n_cell_tot << ")" << std::endl;
        std::cout << "lev_pm_cell_     = " << lev_pm_cell_ << std::endl;
        std::cout << "width_pm_cell_   = " << width_pm_cell_ << std::endl;
        std::cout << "pos_root_cell_   = " << pos_root_cell_ << std::endl;
        std::cout << "pos_root_cell_.getFullLength() = " << pos_root_cell_.getFullLength() << std::endl;
        std::cout << "length_ " << length_ << std::endl;
        std::cout << "center_ " << center_ << std::endl;
    }
    PARTICLE_SIMULATOR_PRINT_MSG("***** END repositionPosRootCellToEmbedParticleMeshImpl", 0);
#endif
}

}  // namespace ParticleSimulator
