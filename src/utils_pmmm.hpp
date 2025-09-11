#pragma once
// This file is included at the end of ps_defs.hpp.
// Classes and functions in this file are used for Particle Mesh Multipole mode.
namespace ParticleSimulator {

template <class Tepi, class Tforce>
class IParticleInfo {
   public:
    S32 n_epi;
    S32 adr_epi_first;
    Tepi* epi_first;
    Tforce* force_first;

    IParticleInfo() {
        n_epi = 0;
        adr_epi_first = -1;
        epi_first = nullptr;
        force_first = nullptr;
    }

    void dump(std::ostream& fout = std::cerr) {
        fout << "n_epi= " << n_epi << " adr_epi_first= " << adr_epi_first << " epi_first= " << epi_first << " force_first= " << force_first
             << std::endl;
    }

    friend std::ostream& operator<<(std::ostream& c, const IParticleInfo& info) {
        c << info.n_epi << " " << info.adr_epi_first << " " << info.epi_first << " " << info.force_first;
        return c;
    }
};

//---- GetWidthOfParticleMeshCell
inline F64vec GetWidthOfParticleMeshCell(const F64ort& pos_unit_cell, const S32vec& n_cell) {
    #if 1
    return pos_unit_cell.getFullLength() / n_cell;
    #else
    F64vec width;
    width.x = pos_unit_cell.getFullLength().x / n_cell.x;
    width.y = pos_unit_cell.getFullLength().y / n_cell.y;
    width.z = pos_unit_cell.getFullLength().z / n_cell.z;
    return width;
    #endif
}

//---- GetCenterOfParticleMeshCell
inline F64vec GetCenterOfParticleMeshCell(const F64ort& pos_unit_cell, const F64vec& width_cell, const S32vec& idx) {
    F64vec pos;
    pos.x = pos_unit_cell.low_.x + (idx.x + 0.5) * width_cell.x;
    pos.y = pos_unit_cell.low_.y + (idx.y + 0.5) * width_cell.y;
    pos.z = pos_unit_cell.low_.z + (idx.z + 0.5) * width_cell.z;
    return pos;
}

inline F64vec GetCenterOfParticleMeshCell(const F64ort& pos_unit_cell, const S32vec& n_cell, const S32vec& idx) {
    F64vec width = GetWidthOfParticleMeshCell(pos_unit_cell, n_cell);
    return GetCenterOfParticleMeshCell(pos_unit_cell, width, idx);
}

//---- GetCellIDMeasuredInUnitCell
#if 1
inline S64vec GetCellIDMeasuredInUnitCell(const F64vec& pos_unit_cell, const F64vec& width_cell, const F64vec& pos) {
    // This function returns a cell index of a PM cell to which
    // a particle whose position is `pos` belongs.
    // The cell index is calculated so that the index (0,0,0)
    // corresponds to the PM cell located at the lower limits
    // of the unit cell (computational box).
    // std::cout<<"pos_unit_cell= "<<pos_unit_cell<<std::endl;
    // std::cout<<"width_cell= "<<width_cell<<std::endl;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    F64vec idx_f;
    idx_f.x = (pos.x - pos_unit_cell.x) / width_cell.x;
    idx_f.y = (pos.y - pos_unit_cell.y) / width_cell.y;
#if 0
        if ( (fabs(idx_f.x) >= (F64)std::numeric_limits<S32>::max()) ||
             (fabs(idx_f.y) >= (F64)std::numeric_limits<S32>::max()) ) {
            std::cout << "pos = " << pos << std::endl;
        }
#endif
    S64vec idx;
    if (idx_f.x >= 0.0)
        idx.x = (S64)idx_f.x;
    else
        idx.x = ((S64)idx_f.x) - 1;
    if (idx_f.y >= 0.0)
        idx.y = (S64)idx_f.y;
    else
        idx.y = ((S64)idx_f.y) - 1;
    return idx;
#else  // PARTICLE_SIMULATOR_TWO_DIMENSION
    F64vec idx_f;
    idx_f.x = (pos.x - pos_unit_cell.x) / width_cell.x;
    idx_f.y = (pos.y - pos_unit_cell.y) / width_cell.y;
    idx_f.z = (pos.z - pos_unit_cell.z) / width_cell.z;
#if 1
    /*
    if ( (fabs(idx_f.x) >= (F64)std::numeric_limits<S32>::max()) ||
         (fabs(idx_f.y) >= (F64)std::numeric_limits<S32>::max()) ||
         (fabs(idx_f.z) >= (F64)std::numeric_limits<S32>::max()) ) {
        //std::cout << "pos = " << pos
        //          << ", rank = " << Comm::getRank()
        //          << std::endl;
        DebugUtils::fouts(Comm::getThreadNum()) << "pos = " << pos << std::endl;
    }
    */
    // std::cout<<"pos= "<<pos<<std::endl;
    // std::cout<<"idx_f= "<<idx_f<<std::endl;
    // std::cout<<"(F64)std::numeric_limits<S64>::max()= "<<(F64)std::numeric_limits<S64>::max()<<std::endl;
    assert(fabs(idx_f.x) < (F64)std::numeric_limits<S64>::max());
    assert(fabs(idx_f.y) < (F64)std::numeric_limits<S64>::max());
    assert(fabs(idx_f.z) < (F64)std::numeric_limits<S64>::max());
#endif
    S64vec idx;
    if (idx_f.x >= 0.0)
        idx.x = (S64)idx_f.x;
    else
        idx.x = ((S64)idx_f.x) - 1;
    if (idx_f.y >= 0.0)
        idx.y = (S64)idx_f.y;
    else
        idx.y = ((S64)idx_f.y) - 1;
    if (idx_f.z >= 0.0)
        idx.z = (S64)idx_f.z;
    else
        idx.z = ((S64)idx_f.z) - 1;

    //std::cout << "idx_f= " << idx_f << " idx= " << idx << std::endl;
    return idx;
#endif  // PARTICLE_SIMULATOR_TWO_DIMENSION
}
inline S64vec GetCellIDMeasuredInUnitCell(const F64ort& pos_unit_cell, const F64vec& width_cell, const F64vec& pos) {
    return GetCellIDMeasuredInUnitCell(pos_unit_cell.low_, width_cell, pos);
}
#else
inline S32vec GetCellIDMeasuredInUnitCell(const F64vec& pos_unit_cell, const F64vec& width_cell, const F64vec& pos) {
    // This function returns a cell index of a PM cell to which
    // a particle whose position is `pos` belongs.
    // The cell index is calculated so that the index (0,0,0)
    // corresponds to the PM cell located at the lower limits
    // of the unit cell (computational box).
    std::cout << "pos_unit_cell= " << pos_unit_cell << std::endl;
    std::cout << "width_cell= " << width_cell << std::endl;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    F64vec idx_f;
    idx_f.x = (pos.x - pos_unit_cell.x) / width_cell.x;
    idx_f.y = (pos.y - pos_unit_cell.y) / width_cell.y;
#if 0
        if ( (fabs(idx_f.x) >= (F64)std::numeric_limits<S32>::max()) ||
             (fabs(idx_f.y) >= (F64)std::numeric_limits<S32>::max()) ) {
            std::cout << "pos = " << pos << std::endl;
        }
#endif
    assert(fabs(idx_f.x) < (F64)std::numeric_limits<S32>::max());
    assert(fabs(idx_f.y) < (F64)std::numeric_limits<S32>::max());
    S32vec idx;
    if (idx_f.x >= 0.0)
        idx.x = (S32)idx_f.x;
    else
        idx.x = ((S32)idx_f.x) - 1;
    if (idx_f.y >= 0.0)
        idx.y = (S32)idx_f.y;
    else
        idx.y = ((S32)idx_f.y) - 1;
    return idx;
#else  // PARTICLE_SIMULATOR_TWO_DIMENSION
    F64vec idx_f;
    idx_f.x = (pos.x - pos_unit_cell.x) / width_cell.x;
    idx_f.y = (pos.y - pos_unit_cell.y) / width_cell.y;
    idx_f.z = (pos.z - pos_unit_cell.z) / width_cell.z;
#if 0
        if ( (fabs(idx_f.x) >= (F64)std::numeric_limits<S32>::max()) ||
             (fabs(idx_f.y) >= (F64)std::numeric_limits<S32>::max()) ||
             (fabs(idx_f.z) >= (F64)std::numeric_limits<S32>::max()) ) {
            //std::cout << "pos = " << pos
            //          << ", rank = " << Comm::getRank()
            //          << std::endl;
            DebugUtils::fouts(Comm::getThreadNum()) << "pos = " << pos << std::endl;
        }
#endif
    assert(fabs(idx_f.x) < (F64)std::numeric_limits<S32>::max());
    assert(fabs(idx_f.y) < (F64)std::numeric_limits<S32>::max());
    assert(fabs(idx_f.z) < (F64)std::numeric_limits<S32>::max());
    S32vec idx;
    if (idx_f.x >= 0.0)
        idx.x = (S32)idx_f.x;
    else
        idx.x = ((S32)idx_f.x) - 1;
    if (idx_f.y >= 0.0)
        idx.y = (S32)idx_f.y;
    else
        idx.y = ((S32)idx_f.y) - 1;
    if (idx_f.z >= 0.0)
        idx.z = (S32)idx_f.z;
    else
        idx.z = ((S32)idx_f.z) - 1;
    return idx;
#endif  // PARTICLE_SIMULATOR_TWO_DIMENSION
}

inline S32vec GetCellIDMeasuredInUnitCell(const F64ort& pos_unit_cell, const F64vec& width_cell, const F64vec& pos) {
    return GetCellIDMeasuredInUnitCell(pos_unit_cell.low_, width_cell, pos);
}
#endif

//---- GetBoxOfParticleMeshCell
inline F64ort GetBoxOfParticleMeshCell(const S32vec& idx, const F64vec& pos_unit_cell, const F64vec& width_cell) {
    F64ort box;
    box.low_.x = pos_unit_cell.x + width_cell.x * idx.x;
    box.low_.y = pos_unit_cell.y + width_cell.y * idx.y;
    box.low_.z = pos_unit_cell.z + width_cell.z * idx.z;
    box.high_.x = pos_unit_cell.x + width_cell.x * (idx.x + 1);
    box.high_.y = pos_unit_cell.y + width_cell.y * (idx.y + 1);
    box.high_.z = pos_unit_cell.z + width_cell.z * (idx.z + 1);
    return box;
}

}  // namespace ParticleSimulator
