#include <set>
//enum class DistributionType { RANDOM, REGULAR };
//enum class InteractionType { GRAVITY, COULOMB };

template <typename Tp>
void shift_position_impl(Tp &ptcl, const PS::F64vec &shift, const PS::F64ort &pos_root) {
    PS::F64vec len_root = pos_root.getFullLength();
    PS::F64vec pos_new = ptcl.pos + shift;
    while (pos_new.x < pos_root.low_.x) {
        pos_new.x += len_root.x;
    }
    while (pos_new.x >= pos_root.high_.x) {
        pos_new.x -= len_root.x;
    }
    if (pos_new.x == pos_root.high_.x) {
        pos_new.x = pos_root.low_.x;
    }
    while (pos_new.y < pos_root.low_.y) {
        pos_new.y += len_root.y;
    }
    while (pos_new.y >= pos_root.high_.y) {
        pos_new.y -= len_root.y;
    }
    if (pos_new.y == pos_root.high_.y) {
        pos_new.y = pos_root.low_.y;
    }
    while (pos_new.z < pos_root.low_.z) {
        pos_new.z += len_root.z;
    }
    while (pos_new.z >= pos_root.high_.z) {
        pos_new.z -= len_root.z;
    }
    if (pos_new.z == pos_root.high_.z) {
        pos_new.z = pos_root.low_.z;
    }
    ptcl.pos = pos_new;
}

template <typename Tp>
void shift_position(Tp ptcl[], const int n, const PS::F64vec shift, const PS::F64ort &pos_root) {
    // const PS::F64vec len_root = pos_root.getFullLength();
    for (PS::S32 i = 0; i < n; i++) {
        shift_position_impl(ptcl[i], shift, pos_root);
    }
}

template <typename Tp, typename Tparam>
void setup_particle(Tp &ptcl, const Tparam &params, const PS::F64ort &pos_unit_cell, const PS::S32 seed = 0) {
    auto my_rank = PS::Comm::getRank();
    auto n_proc = PS::Comm::getNumberOfProc();
    auto n_ptcl_glb = params.n_ptcl_glb;
    auto [n_disp, n_disp_next] = PS::CalcAdrToSplitData(my_rank, n_proc, n_ptcl_glb);
    auto n_ptcl_loc = n_disp_next - n_disp;
    // std::cerr<<"n_disp= "<<n_disp<<" n_disp_next= "<<n_disp_next<<" n_ptcl_loc= "<<n_ptcl_loc<<std::endl;
    ptcl.setNumberOfParticleLocal(n_ptcl_loc);
    PS::MT::init_genrand(seed);
    for (PS::S32 i = 0; i < n_disp; i++) {
        for (PS::S32 j = 0; j < 3; j++) {
            // auto x = PS::MT::genrand_res53();  // discard the first n_disp random numbers
            PS::MT::genrand_res53();  // discard the first n_disp random numbers
        }
    }
    // std::cerr<<"B) check 1"<<std::endl;
    // PS::F64vec center = (pos_unit_cell.high_ + pos_unit_cell.low_) * 0.5;
    PS::F64vec pos_origin_unit_cell = pos_unit_cell.low_;
    PS::F64vec len_unit_cell = pos_unit_cell.high_ - pos_unit_cell.low_;
    constexpr PS::F64 m_tot = 1.0;
    for (PS::S32 i = 0; i < n_ptcl_loc; i++) {
        // ptcl[i].mass = m_tot / n_ptcl_glb;
        // ptcl[i].charge = ptcl[i].mass;
        ptcl[i].charge = m_tot / n_ptcl_glb;
        ptcl[i].id = i + n_disp;
        ptcl[i].pos[0] = pos_origin_unit_cell.x + len_unit_cell.x * PS::MT::genrand_res53();
        ptcl[i].pos[1] = pos_origin_unit_cell.y + len_unit_cell.y * PS::MT::genrand_res53();
        ptcl[i].pos[2] = pos_origin_unit_cell.z + len_unit_cell.z * PS::MT::genrand_res53();
        ptcl[i].vel = 0.0;
    }
    
//#if defined(COULOMB)
#if 1
    auto n_ptcl_negative = n_ptcl_glb * params.negative_charge_ratio;
    std::set<PS::S32> id_negative;
    std::mt19937 mt(seed);
    std::uniform_int_distribution<PS::S32> dist(0, n_ptcl_glb);
    while (id_negative.size() < n_ptcl_negative) {
        id_negative.insert(dist(mt));
    }
    auto head = id_negative.lower_bound(ptcl[0].id);
    auto end = id_negative.upper_bound(ptcl[n_ptcl_loc - 1].id);
    for (auto it = head; it != end; it++) {
        // ptcl[*it].charge *= -1.0;
        ptcl[*it - ptcl[0].id].charge *= -1.0;
    }
#endif
}
