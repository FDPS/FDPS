#pragma once
// In this file, we implement utility functions used for TreeForForce
// for particle mesh multipole method.
#include <unordered_map>

namespace ParticleSimulator {

template <class Ttc>
inline U32 GetOpenBit(const ReallocatableArray<Ttc> &tc_first, const S32 adr_tc, const ReallocatableArray<F64ort> &target_box, const F64 r_crit_sq) {
    U32 open_bit = 0;
    for (S32 i = 0; i < N_CHILDREN; i++) {
        const F64vec pos = tc_first[adr_tc + i].mom_.getPos();
        bool too_close = false;
        for (S32 k = 0; k < target_box.size(); k++) {
            const F64 dis = target_box[k].getDistanceMinSq(pos);
            too_close = too_close || (dis <= r_crit_sq);
        }
        open_bit |= (too_close << i);
        open_bit |= (tc_first[adr_tc + i].n_ptcl_ > 0) << (i + N_CHILDREN);  // if true, it should be checked
    }
    return open_bit;
}

template <class TSM, class Ttc, class Ttp, class Tep, class Tsp>
inline void MakeList(const ReallocatableArray<Ttc> &tc_first, const S32 adr_tc, const ReallocatableArray<Ttp> &tp_first,
                     const ReallocatableArray<Tep> &ep_first, ReallocatableArray<S32> &adr_ep_list, const ReallocatableArray<Tsp> &sp_first,
                     ReallocatableArray<S32> &adr_sp_list, const ReallocatableArray<F64ort> &target_box, const F64 r_crit_sq, const S32 n_leaf_limit,
                     const S32 lev_leaf_limit, const S32 adr_tree_sp_first) {
    const Ttc *tc_cur = tc_first.getPointer(adr_tc);
    const S32 n_ptcl = tc_cur->n_ptcl_;
    const S32 adr_ptcl = tc_cur->adr_ptcl_;
    const S32 adr_tc_child = tc_cur->adr_tc_;
    if (!(tc_cur->isLeafPMM(n_leaf_limit, lev_leaf_limit))) {  // not leaf
        U32 open_bit = GetOpenBit(tc_first, adr_tc_child, target_box, r_crit_sq * 0.25);
        for (S32 i = 0; i < N_CHILDREN; i++) {
            if (!((open_bit >> (i + N_CHILDREN)) & 0x1))
                continue;
            else if ((open_bit >> i) & 0x1) {  // close
                MakeList<TSM, Ttc, Ttp, Tep, Tsp>(tc_first, adr_tc_child + i, tp_first, ep_first, adr_ep_list, sp_first, adr_sp_list, target_box,
                                                  r_crit_sq * 0.25, n_leaf_limit, lev_leaf_limit, adr_tree_sp_first);
            } else {  // far
                const S32 adr_sp = adr_tree_sp_first + adr_tc_child + i;
                adr_sp_list.push_back(adr_sp);
            }
        }
    } else {  // leaf
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        for (S32 ip = 0; ip < n_ptcl; ip++, cnt_adr_ptcl++) {
            if (GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0) {
                const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                adr_ep_list.pushBackNoCheck(adr_ep);
            } else {
                const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                adr_sp_list.pushBackNoCheck(adr_sp);
            }
        }
    }
}

template <class Ttc>
inline bool IsOpen(const ReallocatableArray<Ttc> &tc, const S32 adr_tc, const ReallocatableArray<F64ort> &target_box, const F64 r_crit_sq) {
    bool too_close = false;
    const F64vec pos = tc[adr_tc].mom_.getPos();
    for (S32 i = 0; i < target_box.size(); i++) {
        const F64 dis = target_box[i].getDistanceMinSq(pos);
        too_close = too_close || (dis <= r_crit_sq);
    }
    return too_close;
}

///////////
// CalcAdrTCFromParticleMeshCellIndex
///////////
template <class Ttc, class Ttp, class Tep, class Tsp, class Tmorton_key>
inline void CalcAdrTCFromParticleMeshCellIndex(S32 adr_tc_level_partition[], const ReallocatableArray<Ttc> &tc, const ReallocatableArray<Ttp> &tp,
                                               const ReallocatableArray<Tep> &ep, const ReallocatableArray<Tsp> &sp, const S32vec &n_cell,
                                               const S32 lev_pm_cell, const Tmorton_key &morton_key,
                                               std::unordered_map<S32, S32> &adr_tc_from_pm_cell_idx) {
#if defined(PS_DEBUG_PRINT_CalcAdrTCFromParticleMeshCellIndex)
    PARTICLE_SIMULATOR_PRINT_MSG("***** START CalcAdrTCFromParticleMeshCellIndex *****", 0);
#endif
    const S32 head = adr_tc_level_partition[lev_pm_cell];
    const S32 next = adr_tc_level_partition[lev_pm_cell + 1];
    for (S32 adr_tc = head; adr_tc < next; adr_tc++) {
        // std::cerr<<"adr_tc= "<<adr_tc<<std::endl;
        const Ttc *tc_tmp = tc.getPointer(adr_tc);
        const S32 n_ptcl = tc_tmp->n_ptcl_;
        // std::cerr<<"n_ptcl= "<<n_ptcl<<std::endl;
        const U32 adr_tp = tc_tmp->adr_ptcl_;
        // std::cerr<<"adr_tp= "<<adr_tp<<std::endl;
        if (n_ptcl == 0)
            continue;
        else {  // PM cell having particles
            // std::cerr<<"tp.size()= "<<tp.size()<<std::endl;
            const KeyT mkey = tp[adr_tp].key_;
            // std::cerr<<"mkey= "<<mkey<<std::endl;
            const S32vec idx = morton_key.getPMCellID(mkey);
            // std::cerr<<"idx "<<idx<<std::endl;
            if (idx.x < 0 || n_cell.x <= idx.x || idx.y < 0 || n_cell.y <= idx.y || idx.z < 0 || n_cell.z <= idx.z) continue;
            const S32 idx_1d = idx.x + n_cell.x * (idx.y + n_cell.y * idx.z);
            adr_tc_from_pm_cell_idx[idx_1d] = adr_tc;

#if PS_DEBUG_PRINT_CalcAdrTCFromParticleMeshCellIndex >= 2

            // for debug
            if (idx_1d == 0) {
                std::cout << "rank = " << Comm::getRank() << " mkey = " << mkey << " idx.x = " << idx.x << " idx.y = " << idx.y
                          << " idx.z = " << idx.z << std::endl;
                for (S32 i = 0; i < n_ptcl; i++) {
                    const U32 adr_ptcl = tp[adr_tp + i].adr_ptcl_;
                    if (GetMSB(adr_ptcl) == 0) {  // ep
                        std::cout << "EPJ: (" << Comm::getRank() << ") i = " << i << " " << ep[adr_ptcl].pos.x << " " << ep[adr_ptcl].pos.y << " "
                                  << ep[adr_ptcl].pos.z << " " << ep[adr_ptcl].mass << std::endl;
                    } else {  // sp
                        const U32 adr = ClearMSB(adr_ptcl);
                        std::cout << "SPJ: (" << Comm::getRank() << ") i = " << i << " " << sp[adr].getPos().x << " " << sp[adr].getPos().y << " "
                                  << sp[adr].getPos().z << " " << sp[adr].getCharge() << std::endl;
                    }
                }
            }
#endif
        }
    }
#if defined(PS_DEBUG_PRINT_CalcAdrTCFromParticleMeshCellIndex)
    PARTICLE_SIMULATOR_PRINT_MSG("END CalcAdrTCFromParticleMeshCellIndex", 0);
#endif
}

///////////
// CalcTotalChargeAndDispersionOfParticleMeshCell*
///////////
template <class Ttc, class Ttp, class Tep>
inline void CalcTotalChargeAndDispersionOfParticleMeshCellFull(const S32 adr_tc, const ReallocatableArray<Ttc> &tc, const ReallocatableArray<Ttp> &tp,
                                                               const ReallocatableArray<Tep> &ep, const F64vec &center, const F64ort &pos_my_domain,
                                                               F64 &msum, F64 &quad0) {
    const Ttc *tc_tmp = tc.getPointer(adr_tc);
    const S32 n_ptcl = tc_tmp->n_ptcl_;
    const U32 adr_ptcl = tc_tmp->adr_ptcl_;
    if (n_ptcl > 0) {  // PM cell having particles
        ReallocatableArray<S32> adr_ep_list;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        // Make a list of array index of ep
        S32 cnt_adr_ptcl = adr_ptcl;
        for (S32 ip = 0; ip < n_ptcl; ip++, cnt_adr_ptcl++) {
            if (GetMSB(tp[cnt_adr_ptcl].adr_ptcl_) == 0) {  // EPJ
                const S32 adr_ep = tp[cnt_adr_ptcl].adr_ptcl_;
                adr_ep_list.pushBackNoCheck(adr_ep);
            }
        }
        // Calculate total charge and dispersion using the local EPJs
        msum = quad0 = 0;
        for (S32 i = 0; i < adr_ep_list.size(); i++) {
            const S32 adr = adr_ep_list[i];
            const F64 charge = ep[adr].getCharge();
            const F64vec pos = ep[adr].getPos();
            if (pos_my_domain.contains(pos)) {  // check if the pos is in the local domain or not
                const F64vec dr = pos - center;
                msum += charge;
                quad0 += charge * (dr * dr);
            }
        }
        // Free memory
        adr_ep_list.freeMem();
    }
}

template <class Ttc, class Ttp, class Tep>
inline void CalcTotalChargeAndDispersionOfParticleMeshCellStripe(const S32 adr_tc, const ReallocatableArray<Ttc> &tc,
                                                                 const ReallocatableArray<Ttp> &tp, const ReallocatableArray<Tep> &ep,
                                                                 const F64vec &center, const F64ort &pos_my_domain, F64 &msum, F64 &quad0,
                                                                 const S32 n_thread, const S32 thread_id) {
    const Ttc *tc_tmp = tc.getPointer(adr_tc);
    const S32 n_ptcl = tc_tmp->n_ptcl_;
    const U32 adr_ptcl = tc_tmp->adr_ptcl_;
    if (n_ptcl > 0) {  // PM cell having particles
        const S32 cap = (n_ptcl / n_thread) + 1;
        ReallocatableArray<S32> adr_ep_list;
        adr_ep_list.reserveEmptyAreaAtLeast(cap);
        // Calculate the range of array index for thread_id
        const S32 head = (n_ptcl / n_thread) * thread_id + std::min(n_ptcl % n_thread, thread_id);
        const S32 tail = (n_ptcl / n_thread) * (thread_id + 1) + std::min(n_ptcl % n_thread, (thread_id + 1));
        // Make a list of array index of ep
        for (S32 i = head; i < tail; i++) {
            const S32 adr_tp = adr_ptcl + i;
            if (GetMSB(tp[adr_tp].adr_ptcl_) == 0) {  // EPJ
                const S32 adr_ep = tp[adr_tp].adr_ptcl_;
                adr_ep_list.pushBackNoCheck(adr_ep);
            }
        }
        // Calculate total charge and dispersion using the local EPJs
        msum = quad0 = 0;
        for (S32 i = 0; i < adr_ep_list.size(); i++) {
            const S32 adr = adr_ep_list[i];
            const F64 charge = ep[adr].getCharge();
            const F64vec pos = ep[adr].getPos();
            if (pos_my_domain.contains(pos)) {  // check if local or not
                const F64vec dr = pos - center;
                msum += charge;
                quad0 += charge * (dr * dr);
            }
        }
        // Free memory
        adr_ep_list.freeMem();
    }
}

///////////
// CalcMultipoleMomentOfParticleMeshCell
///////////
#if 0
// NO TREE WALK VERSION
// In this version, we just use local MM for the tree code to obtain that for the PM cell.
template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, class real_t, class cplx_t>
inline void CalcMultipoleMomentOfParticleMeshCell(const S32vec &idx, const S32 adr_tc, const ReallocatableArray<Ttc> &tc,
                                                  const ReallocatableArray<Ttp> &tp, const ReallocatableArray<Tep> &ep,
                                                  const ReallocatableArray<Tsp> &sp, const S32 n_leaf_limit, const S32 lev_leaf_limit,
                                                  const S32 adr_tree_sp_first, const F64 r_crit_sq_on_root_cell, const F64vec &pos_unit_cell,
                                                  const S32 icut, const S32 lev_pm_cell, const F64vec &width_pm_cell, const F64vec &center,
                                                  const F64ort &pos_my_domain, MultipoleMoment<real_t, cplx_t> &mm, const S32 p_spj2mm) {
    assert(p_spj2mm >= 0);
    MultipoleMoment<real_t, cplx_t> mm_spj;
    mm_spj.alloc(mm.p, true);  // if the second argument is true, the buffer is initialized to zero.
    const auto &tc_mom = tc[adr_tc].mom_;
    mm_spj.buf[0] = tc_mom.getCharge();
    // Set dipole (l=1)
    if (p_spj2mm >= 1) {
        const F64vec dipole = GetMyDipole(tc_mom);
        mm_spj.buf[1] = 0.5 * dipole.y;
        mm_spj.buf[2] = -dipole.z;
        mm_spj.buf[3] = 0.5 * dipole.x;
    }
    // Set quadrupole (l=2)
    if (p_spj2mm >= 2) {
        const F64mat quad = GetMyQuadrupole(tc_mom);
        mm_spj.buf[4] = 0.25 * quad.xy;
        mm_spj.buf[5] = 0.5 * quad.yz;
        mm_spj.buf[6] = 0.25 * (2.0 * quad.zz - quad.xx - quad.yy);
        mm_spj.buf[7] = -0.5 * quad.xz;
        mm_spj.buf[8] = 0.125 * (quad.xx - quad.yy);
    }
    F64vec center_tc = tc[adr_tc].mom_.pos;  // FOR TEST by M.I.
    F64vec center_pc = center;               // FOR TEST by M.I.
    mm.assign_from_MM(mm_spj, center_pc, center_tc);

    for (int i = 0; i < 9; i++) {
        std::cout << "mm.buf[" << i << "] = " << mm.buf[i] << std::endl;
    }

#if 0
    const Ttc *tc_tmp = tc.getPointer(adr_tc);
    const S32 n_ptcl = tc_tmp->n_ptcl_;
    const U32 adr_ptcl = tc_tmp->adr_ptcl_;
    if (n_ptcl > 0) {  // PM cell having particles
        ReallocatableArray<S32> adr_ep_list;
        ReallocatableArray<S32> adr_sp_list;
        // Make lists of array index of ep and sp
#if !defined(PARTICLE_SIMULATOR_USE_EPJ_ONLY_TO_EVAL_MM_IN_PMM)
        F64 r_crit_sq = r_crit_sq_on_root_cell;
        for (S32 lev = 1; lev <= lev_pm_cell; lev++) r_crit_sq *= 0.25;
        if (r_crit_sq >= 0.0) {
            // This case corresponds to \theta > 0.
            // In this case, we use local SPJs to calculate MM if local SPJs
            // satisfies the opening angle criterion to the PM cells separated
            // by distance specified by icut.

            const F64vec pos = tc_tmp->mom_.getPos();
            const F64ort my_box = GetBoxOfParticleMeshCell(idx, pos_unit_cell, width_pm_cell);
            ReallocatableArray<F64ort> target_box;
            ReallocatableArray<F64vec> shift_vec;
            shift_vec.push_back(F64vec(width_pm_cell.x * (icut + 1), 0, 0));   // +x
            shift_vec.push_back(F64vec(-width_pm_cell.x * (icut + 1), 0, 0));  // -x
            shift_vec.push_back(F64vec(0, width_pm_cell.y * (icut + 1), 0));   // +y
            shift_vec.push_back(F64vec(0, -width_pm_cell.y * (icut + 1), 0));  // -y
            shift_vec.push_back(F64vec(0, 0, width_pm_cell.z * (icut + 1)));   // +z
            shift_vec.push_back(F64vec(0, 0, -width_pm_cell.z * (icut + 1)));  // -z
            for (S32 k = 0; k < shift_vec.size(); k++) target_box.push_back(my_box.shift(shift_vec[k]));
            if (IsOpen(tc, adr_tc, target_box, r_crit_sq)) {
                MakeList<TSM, Ttc, Ttp, Tep, Tsp>(tc, adr_tc, tp, ep, adr_ep_list, sp, adr_sp_list, target_box, r_crit_sq, n_leaf_limit,
                                                  lev_leaf_limit, adr_tree_sp_first);
            } else {
                // In this case, we use SPJ corresponding to this tree cell
                // for MM calculation.
                const S32 adr_sp = adr_tree_sp_first + adr_tc;
                adr_sp_list.push_back(adr_sp);
            }
        } else {
            // This case corresponds to \theta = 0.
            // In this case, we use only EPJs to calculate multipole moments.
            S32 cnt_adr_ptcl = adr_ptcl;
            adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
            adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
            for (S32 ip = 0; ip < n_ptcl; ip++, cnt_adr_ptcl++) {
                if (GetMSB(tp[cnt_adr_ptcl].adr_ptcl_) == 0) {
                    const S32 adr_ep = tp[cnt_adr_ptcl].adr_ptcl_;
                    adr_ep_list.pushBackNoCheck(adr_ep);
                } else {
                    const S32 adr_sp = ClearMSB(tp[cnt_adr_ptcl].adr_ptcl_);
                    adr_sp_list.pushBackNoCheck(adr_sp);
                }
            }
        }
#else
        // In this case, we use only ``local" EPJs to calculate multipole moments.
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for (S32 ip = 0; ip < n_ptcl; ip++, cnt_adr_ptcl++) {
            if (GetMSB(tp[cnt_adr_ptcl].adr_ptcl_) == 0) {
                const S32 adr_ep = tp[cnt_adr_ptcl].adr_ptcl_;
                if (pos_my_domain.contains(ep[adr_ep].getPos())) {
                    adr_ep_list.pushBackNoCheck(adr_ep);
                }
            }
        }
#endif
        /// TEST ///
        F64vec center_tree_cell = tc[adr_tc].mom_.pos;  // FOR TEST by M.I.
        F64vec center_pm_cell = center;                 // FOR TEST by M.I.
        F64vec di_tc{0.0};                              // FOR TEST by M.I.
        F64mat quad_tc{0.0};                            // FOR TEST by M.I.
        F64vec di_pc{0.0};                              // FOR TEST by M.I.
        F64mat quad_pc{0.0};                            // FOR TEST by M.I.
        F64vec di_pc_from_tc{0.0};                      // FOR TEST by M.I.
        F64mat quad_pc_from_tc{0.0};                    // FOR TEST by M.I.
        F64vec center_test = center_tree_cell;          // FOR TEST by M.I.
        /// TEST ///

        // Calculate multipole moments
        F64 msum{0.0};
        for (S32 i = 0; i < adr_ep_list.size(); i++) {
            const S32 adr = adr_ep_list[i];
            const F64 charge = ep[adr].getCharge();
            const F64vec pos = ep[adr].getPos();
            msum += charge;
            mm.assign_particle(center, pos, charge);
            {
                /// TEST ///
                /*
                std::cout << " CHECK " << std::endl;
                for (int i2 = 0; i2 < 9; i2++) {
                    std::cout << "mm.buf[" << i2 << "]=" << mm.buf[i2] << std::endl;
                }
                */
                const F64vec dr_tc = pos - center_tree_cell;
                di_tc.x += charge * dr_tc.x;
                di_tc.y += charge * dr_tc.y;
                di_tc.z += charge * dr_tc.z;
                quad_tc.xx += charge * dr_tc.x * dr_tc.x;
                quad_tc.yy += charge * dr_tc.y * dr_tc.y;
                quad_tc.zz += charge * dr_tc.z * dr_tc.z;
                quad_tc.xy += charge * dr_tc.x * dr_tc.y;
                quad_tc.yz += charge * dr_tc.y * dr_tc.z;
                quad_tc.xz += charge * dr_tc.x * dr_tc.z;
                // std::cout<<"quad_tc=\n"<<quad_tc<<std::endl;

                auto dr_pc = pos - center_pm_cell;
                di_pc.x += charge * dr_pc.x;
                di_pc.y += charge * dr_pc.y;
                di_pc.z += charge * dr_pc.z;
                quad_pc.xx += charge * dr_pc.x * dr_pc.x;
                quad_pc.yy += charge * dr_pc.y * dr_pc.y;
                quad_pc.zz += charge * dr_pc.z * dr_pc.z;
                quad_pc.xy += charge * dr_pc.x * dr_pc.y;
                quad_pc.yz += charge * dr_pc.y * dr_pc.z;
                quad_pc.xz += charge * dr_pc.x * dr_pc.z;
                // std::cout<<"quad_pc=\n"<<quad_pc<<std::endl;
                // double R_2_0 = 0.25 * charge * (3.0 * dr_pc.z * dr_pc.z - (dr_pc * dr_pc));
                // std::cout << "R_2_0(PM cell center)=" << R_2_0 << std::endl;

                // std::cout<<"quad_test_tmp=\n"<<quad_test_tmp<<std::endl;
                // double R_2_0_tmp = 0.25 * (2.0 * quad_test_tmp.zz * quad_test_tmp.zz - quad_test_tmp.xx * quad_test_tmp.xx - quad_test_tmp.yy *
                // quad_test_tmp.yy); std::cout << "R_2_0_tmp=" << R_2_0_tmp << std::endl;
                /// TEST ///
            }
        }
        MultipoleMoment<real_t, cplx_t> mm_spj;
        mm_spj.alloc(mm.p);

        for (S32 i = 0; i < adr_sp_list.size(); i++) {
            const S32 adr = adr_sp_list[i];
            Tsp sp_tmp;
            if (adr < adr_tree_sp_first) {
                sp_tmp = sp[adr];
            } else {
                const S32 adr_tc = adr - adr_tree_sp_first;
                sp_tmp.copyFromMoment(tc[adr_tc].mom_);
            }
            const F64 charge = sp_tmp.getCharge();
            const F64vec pos = sp_tmp.getPos();
            msum += charge;
            mm_spj.clear();
            {
                /// TEST ///
                const F64vec dr_tc = pos - center_tree_cell;
                const F64vec dipole = GetMyDipole(sp_tmp);
                const F64mat quad = GetMyQuadrupole(sp_tmp);
                quad_tc.xx += quad.xx + 2.0 * dr_tc.x * dipole.x + charge * dr_tc.x * dr_tc.x;
                quad_tc.yy += quad.yy + 2.0 * dr_tc.y * dipole.y + charge * dr_tc.y * dr_tc.y;
                quad_tc.zz += quad.zz + 2.0 * dr_tc.z * dipole.z + charge * dr_tc.z * dr_tc.z;
                quad_tc.xy += quad.xy + (dr_tc.x * dipole.y + dr_tc.y * dipole.x) + charge * dr_tc.x * dr_tc.y;
                quad_tc.xz += quad.xz + (dr_tc.x * dipole.z + dr_tc.z * dipole.x) + charge * dr_tc.x * dr_tc.z;
                quad_tc.yz += quad.yz + (dr_tc.y * dipole.z + dr_tc.z * dipole.y) + charge * dr_tc.y * dr_tc.z;
                /// TEST ///
            }
            // m>=1  mm_spj[l*(l+1)+m] = REAL(R_l^|m|)
            // m<=-1 mm_spj[l*(l+1)+m] = IMG(R_l^|m|)
            // e.g. In the case of l=1, m=0 => idx=1*2+0=2; m=1 => idx=1*2+1=3; m=-1 => idx=1*2-1=1
            // e.g. In the case of l=2, m=0 => idx=2*3+0=6; m=1 => idx=2*3+1=7; m=-1 => idx=2*3-1=5; m=2 => idx=2*3+2=8; m=-2 => idx=2*3-2=4
            // so, R_l^m = mm_spj[l*(l+1)+m] + i mm_spj[l*(l+1)-m]
            // Set monopole (l=0)
            mm_spj.buf[0] = charge;
            // Set dipole (l=1)
            if (p_spj2mm >= 1) {
                const F64vec dipole = GetMyDipole(sp_tmp);
                mm_spj.buf[1] = 0.5 * dipole.y;
                mm_spj.buf[2] = -dipole.z;
                mm_spj.buf[3] = 0.5 * dipole.x;
            }
            // Set quadrupole (l=2)
            if (p_spj2mm >= 2) {
                const F64mat quad = GetMyQuadrupole(sp_tmp);
                mm_spj.buf[4] = 0.25 * quad.xy;
                mm_spj.buf[5] = 0.5 * quad.yz;
                mm_spj.buf[6] = 0.25 * (2.0 * quad.zz - quad.xx - quad.yy);
                mm_spj.buf[7] = -0.5 * quad.xz;
                mm_spj.buf[8] = 0.125 * (quad.xx - quad.yy);
            }
#ifdef PARTICLE_SIMULATOR_USE_PMMM_EXPERIMENTAL_FEATURE
            // Set multipole
            constexpr S32 p_spj = GetMyMultipoleOrder<Tsp>();
            if (p_spj != -1) {
                constexpr S32 buflen = MultipoleMoment0<p_spj>::length;
                F64 *buf = new F64[buflen];
                GetMyMultipole(sp_tmp, buflen, buf);
                const S32 mmlen = mm_spj.buf.size();
                for (S32 k = 0; k < std::min(buflen, mmlen); k++) {
                    mm_spj.buf[k] = buf[k];
                }
                delete[] buf;
            }
#endif

            // [Notes]
            // (1) For the data format of MultipoleMoment class,
            //     see the descriptions of Appendix A.6 in Nitadori (2014)
            //     [arXiv:1409.5981v2].
            // (2) In the above, we assume that `dipole` and `quad` are
            //     calculated by the TreeForForce class as
            //     dipole = \sum_{k} q_{k}(r_{k}-r_{g.c.}) and
            //     quad = \sum_{k} q_{k}(r_{k}-r_{g.c.})_{i}(r_{k}-r_{g.c.})_{j},
            //     where
            //     q_{k} is the electric charge of particle k,
            //     r_{k} is the position vector of particle k,
            //     r_{g.c.} is the geometric center of a tree cell,
            //     and i, j take one of x,y,z.
            // M2M transformation
            mm.assign_from_MM(mm_spj, center, pos);
        }
        mm_spj.freeMem();

        {
            /// TEST ///
            //std::cout << "**************" << std::endl;
            //std::cout << "adr_ep_list.size()= " << adr_ep_list.size() << std::endl;
            //std::cout << "adr_sp_list.size()= " << adr_sp_list.size() << std::endl;
            //std::cout << "adr_tc= " << adr_tc << " adr_tree_sp_first= " << adr_tree_sp_first << std::endl;
            //std::cout << "idx= " << idx << std::endl;
            //tc[adr_tc].dump(std::cout);

            auto mono_org = tc[adr_tc].mom_.charge;

            const auto dr_pc_from_tc = center_tree_cell - center_pm_cell;
            di_pc_from_tc.x = di_tc.x + mono_org * dr_pc_from_tc.x;
            di_pc_from_tc.y = di_tc.y + mono_org * dr_pc_from_tc.y;
            di_pc_from_tc.z = di_tc.z + mono_org * dr_pc_from_tc.z;
            quad_pc_from_tc.xx = quad_tc.xx + 2.0 * di_tc.x * dr_pc_from_tc.x + mono_org * dr_pc_from_tc.x * dr_pc_from_tc.x;
            quad_pc_from_tc.yy = quad_tc.yy + 2.0 * di_tc.y * dr_pc_from_tc.y + mono_org * dr_pc_from_tc.y * dr_pc_from_tc.y;
            quad_pc_from_tc.zz = quad_tc.zz + 2.0 * di_tc.z * dr_pc_from_tc.z + mono_org * dr_pc_from_tc.z * dr_pc_from_tc.z;
            quad_pc_from_tc.xy = quad_tc.xy + (di_tc.x * dr_pc_from_tc.y + di_tc.y * dr_pc_from_tc.x) + mono_org * dr_pc_from_tc.x * dr_pc_from_tc.y;
            quad_pc_from_tc.xz = quad_tc.xz + (di_tc.x * dr_pc_from_tc.z + di_tc.z * dr_pc_from_tc.x) + mono_org * dr_pc_from_tc.x * dr_pc_from_tc.z;
            quad_pc_from_tc.yz = quad_tc.yz + (di_tc.y * dr_pc_from_tc.z + di_tc.z * dr_pc_from_tc.y) + mono_org * dr_pc_from_tc.y * dr_pc_from_tc.z;
            //std::cout << "di_pc=" << di_pc << std::endl;
            //std::cout << "di_pc_from_tc= " << di_pc_from_tc << std::endl;
            //std::cout << "tc[adr_tc].mom_.quad=\n" << tc[adr_tc].mom_.quadrupole << std::endl;
            //std::cout << "quad_tc=\n" << quad_tc << std::endl;
            //std::cout << "quad_pc=\n" << quad_pc << std::endl;
            //std::cout << "quad_pc_from_tc=\n" << quad_pc_from_tc << std::endl;

            MultipoleMoment<real_t, cplx_t> mm_pc;
            mm_pc.alloc(mm.p);
            mm_pc.clear();
            mm_pc.buf[0] = mono_org;
            mm_pc.buf[1] = 0.5 * di_pc.y;
            mm_pc.buf[2] = -di_pc.z;
            mm_pc.buf[3] = 0.5 * di_pc.x;
            mm_pc.buf[4] = 0.25 * quad_pc.xy;
            mm_pc.buf[5] = 0.5 * quad_pc.yz;
            mm_pc.buf[6] = 0.25 * (2.0 * quad_pc.zz - quad_pc.xx - quad_pc.yy);
            mm_pc.buf[7] = -0.5 * quad_pc.xz;
            mm_pc.buf[8] = 0.125 * (quad_pc.xx - quad_pc.yy);
            //for (int i = 0; i < 9; i++) std::cout << "mm_pc.buf[" << i << "]=" << mm_pc.buf[i] << std::endl;

            MultipoleMoment<real_t, cplx_t> mm_tc;
            mm_tc.alloc(mm.p);
            mm_tc.clear();
            mm_tc.buf[0] = mono_org;
            mm_tc.buf[1] = 0.5 * di_tc.y;
            mm_tc.buf[2] = -di_tc.z;
            mm_tc.buf[3] = 0.5 * di_tc.x;
            mm_tc.buf[4] = 0.25 * quad_tc.xy;
            mm_tc.buf[5] = 0.5 * quad_tc.yz;
            mm_tc.buf[6] = 0.25 * (2.0 * quad_tc.zz - quad_tc.xx - quad_tc.yy);
            mm_tc.buf[7] = -0.5 * quad_tc.xz;
            mm_tc.buf[8] = 0.125 * (quad_tc.xx - quad_tc.yy);
            MultipoleMoment<real_t, cplx_t> mm_pc_from_tc;
            mm_pc_from_tc.alloc(mm.p, true);  // true means that the memory is initialized by zero
            //std::cout << "mm_pc_from_tc initialized" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << "mm_pc_from_tc.buf[" << i << "]=" << mm_pc_from_tc.buf[i] << std::endl;
            mm_pc_from_tc.assign_from_MM(mm_tc, center_pm_cell, center_tree_cell);
            //std::cout << "mm_pc_from_tc assigned" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << "mm_pc_from_tc.buf[" << i << "]=" << mm_pc_from_tc.buf[i] << std::endl;

            MultipoleMoment<real_t, cplx_t> mm_tc_from_pc_from_tc;
            mm_tc_from_pc_from_tc.alloc(mm.p, true);
            mm_tc_from_pc_from_tc.assign_from_MM(mm_pc_from_tc, center_tree_cell, center_pm_cell);
            //std::cout << "mm_tc, mm_tc_from_pc_from_tc" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << mm_tc.buf[i] << ", " << mm_tc_from_pc_from_tc.buf[i] << std::endl;
            

            //std::cout << "ORIGINAL RESULT" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << "mm.buf[" << i << "]=" << mm.buf[i] << std::endl;
            //std::cout << "CHECK" << std::endl;
            //std::cout << "mm.buf,    mm_pc,buf,   mm_pc_from_tc.buf" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << mm.buf[i] << " " << mm_pc.buf[i] << " " << mm_pc_from_tc.buf[i] << std::endl;
            /// TEST ///
        }
        // Free memory
        adr_ep_list.freeMem();
        adr_sp_list.freeMem();
    }
#endif
}
#else
// ORIGINAL IMPLEMENTATION WITH DEBUG PRINT
// not sure if tree walk is necessary
/*
template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, class real_t, class cplx_t>
inline void CalcMultipoleMomentOfParticleMeshCell(const S32vec &idx, const S32 adr_tc, const ReallocatableArray<Ttc> &tc,
                                                  const ReallocatableArray<Ttp> &tp, const ReallocatableArray<Tep> &ep,
                                                  const ReallocatableArray<Tsp> &sp, const S32 n_leaf_limit, const S32 lev_leaf_limit,
                                                  const S32 adr_tree_sp_first, const F64 r_crit_sq_on_root_cell, const F64vec &pos_unit_cell,
                                                  const S32 icut, const S32 lev_pm_cell, const F64vec &width_pm_cell, const F64vec &center,
                                                  const F64ort &pos_my_domain, MultipoleMoment<real_t, cplx_t> &mm, const S32 p_spj2mm) {
                                                  */
template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, class real_t, class cplx_t, template<class, class> class Tmm>
inline void CalcMultipoleMomentOfParticleMeshCell(const S32vec &idx, const S32 adr_tc, const ReallocatableArray<Ttc> &tc,
                                                  const ReallocatableArray<Ttp> &tp, const ReallocatableArray<Tep> &ep,
                                                  const ReallocatableArray<Tsp> &sp, const S32 n_leaf_limit, const S32 lev_leaf_limit,
                                                  const S32 adr_tree_sp_first, const F64 r_crit_sq_on_root_cell, const F64vec &pos_unit_cell,
                                                  const S32 icut, const S32 lev_pm_cell, const F64vec &width_pm_cell, const F64vec &center,
                                                  const F64ort &pos_my_domain, Tmm<real_t, cplx_t> &mm, const S32 p_spj2mm) {
    assert(p_spj2mm >= 0);
    const Ttc *tc_tmp = tc.getPointer(adr_tc);
    const S32 n_ptcl = tc_tmp->n_ptcl_;
    const U32 adr_ptcl = tc_tmp->adr_ptcl_;
    //std::cout<<"check 1"<<std::endl;
    if (n_ptcl > 0) {  // PM cell having particles
        
        ReallocatableArray<S32> adr_ep_list;
        ReallocatableArray<S32> adr_sp_list;
        // Make lists of array index of ep and sp
#if !defined(PARTICLE_SIMULATOR_USE_EPJ_ONLY_TO_EVAL_MM_IN_PMM)
        //std::cout<<"check 2"<<std::endl;
        F64 r_crit_sq = r_crit_sq_on_root_cell;
        for (S32 lev = 1; lev <= lev_pm_cell; lev++) r_crit_sq *= 0.25;
        if (r_crit_sq >= 0.0) {
            // This case corresponds to \theta > 0.
            // In this case, we use local SPJs to calculate MM if local SPJs
            // satisfies the opening angle criterion to the PM cells separated
            // by distance specified by icut.

            //const F64vec pos = tc_tmp->mom_.getPos();
            const F64ort my_box = GetBoxOfParticleMeshCell(idx, pos_unit_cell, width_pm_cell);
            ReallocatableArray<F64ort> target_box;
            ReallocatableArray<F64vec> shift_vec;
            shift_vec.push_back(F64vec(width_pm_cell.x * (icut + 1), 0, 0));   // +x
            shift_vec.push_back(F64vec(-width_pm_cell.x * (icut + 1), 0, 0));  // -x
            shift_vec.push_back(F64vec(0, width_pm_cell.y * (icut + 1), 0));   // +y
            shift_vec.push_back(F64vec(0, -width_pm_cell.y * (icut + 1), 0));  // -y
            shift_vec.push_back(F64vec(0, 0, width_pm_cell.z * (icut + 1)));   // +z
            shift_vec.push_back(F64vec(0, 0, -width_pm_cell.z * (icut + 1)));  // -z
            for (S32 k = 0; k < shift_vec.size(); k++) target_box.push_back(my_box.shift(shift_vec[k]));
            if (IsOpen(tc, adr_tc, target_box, r_crit_sq)) {
                MakeList<TSM, Ttc, Ttp, Tep, Tsp>(tc, adr_tc, tp, ep, adr_ep_list, sp, adr_sp_list, target_box, r_crit_sq, n_leaf_limit,
                                                  lev_leaf_limit, adr_tree_sp_first);
            } else {
                // In this case, we use SPJ corresponding to this tree cell
                // for MM calculation.
                const S32 adr_sp = adr_tree_sp_first + adr_tc;
                adr_sp_list.push_back(adr_sp);
            }
        } else {
            // This case corresponds to \theta = 0.
            S32 cnt_adr_ptcl = adr_ptcl;
            adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
            adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
            for (S32 ip = 0; ip < n_ptcl; ip++, cnt_adr_ptcl++) {
                if (GetMSB(tp[cnt_adr_ptcl].adr_ptcl_) == 0) {
                    const S32 adr_ep = tp[cnt_adr_ptcl].adr_ptcl_;
                    adr_ep_list.pushBackNoCheck(adr_ep);
                } else {
                    const S32 adr_sp = ClearMSB(tp[cnt_adr_ptcl].adr_ptcl_);
                    adr_sp_list.pushBackNoCheck(adr_sp);
                }
            }
        }
#else
        //std::cout<<"check 3"<<std::endl;
        // In this case, we use only ``local" EPJs to calculate multipole moments.
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for (S32 ip = 0; ip < n_ptcl; ip++, cnt_adr_ptcl++) {
            if (GetMSB(tp[cnt_adr_ptcl].adr_ptcl_) == 0) {
                const S32 adr_ep = tp[cnt_adr_ptcl].adr_ptcl_;
                if (pos_my_domain.contains(ep[adr_ep].getPos())) {
                    adr_ep_list.pushBackNoCheck(adr_ep);
                }
            }
        }
#endif
        //std::cout<<"check 4"<<std::endl;
        /// TEST ///
        F64vec center_tree_cell = tc[adr_tc].mom_.pos;  // FOR TEST by M.I.
        F64vec center_pm_cell = center;                 // FOR TEST by M.I.
        F64vec di_tc{0.0};                              // FOR TEST by M.I.
        F64mat quad_tc{0.0};                            // FOR TEST by M.I.
        F64vec di_pc{0.0};                              // FOR TEST by M.I.
        F64mat quad_pc{0.0};                            // FOR TEST by M.I.
        F64vec di_pc_from_tc{0.0};                      // FOR TEST by M.I.
        F64mat quad_pc_from_tc{0.0};                    // FOR TEST by M.I.
        //F64vec center_test = center_tree_cell;          // FOR TEST by M.I.
        /// TEST ///

        // Calculate multipole moments
        F64 msum{0.0};
        for (S32 i = 0; i < adr_ep_list.size(); i++) {
            const S32 adr = adr_ep_list[i];
            const F64 charge = ep[adr].getCharge();
            const F64vec pos = ep[adr].getPos();
            msum += charge;
            mm.assign_particle(center, pos, charge);
            {
                /// TEST ///
                /*
                std::cout << " CHECK " << std::endl;
                for (int i2 = 0; i2 < 9; i2++) {
                    std::cout << "mm.buf[" << i2 << "]=" << mm.buf[i2] << std::endl;
                }
                */
                const F64vec dr_tc = pos - center_tree_cell;
                di_tc.x += charge * dr_tc.x;
                di_tc.y += charge * dr_tc.y;
                di_tc.z += charge * dr_tc.z;
                quad_tc.xx += charge * dr_tc.x * dr_tc.x;
                quad_tc.yy += charge * dr_tc.y * dr_tc.y;
                quad_tc.zz += charge * dr_tc.z * dr_tc.z;
                quad_tc.xy += charge * dr_tc.x * dr_tc.y;
                quad_tc.yz += charge * dr_tc.y * dr_tc.z;
                quad_tc.xz += charge * dr_tc.x * dr_tc.z;
                // std::cout<<"quad_tc=\n"<<quad_tc<<std::endl;

                auto dr_pc = pos - center_pm_cell;
                di_pc.x += charge * dr_pc.x;
                di_pc.y += charge * dr_pc.y;
                di_pc.z += charge * dr_pc.z;
                quad_pc.xx += charge * dr_pc.x * dr_pc.x;
                quad_pc.yy += charge * dr_pc.y * dr_pc.y;
                quad_pc.zz += charge * dr_pc.z * dr_pc.z;
                quad_pc.xy += charge * dr_pc.x * dr_pc.y;
                quad_pc.yz += charge * dr_pc.y * dr_pc.z;
                quad_pc.xz += charge * dr_pc.x * dr_pc.z;
                // std::cout<<"quad_pc=\n"<<quad_pc<<std::endl;
                // double R_2_0 = 0.25 * charge * (3.0 * dr_pc.z * dr_pc.z - (dr_pc * dr_pc));
                // std::cout << "R_2_0(PM cell center)=" << R_2_0 << std::endl;

                // std::cout<<"quad_test_tmp=\n"<<quad_test_tmp<<std::endl;
                // double R_2_0_tmp = 0.25 * (2.0 * quad_test_tmp.zz * quad_test_tmp.zz - quad_test_tmp.xx * quad_test_tmp.xx - quad_test_tmp.yy *
                // quad_test_tmp.yy); std::cout << "R_2_0_tmp=" << R_2_0_tmp << std::endl;
                /// TEST ///
            }
        }
//std::cout<<"check 5"<<std::endl;        
        //MultipoleMoment<real_t, cplx_t> mm_spj;
        Tmm<real_t, cplx_t> mm_spj;
        mm_spj.alloc(mm.p);

//std::cout<<"check 6"<<std::endl;
        for (S32 i = 0; i < adr_sp_list.size(); i++) {
            const S32 adr = adr_sp_list[i];
            Tsp sp_tmp;
            if (adr < adr_tree_sp_first) {
                sp_tmp = sp[adr];
            } else {
                const S32 adr_tc = adr - adr_tree_sp_first;
                sp_tmp.copyFromMoment(tc[adr_tc].mom_);
            }
            const F64 charge = sp_tmp.getCharge();
            const F64vec pos = sp_tmp.getPos();
            msum += charge;
            mm_spj.clear();
            {
                /// TEST ///
                const F64vec dr_tc = pos - center_tree_cell;
                const F64vec di = GetMyDipole(sp_tmp);
                const F64mat quad = GetMyQuadrupole(sp_tmp);
                di_tc.x += di.x + charge * dr_tc.x;
                di_tc.y += di.y + charge * dr_tc.y;
                di_tc.z += di.z + charge * dr_tc.z;
                quad_tc.xx += quad.xx + 2.0 * dr_tc.x * di.x + charge * dr_tc.x * dr_tc.x;
                quad_tc.yy += quad.yy + 2.0 * dr_tc.y * di.y + charge * dr_tc.y * dr_tc.y;
                quad_tc.zz += quad.zz + 2.0 * dr_tc.z * di.z + charge * dr_tc.z * dr_tc.z;
                quad_tc.xy += quad.xy + (dr_tc.x * di.y + dr_tc.y * di.x) + charge * dr_tc.x * dr_tc.y;
                quad_tc.xz += quad.xz + (dr_tc.x * di.z + dr_tc.z * di.x) + charge * dr_tc.x * dr_tc.z;
                quad_tc.yz += quad.yz + (dr_tc.y * di.z + dr_tc.z * di.y) + charge * dr_tc.y * dr_tc.z;

                const F64vec dr_pc = pos - center_pm_cell;
                di_pc.x += di.x + charge * dr_pc.x;
                di_pc.y += di.y + charge * dr_pc.y;
                di_pc.z += di.z + charge * dr_pc.z;                
                quad_pc.xx += quad.xx + 2.0 * dr_pc.x * di.x + charge * dr_pc.x * dr_pc.x;
                quad_pc.yy += quad.yy + 2.0 * dr_pc.y * di.y + charge * dr_pc.y * dr_pc.y;
                quad_pc.zz += quad.zz + 2.0 * dr_pc.z * di.z + charge * dr_pc.z * dr_pc.z;
                quad_pc.xy += quad.xy + (dr_pc.x * di.y + dr_pc.y * di.x) + charge * dr_pc.x * dr_pc.y;
                quad_pc.xz += quad.xz + (dr_pc.x * di.z + dr_pc.z * di.x) + charge * dr_pc.x * dr_pc.z;
                quad_pc.yz += quad.yz + (dr_pc.y * di.z + dr_pc.z * di.y) + charge * dr_pc.y * dr_pc.z;                
                /// TEST ///
            }
            // m>=1  mm_spj[l*(l+1)+m] = REAL(R_l^|m|)
            // m<=-1 mm_spj[l*(l+1)+m] = IMG(R_l^|m|)
            // e.g. In the case of l=1, m=0 => idx=1*2+0=2; m=1 => idx=1*2+1=3; m=-1 => idx=1*2-1=1
            // e.g. In the case of l=2, m=0 => idx=2*3+0=6; m=1 => idx=2*3+1=7; m=-1 => idx=2*3-1=5; m=2 => idx=2*3+2=8; m=-2 => idx=2*3-2=4
            // so, R_l^m = mm_spj[l*(l+1)+m] + i mm_spj[l*(l+1)-m]
            // Set monopole (l=0)
            mm_spj.buf[0] = charge;
            // Set dipole (l=1)
            if (p_spj2mm >= 1) {
                const F64vec dipole = GetMyDipole(sp_tmp);
                mm_spj.buf[1] = 0.5 * dipole.y;
                mm_spj.buf[2] = -dipole.z;
                mm_spj.buf[3] = 0.5 * dipole.x;
            }
            // Set quadrupole (l=2)
            if (p_spj2mm >= 2) {
                const F64mat quad = GetMyQuadrupole(sp_tmp);
                //std::cout << "sp_tmp.getQuadrupole()= " << sp_tmp.getQuadrupole() << std::endl;
                mm_spj.buf[4] = 0.25 * quad.xy;
                mm_spj.buf[5] = 0.5 * quad.yz;
                mm_spj.buf[6] = 0.25 * (2.0 * quad.zz - quad.xx - quad.yy);
                mm_spj.buf[7] = -0.5 * quad.xz;
                mm_spj.buf[8] = 0.125 * (quad.xx - quad.yy);
            }
#ifdef PARTICLE_SIMULATOR_USE_PMMM_EXPERIMENTAL_FEATURE
            // Set multipole
            constexpr S32 p_spj = GetMyMultipoleOrder<Tsp>();
            if (p_spj != -1) {
                constexpr S32 buflen = MultipoleMoment0<p_spj>::length;
                F64 *buf = new F64[buflen];
                GetMyMultipole(sp_tmp, buflen, buf);
                const S32 mmlen = mm_spj.buf.size();
                for (S32 k = 0; k < std::min(buflen, mmlen); k++) {
                    mm_spj.buf[k] = buf[k];
                }
                delete[] buf;
            }
#endif

            // [Notes]
            // (1) For the data format of MultipoleMoment class,
            //     see the descriptions of Appendix A.6 in Nitadori (2014)
            //     [arXiv:1409.5981v2].
            // (2) In the above, we assume that `dipole` and `quad` are
            //     calculated by the TreeForForce class as
            //     dipole = \sum_{k} q_{k}(r_{k}-r_{g.c.}) and
            //     quad = \sum_{k} q_{k}(r_{k}-r_{g.c.})_{i}(r_{k}-r_{g.c.})_{j},
            //     where
            //     q_{k} is the electric charge of particle k,
            //     r_{k} is the position vector of particle k,
            //     r_{g.c.} is the geometric center of a tree cell,
            //     and i, j take one of x,y,z.
            // M2M transformation
            mm.assign_from_MM(mm_spj, center, pos);
        }
        //std::cout<<"check 7"<<std::endl;
        mm_spj.freeMem();
        //std::cout<<"check 8"<<std::endl;
        #if 0
        {
            /// TEST ///
            //std::cout << "**************" << std::endl;
            //std::cout << "adr_ep_list.size()= " << adr_ep_list.size() << std::endl;
            //std::cout << "adr_sp_list.size()= " << adr_sp_list.size() << std::endl;
            //std::cout << "adr_tc= " << adr_tc << " adr_tree_sp_first= " << adr_tree_sp_first << std::endl;
            //std::cout << "idx= " << idx << std::endl;
            //tc[adr_tc].dump(std::cout);

            auto mono_org = tc[adr_tc].mom_.charge;

            const auto dr_pc_from_tc = center_tree_cell - center_pm_cell;
            di_pc_from_tc.x = di_tc.x + mono_org * dr_pc_from_tc.x;
            di_pc_from_tc.y = di_tc.y + mono_org * dr_pc_from_tc.y;
            di_pc_from_tc.z = di_tc.z + mono_org * dr_pc_from_tc.z;
            quad_pc_from_tc.xx = quad_tc.xx + 2.0 * di_tc.x * dr_pc_from_tc.x + mono_org * dr_pc_from_tc.x * dr_pc_from_tc.x;
            quad_pc_from_tc.yy = quad_tc.yy + 2.0 * di_tc.y * dr_pc_from_tc.y + mono_org * dr_pc_from_tc.y * dr_pc_from_tc.y;
            quad_pc_from_tc.zz = quad_tc.zz + 2.0 * di_tc.z * dr_pc_from_tc.z + mono_org * dr_pc_from_tc.z * dr_pc_from_tc.z;
            quad_pc_from_tc.xy = quad_tc.xy + (di_tc.x * dr_pc_from_tc.y + di_tc.y * dr_pc_from_tc.x) + mono_org * dr_pc_from_tc.x * dr_pc_from_tc.y;
            quad_pc_from_tc.xz = quad_tc.xz + (di_tc.x * dr_pc_from_tc.z + di_tc.z * dr_pc_from_tc.x) + mono_org * dr_pc_from_tc.x * dr_pc_from_tc.z;
            quad_pc_from_tc.yz = quad_tc.yz + (di_tc.y * dr_pc_from_tc.z + di_tc.z * dr_pc_from_tc.y) + mono_org * dr_pc_from_tc.y * dr_pc_from_tc.z;
            //std::cout << "di_pc=" << di_pc << std::endl;
            //std::cout << "di_pc_from_tc= " << di_pc_from_tc << std::endl;
            //std::cout << "tc[adr_tc].mom_.quad=\n" << tc[adr_tc].mom_.quadrupole << std::endl;
            //std::cout << "quad_tc=\n" << quad_tc << std::endl;
            //std::cout << "quad_pc=\n" << quad_pc << std::endl;
            //std::cout << "quad_pc_from_tc=\n" << quad_pc_from_tc << std::endl;

            MultipoleMoment<real_t, cplx_t> mm_pc;
            mm_pc.alloc(mm.p, true);
            mm_pc.buf[0] = mono_org;
            mm_pc.buf[1] = 0.5 * di_pc.y;
            mm_pc.buf[2] = -di_pc.z;
            mm_pc.buf[3] = 0.5 * di_pc.x;
            mm_pc.buf[4] = 0.25 * quad_pc.xy;
            mm_pc.buf[5] = 0.5 * quad_pc.yz;
            mm_pc.buf[6] = 0.25 * (2.0 * quad_pc.zz - quad_pc.xx - quad_pc.yy);
            mm_pc.buf[7] = -0.5 * quad_pc.xz;
            mm_pc.buf[8] = 0.125 * (quad_pc.xx - quad_pc.yy);
            //for (int i = 0; i < 9; i++) std::cout << "mm_pc.buf[" << i << "]=" << mm_pc.buf[i] << std::endl;
            
            MultipoleMoment<real_t, cplx_t> mm_tc;
            mm_tc.alloc(mm.p);
            mm_tc.clear();
            mm_tc.buf[0] = mono_org;
            mm_tc.buf[1] = 0.5 * di_tc.y;
            mm_tc.buf[2] = -di_tc.z;
            mm_tc.buf[3] = 0.5 * di_tc.x;
            mm_tc.buf[4] = 0.25 * quad_tc.xy;
            mm_tc.buf[5] = 0.5 * quad_tc.yz;
            mm_tc.buf[6] = 0.25 * (2.0 * quad_tc.zz - quad_tc.xx - quad_tc.yy);
            mm_tc.buf[7] = -0.5 * quad_tc.xz;
            mm_tc.buf[8] = 0.125 * (quad_tc.xx - quad_tc.yy);
            MultipoleMoment<real_t, cplx_t> mm_pc_from_tc;
            mm_pc_from_tc.alloc(mm.p, true);  // true means that the memory is initialized by zero
            //std::cout << "mm_pc_from_tc initialized" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << "mm_pc_from_tc.buf[" << i << "]=" << mm_pc_from_tc.buf[i] << std::endl;
            
            mm_pc_from_tc.assign_from_MM(mm_tc, center_pm_cell, center_tree_cell);
            //std::cout << "mm_pc_from_tc assigned" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << "mm_pc_from_tc.buf[" << i << "]=" << mm_pc_from_tc.buf[i] << std::endl;

            MultipoleMoment<real_t, cplx_t> mm_tc_from_pc_from_tc;
            mm_tc_from_pc_from_tc.alloc(mm.p, true);
            mm_tc_from_pc_from_tc.assign_from_MM(mm_pc_from_tc, center_tree_cell, center_pm_cell);
            //std::cout << "mm_tc, mm_tc_from_pc_from_tc" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << mm_tc.buf[i] << ", " << mm_tc_from_pc_from_tc.buf[i] << std::endl;
            //std::cout << "ORIGINAL RESULT" << std::endl;
            //for (int i = 0; i < 9; i++) std::cout << "mm.buf[" << i << "]=" << mm.buf[i] << std::endl;
            //std::cout << "CHECK" << std::endl;
            //std::cout << "mm.buf,    mm_pc.buf,   mm_pc_from_tc.buf" << std::endl;
            //for (int i = 0; i < 9; i++)  std::cout << mm.buf[i] << " " << mm_pc.buf[i] << " " << mm_pc_from_tc.buf[i] << std::endl;
            /// TEST ///
        }
        #endif
        // Free memory
//std::cout<<"check 9"<<std::endl;
        adr_ep_list.freeMem();
        //std::cout<<"check 10"<<std::endl;
        adr_sp_list.freeMem();
        //std::cout<<"check 11"<<std::endl;
    }
}
#endif

}  // namespace ParticleSimulator
