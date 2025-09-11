/*
 Derived from PMMM (Apache License 2.0)
 Original Copyright (c) 2012-2014 Keigo Nitadori
 Modifications (c) 2025 FDPS team
 Licensed under the Apache License, Version 2.0 (see LICENSE)
*/
#pragma once
#include "../ps_defs.hpp"

namespace ParticleSimulator {

namespace ParticleMeshMultipole {

template <class real_t, class cplx_t, class Tepi, class Tforce>
class Cell_FMM {
   public:
    S32 idx;
    F64vec center;

    S32 n_epi;
    Tepi *epi_first;
    Tforce *force_first;
    // The information of i particles is managed by pairs of the
    // begining address and the number of particles.
    // epi_first and force_first are the pointers to epi_sorted_[]
    // and force_sorted_[] in TreeForForce class, respectively.
    // [TODO]
    //    Fix the problem that epi_first and force_first do not
    //    point correct memory addresses if TreeForForce class
    //    use the MemoryPool class for epi_sorted_ and force_sorted_
    //    and the memory pool inside the MemoryPool class is
    //    reallocated inbetween the PP part and the PM part.

    bool is_mm_defined;
    MultipoleMoment<real_t, cplx_t> mm;
    LocalExpansion<real_t, cplx_t> le;

    Cell_FMM() { this->clear(); }

    ~Cell_FMM() {
        mm.freeMem();
        le.freeMem();
    }

    void clear() {
        idx = -1;
        center = F64vec(0.0);

        n_epi = 0;
        epi_first = nullptr;
        force_first = nullptr;

        is_mm_defined = false;
        mm.freeMem();
        le.freeMem();
    }

    void init(const S32 p) {
        mm.alloc(p);
        le.alloc(p);
    }

    void setIdx(const S32 idx) { this->idx = idx; }

    void setPos(const F64vec &pos) { center = pos; }

    inline void setIParticleInfoToCell(const S32 n_epi, Tepi *epi_first, Tforce *force_first) {
        this->n_epi = n_epi;
        this->epi_first = epi_first;
        this->force_first = force_first;
    }

    void setMultipoleMoment(MultipoleMoment<real_t, cplx_t> &_mm) {
        assert(_mm.p == mm.p);
        is_mm_defined = true;
        const S32 len = (mm.p + 1) * (mm.p + 1);
        for (S32 i = 0; i < len; i++) {
            mm.buf[i] = _mm.buf[i];
        }
    }

    // new
    void do_L2P_full(const bool clear, const bool correct, const F64 msum, const F64 alpha, LocalExpansion<real_t, cplx_t> &le1,
                     Rlm<real_t, cplx_t> &rlm) {
        static const F64 pi = 4.0 * std::atan(1.0);
        const F64 ainv2 = pi / (alpha * alpha);
        const F64 mc = ((2. / 3.) * pi) * msum;
        for (S32 k = 0; k < n_epi; k++) {
            const F64vec pos = epi_first[k].getPos();
            le1.clear();
            le1.assign_from_LE_opt(le, pos, center, rlm);
            F64vec acc;
            F64 pot;
            le1.read_phi_and_grad(pot, acc.x, acc.y, acc.z);
            if (correct) {
                const F64vec dr = pos - center;
                const F64vec acc_corr = (2.0 * mc) * dr;
                const F64 pot_corr = mc * (dr * dr) - (msum * ainv2);
                acc += acc_corr;
                pot += pot_corr;
            }
            // added 2024/08/27
            // acc *= -1.0; //removed 2025/02/06
            // added 2024/08/27
#if 0
            if (clear) force_first[k].clearPMMM(); // removed 2025/08/06
#endif
            // force_first[k].accumulateForcePMMM(acc, pot);
            force_first[k].addForcePMMM(acc, pot);
        }
    }
    void do_L2P_stripe(const bool clear, const bool correct, const F64 msum, const F64 alpha, const S32 n_thread, const S32 thread_id,
                       LocalExpansion<real_t, cplx_t> &le1, Rlm<real_t, cplx_t> &rlm) {
        static const F64 pi = 4.0 * std::atan(1.0);
        const F64 ainv2 = pi / (alpha * alpha);
        const F64 mc = ((2. / 3.) * pi) * msum;
        const S32 head = (n_epi / n_thread) * thread_id + std::min(n_epi % n_thread, thread_id);
        const S32 tail = (n_epi / n_thread) * (thread_id + 1) + std::min(n_epi % n_thread, (thread_id + 1));
        for (S32 k = head; k < tail; k++) {
            const F64vec pos = epi_first[k].getPos();
            le1.clear();
            le1.assign_from_LE_opt(le, pos, center, rlm);
            F64vec acc;
            F64 pot;
            le1.read_phi_and_grad(pot, acc.x, acc.y, acc.z);
            if (correct) {
                const F64vec dr = pos - center;
                const F64vec acc_corr = (2.0 * mc) * dr;
                const F64 pot_corr = mc * (dr * dr) - (msum * ainv2);
                acc += acc_corr;
                pot += pot_corr;
            }
// added 2024/08/27
// acc *= -1.0; //removed 2025/02/06
// added 2024/08/27
#if 0
            if (clear) force_first[k].clearPMMM(); // removed 2025/08/06
#endif
            // force_first[k].accumulateForcePMMM(acc, pot);
            force_first[k].addForcePMMM(acc, pot);
        }
    }

    // old
    void do_L2P(const bool clear) {
        for (S32 k = 0; k < n_epi; k++) {
            LocalExpansion<real_t, cplx_t> le1;
            le1.alloc(1, true);  // allocate and clear
            le1.assign_from_LE(le, epi_first[k].getPos(), center);
            F64 pot, ax, ay, az;
            le1.read_phi_and_grad(pot, ax, ay, az);
#if 0
            if (clear) force_first[k].clearPMMM(); // removed 2025/08/06
#endif
            // force_first[k].accumulateForcePMMM(F64vec(ax, ay, az), pot);
            force_first[k].addForcePMMM(F64vec(ax, ay, az), pot);
            le1.freeMem();
        }
    }
    void do_L2P_corr(const F64 msum, const F64 alpha) {
        const F64 pi = 4.0 * atan(1.0);
        const F64 ainv2 = pi / (alpha * alpha);
        const F64 mc = ((2. / 3.) * pi) * msum;
        for (S32 k = 0; k < n_epi; k++) {
            const F64vec dr = epi_first[k].getPos() - center;
            const F64vec acc_corr = (2.0 * mc) * dr;
            const F64 pot_corr = mc * (dr * dr) - (msum * ainv2);
            // force_first[k].accumulateForcePMMM(acc_corr, pot_corr);
            force_first[k].addForcePMMM(acc_corr, pot_corr);
        }
    }
};

}  // namespace ParticleMeshMultipole
}  // namespace ParticleSimulator
