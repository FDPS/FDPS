#pragma once
/* C++ headers */
#include <complex>
/* FDPS headers */
#include "../ps_defs.hpp"
#include "particle_mesh_multipole_defs.hpp"
#include "cell.hpp"
#include "M2L_engine.hpp"

namespace ParticleSimulator {
namespace ParticleMeshMultipole {

template <class Tforce, class Tepi>
class ParticleMeshMultipole {
   public:
    using class_type = TagClassTypeParticleMeshMultipole;
    using real_t = double;
    using cplx_t = std::complex<real_t>;
    using Cell_t = Cell_FMM<real_t, cplx_t, Tepi, Tforce>;

    const bool is_dummy;

    //        private:
   public:
    // Idenfifier of an instance of ParticleMeshMultipole class
    // (May be useful when debugging)
    std::string id_;

    // MPI communicators
    CommInfo comm_info_;  // MUST be the same communicator used in TreeForForce/DomainInfo classes

    // Variables representing various statuses of ParticleMeshMultipole class
    bool first_call_by_initialize_;

    // Parameters
    Parameters param_;
    S32 p_spj2mm_;
    U32 fftw_planning_rigor_flag_;
    bool use_mpifft_if_possible_;
    S32 fft_size_crit_;

    // Particle information
    S32 n_loc_tot_;
    F64 msum_;
    F64 quad0_;

    // Cell information
    S32 n_cell_loc_;
    S32 n_cell_tot_;
    F64vec width_cell_;
    std::vector<Cell_t> cell_loc_;
    std::unordered_map<S32, S32> adr_cell_loc_from_cell_index_;

    // M2L calculation
    M2L_Engine m2l_engine_;

    // Timing measurement
    TimeProfile time_profile_;

   public:
    F64 wtime_start_GFC_initialize_;
    F64 wtime_end_GFC_initialize_;
    F64 wtime_start_GFC_calc_gf_r_;
    F64 wtime_end_GFC_calc_gf_r_;
    F64 wtime_start_GFC_redist_gf_r_;
    F64 wtime_end_GFC_redist_gf_r_;
    F64 wtime_start_GFC_calc_gf_k_;
    F64 wtime_end_GFC_calc_gf_k_;
    F64 wtime_start_GFC_redist_gf_k_;
    F64 wtime_end_GFC_redist_gf_k_;

    F64 wtime_start_GFC_initialize__copy_;
    F64 wtime_end_GFC_initialize__copy_;
    F64 wtime_start_GFC_initialize__main_;
    F64 wtime_end_GFC_initialize__main_;

    F64 wtime_start_M2L_initialize_;
    F64 wtime_end_M2L_initialize_;
    F64 wtime_start_M2L_preproc_mm_r_comm_;
    F64 wtime_end_M2L_preproc_mm_r_comm_;
    F64 wtime_start_M2L_redist_mm_r_;
    F64 wtime_end_M2L_redist_mm_r_;
    F64 wtime_start_M2L_postproc_mm_r_comm_;
    F64 wtime_end_M2L_postproc_mm_r_comm_;
    F64 wtime_start_M2L_mm_r_to_mm_k_;
    F64 wtime_end_M2L_mm_r_to_mm_k_;
    F64 wtime_start_M2L_gather_mm_k_trans_;
    F64 wtime_end_M2L_gather_mm_k_trans_;
    F64 wtime_start_M2L_transform_;
    F64 wtime_end_M2L_transform_;
    F64 wtime_start_M2L_scatter_le_k_trans_;
    F64 wtime_end_M2L_scatter_le_k_trans_;
    F64 wtime_start_M2L_le_k_to_le_r_;
    F64 wtime_end_M2L_le_k_to_le_r_;
    F64 wtime_start_M2L_preproc_le_r_comm_;
    F64 wtime_end_M2L_preproc_le_r_comm_;
    F64 wtime_start_M2L_redist_le_r_;
    F64 wtime_end_M2L_redist_le_r_;
    F64 wtime_start_M2L_postproc_le_r_comm_;
    F64 wtime_end_M2L_postproc_le_r_comm_;

    F64 wtime_start_M2L_initialize__copy_;
    F64 wtime_end_M2L_initialize__copy_;
    F64 wtime_start_M2L_initialize__main_;
    F64 wtime_end_M2L_initialize__main_;

    F64 wtime_start_calc_force_all_;
    F64 wtime_end_calc_force_all_;
    F64 wtime_start_write_back_;
    F64 wtime_end_write_back_;
    F64 wtime_start_initialize_;
    F64 wtime_end_initialize_;
    F64 wtime_start_set_cell_;
    F64 wtime_end_set_cell_;
    F64 wtime_start_set_ip_info_to_cell_;
    F64 wtime_end_set_ip_info_to_cell_;
    F64 wtime_start_calc_msum_and_quad0_;
    F64 wtime_end_calc_msum_and_quad0_;
    F64 wtime_start_calc_msum_and_quad0__main_;
    F64 wtime_end_calc_msum_and_quad0__main_;
    F64 wtime_start_calc_msum_and_quad0__allreduce_;
    F64 wtime_end_calc_msum_and_quad0__allreduce_;
    F64 wtime_start_calc_multipole_moment_;
    F64 wtime_end_calc_multipole_moment_;
    F64 wtime_start_collect_multipole_moment_;
    F64 wtime_end_collect_multipole_moment_;
    F64 wtime_start_dipole_correction_;
    F64 wtime_end_dipole_correction_;
    F64 wtime_start_L2P_;
    F64 wtime_end_L2P_;

    //        private:
   public:
    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void setParam(const TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
                  const DomainInfo &dinfo) {
        param_.icut = tree.getICut();
        param_.n_cell = tree.getNCell();
        param_.pos_unit_cell = tree.getPosUnitCell();
        param_.bc = dinfo.getBoundaryCondition();
        // Check consistency
        const F64ort pos_root_domain = dinfo.getPosRootDomain();
        // if (Comm::getRank() == 0) {
        //     std::cout << "param_.pos_unit_cell= " << param_.pos_unit_cell << std::endl;
        //     std::cout << "pos_root_domain= " << pos_root_domain << std::endl;
        // }
        if ((param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) &&
            ((param_.pos_unit_cell.low_ != pos_root_domain.low_) || (param_.pos_unit_cell.high_ != pos_root_domain.high_))) {
            PARTICLE_SIMULATOR_PRINT_ERROR("param_ is not consistent with a given DomainInfo.")
            Abort(-1);
        }
    }

    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void setCell(const TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
                 const S32 n_thread) {
        wtime_start_set_cell_ = GetWtime();

        const S32 nx = param_.n_cell.x;
        const S32 ny = param_.n_cell.y;
        const S32 nz = param_.n_cell.z;
        n_cell_tot_ = nx * ny * nz;
        width_cell_ = GetWidthOfParticleMeshCell(param_.pos_unit_cell, param_.n_cell);

        const std::vector<S32> pm_cell_idx_loc = tree.getParticleMeshCellIndexLocal();
        n_cell_loc_ = pm_cell_idx_loc.size();
        cell_loc_.resize(n_cell_loc_);
        adr_cell_loc_from_cell_index_.clear();

        // if(Comm::getRank()==0) std::cout<<"n_cell_loc_="<<n_cell_loc_<<std::endl;

        for (S32 i = 0; i < n_cell_loc_; i++) {
            const S32 idx = pm_cell_idx_loc[i];
            S32vec idx_3d;
            idx_3d.z = idx / (nx * ny);
            idx_3d.y = (idx - (nx * ny) * idx_3d.z) / nx;
            idx_3d.x = idx - (nx * ny) * idx_3d.z - nx * idx_3d.y;
            const F64vec pos = GetCenterOfParticleMeshCell(param_.pos_unit_cell, width_cell_, idx_3d);
            cell_loc_[i].clear();
            cell_loc_[i].init(param_.p);
            cell_loc_[i].setIdx(idx);
            cell_loc_[i].setPos(pos);
            adr_cell_loc_from_cell_index_[idx] = i;
        }
        wtime_end_set_cell_ = GetWtime();
        time_profile_.PMMM__set_cell += (wtime_end_set_cell_ - wtime_start_set_cell_);
#if defined(DEBUG_PRINT_PMMM)
        if (Comm::getRank() == 0) {
            std::cout << "cell_loc_.size()= " << cell_loc_.size() << std::endl;
            for (auto i = 0; i < cell_loc_.size(); i++) {
                std::cout << "i= " << i << " idx= " << cell_loc_[i].idx << " center= " << cell_loc_[i].center << " n_epi= " << cell_loc_[i].n_epi
                          << std::endl;
            }
        }
#endif
    }

    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void setIParticleInfoToCell(
        TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
        const DomainInfo &dinfo, const S32 n_thread) {
        wtime_start_set_ip_info_to_cell_ = GetWtime();

        n_loc_tot_ = tree.getNumberOfEpiSorted();
        PS_OMP(omp parallel for num_threads(n_thread))
        for (S32 i = 0; i < n_cell_loc_; i++) {
            const S32 idx = cell_loc_[i].idx;
            IParticleInfo<Tepi, Tforce> info;
            if (tree.getIParticleInfoOfParticleMeshCell(idx, info) == 0) {
                cell_loc_[i].setIParticleInfoToCell(info.n_epi, info.epi_first, info.force_first);
                // for debug
                // if (comm_info_.getRank() == 3 && (idx == 288 || idx == 226)) {
                //    std::cout << "########" << std::endl;
                //    std::cout << "i = " << i << " idx = " << idx << std::endl;
                //    const S32 nx = param_.n_cell.x;
                //    const S32 ny = param_.n_cell.y;
                //    const S32 nz = param_.n_cell.z;
                //    const F64ort pos_unit_cell = param_.pos_unit_cell;
                //    S32vec idx_3d;
                //    idx_3d.z = idx / (nx * ny);
                //    idx_3d.y = (idx - (nx * ny) * idx_3d.z) / nx;
                //    idx_3d.x = idx - (nx * ny) * idx_3d.z - nx * idx_3d.y;
                //    F64ort box = GetBoxOfParticleMeshCell(idx_3d, pos_unit_cell.low_, width_cell_);
                //    std::cout << "(ix,iy,iz) = " << idx_3d.x << "   " << idx_3d.y << "   " << idx_3d.z << std::endl;
                //    std::cout << "box = " << box << std::endl;
                //    std::cout << "--------" << std::endl;
                //    for (S32 k=0; k < info.n_epi; k++) {
                //        const S64 id = info.epi_first[k].getId();
                //        const F64vec pos = info.epi_first[k].getPos();
                //        std::cout << "k = " << k
                //                  << " id = " << id
                //                  << " pos = " << pos
                //                  << " T/F = " << box.contains(pos)
                //                  << std::endl;
                //    }
                //}
            }
#if defined(DEBUG_PRINT_PMMM)
            if (Comm::getRank() == 0) {
                std::cout << "i= " << i << " cell_loc_[i].idx= " << cell_loc_[i].idx << " center= " << cell_loc_[i].center
                          << " n_epi= " << cell_loc_[i].n_epi << std::endl;
            }
#endif
        }

        wtime_end_set_ip_info_to_cell_ = GetWtime();
        time_profile_.PMMM__set_ip_info_to_cell += (wtime_end_set_ip_info_to_cell_ - wtime_start_set_ip_info_to_cell_);
    }

    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void calcTotalChargeAndDispersion(
        TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
        const DomainInfo &dinfo, const S32 n_thread) {
        F64 wtime_tmp = GetWtime();
        wtime_start_calc_msum_and_quad0_ = wtime_tmp;
        wtime_start_calc_msum_and_quad0__main_ = wtime_tmp;
        // const S32 nx = param_.n_cell.x;
        // const S32 ny = param_.n_cell.y;
        // const S32 nz = param_.n_cell.z;
        const F64ort pos_my_domain = dinfo.getPosDomain(comm_info_.getRank());
        F64 msum_loc{0.0};
        F64 quad0_loc{0.0};
        constexpr S32 branching_coef = 4;
        if (n_cell_loc_ >= branching_coef * n_thread) {
            PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic,1) reduction (+:msum_loc), reduction(+:quad0_loc))
            for (S32 i = 0; i < n_cell_loc_; i++) {
                const F64vec center = cell_loc_[i].center;
                const S32 idx = cell_loc_[i].idx;
                F64 msum, quad0;
                tree.calcTotalChargeAndDispersionOfParticleMeshCellFull(idx, center, pos_my_domain, msum, quad0);
                msum_loc += msum;
                quad0_loc += quad0;
            }
        } else {
            PS_OMP(omp parallel num_threads(n_thread) reduction(+ : msum_loc), reduction(+ : quad0_loc)) {
                const S32 ith = Comm::getThreadNum();
                for (S32 i = 0; i < n_cell_loc_; i++) {
                    const F64vec center = cell_loc_[i].center;
                    const S32 idx = cell_loc_[i].idx;
                    F64 msum, quad0;
                    tree.calcTotalChargeAndDispersionOfParticleMeshCellStripe(idx, center, pos_my_domain, msum, quad0, n_thread, ith);
                    msum_loc += msum;
                    quad0_loc += quad0;
                }
            }
        }
        wtime_tmp = GetWtime();
        wtime_end_calc_msum_and_quad0__main_ = wtime_tmp;
        wtime_start_calc_msum_and_quad0__allreduce_ = wtime_tmp;
        F64vec2 sendbuf, recvbuf;
        sendbuf.x = msum_loc;
        sendbuf.y = quad0_loc;
        recvbuf = comm_info_.getSum(sendbuf);
        msum_ = recvbuf.x;
        quad0_ = recvbuf.y;
        wtime_tmp = GetWtime();
        wtime_end_calc_msum_and_quad0__allreduce_ = wtime_tmp;
        wtime_end_calc_msum_and_quad0_ = wtime_tmp;
        time_profile_.PMMM__calc_msum_and_quad0 += (wtime_end_calc_msum_and_quad0_ - wtime_start_calc_msum_and_quad0_);
#if 0
                if (comm_info_.getRank() == 0) {
                    std::cout << "msum_ = " << msum_ << std::endl;
                    std::cout << "quad0_ = " << quad0_ << std::endl;
                }
#endif
    }

    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void calcMultipoleMoment(
        TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
        const DomainInfo &dinfo, const S32 n_thread) {
        wtime_start_calc_multipole_moment_ = GetWtime();
        // const S32 n_proc = comm_info_.getNumberOfProc();
        // const S32 my_rank = comm_info_.getRank();
        //  MultipoleMoment<real_t, cplx_t> *mm = new MultipoleMoment<real_t, cplx_t>[n_thread];
        MultipoleMoment<real_t, cplx_t> mm[n_thread];
        for (S32 ith = 0; ith < n_thread; ith++) mm[ith].alloc(param_.p);
        // const S32 nx = param_.n_cell.x;
        // const S32 ny = param_.n_cell.y;
        // const S32 nz = param_.n_cell.z;
        const F64ort pos_my_domain = dinfo.getPosDomain(comm_info_.getRank());
        // Calculate multipole moments
        PS_OMP(omp parallel for num_threads(n_thread))
        for (S32 i = 0; i < n_cell_loc_; i++) {
            const S32 ith = comm_info_.getThreadNum();
            const S32 idx = cell_loc_[i].idx;
            const F64vec center = cell_loc_[i].center;
            mm[ith].clear();
            if (tree.template calcMultipoleMomentOfParticleMeshCell<real_t, cplx_t>(idx, center, pos_my_domain, mm[ith], p_spj2mm_) == 0) {
                cell_loc_[i].setMultipoleMoment(mm[ith]);
            }
        }
        for (S32 ith = 0; ith < n_thread; ith++) mm[ith].freeMem();
        // delete[] mm;
        wtime_end_calc_multipole_moment_ = GetWtime();
        time_profile_.PMMM__calc_multipole_moment += (wtime_end_calc_multipole_moment_ - wtime_start_calc_multipole_moment_);
    }

    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void collectMultipoleMoment(
        TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
        const DomainInfo &dinfo, const S32 n_thread) {
        static S64 n_called{0};
        n_called++;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        wtime_start_collect_multipole_moment_ = GetWtime();

        // Calculate convenient variables
        using pair_t = std::pair<S32, S32>;
        // const S32 n_proc = comm_info_.getNumberOfProc();
        // const S32 my_rank = comm_info_.getRank();
        const S32 LEN = (param_.p + 1) * (param_.p + 1);
        // const S32 nx = param_.n_cell.x;
        // const S32 ny = param_.n_cell.y;
        // const S32 nz = param_.n_cell.z;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_COLLECT_MULTIPOLE_MOMENT
        comm_info_.barrier();
        if (Comm::getRank() == 0) std::cout << "ok1 @ collectMultipoleMoment() [PMMM]" << std::endl;
#endif

        // Collect multipole moments necessary to calculate self-energy correction
        std::unordered_map<S32, std::vector<S32>> rank_list = tree.getRankListFromParticleMeshCellIndexLocal();
        std::vector<pair_t> rank_idx_pairs_send;  // send rank & send cell index
        std::vector<pair_t> rank_idx_pairs_recv;  // recv rank & recv cell index
        for (S32 i = 0; i < n_cell_loc_; i++) {
            const S32 idx = cell_loc_[i].idx;
            if (cell_loc_[i].is_mm_defined) {
                const S32 cnt = rank_list[idx].size();
                if (cnt > 1) {
                    for (S32 k = 1; k < cnt; k++) {  // skip k=0 because it is myself.
                        const S32 rank = rank_list[idx][k];
                        pair_t tmp;
                        tmp.first = rank;
                        tmp.second = idx;
                        rank_idx_pairs_send.push_back(tmp);
                    }
                }
            } else {
                const S32 rank = rank_list[idx].front();
                pair_t tmp;
                tmp.first = rank;
                tmp.second = idx;
                rank_idx_pairs_recv.push_back(tmp);
            }
        }
        // Prepare sendbuf
        union buf_t {
            real_t r;  // to store mm
            S32 i;     // to store idx
        };
        CommBuffer<buf_t> sendbuf;
        std::unordered_map<S32, S32> cnt_mp;
        for (size_t i = 0; i < rank_idx_pairs_send.size(); i++) {
            const S32 rank = rank_idx_pairs_send[i].first;
            cnt_mp[rank]++;
        }
        sendbuf.n_comm = cnt_mp.size();
        sendbuf.allocCommInfo();
        S32 ii{0};
        sendbuf.count_tot = 0;
        for (auto itr = cnt_mp.begin(); itr != cnt_mp.end(); ++itr) {
            const S32 rank = itr->first;  // key: rank
            const S32 nc = itr->second;   // val: # of cells
            sendbuf.ranks[ii] = rank;
            sendbuf.counts[ii] = nc * (1 + LEN);
            sendbuf.count_tot += nc * (1 + LEN);
            ii++;
        }
        sendbuf.calcDispls();
        if (sendbuf.count_tot > 0) {
            sendbuf.allocBuffer();
            sendbuf.clearCounts();
            for (size_t i = 0; i < rank_idx_pairs_send.size(); i++) {
                const S32 rank = rank_idx_pairs_send[i].first;
                const S32 idx = rank_idx_pairs_send[i].second;
                S32 adr{-1};
                for (S32 k = 0; k < sendbuf.n_comm; k++) {
                    if (rank == sendbuf.ranks[k]) {
                        adr = k;
                        break;
                    }
                }
                assert(adr != -1);
                // set cell index
                S32 adr_buf = sendbuf.displs[adr] + sendbuf.counts[adr];
                sendbuf.buf[adr_buf++].i = idx;
                sendbuf.counts[adr]++;
                // set mm
                const S32 adr_cell = adr_cell_loc_from_cell_index_[idx];
                for (S32 lm = 0; lm < LEN; lm++) {
                    sendbuf.buf[adr_buf++].r = cell_loc_[adr_cell].mm.buf[lm];
                }
                sendbuf.counts[adr] += LEN;
            }
        }
        // Prepare recvbuf
        CommBuffer<buf_t> recvbuf;
        cnt_mp.clear();
        for (size_t i = 0; i < rank_idx_pairs_recv.size(); i++) {
            const S32 rank = rank_idx_pairs_recv[i].first;
            cnt_mp[rank]++;
        }
        recvbuf.n_comm = cnt_mp.size();
        recvbuf.allocCommInfo();
        ii = 0;
        recvbuf.count_tot = 0;
        for (auto itr = cnt_mp.begin(); itr != cnt_mp.end(); ++itr) {
            const S32 rank = itr->first;  // key: rank
            const S32 nc = itr->second;   // val: # of cells
            recvbuf.ranks[ii] = rank;
            recvbuf.counts[ii] = nc * (1 + LEN);
            recvbuf.count_tot += nc * (1 + LEN);
            ii++;
        }
        recvbuf.calcDispls();
        recvbuf.allocBuffer();
#if 0
                // Check
                //if (n_called == 2) {
                {
                    std::stringstream ss;
                    ss << "mm_comm_table_send_" << std::setfill('0') << std::setw(5) << my_rank << ".txt";
                    const std::string file_name = ss.str();
                    std::ofstream ofs;
                    ofs.open(file_name.c_str(), std::ios::trunc);
                    for (S32 i = 0; i < sendbuf.n_comm; i++) {
                        const S32 rank = sendbuf.ranks[i];
                        const S32 cnt = sendbuf.counts[i];
                        if (cnt > 0) {
                            ofs << my_rank << " ---> " << rank << " : " << cnt << std::endl;
                        }
                    }
                    ofs.close();
                }
                {
                    std::stringstream ss;
                    ss << "mm_comm_table_recv_" << std::setfill('0') << std::setw(5) << my_rank << ".txt";
                    const std::string file_name = ss.str();
                    std::ofstream ofs;
                    ofs.open(file_name.c_str(), std::ios::trunc);
                    for (S32 i = 0; i < recvbuf.n_comm; i++) {
                        const S32 rank = recvbuf.ranks[i];
                        const S32 cnt = recvbuf.counts[i];
                        if (cnt > 0) {
                            ofs << rank << " ---> " << my_rank << " : " << cnt << std::endl;
                        }
                    }
                    ofs.close();
                }
                Finalize(); std::exit(0);
#if 0
                {
                    std::stringstream ss;
                    ss << "mm_comm_info_" << std::setfill('0') << std::setw(5) << my_rank << ".txt";
                    const std::string file_name = ss.str();
                    std::ofstream ofs;
                    ofs.open(file_name.c_str(), std::ios::trunc);
                    ofs << "< sendbuf >" << std::endl;
                    idx_send.dumpCommInfo(ofs);
                    ofs << std::endl;
                    ofs << "< recvbuf >" << std::endl;
                    idx_recv.dumpCommInfo(ofs);
                    ofs << std::endl;
                    ofs.close();
                    Finalize();
                    std::exit(0);
                }
#endif
#endif
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_COLLECT_MULTIPOLE_MOMENT
        comm_info_.barrier();
        if (Comm::getRank() == 0) std::cout << "ok2 @ collectMultipoleMoment() [PMMM]" << std::endl;
#endif
        // Perform MPI comm.
        performComm(sendbuf, recvbuf, comm_info_.getCommunicator());

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_COLLECT_MULTIPOLE_MOMENT
        comm_info_.barrier();
        if (Comm::getRank() == 0) std::cout << "ok3 @ collectMultipoleMoment() [PMMM]" << std::endl;
#endif

        // Copy the contents of recieve buffers to cell_loc_[].mm.buf[]
        if (recvbuf.n_comm > 0) {
            S32 adr_buf{0};
            while (adr_buf < recvbuf.count_tot) {
                const S32 idx = recvbuf.buf[adr_buf++].i;
                const S32 adr_cell = adr_cell_loc_from_cell_index_[idx];
                for (S32 lm = 0; lm < LEN; lm++) {
                    cell_loc_[adr_cell].mm.buf[lm] = recvbuf.buf[adr_buf++].r;
                }
            }
            assert(adr_buf == recvbuf.count_tot);
        }
        wtime_end_collect_multipole_moment_ = GetWtime();
        time_profile_.PMMM__collect_multipole_moment += (wtime_end_collect_multipole_moment_ - wtime_start_collect_multipole_moment_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_COLLECT_MULTIPOLE_MOMENT
        comm_info_.barrier();
        if (Comm::getRank() == 0) std::cout << "ok_last @ collectMultipoleMoment() [PMMM]" << std::endl;
#endif
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
        const F64 now = GetWtime();
        wtime_start_collect_multipole_moment_ = now;
        wtime_end_collect_multipole_moment_ = now;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
    }

    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void calcForceParticleMesh(
        TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
        const DomainInfo &dinfo, const bool clear_force, const S32 n_thread) {
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        constexpr S32 rank_wld_debug_print = 0;
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok1 @ calcForceParticleMesh() [PMMM]" << std::endl;
#endif
        // Get data describing the set of PM cells that covers the domain for each rank
        std::vector<S32ort> pos_pm_domain = tree.getPosParticleMeshDomain();

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok2 @ calcForceParticleMesh() [PMMM]" << std::endl;
#endif

        // Determine the way of parallization
        m2l_engine_.initialize(param_, fftw_planning_rigor_flag_, use_mpifft_if_possible_, fft_size_crit_, comm_info_, n_thread);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok3 @ calcForceParticleMesh() [PMMM]" << std::endl;
#endif

        // Set green function
        m2l_engine_.setGreenFunction();

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok4 @ calcForceParticleMesh() [PMMM]" << std::endl;
#endif

        // Collect information of MM necessary to perform local task
        m2l_engine_.redistMM(n_cell_loc_, cell_loc_);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok5 @ calcForceParticleMesh() [PMMM]" << std::endl;
#endif

        // Perform M2L transform
        m2l_engine_.convolution();

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok6 @ calcForceParticleMesh() [PMMM]" << std::endl;
#endif

        // Collect information of LE necessary to evaluate the potential and
        // its gradient at the positions of local particles.
        m2l_engine_.redistLE(pos_pm_domain, n_cell_loc_, adr_cell_loc_from_cell_index_, cell_loc_);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok7 @ calcForceParticleMesh() [PMMM]" << std::endl;
#endif



        // Calculate total charge and dispersion

        // Collect multipole moment necessary to calculate the dipole correction term

        // Dipole correction
        // [TODO]
        //   - Move the calculation of msum_ & quad0_ here.
        //   - Move the call of collectMultipoleMoment() here.
        wtime_start_dipole_correction_ = GetWtime();
        if (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
            F64 dipole_x, dipole_y, dipole_z;
            dipole_x = dipole_y = dipole_z = 0.0;
            // auto n_cell_loc_dgb = 0;
            PS_OMP(omp parallel for num_threads(n_thread) reduction(+:dipole_x), reduction(+:dipole_y), reduction(+:dipole_z))
            for (S32 i = 0; i < n_cell_loc_; i++) {
                if (cell_loc_[i].is_mm_defined) {
                    dipole_x += cell_loc_[i].mm.buf[3];
                    dipole_y += cell_loc_[i].mm.buf[1];
                    dipole_z += cell_loc_[i].mm.buf[2];
                    // n_cell_loc_dgb++;
                }
            }
            // auto n_cell_glb_dgb = comm_info_.getSum(n_cell_loc_dgb);
            // std::cout<<"n_cell_glb_dgb= "<<n_cell_glb_dgb<<std::endl;
            const F64vec dipole_loc(dipole_x, dipole_y, dipole_z);
            F64vec dipole = comm_info_.getSum(dipole_loc);
            const F64 pi = 4.0 * atan(1.0);
            dipole *= (4. / 3.) * pi;
#if 0
                    DebugUtils::fout() << "msum_  = " << msum_ << std::endl;
                    DebugUtils::fout() << "quad0_ = " << quad0_ << std::endl;
                    DebugUtils::fout() << "dipole = " << dipole << std::endl;
#endif
            /*
                        if (PS::Comm::getRank() == 0) {
                            std::cout << "msum_  = " << msum_ << std::endl;
                            std::cout << "quad0_ = " << quad0_ << std::endl;
                            std::cout << "dipole = " << dipole << std::endl;
                            for (S32 i = 0; i < n_cell_loc_; i++) {
                                if (cell_loc_[i].is_mm_defined) {
                                    std::cout<<"cell_loc_[i].idx= "<<cell_loc_[i].idx<<" n_epi= "<<cell_loc_[i].n_epi<<" center=
               "<<cell_loc_[i].center<<std::endl; std::cout<<"cell_loc_[i].mm= "<<cell_loc_[i].mm.buf[0]<<" "<<cell_loc_[i].mm.buf[1]<<"
               "<<cell_loc_[i].mm.buf[2]<<" "<<cell_loc_[i].mm.buf[3]
                                    <<std::endl;
                                    if(cell_loc_[i].n_epi > 0){
                                        std::cout<<" "<<cell_loc_[i].epi_first[0].id<<" "<<cell_loc_[i].epi_first[0].pos<<std::endl;
                                    }

                                }
                            }
                        }
            */
            PS_OMP(omp parallel for num_threads(n_thread))
            for (S32 i = 0; i < n_cell_loc_; i++) {
                // force correction (B.9)
                cell_loc_[i].le.buf[3] += 2.0 * dipole.x;
                cell_loc_[i].le.buf[1] -= 2.0 * dipole.y;
                cell_loc_[i].le.buf[2] += 1.0 * dipole.z;

                /*
                                        // added 2024/08/24 wrong
                                        const auto qi = tree.getEPISorted(i).getCharge();
                                        std::cout<<"qi= "<<qi<<std::endl;
                                        // const auto mi = tree.getEPISorted(i).mass;
                                        if (qi < 0.0)
                                        {
                                            cell_loc_[i].le.buf[3] *= -1.0;
                                            cell_loc_[i].le.buf[1] *= -1.0;
                                            cell_loc_[i].le.buf[2] *= -1.0;
                                        }
                */
                // std::cerr << "A) cell_loc_[i].le.buf[0]= " << cell_loc_[i].le.buf[0] << std::endl;
                //  potential correction (B.10)
                cell_loc_[i].le.buf[0] += ((2. / 3.) * pi) * quad0_;
                // std::cerr << "B) cell_loc_[i].le.buf[0]= " << cell_loc_[i].le.buf[0] << std::endl;
                //  self energy correction (B.10)
                cell_loc_[i].le.buf[0] -= param_.alpha * (2.0 / sqrt(pi)) * cell_loc_[i].mm.buf[0];
                // std::cerr << "C) cell_loc_[i].le.buf[0]= " << cell_loc_[i].le.buf[0] << std::endl;
            }
        }
        wtime_end_dipole_correction_ = GetWtime();
        time_profile_.PMMM__dipole_correction += (wtime_end_dipole_correction_ - wtime_start_dipole_correction_);

//#if 0  // 20250815

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok8 @ calcForceParticleMesh() [PMMM]" << std::endl;

        // Check le_r after dipole correction
#if 0
                static S32 n_called {0};
                n_called++;
                if (n_called == 2) {
                    std::stringstream ss;
                    ss << "le_r_" << std::setfill('0') << std::setw(5)
                       << Comm::getRank() << ".txt";
                    const std::string filename = ss.str();
                    std::ofstream output_file;
                    output_file.open(filename.c_str(), std::ios::trunc);
                    for (S32 i = 0; i < n_cell_loc_; i++) {
                        const S32 idx_1d = cell_loc_[i].idx;
                        const S32 nx = param_.n_cell.x;
                        const S32 ny = param_.n_cell.y;
                        const S32 nz = param_.n_cell.z;
                        S32vec idx_3d;
                        idx_3d.z = idx_1d / (nx * ny);
                        idx_3d.y = (idx_1d - (nx * ny) * idx_3d.z) / nx;
                        idx_3d.x = idx_1d - (nx * ny) * idx_3d.z - nx * idx_3d.y;
                        const S32 LEN = (param_.p + 1) * (param_.p + 1);
                        for (S32 lm = 0; lm < LEN; lm++) {
                            const S32 idx = lm
                                          + LEN * (idx_3d.x
                                          + nx * (idx_3d.y
                                          + ny * idx_3d.z));
                            const real_t val = cell_loc_[i].le.buf[lm];
                            output_file << idx << "    " << val << std::endl;
                        }
                    }
                    output_file.close();
                    Finalize();
                    std::exit(0);
                }
#endif
#endif
        // To check
        // const S32 n_cell_loc_max = comm_info_.getMaxValue(n_cell_loc_);
        // const S32 n_cell_loc_min = comm_info_.getMinValue(n_cell_loc_);
        // DebugUtils::fout() << "n_cell_loc_ = " << n_cell_loc_ << std::endl;
        // DebugUtils::fout() << "n_cell_loc_max = " << n_cell_loc_max << std::endl;
        // DebugUtils::fout() << "n_cell_loc_min = " << n_cell_loc_min << std::endl;

        wtime_start_L2P_ = GetWtime();
#if 1
        //-----
        // New
        //-----
        const bool correct = (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ);
        std::vector<LocalExpansion<real_t, cplx_t>> le1_tmp(n_thread);
        std::vector<Rlm<real_t, cplx_t>> rlm_tmp(n_thread);
        for (S32 i = 0; i < n_thread; i++) {
            le1_tmp[i].alloc(1, true);
            rlm_tmp[i].alloc(param_.p);
        }
        constexpr S32 branching_coef = 4;
        if (n_cell_loc_ >= branching_coef * n_thread) {
            PS_OMP(omp parallel for num_threads(n_thread) schedule(dynamic, 1))
            for (S32 i = 0; i < n_cell_loc_; i++) {
                const S32 ith = Comm::getThreadNum();
                cell_loc_[i].do_L2P_full(clear_force, correct, msum_, param_.alpha, le1_tmp[ith], rlm_tmp[ith]);
            }
        } else {
            PS_OMP(omp parallel num_threads(n_thread)) {
                const S32 ith = Comm::getThreadNum();
                for (S32 i = 0; i < n_cell_loc_; i++) {
                    cell_loc_[i].do_L2P_stripe(clear_force, correct, msum_, param_.alpha, n_thread, ith, le1_tmp[ith], rlm_tmp[ith]);
                }
            }
        }
        for (S32 i = 0; i < n_thread; i++) {
            le1_tmp[i].freeMem();
            rlm_tmp[i].freeMem();
        }
#else
        //-------
        // Old
        //-------
                PS_OMP(omp parallel for schedule(dynamic, 1))
                for (S32 i = 0; i < n_cell_loc_; i++) {
                    cell_loc_[i].do_L2P(clear_force);
                    if (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                        cell_loc_[i].do_L2P_corr(msum_, param_.alpha);
                    }
                }
#endif
        wtime_end_L2P_ = GetWtime();
        time_profile_.PMMM__L2P += (wtime_end_L2P_ - wtime_start_L2P_);
//#endif 0  // 20250815

        // Copy wall time information
        wtime_start_GFC_initialize_ = m2l_engine_.wtime_start_GFC_initialize_;
        wtime_end_GFC_initialize_ = m2l_engine_.wtime_end_GFC_initialize_;
        wtime_start_GFC_calc_gf_r_ = m2l_engine_.wtime_start_GFC_calc_gf_r_;
        wtime_end_GFC_calc_gf_r_ = m2l_engine_.wtime_end_GFC_calc_gf_r_;
        wtime_start_GFC_redist_gf_r_ = m2l_engine_.wtime_start_GFC_redist_gf_r_;
        wtime_end_GFC_redist_gf_r_ = m2l_engine_.wtime_end_GFC_redist_gf_r_;
        wtime_start_GFC_calc_gf_k_ = m2l_engine_.wtime_start_GFC_calc_gf_k_;
        wtime_end_GFC_calc_gf_k_ = m2l_engine_.wtime_end_GFC_calc_gf_k_;
        wtime_start_GFC_redist_gf_k_ = m2l_engine_.wtime_start_GFC_redist_gf_k_;
        wtime_end_GFC_redist_gf_k_ = m2l_engine_.wtime_end_GFC_redist_gf_k_;

        wtime_start_GFC_initialize__copy_ = m2l_engine_.wtime_start_GFC_initialize__copy_;
        wtime_end_GFC_initialize__copy_ = m2l_engine_.wtime_end_GFC_initialize__copy_;
        wtime_start_GFC_initialize__main_ = m2l_engine_.wtime_start_GFC_initialize__main_;
        wtime_end_GFC_initialize__main_ = m2l_engine_.wtime_end_GFC_initialize__main_;

        wtime_start_M2L_initialize_ = m2l_engine_.wtime_start_initialize_;
        wtime_end_M2L_initialize_ = m2l_engine_.wtime_end_initialize_;
        wtime_start_M2L_preproc_mm_r_comm_ = m2l_engine_.wtime_start_preproc_mm_r_comm_;
        wtime_end_M2L_preproc_mm_r_comm_ = m2l_engine_.wtime_end_preproc_mm_r_comm_;
        wtime_start_M2L_redist_mm_r_ = m2l_engine_.wtime_start_redist_mm_r_;
        wtime_end_M2L_redist_mm_r_ = m2l_engine_.wtime_end_redist_mm_r_;
        wtime_start_M2L_postproc_mm_r_comm_ = m2l_engine_.wtime_start_postproc_mm_r_comm_;
        wtime_end_M2L_postproc_mm_r_comm_ = m2l_engine_.wtime_end_postproc_mm_r_comm_;
        wtime_start_M2L_mm_r_to_mm_k_ = m2l_engine_.wtime_start_mm_r_to_mm_k_;
        wtime_end_M2L_mm_r_to_mm_k_ = m2l_engine_.wtime_end_mm_r_to_mm_k_;
        wtime_start_M2L_gather_mm_k_trans_ = m2l_engine_.wtime_start_gather_mm_k_trans_;
        wtime_end_M2L_gather_mm_k_trans_ = m2l_engine_.wtime_end_gather_mm_k_trans_;
        wtime_start_M2L_transform_ = m2l_engine_.wtime_start_transform_;
        wtime_end_M2L_transform_ = m2l_engine_.wtime_end_transform_;
        wtime_start_M2L_scatter_le_k_trans_ = m2l_engine_.wtime_start_scatter_le_k_trans_;
        wtime_end_M2L_scatter_le_k_trans_ = m2l_engine_.wtime_end_scatter_le_k_trans_;
        wtime_start_M2L_le_k_to_le_r_ = m2l_engine_.wtime_start_le_k_to_le_r_;
        wtime_end_M2L_le_k_to_le_r_ = m2l_engine_.wtime_end_le_k_to_le_r_;
        wtime_start_M2L_preproc_le_r_comm_ = m2l_engine_.wtime_start_preproc_le_r_comm_;
        wtime_end_M2L_preproc_le_r_comm_ = m2l_engine_.wtime_end_preproc_le_r_comm_;
        wtime_start_M2L_redist_le_r_ = m2l_engine_.wtime_start_redist_le_r_;
        wtime_end_M2L_redist_le_r_ = m2l_engine_.wtime_end_redist_le_r_;
        wtime_start_M2L_postproc_le_r_comm_ = m2l_engine_.wtime_start_postproc_le_r_comm_;
        wtime_end_M2L_postproc_le_r_comm_ = m2l_engine_.wtime_end_postproc_le_r_comm_;

        wtime_start_M2L_initialize__copy_ = m2l_engine_.wtime_start_initialize__copy_;
        wtime_end_M2L_initialize__copy_ = m2l_engine_.wtime_end_initialize__copy_;
        wtime_start_M2L_initialize__main_ = m2l_engine_.wtime_start_initialize__main_;
        wtime_end_M2L_initialize__main_ = m2l_engine_.wtime_end_initialize__main_;

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_PARTICLE_MESH)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok_last @ calcForceParticleMesh() [PMMM]" << std::endl;
#endif
    }

   public:
    ParticleMeshMultipole(const std::string id = "", const CommInfo &comm_info = CommInfo()) : is_dummy(false) {
        id_ = id;
        if (comm_info.isCommNull())
            comm_info_.setCommunicator();
        else
            comm_info_ = comm_info;
        first_call_by_initialize_ = true;
        msum_ = 0.0;
        quad0_ = 0.0;
    }
    ~ParticleMeshMultipole() {
        //std::cout<<"PM3 destructor called"<<std::endl;
        }

    void initialize(const S32 p, S32 p_spj2mm=0, const U32 fftw_planning_rigor_flag = FFTW_MEASURE, const bool use_mpifft_if_possible = true,
                    const S32 fft_size_crit = 10000) {
        wtime_start_initialize_ = GetWtime();
        // Check error
        assert(first_call_by_initialize_);
        assert(p >= 0);
        assert(p_spj2mm >= 0);
        if ((fftw_planning_rigor_flag != FFTW_ESTIMATE) && (fftw_planning_rigor_flag != FFTW_MEASURE) && (fftw_planning_rigor_flag != FFTW_PATIENT) &&
            (fftw_planning_rigor_flag != FFTW_EXHAUSTIVE)) {
            if (comm_info_.getRank() == 0) {
                std::stringstream err_msg;
                err_msg << "The value of argument fftw_planning_rigor_flag must be\n"
                        << "one of the following:\n"
                        << "--------------------------\n"
                        << "    FFTW_ESTIMATE   (" << FFTW_ESTIMATE << ")\n"
                        << "    FFTW_MEASURE    (" << FFTW_MEASURE << ")\n"
                        << "    FFTW_PATIENT    (" << FFTW_PATIENT << ")\n"
                        << "    FFTW_EXHAUSTIVE (" << FFTW_EXHAUSTIVE << ")\n"
                        << "--------------------------\n"
                        << "On the other hand, the value you specified for\n"
                        << "fftw_planning_rigor_flag was " << fftw_planning_rigor_flag << ".\n"
                        << "The following are possible causes:\n"
                        << "(1) You just passed an invalid value.\n"
                        << "(2) The value is valid but its position in the argument list\n"
                        << "    is wrong.";
                std::cout << err_msg.str() << std::endl;
            }
            Abort(-1);
        }

        // Make tables for Rlm class
        {
            Rlm<real_t> tmp;
            tmp.make_table(p);  // for P2M & L2P
        }
        // [Note (tag: #144eddb3)]
        //   (1) Tables required by Slm class are created when needed.
        //       (see transform() in M2L_Engine.hpp &
        //        calcGreenFunctionInRealSpace() in green_function.hpp)
        //   (2) In the current implementation, we do not need to care about p_spj2mm
        //       because MomentMultipoleGeometricCenter in tree.hpp is implemented
        //       using fmm_org.hpp, where p is a template parameter and hence necessary
        //       tables are made at the compilation time.

        // Store parameters
        param_.p = p;
        //p_spj2mm_ = p_spj2mm;
	p_spj2mm_ = p;
        fftw_planning_rigor_flag_ = fftw_planning_rigor_flag;
        use_mpifft_if_possible_ = use_mpifft_if_possible;
        fft_size_crit_ = fft_size_crit;

        // Update the flag
        first_call_by_initialize_ = false;

        wtime_end_initialize_ = GetWtime();
    }

    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE>
    void calcForceAll(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
                      const DomainInfo &dinfo, const bool clear_force = true, const S32 n_thread = Comm::getMaxThreads()) {
        wtime_start_calc_force_all_ = GetWtime();
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL)
        constexpr S32 rank_wld_debug_print = 0;
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok1 @ calcForceAll() [PMMM]" << std::endl;
#endif

        setParam(tree, dinfo);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok2 @ calcForceAll() [PMMM]" << std::endl;
#endif

        setCell(tree, n_thread);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok3 @ calcForceAll() [PMMM]" << std::endl;
#endif

        setIParticleInfoToCell(tree, dinfo, n_thread);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok4 @ calcForceAll() [PMMM]" << std::endl;
#endif

        calcTotalChargeAndDispersion(tree, dinfo, n_thread);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok5 @ calcForceAll() [PMMM]" << std::endl;
#endif

        calcMultipoleMoment(tree, dinfo, n_thread);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok6 @ calcForceAll() [PMMM]" << std::endl;
#endif

        collectMultipoleMoment(tree, dinfo, n_thread);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok7 @ calcForceAll() [PMMM]" << std::endl;
#endif

        calcForceParticleMesh(tree, dinfo, clear_force, n_thread);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok_last @ calcForceAll() [PMMM]" << std::endl;
#endif
        wtime_end_calc_force_all_ = GetWtime();
        time_profile_.PMMM__calc_force_all += (wtime_end_calc_force_all_ - wtime_start_calc_force_all_);
    }

    template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, enum CALC_DISTANCE_TYPE CALC_DISTANCE_TYPE, class Tpsys>
    void calcForceAllAndWriteBack(
        TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj, CALC_DISTANCE_TYPE> &tree,
        const DomainInfo &dinfo, Tpsys &psys, const bool clear_force = true, const S32 n_thread = Comm::getMaxThreads()) {
        //comm_info_.barrier();
        //if (Comm::getRank() == 0) std::cout << "ok1 @ calcForceAllAndWriteBack() [PMMM]" << std::endl;
        calcForceAll(tree, dinfo, clear_force, n_thread);
        //comm_info_.barrier();
        //if (Comm::getRank() == 0) std::cout << "ok2 @ calcForceAllAndWriteBack() [PMMM]" << std::endl;
        wtime_start_write_back_ = GetWtime();
#if 1
        tree.copyForceOriginalOrder();
#else
        tree.copyForceOriginalOrder(n_thread);
#endif

        // Comm::barrier();
        // std::cerr<<"check c"<<std::endl;

        PS_OMP(omp parallel for num_threads(n_thread))
        for (S32 i = 0; i < n_loc_tot_; i++) psys[i].copyFromForcePMMM(tree.getForce(i));
        wtime_end_write_back_ = GetWtime();
        time_profile_.PMMM__write_back += (wtime_end_write_back_ - wtime_start_write_back_);
    }

    void clearTimeProfile() {
        m2l_engine_.clearTimeProfile();
        time_profile_.clear();
    }

    void clearWallTime() {
        m2l_engine_.clearWallTime();

        wtime_start_GFC_initialize_ = 0;
        wtime_end_GFC_initialize_ = 0;
        wtime_start_GFC_calc_gf_r_ = 0;
        wtime_end_GFC_calc_gf_r_ = 0;
        wtime_start_GFC_redist_gf_r_ = 0;
        wtime_end_GFC_redist_gf_r_ = 0;
        wtime_start_GFC_calc_gf_k_ = 0;
        wtime_end_GFC_calc_gf_k_ = 0;
        wtime_start_GFC_redist_gf_k_ = 0;
        wtime_end_GFC_redist_gf_k_ = 0;

        wtime_start_GFC_initialize__copy_ = 0;
        wtime_end_GFC_initialize__copy_ = 0;
        wtime_start_GFC_initialize__main_ = 0;
        wtime_end_GFC_initialize__main_ = 0;

        wtime_start_M2L_initialize_ = 0;
        wtime_end_M2L_initialize_ = 0;
        wtime_start_M2L_preproc_mm_r_comm_ = 0;
        wtime_end_M2L_preproc_mm_r_comm_ = 0;
        wtime_start_M2L_redist_mm_r_ = 0;
        wtime_end_M2L_redist_mm_r_ = 0;
        wtime_start_M2L_postproc_mm_r_comm_ = 0;
        wtime_end_M2L_postproc_mm_r_comm_ = 0;
        wtime_start_M2L_mm_r_to_mm_k_ = 0;
        wtime_end_M2L_mm_r_to_mm_k_ = 0;
        wtime_start_M2L_gather_mm_k_trans_ = 0;
        wtime_end_M2L_gather_mm_k_trans_ = 0;
        wtime_start_M2L_transform_ = 0;
        wtime_end_M2L_transform_ = 0;
        wtime_start_M2L_scatter_le_k_trans_ = 0;
        wtime_end_M2L_scatter_le_k_trans_ = 0;
        wtime_start_M2L_le_k_to_le_r_ = 0;
        wtime_end_M2L_le_k_to_le_r_ = 0;
        wtime_start_M2L_preproc_le_r_comm_ = 0;
        wtime_end_M2L_preproc_le_r_comm_ = 0;
        wtime_start_M2L_redist_le_r_ = 0;
        wtime_end_M2L_redist_le_r_ = 0;
        wtime_start_M2L_postproc_le_r_comm_ = 0;
        wtime_end_M2L_postproc_le_r_comm_ = 0;

        wtime_start_M2L_initialize__copy_ = 0;
        wtime_end_M2L_initialize__copy_ = 0;
        wtime_start_M2L_initialize__main_ = 0;
        wtime_end_M2L_initialize__main_ = 0;

        wtime_start_calc_force_all_ = 0;
        wtime_end_calc_force_all_ = 0;
        wtime_start_write_back_ = 0;
        wtime_end_write_back_ = 0;
        wtime_start_initialize_ = 0;
        wtime_end_initialize_ = 0;
        wtime_start_set_cell_ = 0;
        wtime_end_set_cell_ = 0;
        wtime_start_set_ip_info_to_cell_ = 0;
        wtime_end_set_ip_info_to_cell_ = 0;
        wtime_start_calc_msum_and_quad0_ = 0;
        wtime_end_calc_msum_and_quad0_ = 0;
        wtime_start_calc_msum_and_quad0__main_ = 0;
        wtime_end_calc_msum_and_quad0__main_ = 0;
        wtime_start_calc_msum_and_quad0__allreduce_ = 0;
        wtime_end_calc_msum_and_quad0__allreduce_ = 0;
        wtime_start_calc_multipole_moment_ = 0;
        wtime_end_calc_multipole_moment_ = 0;
        wtime_start_collect_multipole_moment_ = 0;
        wtime_end_collect_multipole_moment_ = 0;
        wtime_start_dipole_correction_ = 0;
        wtime_end_dipole_correction_ = 0;
        wtime_start_L2P_ = 0;
        wtime_end_L2P_ = 0;
    }

    void clearCounterAll() {
        clearWallTime();
        clearTimeProfile();
    }

    TimeProfile getTimeProfile() const { return time_profile_ + m2l_engine_.getTimeProfile(); }
};

class ParticleMeshMultipoleDummy {
    // This class has the same interface as ParticleMeshMultipole class,
    // but do nothing. This class may be useful when FDPS users want to
    // make it possible to use the same codes with and without the PMMM
    // feature.
   public:
    using class_type = TagClassTypeParticleMeshMultipole;
    const bool is_dummy;

   private:
    std::string id_;

    CommInfo comm_info_;

    TimeProfile time_profile_;

    F64 wtime_start_GFC_initialize_;
    F64 wtime_end_GFC_initialize_;
    F64 wtime_start_GFC_calc_gf_r_;
    F64 wtime_end_GFC_calc_gf_r_;
    F64 wtime_start_GFC_redist_gf_r_;
    F64 wtime_end_GFC_redist_gf_r_;
    F64 wtime_start_GFC_calc_gf_k_;
    F64 wtime_end_GFC_calc_gf_k_;
    F64 wtime_start_GFC_redist_gf_k_;
    F64 wtime_end_GFC_redist_gf_k_;

    F64 wtime_start_GFC_initialize__copy_;
    F64 wtime_end_GFC_initialize__copy_;
    F64 wtime_start_GFC_initialize__main_;
    F64 wtime_end_GFC_initialize__main_;

    F64 wtime_start_M2L_initialize_;
    F64 wtime_end_M2L_initialize_;
    F64 wtime_start_M2L_preproc_mm_r_comm_;
    F64 wtime_end_M2L_preproc_mm_r_comm_;
    F64 wtime_start_M2L_redist_mm_r_;
    F64 wtime_end_M2L_redist_mm_r_;
    F64 wtime_start_M2L_postproc_mm_r_comm_;
    F64 wtime_end_M2L_postproc_mm_r_comm_;
    F64 wtime_start_M2L_mm_r_to_mm_k_;
    F64 wtime_end_M2L_mm_r_to_mm_k_;
    F64 wtime_start_M2L_gather_mm_k_trans_;
    F64 wtime_end_M2L_gather_mm_k_trans_;
    F64 wtime_start_M2L_transform_;
    F64 wtime_end_M2L_transform_;
    F64 wtime_start_M2L_scatter_le_k_trans_;
    F64 wtime_end_M2L_scatter_le_k_trans_;
    F64 wtime_start_M2L_le_k_to_le_r_;
    F64 wtime_end_M2L_le_k_to_le_r_;
    F64 wtime_start_M2L_preproc_le_r_comm_;
    F64 wtime_end_M2L_preproc_le_r_comm_;
    F64 wtime_start_M2L_redist_le_r_;
    F64 wtime_end_M2L_redist_le_r_;
    F64 wtime_start_M2L_postproc_le_r_comm_;
    F64 wtime_end_M2L_postproc_le_r_comm_;

    F64 wtime_start_M2L_initialize__copy_;
    F64 wtime_end_M2L_initialize__copy_;
    F64 wtime_start_M2L_initialize__main_;
    F64 wtime_end_M2L_initialize__main_;

    F64 wtime_start_calc_force_all_;
    F64 wtime_end_calc_force_all_;
    F64 wtime_start_write_back_;
    F64 wtime_end_write_back_;
    F64 wtime_start_initialize_;
    F64 wtime_end_initialize_;
    F64 wtime_start_set_cell_;
    F64 wtime_end_set_cell_;
    F64 wtime_start_set_ip_info_to_cell_;
    F64 wtime_end_set_ip_info_to_cell_;
    F64 wtime_start_calc_msum_and_quad0_;
    F64 wtime_end_calc_msum_and_quad0_;
    F64 wtime_start_calc_msum_and_quad0__main_;
    F64 wtime_end_calc_msum_and_quad0__main_;
    F64 wtime_start_calc_msum_and_quad0__allreduce_;
    F64 wtime_end_calc_msum_and_quad0__allreduce_;
    F64 wtime_start_calc_multipole_moment_;
    F64 wtime_end_calc_multipole_moment_;
    F64 wtime_start_collect_multipole_moment_;
    F64 wtime_end_collect_multipole_moment_;
    F64 wtime_start_dipole_correction_;
    F64 wtime_end_dipole_correction_;
    F64 wtime_start_L2P_;
    F64 wtime_end_L2P_;

   public:
    ParticleMeshMultipoleDummy(const std::string id = "", const CommInfo &comm_info = CommInfo()) : is_dummy(true) {
        id_ = id;
        if (comm_info.isCommNull())
            comm_info_.setCommunicator();
        else
            comm_info_ = comm_info;
        time_profile_.clear();
    }
    ~ParticleMeshMultipoleDummy() {}

    void initialize(const S32 p, const S32 p_spj2mm=0, const U32 fftw_planning_rigor_flag = FFTW_MEASURE, const bool use_mpifft_if_possible = true,
                    const S32 fft_size_crit = 10000) {}

    template <class Ttree>
    void calcForceAll(Ttree &tree, const DomainInfo &dinfo, const bool clear_force = true) {}

    template <class Ttree, class Tpsys>
    void calcForceAllAndWriteBack(Ttree &tree, const DomainInfo &dinfo, Tpsys &psys, const bool clear_force = true) {}
    void clearTimeProfile() { time_profile_.clear(); }

    void clearWallTime() {
        wtime_start_GFC_initialize_ = 0;
        wtime_end_GFC_initialize_ = 0;
        wtime_start_GFC_calc_gf_r_ = 0;
        wtime_end_GFC_calc_gf_r_ = 0;
        wtime_start_GFC_redist_gf_r_ = 0;
        wtime_end_GFC_redist_gf_r_ = 0;
        wtime_start_GFC_calc_gf_k_ = 0;
        wtime_end_GFC_calc_gf_k_ = 0;
        wtime_start_GFC_redist_gf_k_ = 0;
        wtime_end_GFC_redist_gf_k_ = 0;

        wtime_start_GFC_initialize__copy_ = 0;
        wtime_end_GFC_initialize__copy_ = 0;
        wtime_start_GFC_initialize__main_ = 0;
        wtime_end_GFC_initialize__main_ = 0;

        wtime_start_M2L_initialize_ = 0;
        wtime_end_M2L_initialize_ = 0;
        wtime_start_M2L_preproc_mm_r_comm_ = 0;
        wtime_end_M2L_preproc_mm_r_comm_ = 0;
        wtime_start_M2L_redist_mm_r_ = 0;
        wtime_end_M2L_redist_mm_r_ = 0;
        wtime_start_M2L_postproc_mm_r_comm_ = 0;
        wtime_end_M2L_postproc_mm_r_comm_ = 0;
        wtime_start_M2L_mm_r_to_mm_k_ = 0;
        wtime_end_M2L_mm_r_to_mm_k_ = 0;
        wtime_start_M2L_gather_mm_k_trans_ = 0;
        wtime_end_M2L_gather_mm_k_trans_ = 0;
        wtime_start_M2L_transform_ = 0;
        wtime_end_M2L_transform_ = 0;
        wtime_start_M2L_scatter_le_k_trans_ = 0;
        wtime_end_M2L_scatter_le_k_trans_ = 0;
        wtime_start_M2L_le_k_to_le_r_ = 0;
        wtime_end_M2L_le_k_to_le_r_ = 0;
        wtime_start_M2L_preproc_le_r_comm_ = 0;
        wtime_end_M2L_preproc_le_r_comm_ = 0;
        wtime_start_M2L_redist_le_r_ = 0;
        wtime_end_M2L_redist_le_r_ = 0;
        wtime_start_M2L_postproc_le_r_comm_ = 0;
        wtime_end_M2L_postproc_le_r_comm_ = 0;

        wtime_start_M2L_initialize__copy_ = 0;
        wtime_end_M2L_initialize__copy_ = 0;
        wtime_start_M2L_initialize__main_ = 0;
        wtime_end_M2L_initialize__main_ = 0;

        wtime_start_calc_force_all_ = 0;
        wtime_end_calc_force_all_ = 0;
        wtime_start_write_back_ = 0;
        wtime_end_write_back_ = 0;
        wtime_start_initialize_ = 0;
        wtime_end_initialize_ = 0;
        wtime_start_set_cell_ = 0;
        wtime_end_set_cell_ = 0;
        wtime_start_set_ip_info_to_cell_ = 0;
        wtime_end_set_ip_info_to_cell_ = 0;
        wtime_start_calc_msum_and_quad0_ = 0;
        wtime_end_calc_msum_and_quad0_ = 0;
        wtime_start_calc_msum_and_quad0__main_ = 0;
        wtime_end_calc_msum_and_quad0__main_ = 0;
        wtime_start_calc_msum_and_quad0__allreduce_ = 0;
        wtime_end_calc_msum_and_quad0__allreduce_ = 0;
        wtime_start_calc_multipole_moment_ = 0;
        wtime_end_calc_multipole_moment_ = 0;
        wtime_start_collect_multipole_moment_ = 0;
        wtime_end_collect_multipole_moment_ = 0;
        wtime_start_dipole_correction_ = 0;
        wtime_end_dipole_correction_ = 0;
        wtime_start_L2P_ = 0;
        wtime_end_L2P_ = 0;
    }

    void clearCounterAll() {
        clearWallTime();
        clearTimeProfile();
    }

    TimeProfile getTimeProfile() const { return time_profile_; }
};

}  // namespace ParticleMeshMultipole
}  // namespace ParticleSimulator
