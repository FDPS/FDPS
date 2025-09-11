#pragma once
#include <complex>
#include <unordered_map>
#include "../ps_defs.hpp"
#include "fmm.hpp"
#include "cutoff.hpp"
#include "particle_mesh_multipole_defs.hpp"

namespace ParticleSimulator {
namespace ParticleMeshMultipole {

class GreenFunctionCalculator {
   public:
    using real_t = double;
    using cplx_t = std::complex<real_t>;

    Parameters param_;
    Parameters param_prev_;
    S32 n_thread_;
    S32 n_thread_prev_;

    // Variables for communication
    CommInfo comm_info_;
    S32 n_proc_in_parent_group_;
    S32 rank_in_parent_group_;
    S32 n_group_;
    S32 n_proc_in_my_group_;
    S32 rank_in_my_group_;
    S32 n_proc_min_in_group_;
    std::vector<S32> rank_start_, rank_end_;
    std::vector<S32> lm_start_, lm_end_;
    std::vector<S32> local_0_start_, local_0_end_;  // z of (z,y,x)
    std::vector<S32> local_1_start_, local_1_end_;  // y of (y,z,x)
    S32 n_proc_for_fft_;
    S32 idx_to_my_group_;
    bool is_mpifft_usable_;
    bool use_mpifft_;
    bool use_ompfft_;

    S32 n_proc_tot_gf_r_calc_;
    S32 n_proc_gf_r_calc_[3];
    std::vector<S32> local_start_gf_r_calc_[3], local_end_gf_r_calc_[3];
    // TODO: change to std::vector<S32ort>
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Group parent_group_;
    MPI_Comm parent_comm_;
    MPI_Datatype mpi_real_t_;
    MPI_Datatype mpi_cplx_t_;
    MPI_Group group_all_;
    MPI_Group group_fft_;
    MPI_Group group_int_;  // inter-group
    MPI_Comm comm_all_;
    MPI_Comm comm_fft_;
    MPI_Comm comm_int_;
#endif
    // Variables to store green function
    Buffer<real_t> gf_r_tmp_part_;
    Buffer<real_t> gf_r_tmp_;
    Buffer<cplx_t> gf_k_tmp_;

    // Timing measurement
    TimeProfile time_profile_;
    F64 wtime_start_initialize_;
    F64 wtime_end_initialize_;
    F64 wtime_start_calc_gf_r_;
    F64 wtime_end_calc_gf_r_;
    F64 wtime_start_redist_gf_r_;
    F64 wtime_end_redist_gf_r_;
    F64 wtime_start_calc_gf_k_;
    F64 wtime_end_calc_gf_k_;
    F64 wtime_start_redist_gf_k_;
    F64 wtime_end_redist_gf_k_;

    F64 wtime_start_initialize__copy_;
    F64 wtime_end_initialize__copy_;
    F64 wtime_start_initialize__main_;
    F64 wtime_end_initialize__main_;

    GreenFunctionCalculator() {
        n_thread_ = -1;
        n_thread_prev_ = -1;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        group_all_ = MPI_GROUP_NULL;
        group_fft_ = MPI_GROUP_NULL;
        group_int_ = MPI_GROUP_NULL;
        comm_all_ = MPI_COMM_NULL;
        comm_fft_ = MPI_COMM_NULL;
        comm_int_ = MPI_COMM_NULL;
#endif
    }

    void finalize() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if (group_all_ != MPI_GROUP_NULL) MPI_Group_free(&group_all_);
        if (group_fft_ != MPI_GROUP_NULL) MPI_Group_free(&group_fft_);
        if (group_int_ != MPI_GROUP_NULL) MPI_Group_free(&group_int_);
        if (comm_all_ != MPI_COMM_NULL) MPI_Comm_free(&comm_all_);
        if (comm_fft_ != MPI_COMM_NULL) MPI_Comm_free(&comm_fft_);
        if (comm_int_ != MPI_COMM_NULL) MPI_Comm_free(&comm_int_);
#endif
    }
#if 1
    void initialize(const Parameters &param, const bool use_mpifft_if_possible, const bool use_ompfft, const CommInfo &comm_info, const S32 n_thread,
                    const bool debug_flag = false) {
        F64 wtime_tmp = GetWtime();
        wtime_start_initialize_ = wtime_tmp;
        // Copy parameters, etc.
        wtime_start_initialize__copy_ = wtime_tmp;
        param_prev_ = param_;  // save the previous parameters
        param_ = param;
        use_ompfft_ = use_ompfft;
        comm_info_ = comm_info;
        n_thread_prev_ = n_thread_;
        n_thread_ = n_thread;
        wtime_tmp = GetWtime();
        wtime_end_initialize__copy_ = wtime_tmp;
        // Initialize
        wtime_start_initialize__main_ = wtime_tmp;
        if ((param_ != param_prev_) || (n_thread_ != n_thread_prev_) || (use_mpifft_ != use_mpifft_if_possible)) {
            // [TODO] There is the same problem as member function `initialize`
            //        of M2L_Engine class. Refer to the comment described in
            //        M2L_Engine::initialize() for more details.
            // PS::DebugUtils::fout() << "GFC initialize performed." << std::endl;
            const S32 n_mm_compo = (2 * param_.p + 1) * (2 * param_.p + 1);
            const S32vec n_cell = param_.n_cell;
            const S32 bc = param_.bc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            // Set parent_comm
            parent_comm_ = comm_info_.getCommunicator();
            // Get information of parent_comm
            S32 n_proc, my_rank;
            MPI_Comm_size(parent_comm_, &n_proc);
            MPI_Comm_rank(parent_comm_, &my_rank);
            MPI_Comm_group(parent_comm_, &parent_group_);
            n_proc_in_parent_group_ = n_proc;
            rank_in_parent_group_ = my_rank;
            // Get MPI_Datatype for ptrdiff_t, real_t, cplx_t
            mpi_real_t_ = GetDataType<real_t>();
            mpi_cplx_t_ = GetDataType<cplx_t>();
            // Divide parent_comm into groups
            if (n_proc >= n_mm_compo) {
                n_group_ = n_mm_compo;
                // Set rank_start_, rank_end_
                n_proc_min_in_group_ = n_proc / n_group_;
                S32 n_rem = n_proc - n_proc_min_in_group_ * n_group_;
                rank_start_.resize(n_group_);
                rank_end_.resize(n_group_);
                rank_start_[0] = 0;
                for (S32 i = 0; i < n_group_; i++) {
                    S32 n_proc_in_group = n_proc_min_in_group_;
                    if (i < n_rem) n_proc_in_group++;
                    rank_end_[i] = rank_start_[i] + (n_proc_in_group - 1);
                    if (i < n_group_ - 1) rank_start_[i + 1] = rank_end_[i] + 1;
                }
                if (debug_flag) {
                    if (my_rank == 0) {
                        std::cout << "### Information of groups for G2G ###" << std::endl;
                        std::cout << "(n_proc = " << n_proc << ")" << std::endl;
                        for (S32 i = 0; i < n_group_; i++) {
                            std::cout << "i = " << i << " rank_start = " << rank_start_[i] << " rank_end = " << rank_end_[i] << std::endl;
                        }
                    }
                }
                // Calculate idx_to_my_group, n_proc_in_my_group, rank_in_my_group
                for (S32 i = 0; i < n_group_; i++)
                    if ((rank_start_[i] <= my_rank) && (my_rank <= rank_end_[i])) {
                        idx_to_my_group_ = i;
                        n_proc_in_my_group_ = rank_end_[i] - rank_start_[i] + 1;
                        rank_in_my_group_ = my_rank - rank_start_[i];
                    }
                // Set lm_start_, lm_end_
                lm_start_.resize(n_group_);
                lm_end_.resize(n_group_);
                for (S32 i = 0; i < n_group_; i++) {
                    lm_start_[i] = i;
                    lm_end_[i] = i;
                    // Each group is responsible for a single (l,m)
                }
                // Make group_all & comm_all
                if (group_all_ != MPI_GROUP_NULL) MPI_Group_free(&group_all_);
                if (comm_all_ != MPI_COMM_NULL) MPI_Comm_free(&comm_all_);
                std::vector<S32> ranks;
                S32 rnk_start = rank_start_[idx_to_my_group_];
                S32 rnk_end = rank_end_[idx_to_my_group_];
                for (S32 rnk = rnk_start; rnk <= rnk_end; rnk++) ranks.push_back(rnk);
                MPI_Group_incl(parent_group_, n_proc_in_my_group_, &ranks[0], &group_all_);
                MPI_Comm_create(parent_comm_, group_all_, &comm_all_);
                // Calculate n_proc_for_fft
                S32 n_proc_max_for_fft;
                if (bc == BOUNDARY_CONDITION_OPEN)
                    n_proc_max_for_fft = n_cell.z;
                else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
                    n_proc_max_for_fft = n_cell.z / 2;
                else
                    assert(false);
                if (n_proc_min_in_group_ <= n_proc_max_for_fft)
                    n_proc_for_fft_ = n_proc_min_in_group_;
                else
                    n_proc_for_fft_ = n_proc_max_for_fft;
                // Make group_fft, comm_fft
                if (n_proc_for_fft_ > 1)
                    is_mpifft_usable_ = true;
                else
                    is_mpifft_usable_ = false;
                if (is_mpifft_usable_)
                    use_mpifft_ = use_mpifft_if_possible;  // specified by user.
                else
                    use_mpifft_ = false;
                if (use_mpifft_) {
                    if (!FFTW_Status::is_fftw_mpi_init_called()) {
                        fftw_mpi_init();
                        FFTW_Status::is_fftw_mpi_init_called() = true;
                    }
                    if (group_fft_ != MPI_GROUP_NULL) MPI_Group_free(&group_fft_);
                    if (comm_fft_ != MPI_COMM_NULL) MPI_Comm_free(&comm_fft_);
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        ranks.resize(n_proc_for_fft_);
                        for (S32 i = 0; i < n_proc_for_fft_; i++) ranks[i] = rank_start_[idx_to_my_group_] + i;
                        MPI_Group_incl(parent_group_, n_proc_for_fft_, &ranks[0], &group_fft_);
                        MPI_Comm_create(parent_comm_, group_fft_, &comm_fft_);
                    } else {
                        group_fft_ = MPI_GROUP_NULL;
                        MPI_Comm_create(parent_comm_, MPI_GROUP_EMPTY, &comm_fft_);
                    }
                    // In this case, we need to calculate local_0_start & local_0_end
                    // using FFTW API.
                    // Resize local_?_start, local_?_end
                    local_0_start_.resize(n_proc_for_fft_);
                    local_0_end_.resize(n_proc_for_fft_);
                    local_1_start_.resize(n_proc_for_fft_);
                    local_1_end_.resize(n_proc_for_fft_);
                    // Calculate loca_?_start, local_?_end
                    if (comm_fft_ != MPI_COMM_NULL) {
                        // ptrdiff_t alloc_local, local_n0, local_0_start;
                        ptrdiff_t local_n0, local_0_start;
                        ptrdiff_t local_n1, local_1_start;
                        if (bc == BOUNDARY_CONDITION_OPEN) {
                            // alloc_local = fftw_mpi_local_size_3d_transposed(2 * n_cell.z, 2 * n_cell.y, n_cell.x + 1, comm_fft_, &local_n0,
                            //                                                 &local_0_start, &local_n1, &local_1_start);
                            fftw_mpi_local_size_3d_transposed(2 * n_cell.z, 2 * n_cell.y, n_cell.x + 1, comm_fft_, &local_n0, &local_0_start,
                                                              &local_n1, &local_1_start);
                        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                            // alloc_local = fftw_mpi_local_size_3d_transposed(n_cell.z, n_cell.y, n_cell.x / 2 + 1, comm_fft_, &local_n0,
                            //&local_0_start, &local_n1, &local_1_start);
                            fftw_mpi_local_size_3d_transposed(n_cell.z, n_cell.y, n_cell.x / 2 + 1, comm_fft_, &local_n0, &local_0_start, &local_n1,
                                                              &local_1_start);
                        } else {
                            if (my_rank == 0) std::cout << "This boundary condition is not supported yet." << std::endl;
                            Abort(-1);
                        }
                        S32 sendbuf[4];
                        sendbuf[0] = local_0_start;
                        sendbuf[1] = local_n0;
                        sendbuf[2] = local_1_start;
                        sendbuf[3] = local_n1;
                        S32 *recvbuf = new S32[4 * n_proc_for_fft_];
                        MPI_Gather(sendbuf, 4, MPI_INT, recvbuf, 4, MPI_INT, 0, comm_fft_);
                        if (rank_in_my_group_ == 0) {
                            for (S32 i = 0; i < n_proc_for_fft_; i++) {
                                local_0_start_[i] = recvbuf[4 * i];
                                local_0_end_[i] = local_0_start_[i] + recvbuf[4 * i + 1] - 1;
                                local_1_start_[i] = recvbuf[4 * i + 2];
                                local_1_end_[i] = local_1_start_[i] + recvbuf[4 * i + 3] - 1;
                            }
                        }
                        delete[] recvbuf;
                    }
                    MPI_Bcast(&local_0_start_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                    MPI_Bcast(&local_0_end_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                    MPI_Bcast(&local_1_start_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                    MPI_Bcast(&local_1_end_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                    // Check
                    if (debug_flag) {
                        if (rank_in_parent_group_ == 0) {
                            std::cout << "### Information of decomposition of FFTW ###" << std::endl;
                            std::cout << "(n_cell.z = " << n_cell.z << ")" << std::endl;
                            for (S32 i = 0; i < n_proc_for_fft_; i++) {
                                std::cout << "i = " << i << " local_0_start = " << local_0_start_[i] << " local_0_end = " << local_0_end_[i]
                                          << " local_1_start = " << local_1_start_[i] << " local_1_end = " << local_1_end_[i] << std::endl;
                            }
                        }
                    }
                } else {
                    n_proc_for_fft_ = 1;  // reset
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                    if (!FFTW_Status::is_fftw_init_threads_called()) {
                        fftw_init_threads();
                        FFTW_Status::is_fftw_init_threads_called() = true;
                    }
#endif
                    group_fft_ = MPI_GROUP_NULL;
                    comm_fft_ = MPI_COMM_NULL;
                    local_0_start_.resize(1);
                    local_0_end_.resize(1);
                    if (bc == BOUNDARY_CONDITION_OPEN) {
                        local_0_start_[0] = 0;
                        local_0_end_[0] = 2 * n_cell.z - 1;
                    } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                        local_0_start_[0] = 0;
                        local_0_end_[0] = n_cell.z - 1;
                    }
                }
                // Make group_int, comm_int
                if (group_int_ != MPI_GROUP_NULL) MPI_Group_free(&group_int_);
                if (comm_int_ != MPI_COMM_NULL) MPI_Comm_free(&comm_int_);
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    ranks.resize(n_group_);
                    for (S32 i = 0; i < n_group_; i++) ranks[i] = rank_start_[i] + rank_in_my_group_;
                    MPI_Group_incl(parent_group_, n_group_, &ranks[0], &group_int_);
                    MPI_Comm_create(parent_comm_, group_int_, &comm_int_);
                } else {
                    group_int_ = MPI_GROUP_NULL;
                    MPI_Comm_create(parent_comm_, MPI_GROUP_EMPTY, &comm_int_);
                }
            } else {
                n_group_ = n_proc;
                n_proc_in_my_group_ = 1;
                rank_in_my_group_ = 0;
                n_proc_min_in_group_ = 1;
                n_proc_for_fft_ = 1;
                idx_to_my_group_ = my_rank;
                is_mpifft_usable_ = false;
                use_mpifft_ = false;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                if (!FFTW_Status::is_fftw_init_threads_called()) {
                    fftw_init_threads();
                    FFTW_Status::is_fftw_init_threads_called() = true;
                }
#endif
                // Set rank_start_, rank_end_
                rank_start_.resize(n_group_);
                rank_end_.resize(n_group_);
                for (S32 i = 0; i < n_group_; i++) {
                    rank_start_[i] = i;
                    rank_end_[i] = i;
                }
                // Set lm_start_, lm_end_
                const S32 n_mm_compo_per_proc_min = n_mm_compo / n_proc;
                S32 n_rem = n_mm_compo - n_mm_compo_per_proc_min * n_proc;
                lm_start_.resize(n_proc);
                lm_end_.resize(n_proc);
                lm_start_[0] = 0;
                for (S32 i = 0; i < n_proc; i++) {
                    S32 n_mm_compo_per_proc = n_mm_compo_per_proc_min;
                    if (i < n_rem) n_mm_compo_per_proc++;
                    lm_end_[i] = lm_start_[i] + (n_mm_compo_per_proc - 1);
                    if (i < n_proc - 1) lm_start_[i + 1] = lm_end_[i] + 1;
                }
                if (debug_flag) {
                    if (rank_in_parent_group_ == 0) {
                        std::cout << "### Information of decomposition of (l,m) ###" << std::endl;
                        std::cout << "(n_mm_compo = " << n_mm_compo << ")" << std::endl;
                        for (S32 i = 0; i < n_proc; i++) {
                            std::cout << "i = " << i << " local_0_start = " << lm_start_[i] << " local_0_end = " << lm_end_[i] << std::endl;
                        }
                    }
                }
                // Set local_0_start_, local_0_end_
                group_fft_ = MPI_GROUP_NULL;
                comm_fft_ = MPI_COMM_NULL;
                local_0_start_.resize(n_proc_for_fft_);
                local_0_end_.resize(n_proc_for_fft_);
                if (bc == BOUNDARY_CONDITION_OPEN) {
                    local_0_start_[0] = 0;
                    local_0_end_[0] = 2 * n_cell.z - 1;
                } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                    local_0_start_[0] = 0;
                    local_0_end_[0] = n_cell.z - 1;
                }
                // Note that we do not need to set local_1_start_
                // and local_1_end_ because of use_mpifft = false.
                // (hence, both arrays are not accessed)

                // Make group_int, comm_int
                group_int_ = parent_group_;
                comm_int_ = parent_comm_;
            }
#else  // PARTICLE_SIMULATOR_MPI_PARALLEL
       // Sequential execution/ OpenMP parallelization
            n_proc_in_parent_group_ = 1;
            rank_in_parent_group_ = 0;
            n_group_ = 1;
            n_proc_in_my_group_ = 1;
            rank_in_my_group_ = 0;
            n_proc_min_in_group_ = 1;
            n_proc_for_fft_ = 1;
            idx_to_my_group_ = 0;
            is_mpifft_usable_ = false;
            use_mpifft_ = false;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            if (!FFTW_Status::is_fftw_init_threads_called()) {
                fftw_init_threads();
                FFTW_Status::is_fftw_init_threads_called() = true;
            }
#endif
            // Set rank_start_(n_group_), rank_end_(n_group_) [idx_to_my_group_]
            rank_start_.resize(n_group_);
            rank_end_.resize(n_group_);
            rank_start_[0] = 0;
            rank_end_[0] = 0;
            // Set lm_start_(n_group_), lm_end_(n_group_) [idx_to_my_group_]
            lm_start_.resize(n_group_);
            lm_end_.resize(n_group_);
            lm_start_[0] = 0;
            lm_end_[0] = n_mm_compo - 1;
            // Set local_0_start__, local_0_end_
            local_0_start_.resize(n_proc_for_fft_);
            local_0_end_.resize(n_proc_for_fft_);
            if (bc == BOUNDARY_CONDITION_OPEN) {
                local_0_start_[0] = 0;
                local_0_end_[0] = 2 * n_cell.z - 1;
            } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                local_0_start_[0] = 0;
                local_0_end_[0] = n_cell.z - 1;
            }
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL

            // Set the size information of gf_r_tmp_, gf_k_tmp_
            if (bc == BOUNDARY_CONDITION_OPEN) {
                // Set the size information of gf_r_tmp_
                gf_r_tmp_.start_glb_[0] = 0;
                gf_r_tmp_.start_glb_[1] = 0;
                gf_r_tmp_.start_glb_[2] = 0;
                gf_r_tmp_.start_glb_[3] = 0;
                gf_r_tmp_.end_glb_[0] = n_mm_compo - 1;
                gf_r_tmp_.end_glb_[1] = 2 * n_cell.x - 1;
                gf_r_tmp_.end_glb_[2] = 2 * n_cell.y - 1;
                gf_r_tmp_.end_glb_[3] = 2 * n_cell.z - 1;
                gf_r_tmp_.calcSizeGlb();
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    gf_r_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];
                    gf_r_tmp_.start_loc_[1] = 0;
                    gf_r_tmp_.start_loc_[2] = 0;
                    gf_r_tmp_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                    gf_r_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                    gf_r_tmp_.end_loc_[1] = 2 * n_cell.x - 1;
                    gf_r_tmp_.end_loc_[2] = 2 * n_cell.y - 1;
                    gf_r_tmp_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                    gf_r_tmp_.calcSizeLoc();
                }
                // Set the size information of gf_k_tmp_
                if (use_mpifft_) {
                    // In this case, we set the size information
                    // basend on FFTW_MPI_TRANSPOSED_OUT format.
                    gf_k_tmp_.start_glb_[0] = 0;
                    gf_k_tmp_.start_glb_[1] = 0;
                    gf_k_tmp_.start_glb_[2] = 0;
                    gf_k_tmp_.start_glb_[3] = 0;
                    gf_k_tmp_.end_glb_[0] = n_mm_compo - 1;
                    gf_k_tmp_.end_glb_[1] = 1 + n_cell.x - 1;
                    gf_k_tmp_.end_glb_[2] = 2 * n_cell.z - 1;
                    gf_k_tmp_.end_glb_[3] = 2 * n_cell.y - 1;
                    gf_k_tmp_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        gf_k_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];  // only one (l,m)
                        gf_k_tmp_.start_loc_[1] = 0;
                        gf_k_tmp_.start_loc_[2] = 0;
                        gf_k_tmp_.start_loc_[3] = local_1_start_[rank_in_my_group_];
                        gf_k_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        gf_k_tmp_.end_loc_[1] = n_cell.x;
                        gf_k_tmp_.end_loc_[2] = 2 * n_cell.z - 1;
                        gf_k_tmp_.end_loc_[3] = local_1_end_[rank_in_my_group_];
                        gf_k_tmp_.calcSizeLoc();
                    }
                } else {
                    // In this case, we set the size information
                    // in the same way as gf_r_tmp_.
                    gf_k_tmp_.start_glb_[0] = 0;
                    gf_k_tmp_.start_glb_[1] = 0;
                    gf_k_tmp_.start_glb_[2] = 0;
                    gf_k_tmp_.start_glb_[3] = 0;
                    gf_k_tmp_.end_glb_[0] = n_mm_compo - 1;
                    gf_k_tmp_.end_glb_[1] = 1 + n_cell.x - 1;
                    gf_k_tmp_.end_glb_[2] = 2 * n_cell.y - 1;
                    gf_k_tmp_.end_glb_[3] = 2 * n_cell.z - 1;
                    gf_k_tmp_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        gf_k_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];
                        gf_k_tmp_.start_loc_[1] = 0;
                        gf_k_tmp_.start_loc_[2] = 0;
                        gf_k_tmp_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                        gf_k_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        gf_k_tmp_.end_loc_[1] = n_cell.x;
                        gf_k_tmp_.end_loc_[2] = 2 * n_cell.y - 1;
                        gf_k_tmp_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                        gf_k_tmp_.calcSizeLoc();
                    }
                }
            } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                // Set the size information of gf_r
                gf_r_tmp_.start_glb_[0] = 0;
                gf_r_tmp_.start_glb_[1] = 0;
                gf_r_tmp_.start_glb_[2] = 0;
                gf_r_tmp_.start_glb_[3] = 0;
                gf_r_tmp_.end_glb_[0] = n_mm_compo - 1;
                gf_r_tmp_.end_glb_[1] = n_cell.x - 1;
                gf_r_tmp_.end_glb_[2] = n_cell.y - 1;
                gf_r_tmp_.end_glb_[3] = n_cell.z - 1;
                gf_r_tmp_.calcSizeGlb();
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    gf_r_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];
                    gf_r_tmp_.start_loc_[1] = 0;
                    gf_r_tmp_.start_loc_[2] = 0;
                    gf_r_tmp_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                    gf_r_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                    gf_r_tmp_.end_loc_[1] = n_cell.x - 1;
                    gf_r_tmp_.end_loc_[2] = n_cell.y - 1;
                    gf_r_tmp_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                    gf_r_tmp_.calcSizeLoc();
                }
                // Set the size information of gf_k_tmp_
                if (use_mpifft_) {
                    // In this case, we set the size information
                    // basend on FFTW_MPI_TRANSPOSED_OUT format.
                    gf_k_tmp_.start_glb_[0] = 0;
                    gf_k_tmp_.start_glb_[1] = 0;
                    gf_k_tmp_.start_glb_[2] = 0;
                    gf_k_tmp_.start_glb_[3] = 0;
                    gf_k_tmp_.end_glb_[0] = n_mm_compo - 1;
                    gf_k_tmp_.end_glb_[1] = 1 + n_cell.x / 2 - 1;
                    gf_k_tmp_.end_glb_[2] = n_cell.z - 1;
                    gf_k_tmp_.end_glb_[3] = n_cell.y - 1;
                    gf_k_tmp_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        gf_k_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];  // only one (l,m)
                        gf_k_tmp_.start_loc_[1] = 0;
                        gf_k_tmp_.start_loc_[2] = 0;
                        gf_k_tmp_.start_loc_[3] = local_1_start_[rank_in_my_group_];
                        gf_k_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        gf_k_tmp_.end_loc_[1] = n_cell.x / 2;
                        gf_k_tmp_.end_loc_[2] = n_cell.z - 1;
                        gf_k_tmp_.end_loc_[3] = local_1_end_[rank_in_my_group_];
                        gf_k_tmp_.calcSizeLoc();
                    }
                } else {
                    // In this case, we set the size information
                    // in the same way as gf_r_tmp_.
                    gf_k_tmp_.start_glb_[0] = 0;
                    gf_k_tmp_.start_glb_[1] = 0;
                    gf_k_tmp_.start_glb_[2] = 0;
                    gf_k_tmp_.start_glb_[3] = 0;
                    gf_k_tmp_.end_glb_[0] = n_mm_compo - 1;
                    gf_k_tmp_.end_glb_[1] = 1 + n_cell.x / 2 - 1;
                    gf_k_tmp_.end_glb_[2] = n_cell.y - 1;
                    gf_k_tmp_.end_glb_[3] = n_cell.z - 1;
                    gf_k_tmp_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        gf_k_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];
                        gf_k_tmp_.start_loc_[1] = 0;
                        gf_k_tmp_.start_loc_[2] = 0;
                        gf_k_tmp_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                        gf_k_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        gf_k_tmp_.end_loc_[1] = n_cell.x / 2;
                        gf_k_tmp_.end_loc_[2] = n_cell.y - 1;
                        gf_k_tmp_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                        gf_k_tmp_.calcSizeLoc();
                    }
                }
            }
            // Calculate total sizes
            gf_r_tmp_.calcSizeTot();
            gf_k_tmp_.calcSizeTot();

            // Set local_start_gf_r_calc_ and local_end_gf_r_calc_
            n_proc_gf_r_calc_[0] = 1;
            n_proc_gf_r_calc_[1] = 1;
            n_proc_gf_r_calc_[2] = 1;
            if (n_proc_in_parent_group_ > 1) {
                S32 tmp = n_proc_in_parent_group_, ret;
                ret = getMaximumFactorLessThanOrEqualTo(tmp, gf_r_tmp_.size_glb_[3]);  // z
                if (ret != -1) {
                    n_proc_gf_r_calc_[2] = ret;
                    tmp /= n_proc_gf_r_calc_[2];
                }
                if (tmp > 1) {
                    ret = getMaximumFactorLessThanOrEqualTo(tmp, gf_r_tmp_.size_glb_[2]);  // y
                    if (ret != -1) {
                        n_proc_gf_r_calc_[1] = ret;
                        tmp /= n_proc_gf_r_calc_[1];
                    }
                }
                if (tmp > 1) {
                    ret = getMaximumFactorLessThanOrEqualTo(tmp, gf_r_tmp_.size_glb_[1]);  // x
                    if (ret != -1) {
                        n_proc_gf_r_calc_[0] = ret;
                        tmp /= n_proc_gf_r_calc_[0];
                    }
                }
            }
            n_proc_tot_gf_r_calc_ = n_proc_gf_r_calc_[0] * n_proc_gf_r_calc_[1] * n_proc_gf_r_calc_[2];

#if 0
                    // Check
                    if (rank_in_parent_group_ == 0) {
                        std::cout << "n_proc_tot_gf_r_calc_ = " << n_proc_tot_gf_r_calc_ << std::endl;
                        for (S32 i = 0; i < 3; i++)
                            std::cout << "(i, n_proc_gf_r_calc_[i]) = "
                                      << i << ", " << n_proc_gf_r_calc_[i]
                                      << std::endl; 
                    }
#endif
            for (S32 i = 0; i < 3; i++) {
                local_start_gf_r_calc_[i].resize(n_proc_tot_gf_r_calc_);
                local_end_gf_r_calc_[i].resize(n_proc_tot_gf_r_calc_);
            }
            for (S32 i = 0; i < n_proc_tot_gf_r_calc_; i++) {
                // Calculate idx
                S32 tmp = i, idx[3];
                idx[2] = tmp / (n_proc_gf_r_calc_[0] * n_proc_gf_r_calc_[1]);
                tmp -= idx[2] * (n_proc_gf_r_calc_[0] * n_proc_gf_r_calc_[1]);
                idx[1] = tmp / n_proc_gf_r_calc_[0];
                tmp -= idx[1] * n_proc_gf_r_calc_[0];
                idx[0] = tmp;
                // Calculate local_start_gf_r_calc_[i] and local_end_gf_r_calc_[i];
                for (S32 k = 0; k < 3; k++) {
                    const S32 size = gf_r_tmp_.size_glb_[k + 1];
                    const S32 n_slab_min = size / n_proc_gf_r_calc_[k];
                    const S32 n_rem = size - n_slab_min * n_proc_gf_r_calc_[k];
                    S32 start = gf_r_tmp_.start_glb_[k + 1];
                    for (S32 j = 0; j < idx[k]; j++) {
                        S32 n_slab = n_slab_min;
                        if (j < n_rem) n_slab++;
                        start += n_slab;
                    }
                    S32 n_slab = n_slab_min;
                    if (idx[k] < n_rem) n_slab++;
                    local_start_gf_r_calc_[k][i] = start;
                    local_end_gf_r_calc_[k][i] = start + (n_slab - 1);
                }
            }
#if 0
                    // Check
                    if (rank_in_parent_group_ == 0) {
                        const std::string file_name = "dinfo_gf_r_calc_.txt";
                        std::ofstream output_file;
                        output_file.open(file_name.c_str(), std::ios::trunc);
                        for (S32 i = 0; i < n_proc_tot_gf_r_calc_; i++) {
                            output_file << "----------------" << std::endl;
                            output_file << "Rank = " << i << std::endl;
                            for (S32 k = 0; k < 3; k++) {
                                output_file << local_start_gf_r_calc_[k][i] << "    "
                                            << local_end_gf_r_calc_[k][i] << std::endl;
                            }
                        }
                        output_file.close();
                    }
                    //Finalize();
                    //std::exit(0);
#endif
            // Set the size information of gf_r_tmp_part_
            gf_r_tmp_part_.copySizeInfoGlbOnlyFrom(gf_r_tmp_);
            if (rank_in_parent_group_ < n_proc_tot_gf_r_calc_) {
                gf_r_tmp_part_.start_loc_[0] = 0;
                gf_r_tmp_part_.start_loc_[1] = local_start_gf_r_calc_[0][rank_in_parent_group_];
                gf_r_tmp_part_.start_loc_[2] = local_start_gf_r_calc_[1][rank_in_parent_group_];
                gf_r_tmp_part_.start_loc_[3] = local_start_gf_r_calc_[2][rank_in_parent_group_];
                gf_r_tmp_part_.end_loc_[0] = n_mm_compo - 1;  // all (l,m)
                gf_r_tmp_part_.end_loc_[1] = local_end_gf_r_calc_[0][rank_in_parent_group_];
                gf_r_tmp_part_.end_loc_[2] = local_end_gf_r_calc_[1][rank_in_parent_group_];
                gf_r_tmp_part_.end_loc_[3] = local_end_gf_r_calc_[2][rank_in_parent_group_];
                gf_r_tmp_part_.calcSizeLoc();
                gf_r_tmp_part_.calcSizeTot();
            }
        } else {
            // PS::DebugUtils::fout() << "GFC initialize skipped." << std::endl;
        }
        wtime_tmp = GetWtime();
        wtime_end_initialize__main_ = wtime_tmp;
        wtime_end_initialize_ = wtime_tmp;
        time_profile_.PMMM__GFC__initialize += (wtime_end_initialize_ - wtime_start_initialize_);
    }

#else
    void initialize(const Parameters &param, const bool use_mpifft_if_possible, const bool use_ompfft, const CommInfo &comm_info, const S32 n_thread,
                    const bool debug_flag = false) {
        F64 wtime_tmp = GetWtime();
        wtime_start_initialize_ = wtime_tmp;
        // Copy parameters, etc.
        wtime_start_initialize__copy_ = wtime_tmp;
        param_prev_ = param_;  // save the previous parameters
        param_ = param;
        use_ompfft_ = use_ompfft;
        comm_info_ = comm_info;
        n_thread_prev_ = n_thread_;
        n_thread_ = n_thread;
        wtime_tmp = GetWtime();
        wtime_end_initialize__copy_ = wtime_tmp;
        // Initialize
        wtime_start_initialize__main_ = wtime_tmp;
        if ((param_ != param_prev_) || (n_thread_ != n_thread_prev_) || (use_mpifft_ != use_mpifft_if_possible)) {
            // [TODO] There is the same problem as member function `initialize`
            //        of M2L_Engine class. Refer to the comment described in
            //        M2L_Engine::initialize() for more details.
            // PS::DebugUtils::fout() << "GFC initialize performed." << std::endl;
            const S32 n_mm_compo = (2 * param_.p + 1) * (2 * param_.p + 1);
            const S32vec n_cell = param_.n_cell;
            const S32 bc = param_.bc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            // Set parent_comm
            parent_comm_ = comm_info_.getCommunicator();
            // Get information of parent_comm
            S32 n_proc, my_rank;
            MPI_Comm_size(parent_comm_, &n_proc);
            MPI_Comm_rank(parent_comm_, &my_rank);
            MPI_Comm_group(parent_comm_, &parent_group_);
            n_proc_in_parent_group_ = n_proc;
            rank_in_parent_group_ = my_rank;
            // Get MPI_Datatype for ptrdiff_t, real_t, cplx_t
            mpi_real_t_ = GetDataType<real_t>();
            mpi_cplx_t_ = GetDataType<cplx_t>();
            // Divide parent_comm into groups
            if (n_proc >= n_mm_compo) {
                n_group_ = n_mm_compo;
                // Set rank_start_, rank_end_
                n_proc_min_in_group_ = n_proc / n_group_;
                S32 n_rem = n_proc - n_proc_min_in_group_ * n_group_;
                rank_start_.resize(n_group_);
                rank_end_.resize(n_group_);
                rank_start_[0] = 0;
                for (S32 i = 0; i < n_group_; i++) {
                    S32 n_proc_in_group = n_proc_min_in_group_;
                    if (i < n_rem) n_proc_in_group++;
                    rank_end_[i] = rank_start_[i] + (n_proc_in_group - 1);
                    if (i < n_group_ - 1) rank_start_[i + 1] = rank_end_[i] + 1;
                }
                if (debug_flag) {
                    if (my_rank == 0) {
                        std::cout << "### Information of groups for G2G ###" << std::endl;
                        std::cout << "(n_proc = " << n_proc << ")" << std::endl;
                        for (S32 i = 0; i < n_group_; i++) {
                            std::cout << "i = " << i << " rank_start = " << rank_start_[i] << " rank_end = " << rank_end_[i] << std::endl;
                        }
                    }
                }
                // Calculate idx_to_my_group, n_proc_in_my_group, rank_in_my_group
                for (S32 i = 0; i < n_group_; i++)
                    if ((rank_start_[i] <= my_rank) && (my_rank <= rank_end_[i])) {
                        idx_to_my_group_ = i;
                        n_proc_in_my_group_ = rank_end_[i] - rank_start_[i] + 1;
                        rank_in_my_group_ = my_rank - rank_start_[i];
                    }
                // Set lm_start_, lm_end_
                lm_start_.resize(n_group_);
                lm_end_.resize(n_group_);
                for (S32 i = 0; i < n_group_; i++) {
                    lm_start_[i] = i;
                    lm_end_[i] = i;
                    // Each group is responsible for a single (l,m)
                }
                // Make group_all & comm_all
                if (group_all_ != MPI_GROUP_NULL) MPI_Group_free(&group_all_);
                if (comm_all_ != MPI_COMM_NULL) MPI_Comm_free(&comm_all_);
                std::vector<S32> ranks;
                S32 rnk_start = rank_start_[idx_to_my_group_];
                S32 rnk_end = rank_end_[idx_to_my_group_];
                for (S32 rnk = rnk_start; rnk <= rnk_end; rnk++) ranks.push_back(rnk);
                MPI_Group_incl(parent_group_, n_proc_in_my_group_, &ranks[0], &group_all_);
                MPI_Comm_create(parent_comm_, group_all_, &comm_all_);
                // Calculate n_proc_for_fft
                S32 n_proc_max_for_fft;
                if (bc == BOUNDARY_CONDITION_OPEN)
                    n_proc_max_for_fft = n_cell.z;
                else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
                    n_proc_max_for_fft = n_cell.z / 2;
                else
                    assert(false);
                if (n_proc_min_in_group_ <= n_proc_max_for_fft)
                    n_proc_for_fft_ = n_proc_min_in_group_;
                else
                    n_proc_for_fft_ = n_proc_max_for_fft;
                // Make group_fft, comm_fft
                if (n_proc_for_fft_ > 1)
                    is_mpifft_usable_ = true;
                else
                    is_mpifft_usable_ = false;
                if (is_mpifft_usable_)
                    use_mpifft_ = use_mpifft_if_possible;  // specified by user.
                else
                    use_mpifft_ = false;
                if (use_mpifft_) {
                    if (!FFTW_Status::is_fftw_mpi_init_called()) {
                        fftw_mpi_init();
                        FFTW_Status::is_fftw_mpi_init_called() = true;
                    }
                    if (group_fft_ != MPI_GROUP_NULL) MPI_Group_free(&group_fft_);
                    if (comm_fft_ != MPI_COMM_NULL) MPI_Comm_free(&comm_fft_);
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        ranks.resize(n_proc_for_fft_);
                        for (S32 i = 0; i < n_proc_for_fft_; i++) ranks[i] = rank_start_[idx_to_my_group_] + i;
                        MPI_Group_incl(parent_group_, n_proc_for_fft_, &ranks[0], &group_fft_);
                        MPI_Comm_create(parent_comm_, group_fft_, &comm_fft_);
                    } else {
                        group_fft_ = MPI_GROUP_NULL;
                        MPI_Comm_create(parent_comm_, MPI_GROUP_EMPTY, &comm_fft_);
                    }
                    // In this case, we need to calculate local_0_start & local_0_end
                    // using FFTW API.
                    // Resize local_?_start, local_?_end
                    local_0_start_.resize(n_proc_for_fft_);
                    local_0_end_.resize(n_proc_for_fft_);
                    local_1_start_.resize(n_proc_for_fft_);
                    local_1_end_.resize(n_proc_for_fft_);
                    // Calculate loca_?_start, local_?_end
                    if (comm_fft_ != MPI_COMM_NULL) {
                        ptrdiff_t alloc_local, local_n0, local_0_start;
                        ptrdiff_t local_n1, local_1_start;
                        if (bc == BOUNDARY_CONDITION_OPEN) {
                            alloc_local = fftw_mpi_local_size_3d_transposed(2 * n_cell.z, 2 * n_cell.y, n_cell.x + 1, comm_fft_, &local_n0,
                                                                            &local_0_start, &local_n1, &local_1_start);
                        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                            alloc_local = fftw_mpi_local_size_3d_transposed(n_cell.z, n_cell.y, n_cell.x / 2 + 1, comm_fft_, &local_n0,
                                                                            &local_0_start, &local_n1, &local_1_start);
                        } else {
                            if (my_rank == 0) std::cout << "This boundary condition is not supported yet." << std::endl;
                            Abort(-1);
                        }
                        S32 sendbuf[4];
                        sendbuf[0] = local_0_start;
                        sendbuf[1] = local_n0;
                        sendbuf[2] = local_1_start;
                        sendbuf[3] = local_n1;
                        S32 *recvbuf = new S32[4 * n_proc_for_fft_];
                        MPI_Gather(sendbuf, 4, MPI_INT, recvbuf, 4, MPI_INT, 0, comm_fft_);
                        if (rank_in_my_group_ == 0) {
                            for (S32 i = 0; i < n_proc_for_fft_; i++) {
                                local_0_start_[i] = recvbuf[4 * i];
                                local_0_end_[i] = local_0_start_[i] + recvbuf[4 * i + 1] - 1;
                                local_1_start_[i] = recvbuf[4 * i + 2];
                                local_1_end_[i] = local_1_start_[i] + recvbuf[4 * i + 3] - 1;
                            }
                        }
                        delete[] recvbuf;
                    }
                    MPI_Bcast(&local_0_start_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                    MPI_Bcast(&local_0_end_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                    MPI_Bcast(&local_1_start_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                    MPI_Bcast(&local_1_end_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                    // Check
                    if (debug_flag) {
                        if (rank_in_parent_group_ == 0) {
                            std::cout << "### Information of decomposition of FFTW ###" << std::endl;
                            std::cout << "(n_cell.z = " << n_cell.z << ")" << std::endl;
                            for (S32 i = 0; i < n_proc_for_fft_; i++) {
                                std::cout << "i = " << i << " local_0_start = " << local_0_start_[i] << " local_0_end = " << local_0_end_[i]
                                          << " local_1_start = " << local_1_start_[i] << " local_1_end = " << local_1_end_[i] << std::endl;
                            }
                        }
                    }
                } else {
                    n_proc_for_fft_ = 1;  // reset
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                    if (!FFTW_Status::is_fftw_init_threads_called()) {
                        fftw_init_threads();
                        FFTW_Status::is_fftw_init_threads_called() = true;
                    }
#endif
                    group_fft_ = MPI_GROUP_NULL;
                    comm_fft_ = MPI_COMM_NULL;
                    local_0_start_.resize(1);
                    local_0_end_.resize(1);
                    if (bc == BOUNDARY_CONDITION_OPEN) {
                        local_0_start_[0] = 0;
                        local_0_end_[0] = 2 * n_cell.z - 1;
                    } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                        local_0_start_[0] = 0;
                        local_0_end_[0] = n_cell.z - 1;
                    }
                }
                // Make group_int, comm_int
                if (group_int_ != MPI_GROUP_NULL) MPI_Group_free(&group_int_);
                if (comm_int_ != MPI_COMM_NULL) MPI_Comm_free(&comm_int_);
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    ranks.resize(n_group_);
                    for (S32 i = 0; i < n_group_; i++) ranks[i] = rank_start_[i] + rank_in_my_group_;
                    MPI_Group_incl(parent_group_, n_group_, &ranks[0], &group_int_);
                    MPI_Comm_create(parent_comm_, group_int_, &comm_int_);
                } else {
                    group_int_ = MPI_GROUP_NULL;
                    MPI_Comm_create(parent_comm_, MPI_GROUP_EMPTY, &comm_int_);
                }
            } else {
                n_group_ = n_proc;
                n_proc_in_my_group_ = 1;
                rank_in_my_group_ = 0;
                n_proc_min_in_group_ = 1;
                n_proc_for_fft_ = 1;
                idx_to_my_group_ = my_rank;
                is_mpifft_usable_ = false;
                use_mpifft_ = false;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                if (!FFTW_Status::is_fftw_init_threads_called()) {
                    fftw_init_threads();
                    FFTW_Status::is_fftw_init_threads_called() = true;
                }
#endif
                // Set rank_start_, rank_end_
                rank_start_.resize(n_group_);
                rank_end_.resize(n_group_);
                for (S32 i = 0; i < n_group_; i++) {
                    rank_start_[i] = i;
                    rank_end_[i] = i;
                }
                // Set lm_start_, lm_end_
                const S32 n_mm_compo_per_proc_min = n_mm_compo / n_proc;
                S32 n_rem = n_mm_compo - n_mm_compo_per_proc_min * n_proc;
                lm_start_.resize(n_proc);
                lm_end_.resize(n_proc);
                lm_start_[0] = 0;
                for (S32 i = 0; i < n_proc; i++) {
                    S32 n_mm_compo_per_proc = n_mm_compo_per_proc_min;
                    if (i < n_rem) n_mm_compo_per_proc++;
                    lm_end_[i] = lm_start_[i] + (n_mm_compo_per_proc - 1);
                    if (i < n_proc - 1) lm_start_[i + 1] = lm_end_[i] + 1;
                }
                if (debug_flag) {
                    if (rank_in_parent_group_ == 0) {
                        std::cout << "### Information of decomposition of (l,m) ###" << std::endl;
                        std::cout << "(n_mm_compo = " << n_mm_compo << ")" << std::endl;
                        for (S32 i = 0; i < n_proc; i++) {
                            std::cout << "i = " << i << " local_0_start = " << lm_start_[i] << " local_0_end = " << lm_end_[i] << std::endl;
                        }
                    }
                }
                // Set local_0_start_, local_0_end_
                group_fft_ = MPI_GROUP_NULL;
                comm_fft_ = MPI_COMM_NULL;
                local_0_start_.resize(n_proc_for_fft_);
                local_0_end_.resize(n_proc_for_fft_);
                if (bc == BOUNDARY_CONDITION_OPEN) {
                    local_0_start_[0] = 0;
                    local_0_end_[0] = 2 * n_cell.z - 1;
                } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                    local_0_start_[0] = 0;
                    local_0_end_[0] = n_cell.z - 1;
                }
                // Note that we do not need to set local_1_start_
                // and local_1_end_ because of use_mpifft = false.
                // (hence, both arrays are not accessed)

                // Make group_int, comm_int
                group_int_ = parent_group_;
                comm_int_ = parent_comm_;
            }
#else  // PARTICLE_SIMULATOR_MPI_PARALLEL
       // Sequential execution/ OpenMP parallelization
            n_proc_in_parent_group_ = 1;
            rank_in_parent_group_ = 0;
            n_group_ = 1;
            n_proc_in_my_group_ = 1;
            rank_in_my_group_ = 0;
            n_proc_min_in_group_ = 1;
            n_proc_for_fft_ = 1;
            idx_to_my_group_ = 0;
            is_mpifft_usable_ = false;
            use_mpifft_ = false;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            if (!FFTW_Status::is_fftw_init_threads_called()) {
                fftw_init_threads();
                FFTW_Status::is_fftw_init_threads_called() = true;
            }
#endif
            // Set rank_start_(n_group_), rank_end_(n_group_) [idx_to_my_group_]
            rank_start_.resize(n_group_);
            rank_end_.resize(n_group_);
            rank_start_[0] = 0;
            rank_end_[0] = 0;
            // Set lm_start_(n_group_), lm_end_(n_group_) [idx_to_my_group_]
            lm_start_.resize(n_group_);
            lm_end_.resize(n_group_);
            lm_start_[0] = 0;
            lm_end_[0] = n_mm_compo - 1;
            // Set local_0_start__, local_0_end_
            local_0_start_.resize(n_proc_for_fft_);
            local_0_end_.resize(n_proc_for_fft_);
            if (bc == BOUNDARY_CONDITION_OPEN) {
                local_0_start_[0] = 0;
                local_0_end_[0] = 2 * n_cell.z - 1;
            } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                local_0_start_[0] = 0;
                local_0_end_[0] = n_cell.z - 1;
            }
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL

            // Set the size information of gf_r_tmp_, gf_k_tmp_
            if (bc == BOUNDARY_CONDITION_OPEN) {
                // Set the size information of gf_r_tmp_
                gf_r_tmp_.start_glb_[0] = 0;
                gf_r_tmp_.start_glb_[1] = 0;
                gf_r_tmp_.start_glb_[2] = 0;
                gf_r_tmp_.start_glb_[3] = 0;
                gf_r_tmp_.end_glb_[0] = n_mm_compo - 1;
                gf_r_tmp_.end_glb_[1] = 2 * n_cell.x - 1;
                gf_r_tmp_.end_glb_[2] = 2 * n_cell.y - 1;
                gf_r_tmp_.end_glb_[3] = 2 * n_cell.z - 1;
                gf_r_tmp_.calcSizeGlb();
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    gf_r_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];
                    gf_r_tmp_.start_loc_[1] = 0;
                    gf_r_tmp_.start_loc_[2] = 0;
                    gf_r_tmp_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                    gf_r_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                    gf_r_tmp_.end_loc_[1] = 2 * n_cell.x - 1;
                    gf_r_tmp_.end_loc_[2] = 2 * n_cell.y - 1;
                    gf_r_tmp_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                    gf_r_tmp_.calcSizeLoc();
                }
                // Set the size information of gf_k_tmp_
                if (use_mpifft_) {
                    // In this case, we set the size information
                    // basend on FFTW_MPI_TRANSPOSED_OUT format.
                    gf_k_tmp_.start_glb_[0] = 0;
                    gf_k_tmp_.start_glb_[1] = 0;
                    gf_k_tmp_.start_glb_[2] = 0;
                    gf_k_tmp_.start_glb_[3] = 0;
                    gf_k_tmp_.end_glb_[0] = n_mm_compo - 1;
                    gf_k_tmp_.end_glb_[1] = 1 + n_cell.x - 1;
                    gf_k_tmp_.end_glb_[2] = 2 * n_cell.z - 1;
                    gf_k_tmp_.end_glb_[3] = 2 * n_cell.y - 1;
                    gf_k_tmp_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        gf_k_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];  // only one (l,m)
                        gf_k_tmp_.start_loc_[1] = 0;
                        gf_k_tmp_.start_loc_[2] = 0;
                        gf_k_tmp_.start_loc_[3] = local_1_start_[rank_in_my_group_];
                        gf_k_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        gf_k_tmp_.end_loc_[1] = n_cell.x;
                        gf_k_tmp_.end_loc_[2] = 2 * n_cell.z - 1;
                        gf_k_tmp_.end_loc_[3] = local_1_end_[rank_in_my_group_];
                        gf_k_tmp_.calcSizeLoc();
                    }
                } else {
                    // In this case, we set the size information
                    // in the same way as gf_r_tmp_.
                    gf_k_tmp_.start_glb_[0] = 0;
                    gf_k_tmp_.start_glb_[1] = 0;
                    gf_k_tmp_.start_glb_[2] = 0;
                    gf_k_tmp_.start_glb_[3] = 0;
                    gf_k_tmp_.end_glb_[0] = n_mm_compo - 1;
                    gf_k_tmp_.end_glb_[1] = 1 + n_cell.x - 1;
                    gf_k_tmp_.end_glb_[2] = 2 * n_cell.y - 1;
                    gf_k_tmp_.end_glb_[3] = 2 * n_cell.z - 1;
                    gf_k_tmp_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        gf_k_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];
                        gf_k_tmp_.start_loc_[1] = 0;
                        gf_k_tmp_.start_loc_[2] = 0;
                        gf_k_tmp_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                        gf_k_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        gf_k_tmp_.end_loc_[1] = n_cell.x;
                        gf_k_tmp_.end_loc_[2] = 2 * n_cell.y - 1;
                        gf_k_tmp_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                        gf_k_tmp_.calcSizeLoc();
                    }
                }
            } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                // Set the size information of gf_r
                gf_r_tmp_.start_glb_[0] = 0;
                gf_r_tmp_.start_glb_[1] = 0;
                gf_r_tmp_.start_glb_[2] = 0;
                gf_r_tmp_.start_glb_[3] = 0;
                gf_r_tmp_.end_glb_[0] = n_mm_compo - 1;
                gf_r_tmp_.end_glb_[1] = n_cell.x - 1;
                gf_r_tmp_.end_glb_[2] = n_cell.y - 1;
                gf_r_tmp_.end_glb_[3] = n_cell.z - 1;
                gf_r_tmp_.calcSizeGlb();
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    gf_r_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];
                    gf_r_tmp_.start_loc_[1] = 0;
                    gf_r_tmp_.start_loc_[2] = 0;
                    gf_r_tmp_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                    gf_r_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                    gf_r_tmp_.end_loc_[1] = n_cell.x - 1;
                    gf_r_tmp_.end_loc_[2] = n_cell.y - 1;
                    gf_r_tmp_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                    gf_r_tmp_.calcSizeLoc();
                }
                // Set the size information of gf_k_tmp_
                if (use_mpifft_) {
                    // In this case, we set the size information
                    // basend on FFTW_MPI_TRANSPOSED_OUT format.
                    gf_k_tmp_.start_glb_[0] = 0;
                    gf_k_tmp_.start_glb_[1] = 0;
                    gf_k_tmp_.start_glb_[2] = 0;
                    gf_k_tmp_.start_glb_[3] = 0;
                    gf_k_tmp_.end_glb_[0] = n_mm_compo - 1;
                    gf_k_tmp_.end_glb_[1] = 1 + n_cell.x / 2 - 1;
                    gf_k_tmp_.end_glb_[2] = n_cell.z - 1;
                    gf_k_tmp_.end_glb_[3] = n_cell.y - 1;
                    gf_k_tmp_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        gf_k_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];  // only one (l,m)
                        gf_k_tmp_.start_loc_[1] = 0;
                        gf_k_tmp_.start_loc_[2] = 0;
                        gf_k_tmp_.start_loc_[3] = local_1_start_[rank_in_my_group_];
                        gf_k_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        gf_k_tmp_.end_loc_[1] = n_cell.x / 2;
                        gf_k_tmp_.end_loc_[2] = n_cell.z - 1;
                        gf_k_tmp_.end_loc_[3] = local_1_end_[rank_in_my_group_];
                        gf_k_tmp_.calcSizeLoc();
                    }
                } else {
                    // In this case, we set the size information
                    // in the same way as gf_r_tmp_.
                    gf_k_tmp_.start_glb_[0] = 0;
                    gf_k_tmp_.start_glb_[1] = 0;
                    gf_k_tmp_.start_glb_[2] = 0;
                    gf_k_tmp_.start_glb_[3] = 0;
                    gf_k_tmp_.end_glb_[0] = n_mm_compo - 1;
                    gf_k_tmp_.end_glb_[1] = 1 + n_cell.x / 2 - 1;
                    gf_k_tmp_.end_glb_[2] = n_cell.y - 1;
                    gf_k_tmp_.end_glb_[3] = n_cell.z - 1;
                    gf_k_tmp_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        gf_k_tmp_.start_loc_[0] = lm_start_[idx_to_my_group_];
                        gf_k_tmp_.start_loc_[1] = 0;
                        gf_k_tmp_.start_loc_[2] = 0;
                        gf_k_tmp_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                        gf_k_tmp_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        gf_k_tmp_.end_loc_[1] = n_cell.x / 2;
                        gf_k_tmp_.end_loc_[2] = n_cell.y - 1;
                        gf_k_tmp_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                        gf_k_tmp_.calcSizeLoc();
                    }
                }
            }
            // Calculate total sizes
            gf_r_tmp_.calcSizeTot();
            gf_k_tmp_.calcSizeTot();

            // Set local_start_gf_r_calc_ and local_end_gf_r_calc_
            n_proc_gf_r_calc_[0] = 1;
            n_proc_gf_r_calc_[1] = 1;
            n_proc_gf_r_calc_[2] = 1;
            if (n_proc_in_parent_group_ > 1) {
                S32 tmp = n_proc_in_parent_group_, ret;
                ret = getMaximumFactorLessThanOrEqualTo(tmp, gf_r_tmp_.size_glb_[3]);  // z
                if (ret != -1) {
                    n_proc_gf_r_calc_[2] = ret;
                    tmp /= n_proc_gf_r_calc_[2];
                }
                if (tmp > 1) {
                    ret = getMaximumFactorLessThanOrEqualTo(tmp, gf_r_tmp_.size_glb_[2]);  // y
                    if (ret != -1) {
                        n_proc_gf_r_calc_[1] = ret;
                        tmp /= n_proc_gf_r_calc_[1];
                    }
                }
                if (tmp > 1) {
                    ret = getMaximumFactorLessThanOrEqualTo(tmp, gf_r_tmp_.size_glb_[1]);  // x
                    if (ret != -1) {
                        n_proc_gf_r_calc_[0] = ret;
                        tmp /= n_proc_gf_r_calc_[0];
                    }
                }
            }
            n_proc_tot_gf_r_calc_ = n_proc_gf_r_calc_[0] * n_proc_gf_r_calc_[1] * n_proc_gf_r_calc_[2];

#if 0
                    // Check
                    if (rank_in_parent_group_ == 0) {
                        std::cout << "n_proc_tot_gf_r_calc_ = " << n_proc_tot_gf_r_calc_ << std::endl;
                        for (S32 i = 0; i < 3; i++)
                            std::cout << "(i, n_proc_gf_r_calc_[i]) = "
                                      << i << ", " << n_proc_gf_r_calc_[i]
                                      << std::endl; 
                    }
#endif
            for (S32 i = 0; i < 3; i++) {
                local_start_gf_r_calc_[i].resize(n_proc_tot_gf_r_calc_);
                local_end_gf_r_calc_[i].resize(n_proc_tot_gf_r_calc_);
            }
            for (S32 i = 0; i < n_proc_tot_gf_r_calc_; i++) {
                // Calculate idx
                S32 tmp = i, idx[3];
                idx[2] = tmp / (n_proc_gf_r_calc_[0] * n_proc_gf_r_calc_[1]);
                tmp -= idx[2] * (n_proc_gf_r_calc_[0] * n_proc_gf_r_calc_[1]);
                idx[1] = tmp / n_proc_gf_r_calc_[0];
                tmp -= idx[1] * n_proc_gf_r_calc_[0];
                idx[0] = tmp;
                // Calculate local_start_gf_r_calc_[i] and local_end_gf_r_calc_[i];
                for (S32 k = 0; k < 3; k++) {
                    const S32 size = gf_r_tmp_.size_glb_[k + 1];
                    const S32 n_slab_min = size / n_proc_gf_r_calc_[k];
                    const S32 n_rem = size - n_slab_min * n_proc_gf_r_calc_[k];
                    S32 start = gf_r_tmp_.start_glb_[k + 1];
                    for (S32 j = 0; j < idx[k]; j++) {
                        S32 n_slab = n_slab_min;
                        if (j < n_rem) n_slab++;
                        start += n_slab;
                    }
                    S32 n_slab = n_slab_min;
                    if (idx[k] < n_rem) n_slab++;
                    local_start_gf_r_calc_[k][i] = start;
                    local_end_gf_r_calc_[k][i] = start + (n_slab - 1);
                }
            }
#if 0
                    // Check
                    if (rank_in_parent_group_ == 0) {
                        const std::string file_name = "dinfo_gf_r_calc_.txt";
                        std::ofstream output_file;
                        output_file.open(file_name.c_str(), std::ios::trunc);
                        for (S32 i = 0; i < n_proc_tot_gf_r_calc_; i++) {
                            output_file << "----------------" << std::endl;
                            output_file << "Rank = " << i << std::endl;
                            for (S32 k = 0; k < 3; k++) {
                                output_file << local_start_gf_r_calc_[k][i] << "    "
                                            << local_end_gf_r_calc_[k][i] << std::endl;
                            }
                        }
                        output_file.close();
                    }
                    //Finalize();
                    //std::exit(0);
#endif
            // Set the size information of gf_r_tmp_part_
            gf_r_tmp_part_.copySizeInfoGlbOnlyFrom(gf_r_tmp_);
            if (rank_in_parent_group_ < n_proc_tot_gf_r_calc_) {
                gf_r_tmp_part_.start_loc_[0] = 0;
                gf_r_tmp_part_.start_loc_[1] = local_start_gf_r_calc_[0][rank_in_parent_group_];
                gf_r_tmp_part_.start_loc_[2] = local_start_gf_r_calc_[1][rank_in_parent_group_];
                gf_r_tmp_part_.start_loc_[3] = local_start_gf_r_calc_[2][rank_in_parent_group_];
                gf_r_tmp_part_.end_loc_[0] = n_mm_compo - 1;  // all (l,m)
                gf_r_tmp_part_.end_loc_[1] = local_end_gf_r_calc_[0][rank_in_parent_group_];
                gf_r_tmp_part_.end_loc_[2] = local_end_gf_r_calc_[1][rank_in_parent_group_];
                gf_r_tmp_part_.end_loc_[3] = local_end_gf_r_calc_[2][rank_in_parent_group_];
                gf_r_tmp_part_.calcSizeLoc();
                gf_r_tmp_part_.calcSizeTot();
            }
        } else {
            // PS::DebugUtils::fout() << "GFC initialize skipped." << std::endl;
        }
        wtime_tmp = GetWtime();
        wtime_end_initialize__main_ = wtime_tmp;
        wtime_end_initialize_ = wtime_tmp;
        time_profile_.PMMM__GFC__initialize += (wtime_end_initialize_ - wtime_start_initialize_);
    }
#endif

    void calcGreenFunction(Buffer<cplx_t> &gf_k, const bool is_gf_k_transposed) {
        if (param_ != param_prev_) {
            // Extract parameters
            // const S32 p = param_.p;
            // const S32 bc = param_.bc;
            const F64ort pos_unit_cell = param_.pos_unit_cell;
            const S32vec n_cell = param_.n_cell;
            const F64vec width_cell = GetWidthOfParticleMeshCell(pos_unit_cell, n_cell);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok1 @ calcGreenFunction() [PMMM]" << std::endl;
#endif

            // Calculate gf_r_tmp_part_
            wtime_start_calc_gf_r_ = GetWtime();
            calcGreenFunctionInRealSpace(width_cell);
            wtime_end_calc_gf_r_ = GetWtime();
            time_profile_.PMMM__GFC__calc_gf_r += (wtime_end_calc_gf_r_ - wtime_start_calc_gf_r_);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok2 @ calcGreenFunction() [PMMM]" << std::endl;
            {
                const std::string file_prefix = "gf_r_tmp_part_";
                const S32 file_num = rank_in_parent_group_;
                if (rank_in_parent_group_ < n_proc_tot_gf_r_calc_) {
                    gf_r_tmp_part_.writeBufferToFile(file_prefix, file_num);
                }
            }
#endif

            // Collect gf_r_tmp_ necessary to perform FFT
            wtime_start_redist_gf_r_ = GetWtime();
            redistGFR();
            wtime_end_redist_gf_r_ = GetWtime();
            time_profile_.PMMM__GFC__redist_gf_r += (wtime_end_redist_gf_r_ - wtime_start_redist_gf_r_);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok3 @ calcGreenFunction() [PMMM]" << std::endl;
            {
                const std::string file_prefix = "gf_r_tmp_";
                const S32 file_num = rank_in_parent_group_;
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    gf_r_tmp_.writeBufferToFile(file_prefix, file_num);
                }
            }
#endif

            // Calculate gf_k_ using FFT
            wtime_start_calc_gf_k_ = GetWtime();
            calcGreenFunctionInWavenumberSpace();
            wtime_end_calc_gf_k_ = GetWtime();
            time_profile_.PMMM__GFC__calc_gf_k += (wtime_end_calc_gf_k_ - wtime_start_calc_gf_k_);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok4 @calcGreenFunction() [PMMM]" << std::endl;
            {
                const std::string file_prefix = "gf_k_tmp_";
                const S32 file_num = rank_in_parent_group_;
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    if (use_mpifft_) {
                        gf_k_tmp_.writeBufferToFile(file_prefix, file_num, true);
                    } else {
                        gf_k_tmp_.writeBufferToFile(file_prefix, file_num, false);
                    }
                }
            }
#endif

            // Collect gf_k_
            wtime_start_redist_gf_k_ = GetWtime();
            redistGFK(gf_k, is_gf_k_transposed);
            wtime_end_redist_gf_k_ = GetWtime();
            time_profile_.PMMM__GFC__redist_gf_k += (wtime_end_redist_gf_k_ - wtime_start_redist_gf_k_);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok_last @calcGreenFunction() [PMMM]" << std::endl;
            {
                const std::string file_prefix = "gf_k_new_";
                const S32 file_num = rank_in_parent_group_;
                if (gf_k.size_loc_tot_ > 0) {
                    if (is_gf_k_transposed) {
                        gf_k.writeBufferToFile(file_prefix, file_num, false, true);
                    } else {
                        gf_k.writeBufferToFile(file_prefix, file_num, false, false);
                    }
                }
            }
#endif
        } else {
            const F64 now = GetWtime();
            wtime_start_calc_gf_r_ = now;
            wtime_end_calc_gf_r_ = now;
            wtime_start_redist_gf_r_ = now;
            wtime_end_redist_gf_r_ = now;
            wtime_start_calc_gf_k_ = now;
            wtime_end_calc_gf_k_ = now;
            wtime_start_redist_gf_k_ = now;
            wtime_end_redist_gf_k_ = now;
        }
    }

    void clearTimeProfile() { time_profile_.clear(); }

    void clearWallTime() {
        wtime_start_initialize_ = 0;
        wtime_end_initialize_ = 0;
        wtime_start_calc_gf_r_ = 0;
        wtime_end_calc_gf_r_ = 0;
        wtime_start_redist_gf_r_ = 0;
        wtime_end_redist_gf_r_ = 0;
        wtime_start_calc_gf_k_ = 0;
        wtime_end_calc_gf_k_ = 0;
        wtime_start_redist_gf_k_ = 0;
        wtime_end_redist_gf_k_ = 0;

        wtime_start_initialize__copy_ = 0;
        wtime_end_initialize__copy_ = 0;
        wtime_start_initialize__main_ = 0;
        wtime_end_initialize__main_ = 0;
    }

    void clearCounterAll() {
        clearWallTime();
        clearTimeProfile();
    }

    TimeProfile getTimeProfile() const { return time_profile_; }

   private:
    void calcGreenFunctionInRealSpace(const F64vec cell_length) {
        gf_r_tmp_part_.resize();  // allocate
        if (param_.bc == BOUNDARY_CONDITION_OPEN) {
            calcGreenFunctionInRealSpaceForOpenBC(cell_length);
        } else if (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
            calcGreenFunctionInRealSpaceForPeriodicBC(cell_length);
        }
    }

    void calcGreenFunctionInRealSpaceForOpenBC(const F64vec cell_length) {
        if (rank_in_parent_group_ < n_proc_tot_gf_r_calc_) {
            // Extract parameters
            const S32 p = param_.p;
            const S32 icut = param_.icut;
            // const S32vec n_cell = param_.n_cell;
            //  Local variables
            Slm<real_t> *slm = new Slm<real_t>[n_thread_];
            for (S32 i = 0; i < n_thread_; i++) slm[i].alloc(2 * p);
            {
                Slm<real_t> tmp;
                tmp.make_table(2 * p);
            }
            const S32 size_rspace_loc_tot = gf_r_tmp_part_.size_loc_[1] * gf_r_tmp_part_.size_loc_[2] * gf_r_tmp_part_.size_loc_[3];
            PS_OMP(omp parallel for num_threads(n_thread_))
            for (S32 idx = 0; idx < size_rspace_loc_tot; idx++) {
                // See Notes #c2f9fbef for the reason why we used idx-loop.

                // Get thread ID
                const S32 ith = comm_info_.getThreadNum();

                // Calculate the cell index in the local region
                S32 idx_tmp = idx;
                const S32 k_loc = idx_tmp / (gf_r_tmp_part_.size_loc_[1] * gf_r_tmp_part_.size_loc_[2]);
                idx_tmp -= k_loc * (gf_r_tmp_part_.size_loc_[1] * gf_r_tmp_part_.size_loc_[2]);
                const S32 j_loc = idx_tmp / gf_r_tmp_part_.size_loc_[1];
                idx_tmp -= j_loc * gf_r_tmp_part_.size_loc_[1];
                const S32 i_loc = idx_tmp;

                // Calculate the cell index in the global region
                const S32 i = i_loc + gf_r_tmp_part_.start_loc_[1];
                const S32 j = j_loc + gf_r_tmp_part_.start_loc_[2];
                const S32 k = k_loc + gf_r_tmp_part_.start_loc_[3];

                if ((k == gf_r_tmp_part_.size_glb_[3] / 2) || (j == gf_r_tmp_part_.size_glb_[2] / 2) || (i == gf_r_tmp_part_.size_glb_[1] / 2)) {
                    for (S32 lm = gf_r_tmp_part_.start_loc_[0]; lm <= gf_r_tmp_part_.end_loc_[0]; lm++) {
                        const S32 lm_loc = lm - gf_r_tmp_part_.start_loc_[0];
                        const S32 adr = lm_loc + gf_r_tmp_part_.size_loc_[0] *
                                                     (i_loc + gf_r_tmp_part_.size_loc_[1] * (j_loc + gf_r_tmp_part_.size_loc_[2] * k_loc));
                        gf_r_tmp_part_.buf_[adr] = 0.0;
                    }
                    continue;
                }
                const S32 kk = (k > gf_r_tmp_part_.size_glb_[3] / 2) ? k - gf_r_tmp_part_.size_glb_[3] : k;
                const S32 jj = (j > gf_r_tmp_part_.size_glb_[2] / 2) ? j - gf_r_tmp_part_.size_glb_[2] : j;
                const S32 ii = (i > gf_r_tmp_part_.size_glb_[1] / 2) ? i - gf_r_tmp_part_.size_glb_[1] : i;
                // [TODO]
                //    In the calculation above, we must also use
                //    gf_r_tmp_.start_glb_[] and gf_r_tmp_.end_glb_[]
                //    if gf_r_tmp_.start_glb_[?] != 0.
                if (abs(kk) <= icut && abs(jj) <= icut && abs(ii) <= icut) {
                    for (S32 lm = gf_r_tmp_part_.start_loc_[0]; lm <= gf_r_tmp_part_.end_loc_[0]; lm++) {
                        const S32 lm_loc = lm - gf_r_tmp_part_.start_loc_[0];
                        const S32 adr = lm_loc + gf_r_tmp_part_.size_loc_[0] *
                                                     (i_loc + gf_r_tmp_part_.size_loc_[1] * (j_loc + gf_r_tmp_part_.size_loc_[2] * k_loc));
                        gf_r_tmp_part_.buf_[adr] = 0.0;
                    }
                    continue;
                }
                const F64 dx = cell_length.x * F64(ii);
                const F64 dy = cell_length.y * F64(jj);
                const F64 dz = cell_length.z * F64(kk);
                slm[ith].eval_opt(-dx, dy, dz);  // eval S_l^{-m}
                for (S32 lm = gf_r_tmp_part_.start_loc_[0]; lm <= gf_r_tmp_part_.end_loc_[0]; lm++) {
                    const S32 lm_loc = lm - gf_r_tmp_part_.start_loc_[0];
                    const S32 adr =
                        lm_loc + gf_r_tmp_part_.size_loc_[0] * (i_loc + gf_r_tmp_part_.size_loc_[1] * (j_loc + gf_r_tmp_part_.size_loc_[2] * k_loc));
                    assert(0 <= adr && adr < gf_r_tmp_part_.size_loc_tot_);
                    gf_r_tmp_part_.buf_[adr] = slm[ith].buf[lm];
                }
            }
            // Free memory
            for (S32 i = 0; i < n_thread_; i++) slm[i].freeMem();
            delete[] slm;
        }
    }

    void calcGreenFunctionInRealSpaceForPeriodicBC(const F64vec cell_length) {
        if (rank_in_parent_group_ < n_proc_tot_gf_r_calc_) {
            // Extract parameters
            const S32 p = param_.p;
            const S32 icut = param_.icut;
            const F64vec period_length = param_.pos_unit_cell.getFullLength();
            // const S32vec n_cell = param_.n_cell;
            const F64 alpha = param_.alpha;
            const S32 NMAX = param_.nmax;
            const S32 MMAX = param_.mmax;
            // Local variables
            const F64 bx = period_length.x;
            const F64 by = period_length.y;
            const F64 bz = period_length.z;
            const F64 bx_inv = 1.0 / bx;
            const F64 by_inv = 1.0 / by;
            const F64 bz_inv = 1.0 / bz;
            const F64 bx2_inv = 1.0 / (bx * bx);
            const F64 by2_inv = 1.0 / (by * by);
            const F64 bz2_inv = 1.0 / (bz * bz);
            const F64 vol_inv = 1.0 / (bx * by * bz);
            CutFunc<real_t, cplx_t> *cf;
            cf = new CutFunc<real_t, cplx_t>[n_thread_];
            Slm<real_t> *slm, *rsum, *ksum;
            slm = new Slm<real_t>[n_thread_];
            rsum = new Slm<real_t>[n_thread_];
            ksum = new Slm<real_t>[n_thread_];
            for (S32 i = 0; i < n_thread_; i++) {
                cf[i].init(2 * p);
                slm[i].alloc(2 * p);
                rsum[i].alloc(2 * p);
                ksum[i].alloc(2 * p);
            }
            {
                Slm<real_t> tmp;
                tmp.make_table(2 * p);
            }
            const S32 size_rspace_loc_tot = gf_r_tmp_part_.size_loc_[1] * gf_r_tmp_part_.size_loc_[2] * gf_r_tmp_part_.size_loc_[3];
            PS_OMP(omp parallel for num_threads(n_thread_))
            for (S32 idx = 0; idx < size_rspace_loc_tot; idx++) {
                // [Notes (tag: #c2f9fbef)]
                //     In order to increase the efficiency of OpenMP parallerization,
                //     we use a loop for cell index instead of the usual (i,j,k) triple
                //     loops. By this, we can obtain a high efficiency even in the case
                //     that the length of the outermost loop of triple loops is too short.

                // Get thread ID
                const S32 ith = comm_info_.getThreadNum();

                // Calculate the cell index in the local region
                S32 idx_tmp = idx;
                const S32 k_loc = idx_tmp / (gf_r_tmp_part_.size_loc_[1] * gf_r_tmp_part_.size_loc_[2]);
                idx_tmp -= k_loc * (gf_r_tmp_part_.size_loc_[1] * gf_r_tmp_part_.size_loc_[2]);
                const S32 j_loc = idx_tmp / gf_r_tmp_part_.size_loc_[1];
                idx_tmp -= j_loc * gf_r_tmp_part_.size_loc_[1];
                const S32 i_loc = idx_tmp;

                // Calculate the cell index in the global region
                const S32 i = i_loc + gf_r_tmp_part_.start_loc_[1];
                const S32 j = j_loc + gf_r_tmp_part_.start_loc_[2];
                const S32 k = k_loc + gf_r_tmp_part_.start_loc_[3];

                // real-space sum
                rsum[ith].clear();
                const S32 kk = (k >= gf_r_tmp_part_.size_glb_[3] / 2) ? k - gf_r_tmp_part_.size_glb_[3] : k;
                const S32 jj = (j >= gf_r_tmp_part_.size_glb_[2] / 2) ? j - gf_r_tmp_part_.size_glb_[2] : j;
                const S32 ii = (i >= gf_r_tmp_part_.size_glb_[1] / 2) ? i - gf_r_tmp_part_.size_glb_[1] : i;
                // [TODO]
                //    In the calculation above, we must also use
                //    gf_r_tmp_part_.start_glb_[] and gf_r_tmp_part_.end_glb_[]
                //    if gf_r_tmp_part_.start_glb_[?] != 0.
                S32 nx, ny, nz;
                for (nz = -NMAX; nz <= NMAX; nz++)
                    for (ny = -NMAX; ny <= NMAX; ny++)
                        for (nx = -NMAX; nx <= NMAX; nx++) {
                            const S32 kkk = kk + nz * gf_r_tmp_part_.size_glb_[3];
                            const S32 jjj = jj + ny * gf_r_tmp_part_.size_glb_[2];
                            const S32 iii = ii + nx * gf_r_tmp_part_.size_glb_[1];
                            if (0 == (iii | jjj | kkk)) continue;
                            const F64 dx = cell_length.x * F64(iii);
                            const F64 dy = cell_length.y * F64(jjj);
                            const F64 dz = cell_length.z * F64(kkk);
                            const F64 dr2 = dx * dx + dy * dy + dz * dz;
                            cf[ith].eval_rcut(dr2 * (alpha * alpha));
                            slm[ith].eval_opt(-dx, dy, dz);  // eval S_l^{-m}
                            // near cell correction
                            const bool near = (abs(kkk) <= icut && abs(jjj) <= icut && abs(iii) <= icut);
                            if (near) {
                                for (S32 l = 0; l <= 2 * p; l++) {
                                    cf[ith].rcut[l] -= 1.0;
                                }
                            }

                            for (S32 l = 0; l <= 2 * p; l++) {
                                for (S32 m = 0; m <= l; m++) {
                                    const cplx_t val = cf[ith].rcut[l] * slm[ith].val_at(l, m);
                                    rsum[ith].accum_at(l, m, val);
                                }
                            }
                        }

                // wave-space sum
                ksum[ith].clear();
                S32 mx, my, mz;
                for (mz = -MMAX; mz <= MMAX; mz++)
                    for (my = -MMAX; my <= MMAX; my++)
                        for (mx = -MMAX; mx <= MMAX; mx++) {
                            if (0 == (mx | my | mz)) continue;
                            const F64 dx = cell_length.x * F64(ii);
                            const F64 dy = cell_length.y * F64(jj);
                            const F64 dz = cell_length.z * F64(kk);
                            //---
                            // The following code is the same as that of the original PMMM code.
                            // But, the use of (i,j,k) is probably WRONG.
                            // const F64 dx = cell_length.x * F64(i);
                            // const F64 dy = cell_length.y * F64(j);
                            // const F64 dz = cell_length.z * F64(k);
                            //---
                            const F64 twopi = 8.0 * std::atan(1.0);
                            const F64 theta = twopi * (dx * mx * bx_inv + dy * my * by_inv + dz * mz * bz_inv);
                            const cplx_t phase(cos(theta), sin(theta));
                            slm[ith].eval_opt(-mx * bx_inv, my * by_inv, mz * bz_inv);
                            const F64 kn2 = mx * mx * bx2_inv + my * my * by2_inv + mz * mz * bz2_inv;
                            cf[ith].eval_kcut(kn2, alpha);
                            for (S32 l = 0; l <= 2 * p; l++) {
                                for (S32 m = 0; m <= l; m++) {
                                    const cplx_t val = (cf[ith].kcut[l] * phase) * slm[ith].val_at(l, m) * vol_inv;
                                    ksum[ith].accum_at(l, m, val);
                                }
                            }
                        }
                // store sum
#ifdef PARTICLE_SIMULATOR_PMMM_IGNORE_RSPACE
                rsum[ith].clear();
#endif
#ifdef PARTICLE_SIMULATOR_PMMM_IGNORE_KSPACE
                ksum[ith].clear();
#endif
                for (S32 lm = gf_r_tmp_part_.start_loc_[0]; lm <= gf_r_tmp_part_.end_loc_[0]; lm++) {
                    const S32 lm_loc = lm - gf_r_tmp_part_.start_loc_[0];
                    const S32 adr =
                        lm_loc + gf_r_tmp_part_.size_loc_[0] * (i_loc + gf_r_tmp_part_.size_loc_[1] * (j_loc + gf_r_tmp_part_.size_loc_[2] * k_loc));
                    assert(0 <= adr && adr < gf_r_tmp_part_.size_loc_tot_);
                    gf_r_tmp_part_.buf_[adr] = rsum[ith].buf[lm] + ksum[ith].buf[lm];
                }
#ifdef PARTICLE_SIMULATOR_PMMM_SHOW_RSUM_AND_KSUM
                if (0 == (i | j | k)) {
                    rsum.show();
                    ksum.show();
                };
#endif
            }  // for(idx)

            // Free memory
            for (S32 i = 0; i < n_thread_; i++) {
                cf[i].freeMem();
                slm[i].freeMem();
                rsum[i].freeMem();
                ksum[i].freeMem();
            }
            delete[] cf;
            delete[] slm;
            delete[] rsum;
            delete[] ksum;
        }
    }

    void redistGFR() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // Memory allocation
        gf_r_tmp_.resize();
        // Prepare convenient variables
        // const S32 p = param_.p;
        // Make a recieve buffer
        union buf_t {
            real_t r;  // to store val
            S32 i;     // to store idx
        };
        CommBuffer<buf_t> recvbuf;
        recvbuf.n_comm = n_proc_in_parent_group_;
        recvbuf.allocCommInfo();
        for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
            recvbuf.ranks[i] = i;
            recvbuf.adr_from_rank[i] = i;
        }
        recvbuf.clearCounts();
        recvbuf.count_tot = 0;
        if (gf_r_tmp_.size_loc_tot_ > 0) {
            S32ort my_domain;
            my_domain.low_.x = gf_r_tmp_.start_loc_[1];
            my_domain.high_.x = gf_r_tmp_.end_loc_[1];
            my_domain.low_.y = gf_r_tmp_.start_loc_[2];
            my_domain.high_.y = gf_r_tmp_.end_loc_[2];
            my_domain.low_.z = gf_r_tmp_.start_loc_[3];
            my_domain.high_.z = gf_r_tmp_.end_loc_[3];
            //if (Comm::getRank() == 0) std::cout << "my_domain= " << my_domain << std::endl;
            for (S32 i = 0; i < n_proc_tot_gf_r_calc_; i++) {  // source rank
                S32ort domain;
                domain.low_.x = local_start_gf_r_calc_[0][i];
                domain.high_.x = local_end_gf_r_calc_[0][i];
                domain.low_.y = local_start_gf_r_calc_[1][i];
                domain.high_.y = local_end_gf_r_calc_[1][i];
                domain.low_.z = local_start_gf_r_calc_[2][i];
                domain.high_.z = local_end_gf_r_calc_[2][i];

                //if (Comm::getRank() == 0) std::cout << "i= " << i << " domain= " << domain << std::endl;
                
                S32ort intxn;
                if (my_domain.calcIntersection(domain, intxn)) {
                    const S32vec size = intxn.getLength();
                    const S32 size_tot = size.x * size.y * size.z;
                    const S32 adr = recvbuf.adr_from_rank[i];
                    const S32 cnt = gf_r_tmp_.size_loc_[0];
                    //if (Comm::getRank() == 0) std::cout << "i= " << i << " size= " << size << " size_tot= " << size_tot << " adr= " << adr << " cnt= " << cnt << std::endl;
                    recvbuf.counts[adr] += size_tot * (1 + cnt);
                    recvbuf.count_tot += size_tot * (1 + cnt);
                }
            }
        }
        if (recvbuf.count_tot > 0) {
            recvbuf.allocBuffer();
            recvbuf.calcDispls();
        }

        //Comm::barrier();
        //if (Comm::getRank() == 0) {
        //    std::cout << "recvbuf.dumpCommInfo();" << std::endl;
        //    recvbuf.dumpCommInfo();
        //}

        // Make a send buffer
        CommBuffer<buf_t> sendbuf;
        sendbuf.n_comm = n_proc_in_parent_group_;
        sendbuf.allocCommInfo();
        for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
            const S32 rank = (rank_in_parent_group_ + i) % n_proc_in_parent_group_;
            sendbuf.ranks[i] = rank;
            sendbuf.adr_from_rank[rank] = i;
        }
        sendbuf.clearCounts();
        sendbuf.count_tot = 0;
        if (gf_r_tmp_part_.size_loc_tot_ > 0) {
            for (S32 k = gf_r_tmp_part_.start_loc_[3]; k <= gf_r_tmp_part_.end_loc_[3]; k++)
                for (S32 j = gf_r_tmp_part_.start_loc_[2]; j <= gf_r_tmp_part_.end_loc_[2]; j++)
                    for (S32 i = gf_r_tmp_part_.start_loc_[1]; i <= gf_r_tmp_part_.end_loc_[1]; i++) {  // for(k,j,i)
                        // S32 slab_id;
                        S32 slab_id = 0;
                        for (S32 n = 0; n < n_proc_for_fft_; n++) {
                            if (local_0_start_[n] <= k && k <= local_0_end_[n]) {
                                slab_id = n;
                                break;
                            }
                        }
                        for (S32 n = 0; n < n_group_; n++) {
                            const S32 rank = rank_start_[n] + slab_id;
                            const S32 adr = sendbuf.adr_from_rank[rank];
                            const S32 cnt = lm_end_[n] - lm_start_[n] + 1;
                            sendbuf.counts[adr] += (1 + cnt);
                            sendbuf.count_tot += (1 + cnt);
                        }
                    }  // for(k,j,i)
        }

        if (sendbuf.count_tot > 0) {
            sendbuf.allocBuffer();
            sendbuf.calcDispls();
            sendbuf.clearCounts();
            for (S32 k = gf_r_tmp_part_.start_loc_[3]; k <= gf_r_tmp_part_.end_loc_[3]; k++)
                for (S32 j = gf_r_tmp_part_.start_loc_[2]; j <= gf_r_tmp_part_.end_loc_[2]; j++)
                    for (S32 i = gf_r_tmp_part_.start_loc_[1]; i <= gf_r_tmp_part_.end_loc_[1]; i++) {  // for(k,j,i)
                        // S32 slab_id;
                        S32 slab_id = 0;
                        for (S32 n = 0; n < n_proc_for_fft_; n++) {
                            if (local_0_start_[n] <= k && k <= local_0_end_[n]) {
                                slab_id = n;
                                break;
                            }
                        }
                        for (S32 n = 0; n < n_group_; n++) {
                            const S32 rank = rank_start_[n] + slab_id;
                            const S32 adr = sendbuf.adr_from_rank[rank];
                            // set the cell index
                            S32 adr_dest{sendbuf.displs[adr] + sendbuf.counts[adr]};
                            const S32 idx = i + gf_r_tmp_part_.size_glb_[1] * (j + gf_r_tmp_part_.size_glb_[2] * k);
                            sendbuf.buf[adr_dest++].i = idx;
                            // set the value of the green function
                            const S32 bgn = lm_start_[n];
                            const S32 end = lm_end_[n];
                            for (S32 lm = bgn; lm <= end; lm++) {
                                const S32 lm_loc = lm - gf_r_tmp_part_.start_loc_[0];
                                const S32 i_loc = i - gf_r_tmp_part_.start_loc_[1];
                                const S32 j_loc = j - gf_r_tmp_part_.start_loc_[2];
                                const S32 k_loc = k - gf_r_tmp_part_.start_loc_[3];
                                const S32 adr_src =
                                    lm_loc + gf_r_tmp_part_.size_loc_[0] *
                                                 (i_loc + gf_r_tmp_part_.size_loc_[1] * (j_loc + gf_r_tmp_part_.size_loc_[2] * k_loc));
                                assert(0 <= adr_src && adr_src < gf_r_tmp_part_.size_loc_tot_);
                                sendbuf.buf[adr_dest++].r = gf_r_tmp_part_.buf_[adr_src];
                            }
                            sendbuf.counts[adr] += (1 + end - bgn + 1);
                        }
                    }  // for(k,j,i)
        }
        // Perform MPI comm.
        //Comm::barrier();
        //if (Comm::getRank() == 0) {
        //    std::cout << "sendbuf.dumpCommInfo();" << std::endl;
        //    sendbuf.dumpCommInfo();
        //}

        performComm(sendbuf, recvbuf, parent_comm_);
        // performComm(sendbuf, recvbuf, MPI_COMM_WORLD);

        // Copy from the recieve buffer to gf_r_tmp_
        if (recvbuf.count_tot > 0) {
            for (S32 i = 0; i < recvbuf.n_comm; i++) {
                const S32 rank = recvbuf.ranks[i];
                const S32 adr = recvbuf.adr_from_rank[rank];
                if (recvbuf.counts[adr] > 0) {
                    const S32 head = recvbuf.displs[adr];
                    const S32 tail = head + recvbuf.counts[adr];
                    S32 adr_src{head};
                    while (adr_src < tail) {
                        const S32 idx_1d = recvbuf.buf[adr_src++].i;
                        S32 ix, iy, iz;
                        S32 idx_tmp = idx_1d;
                        iz = idx_tmp / (gf_r_tmp_.size_glb_[1] * gf_r_tmp_.size_glb_[2]);
                        idx_tmp -= iz * (gf_r_tmp_.size_glb_[1] * gf_r_tmp_.size_glb_[2]);
                        iy = idx_tmp / gf_r_tmp_.size_glb_[1];
                        idx_tmp -= iy * gf_r_tmp_.size_glb_[1];
                        ix = idx_tmp;
                        for (S32 lm = gf_r_tmp_.start_loc_[0]; lm <= gf_r_tmp_.end_loc_[0]; lm++) {
                            const S32 lm_loc = lm - gf_r_tmp_.start_loc_[0];
                            const S32 i_loc = ix - gf_r_tmp_.start_loc_[1];
                            const S32 j_loc = iy - gf_r_tmp_.start_loc_[2];
                            const S32 k_loc = iz - gf_r_tmp_.start_loc_[3];
                            const S32 adr_dest =
                                lm_loc + gf_r_tmp_.size_loc_[0] * (i_loc + gf_r_tmp_.size_loc_[1] * (j_loc + gf_r_tmp_.size_loc_[2] * k_loc));
                            assert(0 <= adr_dest && adr_dest < gf_r_tmp_.size_loc_tot_);
                            assert(0 <= adr_src && adr_src < recvbuf.count_tot);
                            gf_r_tmp_.buf_[adr_dest] = recvbuf.buf[adr_src++].r;
                        }
                    }
                    assert(adr_src == tail);
                }
            }
        }

        // Free memory
        gf_r_tmp_part_.free();
#if 0
                // For debug
                {
                    if (gf_r_tmp_.size_loc_tot_ > 0) {
                        std::stringstream ss;
                        ss << "gf_r_tmp_new_" << std::setw(5) << std::setfill('0')
                           << rank_in_parent_group_ << ".txt";
                        const std::string file_name = ss.str();
                        std::ofstream ofs;
                        ofs.open(file_name.c_str(), std::ios::trunc);
                        for (S32 k  = gf_r_tmp_.start_loc_[3]; k  <= gf_r_tmp_.end_loc_[3]; k++)
                        for (S32 j  = gf_r_tmp_.start_loc_[2]; j  <= gf_r_tmp_.end_loc_[2]; j++)
                        for (S32 i  = gf_r_tmp_.start_loc_[1]; i  <= gf_r_tmp_.end_loc_[1]; i++)
                        for (S32 lm = gf_r_tmp_.start_loc_[0]; lm <= gf_r_tmp_.end_loc_[0]; lm++)
                        { // for(lm,i,j,k)
                            const S32 lm_loc = lm - gf_r_tmp_.start_loc_[0];
                            const S32 i_loc  = i  - gf_r_tmp_.start_loc_[1];
                            const S32 j_loc  = j  - gf_r_tmp_.start_loc_[2];
                            const S32 k_loc  = k  - gf_r_tmp_.start_loc_[3];
                            const S32 adr = lm_loc
                                          + gf_r_tmp_.size_loc_[0] * (i_loc
                                          + gf_r_tmp_.size_loc_[1] * (j_loc
                                          + gf_r_tmp_.size_loc_[2] * k_loc));
                            assert(0 <= adr && adr < gf_r_tmp_.size_loc_tot_);
                            ofs << gf_r_tmp_.buf_[adr] << std::endl;
                        } // for(lm,i,j,k)
                        ofs.close();
                    }
                }
#endif
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
        assert(gf_r_tmp_.size_loc_tot_ == gf_r_tmp_part_.size_loc_tot_);
        gf_r_tmp_.moveFrom(gf_r_tmp_part_);
        // [Note (tag: #04dd61ad)]
        //     In the sequantial program, the size information of
        //     gf_r_tmp_ is exactly the same as that of gf_r_tmp_part_.
        //     Hence, we do not need to perform copy.
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
    }

    void calcGreenFunctionInWavenumberSpace() {
        // Memory allocation
        gf_k_tmp_.resize();
        // Extract basic parameters
        // const S32 p = param_.p;
        // const S32vec n_cell = param_.n_cell;
        // Local variables
        fftw_real_t *rbuf = nullptr;
        fftw_cplx_t *kbuf = nullptr;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if (use_mpifft_) {
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION_IN_WAVENUMBER_SPACE
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok1a @ calcGreenFunctionInWavenumberSpace() [PMMM]" << std::endl;
#endif
            if (rank_in_my_group_ < n_proc_for_fft_) {
                // Check if conditions that need to be fulfilled in this IF branch
                // are actually fulfilled.
                assert(gf_r_tmp_.start_loc_[0] == gf_k_tmp_.start_loc_[0]);
                assert(gf_r_tmp_.end_loc_[0] == gf_k_tmp_.end_loc_[0]);
                assert(gf_r_tmp_.start_loc_[0] == gf_k_tmp_.start_loc_[0]);
                assert(gf_r_tmp_.end_loc_[0] == gf_k_tmp_.end_loc_[0]);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION_IN_WAVENUMBER_SPACE
                MPI_Barrier(comm_fft_);
                if (rank_in_my_group_ == 0) DebugUtils::fout() << "ok1a-1 @ calcGreenFunctionInWavenumberSpace() [PMMM]" << std::endl;
#endif

                // Allocate buffer memory
                ptrdiff_t alloc_local;
                ptrdiff_t local_n0, local_0_start;
                ptrdiff_t local_n1, local_1_start;
                alloc_local = fftw_mpi_local_size_3d_transposed(gf_k_tmp_.size_glb_[3], gf_k_tmp_.size_glb_[2], gf_k_tmp_.size_glb_[1], comm_fft_,
                                                                &local_n0, &local_0_start, &local_n1, &local_1_start);
                // Note that an user must indicate the logical size of FFT by
                // the sizes of a COMPLEX array when using fftw_mpi_local_size_3d_transposed.
                const ptrdiff_t size_rbuf = 2 * alloc_local;
                const ptrdiff_t size_kbuf = alloc_local;
                rbuf = fftw_alloc_real(size_rbuf);
                kbuf = fftw_alloc_complex(size_kbuf);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION_IN_WAVENUMBER_SPACE
                MPI_Barrier(comm_fft_);
                if (rank_in_my_group_ == 0) DebugUtils::fout() << "ok1a-2 @ calcGreenFunctionInWavenumberSpace() [PMMM]" << std::endl;
#endif

                // Create plan
                fftw_plan plan_fwd = fftw_mpi_plan_dft_r2c_3d(gf_r_tmp_.size_glb_[3], gf_r_tmp_.size_glb_[2], gf_r_tmp_.size_glb_[1], &rbuf[0],
                                                              &kbuf[0], comm_fft_, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | FFTW_MPI_TRANSPOSED_OUT);
                // Note that fftw_mpi_plan_dft_r2c_3d requires an user to indicate
                // the logical size of FFT by the size of a REAL array.

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION_IN_WAVENUMBER_SPACE
                MPI_Barrier(comm_fft_);
                if (rank_in_my_group_ == 0) DebugUtils::fout() << "ok1a-3 @ calcGreenFunctionInWavenumberSpace() [PMMM]" << std::endl;
#endif

                // Perform FFT
                for (S32 lm = gf_r_tmp_.start_loc_[0]; lm <= gf_r_tmp_.end_loc_[0]; lm++) {
                    const S32 lm_loc = lm - gf_r_tmp_.start_loc_[0];
                    for (S32 k = gf_r_tmp_.start_loc_[3]; k <= gf_r_tmp_.end_loc_[3]; k++)
                        for (S32 j = gf_r_tmp_.start_loc_[2]; j <= gf_r_tmp_.end_loc_[2]; j++)
                            for (S32 i = gf_r_tmp_.start_loc_[1]; i <= gf_r_tmp_.end_loc_[1]; i++) {  // for(k,j,i)
                                const S32 i_loc = i - gf_r_tmp_.start_loc_[1];
                                const S32 j_loc = j - gf_r_tmp_.start_loc_[2];
                                const S32 k_loc = k - gf_r_tmp_.start_loc_[3];
                                const S32 adr_src =
                                    lm_loc + gf_r_tmp_.size_loc_[0] * (i_loc + gf_r_tmp_.size_loc_[1] * (j_loc + gf_r_tmp_.size_loc_[2] * k_loc));
                                const S32 i_buf = i - gf_r_tmp_.start_loc_[1];
                                const S32 j_buf = j - gf_r_tmp_.start_loc_[2];
                                const S32 k_buf = k - gf_r_tmp_.start_loc_[3];
                                const S32 adr_dest = i_buf + (2 * (gf_r_tmp_.size_glb_[1] / 2 + 1)) * (j_buf + gf_r_tmp_.size_glb_[2] * k_buf);
                                // Note that the size of 1st dimension of rbuf is not NX.
                                // (see Section 6.5 of the mamual of FFTW3)
                                rbuf[adr_dest] = gf_r_tmp_.buf_[adr_src];
                            }  // for(k,j,i)

                    // CALL FFTW
                    fftw_execute(plan_fwd);

                    for (S32 j = gf_k_tmp_.start_loc_[3]; j <= gf_k_tmp_.end_loc_[3]; j++)
                        for (S32 k = gf_k_tmp_.start_loc_[2]; k <= gf_k_tmp_.end_loc_[2]; k++)
                            for (S32 i = gf_k_tmp_.start_loc_[1]; i <= gf_k_tmp_.end_loc_[1]; i++) {  // for(k,j,i)
                                const S32 i_buf = i - gf_k_tmp_.start_loc_[1];
                                const S32 k_buf = k - gf_k_tmp_.start_loc_[2];
                                const S32 j_buf = j - gf_k_tmp_.start_loc_[3];
                                const S32 adr_src = i_buf + gf_k_tmp_.size_glb_[1] * (k_buf + gf_k_tmp_.size_glb_[2] * j_buf);
                                const S32 i_loc = i - gf_k_tmp_.start_loc_[1];
                                const S32 k_loc = k - gf_k_tmp_.start_loc_[2];
                                const S32 j_loc = j - gf_k_tmp_.start_loc_[3];
                                const S32 adr_dest =
                                    lm_loc + gf_k_tmp_.size_loc_[0] * (i_loc + gf_k_tmp_.size_loc_[1] * (k_loc + gf_k_tmp_.size_loc_[2] * j_loc));
                                // Note that j and k loops are exchanged because we specified
                                // FFTW_MPI_TRANSPOSED_OUT.
                                assert(0 <= adr_src && adr_src < size_kbuf);
                                assert(0 <= adr_dest && adr_dest < gf_k_tmp_.size_loc_tot_);
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                const cplx_t ctmp(kbuf[adr_src][0], kbuf[adr_src][1]);
                                gf_k_tmp_.buf_[adr_dest] = ctmp;
#else
                                gf_k_tmp_.buf_[adr_dest].real(kbuf[adr_src][0]);
                                gf_k_tmp_.buf_[adr_dest].imag(kbuf[adr_src][1]);
#endif
                                // Note that fftw_complex = double[2]
                            }  // for(k,j,i)
                }
                // Destropy the plan
                fftw_destroy_plan(plan_fwd);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION_IN_WAVENUMBER_SPACE
                MPI_Barrier(comm_fft_);
                if (rank_in_my_group_ == 0) DebugUtils::fout() << "ok1a-last @ calcGreenFunctionInWavenumberSpace() [PMMM]" << std::endl;
#endif
            }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION_IN_WAVENUMBER_SPACE
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok2a @ calcGreenFunctionInWavenumberSpace() [PMMM]" << std::endl;
#endif
        } else {
#endif
            // In this case, we can use only a single MPI process to calculate
            // gf_k_tmp_ for a particular set of (l,m).
            // We use OpenMP parallerization of FFT if available.
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION_IN_WAVENUMBER_SPACE
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok1b @ calcGreenFunctionInWavenumberSpace() [PMMM]" << std::endl;
#endif
            if (rank_in_my_group_ == 0) {
                // Check if conditions that need to be fulfilled in this IF branch
                // are actually fulfilled.
                assert(gf_r_tmp_.start_loc_[0] == gf_k_tmp_.start_loc_[0]);
                assert(gf_r_tmp_.end_loc_[0] == gf_k_tmp_.end_loc_[0]);
                for (S32 i = 1; i <= 3; i++) {
                    assert(gf_r_tmp_.size_loc_[i] == gf_r_tmp_.size_glb_[i]);
                    assert(gf_k_tmp_.size_loc_[i] == gf_k_tmp_.size_glb_[i]);
                }

                // Calculate gf_k
                const S32 fft_size = gf_r_tmp_.size_glb_[3] * gf_r_tmp_.size_glb_[2] * gf_r_tmp_.size_glb_[1];
                if (!use_ompfft_) {
                    // In this case, the size of array is too small to speed up by multithreaded FFT.
                    // Hence, we assign a single FFT to each thread.

                    // Allocate buffer memory
                    const S32 size_rbuf = n_thread_ * fft_size;
                    const S32 size_kbuf = n_thread_ * fft_size;
                    rbuf = fftw_alloc_real(size_rbuf);
                    kbuf = fftw_alloc_complex(size_kbuf);

                    // Create plan
                    fftw_plan *plan_fwd;
                    plan_fwd = (fftw_plan *)malloc((size_t)sizeof(fftw_plan) * n_thread_);
                    for (S32 i = 0; i < n_thread_; i++) {
                        const S32 offset = fft_size * i;
                        plan_fwd[i] = fftw_plan_dft_r2c_3d(gf_r_tmp_.size_glb_[3], gf_r_tmp_.size_glb_[2], gf_r_tmp_.size_glb_[1], &rbuf[offset],
                                                           &kbuf[offset], FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
                    }

                    // Perform FFT
                    PS_OMP(omp parallel for num_threads(n_thread_))
                    for (S32 lm = gf_r_tmp_.start_loc_[0]; lm <= gf_r_tmp_.end_loc_[0]; lm++) {
                        const S32 ith = comm_info_.getThreadNum();
                        const S32 offset = fft_size * ith;

                        const S32 lm_loc = lm - gf_r_tmp_.start_loc_[0];
                        for (S32 k = gf_r_tmp_.start_loc_[3]; k <= gf_r_tmp_.end_loc_[3]; k++)
                            for (S32 j = gf_r_tmp_.start_loc_[2]; j <= gf_r_tmp_.end_loc_[2]; j++)
                                for (S32 i = gf_r_tmp_.start_loc_[1]; i <= gf_r_tmp_.end_loc_[1]; i++) {  // for(k,j,i)
                                    const S32 i_loc = i - gf_r_tmp_.start_loc_[1];
                                    const S32 j_loc = j - gf_r_tmp_.start_loc_[2];
                                    const S32 k_loc = k - gf_r_tmp_.start_loc_[3];
                                    const S32 adr_src =
                                        lm_loc + gf_r_tmp_.size_loc_[0] * (i_loc + gf_r_tmp_.size_loc_[1] * (j_loc + gf_r_tmp_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i + gf_r_tmp_.size_glb_[1] * (j + gf_r_tmp_.size_glb_[2] * k);
#ifdef PARTICLE_SIMULATOR_PMMM_RANGE_CHECK_GREEN_FUNC
                                    assert(0 <= adr_src && adr_src < gf_r_tmp_.size_loc_tot_);
                                    assert(0 <= (adr_dest + offset) && (adr_dest + offset) < size_rbuf);
#endif
                                    rbuf[adr_dest + offset] = gf_r_tmp_.buf_[adr_src];
                                }

                        // CALL FFTW
                        fftw_execute(plan_fwd[ith]);

                        for (S32 k = gf_k_tmp_.start_loc_[3]; k <= gf_k_tmp_.end_loc_[3]; k++)
                            for (S32 j = gf_k_tmp_.start_loc_[2]; j <= gf_k_tmp_.end_loc_[2]; j++)
                                for (S32 i = gf_k_tmp_.start_loc_[1]; i <= gf_k_tmp_.end_loc_[1]; i++) {  // for(k,j,i)
                                    const S32 i_loc = i - gf_k_tmp_.start_loc_[1];
                                    const S32 j_loc = j - gf_k_tmp_.start_loc_[2];
                                    const S32 k_loc = k - gf_k_tmp_.start_loc_[3];
                                    const S32 adr_dest =
                                        lm_loc + gf_k_tmp_.size_loc_[0] * (i_loc + gf_k_tmp_.size_loc_[1] * (j_loc + gf_k_tmp_.size_loc_[2] * k_loc));
                                    const S32 adr_src = i + gf_k_tmp_.size_glb_[1] * (j + gf_k_tmp_.size_glb_[2] * k);
#ifdef PARTICLE_SIMULATOR_PMMM_RANGE_CHECK_GREEN_FUNC
                                    assert(0 <= adr_dest && adr_dest < gf_k_tmp_.size_loc_tot_);
                                    assert(0 <= (adr_src + offset) && (adr_src + offset) < size_kbuf);
#endif
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                    const cplx_t ctmp(kbuf[adr_src + offset][0], kbuf[adr_src + offset][1]);
                                    gf_k_tmp_.buf_[adr_dest] = ctmp;
#else
                                    gf_k_tmp_.buf_[adr_dest].real(kbuf[adr_src + offset][0]);
                                    gf_k_tmp_.buf_[adr_dest].imag(kbuf[adr_src + offset][1]);
#endif
                                    // Note that fftw_complex = double[2]
                                }  // for(k,j,i)
                    }

                    // Destroy the plan
                    for (S32 i = 0; i < n_thread_; i++) fftw_destroy_plan(plan_fwd[i]);
                    // delete[] plan_fwd;
                    free(plan_fwd);

                } else {
                    // In this case, we can expect that the FFT calculation become fast by using
                    // a multithreaded FFT.

                    // Allocate buffer memory
                    const S32 size_rbuf = gf_r_tmp_.size_glb_[3] * gf_r_tmp_.size_glb_[2] * gf_r_tmp_.size_glb_[1];
                    const S32 size_kbuf = gf_k_tmp_.size_glb_[3] * gf_k_tmp_.size_glb_[2] * gf_k_tmp_.size_glb_[1];
                    rbuf = fftw_alloc_real(size_rbuf);
                    kbuf = fftw_alloc_complex(size_kbuf);

                    // Create plan
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                    fftw_plan_with_nthreads(n_thread_);
#endif

                    fftw_plan plan_fwd = fftw_plan_dft_r2c_3d(gf_r_tmp_.size_glb_[3], gf_r_tmp_.size_glb_[2], gf_r_tmp_.size_glb_[1], &rbuf[0],
                                                              &kbuf[0], FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

                    // Perform FFT
                    for (S32 lm = gf_r_tmp_.start_loc_[0]; lm <= gf_r_tmp_.end_loc_[0]; lm++) {
                        const S32 lm_loc = lm - gf_r_tmp_.start_loc_[0];
                        PS_OMP(omp parallel for num_threads(n_thread_))
                        for (S32 k = gf_r_tmp_.start_loc_[3]; k <= gf_r_tmp_.end_loc_[3]; k++)
                            for (S32 j = gf_r_tmp_.start_loc_[2]; j <= gf_r_tmp_.end_loc_[2]; j++)
                                for (S32 i = gf_r_tmp_.start_loc_[1]; i <= gf_r_tmp_.end_loc_[1]; i++) {  // for (k,j,i)
                                    const S32 i_loc = i - gf_r_tmp_.start_loc_[1];
                                    const S32 j_loc = j - gf_r_tmp_.start_loc_[2];
                                    const S32 k_loc = k - gf_r_tmp_.start_loc_[3];
                                    const S32 adr_src =
                                        lm_loc + gf_r_tmp_.size_loc_[0] * (i_loc + gf_r_tmp_.size_loc_[1] * (j_loc + gf_r_tmp_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i + gf_r_tmp_.size_glb_[1] * (j + gf_r_tmp_.size_glb_[2] * k);
#ifdef PARTICLE_SIMULATOR_PMMM_RANGE_CHECK_GREEN_FUNC
                                    assert(0 <= adr_dest && adr_dest < size_rbuf);
                                    assert(0 <= adr_src && adr_src < gf_r_tmp_.size_loc_tot_);
#endif
                                    rbuf[adr_dest] = gf_r_tmp_.buf_[adr_src];
                                }  // for(k,j,i)

                        // CALL FFTW
                        fftw_execute(plan_fwd);

                        PS_OMP(omp parallel for num_threads(n_thread_))
                        for (S32 k = gf_k_tmp_.start_loc_[3]; k <= gf_k_tmp_.end_loc_[3]; k++)
                            for (S32 j = gf_k_tmp_.start_loc_[2]; j <= gf_k_tmp_.end_loc_[2]; j++)
                                for (S32 i = gf_k_tmp_.start_loc_[1]; i <= gf_k_tmp_.end_loc_[1]; i++) {  // for(k,j,i)
                                    const S32 i_loc = i - gf_k_tmp_.start_loc_[1];
                                    const S32 j_loc = j - gf_k_tmp_.start_loc_[2];
                                    const S32 k_loc = k - gf_k_tmp_.start_loc_[3];
                                    const S32 adr_dest =
                                        lm_loc + gf_k_tmp_.size_loc_[0] * (i_loc + gf_k_tmp_.size_loc_[1] * (j_loc + gf_k_tmp_.size_loc_[2] * k_loc));
                                    const S32 adr_src = i + gf_k_tmp_.size_glb_[1] * (j + gf_k_tmp_.size_glb_[2] * k);
#ifdef PARTICLE_SIMULATOR_PMMM_RANGE_CHECK_GREEN_FUNC
                                    assert(0 <= adr_dest && adr_dest < gf_k_tmp_.size_loc_tot_);
                                    assert(0 <= adr_src && adr_src < size_kbuf);
#endif
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                    const cplx_t ctmp(kbuf[adr_src][0], kbuf[adr_src][1]);
                                    gf_k_tmp_.buf_[adr_dest] = ctmp;
#else
                                    gf_k_tmp_.buf_[adr_dest].real(kbuf[adr_src][0]);
                                    gf_k_tmp_.buf_[adr_dest].imag(kbuf[adr_src][1]);
#endif
                                    // Note that fftw_complex = double[2]
                                }
                    }

                    // Destroy the plan
                    fftw_destroy_plan(plan_fwd);
                }  // END of if (rank_in_my_group_ == 0)
            }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_GFC_CALC_GREEN_FUNCTION_IN_WAVENUMBER_SPACE
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok2b @ calcGreenFunctionInWavenumberSpace() [PMMM]" << std::endl;
#endif
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        }
#endif
        // Free
        if (rbuf != nullptr) fftw_free(rbuf);
        if (kbuf != nullptr) fftw_free(kbuf);
        gf_r_tmp_.free();
    }

    void redistGFK(Buffer<cplx_t> &gf_k, const bool is_gf_k_transposed) {
        // Memory allocation
        gf_k.resize();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // Prepare convenient variables
        S32vec size_glb_kspace;
        size_glb_kspace.x = gf_k_tmp_.size_glb_[1];
        if (use_mpifft_) {
            size_glb_kspace.y = gf_k_tmp_.size_glb_[3];
            size_glb_kspace.z = gf_k_tmp_.size_glb_[2];
        } else {
            size_glb_kspace.y = gf_k_tmp_.size_glb_[2];
            size_glb_kspace.z = gf_k_tmp_.size_glb_[3];
        }
        // Share gf_k.size_loc_tot_, gf_k.start_loc_[], and gf_k.end_loc_[]
        S32 *size_loc_req_tbl = new S32[n_proc_in_parent_group_];
        S32 *start_loc_req_tbl = new S32[3 * n_proc_in_parent_group_];
        S32 *end_loc_req_tbl = new S32[3 * n_proc_in_parent_group_];
        MPI_Allgather(&gf_k.size_loc_tot_, 1, MPI_INT, size_loc_req_tbl, 1, MPI_INT, parent_comm_);
        MPI_Allgather(&gf_k.start_loc_[1], 3, MPI_INT, start_loc_req_tbl, 3, MPI_INT, parent_comm_);
        MPI_Allgather(&gf_k.end_loc_[1], 3, MPI_INT, end_loc_req_tbl, 3, MPI_INT, parent_comm_);
        // Make a recieve buffer
        union buf_t {
            cplx_t c;  // to store val
            S32 i;     // to store idx

            buf_t() : i(0) {}
            ~buf_t() {}
        };
        // [Ref.] https://cpprefjp.github.io/lang/cpp11/unrestricted_unions.html
        CommBuffer<buf_t> recvbuf;
        recvbuf.n_comm = n_proc_in_parent_group_;
        recvbuf.allocCommInfo();
        for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
            recvbuf.ranks[i] = i;
            recvbuf.adr_from_rank[i] = i;
        }
        recvbuf.clearCounts();
        recvbuf.count_tot = 0;
        if (gf_k.size_loc_tot_ > 0) {
            for (S32 k = gf_k.start_loc_[3]; k <= gf_k.end_loc_[3]; k++)
                for (S32 j = gf_k.start_loc_[2]; j <= gf_k.end_loc_[2]; j++)
                    for (S32 i = gf_k.start_loc_[1]; i <= gf_k.end_loc_[1]; i++) {  // for(k,j,i)
                        // Compute physical index
                        // S32 ix, iy, iz;
                        // S32 ix = 0, iy = 0, iz = 0;
                        S32 iy = 0, iz = 0;
                        if (is_gf_k_transposed) {
                            // ix = i;
                            iy = k;
                            iz = j;
                        } else {
                            // ix = i;
                            iy = j;
                            iz = k;
                        }
                        // Find which rank in each group is in charge of a cell (ix,iy,iz)
                        // (In other words, find slab id)
                        // S32 slab_id;
                        S32 slab_id = 0;
                        if (use_mpifft_) {
                            // In this case, slab decomposition along y direction
                            for (S32 n = 0; n < n_proc_for_fft_; n++) {
                                if ((local_1_start_[n] <= iy) && (iy <= local_1_end_[n])) {
                                    slab_id = n;
                                    break;
                                }
                            }
                        } else {
                            // In this case, slab decomposition along z direction
                            for (S32 n = 0; n < n_proc_for_fft_; n++) {
                                if ((local_0_start_[n] <= iz) && (iz <= local_0_end_[n])) {
                                    slab_id = n;
                                    break;
                                }
                            }
                        }
                        // Find ranks are responsible for this cell (ix,iy,iz)
                        for (S32 n = 0; n < n_group_; n++) {
                            const S32 rank = rank_start_[n] + slab_id;
                            const S32 cnt = lm_end_[n] - lm_start_[n] + 1;
                            recvbuf.counts[rank] += (1 + cnt);
                            recvbuf.count_tot += (1 + cnt);
                        }
                    }  // for(k,j,i)
        }
        if (recvbuf.count_tot > 0) {
            recvbuf.allocBuffer();
            recvbuf.calcDispls();
        }
        // Make a send buffer
        CommBuffer<buf_t> sendbuf;
        sendbuf.n_comm = n_proc_in_parent_group_;
        sendbuf.allocCommInfo();
        for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
            const S32 rank = (rank_in_parent_group_ + i) % n_proc_in_parent_group_;
            sendbuf.ranks[i] = rank;
            sendbuf.adr_from_rank[rank] = i;
        }
        sendbuf.clearCounts();
        sendbuf.count_tot = 0;
        // calculate the range of a region covered by gf_k_tmp_
        S32ort my_domain;
        if (rank_in_my_group_ < n_proc_for_fft_) {
            // In this case, this process has data of gf_k_tmp_.
            my_domain.low_.x = gf_k_tmp_.start_loc_[1];
            my_domain.high_.x = gf_k_tmp_.end_loc_[1];
            if (use_mpifft_) {
                my_domain.low_.y = gf_k_tmp_.start_loc_[3];
                my_domain.high_.y = gf_k_tmp_.end_loc_[3];
                my_domain.low_.z = gf_k_tmp_.start_loc_[2];
                my_domain.high_.z = gf_k_tmp_.end_loc_[2];
            } else {
                my_domain.low_.y = gf_k_tmp_.start_loc_[2];
                my_domain.high_.y = gf_k_tmp_.end_loc_[2];
                my_domain.low_.z = gf_k_tmp_.start_loc_[3];
                my_domain.high_.z = gf_k_tmp_.end_loc_[3];
            }
            for (S32 n = 0; n < n_proc_in_parent_group_; n++) {  // destination rank
                // Calculate a region required by rank `n`
                if (size_loc_req_tbl[n] > 0) {
                    // In this case, rank `n` requires information.
                    S32ort domain;
                    domain.low_.x = start_loc_req_tbl[3 * n + 0];
                    domain.high_.x = end_loc_req_tbl[3 * n + 0];
                    if (is_gf_k_transposed) {
                        domain.low_.y = start_loc_req_tbl[3 * n + 2];
                        domain.high_.y = end_loc_req_tbl[3 * n + 2];
                        domain.low_.z = start_loc_req_tbl[3 * n + 1];
                        domain.high_.z = end_loc_req_tbl[3 * n + 1];
                    } else {
                        domain.low_.y = start_loc_req_tbl[3 * n + 1];
                        domain.high_.y = end_loc_req_tbl[3 * n + 1];
                        domain.low_.z = start_loc_req_tbl[3 * n + 2];
                        domain.high_.z = end_loc_req_tbl[3 * n + 2];
                    }
                    // Check if there is a intersection between my_domain and domain
                    // and calculate # of cells in the intersection.
                    S32ort intxn;
                    if (my_domain.calcIntersection(domain, intxn)) {
                        const S32vec size = intxn.getLength();
                        const S32 size_tot = size.x * size.y * size.z;
                        const S32 adr = sendbuf.adr_from_rank[n];
                        const S32 cnt = lm_end_[idx_to_my_group_] - lm_start_[idx_to_my_group_] + 1;
                        sendbuf.counts[adr] += (1 + cnt) * size_tot;
                        sendbuf.count_tot += (1 + cnt) * size_tot;
                    }
                }
            }
        }
        if (sendbuf.count_tot > 0) {
            // The program must pass this branch if and only if
            // rank_in_my_group_ < n_proc_for_fft_.
            assert(rank_in_my_group_ < n_proc_for_fft_);
            sendbuf.allocBuffer();
            sendbuf.calcDispls();
            for (S32 n = 0; n < n_proc_in_parent_group_; n++) {
                if (size_loc_req_tbl[n] > 0) {
                    // std::cout << "n = " << n << std::endl;
                    //  Calculate a region required by rank `n`
                    S32ort domain;
                    domain.low_.x = start_loc_req_tbl[3 * n + 0];
                    domain.high_.x = end_loc_req_tbl[3 * n + 0];
                    if (is_gf_k_transposed) {
                        domain.low_.y = start_loc_req_tbl[3 * n + 2];
                        domain.high_.y = end_loc_req_tbl[3 * n + 2];
                        domain.low_.z = start_loc_req_tbl[3 * n + 1];
                        domain.high_.z = end_loc_req_tbl[3 * n + 1];
                    } else {
                        domain.low_.y = start_loc_req_tbl[3 * n + 1];
                        domain.high_.y = end_loc_req_tbl[3 * n + 1];
                        domain.low_.z = start_loc_req_tbl[3 * n + 2];
                        domain.high_.z = end_loc_req_tbl[3 * n + 2];
                    }
                    // Copy the values of green function in the intersection
                    // to send buffer.
                    S32ort intxn;
                    if (my_domain.calcIntersection(domain, intxn)) {
                        // Copy to send buffer
                        const S32 adr = sendbuf.adr_from_rank[n];
                        S32 adr_dest = sendbuf.displs[adr];
                        for (S32 iz = intxn.low_.z; iz <= intxn.high_.z; iz++)
                            for (S32 iy = intxn.low_.y; iy <= intxn.high_.y; iy++)
                                for (S32 ix = intxn.low_.x; ix <= intxn.high_.x; ix++) {  // for(iz,iy,ix)
                                    // std::cout << "(ix,iy,iz) = " << ix << ", " << iy << ", " << iz << std::endl;
                                    //  Copy a real index to send buffer
                                    const S32 idx = ix + size_glb_kspace.x * (iy + size_glb_kspace.y * iz);
                                    sendbuf.buf[adr_dest++].i = idx;
                                    // Copy the values of green function to send buffer
                                    const S32 bgn = lm_start_[idx_to_my_group_];
                                    const S32 end = lm_end_[idx_to_my_group_];
                                    for (S32 lm = bgn; lm <= end; lm++) {
                                        const S32 lm_src = lm - gf_k_tmp_.start_loc_[0];
                                        const S32 i_src = ix - gf_k_tmp_.start_loc_[1];
                                        S32 j_src, k_src;
                                        if (use_mpifft_) {
                                            j_src = iz - gf_k_tmp_.start_loc_[2];
                                            k_src = iy - gf_k_tmp_.start_loc_[3];
                                        } else {
                                            j_src = iy - gf_k_tmp_.start_loc_[2];
                                            k_src = iz - gf_k_tmp_.start_loc_[3];
                                        }
                                        const S32 adr_src = lm_src + gf_k_tmp_.size_loc_[0] *
                                                                         (i_src + gf_k_tmp_.size_loc_[1] * (j_src + gf_k_tmp_.size_loc_[2] * k_src));
                                        assert(0 <= adr_src && adr_src < gf_k_tmp_.size_loc_tot_);
                                        sendbuf.buf[adr_dest++].c = gf_k_tmp_.buf_[adr_src];
                                    }  // for(lm)
                                }  // for(iz,iy,ix)
                    }
                }
            }
        }
        // Perform MPI comm.
        performComm(sendbuf, recvbuf, parent_comm_);
        // Copy from recvbuf to gf_k
        if (recvbuf.count_tot > 0) {
            for (S32 n = 0; n < recvbuf.n_comm; n++) {
                const S32 rank = recvbuf.ranks[n];
                // Find a group
                // S32 idx_to_group;
                S32 idx_to_group = 0;
                for (S32 k = 0; k < n_group_; k++) {
                    if ((rank_start_[k] <= rank) && (rank <= rank_end_[k])) {
                        idx_to_group = k;
                        break;
                    }
                }
                // Copy
                const S32 head = recvbuf.displs[n];
                const S32 tail = head + recvbuf.counts[n];
                S32 adr_src{head};
                while (adr_src < tail) {
                    // (1) Convert recvbuf.buf[].i to (ix,iy,iz).
                    const S32 idx = recvbuf.buf[adr_src++].i;
                    S32 iz, iy, ix;
                    S32 idx_tmp = idx;
                    iz = idx_tmp / (size_glb_kspace.x * size_glb_kspace.y);
                    idx_tmp -= (size_glb_kspace.x * size_glb_kspace.y) * iz;
                    iy = idx_tmp / size_glb_kspace.x;
                    idx_tmp -= size_glb_kspace.x * iy;
                    ix = idx_tmp;
                    // (3) Copy
                    const S32 beg = lm_start_[idx_to_group];
                    const S32 end = lm_end_[idx_to_group];
                    for (S32 lm = beg; lm <= end; lm++) {
                        // (3-1) convert idx_3d to the one-dimensional cell index
                        //       based on the size information of gf_k_trans_
                        const S32 lm_loc = lm - gf_k.start_loc_[0];
                        const S32 i_loc = ix - gf_k.start_loc_[1];
                        S32 j_loc, k_loc;
                        if (is_gf_k_transposed) {
                            j_loc = iz - gf_k.start_loc_[2];
                            k_loc = iy - gf_k.start_loc_[3];
                        } else {
                            j_loc = iy - gf_k.start_loc_[2];
                            k_loc = iz - gf_k.start_loc_[3];
                        }
                        const S32 adr_dest = lm_loc + gf_k.size_loc_[0] * (i_loc + gf_k.size_loc_[1] * (j_loc + gf_k.size_loc_[2] * k_loc));
                        // (3-3) copy recvbuf.buf[] to gf_k.buf_[]
#ifdef PARTICLE_SIMULATOR_PMMM_RANGE_CHECK_GREEN_FUNC
                        assert(0 <= adr_dest && adr_dest < gf_k.size_loc_tot_);
                        assert(0 <= adr_src && adr_src < recvbuf.count_tot);
#endif
                        gf_k.buf_[adr_dest] = recvbuf.buf[adr_src++].c;
                    }
                }
                assert(adr_src == tail);
            }
        }
        // Free
        delete[] size_loc_req_tbl;
        delete[] start_loc_req_tbl;
        delete[] end_loc_req_tbl;
#else  // PARTICLE_SIMULATOR_MPI_PARALLEL
        for (S32 k = gf_k.start_loc_[3]; k <= gf_k.end_loc_[3]; k++)
            for (S32 j = gf_k.start_loc_[2]; j <= gf_k.end_loc_[2]; j++)
                for (S32 i = gf_k.start_loc_[1]; i <= gf_k.end_loc_[1]; i++)
                    for (S32 lm = gf_k.start_loc_[0]; lm <= gf_k.end_loc_[0]; lm++) {  // for(k,j,i,lm)
                        const S32 lm_dest = lm - gf_k.start_loc_[0];
                        const S32 i_dest = i - gf_k.start_loc_[1];
                        const S32 j_dest = j - gf_k.start_loc_[2];
                        const S32 k_dest = k - gf_k.start_loc_[3];
                        const S32 adr_dest = lm_dest + gf_k.size_loc_[0] * (i_dest + gf_k.size_loc_[1] * (j_dest + gf_k.size_loc_[2] * k_dest));

                        const S32 lm_src = lm - gf_k_tmp_.start_loc_[0];
                        const S32 i_src = i - gf_k_tmp_.start_loc_[1];
                        const S32 j_src = j - gf_k_tmp_.start_loc_[2];
                        const S32 k_src = k - gf_k_tmp_.start_loc_[3];
                        const S32 adr_src =
                            lm_src + gf_k_tmp_.size_loc_[0] * (i_src + gf_k_tmp_.size_loc_[1] * (j_src + gf_k_tmp_.size_loc_[2] * k_src));
#ifdef PARTICLE_SIMULATOR_PMMM_RANGE_CHECK_GREEN_FUNC
                        assert(0 <= adr_dest && adr_dest < gf_k.size_loc_tot_);
                        assert(0 <= adr_src && adr_src < gf_k_tmp_.size_loc_tot_);
#endif
                        gf_k.buf_[adr_dest] = gf_k_tmp_.buf_[adr_src];
                    }  // for(k,j,i,lm)
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
        // Free memory
        gf_k_tmp_.free();
    }
};

}  // namespace ParticleMeshMultipole
}  // namespace ParticleSimulator
