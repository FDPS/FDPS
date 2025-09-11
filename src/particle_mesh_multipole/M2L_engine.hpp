#pragma once
#include <complex>
#include <unordered_map>
#include "../ps_defs.hpp"
#include "particle_mesh_multipole_defs.hpp"
#include "green_function.hpp"

#define DEBUG_M2L_ENGINE_INITIALIZE
#define DEBUG_M2L_ENGINE_REDISTMM
#define DEBUG_M2L_ENGINE_CONVOLUTION
#define DEBUG_M2L_ENGINE_REDISTLE
#define NEW_INITIALIZE
namespace ParticleSimulator {
namespace ParticleMeshMultipole {

class M2L_Engine {
   public:
    using real_t = double;
    using cplx_t = std::complex<real_t>;

    Parameters param_;
    Parameters param_prev_;
    S32 n_thread_;
    S32 n_thread_prev_;
    bool first_call_of_convolution_;

    // Variables for communication
    CommInfo comm_info_;
    S32 mode_;
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
    U32 fftw_planning_rigor_flag_;
    S32 fft_size_crit_;

    S32 n_group_trans_;
    S32 n_group_0_trans_, n_group_1_trans_, n_group_2_trans_;
    std::vector<S32> local_0_start_trans_, local_0_end_trans_;  // transform; slowest varying index
    std::vector<S32> local_1_start_trans_, local_1_end_trans_;  // transform
    std::vector<S32> local_2_start_trans_, local_2_end_trans_;  // transform; fastest varying index
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

    // Variables for Green function
    GreenFunctionCalculator gf_calc_;
    Buffer<cplx_t> gf_k_trans_;

    // Variables for M2L and green function
    Buffer<real_t> mm_r_, le_r_;
    Buffer<cplx_t> mm_k_, le_k_;
    Buffer<cplx_t> mm_k_trans_, le_k_trans_;

    // Variables for FFT
    S32 fft_mode_;
    S32 fft_mode_prev_;
    S32 size_fft_[3];
    S32 size_fft_prev_[3];
    fftw_real_t *rbuf_;
    fftw_cplx_t *kbuf_;
    std::vector<fftw_plan> plan_fwd_;
    std::vector<fftw_plan> plan_bkw_;

    // Timing measurement
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

    F64 wtime_start_initialize_;
    F64 wtime_end_initialize_;
    F64 wtime_start_preproc_mm_r_comm_;
    F64 wtime_end_preproc_mm_r_comm_;
    F64 wtime_start_redist_mm_r_;
    F64 wtime_end_redist_mm_r_;
    F64 wtime_start_postproc_mm_r_comm_;
    F64 wtime_end_postproc_mm_r_comm_;
    F64 wtime_start_mm_r_to_mm_k_;
    F64 wtime_end_mm_r_to_mm_k_;
    F64 wtime_start_gather_mm_k_trans_;
    F64 wtime_end_gather_mm_k_trans_;
    F64 wtime_start_transform_;
    F64 wtime_end_transform_;
    F64 wtime_start_scatter_le_k_trans_;
    F64 wtime_end_scatter_le_k_trans_;
    F64 wtime_start_le_k_to_le_r_;
    F64 wtime_end_le_k_to_le_r_;
    F64 wtime_start_preproc_le_r_comm_;
    F64 wtime_end_preproc_le_r_comm_;
    F64 wtime_start_redist_le_r_;
    F64 wtime_end_redist_le_r_;
    F64 wtime_start_postproc_le_r_comm_;
    F64 wtime_end_postproc_le_r_comm_;

    // Timing analysis (detailed information)
    F64 wtime_start_initialize__copy_;
    F64 wtime_end_initialize__copy_;
    F64 wtime_start_initialize__main_;
    F64 wtime_end_initialize__main_;

    M2L_Engine() {
        n_thread_ = -1;
        n_thread_prev_ = -1;
        first_call_of_convolution_ = true;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        group_all_ = MPI_GROUP_NULL;
        group_fft_ = MPI_GROUP_NULL;
        group_int_ = MPI_GROUP_NULL;
        comm_all_ = MPI_COMM_NULL;
        comm_fft_ = MPI_COMM_NULL;
        comm_int_ = MPI_COMM_NULL;
#endif
        fft_mode_ = -1;
        fft_mode_prev_ = -1;
        for (S32 i = 0; i < 3; i++) {
            size_fft_[i] = 0;
            size_fft_prev_[i] = 0;
        }
        rbuf_ = nullptr;
        kbuf_ = nullptr;
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
        if (rbuf_ != nullptr) fftw_free(rbuf_);
        if (kbuf_ != nullptr) fftw_free(kbuf_);
        for (size_t i = 0; i < plan_fwd_.size(); i++) fftw_destroy_plan(plan_fwd_[i]);
        for (size_t i = 0; i < plan_bkw_.size(); i++) fftw_destroy_plan(plan_bkw_[i]);
    }

#if defined(NEW_INITIALIZE)
    // T is real_t or cplx_t
    // if is_mm_r == true, it means that the buffer is mm_r_, othewise le_k_
    template <typename T>
    void setSizeInfo(Buffer<T> &buf, const S32 n_mm_compo, const S32vec &n_cell, const bool use_mpifft_, const bool is_mm_r) {
        // グローバル開始・終了インデックスの設定
        buf.start_glb_[0] = 0;
        buf.start_glb_[1] = 0;
        buf.start_glb_[2] = 0;
        buf.start_glb_[3] = 0;
        const S32 bc = param_.bc;
        if (bc == BOUNDARY_CONDITION_OPEN) {
            buf.end_glb_[0] = n_mm_compo - 1;
            buf.end_glb_[1] = is_mm_r ? (2 * n_cell.x - 1) : (1 + n_cell.x - 1);
            buf.end_glb_[2] = 2 * n_cell.y - 1;
            buf.end_glb_[3] = 2 * n_cell.z - 1;
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
            buf.end_glb_[0] = n_mm_compo - 1;
            buf.end_glb_[1] = is_mm_r ? (n_cell.x - 1) : (1 + n_cell.x / 2 - 1);
            buf.end_glb_[2] = n_cell.y - 1;
            buf.end_glb_[3] = n_cell.z - 1;
        } else {
            assert(false);
        }
        buf.calcSizeGlb();

        // ローカル開始・終了インデックスの設定
        if (rank_in_my_group_ < n_proc_for_fft_) {
            buf.start_loc_[0] = lm_start_[idx_to_my_group_];
            buf.start_loc_[1] = 0;
            buf.start_loc_[2] = 0;
            buf.start_loc_[3] = (use_mpifft_ && !is_mm_r) ? local_1_start_[rank_in_my_group_] : local_0_start_[rank_in_my_group_];

            buf.end_loc_[0] = lm_end_[idx_to_my_group_];
            if (bc == BOUNDARY_CONDITION_OPEN) {
                buf.end_loc_[1] = is_mm_r ? (2 * n_cell.x - 1) : n_cell.x;
                buf.end_loc_[2] = is_mm_r ? (2 * n_cell.y - 1) : (use_mpifft_ ? (2 * n_cell.z - 1) : (2 * n_cell.y - 1));
            } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                buf.end_loc_[1] = is_mm_r ? (n_cell.x - 1) : (n_cell.x / 2);
                buf.end_loc_[2] = n_cell.y - 1;
            }
            buf.end_loc_[3] = (use_mpifft_ && !is_mm_r) ? local_1_end_[rank_in_my_group_] : local_0_end_[rank_in_my_group_];
            buf.calcSizeLoc();
        }
    }

    // local_?_start_trans_ と local_?_end_trans_ を設定する関数

    void setLocalTransIndices() {
        if (rank_in_my_group_ < n_proc_for_fft_) {
            initializeLocalTransIndices();

            if (n_group_ > 1) {
                calculateGroupTransSizes();
                distributeTransIndices();
            } else {
                n_group_trans_ = n_group_;
                for (S32 i = 0; i < n_group_; i++) {
                    local_0_start_trans_[i] = le_k_.start_loc_[3];
                    local_0_end_trans_[i] = le_k_.end_loc_[3];
                    local_1_start_trans_[i] = le_k_.start_loc_[2];
                    local_1_end_trans_[i] = le_k_.end_loc_[2];
                    local_2_start_trans_[i] = le_k_.start_loc_[1];
                    local_2_end_trans_[i] = le_k_.end_loc_[1];
                }
            }
        }
    }

    void initializeLocalTransIndices() {
        local_0_start_trans_.resize(n_group_);
        local_0_end_trans_.resize(n_group_);
        local_1_start_trans_.resize(n_group_);
        local_1_end_trans_.resize(n_group_);
        local_2_start_trans_.resize(n_group_);
        local_2_end_trans_.resize(n_group_);
        n_group_0_trans_ = 1;
        n_group_1_trans_ = 1;
        n_group_2_trans_ = 1;
    }

    void calculateGroupTransSizes() {
        S32 tmp = n_group_, ret;
        ret = getMaximumFactorLessThanOrEqualTo(tmp, le_k_.size_loc_[3]);
        if (ret != -1) {
            n_group_0_trans_ = ret;
            tmp /= n_group_0_trans_;
        }
        if (tmp > 1) {
            ret = getMaximumFactorLessThanOrEqualTo(tmp, le_k_.size_loc_[2]);
            if (ret != -1) {
                n_group_1_trans_ = ret;
                tmp /= n_group_1_trans_;
            }
        }
        if (tmp > 1) {
            ret = getMaximumFactorLessThanOrEqualTo(tmp, le_k_.size_loc_[1]);
            if (ret != -1) {
                n_group_2_trans_ = ret;
                tmp /= n_group_2_trans_;
            }
        }
        n_group_trans_ = n_group_0_trans_ * n_group_1_trans_ * n_group_2_trans_;
    }

    void distributeTransIndices() {
        for (S32 i = 0; i < n_group_; i++) {
            if (i < n_group_trans_) {
                S32 idx[3];
                calculateIndices(i, idx);

                setLocalTransIndexForDimension(idx[0], n_group_0_trans_, le_k_.size_loc_[3], le_k_.start_loc_[3], local_0_start_trans_[i],
                                               local_0_end_trans_[i]);
                setLocalTransIndexForDimension(idx[1], n_group_1_trans_, le_k_.size_loc_[2], le_k_.start_loc_[2], local_1_start_trans_[i],
                                               local_1_end_trans_[i]);
                setLocalTransIndexForDimension(idx[2], n_group_2_trans_, le_k_.size_loc_[1], le_k_.start_loc_[1], local_2_start_trans_[i],
                                               local_2_end_trans_[i]);
            } else {
                // 計算に参加しない場合
                local_0_start_trans_[i] = -1;
                local_0_end_trans_[i] = -1;
                local_1_start_trans_[i] = -1;
                local_1_end_trans_[i] = -1;
                local_2_start_trans_[i] = -1;
                local_2_end_trans_[i] = -1;
            }
        }
    }

    void calculateIndices(S32 idx_flat, S32 idx[3]) {
        S32 tmp = idx_flat;
        idx[0] = tmp / (n_group_1_trans_ * n_group_2_trans_);
        tmp -= idx[0] * (n_group_1_trans_ * n_group_2_trans_);
        idx[1] = tmp / n_group_2_trans_;
        idx[2] = tmp % n_group_2_trans_;
    }

    void setLocalTransIndexForDimension(S32 idx, S32 n_group_trans, S32 size, S32 start_loc, S32 &local_start, S32 &local_end) {
        const S32 n_slab_min = size / n_group_trans;
        const S32 n_rem = size % n_group_trans;
        S32 start = start_loc;

        for (S32 k = 0; k < idx; k++) {
            S32 n_slab = n_slab_min + (k < n_rem ? 1 : 0);
            start += n_slab;
        }
        S32 n_slab = n_slab_min + (idx < n_rem ? 1 : 0);
        local_start = start;
        local_end = start + n_slab - 1;
    }

    void setSizeInfoForTransBuffers() {
        mm_k_trans_.copySizeInfoGlbOnlyFrom(mm_k_);
        le_k_trans_.copySizeInfoGlbOnlyFrom(le_k_);
        const auto n_mm_compo = (param_.p + 1) * (param_.p + 1);
        if (rank_in_my_group_ < n_proc_for_fft_ && idx_to_my_group_ < n_group_trans_) {
            mm_k_trans_.copySizeInfoFrom(mm_k_);
            mm_k_trans_.start_loc_[0] = 0;
            mm_k_trans_.start_loc_[1] = local_2_start_trans_[idx_to_my_group_];
            mm_k_trans_.start_loc_[2] = local_1_start_trans_[idx_to_my_group_];
            mm_k_trans_.start_loc_[3] = local_0_start_trans_[idx_to_my_group_];
            mm_k_trans_.end_loc_[0] = n_mm_compo - 1;
            mm_k_trans_.end_loc_[1] = local_2_end_trans_[idx_to_my_group_];
            mm_k_trans_.end_loc_[2] = local_1_end_trans_[idx_to_my_group_];
            mm_k_trans_.end_loc_[3] = local_0_end_trans_[idx_to_my_group_];
            mm_k_trans_.calcSizeLoc();
            mm_k_trans_.calcSizeTot();

            le_k_trans_.copySizeInfoFrom(mm_k_trans_);
        }
    }

    /*
        void determineOMPFFTUsage() {
            if (fft_size_crit_ < 0) {
                use_ompfft_ = determineWhetherToUseOMPFFT();
            } else {
                const S32 fft_size = n_cell.x * n_cell.y * n_cell.z;
                use_ompfft_ = (fft_size >= fft_size_crit_);
            }
        }
    */

    void initialize_multi_components_on_one_proc_mode() {
        // n_proc < n_mm_compo
        // auto n_proc = Comm::getNumberOfProc();
        auto my_rank = Comm::getRank();
        auto n_mm_compo = (param_.p + 1) * (param_.p + 1);
        auto bc = param_.bc;
        auto n_cell = param_.n_cell;
        mode_ = 2;
        // n_group_ = n_proc;
        n_proc_min_in_group_ = 1;
        n_proc_in_my_group_ = 1;
        rank_in_my_group_ = 0;
        n_proc_for_fft_ = 1;
        idx_to_my_group_ = my_rank;
        is_mpifft_usable_ = false;
        use_mpifft_ = false;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
        if (!FFTW_Status::is_fftw_init_threads_called()) {
            fftw_init_threads();
            FFTW_Status::is_fftw_init_threads_called() = true;
        }
#endif
        // Set rank_start_, rank_end_
        // rank_start_.resize(n_group_);
        // rank_end_.resize(n_group_);
        for (S32 i = 0; i < n_group_; i++) {
            rank_start_[i] = i;
            rank_end_[i] = i;
        }
        // Set lm_start_, lm_end_

#if 1
        for (S32 i = 0; i < n_group_; i++) {
            CalcAdrToSplitData(lm_start_[i], lm_end_[i], i, n_group_, n_mm_compo);
            lm_end_[i]--;
        }
        assert(lm_end_[n_group_ - 1] == n_mm_compo - 1);
#else
        const S32 n_mm_compo_per_proc_min = n_mm_compo / n_group_;
        S32 n_rem = n_mm_compo - n_mm_compo_per_proc_min * n_group_;
        // lm_start_.resize(n_group_);
        // lm_end_.resize(n_group_);
        lm_start_[0] = 0;
        for (S32 i = 0; i < n_group_; i++) {
            S32 n_mm_compo_per_proc = n_mm_compo_per_proc_min;
            if (i < n_rem) n_mm_compo_per_proc++;
            lm_end_[i] = lm_start_[i] + (n_mm_compo_per_proc - 1);
            if (i < n_group_ - 1) lm_start_[i + 1] = lm_end_[i] + 1;
        }
#endif
        // Set local_0_start_, local_0_end_
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        group_fft_ = MPI_GROUP_NULL;
        comm_fft_ = MPI_COMM_NULL;
#endif
        local_0_start_.resize(n_proc_for_fft_);
        local_0_end_.resize(n_proc_for_fft_);
        if (bc == BOUNDARY_CONDITION_OPEN) {
            local_0_start_[0] = 0;
            local_0_end_[0] = 2 * n_cell.z - 1;
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
            local_0_start_[0] = 0;
            local_0_end_[0] = n_cell.z - 1;
        } else {
            assert(false);
        }
        // Note that we do not need to set local_1_start_
        // and local_1_end_ because of use_mpifft = false.
        // (hence, both arrays are not accessed)
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        // Make group_int, comm_int
        group_int_ = parent_group_;
        comm_int_ = parent_comm_;
#endif
    }

    void initialize_single_component_on_multi_procs_mode() {
        // auto n_mm_compo = (param_.p + 1) * (param_.p + 1);
        mode_ = 1;
        // // n_group_ = n_mm_compo;
    }

    void initialize(const Parameters &param, const U32 fftw_planning_rigor_flag, const bool use_mpifft_if_possible, const S32 fft_size_crit,
                    const CommInfo &comm_info, const S32 n_thread, const bool debug_flag = false) {
        F64 wtime_tmp = GetWtime();
        wtime_start_initialize_ = wtime_tmp;
        // Copy parameters, etc.
        wtime_start_initialize__copy_ = wtime_tmp;
        param_prev_ = param_;  // save the previous parameters
        param_ = param;
        fftw_planning_rigor_flag_ = fftw_planning_rigor_flag;
        fft_size_crit_ = fft_size_crit;
        comm_info_ = comm_info;
        n_thread_prev_ = n_thread_;
        n_thread_ = n_thread;
        wtime_tmp = GetWtime();
        wtime_end_initialize__copy_ = wtime_tmp;
        wtime_start_initialize__main_ = wtime_tmp;
        // (Re-)Initialize
        // comm_info_.barrier();
        // if (comm_info_.getRank() == 0) std::cout << "ok1 @ initialize() [M2L_Engine]" << std::endl;
        if ((param_ != param_prev_) || (n_thread_ != n_thread_prev_) || (use_mpifft_ != use_mpifft_if_possible)) {
            // comm_info_.barrier();
            // if (comm_info_.getRank() == 0) std::cout << "ok2 @ initialize() [M2L_Engine]" << std::endl;
            //  [TODO] Support the case that comm_info is changed from the previous call
            //         of this function. To support this, we need to investigate how to
            //         know if two instances (current, previous) of CommInfo class
            //         represent the same MPI_Comm. In other words, we must know if it
            //         is allowed to create instances of MPI_Comm class with different
            //         values from single set of processes.
            // PS::DebugUtils::fout() << "M2L initialize performed." << std::endl;
            const S32 n_mm_compo = (param_.p + 1) * (param_.p + 1);
            const S32vec n_cell = param_.n_cell;
            const S32 bc = param_.bc;
            const S32 n_proc = comm_info_.getNumberOfProc();
            n_group_ = std::min(n_proc, n_mm_compo);
            rank_start_.resize(n_group_);
            rank_end_.resize(n_group_);
            lm_start_.resize(n_group_);
            lm_end_.resize(n_group_);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            // Set parent_comm
            parent_comm_ = comm_info_.getCommunicator();
            // Get information of parent_comm
            // S32 n_proc, my_rank;
            S32 my_rank;
            // MPI_Comm_size(parent_comm_, &n_proc);
            MPI_Comm_rank(parent_comm_, &my_rank);
            MPI_Comm_group(parent_comm_, &parent_group_);
            n_proc_in_parent_group_ = n_proc;
            rank_in_parent_group_ = my_rank;
            // Get MPI_Datatype for ptrdiff_t, real_t, cplx_t
            mpi_real_t_ = GetDataType<real_t>();
            mpi_cplx_t_ = GetDataType<cplx_t>();
            // Divide parent_comm into groups
            if (n_proc >= n_mm_compo) {
                initialize_single_component_on_multi_procs_mode();
                // Set rank_start_, rank_end_
                n_proc_min_in_group_ = n_proc / n_group_;
#if 1
                for (S32 i = 0; i < n_group_; i++) {
                    CalcAdrToSplitData(rank_start_[i], rank_end_[i], i, n_group_, n_proc);
                    rank_end_[i]--;  // rank_end_ is inclusive
                }
                assert(rank_end_[n_group_ - 1] == n_proc - 1);
#else
                // new implementation
                S32 n_rem = n_proc - n_proc_min_in_group_ * n_group_;
                // rank_start_.resize(n_group_);
                // rank_end_.resize(n_group_);
                rank_start_[0] = 0;
                for (S32 i = 0; i < n_group_; i++) {
                    S32 n_proc_in_group = n_proc_min_in_group_;
                    if (i < n_rem) n_proc_in_group++;
                    rank_end_[i] = rank_start_[i] + (n_proc_in_group - 1);
                    if (i < n_group_ - 1) rank_start_[i + 1] = rank_end_[i] + 1;
                }
#endif
                if (debug_flag) {
                    if (my_rank == 0) {
                        std::cout << "### Information of groups for M2L ###" << std::endl;
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
                // lm_start_.resize(n_group_);
                // lm_end_.resize(n_group_);
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
#if 1
                // n_proc_max_for_fft_ = getMaxProcForFFT(bc, n_cell);
                n_proc_max_for_fft = (bc == BOUNDARY_CONDITION_OPEN) ? n_cell.z
                                     : (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
                                         ? n_cell.z / 2
                                         : (assert(false), 0);  // Unreachable, but required to avoid compiler warning
                n_proc_for_fft_ = std::min(n_proc_min_in_group_, n_proc_max_for_fft);
                // Make group_fft, comm_fft
                is_mpifft_usable_ = (n_proc_for_fft_ > 1);
                use_mpifft_ = is_mpifft_usable_ ? use_mpifft_if_possible : false;
#else
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
#endif

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
#if 1
                        if (bc != BOUNDARY_CONDITION_OPEN && bc != BOUNDARY_CONDITION_PERIODIC_XYZ) {
                            if (my_rank == 0) std::cout << "This boundary condition is not supported yet." << std::endl;
                            Abort(-1);
                        }
                        auto n_cell_z_tmp = (bc == BOUNDARY_CONDITION_OPEN) ? 2 * n_cell.z : n_cell.z;
                        auto n_cell_y_tmp = (bc == BOUNDARY_CONDITION_OPEN) ? 2 * n_cell.y : n_cell.y;
                        auto n_cell_x_tmp = (bc == BOUNDARY_CONDITION_OPEN) ? n_cell.x + 1 : n_cell.x / 2 + 1;
                        // alloc_local = fftw_mpi_local_size_3d_transposed(n_cell_z_tmp, n_cell_y_tmp, n_cell_x_tmp, comm_fft_, &local_n0,
                        //                                                 &local_0_start, &local_n1, &local_1_start);
                        fftw_mpi_local_size_3d_transposed(n_cell_z_tmp, n_cell_y_tmp, n_cell_x_tmp, comm_fft_, &local_n0, &local_0_start, &local_n1,
                                                          &local_1_start);
#else
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
#endif
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
                } else {                  // Case of use_mpifft_ == false
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
                    } else {
                        assert(false);
                    }
                }  // End of the if-branch due to the value of use_mpifft_
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
                // n_proc < n_mm_compo
                initialize_multi_components_on_one_proc_mode();
            }
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
        // Sequential execution/ OpenMP parallelization
            n_proc_in_parent_group_ = 1;
            rank_in_parent_group_ = 0;
            initialize_multi_components_on_one_proc_mode();
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
        // Set the size information of mm_r, le_r, mm_k, le_k
        // set the size information of mm_r_
            setSizeInfo(mm_r_, n_mm_compo, n_cell, use_mpifft_, true);
            // le_k_ のサイズ情報を設定
            // set the size information of le_k_
            setSizeInfo(le_k_, n_mm_compo, n_cell, use_mpifft_, false);
            // Set the size information of le_r_
            le_r_.copySizeInfoFrom(mm_r_);
            // Set the size information of mm_k_
            mm_k_.copySizeInfoFrom(le_k_);
            // Calculate total sizes
            mm_r_.calcSizeTot();
            mm_k_.calcSizeTot();
            le_k_.calcSizeTot();
            le_r_.calcSizeTot();

            // Set local_?_start_trans_, local_?_end_trans_
            setLocalTransIndices();
            setSizeInfoForTransBuffers();

            if (fft_size_crit_ < 0) {
                use_ompfft_ = determineWhetherToUseOMPFFT();
            } else {
                const S32 fft_size = n_cell.x * n_cell.y * n_cell.z;
                use_ompfft_ = (fft_size >= fft_size_crit_);
            }

        } else {
            // PS::DebugUtils::fout() << "M2L initialize skipped." << std::endl;
        }
        wtime_tmp = GetWtime();
        wtime_end_initialize__main_ = wtime_tmp;
        wtime_end_initialize_ = wtime_tmp;
        time_profile_.PMMM__M2L__initialize += (wtime_end_initialize_ - wtime_start_initialize_);
    }
#else
    // ORIGINAL IMPLEMENTATION
    void initialize(const Parameters &param, const U32 fftw_planning_rigor_flag, const bool use_mpifft_if_possible, const S32 fft_size_crit,
                    const CommInfo &comm_info, const S32 n_thread, const bool debug_flag = false) {
        F64 wtime_tmp = GetWtime();
        wtime_start_initialize_ = wtime_tmp;
        // Copy parameters, etc.
        wtime_start_initialize__copy_ = wtime_tmp;
        param_prev_ = param_;  // save the previous parameters
        param_ = param;
        fftw_planning_rigor_flag_ = fftw_planning_rigor_flag;
        fft_size_crit_ = fft_size_crit;
        comm_info_ = comm_info;
        n_thread_prev_ = n_thread_;
        n_thread_ = n_thread;
        wtime_tmp = GetWtime();
        wtime_end_initialize__copy_ = wtime_tmp;
        wtime_start_initialize__main_ = wtime_tmp;
        // (Re-)Initialize
        // comm_info_.barrier();
        // if (comm_info_.getRank() == 0) std::cout << "ok1 @ initialize() [M2L_Engine]" << std::endl;
        if ((param_ != param_prev_) || (n_thread_ != n_thread_prev_) || (use_mpifft_ != use_mpifft_if_possible)) {
            // comm_info_.barrier();
            // if (comm_info_.getRank() == 0) std::cout << "ok2 @ initialize() [M2L_Engine]" << std::endl;
            //  [TODO] Support the case that comm_info is changed from the previous call
            //         of this function. To support this, we need to investigate how to
            //         know if two instances (current, previous) of CommInfo class
            //         represent the same MPI_Comm. In other words, we must know if it
            //         is allowed to create instances of MPI_Comm class with different
            //         values from single set of processes.
            // PS::DebugUtils::fout() << "M2L initialize performed." << std::endl;
            const S32 n_mm_compo = (param_.p + 1) * (param_.p + 1);
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
                mode_ = 1;
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
                        std::cout << "### Information of groups for M2L ###" << std::endl;
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
                } else {                  // Case of use_mpifft_ == false
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
                    } else {
                        assert(false);
                    }
                }  // End of the if-branch due to the value of use_mpifft_
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
                mode_ = 2;
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
                } else {
                    assert(false);
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
            mode_ = 0;
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
            } else {
                assert(false);
            }
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL

            // Set the size information of mm_r, le_r, mm_k, le_k
            if (bc == BOUNDARY_CONDITION_OPEN) {
                // Set the size information of mm_r
                mm_r_.start_glb_[0] = 0;
                mm_r_.start_glb_[1] = 0;
                mm_r_.start_glb_[2] = 0;
                mm_r_.start_glb_[3] = 0;
                mm_r_.end_glb_[0] = n_mm_compo - 1;
                mm_r_.end_glb_[1] = 2 * n_cell.x - 1;
                mm_r_.end_glb_[2] = 2 * n_cell.y - 1;
                mm_r_.end_glb_[3] = 2 * n_cell.z - 1;
                mm_r_.calcSizeGlb();
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    mm_r_.start_loc_[0] = lm_start_[idx_to_my_group_];
                    mm_r_.start_loc_[1] = 0;
                    mm_r_.start_loc_[2] = 0;
                    mm_r_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                    mm_r_.end_loc_[0] = lm_end_[idx_to_my_group_];
                    mm_r_.end_loc_[1] = 2 * n_cell.x - 1;
                    mm_r_.end_loc_[2] = 2 * n_cell.y - 1;
                    mm_r_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                    mm_r_.calcSizeLoc();
                }
                // Set the size information of le_k
                if (use_mpifft_) {
                    // In this case, we set the size information
                    // basend on FFTW_MPI_TRANSPOSED_OUT format.
                    le_k_.start_glb_[0] = 0;
                    le_k_.start_glb_[1] = 0;
                    le_k_.start_glb_[2] = 0;
                    le_k_.start_glb_[3] = 0;
                    le_k_.end_glb_[0] = n_mm_compo - 1;
                    le_k_.end_glb_[1] = 1 + n_cell.x - 1;
                    le_k_.end_glb_[2] = 2 * n_cell.z - 1;
                    le_k_.end_glb_[3] = 2 * n_cell.y - 1;
                    le_k_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        le_k_.start_loc_[0] = lm_start_[idx_to_my_group_];  // only one (l,m)
                        le_k_.start_loc_[1] = 0;
                        le_k_.start_loc_[2] = 0;
                        le_k_.start_loc_[3] = local_1_start_[rank_in_my_group_];
                        le_k_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        le_k_.end_loc_[1] = n_cell.x;
                        le_k_.end_loc_[2] = 2 * n_cell.z - 1;
                        le_k_.end_loc_[3] = local_1_end_[rank_in_my_group_];
                        le_k_.calcSizeLoc();
                    }
                } else {
                    // In this case, we set the size information
                    // in the same way as le_r.
                    le_k_.start_glb_[0] = 0;
                    le_k_.start_glb_[1] = 0;
                    le_k_.start_glb_[2] = 0;
                    le_k_.start_glb_[3] = 0;
                    le_k_.end_glb_[0] = n_mm_compo - 1;
                    le_k_.end_glb_[1] = 1 + n_cell.x - 1;
                    le_k_.end_glb_[2] = 2 * n_cell.y - 1;
                    le_k_.end_glb_[3] = 2 * n_cell.z - 1;
                    le_k_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        le_k_.start_loc_[0] = lm_start_[idx_to_my_group_];
                        le_k_.start_loc_[1] = 0;
                        le_k_.start_loc_[2] = 0;
                        le_k_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                        le_k_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        le_k_.end_loc_[1] = n_cell.x;
                        le_k_.end_loc_[2] = 2 * n_cell.y - 1;
                        le_k_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                        le_k_.calcSizeLoc();
                    }
                }
            } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                // Set the size information of mm_r
                mm_r_.start_glb_[0] = 0;
                mm_r_.start_glb_[1] = 0;
                mm_r_.start_glb_[2] = 0;
                mm_r_.start_glb_[3] = 0;
                mm_r_.end_glb_[0] = n_mm_compo - 1;
                mm_r_.end_glb_[1] = n_cell.x - 1;
                mm_r_.end_glb_[2] = n_cell.y - 1;
                mm_r_.end_glb_[3] = n_cell.z - 1;
                mm_r_.calcSizeGlb();
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    mm_r_.start_loc_[0] = lm_start_[idx_to_my_group_];
                    mm_r_.start_loc_[1] = 0;
                    mm_r_.start_loc_[2] = 0;
                    mm_r_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                    mm_r_.end_loc_[0] = lm_end_[idx_to_my_group_];
                    mm_r_.end_loc_[1] = n_cell.x - 1;
                    mm_r_.end_loc_[2] = n_cell.y - 1;
                    mm_r_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                    mm_r_.calcSizeLoc();
                }
                // Set the size information of le_k
                if (use_mpifft_) {
                    // In this case, we set the size information
                    // basend on FFTW_MPI_TRANSPOSED_OUT format.
                    le_k_.start_glb_[0] = 0;
                    le_k_.start_glb_[1] = 0;
                    le_k_.start_glb_[2] = 0;
                    le_k_.start_glb_[3] = 0;
                    le_k_.end_glb_[0] = n_mm_compo - 1;
                    le_k_.end_glb_[1] = 1 + n_cell.x / 2 - 1;
                    le_k_.end_glb_[2] = n_cell.z - 1;
                    le_k_.end_glb_[3] = n_cell.y - 1;
                    le_k_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        le_k_.start_loc_[0] = lm_start_[idx_to_my_group_];  // only one (l,m)
                        le_k_.start_loc_[1] = 0;
                        le_k_.start_loc_[2] = 0;
                        le_k_.start_loc_[3] = local_1_start_[rank_in_my_group_];
                        le_k_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        le_k_.end_loc_[1] = n_cell.x / 2;
                        le_k_.end_loc_[2] = n_cell.z - 1;
                        le_k_.end_loc_[3] = local_1_end_[rank_in_my_group_];
                        le_k_.calcSizeLoc();
                    }
                } else {
                    // In this case, we set the size information
                    // in the same way as le_r.
                    le_k_.start_glb_[0] = 0;
                    le_k_.start_glb_[1] = 0;
                    le_k_.start_glb_[2] = 0;
                    le_k_.start_glb_[3] = 0;
                    le_k_.end_glb_[0] = n_mm_compo - 1;
                    le_k_.end_glb_[1] = 1 + n_cell.x / 2 - 1;
                    le_k_.end_glb_[2] = n_cell.y - 1;
                    le_k_.end_glb_[3] = n_cell.z - 1;
                    le_k_.calcSizeGlb();
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        le_k_.start_loc_[0] = lm_start_[idx_to_my_group_];
                        le_k_.start_loc_[1] = 0;
                        le_k_.start_loc_[2] = 0;
                        le_k_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                        le_k_.end_loc_[0] = lm_end_[idx_to_my_group_];
                        le_k_.end_loc_[1] = n_cell.x / 2;
                        le_k_.end_loc_[2] = n_cell.y - 1;
                        le_k_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                        le_k_.calcSizeLoc();
                    }
                }
            } else {
                assert(false);
            }
            // Set the size information of le_r_
            le_r_.copySizeInfoFrom(mm_r_);
            // Set the size information of mm_k_
            mm_k_.copySizeInfoFrom(le_k_);
            // Calculate total sizes
            mm_r_.calcSizeTot();
            mm_k_.calcSizeTot();
            le_k_.calcSizeTot();
            le_r_.calcSizeTot();

            // Set local_?_start_trans_, local_?_end_trans_
            if (rank_in_my_group_ < n_proc_for_fft_) {
                local_0_start_trans_.resize(n_group_);
                local_0_end_trans_.resize(n_group_);
                local_1_start_trans_.resize(n_group_);
                local_1_end_trans_.resize(n_group_);
                local_2_start_trans_.resize(n_group_);
                local_2_end_trans_.resize(n_group_);
                n_group_0_trans_ = 1;  // initialize
                n_group_1_trans_ = 1;
                n_group_2_trans_ = 1;
                if (n_group_ > 1) {
                    S32 tmp = n_group_, ret;
                    ret = getMaximumFactorLessThanOrEqualTo(tmp, le_k_.size_loc_[3]);
                    if (ret != -1) {
                        n_group_0_trans_ = ret;
                        tmp /= n_group_0_trans_;
                    }
                    if (tmp > 1) {
                        ret = getMaximumFactorLessThanOrEqualTo(tmp, le_k_.size_loc_[2]);
                        if (ret != -1) {
                            n_group_1_trans_ = ret;
                            tmp /= n_group_1_trans_;
                        }
                    }
                    if (tmp > 1) {
                        ret = getMaximumFactorLessThanOrEqualTo(tmp, le_k_.size_loc_[1]);
                        if (ret != -1) {
                            n_group_2_trans_ = ret;
                            tmp /= n_group_2_trans_;
                        }
                    }
                    n_group_trans_ = n_group_0_trans_ * n_group_1_trans_ * n_group_2_trans_;
#if 0
                            // Check
                            if (idx_to_my_group_ == 0) {
                                std::cout << "my_rank       = " << rank_in_parent_group_ << std::endl;
                                std::cout << "n_group_0_trans_ = " << n_group_0_trans_ 
                                          << " (n_cell_0_ = " << le_k_.size_loc_[3] << ")"
                                          << std::endl;
                                std::cout << "n_group_1_trans_ = " << n_group_1_trans_
                                          << " (n_cell_1_ = " << le_k_.size_loc_[2] << ")"
                                          << std::endl;
                                std::cout << "n_group_2_trans_ = " << n_group_2_trans_
                                          << " (n_cell_2_ = " << le_k_.size_loc_[1] << ")"
                                          << std::endl;
                                std::cout << "n_group_trans_   = " << n_group_trans_   << std::endl;
                            }
#endif
                    for (S32 i = 0; i < n_group_; i++) {
                        if (i < n_group_trans_) {
                            // In this case, rank i in comm_int is involved
                            // in the calculation of transform().
                            // (1) Calculate idx
                            S32 tmp = i, idx[3];
                            idx[0] = tmp / (n_group_1_trans_ * n_group_2_trans_);
                            tmp -= idx[0] * (n_group_1_trans_ * n_group_2_trans_);
                            idx[1] = tmp / n_group_2_trans_;
                            tmp -= idx[1] * n_group_2_trans_;
                            idx[2] = tmp;
#if 0
                                    if (idx_to_my_group_ == 0)
                                        std::cout << "idx: " << idx[0] << " " << idx[1] << " " << idx[2] << std::endl;
#endif
                            // (2) local_0_start_trans_, local_0_end_trans_
                            {
                                const S32 size = le_k_.size_loc_[3];
                                const S32 n_slab_min = size / n_group_0_trans_;
                                const S32 n_rem = size - n_slab_min * n_group_0_trans_;
                                S32 start = le_k_.start_loc_[3];
                                for (S32 k = 0; k < idx[0]; k++) {
                                    S32 n_slab = n_slab_min;
                                    if (k < n_rem) n_slab++;
                                    start += n_slab;
                                }
                                S32 n_slab = n_slab_min;
                                if (idx[0] < n_rem) n_slab++;
                                local_0_start_trans_[i] = start;
                                local_0_end_trans_[i] = start + (n_slab - 1);
                            }
                            // (3) local_1_start_trans_, local_1_end_trans_
                            {
                                const S32 size = le_k_.size_loc_[2];
                                const S32 n_slab_min = size / n_group_1_trans_;
                                const S32 n_rem = size - n_slab_min * n_group_1_trans_;
                                S32 start = le_k_.start_loc_[2];
                                for (S32 k = 0; k < idx[1]; k++) {
                                    S32 n_slab = n_slab_min;
                                    if (k < n_rem) n_slab++;
                                    start += n_slab;
                                }
                                S32 n_slab = n_slab_min;
                                if (idx[1] < n_rem) n_slab++;
                                local_1_start_trans_[i] = start;
                                local_1_end_trans_[i] = start + (n_slab - 1);
                            }
                            // (3) local_2_start_trans_, local_2_end_trans_
                            {
                                const S32 size = le_k_.size_loc_[1];
                                const S32 n_slab_min = size / n_group_2_trans_;
                                const S32 n_rem = size - n_slab_min * n_group_2_trans_;
                                S32 start = le_k_.start_loc_[1];
                                for (S32 k = 0; k < idx[2]; k++) {
                                    S32 n_slab = n_slab_min;
                                    if (k < n_rem) n_slab++;
                                    start += n_slab;
                                }
                                S32 n_slab = n_slab_min;
                                if (idx[2] < n_rem) n_slab++;
                                local_2_start_trans_[i] = start;
                                local_2_end_trans_[i] = start + (n_slab - 1);
                            }
                        } else {
                            // In this case, rank i in comm_int is not involved
                            // in the calculation of transform().
                            local_0_start_trans_[i] = -1;
                            local_0_end_trans_[i] = -1;
                            local_1_start_trans_[i] = -1;
                            local_1_end_trans_[i] = -1;
                            local_2_start_trans_[i] = -1;
                            local_2_end_trans_[i] = -1;
                        }
                    }
                } else {
                    n_group_trans_ = n_group_;
                    for (S32 i = 0; i < n_group_; i++) {
                        local_0_start_trans_[i] = le_k_.start_loc_[3];
                        local_0_end_trans_[i] = le_k_.end_loc_[3];
                        local_1_start_trans_[i] = le_k_.start_loc_[2];
                        local_1_end_trans_[i] = le_k_.end_loc_[2];
                        local_2_start_trans_[i] = le_k_.start_loc_[1];
                        local_2_end_trans_[i] = le_k_.end_loc_[1];
                    }
                }
#if 0
                        // Check
                        if (idx_to_my_group_ == 0 && rank_in_my_group_ == 0) {
                            for (S32 i = 0; i < n_group_trans_; i++) {
                                std::cout << "i = " << i << ": "
                                          << local_0_start_trans_[i] << " "
                                          << local_0_end_trans_[i] << " "
                                          << local_1_start_trans_[i] << " "
                                          << local_1_end_trans_[i] << " "
                                          << local_2_start_trans_[i] << " "
                                          << local_2_end_trans_[i] << " "
                                          << std::endl;
                            }
                        }
#endif
            }
            // Set the size information of mm_k_trans_, le_k_trans_
            mm_k_trans_.copySizeInfoGlbOnlyFrom(mm_k_);
            le_k_trans_.copySizeInfoGlbOnlyFrom(le_k_);
            if (rank_in_my_group_ < n_proc_for_fft_ && idx_to_my_group_ < n_group_trans_) {
                // mm_k_trans_
                mm_k_trans_.copySizeInfoFrom(mm_k_);  // overwrite
                mm_k_trans_.start_loc_[0] = 0;
                mm_k_trans_.start_loc_[1] = local_2_start_trans_[idx_to_my_group_];
                mm_k_trans_.start_loc_[2] = local_1_start_trans_[idx_to_my_group_];
                mm_k_trans_.start_loc_[3] = local_0_start_trans_[idx_to_my_group_];
                mm_k_trans_.end_loc_[0] = n_mm_compo - 1;  // all (l,m)
                mm_k_trans_.end_loc_[1] = local_2_end_trans_[idx_to_my_group_];
                mm_k_trans_.end_loc_[2] = local_1_end_trans_[idx_to_my_group_];
                mm_k_trans_.end_loc_[3] = local_0_end_trans_[idx_to_my_group_];
                mm_k_trans_.calcSizeLoc();
                mm_k_trans_.calcSizeTot();
                // le_k_trans_
                le_k_trans_.copySizeInfoFrom(mm_k_trans_);  // overwrite
            }

            // Determine whether to use OpenMP-parallerized FFT or not
            if (fft_size_crit_ < 0) {
                use_ompfft_ = determineWhetherToUseOMPFFT();
            } else {
                const S32 fft_size = n_cell.x * n_cell.y * n_cell.z;
                if (fft_size < fft_size_crit_)
                    use_ompfft_ = false;
                else
                    use_ompfft_ = true;
            }
        } else {
            // PS::DebugUtils::fout() << "M2L initialize skipped." << std::endl;
        }
        wtime_tmp = GetWtime();
        wtime_end_initialize__main_ = wtime_tmp;
        wtime_end_initialize_ = wtime_tmp;
        time_profile_.PMMM__M2L__initialize += (wtime_end_initialize_ - wtime_start_initialize_);
    }
#endif

    void setGreenFunction() {
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_SET_GREEN_FUNCTION
        comm_info_.barrier();
        if (comm_info_.getRank() == 0) std::cout << "ok1 @ setGreenFunction [M2L_Engine]" << std::endl;
#endif

        // Extract parameters
        const S32 p = param_.p;
        const S32 LEN2 = (2 * p + 1) * (2 * p + 1);

        // Set the size information of gf_k_trans_
        gf_k_trans_.copySizeInfoFrom(mm_k_trans_);
        gf_k_trans_.start_glb_[0] = 0;
        gf_k_trans_.end_glb_[0] = LEN2 - 1;
        gf_k_trans_.calcSizeGlb();
        if (rank_in_my_group_ < n_proc_for_fft_ && idx_to_my_group_ < n_group_trans_) {
            gf_k_trans_.start_loc_[0] = 0;
            gf_k_trans_.end_loc_[0] = LEN2 - 1;  // all (l,m)
            gf_k_trans_.calcSizeLoc();
        }
        gf_k_trans_.calcSizeTot();
        // gf_k_trans_.outputSizeInfo();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_SET_GREEN_FUNCTION
        comm_info_.barrier();
        if (comm_info_.getRank() == 0) std::cout << "ok2 @ setGreenFunction [M2L_Engine]" << std::endl;
#endif

        // Initialize green function calculator
        gf_calc_.initialize(param_, use_mpifft_, use_ompfft_, comm_info_, n_thread_);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_SET_GREEN_FUNCTION
        comm_info_.barrier();
        if (comm_info_.getRank() == 0) std::cout << "ok3 @ setGreenFunction [M2L_Engine]" << std::endl;
#endif

        // Calculate green function using the calculator
        gf_calc_.calcGreenFunction(gf_k_trans_, use_mpifft_);

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_SET_GREEN_FUNCTION
        comm_info_.barrier();
        if (comm_info_.getRank() == 0) std::cout << "ok4 @ setGreenFunction [M2L_Engine]" << std::endl;
#endif

        // Copy wall time information
        wtime_start_GFC_initialize_ = gf_calc_.wtime_start_initialize_;
        wtime_end_GFC_initialize_ = gf_calc_.wtime_end_initialize_;
        wtime_start_GFC_calc_gf_r_ = gf_calc_.wtime_start_calc_gf_r_;
        wtime_end_GFC_calc_gf_r_ = gf_calc_.wtime_end_calc_gf_r_;
        wtime_start_GFC_redist_gf_r_ = gf_calc_.wtime_start_redist_gf_r_;
        wtime_end_GFC_redist_gf_r_ = gf_calc_.wtime_end_redist_gf_r_;
        wtime_start_GFC_calc_gf_k_ = gf_calc_.wtime_start_calc_gf_k_;
        wtime_end_GFC_calc_gf_k_ = gf_calc_.wtime_end_calc_gf_k_;
        wtime_start_GFC_redist_gf_k_ = gf_calc_.wtime_start_redist_gf_k_;
        wtime_end_GFC_redist_gf_k_ = gf_calc_.wtime_end_redist_gf_k_;

        wtime_start_GFC_initialize__copy_ = gf_calc_.wtime_start_initialize__copy_;
        wtime_end_GFC_initialize__copy_ = gf_calc_.wtime_end_initialize__copy_;
        wtime_start_GFC_initialize__main_ = gf_calc_.wtime_start_initialize__main_;
        wtime_end_GFC_initialize__main_ = gf_calc_.wtime_end_initialize__main_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_SET_GREEN_FUNCTION
        comm_info_.barrier();
        if (comm_info_.getRank() == 0) std::cout << "ok_last @ setGreenFunction [M2L_Engine]" << std::endl;
#endif
    }

    template <class Cell_t>
    void redistMM(const int n_cell_loc, const std::vector<Cell_t> &cell_loc) {
        // Extract parameters
        // const S32 p = param_.p;
        const S32vec n_cell = param_.n_cell;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // Redistribute MM data
        if (mode_ == 1) {
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok1a @ redistMM [M2L_Engine]" << std::endl;
#endif
            using pair_t = std::pair<S32, real_t>;
            wtime_start_preproc_mm_r_comm_ = GetWtime();
            // Make a send buffer
            CommBuffer<pair_t> sendbuf;
            sendbuf.n_comm = n_proc_in_parent_group_;
            sendbuf.allocCommInfo();
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                const S32 rnk = (rank_in_parent_group_ + i) % n_proc_in_parent_group_;
                sendbuf.ranks[i] = rnk;
                sendbuf.adr_from_rank[rnk] = i;
            }
            sendbuf.clearCounts();
            sendbuf.count_tot = 0;
            for (S32 n = 0; n < n_cell_loc; n++) {
                if (!cell_loc[n].is_mm_defined) continue;
                const S32 idx = cell_loc[n].idx;
                const S32 iz = idx / (n_cell.x * n_cell.y);
                // S32 slab_id;
                S32 slab_id = 0;
                for (S32 k = 0; k < n_proc_for_fft_; k++) {
                    if ((local_0_start_[k] <= iz) && (iz <= local_0_end_[k])) {
                        slab_id = k;
                        break;
                    }
                }
                for (S32 k = 0; k < n_group_; k++) {
                    const S32 rnk = rank_start_[k] + slab_id;
                    const S32 adr = sendbuf.adr_from_rank[rnk];
                    sendbuf.counts[adr]++;
                    sendbuf.count_tot++;
                }
            }
            std::vector<S32> sendcounts;  // used to set recvcounts
            sendcounts.resize(n_proc_in_parent_group_);
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                const S32 rank = sendbuf.ranks[i];
                const S32 cnt = sendbuf.counts[i];
                sendcounts[rank] = cnt;
            }
            if (sendbuf.count_tot > 0) {
                sendbuf.allocBuffer();
                sendbuf.calcDispls();
                sendbuf.clearCounts();
                for (S32 n = 0; n < n_cell_loc; n++) {
                    if (!cell_loc[n].is_mm_defined) continue;
                    const S32 idx = cell_loc[n].idx;
                    const S32 iz = idx / (n_cell.x * n_cell.y);
                    // S32 slab_id;
                    S32 slab_id = 0;
                    for (S32 k = 0; k < n_proc_for_fft_; k++) {
                        if ((local_0_start_[k] <= iz) && (iz <= local_0_end_[k])) {
                            slab_id = k;
                            break;
                        }
                    }
                    for (S32 k = 0; k < n_group_; k++) {
                        const S32 rnk = rank_start_[k] + slab_id;
                        const S32 adr = sendbuf.adr_from_rank[rnk];
                        const S32 adr_buf = sendbuf.displs[adr] + sendbuf.counts[adr];
                        sendbuf.buf[adr_buf].first = idx;
                        sendbuf.buf[adr_buf].second = cell_loc[n].mm.buf[k];
                        sendbuf.counts[adr]++;
                        // NaN check [option]
#ifdef PARTICLE_SIMULATOR_NAN_CHECK_PMMM_M2L_REDIST_MM
                        if (std::isnan(cell_loc[n].mm.buf[k])) {
                            DebugUtils::fout() << "[NaN] n = " << n << ", idx = " << idx << ", iz = " << iz << ", slab_id = " << slab_id
                                               << ", k = " << k << ", rnk = " << rnk << std::endl;
                        }
#endif
                    }
                }
            }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok2a @ redistMM [M2L_Engine]" << std::endl;
            // Output information of sendbuf
            DebugUtils::fout() << "--- Information of sendbuf (start) @ redistMM [M2L_Engine] ---" << std::endl;
            for (S32 i = 0; i < sendbuf.n_comm; i++) {
                DebugUtils::fout() << "i = " << i << ", rnk = " << sendbuf.ranks[i] << ", displ = " << sendbuf.displs[i]
                                   << ", count = " << sendbuf.counts[i] << std::endl;
            }
            DebugUtils::fout() << "--- Information of sendbuf (end) @ redistMM [M2L_Engine] ---" << std::endl;
#endif
            // Make a receive buffer
            CommBuffer<pair_t> recvbuf;
            recvbuf.n_comm = n_proc_in_parent_group_;
            recvbuf.allocCommInfo();
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                recvbuf.ranks[i] = i;
                recvbuf.adr_from_rank[i] = i;
            }
            std::vector<S32> recvcounts;
            recvcounts.resize(n_proc_in_parent_group_);
            MPI_Alltoall(sendcounts.data(), 1, GetDataType<S32>(), recvcounts.data(), 1, GetDataType<S32>(), parent_comm_);
            recvbuf.clearCounts();
            recvbuf.count_tot = 0;
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                recvbuf.counts[i] = recvcounts[i];
                recvbuf.count_tot += recvcounts[i];
            }
            if (recvbuf.count_tot > 0) {
                recvbuf.allocBuffer();
                recvbuf.calcDispls();
            }
            wtime_end_preproc_mm_r_comm_ = GetWtime();
            time_profile_.PMMM__M2L__preproc_mm_r_comm += (wtime_end_preproc_mm_r_comm_ - wtime_start_preproc_mm_r_comm_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok3a @ redistMM [M2L_Engine]" << std::endl;
#endif
            // Perform MPI comm.
            wtime_start_redist_mm_r_ = GetWtime();
            performComm(sendbuf, recvbuf, parent_comm_);
            wtime_end_redist_mm_r_ = GetWtime();
            time_profile_.PMMM__M2L__redist_mm_r += (wtime_end_redist_mm_r_ - wtime_start_redist_mm_r_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok4a @ redistMM [M2L_Engine]" << std::endl;
#endif
            // Copy from recvbuf to mm_r
            wtime_start_postproc_mm_r_comm_ = GetWtime();
            mm_r_.resize();
            mm_r_.clear();
            for (S32 i = 0; i < recvbuf.count_tot; i++) {
                // (1) Convert recvbuf[i].idx to S32vec idx_3d.
                const S32 idx_1d = recvbuf.buf[i].first;  // idx
                S32vec idx_3d;
                idx_3d.z = idx_1d / (n_cell.x * n_cell.y);
                idx_3d.y = (idx_1d - (n_cell.x * n_cell.y) * idx_3d.z) / n_cell.x;
                idx_3d.x = idx_1d - (n_cell.x * n_cell.y) * idx_3d.z - n_cell.x * idx_3d.y;
                // (2) Convert idx_3d to the one-dimensional cell index
                //     based on the size information of mm_r
                const S32 i_loc = idx_3d.x - mm_r_.start_loc_[1];
                const S32 j_loc = idx_3d.y - mm_r_.start_loc_[2];
                const S32 k_loc = idx_3d.z - mm_r_.start_loc_[3];
                const S32 idx = i_loc + mm_r_.size_loc_[1] * (j_loc + mm_r_.size_loc_[2] * k_loc);
                // Note that recvbuf[i].idx is calculated based on
                // the size of PM mesh, n_cell. On the other hand,
                // the cell index of mm_r[] is defined based on
                // its size, which depends on the boundary condition.
                // This is the reason why we have performed the index
                // calculation above.
                assert(0 <= idx && idx < mm_r_.size_loc_tot_);
                mm_r_.buf_[idx] += recvbuf.buf[i].second;
                // Nan check [option]
#ifdef PARTICLE_SIMULATOR_NAN_CHECK_PMMM_M2L_REDIST_MM
                if (std::isnan(recvbuf.buf[i].second)) {
                    // Examine which rank data comes from
                    S32 rank_src{-1};
                    S32 displ{-1}, count{-1};
                    S32 offset{-1};
                    for (S32 k = 0; k < recvbuf.n_comm; k++) {
                        const S32 head = recvbuf.displs[k];
                        const S32 tail = head + recvbuf.counts[k];
                        if (head <= i && i < tail) {
                            rank_src = recvbuf.ranks[k];
                            displ = head;
                            count = recvbuf.counts[k];
                            offset = i - head;
                            break;
                        }
                    }
                    assert(rank_src != -1);
                    // Output information
                    DebugUtils::fout() << "[NaN] i = " << i << ", idx_1d = " << idx_1d << ", idx_3d = " << idx_3d << ", arr_idx = " << idx
                                       << ", rank_src = " << rank_src << ", displ = " << displ << ", count = " << count << ", offset = " << offset
                                       << std::endl;
                }
#endif
            }
            wtime_end_postproc_mm_r_comm_ = GetWtime();
            time_profile_.PMMM__M2L__postproc_mm_r_comm += (wtime_end_postproc_mm_r_comm_ - wtime_start_postproc_mm_r_comm_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok5a @ redistMM [M2L_Engine]" << std::endl;
#endif
        } else if (mode_ == 2) {
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok1b @ redistMM [M2L_Engine]" << std::endl;
#endif
            wtime_start_preproc_mm_r_comm_ = GetWtime();
            // Make a list of # of local PM cells
            S32 *n_cell_loc_tbl = new S32[n_proc_in_parent_group_];
            comm_info_.allGather(&n_cell_loc, 1, n_cell_loc_tbl);
            // Make a send buffer
            union buf_t {
                real_t r;  // to store mm
                S32 i;     // to store idx
            };
            CommBuffer<buf_t> sendbuf;
            sendbuf.n_comm = n_proc_in_parent_group_;
            sendbuf.allocCommInfo();
            sendbuf.count_tot = 0;
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                const S32 n_lm = lm_end_[i] - lm_start_[i] + 1;
                sendbuf.ranks[i] = i;
                sendbuf.counts[i] = n_cell_loc * (1 + n_lm);
                sendbuf.count_tot += n_cell_loc * (1 + n_lm);
            }
            if (sendbuf.count_tot > 0) {
                sendbuf.allocBuffer();
                sendbuf.calcDispls();
                sendbuf.clearCounts();
                for (S32 n = 0; n < n_cell_loc; n++) {
                    const S32 idx = cell_loc[n].idx;
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        // set idx
                        S32 adr_buf = sendbuf.displs[i] + sendbuf.counts[i];
                        sendbuf.buf[adr_buf++].i = idx;
                        sendbuf.counts[i]++;
                        // set mm
                        for (S32 lm = lm_start_[i]; lm <= lm_end_[i]; lm++) {
                            sendbuf.buf[adr_buf++].r = cell_loc[n].mm.buf[lm];
                            sendbuf.counts[i]++;
                        }
                    }
                }
            }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok2b @ redistMM [M2L_Engine]" << std::endl;
#endif
            // Make a receive buffer
            CommBuffer<buf_t> recvbuf;
            recvbuf.n_comm = n_proc_in_parent_group_;
            recvbuf.allocCommInfo();
            recvbuf.count_tot = 0;
            const S32 n_lm = lm_end_[rank_in_parent_group_] - lm_start_[rank_in_parent_group_] + 1;
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                recvbuf.ranks[i] = i;
                recvbuf.counts[i] = n_cell_loc_tbl[i] * (1 + n_lm);
                recvbuf.count_tot += n_cell_loc_tbl[i] * (1 + n_lm);
            }
            if (recvbuf.count_tot > 0) {
                recvbuf.allocBuffer();
                recvbuf.calcDispls();
            }
            wtime_end_preproc_mm_r_comm_ = GetWtime();
            time_profile_.PMMM__M2L__preproc_mm_r_comm += (wtime_end_preproc_mm_r_comm_ - wtime_start_preproc_mm_r_comm_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok3b @ redistMM [M2L_Engine]" << std::endl;
#endif
            // Perform MPI comm.
            wtime_start_redist_mm_r_ = GetWtime();
            performComm(sendbuf, recvbuf, parent_comm_);
            wtime_end_redist_mm_r_ = GetWtime();
            time_profile_.PMMM__M2L__redist_mm_r += (wtime_end_redist_mm_r_ - wtime_start_redist_mm_r_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok4b @ redistMM [M2L_Engine]" << std::endl;
#endif
            // Copy from recvbuf to mm_r
            wtime_start_postproc_mm_r_comm_ = GetWtime();
            mm_r_.resize();
            mm_r_.clear();
            for (S32 n = 0; n < recvbuf.n_comm; n++) {
                const S32 head = recvbuf.displs[n];
                const S32 tail = head + recvbuf.counts[n];
                S32 adr_src{head};
                while (adr_src < tail) {
                    // (1) Convert idx_recv.buf[] to S32vec idx_3d.
                    const S32 idx_1d = recvbuf.buf[adr_src++].i;
                    S32vec idx_3d;
                    idx_3d.z = idx_1d / (n_cell.x * n_cell.y);
                    idx_3d.y = (idx_1d - (n_cell.x * n_cell.y) * idx_3d.z) / n_cell.x;
                    idx_3d.x = idx_1d - (n_cell.x * n_cell.y) * idx_3d.z - n_cell.x * idx_3d.y;
                    // (2) Copy
                    const S32 beg = lm_start_[rank_in_parent_group_];
                    const S32 end = lm_end_[rank_in_parent_group_];
                    for (S32 lm = beg; lm <= end; lm++) {
                        // (2-1) convert idx_3d to the one-dimensional cell index
                        //       based on the size information of mm_r
                        const S32 lm_loc = lm - mm_r_.start_loc_[0];
                        const S32 i_loc = idx_3d.x - mm_r_.start_loc_[1];
                        const S32 j_loc = idx_3d.y - mm_r_.start_loc_[2];
                        const S32 k_loc = idx_3d.z - mm_r_.start_loc_[3];
                        const S32 adr_dest = lm_loc + mm_r_.size_loc_[0] * (i_loc + mm_r_.size_loc_[1] * (j_loc + mm_r_.size_loc_[2] * k_loc));
                        // Note that idx in recvbuf is calculated based on
                        // the size of PM mesh, n_cell. On the other hand,
                        // the cell index of mm_r[] is defined based on
                        // its size, which depends on the boundary condition.
                        // This is the reason why we have performed the index
                        // calculation above.
                        assert(0 <= adr_dest && adr_dest < mm_r_.size_loc_tot_);
                        // (2-2) Copy mm_recv.buf[] to mm_r.buf_[]
                        mm_r_.buf_[adr_dest] = recvbuf.buf[adr_src++].r;
                    }
                }
                assert(adr_src == tail);
            }
            wtime_end_postproc_mm_r_comm_ = GetWtime();
            time_profile_.PMMM__M2L__postproc_mm_r_comm += (wtime_end_postproc_mm_r_comm_ - wtime_start_postproc_mm_r_comm_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
            comm_info_.barrier();
            if (comm_info_.getRank() == 0) std::cout << "ok5b @ redistMM [M2L_Engine]" << std::endl;
#endif
        }
#else  // PARTICLE_SIMULATOR_MPI_PARALLEL
       // Copy from cell_loc to mm_r_
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
        std::cout << "ok1c @ redistMM [M2L_Engine]" << std::endl;
#endif
        const F64 wtime_start = GetWtime();
        wtime_start_preproc_mm_r_comm_ = wtime_start;
        wtime_end_preproc_mm_r_comm_ = wtime_start;
        wtime_start_redist_mm_r_ = wtime_start;
        mm_r_.resize();
        mm_r_.clear();
        for (S32 i = 0; i < n_cell_loc; i++) {
            if (!cell_loc[i].is_mm_defined) continue;
            const S32 idx_1d = cell_loc[i].idx;
            S32vec idx_3d;
            idx_3d.z = idx_1d / (n_cell.x * n_cell.y);
            idx_3d.y = (idx_1d - (n_cell.x * n_cell.y) * idx_3d.z) / n_cell.x;
            idx_3d.x = idx_1d - (n_cell.x * n_cell.y) * idx_3d.z - n_cell.x * idx_3d.y;
            const S32 i_loc = idx_3d.x - mm_r_.start_loc_[1];
            const S32 j_loc = idx_3d.y - mm_r_.start_loc_[2];
            const S32 k_loc = idx_3d.z - mm_r_.start_loc_[3];
            for (S32 lm = 0; lm < cell_loc[i].mm.size(); lm++) {
                const S32 lm_loc = lm - mm_r_.start_loc_[0];
                const S32 idx = lm_loc + mm_r_.size_loc_[0] * (i_loc + mm_r_.size_loc_[1] * (j_loc + mm_r_.size_loc_[2] * k_loc));
                assert(0 <= idx && idx < mm_r_.size_loc_tot_);
                mm_r_.buf_[idx] += cell_loc[i].mm.buf[lm];
            }
        }
        const F64 wtime_end = GetWtime();
        wtime_end_redist_mm_r_ = wtime_end;
        time_profile_.PMMM__M2L__redist_mm_r += (wtime_end_redist_mm_r_ - wtime_start_redist_mm_r_);
        wtime_start_postproc_mm_r_comm_ = wtime_end;
        wtime_end_postproc_mm_r_comm_ = wtime_end;
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
        std::cout << "ok2c @ redistMM [M2L_Engine]" << std::endl;
#endif
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_M2L_REDIST_MM
        comm_info_.barrier();
        if (comm_info_.getRank() == 0) std::cout << "ok_last @ redistMM [M2L_Engine]" << std::endl;
        // Check mm_r
#if 0
                static S32 n_called {0};
                n_called++; 
                if (n_called == 1) {
                    const std::string file_prefix = "mm_r_";
                    const S32 file_num = rank_in_parent_group_;
                    mm_r_.writeBufferToFile(file_prefix, file_num);
                    Finalize();
                    std::exit(0);
                }
#endif
#endif
    }

    void redistMMK() {
        // Extract parameters
        // const S32 p = param_.p;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // Make a send buffer
        CommBuffer<cplx_t> sendbuf;
        sendbuf.n_comm = n_group_;
        sendbuf.allocCommInfo();
        for (S32 i = 0; i < n_group_; i++) {
            const S32 rnk = (idx_to_my_group_ + i) % n_group_;  // rank in comm_int_
            sendbuf.ranks[i] = rnk;
            sendbuf.adr_from_rank[rnk] = i;
        }
        sendbuf.count_tot = 0;
        sendbuf.clearCounts();
        for (S32 i = 0; i < n_group_; i++) {
            if (i < n_group_trans_) {
                const S32 cnt = (lm_end_[idx_to_my_group_] - lm_start_[idx_to_my_group_] + 1) *
                                (local_0_end_trans_[i] - local_0_start_trans_[i] + 1) * (local_1_end_trans_[i] - local_1_start_trans_[i] + 1) *
                                (local_2_end_trans_[i] - local_2_start_trans_[i] + 1);
                // [Note]
                //    cnt = (# of (l,m) that THIS PROCESS HAS)
                //        * (# of cells that a remote process requires)
                const S32 rnk = i;  // rank in comm_int_
                const S32 adr = sendbuf.adr_from_rank[rnk];
                sendbuf.counts[adr] = cnt;
                sendbuf.count_tot += cnt;
            }
        }
        if (sendbuf.count_tot > 0) {
            sendbuf.allocBuffer();
            sendbuf.calcDispls();
            for (S32 n = 0; n < n_group_; n++) {
                if (n < n_group_trans_) {
                    const S32 rnk = n;  // rank in comm_int_
                    const S32 adr = sendbuf.adr_from_rank[rnk];
                    S32 adr_dest = sendbuf.displs[adr];
                    for (S32 k = local_0_start_trans_[n]; k <= local_0_end_trans_[n]; k++)
                        for (S32 j = local_1_start_trans_[n]; j <= local_1_end_trans_[n]; j++)
                            for (S32 i = local_2_start_trans_[n]; i <= local_2_end_trans_[n]; i++)
                                for (S32 lm = lm_start_[idx_to_my_group_]; lm <= lm_end_[idx_to_my_group_]; lm++) {  // for(k,j,i,lm)
                                    const S32 lm_src = lm - mm_k_.start_loc_[0];
                                    const S32 i_src = i - mm_k_.start_loc_[1];
                                    const S32 j_src = j - mm_k_.start_loc_[2];
                                    const S32 k_src = k - mm_k_.start_loc_[3];
                                    const S32 adr_src =
                                        lm_src + mm_k_.size_loc_[0] * (i_src + mm_k_.size_loc_[1] * (j_src + mm_k_.size_loc_[2] * k_src));
                                    sendbuf.buf[adr_dest++] = mm_k_.buf_[adr_src];
                                }  // for(k,j,i,lm)
                }
            }
        }
        // Make a receive buffer
        CommBuffer<cplx_t> recvbuf;
        recvbuf.n_comm = n_group_;
        recvbuf.allocCommInfo();
        for (S32 i = 0; i < n_group_; i++) {
            recvbuf.ranks[i] = i;  // rank in comm_int_
            recvbuf.adr_from_rank[i] = i;
        }
        recvbuf.count_tot = 0;
        recvbuf.clearCounts();
        if (idx_to_my_group_ < n_group_trans_) {
            // In this case, this process receive data.
            const S32 cnt_cell = (local_0_end_trans_[idx_to_my_group_] - local_0_start_trans_[idx_to_my_group_] + 1) *
                                 (local_1_end_trans_[idx_to_my_group_] - local_1_start_trans_[idx_to_my_group_] + 1) *
                                 (local_2_end_trans_[idx_to_my_group_] - local_2_start_trans_[idx_to_my_group_] + 1);
            for (S32 i = 0; i < n_group_; i++) {
                const S32 cnt = (lm_end_[i] - lm_start_[i] + 1) * cnt_cell;
                recvbuf.counts[i] = cnt;
                recvbuf.count_tot += cnt;
            }
        }
        if (recvbuf.count_tot > 0) {
            recvbuf.allocBuffer();
            recvbuf.calcDispls();
        }
        // Perfrom MPI comm.
        performComm(sendbuf, recvbuf, comm_int_);
        // Copy from recvbuf to mm_k_trans_
        mm_k_trans_.resize();
        const S32 k_start = local_0_start_trans_[idx_to_my_group_];
        const S32 k_end = local_0_end_trans_[idx_to_my_group_];
        const S32 j_start = local_1_start_trans_[idx_to_my_group_];
        const S32 j_end = local_1_end_trans_[idx_to_my_group_];
        const S32 i_start = local_2_start_trans_[idx_to_my_group_];
        const S32 i_end = local_2_end_trans_[idx_to_my_group_];
        for (S32 n = 0; n < n_group_; n++) {
            const S32 rnk = recvbuf.ranks[n];
            const S32 adr = recvbuf.adr_from_rank[rnk];
            S32 adr_src = recvbuf.displs[adr];
            const S32 lm_start = lm_start_[n];
            const S32 lm_end = lm_end_[n];
            for (S32 k = k_start; k <= k_end; k++)
                for (S32 j = j_start; j <= j_end; j++)
                    for (S32 i = i_start; i <= i_end; i++)
                        for (S32 lm = lm_start; lm <= lm_end; lm++) {  // for(k,j,i,lm)
                            const S32 lm_dest = lm - mm_k_trans_.start_loc_[0];
                            const S32 i_dest = i - mm_k_trans_.start_loc_[1];
                            const S32 j_dest = j - mm_k_trans_.start_loc_[2];
                            const S32 k_dest = k - mm_k_trans_.start_loc_[3];
                            const S32 adr_dest = lm_dest + mm_k_trans_.size_loc_[0] *
                                                               (i_dest + mm_k_trans_.size_loc_[1] * (j_dest + mm_k_trans_.size_loc_[2] * k_dest));
                            mm_k_trans_.buf_[adr_dest] = recvbuf.buf[adr_src++];
                        }  // for(k,j,i,lm)
        }
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
        // Copy from mm_k_ to mm_k_trans_
        assert(mm_k_.size_loc_tot_ == mm_k_trans_.size_loc_tot_);
        for (S32 i = 0; i < mm_k_trans_.size_loc_tot_; i++) {
            mm_k_trans_.buf_[i] = mm_k_.buf_[i];
        }
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
    }

    void redistLEK() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // Make a send buffer
        CommBuffer<cplx_t> sendbuf;
        sendbuf.n_comm = n_group_;
        sendbuf.allocCommInfo();
        for (S32 i = 0; i < n_group_; i++) {
            const S32 rnk = (idx_to_my_group_ + i) % n_group_;  // rank in comm_int_
            sendbuf.ranks[i] = rnk;
            sendbuf.adr_from_rank[rnk] = i;
        }
        sendbuf.count_tot = 0;
        sendbuf.clearCounts();
        if (idx_to_my_group_ < n_group_trans_) {
            // In this case, this process has le_k_trans_ to be sent.
            const S32 cnt_cell = (local_0_end_trans_[idx_to_my_group_] - local_0_start_trans_[idx_to_my_group_] + 1) *
                                 (local_1_end_trans_[idx_to_my_group_] - local_1_start_trans_[idx_to_my_group_] + 1) *
                                 (local_2_end_trans_[idx_to_my_group_] - local_2_start_trans_[idx_to_my_group_] + 1);
            for (S32 i = 0; i < n_group_; i++) {
                const S32 cnt = (lm_end_[i] - lm_start_[i] + 1) * cnt_cell;
                const S32 rnk = i;  // rank in comm_int_
                const S32 adr = sendbuf.adr_from_rank[rnk];
                sendbuf.counts[adr] = cnt;
                sendbuf.count_tot += cnt;
            }
        }
        if (sendbuf.count_tot > 0) {
            sendbuf.allocBuffer();
            sendbuf.calcDispls();
            const S32 k_start = local_0_start_trans_[idx_to_my_group_];
            const S32 k_end = local_0_end_trans_[idx_to_my_group_];
            const S32 j_start = local_1_start_trans_[idx_to_my_group_];
            const S32 j_end = local_1_end_trans_[idx_to_my_group_];
            const S32 i_start = local_2_start_trans_[idx_to_my_group_];
            const S32 i_end = local_2_end_trans_[idx_to_my_group_];
            for (S32 n = 0; n < n_group_; n++) {
                const S32 rnk = sendbuf.ranks[n];  // rank in comm_int_
                const S32 adr = sendbuf.adr_from_rank[rnk];
                S32 adr_dest = sendbuf.displs[adr];
                const S32 lm_start = lm_start_[rnk];
                const S32 lm_end = lm_end_[rnk];
                for (S32 k = k_start; k <= k_end; k++)
                    for (S32 j = j_start; j <= j_end; j++)
                        for (S32 i = i_start; i <= i_end; i++)
                            for (S32 lm = lm_start; lm <= lm_end; lm++) {  // for(k,j,i,lm)
                                const S32 lm_src = lm - le_k_trans_.start_loc_[0];
                                const S32 i_src = i - le_k_trans_.start_loc_[1];
                                const S32 j_src = j - le_k_trans_.start_loc_[2];
                                const S32 k_src = k - le_k_trans_.start_loc_[3];
                                const S32 adr_src = lm_src + le_k_trans_.size_loc_[0] *
                                                                 (i_src + le_k_trans_.size_loc_[1] * (j_src + le_k_trans_.size_loc_[2] * k_src));
                                sendbuf.buf[adr_dest++] = le_k_trans_.buf_[adr_src];
                            }  // for(k,j,i,lm)
            }
        }
        // Make a receive buffer
        CommBuffer<cplx_t> recvbuf;
        recvbuf.n_comm = n_group_;
        recvbuf.allocCommInfo();
        for (S32 i = 0; i < n_group_; i++) {
            recvbuf.ranks[i] = i;  // rank in comm_int_
            recvbuf.adr_from_rank[i] = i;
        }
        recvbuf.count_tot = 0;
        recvbuf.clearCounts();
        for (S32 i = 0; i < n_group_; i++) {
            if (i < n_group_trans_) {
                const S32 cnt = (lm_end_[idx_to_my_group_] - lm_start_[idx_to_my_group_] + 1) *
                                (local_0_end_trans_[i] - local_0_start_trans_[i] + 1) * (local_1_end_trans_[i] - local_1_start_trans_[i] + 1) *
                                (local_2_end_trans_[i] - local_2_start_trans_[i] + 1);
                // [Note]
                //    cnt = (# of (l,m) that this process is responsible for)
                //        * (# of cells that a remote process has)
                const S32 rnk = i;  // rank in comm_int_
                const S32 adr = recvbuf.adr_from_rank[rnk];
                recvbuf.counts[adr] = cnt;
                recvbuf.count_tot += cnt;
            }
        }
        if (recvbuf.count_tot > 0) {
            recvbuf.allocBuffer();
            recvbuf.calcDispls();
        }
        // Perform MPI comm.
        performComm(sendbuf, recvbuf, comm_int_);
        // Copy from recvbuf to le_k_
        const S32 lm_start = lm_start_[idx_to_my_group_];
        const S32 lm_end = lm_end_[idx_to_my_group_];
        for (S32 n = 0; n < n_group_; n++) {
            const S32 k_start = local_0_start_trans_[n];
            const S32 k_end = local_0_end_trans_[n];
            const S32 j_start = local_1_start_trans_[n];
            const S32 j_end = local_1_end_trans_[n];
            const S32 i_start = local_2_start_trans_[n];
            const S32 i_end = local_2_end_trans_[n];
            S32 adr_src = recvbuf.displs[n];
            for (S32 k = k_start; k <= k_end; k++)
                for (S32 j = j_start; j <= j_end; j++)
                    for (S32 i = i_start; i <= i_end; i++)
                        for (S32 lm = lm_start; lm <= lm_end; lm++) {  // for(k,j,i,lm)
                            const S32 lm_dest = lm - le_k_.start_loc_[0];
                            const S32 i_dest = i - le_k_.start_loc_[1];
                            const S32 j_dest = j - le_k_.start_loc_[2];
                            const S32 k_dest = k - le_k_.start_loc_[3];
                            const S32 adr_dest =
                                lm_dest + le_k_.size_loc_[0] * (i_dest + le_k_.size_loc_[1] * (j_dest + le_k_.size_loc_[2] * k_dest));
                            le_k_.buf_[adr_dest] = recvbuf.buf[adr_src++];
                        }  // for(k,j,i)
        }
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
        // Copy from le_k_trans_ to le_k_
        assert(le_k_.size_loc_tot_ == le_k_trans_.size_loc_tot_);
        for (S32 i = 0; i < le_k_.size_loc_tot_; i++) {
            le_k_.buf_[i] = le_k_trans_.buf_[i];
        }
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
    }

    void convolution() {
        // Extract parameters
        // const S32 p = param_.p;
        // const S32vec n_cell = param_.n_cell;
        // const S32 LEN = (p + 1) * (p + 1);
        // Save the previous FFT configuration
        fft_mode_prev_ = fft_mode_;
        for (S32 i = 0; i < 3; i++) size_fft_prev_[i] = size_fft_[i];
        // Resize buffers
        mm_k_.resize();
        mm_k_trans_.resize();
        le_k_trans_.resize();
        le_k_.resize();
        le_r_.resize();
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
        constexpr S32 rank_wld_debug_print = 0;
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok1 @ convolution() [PMMM]" << std::endl;
#endif
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        if (mode_ == 1 && use_mpifft_) {
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
            comm_info_.barrier();
            if (Comm::getRank() == rank_wld_debug_print) std::cout << "okA @ convolution() [PMMM]" << std::endl;
#endif
            if (rank_in_my_group_ < n_proc_for_fft_) {
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) {
                    std::cout << "okAA @ convolution() [PMMM]" << std::endl;
                }
                std::cout << "rank_in_my_group_= " << rank_in_my_group_ << std::endl;
                std::cout << "n_proc_for_fft_= " << n_proc_for_fft_ << std::endl;
#endif
                wtime_start_mm_r_to_mm_k_ = GetWtime();
                // Check if conditions that need to be fulfilled in this IF branch
                // are actually fulfilled.
                assert(mm_r_.size_loc_[0] == 1);
                assert(le_k_.size_loc_[0] == 1);
                assert(le_r_.size_loc_[0] == 1);
                assert(mm_r_.start_loc_[0] == mm_k_.start_loc_[0]);
                assert(mm_r_.end_loc_[0] == mm_k_.end_loc_[0]);

                // Set FFT configuration
                fft_mode_ = 1;
                size_fft_[0] = mm_r_.size_glb_[1];
                size_fft_[1] = mm_r_.size_glb_[2];
                size_fft_[2] = mm_r_.size_glb_[3];

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok2 @ convolution() [PMMM]" << std::endl;
#endif

                if ((fft_mode_ != fft_mode_prev_) || (size_fft_[0] != size_fft_prev_[0]) || (size_fft_[1] != size_fft_prev_[1]) ||
                    (size_fft_[2] != size_fft_prev_[2])) {
                    // Allocate buffer memory
                    ptrdiff_t alloc_local;
                    ptrdiff_t local_n0, local_0_start;
                    ptrdiff_t local_n1, local_1_start;
                    alloc_local = fftw_mpi_local_size_3d_transposed(mm_k_.size_glb_[3], mm_k_.size_glb_[2], mm_k_.size_glb_[1], comm_fft_, &local_n0,
                                                                    &local_0_start, &local_n1, &local_1_start);
                    // Note that an user must indicate the logical size of FFT by
                    // the sizes of a COMPLEX array when using fftw_mpi_local_size_3d_transposed.
                    const ptrdiff_t size_rbuf = 2 * alloc_local;
                    const ptrdiff_t size_kbuf = alloc_local;
                    if (rbuf_ != nullptr) fftw_free(rbuf_);
                    if (kbuf_ != nullptr) fftw_free(kbuf_);
                    rbuf_ = fftw_alloc_real(size_rbuf);
                    kbuf_ = fftw_alloc_complex(size_kbuf);

                    // Destroy plan if needed
                    for (size_t i = 0; i < plan_fwd_.size(); i++) fftw_destroy_plan(plan_fwd_[i]);
                    plan_fwd_.resize(0);
                    for (size_t i = 0; i < plan_bkw_.size(); i++) fftw_destroy_plan(plan_bkw_[i]);
                    plan_bkw_.resize(0);

                    // Create plan
                    fftw_plan tmp;
                    tmp = fftw_mpi_plan_dft_r2c_3d(mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], &rbuf_[0], &kbuf_[0], comm_fft_,
                                                   fftw_planning_rigor_flag_ | FFTW_DESTROY_INPUT | FFTW_MPI_TRANSPOSED_OUT);
                    plan_fwd_.push_back(tmp);
                    tmp = fftw_mpi_plan_dft_c2r_3d(mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], &kbuf_[0], &rbuf_[0], comm_fft_,
                                                   fftw_planning_rigor_flag_ | FFTW_DESTROY_INPUT | FFTW_MPI_TRANSPOSED_IN);
                    plan_bkw_.push_back(tmp);
                    // Note that fftw_mpi_plan_dft_r2c_3d requires an user to indicate
                    // the logical size of FFT by the size of a REAL array.
                }
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok3 @ convolution() [PMMM]" << std::endl;
#endif
                // Perform FFT (forward multipole moments)
                PS_OMP(omp parallel for num_threads(n_thread_))
                for (S32 k = mm_r_.start_loc_[3]; k <= mm_r_.end_loc_[3]; k++)
                    for (S32 j = mm_r_.start_loc_[2]; j <= mm_r_.end_loc_[2]; j++)
                        for (S32 i = mm_r_.start_loc_[1]; i <= mm_r_.end_loc_[1]; i++) {
                            const S32 i_loc = i - mm_r_.start_loc_[1];
                            const S32 j_loc = j - mm_r_.start_loc_[2];
                            const S32 k_loc = k - mm_r_.start_loc_[3];
                            const S32 adr_src = i_loc + mm_r_.size_loc_[1] * (j_loc + mm_r_.size_loc_[2] * k_loc);
                            const S32 i_buf = i - mm_r_.start_loc_[1];
                            const S32 j_buf = j - mm_r_.start_loc_[2];
                            const S32 k_buf = k - mm_r_.start_loc_[3];
                            const S32 adr_dest = i_buf + (2 * (mm_r_.size_glb_[1] / 2 + 1)) * (j_buf + mm_r_.size_glb_[2] * k_buf);
                            // Note that the size of 1st dimension of rbuf is not NX.
                            // (see Section 6.5 of the mamual of FFTW3)
                            rbuf_[adr_dest] = mm_r_.buf_[adr_src];
                        }

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok4 @ convolution() [PMMM]" << std::endl;
#endif

                // CALL FFTW
                fftw_execute(plan_fwd_[0]);

                PS_OMP(omp parallel for num_threads(n_thread_))
                for (S32 j = mm_k_.start_loc_[3]; j <= mm_k_.end_loc_[3]; j++)
                    for (S32 k = mm_k_.start_loc_[2]; k <= mm_k_.end_loc_[2]; k++)
                        for (S32 i = mm_k_.start_loc_[1]; i <= mm_k_.end_loc_[1]; i++) {
                            const S32 i_buf = i - mm_k_.start_loc_[1];
                            const S32 k_buf = k - mm_k_.start_loc_[2];
                            const S32 j_buf = j - mm_k_.start_loc_[3];
                            const S32 adr_src = i_buf + mm_k_.size_glb_[1] * (k_buf + mm_k_.size_glb_[2] * j_buf);
                            const S32 lm = mm_k_.start_loc_[0];
                            const S32 lm_loc = lm - mm_k_.start_loc_[0];
                            const S32 i_loc = i - mm_k_.start_loc_[1];
                            const S32 k_loc = k - mm_k_.start_loc_[2];
                            const S32 j_loc = j - mm_k_.start_loc_[3];
                            const S32 adr_dest = lm_loc + mm_k_.size_loc_[0] * (i_loc + mm_k_.size_loc_[1] * (k_loc + mm_k_.size_loc_[2] * j_loc));
                    // Note that j and k loops are exchanged because we specified
                    // FFTW_MPI_TRANSPOSED_OUT.

#if defined(PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER)
                            const cplx_t ctmp(kbuf_[adr_src][0], kbuf_[adr_src][1]);
                            mm_k_.buf_[adr_dest] = ctmp;
                    // [Notes (tag: #64d4cd48)]
                    //     Fujitsu C++ compiler installed in K computer
                    //     does not support std::complex<T>::real(T value)
                    //     and std::complex<T>::imag(T value).
#else
                            mm_k_.buf_[adr_dest].real(kbuf_[adr_src][0]);
                            mm_k_.buf_[adr_dest].imag(kbuf_[adr_src][1]);
#endif
                            // Note that fftw_complex = double[2]
                        }

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok5 @ convolution() [PMMM]" << std::endl;
#endif
                wtime_end_mm_r_to_mm_k_ = GetWtime();
                time_profile_.PMMM__M2L__mm_r_to_mm_k += (wtime_end_mm_r_to_mm_k_ - wtime_start_mm_r_to_mm_k_);

                // Gather mm_k_trans_
                wtime_start_gather_mm_k_trans_ = GetWtime();
                redistMMK();
                wtime_end_gather_mm_k_trans_ = GetWtime();
                time_profile_.PMMM__M2L__gather_mm_k_trans += (wtime_end_gather_mm_k_trans_ - wtime_start_gather_mm_k_trans_);

                // M2L transformation
                wtime_start_transform_ = GetWtime();
                transform();
                wtime_end_transform_ = GetWtime();
                time_profile_.PMMM__M2L__transform += (wtime_end_transform_ - wtime_start_transform_);

                // Scatter le_k_trans
                wtime_start_scatter_le_k_trans_ = GetWtime();
                redistLEK();
                wtime_end_scatter_le_k_trans_ = GetWtime();
                time_profile_.PMMM__M2L__scatter_le_k_trans += (wtime_end_scatter_le_k_trans_ - wtime_start_scatter_le_k_trans_);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok6 @ convolution() [PMMM]" << std::endl;
#endif

                // Peform FFT (backward local expansion)
                wtime_start_le_k_to_le_r_ = GetWtime();
                PS_OMP(omp parallel for num_threads(n_thread_))
                for (S32 j = le_k_.start_loc_[3]; j <= le_k_.end_loc_[3]; j++)
                    for (S32 k = le_k_.start_loc_[2]; k <= le_k_.end_loc_[2]; k++)
                        for (S32 i = le_k_.start_loc_[1]; i <= le_k_.end_loc_[1]; i++) {
                            const S32 i_loc = i - le_k_.start_loc_[1];
                            const S32 k_loc = k - le_k_.start_loc_[2];
                            const S32 j_loc = j - le_k_.start_loc_[3];
                            const S32 adr_src = i_loc + le_k_.size_loc_[1] * (k_loc + le_k_.size_loc_[2] * j_loc);
                            const S32 i_buf = i - le_k_.start_loc_[1];
                            const S32 k_buf = k - le_k_.start_loc_[2];
                            const S32 j_buf = j - le_k_.start_loc_[3];
                            const S32 adr_dest = i_buf + le_k_.size_glb_[1] * (k_buf + le_k_.size_glb_[2] * j_buf);
                            kbuf_[adr_dest][0] = le_k_.buf_[adr_src].real();
                            kbuf_[adr_dest][1] = le_k_.buf_[adr_src].imag();
                        }
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok7 @ convolution() [PMMM]" << std::endl;
#endif
                fftw_execute(plan_bkw_[0]);

                const F64 norm = 1.0 / (le_r_.size_glb_[1] * le_r_.size_glb_[2] * le_r_.size_glb_[3]);
                PS_OMP(omp parallel for num_threads(n_thread_))
                for (S32 k = le_r_.start_loc_[3]; k <= le_r_.end_loc_[3]; k++)
                    for (S32 j = le_r_.start_loc_[2]; j <= le_r_.end_loc_[2]; j++)
                        for (S32 i = le_r_.start_loc_[1]; i <= le_r_.end_loc_[1]; i++) {
                            const S32 i_buf = i - le_r_.start_loc_[1];
                            const S32 j_buf = j - le_r_.start_loc_[2];
                            const S32 k_buf = k - le_r_.start_loc_[3];
                            const S32 adr_src = i_buf + (2 * (le_r_.size_glb_[1] / 2 + 1)) * (j_buf + le_r_.size_glb_[2] * k_buf);
                            const S32 i_loc = i - le_r_.start_loc_[1];
                            const S32 j_loc = j - le_r_.start_loc_[2];
                            const S32 k_loc = k - le_r_.start_loc_[3];
                            const S32 adr_dest = i_loc + le_r_.size_loc_[1] * (j_loc + le_r_.size_loc_[2] * k_loc);
                            le_r_.buf_[adr_dest] = norm * rbuf_[adr_src];
                        }
                wtime_end_le_k_to_le_r_ = GetWtime();
                time_profile_.PMMM__M2L__le_k_to_le_r += (wtime_end_le_k_to_le_r_ - wtime_start_le_k_to_le_r_);

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok8 @ convolution() [PMMM]" << std::endl;
#endif
            } else {
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "okAB @ convolution() [PMMM]" << std::endl;
#endif

                const F64 now = GetWtime();
                wtime_start_mm_r_to_mm_k_ = now;
                wtime_end_mm_r_to_mm_k_ = now;
                wtime_start_gather_mm_k_trans_ = now;
                wtime_end_gather_mm_k_trans_ = now;
                wtime_start_transform_ = now;
                wtime_end_transform_ = now;
                wtime_start_scatter_le_k_trans_ = now;
                wtime_end_scatter_le_k_trans_ = now;
                wtime_start_le_k_to_le_r_ = now;
                wtime_end_le_k_to_le_r_ = now;
            }
        } else {
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
        // In this case, we can use only a single MPI process to calculate
        // mm_k for a particular set of (l,m).
        // We use OpenMP parallerization of FFT if available.
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
            comm_info_.barrier();
            if (Comm::getRank() == rank_wld_debug_print) std::cout << "okB @ convolution() [PMMM]" << std::endl;
#endif
            if (rank_in_my_group_ == 0) {
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBA @ convolution() [PMMM]" << std::endl;
#endif
                // Check if conditions that need to be fulfilled in this IF branch
                // are actually fulfilled.
                assert(mm_r_.start_loc_[0] == mm_k_.start_loc_[0]);
                assert(mm_r_.end_loc_[0] == mm_k_.end_loc_[0]);
                for (S32 i = 1; i <= 3; i++) {
                    assert(mm_r_.size_loc_[i] == mm_r_.size_glb_[i]);
                    assert(mm_k_.size_loc_[i] == mm_k_.size_glb_[i]);
                }

                // Calculate mm_k
                const S32 fft_size = mm_r_.size_glb_[3] * mm_r_.size_glb_[2] * mm_r_.size_glb_[1];
                if (!use_ompfft_) {
                    // In this case, the size of array is too small to speed up by multithreaded FFT.
                    // Hence, we assign a single FFT to each thread.
                    wtime_start_mm_r_to_mm_k_ = GetWtime();
                    // Set FFT configuration
                    fft_mode_ = 2;
                    size_fft_[0] = mm_r_.size_glb_[1];
                    size_fft_[1] = mm_r_.size_glb_[2];
                    size_fft_[2] = mm_r_.size_glb_[3];
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                    // comm_info_.barrier();
                    if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBA1 @ convolution() [PMMM]" << std::endl;
#endif
                    if ((fft_mode_ != fft_mode_prev_) || (size_fft_[0] != size_fft_prev_[0]) || (size_fft_[1] != size_fft_prev_[1]) ||
                        (size_fft_[2] != size_fft_prev_[2])) {
                        // Memory allocation
                        const S32 size_rbuf = n_thread_ * fft_size;
                        const S32 size_kbuf = n_thread_ * fft_size;
                        if (rbuf_ != nullptr) fftw_free(rbuf_);
                        if (kbuf_ != nullptr) fftw_free(kbuf_);
                        rbuf_ = fftw_alloc_real(size_rbuf);
                        kbuf_ = fftw_alloc_complex(size_kbuf);

                        // Destroy plan if needed
                        for (size_t i = 0; i < plan_fwd_.size(); i++) fftw_destroy_plan(plan_fwd_[i]);
                        plan_fwd_.resize(n_thread_);
                        for (size_t i = 0; i < plan_bkw_.size(); i++) fftw_destroy_plan(plan_bkw_[i]);
                        plan_bkw_.resize(n_thread_);

                        // Create plans of FFTW
                        for (S32 i = 0; i < n_thread_; i++) {
                            const S32 offset = fft_size * i;
                            plan_fwd_[i] = fftw_plan_dft_r2c_3d(mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], &rbuf_[offset],
                                                                &kbuf_[offset], fftw_planning_rigor_flag_ | FFTW_DESTROY_INPUT);
                            plan_bkw_[i] = fftw_plan_dft_c2r_3d(mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], &kbuf_[offset],
                                                                &rbuf_[offset], fftw_planning_rigor_flag_ | FFTW_DESTROY_INPUT);
                        }
                    }

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                    // comm_info_.barrier();
                    if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBA2 @ convolution() [PMMM]" << std::endl;
#endif

                    // Perform FFT (forward multipole)
                    PS_OMP(omp parallel for num_threads(n_thread_))
                    for (S32 lm = mm_r_.start_loc_[0]; lm <= mm_r_.end_loc_[0]; lm++) {
                        const S32 ith = comm_info_.getThreadNum();
                        const S32 offset = fft_size * ith;

                        S32 lm_loc = lm - mm_r_.start_loc_[0];
                        for (S32 k = mm_r_.start_loc_[3]; k <= mm_r_.end_loc_[3]; k++)
                            for (S32 j = mm_r_.start_loc_[2]; j <= mm_r_.end_loc_[2]; j++)
                                for (S32 i = mm_r_.start_loc_[1]; i <= mm_r_.end_loc_[1]; i++) {
                                    const S32 i_loc = i - mm_r_.start_loc_[1];
                                    const S32 j_loc = j - mm_r_.start_loc_[2];
                                    const S32 k_loc = k - mm_r_.start_loc_[3];
                                    const S32 adr_src =
                                        lm_loc + mm_r_.size_loc_[0] * (i_loc + mm_r_.size_loc_[1] * (j_loc + mm_r_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i + mm_r_.size_glb_[1] * (j + mm_r_.size_glb_[2] * k);
                                    rbuf_[adr_dest + offset] = mm_r_.buf_[adr_src];
                                }

                        fftw_execute(plan_fwd_[ith]);

                        lm_loc = lm - mm_k_.start_loc_[0];
                        for (S32 k = mm_k_.start_loc_[3]; k <= mm_k_.end_loc_[3]; k++)
                            for (S32 j = mm_k_.start_loc_[2]; j <= mm_k_.end_loc_[2]; j++)
                                for (S32 i = mm_k_.start_loc_[1]; i <= mm_k_.end_loc_[1]; i++) {
                                    const S32 i_loc = i - mm_k_.start_loc_[1];
                                    const S32 j_loc = j - mm_k_.start_loc_[2];
                                    const S32 k_loc = k - mm_k_.start_loc_[3];
                                    const S32 adr_src = i + mm_k_.size_glb_[1] * (j + mm_k_.size_glb_[2] * k);
                                    const S32 adr_dest =
                                        lm_loc + mm_k_.size_loc_[0] * (i_loc + mm_k_.size_loc_[1] * (j_loc + mm_k_.size_loc_[2] * k_loc));
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                    const cplx_t ctmp(kbuf_[adr_src + offset][0], kbuf_[adr_src + offset][1]);
                                    mm_k_.buf_[adr_dest] = ctmp;
#else
                                    mm_k_.buf_[adr_dest].real(kbuf_[adr_src + offset][0]);
                                    mm_k_.buf_[adr_dest].imag(kbuf_[adr_src + offset][1]);
#endif
                                }
                    }
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                    // comm_info_.barrier();
                    if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBA3 @ convolution() [PMMM]" << std::endl;
#endif
                    wtime_end_mm_r_to_mm_k_ = GetWtime();
                    time_profile_.PMMM__M2L__mm_r_to_mm_k += (wtime_end_mm_r_to_mm_k_ - wtime_start_mm_r_to_mm_k_);

#if 0
                            // Check mm_k
                            {
                                const std::string file_prefix = "mm_k_";
                                const S32 file_num = rank_in_parent_group_; 
                                mm_k_.writeBufferToFile(file_prefix, file_num, false, true);
                            }
                            Finalize();
                            std::exit(0);
#endif

                    // Gather mm_k_trans_
                    wtime_start_gather_mm_k_trans_ = GetWtime();
                    redistMMK();
                    wtime_end_gather_mm_k_trans_ = GetWtime();
                    time_profile_.PMMM__M2L__gather_mm_k_trans += (wtime_end_gather_mm_k_trans_ - wtime_start_gather_mm_k_trans_);
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                    // comm_info_.barrier();
                    if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBA4 @ convolution() [PMMM]" << std::endl;
#endif
                    // M2L transformation
                    wtime_start_transform_ = GetWtime();
                    transform();
                    wtime_end_transform_ = GetWtime();
                    time_profile_.PMMM__M2L__transform += (wtime_end_transform_ - wtime_start_transform_);
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                    // comm_info_.barrier();
                    if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBA5 @ convolution() [PMMM]" << std::endl;
#endif
#if 0
                            // Check le_k_trans_
                            {
                                const std::string file_prefix = "le_k_trans_";
                                const S32 file_num = rank_in_parent_group_; 
                                le_k_trans_.writeBufferToFile(file_prefix, file_num);
                            }
#endif

                    // Scatter le_k_trans_
                    wtime_start_scatter_le_k_trans_ = GetWtime();
                    redistLEK();
                    wtime_end_scatter_le_k_trans_ = GetWtime();
                    time_profile_.PMMM__M2L__scatter_le_k_trans += (wtime_end_scatter_le_k_trans_ - wtime_start_scatter_le_k_trans_);
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                    // comm_info_.barrier();
                    if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBA6 @ convolution() [PMMM]" << std::endl;
#endif
#if 0
                            // Check le_k_
                            {
                                const std::string file_prefix = "le_k_";
                                const S32 file_num = rank_in_parent_group_; 
                                le_k_.writeBufferToFile(file_prefix, file_num);
                            }
#endif

                    // Peform FFT (backward local expansion)
                    wtime_start_le_k_to_le_r_ = GetWtime();
                    PS_OMP(omp parallel for num_threads(n_thread_))
                    for (S32 lm = le_r_.start_loc_[0]; lm <= le_r_.end_loc_[0]; lm++) {
                        const S32 ith = comm_info_.getThreadNum();
                        const S32 offset = fft_size * ith;
                        const S32 lm_loc = lm - le_r_.start_loc_[0];
                        for (S32 k = le_k_.start_loc_[3]; k <= le_k_.end_loc_[3]; k++)
                            for (S32 j = le_k_.start_loc_[2]; j <= le_k_.end_loc_[2]; j++)
                                for (S32 i = le_k_.start_loc_[1]; i <= le_k_.end_loc_[1]; i++) {
                                    const S32 i_loc = i - le_k_.start_loc_[1];
                                    const S32 j_loc = j - le_k_.start_loc_[2];
                                    const S32 k_loc = k - le_k_.start_loc_[3];
                                    const S32 adr_src =
                                        lm_loc + le_k_.size_loc_[0] * (i_loc + le_k_.size_loc_[1] * (j_loc + le_k_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i + le_k_.size_glb_[1] * (j + le_k_.size_glb_[2] * k);
                                    kbuf_[adr_dest + offset][0] = le_k_.buf_[adr_src].real();
                                    kbuf_[adr_dest + offset][1] = le_k_.buf_[adr_src].imag();
                                }

                        fftw_execute(plan_bkw_[ith]);

                        const F64 norm = 1.0 / (le_r_.size_glb_[1] * le_r_.size_glb_[2] * le_r_.size_glb_[3]);
                        for (S32 k = le_r_.start_loc_[3]; k <= le_r_.end_loc_[3]; k++)
                            for (S32 j = le_r_.start_loc_[2]; j <= le_r_.end_loc_[2]; j++)
                                for (S32 i = le_r_.start_loc_[1]; i <= le_r_.end_loc_[1]; i++) {
                                    const S32 i_loc = i - le_r_.start_loc_[1];
                                    const S32 j_loc = j - le_r_.start_loc_[2];
                                    const S32 k_loc = k - le_r_.start_loc_[3];
                                    const S32 adr_src = i + le_r_.size_glb_[1] * (j + le_r_.size_glb_[2] * k);
                                    const S32 adr_dest =
                                        lm_loc + le_r_.size_loc_[0] * (i_loc + le_r_.size_loc_[1] * (j_loc + le_r_.size_loc_[2] * k_loc));
                                    le_r_.buf_[adr_dest] = norm * rbuf_[adr_src + offset];
                                }
                    }
                    wtime_end_le_k_to_le_r_ = GetWtime();
                    time_profile_.PMMM__M2L__le_k_to_le_r += (wtime_end_le_k_to_le_r_ - wtime_start_le_k_to_le_r_);
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                    // comm_info_.barrier();
                    if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBA7 @ convolution() [PMMM]" << std::endl;
#endif
                } else {
                    // In this case, we can expect that the FFT calculation become fast by using
                    // a multithreaded FFT.
                    wtime_start_mm_r_to_mm_k_ = GetWtime();
                    // Set FFT configuration
                    fft_mode_ = 3;
                    size_fft_[0] = mm_r_.size_glb_[1];
                    size_fft_[1] = mm_r_.size_glb_[2];
                    size_fft_[2] = mm_r_.size_glb_[3];

                    if ((fft_mode_ != fft_mode_prev_) || (size_fft_[0] != size_fft_prev_[0]) || (size_fft_[1] != size_fft_prev_[1]) ||
                        (size_fft_[2] != size_fft_prev_[2])) {
                        // Memory allocation
                        const S32 size_rbuf = mm_r_.size_glb_[3] * mm_r_.size_glb_[2] * mm_r_.size_glb_[1];
                        const S32 size_kbuf = mm_k_.size_glb_[3] * mm_k_.size_glb_[2] * mm_k_.size_glb_[1];
                        if (rbuf_ != nullptr) fftw_free(rbuf_);
                        if (kbuf_ != nullptr) fftw_free(kbuf_);
                        rbuf_ = fftw_alloc_real(size_rbuf);
                        kbuf_ = fftw_alloc_complex(size_kbuf);

                        // Destroy plan if needed
                        for (size_t i = 0; i < plan_fwd_.size(); i++) fftw_destroy_plan(plan_fwd_[i]);
                        plan_fwd_.resize(0);
                        for (size_t i = 0; i < plan_bkw_.size(); i++) fftw_destroy_plan(plan_bkw_[i]);
                        plan_bkw_.resize(0);

                        // Create plan
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                        fftw_plan_with_nthreads(n_thread_);
#endif
                        fftw_plan tmp;
                        tmp = fftw_plan_dft_r2c_3d(mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], &rbuf_[0], &kbuf_[0],
                                                   fftw_planning_rigor_flag_ | FFTW_DESTROY_INPUT);
                        plan_fwd_.push_back(tmp);
                        tmp = fftw_plan_dft_c2r_3d(mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], &kbuf_[0], &rbuf_[0],
                                                   fftw_planning_rigor_flag_ | FFTW_DESTROY_INPUT);
                        plan_bkw_.push_back(tmp);
                    }

                    // Clear
                    // for (S32 i = 0; i < size_rbuf; i++) rbuf[i] = 0;

                    // Perform FFT (forward multipole moments)
                    for (S32 lm = mm_r_.start_loc_[0]; lm <= mm_r_.end_loc_[0]; lm++) {
                        S32 lm_loc = lm - mm_r_.start_loc_[0];
                        PS_OMP(omp parallel for num_threads(n_thread_))
                        for (S32 k = mm_r_.start_loc_[3]; k <= mm_r_.end_loc_[3]; k++)
                            for (S32 j = mm_r_.start_loc_[2]; j <= mm_r_.end_loc_[2]; j++)
                                for (S32 i = mm_r_.start_loc_[1]; i <= mm_r_.end_loc_[1]; i++) {
                                    const S32 i_loc = i - mm_r_.start_loc_[1];
                                    const S32 j_loc = j - mm_r_.start_loc_[2];
                                    const S32 k_loc = k - mm_r_.start_loc_[3];
                                    const S32 adr_src =
                                        lm_loc + mm_r_.size_loc_[0] * (i_loc + mm_r_.size_loc_[1] * (j_loc + mm_r_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i + mm_r_.size_glb_[1] * (j + mm_r_.size_glb_[2] * k);
                                    rbuf_[adr_dest] = mm_r_.buf_[adr_src];
                                }

                        // CALL FFTW
                        fftw_execute(plan_fwd_[0]);

                        lm_loc = lm - mm_k_.start_loc_[0];
                        PS_OMP(omp parallel for num_threads(n_thread_))
                        for (S32 k = mm_k_.start_loc_[3]; k <= mm_k_.end_loc_[3]; k++)
                            for (S32 j = mm_k_.start_loc_[2]; j <= mm_k_.end_loc_[2]; j++)
                                for (S32 i = mm_k_.start_loc_[1]; i <= mm_k_.end_loc_[1]; i++) {
                                    const S32 i_loc = i - mm_k_.start_loc_[1];
                                    const S32 j_loc = j - mm_k_.start_loc_[2];
                                    const S32 k_loc = k - mm_k_.start_loc_[3];
                                    const S32 adr_dest =
                                        lm_loc + mm_k_.size_loc_[0] * (i_loc + mm_k_.size_loc_[1] * (j_loc + mm_k_.size_loc_[2] * k_loc));
                                    const S32 adr_src = i + mm_k_.size_glb_[1] * (j + mm_k_.size_glb_[2] * k);
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                    const cplx_t ctmp(kbuf_[adr_src][0], kbuf_[adr_src][1]);
                                    mm_k_.buf_[adr_dest] = ctmp;
#else
                                    mm_k_.buf_[adr_dest].real(kbuf_[adr_src][0]);
                                    mm_k_.buf_[adr_dest].imag(kbuf_[adr_src][1]);
#endif
                                    // Note that fftw_complex = double[2]
                                }
                    }
                    wtime_end_mm_r_to_mm_k_ = GetWtime();
                    time_profile_.PMMM__M2L__mm_r_to_mm_k += (wtime_end_mm_r_to_mm_k_ - wtime_start_mm_r_to_mm_k_);

#if 0
                            // Check mm_k
                            {
                                const std::string file_prefix = "mm_k_";
                                const S32 file_num = rank_in_parent_group_; 
                                mm_k_.writeBufferToFile(file_prefix, file_num);
                            }
                            Finalize();
                            std::exit(0);
#endif

                    // Gather mm_k_trans_
                    wtime_start_gather_mm_k_trans_ = GetWtime();
                    redistMMK();
                    wtime_end_gather_mm_k_trans_ = GetWtime();
                    time_profile_.PMMM__M2L__gather_mm_k_trans += (wtime_end_gather_mm_k_trans_ - wtime_start_gather_mm_k_trans_);

#if 0
                            // Check mm_k_trans_
                            {
                                const std::string file_prefix = "mm_k_trans_";
                                const S32 file_num = rank_in_parent_group_; 
                                mm_k_trans_.writeBufferToFile(file_prefix, file_num);
                            }
#endif

                    // M2L transformation
                    wtime_start_transform_ = GetWtime();
                    transform();
                    wtime_end_transform_ = GetWtime();
                    time_profile_.PMMM__M2L__transform += (wtime_end_transform_ - wtime_start_transform_);

#if 0
                            // Check le_k_trans_
                            {
                                const std::string file_prefix = "le_k_trans_";
                                const S32 file_num = rank_in_parent_group_; 
                                le_k_trans_.writeBufferToFile(file_prefix, file_num);
                            }
#endif

                    // Scatter le_k_trans
                    wtime_start_scatter_le_k_trans_ = GetWtime();
                    redistLEK();
                    wtime_end_scatter_le_k_trans_ = GetWtime();
                    time_profile_.PMMM__M2L__scatter_le_k_trans += (wtime_end_scatter_le_k_trans_ - wtime_start_scatter_le_k_trans_);

#if 0
                            // Check le_k_
                            {
                                const std::string file_prefix = "le_k_";
                                const S32 file_num = rank_in_parent_group_; 
                                le_k_.writeBufferToFile(file_prefix, file_num);
                            }
#endif

                    // Peform FFT (backward local expansion)
                    wtime_start_le_k_to_le_r_ = GetWtime();
                    for (S32 lm = le_r_.start_loc_[0]; lm <= le_r_.end_loc_[0]; lm++) {
                        const S32 lm_loc = lm - le_r_.start_loc_[0];
                        PS_OMP(omp parallel for)
                        for (S32 k = le_k_.start_loc_[3]; k <= le_k_.end_loc_[3]; k++)
                            for (S32 j = le_k_.start_loc_[2]; j <= le_k_.end_loc_[2]; j++)
                                for (S32 i = le_k_.start_loc_[1]; i <= le_k_.end_loc_[1]; i++) {
                                    const S32 i_loc = i - le_k_.start_loc_[1];
                                    const S32 j_loc = j - le_k_.start_loc_[2];
                                    const S32 k_loc = k - le_k_.start_loc_[3];
                                    const S32 adr_src =
                                        lm_loc + le_k_.size_loc_[0] * (i_loc + le_k_.size_loc_[1] * (j_loc + le_k_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i + le_k_.size_glb_[1] * (j + le_k_.size_glb_[2] * k);
                                    kbuf_[adr_dest][0] = le_k_.buf_[adr_src].real();
                                    kbuf_[adr_dest][1] = le_k_.buf_[adr_src].imag();
                                }

                        fftw_execute(plan_bkw_[0]);

                        const F64 norm = 1.0 / (le_r_.size_glb_[1] * le_r_.size_glb_[2] * le_r_.size_glb_[3]);
                        PS_OMP(omp parallel for)
                        for (S32 k = le_r_.start_loc_[3]; k <= le_r_.end_loc_[3]; k++)
                            for (S32 j = le_r_.start_loc_[2]; j <= le_r_.end_loc_[2]; j++)
                                for (S32 i = le_r_.start_loc_[1]; i <= le_r_.end_loc_[1]; i++) {
                                    const S32 i_loc = i - le_r_.start_loc_[1];
                                    const S32 j_loc = j - le_r_.start_loc_[2];
                                    const S32 k_loc = k - le_r_.start_loc_[3];
                                    const S32 adr_src = i + le_r_.size_glb_[1] * (j + le_r_.size_glb_[2] * k);
                                    const S32 adr_dest =
                                        lm_loc + le_r_.size_loc_[0] * (i_loc + le_r_.size_loc_[1] * (j_loc + le_r_.size_loc_[2] * k_loc));
                                    le_r_.buf_[adr_dest] = norm * rbuf_[adr_src];
                                }
                    }
                    wtime_end_le_k_to_le_r_ = GetWtime();
                    time_profile_.PMMM__M2L__le_k_to_le_r += (wtime_end_le_k_to_le_r_ - wtime_start_le_k_to_le_r_);
                }
            } else {
#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
                // comm_info_.barrier();
                if (Comm::getRank() == rank_wld_debug_print) std::cout << "okBB @ convolution() [PMMM]" << std::endl;
#endif
                // do nothing
                const F64 now = GetWtime();
                wtime_start_mm_r_to_mm_k_ = now;
                wtime_end_mm_r_to_mm_k_ = now;
                wtime_start_gather_mm_k_trans_ = now;
                wtime_end_gather_mm_k_trans_ = now;
                wtime_start_transform_ = now;
                wtime_end_transform_ = now;
                wtime_start_scatter_le_k_trans_ = now;
                wtime_end_scatter_le_k_trans_ = now;
                wtime_start_le_k_to_le_r_ = now;
                wtime_end_le_k_to_le_r_ = now;
            }  // END of if (rank_in_my_group_ == 0)
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        }
#endif

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CONVOLUTION)
        comm_info_.barrier();
        if (Comm::getRank() == rank_wld_debug_print) std::cout << "ok LAST @ convolution() [PMMM]" << std::endl;
#endif

#if 0
                // Check le_r
                {
                    const std::string file_prefix = "le_r_";
                    const S32 file_num = rank_in_parent_group_;
                    le_r_.writeBufferToFile(file_prefix, file_num);
                }
#endif
#if 0
                Finalize();
                std::exit(0);
#endif
    }

    template <class Cell_t>
    void redistLE(const std::vector<S32ort> &pos_pm_domain, const S32 n_cell_loc, std::unordered_map<S32, S32> &adr_cell_loc,
                  std::vector<Cell_t> &cell_loc) {
        // Extract basic parameters
        const S32 n_mm_compo = (param_.p + 1) * (param_.p + 1);
        const S32vec n_cell = param_.n_cell;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        if (mode_ == 1) {
            using pair_t = std::pair<S32, real_t>;
            wtime_start_preproc_le_r_comm_ = GetWtime();
            // Make a send buffer
            CommBuffer<pair_t> sendbuf;
            sendbuf.n_comm = n_proc_in_parent_group_;
            sendbuf.allocCommInfo();
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                const S32 rnk = (rank_in_parent_group_ + i) % n_proc_in_parent_group_;
                sendbuf.ranks[i] = rnk;
                sendbuf.adr_from_rank[rnk] = i;
            }
            sendbuf.clearCounts();
            sendbuf.count_tot = 0;
            if (rank_in_my_group_ < n_proc_for_fft_) {  // this process has LE data
                const S32 ix_start = 0;
                const S32 ix_end = n_cell.x - 1;
                const S32 iy_start = 0;
                const S32 iy_end = n_cell.y - 1;
                const S32 iz_start = local_0_start_[rank_in_my_group_];
                const S32 iz_end = std::min(local_0_end_[rank_in_my_group_], n_cell.z - 1);
                const S32ort pos_self = S32ort(S32vec(ix_start, iy_start, iz_start), S32vec(ix_end, iy_end, iz_end));
                if (pos_self.isValid()) {
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        const S32ort pos_trgt = pos_pm_domain[i];
                        if (pos_self.overlapped(pos_trgt)) {
                            const S32ort pos_intxn = pos_self.getIntersectionNoCheck(pos_trgt);
                            const S32 cnt = pos_intxn.getVolume();
                            const S32 adr = sendbuf.adr_from_rank[i];
                            sendbuf.counts[adr] += cnt;
                            sendbuf.count_tot += cnt;
                        }
                    }
                }
            }
            if (sendbuf.count_tot > 0) {
                sendbuf.allocBuffer();
                sendbuf.calcDispls();
                sendbuf.clearCounts();
                const S32 ix_start = 0;
                const S32 ix_end = n_cell.x - 1;
                const S32 iy_start = 0;
                const S32 iy_end = n_cell.y - 1;
                const S32 iz_start = local_0_start_[rank_in_my_group_];
                const S32 iz_end = std::min(local_0_end_[rank_in_my_group_], n_cell.z - 1);
                const S32ort pos_self = S32ort(S32vec(ix_start, iy_start, iz_start), S32vec(ix_end, iy_end, iz_end));
                if (pos_self.isValid()) {
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        const S32ort pos_trgt = pos_pm_domain[i];
                        if (pos_self.overlapped(pos_trgt)) {
                            const S32ort pos_intxn = pos_self.getIntersectionNoCheck(pos_trgt);
                            const S32 adr = sendbuf.adr_from_rank[i];
                            for (S32 iz = pos_intxn.low_.z; iz <= pos_intxn.high_.z; iz++) {
                                for (S32 iy = pos_intxn.low_.y; iy <= pos_intxn.high_.y; iy++) {
                                    for (S32 ix = pos_intxn.low_.x; ix <= pos_intxn.high_.x; ix++) {
                                        const S32 idx = ix + n_cell.x * (iy + n_cell.y * iz);
                                        const S32 i_loc = ix - le_r_.start_loc_[1];
                                        const S32 j_loc = iy - le_r_.start_loc_[2];
                                        const S32 k_loc = iz - le_r_.start_loc_[3];
                                        const S32 adr_src = i_loc + le_r_.size_loc_[1] * (j_loc + le_r_.size_loc_[2] * k_loc);
                                        const S32 adr_dest = sendbuf.displs[adr] + sendbuf.counts[adr];
                                        sendbuf.buf[adr_dest].first = idx;
                                        sendbuf.buf[adr_dest].second = le_r_.buf_[adr_src];
                                        sendbuf.counts[adr]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Make a receive buffer
            CommBuffer<pair_t> recvbuf;
            recvbuf.n_comm = n_proc_in_parent_group_;
            recvbuf.allocCommInfo();
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                recvbuf.ranks[i] = i;
                recvbuf.adr_from_rank[i] = i;
            }
            recvbuf.clearCounts();
            recvbuf.count_tot = 0;
            for (S32 n = 0; n < n_cell_loc; n++) {
                const S32 idx = cell_loc[n].idx;
                const S32 iz = idx / (n_cell.x * n_cell.y);
                // int slab_id;
                int slab_id = 0;
                for (S32 k = 0; k < n_proc_for_fft_; k++) {
                    if ((local_0_start_[k] <= iz) && (iz <= local_0_end_[k])) {
                        slab_id = k;
                        break;
                    }
                }
                for (S32 k = 0; k < n_group_; k++) {
                    const S32 rnk = rank_start_[k] + slab_id;
                    recvbuf.counts[rnk]++;
                    recvbuf.count_tot++;
                }
            }
            if (recvbuf.count_tot > 0) {
                recvbuf.allocBuffer();
                recvbuf.calcDispls();
            }
            wtime_end_preproc_le_r_comm_ = GetWtime();
            time_profile_.PMMM__M2L__preproc_le_r_comm += (wtime_end_preproc_le_r_comm_ - wtime_start_preproc_le_r_comm_);
            // Perform MPI comm.
            wtime_start_redist_le_r_ = GetWtime();
            performComm(sendbuf, recvbuf, parent_comm_);
            wtime_end_redist_le_r_ = GetWtime();
            time_profile_.PMMM__M2L__redist_le_r += (wtime_end_redist_le_r_ - wtime_start_redist_le_r_);
            // Copy from recvbuf to cell_loc
            wtime_start_postproc_le_r_comm_ = GetWtime();
            S32 *cnts = nullptr;
            if (n_cell_loc > 0) {
                cnts = new S32[n_cell_loc];
                for (S32 i = 0; i < n_cell_loc; i++) cnts[i] = 0;
            }
            for (S32 i = 0; i < recvbuf.count_tot; i++) {
                const S32 idx = recvbuf.buf[i].first;  // idx
                const S32 adr = adr_cell_loc[idx];
                // [Note (tag: #cee527db)]
                //     In C++11, we cannot use operator [] to access
                //     const std::unordered_map<,> and have to use member
                //     function `at` as follows:
                //
                //     adr = adr_cell_loc.at(idx);
                //
                //     However, Fujitsu C++ compiler does not support `at`.
                //     Hence, we had to remove a const qualifier from
                //     adr_cell_loc and used operator[].
                const real_t val = recvbuf.buf[i].second;  // le
                cell_loc[adr].le.buf[cnts[adr]++] = val;
            }
#ifdef DEBUG_M2L_ENGINE_REDISTLE
            for (S32 n = 0; n < n_cell_loc; n++) {
                const S32 size = cnts[n];
                if (size > 0 && size != n_mm_compo) {
                    std::cout << "cell_loc[].le is something wrong:" << std::endl;
                    std::cout << " my_rank = " << rank_in_parent_group_ << " adr = " << n << " idx = " << cell_loc[n].idx << std::endl;
                    assert(false);
                }
            }
#endif
            if (cnts != nullptr) delete[] cnts;
            wtime_end_postproc_le_r_comm_ = GetWtime();
            time_profile_.PMMM__M2L__postproc_le_r_comm += (wtime_end_postproc_le_r_comm_ - wtime_start_postproc_le_r_comm_);
        } else if (mode_ == 2) {
            wtime_start_preproc_le_r_comm_ = GetWtime();
            // Make a send buffer
            union buf_t {
                real_t r;  // to store le
                S32 i;     // to store idx
            };
            CommBuffer<buf_t> sendbuf;
            sendbuf.n_comm = n_proc_in_parent_group_;
            sendbuf.allocCommInfo();
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                sendbuf.ranks[i] = i;
                sendbuf.adr_from_rank[i] = i;
            }
            sendbuf.clearCounts();
            sendbuf.count_tot = 0;
            const S32 n_lm = le_r_.size_loc_[0];
            const S32 ix_start = le_r_.start_loc_[1];
            const S32 ix_end = std::min(le_r_.end_loc_[1], n_cell.x - 1);
            const S32 iy_start = le_r_.start_loc_[2];
            const S32 iy_end = std::min(le_r_.end_loc_[2], n_cell.y - 1);
            const S32 iz_start = le_r_.start_loc_[3];
            const S32 iz_end = std::min(le_r_.end_loc_[3], n_cell.z - 1);
            const S32ort pos_self = S32ort(S32vec(ix_start, iy_start, iz_start), S32vec(ix_end, iy_end, iz_end));
            if (pos_self.isValid()) {
                for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                    const S32ort pos_trgt = pos_pm_domain[i];
                    if (pos_self.overlapped(pos_trgt)) {
                        const S32ort pos_intxn = pos_self.getIntersectionNoCheck(pos_trgt);
                        const S32 cnt = pos_intxn.getVolume();
                        const S32 adr = sendbuf.adr_from_rank[i];
                        sendbuf.counts[adr] += cnt * (1 + n_lm);
                        sendbuf.count_tot += cnt * (1 + n_lm);
                    }
                }
            }
            if (sendbuf.count_tot > 0) {
                sendbuf.allocBuffer();
                sendbuf.calcDispls();
                sendbuf.clearCounts();
                for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                    const S32ort pos_trgt = pos_pm_domain[i];
                    if (pos_self.overlapped(pos_trgt)) {
                        const S32ort pos_intxn = pos_self.getIntersectionNoCheck(pos_trgt);
                        const S32 adr = sendbuf.adr_from_rank[i];
                        for (S32 iz = pos_intxn.low_.z; iz <= pos_intxn.high_.z; iz++) {
                            for (S32 iy = pos_intxn.low_.y; iy <= pos_intxn.high_.y; iy++) {
                                for (S32 ix = pos_intxn.low_.x; ix <= pos_intxn.high_.x; ix++) {
                                    const S32 idx = ix + n_cell.x * (iy + n_cell.y * iz);
                                    // set idx
                                    S32 adr_dest = sendbuf.displs[adr] + sendbuf.counts[adr];
                                    sendbuf.buf[adr_dest++].i = idx;
                                    sendbuf.counts[adr]++;
                                    // set le
                                    for (S32 lm = le_r_.start_loc_[0]; lm <= le_r_.end_loc_[0]; lm++) {
                                        const S32 lm_loc = lm - le_r_.start_loc_[0];
                                        const S32 i_loc = ix - le_r_.start_loc_[1];
                                        const S32 j_loc = iy - le_r_.start_loc_[2];
                                        const S32 k_loc = iz - le_r_.start_loc_[3];
                                        const S32 adr_src =
                                            lm_loc + le_r_.size_loc_[0] * (i_loc + le_r_.size_loc_[1] * (j_loc + le_r_.size_loc_[2] * k_loc));
                                        sendbuf.buf[adr_dest++].r = le_r_.buf_[adr_src];
                                        sendbuf.counts[adr]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Make a receive buffer
            CommBuffer<buf_t> recvbuf;
            recvbuf.n_comm = n_proc_in_parent_group_;
            recvbuf.allocCommInfo();
            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                recvbuf.ranks[i] = i;
                recvbuf.adr_from_rank[i] = i;
            }
            recvbuf.clearCounts();
            recvbuf.count_tot = 0;
            for (S32 n = 0; n < n_cell_loc; n++) {
                // const S32 idx = cell_loc[n].idx;
                for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                    const S32 n_lm = lm_end_[i] - lm_start_[i] + 1;
                    recvbuf.counts[i] += (1 + n_lm);
                    recvbuf.count_tot += (1 + n_lm);
                }
            }
            if (recvbuf.count_tot > 0) {
                recvbuf.allocBuffer();
                recvbuf.calcDispls();
            }
            wtime_end_preproc_le_r_comm_ = GetWtime();
            time_profile_.PMMM__M2L__preproc_le_r_comm += (wtime_end_preproc_le_r_comm_ - wtime_start_preproc_le_r_comm_);
            // Perform MPI comm.
            wtime_start_redist_le_r_ = GetWtime();
            performComm(sendbuf, recvbuf, parent_comm_);
            wtime_end_redist_le_r_ = GetWtime();
            time_profile_.PMMM__M2L__redist_le_r += (wtime_end_redist_le_r_ - wtime_start_redist_le_r_);
            // Copy receive buffers to cell_loc[]
            wtime_start_postproc_le_r_comm_ = GetWtime();
            for (S32 n = 0; n < recvbuf.n_comm; n++) {
                const S32 rank = recvbuf.ranks[n];
                const S32 head = recvbuf.displs[n];
                const S32 tail = head + recvbuf.counts[n];
                S32 adr_src{head};
                while (adr_src < tail) {
                    const S32 idx = recvbuf.buf[adr_src++].i;
                    const S32 adr = adr_cell_loc[idx];  // see Note #cee527db.
                    const S32 beg = lm_start_[rank];
                    const S32 end = lm_end_[rank];
                    for (S32 lm = beg; lm <= end; lm++) {
                        cell_loc[adr].le.buf[lm] = recvbuf.buf[adr_src++].r;
                    }
                }
                assert(adr_src == tail);
            }
            wtime_end_postproc_le_r_comm_ = GetWtime();
            time_profile_.PMMM__M2L__postproc_le_r_comm += (wtime_end_postproc_le_r_comm_ - wtime_start_postproc_le_r_comm_);
        }
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
        // Copy from le_r_[] to cell_loc[]
        const F64 wtime_start = GetWtime();
        wtime_start_preproc_le_r_comm_ = wtime_start;
        wtime_end_preproc_le_r_comm_ = wtime_start;
        wtime_start_redist_le_r_ = wtime_start;
        for (S32 k = le_r_.start_loc_[3]; k <= le_r_.end_loc_[3]; k++) {
            if (k < n_cell.z) {
                for (S32 j = le_r_.start_loc_[2]; j <= le_r_.end_loc_[2]; j++) {
                    if (j < n_cell.y) {
                        for (S32 i = le_r_.start_loc_[1]; i <= le_r_.end_loc_[1]; i++) {
                            if (i < n_cell.x) {
                                for (S32 lm = le_r_.start_loc_[0]; lm <= le_r_.end_loc_[0]; lm++) {
                                    const S32 lm_loc = lm - le_r_.start_loc_[0];
                                    const S32 i_loc = i - le_r_.start_loc_[1];
                                    const S32 j_loc = j - le_r_.start_loc_[2];
                                    const S32 k_loc = k - le_r_.start_loc_[3];
                                    const S32 adr_src =
                                        lm_loc + le_r_.size_loc_[0] * (i_loc + le_r_.size_loc_[1] * (j_loc + le_r_.size_loc_[2] * k_loc));
                                    const S32 idx = i + n_cell.x * (j + n_cell.y * k);
                                    const S32 adr_dest = adr_cell_loc[idx];  // see Note #cee527db.
                                    cell_loc[adr_dest].le.buf[lm] = le_r_.buf_[adr_src];
                                }
                            }
                        }
                    }
                }
            }
        }
        const F64 wtime_end = GetWtime();
        wtime_end_redist_le_r_ = wtime_end;
        time_profile_.PMMM__M2L__redist_le_r += (wtime_end_redist_le_r_ - wtime_start_redist_le_r_);
        wtime_start_postproc_le_r_comm_ = wtime_end;
        wtime_end_postproc_le_r_comm_ = wtime_end;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL

#if 0
                // Check le_r
                {
                    std::stringstream ss;
                    ss << "le_r_" << std::setfill('0') << std::setw(5)
                       << rank_in_parent_group_ << ".txt";
                    const std::string filename = ss.str();
                    std::ofstream output_file;
                    output_file.open(filename.c_str(), std::ios::trunc);
                    for (S32 i = 0; i < n_cell_loc; i++) {
                        const S32 idx_1d = cell_loc[i].idx;
                        S32vec idx_3d;
                        idx_3d.z = idx_1d / (n_cell.x * n_cell.y);
                        idx_3d.y = (idx_1d - (n_cell.x * n_cell.y) * idx_3d.z) / n_cell.x;
                        idx_3d.x = idx_1d - (n_cell.x * n_cell.y) * idx_3d.z - n_cell.x * idx_3d.y; 
                        for (S32 lm = 0; lm < cell_loc[i].le.size(); lm++) {
                            const S32 idx = lm 
                                          + n_mm_compo * (idx_3d.x
                                          + n_cell.x * (idx_3d.y
                                          + n_cell.y * idx_3d.z));
                            const real_t val = cell_loc[i].le.buf[lm];
                            output_file << idx << "    " << val << std::endl;
                        }
                    }
                    output_file.close();
                }
                Finalize();
                std::exit(0);
#endif
    }

    void clearTimeProfile() {
        time_profile_.clear();
        gf_calc_.clearTimeProfile();
    }

    void clearWallTime() {
        gf_calc_.clearWallTime();

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

        wtime_start_initialize_ = 0;
        wtime_end_initialize_ = 0;
        wtime_start_preproc_mm_r_comm_ = 0;
        wtime_end_preproc_mm_r_comm_ = 0;
        wtime_start_redist_mm_r_ = 0;
        wtime_end_redist_mm_r_ = 0;
        wtime_start_postproc_mm_r_comm_ = 0;
        wtime_end_postproc_mm_r_comm_ = 0;
        wtime_start_mm_r_to_mm_k_ = 0;
        wtime_end_mm_r_to_mm_k_ = 0;
        wtime_start_gather_mm_k_trans_ = 0;
        wtime_end_gather_mm_k_trans_ = 0;
        wtime_start_transform_ = 0;
        wtime_end_transform_ = 0;
        wtime_start_scatter_le_k_trans_ = 0;
        wtime_end_scatter_le_k_trans_ = 0;
        wtime_start_le_k_to_le_r_ = 0;
        wtime_end_le_k_to_le_r_ = 0;
        wtime_start_preproc_le_r_comm_ = 0;
        wtime_end_preproc_le_r_comm_ = 0;
        wtime_start_redist_le_r_ = 0;
        wtime_end_redist_le_r_ = 0;
        wtime_start_postproc_le_r_comm_ = 0;
        wtime_end_postproc_le_r_comm_ = 0;

        wtime_start_initialize__copy_ = 0;
        wtime_end_initialize__copy_ = 0;
        wtime_start_initialize__main_ = 0;
        wtime_end_initialize__main_ = 0;
    }

    void clearCounterAll() {
        clearWallTime();
        clearTimeProfile();
    }

    TimeProfile getTimeProfile() const { return time_profile_ + gf_calc_.getTimeProfile(); }

   private:
    bool determineWhetherToUseOMPFFT() {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        if (n_thread_ == 1) return false;

        bool ret;
        if (rank_in_parent_group_ == 0) {
            // Extract parameters
            const S32 p = param_.p;
            const S32 NZ = param_.n_cell.z;
            const S32 NY = param_.n_cell.y;
            const S32 NX = param_.n_cell.x;

            // First, meausre time to perform OpenMP-parallerized FFTs
            F64 etime_ompfft;
            {
                // Allocate
                fftw_real_t *rbuf = fftw_alloc_real(NZ * NY * NX);
                fftw_cplx_t *kbuf = fftw_alloc_complex(NZ * NY * (1 + NX / 2));
                // Create plans
                fftw_plan_with_nthreads(n_thread_);
                fftw_plan plan_fwd = fftw_plan_dft_r2c_3d(NZ, NY, NX, rbuf, kbuf, FFTW_ESTIMATE);
                fftw_plan plan_bkw = fftw_plan_dft_c2r_3d(NZ, NY, NX, kbuf, rbuf, FFTW_ESTIMATE);
                // Make input data
                const long seed = 19810614;
                srand48(seed);
                for (S32 k = 0; k < NZ; k++)
                    for (S32 j = 0; j < NY; j++)
                        for (S32 i = 0; i < NX; i++) {
                            const S32 idx = i + NX * (j + NY * k);
                            rbuf[idx] = 2.0 * drand48() - 1.0;
                        }
                // Measure the performance
                const S32 n_FFTs = (p + 1) * (p + 1);
                const F64 time_offset = GetWtime();
                for (S32 n = 0; n < n_FFTs; n++) {
                    fftw_execute(plan_fwd);
                    fftw_execute(plan_bkw);
                }
                etime_ompfft = GetWtime() - time_offset;
                // Destroy plans
                fftw_destroy_plan(plan_fwd);
                fftw_destroy_plan(plan_bkw);
                // Free
                fftw_free(rbuf);
                fftw_free(kbuf);
            }

            // Next, measure time to perform multiple FFTs simultaneously
            F64 etime_concurrent;
            {
                // Allocate
                fftw_real_t **rbuf;
                fftw_cplx_t **kbuf;
                rbuf = (fftw_real_t **)malloc((size_t)sizeof(fftw_real_t *) * n_thread_);
                kbuf = (fftw_cplx_t **)malloc((size_t)sizeof(fftw_cplx_t *) * n_thread_);
                for (S32 i = 0; i < n_thread_; i++) {
                    rbuf[i] = fftw_alloc_real(NZ * NY * NX);
                    kbuf[i] = fftw_alloc_complex(NZ * NY * (1 + NX / 2));
                }
                // Create plans
                fftw_plan *plan_fwd, *plan_bkw;
                plan_fwd = (fftw_plan *)malloc((size_t)sizeof(fftw_plan) * n_thread_);
                plan_bkw = (fftw_plan *)malloc((size_t)sizeof(fftw_plan) * n_thread_);
                fftw_plan_with_nthreads(1);
                for (S32 i = 0; i < n_thread_; i++) {
                    plan_fwd[i] = fftw_plan_dft_r2c_3d(NZ, NY, NX, rbuf[i], kbuf[i], FFTW_ESTIMATE);
                    plan_bkw[i] = fftw_plan_dft_c2r_3d(NZ, NY, NX, kbuf[i], rbuf[i], FFTW_ESTIMATE);
                }
                // Make input data
                const long seed = 19810614;
                srand48(seed);
                for (S32 n = 0; n < n_thread_; n++) {
                    for (S32 k = 0; k < NZ; k++)
                        for (S32 j = 0; j < NY; j++)
                            for (S32 i = 0; i < NX; i++) {
                                const S32 idx = i + NX * (j + NY * k);
                                rbuf[n][idx] = 2.0 * drand48() - 1.0;
                            }
                }

                // Measure the performance
                const S32 n_FFTs = (p + 1) * (p + 1);
                const F64 time_offset = GetWtime();
                PS_OMP(omp parallel for num_threads(n_thread_))
                for (S32 n = 0; n < n_FFTs; n++) {
                    const S32 ith = comm_info_.getThreadNum();
                    fftw_execute(plan_fwd[ith]);
                    fftw_execute(plan_bkw[ith]);
                }
                etime_concurrent = GetWtime() - time_offset;
                // Destroy plans
                for (S32 i = 0; i < n_thread_; i++) {
                    fftw_destroy_plan(plan_fwd[i]);
                    fftw_destroy_plan(plan_bkw[i]);
                }
                free(plan_fwd);
                free(plan_bkw);
                // Free
                for (S32 i = 0; i < n_thread_; i++) {
                    fftw_free(rbuf[i]);
                    fftw_free(kbuf[i]);
                }
                free(rbuf);
                free(kbuf);
            }

            // Choose the fastest way
            if (etime_ompfft < etime_concurrent)
                ret = true;
            else
                ret = false;
            // Check
            std::cout << "etime_ompfft     = " << etime_ompfft << std::endl;
            std::cout << "etime_concurrnet = " << etime_concurrent << std::endl;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Bcast(&ret, 1, MPI_C_BOOL, 0, parent_comm_);
#endif
        return ret;
#else   // PARTICLE_SIMULATOR_THREAD_PARALLEL
        return false;
#endif  // PARTICLE_SIMULATOR_THREAD_PARALLEL
    }

    void transform() {
        if (idx_to_my_group_ < n_group_trans_) {
            // Extract parameters
            const S32 p = param_.p;

            // Local buffers
            Slm<cplx_t> *slm;
            MultipoleMoment<cplx_t> *mm;
            LocalExpansion<cplx_t> *le;
            slm = new Slm<cplx_t>[n_thread_];
            mm = new MultipoleMoment<cplx_t>[n_thread_];
            le = new LocalExpansion<cplx_t>[n_thread_];
            {
                Slm<cplx_t> tmp;
                tmp.make_table(2 * p);
            }

            const S32 n_cell_loc = le_k_trans_.size_loc_[1] * le_k_trans_.size_loc_[2] * le_k_trans_.size_loc_[3];
            PS_OMP(omp parallel for num_threads(n_thread_))
            for (S32 n = 0; n < n_cell_loc; n++) {
                const S32 ith = comm_info_.getThreadNum();
                S32 tmp = n;  // 1D cell id
                const S32 k_loc = tmp / (le_k_trans_.size_loc_[1] * le_k_trans_.size_loc_[2]);
                tmp -= k_loc * (le_k_trans_.size_loc_[1] * le_k_trans_.size_loc_[2]);
                const S32 j_loc = tmp / le_k_trans_.size_loc_[1];
                tmp -= j_loc * le_k_trans_.size_loc_[1];
                const S32 i_loc = tmp;
                S32 adr;
                // Set mm
                adr = mm_k_trans_.start_glb_[0] +
                      mm_k_trans_.size_loc_[0] * (i_loc + mm_k_trans_.size_loc_[1] * (j_loc + mm_k_trans_.size_loc_[2] * k_loc));
                mm[ith].reinterpret(p, &mm_k_trans_.buf_[adr]);
                // Set slm
                adr = gf_k_trans_.start_glb_[0] +
                      gf_k_trans_.size_loc_[0] * (i_loc + gf_k_trans_.size_loc_[1] * (j_loc + gf_k_trans_.size_loc_[2] * k_loc));
                slm[ith].reinterpret(2 * p, &gf_k_trans_.buf_[adr]);
                // Set le
                adr = le_k_trans_.start_loc_[0] +
                      le_k_trans_.size_loc_[0] * (i_loc + le_k_trans_.size_loc_[1] * (j_loc + le_k_trans_.size_loc_[2] * k_loc));
                le[ith].reinterpret(p, &le_k_trans_.buf_[adr]);
                // do M2L
                slm[ith].transform_M2L(mm[ith], le[ith], false);
            }
            // Note that correspondence relation between indices (k,j) and
            // cell indices (iz,iy) depends on the value of use_mpifft_.
            // If use_mpifft_ = true, k and array index of `3` represent
            // y direction, while j and array index of `2` represent z
            // direction.

            // Free memory
            delete[] slm;
            delete[] mm;
            delete[] le;
        }
    }
};
}  // namespace ParticleMeshMultipole
}  // namespace ParticleSimulator
