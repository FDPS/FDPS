#pragma once

#ifdef PARTICLE_SIMULATOR_USE_SAMPLE_SORT
#include "samplesortlib.hpp"
#endif

#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif

constexpr bool DEBUG_FLAG_DD = false;

namespace ParticleSimulator {
template <typename T>
struct IsParticleSystem;

inline void CommunicatorSplit(S32 np[], CommInfo &comm_info_src, CommInfo comm_info_1d[], CommInfo comm_info_sub[]) {
    int rank_1d[DIMENSION_LIMIT];
    int rank_sub[DIMENSION_LIMIT];
    int rank_glb = comm_info_src.getRank();
    for (S32 d = DIMENSION_LIMIT - 1; d >= 0; d--) {
        rank_1d[d] = rank_glb % np[d];
        rank_glb /= np[d];
        comm_info_sub[d] = comm_info_src.split(rank_1d[d], rank_glb);
        rank_sub[d] = comm_info_sub[d].getRank();
        comm_info_1d[d] = comm_info_src.split(rank_sub[d], rank_glb);
    }
}

class DomainInfo {
   private:
    // bool flag_set_pos_root_domain_ = false;
    bool flag_set_pos_root_domain_[3];
    CommInfo comm_info_;
    TimeProfile time_profile_;
    ReallocatableArray<F64vec> pos_sample_tot_;
    ReallocatableArray<F64vec> pos_sample_loc_;

    ReallocatableArray<F64ort> pos_domain_;
    ReallocatableArray<F64ort> pos_domain_temp_;

    ReallocatableArray<S32> n_smp_array_;
    ReallocatableArray<S32> n_smp_disp_array_;

    F32 coef_ema_;
    S32 target_number_of_sample_particle_;
    S32 number_of_sample_particle_tot_;
    S32 number_of_sample_particle_loc_;
    S32 n_domain_[DIMENSION_LIMIT];  // in 2-dim, n_domain_[2] is always 1.

    F64ort pos_root_domain_;

    bool first_call_by_initialize;
    bool first_call_by_decomposeDomain;

    S32 boundary_condition_;
    bool periodic_axis_[DIMENSION_LIMIT];  // in 2-dim, periodic_axis_[2] is always false.

    CommInfo comm_info_1d_[DIMENSION_LIMIT];
    CommInfo comm_info_sub_[DIMENSION_LIMIT];

#if 1
    template <int cid>  // cid is the dimension index (0: x, 1: y, 2: z)
    void calculateBoundaryOfDomain1D(const S32 np, const F64vec pos_sample[], const S32 istart, const S32 iend, F64 &xlow, F64 &xhigh) {
        if (istart == 0)
            xlow = pos_root_domain_.low_[cid];
        else
            xlow = 0.5 * (pos_sample[istart - 1][cid] + pos_sample[istart][cid]);
        if (iend == np - 1)
            xhigh = pos_root_domain_.high_[cid];
        else
            xhigh = 0.5 * (pos_sample[iend][cid] + pos_sample[iend + 1][cid]);
    }
#else
    void calculateBoundaryOfDomainX(const S32 np, const F64vec pos_sample[], const S32 istart, const S32 iend, F64 &xlow, F64 &xhigh) {
        if (istart == 0)
            xlow = pos_root_domain_.low_.x;
        else
            xlow = 0.5 * (pos_sample[istart - 1].x + pos_sample[istart].x);
        if (iend == np - 1)
            xhigh = pos_root_domain_.high_.x;
        else
            xhigh = 0.5 * (pos_sample[iend].x + pos_sample[iend + 1].x);
    }
    void calculateBoundaryOfDomainY(const S32 np, const F64vec pos_sample[], const S32 istart, const S32 iend, F64 &xlow, F64 &xhigh) {
        if (istart == 0)
            xlow = pos_root_domain_.low_.y;
        else
            xlow = 0.5 * (pos_sample[istart - 1].y + pos_sample[istart].y);
        if (iend == np - 1)
            xhigh = pos_root_domain_.high_.y;
        else
            xhigh = 0.5 * (pos_sample[iend].y + pos_sample[iend + 1].y);
    }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    void calculateBoundaryOfDomainZ(const S32 np, const F64vec pos_sample[], const S32 istart, const S32 iend, F64 &xlow, F64 &xhigh) {
        if (istart == 0)
            xlow = pos_root_domain_.low_.z;
        else
            xlow = 0.5 * (pos_sample[istart - 1].z + pos_sample[istart].z);
        if (iend == np - 1)
            xhigh = pos_root_domain_.high_.z;
        else
            xhigh = 0.5 * (pos_sample[iend].z + pos_sample[iend + 1].z);
    }
#endif
#endif

   public:
    void initialize_impl(const F64 coef_ema) {
        const S32 n_proc = comm_info_.getNumberOfProc();
        const S32 my_rank = comm_info_.getRank();
        n_smp_array_.resizeNoInitialize(n_proc);
        n_smp_disp_array_.resizeNoInitialize(n_proc + 1);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        pos_domain_.resizeNoInitialize(n_proc);
        pos_domain_temp_.resizeNoInitialize(n_proc);
#else
        pos_domain_.resizeNoInitialize(1);
        pos_domain_temp_.resizeNoInitialize(1);
#endif
        coef_ema_ = coef_ema;
        target_number_of_sample_particle_ = 0;
        number_of_sample_particle_tot_ = 0;
        number_of_sample_particle_loc_ = 0;
        S32 rank_tmp[DIMENSION_LIMIT];
        SetNumberOfDomainMultiDimension<DIMENSION>(n_proc, my_rank, n_domain_, rank_tmp);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        CommunicatorSplit(n_domain_, comm_info_, comm_info_1d_, comm_info_sub_);
#endif
    }

    void setCommInfo(const CommInfo &c, const F64 coef_ema = FDPS_DFLT_VAL_COEF_EMA) {
        comm_info_ = c;
        initialize_impl(coef_ema);
    }
    CommInfo getCommInfo() const { return comm_info_; }
#ifdef TEST_VARIADIC_TEMPLATE
    void DomainInfoDummyFunc() {}
#endif

    TimeProfile getTimeProfile() const { return time_profile_; }
    void clearTimeProfile() { time_profile_.clear(); }
    void dumpTimeProfile(std::ostream &fout = std::cout, const S32 level = 0) const { time_profile_.dump_domain_info(fout, level); }

    DomainInfo() {
        first_call_by_initialize = true;
        first_call_by_decomposeDomain = true;
        pos_sample_tot_.setAllocMode(MemoryAllocMode::Pool);
        pos_sample_loc_.setAllocMode(MemoryAllocMode::Pool);
        periodic_axis_[0] = periodic_axis_[1] = false;
        flag_set_pos_root_domain_[0] = flag_set_pos_root_domain_[1] = false;
        pos_root_domain_.low_.x = -LARGE_FLOAT;
        pos_root_domain_.high_.x = LARGE_FLOAT;
        pos_root_domain_.low_.y = -LARGE_FLOAT;
        pos_root_domain_.high_.y = LARGE_FLOAT;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        periodic_axis_[2] = false;
        flag_set_pos_root_domain_[2] = false;
        pos_root_domain_.low_.z = -LARGE_FLOAT;
        pos_root_domain_.high_.z = LARGE_FLOAT;
#endif
        boundary_condition_ = BOUNDARY_CONDITION_OPEN;
        comm_info_.setCommunicator();
    }
    ~DomainInfo() {
        n_smp_array_.freeMem();
        n_smp_disp_array_.freeMem();
        pos_domain_.freeMem();
        pos_domain_temp_.freeMem();
        pos_sample_tot_.freeMem();
        pos_sample_loc_.freeMem();
    }

    void initialize(const F32 coef_ema = FDPS_DFLT_VAL_COEF_EMA, MemoryAllocMode alloc_mode = MemoryAllocMode::Default) {
        if (coef_ema < 0.0 || coef_ema > 1.0) {
            PARTICLE_SIMULATOR_PRINT_ERROR("The smoothing factor of an exponential moving average is must between 0 and 1.");
            std::cerr << "The smoothing factor of an exponential moving average is must between 0 and 1." << std::endl;
            Abort(-1);
        }
        assert(first_call_by_initialize);
        first_call_by_initialize = false;
        initialize_impl(coef_ema);
    }

    void setNumberOfDomainMultiDimension(const S32 nx, const S32 ny, const S32 nz = 1) {
        auto n_proc = comm_info_.getNumberOfProc();
        if (n_proc != nx * ny * nz) {
            PARTICLE_SIMULATOR_PRINT_ERROR("devided number of domains is not consistent with total processe number");
            Abort(-1);
        }
        n_domain_[0] = nx;
        n_domain_[1] = ny;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        n_domain_[2] = 1;
#else
        n_domain_[2] = nz;
#endif

        Comm::setNumberOfProcMultiDim(0, nx);
        Comm::setNumberOfProcMultiDim(1, ny);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        CommunicatorSplit(n_domain_, comm_info_, comm_info_1d_, comm_info_sub_);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        Comm::setNumberOfProcMultiDim(2, nz);
#endif  // PARTICLE_SIMULATOR_TWO_DIMENSION
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
    }

    void setDomain(const S32 nx, const S32 ny, const S32 nz = 1) { setNumberOfDomainMultiDimension(nx, ny, nz); }

    //////////////////////////
    // collectSampleParticle
    template <class Tpsys>
    void collectSampleParticle(Tpsys &psys, const bool clear, const F64 weight) {
        F64 time_offset = GetWtime();
        if (clear) {
            number_of_sample_particle_loc_ = 0;
        }

        const auto number_of_sample_particle = psys.getNumberOfSampleParticleLocal(weight);
        pos_sample_loc_.increaseSize(number_of_sample_particle);
        // psys.getSampleParticle2(number_of_sample_particle, &pos_sample_loc_[number_of_sample_particle_loc_]);
        psys.getSampleParticle2(number_of_sample_particle, pos_sample_loc_, number_of_sample_particle_loc_);

        number_of_sample_particle_loc_ += number_of_sample_particle;
        pos_sample_loc_.resizeNoInitialize(number_of_sample_particle_loc_);

        //time_profile_.decompose_domain_all__collect_sample_particle += GetWtime() - time_offset;
        time_profile_.select_sample_particle += GetWtime() - time_offset;
        time_profile_.decompose_domain_all += GetWtime() - time_offset;
        return;
    }

    template <class Tpsys>
    void collectSampleParticle(Tpsys &psys, const bool clear) {
        const F64 wgh = psys.getNumberOfParticleLocal();
        collectSampleParticle(psys, clear, wgh);
    }

    template <class Tpsys>
    void collectSampleParticle(Tpsys &psys) {
        const F64 wgh = psys.getNumberOfParticleLocal();
        const bool clear = true;
        collectSampleParticle(psys, clear, wgh);
    }
    // collectSampleParticle
    //////////////////////////

    //////////////////
    // decomposeDomain
    void decomposeDomain() {
        F64 time_offset = GetWtime();
        // ****** collect sample particles to process 0. ******
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        auto nproc = comm_info_.getNumberOfProc();
        auto myrank = comm_info_.getRank();
        auto t0 = GetWtime();
        comm_info_.allGather(&number_of_sample_particle_loc_, 1, &n_smp_array_[0]);
        time_profile_.decompose_domain_communication += GetWtime() - t0;
        n_smp_disp_array_[0] = 0;
        for (S32 i = 0; i < nproc; i++) {
            n_smp_disp_array_[i + 1] = n_smp_disp_array_[i] + n_smp_array_[i];
        }
        pos_sample_tot_.resizeNoInitialize(n_smp_disp_array_[nproc]);

        t0 = GetWtime();
        comm_info_.allGatherV(pos_sample_loc_.getPointer(), number_of_sample_particle_loc_, pos_sample_tot_.getPointer(), &n_smp_array_[0],
                              &n_smp_disp_array_[0]);
        time_profile_.decompose_domain_communication += GetWtime() - t0;
        // time_profile_.decompose_domain__gather_particle += t1 - time_offset;

        number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];

        // ****************************************************
        // *** decompose domain *******************************
        bool decompose_in_all_processes = true;
        if (myrank == 0 || decompose_in_all_processes) {
            // S32 * istart = new S32[nproc];
            // S32 * iend   = new S32[nproc];
            S32 istart[nproc];
            S32 iend[nproc];
            // --- x direction --------------------------
            // std::sort(pos_sample_tot_, pos_sample_tot_+number_of_sample_particle_tot_, Cmpvec(&F64vec::x));
            t0 = GetWtime();
#if defined(PARTICLE_SIMULATOR_USE_SAMPLE_SORT)
            SampleSortLib::samplesort_bodies(pos_sample_tot_.getPointer(), number_of_sample_particle_tot_,
                                             [](const F64vec &v) -> auto { return v.x; });
#else
            std::sort(pos_sample_tot_.getPointer(), pos_sample_tot_.getPointer(number_of_sample_particle_tot_), Cmpvec(&F64vec::x));
#endif
            time_profile_.decompose_domain_sort_sample += GetWtime() - t0;
            for (S32 i = 0; i < nproc; i++) {
                istart[i] = ((S64)(i) * (S64)(number_of_sample_particle_tot_)) / (S64)(nproc);
                if (i > 0) iend[i - 1] = istart[i] - 1;
            }
            iend[nproc - 1] = number_of_sample_particle_tot_ - 1;
            for (S32 ix = 0; ix < n_domain_[0]; ix++) {
                S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                F64 x0 = 0.0;
                F64 x1 = 0.0;
                // calculateBoundaryOfDomain(number_of_sample_particle_tot_, pos_sample_tot_, 0, istart[ix0], iend[ix1-1], x0, x1);
                // calculateBoundaryOfDomainX(number_of_sample_particle_tot_, pos_sample_tot_.getPointer(), istart[ix0], iend[ix1 - 1], x0, x1);
                calculateBoundaryOfDomain1D<0>(number_of_sample_particle_tot_, pos_sample_tot_.getPointer(), istart[ix0], iend[ix1 - 1], x0, x1);
                for (S32 i = ix0; i < ix1; i++) {
                    pos_domain_temp_[i].low_.x = x0;
                    pos_domain_temp_[i].high_.x = x1;
                }
            }

            // ------------------------------------------
            // --- y direction --------------------------

            PS_OMP_PARALLEL_FOR_G
            for (S32 ix = 0; ix < n_domain_[0]; ix++) {
                S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                // sortCoordinateOfSampleParticle(pos_sample_tot_, istart[ix0], iend[ix1-1], 1);
                // std::sort(pos_sample_tot_+istart[ix0], pos_sample_tot_+(iend[ix1-1]+1), Cmpvec(&F64vec::y));
                // std::sort(&pos_sample_tot_[istart[ix0]], &pos_sample_tot_[(iend[ix1-1]+1)], Cmpvec(&F64vec::y));
                t0 = GetWtime();
                std::sort(pos_sample_tot_.getPointer(istart[ix0]), pos_sample_tot_.getPointer((iend[ix1 - 1] + 1)), Cmpvec(&F64vec::y));
                time_profile_.decompose_domain_sort_sample += GetWtime() - t0;
                /*
                                      for(S32 i=istart[ix0]; i<iend[ix1-1]+1; i++){
                                        if(pos_sample_tot_[i+1].y < pos_sample_tot_[i].y){
                                            std::cout<<"y sort is wrong: i="<<i<<std::endl;
                                            std::cout<<"pos_sample_tot_[i+1]="<<pos_sample_tot_[i+1]<<std::endl;
                                            std::cout<<"pos_sample_tot_[i]="<<pos_sample_tot_[i]<<std::endl;
                                        }
                                    }
                */
                S32 number_of_sample_particle_tot_y = iend[ix1 - 1] - istart[ix0] + 1;
                for (S32 iy = 0; iy < n_domain_[1]; iy++) {
                    S32 iy0 = ix0 + iy * n_domain_[2];
                    S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                    // F64 y0, y1;
                    F64 y0 = 0.0;
                    F64 y1 = 0.0;
                    // std::cout<<"ix="<<ix<<" ix0="<<ix0<<" ix1="<<ix1<<" iy="<<iy<<" iy0="<<iy0<<" iy1="<<iy1<<" y0="<<y0<<" y1="<<y1<<std::endl;
                    // calculateBoundaryOfDomain(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0], 1, istart[iy0]-istart[ix0],
                    // iend[iy1-1]-istart[ix0], y0, y1); calculateBoundaryOfDomainY(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0],
                    // istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1); calculateBoundaryOfDomainY(number_of_sample_particle_tot_y,
                    // &pos_sample_tot_[istart[ix0]], istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
                    // calculateBoundaryOfDomainY(number_of_sample_particle_tot_y, pos_sample_tot_.getPointer(istart[ix0]), istart[iy0] - istart[ix0],
                    // iend[iy1 - 1] - istart[ix0], y0, y1);
                    calculateBoundaryOfDomain1D<1>(number_of_sample_particle_tot_y, pos_sample_tot_.getPointer(istart[ix0]),
                                                   istart[iy0] - istart[ix0], iend[iy1 - 1] - istart[ix0], y0, y1);
                    for (S32 i = iy0; i < iy1; i++) {
                        pos_domain_temp_[i].low_.y = y0;
                        pos_domain_temp_[i].high_.y = y1;
                    }
                }
            }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            // ------------------------------------------
            // --- z direction --------------------------
            PS_OMP_PARALLEL_FOR_G
            for (S32 ix = 0; ix < n_domain_[0]; ix++) {
                S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                for (S32 iy = 0; iy < n_domain_[1]; iy++) {
                    S32 iy0 = ix0 + iy * n_domain_[2];
                    S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                    // sortCoordinateOfSampleParticle(pos_sample_tot_, istart[iy0], iend[iy1-1], 2);
                    // std::sort(pos_sample_tot_+istart[iy0], pos_sample_tot_+(iend[iy1-1]+1), Cmpvec(&F64vec::z));
                    // std::sort(&pos_sample_tot_[istart[iy0]], &pos_sample_tot_[(iend[iy1-1]+1)], Cmpvec(&F64vec::z));
                    t0 = GetWtime();
                    std::sort(pos_sample_tot_.getPointer(istart[iy0]), pos_sample_tot_.getPointer((iend[iy1 - 1] + 1)), Cmpvec(&F64vec::z));
                    time_profile_.decompose_domain_sort_sample += GetWtime() - t0;
                    S32 number_of_sample_particle_tot_z = iend[iy1 - 1] - istart[iy0] + 1;
                    for (S32 iz = 0; iz < n_domain_[2]; iz++) {
                        S32 iz0 = iy0 + iz;
                        // F64 z0, z1;
                        F64 z0 = 0.0;
                        F64 z1 = 0.0;
                        // calculateBoundaryOfDomain(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], 2, istart[iz0]-istart[iy0],
                        // iend[iz0]-istart[iy0], z0, z1); calculateBoundaryOfDomainZ(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0],
                        // istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
                        // calculateBoundaryOfDomainZ(number_of_sample_particle_tot_z, &pos_sample_tot_[istart[iy0]], istart[iz0] - istart[iy0],
                        // iend[iz0] - istart[iy0], z0, z1);
                        calculateBoundaryOfDomain1D<2>(number_of_sample_particle_tot_z, &pos_sample_tot_[istart[iy0]], istart[iz0] - istart[iy0],
                                                       iend[iz0] - istart[iy0], z0, z1);
                        pos_domain_temp_[iz0].low_.z = z0;
                        pos_domain_temp_[iz0].high_.z = z1;
                    }
                }
            }
#endif  // PARTICLE_SIMULATOR_TWO_DIMENSION
        // ------------------------------------------
        // --- process first ------------------------
            if (first_call_by_decomposeDomain) {
                for (S32 i = 0; i < nproc; i++) {
                    pos_domain_[i].low_ = pos_domain_temp_[i].low_;
                    pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                }
            } else {
                for (S32 i = 0; i < nproc; i++) {
                    pos_domain_[i].low_ = (F64)coef_ema_ * pos_domain_temp_[i].low_ + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                    pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                }
            }
            // ------------------------------------------
            // delete [] istart;
            // delete [] iend;
        }
        // F64 wtmp = GetWtime();
        //  ****************************************************
        //  *** broad cast pos_domain_ *************************
        // comm_info_.barrier();
        // F64 wtmp2 = GetWtime();
        t0 = GetWtime();
        if (!decompose_in_all_processes) {
            comm_info_.broadcast(&pos_domain_[0], nproc);
        }
        if (first_call_by_decomposeDomain) {
            first_call_by_decomposeDomain = false;
            //     comm_info_.broadcast(&first_call_by_decomposeDomain, 1);
        }
        // F64 wtmp3 = GetWtime();
        time_profile_.decompose_domain_communication += GetWtime() - t0;

        // std::cout<<"end of bcast: "<<"time: "<<GetWtime() - Tbegin<<std::endl;
        // Comm::broadcast(pos_domain_, nproc);
        //  ****************************************************
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
        pos_domain_[0] = pos_root_domain_;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout << "pos_root_domain_=" << pos_root_domain_ << std::endl;
        std::cout << "pos_domain_[comm_info_.getRank()]=" << pos_domain_[comm_info_.getRank()] << std::endl;
#endif

        pos_sample_tot_.freeMem();
        pos_sample_loc_.freeMem();

#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        // time_profile_.decompose_domain__barrier += wtmp2 - wtmp;
        // time_profile_.decompose_domain__exchange_pos_domain += wtmp3 - wtmp2;
        // time_profile_.decompose_domain__postprocess += GetWtime() - wtmp3;
#endif
        F64 time_end = GetWtime();
        time_profile_.decompose_domain += time_end - time_offset;
        time_profile_.decompose_domain_all += time_end - time_offset;
        //time_profile_.decompose_domain_all__decompose_domain += time_end - time_offset;
    }
    // decomposeDomain
    //////////////////

    ////////////////////
    // decomposeDomainAll

    template <class Tpsys>
    void decomposeDomainAll(Tpsys &psys, const F64 wgh) {
        const bool clear = true;
        collectSampleParticle(psys, clear, wgh);
        decomposeDomain();
    }

    template <class Tpsys>
    void decomposeDomainAll(Tpsys &psys) {
        const F64 wgh = psys.getNumberOfParticleLocal();
        decomposeDomainAll(psys, wgh);
    }

    // not public
    void decomposeDomainAll(const F64ort &pos_my_domain) {
        // This API is added by D.N. to enable restart run with FDPS
        Comm::allGather(&pos_my_domain, 1, &pos_domain_[0]);
        first_call_by_decomposeDomain = true;
    }

    // decomposeDomainAll
    ////////////////////

    void getRootDomain(FILE *fp) {
        fprintf(fp, "%+e %+e %+e\n", pos_root_domain_.low_[0], pos_root_domain_.low_[1], pos_root_domain_.low_[2]);
        fprintf(fp, "%+e %+e %+e\n", pos_root_domain_.high_[0], pos_root_domain_.high_[1], pos_root_domain_.high_[2]);
        return;
    }

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    void getSampleParticleLocal(FILE *fp) {
        for (S32 i = 0; i < number_of_sample_particle_loc_; i++) {
            fprintf(fp, "%+e %+e %+e\n", pos_sample_loc_[i].x, pos_sample_loc_[i].y, pos_sample_loc_[i].z);
        }
        return;
    }

    void getSampleParticleTotal(FILE *fp) {
        for (S32 i = 0; i < number_of_sample_particle_tot_; i++) {
            fprintf(fp, "%+e %+e %+e\n", pos_sample_tot_[i].x, pos_sample_tot_[i].y, pos_sample_tot_[i].z);
        }
        return;
    }
#endif

    void getPosDomainTotal(FILE *fp) {
        auto nproc = comm_info_.getNumberOfProc();
        for (S32 i = 0; i < nproc; i++) {
            for (S32 k = 0; k < DIMENSION; k++) fprintf(fp, "%+e ", pos_domain_[i].low_[k]);
            for (S32 k = 0; k < DIMENSION; k++) fprintf(fp, "%+e ", pos_domain_[i].high_[k]);
            fprintf(fp, "\n");
        }
    }

    S32 *getPointerOfNDomain() { return n_domain_; };

    F64ort *getPointerOfPosDomain() { return &pos_domain_[0]; };

    // A. Tanikawa need this method for Particle Mesh...
    S32 getNDomain(const S32 dim) const { return n_domain_[dim]; }
    S32 getRank1D(const S32 rank_glb, const S32 dim) const {
        S32 rank_1d;
        if (dim == 0)
            rank_1d = rank_glb / (n_domain_[1] * n_domain_[2]);
        else if (dim == 1)
            rank_1d = (rank_glb / n_domain_[2]) % n_domain_[1];
        else if (dim == 2)
            rank_1d = rank_glb % n_domain_[2];
        return rank_1d;
    }

    /* AT_DEBUG
    // for DEBUG
    F32vec & getPosSample(const S32 id=0){return pos_sample_tot_[id];}
    // for DEBUG // A. Tanikawa would not like to delete this method...
    F32ort & getPosDomain(const S32 id=0) const {return pos_domain_[id];}
    // for DEBUG
    void setPosDomain(const S32 id, const F32ort & pos){ pos_domain_[id] = pos;}
    */
    // for DEBUG
    F64vec &getPosSample(const S32 id = 0) { return pos_sample_tot_[id]; }
    // for DEBUG // A. Tanikawa would not like to delete this method...
    // F64ort & getPosDomain(const S32 id=0) const {return pos_domain_[id];}
    F64ort getPosDomain(const S32 id = 0) const { return pos_domain_[id]; }
    F64ort getPosDomain(const F64ort &box, const S32 id = 0) const {
        F64ort ret = pos_domain_[id];
        if (periodic_axis_[0] == false) {
            if (getRank1D(id, 0) == 0) {
                ret.low_.x = box.low_.x;
            }
            if (getRank1D(id, 0) == n_domain_[0] - 1) {
                ret.high_.x = box.high_.x;
            }
        }
        if (periodic_axis_[1] == false) {
            if (getRank1D(id, 1) == 0) {
                ret.low_.y = box.low_.y;
            }
            if (getRank1D(id, 1) == n_domain_[1] - 1) {
                ret.high_.y = box.high_.y;
            }
        }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        if (periodic_axis_[2] == false) {
            if (getRank1D(id, 2) == 0) {
                ret.low_.z = box.low_.z;
            }
            if (getRank1D(id, 2) == n_domain_[2] - 1) {
                ret.high_.z = box.high_.z;
            }
        }
#endif
        return ret;
    }
    // for DEBUG
    void setPosDomain(const S32 id, const F64ort &pos) { pos_domain_[id] = pos; }

    void setBoundaryCondition(enum BOUNDARY_CONDITION bc) {
        boundary_condition_ = bc;
        if (DIMENSION == 2 && (bc == BOUNDARY_CONDITION_PERIODIC_XYZ || bc == BOUNDARY_CONDITION_PERIODIC_XZ ||
                               bc == BOUNDARY_CONDITION_PERIODIC_YZ || bc == BOUNDARY_CONDITION_PERIODIC_Z)) {
            throw "PS_ERROR: in setBoundaryCondition(enum BOUNDARY_CONDITION) \n boundary condition is incompatible with DIMENSION";
        }
        if (bc == BOUNDARY_CONDITION_OPEN)
            periodic_axis_[0] = periodic_axis_[1] = periodic_axis_[2] = false;
        else if (bc == BOUNDARY_CONDITION_PERIODIC_X) {
            periodic_axis_[0] = true;
            periodic_axis_[1] = periodic_axis_[2] = false;
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_Y) {
            periodic_axis_[1] = true;
            periodic_axis_[0] = periodic_axis_[2] = false;
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_Z) {
            periodic_axis_[2] = true;
            periodic_axis_[0] = periodic_axis_[1] = false;
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XY) {
            periodic_axis_[0] = periodic_axis_[1] = true;
            periodic_axis_[2] = false;
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XZ) {
            periodic_axis_[0] = periodic_axis_[2] = true;
            periodic_axis_[1] = false;
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_YZ) {
            periodic_axis_[1] = periodic_axis_[2] = true;
            periodic_axis_[0] = false;
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
            periodic_axis_[0] = periodic_axis_[1] = periodic_axis_[2] = true;
    }

    S32 getBoundaryCondition() const { return boundary_condition_; }

    void setPosRootDomain1D(const F64 l, const F64 h, const S32 cid) {
        pos_root_domain_.low_[cid] = l;
        pos_root_domain_.high_[cid] = h;
    }

    void setPosRootDomainX(const F64 l, const F64 h) {
        setPosRootDomain1D(l, h, 0);
        flag_set_pos_root_domain_[0] = true;
    }
    void setPosRootDomainY(const F64 l, const F64 h) {
        setPosRootDomain1D(l, h, 1);
        flag_set_pos_root_domain_[1] = true;
    }
    void setPosRootDomainZ(const F64 l, const F64 h) {
        setPosRootDomain1D(l, h, 2);
        flag_set_pos_root_domain_[2] = true;
    }
    void setPosRootDomain(const F64vec &low, const F64vec &high) {
        for (S32 k = 0; k < DIMENSION; k++) {
            if (low[k] > high[k]) {
                PARTICLE_SIMULATOR_PRINT_ERROR("The coodinate of the root domain is inconsistent.");
                std::cerr << "The coordinate of the low vertex of the rood domain=" << low << std::endl;
                std::cerr << "The coordinate of the high vertex of the rood domain=" << high << std::endl;
                Abort(-1);
            }
        }
        for (S32 i = 0; i < DIMENSION; i++) {
            // if( periodic_axis_[i] == false ) continue;
            if (low[i] < high[i]) {
                pos_root_domain_.low_[i] = low[i];
                pos_root_domain_.high_[i] = high[i];
                flag_set_pos_root_domain_[i] = true;
            }
        }
    }
    void setPosRootDomain(const F64ort &box) { setPosRootDomain(box.low_, box.high_); }
    bool getFlagSetPosRootDomain(int cid) const { return flag_set_pos_root_domain_[cid]; }

    const F64ort getPosRootDomain() const { return pos_root_domain_; }
    const F64vec getLenRootDomain() const { return pos_root_domain_.getFullLength(); }

    void getPeriodicAxis(bool pa[]) const {
        for (S32 i = 0; i < DIMENSION; i++) pa[i] = periodic_axis_[i];
    }

    template <class Tpsys>
    bool checkCollectSampleParticleSubset(Tpsys &psys);
    template <class Tpsys>
    bool checkCollectSampleParticleAverage(Tpsys &psys);
    template <class Tpsys>
    bool checkDecomposeDomain(Tpsys &psys);

    S32 whereToGoOpen(const F64 x,
                      const F64 coord_sep[],  // size is n_proc+1
                      const S32 my_rank, const S32 shift = 1) {
        if (coord_sep[my_rank] <= x && x < coord_sep[my_rank + 1]) {
            return my_rank;
        }
        auto rank = my_rank;
        if (x < coord_sep[rank]) {
            do {
                rank -= shift;
            } while (x < coord_sep[rank]);
        } else if (x >= coord_sep[rank + 1]) {
            do {
                rank += shift;
            } while (x >= coord_sep[rank + 1]);
        }
        return rank;
    }

    S32 whereToGoPeriodic(const F64 x,
                          const F64 coord_sep[],  // size is n_proc+1
                          const S32 my_rank, const F64 len_peri, const S32 n_proc, const S32 shift = 1) {
        if (coord_sep[my_rank] <= x && x < coord_sep[my_rank + 1]) {
            return my_rank;
        }
        auto rank = my_rank;
        auto x_tmp = (x < coord_sep[rank]) ? (coord_sep[rank] - x < x + len_peri - coord_sep[rank + 1] ? x : x + len_peri)
                                           : (x - coord_sep[rank + 1] < coord_sep[rank] - (x - len_peri) ? x : x - len_peri);
        if (x_tmp < coord_sep[rank]) {
            if (x_tmp != x)
                rank = n_proc - 1;
            else
                rank -= shift;
            while (x < coord_sep[rank]) {
                rank -= shift;
            };
        } else if (x_tmp >= coord_sep[rank + 1]) {
            if (x_tmp != x)
                rank = 0;
            else
                rank += shift;
            while (x >= coord_sep[rank + 1]) {
                rank += shift;
            };
        }
        return rank;
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL

    void determineRoughPartition(ReallocatableArray<F64> &coord,  // in. not const, becuase the order is temporally changed in Randomsampling.
                                 ReallocatableArray<F64> &coord_partition,  // out
                                 const S32 smp_frac,                        // in
                                 CommInfo &comm_info, bool print = false) {
        const auto my_rank = comm_info.getRank();
        const auto n_proc = comm_info.getNumberOfProc();
        auto n_coord = coord.size();
        auto n_coord_smp = ((n_coord / smp_frac) >= 1) ? n_coord / smp_frac : 1;
        if (n_coord <= 0) n_coord_smp = 0;

        // comm_info.barrier();
        // if(comm_info.getRank()==0) std::cerr<<"check aaa"<<std::endl;

        ReallocatableArray<F64> coord_smp(n_coord_smp, n_coord_smp, MemoryAllocMode::Pool);
        RandomSampling(coord, coord_smp, [](F64 &src, F64 &dst) { dst = src; });

        // comm_info.barrier();
        // if(comm_info.getRank()==0) std::cerr<<"check bbb"<<std::endl;

        ReallocatableArray<S32> n_coord_smp_ar(n_proc, n_proc, MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_disp_coord_smp_ar(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);

        // comm_info.barrier();
        // if(comm_info.getRank()==0) std::cerr<<"check ccc"<<std::endl;

        comm_info.gather(&n_coord_smp, 1, n_coord_smp_ar.getPointer());
        if (my_rank == 0) {
            n_disp_coord_smp_ar[0] = 0;
            for (S32 i = 0; i < n_proc; i++) {
                n_disp_coord_smp_ar[i + 1] = n_disp_coord_smp_ar[i] + n_coord_smp_ar[i];
            }
        }

        // comm_info.barrier();
        // if(comm_info.getRank()==0) std::cerr<<"check ddd"<<std::endl;

        const auto n_coord_smp_glb = (my_rank == 0) ? n_disp_coord_smp_ar[n_proc] : 0;
        ReallocatableArray<F64> coord_smp_glb(n_coord_smp_glb, n_coord_smp_glb, MemoryAllocMode::Pool);
        comm_info.gatherV(coord_smp.getPointer(), n_coord_smp, coord_smp_glb.getPointer(), n_coord_smp_ar.getPointer(),
                          n_disp_coord_smp_ar.getPointer());
        if (my_rank == 0) {
            std::sort(coord_smp_glb.frontP(), coord_smp_glb.endP());
            for (S32 i = 1; i < n_proc; i++) {
                const auto id = (n_coord_smp_glb / n_proc) * i;
                if (id == 0 || id == n_coord_smp_glb - 1)
                    coord_partition[i] = coord_smp_glb[id];
                else
                    coord_partition[i] = (coord_smp_glb[id] + coord_smp_glb[id + 1]) * 0.5;
            }
        }

        // comm_info.barrier();
        // if(comm_info.getRank()==0) std::cerr<<"check eee"<<std::endl;

        comm_info.broadcast(coord_partition.getPointer(), n_proc + 1, 0);

        // comm_info.barrier();
        // if(comm_info.getRank()==0) std::cerr<<"check fff"<<std::endl;
    }

    template <typename T, typename Tgetter>
    void exchangeSample(ReallocatableArray<T> &smp,                // in and out
                        ReallocatableArray<F64> &coord_partition,  // in
                        const S32 n_proc,                          // in
                        const S32 my_rank,                         // in
                        const F64 len_peri,                        // in
                        const bool flag_periodic,                  // in
                        CommInfo comm_info,                        // in
                        Tgetter getter) {
        const auto n_smp = smp.size();
        ReallocatableArray<S32> n_send(n_proc, n_proc, MemoryAllocMode::Pool);
        ReallocatableArray<S32> rank_send(n_smp, n_smp, MemoryAllocMode::Pool);
        ReallocatableArray<MPI_Request> req_recv(n_proc, n_proc, MemoryAllocMode::Pool);
        ReallocatableArray<MPI_Status> status(n_proc, n_proc, MemoryAllocMode::Pool);
        for (S32 i = 0; i < n_proc; i++) n_send[i] = 0;
        S32 dnp_max = 0;
        S32 i_loc = n_smp - 1;

        if (!flag_periodic) {
            for (S32 i = n_smp - 1; i >= 0; i--) {
                const auto rank = whereToGoOpen(getter(smp[i]), coord_partition.getPointer(), my_rank);
                rank_send[i] = rank;
                n_send[rank]++;
                if (rank != my_rank) {
                    std::swap(smp[i], smp[i_loc]);
                    rank_send[i_loc] = rank;
                    i_loc--;
                }
                const auto dnp = std::abs(rank - my_rank);
                dnp_max = std::max(dnp_max, dnp);
            }
        } else {
            for (S32 i = n_smp - 1; i >= 0; i--) {
                const auto rank = whereToGoPeriodic(getter(smp[i]), coord_partition.getPointer(), my_rank, len_peri, n_proc);
                rank_send[i] = rank;
                n_send[rank]++;
                if (rank != my_rank) {
                    std::swap(smp[i], smp[i_loc]);
                    rank_send[i_loc] = rank;
                    i_loc--;
                }
                const auto dnp = std::abs(rank - my_rank);
                dnp_max = std::max(dnp_max, dnp);
            }
        }
        const auto n_remain = i_loc + 1;

        if (DEBUG_FLAG_DD) return;

        ReallocatableArray<S32> n_disp_send(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
        n_disp_send[0] = 0;
        for (S32 i = 0; i < n_proc; i++) {
            n_disp_send[i + 1] = n_disp_send[i] + n_send[i];
            n_send[i] = 0;  // reuse
        }

        ReallocatableArray<T> smp_send(n_disp_send[n_proc], n_disp_send[n_proc], MemoryAllocMode::Pool);
        for (S32 i = n_remain; i < n_smp; i++) {
            const auto rank = rank_send[i];
            const auto loc = n_send[rank] + n_disp_send[rank];
            smp_send[loc] = smp[i];
            n_send[rank]++;
        }

        const auto dnp_max_glb = comm_info.getMaxValue(dnp_max);
        // std::cerr<<"my_rank= "<<my_rank
        //          <<" dnp_max_glb= "<<dnp_max_glb
        //          <<std::endl;
        S32 dnp_head = -dnp_max_glb;
        S32 dnp_tail = dnp_max_glb;
        if (dnp_max_glb * 2 + 1 > n_proc) {
            dnp_head = -n_proc / 2;
            dnp_tail = (n_proc % 2 == 0) ? (n_proc / 2 - 1) : n_proc / 2;
        }
        // std::cerr<<"Comm::getRank()= "<<Comm::getRank()
        //          <<" dnp_head= "<<dnp_head
        //          <<" dnp_tail= "<<dnp_tail
        //          <<std::endl;

        // comm_info.barrier();
        // std::cerr<<"check dd"<<std::endl;
#if 1
        // if(2*dnp_max_glb < n_proc){
        if (1) {
            ReallocatableArray<S32> rank_comm(n_proc, 0, MemoryAllocMode::Pool);
            ReallocatableArray<S32> n_send_tmp(n_proc, 0, MemoryAllocMode::Pool);
            ReallocatableArray<S32> n_recv_tmp(n_proc, 0, MemoryAllocMode::Pool);
            for (S32 i = dnp_head; i <= dnp_tail; i++) {
                auto rank = (my_rank + n_proc + i) % n_proc;
                if (my_rank == rank) continue;
                rank_comm.push_back(rank);
                n_send_tmp.push_back(n_send[rank]);
                n_recv_tmp.push_back(0);
            }

            // comm_info.barrier();
            // std::cerr<<"check ee"<<std::endl;

            comm_info.sendIrecv(n_send_tmp.getPointer(), rank_comm.getPointer(), 1, n_send_tmp.size(), n_recv_tmp.getPointer(),
                                rank_comm.getPointer(), 1, n_recv_tmp.size());

            // comm_info.barrier();
            // std::cerr<<"check ff"<<std::endl;

            ReallocatableArray<S32> rank_send_tmp(n_proc, 0, MemoryAllocMode::Pool);
            ReallocatableArray<S32> rank_recv_tmp(n_proc, 0, MemoryAllocMode::Pool);
            ReallocatableArray<S32> n_send_tmp2(n_proc, 0, MemoryAllocMode::Pool);
            ReallocatableArray<S32> n_recv_tmp2(n_proc, 0, MemoryAllocMode::Pool);
            for (S32 i = 0; i < rank_comm.size(); i++) {
                const auto rank = rank_comm[i];
                if (n_send_tmp[i] > 0) {
                    rank_send_tmp.push_back(rank);
                    n_send_tmp2.push_back(n_send_tmp[i]);
                }
                if (n_recv_tmp[i] > 0) {
                    rank_recv_tmp.push_back(rank);
                    n_recv_tmp2.push_back(n_recv_tmp[i]);
                }
            }

            // comm_info.barrier();
            // std::cerr<<"check gg"<<std::endl;

            ReallocatableArray<S32> n_disp_send_tmp2(n_send_tmp2.size() + 1, n_send_tmp2.size() + 1, MemoryAllocMode::Pool);
            ReallocatableArray<S32> n_disp_recv_tmp2(n_recv_tmp2.size() + 1, n_recv_tmp2.size() + 1, MemoryAllocMode::Pool);
            for (S32 i = 0; i < rank_send_tmp.size(); i++) {
                const auto rank = rank_send_tmp[i];
                n_disp_send_tmp2[i] = n_disp_send[rank];
            }
            // Note that n_disp_send_tmp2[i] must be the beginning array index
            // of the send data for RANK rank_send_tmp[i]. Since rank_send_tmp[]
            // is not necessarily in ascending order (see above), we cannot use
            // function PrefixSum to calculate it as below:
            //     PrefixSum(n_send_tmp2, n_disp_send_tmp2); // WRONG
            PrefixSum(n_recv_tmp2, n_disp_recv_tmp2);
            /*
            comm_info.barrier();
            std::cerr<<"check hh"<<std::endl;
            std::cerr<<"rank_send_tmp.size()= "<<rank_send_tmp.size()<<std::endl;
            for(int i=0; i<rank_send_tmp.size(); i++){
                std::cerr<<"rank_send_tmp[i]= "<<rank_send_tmp[i]<<std::endl;
            }
            std::cerr<<"rank_recv_tmp.size()= "<<rank_recv_tmp.size()<<std::endl;
            for(int i=0; i<rank_recv_tmp.size(); i++){
                std::cerr<<"rank_recv_tmp[i]= "<<rank_recv_tmp[i]<<std::endl;
            }
            */
            smp.resizeNoInitialize(n_disp_recv_tmp2.back() + n_remain);
            comm_info.sendIrecvV(smp.getPointer(), rank_send_tmp.getPointer(), n_send_tmp2.getPointer(), n_disp_send_tmp2.getPointer(),
                                 n_send_tmp2.size(), smp.getPointer(n_remain), rank_recv_tmp.getPointer(), n_recv_tmp2.getPointer(),
                                 n_disp_recv_tmp2.getPointer(), n_recv_tmp2.size());

            // comm_info.barrier();
            // std::cerr<<"check ii"<<std::endl;

        } else {
            ReallocatableArray<S32> n_recv(n_proc, n_proc, MemoryAllocMode::Pool);
            for (S32 i = 0; i < n_proc; i++) {
                n_recv[i] = 0;
            }
            comm_info.allToAll(n_send.getPointer(), 1, n_recv.getPointer());
            ReallocatableArray<S32> n_disp_recv(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
            n_disp_recv[0] = 0;
            for (S32 i = 0; i < n_proc; i++) {
                n_disp_recv[i + 1] = n_disp_recv[i] + n_recv[i];
            }
            smp.resizeNoInitialize(n_disp_recv[n_proc] + n_remain);
            comm_info.allToAllV(smp_send.getPointer(), n_send.getPointer(), n_disp_send.getPointer(), smp.getPointer(n_remain), n_recv.getPointer(),
                                n_disp_recv.getPointer());
        }

#else
        ReallocatableArray<S32> n_recv(n_proc, n_proc, MemoryAllocMode::Pool);
        for (S32 i = 0; i < n_proc; i++) {
            n_recv[i] = 0;
        }
        // if(2*dnp_max_glb < n_proc){
        if (1) {
            S32 n_cnt = 0;
            for (S32 i = dnp_head; i <= dnp_tail; i++) {
                auto rank = (my_rank + n_proc + i) % n_proc;
                if (my_rank == rank) continue;
                MPI_Irecv(n_recv.getPointer(rank), 1, GetDataType<S32>(), rank, 0, comm_info.getCommunicator(), req_recv.getPointer(n_cnt));
                n_cnt++;
            }
            for (S32 i = dnp_head; i <= dnp_tail; i++) {
                auto rank = (my_rank + n_proc + i) % n_proc;
                if (my_rank == rank) continue;
                MPI_Send(n_send.getPointer(rank), 1, GetDataType<S32>(), rank, 0, comm_info.getCommunicator());
            }
            MPI_Waitall(n_cnt, req_recv.getPointer(), status.getPointer());
        } else {
            // use all to all
            comm_info.allToAll(n_send.getPointer(), 1, n_recv.getPointer());
        }

        ReallocatableArray<S32> n_disp_recv(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
        n_disp_recv[0] = 0;
        for (S32 i = 0; i < n_proc; i++) n_disp_recv[i + 1] = n_disp_recv[i] + n_recv[i];
        smp.resizeNoInitialize(n_disp_recv[n_proc] + n_remain);
        // if(2*dnp_max_glb < n_proc){
        if (1) {
            S32 n_cnt = 0;
            for (S32 i = dnp_head; i <= dnp_tail; i++) {
                const auto rank = (my_rank + n_proc + i) % n_proc;
                const auto adr = n_disp_recv[rank] + n_remain;
                MPI_Irecv(smp.getPointer(adr), n_recv[rank], GetDataType<T>(), rank, 1, comm_info.getCommunicator(), req_recv.getPointer(n_cnt));
                n_cnt++;
            }
            for (S32 i = dnp_head; i <= dnp_tail; i++) {
                const auto rank = (my_rank + n_proc + i) % n_proc;
                const auto adr = n_disp_send[rank];
                MPI_Send(smp_send.getPointer(adr), n_send[rank], GetDataType<T>(), rank, 1, comm_info.getCommunicator());
            }
            MPI_Waitall(n_cnt, req_recv.getPointer(), status.getPointer());
        } else {
            // alltoallv
            comm_info.allToAllV(smp_send.getPointer(), n_send.getPointer(), n_disp_send.getPointer(), smp.getPointer(n_remain), n_recv.getPointer(),
                                n_disp_recv.getPointer());
        }
#endif
        // comm_info.barrier();
        // std::cerr<<"check zz"<<std::endl;
    }

    void determinePartition(ReallocatableArray<F64> &coord,                // in
                            ReallocatableArray<F64> &coord_partition_tmp,  // in
                            ReallocatableArray<F64> &coord_partition,      // out
                            CommInfo &comm_info,                           // in
                            bool print = false) {
        const auto my_rank = comm_info.getRank();
        const auto n_proc = comm_info.getNumberOfProc();
        const auto n_smp_loc = coord.size();
        ReallocatableArray<S32> n_smp(n_proc, n_proc, MemoryAllocMode::Pool);
        comm_info.allGather(&n_smp_loc, 1, n_smp.getPointer());
        // std::cerr<<"my_rank= "<<my_rank<<std::endl;
        // for(S32 i=0; i<n_proc; i++){
        //     std::cerr<<"i= "<<i<<" n_smp[i]= "<<n_smp[i]<<std::endl;
        // }

        ReallocatableArray<S32> n_disp_smp(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
        n_disp_smp[0] = 0;
        for (S32 i = 0; i < n_proc; i++) {
            n_disp_smp[i + 1] = n_disp_smp[i] + n_smp[i];
        }
        const auto n_tot = n_disp_smp[n_proc];
        const auto n_ave = n_tot / n_proc;
        const auto n_amari = n_tot % n_proc;
        // std::cerr<<"n_tot= "<<n_tot
        //          <<" n_ave= "<<n_ave
        //          <<" n_amari= "<<n_amari
        //          <<std::endl;
        ReallocatableArray<F64> coord_partition_send(n_proc - 1, 0, MemoryAllocMode::Pool);  // whitout edge
        ReallocatableArray<S32> n_partition_recv(n_proc, n_proc, MemoryAllocMode::Pool);
        for (S32 i = 0; i < n_proc; i++) n_partition_recv[i] = 0;
        S32 i_domain = 0;
        S32 offset = n_disp_smp[my_rank];

        for (S32 i = 1; i < n_proc; i++) {
            const auto index_partition = n_ave * i + ((i < n_amari) ? 1 : 0);
            while (!(n_disp_smp[i_domain] <= index_partition && index_partition < n_disp_smp[i_domain + 1])) {
                i_domain++;
            }
            if (my_rank == i_domain) {
                const auto head = n_disp_smp[i_domain];
                const auto tail = n_disp_smp[i_domain + 1] - 1;
                if (index_partition == head) {
                    coord_partition_send.push_back(coord_partition_tmp[i_domain]);
                } else if (index_partition == tail) {
                    coord_partition_send.push_back(coord_partition_tmp[i_domain + 1]);
                } else {
                    coord_partition_send.push_back((coord[index_partition - offset] + coord[index_partition - offset + 1]) * 0.5);
                }
            }
            n_partition_recv[i_domain]++;
        }
        // n_partition_recv[my_rank] = coord_partition_send.size();

        auto n_partition_send = n_partition_recv[my_rank];
        ReallocatableArray<S32> n_disp_partition_recv(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
        n_disp_partition_recv[0] = 0;
        for (S32 i = 0; i < n_proc; i++) {
            n_disp_partition_recv[i + 1] = n_disp_partition_recv[i] + n_partition_recv[i];
        }
        comm_info.gatherV(coord_partition_send.getPointer(), n_partition_send, coord_partition.getPointer(1), n_partition_recv.getPointer(),
                          n_disp_partition_recv.getPointer());
        // NOTE: in the 3rd argument, coord_partition[1] means we do not calculate edge coordinate.
    }

    template <typename T0, typename T1, typename Tgetter>
    void gatherSampleRandomly(ReallocatableArray<T0> &smp,      // in
                              ReallocatableArray<T1> &smp_gat,  // out
                              const S32 n_smp_gat,              // in
                              CommInfo comm_info,               // in
                              Tgetter getter) {
        // comm_info.barrier();
        // if(comm_info.getRank()==0){std::cerr<<"CHECK E0"<<std::endl;}

        const auto my_rank = comm_info.getRank();
        const auto n_proc = comm_info.getNumberOfProc();
        const auto n_smp = smp.size();
        ReallocatableArray<S32> n_smp_ar(n_proc, n_proc, MemoryAllocMode::Pool);
        ReallocatableArray<S32> n_disp_smp_ar(n_proc + 1, n_proc + 1, MemoryAllocMode::Pool);
        comm_info.allGather(&n_smp, 1, n_smp_ar.getPointer());
        n_disp_smp_ar[0] = 0;
        for (S32 i = 0; i < n_proc; i++) {
            n_disp_smp_ar[i + 1] = n_disp_smp_ar[i] + n_smp_ar[i];
        }
        const auto n_smp_glb = n_disp_smp_ar[n_proc];
        std::mt19937 mt(1);
        std::uniform_int_distribution<S32> dist(0, n_smp_glb - 1);
        ReallocatableArray<S32> index(n_smp_gat, n_smp_gat, MemoryAllocMode::Pool);
        for (S32 i = 0; i < n_smp_gat; i++) {
            S32 j = dist(mt);
            index[i] = j;
        }
        S32 n_smp_send = 0;
        for (S32 i = 0; i < n_smp_gat; i++) {
            if (n_disp_smp_ar[my_rank] <= index[i] && n_disp_smp_ar[my_rank + 1] > index[i]) {
                n_smp_send++;
            }
        }

        // comm_info.barrier();
        // if(comm_info.getRank()==0){std::cerr<<"CHECK E1"<<std::endl;}

        n_smp_send = std::min(n_smp_send, n_smp);
        const auto n_smp_send_glb = comm_info.getSum(n_smp_send);
        assert(n_smp_send_glb == n_smp_gat);
        ReallocatableArray<T1> smp_send(n_smp_send, n_smp_send, MemoryAllocMode::Pool);
        RandomSampling(smp, smp_send, getter);

        // comm_info.barrier();
        // if(comm_info.getRank()==0){std::cerr<<"CHECK E2"<<std::endl;}

        // reuse n_smp_ar and n_disp_smp_ar
        comm_info.allGather(&n_smp_send, 1, n_smp_ar.getPointer());
        n_disp_smp_ar[0] = 0;
        for (S32 i = 0; i < n_proc; i++) {
            n_disp_smp_ar[i + 1] = n_disp_smp_ar[i] + n_smp_ar[i];
        }

        // comm_info.barrier();
        // if(comm_info.getRank()==0){std::cerr<<"CHECK E3"<<std::endl;}

        comm_info.gatherV(smp_send.getPointer(), n_smp_send, smp_gat.getPointer(), n_smp_ar.getPointer(), n_disp_smp_ar.getPointer());

        // comm_info.barrier();
        // if(comm_info.getRank()==0){std::cerr<<"CHECK E4"<<std::endl;}
    }

    void setEdgeBoundary(F64 &x_low, F64 &x_high, const bool flag_periodic, const S32 cid) {
        if (!flag_periodic) {
            x_low = -LARGE_DOUBLE;
            x_high = LARGE_DOUBLE;
        } else {
            x_low = pos_root_domain_.low_[cid];
            x_high = pos_root_domain_.high_[cid];
        }
    }
#endif

    void decomposeDomainMultiStep2(const bool flag_use_old_domain = false) {
        F64 time_offset = GetWtime();

        const auto n_proc = comm_info_.getNumberOfProc();

        if (n_proc == 1) {
            pos_domain_[0] = pos_root_domain_;
            pos_sample_tot_.freeMem();
            pos_sample_loc_.freeMem();
            time_profile_.decompose_domain = GetWtime() - time_offset;
            return;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const auto full_len_root_domain = pos_root_domain_.high_ - pos_root_domain_.low_;
        const auto my_rank = comm_info_.getRank();
        const auto my_rank_x = comm_info_1d_[0].getRank();
        const auto n_proc_x = comm_info_1d_[0].getNumberOfProc();
        const auto my_rank_y = comm_info_1d_[1].getRank();
        const auto n_proc_y = comm_info_1d_[1].getNumberOfProc();
        const auto my_rank_z = comm_info_1d_[2].getRank();
        const auto n_proc_z = comm_info_1d_[2].getNumberOfProc();

        ////////////////
        /// X DIRECTION
        // Comm::barrier(); if(Comm::getRank()==0) std::cerr<<"start x"<<std::endl;
        // const auto n_smp_sub_x = Comm::getSum(pos_sample_loc_.size(), comm_sub_[0]) / (n_proc_y*n_proc_z) ;
        const auto n_smp_sub_x = comm_info_sub_[0].getSum(pos_sample_loc_.size()) / (n_proc_y * n_proc_z);
        ReallocatableArray<F64> x_smp_sub_in_slab(n_smp_sub_x, n_smp_sub_x, MemoryAllocMode::Pool);

        gatherSampleRandomly(pos_sample_loc_, x_smp_sub_in_slab, n_smp_sub_x, comm_info_sub_[0], [](F64vec &src, F64 &dst) { dst = src.x; });

        // calculate x_partition
        ReallocatableArray<F64> x_partition(n_proc_x + 1, n_proc_x + 1, MemoryAllocMode::Pool);

        if (comm_info_sub_[0].getRank() == 0) {
            ReallocatableArray<F64> x_partition_tmp(n_proc_x + 1, n_proc_x + 1, MemoryAllocMode::Pool);
            if (flag_use_old_domain) {
                // reuse old domain
                x_partition_tmp[0] = pos_domain_[0].low_.x;
                for (S32 i = 1; i < n_proc_x + 1; i++) {
                    const auto rank_src = (i - 1) * n_proc_y * n_proc_z;
                    x_partition_tmp[i] = pos_domain_[rank_src].high_.x;
                }
            } else {
                // rougly calc postion of partition in x direction
                determineRoughPartition(x_smp_sub_in_slab, x_partition_tmp, 2 * n_proc_x, comm_info_1d_[0], false);
                setEdgeBoundary(x_partition_tmp[0], x_partition_tmp[n_proc_x], periodic_axis_[0], 0);
            }
            exchangeSample(x_smp_sub_in_slab, x_partition_tmp, n_proc_x, my_rank_x, full_len_root_domain[0], periodic_axis_[0], comm_info_1d_[0],
                           [](const F64 &src) -> F64 { return src; });
            std::sort(x_smp_sub_in_slab.frontP(), x_smp_sub_in_slab.endP());
#ifdef PS_DEBUG_DD
            for (S32 i = 1; i < x_smp_sub_in_slab.size(); i++) {
                assert(x_partition_tmp[my_rank_x] <= x_smp_sub_in_slab[i]);
                assert(x_partition_tmp[my_rank_x + 1] > x_smp_sub_in_slab[i]);
                assert(x_smp_sub_in_slab[i - 1] < x_smp_sub_in_slab[i]);
            }
#endif
            determinePartition(x_smp_sub_in_slab, x_partition_tmp, x_partition, comm_info_1d_[0]);
            setEdgeBoundary(x_partition[0], x_partition[n_proc_x], periodic_axis_[0], 0);
        }

        ////////////////////////////////////////
        // exchange samples in x directino for determination of y partition
        // broadcast x_partition to all processes
        // Comm::broadcast(x_partition.getPointer(), x_partition.size());
        comm_info_.broadcast(x_partition.getPointer(), x_partition.size());
        // exchange samples so that each particles belong to own apropriate domain.
        exchangeSample(pos_sample_loc_, x_partition, n_proc_x, my_rank_x, full_len_root_domain[0], periodic_axis_[0], comm_info_1d_[0],
                       [](const F64vec &src) -> F64 { return src.x; });
#ifdef PS_DEBUG_DD
        for (S32 i = 0; i < pos_sample_loc_.size(); i++) {
            assert(pos_sample_loc_[i].x < x_partition[my_rank_x + 1]);
            assert(pos_sample_loc_[i].x >= x_partition[my_rank_x]);
        }
#endif

        // comm_info_.barrier();
        // if(comm_info_.getRank()==0){std::cerr<<"CHECK D"<<std::endl;}

        ////////////////
        /// Y DIRECTION
        // Comm::barrier(); if(Comm::getRank()==0) std::cerr<<"start y"<<std::endl;
        // const auto n_smp_sub_y = Comm::getSum(pos_sample_loc_.size(), comm_1d_[2])/n_proc_z;
        const auto n_smp_sub_y = comm_info_1d_[2].getSum(pos_sample_loc_.size()) / n_proc_z;

        ReallocatableArray<F64> y_smp_sub_in_rod(n_smp_sub_y, n_smp_sub_y, MemoryAllocMode::Pool);

        // comm_info_1d_[2].barrier();
        // if(comm_info_1d_[2].getRank()==0){std::cerr<<"CHECK E "<<Comm::getRank()<<std::endl;}

        gatherSampleRandomly(pos_sample_loc_, y_smp_sub_in_rod, n_smp_sub_y, comm_info_1d_[2], [](F64vec &src, F64 &dst) { dst = src.y; });

        // comm_info_1d_[2].barrier();
        // if(comm_info_1d_[2].getRank()==0){std::cerr<<"CHECK F "<<Comm::getRank()<<std::endl;}

        // calc y_partition
        ReallocatableArray<F64> y_partition(n_proc_y + 1, n_proc_y + 1, MemoryAllocMode::Pool);

        if (my_rank_z == 0) {
            ReallocatableArray<F64> y_partition_tmp(n_proc_y + 1, n_proc_y + 1, MemoryAllocMode::Pool);
            if (flag_use_old_domain) {
                // reuse old domain
                y_partition_tmp[0] = pos_domain_[my_rank_x * n_proc_y * n_proc_z].low_.y;
                for (S32 i = 1; i < n_proc_y + 1; i++) {
                    const auto rank_src = my_rank_x * n_proc_y * n_proc_z + (i - 1) * n_proc_z;
                    y_partition_tmp[i] = pos_domain_[rank_src].high_.y;
                }
            } else {
                // rougly calc postion of partition in y direction
                determineRoughPartition(y_smp_sub_in_rod, y_partition_tmp, 2 * n_proc_y, comm_info_1d_[1], true);
                setEdgeBoundary(y_partition_tmp[0], y_partition_tmp[n_proc_y], periodic_axis_[1], 1);
            }
            exchangeSample(y_smp_sub_in_rod, y_partition_tmp, n_proc_y, my_rank_y, full_len_root_domain[1], periodic_axis_[1], comm_info_1d_[1],
                           [](const F64 &src) -> F64 { return src; });

            // comm_info_1d_[1].barrier();
            // if(comm_info_1d_[1].getRank()==0){std::cerr<<"CHECK E5"<<std::endl;}

#ifdef PS_DEBUG_DD
            for (S32 i = 1; i < y_smp_sub_in_rod.size(); i++) {
                assert(y_partition_tmp[my_rank_y] <= y_smp_sub_in_rod[i]);
                assert(y_partition_tmp[my_rank_y + 1] > y_smp_sub_in_rod[i]);
            }
#endif
            std::sort(y_smp_sub_in_rod.frontP(), y_smp_sub_in_rod.endP());
#ifdef PS_DEBUG_DD
            for (S32 i = 1; i < y_smp_sub_in_rod.size(); i++) {
                assert(y_smp_sub_in_rod[i - 1] < y_smp_sub_in_rod[i]);
            }
#endif

            // comm_info_1d_[1].barrier();
            // if(comm_info_1d_[1].getRank()==0){std::cerr<<"CHECK E6"<<std::endl;}

            // determinePartition(y_smp_sub_in_rod, y_partition_tmp, y_partition, comm_1d_[1], true);
            determinePartition(y_smp_sub_in_rod, y_partition_tmp, y_partition, comm_info_1d_[1], true);

            // comm_info_1d_[1].barrier();
            // if(comm_info_1d_[1].getRank()==0){std::cerr<<"CHECK E7"<<std::endl;}

            setEdgeBoundary(y_partition[0], y_partition[n_proc_y], periodic_axis_[1], 1);

            // comm_info_1d_[1].barrier();
            // if(comm_info_1d_[1].getRank()==0){std::cerr<<"CHECK E8"<<std::endl;}
        }

        // comm_info_.barrier();
        // if(comm_info_.getRank()==0){std::cerr<<"CHECK H"<<std::endl;}

        ////////////////////////////////////////
        // exchange samples in y directino for determination of z partition
        // broadcast y_partition to processes with the same idx (y-z plane)
        comm_info_sub_[0].broadcast(y_partition.getPointer(), y_partition.size(), 0);
        exchangeSample(pos_sample_loc_, y_partition, n_proc_y, my_rank_y, full_len_root_domain[1], periodic_axis_[1], comm_info_1d_[1],
                       [](const F64vec &src) -> F64 { return src.y; });
#ifdef PS_DEBUG_DD
        for (S32 i = 0; i < pos_sample_loc_.size(); i++) {
            assert(pos_sample_loc_[i].y < y_partition[my_rank_y + 1]);
            assert(pos_sample_loc_[i].y >= y_partition[my_rank_y]);
        }
#endif

        // comm_info_.barrier();
        // if(comm_info_.getRank()==0){std::cerr<<"CHECK G"<<std::endl;}

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        ////////////////
        /// Z DIRECTION
        S32 n_smp_sub_z = pos_sample_loc_.size();
        ReallocatableArray<F64> z_smp_sub(n_smp_sub_z, n_smp_sub_z, MemoryAllocMode::Pool);
        for (S32 i = 0; i < n_smp_sub_z; i++) z_smp_sub[i] = pos_sample_loc_[i].z;
        ReallocatableArray<F64> z_partition(n_proc_z + 1, n_proc_z + 1, MemoryAllocMode::Pool);
        // rougly calc postion of partition in z direction
        ReallocatableArray<F64> z_partition_tmp(n_proc_z + 1, n_proc_z + 1, MemoryAllocMode::Pool);

        // comm_info_.barrier();
        // if(comm_info_.getRank()==0){std::cerr<<"CHECK H"<<std::endl;}

        if (flag_use_old_domain) {
            // reuse old domain
            z_partition_tmp[0] = pos_domain_[my_rank_x * n_proc_y * n_proc_z + my_rank_y * n_proc_z].low_.z;
            for (S32 i = 1; i < n_proc_z + 1; i++) {
                const auto rank_src = my_rank_x * n_proc_y * n_proc_z + my_rank_y * n_proc_z + (i - 1);
                z_partition_tmp[i] = pos_domain_[rank_src].high_.z;
            }
        } else {
            // rougly calc postion of partition in y direction
            determineRoughPartition(z_smp_sub, z_partition_tmp, 2 * n_proc_z, comm_info_1d_[2], true);
            setEdgeBoundary(z_partition_tmp[0], z_partition_tmp[n_proc_z], periodic_axis_[2], 2);
        }
        exchangeSample(z_smp_sub, z_partition_tmp, n_proc_z, my_rank_z, full_len_root_domain[2], periodic_axis_[2], comm_info_1d_[2],
                       [](const F64 &src) -> F64 { return src; });

        // comm_info_.barrier();
        // if(comm_info_.getRank()==0){std::cerr<<"CHECK I"<<std::endl;}

#ifdef PS_DEBUG_DD
        for (S32 i = 1; i < z_smp_sub.size(); i++) {
            assert(z_partition_tmp[my_rank_z] <= z_smp_sub[i]);
            assert(z_partition_tmp[my_rank_z + 1] > z_smp_sub[i]);
        }
#endif
        std::sort(z_smp_sub.frontP(), z_smp_sub.endP());
#ifdef PS_DEBUG_DD
        for (S32 i = 1; i < z_smp_sub.size(); i++) {
            assert(z_partition_tmp[my_rank_z] <= z_smp_sub[i]);
            assert(z_partition_tmp[my_rank_z + 1] > z_smp_sub[i]);
            assert(z_smp_sub[i - 1] < z_smp_sub[i]);
        }
#endif

        // comm_info_.barrier();
        // if(comm_info_.getRank()==0){std::cerr<<"CHECK J"<<std::endl;}

        determinePartition(z_smp_sub, z_partition_tmp, z_partition, comm_info_1d_[2], true);
        setEdgeBoundary(z_partition[0], z_partition[n_proc_z], periodic_axis_[2], 2);
#endif  // PARTICLE_SIMULATOR_TWO_DIMENSION

        // comm_info_.barrier();
        // if(comm_info_.getRank()==0){std::cerr<<"CHECK K"<<std::endl;}

        ////////////////
        // braod cast coordinates of partition
        ReallocatableArray<F64> x_partition_recv(n_proc_x + 1, n_proc_x + 1, MemoryAllocMode::Pool);
        if (my_rank == 0) {
            for (S32 i = 0; i < n_proc_x + 1; i++) {
                x_partition_recv[i] = x_partition[i];
            }
        }
        // Comm::broadcast(x_partition_recv.getPointer(), n_proc_x+1, 0);
        comm_info_.broadcast(x_partition_recv.getPointer(), n_proc_x + 1, 0);

        ReallocatableArray<F64> y_partition_recv((n_proc_y + 1) * n_proc_x, (n_proc_y + 1) * n_proc_x, MemoryAllocMode::Pool);
        comm_info_1d_[0].gather(y_partition.getPointer(), n_proc_y + 1, y_partition_recv.getPointer());
        // Comm::broadcast(y_partition_recv.getPointer(), (n_proc_y+1)*n_proc_x, 0);
        comm_info_.broadcast(y_partition_recv.getPointer(), (n_proc_y + 1) * n_proc_x, 0);

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        ReallocatableArray<F64> z_partition_recv((n_proc_z + 1) * n_proc_y * n_proc_x, (n_proc_z + 1) * n_proc_y * n_proc_x, MemoryAllocMode::Pool);
        comm_info_sub_[2].gather(z_partition.getPointer(), n_proc_z + 1, z_partition_recv.getPointer());
        // Comm::broadcast(z_partition_recv.getPointer(), (n_proc_z+1)*n_proc_y*n_proc_x, 0);
        comm_info_.broadcast(z_partition_recv.getPointer(), (n_proc_z + 1) * n_proc_y * n_proc_x, 0);
#endif

        ////////////////
        // set pos_domain_
        for (S32 i = 0; i < n_proc_x; i++) {
            for (S32 j = 0; j < n_proc_y; j++) {
                for (S32 k = 0; k < n_proc_z; k++) {
                    const auto rank = i * n_proc_y * n_proc_z + j * n_proc_z + k;
                    pos_domain_temp_[rank].low_.x = x_partition_recv[i];
                    pos_domain_temp_[rank].high_.x = x_partition_recv[i + 1];
                    pos_domain_temp_[rank].low_.y = y_partition_recv[i * (n_proc_y + 1) + j];
                    pos_domain_temp_[rank].high_.y = y_partition_recv[i * (n_proc_y + 1) + j + 1];
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                    pos_domain_temp_[rank].low_.z = z_partition_recv[i * (n_proc_z + 1) * n_proc_y + j * (n_proc_z + 1) + k];
                    pos_domain_temp_[rank].high_.z = z_partition_recv[i * (n_proc_z + 1) * n_proc_y + j * (n_proc_z + 1) + k + 1];
#endif
                }
            }
        }

        // set pos_domain by using exponential moving average
        if (first_call_by_decomposeDomain) {
            first_call_by_decomposeDomain = false;
            for (S32 i = 0; i < n_proc; i++) {
                pos_domain_[i].low_ = pos_domain_temp_[i].low_;
                pos_domain_[i].high_ = pos_domain_temp_[i].high_;
            }
        } else {
            for (S32 i = 0; i < n_proc; i++) {
                pos_domain_[i].low_ = (F64)coef_ema_ * pos_domain_temp_[i].low_ + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
            }
        }

        // if(my_rank==1){
        //     for(S32 i=0; i<n_proc; i++){
        //         std::cerr<<"i= "<<i<<" pos_domain_temp_[i]= "<<pos_domain_temp_[i]<<std::endl;
        //     }
        // }

        /*
        ReallocatableArray<S32vec, 1> rank_3d(n_proc, n_proc, MemoryAllocMode::Pool);
        for(S32 i=0; i<n_proc; i++){
            rank_3d[i].x = i / (n_proc_y*n_proc_z);
            rank_3d[i].y = (i / n_proc_z) % n_proc_y ;
            rank_3d[i].z = i % n_proc_z;
        }

        for(S32 i=0; i<n_proc; i++){
            S32 rank_x = rank_3d[i].x;
            assert(pos_domain_temp_[i].low_.x  == x_partition_recv[rank_x]);
            assert(pos_domain_temp_[i].high_.x == x_partition_recv[rank_x+1]);
        }

        for(S32 i=0; i<n_proc; i++){
            S32 rank_x = rank_3d[i].x;
            S32 rank_y = rank_3d[i].y;
            assert(pos_domain_temp_[i].low_.y  == y_partition_recv[(rank_x*(n_proc_y+1))+rank_y]);
            assert(pos_domain_temp_[i].high_.y == y_partition_recv[(rank_x*(n_proc_y+1))+rank_y+1]);
        }
        */
        /*
        for(S32 i=0; i<n_proc; i++){
            S32 rank_x = rank_3d[i].x;
            S32 rank_y = rank_3d[i].y;
            S32 rank_z = rank_3d[i].z;
            assert(pos_domain_temp_[i].low_.z  == z_partition_recv[(rank_x*(n_proc_y+1))+rank_z]);
            assert(pos_domain_temp_[i].high_.z == z_partition_recv[(rank_x*(n_proc_y+1))+rank_z+1]);
        }
        */
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL

        // DEBUG_FLAG_DD = true;

#ifdef DEBUG_DD
#endif  // #ifdef DEBUG_DD

        pos_sample_tot_.freeMem();
        pos_sample_loc_.freeMem();
        auto time_end = GetWtime();
        time_profile_.decompose_domain += time_end - time_offset;
        time_profile_.decompose_domain_all += time_end - time_offset;

        // comm_info_.barrier();
        // if(comm_info_.getRank()==0){std::cerr<<"CHECK Z"<<std::endl;}
    }

    template <class Tpsys>
    void decomposeDomainAll2(Tpsys &psys, const F64 wgh) {
        const bool clear = true;
        collectSampleParticle(psys, clear, wgh);
        decomposeDomainMultiStep2();
    }

    template <class Tpsys>
    void decomposeDomainAll2(Tpsys &psys) {
        const F64 wgh = psys.getNumberOfParticleLocal();
        decomposeDomainAll2(psys, wgh);
    }

#if 0
// VAMODE 
// under construction    
    F64 getWeight() { return 1.0; }
    template <typename Head, typename... Tail>
    F64 getWeight(Head &head, Tail &&...tail) {
        if constexpr (IsParticleSystem<Head>::value) {
            return getWeight();
        } else {
            return head;
        }
    }

    template <typename Head, typename... Tail>
    void collectSampleParticleImplVA(bool &clear, Head &head, Tail &&...tail) {
        if constexpr (IsParticleSystem<Head>::value) {
            auto wgh = getWeight(tail...);
            collectSampleParticle(head, clear, wgh);
            clear = false;
        }
        if constexpr (sizeof...(tail) > 0) {
            collectSampleParticleImplVA(clear, tail...);
        }
    }

    template <typename... Args>
    void collectSampleParticleVA(Args &&...args) {
        bool clear = true;
        collectSampleParticleImplVA(clear, args...);
    }
    /*
    template <typename... Args>
    void decomposeDomainAllVA(Args &&...args) {
        auto t0 = GetWtime();
        collectSampleParticleVA(args...);
        // time_profile_.collect_sample_particle += GetWtime() - t0;
        decomposeDomain();
        time_profile_.decompose_domain_all += GetWtime() - t0;
    }
    */

    template <typename... Args>
    void decomposeDomainAll(Args &&...args) {
        auto t0 = GetWtime();
        collectSampleParticleVA(args...);
        // time_profile_.collect_sample_particle += GetWtime() - t0;
        decomposeDomain();
        time_profile_.decompose_domain_all += GetWtime() - t0;
    }
#endif
};

template <typename T>
struct IsDomainInfo : std::false_type {};
template <>
struct IsDomainInfo<DomainInfo> : std::true_type {};
template <>
struct IsDomainInfo<const DomainInfo> : std::true_type {};

}  // namespace ParticleSimulator
