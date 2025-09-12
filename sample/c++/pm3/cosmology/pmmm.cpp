#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <float.h>
#include <cstdio>
#include <cstdlib>
#include <particle_simulator.hpp>
#include <particle_mesh_multipole.hpp>
#include "pmmm.hpp"
#include "run_param.hpp"
#include "cosmology.hpp"
#include "mpi.h"

// static PS::F64 EPS_FOR_PP = 0.1 / 32.0;  // for N=32^3
static PS::F64 EPS_FOR_PP = 0.1 / 128.0;  // for N=128^3
// static constexpr PS::S64 N_ADD_LOCAL = 100;

// for debugging
template <typename Tptcl>
void AddCloseParticles(Tptcl &ptcl) {
    PS::S64 n_loc = ptcl.getNumberOfParticleLocal();
    // const auto n_add = std::min(n_loc, N_ADD_LOCAL);
    // for (int i = 0; i < n_add; i++) {
    for (int i = 0; i < n_loc; i++) {
        auto p = ptcl[i];
        p.mass = 0.0;
        p.pos.x += 1e-8;
        p.pos.y += 1e-8;
        p.pos.z += 1e-8;
        ptcl.addOneParticle(p);
    }
}

template <class Tp>
void setup_random_particle(Tp &ptcl, run_param &this_run, const PS::S32 npart_total) {
    PS::S32 rank = PS::Comm::getRank();
    PS::S64 npart_local = (rank == 0) ? npart_total : 0;

    ptcl.setNumberOfParticleLocal(npart_local);

    this_run.npart_total = npart_total;
    this_run.mpi_nproc = PS::Comm::getNumberOfProc();
    this_run.mpi_rank = rank;

    PS::MT::init_genrand(0);

    for (PS::S32 i = 0; i < npart_local; i++) {
        ptcl[i].mass = 3.0 * this_run.cosm.omegam / (8.0 * AC::PI * (PS::F32)npart_total);
        ptcl[i].pos[0] = PS::MT::genrand_res53();
        ptcl[i].pos[1] = PS::MT::genrand_res53();
        ptcl[i].pos[2] = PS::MT::genrand_res53();
        ptcl[i].vel[0] = 0.0;
        ptcl[i].vel[1] = 0.0;
        ptcl[i].vel[2] = 0.0;
        // ptcl[i].eps = 0.1/pow(npart_total, 1.0/3.0);
        ptcl[i].eps = EPS_FOR_PP;
    }
}

template <class Tp>
void read_SB_particle(Tp &ptcl, run_param &this_run, const char *input_file) {
    ptcl.readParticleBinary(input_file);
    // AddCloseParticles(ptcl);
    this_run.npart_total = ptcl.getNumberOfParticleGlobal();
    this_run.npart_local = ptcl.getNumberOfParticleLocal();
    this_run.mpi_nproc = PS::Comm::getNumberOfProc();
    this_run.mpi_rank = PS::Comm::getRank();
    for (PS::S32 i = 0; i < this_run.npart_local; i++) {
        // ptcl[i].eps = 0.1/pow(this_run.npart_total, 1.0/3.0);
        ptcl[i].eps = EPS_FOR_PP;
    }
}

template <class Tptcl>
void read_param_file(Tptcl &ptcl, run_param &this_run, const char *input_file) {
    if (PS::Comm::getRank() == 0) std::cerr << "input_file=" << input_file << std::endl;
    FILE *param_fp;
    param_fp = fopen(input_file, "r");
    if (param_fp == NULL) {
        fprintf(stderr, "File %s not found in input_params.\n", input_file);
        PS::Abort();
        exit(1);
    }
    int ret = fscanf(param_fp, "%d", &this_run.mode);
    if (ret == EOF) {
        fprintf(stderr, "Input error of mode in input_params()\n");
        PS::Abort();
        exit(1);
    }
    ret = fscanf(param_fp, "%lf", &EPS_FOR_PP);
    ret = fscanf(param_fp, "%lf", &this_run.theta);
    ret = fscanf(param_fp, "%f", &this_run.zend);
    if (PS::Comm::getRank() == 0) {
        std::cerr << "this_run.theta=" << this_run.theta << std::endl;
        std::cerr << "this_run.zend=" << this_run.zend << std::endl;
    }
    if (this_run.mode == run_param::SANTABARBARA) {
        this_run.znow = 63.0;
        this_run.cosm.omegam = 1.0;
        this_run.cosm.omegav = this_run.cosm.omegab = this_run.cosm.omeganu = 0.0;
        // FPtreepm::H0 = 50.0;
        // FPtreepm::Lbnd = 64.0;
        FP::H0 = 50.0;
        FP::Lbnd = 64.0;
        char snap_name[1024];
        ret = fscanf(param_fp, "%s", snap_name);
        if (PS::Comm::getRank() == 0) std::cerr << "snap_name:" << snap_name << std::endl;
        read_SB_particle(ptcl, this_run, snap_name);
    } else if (this_run.mode == run_param::READ_FILE) {
        char snap_name[1024];
        ret = fscanf(param_fp, "%s", snap_name);
        if (PS::Comm::getRank() == 0) std::cerr << "snap_name:" << snap_name << std::endl;
        ret = fscanf(param_fp, "%f", &this_run.znow);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegam);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegav);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegab);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omeganu);
        // ret = fscanf(param_fp, "%lf", &FPtreepm::H0);
        ret = fscanf(param_fp, "%lf", &FP::H0);
        // ret = fscanf(param_fp, "%lf", &FPtreepm::Lbnd);
        ret = fscanf(param_fp, "%lf", &FP::Lbnd);
        read_SB_particle(ptcl, this_run, snap_name);
    } else if (this_run.mode == run_param::RANDOM) {
        ret = fscanf(param_fp, "%f", &this_run.znow);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegam);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegav);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegab);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omeganu);
        // ret = fscanf(param_fp, "%lld", &this_run.npart_total);
        ret = fscanf(param_fp, "%ld", &this_run.npart_total);
        setup_random_particle(ptcl, this_run, this_run.npart_total);
    } else {
        PS::Abort();
        exit(1);
    }
    this_run.anow = 1.0 / (1.0 + this_run.znow);
    this_run.tnow = this_run.cosm.atotime(this_run.anow);
    this_run.update_expansion(this_run.tnow);
    this_run.input_params(param_fp);
}

int main(int argc, char **argv) {
    PS::Initialize(argc, argv);
    // PS::PM::ParticleMesh pm;
    // PS::ParticleSystem<FPtreepm> ptcl;
    PS::ParticleSystem<FP> ptcl_pm3;
    ptcl_pm3.initialize();
    PS::DomainInfo domain_info;
    domain_info.initialize();
    run_param this_run;

    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "read file" << std::endl;

    read_param_file(ptcl_pm3, this_run, argv[1]);

    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "this_run.npart_local=" << this_run.npart_local << std::endl;
    PS::Comm::barrier();

    domain_info.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    domain_info.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), PS::F64vec(1.0, 1.0, 1.0));
    domain_info.decomposeDomainAll(ptcl_pm3);
    ptcl_pm3.adjustPositionIntoRootDomain(domain_info);
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "start particle exchange" << std::endl;
    PS::Comm::barrier();
    ptcl_pm3.exchangeParticle(domain_info);
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "finish particle exchange" << std::endl;
    PS::Comm::barrier();
    this_run.npart_local = ptcl_pm3.getNumberOfParticleLocal();
    // const auto n_ptcl_loc = ptcl_pm3.getNumberOfParticleLocal();

    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "tree constructor" << std::endl;
    FP::eps = EPS_FOR_PP;
    //PS::TreeForForce<PS::SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Result, EP, EP, PS::MomentQuadrupoleGeometricCenter,
//                     PS::MomentQuadrupoleGeometricCenter, PS::SPJQuadrupoleGeometricCenter>
        //treepm3;
    PS::TreeForForce<PS::SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Result, EPI, EPJ, PS::MomentMonopole, PS::MomentMonopole, PS::SPJMonopole> treepm3;
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "tree initialize" << std::endl;

    // const PS::S32 n_leaf_limit = 1;
    const PS::U32 n_leaf_limit = 8;
    const PS::U32 n_grp_limit = 256;
    treepm3.initialize(ptcl_pm3.getNumberOfParticleLocal(), this_run.theta, n_leaf_limit, n_grp_limit);
    PS::S32 n_cell = 8;
    PS::S32 i_cut = 2;
    treepm3.setParamPMMM(n_cell, i_cut);

    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "start tree force" << std::endl;
    treepm3.calcForceAllAndWriteBack(CalcGravity<EPJ>(), CalcGravity<PS::SPJMonopole>(), ptcl_pm3, domain_info);
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "finish tree force" << std::endl;


    // calc PMM force
    auto p_exp = 2;
    PS::PMMM::ParticleMeshMultipole<Result, EPI> pm3;
    pm3.initialize(p_exp);
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "start pmmm force" << std::endl;
    pm3.calcForceAllAndWriteBack(treepm3, domain_info, ptcl_pm3);

    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "finish pmmm force" << std::endl;
    // PMMM
    //////////////////////

    PS::F64 dtime = calc_dtime(ptcl_pm3, this_run);
    PS::F64 dtime_prev, dtime_mid;
    drift_ptcl(ptcl_pm3, domain_info, 0.5 * dtime);
    domain_info.decomposeDomainAll(ptcl_pm3);
    ptcl_pm3.adjustPositionIntoRootDomain(domain_info);
    ptcl_pm3.exchangeParticle(domain_info);
    this_run.npart_local = ptcl_pm3.getNumberOfParticleLocal();
    Map2D map_2d;
    PS::F64 res = 2.0;
    map_2d.initialize(ptcl_pm3, this_run, res);
    PS::F64 wtime_tot, wtime_tree, wtime_pm3, wtime_misc;
    while (this_run.znow > this_run.zend) {
        auto n_int_ep_ep_glb = treepm3.getNumberOfInteractionEPEPGlobal();
        auto n_int_ep_sp_glb = treepm3.getNumberOfInteractionEPSPGlobal();
        auto n_ptcl_glb = ptcl_pm3.getNumberOfParticleGlobal();
        auto n_ipg = treepm3.getNumberOfIPGGlobal();
        if (PS::Comm::getRank() == 0) {
            std::cout << "====================================================" << std::endl;
            std::cout << "this_run.step=" << this_run.step << " this_run.znow=" << this_run.znow << " this_run.zend=" << this_run.zend
                      << "n_ptcl_glb= " << n_ptcl_glb << std::endl;
            std::cout << "n_int_ep_ep_glb=" << n_int_ep_ep_glb << " n_int_ep_sp_glb=" << n_int_ep_sp_glb << std::endl;
            std::cout << "<n_epi>= " << n_ptcl_glb / n_ipg << " <n_epj>= " << n_int_ep_ep_glb / n_ptcl_glb
                      << " <n_spj>= " << n_int_ep_sp_glb / n_ptcl_glb << std::endl;
            std::cout << "wtime_tot=" << wtime_tot << " wtime_tree=" << wtime_tree << " wtime_pm3=" << wtime_pm3 << " wtime_misc=" << wtime_misc
                      << std::endl;
            treepm3.dumpTimeProfile();
        }
        treepm3.clearCounterAll();
        PS::Comm::barrier();
        auto t0 = PS::GetWtime();
        //treepm3.calcForceAllAndWriteBack(CalcGravity<EP>(), CalcGravity<PS::SPJQuadrupoleGeometricCenter>(), ptcl_pm3, domain_info);
        treepm3.calcForceAllAndWriteBack(CalcGravity<EPJ>(), CalcGravity<PS::SPJMonopole>(), ptcl_pm3, domain_info);
        auto t1 = PS::GetWtime();
        // pm3.calcForceAllAndWriteBack(treepm3, domain_info, ptcl_pm3, flag_clear_force_pmm, n_thread_max);
        pm3.calcForceAllAndWriteBack(treepm3, domain_info, ptcl_pm3);
        auto t2 = PS::GetWtime();
        this_run.tnow += 0.5 * dtime;
        this_run.update_expansion(this_run.tnow);
        kick_ptcl(ptcl_pm3, dtime, this_run);
        this_run.tnow += 0.5 * dtime;
        this_run.update_expansion(this_run.tnow);
        dtime_prev = dtime;
        dtime = calc_dtime(ptcl_pm3, this_run);
        dtime_mid = 0.5 * (dtime_prev + dtime);
        drift_ptcl(ptcl_pm3, domain_info, dtime_mid);
        domain_info.decomposeDomainAll(ptcl_pm3);
        ptcl_pm3.adjustPositionIntoRootDomain(domain_info);
        ptcl_pm3.exchangeParticle(domain_info);
        auto t3 = PS::GetWtime();
        this_run.npart_local = ptcl_pm3.getNumberOfParticleLocal();
        output_data_in_run(ptcl_pm3, this_run, map_2d);
        this_run.step++;
        wtime_tot = t3 - t0;
        wtime_tree = t1 - t0;
        wtime_pm3 = t2 - t1;
        wtime_misc = t3 - t2;
    }

    PS::Finalize();
    return (0);
}
