#include <float.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <particle_simulator.hpp>
#include <particle_mesh_multipole.hpp>
#include "mpi.h"
#include "ic.hpp"
#if defined(USE_PM)
#include <particle_mesh.hpp>
#include <param_fdps.h>
#endif
#include "user_defined.hpp"
#include "options.hpp"

template <typename Tpsys>
void DumpParticles(Tpsys &psys) {
    for (int i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        std::cout << "id= " << psys[i].id << " charge= " << psys[i].charge << " pos= " << psys[i].pos << " pot= " << psys[i].pot
                  << " grad= " << psys[i].grad << std::endl;
    }
}

template <typename Tpsys>
void ClearForces(Tpsys &psys) {
    for (int i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].pot = 0.0;
        psys[i].grad = 0.0;
    }
}

template <typename Tforce, typename Tpsys>
void CopyForces(Tforce &force, Tpsys &psys) {
    for (int i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        force[i].pot = psys[i].pot;
        force[i].grad = psys[i].grad;
    }
}

int main(int argc, char **argv) {
    PS::Initialize(argc, argv, 100000000);
    PS::F64ort POS_UNIT_CELL(PS::F64vec(0.0), PS::F64vec(1.0));  // computational domain
    constexpr PS::S32 n_img_ewald_r = 3;                         // number of images in real space for Ewald summation
    constexpr PS::S32 n_img_ewald_k = 5;                         // number of images in reciprocal space for Ewald summation
    const auto my_rank = PS::Comm::getRank();
    const auto n_proc = PS::Comm::getNumberOfProc();
    Params params;
    select_options(params, argc, argv);

    PS::ParticleSystem<FP> ptcl;
    ptcl.setNumberOfParticleLocal(0);
    ptcl.initialize();
    if (params.read_file) {
        if (my_rank == 0) {
            std::ifstream ifs(params.input_file);
            if (!ifs) {
                std::cerr << "Error: Cannot open input file: " << params.input_file << std::endl;
                exit(1);
            }
            PS::S64 n;
            ifs >> n;
            ptcl.setNumberOfParticleLocal(n);
            PS::F64vec pos;
            PS::F64 charge;
            for (int i = 0; i < n; i++) {
                ifs >> charge >> pos.x >> pos.y >> pos.z;
                FP fp;
                fp.setPos(pos);
                fp.charge = charge;

                ptcl[i].pos = pos;
                ptcl[i].charge = charge;
                ptcl[i].id = i;
            }
            ifs.close();
            POS_UNIT_CELL = PS::F64ort(PS::F64vec(0.0), PS::F64vec(1.0));
        }
    } else {
        if (my_rank == 0) {
            std::cerr << "Generating particles..." << std::endl;
        }
        setup_particle(ptcl, params, POS_UNIT_CELL);
    }
    /*
        // for debugging
        for (int i = 0; i < ptcl.getNumberOfParticleLocal(); i++) {
            if (ptcl[i].id == 0) {
                // charge = -1.0;
                ptcl[i].charge = -100.0;
                ptcl[i].pos = PS::F64vec(0.5);
            }
        }
    */
    PS::DomainInfo domain_info;
    domain_info.initialize();
    domain_info.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    domain_info.setPosRootDomain(POS_UNIT_CELL);
    domain_info.decomposeDomainAll(ptcl);
    ptcl.adjustPositionIntoRootDomain(domain_info);
    ptcl.exchangeParticle(domain_info);
    PS::S32 n_ptcl_loc = ptcl.getNumberOfParticleLocal();
    PS::F64 t0, t1, t2, t3;  // for time measurement
#if 1
    // For Ewald
    if (my_rank == 0) std::cerr << "Ewald" << std::endl;
    EwaldEngine::setQTot(ptcl);
    EwaldEngine::setAlpha(2.4);
    EwaldEngine::setUnitCell(POS_UNIT_CELL);
    EwaldEngine::setNImgR(n_img_ewald_r);
    EwaldEngine::setNImgK(n_img_ewald_k);
    EP::r_search = n_img_ewald_r * POS_UNIT_CELL.getFullLength().x;  // Set the search radius for Ewald particles
    PS::TreeForForce<PS::SEARCH_MODE_LONG_CUTOFF, Force, EP, EP, PS::MomentMonopoleGeometricCenter, PS::MomentMonopoleGeometricCenter,
                     PS::SPJMonopoleGeometricCenter>
        tree_ewald;
    tree_ewald.initialize(n_ptcl_loc, 0.0);  // theta is 0.0
    for (int loop = 0; loop < 2; loop++) {
        ClearForces(ptcl);
        t0 = PS::GetWtime();
        tree_ewald.calcForceAllAndWriteBack(EwaldEngine::calcForceSR<EP, EP, Force>,
                                            EwaldEngine::calcForceSR<EP, PS::SPJMonopoleGeometricCenter, Force>, ptcl, domain_info);
        t1 = PS::GetWtime();
        EwaldEngine::energyCorrection(ptcl);
        t2 = PS::GetWtime();
        EwaldEngine::calcForceLR(ptcl);
        t3 = PS::GetWtime();
    }
    if (my_rank == 0) {
        std::cerr << "***************" << std::endl;
        std::cerr << "Ewald: " << t3 - t0 << std::endl;
        std::cerr << "Ewald Short Force: " << t1 - t0 << std::endl;
        std::cerr << "Ewald Energy Correction: " << t2 - t1 << std::endl;
        std::cerr << "Ewald Long Force: " << t3 - t2 << std::endl;
    }
    //DumpParticles(ptcl);
    std::vector<Force> force_ewald(n_ptcl_loc);
    CopyForces(force_ewald, ptcl);
    ClearForces(ptcl);
#endif

    // PM3
    if (my_rank == 0) std::cerr << "PM3" << std::endl;
    constexpr PS::S32 P_MAX = 7;
    std::vector<Force> force_pm3[P_MAX];
    PS::TreeForForce<PS::SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Force, EP, EP, PS::MomentQuadrupoleGeometricCenter,
                     PS::MomentQuadrupoleGeometricCenter, PS::SPJQuadrupoleGeometricCenter>
        tree_pm3;
    tree_pm3.initialize(n_ptcl_loc, params.theta);
    tree_pm3.setParamPMMM(params.n_cell_1d, params.i_cut);
    for (PS::S32 p = 0; p < P_MAX; p++) {
        PS::PMMM::ParticleMeshMultipole<Force, EP> pm3;
        pm3.initialize(p);
        for (int loop = 0; loop < 2; loop++) {
            ClearForces(ptcl);
            t0 = PS::GetWtime();
            tree_pm3.calcForceAllAndWriteBack(calc_p0<EP, EP, Force>, calc_p2<EP, PS::SPJQuadrupoleGeometricCenter, Force>, ptcl, domain_info);
            t1 = PS::GetWtime();
            pm3.calcForceAllAndWriteBack(tree_pm3, domain_info, ptcl);
            t2 = PS::GetWtime();
        }
        // DumpParticles(ptcl);
        if (my_rank == 0) {
            std::cerr << "***************" << std::endl;
            std::cerr << "PM3 (p= " << p << "): " << t2 - t0 << std::endl;
            std::cerr << "PM3 Tree Force: " << t1 - t0 << std::endl;
            std::cerr << "PM3 Mesh Force: " << t2 - t1 << std::endl;
        }
        force_pm3[p].resize(ptcl.getNumberOfParticleLocal());
        CopyForces(force_pm3[p], ptcl);
        ClearForces(ptcl);
    }

    // 粒子の個数を取得
    PS::S32 n_ptcl_array[n_proc];
    PS::S32 n_disp_ptcl_array[n_proc + 1];
    PS::Comm::allGather(&n_ptcl_loc, 1, &n_ptcl_array[0]);
    n_disp_ptcl_array[0] = 0;
    for (int i = 0; i < n_proc; i++) {
        n_disp_ptcl_array[i + 1] = n_disp_ptcl_array[i] + n_ptcl_array[i];
    }

    std::vector<PS::F64> dgrad_pm3_loc[P_MAX];
    std::vector<PS::F64> dpot_pm3_loc[P_MAX];
    std::vector<PS::F64> dgrad_pm3_glb[P_MAX];
    std::vector<PS::F64> dpot_pm3_glb[P_MAX];
    for (int p = 0; p < P_MAX; p++) {
        dgrad_pm3_loc[p].resize(n_ptcl_loc);
        dpot_pm3_loc[p].resize(n_ptcl_loc);
        dgrad_pm3_glb[p].resize(params.n_ptcl_glb);
        dpot_pm3_glb[p].resize(params.n_ptcl_glb);
        for (int i = 0; i < n_ptcl_loc; i++) {
            auto dgrad = sqrt(((force_pm3[p][i].grad - force_ewald[i].grad) * (force_pm3[p][i].grad - force_ewald[i].grad)) /
                              (force_ewald[i].grad * force_ewald[i].grad));
            auto dpot = fabs((force_pm3[p][i].pot - force_ewald[i].pot) / force_ewald[i].pot);
            dgrad_pm3_loc[p][i] = dgrad;
            dpot_pm3_loc[p][i] = dpot;
        }
        PS::Comm::gatherV(&dgrad_pm3_loc[p][0], n_ptcl_loc, &dgrad_pm3_glb[p][0], &n_ptcl_array[0], &n_disp_ptcl_array[0]);
        PS::Comm::gatherV(&dpot_pm3_loc[p][0], n_ptcl_loc, &dpot_pm3_glb[p][0], &n_ptcl_array[0], &n_disp_ptcl_array[0]);
        std::sort(&dpot_pm3_glb[p][0], &dpot_pm3_glb[p][params.n_ptcl_glb]);
        std::sort(&dgrad_pm3_glb[p][0], &dgrad_pm3_glb[p][params.n_ptcl_glb]);
    }

#if defined(USE_PM)
    if (my_rank == 0) std::cerr << "PM" << std::endl;
    std::vector<Force> force_pm(ptcl.getNumberOfParticleLocal());
    if (n_proc < 2) {
        std::cout << "PM module requires at least 2 processes" << std::endl;
    } else {
        EP::r_search = 3.0 / SIZE_OF_MESH;  // Set the search radius for Particle Mesh
        PS::TreeForForceLong<Force, EP, EP>::MonopoleWithCutoff tree_pm;
        tree_pm.initialize(n_ptcl_loc, params.theta);
        PS::PM::ParticleMesh pm;
        for (int loop = 0; loop < 2; loop++) {
            ClearForces(ptcl);
            t0 = PS::GetWtime();
            tree_pm.calcForceAllAndWriteBack(calc_pp_force_for_PM<EP, EP, Force>(), calc_pp_force_for_PM<EP, PS::SPJMonopole, Force>(), ptcl,
                                             domain_info);
            t1 = PS::GetWtime();
            pm.calcForceAllAndWriteBack(ptcl, domain_info);
            t2 = PS::GetWtime();
        }
        if (my_rank == 0) {
            std::cerr << "***************" << std::endl;
            std::cerr << "PM: " << t2 - t0 << std::endl;
            std::cerr << "PM Tree Force: " << t1 - t0 << std::endl;
            std::cerr << "PM Mesh Force: " << t2 - t1 << std::endl;
        }
        // std::cout << "PM"<<std::endl;
        // DumpParticles(ptcl);
        force_pm.resize(ptcl.getNumberOfParticleLocal());
        CopyForces(force_pm, ptcl);
        ClearForces(ptcl);
    }
    std::vector<PS::F64> dgrad_pm_loc(n_ptcl_loc);
    std::vector<PS::F64> dpot_pm_loc(n_ptcl_loc);
    std::vector<PS::F64> dgrad_pm_glb(params.n_ptcl_glb);
    std::vector<PS::F64> dpot_pm_glb(params.n_ptcl_glb);
    for (int i = 0; i < n_ptcl_loc; i++) {
        auto dgrad =
            sqrt(((force_pm[i].grad - force_ewald[i].grad) * (force_pm[i].grad - force_ewald[i].grad)) / (force_ewald[i].grad * force_ewald[i].grad));
        auto dpot = fabs((force_pm[i].pot - force_ewald[i].pot) / force_ewald[i].pot);
        dgrad_pm_loc[i] = dgrad;
        dpot_pm_loc[i] = dpot;
    }
    PS::Comm::gatherV(&dgrad_pm_loc[0], n_ptcl_loc, &dgrad_pm_glb[0], &n_ptcl_array[0], &n_disp_ptcl_array[0]);
    PS::Comm::gatherV(&dpot_pm_loc[0], n_ptcl_loc, &dpot_pm_glb[0], &n_ptcl_array[0], &n_disp_ptcl_array[0]);
    std::sort(&dpot_pm_glb[0], &dpot_pm_glb[params.n_ptcl_glb]);
    std::sort(&dgrad_pm_glb[0], &dgrad_pm_glb[params.n_ptcl_glb]);
#endif

    // output results
    if (my_rank == 0) {
        std::ofstream fout;
        fout.open("result.dat");
        if (!fout) {
            std::cerr << "cannot open the file: result.dat" << std::endl;
            return 1;
        }
        for (PS::S32 i = 0; i < params.n_ptcl_glb; i++) {
            auto frac = PS::F64(i + 1) / PS::F64(params.n_ptcl_glb);
            fout << frac;
            for (PS::S32 p = 0; p < P_MAX; p++) {
                fout << "   " << dgrad_pm3_glb[p][i] << "   " << dpot_pm3_glb[p][i];
            }
#if defined(USE_PM)
            fout << "   " << dgrad_pm_glb[i];
#endif
            fout << std::endl;
        }
    }
    PS::Finalize();
    return 0;
}
