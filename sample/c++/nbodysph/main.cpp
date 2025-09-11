// #define SANITY_CHECK_REALLOCATABLE_ARRAY
#include <sys/stat.h>
#include "header.h"

template <typename F, typename... Objs>
void InvokeOnEach(F&& f, Objs&&... objs) {
    (std::invoke(std::forward<F>(f), std::forward<Objs>(objs)), ...);
}

void makeOutputDirectory(const char* dir_name) {
    struct stat st;
    PS::S32 ret;
    if (PS::Comm::getRank() == 0) {
        if (stat(dir_name, &st) != 0) {
            ret = mkdir(dir_name, 0777);
        } else {
            ret = 0;  // the directory named dir_name already exists.
        }
    }
    PS::Comm::broadcast(&ret, 1);
    if (ret == 0) {
        if (PS::Comm::getRank() == 0) fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
        if (PS::Comm::getRank() == 0) fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
        PS::Abort();
    }
}

int main(int argc, char* argv[]) {
    //////////////////
    // Create vars.
    //////////////////
    PS::Initialize(argc, argv);
    makeOutputDirectory("result");
    PS::ParticleSystem<RealPtcl> sph_system;
    sph_system.initialize();
    PS::DomainInfo dinfo;
    dinfo.initialize();

    PS::F64 dt, end_time;
    //////////////////
    // Disp. Info
    //////////////////
    DisplayInfo();
    //////////////////
    // Setup Initial
    //////////////////
    SetupIC(sph_system, &end_time, dinfo);
    Initialize(sph_system);
    // Dom. info
    dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
    dinfo.decomposeDomainAll(sph_system);
    sph_system.exchangeParticle(dinfo);
    PS::S32 n_loc = sph_system.getNumberOfParticleLocal();
    // plant tree
    PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather dens_tree;
    PS::TreeForForceShort<RESULT::Drvt, EPI::Drvt, EPJ::Drvt>::Gather drvt_tree;
    PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
    PS::TreeForForceLong<RESULT::Grav, EPI::Grav, EPJ::Grav>::Monopole grav_tree;

    dens_tree.initialize(n_loc);
    drvt_tree.initialize(n_loc);
    hydr_tree.initialize(n_loc);
    grav_tree.initialize(n_loc);

    // dens_tree.setExchangeLETMode(PS::EXCHANGE_LET_P2P_FAST);
    // drvt_tree.setExchangeLETMode(PS::EXCHANGE_LET_P2P_FAST);
    // hydr_tree.setExchangeLETMode(PS::EXCHANGE_LET_P2P_FAST);
    // grav_tree.setExchangeLETMode(PS::EXCHANGE_LET_P2P_FAST);

    for (int loop = 0; loop <= 5; ++loop) {
        dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
    }
    for (PS::S32 i = 0; i < n_loc; ++i) {
        sph_system[i].setPressure();
    }
    drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
    hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
    grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

    dt = getTimeStepGlobal(sph_system);

    PS::F64 time = 0.0;
    std::cout << std::scientific << std::setprecision(16) << "time = " << time << ", dt = " << dt << std::endl;

    PS::S32 step = 0;
    for (; time < end_time; time += dt, ++step) {
        InitialKick(sph_system, dt);
        FullDrift(sph_system, dt);
        sph_system.adjustPositionIntoRootDomain(dinfo);
        Predict(sph_system, dt);
        dinfo.decomposeDomainAll(sph_system);
        sph_system.exchangeParticle(dinfo);
        n_loc = sph_system.getNumberOfParticleLocal();
        for (int loop = 0; loop <= 2; ++loop) {
            dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
        }
        for (PS::S32 i = 0; i < n_loc; ++i) {
            sph_system[i].setPressure();
        }
        drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
        hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
        grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
        PS::F64 mass0, eng0;  // initial total mass and energy
        PS::F64vec mom0;      // initial total momentum
        CalculateConservativeVariables(mass0, eng0, mom0, sph_system);
        dt = getTimeStepGlobal(sph_system);

        FinalKick(sph_system, dt);
        if (step % PARAM::OUTPUT_INTERVAL == 0) {
            FileHeader header;
            header.time = time;
            header.Nbody = sph_system.getNumberOfParticleGlobal();
            char filename[256];
            sprintf(filename, "result/%04d.dat", step);
            sph_system.writeParticleAscii(filename, header);
            if (PS::Comm::getRank() == 0) {
                std::cout << "//================================" << std::endl;
                std::cout << "output " << filename << "." << std::endl;
                std::cout << "//================================" << std::endl;
            }
            /*
                        if (PS::Comm::getRank() == 0) {
                            std::cout << "dinfo" << std::endl;
                            dinfo.dumpTimeProfile();
                            std::cout << "sph_system" << std::endl;
                            sph_system.dumpTimeProfile();
                            std::cout << "drvt_tree" << std::endl;
                            drvt_tree.dumpTimeProfile();
                            std::cout << "hydr_tree" << std::endl;
                            hydr_tree.dumpTimeProfile();
                            std::cout << "grav_tree" << std::endl;
                            grav_tree.dumpTimeProfile();
                            InvokeOnEach([&](auto& obj) { obj.clearTimeProfile(); }, dinfo, sph_system, drvt_tree, hydr_tree, grav_tree);
                        }
            */
        }
        PS::F64 mass, eng;
        PS::F64vec mom;
        CalculateConservativeVariables(mass, eng, mom, sph_system);
        if (PS::Comm::getRank() == 0) {
            std::cout << "//================================" << std::endl;
            std::cout << std::scientific << std::setprecision(16) << "time = " << time << ", dt = " << dt << std::endl;
            std::cout << "step = " << step << std::endl;
            printf("mass-mass0= %.10e, mom-mom0= %.10e %.10e %.10e, (eng-eng0)/eng0= %.10e\n", mass - mass0, mom.x - mom0.x, mom.y - mom0.y,
                   mom.z - mom0.z, (eng - eng0) / eng0);
            std::cout << "//================================" << std::endl;
        }
    }

    PS::Finalize();
    return 0;
}
