// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>
// Include the header file of Phantom-GRAPE library
#if defined(ENABLE_PHANTOM_GRAPE_X86)
#include <gp5util.h>
#endif
// Include user-defined headers
#include "macro_defs.hpp"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "user_defined.hpp"
#include "ic.hpp"
#include "leapfrog.hpp"

void makeOutputDirectory(char * dir_name) {
    struct stat st;
    PS::S32 ret;
    if (PS::Comm::getRank() == 0) {
        if (stat(dir_name, &st) != 0) {
            ret = mkdir(dir_name, 0777);
        } else {
            ret = 0; // the directory named dir_name already exists.
        }
    } 
    PS::Comm::broadcast(&ret, 1);
    if (ret == 0) {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
        PS::Abort();
    }
}

void calcDensity(PS::ParticleSystem<FP_sph> & psys,
                 PS::DomainInfo & dinfo, 
                 PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree) {
#if defined(ENABLE_VARIABLE_SMOOTHING_LENGTH)
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    const PS::S64 n_glb = psys.getNumberOfParticleGlobal();
    // Determine the density and the smoothing length so that Eq.(6) in Springel (2005) 
    // holds within a specified accuracy.
    SCF_smth = 1.25;
    PS::S32 iter = 0;
    for (;;) {
        iter++;
        // Compute density, etc.
        tree.calcForceAllAndWriteBack(CalcDensity(), psys, dinfo);
        // Check convergence
        PS::S32 n_compl_loc = 0;
        for (PS::S32 i = 0; i < n_loc; i++) {
            if (psys[i].flag == 1) n_compl_loc++;
        }
        const PS::S64 n_compl = PS::Comm::getSum(n_compl_loc);
        if (n_compl == n_glb) break;
    }
    // Reset SCF_smth
    SCF_smth = 1.0;
#else
    SCF_smth = 1.0;
    tree.calcForceAllAndWriteBack(CalcDensity(), psys, dinfo);
#endif
}

void setEntropy(PS::ParticleSystem<FP_sph>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].setEntropy();
    }
}

void setPressure(PS::ParticleSystem<FP_sph>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].setPressure();
    }
}

void setCircularVelocity(PS::ParticleSystem<FP_sph>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        const PS::F64vec acc = psys[i].acc_grav + psys[i].acc_hydro;
        const PS::F64 r = std::sqrt(psys[i].pos * psys[i].pos);
        const PS::F64 R = std::sqrt(psys[i].pos.x * psys[i].pos.x 
                                   +psys[i].pos.y * psys[i].pos.y);
        const PS::F64 phi = std::atan2(psys[i].pos.y, psys[i].pos.x);
        const PS::F64 theta = std::atan2(R, psys[i].pos.z);
        PS::F64vec base_vect_r;
        base_vect_r.x = std::sin(theta) * std::cos(phi);
        base_vect_r.y = std::sin(theta) * std::sin(phi);
        base_vect_r.z = std::cos(theta);
        const PS::F64 vel_circ = std::sqrt(r * std::abs(acc * base_vect_r));
        psys[i].vel.x = - vel_circ * std::sin(phi);
        psys[i].vel.y =   vel_circ * std::cos(phi);
        psys[i].vel.z = 0.0;
    }
#if 0
   // [for debug] 
   std::string filename;
   std::ostringstream ss;
   std::ofstream output_file;
   // Output the velocity field
   ss << "velc_fld" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
   filename = ss.str();
   output_file.open(filename.c_str(),std::ios::trunc);
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   output_file << std::setprecision(15) << std::showpos;
   for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
       output_file << psys[i].pos.x << " " 
                   << psys[i].pos.y << " "
                   << psys[i].vel.x << " " 
                   << psys[i].vel.y << " "
                   << std::endl;
   }
   output_file.close();
   // Output the rotation curve
   ss.str(""); // clear
   ss << "rot_curve" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
   filename = ss.str();
   output_file.open(filename.c_str(),std::ios::trunc);
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   output_file << std::setprecision(15) << std::showpos;
   for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
       const PS::F64 r = std::sqrt(psys[i].pos * psys[i].pos);
       const PS::F64 R = std::sqrt(psys[i].pos.x * psys[i].pos.x
                                  +psys[i].pos.y * psys[i].pos.y);
       const PS::F64 v = std::sqrt(psys[i].vel * psys[i].vel);
       const PS::F64 T = 2.0 * math_const::pi * r / v;
       output_file << R << " " << v << " " << T << std::endl;
   }
   output_file.close();
   //PS::Finalize();
   //std::exit(0);
#endif
}

void setupIC(PS::ParticleSystem<FP_nbody>& psys_nbody,
             PS::ParticleSystem<FP_sph>& psys_sph,
             PS::DomainInfo& dinfo,
             PS::F64 & time_dump,
             PS::F64 & dt_dump,
             PS::F64 & time_end) {
    // Make an initial condition at MPI rank 0.
    PS::BOUNDARY_CONDITION bc;
    PS::F64ort pos_root_domain;
    if (PS::Comm::getRank() == 0) {
#if (INITIAL_CONDITION == 0)
        GalaxyIC(psys_nbody, psys_sph,
                 bc, pos_root_domain,
                 time_dump, dt_dump, time_end);
#elif (INITIAL_CONDITION == 1)
        ColdCollapseTestIC(psys_nbody, psys_sph,
                           bc, pos_root_domain,
                           time_dump, dt_dump, time_end);
#elif (INITIAL_CONDITION == 2)
        EvrardTestIC(psys_nbody, psys_sph,
                     bc, pos_root_domain,
                     time_dump, dt_dump, time_end, 1);
#elif (INITIAL_CONDITION == 3)
        MakeGlassIC(psys_nbody, psys_sph,
                    bc, pos_root_domain,
                    time_dump, dt_dump, time_end);
#else
#error Invalid IC number is specified.
#endif

        // Check the distribution of particles
        const std::string filename = "IC.txt";
        std::ofstream output_file;
        output_file.open(filename.c_str(),std::ios::trunc);
        output_file.setf(std::ios_base::scientific,
                         std::ios_base::floatfield);
        output_file << std::setprecision(15) << std::showpos;
        for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
            output_file << psys_nbody[i].pos.x << " "
                        << psys_nbody[i].pos.y << " "
                        << psys_nbody[i].pos.z << " "
                        << std::endl;
        }
        output_file << std::endl << std::endl; // separator
        for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
           output_file << psys_sph[i].pos.x << " "
                       << psys_sph[i].pos.y << " "
                       << psys_sph[i].pos.z << " "
                       << std::endl;
        }
        output_file.close();
    }
    else {
        psys_nbody.setNumberOfParticleLocal(0);
        psys_sph.setNumberOfParticleLocal(0);
    }

    // Broadcast 
    PS::Comm::broadcast(&bc, 1, 0);
    PS::Comm::broadcast(&pos_root_domain, 1, 0);
    PS::Comm::broadcast(&eps_grav, 1, 0);
    PS::Comm::broadcast(&dt_dump, 1, 0);
    PS::Comm::broadcast(&time_dump, 1, 0);
    PS::Comm::broadcast(&time_end, 1, 0);
    PS::Comm::broadcast(&dt_max, 1, 0);

    // Set the boundary condition and the size of the computational domain if needed.
    dinfo.setBoundaryCondition(bc);
    if (bc != PS::BOUNDARY_CONDITION_OPEN) {
        dinfo.setPosRootDomain(pos_root_domain.low_,
                               pos_root_domain.high_);
    }

    // Compute the average mass of SPH particles
    PS::F64 m_sum_loc = 0.0; 
    for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++)
        m_sum_loc += psys_sph[i].getCharge();
    mass_avg = PS::Comm::getSum(m_sum_loc) / psys_sph.getNumberOfParticleGlobal();

    // Output information to stdout
    if (PS::Comm::getRank() == 0)
        std::cout << "setupIC() is completed." << std::endl;
    //PS::Finalize();
    //std::exit(0);
}


PS::F64 getTimeStep(const PS::ParticleSystem<FP_nbody>& psys_nbody,
                    const PS::ParticleSystem<FP_sph>& psys_sph) {
   PS::F64 dt = DBL_MAX; 
   if (dt_max > 0.0) dt = dt_max;

   // Timescale for N-body system
#if defined(ENABLE_GRAVITY_INTERACT)
   for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
      const PS::F64 acc = std::sqrt(psys_nbody[i].acc * psys_nbody[i].acc);
      if (acc > 0.0)
          dt = std::min(dt, CFL_dyn * std::sqrt(eps_grav / acc));
   }
#endif

   // Timescale for SPH system
   for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
#if defined(ENABLE_GRAVITY_INTERACT)
      const PS::F64 acc = std::sqrt((psys_sph[i].acc_grav + psys_sph[i].acc_hydro)
                                   *(psys_sph[i].acc_grav + psys_sph[i].acc_hydro));
      if (acc > 0.0)
          dt = std::min(dt, CFL_dyn * std::sqrt(eps_grav / acc));
#endif
#if defined(ENABLE_HYDRO_INTERACT)
      dt = std::min(dt, psys_sph[i].dt);
#endif
   }
   return PS::Comm::getMinValue(dt);
}

void checkConservativeVariables(const PS::ParticleSystem<FP_nbody>& psys_nbody,
                                const PS::ParticleSystem<FP_sph>& psys_sph) {
    PS::F64    ekin_loc = 0.0;
    PS::F64    epot_loc = 0.0;
    PS::F64    eth_loc  = 0.0; 
    PS::F64vec mom_loc  = 0.0; 
    for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        ekin_loc += 0.5 * psys_nbody[i].mass * psys_nbody[i].vel * psys_nbody[i].vel;
        epot_loc += 0.5 * psys_nbody[i].mass * (psys_nbody[i].pot + psys_nbody[i].mass / eps_grav);
        mom_loc  += psys_nbody[i].mass * psys_nbody[i].vel;
    }
    for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
        ekin_loc += 0.5 * psys_sph[i].mass * psys_sph[i].vel * psys_sph[i].vel;
        epot_loc += 0.5 * psys_sph[i].mass * (psys_sph[i].pot_grav + psys_sph[i].mass / eps_grav);
        eth_loc  += psys_sph[i].mass * psys_sph[i].eng;
        mom_loc  += psys_sph[i].mass * psys_sph[i].vel;
    }
    PS::F64 ekin    = PS::Comm::getSum(ekin_loc);
    PS::F64 epot    = PS::Comm::getSum(epot_loc);
    PS::F64 eth     = PS::Comm::getSum(eth_loc);
    PS::F64vec mom  = PS::Comm::getSum(mom_loc);

    static bool is_initialized = false;
    static PS::F64 emech_ini, etot_ini;
    if (is_initialized == false) {
        emech_ini = ekin + epot;
        etot_ini  = ekin + epot + eth;
        is_initialized = true;
    }

    if (PS::Comm::getRank() == 0){
        const PS::F64 emech = ekin + epot;
        const PS::F64 etot  = ekin + epot + eth;
        const PS::F64 relerr_mech = std::fabs((emech - emech_ini)/emech_ini);
        const PS::F64 relerr_tot  = std::fabs((etot  - etot_ini)/etot_ini);
        std::cout << "-------------------------" << std::endl;
        std::cout << "E_kin  = " << ekin  << std::endl;
        std::cout << "E_pot  = " << epot  << std::endl;
        std::cout << "E_th   = " << eth   << std::endl;
        std::cout << "E_mech = " << emech << " (" << relerr_mech << ")" << std::endl;
        std::cout << "E_tot  = " << etot  << " (" << relerr_tot  << ")" << std::endl;
        std::cout << "Mom (x) = " << mom.x << std::endl;
        std::cout << "Mom (y) = " << mom.y << std::endl;
        std::cout << "Mom (z) = " << mom.z << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
}

void checkDensityFluctuation(PS::ParticleSystem<FP_sph>& psys) {
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    const PS::S64 n_glb = psys.getNumberOfParticleGlobal();
    // Compute the average of density
    PS::F64 tmp = 0.0;
    for (PS::S32 i = 0; i < n_loc; i++) tmp += psys[i].dens;
    const PS::F64 dens_avg = PS::Comm::getSum(tmp) / n_glb;
    // Compute the dispersion of density
    tmp = 0.0;
    for (PS::S32 i = 0; i < n_loc; i++) 
        tmp += std::pow(psys[i].dens - dens_avg, 2.0);
    const PS::F64 dens_disp = std::sqrt(PS::Comm::getSum(tmp)/n_glb);
    // Output the status 
    const PS::F64 fluc_str = dens_disp / dens_avg;
    if (PS::Comm::getRank() == 0) {
        std::cout << "---------------------------" << std::endl;
        std::cout << "avg.       = " << dens_avg << std::endl;
        std::cout << "disp.      = " << dens_disp << std::endl;
        std::cout << "disp./avg. = " << fluc_str << std::endl;
    }
    // Output data and end the simulation if the fluctuation is small
    const PS::F64 eps = 3.05e-3;
    if (fluc_str < eps) {
        char filename[256];
        sprintf(filename, "result/glass_data.txt");
        psys.writeParticleBinary(filename, &FP_sph::writeBinaryPos);
        if (PS::Comm::getRank() == 0) {
            std::cout << "A glass-like distribution is obtained." << std::endl;
            std::cout << "The particle position data is output as file " 
                      << filename << std::endl;
        }
        PS::Finalize();
        std::exit(0);
    }

}


int main(int argc, char* argv[]){
    // Configure stdout & stderr
    std::cout << std::setprecision(15);
    std::cerr << std::setprecision(15);

    // Initialize FDPS
    PS::Initialize(argc, argv);

    // Make a directory
    char dir_name[1024];
    sprintf(dir_name,"./result");
    makeOutputDirectory(dir_name);

    // Display # of MPI processes and threads
    PS::S32 nprocs = PS::Comm::getNumberOfProc();
    PS::S32 nthrds = PS::Comm::getNumberOfThread();
    if (PS::Comm::getRank() == 0) {
        std::cout << "===========================================" << std::endl
                  << " This is a sample program of "               << std::endl
                  << " Nbody + SPH particle system on FDPS!"       << std::endl
                  << " # of processes is " << nprocs               << std::endl
                  << " # of thread is    " << nthrds               << std::endl
                  << "===========================================" << std::endl;
    }

    // Make instances of ParticleSystem and initialize them
    PS::ParticleSystem<FP_nbody> psys_nbody;
    PS::ParticleSystem<FP_sph> psys_sph;
    psys_nbody.initialize();
    psys_sph.initialize();

    // Make an instance of DomainInfo and initialize it
    PS::DomainInfo dinfo;
    dinfo.initialize();

    // Define local variables
    PS::F64 dt, time_dump, dt_dump, time_end;

    // Make an initial condition and initialize the particle systems
    setupIC(psys_nbody, psys_sph, dinfo, time_dump, dt_dump, time_end);

    // Perform domain decomposition 
    dinfo.collectSampleParticle(psys_nbody);
    dinfo.collectSampleParticle(psys_sph,false);
    dinfo.decomposeDomain();
    
    // Perform particle exchange
    psys_nbody.exchangeParticle(dinfo);
    psys_sph.exchangeParticle(dinfo);

    // Make tree structures
    const PS::S64 numPtclSPH = std::max(psys_sph.getNumberOfParticleLocal(),1);
    const PS::S64 numPtclAll = psys_nbody.getNumberOfParticleLocal() + numPtclSPH;

    const PS::F32 theta_grav = 0.5;
    PS::TreeForForceLong<Force_grav, EP_grav, EP_grav>::Monopole tree_grav;
    tree_grav.initialize(3 * numPtclAll, theta_grav);

    PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather tree_dens;
    tree_dens.initialize(3 * numPtclSPH);

    PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry tree_hydro;
    tree_hydro.initialize(3 * numPtclSPH);

#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_open();
    g5_set_eps_to_all(eps_grav);
#endif

    // Peform force calculations 
    //- Gravity calculations
#if defined(ENABLE_GRAVITY_INTERACT)
    tree_grav.setParticleLocalTree(psys_nbody);
    tree_grav.setParticleLocalTree(psys_sph,false);
    tree_grav.calcForceMakingTree(CalcGravity<EP_grav>,
                                  CalcGravity<PS::SPJMonopole>,
                                  dinfo);
    for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        psys_nbody[i].copyFromForce(tree_grav.getForce(i));
    }
    const PS::S32 offset = psys_nbody.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
        psys_sph[i].copyFromForce(tree_grav.getForce(i+offset));
    }
#endif

    //- SPH calculations
#if defined(ENABLE_HYDRO_INTERACT)
    calcDensity(psys_sph, dinfo, tree_dens);
#if defined(USE_ENTROPY)
    setEntropy(psys_sph);
#endif
    setPressure(psys_sph);
    tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_sph, dinfo);
#endif

    // Set the initial velocity of gas particles
#if defined(SET_CIRCULAR_VELOCITY)
    setCircularVelocity(psys_sph);
#endif

    // Get timestep
    dt = getTimeStep(psys_nbody, psys_sph);

    // Calculate energies 
    checkConservativeVariables(psys_nbody, psys_sph);

    // Main loop for time integration
    PS::S32 nstep = 1;
    for (PS::F64 time = 0.0; time < time_end; time += dt, nstep++){
        if (PS::Comm::getRank() == 0) {
            std::cout << "nstep = " << nstep 
                      << " dt = " << dt 
                      << " time = " << time 
                      << " time_end = " << time_end
                      << std::endl;
        }

        // Leap frog: Initial Kick & Full Drift
        InitialKick(psys_nbody, dt);
        InitialKick(psys_sph, dt);
        FullDrift(psys_nbody, dt);
        FullDrift(psys_sph, dt);
        if (dinfo.getBoundaryCondition() != PS::BOUNDARY_CONDITION_OPEN) {
            psys_nbody.adjustPositionIntoRootDomain(dinfo);
            psys_sph.adjustPositionIntoRootDomain(dinfo);
        }

        // Leap frog: Predict
        Predict(psys_sph, dt);

        // Perform domain decomposition again
        dinfo.collectSampleParticle(psys_nbody);
        dinfo.collectSampleParticle(psys_sph,false);
        dinfo.decomposeDomain();

        // Exchange the particles between the (MPI) processes
        psys_nbody.exchangeParticle(dinfo);
        psys_sph.exchangeParticle(dinfo);

        // Peform force calculations
        PS::F64 t_start; 
        //- Gravity calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
#if defined(ENABLE_GRAVITY_INTERACT)
        tree_grav.setParticleLocalTree(psys_nbody);
        tree_grav.setParticleLocalTree(psys_sph,false);
        tree_grav.calcForceMakingTree(CalcGravity<EP_grav>,
                                      CalcGravity<PS::SPJMonopole>,
                                      dinfo);
        for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
            psys_nbody[i].copyFromForce(tree_grav.getForce(i));
        }
        const PS::S32 offset = psys_nbody.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_sph.getNumberOfParticleLocal(); i++) {
            psys_sph[i].copyFromForce(tree_grav.getForce(i+offset));
        }
#endif
        PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) std::cout << "t_grav = " << (PS::GetWtime() - t_start) << std::endl;
        //- SPH calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
#if defined(ENABLE_HYDRO_INTERACT)
        calcDensity(psys_sph, dinfo, tree_dens);
        setPressure(psys_sph);
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_sph, dinfo);
#endif
        PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) std::cout << "t_hydro = " << (PS::GetWtime() - t_start) << std::endl;

        // Get a new timestep
        dt = getTimeStep(psys_nbody, psys_sph);

        // Leap frog: Final Kick
        FinalKick(psys_nbody, dt);
        FinalKick(psys_sph, dt);

        // Output result files
        if (time > time_dump){
           FileHeader header_nbody, header_sph;
           header_nbody.time    = time;
           header_nbody.numPtcl = psys_nbody.getNumberOfParticleGlobal();
           header_sph.time      = time;
           header_sph.numPtcl   = psys_sph.getNumberOfParticleGlobal();
           char filename[256];
           static int ndump = 1;
           sprintf(filename, "result/nbody%05d.txt", ndump);
           psys_nbody.writeParticleAscii(filename, header_nbody);
           sprintf(filename, "result/sph%05d.txt", ndump);
           psys_sph.writeParticleAscii(filename, header_sph);
           if (PS::Comm::getRank() == 0){
              std::cout << "============================================" << std::endl;
              std::cout << "output " << filename << " at time = " << time << std::endl;
              std::cout << "============================================" << std::endl;
           }
           time_dump += dt_dump;
           ndump++;
        }

        // Calculate energies
        checkConservativeVariables(psys_nbody, psys_sph);

        // Check the amplitude of density fluctuation
#if defined(CHECK_DENSITY_FLUCTUATION)
        if (nstep % 100 == 0) 
            checkDensityFluctuation(psys_sph);
#endif

    }

#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_close();
#endif
    // Finalize FDPS
    PS::Finalize();
    return 0;
}

