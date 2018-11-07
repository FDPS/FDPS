/* Headers */
#include "mathematical_constants.h"
#include "physical_constants.h"

/* Class definitions */
// The following two classes are used in function readTipsyFile,
// which is used to read particle data created by MAGI.
class MAGI_Tipsy_header {
public:
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
};

class MAGI_Tipsy_particle {
public:
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    int idx;
};

void readTipsyFile(std::string& input_file_name,
                   PS::ParticleSystem<FP_nbody>& psys) {
    // This function is used to read particle data created by
    // MAGI (https://bitbucket.org/ymiki/magi). The particle
    // data must be in the TIPSY format.

    // Read buffers
    MAGI_Tipsy_header header;
    std::vector<MAGI_Tipsy_particle> ptcl;

    // Read particle data
    std::ifstream input_file;
    input_file.open(input_file_name.c_str(), std::ios::in | std::ios::binary);
    if (input_file) {
        input_file.read((char *)&header, sizeof(MAGI_Tipsy_header));
        ptcl.resize(header.nbodies);
        input_file.read((char *)&ptcl[0], sizeof(MAGI_Tipsy_particle)*header.nbodies);
    }
    else {
        std::cout << "cannot open the file " << input_file_name << std::endl;
        PS::Abort(-1);
        std::exit(1);
    }
    input_file.close();

    // Copy particle data
    psys.setNumberOfParticleLocal(header.nbodies);
    for (PS::S32 i=0; i<header.nbodies; i++) {
        psys[i].mass  = ptcl[i].mass;
        psys[i].pos.x = ptcl[i].pos[0];
        psys[i].pos.y = ptcl[i].pos[1];
        psys[i].pos.z = ptcl[i].pos[2];
        psys[i].vel.x = ptcl[i].vel[0];
        psys[i].vel.y = ptcl[i].vel[1];
        psys[i].vel.z = ptcl[i].vel[2];
    }

}

PS::F64 getPosCellCenter(const PS::F64 pos_left_bound,
                         const PS::F64 pos_right_bound,
                         const PS::S32 number_of_cells,
                         const PS::S32 i) {
    assert(0 <= i && i < number_of_cells);
    assert(pos_left_bound < pos_right_bound);
    const PS::F64 dx = (pos_right_bound - pos_left_bound) / number_of_cells;
    if ( number_of_cells % 2 == 0) {
        // # of cells is even.
        if (i < number_of_cells/2) {
            return pos_left_bound + dx * (i + 0.5);
        }
        else {
            return pos_right_bound - dx * (number_of_cells - i - 0.5);
        }
    } 
    else {
        // # of cells is odd.
        const PS::F64 center = 0.5 * (pos_left_bound + pos_right_bound);
        if (i < number_of_cells/2) {
            return center - dx * (number_of_cells/2 - i);
        }
        else {
            return center + dx * (i - number_of_cells/2);
        }
    }
}

void ColdCollapseTestIC(PS::ParticleSystem<FP_nbody>& psys_nbody,
                        PS::ParticleSystem<FP_sph>& psys_sph,
                        PS::BOUNDARY_CONDITION& bc,
                        PS::F64ort& pos_root_domain,
                        PS::F64 & time_dump,
                        PS::F64 & dt_dump,
                        PS::F64 & time_end) {
    // Model parameters
    const PS::F64 M_nbody = 0.5;
    const PS::F64 M_sph   = 0.5;
    const PS::S64 N_nbody = 512;
    const PS::S64 N_sph   = 512;
    const PS::F64 R_nbody = 3.0;
    const PS::F64 R_sph   = 3.0;
    const PS::F64 E_bind_nbody = - 3.0 * M_nbody * M_nbody / 5.0 / R_nbody;
    const PS::F64 E_bind_sph   = - 3.0 * M_sph * M_sph / 5.0 / R_sph;
    const PS::F64 virial_ratio = 0.5;
    // Initialize pseudorandom number generator
    PS::MTTS mt;
    mt.init_genrand(0);
    // Place Nbody particles
    psys_nbody.setNumberOfParticleLocal(N_nbody);
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_nbody[i].id = (i + 1);
        psys_nbody[i].mass = M_nbody / N_nbody;
        for(;;) {
            psys_nbody[i].pos.x = (2.0 * mt.genrand_res53() - 1.0) * R_nbody;
            psys_nbody[i].pos.y = (2.0 * mt.genrand_res53() - 1.0) * R_nbody;
            psys_nbody[i].pos.z = (2.0 * mt.genrand_res53() - 1.0) * R_nbody;
            const PS::F64 r2 = psys_nbody[i].pos * psys_nbody[i].pos;
            if (r2 <= R_nbody * R_nbody) break;
        }
        psys_nbody[i].vel = 0.0;
        psys_nbody[i].acc = 0.0;
        psys_nbody[i].pot = 0.0;
    }
    // Place SPH particles
    psys_sph.setNumberOfParticleLocal(N_sph);
    for (PS::S32 i = 0; i < N_sph; i++) {
        psys_sph[i].id = (i + 1) + N_nbody;
        psys_sph[i].mass = M_sph / N_sph;
        for(;;) {
            psys_sph[i].pos.x = (2.0 * mt.genrand_res53() - 1.0) * R_sph;
            psys_sph[i].pos.y = (2.0 * mt.genrand_res53() - 1.0) * R_sph;
            psys_sph[i].pos.z = (2.0 * mt.genrand_res53() - 1.0) * R_sph;
            const PS::F64 r2 = psys_sph[i].pos * psys_sph[i].pos;
            if (r2 <= R_sph * R_sph) break;
        }
        psys_sph[i].vel       = 0.0;
        psys_sph[i].acc_grav  = 0.0;
        psys_sph[i].pot_grav  = 0.0;
        psys_sph[i].acc_hydro = 0.0;
        psys_sph[i].eng = virial_ratio * std::abs(E_bind_sph) / M_sph;  // specific thermal energy
        psys_sph[i].smth = R_sph * std::pow((PS::F64)N_neighbor/(PS::F64)N_sph, 1.0/3.0); // smoothing length
        // Note that entropy is determined in setEntropy().
    }
    // Set boundary condition
    bc = PS::BOUNDARY_CONDITION_OPEN;
    // Set gravitational softening
    eps_grav = 0.01 * R_nbody;
    // Set I/O intervals
    const PS::F64 dens_nbody = 3.0 * M_nbody / (4.0 * math_const::pi * R_nbody * R_nbody * R_nbody);
    const PS::F64 dens_sph   = 3.0 * M_sph / (4.0 * math_const::pi * R_sph * R_sph * R_sph);
    const PS::F64 unit_dens = dens_nbody + dens_sph;
    const PS::F64 unit_time = std::sqrt(3.0 * math_const::pi/ (32.0 * unit_dens));
    dt_dump   = 0.1 * unit_time;
    time_dump = dt_dump;
    time_end  = 2.0 * unit_time;
    if (PS::Comm::getRank() == 0) {
        std::cout << "unit_dens = " << unit_dens << std::endl;
        std::cout << "unit_time = " << unit_time << std::endl;
    }
    // Set maximum timestep
    dt_max = 1.0e-3;

}

void EvrardTestIC(PS::ParticleSystem<FP_nbody>& psys_nbody,
                  PS::ParticleSystem<FP_sph>& psys_sph,
                  PS::BOUNDARY_CONDITION& bc,
                  PS::F64ort& pos_root_domain,
                  PS::F64 & time_dump,
                  PS::F64 & dt_dump,
                  PS::F64 & time_end,
                  PS::S32 gen_mode) {
    // Place Nbody particles
    // here, we place only one dummy, mass-less particle.
    const PS::S64 N_nbody = 1;
    psys_nbody.setNumberOfParticleLocal(N_nbody);
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_nbody[i].id = i;
        psys_nbody[i].mass = 0.0;
        psys_nbody[i].pos  = 0.0;
        psys_nbody[i].vel  = 0.0;
        psys_nbody[i].acc  = 0.0;
        psys_nbody[i].pot  = 0.0;
    }
    // Place SPH particles
    PS::S64 N_sph = 0;
    const PS::F64 M_sph = 1.0;
    const PS::F64 radius_of_sphere = 1.0;
    if (gen_mode == 0) {
        // In this mode, we create an initial distribution of particles
        // by rescaling the positions of particles which are placed in a grid.
        // (1) Count # of particles in a sphere of radius 1
        const PS::S64 N_ptcl_per_side = 64;
        PS::F64 x,y,z;
        // x-loop
        for (PS::S32 i = 0; i < N_ptcl_per_side; i++) {
            x = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, i);
            // y-loop
            for (PS::S32 j = 0; j < N_ptcl_per_side; j++) {
                y = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, j);
                // z-loop
                for (PS::S32 k = 0; k < N_ptcl_per_side; k++) {
                    z = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, k);
                    const PS::F64 r = std::sqrt(x*x + y*y + z*z);
                    if (r <= radius_of_sphere) N_sph++;
                }
            }
        }
        assert(N_sph > 0);
        psys_sph.setNumberOfParticleLocal(N_sph);
        // (2) Actually place particles
        PS::S64 id = 0;
        // x-loop
        for (PS::S32 i = 0; i < N_ptcl_per_side; i++) {
            x = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, i);
            // y-loop
            for (PS::S32 j = 0; j < N_ptcl_per_side; j++) {
                y = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, j);
                // z-loop
                for (PS::S32 k = 0; k < N_ptcl_per_side; k++) {
                    z = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, k);
                    const PS::F64 r = std::sqrt(x*x + y*y + z*z);
                    if (r <= radius_of_sphere) {
                        const PS::F64 r_new = radius_of_sphere * std::pow(r/radius_of_sphere,1.5);
                        psys_sph[id].id = id + N_nbody;
                        psys_sph[id].mass = M_sph / N_sph;
                        psys_sph[id].pos.x = (r_new / r) * x;
                        psys_sph[id].pos.y = (r_new / r) * y;
                        psys_sph[id].pos.z = (r_new / r) * z;
                        psys_sph[id].vel       = 0.0;
                        psys_sph[id].acc_grav  = 0.0;
                        psys_sph[id].pot_grav  = 0.0;
                        psys_sph[id].acc_hydro = 0.0;
                        psys_sph[id].eng = 0.05 * M_sph / radius_of_sphere;
                        psys_sph[id].smth = radius_of_sphere * std::pow((PS::F64)N_neighbor/(PS::F64)N_sph, 1.0/3.0);
                        // Note that entropy is determined in setEntropy().
                        id++;
                    }
                }
            }
        }
    }
    else if (gen_mode == 1) {
        // In this mode, we set an initial distribution of particles by reading a file.
        std::vector<PS::F64vec> ptcl;
        // Read position data
        const std::string filename = "result/glass_data.txt";
        std::ifstream data_file; 
        data_file.open(filename.c_str(), std::ios::in | std::ios::binary);
        if (data_file) {
           while(true) {
              PS::F64vec pos;
              data_file.read((char *)&pos, sizeof(PS::F64vec));
              if (data_file.eof()) break;
              ptcl.push_back(pos);
           }
        }
        data_file.close();
        const PS::S64 n_ptcl_in = ptcl.size();
        if (PS::Comm::getRank() == 0)
            std::cout << n_ptcl_in << " particles are read." << std::endl;
        // Count # of particles in a sphere of radius 1
        N_sph = 0;
        for (PS::S32 i = 0; i < n_ptcl_in; i++) {
            const PS::F64 r2 = ptcl[i] * ptcl[i];
            if (r2 < 1.0) N_sph++;
        }
        if (PS::Comm::getRank() == 0)
            std::cout << N_sph << " particles of them will be used to make a Evrard sphere." << std::endl;
        // Place SPH particles
        psys_sph.setNumberOfParticleLocal(N_sph);
        PS::S64 j = -1;
        for (PS::S64 i = 0; i < N_sph; i++) {
            psys_sph[i].id = (i + 1);
            psys_sph[i].mass = M_sph / N_sph;
            PS::F64 r2;
            for(;;) {
                j++;
                r2 = ptcl[j] * ptcl[j];
                if (r2 < 1.0) break;
            }
            const PS::F64 r = std::sqrt(r2);
            if (r > 0.0) {
                const PS::F64 r_new = radius_of_sphere * std::pow(r,1.5);
                psys_sph[i].pos = (r_new / r) * ptcl[j];
            } else {
                psys_sph[i].pos = ptcl[j];
            }
            psys_sph[i].vel       = 0.0;
            psys_sph[i].acc_grav  = 0.0;
            psys_sph[i].pot_grav  = 0.0;
            psys_sph[i].acc_hydro = 0.0;
            psys_sph[i].eng = 0.05 * M_sph / radius_of_sphere;
            psys_sph[i].smth = radius_of_sphere * std::pow((PS::F64)N_neighbor/(PS::F64)N_sph, 1.0/3.0);
            // Note that entropy is determined in setEntropy().
        }
    }
    else {
        if (PS::Comm::getRank() == 0) std::cout << "Given gen_mode is not supported." << std::endl;
        PS::Abort();
        std::exit(-1);
    }
    // Set boundary condition
    bc = PS::BOUNDARY_CONDITION_OPEN;
    // Set gravitational softening
    eps_grav = 1.0e-4 * radius_of_sphere;
    // Set I/O intervals
    const PS::F64 unit_dens = 3.0 * M_sph / (4.0 * math_const::pi * std::pow(radius_of_sphere, 3.0));
    const PS::F64 unit_time = std::sqrt(math_const::pi * math_const::pi / 8.0)
                            * std::pow(radius_of_sphere, 1.5)
                            / std::sqrt(M_sph);
    dt_dump = 0.1 * unit_time;
    time_dump = dt_dump;
    time_end = unit_time;
    if (PS::Comm::getRank() == 0) {
        std::cout << "unit_dens = " << unit_dens << std::endl;
        std::cout << "unit_time = " << unit_time << std::endl;
    }
    // Set maximum timestep
    dt_max = 1.0e-3;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for the Evrard test is made." << std::endl;
        std::cout << "(N_sph = " << N_sph << ")" << std::endl;
    }
}

void MakeGlassIC(PS::ParticleSystem<FP_nbody>& psys_nbody,
                 PS::ParticleSystem<FP_sph>& psys_sph,
                 PS::BOUNDARY_CONDITION& bc,
                 PS::F64ort& pos_root_domain,
                 PS::F64 & time_dump,
                 PS::F64 & dt_dump,
                 PS::F64 & time_end) {
    // Model parameters
    const PS::S64 N_nbody = 1; // dummy value
    const PS::S64 N_sph   = std::pow(2,18);
    // Initialize pseudorandom number generator
    PS::MTTS mt;
    mt.init_genrand(0);
    // Place Nbody particles
    psys_nbody.setNumberOfParticleLocal(N_nbody);
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_nbody[i].id = (i + 1);
        psys_nbody[i].mass  = 0.0;
        psys_nbody[i].pos   = 0.0;
        psys_nbody[i].vel   = 0.0;
    }
    // Place SPH particles
    psys_sph.setNumberOfParticleLocal(N_sph);
    for (PS::S32 i = 0; i < N_sph; i++) {
        psys_sph[i].id = (i + 1) + N_nbody;
        psys_sph[i].mass = 8.0/N_sph;
        PS::F64 x,y,z;
        for(;;) {
            x = (2.0 * mt.genrand_res53() - 1.0);
            y = (2.0 * mt.genrand_res53() - 1.0);
            z = (2.0 * mt.genrand_res53() - 1.0);
            if ((-1.0 <= x) && (x < 1.0) &&
                (-1.0 <= y) && (y < 1.0) &&
                (-1.0 <= z) && (z < 1.0)) break;
        }
        psys_sph[i].pos.x = x;
        psys_sph[i].pos.y = y;
        psys_sph[i].pos.z = z;
        psys_sph[i].vel = 0.0;
        psys_sph[i].eng = 1.0;  
        psys_sph[i].smth = 2.0 * std::pow((PS::F64)N_neighbor/(PS::F64)N_sph, 1.0/3.0); // smoothing length
        // [Notes]
        //   (1) The value of the specific thermal energy is chosen
        //       so that the sound speed is nearly equal to 1.
        //   (2) the value of the entropy function is determined in setEntropy().
    }
    // Set boundary condition
    bc = PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
    pos_root_domain.low_  = (PS::F64vec)(-1.0, -1.0, -1.0);
    pos_root_domain.high_ = (PS::F64vec)( 1.0,  1.0,  1.0);
    // Set gravitational softening
    eps_grav = 0.01;
    // Set I/O intervals
    const PS::F64 tcross = 2.0 * std::sqrt(3.0); 
    dt_dump   = 4.0 * tcross;
    time_dump = dt_dump;
    time_end  = 64.0 * tcross;
    if (PS::Comm::getRank() == 0) {
        std::cout << "The sound crossing time = " << tcross << std::endl;
    }
    // Set maximum timestep
    dt_max = 1.0e-2;

}


void GalaxyIC(PS::ParticleSystem<FP_nbody>& psys_nbody,
              PS::ParticleSystem<FP_sph>& psys_sph,
              PS::BOUNDARY_CONDITION& bc,
              PS::F64ort& pos_root_domain,
              PS::F64 & time_dump,
              PS::F64 & dt_dump,
              PS::F64 & time_end) {
    // Define the code units of MAGI
    //    [Important]
    //    (1) The values MUST BE consistent with "the computational units"
    //        written in the file doc/unit.txt, which is output by MAGI
    //        when we create a particle data with MAGI.
    //    (2) The MAGI's code units are DIFFERENT for unit systems
    //        a user choose. For detail, read Section "Unit systems in
    //        inc/constants.[c h]" in https://bitbucket.org/ymiki/magi.
    //    (3) In this sample code, "Galactic scale" unit is adopted.
    //        It is consistent with ./magi_data/cfg/Galaxy.cfg.
    const PS::F64 magi_unit_mass = 1.0e8 * phys_const::Msolar;
    const PS::F64 magi_unit_leng = phys_const::kpc;
    const PS::F64 magi_unit_time = 1.0e2 * phys_const::Myr;
    const PS::F64 magi_unit_velc = magi_unit_leng/magi_unit_time;
    // Initialize pseudorandom number generator
    PS::MTTS mt;
    mt.init_genrand(0);
    // Place Nbody particles
    std::string filename = "./magi_data/dat/Galaxy.tipsy";
    readTipsyFile(filename, psys_nbody);
    const PS::S64 N_nbody = psys_nbody.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_nbody[i].id = i;
        psys_nbody[i].acc  = 0.0;
        psys_nbody[i].pot  = 0.0;
    }
    // Place SPH particles to form an exponential-disk
    const PS::S64 N_sph = std::pow(2,18);
    const PS::F64 M_gas = 1.0e10 * phys_const::Msolar;
    const PS::F64 Rs = 7.0 * phys_const::kpc; // scale radius
    const PS::F64 Rt = 12.5 * phys_const::kpc; // truncation radius
    const PS::F64 zd = 400.0 * phys_const::pc; // scale height
    const PS::F64 zt = 1.0 * phys_const::kpc; // truncation height
    const PS::F64 temp = 1.0e4; // gas temperature
    const PS::F64 mu = 0.5; // mean molecular weight relative to the mass of hydrogen
    psys_sph.setNumberOfParticleLocal(N_sph);
    PS::S64 id = 0;
    for (PS::S32 i = 0; i < N_sph; i++) {
        // First make a uniform disk with a finite thickness
        PS::F64 x,y,z,R2;
        for(;;) {
            x = (2.0 * mt.genrand_res53() - 1.0) * Rt;
            y = (2.0 * mt.genrand_res53() - 1.0) * Rt;
            z = (2.0 * mt.genrand_res53() - 1.0) * zt;
            R2 = (x * x + y * y);
            if ((R2 < Rt * Rt) && (std::abs(z) < zt)) break;
        }
        const PS::F64 R = std::sqrt(R2);
        // Then re-scale particle position to generate an exponential disk
        PS::F64 R_low = 0.0, R_high = Rt, R_new;
        const PS::F64 eps = 1.0e-6;
        const PS::F64 val_trgt = (R/Rt)*(R/Rt);
        for(;;) {
            R_new = 0.5 * (R_low + R_high);
            const PS::F64 val = (1.0 - (R_new/Rs + 1.0)*std::exp(-R_new/Rs)) 
                              / (1.0 - (   Rt/Rs + 1.0)*std::exp(-   Rt/Rs));
            if (val < val_trgt) R_low  = R_new;
            if (val > val_trgt) R_high = R_new;
            const PS::F64 reldiff = 2.0 * std::abs(R_low - R_high)/(R_low + R_high);
            if (reldiff < eps) {
                R_new = 0.5 * (R_low + R_high);
                break;
            }
        }
        PS::F64 z_new;
        if (z >= 0.0) z_new = - zd * std::log(1.0 - (z/zt) * (1.0 - std::exp(-zt/zd)));
        else z_new = zd * std::log(1.0 + (z/zt) * (1.0 - std::exp(-zt/zd)));
        // Set  
        psys_sph[id].id = id + N_nbody;
        psys_sph[id].mass = M_gas / N_sph;
        psys_sph[id].pos.x = (R_new / R) * x;
        psys_sph[id].pos.y = (R_new / R) * y;
        psys_sph[id].pos.z = z_new;
        psys_sph[id].vel       = 0.0;
        psys_sph[id].acc_grav  = 0.0;
        psys_sph[id].pot_grav  = 0.0;
        psys_sph[id].acc_hydro = 0.0;
        psys_sph[id].eng = (phys_const::kBoltz * temp)/((specific_heat_ratio - 1.0) * mu * phys_const::Mhydrogen);
        psys_sph[id].smth = std::pow(Rt*Rt*zt, 1.0/3.0) * std::pow((PS::F64)N_neighbor/(PS::F64)N_sph, 1.0/3.0);
        // Note that entropy is determined in setEntropy().
        id++;
    }
    // Unit convertion (MAGI unit, CGS unit -> G=M=R=1 system)
    PS::F64 unit_mass = 0.0, unit_leng = 0.0;
    PS::F64 unit_time, unit_velc, unit_eng;
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_nbody[i].mass *= magi_unit_mass;
        psys_nbody[i].pos  *= magi_unit_leng;
        psys_nbody[i].vel  *= magi_unit_velc;
        unit_mass += psys_nbody[i].mass;
        const PS::F64 r = std::sqrt(psys_nbody[i].pos * psys_nbody[i].pos);
        if (r > unit_leng) unit_leng = r;
    }
    std::cout << "Total mass in N-body particles = "
              << unit_mass/phys_const::Msolar << " [Msolar]" << std::endl;
    for (PS::S32 i = 0; i< N_sph; i++) {
        unit_mass += psys_sph[i].mass;
        const PS::F64 r = std::sqrt(psys_sph[i].pos * psys_sph[i].pos);
        if (r > unit_leng) unit_leng = r;
    }
    unit_time = std::sqrt(std::pow(unit_leng, 3.0)/(phys_const::Ggrav * unit_mass));
    unit_velc = unit_leng / unit_time;
    unit_eng  = unit_mass * unit_velc * unit_velc;
    if (PS::Comm::getRank() == 0) {
        std::cout << "unit_mass = " << unit_mass << " [g]" 
                  << " = " << unit_mass/phys_const::Msolar << " [Msolar]" 
                  << std::endl;
        std::cout << "unit_leng = " << unit_leng << " [cm]" 
                  << " = " << unit_leng/phys_const::kpc << " [kpc]"
                  << std::endl;
        std::cout << "unit_time = " << unit_time << " [s]" 
                  << " = " << unit_time/phys_const::Gyr << " [Gyr]"
                  << std::endl;
        std::cout << "unit_velc = " << unit_velc << " [cm/s]"
                  << " = " << unit_velc/phys_const::km << " [km/s]"
                  << std::endl;
    }
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_nbody[i].mass /= unit_mass;
        psys_nbody[i].pos  /= unit_leng;
        psys_nbody[i].vel  /= unit_velc;
    }
    for (PS::S32 i = 0; i < N_sph; i++) {
        psys_sph[i].mass /= unit_mass;
        psys_sph[i].pos  /= unit_leng;
        psys_sph[i].vel  /= unit_velc;
        psys_sph[i].smth /= unit_leng;
        psys_sph[i].eng  /= (unit_eng/unit_mass);
    }
    // Set boundary condition
    bc = PS::BOUNDARY_CONDITION_OPEN;
    // Set the other parameters
    eps_grav = 1.0e-3;
    // Set I/O intervals
    dt_dump = 0.01;
    time_dump = dt_dump;
    time_end = 1.0;
    // Set maximum timestep
    dt_max = 1.0e-4;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for isolated galaxy simulation is made." << std::endl;
        std::cout << "N_nbody = " << N_nbody << ", N_sph = " << N_sph << std::endl;
    }


#if 1
    // [for DEBUG]
    // Compute the surface gas density and output it.
    // Also check the distribution of particles w.r.t. the z coordinate.
    {
        // Set the resolution of bins, etc.
        const PS::S32 nbin = 64;
        const PS::F64 safety = 1.001;
        // Compute the maxium clyndrical radius of gas particles
        PS::F64 Rmax = 0.0;
        for (PS::S32 i = 0; i < N_sph; i++) {
            const PS::F64 R = std::sqrt(psys_sph[i].pos.x * psys_sph[i].pos.x
                                       +psys_sph[i].pos.y * psys_sph[i].pos.y);
            if (R > Rmax) Rmax = R;
        }
        // Compute the surface density
        const PS::F64 dR = (safety * Rmax)/nbin;
        PS::F64 Sigma[nbin] = {};
        for (PS::S32 i = 0; i < N_sph; i++) {
            const PS::F64 R = std::sqrt(psys_sph[i].pos.x * psys_sph[i].pos.x
                                       +psys_sph[i].pos.y * psys_sph[i].pos.y);
            const PS::S32 indx = R/dR;
            Sigma[indx] += psys_sph[i].mass;
        }
        for (PS::S32 n = 0; n < nbin; n++) 
            Sigma[n] /= 2.0 * math_const::pi * ((n+1)*(n+1)-n*n)*dR*dR;
        // Output the surface density
        std::string filename = "Sigma.txt";
        std::ofstream output_file;
        output_file.open(filename.c_str(),std::ios::trunc);
        output_file.setf(std::ios_base::scientific,
                         std::ios_base::floatfield);
        output_file << std::setprecision(15) << std::showpos;
        for (PS::S32 n = 0; n < nbin; n++) {
            output_file << (n + 0.5) * dR << " " << Sigma[n] << " " << std::endl;
        }
        output_file.close();
        // Compute the distribution function 
        const PS::F64 zmax = zt/unit_leng; 
        const PS::F64 dz = (safety * 2.0 * zmax)/nbin;
        PS::S32 dist_func[nbin] = {};
        for (PS::S32 i = 0; i < N_sph; i++) {
            const PS::S32 indx = (psys_sph[i].pos.z + safety * zmax)/dz;
            dist_func[indx]++;
        }
        // Output the distribution function
        filename = "dist_func.txt";
        output_file.open(filename.c_str(),std::ios::trunc);
        for (PS::S32 n = 0; n < nbin; n++) {
            const PS::F64 z = (n + 0.5) * dz - safety * zmax;
            output_file << z << " " << dist_func[n] << std::endl;
        }
        output_file.close();
    }
#endif


}

