#pragma once

#include <sys/times.h>
#include <sys/time.h>
#include <string>
#include <mpi.h>
#include "run_param.hpp"
#include "cosmology.hpp"

#define TINY (1.0e-30)

template <typename T>
auto SQR(const T &x) {
    return x * x;
}

template <typename T>
auto CUBE(const T &x) {
    return x * x * x;
}

#define MAXPATHLEN (1024)

void make_directory(char *directory_name) {
    // Set the absolute PATH of the directory
    static char cwd_path[MAXPATHLEN];
    if (getcwd(cwd_path, sizeof(cwd_path)) == NULL) {
        fprintf(stderr, "Error in make_directory");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    strcat(cwd_path, "/");
    strcat(cwd_path, directory_name);

    // Make the directory at rank 0
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0) {
        mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
        mkdir(cwd_path, mode);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

class Map2D {
   private:
    PS::F64 *dens_loc;
    PS::F64 *dens_glb;
    PS::S32 n_mesh;
    PS::F64 l_mesh;
    PS::S32 snp_id;
    std::string file_name_prefix;
    std::string dir_name_prefix;

   public:
    template <class Tptcl>
    void initialize(const Tptcl &ptcl, const run_param &this_run, const PS::F64 res = 1.0) {
        // PS::S32 n_loc = ptcl.getNumberOfParticleLocal();
        PS::S32 n_tot = ptcl.getNumberOfParticleGlobal();
        n_mesh = cbrt((PS::F64)n_tot * 1.0000001) * res;
        l_mesh = 1.0 / n_mesh;
        if (PS::Comm::getRank() == 0) std::cerr << "n_mesh=" << n_mesh << " l_mesh=" << l_mesh << std::endl;
        dens_loc = new PS::F64[n_mesh * n_mesh];
        dens_glb = new PS::F64[n_mesh * n_mesh];
        snp_id = 0;
        dir_name_prefix = this_run.model_name;
        dir_name_prefix += "_";
        file_name_prefix = "map_2d";
        file_name_prefix += "_";
        if (PS::Comm::getRank() == 0) std::cerr << "dir_name_prefix=" << dir_name_prefix << " file_name_prefix=" << file_name_prefix << std::endl;
    }

    template <class Tptcl>
    void write(const Tptcl &ptcl, const run_param &this_run) {
        for (PS::S32 i = 0; i < n_mesh * n_mesh; i++) {
            dens_loc[i] = dens_glb[i] = 0.0;
        }
        PS::S32 n_loc = ptcl.getNumberOfParticleLocal();
        PS::F64 m_tot_loc = 0.0;
        for (PS::S32 i = 0; i < n_loc; i++) {
            PS::S32 idx = (PS::S32)(ptcl[i].pos.x / l_mesh);
            PS::S32 idy = (PS::S32)(ptcl[i].pos.y / l_mesh);
            dens_loc[idx * n_mesh + idy] += ptcl[i].mass;
            m_tot_loc += ptcl[i].mass;
        }
        MPI_Allreduce((PS::F64 *)&dens_loc[0], (PS::F64 *)&dens_glb[0], n_mesh * n_mesh, PS::GetDataType<PS::F64>(), MPI_SUM, MPI_COMM_WORLD);
        PS::F64 m_tot_glb = 0.0;
        MPI_Allreduce(&m_tot_loc, &m_tot_glb, 1, PS::GetDataType<PS::F64>(), MPI_SUM, MPI_COMM_WORLD);
        for (PS::S32 i = 0; i < n_mesh * n_mesh; i++) {
            dens_glb[i] /= (m_tot_glb * l_mesh * l_mesh);
        }
        if (PS::Comm::getRank() == 0) {
            std::ostringstream snp_id_str;
            snp_id_str << snp_id;
            std::ostringstream input;
            input << dir_name_prefix << snp_id_str.str() << "/" << file_name_prefix << snp_id_str.str();
            std::cerr << "input=" << input.str() << std::endl;
            std::ofstream fout;
            fout.open(input.str());
            for (PS::S32 ix = 0; ix < n_mesh; ix++) {
                fout << ix * l_mesh << "  " << "0.0" << "  " << dens_glb[ix * n_mesh] << std::endl;
                for (PS::S32 iy = 1; iy < n_mesh; iy++) {
                    fout << ix * l_mesh << "  " << iy * l_mesh << "  " << dens_glb[ix * n_mesh + iy] << std::endl;
                }
                fout << std::endl;
            }
            fout.close();
        }
        snp_id++;
    }
};

template <typename Tfp>
void output_data(PS::ParticleSystem<Tfp> &ptcl, run_param &this_run, char *filename) {
    FILE *output_fp = fopen(filename, "w");

    if (output_fp == NULL) {
        fprintf(stderr, "File %s cannot be written.", filename);
        exit(EXIT_FAILURE);
    }

    this_run.mpi_rank = PS::Comm::getRank();
    this_run.mpi_nproc = PS::Comm::getNumberOfProc();
    this_run.npart_local = ptcl.getNumberOfParticleLocal();
    this_run.npart_total = ptcl.getNumberOfParticleGlobal();
    this_run.write_header(output_fp);

    for (PS::S64 i = 0; i < this_run.npart_local; i++) {
        ptcl[i].writeParticleBinary(output_fp);
    }

    fclose(output_fp);
}

template <typename Tfp>
void output_data_in_run(PS::ParticleSystem<Tfp> &ptcl, run_param &this_run) {
    static char filename[256], directory_name[256];

    sprintf(directory_name, "%s_%d", this_run.model_name, this_run.output_indx);
    sprintf(filename, "%s/%s_%d-%d", directory_name, this_run.model_name, this_run.output_indx, this_run.mpi_rank);

    make_directory(directory_name);

    if (this_run.znow < this_run.output_timing[this_run.output_indx] + 0.001) {
        output_data(ptcl, this_run, filename);
        this_run.output_indx++;
    }
}

template <typename Tfp>
void output_data_in_run(PS::ParticleSystem<Tfp> &ptcl, run_param &this_run, Map2D &map) {
    // static char filename[1024], directory_name[1024];
    // sprintf(directory_name, "%s_%d", this_run.model_name, this_run.output_indx);
    // sprintf(filename, "%s/%s_%d-%d", directory_name, this_run.model_name, this_run.output_indx, this_run.mpi_rank);
    std::string d_name = std::string(this_run.model_name) + "_" + std::to_string(this_run.output_indx);
    std::string f_name = d_name + "/" + this_run.model_name + "_" + std::to_string(this_run.output_indx) + "-" + std::to_string(this_run.mpi_rank);
    char *directory_name = d_name.data();
    char *filename = f_name.data();
    make_directory(directory_name);
    if (this_run.znow < this_run.output_timing[this_run.output_indx] + 0.001) {
        output_data(ptcl, this_run, filename);
        map.write(ptcl, this_run);
        this_run.output_indx++;
    }
}

template <typename Tfp>
void drift_ptcl(PS::ParticleSystem<Tfp> &ptcl, PS::DomainInfo &dinfo, const PS::F64 dtime) {
    PS::S32 npart_local = ptcl.getNumberOfParticleLocal();
    for (PS::S64 i = 0; i < npart_local; i++) ptcl[i].pos += ptcl[i].vel * dtime;

    ptcl.adjustPositionIntoRootDomain(dinfo);
}
template <typename Tfp>
void reverse_ptcl_acc(PS::ParticleSystem<Tfp> &ptcl) {
    PS::S64 npart_local = ptcl.getNumberOfParticleLocal();

    for (PS::S64 i = 0; i < npart_local; i++) {
        ptcl[i].acc *= -1.0;
        ptcl[i].acc_pm *= -1.0;
    }
}

template <typename Tfp>
void kick_ptcl(PS::ParticleSystem<Tfp> &ptcl, const PS::F64 dtime, run_param &this_run) {
    PS::S64 npart_local = ptcl.getNumberOfParticleLocal();

    PS::F64 om = this_run.cosm.omegam;
    PS::F64 ov = this_run.cosm.omegav;

    PS::F64 anow = this_run.cosm.timetoa(this_run.tnow);

    PS::F64 at = sqrt(1.e0 + om * (1.e0 / anow - 1.e0) + ov * (SQR(anow) - 1.e0)) / anow;
    PS::F64 bt = 1.0 / CUBE(anow);

    PS::F64 atdt1 = 1.0 + at * dtime;
    PS::F64 vfact = (2.0 - atdt1) / atdt1;
    PS::F64 afact = bt * dtime / atdt1;

    for (PS::S64 i = 0; i < npart_local; i++) {
        // ptcl[i].vel = vfact * ptcl[i].vel + afact * (ptcl[i].force + ptcl[i].force_pmm);
        ptcl[i].vel = vfact * ptcl[i].vel + afact * ptcl[i].acc;
    }
}

template <typename Tfp>
PS::F64 calc_dtime(PS::ParticleSystem<Tfp> &ptcl, run_param &this_run) {
    PS::F64 dtime;

    dtime = DBL_MAX;
    for (PS::S64 i = 0; i < ptcl.getNumberOfParticleLocal(); i++) {
        dtime = fmin(dtime, ptcl[i].calcDtime(this_run));
    }

    if (this_run.noutput >= 1) {
        COSM::REAL zred_next;
        COSM::REAL zred_next_output;
        COSM::REAL time_next_output;

        if (this_run.znow < this_run.output_timing[this_run.output_indx] && this_run.output_indx + 1 < this_run.noutput) {
            zred_next_output = this_run.output_timing[this_run.output_indx + 1];
        } else {
            zred_next_output = this_run.output_timing[this_run.output_indx];
        }
        zred_next = this_run.cosm.timetoz(this_run.tnow + dtime);

        if (zred_next < zred_next_output) {
            time_next_output = this_run.cosm.ztotime(zred_next_output / 1.0001);
            dtime = time_next_output - this_run.tnow;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &dtime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    return dtime;
}

struct Result {
    PS::F64vec acc;
    PS::F64 pot;
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
    void addForcePMMM(const PS::F64vec &f, const PS::F64 p) {
        this->acc += f;
        this->pot -= p;
    }
};

class FP {
    template <class T>
    T reverseEndian(T value) {
        char *first = reinterpret_cast<char *>(&value);
        char *last = first + sizeof(T);
        std::reverse(first, last);
        return value;
    }

   public:
    PS::S64 id;
    PS::F64 mass;
    static inline PS::F64 eps = 0.0;
    static inline PS::F64 Lbnd = 0.0;
    static inline PS::F64 H0 = 0.0;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64vec getPos() const { return pos; }
    PS::F64 getCharge() const { return mass; }
    void copyFromForce(const Result &f) {
        this->acc = f.acc;
        this->pot = f.pot;
    }
    void copyFromForcePMMM(const Result &f) {
        this->acc = f.acc;
        this->pot = f.pot;
    }
    void setPos(const PS::F64vec &pos_new) { this->pos = pos_new; }
    void clear() {
        this->acc = 0.0;
        this->pot = 0.0;
    }
    PS::F64 calcDtime(run_param &this_run) {
        PS::F64 dtime_v, dtime_a, dtime;
        PS::F64 vnorm, anorm;
        vnorm = sqrt(SQR(this->vel)) + TINY;
        anorm = sqrt(SQR(this->acc)) + TINY;
        dtime_v = this->eps / vnorm;
        dtime_a = sqrt(this->eps / anorm) * CUBE(this_run.anow);
        dtime = fmin(0.5 * dtime_v, dtime_a);
        return dtime;
    }
    void writeParticleBinary(FILE *fp) {
        PS::F32 x = pos[0];
        PS::F32 y = pos[1];
        PS::F32 z = pos[2];
        PS::F32 vx = vel[0];
        PS::F32 vy = vel[1];
        PS::F32 vz = vel[2];
        PS::S32 i = id;
        PS::S32 m = mass;
        fwrite(&x, sizeof(PS::F32), 1, fp);
        fwrite(&vx, sizeof(PS::F32), 1, fp);
        fwrite(&y, sizeof(PS::F32), 1, fp);
        fwrite(&vy, sizeof(PS::F32), 1, fp);
        fwrite(&z, sizeof(PS::F32), 1, fp);
        fwrite(&vz, sizeof(PS::F32), 1, fp);
        fwrite(&m, sizeof(PS::F32), 1, fp);
        fwrite(&i, sizeof(PS::F32), 1, fp);
    }

    // for API of FDPS
    // in snapshot, L unit is Mpc/h, M unit is Msun, v unit is km/s
    void readBinary(FILE *fp) {
        static PS::S32 ONE = 1;
        static bool is_little_endian = *reinterpret_cast<char *>(&ONE) == ONE;
        static const PS::F64 Mpc_m = 3.08567e22;   // unit is m
        static const PS::F64 Mpc_km = 3.08567e19;  // unit is km
        static const PS::F64 Msun_kg = 1.9884e30;  // unit is kg
        static const PS::F64 G = 6.67428e-11;      // m^3*kg^-1*s^-2
        static const PS::F64 Cl = 1.0 / FP::Lbnd;
        static const PS::F64 Cv = 1.0 / (FP::Lbnd * FP::H0);
        static const PS::F64 Cm = 1.0 / (pow(Mpc_m * FP::Lbnd, 3.0) / pow(Mpc_km / FP::H0, 2.0) / G / Msun_kg);
        PS::F32 x, y, z, vx, vy, vz, m;
        PS::S32 i;

        auto ret = fread(&x, 4, 1, fp);
        (void)ret;
        ret = fread(&vx, 4, 1, fp);
        ret = fread(&y, 4, 1, fp);
        ret = fread(&vy, 4, 1, fp);
        ret = fread(&z, 4, 1, fp);
        ret = fread(&vz, 4, 1, fp);
        ret = fread(&m, 4, 1, fp);
        ret = fread(&i, 4, 1, fp);
        if (is_little_endian) {
            pos.x = x * Cl;
            pos.y = y * Cl;
            pos.z = z * Cl;
            vel.x = vx * Cv;
            vel.y = vy * Cv;
            vel.z = vz * Cv;
            mass = m * Cm;
            id = i;
        } else {
            pos.x = reverseEndian(x) * Cl;
            pos.y = reverseEndian(y) * Cl;
            pos.z = reverseEndian(z) * Cl;
            vel.x = reverseEndian(vx) * Cv;
            vel.y = reverseEndian(vy) * Cv;
            vel.z = reverseEndian(vz) * Cv;
            mass = reverseEndian(m) * Cm;
            id = reverseEndian(i);
        }
    }

    // for API of FDPS
    void writeBinary(FILE *fp) {
        static const PS::F64 Mpc_m = 3.08567e22;   // unit is m
        static const PS::F64 Mpc_km = 3.08567e19;  // unit is km
        static const PS::F64 Msun_kg = 1.9884e30;  // unit is kg
        static const PS::F64 G = 6.67428e-11;      // m^3*kg^-1*s^-2
        static const PS::F64 Cl = FP::Lbnd;
        static const PS::F64 Cv = (FP::Lbnd * FP::H0);
        static const PS::F64 Cm = (pow(Mpc_m * FP::Lbnd, 3.0) / pow(Mpc_km / FP::H0, 2.0) / G / Msun_kg);
        PS::F32vec x = pos * Cl;
        PS::F32vec v = vel * Cv;
        PS::F32 m = mass * Cm;
        PS::S32 i = id;
        fwrite(&x.x, sizeof(PS::F32), 1, fp);
        fwrite(&v.x, sizeof(PS::F32), 1, fp);
        fwrite(&x.y, sizeof(PS::F32), 1, fp);
        fwrite(&v.y, sizeof(PS::F32), 1, fp);
        fwrite(&x.z, sizeof(PS::F32), 1, fp);
        fwrite(&v.z, sizeof(PS::F32), 1, fp);
        fwrite(&m, sizeof(PS::F32), 1, fp);
        fwrite(&i, sizeof(PS::S32), 1, fp);
    }
    void writeAscii(FILE *fp) {
        fprintf(fp, "%ld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x,
                this->vel.y, this->vel.z, this->acc.x, this->acc.y, this->acc.z, this->pot);
    }
};

struct EPI {
    PS::S64 id;
    PS::F64vec pos;
    PS::F64vec getPos() const { return this->pos; }
    void copyFromFP(const FP &fp) {
        this->id = fp.id;
        this->pos = fp.pos;
    }
};

struct EPJ {
    PS::S64 id;
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64vec getPos() const { return this->pos; }
    void copyFromFP(const FP &fp) {
        this->id = fp.id;
        this->mass = fp.mass;
        this->pos = fp.pos;
    }
    PS::F64 getCharge() const { return this->mass; }
    void setPos(const PS::F64vec &pos_new) { this->pos = pos_new; }
};


#if defined(ENABLE_PIKG)
struct Epi {
    PS::F32vec pos;
};
struct Epj {
    PS::F32vec pos;
    PS::F32 mass;
};
struct Force {
    PS::F32vec acc;
    PS::F32 pot;
};

#include "kernel_pikg.hpp"

template <class TParticleJ>
struct CalcGravity {
    void operator()(const EPI *ep_i, const PS::S32 n_ip, const TParticleJ *ep_j, const PS::S32 n_jp, Result *force) {
        Epi epi[n_ip];
        Force f[n_ip];
        for (int i = 0; i < n_ip; i++) {
            epi[i].pos = (PS::F32vec)(ep_i[i].getPos() - ep_i[0].getPos());
            f[i].acc = force[i].acc;
            f[i].pot = force[i].pot;
        }
        Epj epj[n_jp];
        for (int i = 0; i < n_jp; i++) {
            epj[i].pos = (PS::F32vec)(ep_j[i].getPos() - ep_i[0].getPos());
            epj[i].mass = ep_j[i].getCharge();
        }
        CalcGravityEpEp(FP::eps * FP::eps)(epi, n_ip, epj, n_jp, f);
        for (int i = 0; i < n_ip; i++) {
            force[i].acc = f[i].acc;
            force[i].pot = f[i].pot;
        }
    }
};

#else  // ENABLE_PIKG

template <class TPJ>
struct CalcGravity {
    void operator()(EPI *ep_i, const PS::S32 n_ip, TPJ *ep_j, const PS::S32 n_jp, Result *result) {
        PS::F64 eps2 = FP::eps * FP::eps;
        for (PS::S32 i = 0; i < n_ip; i++) {
            PS::F64vec xi = ep_i[i].getPos();
            PS::F64vec Mgrad = 0.0;
            PS::F64 poti = 0.0;
            for (PS::S32 j = 0; j < n_jp; j++) {
                PS::F64vec rij = xi - ep_j[j].getPos();
                if (rij.x == 0.0 && rij.y == 0.0 && rij.z == 0.0) continue;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0 / sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].getCharge();
                r3_inv *= r_inv;
                Mgrad += r3_inv * rij;
                poti += r_inv;
            }
            result[i].acc -= Mgrad;
            result[i].pot -= poti;
        }
    }
};

#endif  // ENABLE_PIKG
