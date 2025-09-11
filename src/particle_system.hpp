#pragma once

#include <cassert>
#include <fstream>

#include "MT.hpp"
#include "ps_defs.hpp"

namespace ParticleSimulator {
template <class Tptcl>
class ParticleSystem {
   private:
    CountT n_ptcl_send_;
    CountT n_ptcl_recv_;
    TimeProfile time_profile_;
    const S32 n_smp_ave_ = 30;
    ReallocatableArray<Tptcl> ptcl_;
    ReallocatableArray<S32> idx_remove_ptcl_;  // add 2016/09/14
    S32 n_smp_ptcl_tot_;
    bool first_call_by_initialize;
    bool first_call_by_setAverageTargetNumberOfSampleParticlePerProcess;
    bool first_call_by_DomainInfo_collect_sample_particle;
    CountT n_proc_exch_;
    PS_INLINE bool determineWhetherParticleIsInDomain(const F64vec& pos, const F64ort& domain) {
        bool ret = true;
        ret = ret && (domain.low_.x <= pos.x) && (pos.x < domain.high_.x);
        ret = ret && (domain.low_.y <= pos.y) && (pos.y < domain.high_.y);
#if !defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        ret = ret && (domain.low_.z <= pos.z) && (pos.z < domain.high_.z);
#endif
        return ret;
    }

    PS_INLINE S32 searchWhichDomainParticleGoTo(const F64vec& pos, const S32 n_domain[DIMENSION_LIMIT], const F64ort domain[], const S32 rank_start,
                                                S32vec& dnp) {
        if (domain[rank_start].contained(pos)) {
            return rank_start;
        }
        const S32vec rank_vec_start = GetRankVec(n_domain, rank_start);
        S32vec rank_vec = rank_vec_start;
        S32 rank_new = rank_start;

        rank_vec.x = GetAppropriateRank1D(pos, n_domain, domain, 0, rank_new);
        // even if open boundary condition, this function (GetDistanceMin1DPeriImpl) works correctly but, the cost is rather high.
        dnp.x = GetDistanceMin1DPeriImpl(rank_vec.x, rank_vec_start.x, n_domain[0]);  // logical coordinate distance in x-direction
        rank_new = GetRankGlb(rank_vec, n_domain);

        rank_vec.y = GetAppropriateRank1D(pos, n_domain, domain, 1, rank_new);
        dnp.y = GetDistanceMin1DPeriImpl(rank_vec.y, rank_vec_start.y, n_domain[1]);
        rank_new = GetRankGlb(rank_vec, n_domain);

#if !defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        rank_vec.z = GetAppropriateRank1D(pos, n_domain, domain, 2, rank_new);
        dnp.z = GetDistanceMin1DPeriImpl(rank_vec.z, rank_vec_start.z, n_domain[2]);
        rank_new = GetRankGlb(rank_vec, n_domain);
#endif
        return rank_new;
    }

   public:
    CommInfo comm_info_;
    void setCommInfo(const CommInfo& c) { comm_info_ = c; }
#ifdef TEST_VARIADIC_TEMPLATE
    void ParticleSystemDummyFunc() {}
#endif

    TimeProfile getTimeProfile() const { return time_profile_; }
    void clearTimeProfile() { time_profile_.clear(); }
    void dumpTimeProfile(std::ostream& fout = std::cout, const S32 level = 0) const { time_profile_.dump_particle_system(fout, level); }

    ParticleSystem() {
        first_call_by_setAverageTargetNumberOfSampleParticlePerProcess = true;
        first_call_by_initialize = true;
        first_call_by_DomainInfo_collect_sample_particle = true;
        comm_info_.setCommunicator();
    }

    void initialize(const S32 cap = 10000) {
        assert(first_call_by_initialize);
        first_call_by_initialize = false;
        // n_smp_ptcl_tot_ = n_smp_ave_ * Comm::getNumberOfProc();
        // std::cerr << "n_smp_ave_= " << n_smp_ave_ << std::endl;
        // std::cerr << "comm_info_.getNumberOfProc()= " << comm_info_.getNumberOfProc() << std::endl;
        n_smp_ptcl_tot_ = n_smp_ave_ * comm_info_.getNumberOfProc();
        n_ptcl_send_ = n_ptcl_recv_ = 0;
        // MT::init_genrand(Comm::getRank()); // 2017.11.01
        MT::init_genrand(comm_info_.getRank());
        createParticle(cap);
    }

    void setAverageTargetNumberOfSampleParticlePerProcess(const S32& nsampleperprocess) {
        assert(first_call_by_setAverageTargetNumberOfSampleParticlePerProcess);
        first_call_by_setAverageTargetNumberOfSampleParticlePerProcess = false;
        // n_smp_ptcl_tot_ = nsampleperprocess * Comm::getNumberOfProc();
        n_smp_ptcl_tot_ = nsampleperprocess * comm_info_.getNumberOfProc();
    }

    S32 getTargetNumberOfSampleParticle() { return n_smp_ptcl_tot_; }

    bool getFirstCallByDomainInfoCollectSampleParticle() {
        if (first_call_by_DomainInfo_collect_sample_particle) {
            first_call_by_DomainInfo_collect_sample_particle = false;
            return true;
        } else {
            return false;
        }
    }

    void createParticle(const S32 n_limit, bool clear = true) {
        ptcl_.reserve(n_limit);
        ptcl_.resizeNoInitialize(0);
    }

    // void setNumberOfParticleLocal(const S32 n){ n_ptcl_ = n; }
    void setNumberOfParticleLocal(const S32 n) {
        // 15/02/20 Hosono bug(?) fix
        ptcl_.resizeNoInitialize(n);
        // ptcl_.dump();
    }
    ////////////////
    // 05/01/30 Hosono From
    ////////////////
    // dummy class for the case if User does NOT define the file header.
    // TO BE PRIVATE
    struct DummyHeader {
        void writeAscii(FILE* fp) const {}
        S32 readAscii(FILE* fp) { return -1; }
        void writeBinary(FILE* fp) const {}
        S32 readBinary(FILE* fp) { return -1; }
    };

    // fp must points to the head address of the first partticle
    S32 countPtclAscii(FILE* fp) {
        S32 n = 0;
        for (int c; (c = getc(fp)) != EOF; n += '\n' == c ? 1 : 0) {
        }
        return n;
    }
    S32 countPtclBinary(FILE* fp, void (Tptcl::*pFuncPtcl)(FILE*)) {
        S32 n = 0;
        long end_of_header = ftell(fp);
        Tptcl ptcl_tmp;
        (ptcl_tmp.*pFuncPtcl)(fp);
        long size_ptcl = ftell(fp) - end_of_header;
        if (size_ptcl > 0) {
            fseek(fp, 0, SEEK_END);
            n = (ftell(fp) - end_of_header) / size_ptcl;
        }
        return n;
    }
    // 2016_11_03 modified IO functions to handle multiple file formats
    // read
    template <typename Theader, typename TfuncPtcl, typename TfuncHeader>
    void readParticleImpl(const char* const filename, const char* const format, Theader* header, TfuncPtcl pFuncPtcl, TfuncHeader pFuncHead,
                          const char* const open_format) {
        if (format == NULL) {  // Read from single file
            // if(comm_info_.getRank() == 0){
            FILE* fp = fopen(filename, open_format);
            if (fp == NULL) {
                PARTICLE_SIMULATOR_PRINT_ERROR("can not open input file ");
                std::cerr << "filename: " << filename << std::endl;
                Abort(-1);
            }
            S32 n_ptcl_ = (header->*pFuncHead)(fp);
            if (n_ptcl_ < 0) {  // User does NOT return # of ptcl
                if (strcmp(open_format, "rb") == 0) {
                    n_ptcl_ = countPtclBinary(fp, pFuncPtcl);
                } else {
                    // count # of lines
                    // KN
                    n_ptcl_ = countPtclAscii(fp);
                }
                fseek(fp, 0, SEEK_SET);
                (header->*pFuncHead)(fp);
            }
            if (comm_info_.getRank() == 0) {
                // Inform the # of ptcl for each process.
                const S32 n_proc = comm_info_.getNumberOfProc();
                S32* n_ptcl = new S32[n_proc];
                for (S32 i = 0; i < n_proc; ++i) {
                    n_ptcl[i] = n_ptcl_ / n_proc;
                }
                n_ptcl[0] += n_ptcl_ % n_proc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                comm_info_.scatter(n_ptcl, 1, &n_ptcl_);
#endif
                // allocate ptcl.
                // First of all, Rank 0 reads its own particle.
                this->createParticle(n_ptcl_ << 2);  // Magic shift
                ptcl_.resizeNoInitialize(n_ptcl_);
                for (int i = 0; i < n_ptcl_; ++i) {
                    (ptcl_[i].*pFuncPtcl)(fp);
                }
// Read remaining data to buffer and send them to appropriate process.
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                for (S32 rank = 1; rank < n_proc; ++rank) {
                    Tptcl* buffer = new Tptcl[n_ptcl[rank]];
                    for (int i = 0; i < n_ptcl[rank]; ++i) {
                        (buffer[i].*pFuncPtcl)(fp);
                    }
                    MPI_Send(buffer, n_ptcl[rank], GetDataType<Tptcl>(), rank, 0, comm_info_.getCommunicator());
                    delete[] buffer;
                }
#endif
                // End.
                delete[] n_ptcl;
                fclose(fp);
            } else {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                // Receive the # of ptcl from Rank 0
                S32 n_ptcl_loc;
                S32* n_ptcl = new S32[comm_info_.getNumberOfProc()];
                comm_info_.scatter(n_ptcl, 1, &n_ptcl_loc);
                delete[] n_ptcl;
                // allocate ptcl.
                this->createParticle(n_ptcl_loc << 2);  // Magic shift
                ptcl_.resizeNoInitialize(n_ptcl_loc);
                MPI_Status stat;
                MPI_Recv(ptcl_.getPointer(), ptcl_.size(), GetDataType<Tptcl>(), 0, 0, comm_info_.getCommunicator(), &stat);
#endif
            }
        } else {  // Read from multiple file
            char input[256];
            sprintf(input, format, filename, comm_info_.getNumberOfProc(), comm_info_.getRank());
            FILE* fp = fopen(input, open_format);
            if (fp == NULL) {
                PARTICLE_SIMULATOR_PRINT_ERROR("can not open input file");
                std::cerr << "filename: " << input << std::endl;
                Abort(-1);
            }
            S32 n_ptcl_ = (header->*pFuncHead)(fp);
            if (n_ptcl_ < 0) {  // User does NOT return # of ptcl
                if (strcmp(open_format, "rb") == 0) {
                    n_ptcl_ = countPtclBinary(fp, pFuncPtcl);
                } else {
                    n_ptcl_ = countPtclAscii(fp);
                }
                this->createParticle(n_ptcl_ << 2);  // Magic shift
                ptcl_.resizeNoInitialize(n_ptcl_);
                fseek(fp, 0, SEEK_SET);
                // S32 tmp = (header->*pFuncHead)(fp);
                (header->*pFuncHead)(fp);
                for (S32 i = 0; i < ptcl_.size(); ++i) {
                    (ptcl_[i].*pFuncPtcl)(fp);
                }
                fclose(fp);
            } else {
                // User returns # of ptcl.
                this->createParticle(n_ptcl_ << 2);  // Magic shift
                ptcl_.resizeNoInitialize(n_ptcl_);
                for (S32 i = 0; i < n_ptcl_; ++i) {
                    (ptcl_[i].*pFuncPtcl)(fp);
                }
                fclose(fp);
            }
        }
    }

    template <typename Tchar0, typename Tchar1, typename Theader, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr>
    void readParticleAscii(Tchar0 filename, Tchar1 format, Theader& header) {
        readParticleImpl(filename, format, &header, &Tptcl::readAscii, &Theader::readAscii, "r");
    }
    template <typename Tchar0, typename Theader, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr>
    void readParticleAscii(Tchar0 filename, Theader& header) {
        readParticleImpl(filename, NULL, &header, &Tptcl::readAscii, &Theader::readAscii, "r");
    }
    template <typename Tchar0, typename Tchar1, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr>
    void readParticleAscii(Tchar0 filename, Tchar1 format) {
        readParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::readAscii, &DummyHeader::readAscii, "r");
    }
    template <typename Tchar0, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr>
    void readParticleAscii(Tchar0 filename) {
        readParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::readAscii, &DummyHeader::readAscii, "r");
    }

    template <typename Tchar0, typename Tchar1, typename Theader, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void readParticleAscii(Tchar0 filename, Tchar1 format, Theader& header, Tfunc pFunc) {
        readParticleImpl(filename, format, &header, pFunc, &Theader::readAscii, "r");
    }
    template <typename Tchar0, typename Theader, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void readParticleAscii(Tchar0 filename, Theader& header, Tfunc pFunc) {
        readParticleImpl(filename, NULL, &header, pFunc, &Theader::readAscii, "r");
    }
    template <typename Tchar0, typename Tchar1, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void readParticleAscii(Tchar0 filename, Tchar1 format, Tfunc pFunc) {
        readParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::readAscii, "r");
    }
    template <typename Tchar0, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void readParticleAscii(Tchar0 filename, Tfunc pFunc) {
        readParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::readAscii, "r");
    }

    template <typename Tchar0, typename Tchar1, typename Theader, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr>
    void readParticleBinary(Tchar0 filename, Tchar1 format, Theader& header) {
        readParticleImpl(filename, format, &header, &Tptcl::readBinary, &Theader::readBinary, "rb");
    }
    template <typename Tchar0, typename Theader, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr>
    void readParticleBinary(Tchar0 filename, Theader& header) {
        readParticleImpl(filename, NULL, &header, &Tptcl::readBinary, &Theader::readBinary, "rb");
    }
    template <typename Tchar0, typename Tchar1, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr>
    void readParticleBinary(Tchar0 filename, Tchar1 format) {
        readParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::readBinary, &DummyHeader::readBinary, "rb");
    }
    template <typename Tchar0, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr>
    void readParticleBinary(Tchar0 filename) {
        readParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::readBinary, &DummyHeader::readBinary, "rb");
    }

    template <typename Tchar0, typename Tchar1, typename Theader, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void readParticleBinary(Tchar0 filename, Tchar1 format, Theader& header, Tfunc pFunc) {
        readParticleImpl(filename, format, &header, pFunc, &Theader::readBinary, "rb");
    }
    template <typename Tchar0, typename Theader, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void readParticleBinary(Tchar0 filename, Theader& header, Tfunc pFunc) {
        readParticleImpl(filename, NULL, &header, pFunc, &Theader::readBinary, "rb");
    }
    template <typename Tchar0, typename Tchar1, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void readParticleBinary(Tchar0 filename, Tchar1 format, Tfunc pFunc) {
        readParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::readBinary, "rb");
    }
    template <typename Tchar0, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void readParticleBinary(Tchar0 filename, Tfunc pFunc) {
        readParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::readBinary, "rb");
    }

    template <typename Theader, typename TfuncPtcl, typename TfuncHeader>
    void writeParticleImpl(const char* const filename, const char* const format, Theader* header, TfuncPtcl pFuncPtcl, TfuncHeader pFuncHead,
                           const char* const open_format) const {
        if (format == NULL) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            // declare local # of ptcl.
            const S32 n_ptcl_ = ptcl_.size();
            // get # of process.
            const S32 n_proc = comm_info_.getNumberOfProc();
            // get # of ptcls in each process.
            S32* n_ptcl = new S32[n_proc];
            // Gather # of particles.
            comm_info_.allGather(const_cast<S32*>(&n_ptcl_), 1, n_ptcl);
            // set displacement
            S32* n_ptcl_displs = new S32[n_proc + 1];
            n_ptcl_displs[0] = 0;
            for (S32 i = 0; i < n_proc; ++i) {
                n_ptcl_displs[i + 1] = n_ptcl_displs[i] + n_ptcl[i];
            }
            const S32 n_tot = n_ptcl_displs[n_proc];
            Tptcl* ptcl = new Tptcl[n_tot];
            // gather data
            comm_info_.gatherV(ptcl_.getPointer(), n_ptcl_, ptcl, n_ptcl, n_ptcl_displs);
            if (comm_info_.getRank() == 0) {
                FILE* fp = fopen(filename, open_format);
                if (fp == NULL) {
                    PARTICLE_SIMULATOR_PRINT_ERROR("can not open output file");
                    std::cerr << "output file: " << filename << std::endl;
                    Abort(-1);
                }
                (header->*pFuncHead)(fp);
                for (S32 i = 0; i < n_tot; ++i) {
                    (ptcl[i].*pFuncPtcl)(fp);
                }
                fclose(fp);
            }
            delete[] n_ptcl;
            delete[] n_ptcl_displs;
            delete[] ptcl;
#else
            const S32 n_tot = ptcl_.size();
            if (comm_info_.getRank() == 0) {
                FILE* fp = fopen(filename, open_format);
                (header->*pFuncHead)(fp);
                for (S32 i = 0; i < n_tot; ++i) {
                    (ptcl_[i].*pFuncPtcl)(fp);
                }
                fclose(fp);
            }
#endif
        } else {
            char output[256];
            sprintf(output, format, filename, comm_info_.getNumberOfProc(), comm_info_.getRank());
            FILE* fp = fopen(output, open_format);
            if (fp == NULL) {
                PARTICLE_SIMULATOR_PRINT_ERROR("can not open output file");
                std::cerr << "output file: " << output << std::endl;
                Abort(-1);
            }
            (header->*pFuncHead)(fp);
            for (S32 i = 0; i < ptcl_.size(); ++i) {
                (ptcl_[i].*pFuncPtcl)(fp);
            }
            fclose(fp);
        }
    }

    // write
    template <typename Tchar0, typename Tchar1, typename Theader, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr>
    void writeParticleAscii(Tchar0 filename, Tchar1 format, Theader& header) const {
        writeParticleImpl(const_cast<const char*>(filename), const_cast<const char*>(format), &header, &Tptcl::writeAscii, &Theader::writeAscii, "w");
    }
    template <typename Tchar0, typename Theader, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr>
    void writeParticleAscii(Tchar0 filename, Theader& header) const {
        writeParticleImpl(filename, NULL, &header, &Tptcl::writeAscii, &Theader::writeAscii, "w");
    }
    template <typename Tchar0, typename Tchar1, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr>
    void writeParticleAscii(Tchar0 filename, Tchar1 format) const {
        writeParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::writeAscii, &DummyHeader::writeAscii, "w");
    }
    template <typename Tchar0, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr>
    void writeParticleAscii(Tchar0 filename) const {
        writeParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::writeAscii, &DummyHeader::writeAscii, "w");
    }

    template <typename Tchar0, typename Tchar1, typename Theader, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void writeParticleAscii(Tchar0 filename, Tchar1 format, Theader& header, Tfunc pFunc) const {
        writeParticleImpl(filename, format, &header, pFunc, &Theader::writeAscii, "w");
    }
    template <typename Tchar0, typename Theader, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void writeParticleAscii(Tchar0 filename, Theader& header, Tfunc pFunc) const {
        writeParticleImpl(filename, NULL, &header, pFunc, &Theader::writeAscii, "w");
    }
    template <typename Tchar0, typename Tchar1, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void writeParticleAscii(Tchar0 filename, Tchar1 format, Tfunc pFunc) const {
        writeParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::writeAscii, "w");
    }
    template <typename Tchar0, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void writeParticleAscii(Tchar0 filename, Tfunc pFunc) const {
        writeParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::writeAscii, "w");
    }

    template <typename Tchar0, typename Tchar1, typename Theader, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr>
    void writeParticleBinary(Tchar0 filename, Tchar1 format, Theader& header) const {
        writeParticleImpl(filename, format, &header, &Tptcl::writeBinary, &Theader::writeBinary, "wb");
    }
    template <typename Tchar0, typename Theader, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr>
    void writeParticleBinary(Tchar0 filename, Theader& header) const {
        writeParticleImpl(filename, NULL, &header, &Tptcl::writeBinary, &Theader::writeBinary, "wb");
    }
    template <typename Tchar0, typename Tchar1, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr>
    void writeParticleBinary(Tchar0 filename, Tchar1 format) const {
        writeParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
    }
    template <typename Tchar0, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr>
    void writeParticleBinary(Tchar0 filename) const {
        writeParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
    }

    template <typename Tchar0, typename Tchar1, typename Theader, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void writeParticleBinary(Tchar0 filename, Tchar1 format, Theader& header, Tfunc pFunc) const {
        writeParticleImpl(filename, format, &header, pFunc, &Theader::writeBinary, "wb");
    }
    template <typename Tchar0, typename Theader, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_class<Theader>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void writeParticleBinary(Tchar0 filename, Theader& header, Tfunc pFunc) const {
        writeParticleImpl(filename, NULL, &header, pFunc, &Theader::writeBinary, "wb");
    }
    template <typename Tchar0, typename Tchar1, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<my_is_char<Tchar1>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void writeParticleBinary(Tchar0 filename, Tchar1 format, Tfunc pFunc) const {
        writeParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::writeBinary, "wb");
    }
    template <typename Tchar0, typename Tfunc, typename std::enable_if<my_is_char<Tchar0>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<Tfunc>::value>::type* = nullptr>
    void writeParticleBinary(Tchar0 filename, Tfunc pFunc) const {
        writeParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::writeBinary, "wb");
    }

    Tptcl& operator[](const S32 id) { return ptcl_[id]; }
    const Tptcl& operator[](const S32 id) const { return ptcl_[id]; }
    Tptcl& getParticle(const S32 id = 0) { return ptcl_[id]; }
    const Tptcl& getParticle(const S32 id = 0) const { return ptcl_[id]; }
    Tptcl* getParticlePointer(const S32 id = 0) const { return ptcl_ + id; }
    // S32 getNumberOfParticleLocal() const {return n_ptcl_;}
    S32 getNumberOfParticleLocal() const { return ptcl_.size(); }
    ////////////////
    // 05/02/04 Hosono From
    // 11th jul MI modified
    ////////////////
    S64 getNumberOfParticleGlobal() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        const S64 n_ptcl_loc = ptcl_.size();
        // return Comm::getSum(n_ptcl_loc);
        return comm_info_.getSum(n_ptcl_loc);
#else
        return ptcl_.size();
#endif
    }
    ////////////////
    // 05/02/04 Hosono To
    ////////////////

    F64 getHalfLength(const F64vec& center = 0.0) {
        F64 hl_max_loc = (ptcl_[0].getPos() - center).applyEach(Abs<F64>()).getMax();
        const S32 n_ptcl = ptcl_.size();
        for (size_t i = 1; i < n_ptcl; i++) {
            F64 hl_tmp = (ptcl_[i].getPos() - center).applyEach(Abs<F64>()).getMax();
            hl_max_loc = (hl_max_loc > hl_tmp) ? hl_max_loc : hl_tmp;
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // const F64 hl_max_glb = Comm::getMaxValue(hl_max_loc);
        const F64 hl_max_glb = comm_info_.getMaxValue(hl_max_loc);
        return hl_max_glb;
#else
        return hl_max_loc;
#endif
    }

    // *******************************************************************
    // ************** This can be replaced with MT method. ***************
    // *******************************************************************
    inline F64 frand() {
        //                return (double) rand() / ((double) RAND_MAX + 1.);
        // return genrand_res53();
        return MT::genrand_res53();
    }
    // *******************************************************************
    // *******************************************************************

    inline S32 getUniformDistributionFromArg1ToArg2(S32 arg1, S32 arg2) {
        S32 random_number;

        random_number = (S32)((arg2 - arg1 + 1) * frand()) + arg1;

        return random_number;
    }

    /* AT_DEBU
    void getSampleParticle(S32 & number_of_sample_particle,
                           F32vec pos_sample[],
                           const F32 weight=1.0) {
    */

    S32 getNumberOfSampleParticleLocal(const F64 weight) {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        S64 nloc = (S64)ptcl_.size();
        // const auto weight_all = Comm::getSum(weight);
        const auto weight_all = comm_info_.getSum(weight);
        S32 number_of_sample_particle = (weight * n_smp_ptcl_tot_) / weight_all;
#if 0
            // modified to limit # of sample particles by M.I. 
            const F32 coef_limitter = 0.2;
            //S64 nglb = Comm::getSum( (S64)nloc );
            S64 nglb = comm_info_.getSum( (S64)nloc );
            S32 number_of_sample_particle_limit = ((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter)); // lower limit
            number_of_sample_particle = (number_of_sample_particle > number_of_sample_particle_limit) ? number_of_sample_particle : number_of_sample_particle_limit;
#endif
        number_of_sample_particle = (number_of_sample_particle < nloc) ? number_of_sample_particle : nloc;

#if defined(PARTICLE_SIMULATOR_DEBUG_PRINT)
        std::cout << "weight=" << weight << " weight_all=" << weight_all << std::endl;
        std::cout << "n_smp_ptcl_tot_=" << n_smp_ptcl_tot_ << std::endl;
        std::cout << "(weight * n_smp_ptcl_tot_) / weight_all=" << (weight * n_smp_ptcl_tot_) / weight_all << std::endl;
        std::cout << "((S64)nloc * n_smp_ptcl_tot_)=" << ((S64)nloc * n_smp_ptcl_tot_) << std::endl;
        std::cout << "number_of_sample_particle(final)=" << number_of_sample_particle << std::endl;
#endif
        return number_of_sample_particle;
#else
        return 0;
#endif
    }

    void getSampleParticle2(const S32& number_of_sample_particle, ReallocatableArray<F64vec>& pos_sample, const S32 offset_pos_sample) {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        RandomSampling(ptcl_, pos_sample, [](const Tptcl& src, F64vec& dst) { dst = src.getPos(); }, 0, offset_pos_sample);
        return;
#endif
    }

    /*
    void getSampleParticle(S32 & number_of_sample_particle,
                           F64vec pos_sample[],
                           const F64 weight) {

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        S64 nloc = (S64)ptcl_.size();
        //F64 weight_all = Comm::getSum(weight);
        F64 weight_all = comm_info_.getSum(weight);
        number_of_sample_particle = (weight * n_smp_ptcl_tot_) / weight_all;
        //std::cerr<<"number_of_sample_particle= "<<number_of_sample_particle<<std::endl;
#if 0
// modified to limit # of sample particles by M.I.
        const F32 coef_limitter = 0.2;
        //S64 nglb = Comm::getSum( (S64)nloc );
        S64 nglb = comm_info_.getSum( (S64)nloc );
        S32 number_of_sample_particle_limit = ((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter)); // lower limit
        number_of_sample_particle = (number_of_sample_particle > number_of_sample_particle_limit) ? number_of_sample_particle :
number_of_sample_particle_limit; #endif number_of_sample_particle = (number_of_sample_particle < nloc) ? number_of_sample_particle : nloc;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"weight="<<weight<<" weight_all="<<weight_all<<std::endl;
        std::cout<<"n_smp_ptcl_tot_="<<n_smp_ptcl_tot_<<std::endl;
        std::cout<<"(weight * n_smp_ptcl_tot_) / weight_all="<<(weight * n_smp_ptcl_tot_) / weight_all<<std::endl;
        std::cout<<"((S64)nloc * n_smp_ptcl_tot_)="<< ((S64)nloc * n_smp_ptcl_tot_)<<std::endl;
        std::cout<<"number_of_sample_particle(final)="<<number_of_sample_particle<<std::endl;
#endif
        S32 *record = new S32[number_of_sample_particle];
        for(S32 i = 0; i < number_of_sample_particle; i++) {
            S32 j = getUniformDistributionFromArg1ToArg2(i, nloc-1);
            Tptcl hold = ptcl_[j];
            ptcl_[j]   = ptcl_[i];
            ptcl_[i]   = hold;
            record[i]  = j;
        }
        for(S32 i = 0; i < number_of_sample_particle; i++) {
            pos_sample[i] = ptcl_[i].getPos();
        }
        for(S32 i = number_of_sample_particle - 1; i >= 0; i--) {
            S32 j = record[i];
            Tptcl hold = ptcl_[j];
            ptcl_[j]   = ptcl_[i];
            ptcl_[i]   = hold;
        }
        delete [] record;
        return;
#endif
    }
    */

    template <class Tdinfo>
    void exchangeParticle(Tdinfo& dinfo, const bool flag_serialize = false) {
        // std::cerr<<"CHECK -a"<<std::endl;
        if (flag_serialize == true) {
            PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
            Abort(-1);
        }
        //	    Comm::barrier(); //for time measurement only; should not be here.
        F64 time_offset = GetWtime();
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
        const S32 n_loc = ptcl_.size();
        const S32 my_rank = comm_info_.getRank();
        const S32 n_proc = comm_info_.getNumberOfProc();
        const S32* n_domain = dinfo.getPointerOfNDomain();
        const F64ort* pos_domain = dinfo.getPointerOfPosDomain();
        S32 n_send[n_proc], n_disp_send[n_proc + 1], n_recv[n_proc], n_disp_recv[n_proc + 1];
        for (S32 i = 0; i < n_proc; i++) {
            n_send[i] = n_disp_send[i] = n_recv[i] = n_disp_recv[i] = 0;
        }
        n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
        MPI_Request req_recv[n_proc];
        MPI_Status status[n_proc];
        S32 rank_send[n_loc];
        ReallocatableArray<Tptcl> ptcl_send(0, 0, MemoryAllocMode::Pool);
        //F64 time_offset_inner = GetWtime();
        // F64ort pos_root_domain = dinfo.getPosRootDomain();
        // F64vec len_peri = pos_root_domain.high_ - pos_root_domain.low_;
        bool pa[DIMENSION_LIMIT];
        dinfo.getPeriodicAxis(pa);
        S32vec dnp_max_loc(0);
#if 1
        //  parallel version
        const S32 n_thread = Comm::getNumberOfThread();
        // std::vector<ReallocatableArray<S32>> n_send_per_thread(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Stack));
        // std::vector<ReallocatableArray<S32>> n_disp_send_per_thread(n_thread+1, ReallocatableArray<S32>(MemoryAllocMode::Stack));
        ReallocatableArray<S32> n_send_per_thread[n_thread];
        ReallocatableArray<S32> n_disp_send_per_thread[n_thread + 1];

        ReallocatableArray<S32> adr_ptcl_send(n_loc, 0, MemoryAllocMode::Pool);
        ReallocatableArray<S32> adr_ptcl_remain(n_loc, 0, MemoryAllocMode::Pool);
        S32vec dnp_max_loc_per_thread[n_thread];

        // std::vector<ReallocatableArray<S32>> adr_ptcl_send_per_thread(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Stack));
        // std::vector<ReallocatableArray<S32>> adr_ptcl_remain_per_thread(n_thread, ReallocatableArray<S32>(MemoryAllocMode::Stack));
        ReallocatableArray<S32> adr_ptcl_send_per_thread[n_thread];
        ReallocatableArray<S32> adr_ptcl_remain_per_thread[n_thread];
        for (S32 i = 0; i < n_thread; i++) {
            dnp_max_loc_per_thread[i] = S32vec(0);
            n_send_per_thread[i].resizeNoInitialize(n_proc);
            adr_ptcl_send_per_thread[i].reserve(n_loc / n_thread + 10);
            adr_ptcl_remain_per_thread[i].reserve(n_loc / n_thread + 10);
            n_disp_send_per_thread[i].resizeNoInitialize(n_proc);
        }
        n_disp_send_per_thread[n_thread].resizeNoInitialize(n_proc);
        auto t0 = GetWtime();
        PS_OMP_PARALLEL {
            const auto ith = Comm::getThreadNum();
            adr_ptcl_send_per_thread[ith].clearSize();
            adr_ptcl_remain_per_thread[ith].clearSize();
            for (auto i = 0; i < n_proc; i++) {
                n_send_per_thread[ith][i] = 0;
                n_disp_send_per_thread[ith][i] = 0;
            }
            S32 ip_head, ip_end;
            CalcAdrToSplitData(ip_head, ip_end, ith, n_thread, n_loc);
            S32vec dnp_max_tmp = 0;
            for (S32 ip = ip_head; ip < ip_end; ip++) {
                S32vec dnp(0);
                S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain, comm_info_.getRank(), dnp);
                rank_send[ip] = srank;
                n_send_per_thread[ith][srank]++;
                if (srank == my_rank) {
                    adr_ptcl_remain_per_thread[ith].push_back(ip);
                    continue;
                } else {
                    adr_ptcl_send_per_thread[ith].push_back(ip);
                }
                dnp_max_tmp = ApplyEach([](const F64& l, const F64& r) { return std::max(l, r); }, dnp_max_tmp, dnp);
            }
            dnp_max_loc_per_thread[ith] = dnp_max_tmp;
            n_send_per_thread[ith][my_rank] = 0;
            PS_OMP_BARRIER
            PackDataInOmp(&adr_ptcl_send_per_thread[0], adr_ptcl_send, n_thread, ith, 0);
            PackDataInOmp(&adr_ptcl_remain_per_thread[0], adr_ptcl_remain, n_thread, ith, 0);
            PS_OMP_FOR
            for (auto i = 0; i < n_proc; i++) {
                n_send[i] = 0;
                n_disp_send_per_thread[0][i] = 0;
                for (auto j = 0; j < n_thread; j++) {
                    n_send[i] += n_send_per_thread[j][i];
                    n_disp_send_per_thread[j + 1][i] = n_disp_send_per_thread[j][i] + n_send_per_thread[j][i];
                    n_send_per_thread[j][i] = 0;  // set zero to reuse this buffer
                }
            }
            PS_OMP_BARRIER
            PS_OMP_SINGLE {
                dnp_max_loc = S32vec(0);
                for (auto i = 0; i < n_thread; i++) {
                    dnp_max_loc = ApplyEach([](const S32& l, const S32& r) { return std::max(l, r); }, dnp_max_loc, dnp_max_loc_per_thread[i]);
                }
                n_disp_send[0] = 0;
                for (S32 i = 0; i < n_proc; i++) {
                    n_disp_send[i + 1] = n_disp_send[i] + n_send[i];
                }
                ptcl_send.resizeNoInitialize(n_disp_send[n_proc]);
            }
            PS_OMP_BARRIER
            for (auto ip = 0; ip < adr_ptcl_send_per_thread[ith].size(); ip++) {
                const auto adr = adr_ptcl_send_per_thread[ith][ip];
                const auto srank = rank_send[adr];
                const auto adr_dst = n_disp_send[srank] + n_disp_send_per_thread[ith][srank] + n_send_per_thread[ith][srank];
                ptcl_send[adr_dst] = ptcl_[adr];
                n_send_per_thread[ith][srank]++;
            }
        }  // end of OMP scope
        const auto n_remain = adr_ptcl_remain.size();
        const auto n_sml = std::min(adr_ptcl_send.size(), n_remain);
        for (auto i = 0; i < n_sml; i++) {
            const auto adr_remain = adr_ptcl_remain[n_remain - 1 - i];
            const auto adr_send = adr_ptcl_send[i];
            if (adr_remain < adr_send) break;
            std::swap(ptcl_[adr_send], ptcl_[adr_remain]);
        }
#else  // #if 0
       // non-parallel version
        S32 iloc = n_loc - 1;
        for (auto ip = n_loc - 1; ip >= 0; ip--) {
            if (dinfo.getPosRootDomain().notOverlapped(ptcl_[ip].getPos())) {
                PARTICLE_SIMULATOR_PRINT_ERROR("A particle is out of root domain");
                std::cerr << "position of the particle=" << ptcl_[ip].getPos() << std::endl;
                std::cerr << "position of the root domain=" << dinfo.getPosRootDomain() << std::endl;
                Abort(-1);
            }
            S32vec dnp;
            // S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain, comm_info_.getRank(), len_peri, pa, dnp,
            // dis_domain);
            S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain, comm_info_.getRank(), dnp);
            rank_send[ip] = my_rank;
            if (srank == my_rank) continue;
            std::swap(ptcl_[ip], ptcl_[iloc]);
            rank_send[iloc] = srank;
            iloc--;
            dnp_max_loc = ApplyEach([](const F64& l, const F64& r) { return std::max(l, r); }, dnp_max_loc, dnp);
            n_send[srank]++;
        }
        const auto n_remain = iloc + 1;
        n_disp_send[0] = 0;
        for (S32 i = 0; i < n_proc; i++) {
            n_disp_send[i + 1] += n_disp_send[i] + n_send[i];
        }
        ptcl_send.resizeNoInitialize(n_disp_send[n_proc]);
        // ****************************************************
        // *** align send particles on ptcl_send_ *************
        for (S32 i = 0; i < n_proc; i++) n_send[i] = 0;
        for (S32 ip = n_remain; ip < n_loc; ip++) {
            const auto srank = rank_send[ip];
            S32 jloc = n_send[srank] + n_disp_send[srank];
            ptcl_send[jloc] = ptcl_[ip];
            n_send[srank]++;
        }
#endif
        ptcl_.resizeNoInitialize(n_remain);
        time_profile_.exchange_particle_search += GetWtime() - t0;
        t0 = GetWtime();
        // ****************************************************
        // *** receive the number of receive particles ********
        //time_offset_inner = GetWtime();
        //F64 etime1;
        //F64 etime2;
        //F64 etime3;
        //	    Comm::barrier();
        //etime1 = GetWtime();
        auto dnp_max_glb = comm_info_.getMaxValue(dnp_max_loc);
        //etime2 = GetWtime();

        /*
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        if(2*dnp_max_glb[0]*2*dnp_max_glb[1] < nproc/4){
#else
        if(2*dnp_max_glb[0]*2*dnp_max_glb[1]*2*dnp_max_glb[2] < nproc/8){
#endif
        */
        if (1) {
            exchangeNumberOfPtcl(n_send, n_recv, n_domain, dnp_max_glb);
        } else {
            // comm_info_.allToAll(n_send.getPointer(), 1, n_recv.getPointer());
            comm_info_.allToAll(n_send, 1, n_recv);
        }
        n_disp_recv[0] = 0;
        for (auto i = 0; i < n_proc; i++) {
            n_disp_recv[i + 1] = n_disp_recv[i] + n_recv[i];
        }
        ptcl_.resizeNoInitialize(n_disp_recv[n_proc] + n_remain);
        /*
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        if(2*dnp_max_glb[0]*2*dnp_max_glb[1] < nproc/4){
#else
        if(2*dnp_max_glb[0]*2*dnp_max_glb[1]*2*dnp_max_glb[2] < nproc/8){
#endif
        */
        if (1) {
            S32 n_cnt = 0;
            for (S32 ib = 1; ib < n_proc; ib++) {
                S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                if (n_recv[idrecv] > 0) {
                    S32 adrrecv = n_disp_recv[idrecv] + n_remain;
                    S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                    MPI_Irecv(ptcl_.getPointer(adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv, comm_info_.getCommunicator(),
                              &req_recv[n_cnt]);
                    n_cnt++;
                }
            }
            //etime3 = GetWtime();
            for (S32 ib = 1; ib < n_proc; ib++) {
                S32 idsend = (ib + my_rank) % n_proc;
                if (n_send[idsend] > 0) {
                    S32 adrsend = n_disp_send[idsend];
                    S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                    MPI_Send(ptcl_send.getPointer(adrsend), n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend, comm_info_.getCommunicator());
                }
            }
            MPI_Waitall(n_cnt, &req_recv[0], &status[0]);

        } else {
            comm_info_.allToAllV(ptcl_send.getPointer(), n_send, n_disp_send, ptcl_.getPointer(n_remain), n_recv, n_disp_recv);
        }
        n_ptcl_send_ += n_disp_send[n_proc];
        n_ptcl_recv_ += n_disp_recv[n_proc];
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        n_proc_exch_ += 2 * dnp_max_glb[0] * 2 * dnp_max_glb[1];
#else
        n_proc_exch_ += 2 * dnp_max_glb[0] * 2 * dnp_max_glb[1] * 2 * dnp_max_glb[2];
#endif
        auto wt = GetWtime();
        time_profile_.exchange_particle_communication += wt - t0;
        ptcl_send.freeMem();
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
        n_ptcl_send_ = 0;
        n_ptcl_recv_ = 0;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
        time_profile_.exchange_particle += GetWtime() - time_offset;
    }

    template <class Tcomp>
    void sortParticle(Tcomp comp) {
        const S32 n = ptcl_.size();
        std::sort(ptcl_.getPointer(), ptcl_.getPointer(n), comp);
    }

    template <class Tcomp>
    void sortParticleStable(Tcomp comp) {
        const S32 n = ptcl_.size();
        std::stable_sort(ptcl_.getPointer(), ptcl_.getPointer(n), comp);
    }

    // for DEBUG functions
    template <class Treal, class Tvec>
    void calcCMDirect(Treal& mass_cm, Tvec& pos_cm) {
        mass_cm = 0.0;
        pos_cm = 0.0;
        const S32 n_ptcl = ptcl_.size();
        for (S32 i = 0; i < n_ptcl; i++) {
            mass_cm += ptcl_[i].mass;
            pos_cm += ptcl_[i].mass * ptcl_[i].pos;
        }
        pos_cm /= mass_cm;
    }

    bool checkExchangeParticleAllParticleInside(DomainInfo& dinfo);
    bool checkExchangeParticleSumOfNumberOfParticle(DomainInfo& dinfo, S32 ntot_init);

    void adjustPositionIntoRootDomain(const DomainInfo& dinfo) {
        const F64ort pos_root = dinfo.getPosRootDomain();
        const F64vec len_root = pos_root.getFullLength();
        const S32 n = ptcl_.size();
        PS_OMP_PARALLEL_FOR
        for (S32 i = 0; i < n; i++) {
            auto pos_new = ptcl_[i].getPos();
            while (pos_new.x < pos_root.low_.x) {
                pos_new.x += len_root.x;
            }
            while (pos_new.x >= pos_root.high_.x) {
                pos_new.x -= len_root.x;
            }
            if (pos_new.x == pos_root.high_.x) {
                pos_new.x = pos_root.low_.x;
            }
            while (pos_new.y < pos_root.low_.y) {
                pos_new.y += len_root.y;
            }
            while (pos_new.y >= pos_root.high_.y) {
                pos_new.y -= len_root.y;
            }
            if (pos_new.y == pos_root.high_.y) {
                pos_new.y = pos_root.low_.y;
            }
#if !defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
            while (pos_new.z < pos_root.low_.z) {
                pos_new.z += len_root.z;
            }
            while (pos_new.z >= pos_root.high_.z) {
                pos_new.z -= len_root.z;
            }
            if (pos_new.z == pos_root.high_.z) {
                pos_new.z = pos_root.low_.z;
            }
#endif
            ptcl_[i].setPos(pos_new);
        }
    }

    size_t getMemSizeUsed() const { return ptcl_.getMemSize(); }
    size_t getUsedMemorySize() const { return getMemSizeUsed(); }
    CountT getNumberOfParticleSendLocal() const { return (CountT)n_ptcl_send_; }
    CountT getNumberOfParticleRecvLocal() const { return (CountT)n_ptcl_recv_; }
    CountT getNumberOfParticleSendGlobal() const { return comm_info_.getSum((CountT)n_ptcl_send_); }
    CountT getNumberOfParticleRecvGlobal() const { return comm_info_.getSum((CountT)n_ptcl_recv_); }
    CountT getNumberOfProcExchangeLocal() const { return (CountT)n_proc_exch_; }
    CountT getNumberOfProcExchangeGlobal() const { return comm_info_.getSum((CountT)n_proc_exch_); }
    void clearCounterAll() {
        n_ptcl_send_ = n_ptcl_recv_ = 0;
        n_proc_exch_ = 0;
        time_profile_.clear();
    }

    //////////
    // add and remove particles
    void addOneParticle(const Tptcl& fp) { ptcl_.push_back(fp); }
    void removeParticle(const S32* idx, const S32 n_remove) {
        idx_remove_ptcl_.resizeNoInitialize(n_remove);
        for (S32 i = 0; i < n_remove; i++) {
            idx_remove_ptcl_[i] = idx[i];
        }
        std::sort(idx_remove_ptcl_.getPointer(), idx_remove_ptcl_.getPointer(n_remove));
        S32* ptr_end = std::unique(idx_remove_ptcl_.getPointer(), idx_remove_ptcl_.getPointer(n_remove));
        const S32 n_remove_tmp = ptr_end - idx_remove_ptcl_.getPointer();
        const S32 n_prev = ptcl_.size();
        S32 i_loc = n_prev - 1;
        for (S32 i = n_remove_tmp - 1; i >= 0; i--) {
            std::swap(ptcl_[idx_remove_ptcl_[i]], ptcl_[i_loc]);
            i_loc--;
        }
        ptcl_.resizeNoInitialize(i_loc + 1);
        /*
        // original
        const S32 n_prev = ptcl_.size();
        flag_remove.resizeNoInitialize(n_prev);
        for(S32 i=0; i<n_prev; i++){
            flag_remove[i] = false;
        }
        for(S32 i=0; i<n_remove; i++){
            S32 idx_tmp = idx[i];
            flag_remove[idx_tmp] = true;
        }
        S32 i_loc = n_prev-1;
        for(S32 i=n_prev-1; i>=0; i--){
            if(flag_remove[i] == true){
                std::swap(ptcl_[i], ptcl_[i_loc]);
                i_loc--;
            }
        }
        ptcl_.resizeNoInitialize(i_loc+1);
        */
    }

    void dumpMemSizeUsed(std::ostream& fout) {
        F64 sum_loc, sum_max;
        sum_loc = (double)(ptcl_.getMemSize() + idx_remove_ptcl_.getMemSize()) * 1e-9;
        //+ptcl_send_.getMemSize()
        //+ptcl_recv_.getMemSize()) * 1e-9;
        // if (Comm::getRank() == 0) {
        if (comm_info_.getRank() == 0) {
            fout << "ptcl_.getMemSize()= " << ptcl_.getMemSize() << std::endl;
            fout << "    ptcl_.size()= " << ptcl_.size() << std::endl;
            fout << "    sizeof(Tptcl)= " << sizeof(Tptcl) << std::endl;
            fout << "idx_remove_ptcl_.getMemSize()= " << idx_remove_ptcl_.getMemSize() << std::endl;
            fout << "sum[GB]= " << sum_loc << std::endl;
        }
        // sum_max = Comm::getMaxValue(sum_loc);
        sum_max = comm_info_.getMaxValue(sum_loc);
        // if (Comm::getRank() == 0) {
        if (comm_info_.getRank() == 0) {
            fout << "sum[GB](psys,max.) = " << sum_max << std::endl;
        }
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    void exchangeNumberOfPtcl(const S32 n_send[],  // const
                              S32 n_recv[], const S32 n_domain[], const S32vec& dnp_max) {
        const auto my_rank = comm_info_.getRank();
        const auto n_proc = comm_info_.getNumberOfProc();
        // ReallocatableArray<MPI_Request> req_recv(n_proc, n_proc, MemoryAllocMode::Stack);
        MPI_Request req_recv[n_proc];
        // ReallocatableArray<MPI_Status> status(n_proc, n_proc, MemoryAllocMode::Stack);
        MPI_Status status[n_proc];

        // auto rank_vec_org = GetRankVec(n_domain, my_rank, DIMENSION);
        auto rank_vec_org = GetRankVec(n_domain, my_rank);
        auto dnp_head = -dnp_max;
        auto dnp_tail = dnp_max;
        // if (Comm::getRank()==0){
        // 	std::cerr << "exchangenop" << dnp_max.x << " " << dnp_max.y
        //          << " " <<  dnp_max.z << "\n" ;
        // }
        for (auto k = 0; k < DIMENSION; k++) {
            if (dnp_max[k] * 2 + 1 > n_domain[k]) {
                dnp_head[k] = -n_domain[k] / 2;
                dnp_tail[k] = (n_domain[k] % 2 == 0) ? (n_domain[k] / 2 - 1) : n_domain[k] / 2;
            }
        }
#if 1
        ReallocatableArray<S32> rank_comm(n_proc, 0, MemoryAllocMode::Pool);
        S32vec rank_vec_tmp = 0;
        for (auto i = dnp_head.x; i <= dnp_tail.x; i++) {
            rank_vec_tmp.x = (rank_vec_org.x + n_domain[0] + i) % n_domain[0];
            for (auto j = dnp_head.y; j <= dnp_tail.y; j++) {
                rank_vec_tmp.y = (rank_vec_org.y + n_domain[1] + j) % n_domain[1];
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
                const auto rank_glb = GetRankGlb(rank_vec_tmp, n_domain);
                MPI_Irecv(&n_recv[rank_glb], 1, GetDataType<S32>(), rank_glb, 0, comm_info_.getCommunicator(), &req_recv[rank_comm.size()]);
                rank_comm.pushBackNoCheck(rank_glb);
#else
                for (auto k = dnp_head.z; k <= dnp_tail.z; k++) {
                    rank_vec_tmp.z = (rank_vec_org.z + n_domain[2] + k) % n_domain[2];
                    const auto rank_glb = GetRankGlb(rank_vec_tmp, n_domain);
                    MPI_Irecv(&n_recv[rank_glb], 1, GetDataType<S32>(), rank_glb, 0, comm_info_.getCommunicator(), &req_recv[rank_comm.size()]);
                    rank_comm.pushBackNoCheck(rank_glb);
                }
#endif
            }
        }
        for (auto i = 0; i < rank_comm.size(); i++) {
            MPI_Send(&n_send[rank_comm[i]], 1, GetDataType<S32>(), rank_comm[i], 0, comm_info_.getCommunicator());
        }
        MPI_Waitall(rank_comm.size(), &req_recv[0], &status[0]);
        rank_comm.freeMem();
#else
        const auto rank_x_org = my_rank / (n_domain[1] * n_domain[2]);
        const auto rank_y_org = (my_rank / n_domain[2]) % n_domain[1];
        const auto rank_z_org = my_rank % n_domain[2];
        S32 n_cnt = 0;
        S32 dnp_x_head = -dnp_max[0];
        S32 dnp_x_tail = dnp_max[0];
        if (dnp_max[0] * 2 + 1 > n_domain[0]) {
            dnp_x_head = -n_domain[0] / 2;
            dnp_x_tail = (n_domain[0] % 2 == 0) ? (n_domain[0] / 2 - 1) : n_domain[0] / 2;
        }
        S32 dnp_y_head = -dnp_max[1];
        S32 dnp_y_tail = dnp_max[1];
        if (dnp_max[1] * 2 + 1 > n_domain[1]) {
            dnp_y_head = -n_domain[1] / 2;
            dnp_y_tail = (n_domain[1] % 2 == 0) ? (n_domain[1] / 2 - 1) : n_domain[1] / 2;
        }
        S32 dnp_z_head = -dnp_max[2];
        S32 dnp_z_tail = dnp_max[2];
        if (dnp_max[2] * 2 + 1 > n_domain[2]) {
            dnp_z_head = -n_domain[2] / 2;
            dnp_z_tail = (n_domain[2] % 2 == 0) ? (n_domain[2] / 2 - 1) : n_domain[2] / 2;
        }
        for (S32 i = dnp_x_head; i <= dnp_x_tail; i++) {
            auto rank_x = (rank_x_org + n_domain[0] + i) % n_domain[0];
            for (S32 j = dnp_y_head; j <= dnp_y_tail; j++) {
                auto rank_y = (rank_y_org + n_domain[1] + j) % n_domain[1];
                for (S32 k = dnp_z_head; k <= dnp_z_tail; k++) {
                    auto rank_z = (rank_z_org + n_domain[2] + k) % n_domain[2];
                    const auto rank = rank_x * (n_domain[1] * n_domain[2]) + rank_y * n_domain[2] + rank_z;
                    MPI_Irecv(&n_recv[rank], 1, GetDataType<S32>(), rank, 0, comm_info_.getCommunicator(), &req_recv[n_cnt]);
                    n_cnt++;
                }
            }
        }
        n_cnt = 0;
        for (S32 i = dnp_x_head; i <= dnp_x_tail; i++) {
            auto rank_x = (rank_x_org + n_domain[0] + i) % n_domain[0];
            for (S32 j = dnp_y_head; j <= dnp_y_tail; j++) {
                auto rank_y = (rank_y_org + n_domain[1] + j) % n_domain[1];
                for (S32 k = dnp_z_head; k <= dnp_z_tail; k++) {
                    auto rank_z = (rank_z_org + n_domain[2] + k) % n_domain[2];
                    const auto rank = rank_x * (n_domain[1] * n_domain[2]) + rank_y * n_domain[2] + rank_z;
                    MPI_Send(&n_send[rank], 1, GetDataType<S32>(), rank, 0, comm_info_.getCommunicator());
                    n_cnt++;
                }
            }
        }
        MPI_Waitall(n_cnt, &req_recv[0], &status[0]);
#endif
    }
#endif

#if 0
        template<typename Tfunc_ep_ep, typename Tforce>
        void calcForceDirectParallelSimple(Tfunc_ep_ep pfunc_ep_ep,
                                           Tforce force[],
                                           const bool clear=true){
            // Simply recv all particles from all processes and
            // calculate interactions.  In the case of PERIODIC
            // BOUNDARY, must modify force kernel to calculate
            // interactions from not only real particles but also image particles.
            auto n_loc = ptcl_.size();
            if(clear){
                for(auto i=0; i<n_loc; i++) force[i].clear();
            }
            //Tptcl pi[n_loc];
            //Tptcl pj_send[n_loc];
            ReallocatableArray<Tptcl> pi(n_loc, n_loc, MemoryAllocMode::Pool);
            ReallocatableArray<Tptcl> pj_send(n_loc, n_loc, MemoryAllocMode::Pool);
            for(auto i=0; i<n_loc; i++){
                pi[i] = ptcl_[i];
                pj_send[i] = ptcl_[i];
            }
            const auto n_proc  = comm_info_.getNumberOfProc();
            const auto my_rank = comm_info_.getRank();
            for(auto i=0; i<n_proc; i++){
                const S32 n_image = 1;
                const auto rank_send = (my_rank + n_proc + i) % n_proc;
                const auto rank_recv = (my_rank + n_proc - i) % n_proc;
                S32 n_recv = 0;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                MPI_Status stat;
                MPI_Sendrecv(&n_loc, 1, GetDataType<S32>(), rank_send, 0, &n_recv, 1, GetDataType<S32>(), rank_recv, 0, comm_info_.getCommunicator(), &stat);
                //Tptcl pj_recv[n_recv*n_image];
                ReallocatableArray<Tptcl> pj_recv(n_recv*n_image, n_recv*n_image, MemoryAllocMode::Pool);
                MPI_Sendrecv(pj_send.getPointer(), n_loc, GetDataType<Tptcl>(), rank_send, 1, pj_recv.getPointer(), n_recv, GetDataType<Tptcl>(), rank_recv, 1, comm_info_.getCommunicator(), &stat);
#else
                n_recv = n_loc;
                Tptcl pj_recv[n_recv*n_image];
                for(S32 j=0; j<n_recv; j++){
                    pj_recv[j] = pj_send[j];
                }
#endif
                if(Comm::getRank()==0){
                    PS_PRINT_VARIABLE3(i, n_image, n_recv);
                }
                for(S32 j=0; j<n_image; j++){
                    for(S32 k=0; k<n_recv; k++){
                        pj_recv[j*n_image+k] = pj_recv[k];
                    }
                }
                pfunc_ep_ep(pi.getPointer(), n_loc, pj_recv.getPointer(), n_recv*n_image, force);
            }
        }
#endif
#if 0
        template<typename Tfunc_ep_ep, typename Tforce>
        void calcForceDirectParallel(Tfunc_ep_ep pfunc_ep_ep,
                                     Tforce force[],
                                     const DomainInfo & dinfo,
                                     const bool clear=true){
	    // Simply recv both real and image particles from all processes.
	    // Thus must not modify force kernel.
            auto n_loc = ptcl_.size();
            if(clear){
                for(auto i=0; i<n_loc; i++) force[i].clear();
            }
	    Tptcl pi[n_loc];
	    Tptcl pj_send[n_loc];
            for(auto i=0; i<n_loc; i++){
                pi[i] = ptcl_[i];
                pj_send[i] = ptcl_[i];
            }
            const auto n_proc  = comm_info_.getNumberOfProc();
            const auto my_rank = comm_info_.getRank();

            F64ort pos_root_domain = dinfo.getPosRootDomain();
            F64vec domain_size = pos_root_domain.getFullLength();
            bool pa[DIMENSION_LIMIT];
            dinfo.getPeriodicAxis(pa);
            F64ort my_box;
            for(S32 i=0; i<n_loc; i++){
                my_box.merge(ptcl_[i].getPos(), GetMyRSearch(ptcl_[i]));
            }
	    F64ort box_ar[n_proc];

            //Comm::allGather(&my_box, 1, box_ar.getPointer());
            comm_info_.allGather(&my_box, 1, box_ar);
	    //std::cerr<<"CHECK C"<<std::endl;
            for(S32 i=0; i<n_proc; i++){
                ReallocatableArray<F64vec> shift;
                CalcNumberAndShiftOfImageDomain(shift, domain_size, my_box, box_ar[i], pa);
                const S32 n_image = shift.size();
                auto rank_send = (my_rank + n_proc + i) % n_proc;
                auto rank_recv = (my_rank + n_proc - i) % n_proc;
                S32 n_recv = 0;
		//std::cerr<<"CHECK D"<<std::endl;
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                MPI_Status stat;
                //MPI_Sendrecv(&n_loc, 1, GetDataType<S32>(), rank_send, 0, &n_recv, 1, GetDataType<S32>(), rank_recv, 0, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv(&n_loc, 1, GetDataType<S32>(), rank_send, 0, &n_recv, 1, GetDataType<S32>(), rank_recv, 0, comm_info_.getCommunicator(), &stat);
		Tptcl pj_recv[n_recv*n_image];
                //MPI_Sendrecv(my_pj.getPointer(), n_loc, GetDataType<Tptcl>(), rank_send, 1, pj.getPointer(), n_recv, GetDataType<Tptcl>(), rank_recv, 1, MPI_COMM_WORLD, &stat);
                MPI_Sendrecv(pj_send, n_loc, GetDataType<Tptcl>(), rank_send, 1, pj_recv, n_recv, GetDataType<Tptcl>(), rank_recv, 1, comm_info_.getCommunicator(), &stat);
#else
		n_recv = n_loc;
		Tptcl pj_recv[n_recv*n_image];
                for(S32 j=0; j<n_recv; j++){
                    pj_recv[j] = pj_send[j];
                }
#endif
		//std::cerr<<"CHECK E"<<std::endl;
		if(Comm::getRank()==0){
		    PS_PRINT_VARIABLE3(i, n_image, n_recv);
		}
                for(S32 j=0; j<n_image; j++){
                    for(S32 k=0; k<n_recv; k++){
                        pj_recv[j*n_image+k] = pj_recv[k];
                        const F64vec pos_new = pj_recv[k].getPos() + shift[j];
                        pj_recv[j*n_image+k].setPos(pos_new);
                        //pj[j*n_image+k].pos = pos_new;
                    }
                }
                pfunc_ep_ep(pi, n_loc, pj_recv, n_recv*n_image, force);
            }

        }
#endif
};

template <typename T>
struct IsParticleSystem : std::false_type {};
template <typename T>
struct IsParticleSystem<ParticleSystem<T> > : std::true_type {};
template <typename T>
struct IsParticleSystem<const ParticleSystem<T> > : std::true_type {};
template <typename T>
struct IsParticleSystem<ParticleSystem<T>&> : std::true_type {};
template <typename T>
struct IsParticleSystem<const ParticleSystem<T>&> : std::true_type {};
template <typename T>
struct IsParticleSystem<ParticleSystem<T>&&> : std::true_type {};
template <typename T>
struct IsParticleSystem<const ParticleSystem<T>&&> : std::true_type {};
}  // namespace ParticleSimulator
