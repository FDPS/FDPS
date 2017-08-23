#pragma once

#include<cassert>
#include<fstream>

#include "MT.hpp"
#include"ps_defs.hpp"

namespace ParticleSimulator{
    template<class Tptcl>
    class ParticleSystem{
    private:
        CountT n_ptcl_send_;
        CountT n_ptcl_recv_;
        TimeProfile time_profile_;
        static const S32 n_smp_ave_ = 30;
        ReallocatableArray<Tptcl> ptcl_;
	ReallocatableArray<S32> idx_remove_ptcl_; // add 2016/09/14
        ReallocatableArray<Tptcl> ptcl_send_;
        ReallocatableArray<Tptcl> ptcl_recv_;
        S32 n_smp_ptcl_tot_;
        bool first_call_by_initialize;
        bool first_call_by_setAverageTargetNumberOfSampleParticlePerProcess;
        bool first_call_by_DomainInfo_collect_sample_particle;
        inline bool determineWhetherParticleIsInDomain(const F64vec & pos,
                                                       const F64ort & domain) {
            bool ret = true;
	    ret *= (domain.low_.x <= pos.x) * (pos.x < domain.high_.x);
	    ret *= (domain.low_.y <= pos.y) * (pos.y < domain.high_.y);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	    ret *= (domain.low_.z <= pos.z) * (pos.z < domain.high_.z);
#endif
            //for(S32 k = 0; k < DIMENSION; k++) ret *= (domain.low_[k] <= pos[k]) * (pos[k] < domain.high_[k]);
            return ret;
        }

        S32 searchWhichDomainParticleGoTo(const F64vec & pos,
                                          const S32 n_domain [],
                                          const F64ort domain []) {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            S32 idomain = 0;
            const S32 ny = n_domain[1];
            while(domain[idomain].high_.x <= pos.x)
                idomain += ny;
            while(domain[idomain].high_.y <= pos.y)
                idomain++;
            return idomain;
#else
            S32 idomain = 0;
            const S32 nynz = n_domain[1] * n_domain[2];
            while(domain[idomain].high_.x <= pos.x)
                idomain += nynz;
            const S32 nz   = n_domain[2];
            while(domain[idomain].high_.y <= pos.y)
                idomain += nz;
            while(domain[idomain].high_.z <= pos.z)
                idomain++;            
            return idomain;
#endif
        }

//        template<class Tvec, class Tdinfo>
//        S32 whereToGo(const Tvec & pos, const Tdinfo & dinfo);

    public:
        TimeProfile getTimeProfile() const {
            return time_profile_;
        }
        void clearTimeProfile(){
            time_profile_.clear();
        }
        ParticleSystem() {
            first_call_by_setAverageTargetNumberOfSampleParticlePerProcess = true;
            first_call_by_initialize = true;
            first_call_by_DomainInfo_collect_sample_particle = true;
            //n_smp_ptcl_tot_ = n_smp_ave_ * Comm::getNumberOfProc();
        }
	
        void initialize() {
            assert(first_call_by_initialize);
            first_call_by_initialize = false;
	    n_smp_ptcl_tot_ = n_smp_ave_ * Comm::getNumberOfProc();
            n_ptcl_send_ = n_ptcl_recv_ = 0;
	    //first_call_by_DomainInfo_collect_sample_particle = true;
	    //n_smp_ptcl_tot_ = n_smp_ave_ * Comm::getNumberOfProc();
        }

        void setAverageTargetNumberOfSampleParticlePerProcess(const S32 &nsampleperprocess) {
            assert(first_call_by_setAverageTargetNumberOfSampleParticlePerProcess);
            first_call_by_setAverageTargetNumberOfSampleParticlePerProcess = false;
            n_smp_ptcl_tot_ = nsampleperprocess * Comm::getNumberOfProc();
        }

        S32 getTargetNumberOfSampleParticle() {
            return n_smp_ptcl_tot_;
        }

        bool getFirstCallByDomainInfoCollectSampleParticle() {
            if(first_call_by_DomainInfo_collect_sample_particle) {
                first_call_by_DomainInfo_collect_sample_particle = false;
                return true;
            } else {
                return false;
            }
        }

        void createParticle(const S32 n_limit, bool clear=true){
            //n_ptcl_limit_ = n_limit;
            //ptcl_ = new Tptcl[n_ptcl_limit_];
            ptcl_.reserve(n_limit);
            ptcl_.resizeNoInitialize(0);
        }

	
        //void setNumberOfParticleLocal(const S32 n){ n_ptcl_ = n; }
        void setNumberOfParticleLocal(const S32 n){
            //15/02/20 Hosono bug(?) fix.
            //ptcl_.reserve(n*3+1000);
            ptcl_.resizeNoInitialize(n);
        }
        ////////////////
        // 05/01/30 Hosono From
        ////////////////
        //dummy class for the case if User does NOT define the file header.
        //TO BE PRIVATE
        struct DummyHeader{
            void writeAscii(FILE* fp) const{
            }
            S32 readAscii (FILE* fp){
                return -1;
            }
            void writeBinary(FILE* fp) const{
            }
            S32 readBinary (FILE* fp){
                return -1;
            }
        };
        //2016_11_03 modified IO functions to handle multiple file formats
        //read
        template <class Theader>
        void readParticleImpl(const char * const filename,
                              const char * const format,
                              Theader * const header,
                              void (Tptcl::*pFuncPtcl)(FILE*),
                              S32 (Theader::*pFuncHead)(FILE*),
                              const char * open_format){

            if(format == NULL){//Read from single file
                if(Comm::getRank() == 0){
                    FILE* fp = fopen(filename, open_format);
                    if(fp == NULL){
                        PARTICLE_SIMULATOR_PRINT_ERROR("can not open input file ");
                        std::cerr<<"filename: "<<filename<<std::endl;
                        Abort(-1);
                    }
                    //S32 n_ptcl_ = header->readAscii(fp);
                    S32 n_ptcl_ = (header->*pFuncHead)(fp);
                    while('\n' == getc(fp));
                    fseek(fp, -1, SEEK_CUR);
                    if(n_ptcl_ < 0){//User does NOT return # of ptcl
                        //count # of lines
                        n_ptcl_ = 0;
                        //KN
                        for(int c ; (c = getc(fp)) != EOF ; n_ptcl_ += '\n' == c ? 1 : 0){}
                        fclose(fp);
                        fp = fopen(filename, "r");
                        //header->readAscii(fp);
                        (header->*pFuncHead)(fp);
                        while('\n' == getc(fp));
			fseek(fp, -1, SEEK_CUR);
                    }
                    //Inform the # of ptcl for each process.
                    const S32 n_proc = Comm::getNumberOfProc();
                    S32 *n_ptcl = new S32[n_proc];
                    for(S32 i = 0 ; i < n_proc ; ++ i){
                        n_ptcl[i] = n_ptcl_ / n_proc;
                    }
                    n_ptcl[0] += n_ptcl_ % n_proc;
                    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    MPI::COMM_WORLD.Scatter(n_ptcl, 1, GetDataType<S32>(), &n_ptcl_, 1, GetDataType<S32>(), 0);
                    #endif
                    //allocate ptcl.
                    //First of all, Rank 0 reads its own particle.
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    for(int i = 0 ; i < n_ptcl_ ; ++ i){
                        //ptcl_[i].readAscii(fp);
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    //Read remaining data to buffer and send them to appropriate process.
                    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    for(S32 rank = 1 ; rank < n_proc ; ++ rank){
                        Tptcl * buffer = new Tptcl[n_ptcl[rank]];
                        for(int i = 0 ; i < n_ptcl[rank] ; ++ i){
                            //buffer[i].readAscii(fp);
                            (buffer[i].*pFuncPtcl)(fp);
                        }
                        MPI::COMM_WORLD.Send(buffer, n_ptcl[rank], GetDataType<Tptcl>(), rank, 0);
                        delete [] buffer;
                    }
                    #endif
                    //End.
                    delete [] n_ptcl;
                    fclose(fp);
                }else{
                    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    //Receive the # of ptcl from Rank 0
                    S32 n_ptcl_loc;
                    S32 *n_ptcl = new S32[Comm::getNumberOfProc()];
                    MPI::COMM_WORLD.Scatter(n_ptcl, 1, GetDataType<S32>(), &n_ptcl_loc, 1, GetDataType<S32>(), 0);
                    delete [] n_ptcl;
                    //allocate ptcl.
                    this->createParticle(n_ptcl_loc << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_loc);
                    MPI::COMM_WORLD.Recv(ptcl_.getPointer(), ptcl_.size(), GetDataType<Tptcl>(), 0, 0);
                    #endif
                }
            }else{//Read from multiple file
                char input[256];
                sprintf(input, format, filename, Comm::getNumberOfProc(), Comm::getRank());
                FILE* fp = fopen(input, "r");
                if(fp == NULL){
                    PARTICLE_SIMULATOR_PRINT_ERROR("can not open input file");
                    std::cerr<<"filename: "<<input<<std::endl;
                    Abort(-1);
                }
                //S32 n_ptcl_ = header->readAscii(fp);
                S32 n_ptcl_ = (header->*pFuncHead)(fp);
                while('\n' == getc(fp));
                fseek(fp, -1, SEEK_CUR);
                if(n_ptcl_ >= 0){
                    //User returns # of ptcl.
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    for(S32 i = 0 ; i < n_ptcl_ ; ++ i){
                        //ptcl_[i].readAscii(fp);
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }else{//User does NOT return # of ptcl
                    //count # of lines
                    n_ptcl_ = 0;
                    for(int c ; (c = getc(fp)) != EOF ; n_ptcl_ += c == '\n' ? 1 : 0){}
                    fclose(fp);
                    //
                    FILE* fp = fopen(input, "r");
                    //header->readAscii(fp);
                    (header->*pFuncHead)(fp);
                    while('\n' == getc(fp));
		    fseek(fp, -1, SEEK_CUR);
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    for(S32 i = 0 ; i < ptcl_.size() ; ++ i){
                        //ptcl_[i].readAscii(fp);
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }
            }
        }

        template <class Theader>
        void readParticleAscii(const char * const filename, const char * const format, Theader& header){
            readParticleImpl(filename, format, &header, &Tptcl::readAscii, &Theader::readAscii, "r");
        }
        template <class Theader>
        void readParticleAscii(const char * const filename, Theader& header){
            readParticleImpl(filename, NULL, &header, &Tptcl::readAscii, &Theader::readAscii, "r");
        }
        void readParticleAscii(const char * const filename, const char * const format){
            readParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::readAscii, &DummyHeader::readAscii, "r");
        }
        void readParticleAscii(const char * const filename){
            readParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::readAscii, &DummyHeader::readAscii, "r");
        }

        template <class Theader>
        void readParticleAscii(const char * const filename, const char * const format, Theader& header, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl(filename, format, &header, pFunc, &Theader::readAscii, "r");
        }
        template <class Theader>
        void readParticleAscii(const char * const filename, Theader& header, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl(filename, NULL, &header, pFunc, &Theader::readAscii, "r");
        }
        void readParticleAscii(const char * const filename, const char * const format, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::readAscii, "r");
        }
        void readParticleAscii(const char * const filename, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::readAscii, "r");
        }

        template <class Theader>
        void readParticleBinary(const char * const filename, const char * const format, Theader& header){
            readParticleImpl(filename, format, &header, &Tptcl::readBinary, &Theader::readBinary, "rb");
        }
        template <class Theader>
        void readParticleBinary(const char * const filename, Theader& header){
            readParticleImpl(filename, NULL, &header, &Tptcl::readBinary, &Theader::readBinary, "rb");
        }
        void readParticleBinary(const char * const filename, const char * const format){
            readParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::readBinary, &DummyHeader::readBinary, "rb");
        }
        void readParticleBinary(const char * const filename){
            readParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::readBinary, &DummyHeader::readBinary, "rb");
        }
        
        template <class Theader>
        void readParticleBinary(const char * const filename, const char * const format, Theader& header, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl(filename, format, &header, pFunc, &Theader::readBinary, "rb");
        }
        template <class Theader>
        void readParticleBinary(const char * const filename, Theader& header, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl(filename, NULL, &header, pFunc, &Theader::readBinary, "rb");
        }
        void readParticleBinary(const char * const filename, const char * const format, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::readBinary, "rb");
        }
        void readParticleBinary(const char * const filename, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::readBinary, "rb");
        }

        template <class Theader>
        void writeParticleImpl(const char * const filename,
                               const char * const format,
                               const Theader * const header,
                               void (Tptcl::*pFuncPtcl)(FILE*)const,
                               void (Theader::*pFuncHead)(FILE*)const,
                               const char * open_format){
            if(format == NULL){
                #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                //declare local # of ptcl.
                const S32 n_ptcl_ = ptcl_.size();
                //get # of process.
                const S32 n_proc = Comm::getNumberOfProc();
                //get # of ptcls in each process.
                S32 *n_ptcl = new S32[n_proc];
                //Gather # of particles.
                MPI::COMM_WORLD.Allgather(&n_ptcl_, 1, GetDataType<S32>(), n_ptcl, 1, GetDataType<S32>());
                //set displacement
                S32 *n_ptcl_displs = new S32[n_proc+1];
                n_ptcl_displs[0] = 0;
                for(S32 i = 0 ; i < n_proc ; ++ i){
                    n_ptcl_displs[i+1] = n_ptcl_displs[i] + n_ptcl[i];
                }
                const S32 n_tot = n_ptcl_displs[n_proc];
                Tptcl *ptcl = new Tptcl[n_tot];
                //gather data
                MPI::COMM_WORLD.Gatherv(ptcl_.getPointer(), n_ptcl_, GetDataType<Tptcl>(), ptcl, n_ptcl, n_ptcl_displs, GetDataType<Tptcl>(), 0);
                if(Comm::getRank() == 0){
                    FILE* fp = fopen(filename, open_format);
                    if(fp == NULL){
                        PARTICLE_SIMULATOR_PRINT_ERROR("can not open output file");
                        std::cerr<<"output file: "<<filename<<std::endl;
                        Abort(-1);
                    }
                    //header->writeAscii(fp);
                    (header->*pFuncHead)(fp);
                    for(S32 i = 0 ; i < n_tot ; ++ i){
                        //ptcl[i].writeAscii(fp);
                        (ptcl[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }
                delete [] n_ptcl;
                delete [] n_ptcl_displs;
                delete [] ptcl;
                #else
                const S32 n_tot = ptcl_.size();
                if(Comm::getRank() == 0){
                    FILE* fp = fopen(filename, open_format);
                    //header->writeAscii(fp);
                    (header->*pFuncHead)(fp);
                    for(S32 i = 0 ; i < n_tot ; ++ i){
                        //ptcl_[i].writeAscii(fp);
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }
                #endif
            }else{
                char output[256];
                sprintf(output, format, filename, Comm::getNumberOfProc(), Comm::getRank());
                FILE* fp = fopen(output, open_format);
                if(fp == NULL){
                    PARTICLE_SIMULATOR_PRINT_ERROR("can not open output file");
		    std::cerr<<"output file: "<<output<<std::endl;
                    Abort(-1);
                }
                //header->writeAscii(fp);
                (header->*pFuncHead)(fp);
                for(S32 i = 0 ; i < ptcl_.size() ; ++ i){
                    //ptcl_[i].writeAscii(fp);
                    (ptcl_[i].*pFuncPtcl)(fp);
                }
                fclose(fp);
            }
        }
        //write
        template <class Theader>
        void writeParticleAscii(const char * const filename, const char * const format, const Theader& header){
            writeParticleImpl(filename, format, &header, &Tptcl::writeAscii, &Theader::writeAscii, "w");
        }
        template <class Theader>
        void writeParticleAscii(const char * const filename, const Theader& header){
            writeParticleImpl(filename, NULL, &header, &Tptcl::writeAscii, &Theader::writeAscii, "w");
        }
        void writeParticleAscii(const char * const filename, const char * format){
            writeParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::writeAscii, &DummyHeader::writeAscii, "w");
        }
        void writeParticleAscii(const char * const filename){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::writeAscii, &DummyHeader::writeAscii, "w");
        }

        template <class Theader>
        void writeParticleAscii(const char * const filename, const char * const format, const Theader& header, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl(filename, format, &header, pFunc, &Theader::writeAscii, "w");
        }
        template <class Theader>
        void writeParticleAscii(const char * const filename, const Theader& header, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl(filename, NULL, &header, pFunc, &Theader::writeAscii, "w");
        }
        void writeParticleAscii(const char * const filename, const char * format, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::writeAscii, "w");
        }
        void writeParticleAscii(const char * const filename, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::writeAscii, "w");
        }


        template <class Theader>
        void writeParticleBinary(const char * const filename, const char * const format, const Theader& header){
            writeParticleImpl(filename, format, &header, &Tptcl::writeBinary, &Theader::writeBinary, "wb");
        }
        template <class Theader>
        void writeParticleBinary(const char * const filename, const Theader& header){
            writeParticleImpl(filename, NULL, &header, &Tptcl::writeBinary, &Theader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename, const char * format){
            writeParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
        }

        template <class Theader>
        void writeParticleBinary(const char * const filename, const char * const format, const Theader& header, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl(filename, format, &header, pFunc, &Theader::writeBinary, "wb");
        }
        template <class Theader>
        void writeParticleBinary(const char * const filename, const Theader& header, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl(filename, NULL, &header, pFunc, &Theader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename, const char * format, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
        }
        

        
	
        Tptcl & operator [] (const S32 id) {return ptcl_[id];}
        const Tptcl & operator [] (const S32 id) const {return ptcl_[id];}
        Tptcl & getParticle(const S32 id=0) {return ptcl_[id];}
        const Tptcl & getParticle(const S32 id=0) const {return ptcl_[id];}
        Tptcl * getParticlePointer(const S32 id=0) const {return ptcl_+id;}
        //S32 getNumberOfParticleLocal() const {return n_ptcl_;}
        S32 getNumberOfParticleLocal() const {return ptcl_.size();}
        ////////////////
        // 05/02/04 Hosono From
        // 11th jul MI modified
        ////////////////
        S64 getNumberOfParticleGlobal() const {
/*
            #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            //get # of process.
            const S32 n_proc = Comm::getNumberOfProc();
            //get # of ptcls in each process.
            S64 *n_ptcl = new S64[n_proc];
            //Gather # of particles.
            S64 n_ptcl_ = ptcl_.size();
            MPI::COMM_WORLD.Allgather(&n_ptcl_, 1, GetDataType<S64>(), n_ptcl, 1, GetDataType<S64>());
            //set displacement
            S32 *n_ptcl_displs = new S32[n_proc+1];
            n_ptcl_displs[0] = 0;
            for(S32 i = 0 ; i < n_proc ; ++ i){
                n_ptcl_displs[i+1] = n_ptcl_displs[i] + n_ptcl[i];
            }
            return n_ptcl_displs[n_proc];
            #else
            return ptcl_.size();
            #endif
*/
            #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const S64 n_ptcl_loc = ptcl_.size();
            return Comm::getSum(n_ptcl_loc);
            #else
            return ptcl_.size();
            #endif
        }
        ////////////////
        // 05/02/04 Hosono To
        ////////////////

        F64 getHalfLength(const F64vec & center=0.0){
            F64 hl_max_loc = (ptcl_[0].getPos() - center).applyEach(Abs<F64>()).getMax();
            const S32 n_ptcl = ptcl_.size();
            for(size_t i=1; i<n_ptcl; i++){
                F64 hl_tmp = (ptcl_[i].getPos() - center).applyEach(Abs<F64>()).getMax();
                hl_max_loc = (hl_max_loc > hl_tmp) ? hl_max_loc : hl_tmp;
            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            F64 hl_max_glb;
            MPI::COMM_WORLD.Allreduce(&hl_max_loc, &hl_max_glb, 1, GetDataType<F64>(), MPI::MAX);
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
            //return genrand_res53();
            return MT::genrand_res53();
        }
// *******************************************************************        
// *******************************************************************        

        inline S32 getUniformDistributionFromArg1ToArg2(S32 arg1,
                                                        S32 arg2) {
            S32 random_number;
            
            random_number = (S32)((arg2 - arg1 + 1) * frand()) + arg1;
            
            return random_number;
        }

        /* AT_DEBU
        void getSampleParticle(S32 & number_of_sample_particle,
                               F32vec pos_sample[],
                               const F32 weight=1.0) {
        */
        void getSampleParticle(S32 & number_of_sample_particle,
                               F64vec pos_sample[],
                               const F32 weight) {

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL

            S64 nloc = (S64)ptcl_.size();

            F32 weight_all = Comm::getSum(weight);
            number_of_sample_particle = (weight * n_smp_ptcl_tot_) / weight_all;
#if 0
// modified to limit # of sample particles by M.I. 
            const F32 coef_limitter = 0.2;
            S64 nglb = Comm::getSum( (S64)nloc );
            S32 number_of_sample_particle_limit = ((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter)); // lower limit
            number_of_sample_particle = (number_of_sample_particle > number_of_sample_particle_limit) ? number_of_sample_particle : number_of_sample_particle_limit;
#endif
            number_of_sample_particle = (number_of_sample_particle < nloc) ? number_of_sample_particle : nloc;
            //std::cout<<"number_of_sample_particle= "<<number_of_sample_particle<<std::endl;
            //std::cout<<"weight= "<<weight<<" weight_all= "<<weight_all<<std::endl;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            std::cout<<"weight="<<weight<<" weight_all="<<weight_all<<std::endl;
            std::cout<<"n_smp_ptcl_tot_="<<n_smp_ptcl_tot_<<std::endl;
            //std::cout<<"nglb="<<nglb<<" nloc="<<nloc<<std::endl;
            std::cout<<"(weight * n_smp_ptcl_tot_) / weight_all="<<(weight * n_smp_ptcl_tot_) / weight_all<<std::endl;
            std::cout<<"((S64)nloc * n_smp_ptcl_tot_)="<< ((S64)nloc * n_smp_ptcl_tot_)<<std::endl;
            //std::cout<<"((F32)nglb * (1.0 + coef_limitter))="<<((F32)nglb * (1.0 + coef_limitter))<<std::endl;
            //std::cout<<"((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter))="<<((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter))<<std::endl;
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



#if 1
      // new version (with switch)
      // for search, use tree with level 3(x, y, z) 
      // must be consistend with geometry of domains.
        template<class Tdinfo>
        void exchangeParticle(Tdinfo & dinfo) {
            F64 time_offset = GetWtime();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const S32 nloc  = ptcl_.size();
            const S32 rank  = MPI::COMM_WORLD.Get_rank();
            const S32 nproc = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();

            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(rank);

            S32 * nsend  = new S32[nproc];
            S32 * nsend_disp  = new S32[nproc+1];
            S32 * nrecv  = new S32[nproc];
            S32 * nrecv_disp  = new S32[nproc+1];
            MPI::Request * req_send = new MPI::Request[nproc];
            MPI::Request * req_recv = new MPI::Request[nproc];
            for(S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
            nsend_disp[nproc] = nrecv_disp[nproc] = 0;
            F64 time_offset_inner = GetWtime();
            for(S32 ip = 0; ip < nloc; ip++) {
                if( dinfo.getPosRootDomain().notOverlapped(ptcl_[ip].getPos()) ){
                    PARTICLE_SIMULATOR_PRINT_ERROR("A particle is out of root domain");
                    std::cerr<<"position of the particle="<<ptcl_[ip].getPos()<<std::endl;
                    std::cerr<<"position of the root domain="<<dinfo.getPosRootDomain()<<std::endl;
                    Abort(-1);
                }
                if(!determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    nsend[srank]++;
                }
            }

            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] += nsend_disp[i] + nsend[i];
            }
            ptcl_send_.resizeNoInitialize( nsend_disp[nproc] );
            // ****************************************************
            // *** align send particles on ptcl_send_ *************
            for(S32 i = 0; i < nproc; i++) nsend[i] = 0;
            S32 iloc = 0;
            for(S32 ip = 0; ip < nloc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    iloc++;
                } else {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    S32 jloc = nsend[srank] + nsend_disp[srank];
                    ptcl_send_[jloc] = ptcl_[ip];
                    nsend[srank]++;
                }
            }
            ptcl_.resizeNoInitialize(iloc);
            //time_profile_.exchange_particle__find_particle = GetWtime() - time_offset_inner;
            time_profile_.exchange_particle__find_particle += GetWtime() - time_offset_inner;

            // ****************************************************
            // *** receive the number of receive particles ********
            time_offset_inner = GetWtime();
            Comm::allToAll(nsend, 1, nrecv);

            nrecv_disp[0] = 0;
            for(S32 i=0; i<nproc; i++){
                nrecv_disp[i+1] = nrecv_disp[i] + nrecv[i];
            }
            ptcl_recv_.resizeNoInitialize( nrecv_disp[nproc] );

            const S32 n_proc_comm_limit = 500;
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0) n_proc_send++;
                if(nrecv[i] > 0) n_proc_recv++;
            }
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;

            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idsend = (ib + rank) % nproc;
                    if(nsend[idsend] > 0) {
                        S32 adrsend = nsend_disp[idsend];
                        S32 tagsend = (rank < idsend) ? rank : idsend;
                        req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    }
                    S32 idrecv = (nproc + rank - ib) % nproc;
                    if(nrecv[idrecv] > 0) {
                        S32 adrrecv = nrecv_disp[idrecv];
                        S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                        req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    }
                }
                MPI::Request::Waitall(n_proc_send, req_send);
                MPI::Request::Waitall(n_proc_recv, req_recv);
            }
            else{
                Comm::allToAllV(ptcl_send_.getPointer(), nsend, nsend_disp,
                                ptcl_recv_.getPointer(), nrecv, nrecv_disp);
            }

            const S32 nrecv_tot = nrecv_disp[nproc];
            ptcl_.reserveEmptyAreaAtLeast( nrecv_tot );
            for(S32 ip = 0; ip < nrecv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
            // **************************************************** 
            n_ptcl_send_ += nsend_disp[nproc];
            n_ptcl_recv_ += nrecv_disp[nproc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - time_offset_inner;

            delete [] nsend;
            delete [] nsend_disp;
            delete [] nrecv;
            delete [] nrecv_disp;
            delete [] req_send;
            delete [] req_recv;
#else
            n_ptcl_send_ = 0;
            n_ptcl_recv_ = 0;
#endif
            time_profile_.exchange_particle += GetWtime() - time_offset;
        }
#else

        template<class Tdinfo>
        void exchangeParticle(Tdinfo & dinfo) {
            F64 time_offset = GetWtime();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            //const S32 nloc  = n_ptcl_;
            const S32 nloc  = ptcl_.size();
            const S32 rank  = MPI::COMM_WORLD.Get_rank();
            const S32 nproc = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();

            /* AT_DEBUG
            F32ort * pos_domain = dinfo.getPointerOfPosDomain();
            F32ort thisdomain = dinfo.getPosDomain(rank);
            */
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(rank);

            S32 nsendtot = 0;
            S32 *nsend0  = new S32[nproc];
            S32 *nsend1  = new S32[nproc];
            S32 nrecvtot = 0;
            S32 *nrecv0  = new S32[nproc];
            S32 *nrecv1  = new S32[nproc];
            for(S32 i = 0; i < nproc; i++) {
                nsend0[i] = 0;
                nsend1[i] = 0;
                nrecv0[i] = 0;
                nrecv1[i] = 0;
            }
            MPI::Request * req_send = new MPI::Request[nproc];
            MPI::Request * req_recv = new MPI::Request[nproc];
            // *** count the number of send particles preliminary *
            for(S32 ip = 0; ip < nloc; ip++) {
                if( dinfo.getPosRootDomain().notOverlapped(ptcl_[ip].getPos()) ){
                    PARTICLE_SIMULATOR_PRINT_ERROR("A particle is out of root domain");
                    std::cerr<<"position of the particle="<<ptcl_[ip].getPos()<<std::endl;
                    std::cerr<<"position of the root domain="<<dinfo.getPosRootDomain()<<std::endl;
                    Abort(-1);
                }
                if(!determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    nsend1[srank]++;
                    nsendtot++;
                }
            }
            nsend0[0] = 0;
            for(S32 i = 1; i < nproc; i++) {
                nsend0[i] += nsend0[i-1] + nsend1[i-1];
            }
            //ptcl_send_ = new Tptcl[nsendtot];
            ptcl_send_.resizeNoInitialize(nsendtot);
            // ****************************************************
            // *** align send particles on ptcl_send_ *************
            //std::cerr<<"check 0"<<std::endl;
            for(S32 i = 0; i < nproc; i++)
                nsend1[i] = 0;
            S32 iloc = 0;
            for(S32 ip = 0; ip < nloc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    iloc++;
                } else {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    S32 jloc = nsend0[srank] + nsend1[srank];
                    ptcl_send_[jloc] = ptcl_[ip];
                    nsend1[srank]++;
                }
            }           
            //n_ptcl_ = iloc;
	    //std::cerr<<"check 1"<<std::endl;
	    ptcl_.resizeNoInitialize(iloc);
	    //std::cerr<<"ptcl_.size()="<<ptcl_.size()<<std::endl;
/*
            if(rank == 0) {
                char filename[1024];
                FILE *fp;
                sprintf(filename, "out/send_%04d_%04d.txt", rank, rank);
                fp = fopen(filename, "w");
                for(S32 j = 0; j < n_ptcl_; j++)
                    fprintf(fp, "%+e %+e %+e\n",
                            ptcl_[j].getPos()[0],
                            ptcl_[j].getPos()[1],
                            ptcl_[j].getPos()[2]);
                fclose(fp);
                for(S32 i = 1; i < nproc; i++) {
                    S32 srank = (rank + i) % nproc;
                    sprintf(filename, "out/send_%04d_%04d.txt", rank, srank);
                    fp = fopen(filename, "w");
                    S32 next = (srank + 1 < nproc) ? nsend0[srank+1] : nsendtot;
                    for(S32 j = nsend0[srank]; j < next; j++)
                        fprintf(fp, "%+e %+e %+e\n",
                                ptcl_send_[j].getPos()[0],
                                ptcl_send_[j].getPos()[1],
                                ptcl_send_[j].getPos()[2]);
                    fclose(fp);
                }
            }
*/
            // ****************************************************
            // *** receive the number of receive particles ********
            MPI::COMM_WORLD.Alltoall(nsend1, 1, GetDataType<S32>(), nrecv1, 1, GetDataType<S32>());
	    //std::cerr<<"check 2"<<std::endl;
            for(S32 i = 0; i < nproc; i++)
                nrecvtot += nrecv1[i];
            //assert(n_ptcl_ + nrecvtot <= n_ptcl_limit_);
	    //ptcl_.reserve(n_ptcl_ + nrecvtot);
	    //ptcl_.resizeNoInitialize(n_ptcl_ + nrecvtot);
            //ptcl_recv_ = new Tptcl[nrecvtot];
	    ptcl_recv_.resizeNoInitialize(nrecvtot);
	    //std::cerr<<"check 3"<<std::endl;
/*
            {
                char filename[1024];
                FILE *fp;
                sprintf(filename, "out/send_%04d.txt", rank);
                fp = fopen(filename, "w");
                fprintf(fp, "%d", rank);
                for(S32 i = 0; i < nproc; i++) {
                    fprintf(fp, "%6d", nsend1[i]);
                }
                fprintf(fp, "\n");
                fclose(fp);
                sprintf(filename, "out/recv_%04d.txt", rank);
                fp = fopen(filename, "w");
                fprintf(fp, "%d", rank);
                for(S32 i = 0; i < nproc; i++) {
                    fprintf(fp, "%6d", nrecv[i]);
                }
                fprintf(fp, "\n");
                fclose(fp);
            }
*/
            // ****************************************************
            // *** send and receive particles *********************
            nrecv0[0] = 0;
            for(S32 i = 1; i < nproc; i++) {
                nrecv0[i] += nrecv0[i-1] + nrecv1[i-1];
            }
            S32 nsendnode = 0;
            S32 nrecvnode = 0;
	    //std::cerr<<"check 4"<<std::endl;
            for(S32 ib = 1; ib < nproc; ib++) {
                S32 idsend = (ib + rank) % nproc;
                if(nsend1[idsend] > 0) {
                    S32 adrsend = nsend0[idsend];                    
                    S32 tagsend = (rank < idsend) ? rank : idsend;                    
                    //req_send[nsendnode] = MPI::COMM_WORLD.Isend(ptcl_send_+adrsend, nsend1[idsend]*sizeof(Tptcl), MPI::BYTE, idsend, tagsend);
		    req_send[nsendnode] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend1[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    nsendnode++;
                }                
                S32 idrecv = (nproc + rank - ib) % nproc;
                if(nrecv1[idrecv] > 0) {
                    S32 adrrecv = nrecv0[idrecv];
                    S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                    //req_recv[nrecvnode] = MPI::COMM_WORLD.Irecv(ptcl_recv_+adrrecv, nrecv1[idrecv]*sizeof(Tptcl), MPI::BYTE, idrecv, tagrecv);
		    req_recv[nrecvnode] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv1[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    nrecvnode++;
                }
            }
            MPI::Request::Waitall(nsendnode, req_send);
            MPI::Request::Waitall(nrecvnode, req_recv);
            //std::cerr<<"check 5"<<std::endl;
            // ****************************************************            
            // *** align particles ********************************
            /*
              for(S32 ip = 0; ip < nrecvtot; ip++) {
              ptcl_[n_ptcl_] = ptcl_recv_[ip];
              n_ptcl_++;
            }
	    */
            ptcl_.reserve( ptcl_.size()+nrecvtot );
            //ptcl_.dump("dump ptcl");
            for(S32 ip = 0; ip < nrecvtot; ip++) {
                ptcl_.pushBackNoCheck(ptcl_recv_[ip]);
            }
	    //std::cerr<<"ptcl_.size()="<<ptcl_.size()<<std::endl;
            // ****************************************************            

            delete [] nsend0;
            delete [] nsend1;
            delete [] nrecv0;
            delete [] nrecv1;
            delete [] req_send;
            delete [] req_recv;
            //delete [] ptcl_send_;
            //delete [] ptcl_recv_;
#endif
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            std::cout<<"ptcl_.size()="<<ptcl_.size()<<std::endl;
            if(ptcl_.size() > 0){
                std::cout<<"ptcl_[0].getPos()="<<ptcl_[0].getPos()<<std::endl;
                std::cout<<"ptcl_[0].getRSearch()="<<ptcl_[0].getRSearch()<<std::endl;
            }
#endif
            //time_profile_.exchange_particle = GetWtime() - time_offset;
            time_profile_.exchange_particle += GetWtime() - time_offset;
        }
#endif

        // for DEBUG functions
        template<class Treal, class Tvec>
        void calcCMDirect(Treal & mass_cm, Tvec & pos_cm){
            mass_cm = 0.0;
            pos_cm = 0.0;
            const S32 n_ptcl = ptcl_.size();
            for(S32 i=0; i<n_ptcl; i++){
                mass_cm += ptcl_[i].mass;
                pos_cm += ptcl_[i].mass * ptcl_[i].pos;
            }
            pos_cm /= mass_cm;
        }

        bool checkExchangeParticleAllParticleInside(DomainInfo & dinfo);
        bool checkExchangeParticleSumOfNumberOfParticle(DomainInfo & dinfo,
                                                        S32 ntot_init);

        void adjustPositionIntoRootDomain(const DomainInfo & dinfo){
            const F64ort pos_root = dinfo.getPosRootDomain();
            const F64vec len_root = pos_root.getFullLength();
            const S32 n = ptcl_.size();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n; i++){
                F64vec pos_new = ptcl_[i].getPos() ;
                //if( pos_root.notOverlapped(pos_new) ){
                    while(pos_new.x < pos_root.low_.x){
                        pos_new.x += len_root.x;
                    }
                    while(pos_new.x > pos_root.high_.x){
                        pos_new.x -= len_root.x;
                    }
                    if(pos_new.x == pos_root.high_.x){
                        pos_new.x = pos_root.low_.x;
                    }
                    while(pos_new.y < pos_root.low_.y){
                        pos_new.y += len_root.y;
                    }
                    while(pos_new.y >= pos_root.high_.y){
                        pos_new.y -= len_root.y;
                    }
                    if(pos_new.y == pos_root.high_.y){
                        pos_new.y = pos_root.low_.y;
                    }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                    while(pos_new.z < pos_root.low_.z){
                        pos_new.z += len_root.z;
                    }
                    while(pos_new.z >= pos_root.high_.z){
                        pos_new.z -= len_root.z;
                    }
                    if(pos_new.z == pos_root.high_.z){
                        pos_new.z = pos_root.low_.z;
                    }
#endif
                    //}
                ptcl_[i].setPos(pos_new);
            }
        }

        size_t getMemSizeUsed() const {
            return ptcl_.getMemSize() + ptcl_send_.getMemSize() + ptcl_recv_.getMemSize();
        }
	/*
        CountT getNumberOfParticleSendLocal() const { return (CountT)ptcl_send_.size(); }
        CountT getNumberOfParticleRecvLocal() const { return (CountT)ptcl_recv_.size(); }
        CountT getNumberOfParticleSendGlobal() const { return Comm::getSum((CountT)ptcl_send_.size()); }
        CountT getNumberOfParticleRecvGlobal() const { return Comm::getSum((CountT)ptcl_recv_.size()); }
	*/
        CountT getNumberOfParticleSendLocal() const { return (CountT)n_ptcl_send_; }
        CountT getNumberOfParticleRecvLocal() const { return (CountT)n_ptcl_recv_; }
        CountT getNumberOfParticleSendGlobal() const { return Comm::getSum((CountT)n_ptcl_send_); }
        CountT getNumberOfParticleRecvGlobal() const { return Comm::getSum((CountT)n_ptcl_recv_); }
        void clearCounterAll(){
            n_ptcl_send_ = n_ptcl_recv_ = 0;
            time_profile_.clear();
        }

	//////////
	// add and remove particles
	void addOneParticle(const Tptcl & fp){
	    ptcl_.push_back(fp);
	}
	void removeParticle(const S32 * idx, const S32 n_remove){
	    idx_remove_ptcl_.resizeNoInitialize(n_remove);
	    for(S32 i=0; i<n_remove; i++){
		idx_remove_ptcl_[i] = idx[i];
	    }
	    std::sort(idx_remove_ptcl_.getPointer(), idx_remove_ptcl_.getPointer(n_remove));
	    S32 * ptr_end = std::unique(idx_remove_ptcl_.getPointer(), idx_remove_ptcl_.getPointer(n_remove));
	    const S32 n_remove_tmp = ptr_end - idx_remove_ptcl_.getPointer();
	    const S32 n_prev = ptcl_.size();
	    S32 i_loc = n_prev-1;
	    for(S32 i=n_remove_tmp-1; i>=0; i--){
		std::swap(ptcl_[idx_remove_ptcl_[i]], ptcl_[i_loc]);
		i_loc--;
	    }
	    ptcl_.resizeNoInitialize(i_loc+1);
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

        void dumpMemSizeUsed(std::ostream & fout){
            F64 sum_loc,sum_max;
            sum_loc = (double)(ptcl_.getMemSize()
                              +idx_remove_ptcl_.getMemSize()
                              +ptcl_send_.getMemSize()
                              +ptcl_recv_.getMemSize()) * 1e-9;
            if (Comm::getRank() == 0) {
                fout<<"ptcl_.getMemSize()= "<<ptcl_.getMemSize()<<std::endl;
                fout<<"    ptcl_.size()= "<<ptcl_.size()<<std::endl;
                fout<<"    sizeof(Tptcl)= "<<sizeof(Tptcl)<<std::endl;
                fout<<"idx_remove_ptcl_.getMemSize()= "<<idx_remove_ptcl_.getMemSize()<<std::endl;
                fout<<"ptcl_send_.getMemSize()= "<<ptcl_send_.getMemSize()<<std::endl;
                fout<<"ptcl_recv_.getMemSize()= "<<ptcl_recv_.getMemSize()<<std::endl;
                fout<<"sum[GB]= " << sum_loc << std::endl;
            }
            //MPI::COMM_WORLD.Allreduce(&sum_loc,&sum,1,MPI::DOUBLE,MPI::MAX);
            sum_max = Comm::getMaxValue(sum_loc);
            if (Comm::getRank() == 0) {
               fout << "sum[GB](psys,max.) = " << sum_max << std::endl;
            }
        }
        
        void freeCommunicationBuffer(){
            ptcl_send_.freeMem();
            ptcl_recv_.freeMem();
        }
        
    };
}
