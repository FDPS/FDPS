#pragma once

#include<iostream>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include<mpi.h>
#endif

namespace  ParticleSimulator{

    template<S32 DIM>
    void SetNumberOfDomainMultiDimension(S32 np[], S32 rank[]){
        for(S32 i=0; i<DIMENSION_LIMIT; i++){
            np[i] = 1;
            rank[i] = 1;
        }
        std::vector<S32> npv;
        npv.resize(DIM);
        S32 np_tmp = Comm::getNumberOfProc();
        for(S32 d=DIM, cid=0; cid<DIM-1; d--, cid++){
            S32 tmp = (S32)pow(np_tmp+0.000001, (1.0/d)*1.000001 );
            while(np_tmp%tmp){
                tmp--;
            }
            npv[cid] = tmp;
            np_tmp /= npv[cid];
        }
        npv[DIM-1] = np_tmp;
        S32 rank_tmp = Comm::getRank();
        std::sort(npv.begin(), npv.end(), std::greater<S32>());
        for(S32 i=DIM-1; i>=0; i--){
            np[i] = npv[i];
            rank[i] = rank_tmp % np[i];
            rank_tmp /= np[i];
        }
    }

    class DomainInfo{
    private:

        F64vec * pos_sample_tot_;
        F64vec * pos_sample_loc_;
        //ReallocatableArray<F64vec> pos_sample_tot_;
        //ReallocatableArray<F64vec> pos_sample_loc_;
	
        F64ort * pos_domain_;
        F64ort * pos_domain_temp_;

        F32 coef_ema_;
        S32 target_number_of_sample_particle_;
        S32 number_of_sample_particle_tot_;
        S32 number_of_sample_particle_loc_;
        S32 n_domain_[DIMENSION_LIMIT]; // in 2-dim, n_domain_[2] is always 1.

        F64ort pos_root_domain_;
                
        bool first_call_by_initialize;
        bool first_call_by_decomposeDomain;

        S32 boundary_condition_;
        bool periodic_axis_[DIMENSION_LIMIT]; // in 2-dim, periodic_axis_[2] is always false.

        void sortCoordinateOfSampleParticle(F64vec pos[],
                                            S32 lo,
                                            S32 up,
                                            const S32 &cid) {

            S32 i, j;
            F64vec tpos;

            while (up > lo) {
                i = lo;
                j = up;
                tpos = pos[lo];
                /*** Split file in two ***/
                while (i < j) {
                    for (; pos[j][cid] > tpos[cid]; j--);
                    for (pos[i] = pos[j]; i < j && pos[i][cid] <= tpos[cid]; i++);
                    pos[j] = pos[i];
                }
                pos[i] = tpos;
                /*** Sort recursively, the smallest first ***/
                if (i - lo < up - i) {
                    sortCoordinateOfSampleParticle(pos, lo, i-1, cid);
                    lo = i + 1;
                } else {
                    sortCoordinateOfSampleParticle(pos, i+1, up, cid); 
                    up = i - 1;
                }
            }
        }

        void calculateBoundaryOfDomain(const S32 &np,
                                       const F64vec pos_sample[],
                                       const S32 cid,
                                       const S32 &istart,
                                       const S32 &iend,
                                       F64 & xlow,
                                       F64 & xhigh) {
            if(istart == 0) {
                xlow  = pos_root_domain_.low_[cid];
            } else {
                xlow  = 0.5 * (pos_sample[istart-1][cid] + pos_sample[istart][cid]);
            }
            if(iend == np - 1) {
                xhigh = pos_root_domain_.high_[cid];
            } else {
                xhigh = 0.5 * (pos_sample[iend][cid] + pos_sample[iend+1][cid]);
            }
        }

    public:
        DomainInfo() {
            first_call_by_initialize = true;
            first_call_by_decomposeDomain = true;
            for(S32 k = 0; k < DIMENSION; k++) {
                pos_root_domain_.low_[k]  = - std::numeric_limits<float>::max()*0.0625;
                pos_root_domain_.high_[k] = + std::numeric_limits<float>::max()*0.0625;
                periodic_axis_[k] = false;
            }
            boundary_condition_ = BOUNDARY_CONDITION_OPEN;
        }

        void initialize(const F32 coef_ema = 1.0){
            assert(first_call_by_initialize);
            first_call_by_initialize = false;
            pos_sample_tot_ = NULL;
            pos_sample_loc_ = NULL;
	    

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_      = new F64ort[Comm::getNumberOfProc()];
            pos_domain_temp_ = new F64ort[Comm::getNumberOfProc()];
#else
            pos_domain_     = new F64ort[1];
            pos_domain_temp_= new F64ort[1];
#endif

            coef_ema_ = coef_ema;
            target_number_of_sample_particle_ = 0;
            number_of_sample_particle_tot_ = 0;
            number_of_sample_particle_loc_ = 0;

            S32 rank_tmp[DIMENSION];
            SetNumberOfDomainMultiDimension<DIMENSION>(n_domain_, rank_tmp);
        }

        void setNumberOfDomainMultiDimension(const S32 nx, const S32 ny, const S32 nz=1){
            S32 n_proc = Comm::getNumberOfProc();
            if(n_proc != nx*ny*nz){
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
        }

        void setDomain(const S32 nx, const S32 ny, const S32 nz=1){
            setNumberOfDomainMultiDimension(nx, ny, nz);
        }

/*
        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys,
                                   const F32 weight,
                                   const bool clear = true) {
*/
        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys,
                                   const bool clear,
                                   const F32 weight) {
            if(psys.getFirstCallByDomainInfoCollectSampleParticle()) {
                F64vec *temp_loc = new F64vec[target_number_of_sample_particle_];
                for(S32 i = 0; i < number_of_sample_particle_loc_; i++)
                    temp_loc[i] = pos_sample_loc_[i];                
                
                target_number_of_sample_particle_ += psys.getTargetNumberOfSampleParticle();
                delete [] pos_sample_tot_;
                delete [] pos_sample_loc_;
                
                pos_sample_tot_ = new F64vec[target_number_of_sample_particle_];
                pos_sample_loc_ = new F64vec[target_number_of_sample_particle_];
                for(S32 i = 0; i < number_of_sample_particle_loc_; i++)
                    pos_sample_loc_[i] = temp_loc[i];                                
                delete [] temp_loc;
            }
            
            if(clear) {
                number_of_sample_particle_loc_ = 0;
            }
            
            S32 number_of_sample_particle = 0;

            psys.getSampleParticle(number_of_sample_particle, &pos_sample_loc_[number_of_sample_particle_loc_], weight);
            number_of_sample_particle_loc_ += number_of_sample_particle;
            return;
        }

        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys,
                                   const bool clear){
            const F32 wgh = psys.getNumberOfParticleLocal();
            collectSampleParticle(psys, clear, wgh);
        }

        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys){
            const F32 wgh = psys.getNumberOfParticleLocal();
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
        }

        void decomposeDomain() {
            // ****** collect sample particles to process 0. ****** 
            // ****** Here, Gatherv could be used *****************
            //S32 myrank = MPI::COMM_WORLD.Get_rank();
            //S32 nproc  = MPI::COMM_WORLD.Get_size();
            S32 nproc  = Comm::getNumberOfProc();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            S32 myrank = Comm::getRank();
            if(myrank != 0) {
                MPI::COMM_WORLD.Send(&number_of_sample_particle_loc_, 1, GetDataType<S32>(), 0, myrank*2);
                //MPI::COMM_WORLD.Send((F64 *)&pos_sample_loc_[0][0], number_of_sample_particle_loc_*DIMENSION, GetDataType<F64>(), 0, myrank*2+1);
                MPI::COMM_WORLD.Send(pos_sample_loc_, number_of_sample_particle_loc_, GetDataType<F64vec>(), 0, myrank*2+1);
            } else {
                number_of_sample_particle_tot_ = number_of_sample_particle_loc_;
                for(S32 i = 0; i < number_of_sample_particle_loc_; i++) {                    
                    pos_sample_tot_[i] = pos_sample_loc_[i];
                }                
                for(S32 i = 1; i < nproc; i++) {
                    S32 nreceive;
                    MPI::COMM_WORLD.Recv(&nreceive, 1, GetDataType<S32>(), i, i*2);
                    
                    //MPI::COMM_WORLD.Recv((F64 *)(&pos_sample_tot_[0][0]+number_of_sample_particle_tot_*DIMENSION), nreceive*DIMENSION, GetDataType<F64>(), i, i*2+1);
                    MPI::COMM_WORLD.Recv(pos_sample_tot_+number_of_sample_particle_tot_, nreceive, GetDataType<F64vec>(), i, i*2+1);
                    number_of_sample_particle_tot_ += nreceive;
                }
            }
            // ****************************************************
            // *** decompose domain *******************************
            if(myrank == 0) {
                S32 * istart = new S32[nproc];
                S32 * iend   = new S32[nproc];
                // --- x direction --------------------------
                sortCoordinateOfSampleParticle(pos_sample_tot_, 0, number_of_sample_particle_tot_-1, 0);
                for(S32 i = 0; i < nproc; i++) {
                    istart[i] = (i * number_of_sample_particle_tot_) / nproc;
                    if(i > 0)
                        iend[i-1] = istart[i] - 1;
                }
                iend[nproc-1] = number_of_sample_particle_tot_ - 1;
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    
                    F64 x0, x1;
                    
                    calculateBoundaryOfDomain(number_of_sample_particle_tot_, pos_sample_tot_, 0, istart[ix0], iend[ix1-1], x0, x1);
                    for(S32 i = ix0; i < ix1; i++) {
                        pos_domain_temp_[i].low_[0]  = x0;
                        pos_domain_temp_[i].high_[0] = x1;
                    }
                }
                // ------------------------------------------
                // --- y direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    sortCoordinateOfSampleParticle(pos_sample_tot_, istart[ix0], iend[ix1-1], 1);
                    S32 number_of_sample_particle_tot_y = iend[ix1-1] - istart[ix0] + 1;
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        F64 y0, y1;
                        
                        calculateBoundaryOfDomain(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0], 1, istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
                        for(S32 i = iy0; i < iy1; i++) {
                            pos_domain_temp_[i].low_[1]  = y0;
                            pos_domain_temp_[i].high_[1] = y1;
                        }
                    }
                }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- z direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        sortCoordinateOfSampleParticle(pos_sample_tot_, istart[iy0], iend[iy1-1], 2);
                        S32 number_of_sample_particle_tot_z = iend[iy1-1] - istart[iy0] + 1;
                        for(S32 iz = 0; iz < n_domain_[2]; iz++) {
                            S32 iz0 = iy0 + iz;
                            F64 z0, z1;
                            
                            calculateBoundaryOfDomain(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], 2, istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
                            pos_domain_temp_[iz0].low_[2]  = z0;
                            pos_domain_temp_[iz0].high_[2] = z1;
                        }
                    }
                }
#endif // PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- process first ------------------------
                if(first_call_by_decomposeDomain) {
                    first_call_by_decomposeDomain = false;
                    for(S32 i = 0; i < nproc; i++) {
                        pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                        pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                    }
                } else {
                    for(S32 i = 0; i < nproc; i++) {
                        pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                        pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                        
                    }
                }
                // ------------------------------------------
                delete [] istart;
                delete [] iend;
            }
            // ****************************************************
            // *** broad cast pos_domain_ *************************
            MPI::COMM_WORLD.Bcast(pos_domain_, nproc, GetDataType<F64ort>(), 0);
            //Comm::broadcast(pos_domain_, nproc);
            // ****************************************************
#else       // PARTICLE_SIMULATOR_MPI_PARALLEL
            
            pos_domain_[0] = pos_root_domain_;
            
#endif     // PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"pos_root_domain_="<<pos_root_domain_<<std::endl;
            std::cout<<"pos_domain_[Comm::getRank()]="<<pos_domain_[Comm::getRank()]<<std::endl;
            if(Comm::getRank()==0){
                for(S32 i=0; i<Comm::getNumberOfProc(); i++){
                    std::cout<<"pos_domain_["<<i<<"]="<<pos_domain_[i]<<std::endl;
                }
            }
#endif
        }

        template<class Tpsys>
        void decomposeDomainAll(Tpsys & psys,
                                const F32 wgh){
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
            decomposeDomain();
        }

        template<class Tpsys>
        void decomposeDomainAll(Tpsys & psys){
            const F32 wgh = psys.getNumberOfParticleLocal();
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
            decomposeDomain();
        }

        void getRootDomain(FILE *fp) {
            fprintf(fp, "%+e %+e %+e\n",
                pos_root_domain_.low_[0],
                pos_root_domain_.low_[1],
                pos_root_domain_.low_[2]);
            fprintf(fp, "%+e %+e %+e\n",
                pos_root_domain_.high_[0],
                pos_root_domain_.high_[1],
                pos_root_domain_.high_[2]);

            return;
        }


#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        void getSampleParticleLocal(FILE *fp) {
            for(S32 i = 0; i < number_of_sample_particle_loc_; i++) {
                fprintf(fp, "%+e %+e %+e\n", pos_sample_loc_[i][0], pos_sample_loc_[i][1], pos_sample_loc_[i][2]);
            }            
            return;
        }

        void getSampleParticleTotal(FILE *fp) {
            for(S32 i = 0; i < number_of_sample_particle_tot_; i++) {
                fprintf(fp, "%+e %+e %+e\n", pos_sample_tot_[i][0], pos_sample_tot_[i][1], pos_sample_tot_[i][2]);
            }            
            return;
        }
#endif
	
        void getPosDomainTotal(FILE *fp) {
            //S32 nproc = MPI::COMM_WORLD.Get_size();
            S32 nproc = Comm::getNumberOfProc();
            for(S32 i = 0; i < nproc; i++) {
                for(S32 k = 0; k < DIMENSION; k++)
                    fprintf(fp, "%+e ", pos_domain_[i].low_[k]);
                for(S32 k = 0; k < DIMENSION; k++)
                    fprintf(fp, "%+e ", pos_domain_[i].high_[k]);
                fprintf(fp, "\n");
            }
        }

        S32 * getPointerOfNDomain(){return n_domain_;};

        /* AT_DEBUG
        F32ort * getPointerOfPosDomain(){return pos_domain_;};
        */
        F64ort * getPointerOfPosDomain(){return pos_domain_;};

        // A. Tanikawa need this method for Particle Mesh...
        S32 getNDomain(const S32 dim) const {return n_domain_[dim];};

        /* AT_DEBUG
        // for DEBUG
        F32vec & getPosSample(const S32 id=0){return pos_sample_tot_[id];}
        // for DEBUG // A. Tanikawa would not like to delete this method...
        F32ort & getPosDomain(const S32 id=0) const {return pos_domain_[id];}
        // for DEBUG
        void setPosDomain(const S32 id, const F32ort & pos){ pos_domain_[id] = pos;}
        */
        // for DEBUG
        F64vec & getPosSample(const S32 id=0){return pos_sample_tot_[id];}
        // for DEBUG // A. Tanikawa would not like to delete this method...
        F64ort & getPosDomain(const S32 id=0) const {return pos_domain_[id];}
        // for DEBUG
        void setPosDomain(const S32 id, const F64ort & pos){ pos_domain_[id] = pos;}

        void setBoundaryCondition(enum BOUNDARY_CONDITION bc){
            boundary_condition_ = bc;
            if(DIMENSION == 2 && 
               (bc == BOUNDARY_CONDITION_PERIODIC_XYZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_XZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_YZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_Z ) ){
                throw "PS_ERROR: in setBoundaryCondition(enum BOUNDARY_CONDITION) \n boundary condition is incompatible with DIMENSION";
            }
            if(bc == BOUNDARY_CONDITION_PERIODIC_X) periodic_axis_[0] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_Y) periodic_axis_[1] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_Z) periodic_axis_[2] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XY) periodic_axis_[0] = periodic_axis_[1] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XZ) periodic_axis_[0] = periodic_axis_[2] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_YZ) periodic_axis_[1] = periodic_axis_[2] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XYZ) periodic_axis_[0] = periodic_axis_[1] = periodic_axis_[2] = true;
        }

        S32 getBoundaryCondition() const { return boundary_condition_; }

        /* AT_DEBUG
        void setPosRootDomain(const F32vec & low, const F32vec & high){
        */
        void setPosRootDomain(const F64vec & low, const F64vec & high){
            for(S32 i=0; i<DIMENSION; i++){
                //std::cerr<<"low[i]="<<low[i]<<std::endl;
                //std::cerr<<"high[i]="<<high[i]<<std::endl;
                //std::cerr<<"periodic_axis_[i]="<<periodic_axis_[i]<<std::endl;
                if( periodic_axis_[i] == false ) continue;
                if(low[i] < high[i]){
                    pos_root_domain_.low_[i] = low[i];
                    pos_root_domain_.high_[i] = high[i];
                }
            }
        }

        /* AT_DEBUG
        F32ort getPosRootDomain() const { return pos_root_domain_; }
        */
        const F64ort getPosRootDomain() const { return pos_root_domain_; }

        void getPeriodicAxis(bool pa[]) const {
            for(S32 i=0; i<DIMENSION; i++) pa[i] = periodic_axis_[i];
        }

        template<class Tpsys>
        bool checkCollectSampleParticleSubset(Tpsys & psys);
        template<class Tpsys>
        bool checkCollectSampleParticleAverage(Tpsys & psys);
        template<class Tpsys>
        bool checkDecomposeDomain(Tpsys & psys);

    };
}
