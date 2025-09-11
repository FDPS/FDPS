#pragma once
#include "../ps_defs.hpp"
#include "multidimensional_array.hpp"
#include "fmm.hpp"
#include "cell.hpp"

namespace ParticleSimulator {

    static void ForceKernelOpen(Cell &ci, const Cell &cj){
        const int ni = ci.plist.size();
        const int nj = cj.plist.size();
        for(int i=0; i<ni; i++){
            Particle &pi = *ci.plist[i];
            for(int j=0; j<nj; j++){
                const Particle &pj = *cj.plist[j];
                if(&pi == &pj) continue;
    
                const F64vec  dr  = pj.pos - pi.pos;
                const F64     r2  = dr*dr;
                const F64     ri2 = 1.0 / r2;
                const F64     ri  = sqrt(ri2);
                const F64     ri3 = ri * ri2;
                pi.phi_direct += pj.mass * ri;
                pi.acc_direct += (pj.mass * ri3) * dr;
            }
        }
    }

    static void ForceKernelPeriodicXYZ(Cell &ci, const Cell &cj){
        const int ni = ci.plist.size();
        const int nj = cj.plist.size();
        for(int i=0; i<ni; i++){
            Particle &pi = *ci.plist[i];
            for(int j=0; j<nj; j++){
                const Particle &pj = *cj.plist[j];
                if(&pi == &pj) continue;
    
                const F64vec dr  = minimum_image(pj.pos - pi.pos);
                const F64    r2  = dr*dr;
                const F64    ri2 = 1.0 / r2;
                const F64    ri  = sqrt(ri2);
                const F64    ri3 = ri * ri2;
                pi.phi_direct += pj.mass * ri;
                pi.acc_direct += (pj.mass * ri3) * dr;
                // puts("evaluated PP");
            }
        }
    }
    
    template <int PFMM, int ICUT, int DIM>
    static void calcForceParticleParticleOpenImpl(const int NX,
                                                  const int NY,
                                                  const int NZ,
                                                  MultidimensionalArray< Cell_FMM<PFMM>, DIM> & cell)
    {
        int i, j, k;
        for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
        {
            int ii, jj, kk;
            for(kk=k-ICUT; kk<=k+ICUT; kk++) for(jj=j-ICUT; jj<=j+ICUT; jj++) for(ii=i-ICUT; ii<=i+ICUT; ii++)
            {
                if(kk < 0 || kk >= NZ) continue;
                if(jj < 0 || jj >= NY) continue;
                if(ii < 0 || ii >= NX) continue;
                ForceKernelOpen(cell(k,j,i), cell(kk,jj,ii));
            }
        }
    }

    template <int PFMM, int ICUT, int DIM>
    static void calcForceParticleParticlePeriodicXYZImpl(const int NX, 
                                                         const int NY,
                                                         const int NZ,
                                                         MultidimensionalArray< Cell_FMM<PFMM>, DIM> & cell)
    {
        int i, j, k;
        for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
        {
            int ii, jj, kk;
            for(kk=k-ICUT; kk<=k+ICUT; kk++) for(jj=j-ICUT; jj<=j+ICUT; jj++) for(ii=i-ICUT; ii<=i+ICUT; ii++)
            {
                const int iii = (ii+NX)%NX;
                const int jjj = (jj+NY)%NY;
                const int kkk = (kk+NZ)%NZ;
                ForceKernelPeriodicXYZ(cell(k,j,i), cell(kkk,jjj,iii));
            }
        }
    }

    template <int PFMM, int ICUT, int DIM>
    void calcForceParticleParticle(const int NX,
                                   const int NY,
                                   const int NZ,
                                   MultidimensionalArray< Cell_FMM<PFMM>, DIM> & cell,
                                   const int bc) {
        if (bc == BOUNDARY_CONDITION_OPEN) {
            calcForceParticleParticleOpenImpl<PFMM, ICUT, DIM>(NX, NY, NZ, cell);
        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
            calcForceParticleParticlePeriodicXYZImpl<PFMM, ICUT, DIM>(NX, NY, NZ, cell);
        }
    }

    
} // END of namespace of ParticleSimulator
