/*
 Derived from PMMM (Apache License 2.0)
 Original Copyright (c) 2012-2014 Keigo Nitadori
 Modifications (c) 2025 FDPS team
 Licensed under the Apache License, Version 2.0 (see LICENSE)
*/
#pragma once
#include <cmath>
#include <complex>
#include <vector>

namespace ParticleSimulator {

    template <class real_t=double,
              class cplx_t=std::complex<real_t>>
    class CutFunc {
    // This class is only used in the Particle Mesh Multipole module.
    public:
        int LMAX;
        std::vector<real_t> gammainv;
        std::vector<real_t> rcut;
        std::vector<cplx_t> kcut;
    
        CutFunc() {}
    
        void init(int _LMAX){
            LMAX = _LMAX;
            gammainv.reserve(1+LMAX);
            rcut.reserve(1+LMAX);
            kcut.reserve(1+LMAX);
            double val = sqrt(4.0*atan(1.0));
            for(int ell=0; ell<=LMAX; ell++){
                gammainv[ell] = 1.0 / val;
                val *= (0.5 + ell);
            }
        }
    
        void eval_rcut(const double x){
            const double xexp = exp(-x);
            double xpow = sqrt(x);
            double val  = sqrt(4.0*atan(1.0)) * erfc(xpow);
            for(int ell=0; ell<=LMAX; ell++){
                rcut[ell] = val * gammainv[ell];
                val = (0.5 + ell)*val + xpow * xexp;
                xpow *= x;
            }
        }
    
        void eval_kcut(const double kn2, const double kappa){
            const double pi = 4.0 * atan(1.0);
            const double gauss = exp(-kn2 * ((pi/kappa)*(pi/kappa)));
    
            const cplx_t ratio(0.0, -pi * kn2);
            cplx_t coef = 1.0 / sqrt(pi * kn2);
            for(int ell=0; ell<=LMAX; ell++){
                kcut[ell] = coef * (gauss * gammainv[ell]);
                coef *= ratio;
            }
        }

        void freeMem() {
            // This function explicitly free memories.
            std::vector<real_t>().swap(gammainv);
            std::vector<real_t>().swap(rcut);
            std::vector<cplx_t>().swap(kcut);
        }
    };

} // END of namespace of ParticleSimulator
